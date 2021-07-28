/***** Electron drift model for TPC *****
 * Implementation of TPC physics simulated in Heed/Garfield/Degrad.
 * Based off method described in https://arxiv.org/pdf/1806.05880.pdf and
 * https://github.com/lennertdekeukeleere/Geant4GarfieldDegradInterface/tree/master/ALICE.
 * Uses G4VFastSimulation to implement manual definition of electron drift.
 ***** AC Postuma 2021 *****/

#include "A2HeedModel.hh"

#include "G4VPhysicalVolume.hh" //volumes
#include "G4Electron.hh" //particle the model applies to
//#include "G4Gamma.hh" //another relevant (?) particle
#include "G4SystemOfUnits.hh" //units
#include "A2DetectorConstruction.hh" //detector construction
#include "G4RunManager.hh" //run
//#include <stdio.h>
#include "A2Target.hh" //targets
#include "A2SD.hh" //sensitive detector
//#include "DriftLineTrajectory.hh"
#include "G4TrackingManager.hh" //tracking of particles
#include "G4EventManager.hh" //events
#include "G4TransportationManager.hh" //particle transport
#include "G4VVisManager.hh" //track visualization
#include "CLHEP/Random/RandGauss.h" //random generation

using namespace CLHEP;

/**** Constructor *****/
A2HeedModel::A2HeedModel(G4String modelName, G4Region* actVol, A2Target* target, A2SD* anode)
: G4VFastSimulationModel(modelName, actVol), fA2Target(target), fA2SD(anode) { //fast simulation implements user defined physics response
	//initiate pointers
	fFakeStep = new G4Step(); //step used to call hit in SD
	fFakePreStepPoint  = fFakeStep->GetPreStepPoint(); //step point
  	fFakePostStepPoint = fFakeStep->GetPostStepPoint(); //step point
  	fTouchableHandle   = new G4TouchableHistory(); //touchable for step
  	fpNavigator        = new G4Navigator(); //navigator to find SD
	fNaviSetup = false; //if setup has already been done
}

/***** Destructor *****/
A2HeedModel::~A2HeedModel(){
	//remove objects requiring manual deletion
	delete fFakeStep;
  	delete fpNavigator;
}

/***** Called in SteppingAction: checks particle type if model is applicable ******/
G4bool A2HeedModel::IsApplicable(const G4ParticleDefinition& particleType){
	G4String particleName = particleType.GetParticleName();
	if(particleName=="e-")return true; //only applicable to electrons
	return false;
}

/***** Called in SteppingAction: conditions in which to trigger model *****/
G4bool A2HeedModel::ModelTrigger(const G4FastTrack& fastTrack){
	G4double ekin = fastTrack.GetPrimaryTrack()->GetKineticEnergy()/keV;
	if (ekin <=1)return true; //trigger for kinetic energy below 1 keV
	return false;
}

/***** This function contains the main operation of the model *****/
void A2HeedModel::DoIt(const G4FastTrack& fastTrack, G4FastStep& fastStep){
	//get all relevant step data from the track
	G4ThreeVector direction = fastTrack.GetPrimaryTrack()->GetMomentumDirection();
	G4ThreeVector worldPosition = fastTrack.GetPrimaryTrack()->GetPosition()/mm;
	G4double ekin = fastTrack.GetPrimaryTrack()->GetKineticEnergy()/keV;
	G4double time = fastTrack.GetPrimaryTrack()->GetGlobalTime();
	G4String particleName = fastTrack.GetPrimaryTrack()->GetParticleDefinition()->GetParticleName();
	//pass step data to transportation 
	Transport(fastStep, fastTrack, particleName, ekin, time, worldPosition.x(), worldPosition.y(), worldPosition.z(), direction.x(), direction.y(), direction.z());
}

//take just the delta electron run instructions
void A2HeedModel::Transport(G4FastStep& fastStep,const G4FastTrack& fastTrack, G4String particleName, double ekin_keV, double t, double x_mm, double y_mm, double z_mm, double dx, double dy, double dz){
	//G4cout<<"Transporting delta electron of energy "<< ekin_keV <<" keV"<<G4endl; //debugging message
	/****transport each electron to anode ****/
	G4double z_pos = -115.5; //set final z position to anode z position
	G4double pathLength = z_pos - z_mm; //total length to final z position in mm

	//use Gaussians as defined by Fabian Metzger in his TPC work
	//set gas parameters as calculated from Garfield
	//currently set only for helium-4 at 30 bar and 2kV/cm
	G4double drift_vel = 2633; //drift velocity in mm/ms
	G4double trans_diff = 0.0219; //transverse diffusion mm^1/2
	G4double long_diff = 0.0353; //longitudinal diffusion mm^1/2

	//set means and sigmas of Gaussian distribtuions
	G4double mean_x=x_mm; //mean for Gaussian calc of x pos
	G4double mean_y=y_mm; //mean for Gaussian calc of y pos
	G4double mean_t = abs(pathLength/drift_vel); //mean for Gaussian calc of time: comes out in ms
	G4double sigma_diff = trans_diff*sqrt(abs(pathLength)); //sigma for Gaussian calc of x,y positions: close enough to mm
	G4double sigma_time = long_diff/drift_vel*sqrt(abs(pathLength)); //comes out close enough to ms

	//use random number generation to get values for positions and times
	RanluxEngine aRandEngine; //make an engine for random number generation
	RandGauss gaussian(aRandEngine);
	G4double x_pos = gaussian.shoot(mean_x,sigma_diff); //calculate an x position: mm
	G4double y_pos = gaussian.shoot(mean_y,sigma_diff); //calc y mm
	G4double time = gaussian.shoot(mean_t,sigma_time); //calc a time ms

	//combine position data into a vector
	G4ThreeVector position = G4ThreeVector(x_pos*mm,y_pos*mm,z_pos*mm);
	
	/**** set final track parameters *****/
	fastStep.SetPrimaryTrackFinalProperTime(time);
	fastStep.SetPrimaryTrackPathLength(pathLength*mm); //travel calculated distance
	fastStep.SetPrimaryTrackFinalPosition(position); //final calculated position
	fastStep.SetTotalEnergyDeposited(ekin_keV*keV); //deposit all energy
	
	/**** kill step and call hit ****/
	ProcessHit(fastTrack,position,ekin_keV,time);
	fastStep.KillPrimaryTrack();
	
	/**** troubleshooting ****/
	//G4cout<<position;
	//G4double radius = sqrt(x_pos*x_pos + y_pos+y_pos); //calculate radius - will it hit the anode?
	//if(radius <=50)G4cout<<"Anode Hit!"<<G4endl;//if it hits the anode radius
	//else G4cout<<"Electron misses anode"<<G4endl;
	//	G4cout<<"Anode Hit!"<<G4endl;
	//	fSensitive->Hit(fastStep);

	//} else { 
	//	G4cout<<"Electron misses anode"<<G4endl;
	//}
	//fastStep.DumpInfo(); //print all step information
}

/***** Process secondary electrons created along drift track *****/
void A2HeedModel::ProcessSecondaries(){
	//do nothing for now
	//later deal with secondary electrons created
}


/***** Call a hit in the anode for each electron that reaches it *****/
void A2HeedModel::ProcessHit(const G4FastTrack& fastTrack,G4ThreeVector position, G4double ekin_keV, G4double drift_time){
/**** set up touchable in current volume ****/
	if (!fNaviSetup) {
		fpNavigator->SetWorldVolume(G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking()->GetWorldVolume());
		fpNavigator->LocateGlobalPointAndUpdateTouchableHandle(position,G4ThreeVector(0.,0.,0.),fTouchableHandle,true);
		fNaviSetup = true;
    	} else {
      		fpNavigator->
        	LocateGlobalPointAndUpdateTouchableHandle(position,G4ThreeVector(0.,0.,0.),fTouchableHandle);
     	}

  	/**** fill G4Step with information necessary for the sensitive detector ****/
  	//set track to step
  	fFakeStep->SetTrack(const_cast<G4Track*>(fastTrack.GetPrimaryTrack()));
	//set touchable for step position
  	fFakePreStepPoint->SetTouchableHandle(fTouchableHandle);
	fFakePreStepPoint->SetPosition(position); //try to get it to work properly
  	//set total energy deposit
  	fFakeStep->SetTotalEnergyDeposit(ekin_keV);
  	//set the time of hit
  	fFakeStep->GetPreStepPoint()->SetGlobalTime(drift_time*ms); //time hist will not print due to large value - why???
  	//debugging statement: print relevant parameters here
  	//G4cout<<fFakeStep->GetPreStepPoint()->GetTouchable()->GetVolume()->GetCopyNo()<<G4endl;
	//if (position.x()*position.x()+position.y()*position.y() > 2500) G4cout<<fFakeStep->GetPreStepPoint()->GetPosition()/mm<<G4endl;

	/**** call hit in sensitive detector ****/
  	G4VPhysicalVolume* fCurrentVolume = fFakeStep->GetPreStepPoint()->GetPhysicalVolume();
	G4VSensitiveDetector* fSensitive;
	if( fCurrentVolume != 0 ) {
	fSensitive = fCurrentVolume->GetLogicalVolume()->GetSensitiveDetector();
	if( fSensitive != 0 ) {
		fSensitive->Hit(fFakeStep);
		}
	}
}

/***** Create proper response in detector *****/
void A2HeedModel::GenerateDetectorResponse(){
	//do nothing for now
	//later simulate a proper anode signal
}

/***** Reimplement so that class works *****/
void A2HeedModel::ProcessEvent(){
	//reimplement from G4VFastSimulation
}
void A2HeedModel::Reset(){
	//reimplement from G4VFastSimulation
}

