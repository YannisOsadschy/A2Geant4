/***** Electron drift model for TPC *****
 * Implementation of TPC physics simulated in Heed/Garfield/Degrad.
 * Based off method described in https://arxiv.org/pdf/1806.05880.pdf and
 * https://github.com/lennertdekeukeleere/Geant4GarfieldDegradInterface/tree/master/ALICE.
 * Uses G4VFastSimulation to implement manual definition of electron drift.
 * Gaussians based on work by Fabian Metzger, presentation shared by PNPI.
 ***** AC Postuma 2021 *****/

#include "A2DriftModel.hh"

#include "A2DetectorConstruction.hh" //detector construction
#include "A2DriftandHitLogic.hh"
#include "A2UserTrackInformation.hh"

#include "G4VPhysicalVolume.hh" //volumes
#include "G4Electron.hh" //particle the model applies to
#include "G4SystemOfUnits.hh" //units
#include "G4RunManager.hh" //run
#include "G4TrackingManager.hh" //tracking of particles
#include "G4EventManager.hh" //events
#include "G4VVisManager.hh" //track visualization

#include <iostream>

using namespace CLHEP;

/**** Constructor *****/
A2DriftModel::A2DriftModel(G4String modelName, G4Region* actVol)
: G4VFastSimulationModel(modelName, actVol), fDrifter(actVol){ //fast simulation implements user defined physics response   
}

/***** Destructor *****/
A2DriftModel::~A2DriftModel(){

}

/***** Called in SteppingAction: checks particle type if model is applicable ******/
G4bool A2DriftModel::IsApplicable(const G4ParticleDefinition& particleType){
	G4String particleName = particleType.GetParticleName();
	if(particleName=="e-")return true; //only applicable to electrons
	return false;
}

/***** Called in SteppingAction: conditions in which to trigger model *****/
G4bool A2DriftModel::ModelTrigger(const G4FastTrack& fastTrack){
	G4double ekin = fastTrack.GetPrimaryTrack()->GetKineticEnergy()/keV;
	if (ekin <=1)return true; //trigger for kinetic energy below 1 keV
	return false;
}

/***** This function contains the main operation of the model *****/
void A2DriftModel::DoIt(const G4FastTrack& fastTrack, G4FastStep& fastStep){
	//get all relevant step data from the track
	G4ThreeVector direction = fastTrack.GetPrimaryTrack()->GetMomentumDirection();
	G4ThreeVector worldPosition = fastTrack.GetPrimaryTrack()->GetPosition()/mm;
	G4double ekin = fastTrack.GetPrimaryTrack()->GetKineticEnergy()/keV;
	G4double time = fastTrack.GetPrimaryTrack()->GetGlobalTime();
	G4String particleName = fastTrack.GetPrimaryTrack()->GetParticleDefinition()->GetParticleName();
	//make sure the electron is in a position to be transported
	G4bool position = true; //it normally will be
	if (worldPosition.z() < -115.5) position = false; //electron already past anode
	G4double radius = sqrt(worldPosition.x()*worldPosition.x()+worldPosition.y()*worldPosition.y());
	if ((worldPosition.z() > 115.5)&&(radius<5)) position=false; //electron behind cathode
	if (position == true){
		Transport(fastStep, fastTrack, particleName, ekin, time, worldPosition.x(), worldPosition.y(), 
                worldPosition.z());
	}
	}

//Generate electron position, time when reaching anode
void A2DriftModel::Transport(G4FastStep& fastStep,const G4FastTrack& fastTrack, G4String particleName, 
                            double ekin_keV, double t, double x_mm, double y_mm, double z_mm){
	A2DriftandHitLogic::TransportValues transportValues 
    = fDrifter.GetTransportValues(particleName, ekin_keV, t, x_mm, y_mm, z_mm);
    fastStep.SetPrimaryTrackFinalProperTime(transportValues.time);
	fastStep.SetPrimaryTrackPathLength(transportValues.pathLength*mm); //travel calculated distance
	fastStep.SetPrimaryTrackFinalPosition(transportValues.position); //final calculated position
	fastStep.SetTotalEnergyDeposited(transportValues.eKin_keV*keV); //deposit all energy
	/**** kill step and call hit ****/
    fDrifter.ProcessHit(transportValues.position,transportValues.eKin_keV,transportValues.time, 
                    static_cast<A2UserTrackInformation*>((fastTrack.GetPrimaryTrack())->GetUserInformation())->GetTrackID(), 
                    static_cast<A2UserTrackInformation*>((fastTrack.GetPrimaryTrack())->GetUserInformation())->GetPartID(), 
                    fastTrack.GetPrimaryTrack()->GetDynamicParticle()->GetCharge());
	fastStep.KillPrimaryTrack();
}

/***** Reimplement so that class works *****/
void A2DriftModel::ProcessEvent(){
	//reimplement from G4VFastSimulation
}
void A2DriftModel::Reset(){
	//reimplement from G4VFastSimulation
}
