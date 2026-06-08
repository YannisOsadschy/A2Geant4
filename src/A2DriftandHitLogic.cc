#include "G4Region.hh"
#include "G4VPhysicalVolume.hh" //volumes
#include "G4VSensitiveDetector.hh"
#include "G4TransportationManager.hh" //particle transport

#include "TSpline.h"

#include "A2DriftandHitLogic.hh"
#include "A2UserTrackInformation.hh"

#include <math.h>
using namespace CLHEP;

/**** Constructor *****/
A2DriftandHitLogic::A2DriftandHitLogic(G4Region* actVol)
{ //fast simulation implements user defined physics response
	//initiate pointers
	fFakeStep = new G4Step(); //step used to call hit in SD
    fFakeTrack = new G4Track(); //dummy Track just to hold A2UserTrackinformation, which is required for A2SD to work properly
	fFakePreStepPoint  = fFakeStep->GetPreStepPoint(); //step point
  	fFakePostStepPoint = fFakeStep->GetPostStepPoint(); //step point
  	fTouchableHandle   = new G4TouchableHistory(); //touchable for step
  	fpNavigator        = new G4Navigator(); //navigator to find SD
	fNaviSetup = false; //if setup has already been done
	//read in data to set drift velocity, diffusion coefficients
	SetConstants(actVol);
}

A2DriftandHitLogic::~A2DriftandHitLogic()
{
    //remove objects requiring manual deletion
	delete fFakeStep;
  	delete fpNavigator;
}

//Generate electron position, time when reaching anode
A2DriftandHitLogic::TransportValues A2DriftandHitLogic::GetTransportValues(G4String particleName, double ekin_keV, double t, double x_mm,
                                            double y_mm, double z_mm)
{
	//G4cout<<"Transporting delta electron of energy "<< ekin_keV <<" keV"<<G4endl; //debugging message
	/****transport each electron to anode ****/
	G4double z_pos = -115.5; //set final z position to anode z position
	G4double pathLength = z_pos - z_mm; //total length to final z position in mm

	//use Gaussians as defined by Fabian Metzger in his TPC work
	//set means and sigmas of Gaussian distribtuions
	G4double mean_x=x_mm; //mean for Gaussian calc of x pos
	G4double mean_y=y_mm; //mean for Gaussian calc of y pos
	G4double mean_t = abs(pathLength/drift_vel); //mean for Gaussian calc of time: comes out in ms
	G4double sigma_diff = trans_diff*sqrt(abs(pathLength)); //sigma for Gaussian calc of x,y positions: close enough to mm
	G4double sigma_time = long_diff/drift_vel*sqrt(abs(pathLength)); //comes out close enough to ms

	//use random number generation to get values for positions and times
	G4double x_pos = fGaussian.shoot(mean_x,sigma_diff); //calculate an x position: mm
	G4double y_pos = fGaussian.shoot(mean_y,sigma_diff); //calc y mm
	G4double time = fGaussian.shoot(mean_t,sigma_time); //calc a time ms

	//combine position data into a vector
	G4ThreeVector position = G4ThreeVector(x_pos*mm,y_pos*mm,z_pos*mm);

	//let's see what's up here
	//G4cout<<"("<<x_mm<<","<<y_mm<<","<<z_mm<<") to ("<<x_pos<<","<<y_pos<<","<<z_pos<<")"<<G4endl;
	//G4cout<<"r="<<sqrt(x_mm*x_mm+y_mm*y_mm)<<"mm to r="<<sqrt(x_pos*x_pos+y_pos*y_pos)<<"mm in "<<time<<" ms"<<G4endl;

    TransportValues transportValues;
    transportValues.time = time;
    transportValues.pathLength = pathLength;
    transportValues.eKin_keV = ekin_keV;
    transportValues.position = position;
    
    return transportValues;
}

/***** Call a hit in the anode for each electron that reaches it *****/
void A2DriftandHitLogic::ProcessHit(G4ThreeVector position, G4double ekin_keV, G4double drift_time, G4int trackID, G4int partID, G4double charge){
/**** set up touchable in current volume ****/
    if (!fNaviSetup) {
		fpNavigator->SetWorldVolume(G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking()->GetWorldVolume());
		fpNavigator->LocateGlobalPointAndUpdateTouchableHandle(position,G4ThreeVector(0.,0.,0.),fTouchableHandle,true);
		fNaviSetup = true;
    	} else {
      		fpNavigator->
        	LocateGlobalPointAndUpdateTouchableHandle(position,G4ThreeVector(0.,0.,0.),fTouchableHandle);
     	}
    A2UserTrackInformation* anInfo = new A2UserTrackInformation();
    anInfo->SetTrackID(trackID);
    anInfo->SetPartID(partID);
    fFakeTrack->SetUserInformation(anInfo);
    const_cast<G4DynamicParticle*>(fFakeTrack->GetDynamicParticle())->SetCharge(charge);
  	/**** fill G4Step with information necessary for the sensitive detector ****/
  	fFakeStep->SetTrack(fFakeTrack);
	//set touchable for step position
  	fFakePreStepPoint->SetTouchableHandle(fTouchableHandle);
	fFakePreStepPoint->SetPosition(position); //try to get it to work properly
  	//G4cout<<position<<" "<<sqrt(position.x()*position.x()+position.y()*position.y())<<G4endl;
	//set total energy deposit
  	fFakeStep->SetTotalEnergyDeposit(ekin_keV);
  	//set the time of hit
  	fFakeStep->GetPreStepPoint()->SetGlobalTime(drift_time*ms); 

	/**** call hit in sensitive detector ****/
  	G4VPhysicalVolume* fCurrentVolume = fFakeStep->GetPreStepPoint()->GetPhysicalVolume();
	G4VSensitiveDetector* fSensitive;
	if( fCurrentVolume != 0 ) { //if you can find volume
	fSensitive = fCurrentVolume->GetLogicalVolume()->GetSensitiveDetector();
	if( fSensitive != 0 ) { //if volume has sensitive detector
		fSensitive->Hit(fFakeStep); //call hit for the fake step
		//G4cout<<"Calling new hit"<<G4endl;
		}
	}


//interface to manually start the sample drift process hit chain, for edep ion pairs in a TPC
void A2DriftandHitLogic::SampleEdep(const G4Step* aStep)
{   
    
    G4double meanNElec = aStep->GetTotalEnergyDeposit() / fWorkFunction;
    G4double nElec = std::round(fPoisson.shoot(meanNElec));
    
    G4ThreeVector unitDirection = aStep->GetDeltaPosition().unit();
    G4ThreeVector preStepPosition = aStep->GetPreStepPoint()->GetPosition();
    G4double stepLength = aStep->GetStepLength();
    G4double deltaTime = aStep->GetDeltaTime();
    G4double preStepTime = aStep->GetPreStepPoint()->GetGlobalTime();
    G4int trackID = static_cast<A2UserTrackInformation*>(aStep->GetTrack()->GetUserInformation())->GetTrackID();
    G4int parentID = static_cast<A2UserTrackInformation*>(aStep->GetTrack()->GetUserInformation())->GetPartID();
    for (G4int i = 0; i < nElec; ++i)
    {
        G4double randFlat = fFlat.shoot();
        G4ThreeVector positionElec = preStepPosition + randFlat * stepLength * unitDirection;
        G4double timeElec = preStepTime + randFlat * deltaTime;
        G4double eKin_keV = 0.001;
        TransportValues transportValues = GetTransportValues("e-", eKin_keV, timeElec, positionElec.x(),
                                            positionElec.y(), positionElec.z());
        ProcessHit(transportValues.position, transportValues.eKin_keV, transportValues.time, trackID, parentID, -1);
        
    }
}

























/**** Assign gas parameters depending on isotope, pressure of helium ****/
void A2DriftandHitLogic::SetConstants(G4Region *gasRegion){
	G4String name=gasRegion->GetRootLogicalVolumeIterator()[0]->GetMaterial()->GetName();
	G4double pressure = gasRegion->GetRootLogicalVolumeIterator()[0]->GetMaterial()->GetPressure()/bar;
	G4double p_bar[6]={5,10,15,20,25,30}; //supported pressures
	drift_vel=trans_diff=long_diff=0; //initialize
	//set the correct set of constants for helium isotope
	if (name.contains("3")){ //helium-3
		G4double v[6]={7825,5445,4498,3940,3540,3235};
		G4double dl[6]={0.0421,0.0310,0.0259,0.0236,0.0215,0.0202};
		G4double dt[6]={0.0591,0.0447,0.0370,0.0329,0.0294,0.0275};
	for (G4int i=0; i>6; i++){
		if (pressure == p_bar[i]){ //if pressure at exact point
			drift_vel = v[i]; //use exact values
			trans_diff = dt[i];
			long_diff = dl[i];
			}
		}
	if (drift_vel ==0){ //if none of the exact values are found
		//make some splines and use those
		TSpline3* v_spline = new TSpline3("v_spline",p_bar,v,6);
		TSpline3* dt_spline = new TSpline3("dt_spline",p_bar,dt,6);
		TSpline3* dl_spline = new TSpline3("dl_spline",p_bar,dl,6);
		drift_vel = v_spline->Eval(pressure);
		trans_diff=dt_spline->Eval(pressure);
		long_diff=dl_spline->Eval(pressure);
		}	
	} else { //helium-4
		G4double v[6]={7680,5297,4367,3830,3437,3145};
		G4double dl[6]={0.0433,0.0321,0.0266,0.0243,0.0220,0.0209};
		G4double dt[6]={0.0601,0.0462,0.0381,0.0332,0.0307,0.0289};
		for (G4int i=0; i>6; i++){
			if (pressure == p_bar[i]){ //if pressure at exact point
				drift_vel = v[i]; //use exact values
				trans_diff = dt[i];
				long_diff = dl[i];
			}
		}
	if (drift_vel ==0){ //if no exact value found
		//make some splines and use those
		TSpline3* v_spline = new TSpline3("v_spline",p_bar,v,6);
		TSpline3* dt_spline = new TSpline3("dt_spline",p_bar,dt,6);
		TSpline3* dl_spline = new TSpline3("dl_spline",p_bar,dl,6);
		drift_vel = v_spline->Eval(pressure);
		trans_diff=dt_spline->Eval(pressure);
		long_diff=dl_spline->Eval(pressure);
		}		
	}
	G4cout<<name<<" "<<pressure<<" "<<drift_vel<<" "<<trans_diff<<" "<<long_diff<<G4endl;
}



