#include "G4Region.hh"
#include "G4VPhysicalVolume.hh" //volumes
#include "G4VSensitiveDetector.hh"
#include "G4TransportationManager.hh" //particle transport

#include "TSpline.h"

#include "A2DriftandHitLogic.hh"
#include "A2UserTrackInformation.hh"
#include "TimeDebugger.hh"
#include "A2UserRegionInformation.hh"

#include "Randomize.hh"

#include <math.h>
#include <chrono>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <array>
#include <stdexcept>

#include "TH3D.h"
#include "TH2D.h"

using namespace CLHEP;

struct Point
{
    double p,T,E,v,dL,dT;
};

enum Parameter
{
    T, p, E
};





/**** Constructor *****/
A2DriftandHitLogic::A2DriftandHitLogic(G4Region* actVol)
{ //fast simulation implements user defined physics response
	//initiate pointers
	fFakeStep = new G4Step(); //step used to call hit in SD
    fFakeTrack = new G4Track(); //dummy Track just to hold A2UserTrackinformation, which is required for A2SD to work properly
    fFakeTrack->SetUserInformation(new A2UserTrackInformation());
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
    delete fFakeTrack;
  	delete fpNavigator;
}

//Generate electron position, time when reaching anode
A2DriftandHitLogic::TransportValues A2DriftandHitLogic::GetTransportValues(G4String particleName, double ekin_keV, double t, double x_mm,
                                            double y_mm, double z_mm)
{
    //auto t0 = std::chrono::high_resolution_clock::now();
    //G4cout<<"Transporting delta electron of energy "<< ekin_keV <<" keV"<<G4endl; //debugging message
	/****transport each electron to anode ****/
	G4double z_pos = -fAnodeCathodeDistance/2; //set final z position to anode z position
	G4double pathLength = z_pos - z_mm; //total length to final z position in mm

	//use Gaussians as defined by Fabian Metzger in his TPC work
	//set means and sigmas of Gaussian distribtuions
	G4double mean_x=x_mm; //mean for Gaussian calc of x pos
	G4double mean_y=y_mm; //mean for Gaussian calc of y pos
	G4double mean_t = abs(pathLength/fDriftVel); //mean for Gaussian calc of time: comes out in ms
	G4double sigma_diff = fTransDiff*sqrt(abs(pathLength)); //sigma for Gaussian calc of x,y positions: close enough to mm
	G4double sigma_time = fLongDiff/fDriftVel*sqrt(abs(pathLength)); //comes out close enough to ms

	//use random number generation to get values for positions and times
	G4double x_pos = RandGauss::shoot(mean_x,sigma_diff); //calculate an x position: mm
	G4double y_pos = RandGauss::shoot(mean_y,sigma_diff); //calc y mm
	G4double time = RandGauss::shoot(mean_t,sigma_time); //calc a time ms

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
    //auto t1 = std::chrono::high_resolution_clock::now();
    //TimeDebugger::getTransportValuesTime +=std::chrono::duration<double>(t1-t0).count();
    return transportValues;
}

/***** Call a hit in the anode for each electron that reaches it *****/
void A2DriftandHitLogic::ProcessHit(G4ThreeVector position, G4double ekin_keV, G4double drift_time, G4int trackID, G4int partID, G4double charge){
/**** set up touchable in current volume ****/
    //auto t0 = std::chrono::high_resolution_clock::now();
    if (!fNaviSetup) {
		fpNavigator->SetWorldVolume(G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking()->GetWorldVolume());
		fpNavigator->LocateGlobalPointAndUpdateTouchableHandle(position,G4ThreeVector(0.,0.,0.),fTouchableHandle,true);
		fNaviSetup = true;
	} else {
  		fpNavigator->
    	LocateGlobalPointAndUpdateTouchableHandle(position,G4ThreeVector(0.,0.,0.),fTouchableHandle);
 	}
    //auto tA = std::chrono::high_resolution_clock::now();
    //TimeDebugger::navigatorTime += std::chrono::duration<double>(tA-t0).count();
    A2UserTrackInformation* anInfo = static_cast<A2UserTrackInformation*>(fFakeTrack->GetUserInformation());
    anInfo->SetTrackID(trackID);
    anInfo->SetPartID(partID);
    anInfo->SetHasDriftParametersTPC(true);
    //TimeDebugger::A2SDIfTrue += std::chrono::duration<double>(tB-tA).count();
    //TPC driftParameters are currently passed as the following: A2DriftandHitLogic->A2UserTrackInformation->A2Hit->A2CBOutput for every track. 
    //To save memory it would better to move the determination of the Drift parameters from A2DriftandHitLogic to the A2TPC 
    //and pass the parameters to A2CBOutput via the A2DetectorConstruction
    anInfo->SetDriftParametersTPC(fPressure,fEfield,fDriftVel,fLongDiff,fTransDiff);
    //G4cout<<"vel"<<fDriftVel<<G4endl;
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
    //auto tB = std::chrono::high_resolution_clock::now();
    //TimeDebugger::inBetweenStuffTime += std::chrono::duration<double>(tB-tA).count();
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
    //auto t1 = std::chrono::high_resolution_clock::now();
    //TimeDebugger::processHitTime += std::chrono::duration<double>(t1-t0).count();
    //TimeDebugger::sdStuffTime += std::chrono::duration<double>(t1-tB).count();
}

//interface to manually start the sample drift process hit chain, for edep ion pairs in a TPC
void A2DriftandHitLogic::SampleEdep(const G4Step* aStep)
{   
    //auto t0 = std::chrono::high_resolution_clock::now();
    G4double meanNElec = aStep->GetTotalEnergyDeposit() / fWorkFunction;
    G4double nElec = std::round(RandPoisson::shoot(meanNElec));
    
    G4ThreeVector unitDirection = aStep->GetDeltaPosition().unit();
    G4ThreeVector preStepPosition = aStep->GetPreStepPoint()->GetPosition();
    G4double stepLength = aStep->GetStepLength();
    G4double deltaTime = aStep->GetDeltaTime();
    G4double preStepTime = aStep->GetPreStepPoint()->GetGlobalTime();
    G4int trackID = static_cast<A2UserTrackInformation*>(aStep->GetTrack()->GetUserInformation())->GetTrackID();
    G4int parentID = static_cast<A2UserTrackInformation*>(aStep->GetTrack()->GetUserInformation())->GetPartID();
    for (G4int i = 0; i < nElec; ++i)
    {
        G4double randFlat = RandFlat::shoot();
        G4ThreeVector positionElec = preStepPosition + randFlat * stepLength * unitDirection;
        G4double timeElec = preStepTime + randFlat * deltaTime;
        G4double eKin_keV = 0.000001;   //small value != 0 otherwise senitive detector does not work
        if (!fSimulateElectronLoss || DoesElectronSurvive(positionElec.z()))
	    {
            TransportValues transportValues = GetTransportValues("e-", eKin_keV, timeElec, positionElec.x(),
                                                positionElec.y(), positionElec.z());
            ProcessHit(transportValues.position, transportValues.eKin_keV, transportValues.time, trackID, parentID, -1);
        }
    //auto t1 = std::chrono::high_resolution_clock::now();
    //TimeDebugger::sampleEdepTime += std::chrono::duration<double>(t1-t0).count();
    }
}

G4bool A2DriftandHitLogic::DoesElectronSurvive(G4double zTrue) const
{
	G4double zDistance=zTrue+fAnodeCathodeDistance/2;
	G4double lambdaAttachement = fAnodeCathodeDistance*fLambdaAttachementFactor ; //free path for electron attachement
	G4double fAttachmentSurvivalProbability = exp(-zDistance/lambdaAttachement);
	G4double uniformX1 = RandFlat::shoot();  //losses along the path
	G4double uniformX2 = RandFlat::shoot();  //detector losses
	if (fAttachmentSurvivalProbability > uniformX1 && fDetectionSurvivalProbability > uniformX2)
	{
		return true;
	}
    else
    {
        return false;
    }
}


/**** Assign gas parameters depending on isotope, fPressure of helium ****/
void A2DriftandHitLogic::SetConstants(G4Region *gasRegion){
	G4String name=gasRegion->GetRootLogicalVolumeIterator()[0]->GetMaterial()->GetName();
	fPressure = gasRegion->GetRootLogicalVolumeIterator()[0]->GetMaterial()->GetPressure()/bar;
    A2UserRegionInformation* regionInfo = static_cast<A2UserRegionInformation*>(gasRegion->GetUserInformation());
    fEfield = regionInfo->GetEfield()/(volt/mm);
    fTemperature = regionInfo->GetTemperature();
    fSimulateElectronLoss = regionInfo->GetSimulateElectronLoss();
    fLambdaAttachementFactor = regionInfo->GetLambdaAttachementFactor();
    fDetectionSurvivalProbability= regionInfo->GetDetectionSurvivalProbability();
    fAnodeCathodeDistance = regionInfo->GetAnodeCathodeDistance();
    std::string inputString;
	if (name == "He3ActiveGas")
    { //helium-3
        fWorkFunction = 42.7e-6;
        inputString = "data/drift_cali_He3.tsv";
	} 
    else if (name=="He4ActiveGas")
    { //helium-4
        fWorkFunction = 42.7e-6;
		inputString = "data/drift_cali_He4.tsv";
	}
    else if (name=="D2GasPure")
    { //deuterium
        fWorkFunction = 36.5e-6 ;
		inputString = "data/drift_cali_D2.tsv";
	}


    double eps = 1e-8; //to supress numrical floating point errors
    std::ifstream file(inputString.c_str());

    if(!file.is_open())
    {
        throw std::runtime_error("failed to open drift data input-file");
    }

    std::array<double,3> selectedValues = {fTemperature,fPressure,fEfield};

    std::string header;
    std::getline(file,header);

    std::vector<Point> points;
    Point point;
    bool interpolBool = true;
    while (file >> point.T >> point.p >> point.E >> point.v >> point.dL >> point.dT)
    {
        if (std::abs(point.T-selectedValues[T])<eps && std::abs(point.p-selectedValues[p])<eps && std::abs(point.E-selectedValues[E])<eps  )
        //if (false)
        {
            fDriftVel = point.v;
            fLongDiff = point.dL;
            fTransDiff = point.dT;
            interpolBool = false;
            break;
        }
        points.push_back(point);
    }

    if(points.empty())
    {
        throw std::runtime_error("no calibration points found");
    }

    if (interpolBool)
    {
        InterpolateDriftConstants(points, selectedValues, eps);
    }
//G4cout<<name<<" "<<fPressure<<" "<<fDriftVel<<" "<<fTransDiff<<" "<<fLongDiff<<G4endl;
}

void A2DriftandHitLogic::InterpolateDriftConstants(std::vector<A2DriftandHitLogic::Point>& points, std::array<double,3>& selectedValues, double eps)
{
    //v, dL, dT will be calculated by interpolation on TH3 bins
    std::array<std::vector<double>,3> uniqueValueVectors;
    std::array<std::string,3> uniqueValueNames = {"T", "p", "E"};

    for (Point& point : points)
    {
        uniqueValueVectors[T].push_back(point.T);
        uniqueValueVectors[p].push_back(point.p);
        uniqueValueVectors[E].push_back(point.E);
    }
    for (size_t i=0; i<uniqueValueVectors.size(); ++i)
    {
        std::sort(uniqueValueVectors[i].begin(), uniqueValueVectors[i].end());
        std::vector<double> uniqueValueVectorsTmp;
        uniqueValueVectorsTmp.push_back(uniqueValueVectors[i][0]);
        for (size_t j=0; j<uniqueValueVectors[i].size(); ++j)
        {
            if (std::abs(uniqueValueVectors[i][j]-uniqueValueVectorsTmp[uniqueValueVectorsTmp.size()-1]) > eps)
            {
                uniqueValueVectorsTmp.push_back(uniqueValueVectors[i][j]);
            }
        }
        uniqueValueVectors[i] = uniqueValueVectorsTmp;
    }

    /////////////////////////////////////////////////////////////////////////////////////////////
    // a few checks to ensure the calibration data is in the from of a regular grid (i.e. all combinations (T,E,p) 
    // for given possible values {T_i}, {E_i}, {p_i} appear excatly onece and the possible values are equally spaced)
    //exactly once and the posiible values {T_i}, {E_i}, {p_i} are all equally spaced E)
    if (uniqueValueVectors[T].size()*uniqueValueVectors[p].size()*uniqueValueVectors[E].size() < points.size())
    {
        throw std::runtime_error("at least one combination of (T,p,E) occurs at least twice");
    }
    else if (uniqueValueVectors[T].size()*uniqueValueVectors[p].size()*uniqueValueVectors[E].size() > points.size())
    {
        throw std::runtime_error("not all combinations of (T,p,E) occur");
    }
    //check for duplicate points (part1):
    //from here nPoints=nT*np*nE. The only possible remaining faulty states have the same number of dublicate points 
    //as the number of missing possible combinations. Therefore if it can be rulled out, that there are duplicates, 
    //it can be concluded that every combination (T,p,E) exist exactly once in points.
    //Unfortunately i have no idea how to implemt such a check without nesting iterations over all points. Therefore I implement
    //this check further down in the code, where the programm iterates over all calibration points anyway.


    //check wheter the grid is regular with respect to T,p and E
    for (size_t i=0; i<uniqueValueVectors.size(); ++i)
    {
        if (uniqueValueVectors[i].size()>1)
        {

            double d0 = uniqueValueVectors[i][1] - uniqueValueVectors[i][0];
            bool equallySpaced = true;
            for (size_t j=0; j<uniqueValueVectors[i].size()-1; ++j)
            {    
                if (std::abs(d0 - (uniqueValueVectors[i][j+1]-uniqueValueVectors[i][j])) > eps)
                {
                    equallySpaced = false;
                    break;
                }
            }
            if (!equallySpaced)
            {
                throw std::runtime_error(uniqueValueNames[i] + " is not equally spaced") ;
            }
        }
    }   
    //creating the histogramms under the assumption that at this point grid is regular
    std::array<double, uniqueValueVectors.size()> binWidths;
    for (size_t i=0; i<uniqueValueVectors.size(); ++i)
    {
        binWidths[i] = 1;  //default bin width if there is just one unique Value
        if (uniqueValueVectors[i].size()>1)
        {
            binWidths[i] = uniqueValueVectors[i][1] - uniqueValueVectors[i][0];
        }
    }

    TH3D vGrid("vGrid", "vGrid",
        uniqueValueVectors[T].size(), uniqueValueVectors[T][0]-binWidths[T]/2, uniqueValueVectors[T][uniqueValueVectors[T].size()-1]+binWidths[T]/2,
        uniqueValueVectors[p].size(), uniqueValueVectors[p][0]-binWidths[p]/2, uniqueValueVectors[p][uniqueValueVectors[p].size()-1]+binWidths[p]/2,
        uniqueValueVectors[E].size(), uniqueValueVectors[E][0]-binWidths[E]/2, uniqueValueVectors[E][uniqueValueVectors[E].size()-1]+binWidths[E]/2    );
    
    TH3D dLGrid = vGrid;
    dLGrid.SetTitle("dLGrid");
    dLGrid.SetTitle("dLGrid");

    TH3D dTGrid = vGrid;
    dTGrid.SetTitle("dTGrid");
    dTGrid.SetTitle("dTGrid");
    
    std::vector<int> uniqueIDPoints; //necessary for duplicate check
    for (Point& point : points)
    {
        //grid layout is equal by design for v, dL and dT, therfore only v indeces have to be found
        int iT = vGrid.GetXaxis()->FindBin(point.T);
        if (abs(vGrid.GetXaxis()->GetBinCenter(iT)-point.T)<eps)
        {
            throw std::runtime_error("A bins center of the interpolation hist does not equal its respective calibration point [fTemperature]");
        }
        int ip = vGrid.GetYaxis()->FindBin(point.p);
        if (abs(vGrid.GetYaxis()->GetBinCenter(ip)-point.p)<eps)
        {
            throw std::runtime_error("A bins center of the interpolation hist does not equal its respective calibration point [fPressure]");
        }
        int iE = vGrid.GetZaxis()->FindBin(point.E);
        if (abs(vGrid.GetZaxis()->GetBinCenter(iE)-point.E)>eps)
        {
            throw std::runtime_error("A bins center of the interpolation hist does not equal its respective calibration point [efield]");
        }
        vGrid.SetBinContent(iT,ip,iE,point.v);
        dLGrid.SetBinContent(iT,ip,iE,point.dL);
        dTGrid.SetBinContent(iT,ip,iE,point.dT);
        ///////////////////////////////////////////////////////////////
        //check for duplicate points (part2)
        int nT = uniqueValueVectors[T].size();
        int np = uniqueValueVectors[p].size();
        uniqueIDPoints.push_back((iT-1) + (ip-1)*nT + (iE-1)*nT*np);  //root bin indices start at 1
    }
    //check for duplicate points (part3)
    std::sort(uniqueIDPoints.begin(), uniqueIDPoints.end());
    for (size_t i=0; i<uniqueIDPoints.size()-1; ++i)
    {
        if (uniqueIDPoints[i] == uniqueIDPoints[i+1])
        {
            throw std::runtime_error("there is at least one pair of a duplicate and missing combination");
        }
    }
    
    ////////////////////////////////////////////////////////
    //it needs to be checked wheter the selected values of T, p and E lie in the range covered by the calibration data. If the calibration data is in form of a regular rectangular grid the following test should be sufficient: 
    for (size_t i=0; i<uniqueValueVectors.size(); ++i)
    {
        if (uniqueValueVectors[i].size()==1)
        {
            if (std::abs(uniqueValueVectors[i][0] - selectedValues[i]) > eps)
            {
                throw std::runtime_error("only " + uniqueValueNames[i] + "=" + std::to_string(uniqueValueVectors[i][0]) +" is supported");
            }
        }
        else
        {
            if (selectedValues[i]<uniqueValueVectors[i][0] || selectedValues[i]>uniqueValueVectors[i][uniqueValueVectors[i].size()-1])
            {
                throw std::runtime_error(uniqueValueNames[i] + " lies not in the supported range: [" + 
                std::to_string(uniqueValueVectors[i][0]) + "," + std::to_string(uniqueValueVectors[i][uniqueValueVectors[i].size()-1]) + "]");
            }
        }
    }

    std::array<int,3> identifier = {0,0,0};    //0: size()=1  ,  1: size()>1
    if (vGrid.GetNbinsX() > 1)
    {
        identifier[T] = 1;
    }
    if (vGrid.GetNbinsY() > 1)
    {
        identifier[p] = 1;
    }
    if (vGrid.GetNbinsZ() > 1)
    {
        identifier[E] = 1;
    }
    switch (int dim=identifier[0]+identifier[1]+identifier[2])
    {
        case 0:
        {
            //no interpolation
            throw std::runtime_error("0d interpolation not yet implemented");
            break;
        }
        case 1:
        {
            //1d interpolation
            throw std::runtime_error("1d interpolation not yet implemented");
            break;
        }
        case 2:
        {
            //2d interpolation
            std::array<std::string,3> projectionKey = {"yz","xz","xy"};
            std::vector<int> remainingVariableIndex = {0,1,2};
            for (size_t i=0; i<identifier.size(); ++i)
            {
                if (identifier[i] == 0)
                {
                    TH2D* v2dGrid = (TH2D*)vGrid.Project3D(projectionKey[i].c_str());
                    TH2D* dL2dGrid = (TH2D*)dLGrid.Project3D(projectionKey[i].c_str());
                    TH2D* dT2dGrid = (TH2D*)dTGrid.Project3D(projectionKey[i].c_str());
                    remainingVariableIndex.erase(remainingVariableIndex.begin()+i);
                    fDriftVel = v2dGrid->Interpolate(selectedValues[remainingVariableIndex[0]], selectedValues[remainingVariableIndex[1]]);
                    fLongDiff = dL2dGrid->Interpolate(selectedValues[remainingVariableIndex[0]], selectedValues[remainingVariableIndex[1]]);
                    fTransDiff = dT2dGrid->Interpolate(selectedValues[remainingVariableIndex[0]], selectedValues[remainingVariableIndex[1]]);
                    delete v2dGrid;
                    delete dL2dGrid;
                    delete dT2dGrid;
                }
            }
            break;
        }
        case 3:
        {
            //3d interpolation
            fDriftVel = vGrid.Interpolate(selectedValues[T], selectedValues[p], selectedValues[E]);
            fLongDiff = dLGrid.Interpolate(selectedValues[T], selectedValues[p], selectedValues[E]);
            fTransDiff = dTGrid.Interpolate(selectedValues[T], selectedValues[p], selectedValues[E]);
            break;
        }
    }
}