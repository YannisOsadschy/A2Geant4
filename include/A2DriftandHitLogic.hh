//this class serves the purpose of having a common implementation of the parmetric drift and 
//hit detection logic used by the two different electron sources
//in the simulation of a TPC. Firstly there are Delta elctrons generated and tracked by Geant itself, 
//secondly in A2UserSteppingAction sampled electrons are manually created
//to deal with the large part of ionization, which was just continously accounted for by deposited energy, 
//but with no real tracks. Since there is no real need to create track
//objects for each and every sampled electron, which come with overhaead. 
//The drifting and hit detection logic was moved from the A2Driftmodell class, which has now solely the task of
//dealing with delta electrons, which have real track objects. 

#ifndef A2DriftandHitLogic_h
#define A2DriftandHitLogic_h 1

#include "G4ThreeVector.hh"
#include "Randomize.hh"
#include <array>
#include <vector>

class A2DriftandHitLogic
{
    public:
        A2DriftandHitLogic(G4Region*);
        ~A2DriftandHitLogic();
        struct TransportValues
        {
            double time;
            double pathLength;      //mm
            double eKin_keV;        //keV
            G4ThreeVector position; //mm
        };
        virtual TransportValues GetTransportValues(G4String, G4double, G4double, G4double, G4double, 
                                                G4double); 
                                                //move the electron through the active volume
        void ProcessHit(G4ThreeVector, G4double, G4double, G4int, G4int, G4double);
        void SampleEdep(const G4Step* aStep);

        //void SetSimulateElectronLoss(G4bool simulateElectronLoss) {fSimulateElectronLoss=simulateElectronLoss;}
        G4bool GetSimulateElectronLoss() const {return fSimulateElectronLoss;}
        G4bool DoesElectronSurvive(G4double) const;

    private:
        struct Point
        {
            double p,T,E,v,dL,dT;
        };

    	void SetConstants(G4Region*);
        void InterpolateDriftConstants(std::vector<Point>&, std::array<double,3>&, double);
        G4double fPressure;
        G4double fEfield;
        G4double fTemperature;       
        G4double fDriftVel; //drift velocity
    	G4double fLongDiff; //longitudinal diffusion
    	G4double fTransDiff; //transverse diffusion
    	//G4double fHePressure; //from TPC file -  to pick correct data
    	//G4int fHeIsotope;

    	G4Step* fFakeStep;
        G4Track* fFakeTrack;
    	G4TouchableHandle fTouchableHandle;
    	G4Navigator* fpNavigator;
    	G4bool fNaviSetup;
    	G4StepPoint* fFakePreStepPoint;
        G4StepPoint* fFakePostStepPoint;
        G4double fWorkFunction;

        G4bool fSimulateElectronLoss;
        G4double fLambdaAttachementFactor;
        G4double fDetectionSurvivalProbability;
        G4double fAnodeCathodeDistance;
};

#endif