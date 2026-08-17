#ifndef A2UserRegionInformation_h
#define A2UserRegionInformation_h 1

#include "G4VUserRegionInformation.hh"


//this class is solely used as a container to pass the values of the constant Efield and temperature in the TPC (implemented in A2TPC) to the A2DriftandHitLogic class for parametric drift via a region object.
class A2UserRegionInformation : public G4VUserRegionInformation
{
    public:
        void SetEfield(double Efield) {fEfield = Efield;}
        double GetEfield() const {return fEfield;}
        void SetTemperature(double temperature) {fTemperature = temperature;}
        double GetTemperature() const {return fTemperature;}
        void SetSimulateElectronLoss(bool simulateElectronLoss){fSimulateElectronLoss=simulateElectronLoss;}
	    void SetLambdaAttachementFactor(double lambdaAttachementFactor){fLambdaAttachementFactor=lambdaAttachementFactor;}
	    void SetDetectionSurvivalProbability(double detectionSurvivalProbability){fDetectionSurvivalProbability=detectionSurvivalProbability;}
        void SetAnodeCathodeDistance(double anodeCathodeDistance){fAnodeCathodeDistance=anodeCathodeDistance;}
        bool GetSimulateElectronLoss() const{return fSimulateElectronLoss;}
        double GetLambdaAttachementFactor() const{return fLambdaAttachementFactor;}
        double GetDetectionSurvivalProbability() const{return fDetectionSurvivalProbability;}
        double GetAnodeCathodeDistance() const{return fAnodeCathodeDistance;}

        void Print() const override {}

    private:
        double fEfield;
        double fTemperature;
        bool fSimulateElectronLoss;
        double fLambdaAttachementFactor;
        double fDetectionSurvivalProbability;
        double fAnodeCathodeDistance;


};

#endif
