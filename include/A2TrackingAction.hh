// A2TrackingAction
// Author: Dominik Werthmueller, 2017

#ifndef A2TrackingAction_h
#define A2TrackingAction_h 1

#include "G4UserTrackingAction.hh"

#include "A2TrueData.hh"

#include "A2EventAction.hh"

class A2PrimaryGeneratorAction;

class A2TrackingAction : public G4UserTrackingAction
{

private:
    A2PrimaryGeneratorAction* fPGA;     // pointer to generator

    TrackData currentTrackData; //is used for persistent dataTree

    A2EventAction* fEventAction; 

public:
    A2TrackingAction(A2EventAction*);
    virtual ~A2TrackingAction();

    virtual void PreUserTrackingAction(const G4Track*);
    virtual void PostUserTrackingAction(const G4Track*);

    TrackData& GetCurrentTrackData();
};

#endif

