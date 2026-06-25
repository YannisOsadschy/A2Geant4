
#include "A2SteppingActionMessenger.hh"
#include "A2SteppingAction.hh"

#include "A2SteppingAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithABool.hh"
#include "globals.hh"



A2SteppingActionMessenger::A2SteppingActionMessenger(A2SteppingAction* SteAct)
:fSteppingAction(SteAct)
{
  fSteppingDir = new G4UIdirectory("/A2/stepping/");
  fSteppingDir->SetGuidance("stepping control");

  fSampleElectronsCmd = new G4UIcmdWithABool("/A2/stepping/sampleElectrons", this);
  fSampleElectronsCmd->SetGuidance("Should electons of ion pairs be sampled depending of the deposited energy?");
  fSampleElectronsCmd->SetParameterName("sampleElec", true);
  fSampleElectronsCmd->SetDefaultValue(true);
  fSampleElectronsCmd->AvailableForStates(G4State_Idle);
}



A2SteppingActionMessenger::~A2SteppingActionMessenger()
{
  delete fSteppingDir;
  delete fSampleElectronsCmd;
}



void A2SteppingActionMessenger::SetNewValue(
                                        G4UIcommand* command,G4String newValue)
{ 
  if (command == fSampleElectronsCmd)
    fSteppingAction->SetSampleElectrons(fSampleElectronsCmd->GetNewBoolValue(newValue));
}



