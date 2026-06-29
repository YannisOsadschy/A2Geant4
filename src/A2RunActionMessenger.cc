#include "A2RunActionMessenger.hh"
#include "A2RunAction.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithABool.hh"
#include "globals.hh"



A2RunActionMessenger::A2RunActionMessenger(A2RunAction* SteAct)
:fRunAction(SteAct)
{
  fRunDir = new G4UIdirectory("/A2/run/");
  fRunDir->SetGuidance("run control");

  fTrueDataCollectionCmd = new G4UIcmdWithABool("/A2/run/collectTrueData", this);
  fTrueDataCollectionCmd->SetGuidance("Should true data be collected");
  fTrueDataCollectionCmd->SetParameterName("collect", false);
  fTrueDataCollectionCmd->SetDefaultValue(false);
  fTrueDataCollectionCmd->AvailableForStates(G4State_Idle);
}



A2RunActionMessenger::~A2RunActionMessenger()
{
  delete fRunDir;
  delete fTrueDataCollectionCmd;
}



void A2RunActionMessenger::SetNewValue(
                                        G4UIcommand* command,G4String newValue)
{ 
  if (command == fTrueDataCollectionCmd)
    fRunAction->SetCollectRunData(fTrueDataCollectionCmd->GetNewBoolValue(newValue));
}



