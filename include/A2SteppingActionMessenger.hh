
#ifndef A2SteppingActionMessenger_h
#define A2SteppingActionMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class A2SteppingAction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithoutParameter;
class G4UIcmdWithABool;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class A2SteppingActionMessenger: public G4UImessenger
{
  public:
    A2SteppingActionMessenger(A2SteppingAction*);
   ~A2SteppingActionMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    A2SteppingAction*  fSteppingAction;
    G4UIdirectory*     fSteppingDir;     
    G4UIcmdWithABool*  fSampleElectronsCmd;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
