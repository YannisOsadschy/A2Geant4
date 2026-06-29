
#ifndef A2RunActionMessenger_h
#define A2RunActionMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class A2RunAction;
class G4UIdirectory;
class G4UIcmdWithABool;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class A2RunActionMessenger: public G4UImessenger
{
  public:
    A2RunActionMessenger(A2RunAction*);
   ~A2RunActionMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    A2RunAction*  fRunAction;
    G4UIdirectory*     fRunDir;     
    G4UIcmdWithABool*  fTrueDataCollectionCmd;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
