
#ifndef A2RunAction_h
#define A2RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"
#include "A2EventAction.hh"

#include "A2TrueData.hh"
#include "TimeDebugger.hh"

class G4Run;
class A2RunActionMessenger;

class A2RunAction : public G4UserRunAction
{
    public:
        A2RunAction();
        ~A2RunAction();

    public:
        void BeginOfRunAction(const G4Run*);
        void EndOfRunAction(const G4Run*);

        void inline SetCollectRunData(G4bool collectRunData ) {fCollectRunData = collectRunData;}
        G4bool inline GetCollectRunData() {return fCollectRunData;}
    private:
        A2EventAction *fEventAction;

        RunData currentRunData;
        G4bool fCollectRunData;

        TimeDebugger fTimeDebugger;

        A2RunActionMessenger* fRunMessenger;



    public:
        RunData& GetCurrentRunData();
};

#endif

