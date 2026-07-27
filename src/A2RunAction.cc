
#include "A2RunAction.hh"
#include "A2RunActionMessenger.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"

#include "A2TrueDataAnalyser.hh"
#include "TimeDebugger.hh"
#include "G4ProductionCutsTable.hh"
#include "CLHEP/Units/SystemOfUnits.h" //units

using namespace CLHEP; //units


A2RunAction::A2RunAction()
{
  fEventAction=NULL;
  fCollectRunData=false;
  fRunMessenger=new A2RunActionMessenger(this);
}



A2RunAction::~A2RunAction()
{}



void A2RunAction::BeginOfRunAction(const G4Run* aRun)
{ 
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;

  //inform the runManager to save random number seed
  //G4RunManager::GetRunManager()->SetRandomNumberStore(true);

  //Open output file
  fEventAction=  const_cast<A2EventAction*>(static_cast<const A2EventAction*>(G4RunManager::GetRunManager()->GetUserEventAction()));
  fEventAction->PrepareOutput();
}




void A2RunAction::EndOfRunAction(const G4Run* aRun)
{
  G4int NbOfEvents = aRun->GetNumberOfEvent();
  if (NbOfEvents == 0) return;

  fEventAction->CloseOutput();
  if(fCollectRunData)
  {
    A2TrueDataAnalyser trueDataAnalyser = A2TrueDataAnalyser(std::move(currentRunData));
    trueDataAnalyser.VisualizeTree();
    trueDataAnalyser.MakePrimaryTrackInfoHists(0);
    //trueDataAnalyser.StepLengthPlots(0);
    //trueDataAnalyser.MakeEKinHists(1);
    //trueDataAnalyser.MakeEdepEKinHists(0);


    //debugging productioncut length to energy conversion
    //quick and dirty chatgpt code:
    /*
    auto* pct = G4ProductionCutsTable::GetProductionCutsTable();

    auto* ranges =
        pct->GetRangeCutsVector(G4ProductionCuts::GetIndex("e-"));

    auto* energies =
        pct->GetEnergyCutsVector(G4ProductionCuts::GetIndex("e-"));

    for (size_t i = 0; i < pct->GetTableSize(); ++i)
    {
        auto* couple = pct->GetMaterialCutsCouple(i);

        G4cout
            << i
            << ";material = " << couple->GetMaterial()->GetName()
            << ";range = " << (*ranges)[i]/mm
            << "mm"
            << ";energy = " << (*energies)[i]/keV
            << "keV"
            << G4endl;
    }
    */

  }

  
  
  /*
  G4cout<<"sampleEdepTime "<<TimeDebugger::sampleEdepTime<<G4endl;
  G4cout<<"getTransportValuesTime "<<TimeDebugger::getTransportValuesTime<<G4endl;
  G4cout<<"processHitTime "<<TimeDebugger::processHitTime<<G4endl;
  G4cout<<"navigatorTime "<<TimeDebugger::navigatorTime<<G4endl;
  G4cout<<"inBetweenStuffTime "<<TimeDebugger::inBetweenStuffTime<<G4endl;
  G4cout<<"sdStuffTime "<<TimeDebugger::sdStuffTime<<G4endl;
  G4cout<<"A2SDProcessHitsTime "<<TimeDebugger::A2SDProcessHitsTime<<G4endl;
  G4cout<<"A2SDSetupTime "<<TimeDebugger::A2SDSetup<<G4endl;
  G4cout<<"A2SDIfTrueTime "<<TimeDebugger::A2SDIfTrue<<G4endl;
  G4cout<<"A2SDIfFalseTime "<<TimeDebugger::A2SDIfFalse<<G4endl;
  G4cout<<"falseBlockEdepQdepTime "<<TimeDebugger::falseBlockEdepQdepTime<<G4endl;
  G4cout<<"falseBlockTPCBlockTime "<<TimeDebugger::falseBlockTPCBlockTime<<G4endl;
  G4cout<<"line1Time "<<TimeDebugger::line1<<G4endl;
  G4cout<<"line2Time "<<TimeDebugger::line2<<G4endl;
  G4cout<<"line3Time "<<TimeDebugger::line3<<G4endl;
  */

}

RunData& A2RunAction::GetCurrentRunData()
{
    return currentRunData;
}

