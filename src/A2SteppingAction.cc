

#include "A2SteppingAction.hh"

#include "A2DetectorConstruction.hh"
#include "A2EventAction.hh"

#include "G4Cerenkov.hh"
#include "G4Scintillation.hh"
#include "G4OpBoundaryProcess.hh"

#include "G4Track.hh"
#include "G4Gamma.hh"
#include "G4Proton.hh"
#include "G4OpticalPhoton.hh"
#include "G4SDManager.hh"
#include "G4SteppingManager.hh"
#include "G4RunManager.hh"
#include "CLHEP/Units/SystemOfUnits.h"
#include "G4FastSimulationManagerProcess.hh"
#include "A2DriftModel.hh"

#include "A2SteppingActionMessenger.hh"

#include "A2TrueData.hh"
#include "A2DriftandHitLogic.hh"

#include "A2UserRegionInformation.hh"



/*
#include "G4Version.hh"

#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4Track.hh"

#include "G4He3.hh"
#include "G4Alpha.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4Isotope.hh"
#include "G4IonisParamMat.hh"

#include "G4Region.hh"
#include "G4ProductionCuts.hh"

#include "G4EmCalculator.hh"
#include "G4NistManager.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "G4ios.hh"
*/

using namespace CLHEP;

////#include "G4RunManager.hh"



A2SteppingAction::A2SteppingAction(A2DetectorConstruction* det, A2EventAction* evt, A2TrackingAction* trc)
{
    detector = det;
    eventaction = evt;
    fTrackingAction = trc;
    fRegion = NULL;
    fFSManager = NULL;
    fDrifter = nullptr;
    fSampleElectrons = true;
    fSteppingMessenger = new A2SteppingActionMessenger(this);
}



A2SteppingAction::~A2SteppingAction()
{ }



void A2SteppingAction::UserSteppingAction(const G4Step* aStep)
{ 
/*  
    // ============================================================
// One-time He3 stopping/material diagnostic
//
// IMPORTANT:
// This code only PRINTS information.
// It does NOT modify particles, materials or physics.
// ============================================================

static G4bool he3DiagnosticPrinted = false;


{
    G4Track* track = aStep->GetTrack();
    if (!he3DiagnosticPrinted &&
    track->GetDefinition() == G4He3::He3Definition())
    {
        const G4StepPoint* pre = aStep->GetPreStepPoint();

        const G4Material* material = pre->GetMaterial();

        const G4Region* region = nullptr;

        if (pre->GetPhysicalVolume() &&
            pre->GetPhysicalVolume()->GetLogicalVolume())
        {
            region =
                pre->GetPhysicalVolume()
                ->GetLogicalVolume()
                ->GetRegion();
        }

        // Only print once when He3 is really inside ActiveGas
        if (region && region->GetName() == "ActiveGas")
        {
            he3DiagnosticPrinted = true;

            G4ParticleDefinition* he3 =
                G4He3::He3Definition();

            G4ParticleDefinition* alpha =
                G4Alpha::AlphaDefinition();


            G4cout << G4endl;
            G4cout << "============================================================" << G4endl;
            G4cout << "         A2 He3 STOPPING DIAGNOSTIC" << G4endl;
            G4cout << "============================================================" << G4endl;


            // ------------------------------------------------------------
            // Geant4 version
            // ------------------------------------------------------------

            G4cout << G4endl;
            G4cout << "===== GEANT4 =====" << G4endl;

            G4cout
                << "Version = "
                << G4Version
                << G4endl;


            // ------------------------------------------------------------
            // Actual track
            // ------------------------------------------------------------

            G4cout << G4endl;
            G4cout << "===== ACTUAL TRACK =====" << G4endl;

            G4cout
                << "Particle       = "
                << track->GetDefinition()->GetParticleName()
                << G4endl;

            G4cout
                << "PDG            = "
                << track->GetDefinition()->GetPDGEncoding()
                << G4endl;

            G4cout
                << "Track ID       = "
                << track->GetTrackID()
                << G4endl;

            G4cout
                << "Parent ID      = "
                << track->GetParentID()
                << G4endl;

            G4cout
                << "Kinetic energy = "
                << track->GetKineticEnergy()/MeV
                << " MeV"
                << G4endl;

            G4cout
                << "Charge         = "
                << track->GetDefinition()->GetPDGCharge()/eplus
                << " e"
                << G4endl;

            G4cout
                << "Mass           = "
                << track->GetDefinition()->GetPDGMass()/MeV
                << " MeV"
                << G4endl;


            // ------------------------------------------------------------
            // Particle definition
            // ------------------------------------------------------------

            G4cout << G4endl;
            G4cout << "===== HE3 PARTICLE DEFINITION =====" << G4endl;

            G4cout
                << "Name   = "
                << he3->GetParticleName()
                << G4endl;

            G4cout
                << "PDG    = "
                << he3->GetPDGEncoding()
                << G4endl;

            G4cout
                << "Charge = "
                << he3->GetPDGCharge()/eplus
                << " e"
                << G4endl;

            G4cout
                << "Mass   = "
                << he3->GetPDGMass()/MeV
                << " MeV"
                << G4endl;


            // ------------------------------------------------------------
            // Physics processes
            // ------------------------------------------------------------

            G4cout << G4endl;
            G4cout << "===== HE3 PROCESS LIST =====" << G4endl;

            if (he3->GetProcessManager())
            {
                he3->GetProcessManager()->DumpInfo();
            }
            else
            {
                G4cout << "ERROR: No process manager!" << G4endl;
            }


            // ------------------------------------------------------------
            // Geometry / region
            // ------------------------------------------------------------

            G4cout << G4endl;
            G4cout << "===== GEOMETRY / REGION =====" << G4endl;

            if (pre->GetPhysicalVolume())
            {
                G4cout
                    << "Physical volume = "
                    << pre->GetPhysicalVolume()->GetName()
                    << G4endl;
            }

            G4cout
                << "Region          = "
                << region->GetName()
                << G4endl;

            G4cout
                << "Material        = "
                << material->GetName()
                << G4endl;


            // ------------------------------------------------------------
            // Material macroscopic properties
            // ------------------------------------------------------------

            G4cout << G4endl;
            G4cout << "===== MATERIAL MACROSCOPIC PROPERTIES =====" << G4endl;

            G4cout
                << "Density       = "
                << material->GetDensity()/(g/cm3)
                << " g/cm3"
                << G4endl;

            G4cout
                << "Temperature   = "
                << material->GetTemperature()/kelvin
                << " K"
                << G4endl;

            G4cout
                << "Pressure      = "
                << material->GetPressure()/bar
                << " bar"
                << G4endl;

            G4cout
                << "State enum    = "
                << static_cast<G4int>(material->GetState())
                << G4endl;

            G4cout
                << "N elements    = "
                << material->GetNumberOfElements()
                << G4endl;

            G4cout
                << "Total atoms   = "
                << material->GetTotNbOfAtomsPerVolume()*cm3
                << " atoms/cm3"
                << G4endl;

            G4cout
                << "Total electrons = "
                << material->GetTotNbOfElectPerVolume()*cm3
                << " electrons/cm3"
                << G4endl;

            if (material->GetBaseMaterial())
            {
                G4cout
                    << "Base material = "
                    << material->GetBaseMaterial()->GetName()
                    << G4endl;
            }
            else
            {
                G4cout
                    << "Base material = NONE (custom material)"
                    << G4endl;
            }


            // ------------------------------------------------------------
            // Ionisation parameters
            // ------------------------------------------------------------

            if (material->GetIonisation())
            {
                G4cout
                    << "Mean excitation energy = "
                    << material->GetIonisation()
                            ->GetMeanExcitationEnergy()/eV
                    << " eV"
                    << G4endl;
            }


            // ------------------------------------------------------------
            // Element and isotope composition
            // ------------------------------------------------------------

            G4cout << G4endl;
            G4cout << "===== ELEMENT / ISOTOPE COMPOSITION =====" << G4endl;

            const G4double* fractions =
                material->GetFractionVector();

            const G4double* atomDensities =
                material->GetVecNbOfAtomsPerVolume();

            for (size_t i = 0;
                i < material->GetNumberOfElements();
                ++i)
            {
                const G4Element* element =
                    material->GetElement(i);

                G4cout << G4endl;

                G4cout
                    << "Element[" << i << "]"
                    << " name=" << element->GetName()
                    << " symbol=" << element->GetSymbol()
                    << G4endl;

                G4cout
                    << "  Z             = "
                    << element->GetZ()
                    << G4endl;

                G4cout
                    << "  A             = "
                    << element->GetA()/(g/mole)
                    << " g/mol"
                    << G4endl;

                if (fractions)
                {
                    G4cout
                        << "  mass fraction = "
                        << fractions[i]
                        << G4endl;
                }

                if (atomDensities)
                {
                    G4cout
                        << "  atom density  = "
                        << atomDensities[i]*cm3
                        << " atoms/cm3"
                        << G4endl;
                }


                // Isotopes contained in this element
                const size_t nIso =
                    element->GetNumberOfIsotopes();

                G4cout
                    << "  N isotopes    = "
                    << nIso
                    << G4endl;

                if (nIso > 0)
                {
                    G4IsotopeVector* isotopes =
                        element->GetIsotopeVector();

                    G4double* abundances =
                        element->GetRelativeAbundanceVector();

                    for (size_t j = 0; j < nIso; ++j)
                    {
                        const G4Isotope* iso =
                            (*isotopes)[j];

                        G4cout
                            << "    isotope[" << j << "] "
                            << iso->GetName()
                            << " Z=" << iso->GetZ()
                            << " N=" << iso->GetN()
                            << " A=" << iso->GetA()/(g/mole)
                            << " g/mol";

                        if (abundances)
                        {
                            G4cout
                                << " abundance="
                                << abundances[j];
                        }

                        G4cout << G4endl;
                    }
                }
            }


            // ------------------------------------------------------------
            // Independent ideal-gas number-density check
            //
            // For PURE MONATOMIC He3:
            //
            //     n = p / (k_B T)
            //
            // This should agree with GetTotNbOfAtomsPerVolume().
            // ------------------------------------------------------------

            G4cout << G4endl;
            G4cout << "===== IDEAL GAS CONSISTENCY CHECK =====" << G4endl;

            const G4double p =
                material->GetPressure();

            const G4double T =
                material->GetTemperature();

            if (p > 0.0 && T > 0.0)
            {
                const G4double nIdeal =
                    p / (k_Boltzmann * T);

                const G4double nG4 =
                    material->GetTotNbOfAtomsPerVolume();

                G4cout
                    << "n from p/(kB*T) = "
                    << nIdeal*cm3
                    << " atoms/cm3"
                    << G4endl;

                G4cout
                    << "n stored by G4  = "
                    << nG4*cm3
                    << " atoms/cm3"
                    << G4endl;

                G4cout
                    << "ratio G4/ideal  = "
                    << nG4/nIdeal
                    << G4endl;

                G4cout
                    << "(For pure monatomic He3 this ratio should be ~1.)"
                    << G4endl;
            }


            // ------------------------------------------------------------
            // Production cuts
            // ------------------------------------------------------------

            G4cout << G4endl;
            G4cout << "===== PRODUCTION CUTS IN ActiveGas =====" << G4endl;

            const G4ProductionCuts* cuts =
                region->GetProductionCuts();

            if (cuts)
            {
                const G4int iGamma =
                    G4ProductionCuts::GetIndex("gamma");

                const G4int iElectron =
                    G4ProductionCuts::GetIndex("e-");

                const G4int iPositron =
                    G4ProductionCuts::GetIndex("e+");

                const G4int iProton =
                    G4ProductionCuts::GetIndex("proton");

                G4cout
                    << "gamma cut  = "
                    << cuts->GetProductionCut(iGamma)/mm
                    << " mm"
                    << G4endl;

                G4cout
                    << "e- cut     = "
                    << cuts->GetProductionCut(iElectron)/mm
                    << " mm"
                    << G4endl;

                G4cout
                    << "e+ cut     = "
                    << cuts->GetProductionCut(iPositron)/mm
                    << " mm"
                    << G4endl;

                G4cout
                    << "proton cut = "
                    << cuts->GetProductionCut(iProton)/mm
                    << " mm"
                    << G4endl;
            }
            else
            {
                G4cout
                    << "No production cuts attached to region."
                    << G4endl;
            }


            // ------------------------------------------------------------
            // Geant4 stopping / range values
            // ------------------------------------------------------------

            G4cout << G4endl;
            G4cout << "===== GEANT4 INTERNAL STOPPING / RANGE =====" << G4endl;

            G4EmCalculator emCalc;

            G4cout
                << " E[MeV]"
                << "  model-total-dEdx[MeV/mm]"
                << "  table-dEdx[MeV/mm]"
                << "  restricted-range[mm]"
                << "  CSDA-range[mm]"
                << G4endl;

            for (G4double E = 1.0;
                E <= 10.0;
                E += 1.0)
            {
                const G4double kineticEnergy =
                    E*MeV;

                // Direct evaluation of total electronic stopping
                const G4double modelDEDX =
                    emCalc.ComputeTotalDEDX(
                        kineticEnergy,
                        he3,
                        material
                    );

                // Value from Geant4's energy-loss tables for
                // this material/region
                const G4double tableDEDX =
                    emCalc.GetDEDX(
                        kineticEnergy,
                        he3,
                        material,
                        region
                    );

                const G4double restrictedRange =
                    emCalc.GetRangeFromRestricteDEDX(
                        kineticEnergy,
                        he3,
                        material,
                        region
                    );

                const G4double csdaRange =
                    emCalc.GetCSDARange(
                        kineticEnergy,
                        he3,
                        material,
                        region
                    );

                G4cout
                    << E
                    << "  "
                    << modelDEDX/(MeV/mm)
                    << "  "
                    << tableDEDX/(MeV/mm)
                    << "  "
                    << restrictedRange/mm
                    << "  "
                    << csdaRange/mm
                    << G4endl;
            }


            // ------------------------------------------------------------
            // Direct He3 vs alpha comparison in EXACT SAME MATERIAL
            // ------------------------------------------------------------

            G4cout << G4endl;
            G4cout << "===== HE3 vs ALPHA IN SAME MATERIAL =====" << G4endl;

            for (G4double E :
                {2.0, 4.0, 6.0, 8.0, 10.0})
            {
                const G4double he3DEDX =
                    emCalc.ComputeTotalDEDX(
                        E*MeV,
                        he3,
                        material
                    );

                const G4double alphaDEDX =
                    emCalc.ComputeTotalDEDX(
                        E*MeV,
                        alpha,
                        material
                    );

                G4cout
                    << "E=" << E << " MeV"
                    << "  He3=" << he3DEDX/(MeV/mm)
                    << " MeV/mm"
                    << "  alpha=" << alphaDEDX/(MeV/mm)
                    << " MeV/mm"
                    << G4endl;
            }


            // ------------------------------------------------------------
            // NIST helium cross-check
            //
            // Compare stopping per atom by scaling G4_He to the same
            // atom density as our custom He3 material.
            // This does NOT modify the simulation.
            // ------------------------------------------------------------

            G4cout << G4endl;
            G4cout << "===== NIST G4_He CROSS-CHECK =====" << G4endl;

            G4Material* nistHe =
                G4NistManager::Instance()
                    ->FindOrBuildMaterial("G4_He");

            if (nistHe)
            {
                const G4double nCustom =
                    material->GetTotNbOfAtomsPerVolume();

                const G4double nNist =
                    nistHe->GetTotNbOfAtomsPerVolume();

                G4cout
                    << "NIST G4_He density = "
                    << nistHe->GetDensity()/(g/cm3)
                    << " g/cm3"
                    << G4endl;

                G4cout
                    << "NIST atom density  = "
                    << nNist*cm3
                    << " atoms/cm3"
                    << G4endl;

                for (G4double E :
                    {2.0, 4.0, 6.0, 8.0, 10.0})
                {
                    const G4double customDEDX =
                        emCalc.ComputeTotalDEDX(
                            E*MeV,
                            he3,
                            material
                        );

                    const G4double nistDEDX =
                        emCalc.ComputeTotalDEDX(
                            E*MeV,
                            he3,
                            nistHe
                        );

                    const G4double nistSameAtomDensity =
                        nistDEDX * nCustom/nNist;

                    G4cout
                        << "E=" << E << " MeV"
                        << " custom-He3="
                        << customDEDX/(MeV/mm)
                        << " MeV/mm"
                        << "  NIST-He(same atom density)="
                        << nistSameAtomDensity/(MeV/mm)
                        << " MeV/mm"
                        << G4endl;
                }
            }


            G4cout << G4endl;
            G4cout << "============================================================" << G4endl;
            G4cout << "       END A2 He3 STOPPING DIAGNOSTIC" << G4endl;
            G4cout << "============================================================" << G4endl;
            G4cout << G4endl;
        }
    }

}
*/






    G4Track* track = aStep->GetTrack();

    G4StepPoint* startPoint = aStep->GetPreStepPoint();
    G4StepPoint* endPoint   = aStep->GetPostStepPoint();

    G4String particleName = track->GetDynamicParticle()->GetParticleDefinition()->GetParticleName();

    //G4cout << particleName << "\t" << pds->GetProcessName() << "\t" << startPoint->GetPhysicalVolume()->GetName()<< " to " << endPoint->GetPhysicalVolume()->GetName();
    //G4cout << "\t" << startPoint->GetPhysicalVolume()->GetLogicalVolume()->GetMaterial()->GetMaterialPropertiesTable()->GetConstProperty("SCINTILLATIONYIELD");
    //G4cout << G4endl;
    
    //Need to add some stuff here for sure! Follow AHeT code
    if (particleName == "e-") 
    {
        //propogate through field
        //and see if you hit the anode
    //G4https://root.cern/manual/histograms/histo-trial.png 
    }

    if (particleName == "opticalphoton") 
    {
        G4OpBoundaryProcess* fOpProcess;

        // Retrieve the status of the photon
        G4OpBoundaryProcessStatus theStatus = Undefined;

        G4ProcessManager* OpManager =
                G4OpticalPhoton::OpticalPhoton()->GetProcessManager();

        if (OpManager) 
        {
            G4int MAXofPostStepLoops =
                    OpManager->GetPostStepProcessVector()->entries();
            G4ProcessVector* fPostStepDoItVector =
                    OpManager->GetPostStepProcessVector(typeDoIt);

            for ( G4int i=0; i<MAXofPostStepLoops; i++) 
            {
                G4VProcess* fCurrentProcess = (*fPostStepDoItVector)[i];
                fOpProcess = dynamic_cast<G4OpBoundaryProcess*>(fCurrentProcess);
                if (fOpProcess) { theStatus = fOpProcess->GetStatus(); break;}
            }
        }

        // Detected by a detector
        if (theStatus == Detection) 
        {
            //G4cout << "Detected in " << startPoint->GetPhysicalVolume()->GetName() << "\t" << aStep->GetTotalEnergyDeposit()/eV << G4endl;

            // Check if the photon hits the detector and process the hit if it does
            if ( startPoint->GetPhysicalVolume()->GetLogicalVolume()->GetName() == "LogicSiPMT" ) 
            {

                G4SDManager* SDman = G4SDManager::GetSDMpointer();
                A2SD* AHe3SD = (A2SD*)SDman->FindSensitiveDetector("AHe3SD");

                //if (AHe3SD) AHe3SD->ProcessHits_AHe3(aStep, NULL);

                // Stop Tracking when it hits the detector's surface
                //ResetCounters();
                track->SetTrackStatus(fStopAndKill);
            }
        }
    }

    //   G4VPhysicalVolume* volume = track->GetVolume();
    
    //   // collect energy and track length step by step
    G4double edep = aStep->GetTotalEnergyDeposit();
    
    //   G4double stepl = 0.;

    /*
    //force electrons to drift in Time Projection chamber
    if(track->GetDefinition()->GetParticleName()==G4String("e-")){ //only applies to electrons
        if(fpSteppingManager->GetfCurrentVolume()->GetName()=="HELIUM" && track->GetTrackStatus()!=fStopAndKill){ //check for active electron inside active volume 
            //G4cout<<"Drift electron..."<<G4endl; //debugging progress message
            //get the fast simulation manager for the active volume
            if (!fRegion) fRegion = fpSteppingManager->GetfCurrentVolume()->GetLogicalVolume()->GetRegion();
            if (!fFSManager) fFSManager = fRegion->GetFastSimulationManager();
            //breaks for zero eenergy: give just a bit
            if(track->GetKineticEnergy()==0)track->SetKineticEnergy(1*eV);
            //G4cout<<"Checking for trigger..."<<G4endl; //debugging progress message
            if(fFSManager->PostStepGetFastSimulationManagerTrigger(*track)){ //check if Heed process applies
                G4VParticleChange* fastStep = fFSManager->InvokePostStepDoIt(); //force heed process to occur
                track->SetTrackStatus(fStopAndKill); //stop particle after running Heed process
                //do nothing: currently test-running simulation without this
            }
        }	
    }
    else if(aStep->GetPreStepPoint()->GetGlobalTime()>2*ms)track->SetTrackStatus(fStopAndKill);
    */
    if(aStep->GetPreStepPoint()->GetGlobalTime()>2*ms)track->SetTrackStatus(fStopAndKill);
    //stop tracking after the trigger time
    //   if(track->GetDefinition()->GetParticleName()==G4String("pi0"))
    //     {G4cout<<"Got a pi0 "<<aStep->GetPreStepPoint()->GetGlobalTime()/ns<<" "<<track->GetKineticEnergy()/MeV<<" "<< fpSteppingManager->GetfCurrentVolume()->GetName()<<G4endl;track->SetTrackStatus(fStopAndKill);}
    //  if(track->GetDefinition()->GetParticleName()==G4String("pi+"))
    //    {G4cout<<"Got a pi+ "<<aStep->GetPreStepPoint()->GetGlobalTime()/ns<<" "<<track->GetKineticEnergy()/MeV<<" "<< fpSteppingManager->GetfCurrentVolume()->GetName()<<G4endl;track->SetTrackStatus(fStopAndKill);}
    //  if(track->GetDefinition()->GetParticleName()==G4String("mu+"))
    //    {G4cout<<"Got a mu+ "<<aStep->GetPreStepPoint()->GetGlobalTime()/ns<<" "<<track->GetKineticEnergy()/MeV<<" "<< fpSteppingManager->GetfCurrentVolume()->GetName()<<G4endl;track->SetTrackStatus(fStopAndKill);}
    //  if(track->GetDefinition()->GetParticleName()==G4String("pi0"))
    //     {track->SetTrackStatus(fStopAndKill);}
    //  if(track->GetDefinition()->GetParticleName()==G4String("pi+"))
    //    {track->SetTrackStatus(fStopAndKill);}
    //  if(track->GetDefinition()->GetParticleName()==G4String("pi-"))
    //    {track->SetTrackStatus(fStopAndKill);}
    //  if(track->GetDefinition()->GetParticleName()==G4String("mu+"))
    //    {track->SetTrackStatus(fStopAndKill);}
    //if (track->GetDefinition()->GetPDGCharge() != 0.)
    // stepl = aStep->GetStepLength();
    //if(track->GetDefinition()->GetParticleName()==G4Gamma::Gamma()->GetParticleName()){

    //G4cout<<track->GetDefinition()->GetParticleName()<< track->GetTrackID()<<" process " <<fpSteppingManager->GetfCurrentProcess()->GetProcessName()<<" in "<< fpSteppingManager->GetfCurrentVolume()->GetName()<<G4endl;
    // }
    //  if(fpSteppingManager->GetfCurrentVolume()->GetName()==G4String("ANOIP"))G4cout<<"OK "<<track->GetDefinition()->GetParticleName()<<G4endl;
    //if(track->GetDefinition()->GetParticleName()==G4Proton::Proton()->GetParticleName()&&!(fpSteppingManager->GetfCurrentProcess()->GetProcessName()==G4String("msc"))&&!(fpSteppingManager->GetfCurrentProcess()->GetProcessName()==G4String("hIoni"))){
    // if(fpSteppingManager->GetfCurrentProcess()->GetProcessName()==G4String("msc"))
    // if(track->GetTrackID()==1&&fpSteppingManager->GetfCurrentVolume()->GetName()!=G4String("World")) G4cout<<track->GetDefinition()->GetParticleName()<< track->GetTrackID()<< " "<<track->GetParentID()<< " "<<track->GetKineticEnergy()<<" process " <<fpSteppingManager->GetfCurrentProcess()->GetProcessName()<<" in "<< fpSteppingManager->GetfCurrentVolume()->GetName()<<G4endl;
        //G4cout<<"Secondaries "<<aStep->GetSecondary()->size()<<" "<<aStep->GetfSecondary()<<G4endl;
    //}
    //G4cout <<" STEPPING ACTION "<<eventaction->GetNEvent()
    //  if(eventaction->GetNEvent()==1317){
    // G4cout<<track->GetDefinition()->GetParticleName()<< track->GetTrackID()<<" process " <<fpSteppingManager->GetfCurrentProcess()->GetProcessName()<<" in "<< fpSteppingManager->GetfCurrentVolume()->GetName()<< " "<<track->GetKineticEnergy()/MeV<<G4endl;
    // }

    //bug in phot process, can't get rid of gamma with energy 1.2E-5MeV
    //goes into infinite loop!
    if(track->GetDefinition()->GetParticleName()==G4Gamma::Gamma()->GetParticleName()&&track->GetKineticEnergy()/MeV<1E-4&&fpSteppingManager->GetfCurrentProcess()->GetProcessName()==G4String("phot"))track->SetTrackStatus(fStopAndKill);

    StepData stepData;
    stepData.edep = aStep->GetTotalEnergyDeposit();
    stepData.preKinEnergy = aStep->GetPreStepPoint()->GetKineticEnergy();
    stepData.postKinEnergy = aStep->GetPostStepPoint()->GetKineticEnergy();
    stepData.stepLength = aStep->GetStepLength();
    //stepData.volumeName = aStep->GetPreStepPoint()->GetPhysicalVolume()->;
    const G4TouchableHandle& touchableI = aStep->GetPreStepPoint()->GetTouchableHandle();
    const G4ThreeVector& iPosition = touchableI->GetHistory()->GetTopTransform().TransformPoint(aStep->GetPreStepPoint()->GetPosition());
    stepData.iX = iPosition.x(); 
    stepData.iY = iPosition.y(); 
    stepData.iZ = iPosition.z();

    const G4TouchableHandle& touchableF = aStep->GetPostStepPoint()->GetTouchableHandle();
    const G4ThreeVector& fPosition = touchableF->GetHistory()->GetTopTransform().TransformPoint(aStep->GetPostStepPoint()->GetPosition());
    stepData.fX = fPosition.x(); 
    stepData.fY = fPosition.y(); 
    stepData.fZ = fPosition.z();

    for (const auto& childTrack : *(aStep->GetSecondary()))
    {
        stepData.secondariesTrackID.push_back(childTrack->GetTrackID());
    }
    fTrackingAction->GetCurrentTrackData().steps.push_back(stepData);
    
    if (fSampleElectrons)
    {
        G4StepPoint* preStep = aStep->GetPreStepPoint();
        G4LogicalVolume* logicalVolume = preStep->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
        G4Region* region =logicalVolume->GetRegion();
        G4String particleName = aStep->GetTrack()->GetParticleDefinition()->GetParticleName();

        if (region && region->GetName() == "ActiveGas")
        {
            A2UserRegionInformation* regionInfo = static_cast<A2UserRegionInformation*>(region->GetUserInformation());
            const G4double z =preStep->GetPosition().z();
            const G4double anodeZ = regionInfo->GetAnodeZ();
            const G4double cathodeZ = regionInfo->GetCathodeZ();
            const bool insideDriftVolume = z>anodeZ && z<cathodeZ;
            const G4String materialName = logicalVolume->GetMaterial()->GetName();
            const bool isTPCGas = materialName == "He3ActiveGas" || materialName == "He4ActiveGas" || materialName == "D2GasPure";
            const bool isTPCSignalParticle = particleName == "alpha" || particleName == "He3" || particleName == "deuteron";
            if (insideDriftVolume && isTPCGas && isTPCSignalParticle)
            {
                if (!fDrifter)
                {
                    fDrifter =
                        new A2DriftandHitLogic(region);
                }

                fDrifter->SampleEdep(aStep);
            }
        }
    }
}

void A2SteppingAction::SetSampleElectrons(bool sampleElectrons)
{
    fSampleElectrons = sampleElectrons;
}