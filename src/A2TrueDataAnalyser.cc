#include "A2TrueDataAnalyser.hh"
#include <iostream>
#include <fstream>
#include <unordered_map>

#include "TTree.h"
#include "TFile.h"



A2TrueDataAnalyser::A2TrueDataAnalyser(RunData&& runData)
    :fRunData(std::move(runData))
{}

A2TrueDataAnalyser::~A2TrueDataAnalyser()
{}

std::unordered_map<int, std::size_t> A2TrueDataAnalyser::MakeTrackLookUpMap(const EventData& event) const
{
    std::unordered_map<int, std::size_t> trackLookUpMap;
    trackLookUpMap.reserve(event.tracks.size()); 
    for (std::size_t i = 0; i < event.tracks.size(); ++i)
    {
        trackLookUpMap[event.tracks[i].trackID] = i;
    }
    return trackLookUpMap;
}

std::unordered_map<int, std::vector<int>> A2TrueDataAnalyser::MakeChildTrackMap(const EventData& event) const
{
    std::unordered_map<int, std::vector<int>> childTrackMap;
    childTrackMap.reserve(event.tracks.size());
    for (const auto& track : event.tracks)
    {
        childTrackMap[track.parentTrackID].push_back(track.trackID);
    }
    return childTrackMap;
}

double A2TrueDataAnalyser::GetEdepTrack(const TrackData& track) const
{
    double eDepTrack = 0;
    for (const auto& step : track.steps)
    {
        eDepTrack += step.edep;
    }
    return eDepTrack;
}


void A2TrueDataAnalyser::VisualizeTree(bool all) const
{   
    std::ofstream out;
    out.open("DataTreeVisualisation.txt");
    for (std::size_t i = 0; i<fRunData.events.size(); ++i)
    {   
        if (!all && i>10)
        {
            break;
        }
        const EventData& event = fRunData.events[i];
        //creating maps
        
        std::unordered_map<int, std::size_t> trackLookUpMap = MakeTrackLookUpMap(event);
        std::unordered_map<int, std::vector<int>> childTrackMap = MakeChildTrackMap(event);

        //output
        out << "(EventID: " << event.eventID << ")\n";
        const std::vector<TrackData>& tracks = event.tracks;
        auto writeTreeFunc = [&trackLookUpMap, &tracks, &out](int trackID, std::size_t depth)
            {
                out << "\t" << "\t";
                for (std::size_t i=0 ; i<depth ; ++i )
                {
                    out << "\t";
                };
                const TrackData& track = tracks[trackLookUpMap.at(trackID)];
                out << "(TrackID=" << trackID 
                << "; Type=" << track.PDGE 
                << "; T=" << track.kinEnergy 
                << "; length=" << track.trackLength
                << "; Nsteps=" << track.steps.size()
                << ")\n";
            };
            
        for (const auto& track : event.tracks)
        {
            if (track.parentTrackID == 0)
            {   
                ParseTrackTree( track.trackID,
                                childTrackMap,
                                writeTreeFunc );
            }
        }
    }
    out.close();
}

void A2TrueDataAnalyser::MakeEKinHists(std::size_t chosenLayer) const
{
    std::vector<double> selectedEKins;
    std::vector<int> selectedTrackIDs;
    std::vector<int> selectedParentTrackIDs;
    std::vector<int> selectedPDGEs;
    for (auto& event : fRunData.events)
    {
        std::unordered_map<int, std::size_t> trackLookUpMap = MakeTrackLookUpMap(event);
        std::unordered_map<int, std::vector<int>> childTrackMap = MakeChildTrackMap(event);
        const std::vector<TrackData>& tracks = event.tracks;

        auto fillEnergies = [   &trackLookUpMap, &tracks, &chosenLayer, &selectedEKins, 
                                &selectedTrackIDs, &selectedParentTrackIDs, &selectedPDGEs](int trackID, std::size_t depth)
        {
            if (chosenLayer == allLayers || depth==chosenLayer)
            {
                selectedEKins.push_back(tracks[trackLookUpMap.at(trackID)].kinEnergy);
                selectedParentTrackIDs.push_back(tracks[trackLookUpMap.at(trackID)].parentTrackID);
                selectedPDGEs.push_back(tracks[trackLookUpMap.at(trackID)].PDGE);
                selectedTrackIDs.push_back(trackID);
            }
        };

        for (const auto& track : event.tracks)
        {
            if (track.parentTrackID == 0)
            {   
                ParseTrackTree( track.trackID,
                                childTrackMap,
                                fillEnergies,
                                chosenLayer );  //chosenLayer passed as maxDepth, to stop recursion after reaching the chosen layer
            }
        }
    }
    TFile file("Ekin.root", "Recreate");
    double eKin;
    int ID;
    int parentID;
    int PDGE;
    TTree tree("tree", "EKin Data");
    tree.Branch("EKin", &eKin);
    tree.Branch("ID", &ID);
    tree.Branch("ParentID", &parentID);
    tree.Branch("PDGEncoding", &PDGE);
    for (size_t i = 0; i < selectedEKins.size(); ++i)
    {
        eKin = selectedEKins[i];
        ID = selectedTrackIDs[i];
        parentID = selectedParentTrackIDs[i];
        PDGE = selectedPDGEs[i];

        tree.Fill();
    }
    tree.Write();
    file.Close();
}

void A2TrueDataAnalyser::MakeEdepEKinHists(std::size_t chosenLayer) const
{   
    std::vector<double> edeps;
    std::vector<double> sumChildrenEKins;
    std::vector<double> eKins;
    std::vector<int> NSteps;
    std::vector<int> NSecondaries;
    std::vector<int> NLeafs;

    for (auto& event : fRunData.events)
    {   
        std::unordered_map<int, std::size_t> trackLookUpMap = MakeTrackLookUpMap(event);
        std::unordered_map<int, std::vector<int>> childTrackMap = MakeChildTrackMap(event);
        const std::vector<TrackData>& tracks = event.tracks;

        auto fillEnergyComparissons = [ &trackLookUpMap, &childTrackMap, &tracks, 
                                        &chosenLayer, &sumChildrenEKins, &eKins, &edeps, &NSteps, &NSecondaries, this ]
                                            (int trackID, std::size_t depth)
        {
            
            if (chosenLayer==allLayers || depth==chosenLayer)
            {   
                edeps.push_back(GetEdepTrack(tracks[trackLookUpMap.at(trackID)]));
                eKins.push_back(tracks[trackLookUpMap.at(trackID)].kinEnergy);
                NSteps.push_back(tracks[trackLookUpMap.at(trackID)].steps.size());
                auto it = childTrackMap.find(trackID);
                if (it != childTrackMap.end())
                {   
                    NSecondaries.push_back(it->second.size());
                    double sumChildrenEKin=0;
                    for (auto& childTrackID : it->second)
                    {
                        sumChildrenEKin += tracks[trackLookUpMap.at(childTrackID)].kinEnergy; 
                    }
                    sumChildrenEKins.push_back(sumChildrenEKin);

                }
                
            }
        }; 

        for (const auto& track : event.tracks)
        {
            if (track.parentTrackID == 0)
            {   
                ParseTrackTree( track.trackID,
                                childTrackMap,
                                fillEnergyComparissons,
                                chosenLayer );  //chosenLayer passed as maxDepth, to stop recursion after reaching the chosen layer
            }
        }
        int NLeaf = 0;
        auto getNumberofLeafs = [&childTrackMap, &NLeaf](int trackID, std::size_t depth)
        {
            auto it = childTrackMap.find(trackID);
            if (it == childTrackMap.end())
            {
                ++NLeaf;
            }
        };

        for (const auto& track : event.tracks)
        {
            if (track.parentTrackID == 0)
            {   
                ParseTrackTree( track.trackID,
                                childTrackMap,
                                getNumberofLeafs);
            }
        }
        NLeafs.push_back(NLeaf);

    }
    TFile file("edepEkin.root", "Recreate");
    double edep;
    double sumChildrenEKin;
    double eKin;
    double NStep;
    double NSecondary;
    double NLeaf;
    TTree tree("tree", "energy comparisson");
    tree.Branch("EKinPrim", &eKin);
    tree.Branch("EDepPrimSum", &edep);
    tree.Branch("EKinSecSum", &sumChildrenEKin);
    tree.Branch("NSecondaries", &NSecondary);
    tree.Branch("NLeaf", &NLeaf);
    tree.Branch("NStep", &NStep);
    for (size_t i = 0; i < eKins.size(); ++i)
    {
        sumChildrenEKin = sumChildrenEKins[i];
        eKin = eKins[i];
        edep = edeps[i];
        NStep = NSteps[i];
        NSecondary = NSecondaries[i];
        NLeaf = NLeafs[i];
        tree.Fill();
    }
    tree.Write();
    file.Close();
}

void A2TrueDataAnalyser::MakePrimaryTrackInfoHists(std::size_t chosenLayer) const
{   
    int PDGE;
    std::string iVolumeName;
    std::string fVolumeName;
    double iX;
    double iY;
    double iZ;
    double fX;
    double fY;
    double fZ;
    double kinEnergy;
    double trackLength;
    double iGlobalTime;
    double fGlobalTime;
    std::vector<double> iXSteps;
    std::vector<double> iYSteps;
    std::vector<double> iZSteps;
    std::vector<double> fXSteps;
    std::vector<double> fYSteps;
    std::vector<double> fZSteps;

    TFile file("tracks.root", "Recreate");
    TTree tree("tree", "energy comparisson");
    tree.Branch("PDGE", &PDGE);
    tree.Branch("initialVolumeName", &iVolumeName);
    tree.Branch("finalVolumeName", &fVolumeName);
    tree.Branch("iX", &iX);
    tree.Branch("iY", &iY);
    tree.Branch("iZ", &iZ);
    tree.Branch("fX", &fX);
    tree.Branch("fY", &fY);
    tree.Branch("fZ", &fZ);
    tree.Branch("Ekin", &kinEnergy);
    tree.Branch("trackLength", &trackLength);
    tree.Branch("iGlobalTime", &iGlobalTime);
    tree.Branch("fGlobalTime", &fGlobalTime);
    
    tree.Branch("iXSteps", &iXSteps);
    tree.Branch("iYSteps", &iYSteps);
    tree.Branch("iZSteps", &iZSteps);
    tree.Branch("fXSteps", &fXSteps);
    tree.Branch("fYSteps", &fYSteps);
    tree.Branch("fZSteps", &fZSteps);

    
    for (auto& event : fRunData.events)
    {   
        std::unordered_map<int, std::size_t> trackLookUpMap = MakeTrackLookUpMap(event);
        std::unordered_map<int, std::vector<int>> childTrackMap = MakeChildTrackMap(event);
        const std::vector<TrackData>& tracks = event.tracks;

        auto fillTrackInfo = [ &trackLookUpMap, &childTrackMap, &tracks, 
                               &chosenLayer, &PDGE, &iVolumeName, &fVolumeName,
                               &iX, &iY, &iZ, &fX, &fY, &fZ, &kinEnergy,
                               &trackLength, &iGlobalTime, &fGlobalTime, &iXSteps, &iYSteps, &iZSteps, &fXSteps, &fYSteps, &fZSteps, &tree, this ]
                                            (int trackID, std::size_t depth)
        {
            
            if (chosenLayer==allLayers || depth==chosenLayer)
            //if ((chosenLayer==allLayers || depth==chosenLayer) && depth != 0) //temporary no primaries
            {   
                const TrackData& track = tracks[trackLookUpMap.at(trackID)];
                PDGE = track.PDGE;
                iVolumeName = track.iVolumeName;
                fVolumeName = track.fVolumeName;
                iX = track.iX;
                iY = track.iY;
                iZ = track.iZ;
                fX = track.fX;
                fY = track.fY;
                fZ = track.fZ;
                kinEnergy = track.kinEnergy;
                trackLength = track.trackLength;
                iGlobalTime = track.iGlobalTime;
                fGlobalTime = track.fGlobalTime;
                
                for (auto& step : track.steps)
                {
                    iXSteps.push_back(step.iX);
                    iYSteps.push_back(step.iY);
                    iZSteps.push_back(step.iZ);
                    fXSteps.push_back(step.fX);
                    fYSteps.push_back(step.fY);
                    fZSteps.push_back(step.fZ);
                }
                
                tree.Fill();
                
                iXSteps.clear();
                iYSteps.clear();
                iZSteps.clear();
                fXSteps.clear();
                fYSteps.clear();
                fZSteps.clear();
            }
        }; 

        for (const auto& track : event.tracks)
        {
            if (track.parentTrackID == 0)
            {   
                ParseTrackTree( track.trackID,
                                childTrackMap,
                                fillTrackInfo,
                                chosenLayer );  //chosenLayer passed as maxDepth, to stop recursion after reaching the chosen layer
            }
        }
    }
    tree.Write();
    file.Close();
}


void A2TrueDataAnalyser::StepLengthPlots(std::size_t chosenLayer) const
{   
    double stepLength;
    double preStepKinEnergy;
    double accumulatedLength;
    TFile file("steps.root", "Recreate");
    TTree tree("tree", "stepLength");
    tree.Branch("stepLength", &stepLength);
    tree.Branch("preStepKinEnergy", &preStepKinEnergy);
    tree.Branch("addedLength", &accumulatedLength);
    const EventData& event = fRunData.events[0]; //only first event 
    std::unordered_map<int, std::size_t> trackLookUpMap = MakeTrackLookUpMap(event);
    std::unordered_map<int, std::vector<int>> childTrackMap = MakeChildTrackMap(event);
    const std::vector<TrackData>& tracks = event.tracks;
    auto fillStepInfo = [ &trackLookUpMap, &childTrackMap, &tracks, 
                               &chosenLayer, &stepLength, &preStepKinEnergy, &accumulatedLength, &tree, this ]
                                            (int trackID, std::size_t depth)
    {  
        if (chosenLayer==allLayers || depth==chosenLayer)
        //if ((chosenLayer==allLayers || depth==chosenLayer) && depth != 0) //temporary no primaries
        {   
            const TrackData& track = tracks[trackLookUpMap.at(trackID)];
            accumulatedLength=0;
            for (auto& step : track.steps)
            {
               stepLength = step.stepLength;
               accumulatedLength += step.stepLength; 
               preStepKinEnergy = step.preKinEnergy;
               tree.Fill();
            }
        }
    };
    for (const auto& track : event.tracks)
    {
        if (track.parentTrackID == 0)
        {   
            ParseTrackTree( track.trackID,
                            childTrackMap,
                            fillStepInfo,
                            chosenLayer );  //chosenLayer passed as maxDepth, to stop recursion after reaching the chosen layer
        }
    }
    tree.Write();
    file.Close();
}









void A2TrueDataAnalyser::StoppingPower(std::size_t chosenLayer) const
{   
    int PDGE;
    double trackEkin;
    double trackLength;
    std::vector<double> stepLength;
    std::vector<double> accumulatedStepLength;
    std::vector<double> preKinEnergy;
    std::vector<double> postKinEnergy;

    TFile file("steps.root", "Recreate");
    TTree tree("tree", "energy comparisson");
    tree.Branch("PDGE", &PDGE);
    tree.Branch("Ekin", &trackEkin);
    tree.Branch("trackLength", &trackLength);
    tree.Branch("stepLength", &stepLength);
    tree.Branch("accumulatedStepLength", &accumulatedStepLength);
    tree.Branch("preKinEnergy", &preKinEnergy);
    tree.Branch("postKinEnergy", &postKinEnergy);

    
    for (auto& event : fRunData.events)
    {   
        std::unordered_map<int, std::size_t> trackLookUpMap = MakeTrackLookUpMap(event);
        std::unordered_map<int, std::vector<int>> childTrackMap = MakeChildTrackMap(event);
        const std::vector<TrackData>& tracks = event.tracks;

        auto fillTrackInfo = [ &trackLookUpMap, &childTrackMap, &tracks, 
                               &chosenLayer, &PDGE, &stepLength, &accumulatedStepLength, &preKinEnergy, &postKinEnergy, &trackEkin, &trackLength, &tree, this ]
                                            (int trackID, std::size_t depth)
        {
            
            if (chosenLayer==allLayers || depth==chosenLayer)
            //if ((chosenLayer==allLayers || depth==chosenLayer) && depth != 0) //temporary no primaries
            {   
                const TrackData& track = tracks[trackLookUpMap.at(trackID)];
                PDGE = track.PDGE;
                trackEkin = track.kinEnergy;
                trackLength = track.trackLength;
                double sum = 0;
                for (auto& step : track.steps)
                {
                    stepLength.push_back(step.stepLength);
                    accumulatedStepLength.push_back(sum + 0.5*step.stepLength);
                    sum += step.stepLength;
                    preKinEnergy.push_back(step.preKinEnergy);
                    postKinEnergy.push_back(step.postKinEnergy);
                }
                tree.Fill();
                stepLength.clear();
                preKinEnergy.clear();
                postKinEnergy.clear();
                accumulatedStepLength.clear();
            }
        }; 

        for (const auto& track : event.tracks)
        {
            if (track.parentTrackID == 0)
            {   
                ParseTrackTree( track.trackID,
                                childTrackMap,
                                fillTrackInfo,
                                chosenLayer );  //chosenLayer passed as maxDepth, to stop recursion after reaching the chosen layer
            }
        }
    }
    tree.Write();
    file.Close();
}  
