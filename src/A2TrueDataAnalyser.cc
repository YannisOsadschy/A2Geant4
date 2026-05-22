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
        if (!all && i>0)
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
                out << "(TrackID=" << trackID <<"; Type=" << tracks[trackLookUpMap.at(trackID)].PDGE << "; T=" << tracks[trackLookUpMap.at(trackID)].kinEnergy << ")\n";
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

void A2TrueDataAnalyser::MakeEKinHists(int chosenDepth) const
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

        auto fillEnergies = [   &trackLookUpMap, &tracks, &chosenDepth, &selectedEKins, 
                                &selectedTrackIDs, &selectedParentTrackIDs, &selectedPDGEs](int trackID, std::size_t depth)
        {
            if (chosenDepth==-1 || depth==chosenDepth)
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
                                chosenDepth );
            }
        }
    }
    TFile file("data1.root", "Recreate");
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

void A2TrueDataAnalyser::MakeEdepEKinHists(int chosenDepth) const
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
                                        &chosenDepth, &sumChildrenEKins, &eKins, &edeps, &NSteps, &NSecondaries, this ]
                                            (int trackID, std::size_t depth)
        {
            
            if (chosenDepth==-1 || depth==chosenDepth)
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
                                chosenDepth );
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
    TFile file("data2.root", "Recreate");
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

void A2TrueDataAnalyser::MakePrimaryTrackLengthHists() const
{
    std::vector<double> trackLengths;
    std::vector<double> eKins;
    for (auto& event : fRunData.events)
    {
        for (const auto& track : event.tracks)
        {
            if (track.parentTrackID == 0)
            {   
                trackLengths.push_back(track.trackLength);
                eKins.push_back(track.kinEnergy);
            }
        }
    }
    TFile file("data3.root", "Recreate");
    double eKin;
    double trackLength;
    TTree tree("tree", "trackLength");
    tree.Branch("EKin", &eKin);
    tree.Branch("trackLength", &trackLength);
    for (size_t i = 0; i < eKins.size(); ++i)
    {
        eKin = eKins[i];
        trackLength = trackLengths[i];

        tree.Fill();
    }
    tree.Write();
    file.Close();

}


