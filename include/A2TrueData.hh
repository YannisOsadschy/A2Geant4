#ifndef A2TrueData_h
#define A2TrueData_h

#include <vector>
#include <string>

#include "G4ThreeVector.hh"
#include "globals.hh"

//containers to save data collected in the various action classes during a run
//the containers can be analysed in the A2TrueDataAnalyser class 
//This interface was implemented as modular tool to get a brief and temporary overview of events, tracks, and steps, primarily for debugging of simulated processes
//the A2TrueData::structures and A2TrueDataAnalyser are not designed for detector data (see A2CBOutput), just for internal simulation data for testing
//the A2TrueData::structures and the public methods of A2TrueDataAnalyser can be edited as required.
//data collection can be turned on or off with the command: "/A2/run/collectTrueData true/false" default:false
struct StepData
{
    ///////////////////////////////////////////////////////////////
    double edep;
    double preKinEnergy;
    double postKinEnergy;
    ///////////////////////////////////////////////////////////////
    std::vector<int> secondariesTrackID; //do not remove
};

struct TrackData
{
    int trackID; //do not remove
    int parentTrackID; //do not remove
    //////////////////////////////////////////////////////////////
    int PDGE;
    double kinEnergy;
    double trackLength;
    std::string iVolumeName;    //intital volume
    std::string fVolumeName;    //final volume
    double iX, iY, iZ;          //initial local position
    double fX, fY, fZ;          //final local position
    //////////////////////////////////////////////////////////////
    std::vector<StepData> steps; //do not remove
};

struct EventData
{
    int eventID; //do not remove
    //////////////////////////////////////////////////////////////

    
    //////////////////////////////////////////////////////////////
    std::vector<TrackData> tracks; //do not remove
};

struct RunData
{
    std::vector<EventData> events; //do not remove
};

#endif