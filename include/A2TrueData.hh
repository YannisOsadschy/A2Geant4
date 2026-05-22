#ifndef A2TrueData_h
#define A2TrueData_h

#include <vector>

#include "G4ThreeVector.hh"
#include "globals.hh"


struct StepData
{
    //std::string volumeName;
    double edep;
    double preKinEnergy;
    double postKinEnergy;
    std::vector<int> secondariesTrackID;
};

struct TrackData
{
    int trackID;
    int parentTrackID;
    int PDGE;
    double kinEnergy;
    double trackLength;
    std::vector<StepData> steps;
};

struct EventData
{
    int eventID;

    std::vector<int> primaryParticlesID;
    
    std::vector<TrackData> tracks;
};

struct RunData
{
    std::vector<EventData> events;
};


#endif