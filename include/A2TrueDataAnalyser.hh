#ifndef A2TrueDataAnalyser_h
#define A2TrueDataAnalyser_h

#include <vector>
#include <unordered_map>

#include "A2TrueData.hh"



class A2TrueDataAnalyser
{
    public:
        explicit A2TrueDataAnalyser(RunData&& runData);    
        ~A2TrueDataAnalyser();
        void VisualizeTree(bool all = 0) const;
        void MakeEKinHists(std::size_t chosenLayer = allLayers) const;  
        // default: all layers are selected, 0:only primaries, 1:only secondaries, ...
        void MakeEdepEKinHists(std::size_t chosenLayer = allLayers) const; 
        // default: all layers are selected, 0:only primaries, 1:only secondaries, ...
        void MakePrimaryTrackInfoHists(std::size_t chosenLayer = allLayers) const;
        // default: all layers are selected, 0:only primaries, 1:only secondaries, ...
        void StepLengthPlots(std::size_t chosenLayer = allLayers) const;
        // default: all layers are selected, 0:only primaries, 1:only secondaries, ...
        void StoppingPower(std::size_t chosenLayer = allLayers) const;
    
    private:
        static constexpr std::size_t allLayers = std::numeric_limits<std::size_t>::max();      
        //used to select all layers if a function requires the selction of a specific layer
        static constexpr std::size_t noDepthLimit = std::numeric_limits<std::size_t>::max();   
        //used to to select no premature stopping of recursion if a function requires a max recursion depth
        const RunData fRunData;   //copy instead of reference for later seperation of sim and ana
        std::unordered_map<int, std::size_t> MakeTrackLookUpMap(const EventData& event) const;
        std::unordered_map<int, std::vector<int>> MakeChildTrackMap(const EventData& event) const;

        double GetEdepTrack(const TrackData& track) const;

        void WriteTreeSegment ( int trackID,
                                const std::vector<TrackData>& tracks,
                                const std::unordered_map<int,std::size_t>& trackLookUpMap,
                                const std::unordered_map<int, std::vector<int>>& childTrackMap,
                                std::ofstream& out,
                                std::size_t depth = 0 ) const; //used in VisualizeTree()
        
        //needs to be implemented in the header due to the template Datatype
        template<typename Func>
        void ParseTrackTree( int trackID,
                            const std::unordered_map<int, std::vector<int>>& childTrackMap,
                            const Func& func,
                            std::size_t maxDepth = noDepthLimit,   //terminationDepth default: no temernationDepth
                            std::size_t depth = 0 ) const          //currentDepth
                                                                 
            {
                if ( depth<=maxDepth )
                {
                    func(trackID,depth);

                    auto it = childTrackMap.find(trackID);
                    if (it != childTrackMap.end())
                    {
                        for (const int& childID : it->second)
                        {
                            ParseTrackTree ( childID,
                                            childTrackMap,
                                            func,
                                            maxDepth,
                                            depth+1 );
                        }
                    }
                }  
            }
};

#endif