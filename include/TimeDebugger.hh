#ifndef TimeDebugger_h
#define TimeDebugger_h 1

class TimeDebugger
{
    public:
        inline static double sampleEdepTime = 0;
        inline static double getTransportValuesTime = 0;
        inline static double processHitTime = 0;
        inline static double navigatorTime = 0;
        inline static double inBetweenStuffTime = 0;
        inline static double sdStuffTime = 0;
        inline static double A2SDProcessHitsTime = 0;
        inline static double A2SDSetup = 0;
        inline static double A2SDIfTrue = 0;
        inline static double A2SDIfFalse = 0;
        inline static double falseBlockEdepQdepTime = 0;
        inline static double falseBlockTPCBlockTime = 0;
        inline static double line1 = 0;
        inline static double line2 = 0;
        inline static double line3 = 0;
};

#endif
