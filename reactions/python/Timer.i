%{
#include "Timer.h"
%}

class Timer {
public:
    Timer();

    void Start();
    float CpuTimeElapsed();
    int WallTimeElapsed();
    void Report(const char* s = "");
    void ReportAndReset(const char* s = "");
};
