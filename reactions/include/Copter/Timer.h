#ifndef TIMER_H
#define TIMER_H

#include <ctime>


/* Timer: A simple utility class for timing operations. */

class Timer {
public:
    Timer();

    void Start();

    float CpuTimeElapsed();
    int WallTimeElapsed();

    void Report(const char* s = "");
    void ReportAndReset(const char* s = "");

    clock_t clock_start;
    time_t time_start;
};

#endif // TIMER_H
