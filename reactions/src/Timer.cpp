#if HAVE_CONFIG_H
# include <config.h>
#endif

#include "Timer.h"
#include "Common.h"

Timer::Timer() {
    Start();
}

void Timer::Start() {
    clock_start = clock();
    time_start = time(NULL);
}

float Timer::CpuTimeElapsed() {
    return (float)(clock() - clock_start)/(float)CLOCKS_PER_SEC;
}

int Timer::WallTimeElapsed() {
    return time(NULL) - time_start;
}

void Timer::Report(const char* s) {
    if(s && *s != '\0')
        info("%s: time elapsed = %d seconds\n", s, WallTimeElapsed());
    else
        info("time elapsed = %d seconds\n", WallTimeElapsed());
}

void Timer::ReportAndReset(const char* s) {
    Report(s);
    Start();
}
