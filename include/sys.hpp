#ifndef SYS_HPP
#define SYS_HPP

#include <stdio.h>
#include <sys/resource.h>
#include <sys/time.h>

double cputime() {
    rusage r;
    getrusage(RUSAGE_SELF, &r);
    return r.ru_utime.tv_sec + r.ru_stime.tv_sec + 1e-6 * (r.ru_utime.tv_usec + r.ru_stime.tv_usec);
}

long peakrss() {
    rusage r;
    getrusage(RUSAGE_SELF, &r);
#ifdef __linux__
    return r.ru_maxrss * 1024;
#else
    return r.ru_maxrss;
#endif
}

double realtime() {
    timeval tp;
    gettimeofday(&tp, NULL);
    return tp.tv_sec + tp.tv_usec * 1e-6;
}

#endif
