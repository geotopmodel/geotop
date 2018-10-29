//
// Created by elisa on 29/10/18.
//
#ifndef GEOTOP_TIMER_H
#define GEOTOP_TIMER_H

#include <iostream>
#include <chrono>
#include <map>

using namespace std::chrono;

class Timer {
public:
    std::map<std::string, std::pair<int,duration<double>>> times;
    ~ Timer() {
        for (const auto& x:times){
            std::cerr << " --------------------------------------------------- " << std::endl;
            std::cerr << x.first << "\t\t" /** function name **/
                      << x.second.first << "\t\t" /** n° of times the function is called **/
                      << x.second.second.count() << std::endl; /** time passed inside the function **/
        }
    }
    class ScopedTimer;
};

/** the default timer **/
extern Timer geotimer;

class Timer::ScopedTimer {
public:
    high_resolution_clock::time_point t;
    const std::string _s;
    ScopedTimer(const std::string &s): _s{s}{
        geotimer.times[s].first += 1;
        t = high_resolution_clock::now();
    }

    ~ScopedTimer(){
        auto t_end = high_resolution_clock::now();
        geotimer.times[_s].second += duration_cast<duration<double>>(t_end-t);
    }
};


#ifdef MUTE_GEOTIMER
#define GEOTIMER_PREFIX(dummy)
#else
# define GEOTIMER_PREFIX(string) Timer::ScopedTimer __geotimer_prefix__ {string}
#endif

#endif // GEOTOP_TIMER_H
