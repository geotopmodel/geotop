//
// Created by elisa on 29/10/18.
//
#ifndef GEOTOP_TIMER_H
#define GEOTOP_TIMER_H

#include <chrono>
#include <iostream>
#include <map>

using namespace std::chrono;

/**
 * Utility structure to collect all the measurements
 */
struct ClockMeasurements {
  unsigned int number_of_calls{};
  double elapsed_time{};
};

/**
 *  Simple class timer. By default, a summary of the time a function
 *  has been called and the cpu time spent on it is reported.
 *  TODO: The class is NOT thread safe.
 *
 * The design of output produced by the print_summary() function is
 * inspired from the one of the deal.II library
 * (github.com/dealii/dealii).
 */
class Timer {
  void print_summary();
  high_resolution_clock::time_point t_start;
  std::map<std::string, ClockMeasurements> times;
 public:
  Timer() : t_start{high_resolution_clock::now()} {}
  ~Timer() {
#   ifndef MUTE_GEOTIMER
      print_summary();
#   endif
  }
  class ScopedTimer;
};

/** the default timer **/
extern Timer geotimer;


/** 
 * utility class. The constructor takes a string and a
 * timer. Explointing the RAII pattern, the number of function calls
 * and the elapsed time passed inside the scope labeled with a chosen
 * string are updated automatically.
 */
class Timer::ScopedTimer {
 public:
  const high_resolution_clock::time_point t;
  const std::string _s;
  Timer& _timer;
  ScopedTimer(const std::string& s, Timer& t = geotimer)
    : t{high_resolution_clock::now()}, _s{s}, _timer(t) {
    _timer.times[s].number_of_calls += 1;
  }

  ~ScopedTimer() {
    auto t_end = high_resolution_clock::now();
    _timer.times[_s].elapsed_time +=
        duration_cast<duration<double>>(t_end - t).count();
  }
};

#ifdef MUTE_GEOTIMER
#  define GEOTIMER_PREFIX(dummy)
#else
#  define GEOTIMER_PREFIX(string)                                              \
    Timer::ScopedTimer __geotimer_prefix__ { string }
#endif

#endif  // GEOTOP_TIMER_H
