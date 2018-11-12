//
// Created by elisa on 29/10/18.
//

#include "timer.h"
#include <algorithm>
#include <iomanip>

// the default timer
Timer geotimer;

void Timer::print_summary() {
  if (times.size() == 0 ) return; // do not print empty table
  // compute the total elapsed time in seconds
  const auto t_end = high_resolution_clock::now();
  const auto t_tot = duration_cast<duration<double>>(t_end - t_start).count();

  // get the maximum length among the section (function) names
  unsigned int max_length{0};
  for (const auto& x : times) {
    max_length = std::max(max_length, (unsigned int)x.first.length());
  }

  // 32 is the default width until | character
  max_length = std::max(max_length + 1, 32u);
  const std::string extra_dash = std::string(max_length - 32, '-');
  const std::string extra_space = std::string(max_length - 32, ' ');

  // generate a nice table
  std::cerr << "\n\n"
            << "+---------------------------------------------" << extra_dash
            << "+------------"
            << "+------------+\n"
            << "| Total CPU time elapsed since start          " << extra_space
            << "|";
  std::cerr << std::setw(10) << std::setprecision(3) << std::right;
  std::cerr << t_tot << "s |            |\n";
  std::cerr << "|                                             " << extra_space
            << "|            "
            << "|            |\n";
  std::cerr << "| Section                         " << extra_space
            << "| no. calls |";
  std::cerr << std::setw(10);
  std::cerr << std::setprecision(3);
  std::cerr << "  CPU time "
            << " | % of total |\n";
  std::cerr << "+---------------------------------" << extra_dash
            << "+-----------+------------"
            << "+------------+";
  for (const auto& i : times) {
    std::string name_out = i.first;

    // resize the array so that it is always of the same size
    unsigned int pos_non_space = name_out.find_first_not_of(' ');
    name_out.erase(0, pos_non_space);
    name_out.resize(max_length, ' ');
    std::cerr << std::endl;
    std::cerr << "| " << name_out;
    std::cerr << "| ";
    std::cerr << std::setw(9);
    std::cerr << i.second.number_of_calls << " |";
    std::cerr << std::setw(10);
    std::cerr << std::fixed << std::setprecision(2);
    std::cerr << i.second.elapsed_time << "s |";
    std::cerr << std::setw(10);

    // do not print too small percentage
    const double fraction = i.second.elapsed_time / t_tot;
    if (fraction > 0.001) {
      std::cerr << std::fixed << std::setprecision(1);
      std::cerr << fraction * 100.;
    } else
      std::cerr << 0.0;

    std::cerr << "% |";
  }
  std::cerr << std::endl
            << "+---------------------------------" << extra_dash
            << "+-----------+"
            << "------------+------------+\n"
            << std::endl;
}
