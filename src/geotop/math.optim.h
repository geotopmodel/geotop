#ifndef __math_optim_h_
#define __math_optim_h_

#include "config.h"

#include <cmath>
#include <fstream>

#define pow_2(a) ((a)*(a))

template <class T>
T power(const T a, const T b){
#ifdef MATH_OPTIM
//    std::ofstream myfile{"analyze_power.txt", std::ios_base::app};
//    myfile << "(" << a << "," << b << ")\t\t" << std::endl;
    return std::exp(b*std::log(a));
#else
//    std::ofstream file{"analyze_pow.txt", std::ios_base::app};
//    file << "(" << a << "," << b << ")\t\t" << std::endl;
    return std::pow(a,b);
#endif
}

#endif
