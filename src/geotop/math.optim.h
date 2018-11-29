#ifndef __math_optim_h_
#define __math_optim_h_

#include "config.h"

#include <cmath>

#define pow_2(a) ((a)*(a))

template <class T>
T power(const T a, const T b){
#ifdef MATH_OPTIM
    return std::exp(b*std::log(a));
//        return std::exp(b*std::log(static_cast<long double>(a)));
#else
    return std::pow(a,b);
#endif
}

#endif
