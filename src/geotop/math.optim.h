#ifndef _GEOTOP_MATH_OPTIM_H
#define _GEOTOP_MATH_OPTIM_H

#include "config.h"

#include <cmath>

#define pow_2(a) ((a)*(a))

template <class T>
T power(const T a, const T b){
#ifdef MATH_OPTIM
    return std::exp(b*std::log(a));
#else
    return std::pow(a,b);
#endif
}

#endif
