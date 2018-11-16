#ifndef __math_optim_h_
#define __math_optim_h_

#include <cmath>

#define pow_2(a) ((a)*(a))
#define pow_3(a) ((a)*(pow_2(a)))
#define pow_4(a) ((pow_2(a))*(pow_2(a)))
#define pow_6(a) ((pow_3(a))*(pow_3(a)))

template <class T>
T power(const T a, const T b){
#ifdef MATH_OPTIM
  return std::exp(b*std::log(a));
#else
  return std::pow(a,b);
#endif
}

#endif 
