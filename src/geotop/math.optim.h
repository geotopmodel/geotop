#ifndef __math_optim_h_
#define __math_optim_h_

#include <cmath>

template <class T>
T pow_2(const T a){
#ifdef MATH_OPTIM
return a*a;
#else
return std::pow(a,2);
#endif
}
//#define pow_2(a) a*a

template <class T>
T power(const T a, const T b){
#ifdef MATH_OPTIM
  return std::exp(b*std::log(a));
#else
  return std::pow(a,b);
#endif
}

#endif 
