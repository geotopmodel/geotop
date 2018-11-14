#include <cmath>

double power(double a, double b){
#ifdef MATH_OPTIM
  return exp(b*log(a));
#else
  return std::pow(a,b);
#endif
}
