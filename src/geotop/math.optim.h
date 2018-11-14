#include <cmath>

template <class T>
T power(const T a, const T b){
#ifdef MATH_OPTIM
    return std::exp(b*std::log(a));
#else
    return std::pow(a,b);
#endif
}