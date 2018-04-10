/**
 * @brief Soil State Variables implementation
 * @date September 2014
 */

#include "soilstatevar.h"

SoilState::SoilState(size_t total_pixel, size_t layers):
  P{layers + 1, total_pixel + 1, 0.},
  thi{layers + 1, total_pixel + 1, 0.},
  T{layers + 1, total_pixel + 1}
{}
