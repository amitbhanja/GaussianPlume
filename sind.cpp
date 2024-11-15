//
// File: sind.cpp
//
// MATLAB Coder version            : 24.2
// C/C++ source code generated on  : 29-Oct-2024 13:42:30
//

// Include Files
#include "sind.h"
#include <cmath>

// Function Definitions
//
// Arguments    : double x[100]
// Return Type  : void
//
namespace coder {
void b_sind(double x[100])
{
  for (int k{0}; k < 100; k++) {
    double b_x;
    b_x = x[k];
    if (std::isinf(b_x) || std::isnan(b_x)) {
      x[k] = rtNaN;
    } else {
      double absx;
      signed char n;
      b_x = std::fmod(b_x, 360.0);
      absx = std::abs(b_x);
      if (absx > 180.0) {
        if (b_x > 0.0) {
          b_x -= 360.0;
        } else {
          b_x += 360.0;
        }
        absx = std::abs(b_x);
      }
      if (absx <= 45.0) {
        b_x *= 0.017453292519943295;
        n = 0;
      } else if (absx <= 135.0) {
        if (b_x > 0.0) {
          b_x = 0.017453292519943295 * (b_x - 90.0);
          n = 1;
        } else {
          b_x = 0.017453292519943295 * (b_x + 90.0);
          n = -1;
        }
      } else if (b_x > 0.0) {
        b_x = 0.017453292519943295 * (b_x - 180.0);
        n = 2;
      } else {
        b_x = 0.017453292519943295 * (b_x + 180.0);
        n = -2;
      }
      if (n == 0) {
        x[k] = std::sin(b_x);
      } else if (n == 1) {
        x[k] = std::cos(b_x);
      } else if (n == -1) {
        x[k] = -std::cos(b_x);
      } else {
        x[k] = -std::sin(b_x);
      }
    }
  }
}

//
// Arguments    : double &x
// Return Type  : void
//
void b_sind(double &x)
{
  if (std::isinf(x) || std::isnan(x)) {
    x = rtNaN;
  } else {
    double absx;
    signed char n;
    if (std::isnan(x) || std::isinf(x)) {
      x = rtNaN;
    } else {
      x = std::fmod(x, 360.0);
    }
    absx = std::abs(x);
    if (absx > 180.0) {
      if (x > 0.0) {
        x -= 360.0;
      } else {
        x += 360.0;
      }
      absx = std::abs(x);
    }
    if (absx <= 45.0) {
      x *= 0.017453292519943295;
      n = 0;
    } else if (absx <= 135.0) {
      if (x > 0.0) {
        x = 0.017453292519943295 * (x - 90.0);
        n = 1;
      } else {
        x = 0.017453292519943295 * (x + 90.0);
        n = -1;
      }
    } else if (x > 0.0) {
      x = 0.017453292519943295 * (x - 180.0);
      n = 2;
    } else {
      x = 0.017453292519943295 * (x + 180.0);
      n = -2;
    }
    if (n == 0) {
      x = std::sin(x);
    } else if (n == 1) {
      x = std::cos(x);
    } else if (n == -1) {
      x = -std::cos(x);
    } else {
      x = -std::sin(x);
    }
  }
}

} // namespace coder

//
// File trailer for sind.cpp
//
// [EOF]
//
