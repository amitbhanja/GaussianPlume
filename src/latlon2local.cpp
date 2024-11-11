//
// File: latlon2local.cpp
//
// MATLAB Coder version            : 24.2
// C/C++ source code generated on  : 29-Oct-2024 13:42:30
//

// Include Files
#include "latlon2local.h"
#include "cosd.h"
#include "sind.h"
#include <cmath>
#include <emmintrin.h>

// Function Definitions
//
// Arguments    : const double lat[100]
//                const double lon[100]
//                const double origin[3]
//                double xEast[100]
//                double yNorth[100]
// Return Type  : void
//
namespace coder {
void latlon2local(const double lat[100], const double lon[100],
                  const double origin[3], double xEast[100], double yNorth[100])
{
  __m128d r;
  __m128d r1;
  __m128d r2;
  __m128d r3;
  double c2[100];
  double s2[100];
  double w2[100];
  double b_p1;
  double c1;
  double cosLambda;
  double cosPhi;
  double d;
  double d1;
  double p1;
  double q1;
  double s1;
  double sinLambda;
  double sinPhi;
  s1 = origin[0];
  b_sind(s1);
  c1 = origin[0];
  b_cosd(c1);
  d = origin[1];
  b_cosd(d);
  p1 = c1 * d;
  for (int i{0}; i < 100; i++) {
    d = lat[i];
    s2[i] = d;
    c2[i] = d;
    yNorth[i] = lon[i];
  }
  b_sind(s2);
  b_cosd(c2);
  b_cosd(yNorth);
  d = origin[1];
  b_sind(d);
  c1 *= d;
  for (int i{0}; i <= 98; i += 2) {
    r = _mm_loadu_pd(&c2[i]);
    r1 = _mm_loadu_pd(&yNorth[i]);
    _mm_storeu_pd(&yNorth[i], _mm_mul_pd(r, r1));
    _mm_storeu_pd(&w2[i], _mm_loadu_pd(&lon[i]));
  }
  double w1;
  b_sind(w2);
  w1 = 1.0 / std::sqrt(1.0 - 0.0066943799901413165 * (s1 * s1));
  b_p1 = p1 * w1;
  d = 10.0 * p1;
  q1 = c1 * w1;
  d1 = 10.0 * c1;
  cosPhi = origin[0];
  b_cosd(cosPhi);
  sinPhi = origin[0];
  b_sind(sinPhi);
  cosLambda = origin[1];
  b_cosd(cosLambda);
  sinLambda = origin[1];
  b_sind(sinLambda);
  p1 = s1 * w1;
  c1 = 10.0 * s1;
  r = _mm_set1_pd(1.0);
  r1 = _mm_set1_pd(6.378137E+6);
  r2 = _mm_set1_pd(10.0);
  r3 = _mm_set1_pd(cosLambda);
  for (int i{0}; i <= 98; i += 2) {
    __m128d r4;
    __m128d r5;
    __m128d r6;
    __m128d r7;
    r4 = _mm_loadu_pd(&c2[i]);
    r5 = _mm_loadu_pd(&w2[i]);
    r4 = _mm_mul_pd(r4, r5);
    r5 = _mm_loadu_pd(&s2[i]);
    r6 = _mm_div_pd(r, _mm_sqrt_pd(_mm_sub_pd(
                           r, _mm_mul_pd(_mm_set1_pd(0.0066943799901413165),
                                         _mm_mul_pd(r5, r5)))));
    _mm_storeu_pd(&w2[i], r6);
    r7 = _mm_loadu_pd(&yNorth[i]);
    r7 = _mm_add_pd(
        _mm_mul_pd(r1, _mm_sub_pd(_mm_mul_pd(r7, r6), _mm_set1_pd(b_p1))),
        _mm_sub_pd(_mm_mul_pd(r2, r7), _mm_set1_pd(d)));
    r4 = _mm_add_pd(
        _mm_mul_pd(r1, _mm_sub_pd(_mm_mul_pd(r4, r6), _mm_set1_pd(q1))),
        _mm_sub_pd(_mm_mul_pd(r2, r4), _mm_set1_pd(d1)));
    _mm_storeu_pd(&c2[i], r4);
    _mm_storeu_pd(&xEast[i], _mm_add_pd(_mm_mul_pd(_mm_set1_pd(-sinLambda), r7),
                                        _mm_mul_pd(r3, r4)));
    _mm_storeu_pd(
        &yNorth[i],
        _mm_add_pd(
            _mm_mul_pd(_mm_set1_pd(-sinPhi),
                       _mm_add_pd(_mm_mul_pd(r3, r7),
                                  _mm_mul_pd(_mm_set1_pd(sinLambda), r4))),
            _mm_mul_pd(
                _mm_set1_pd(cosPhi),
                _mm_add_pd(
                    _mm_mul_pd(_mm_set1_pd(6.3354393272928195E+6),
                               _mm_sub_pd(_mm_mul_pd(r5, r6), _mm_set1_pd(p1))),
                    _mm_sub_pd(_mm_mul_pd(r2, r5), _mm_set1_pd(c1))))));
  }
}

} // namespace coder

//
// File trailer for latlon2local.cpp
//
// [EOF]
//
