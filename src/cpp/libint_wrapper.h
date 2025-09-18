#pragma once
#include <cstddef>

/* Number of cartesian function for a given angular momentum */
inline int nint(int am) {
  return (am+1)*(am+2)/2;
}

extern "C" {
    void dboysfun12(const double* x, double* vals);
}