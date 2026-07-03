#include "../include/algo.h"

float det(const float src[3][3])
{
  return (
    src[0][0] * (src[1][1] * src[2][2] - src[2][1] * src[1][2])
    + src[0][1] * (src[2][0] * src[1][2] - src[1][0] * src[2][2])
    + src[0][2] * (src[1][0] * src[2][1] - src[2][0] * src[1][1])
  );
}

double det_to_d(const float src[3][3])
{
  return (
    src[0][0] * ((double)src[1][1] * src[2][2] - (double)src[2][1] * src[1][2])
    + src[0][1] * ((double)src[2][0] * src[1][2] - (double)src[1][0] * src[2][2])
    + src[0][2] * ((double)src[1][0] * src[2][1] - (double)src[2][0] * src[1][1])
  );
}

double det_d(const double src[3][3])
{
  return (
    src[0][0] * (src[1][1] * src[2][2] - src[2][1] * src[1][2])
    + src[0][1] * (src[2][0] * src[1][2] - src[1][0] * src[2][2])
    + src[0][2] * (src[1][0] * src[2][1] - src[2][0] * src[1][1])
  );
}
