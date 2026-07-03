#include "../include/algo.h"
#include <string.h>
#include <stdio.h>

void transposed(float result[3][3], const float src[3][3])
{
  result[0][0] = src[0][0], result[1][1] = src[1][1], result[2][2] = src[2][2];
  result[0][1] = src[1][0], result[0][2] = src[2][0], result[1][2] = src[2][1];
  result[1][0] = src[0][1], result[2][0] = src[0][2], result[2][1] = src[1][2];
}

void inverse(float result[3][3], const float src[3][3])
{
  float elem11 = src[1][1] * src[2][2] - src[2][1] * src[1][2];
  float elem12 = src[2][0] * src[1][2] - src[1][0] * src[2][2];
  float elem13 = src[1][0] * src[2][1] - src[2][0] * src[1][1];
  float det = src[0][0] * elem11 + src[0][1] * elem12 + src[0][2] * elem13;
  
  if (!det)
  {
    puts("The matrix was singular! Result matrix inchanged.");
    return;
  }
  
  float inv_det = 1.f / det;
  
  result[0][0] = elem11 * inv_det, result[1][0] = elem12 * inv_det, result[2][0] = elem13 * inv_det;
  result[1][1] = (src[0][0] * src[2][2] - src[2][0] * src[0][2]) * inv_det;
  result[2][2] = (src[0][0] * src[1][1] - src[1][0] * src[0][1]) * inv_det;
  result[0][1] = (src[2][1] * src[0][2] - src[0][1] * src[2][2]) * inv_det;
  result[0][2] = (src[0][1] * src[1][2] - src[1][1] * src[0][2]) * inv_det;
  result[1][2] = (src[1][0] * src[0][2] - src[0][0] * src[1][2]) * inv_det;
  result[2][1] = (src[2][0] * src[0][1] - src[0][0] * src[2][1]) * inv_det;
}

void scaled_m(float result[3][3], const float m[3][3], float scale)
{
  result[0][0] = m[0][0] * scale, result[1][1] = m[1][1] * scale, result[2][2] = m[2][2] * scale;
  result[0][1] = m[0][1] * scale, result[0][2] = m[0][2] * scale, result[1][2] = m[1][2] * scale;
  result[1][0] = m[1][0] * scale, result[2][0] = m[2][0] * scale, result[2][1] = m[2][1] * scale;
}

void triangulated(float result[3][3], const float src[3][3])
{
  memcpy(result, src, sizeof(float) * 3);
  float denom = 1.f / result[0][0];
  float div = -src[1][0] * denom;
  result[1][0] = result[2][0] = 0.f;
  result[1][1] = src[1][1] + result[0][1] * div, result[1][2] = src[1][2] + result[0][2] * div;
  div = -src[2][0] * denom;
  result[2][1] = src[2][1] + result[0][1] * div, result[2][2] = src[2][2] + result[0][2] * div;
  div = -result[2][1] / result[1][1];
  result[2][1] = 0.f, result[2][2] += result[1][2] * div;
}
