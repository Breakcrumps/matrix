#include "../include/algo.h"
#include <string.h>

static void triangulated_b(float result[restrict 3][3], const float src[restrict 3][3], float new_b[restrict 3], const float b[restrict 3])
{
  memcpy(result, src, sizeof(float) * 3);
  new_b[0] = b[0];
  float denom = 1.f / result[0][0];
  float div = -src[1][0] * denom;
  result[1][0] = result[2][0] = 0.f;
  result[1][1] = src[1][1] + result[0][1] * div, result[1][2] = src[1][2] + result[0][2] * div;
  new_b[1] = b[1] + new_b[0] * div;
  div = -src[2][0] * denom;
  result[2][1] = src[2][1] + result[0][1] * div, result[2][2] = src[2][2] + result[0][2] * div;
  new_b[2] = b[2] + new_b[0] * div;
  div = -result[2][1] / result[1][1];
  result[2][1] = 0.f, result[2][2] += result[1][2] * div;
  new_b[2] += new_b[1] * div;
}

void solve_gauss(float result[3], const float src[3][3], const float b[3])
{
  float tri[3][3], tri_b[3];
  triangulated_b(tri, src, tri_b, b);
  result[2] = tri_b[2] / tri[2][2];
  result[1] = (tri_b[1] - result[2] * tri[1][2]) / tri[1][1];
  result[0] = (tri_b[0] - result[1] * tri[0][1] - result[2] * tri[0][2]) / tri[0][0];
}

void solve_lu(float result[3], const float l[3][3], const float u[3][3], const float b[3])
{
  float y[3];
  y[0] = b[0], y[1] = b[1] - l[1][0] * b[0], y[2] = b[2] - l[2][0] * b[0] - l[2][1] * y[1];
  result[2] = y[2] / u[2][2], result[1] = (y[1] - u[1][2] * result[2]) / u[1][1];
  result[0] = (y[0] - u[0][1] * result[1] - u[0][2] * result[2]) / u[0][0];
}

void solve_cholesky(float result[3], const float l[3][3], const float b[3])
{
  float y[3];
  y[0] = b[0] / l[0][0], y[1] = (b[1] - l[1][0] * y[0]) / l[1][1];
  y[2] = (b[2] - l[2][0] * y[0] - l[2][1] * y[1]) / l[2][2];
  result[2] = y[2] / l[2][2], result[1] = (y[1] - l[2][1] * result[2]) / l[1][1];
  result[0] = (y[0] - l[2][0] * result[2] - l[1][0] * result[1]) / l[0][0];
}
