#include "../include/algo.h"
#include <string.h>

void outer_product(float result[3][3], const float left[3], const float right[3])
{
  result[0][0] = left[0] * right[0], result[1][1] = left[1] * right[1], result[2][2] = left[2] * right[2];
  result[1][0] = left[1] * right[0], result[2][0] = left[2] * right[0], result[2][1] = left[2] * right[1];
  result[0][1] = left[0] * right[1], result[0][2] = left[0] * right[2], result[1][2] = left[1] * right[2];
}

void cross_product(float result[3], const float left[3], const float right[3])
{
  result[0] = left[1] * right[2] - left[2] * right[1];
  result[1] = left[2] * right[0] - left[0] * right[2];
  result[2] = left[0] * right[1] - left[1] * right[0];
}

void cross_product_d(double result[3], const double left[3], const double right[3])
{
  result[0] = left[1] * right[2] - left[2] * right[1];
  result[1] = left[2] * right[0] - left[0] * right[2];
  result[2] = left[0] * right[1] - left[1] * right[0];
}

void mult_m_v(float result[3], const float m[3][3], const float v[3])
{
  result[0] = m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2];
  result[1] = m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2];
  result[2] = m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2];
}

void mult_m_v_d(double result[3], const double m[3][3], const double v[3])
{
  result[0] = m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2];
  result[1] = m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2];
  result[2] = m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2];
}

void mult_m_m_3x3(float result[restrict 3][3], const float left[restrict 3][3], float const right[restrict 3][3])
{
  result[0][0] = left[0][0] * right[0][0] + left[0][1] * right[1][0] + left[0][2] * right[2][0];
  result[0][1] = left[0][0] * right[0][1] + left[0][1] * right[1][1] + left[0][2] * right[2][1];
  result[0][2] = left[0][0] * right[0][2] + left[0][1] * right[1][2] + left[0][2] * right[2][2];
  result[1][0] = left[1][0] * right[0][0] + left[1][1] * right[1][0] + left[1][2] * right[2][0];
  result[2][0] = left[2][0] * right[0][0] + left[2][1] * right[1][0] + left[2][2] * right[2][0];
  result[1][1] = left[1][0] * right[0][1] + left[1][1] * right[1][1] + left[1][2] * right[2][1];
  result[2][1] = left[2][0] * right[0][1] + left[2][1] * right[1][1] + left[2][2] * right[2][1];
  result[1][2] = left[1][0] * right[0][2] + left[1][1] * right[1][2] + left[1][2] * right[2][2];
  result[2][2] = left[2][0] * right[0][2] + left[2][1] * right[1][2] + left[2][2] * right[2][2];
}

void mult_m_m(
  int left_rows, int left_cols, int right_cols,
  float result[restrict left_rows][right_cols],
  const float left[restrict left_rows][left_cols],
  const float right[restrict left_cols][right_cols]
)
{
  memset(result, 0, left_rows * right_cols * sizeof(float));

  for (int i = 0; i < left_rows; i++)
  {
    for (int k = 0; k < left_cols; k++)
    {
      float r_left = left[i][k];
      for (int j = 0; j < right_cols; j++)
      {
        result[i][j] += r_left * right[k][j];
      }
    }
  }
}

void add_m_m(float result[restrict 3][3], const float left[restrict 3][3], const float right[restrict 3][3])
{
  result[0][0] = left[0][0] + right[0][0], result[1][1] = left[1][1] + right[1][1], result[2][2] = left[2][2] + right[2][2];
  result[1][0] = left[1][0] + right[1][0], result[2][0] = left[2][0] + right[2][0], result[2][1] = left[2][1] + right[2][1];
  result[0][1] = left[0][1] + right[0][1], result[0][2] = left[0][2] + right[0][2], result[1][2] = left[1][2] + right[1][2];
}
