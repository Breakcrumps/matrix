#include "../include/algo.h"

void rotation_mat_x(float result[3][3], float angle)
{
  float cos = cosf(angle), sin = sinf(angle);
  result[0][0] = 1.f;
  result[0][1] = result[0][2] = result[1][0] = result[2][0] = 0.f;
  result[1][1] = result[2][2] = cos;
  result[1][2] = -(result[2][1] = sin);
}

void rotation_mat_y(float result[3][3], float angle)
{
  float cos = cosf(angle), sin = sinf(angle);
  result[1][1] = 1.f;
  result[1][0] = result[1][2] = result[0][1] = result[2][1] = 0.f;
  result[0][0] = result[2][2] = cos;
  result[2][0] = -(result[0][2] = sin);
}

void rotation_mat_z(float result[3][3], float angle)
{
  float cos = cosf(angle), sin = sinf(angle);
  result[2][2] = 1.f;
  result[2][1] = result[2][0] = result[1][2] = result[0][2] = 0.f;
  result[0][0] = result[1][1] = cos;
  result[0][1] = -(result[1][0] = sin);
}

void rotation_mat(float result[3][3], const float axis[3], float angle)
{
  float outer[3][3], cos = cosf(angle), sin = sinf(angle), one_min_cos= 1.f - cos;
  outer[0][0] = axis[0] * axis[0], outer[1][1] = axis[1] * axis[1], outer[2][2] = axis[2] * axis[2];
  outer[1][0] = axis[1] * axis[0], outer[2][0] = axis[2] * axis[0], outer[2][1] = axis[2] * axis[1];
  outer[0][1] = outer[1][0], outer[0][2] = outer[2][0], outer[1][2] = outer[2][1];
  result[0][0] = cos + outer[0][0] * one_min_cos, result[1][1] = cos + outer[1][1] * one_min_cos;
  result[2][2] = cos + outer[2][2] * one_min_cos;
  result[1][0] = outer[1][0] * one_min_cos + axis[2] * sin, result[2][0] = outer[2][0] * one_min_cos - axis[1] * sin;
  result[2][1] = outer[2][1] * one_min_cos + axis[0] * sin;
  result[0][1] = outer[0][1] * one_min_cos - axis[2] * sin, result[0][2] = outer[0][2] * one_min_cos + axis[1] * sin;
  result[1][2] = outer[1][2] * one_min_cos - axis[0] * sin;
}

void cross_product_mat(float result[3][3], const float v[3])
{
  result[0][0] = result[1][1] = result[2][2] = 0.f;
  result[0][1] = -(result[1][0] = v[2]);
  result[2][0] = -(result[0][2] = v[1]);
  result[1][2] = -(result[2][1] = v[0]);
}

void euler_to_mat_zyx(float res[3][3], float pitch, float yaw, float roll)
{
  float cp = cosf(pitch), sp = sinf(pitch);
  float cy = cosf(yaw), sy = sinf(yaw);
  float cr = cosf(roll), sr = sinf(roll);
  res[0][0] = cy * cr, res[1][0] = cy * sr, res[2][0] = -sy, res[2][1] = sp * cy, res[2][2] = cp * cy;
  res[0][1] = sp * sy * cr - cp * sr, res[0][2] = cp * sy * cr + sp * sr;
  res[1][1] = sp * sy * sr + cp * cr, res[1][2] = cp * sy * sr - sp * cr;
}
