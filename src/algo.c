#include "../include/algo.h"
#include <stdio.h>
#include <string.h>

#define TWO_THIRDS_PI 2.0943951023931954923

float det(float src[3][3])
{
  return (
    src[0][0] * (src[1][1] * src[2][2] - src[2][1] * src[1][2])
    + src[0][1] * (src[2][0] * src[1][2] - src[1][0] * src[2][2])
    + src[0][2] * (src[1][0] * src[2][1] - src[2][0] * src[1][1])
  );
}

double det_to_d(float src[3][3])
{
  return (
    src[0][0] * ((double)src[1][1] * src[2][2] - (double)src[2][1] * src[1][2])
    + src[0][1] * ((double)src[2][0] * src[1][2] - (double)src[1][0] * src[2][2])
    + src[0][2] * ((double)src[1][0] * src[2][1] - (double)src[2][0] * src[1][1])
  );
}

double det_d(double src[3][3])
{
  return (
    src[0][0] * (src[1][1] * src[2][2] - src[2][1] * src[1][2])
    + src[0][1] * (src[2][0] * src[1][2] - src[1][0] * src[2][2])
    + src[0][2] * (src[1][0] * src[2][1] - src[2][0] * src[1][1])
  );
}

int eigen_values(float lambda[3], float src[3][3])
{
  double c2 = (double)src[0][0] + (double)src[1][1] + src[2][2];
  double c1 = (
    (double)src[1][1] * src[2][2] - (double)src[2][1] * src[1][2]
    + (double)src[0][0] * src[2][2] - (double)src[2][0] * src[0][2]
    + (double)src[0][0] * src[1][1] - (double)src[1][0] * src[0][1]
  );
  double c0 = det_to_d(src);

  double p = c1 - c2 * c2 / 3.;

  if (fabs(p) < 1e-9) 
  {
    lambda[0] = lambda[1] = lambda[2] = c2 / 3.;
    return 3;
  }

  double q = (9. * c2 * c1 - 2. * c2 * c2 * c2) / 27. - c0;
  double d = -4. * p * p * p - 27. * q * q;

  if (d < -1e-9)
  {
    double sqrt_minus_d = sqrt(-d / 27.);
    double u = cbrt(-q / 2. + sqrt_minus_d);
    double v = cbrt(-q / 2. - sqrt_minus_d);
    lambda[0] = (u + v) + c2 / 3.;
    return 1;
  }
  
  double acos_arg = -q / (2. * sqrt(-p * p * p / 27.));
  if (acos_arg > 1.)
    acos_arg = 1.;
  else if (acos_arg < -1.)
    acos_arg = -1.;

  double phi = acos(acos_arg) / 3., r = 2. * sqrt(-p / 3.);
  c2 /= 3.;

  if (d < 1e-9)
  {
    lambda[0] = (3. * q) / p + c2;
    lambda[1] = lambda[2] = (-3. * q) / (2. * p) + c2;
  }
  else if (d > 0)
  {
    lambda[0] = r * cos(phi) + c2;
    lambda[1] = r * cos(phi + TWO_THIRDS_PI) + c2;
    lambda[2] = r * cos(phi - TWO_THIRDS_PI) + c2;
  }

  return 3;
}

int eigen_values_to_d(double lambda[3], float src[3][3])
{
  double c2 = (double)src[0][0] + (double)src[1][1] + src[2][2];
  double c1 = (
    (double)src[1][1] * src[2][2] - (double)src[2][1] * src[1][2]
    + (double)src[0][0] * src[2][2] - (double)src[2][0] * src[0][2]
    + (double)src[0][0] * src[1][1] - (double)src[1][0] * src[0][1]
  );
  double c0 = det_to_d(src);

  double p = c1 - c2 * c2 / 3.;

  if (fabs(p) < 1e-9) 
  {
    lambda[0] = lambda[1] = lambda[2] = c2 / 3.;
    return 3;
  }

  double q = (9. * c2 * c1 - 2. * c2 * c2 * c2) / 27. - c0;
  double d = -4. * p * p * p - 27. * q * q;

  if (d < -1e-9)
  {
    double sqrt_minus_d = sqrt(-d / 27.);
    double u = cbrt(-q / 2. + sqrt_minus_d);
    double v = cbrt(-q / 2. - sqrt_minus_d);
    lambda[0] = (u + v) + c2 / 3.;
    return 1;
  }
  
  double acos_arg = -q / (2. * sqrt(-p * p * p / 27.));
  if (acos_arg > 1.)
    acos_arg = 1.;
  else if (acos_arg < -1.)
    acos_arg = -1.;

  double phi = acos(acos_arg) / 3., r = 2. * sqrt(-p / 3.);
  c2 /= 3.;

  if (d < 1e-9)
  {
    lambda[0] = (3. * q) / p + c2;
    lambda[1] = lambda[2] = (-3. * q) / (2. * p) + c2;
  }
  else if (d > 0)
  {
    lambda[0] = r * cos(phi) + c2;
    lambda[1] = r * cos(phi + TWO_THIRDS_PI) + c2;
    lambda[2] = r * cos(phi - TWO_THIRDS_PI) + c2;
  }

  return 3;
}

int eigen_values_d(double lambda[3], double src[3][3])
{
  double c2 = src[0][0] + src[1][1] + src[2][2];
  double c1 = (
    src[1][1] * src[2][2] - src[2][1] * src[1][2]
    + src[0][0] * src[2][2] - src[2][0] * src[0][2]
    + src[0][0] * src[1][1] - src[1][0] * src[0][1]
  );
  double c0 = det_d(src);

  double p = c1 - c2 * c2 / 3.;

  if (fabs(p) < 1e-9) 
  {
    lambda[0] = lambda[1] = lambda[2] = c2 / 3.;
    return 3;
  }

  double q = (9. * c2 * c1 - 2. * c2 * c2 * c2) / 27. - c0;
  double d = -4. * p * p * p - 27. * q * q;

  if (d < -1e-9)
  {
    double sqrt_minus_d = sqrt(-d / 27.);
    double u = cbrt(-q / 2. + sqrt_minus_d);
    double v = cbrt(-q / 2. - sqrt_minus_d);
    lambda[0] = (u + v) + c2 / 3.;
    return 1;
  }
  
  double acos_arg = -q / (2. * sqrt(-p * p * p / 27.));
  if (acos_arg > 1.)
    acos_arg = 1.;
  else if (acos_arg < -1.)
    acos_arg = -1.;

  double phi = acos(acos_arg) / 3., r = 2. * sqrt(-p / 3.);
  c2 /= 3.;

  if (d < 1e-9)
  {
    lambda[0] = (3. * q) / p + c2;
    lambda[1] = lambda[2] = (-3. * q) / (2. * p) + c2;
  }
  else if (d > 0)
  {
    lambda[0] = r * cos(phi) + c2;
    lambda[1] = r * cos(phi + TWO_THIRDS_PI) + c2;
    lambda[2] = r * cos(phi - TWO_THIRDS_PI) + c2;
  }

  return 3;
}

void transposed(float result[3][3], float src[3][3])
{
  float copy[3][3];
  copy[0][0] = src[0][0], copy[1][1] = src[1][1], copy[2][2] = src[2][2];
  copy[0][1] = src[1][0], copy[0][2] = src[2][0], copy[1][2] = src[2][1];
  copy[1][0] = src[0][1], copy[2][0] = src[0][2], copy[2][1] = src[1][2];
  memcpy(result, copy, sizeof(float) * 9);
}

void inverse(float result[3][3], float src[3][3])
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
  
  float inv_det = 1.f / det, copy[3][3];
  
  copy[0][0] = elem11 * inv_det, copy[1][0] = elem12 * inv_det, copy[2][0] = elem13 * inv_det;
  copy[1][1] = (src[0][0] * src[2][2] - src[2][0] * src[0][2]) * inv_det;
  copy[2][2] = (src[0][0] * src[1][1] - src[1][0] * src[0][1]) * inv_det;
  copy[0][1] = (src[2][1] * src[0][2] - src[0][1] * src[2][2]) * inv_det;
  copy[0][2] = (src[0][1] * src[1][2] - src[1][1] * src[0][2]) * inv_det;
  copy[1][2] = (src[1][0] * src[0][2] - src[0][0] * src[1][2]) * inv_det;
  copy[2][1] = (src[2][0] * src[0][1] - src[0][0] * src[2][1]) * inv_det;
  memcpy(result, copy, sizeof(float) * 9);
}

void scaled_m(float result[3][3], float m[3][3], float scale)
{
  result[0][0] = m[0][0] * scale, result[1][1] = m[1][1] * scale, result[2][2] = m[2][2] * scale;
  result[0][1] = m[0][1] * scale, result[0][2] = m[0][2] * scale, result[1][2] = m[1][2] * scale;
  result[1][0] = m[1][0] * scale, result[2][0] = m[2][0] * scale, result[2][1] = m[2][1] * scale;
}

void triangulated(float result[3][3], float src[3][3])
{
  float copy[3][3];
  memcpy(copy, src, sizeof(float) * 3);
  float denom = 1.f / copy[0][0];
  float div = -src[1][0] * denom;
  copy[1][0] = copy[2][0] = 0.f;
  copy[1][1] = src[1][1] + copy[0][1] * div, copy[1][2] = src[1][2] + copy[0][2] * div;
  div = -src[2][0] * denom;
  copy[2][1] = src[2][1] + copy[0][1] * div, copy[2][2] = src[2][2] + copy[0][2] * div;
  div = -copy[2][1] / copy[1][1];
  copy[2][1] = 0.f;
  copy[2][2] = copy[2][2] + copy[1][2] * div;
  memcpy(result, copy, sizeof(float) * 9);
}

static void triangulated_b(float result[3][3], float src[3][3], float new_b[3], float b[3])
{
  float result_copy[3][3], new_b_copy[3];
  memcpy(result_copy, src, sizeof(float) * 3);
  new_b_copy[0] = b[0];
  float denom = 1.f / result_copy[0][0];
  float div = -src[1][0] * denom;
  result_copy[1][0] = result_copy[2][0] = 0.f;
  result_copy[1][1] = src[1][1] + result_copy[0][1] * div, result_copy[1][2] = src[1][2] + result_copy[0][2] * div;
  new_b_copy[1] = b[1] + new_b_copy[0] * div;
  div = -src[2][0] * denom;
  result_copy[2][1] = src[2][1] + result_copy[0][1] * div, result_copy[2][2] = src[2][2] + result_copy[0][2] * div;
  new_b_copy[2] = b[2] + new_b_copy[0] * div;
  div = -result_copy[2][1] / result_copy[1][1];
  result_copy[2][1] = 0.f;
  result_copy[2][2] = result_copy[2][2] + result_copy[1][2] * div;
  new_b_copy[2] = new_b_copy[2] + new_b_copy[1] * div;
  memcpy(result, result_copy, sizeof(float) * 9);
  memcpy(new_b, new_b_copy, sizeof(float) * 3);
}

void outer_product(float result[3][3], float left[3], float right[3])
{
  result[0][0] = left[0] * right[0], result[1][1] = left[1] * right[1], result[2][2] = left[2] * right[2];
  result[1][0] = left[1] * right[0], result[2][0] = left[2] * right[0], result[2][1] = left[2] * right[1];
  result[0][1] = left[0] * right[1], result[0][2] = left[0] * right[2], result[1][2] = left[1] * right[2];
}

void cross_product(float result[3], float left[3], float right[3])
{
  float copy[3];
  copy[0] = left[1] * right[2] - left[2] * right[1];
  copy[1] = left[2] * right[0] - left[0] * right[2];
  copy[2] = left[0] * right[1] - left[1] * right[0];
  memcpy(result, copy, sizeof(float) * 3);
}

void cross_product_d(double result[3], double left[3], double right[3])
{
  double copy[3];
  copy[0] = left[1] * right[2] - left[2] * right[1];
  copy[1] = left[2] * right[0] - left[0] * right[2];
  copy[2] = left[0] * right[1] - left[1] * right[0];
  memcpy(result, copy, sizeof(double) * 3);
}

void mult_m_v(float result[3], float m[3][3], float v[3])
{
  float copy[3];
  copy[0] = m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2];
  copy[1] = m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2];
  copy[2] = m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2];
  memcpy(result, copy, sizeof(float) * 3);
}

void mult_m_v_d(double result[3], double m[3][3], double v[3])
{
  double copy[3];
  copy[0] = m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2];
  copy[1] = m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2];
  copy[2] = m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2];
  memcpy(result, copy, sizeof(double) * 3);
}

void mult_m_m(float result[3][3], float left[3][3], float right[3][3])
{
  float copy[3][3];
  copy[0][0] = left[0][0] * right[0][0] + left[0][1] * right[1][0] + left[0][2] * right[2][0];
  copy[0][1] = left[0][0] * right[0][1] + left[0][1] * right[1][1] + left[0][2] * right[2][1];
  copy[0][2] = left[0][0] * right[0][2] + left[0][1] * right[1][2] + left[0][2] * right[2][2];
  copy[1][0] = left[1][0] * right[0][0] + left[1][1] * right[1][0] + left[1][2] * right[2][0];
  copy[2][0] = left[2][0] * right[0][0] + left[2][1] * right[1][0] + left[2][2] * right[2][0];
  copy[1][1] = left[1][0] * right[0][1] + left[1][1] * right[1][1] + left[1][2] * right[2][1];
  copy[2][1] = left[2][0] * right[0][1] + left[2][1] * right[1][1] + left[2][2] * right[2][1];
  copy[1][2] = left[1][0] * right[0][2] + left[1][1] * right[1][2] + left[1][2] * right[2][2];
  copy[2][2] = left[2][0] * right[0][2] + left[2][1] * right[1][2] + left[2][2] * right[2][2];
  memcpy(result, copy, sizeof(float) * 9);
}

void add_m_m(float result[3][3], float left[3][3], float right[3][3])
{
  result[0][0] = left[0][0] + right[0][0], result[1][1] = left[1][1] + right[1][1], result[2][2] = left[2][2] + right[2][2];
  result[1][0] = left[1][0] + right[1][0], result[2][0] = left[2][0] + right[2][0], result[2][1] = left[2][1] + right[2][1];
  result[0][1] = left[0][1] + right[0][1], result[0][2] = left[0][2] + right[0][2], result[1][2] = left[1][2] + right[1][2];
}

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

void rotation_mat(float result[3][3], float axis[3], float angle)
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

void cross_product_mat(float result[3][3], float v[3])
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

void lu(float l[3][3], float u[3][3], float src[3][3])
{
  l[0][1] = l[0][2] = l[1][2] = u[1][0] = u[2][0] = u[2][1] = 0.f;
  l[0][0] = l[1][1] = l[2][2] = 1.f;
  u[0][0] = src[0][0], u[0][1] = src[0][1], u[0][2] = src[0][2];
  float div = 1.f / u[0][0];
  l[1][0] = src[1][0] * div, l[2][0] = src[2][0] * div;
  u[1][1] = src[1][1] - u[0][1] * l[1][0], u[1][2] = src[1][2] - l[1][0] * u[0][2];
  l[2][1] = (src[2][1] - l[2][0] * u[0][1]) / u[1][1];
  u[2][2] = src[2][2] - l[2][0] * u[0][2] - l[2][1] * u[1][2];
}

void cholesky(float l[3][3], float src[3][3])
{
  l[0][1] = l[0][2] = l[1][2] = 0.f;
  l[0][0] = sqrtf(src[0][0]);
  float div = 1.f / l[0][0];
  l[1][0] = src[0][1] * div, l[2][0] = src[0][2] * div;
  l[1][1] = sqrtf(src[1][1] - l[1][0] * l[1][0]), l[2][1] = (src[1][2] - l[2][0] * l[1][0]) / l[1][1];
  l[2][2] = sqrtf(src[2][2] - l[2][0] * l[2][0] - l[2][1] * l[2][1]);
}

static void kernel_v_d(double result[3], double src[3][3])
{
  double cross01[3], cross02[3], cross12[3];
  cross_product_d(cross01, src[0], src[1]);
  cross_product_d(cross02, src[0], src[2]);
  cross_product_d(cross12, src[1], src[2]);
  double len_sq1 = cross01[0] * cross01[0] + cross01[1] * cross01[1] + cross01[2] * cross01[2];
  double len_sq2 = cross02[0] * cross02[0] + cross02[1] * cross02[1] + cross02[2] * cross02[2];
  double len_sq3 = cross12[0] * cross12[0] + cross12[1] * cross12[1] + cross12[2] * cross12[2];
  double *best_cross = len_sq1 >= len_sq2 && len_sq1 >= len_sq3 ? cross01 : len_sq2 >= len_sq3 ? cross02 : cross12;
  memcpy(result, best_cross, sizeof(double) * 3);
}

void svd(float u[3][3], float s[3][3], float v_t[3][3], float src[3][3])
{  
  double ata[3][3], src_d[3][3];
  src_d[0][0] = src[0][0], src_d[1][1] = src[1][1], src_d[2][2] = src[2][2];
  src_d[1][0] = src[1][0], src_d[2][0] = src[2][0], src_d[2][1] = src[2][1];
  src_d[0][1] = src[0][1], src_d[0][2] = src[0][2], src_d[1][2] = src[1][2];
  
  ata[0][0] = src_d[0][0] * src_d[0][0] + src_d[1][0] * src_d[1][0] + src_d[2][0] * src_d[2][0];
  ata[0][1] = src_d[0][0] * src_d[0][1] + src_d[1][0] * src_d[1][1] + src_d[2][0] * src_d[2][1];
  ata[0][2] = src_d[0][0] * src_d[0][2] + src_d[1][0] * src_d[1][2] + src_d[2][0] * src_d[2][2];
  ata[1][1] = src_d[0][1] * src_d[0][1] + src_d[1][1] * src_d[1][1] + src_d[2][1] * src_d[2][1];
  ata[1][2] = src_d[0][1] * src_d[0][2] + src_d[1][1] * src_d[1][2] + src_d[2][1] * src_d[2][2];
  ata[2][2] = src_d[0][2] * src_d[0][2] + src_d[1][2] * src_d[1][2] + src_d[2][2] * src_d[2][2];
  ata[1][0] = ata[0][1], ata[2][0] = ata[0][2], ata[2][1] = ata[1][2];

  double sigmas[3];
  eigen_values_d(sigmas, ata);

  if (sigmas[0] < sigmas[1]) { double temp = sigmas[0]; sigmas[0] = sigmas[1]; sigmas[1] = temp; }
  if (sigmas[0] < sigmas[2]) { double temp = sigmas[0]; sigmas[0] = sigmas[2]; sigmas[2] = temp; }
  if (sigmas[1] < sigmas[2]) { double temp = sigmas[1]; sigmas[1] = sigmas[2]; sigmas[2] = temp; }

  if (fabs(sigmas[0]) < 1e-6)
  {
    memset(u, 0, sizeof(float) * 9);
    memset(s, 0, sizeof(float) * 9);
    memset(v_t, 0, sizeof(float) * 9);
    return;
  }

  double ata00 = ata[0][0], ata11 = ata[1][1], ata22 = ata[2][2];
  double v1[3], v2[3], v3[3];

  ata[0][0] -= sigmas[0], ata[1][1] -= sigmas[0], ata[2][2] -= sigmas[0];
  kernel_v_d(v1, ata);

  if (len_sq_d(v1) < 1e-6)
  {
    v1[0] = 1., v1[1] = 0., v1[2] = 0.;
    v2[0] = 0., v2[1] = 1., v2[2] = 0.;
    v3[0] = 0., v3[1] = 0., v3[2] = 1.;
  }
  else
  {
    normalise_d(v1);
  
    if (fabs(sigmas[0] - sigmas[1]) > 1e-6)
    {
      ata[0][0] = ata00 - sigmas[1], ata[1][1] = ata11 - sigmas[1], ata[2][2] = ata22 - sigmas[1];
      kernel_v_d(v2, ata);
    }
    else
    {
      double t[3] = { (fabs(v1[0]) > .5) ? 0. : 1., 0., 0. };
      cross_product_d(v2, v1, t);
    }
    normalise_d(v2);
  
    cross_product_d(v3, v1, v2);
    normalise_d(v3);
  }

  v_t[0][0] = v1[0], v_t[0][1] = v1[1], v_t[0][2] = v1[2];
  v_t[1][0] = v2[0], v_t[1][1] = v2[1], v_t[1][2] = v2[2];
  v_t[2][0] = v3[0], v_t[2][1] = v3[1], v_t[2][2] = v3[2];

  s[0][0] = sigmas[0] = sqrt(sigmas[0]);
  s[1][0] = s[2][0] = s[2][1] = s[0][1] = s[0][2] = s[1][2] = 0.f;

  double u1[3], u2[3], u3[3];
  mult_m_v_d(u1, src_d, v1);
  scale_v_d(u1, 1. / sigmas[0]);

  if (fabs(sigmas[1]) > 1e-6)
  {
    s[1][1] = sigmas[1] = sqrt(sigmas[1] > 1e-6 ? sigmas[1] : 0.);
    mult_m_v_d(u2, src_d, v2);
    scale_v_d(u2, 1. / sigmas[1]);
  
    if (fabs(sigmas[2]) > 1e-6)
    {
      s[2][2] = sigmas[2] = sqrt(sigmas[2] > 1e-6 ? sigmas[2] : 0.);
      mult_m_v_d(u3, src_d, v3);
      scale_v_d(u3, 1. / sigmas[2]);
    }
    else
    {
      cross_product_d(u3, u1, u2);
      scale_v_d(u3, 1. / len_d(u3));
      s[2][2] = sqrt(sigmas[2] > 1e-6 ? sigmas[2] : 0.);
    }
  }
  else
  {
    double t[3] = { (fabs(u1[0]) > .5) ? 0. : 1., 0., 0. };
    cross_product_d(u2, u1, t);
    scale_v_d(u2, 1. / len_d(u2));

    cross_product_d(u3, u1, u2);
    scale_v_d(u3, 1. / len_d(u3));
    s[1][1] = 0.f;
    s[2][2] = 0.f;
  }

  u[0][0] = u1[0], u[1][0] = u1[1], u[2][0] = u1[2];
  u[0][1] = u2[0], u[1][1] = u2[1], u[2][1] = u2[2];
  u[0][2] = u3[0], u[1][2] = u3[1], u[2][2] = u3[2];
}

void solve_gauss(float result[3], float src[3][3], float b[3])
{
  float tri[3][3], tri_b[3];
  triangulated_b(tri, src, tri_b, b);
  result[2] = tri_b[2] / tri[2][2];
  result[1] = (tri_b[1] - result[2] * tri[1][2]) / tri[1][1];
  result[0] = (tri_b[0] - result[1] * tri[0][1] - result[2] * tri[0][2]) / tri[0][0];
}

void solve_lu(float result[3], float l[3][3], float u[3][3], float b[3])
{
  float y[3];
  y[0] = b[0], y[1] = b[1] - l[1][0] * b[0], y[2] = b[2] - l[2][0] * b[0] - l[2][1] * y[1];
  result[2] = y[2] / u[2][2], result[1] = (y[1] - u[1][2] * result[2]) / u[1][1];
  result[0] = (y[0] - u[0][1] * result[1] - u[0][2] * result[2]) / u[0][0];
}

void solve_cholesky(float result[3], float l[3][3], float b[3])
{
  float y[3];
  y[0] = b[0] / l[0][0], y[1] = (b[1] - l[1][0] * y[0]) / l[1][1];
  y[2] = (b[2] - l[2][0] * y[0] - l[2][1] * y[1]) / l[2][2];
  result[2] = y[2] / l[2][2], result[1] = (y[1] - l[2][1] * result[2]) / l[1][1];
  result[0] = (y[0] - l[2][0] * result[2] - l[1][0] * result[1]) / l[0][0];
}
