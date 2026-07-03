#include "../include/algo.h"
#include <stdio.h>
#include <string.h>

#define TWO_THIRDS_PI 2.0943951023931954923
#define EPS_D 1e-14
#define SQRT3 1.73205080756887729
#define ONE_THIRD 0.333333333333333314829616256247390992939472198486328125

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

int eigen_values(float lambda[3], const float src[3][3])
{
  double c2 = (double)src[0][0] + (double)src[1][1] + src[2][2];
  double c1 = (
    (double)src[1][1] * src[2][2] - (double)src[2][1] * src[1][2]
    + (double)src[0][0] * src[2][2] - (double)src[2][0] * src[0][2]
    + (double)src[0][0] * src[1][1]  - (double)src[1][0] * src[0][1]
  );
  double c0 = det_to_d(src);

  double p = c1 - c2 * c2 * ONE_THIRD;
  double q = (9. * c2 * c1 - 2. * c2 * c2 * c2) / 27. - c0;

  if (fabs(p) < EPS_D && fabs(q) < EPS_D) 
  {
    lambda[0] = lambda[1] = lambda[2] = c2 * ONE_THIRD;
    return 3;
  }

  double d = -4. * p * p * p - 27. * q * q;

  if (d < -EPS_D)
  {
    double sqrt_minus_d = sqrt(-d / 108.);
    double u = cbrt(-q / 2. + sqrt_minus_d);
    double v = (fabs(u) > EPS_D) ? (-p / (3.0 * u)) : 0.0;
    lambda[0] = (u + v) + c2 * ONE_THIRD;
    return 1;
  }
  
  double acos_arg = -q / (2. * sqrt(-p * p * p / 27.));
  if (acos_arg > 1.)
    acos_arg = 1.;
  else if (acos_arg < -1.)
    acos_arg = -1.;

  double phi = acos(acos_arg) * ONE_THIRD, r = 2. * sqrt(-p * ONE_THIRD);
  c2 *= ONE_THIRD;

  if (d < EPS_D)
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

int eigen_values_to_d(double lambda[3], const float src[3][3])
{
  double c2 = (double)src[0][0] + (double)src[1][1] + src[2][2];
  double c1 = (
    (double)src[1][1] * src[2][2] - (double)src[2][1] * src[1][2]
    + (double)src[0][0] * src[2][2] - (double)src[2][0] * src[0][2]
    + (double)src[0][0] * src[1][1] - (double)src[1][0] * src[0][1]
  );
  double c0 = det_to_d(src);

  double p = c1 - c2 * c2 * ONE_THIRD;
  double q = (9. * c2 * c1 - 2. * c2 * c2 * c2) / 27. - c0;

  if (fabs(p) < EPS_D && fabs(q) < EPS_D) 
  {
    lambda[0] = lambda[1] = lambda[2] = c2 * ONE_THIRD;
    return 3;
  }

  double d = -4. * p * p * p - 27. * q * q;

  if (d < -EPS_D)
  {
    double sqrt_minus_d = sqrt(-d / 108.);
    double u = cbrt(-q / 2. + sqrt_minus_d);
    double v = (fabs(u) > EPS_D) ? (-p / (3.0 * u)) : 0.0;
    lambda[0] = (u + v) + c2 * ONE_THIRD;
    return 1;
  }
  
  double acos_arg = -q / (2. * sqrt(-p * p * p / 27.));
  if (acos_arg > 1.)
    acos_arg = 1.;
  else if (acos_arg < -1.)
    acos_arg = -1.;

  double phi = acos(acos_arg) * ONE_THIRD, r = 2. * sqrt(-p * ONE_THIRD);
  c2 *= ONE_THIRD;

  if (d < EPS_D)
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

int eigen_values_d(double lambda[3], const double src[3][3])
{
  double c2 = src[0][0] + src[1][1] + src[2][2];
  double c1 = (
    src[1][1] * src[2][2] - src[2][1] * src[1][2]
    + src[0][0] * src[2][2] - src[2][0] * src[0][2]
    + src[0][0] * src[1][1] - src[1][0] * src[0][1]
  );
  double c0 = det_d(src);

  double p = c1 - c2 * c2 * ONE_THIRD;
  double q = (9. * c2 * c1 - 2. * c2 * c2 * c2) / 27. - c0;

  if (fabs(p) < EPS_D && fabs(q) < EPS_D) 
  {
    lambda[0] = lambda[1] = lambda[2] = c2 * ONE_THIRD;
    return 3;
  }

  double d = -4. * p * p * p - 27. * q * q;

  if (d < -EPS_D)
  {
    double sqrt_minus_d = sqrt(-d / 108.);
    double u = cbrt(-q / 2. + sqrt_minus_d);
    double v = (fabs(u) > EPS_D) ? (-p / (3.0 * u)) : 0.0;
    lambda[0] = (u + v) + c2 * ONE_THIRD;
    return 1;
  }
  
  double acos_arg = -q / (2. * sqrt(-p * p * p / 27.));
  if (acos_arg > 1.)
    acos_arg = 1.;
  else if (acos_arg < -1.)
    acos_arg = -1.;

  double phi = acos(acos_arg) * ONE_THIRD, r = 2. * sqrt(-p * ONE_THIRD);
  c2 *= ONE_THIRD;

  if (d < EPS_D)
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

int eigen_values_symm_d(double lambda[3], const double src[3][3])
{
  double shift = (src[0][0] + src[1][1] + src[2][2]) * ONE_THIRD;

  double a00 = src[0][0] - shift, a11 = src[1][1] - shift, a22 = src[2][2] - shift;
  double a01 = src[0][1], a02 = src[0][2], a12 = src[1][2];

  double p = a00 * a11 + a11 * a22 + a22 * a00 - (a01 * a01 + a02 * a02 + a12 * a12);

  if (p >= 0.) 
  {
    lambda[0] = lambda[1] = lambda[2] = shift;
    return 3;
  }

  double q = a00 * a12 * a12 + a11 * a02 * a02 + a22 * a01 * a01 - a00 * a11 * a22 - 2.0 * a01 * a02 * a12;

  double minus_one_third_p = -p * ONE_THIRD;
  double sqrt_minus_p_third = sqrt(minus_one_third_p);
  double r = 2.0 * sqrt_minus_p_third;
  double denom = r * (minus_one_third_p);

  double acos_arg = -q / denom;

  if (acos_arg > 1.)
    acos_arg = 1.;
  else if (acos_arg < -1.)
    acos_arg = -1.;

  double phi = acos(acos_arg) * ONE_THIRD;

  double sin_phi = sin(phi), cos_phi = cos(phi);

  double sin_phi_sqrt3 = sin_phi * SQRT3;

  lambda[0] = r * cos_phi + shift;
  lambda[1] = -.5 * r * (cos_phi + sin_phi_sqrt3) + shift;
  lambda[2] = .5 * r * (sin_phi_sqrt3 - cos_phi) + shift;

  return 3;
}

void transposed(float result[restrict 3][3], const float src[restrict 3][3])
{
  result[0][0] = src[0][0], result[1][1] = src[1][1], result[2][2] = src[2][2];
  result[0][1] = src[1][0], result[0][2] = src[2][0], result[1][2] = src[2][1];
  result[1][0] = src[0][1], result[2][0] = src[0][2], result[2][1] = src[1][2];
}

void inverse(float result[restrict 3][3], const float src[restrict 3][3])
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

void scaled_m(float result[restrict 3][3], const float m[restrict 3][3], float scale)
{
  result[0][0] = m[0][0] * scale, result[1][1] = m[1][1] * scale, result[2][2] = m[2][2] * scale;
  result[0][1] = m[0][1] * scale, result[0][2] = m[0][2] * scale, result[1][2] = m[1][2] * scale;
  result[1][0] = m[1][0] * scale, result[2][0] = m[2][0] * scale, result[2][1] = m[2][1] * scale;
}

void triangulated(float result[restrict 3][3], const float src[restrict 3][3])
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

void lu(float l[restrict 3][3], float u[restrict 3][3], const float src[restrict 3][3])
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

void cholesky(float l[3][3], const float src[3][3])
{
  l[0][1] = l[0][2] = l[1][2] = 0.f;
  l[0][0] = sqrtf(src[0][0]);
  float div = 1.f / l[0][0];
  l[1][0] = src[0][1] * div, l[2][0] = src[0][2] * div;
  l[1][1] = sqrtf(src[1][1] - l[1][0] * l[1][0]), l[2][1] = (src[1][2] - l[2][0] * l[1][0]) / l[1][1];
  l[2][2] = sqrtf(src[2][2] - l[2][0] * l[2][0] - l[2][1] * l[2][1]);
}

static void kernel_v_d(double result[3], const double src[3][3])
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

void svd(float u[3][3], float s[3][3], float v_t[3][3], const float src[3][3])
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
  eigen_values_symm_d(sigmas, ata);

  if (sigmas[0] < sigmas[1]) { double temp = sigmas[0]; sigmas[0] = sigmas[1]; sigmas[1] = temp; }
  if (sigmas[0] < sigmas[2]) { double temp = sigmas[0]; sigmas[0] = sigmas[2]; sigmas[2] = temp; }
  if (sigmas[1] < sigmas[2]) { double temp = sigmas[1]; sigmas[1] = sigmas[2]; sigmas[2] = temp; }

  if (sigmas[0] < EPS_D)
  {
    memset(u, 0, sizeof(float) * 9);
    memset(s, 0, sizeof(float) * 9);
    memset(v_t, 0, sizeof(float) * 9);
    return;
  }

  double ata00 = ata[0][0], ata11 = ata[1][1], ata22 = ata[2][2];
  double v1[3], v2[3], v3[3];

  if (sigmas[0] - sigmas[2] < EPS_D)
  {
    v1[0] = 1., v1[1] = 0., v1[2] = 0.;
    v2[0] = 0., v2[1] = 1., v2[2] = 0.;
    v3[0] = 0., v3[1] = 0., v3[2] = 1.;
  }
  else if (sigmas[0] - sigmas[1] < EPS_D)
  {
    ata[0][0] -= sigmas[2], ata[1][1] -= sigmas[2], ata[2][2] -= sigmas[2];
    kernel_v_d(v3, ata);
    normalise_d(v3);

    double t[3] = { 0., 0., 0. };
    if (fabs(v3[0]) < .7)
      t[0] = 1.;
    else
      t[1] = 1.;
    cross_product_d(v2, v3, t);
    normalise_d(v2);
    
    cross_product_d(v1, v2, v3);
    normalise_d(v1);
  }
  else
  {
    ata[0][0] -= sigmas[0], ata[1][1] -= sigmas[0], ata[2][2] -= sigmas[0];
    kernel_v_d(v1, ata);
    
    normalise_d(v1);
  
    if (sigmas[0] - sigmas[1] > EPS_D)
    {
      ata[0][0] = ata00 - sigmas[1], ata[1][1] = ata11 - sigmas[1], ata[2][2] = ata22 - sigmas[1];
      kernel_v_d(v2, ata);
    }
    else
    {
      double t[3] = { 0., 0., 0. };
      if (fabs(v1[0]) < .7)
        t[0] = 1.;
      else
        t[1] = 1.;
      cross_product_d(v2, v1, t);
    }
    normalise_d(v2);
  
    cross_product_d(v3, v1, v2);
    normalise_d(v3);
  }

  v_t[0][0] = v1[0], v_t[0][1] = v1[1], v_t[0][2] = v1[2];
  v_t[1][0] = v2[0], v_t[1][1] = v2[1], v_t[1][2] = v2[2];
  v_t[2][0] = v3[0], v_t[2][1] = v3[1], v_t[2][2] = v3[2];

  s[0][0] = sigmas[0] = sqrt(sigmas[0]), s[1][1] = sigmas[1] = sqrt(sigmas[1] > EPS_D ? sigmas[1] : 0.);
  s[2][2] = sigmas[2] = sqrt(sigmas[2] > EPS_D ? sigmas[2] : 0.);
  s[1][0] = s[2][0] = s[2][1] = s[0][1] = s[0][2] = s[1][2] = 0.f;

  double u1[3], u2[3], u3[3];
  mult_m_v_d(u1, src_d, v1);
  scaled_v_d(u1, 1. / sigmas[0]);

  if (sigmas[1] > EPS_D)
  {
    mult_m_v_d(u2, src_d, v2);
    scaled_v_d(u2, 1. / sigmas[1]);
  
    if (sigmas[2] > EPS_D)
    {
      mult_m_v_d(u3, src_d, v3);
      scaled_v_d(u3, 1. / sigmas[2]);
    }
    else
    {
      cross_product_d(u3, u1, u2);
      scaled_v_d(u3, 1. / len_d(u3));
    }
  }
  else
  {
    double t[3] = { 0., 0., 0. };
    if (fabs(u1[0]) < .7)
      t[0] = 1.;
    else
      t[1] = 1.;
    cross_product_d(u2, u1, t);
    scaled_v_d(u2, 1. / len_d(u2));

    cross_product_d(u3, u1, u2);
    scaled_v_d(u3, 1. / len_d(u3));
  }

  u[0][0] = u1[0], u[1][0] = u1[1], u[2][0] = u1[2];
  u[0][1] = u2[0], u[1][1] = u2[1], u[2][1] = u2[2];
  u[0][2] = u3[0], u[1][2] = u3[1], u[2][2] = u3[2];
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
