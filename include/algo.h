#ifndef ALGO_H
#define ALGO_H
#include <math.h>

static inline float len(const float v[3])
{
  return sqrtf(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

static inline double len_d(const double v[3])
{
  return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

static inline double len_sq_d(const double v[3])
{
  return v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
}

static inline float tr(const float src[3][3])
{
  return src[0][0] + src[1][1] + src[2][2];
}

static inline void normalise_d(double v[3])
{
  double div = 1. / len_d(v);
  v[0] *= div, v[1] *= div, v[2] *= div;
}

static inline void scale_v_d(double v[3], double scale)
{
  v[0] *= scale, v[1] *= scale, v[2] *= scale;
}

float det(const float src[3][3]);
double det_to_d(const float src[3][3]);
double det_d(const double src[3][3]);

int eigen_values(float lambda[3], const float src[3][3]);
int eigen_values_to_d(double lambda[3], const float src[3][3]);
int eigen_values_d(double lambda[3], const double src[3][3]);
int eigen_values_symm_d(double lambda[3], const double src[3][3]);

void transposed(float result[restrict 3][3], const float src[restrict 3][3]);
void inverse(float result[restrict 3][3], const float src[restrict 3][3]);
void scaled_m(float result[restrict 3][3], const float m[restrict 3][3], float scale);
void triangulated(float result[restrict 3][3], const float src[restrict 3][3]);

void mult_m_v(float result[3], const float m[3][3], const float v[3]);
void mult_m_v_d(double result[3], const double m[3][3], const double v[3]);
void outer_product(float result[3][3], const float left[3], const float right[3]);
void cross_product(float result[3], const float left[3], const float right[3]);
void cross_product_d(double result[3], const double left[3], const double right[3]);
void mult_m_m_3x3(float result[restrict 3][3], const float left[restrict 3][3], const float right[restrict 3][3]);
void mult_m_m(
  int left_rows, int left_cols, int right_cols,
  float result[restrict left_rows][right_cols],
  const float left[restrict left_rows][left_cols],
  const float right[restrict left_cols][right_cols]
);
void add_m_m(float result[restrict 3][3], const float left[restrict 3][3], const float right[restrict 3][3]);

void rotation_mat_x(float result[3][3], float angle);
void rotation_mat_y(float result[3][3], float angle);
void rotation_mat_z(float result[3][3], float angle);
void rotation_mat(float result[3][3], const float axis[3], float angle);
void euler_to_mat_zyx(float result[3][3], float pitch, float yaw, float roll);
void cross_product_mat(float result[3][3], const float v[3]);

void lu(float l[restrict 3][3], float u[restrict 3][3], const float src[restrict 3][3]);
void cholesky(float l[3][3], const float src[3][3]);
void svd(float u[3][3], float s[3][3], float v_t[3][3], const float src[3][3]);

void solve_gauss(float result[3], const float src[3][3], const float b[3]);
void solve_lu(float result[3], const float l[3][3], const float u[3][3], const float b[3]);
void solve_cholesky(float result[3], const float l[3][3], const float b[3]);

#endif
