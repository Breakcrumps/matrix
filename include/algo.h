#ifndef ALGO_H
#define ALGO_H
#include <math.h>

static inline float len(float v[3])
{
  return sqrtf(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

static inline double len_d(double v[3])
{
  return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

static inline double len_sq_d(double v[3])
{
  return v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
}

static inline float tr(float src[3][3])
{
  return src[0][0] + src[1][1] + src[2][2];
}

static inline void normalise_d(double v[3])
{
  double len = len_d(v);
  
  if (len < 1e-9)
    return;
  
  double div = 1. / len;
  v[0] *= div, v[1] *= div, v[2] *= div;
}

static inline void scale_v_d(double v[3], double scale)
{
  v[0] *= scale, v[1] *= scale, v[2] *= scale;
}

float det(float src[3][3]);
double det_to_d(float src[3][3]);
double det_d(double src[3][3]);
int eigen_values(float lambda[3], float src[3][3]);
int eigen_values_to_d(double lambda[3], float src[3][3]);
int eigen_values_d(double lambda[3], double src[3][3]);

void transposed(float result[3][3], float src[3][3]);
void inverse(float result[3][3], float src[3][3]);
void scaled_m(float result[3][3], float m[3][3], float scale);
void triangulated(float result[3][3], float src[3][3]);

void mult_m_v(float result[3], float m[3][3], float v[3]);
void mult_m_v_d(double result[3], double m[3][3], double v[3]);
void outer_product(float result[3][3], float left[3], float right[3]);
void cross_product(float result[3], float left[3], float right[3]);
void cross_product_d(double result[3], double left[3], double right[3]);
void mult_m_m(float result[3][3], float left[3][3], float right[3][3]);
void add_m_m(float result[3][3], float left[3][3], float right[3][3]);

void rotation_mat_x(float result[3][3], float angle);
void rotation_mat_y(float result[3][3], float angle);
void rotation_mat_z(float result[3][3], float angle);
void rotation_mat(float result[3][3], float axis[3], float angle);
void euler_to_mat_zyx(float result[3][3], float pitch, float yaw, float roll);
void cross_product_mat(float result[3][3], float v[3]);

void lu(float l[3][3], float u[3][3], float src[3][3]);
void cholesky(float l[3][3], float src[3][3]);
void svd(float u[3][3], float s[3][3], float v_t[3][3], float src[3][3]);

void solve_gauss(float result[3], float src[3][3], float b[3]);
void solve_lu(float result[3], float l[3][3], float u[3][3], float b[3]);
void solve_cholesky(float result[3], float l[3][3], float b[3]);

#endif
