#include "../include/algo.h"
#include <string.h>

#define EPS_D 1e-14

void lu(float l[3][3], float u[3][3], const float src[3][3])
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
