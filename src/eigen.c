#include "../include/algo.h"

#define TWO_THIRDS_PI 2.0943951023931954923
#define EPS_D 1e-14
#define SQRT3 1.73205080756887729
#define ONE_THIRD 0.333333333333333314829616256247390992939472198486328125

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
