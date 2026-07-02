#include "../include/io.h"
#include <stdio.h>

void read_v(float v[3])
{
  scanf("%f%f%f", &v[0], &v[1], &v[2]);
}

void read_m(float m[3][3])
{
  scanf(
    "%f%f%f%f%f%f%f%f%f",
    &m[0][0], &m[0][1], &m[0][2],
    &m[1][0], &m[1][1], &m[1][2],
    &m[2][0], &m[2][1], &m[2][2]
  );
}

void read_m_d(double m[3][3])
{
  scanf(
    "%lf%lf%lf%lf%lf%lf%lf%lf%lf",
    &m[0][0], &m[0][1], &m[0][2],
    &m[1][0], &m[1][1], &m[1][2],
    &m[2][0], &m[2][1], &m[2][2]
  );
}

void print_v(float v[3])
{
  printf("%f %f %f\n", v[0], v[1], v[2]);
}

void print_m(float m[3][3])
{
  printf(
    "%f %f %f\n%f %f %f\n%f %f %f\n",
    m[0][0], m[0][1], m[0][2],
    m[1][0], m[1][1], m[1][2],
    m[2][0], m[2][1], m[2][2]
  );
}

void print_v_msg(float v[3], char pre_msg[])
{
  puts(pre_msg);
  printf("%f %f %f\n", v[0], v[1], v[2]);
}

void print_v_msg_d(double v[3], char pre_msg[])
{
  puts(pre_msg);
  printf("%lf %lf %lf\n", v[0], v[1], v[2]);
}

void print_m_msg(float m[3][3], char pre_msg[])
{
  puts(pre_msg);
  print_m(m);
}
