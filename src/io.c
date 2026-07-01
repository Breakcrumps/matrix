#include "../include/io.h"
#include <stdio.h>

void read_v(float v[3])
{
  for (int i = 0; i < 3; i++)
  {
    scanf("%f", &v[i]);
  }
}

void read_m(float m[3][3])
{
  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      scanf("%f", &m[i][j]);
    }
  }
}

void print_v(float v[3])
{
  for (int i = 0; i < 3; i++)
  {
    printf("%f ", v[i]);
  }

  putchar('\n');
}

void print_m(float m[3][3])
{
  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      printf("%f ", m[i][j]);
    }

    putchar('\n');
  }
}

void print_v_msg(float v[3], char pre_msg[])
{
  puts(pre_msg);
  print_v(v);
}

void print_m_msg(float m[3][3], char pre_msg[])
{
  puts(pre_msg);
  print_m(m);
}
