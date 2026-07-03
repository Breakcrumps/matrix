#include <stdio.h>
#include "include/io.h"
#include "include/algo.h"

int main()
{
  float u[3][3], s[3][3], v_t[3][3], m[3][3];
  read_m_msg(m, "Input matrix:");
  svd(u, s, v_t, m);
  print_m_msg(u, "U:");
  print_m_msg(s, "Sigma:");
  print_m_msg(v_t, "V_T:");
  return 0;
}
