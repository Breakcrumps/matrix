#include <stdio.h>
#include "include/io.h"
#include "include/algo.h"

int main()
{
  float m[3][3], u[3][3], s[3][3], v_t[3][3];
  read_m(m);
  svd(u, s, v_t, m);
  print_m(u);
  print_m(s);
  print_m(v_t);

  return 0;
}
