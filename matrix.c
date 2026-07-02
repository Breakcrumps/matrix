#include <stdio.h>
#include "include/io.h"
#include "include/algo.h"

int main()
{
  double m[3][3], eigen[3];
  read_m_d(m);
  eigen_values_d(eigen, m);
  print_v_msg_d(eigen, "Eigen: ");

  return 0;
}
