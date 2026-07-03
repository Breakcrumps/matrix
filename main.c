#include <stdio.h>
#include "include/io.h"
#include "include/algo.h"

int main()
{
  float m[3][3], t[3][3];
  read_m_msg(m, "Input matrix:");
  transposed(t, m);
  print_m_msg(t, "Transposed:");
  read_m_msg(m, "Input matrix:");
  transposed(t, m);
  print_m_msg(t, "Transposed:");
  read_m_msg(m, "Input matrix:");
  transposed(t, m);
  print_m_msg(t, "Transposed:");
  read_m_msg(m, "Input matrix:");
  transposed(t, m);
  print_m_msg(t, "Transposed:");
  read_m_msg(m, "Input matrix:");
  transposed(t, m);
  print_m_msg(t, "Transposed:");
  read_m_msg(m, "Input matrix:");
  transposed(t, m);
  print_m_msg(t, "Transposed:");
  read_m_msg(m, "Input matrix:");
  transposed(t, m);
  print_m_msg(t, "Transposed:");
  read_m_msg(m, "Input matrix:");
  transposed(t, m);
  print_m_msg(t, "Transposed:");
  return 0;
}
