#include "helper.h"

double getRandom()
{
  return (5 - 2) * rand() / (double)RAND_MAX + 2;
}

void setMatrix(int n, double *a)
{
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
    {
      a[i * j] = getRandom();
    }
  }
}