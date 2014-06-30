#include "MCGeom.h"

/* Geometric center */

void PF_Center(int n, double (*xyz)[3], double *pcen) {
  int i, j;
  
  for (j = 0; j < 3; j++)
    pcen[j] = 0.0;
    
  for (i = 0; i < n; i++) {
    for (j = 0; j < 3; j++)
      pcen[j] += xyz[i][j];
  }
    
  for (j = 0; j < 3; j++)
    pcen[j] /= n;
}

  
