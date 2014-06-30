#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "MCGeom.h"

/*
   Reflect a point 'p' about an arbitrary plane described by normal 'n' and
   point 'v'
   Find distance 'D'  of 'p' to the plane and find reflected point by
   the formula {p}-2D{n}
*/

#ifdef __cplusplus
extern "C" {
#endif

void P_ReflectGen(double p[], double v[], double n[], double rp[]) {
  double D, vecvp[3];
  int i;

  /* Assume 'n' is normalized, otherwise uncomment next line */
  /* normVt(n,n); */
  
  VDiff3(p,v,vecvp);
  D = fabs(VDot3(vecvp,n));
  
  for (i = 0; i < 3; i++)
    rp[i] = p[i] - 2*D*n[i];
}

#ifdef __cplusplus
}
#endif
