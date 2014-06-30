#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "MCGeom.h"

#ifdef __cplusplus
extern "C" {
#endif

  /* Find projection of point to an infinite line passing through
     lxyz[0] and lxyz[1]. Return the square of the perpendicular
     distance from point to the line */

void P_Prj2Line(double pxyz[], double lxyz[][3], double prxyz[], double *distsqr) {
  double lvec[3], pvec[3], dp;
  int i;
  
  VDiff3(lxyz[1],lxyz[0],lvec);
  VNormalize3(lvec);
  
  VDiff3(pxyz,lxyz[0],pvec);
  
  dp = VDot3(pvec,lvec);
  
  for (i = 0; i < 3; i++)
    prxyz[i] = lxyz[0][i] + dp*lvec[i];
  
  *distsqr = ((pxyz[0]-prxyz[0])*(pxyz[0]-prxyz[0]) +
	      (pxyz[1]-prxyz[1])*(pxyz[1]-prxyz[1]) +
	      (pxyz[2]-prxyz[2])*(pxyz[2]-prxyz[2]));
  
}

#ifdef __cplusplus
}
#endif
