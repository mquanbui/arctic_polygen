#include <stdio.h>
#include <math.h>
#include "MCGeom.h"


#ifdef __cplusplus
extern "C" {
#endif

void Tri_UnitNormal(double (*xyz)[3], double *normal) {
  double vec01[3], vec02[3];
  
  VDiff3(xyz[1],xyz[0],vec01);
  VDiff3(xyz[2],xyz[0],vec02);

  VCross3(vec01,vec02,normal);
  VNormalize3(normal);
}

#ifdef __cplusplus
}
#endif 
