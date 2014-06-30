#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "MCGeom.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Assume that the corner is a trivalent corner */

double Corner_CondNum3(double xyz0[3], double xyz1[3], double xyz2[3], double xyz3[3]) {
  int i, j, iv0, iv1, iv2, iv3, gen;
  double vec01[3], vec02[3], vec03[3], cpvec[3], a, b, vol6, condnum;

  VDiff3(xyz1,xyz0,vec01);
  VDiff3(xyz2,xyz0,vec02);
  VDiff3(xyz3,xyz0,vec03);
    
  a = VLenSqr3(vec01) + VLenSqr3(vec02) + VLenSqr3(vec03);
  VCross3(vec02,vec03,cpvec);
  b = VLenSqr3(cpvec);
  VCross3(vec03,vec01,cpvec);
  b += VLenSqr3(cpvec);
  VCross3(vec01,vec02,cpvec);
  b += VLenSqr3(cpvec);
  vol6 = VDot3(cpvec,vec03);
    
  condnum = sqrt(a*b)/vol6;
  return condnum;
}

#ifdef __cplusplus
}
#endif
