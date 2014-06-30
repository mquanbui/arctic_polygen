#include <math.h>
#include "MCGeom.h"


#ifdef __cplusplus
extern "C" {
#endif

void P_Prj2TriPlane(double vxyz[3], double (*fxyz)[3], double prj[3], double *dist) {
  int i;
  double vec01[3], vec02[3], vec0v[3], normal[3];

  VDiff3(fxyz[1],fxyz[0],vec01);
  VDiff3(fxyz[2],fxyz[0],vec02);
  VCross3(vec01,vec02,normal);
  VNormalize3(normal);

  VDiff3(vxyz,fxyz[0],vec0v);
  *dist = VDot3(vec0v,normal);

  for (i = 0; i < 3; i++)
    prj[i] = vxyz[i] - (*dist)*normal[i];
} 


#ifdef __cplusplus
}
#endif
