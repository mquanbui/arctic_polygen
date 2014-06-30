#include <math.h>
#include "MCGeom.h"

#ifdef __cplusplus
extern "C" {
#endif

double P_MinDist2Quad(double vxyz[], double fxyz[][3], double eps) {
  double dist21, dist22, mindist2, t1xyz[3][3], t2xyz[3][3];
  int t1ind[3], t2ind[3];
  
  Quad_BestTriangles(fxyz,t1ind,t1xyz,t2ind,t2xyz);
  dist21 = P_MinDist2Tri(vxyz, t1xyz, eps);
  dist22 = P_MinDist2Tri(vxyz, t2xyz, eps);

  mindist2 = dist21 < dist22 ? dist21 : dist22;

  return mindist2;
}

#ifdef __cplusplus
}
#endif
