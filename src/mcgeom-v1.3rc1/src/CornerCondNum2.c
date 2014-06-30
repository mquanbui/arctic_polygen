#include <math.h>
#include "MCGeom.h"

#ifdef __cplusplus
extern "C" {
#endif

double Corner_CondNum2(double xyz0[3], double xyz1[3], double xyz2[3], double refnormal[3]) {
  int j;
  double len01_sqr=0.0, len02_sqr=0.0, vec01[3], vec02[3];
  double condnum, areavec[3], area;
  
  VDiff3(xyz2,xyz0,vec02);
  len02_sqr = VLenSqr3(vec02);
    
  VDiff3(xyz1,xyz0,vec01);
  len01_sqr = VLenSqr3(vec01);

  VCross3(vec01,vec02,areavec);
  
  /* To get the signed area we have to dot area vector with reference normal */
  /* Reference normal - local normal, normal to local tangent plane or in    */
  /* case of 2D meshes, normal to plane in which mesh is                     */
  /* area of parallelogram = 2 times area of triangle */;

  area = VDot3(areavec,refnormal);

  condnum = (len01_sqr + len02_sqr)/area;
  return condnum;
}

#ifdef __cplusplus
}
#endif
