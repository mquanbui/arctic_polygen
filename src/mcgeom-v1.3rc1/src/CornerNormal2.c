#include <math.h>
#include "MCGeom.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Return the normal at a corner defined by the vectors (xyz1-xyz0) and
   (xyz2-xyz0). The normal is not normalized */

void CornerNormal2(double xyz0[3], double xyz1[3], double xyz2[3], double *normal) {
  double vec01[3], vec02[3];
  
  VDiff3(xyz2,xyz0,vec02);
    
  VDiff3(xyz1,xyz0,vec01);

  VCross3(vec01,vec02,normal);
}

#ifdef __cplusplus
}
#endif
