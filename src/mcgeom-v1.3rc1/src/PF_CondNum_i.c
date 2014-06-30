#include <math.h>
#include "MCGeom.h"

#ifdef __cplusplus
extern "C" {
#endif

void PF_CondNum_i(double (*fxyz)[3], int n, int i, double *condnum) {
  int j;
  double len01_sqr=0.0, len02_sqr=0.0, vec01[3], vec02[3], areavec[3], area;
  
  VDiff3(fxyz[(i+n-1)%n],fxyz[i],vec02);
  len02_sqr = VLenSqr3(vec02);
    
  VDiff3(fxyz[(i+1)%n],fxyz[i],vec01);
  len01_sqr = VLenSqr3(vec01);

  VCross3(vec01,vec02,areavec);
  area = VLen3(areavec);
  *condnum = (len01_sqr + len02_sqr)/area;
}

#ifdef __cplusplus
}
#endif
