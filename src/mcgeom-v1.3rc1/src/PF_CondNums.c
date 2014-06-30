#include <math.h>
#include "MCGeom.h"

#ifdef __cplusplus
extern "C" {
#endif

void PF_CondNums(double (*fxyz)[3], int n, double *condnums) {
  int i, j;
  double len01_sqr=0.0, len02_sqr=0.0, vec01[3], vec02[3], areavec[3], area;
  
  for (i = 0; i < n; i++) {
    if (i == 0) {
      VDiff3(fxyz[n-1],fxyz[0],vec02);
      len02_sqr = VLenSqr3(vec02);
    }
    else {
      for (j = 0; j < 3; j++) 
	vec02[j] = -vec01[j];
      len02_sqr = len01_sqr;
    }
    
    VDiff3(fxyz[(i+1)%n],fxyz[i],vec01);
    len01_sqr = VLenSqr3(vec01);

    VCross3(vec01,vec02,areavec);
    area = VLen3(areavec);
    condnums[i] = (len01_sqr + len02_sqr)/area;
  }
}

#ifdef __cplusplus
}
#endif
