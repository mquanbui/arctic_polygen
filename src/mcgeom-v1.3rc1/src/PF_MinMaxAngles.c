#include <math.h>
#include "MCGeom.h"

#ifdef __cplusplus
extern "C" {
#endif

void PF_MinMaxAngles(double (*pxyz)[3], int n, double *minang, double *maxang) {
  int i, j;
  double len01_sqr=0.0, len02_sqr=0.0, vec01[3], vec02[3], dp, mincos, maxcos;

  mincos = 99.0;
  maxcos = -99.0;

  for (i = 0; i < n; i++) {
    if (i == 0) {
      VDiff3(pxyz[n-1],pxyz[0],vec02);
      len02_sqr = VLenSqr3(vec02);
    }
    else {
      for (j = 0; j < 3; j++) 
	vec02[j] = -vec01[j];
      len02_sqr = len01_sqr;
    }
    
    VDiff3(pxyz[(i+1)%n],pxyz[i],vec01);
    len01_sqr = VLenSqr3(vec01);

    dp = VDot3(vec01,vec02);
    dp = dp/sqrt(len01_sqr*len02_sqr);
    if (dp > mincos)
      mincos = dp;
    if (dp < maxcos)
      maxcos = dp;
  }


  *minang = acos(maxcos);
  *maxang = acos(mincos);
}

#ifdef __cplusplus
}
#endif
