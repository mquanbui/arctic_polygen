#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "MCGeom.h"

#ifdef __cplusplus
extern "C" {
#endif

double TriCotAngle(double txyz[][3], int idx) {
  double vec01[3], vec02[3], vec3[3];
  double dp, cpmag, cotang;

  VDiff3(txyz[(idx+1)%3],txyz[idx],vec01);
  VDiff3(txyz[(idx+2)%3],txyz[idx],vec02);

  dp = VDot3(vec01,vec02);  /* len01*len02*cos(ang201) */
  VCross3(vec01,vec02,vec3);
  cpmag = VLen3(vec3); /* len01*len02*sin(ang201) */
  cotang = dp/cpmag;  

  return cotang;
}


#ifdef __cplusplus
	   }
#endif















