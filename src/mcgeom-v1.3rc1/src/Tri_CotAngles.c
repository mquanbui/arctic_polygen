#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "MCGeom.h"

#ifdef __cplusplus
extern "C" {
#endif

void Tri_CotAngles(double txyz[][3], double cotangs[]) {
  int i;
  double vec01[3], vec02[3], vec12[3], vec10[3], vec3[3];
  double dp, cpmag;

  VDiff3(txyz[1],txyz[0],vec01);
  VDiff3(txyz[2],txyz[0],vec02);
  VDiff3(txyz[2],txyz[1],vec12);

  dp = VDot3(vec01,vec02);  /* len01*len02*cos(ang201) */
  VCross3(vec01,vec02,vec3);
  cpmag = VLen3(vec3); /* len01*len02*sin(ang201) */
  cotangs[0] = dp/cpmag;  

  for (i = 0; i < 3; i++)
    vec10[i] = -vec01[i];

  dp = VDot3(vec12,vec10);  /* len12*len10*cos(ang012) */
  VCross3(vec12,vec10,vec3);
  cpmag = VLen3(vec3); /* len12*len10*sin(ang012) */
  cotangs[1] = dp/cpmag;

  dp = VDot3(vec02,vec12); /* len02*len12*cos(ang120) */
  VCross3(vec02,vec12,vec3);
  cpmag = VLen3(vec3); /* len02*len12*sin(ang120) */
  cotangs[2] = dp/cpmag;
}


#ifdef __cplusplus
	   }
#endif















