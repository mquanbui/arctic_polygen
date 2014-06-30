#include <math.h>
#include "MCGeom.h"


void Tet_Angles(double (*rxyz)[3], double *angles) {
  int i, j0, j1, j2;
  double vec01[3], vec02[3], n0[3], n1[3], dp;

  static int vtmpl[6][2][3] = {{{0,1,2},{1,0,3}},{{1,2,0},{2,1,3}},
			       {{2,0,1},{0,2,3}},{{0,3,1},{3,0,2}},
			       {{1,3,2},{1,0,3}},{{2,3,0},{3,2,1}}};
  

  /* Assume that vertices of tet are ordered so that they obey the
     right-hand rule. In other words, the normal of the face formed by
     the first three vertices points towards the fourth vertex */


  for (i = 0; i < 6; i++) {
    j0 = vtmpl[i][0][0];
    j1 = vtmpl[i][0][1];
    j2 = vtmpl[i][0][2];

    VDiff3(rxyz[j1],rxyz[j0],vec01);
    VDiff3(rxyz[j2],rxyz[j0],vec02);
    VCross3(vec01,vec02,n0);

    j0 = vtmpl[i][1][0];
    j1 = vtmpl[i][1][1];
    j2 = vtmpl[i][1][2];

    VDiff3(rxyz[j1],rxyz[j0],vec01);
    VDiff3(rxyz[j2],rxyz[j0],vec02);
    VCross3(vec01,vec02,n1);

    dp = VDot3(n0,n1);

    dp = sqrt((dp*dp)/(VLenSqr3(n0)*VLenSqr3(n1)));
    angles[i] = acos(dp);
  }

}
