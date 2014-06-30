#include <math.h>
#include "MCGeom.h"


void Hex_CosAngles(double (*rxyz)[3], double *cosangs) {
  int i, j0, j1, j2;
  double vec01[3], vec02[3], n0[3], n1[3], dp;

  static int vtmpl[12][2][4] = {{{1,0,4,5},{0,1,2,3}},{{2,1,5,6},{1,2,3,0}},
				{{3,2,6,7},{2,3,0,1}},{{0,3,7,4},{3,0,1,2}},
				{{5,4,7,6},{4,5,1,0}},{{6,5,4,7},{5,6,2,1}},
				{{7,6,5,4},{6,7,3,2}},{{4,7,6,5},{7,4,0,3}},
				{{0,4,5,1},{0,3,7,4}},{{1,5,6,2},{1,0,4,5}},
				{{2,6,7,3},{2,1,5,6}},{{3,7,4,0},{3,2,6,7}}};
  

  /* Assume that vertices hex are ordered such that the 'bottom' face
     vertices are listed first ordered so that the normal points into
     the hex and the opposite face vertices are listed in the same
     order (such that the normal points out of the hex */

  /* This computation is just an approximation since quad faces are
     bi-linear and therefore the dihedral angle between them can vary
     along the common edge */

  /*                    6                5
                        *----------------*
                       /|               /|
                      / |              / |
                     /7 |            4/  |
                     *---------------*   |   
                     |  *------------|---*1
		     | /2            |  /
                     |/              | /
                     |3              |/
                     *---------------*0

  */



  for (i = 0; i < 12; i++) {
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

    cosangs[i] = sqrt((dp*dp)/(VLenSqr3(n0)*VLenSqr3(n1)));
  }

}
