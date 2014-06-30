#include <math.h>
#include "MCGeom.h"

/* Subdivide quad into two triangles such that the dihedral angle
   between the triangles is maximum. In most cases this ensures that
   the patch of surface represented by the is approximated accurately
   by the two triangles */

#ifdef __cplusplus
extern "C" {
#endif

void Quad_BestTriangles(double fxyz[][3], int t1ind[], double t1xyz[][3], int t2ind[], double t2xyz[][3]) {
  double vec0[3], vec1[3], n0[3], n1[3], dp13, dp20;
  int sgn, i;
  
  /* Dot Product of normals when quad is divided along 1,3 into two
     triangles, 0,1,3 and 1,2,3 */
  VDiff3(fxyz[1],fxyz[0],vec0);
  VDiff3(fxyz[3],fxyz[0],vec1);
  VCross3(vec0,vec1,n0);

  VDiff3(fxyz[3],fxyz[2],vec0);
  VDiff3(fxyz[1],fxyz[2],vec1);
  VCross3(vec0,vec1,n1);
  
  dp13 = VDot3(n0,n1);
  sgn = dp13 > 0 ? 1 : -1;
  dp13 = sgn*dp13*dp13/(VLenSqr3(n0)*VLenSqr3(n1));
  
  /* Dot Product of normals when quad is divided along 2,0 into two
     triangles, 0,1,2 and 0,2,3 */
  VDiff3(fxyz[2],fxyz[1],vec0);
  VDiff3(fxyz[0],fxyz[1],vec1);
  VCross3(vec0,vec1,n0);

  VDiff3(fxyz[0],fxyz[3],vec0);
  VDiff3(fxyz[2],fxyz[3],vec1);
  VCross3(vec0,vec1,n1);
  
  dp20 = VDot3(n0,n1);
  sgn = dp20 > 0 ? 1 : -1;
  dp20 = sgn*dp20*dp20/(VLenSqr3(n0)*VLenSqr3(n1));
  
  if (dp13 < dp20) {
    /* dihedral angle along 1,3 is greater than along 2,0 */
    t1ind[0] = 0; t1ind[1] = 1; t1ind[2] = 3; 
    for (i = 0; i < 3; i++) t1xyz[0][i] = fxyz[0][i];
    for (i = 0; i < 3; i++) t1xyz[1][i] = fxyz[1][i];
    for (i = 0; i < 3; i++) t1xyz[2][i] = fxyz[3][i];

    t2ind[0] = 1; t2ind[1] = 2; t2ind[2] = 3;
    for (i = 0; i < 3; i++) t2xyz[0][i] = fxyz[1][i];
    for (i = 0; i < 3; i++) t2xyz[1][i] = fxyz[2][i];
    for (i = 0; i < 3; i++) t2xyz[2][i] = fxyz[3][i];
  }
  else {
    /* dihedral angle along 2,0 is greater than along 1,3 */
    t1ind[0] = 0; t1ind[1] = 1; t1ind[2] = 2;
    for (i = 0; i < 3; i++) t1xyz[0][i] = fxyz[0][i];
    for (i = 0; i < 3; i++) t1xyz[1][i] = fxyz[1][i];
    for (i = 0; i < 3; i++) t1xyz[2][i] = fxyz[2][i];

    t2ind[0] = 0; t2ind[1] = 2; t2ind[2] = 3;
    for (i = 0; i < 3; i++) t2xyz[0][i] = fxyz[0][i];
    for (i = 0; i < 3; i++) t2xyz[1][i] = fxyz[2][i];
    for (i = 0; i < 3; i++) t2xyz[2][i] = fxyz[3][i];
  }
}

#ifdef __cplusplus
}
#endif
