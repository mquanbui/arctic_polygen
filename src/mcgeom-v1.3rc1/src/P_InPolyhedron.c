#include <stdio.h>
#include <stdlib.h>
#include "MCGeom.h"

/* Test if point is inside a 3D polyhedral region - faces of polyhedron
   can be non-planar but we assume that they can be broken into planar 
   triangular subfaces by connecting each edge to a central point

   The algorithm checks if the tetrahedron formed by the point and each
   of the polyhedrons triangular subfaces has +ve volume

   The variable 'flag' has no effect here. Points on the boundary are
   flagged as being inside 

   Also, we are not using 'tol' - only relevant if we are trying to
   determine if points are on the boundary
*/
 

int P_InPolyhedron(double *xyz, int nf, int *nfp, double (**pxyz)[3], 
		   double tol, int flag) {
  int i, j, k, jp1;
  double txyz[4][3];
  double MACH_EPS=1.0e-15;

  VCopy3(txyz[3],xyz);

  for (i = 0; i < nf; i++) {

    for (k = 0; k < 3; k++) txyz[2][k] = 0.0;
    for (j = 0; j < nfp[i]; j++)
      for (k = 0; k < 3; k++) txyz[2][k] += pxyz[i][j][k];
    for (k = 0; k < 3; k++) txyz[2][k] /= nfp[i];

    for (j = 0; j < nfp[i]; j++) {
      jp1 = (j+1)%nfp[i];
      VCopy3(txyz[1],pxyz[i][j]);
      VCopy3(txyz[0],pxyz[i][jp1]);
      if (Tet_Volume(txyz) < 0.0)
	return 0;      
    }
  }

  return 1;

} 
