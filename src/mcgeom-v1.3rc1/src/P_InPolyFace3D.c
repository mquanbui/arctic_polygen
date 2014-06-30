#include <stdio.h>
#include <stdlib.h>
#include "MCGeom.h"

/* Test if point is inside a 3D (but planar) polygonal face
   Face has to be convex

   The algorithm checks if the triangle formed by the point and each
   of the polygon's edges has +ve area and the same normal as the
   normal of the original face

   The variable 'flag' has no effect here. Points on the boundary are
   flagged as being inside 


   Also, we are not using 'tol' - only relevant if we are trying to
   determine if points are on the boundary
*/


int P_InPolyFace3D(double *xyz, int np, double (*pxyz)[3], double tol, 
		   int flag) {
  int i, ip1, im1, done;
  double vec1[3], vec2[3], refnormal[3], normal[3], dp, len2;
  double MACH_EPS=1.0e-15;


  /* Get a normal to the face - must make sure it is not degenerate */
  
  done = 0;
  i = 0;
  while (!done) {
    ip1 = (i+1)%np;
    im1 = (i-1+np)%np;
    VDiff3(pxyz[ip1],pxyz[i],vec1);
    VDiff3(pxyz[im1],pxyz[i],vec2);
    VCross3(vec1,vec2,refnormal);
    len2 = VLenSqr3(refnormal);
    if (len2 > MACH_EPS || i == np-1) done = 1;
    i++;
  }

  
  /* Test the triangles formed by the point xyz and each edge of the polygon */

  for (i = 0; i < np; i++) {

    ip1 = (i+1)%np;

    VDiff3(pxyz[i],xyz,vec1);
    VDiff3(pxyz[ip1],xyz,vec2);
    VCross3(vec1,vec2,normal);
    dp = VDot3(normal,refnormal);

    if (dp < 0)
      return 0;

  }

  return 1;

} 
