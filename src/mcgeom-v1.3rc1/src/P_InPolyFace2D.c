#include <stdio.h>
#include <stdlib.h>
#include "MCGeom.h"

/* For now assume the polygon is in 2D. We will use the Jordan curve
   theorem to test in/out. 

   Shoot a ray in the +x direction and count the line crossings. The
   ray is said to have crossed the line if one point is strictly above
   the ray and another point is on or below the ray. If the number of
   crossings is odd, the point is inside; if it is even, the point is
   outside.

   if flag is 0, we don't care how boundary points are classified (fast)

   If flag is set to 1, we want points on the boundary to be
   consistently classified as being inside.  The above algorithm can
   give inconsistent results for points on the boundary. So we do
   additional geometric checks to see if a point classified as
   exterior is really on the boundary */


int P_InPolyFace2D(double *xyz, int np, double (*pxyz)[3], double tol, 
		   int flag) {
  int i, ip1, c;
  double x, y, minx, miny, maxx, maxy, elen, plen, dist2;
  double vec1[3], evec[3], prxyz[3], lxyz[2][3];

  /* So eliminate obviously external points by a bounding box check */

  /*
  minx = miny = 1.0e+16; maxx = maxy = -1.0e-16;
  for (i = 0; i < np; i++) {
    if (pxyz[i][0] < minx) minx = pxyz[i][0];
    if (pxyz[i][0] > maxx) maxx = pxyz[i][0];
    if (pxyz[i][1] < miny) miny = pxyz[i][1];
    if (pxyz[i][1] > maxy) maxy = pxyz[i][1];
  }

  if (xyz[0] < minx || xyz[0] > maxx || xyz[1] < miny || xyz[1] > maxy)
    return 0;
  */

  /* Basic test - will work for strictly interior and exterior points */

  x = xyz[0]; y = xyz[1];

  for (i = 0, c = 0; i < np; i++) {
    ip1 = (i+1)%np;
    if (((pxyz[i][1] > y && pxyz[ip1][1] <= y) ||
	 (pxyz[ip1][1] > y && pxyz[i][1] <= y)) &&
	(x <= (pxyz[i][0] + (y-pxyz[i][1])*(pxyz[ip1][0]-pxyz[i][0])/(pxyz[ip1][1]-pxyz[i][1]))))
      c = !c;
  }

  /* If we don't need consistent classification of points on the
     boundary, we can quit here.  If the point is classified as
     inside, it is definitely inside or on the boundary - no way it
     can be outside and be classified as inside */

  if (!flag || c == 1)
    return c;


  /* If the point is classified as outside but we need all boundary
     points to be consistently classified as being inside, we need
     additional checks */

  for (i = 0; i < np; i++) {
    ip1 = (i+1)%np;

     VCopy3(lxyz[0],pxyz[i]);
     VCopy3(lxyz[1],pxyz[(i+1)%np]);

     if (P_OnLineSeg(xyz,lxyz,tol))
       return 1;
  }

  return 0;

} 
