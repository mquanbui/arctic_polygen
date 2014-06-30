#include <stdio.h>
#include <stdlib.h>
#include "MCGeom.h"

/* Check if a point is on a line segment (to within tolerance) */ 
/* We could do a cross product check to see if point is within
   tolerance of line and then do a dot product to see if point lies
   withing line segment */

int P_OnLineSeg(double *xyz, double (*lxyz)[3], double tol) {
  double lbb[2][3], elen, plen, dist2;
  double vec1[3], evec[3], prxyz[3];

  /* So eliminate obviously far away points by a bounding box check */

  /*
  for (int i = 0; i < 3; i++) {
    if (lxyz[0][i] < lxyz[1][i]) {
      lbb[0][i] = lxyz[0][i]; lbb[1][i] = lxyz[1][i];
    }
    else {
      lbb[0][i] = lxyz[1][i]; lbb[1][i] = lxyz[0][i];
    }
  }
    
  if (xyz[0] < lbb[0][0] || xyz[0] > lbb[1][0] || 
      xyz[1] < lbb[0][1] || xyz[1] > lbb[1][1] || 
      xyz[2] < lbb[0][2] || xyz[2] > lbb[1][2])
    return 0;
  */


  /* Vector from point 1 of line segment to point */

  VDiff3(xyz,lxyz[0],vec1);


  /* line segment vector and its length */

  VDiff3(lxyz[1],lxyz[0],evec);
  elen = VLen3(evec);


  /* Normalize */

  evec[0] /= elen;
  evec[1] /= elen;
  evec[2] /= elen;


  /* Projection of vec1 onto evec (normalized) OR distance b/w P0
     and projection of point onto the line */

  plen = VDot3(vec1,evec);


  /* Check if projection of point is inside the segment - true only the
     (signed) distance from the projected point to pnt 1 of the segment
     is +ve but less than length of segment (accounting for tolerance) */

  if (plen < -tol || plen > elen+tol)
    return 0;


  /* Projected point */

  prxyz[0] = lxyz[0][0] + plen*evec[0];
  prxyz[1] = lxyz[0][1] + plen*evec[1];
  prxyz[2] = lxyz[0][2] + plen*evec[2];
    

  /* distance of point to its projection (and therefore, the line)
     must be below tolerance */

  dist2 = ((prxyz[0]-xyz[0])*(prxyz[0]-xyz[0]) +
	   (prxyz[1]-xyz[1])*(prxyz[1]-xyz[1]) +
	   (prxyz[2]-xyz[2])*(prxyz[2]-xyz[2]));

  if (dist2 <= tol*tol) 
    return 1;
  else
    return 0;
} 
