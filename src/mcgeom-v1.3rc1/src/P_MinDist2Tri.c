#include <math.h>
#include "MCGeom.h"


#ifdef __cplusplus
extern "C" {
#endif

double P_MinDist2Tri(double vxyz[], double fxyz[][3], double eps) {
  double mindist2, dist, dist2, pxyz[3], normal[3], avec[3], lxyz[2][3];
  double vec0[3], vec1[3], dp, dummy;
  int i, j, intri, inedge, sgn;

  mindist2 = 1.0e+14;

  /* Project point on plane of triangle */
  P_Prj2TriPlane(vxyz, fxyz, pxyz, &dist);
  dist2 = dist*dist;

  /* If the projection is inside the triangle we found the minimum distance */
  intri = 1;
  
  VDiff3(fxyz[1],fxyz[0],vec0);
  VDiff3(fxyz[2],fxyz[0],vec1);
  VCross3(vec0,vec1,normal);

  for (i = 0; i < 3; i++) {
    VDiff3(fxyz[(i+1)%3],fxyz[i],vec0);
    VDiff3(pxyz,fxyz[i],vec1);
    VCross3(vec0,vec1,avec);
    if (VDot3(normal,avec) < 0.0) {
      intri = 0;
      break;
    }
  }      

  if (intri)
    mindist2 = dist2;
  else {
    inedge = 0;

    /* Projection of point is outside triangle - find closest point on
       triangle boundary */

    /* project point onto edges */

    for (i = 0; i < 3; i++) {
      for (j = 0; j < 3; j++) {
	lxyz[0][j] = fxyz[i][j];
	lxyz[1][j] = fxyz[(i+1)%3][j];
      }
      P_Prj2Line(vxyz, lxyz, pxyz, &dummy);

      VDiff3(lxyz[0],pxyz,vec0);
      VDiff3(lxyz[1],pxyz,vec1);
      dp = VDot3(vec0,vec1);
      sgn = dp > 0 ? 1 : -1;
      if (fabs(dp) < eps*eps) {
	/* projection is coincident with one of the vertices but it
           may not be the closest vertex of the triangle */
	if (dist2 < mindist2)
	  mindist2 = dist2;
	continue;
      }

      dp = sgn*dp*dp/(VLenSqr3(vec0)*VLenSqr3(vec1));
      if (fabs(dp-(-1.0)) < eps) {
	if (dist2 < mindist2)
	  mindist2 = dist2;
      }
    }

    /* Also check for the closest vertex */

    for (i = 0; i < 3; i++) {
      dist2 = ((fxyz[i][0]-vxyz[0])*(fxyz[i][0]-vxyz[0]) +
	       (fxyz[i][1]-vxyz[1])*(fxyz[i][1]-vxyz[1]) +
	       (fxyz[i][2]-vxyz[2])*(fxyz[i][2]-vxyz[2]));
      if (dist2 < mindist2) 
	mindist2 = dist2;
    }
  }

  return mindist2;
}
  
#ifdef __cplusplus
}
#endif
