#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "MCGeom.h"

#ifdef __cplusplus
extern "C" {
#endif

  /* Given an arbitrary point, find the closest point on a line
     segment. If the projection of the point to the line is outside
     the line segment, the closest point is a vertex. Return the
     square of the distance from the point to the closest point on the
     line */

double P_ClosestPntOnLine(double *pxyz, double *lxyz0, double *lxyz1, 
			  double eps, double *cxyz) {
  int i;
  double vec1[3], vec2[3], vec3[3];
  double d2, dp; 

  /* Vector from point 0 of line segment to point */
  for (i = 0; i < 3; i++)
    vec1[i] = pxyz[i]-lxyz0[i];

  /* Unit vector along line segment (d2 is length of line segment) */
  for (i = 0; i < 3; i++)
    vec2[i] = lxyz1[i]-lxyz0[i];

  d2 = VLen3(vec2);

  for (i = 0; i < 3; i++)
    vec2[i] /= d2;


  /* Projected length of vec1 along vec2 */

  dp = VDot3(vec1,vec2);


  if (dp < eps) { 
    /* projection of point is on vertex 0 or outside line segment
       beyond vertex 0 - so closest point is vertex 0 */

    for (i = 0; i < 3; i++)
      cxyz[i] = lxyz0[i];
    return VDot3(vec1,vec1);
  }
  else if (dp > (d2-eps)) { 
    /* projection of point is on vertex 1 or outside line segment
       beyond vertex 1 - so closest point is vertex 1 */

    for (i = 0; i < 3; i++)
      cxyz[i] = lxyz1[i];
    
    for (i = 0; i < 3; i++)
      vec3[i] = pxyz[i]-lxyz1[i];
    return VDot3(vec3,vec3);
  }
  else {
    /* Projection of point is inside line segment */

    for (i = 0; i < 3; i++)
      cxyz[i] = lxyz0[i] + dp*vec2[i];
    return ((pxyz[0]-cxyz[0])*(pxyz[0]-cxyz[0])+
	    (pxyz[1]-cxyz[1])*(pxyz[1]-cxyz[1])+
	    (pxyz[2]-cxyz[2])*(pxyz[2]-cxyz[2]));
  }

}

#ifdef __cplusplus
}
#endif

