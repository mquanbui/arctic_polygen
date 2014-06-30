#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "MCGeom.h"

/*
   Rotate a point 'p' by angle 'theta' around an arbitrarily oriented
   axis 'ra' at point 'rp'.
   Return the rotated point, 'p1'.  
   Positive angles are anticlockwise looking down the axis towards the origin.

   Paul Borke, Ronald Goldman,
   http://astronomy.swin.edu.au/~pbourke/geometry/rotate/source1.c
   with modification to add 'rotation point rp'
*/

#ifdef __cplusplus
extern "C" {
#endif

void P_RotateGen(double p[], double theta, double rp[], double ra[], 
		   double p1[]) {
   double costheta,sintheta;

   /* Assume ra is normalized */
   /* normVt(ra,ra); */

   costheta = cos(theta);
   sintheta = sin(theta);

   p[0] -= rp[0]; p[1] -= rp[1]; p[2] -= rp[2];

   p1[0] = p1[1] = p1[2] = 0.0;

   p1[0] += (costheta + (1 - costheta) * ra[0] * ra[0]) * p[0];
   p1[0] += ((1 - costheta) * ra[0] * ra[1] - ra[2] * sintheta) * p[1];
   p1[0] += ((1 - costheta) * ra[0] * ra[2] + ra[1] * sintheta) * p[2];

   p1[1] += ((1 - costheta) * ra[0] * ra[1] + ra[2] * sintheta) * p[0];
   p1[1] += (costheta + (1 - costheta) * ra[1] * ra[1]) * p[1];
   p1[1] += ((1 - costheta) * ra[1] * ra[2] - ra[0] * sintheta) * p[2];

   p1[2] += ((1 - costheta) * ra[0] * ra[2] - ra[1] * sintheta) * p[0];
   p1[2] += ((1 - costheta) * ra[1] * ra[2] + ra[0] * sintheta) * p[1];
   p1[2] += (costheta + (1 - costheta) * ra[2] * ra[2]) * p[2];

   p[0] += rp[0]; p[1] += rp[1]; p[2] += rp[2];
   p1[0] += rp[0]; p1[1] += rp[1]; p1[2] += rp[2];
}


#ifdef __cplusplus
}
#endif
