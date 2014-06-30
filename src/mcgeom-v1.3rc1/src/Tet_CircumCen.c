#include <stdio.h>
#include <math.h>
#include "MCGeom.h"


#ifdef __cplusplus
extern "C" {
#endif
 
#define det2(a11,a12,a21,a22) (a11*a22 - a21*a12)
#define det3(a11,a12,a13,a21,a22,a23,a31,a32,a33) \
 ((a31*((a12*a23)-(a22*a13)))-(a32*((a11*a23)-(a13*a21)))+(a33*((a11*a22)-\
 (a21*a12))))
 
/**************************************************************************/
 
void Tet_CircumCen(double (*xyz)[3], double *cen) {
 
/*  Computes the circumcenter of tetrahedron */	

  double dxyz[3][3], drhs[3], den, num;

 
  dxyz[0][0] = xyz[1][0]-xyz[0][0];
  dxyz[0][1] = xyz[1][1]-xyz[0][1];
  dxyz[0][2] = xyz[1][2]-xyz[0][2];

  dxyz[1][0] = xyz[2][0]-xyz[0][0];
  dxyz[1][1] = xyz[2][1]-xyz[0][1];
  dxyz[1][2] = xyz[2][2]-xyz[0][2];

  dxyz[2][0] = xyz[3][0]-xyz[0][0];
  dxyz[2][1] = xyz[3][1]-xyz[0][1];
  dxyz[2][2] = xyz[3][2]-xyz[0][2];

  den = det3(dxyz[0][0],dxyz[0][1],dxyz[0][2],
	     dxyz[1][0],dxyz[1][1],dxyz[1][2],
	     dxyz[2][0],dxyz[2][1],dxyz[2][2]);

  if (fabs(den) <= 1e-30) {
    cen[0] = cen[1] = cen[2] = 1e+30;
    fprintf(stderr,"Coplanar points; Cannot determine circumsphere!\n");
  }

  drhs[0] = 0.5*(dxyz[0][0]*dxyz[0][0] + dxyz[0][1]*dxyz[0][1] + dxyz[0][2]*dxyz[0][2]);
  drhs[1] = 0.5*(dxyz[1][0]*dxyz[1][0] + dxyz[1][1]*dxyz[1][1] + dxyz[1][2]*dxyz[1][2]);
  drhs[2] = 0.5*(dxyz[2][0]*dxyz[2][0] + dxyz[2][1]*dxyz[2][1] + dxyz[2][2]*dxyz[2][2]);
  

  num = det3(drhs[0],dxyz[0][1],dxyz[0][2],
	     drhs[1],dxyz[1][1],dxyz[1][2],
	     drhs[2],dxyz[2][1],dxyz[2][2]);

  cen[0] = xyz[0][0] + num/den;

  num = det3(dxyz[0][0],drhs[0],dxyz[0][2],
	     dxyz[1][0],drhs[1],dxyz[1][2],
	     dxyz[2][0],drhs[2],dxyz[2][2]);

  cen[1] = xyz[0][1] + num/den;

  num = det3(dxyz[0][0],dxyz[0][1],drhs[0],
	     dxyz[1][0],dxyz[1][1],drhs[1],
	     dxyz[2][0],dxyz[2][1],drhs[2]);

  cen[2] = xyz[0][2] + num/den;

}
 
 
#ifdef __cplusplus
}
#endif 
