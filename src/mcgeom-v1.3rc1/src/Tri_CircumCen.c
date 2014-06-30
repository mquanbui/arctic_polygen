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
 
void Tri_CircumCen(double (*xyz)[3], double *cen) {
 
/*  Computes the circumcenter of three dimensional triangle */	

  double A[3][3],B[3];
  double deta, dist[3];
 
 
  A[0][0] =  xyz[2][0] - xyz[0][0];
  A[0][1] =  xyz[2][1] - xyz[0][1];
  A[0][2] =  xyz[2][2] - xyz[0][2];
 
  A[1][0] =  xyz[1][0] - xyz[0][0];
  A[1][1] =  xyz[1][1] - xyz[0][1];
  A[1][2] =  xyz[1][2] - xyz[0][2];
 
  A[2][0] =  det2(A[0][1],A[0][2],A[1][1],A[1][2]);
  A[2][1] = -det2(A[0][0],A[0][2],A[1][0],A[1][2]);
  A[2][2] =  det2(A[0][0],A[0][1],A[1][0],A[1][1]);
 
 
  B[0] = 0.5*(A[0][0]*(xyz[0][0]+xyz[2][0])+
	      A[0][1]*(xyz[0][1]+xyz[2][1])+
	      A[0][2]*(xyz[0][2]+xyz[2][2]));
 
 
  B[1] = 0.5*(A[1][0]*(xyz[0][0]+xyz[1][0])+
	      A[1][1]*(xyz[0][1]+xyz[1][1])+
	      A[1][2]*(xyz[0][2]+xyz[1][2]));
 
 
  B[2] = (A[2][0]*(xyz[0][0]) +
	  A[2][1]*(xyz[0][1]) +
	  A[2][2]*(xyz[0][2]));
 
 
  /* Now solve the system */
  deta = det3(A[0][0],A[0][1],A[0][2],
	      A[1][0],A[1][1],A[1][2],
	      A[2][0],A[2][1],A[2][2])+1e-30;
 
  if (fabs(deta) < 1.0e-30) {
    fprintf(stderr,"ERROR: Singular matrix--Compute Triangle center--%lf.\n",deta);
  }
 
 
  cen[0] = det3(B[0],A[0][1],A[0][2],
	       B[1],A[1][1],A[1][2],
	       B[2],A[2][1],A[2][2])/deta;
 
 
  cen[1] = det3(A[0][0],B[0],A[0][2],
	       A[1][0],B[1],A[1][2],
	       A[2][0],B[2],A[2][2])/deta;
 
  cen[2] = det3(A[0][0],A[0][1],B[0],
	       A[1][0],A[1][1],B[1],
	       A[2][0],A[2][1],B[2])/deta;
 
}
 
 
#ifdef __cplusplus
}
#endif 
