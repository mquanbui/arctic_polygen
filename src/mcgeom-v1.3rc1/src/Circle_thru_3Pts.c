#include <math.h>
#include "MCGeom.h"

void Circle_thru_3Pts(double *xyz0, double *xyz1, double *xyz2,
		       double cen[3], double *R) {

  double A[2][2], B[2], den;

  A[0][0] = 2*(xyz0[0]-xyz1[0]); 
  A[0][1] = 2*(xyz0[1]-xyz1[1]);
  A[1][0] = 2*(xyz0[0]-xyz2[0]);
  A[1][1] = 2*(xyz0[1]-xyz2[1]);
  B[0] = (xyz0[0]*xyz0[0] + xyz0[1]*xyz0[1] -
	  xyz1[0]*xyz1[0] - xyz1[1]*xyz1[1]);
  B[1] = (xyz0[0]*xyz0[0] + xyz0[1]*xyz0[1] -
	  xyz2[0]*xyz2[0] - xyz2[1]*xyz2[1]);
  
  den = A[1][0]*A[0][1]-A[0][0]*A[1][1];
  
  cen[0] = (A[0][1]*B[1]-A[1][1]*B[0])/den;
  cen[1] = (A[1][0]*B[0]-A[0][0]*B[1])/den;
  cen[2] = 0.0;

  *R = sqrt((xyz0[0]-cen[0])*(xyz0[0]-cen[0])+(xyz0[1]-cen[1])*(xyz0[1]-cen[1]));

}
