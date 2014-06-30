#include <stdio.h>
#include <math.h>
#include "MCGeom.h"


#ifdef __cplusplus
extern "C" {
#endif

/**************************************************************************/
/* area (a,b,c) = length of AxB, A= b-a, B = c-a */
 
double Tri_Area(double (*xyz)[3]) {

  double A0, A1, A2, B0,B1,B2, C0,C1,C2;
  double D1,D2,D3,area,area2,L0,L1,L2;
 
  A0 = xyz[1][0] - xyz[0][0];
  A1 = xyz[1][1] - xyz[0][1];
  A2 = xyz[1][2] - xyz[0][2];
 
  B0 = xyz[2][0] - xyz[0][0];
  B1 = xyz[2][1] - xyz[0][1];
  B2 = xyz[2][2] - xyz[0][2];
 
  D1 = A1*B2-A2*B1;
  D1 = D1*D1;
 
  D2 = A2*B0 - A0*B2;
  D2 = D2*D2;
 
  D3 = A0*B1 - A1*B0;
  D3 = D3*D3;

  area = sqrt(D1+D2+D3)/2.0;

  return area;
}
 
#ifdef __cplusplus
}
#endif 
