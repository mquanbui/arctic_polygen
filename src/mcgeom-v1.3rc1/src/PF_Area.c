#include <math.h>
#include "MCGeom.h"

/* This works only for 2D. */

double PF_Area(int n, double (*xyz)[3]) {
  int i, j;
  double pcen[3], txyz[3][3], area = 0.0;
  
  if (n < 3)
    return 0;
  
  /* This is based on Green's Theorem in the Plane - Works for all
     3D polygons 
       
     Area = 0.5*Sum_over_i(a_i);
     a_i = x(i)*y(i+1)-x(i+1)*y(i);
       
     However, if the coordinates are very large, then a*b-c*d can
     result in roundoff error. To improve accuracy, we will make
     all coordinates relative to x0,y0. But what if xi is very close 
     to x0? Then xi-x0 will also generate high error. Which is better?

*/

  if (fabs(xyz[0][0]) > 1000.0 || fabs(xyz[0][1]) > 1000.0) { 
    for (i = 0; i < n; i++)
      area += ((xyz[i][0]-xyz[0][0])*(xyz[(i+1)%n][1]-xyz[0][1]) - 
               (xyz[(i+1)%n][0]-xyz[0][0])*(xyz[i][1]-xyz[0][1]));
  }
  else {
    for (i = 0; i < n; i++)
      area += ((xyz[i][0])*(xyz[(i+1)%n][1]) - 
               (xyz[(i+1)%n][0])*(xyz[i][1]));
  }

  area = 0.5*area;
    
  return area;

}

  
