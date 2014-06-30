#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "MCGeom.h"

/* Intersect a line segment and a circle. Return 1 if they intersect,
   0 if they don't. Return the number of intersection points in nx,
   intersection points in the xpt array */

int X_LineCircle(double (*La)[3], double X0[3], double R, double ptol, int *nx, double (*xpt)[3]) {
  int i, j, k, soln, sgn;
  double lla[2][3], aq, bq, cq;
  double xa, xb, ya, yb, m, c;
  double MACH_EPS=1.0e-16, disc, expr;
#ifdef DEBUG
  double pt[3];
#endif

  if (fabs(La[0][2]-La[1][2]) > ptol || fabs(La[0][2]-X0[2]) > ptol) {
    fprintf(stderr,
	    "Line-circle intersection only implemented for the z=const case\n");
    return 0;
  }

  /* Translate line segment coordinates to coordinate system centered
     at the circle center */

  for (k = 0; k < 3; k++) {
    lla[0][k] = La[0][k]-X0[k];
    lla[1][k] = La[1][k]-X0[k];
  }

  *nx = 0; 
  soln = 1;

  if (fabs(lla[0][0]-lla[1][0]) > MACH_EPS) {

    /* If the line is in general orientation or horizontal do this */

    m = (lla[1][1]-lla[0][1])/(lla[1][0]-lla[0][0]);
    c = lla[0][1]-m*lla[0][0];

    aq = 1+m*m;   /* Cannot be zero */
    bq = 2*m*c;
    cq = c*c-R*R;

    sgn = (bq >= 0) ? 1 : -1;
    
    disc = bq*bq-4*aq*cq;
    if (disc >= 0.0) {
      
      expr = -bq-sgn*sqrt(disc);

      if (fabs(expr) > MACH_EPS) {      
	xa = (2*cq)/expr;
	ya = m*xa + c;
	xb = expr/(2*aq);  /* in this case, aq is never 0 */
	yb = m*xb + c;
      }
      else {	
	xa = xb = 0.0;
	yb = c;
      }

    }
    else
      soln = 0;
  }
  else {

    /* If the line is vertical, do this */

    xa = xb = lla[0][0];
    if (R*R-xa*xa >= 0.0) {
      ya =  sqrt(R*R-xa*xa);
      yb = -sqrt(R*R-xb*xb);
    }
    else
      soln = 0;
  }

  /* Check which of the two solutions are valid (if any) */

  if (soln) {
    if ((lla[0][0] <= xa && xa <= lla[1][0]) ||
	(lla[1][0] <= xa && xa <= lla[0][0])) {
      if ((lla[0][1] <= ya && ya <= lla[1][1]) ||
	  (lla[1][1] <= ya && ya <= lla[0][1])) {

	xpt[*nx][0] = xa + X0[0]; 
	xpt[*nx][1] = ya + X0[1];
	xpt[*nx][2] = X0[2];
	(*nx)++;

      }
    }

    /* Process second point only if it is distinguishable from the first */

    if ((xa-xb)*(xa-xb)+(ya-yb)*(ya-yb) > ptol*ptol) {
      if ((lla[0][0] <= xb && xb <= lla[1][0]) ||
	  (lla[1][0] <= xb && xb <= lla[0][0])) {
	if ((lla[0][1] <= yb && yb <= lla[1][1]) ||
	    (lla[1][1] <= yb && yb <= lla[0][1])) {
	  
	  xpt[*nx][0] = xb + X0[0]; 
	  xpt[*nx][1] = yb + X0[1];
	  xpt[*nx][2] = X0[2];
	  (*nx)++;
	  
	}
      }
    }

  }

  
  return (*nx ? 1 : 0);
}
