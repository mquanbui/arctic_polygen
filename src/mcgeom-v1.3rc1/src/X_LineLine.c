#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "MCGeom.h"

/* Intersect two line segments (including collinear cases). Return 1
   if they intersect, 0 if they don't. Return the number of
   intersection points in nx, intersection points in the xpt array */

int X_LineLine(double (*La)[3], double (*Lb)[3], double ptol, int *nx, double (*xpt)[3]) {
  int i, j, k, iop[4], xtrachk, imax, itmp, coplanar=0, W1W2_len_calc=0;
  double LLa[2][3], LLb[2][3];
  double A[3], B[3], W1[3], W2[3], W3[3], V1[3], V2[3], U1[3], U2[3], X1[3];
  double A_len, A_len2, B_len2, W1_len, W2_len;
  double rtmp, s1, s2, ldist[4];
  double sa, sb, dpu, dpv;
  double dist, dist2, ptol2;
  double MACH_EPS=1.0e-16;
#ifdef DEBUG
  double pt[3];
#endif
  

  ptol2 = ptol*ptol;

  xtrachk = 0;

  VCopy3(LLa[0],La[0]);
  VCopy3(LLa[1],La[1]);
  VCopy3(LLb[0],Lb[0]);
  VCopy3(LLb[1],Lb[1]);
  

  *nx = 0; 


  VDiff3(LLa[1],LLa[0],A);

  VDiff3(LLb[1],LLb[0],B);



  /* We are going to make the shorter of the two segments, segment A. */

  A_len2 = 0.0;
  B_len2 = 0.0;

  for (k = 0; k < 3; k++) {
    A_len2 = A_len2 + A[k]*A[k];
    B_len2 = B_len2 + B[k]*B[k];
  }

  if (A_len2 > 1.2*B_len2) {

    /* Switch the segments around */

    VCopy3(LLa[0],Lb[0]);
    VCopy3(LLa[1],Lb[1]);
    VCopy3(LLb[0],La[0]);
    VCopy3(LLb[1],La[1]);

    VDiff3(LLa[1],LLa[0],A);

    VDiff3(LLb[1],LLb[0],B);

    rtmp = A_len2;
    A_len2 = B_len2;
    B_len2 = rtmp;

  }


  VDiff3(LLb[0],LLa[0],W1);
  
  VDiff3(LLb[1],LLa[1],W2);
  

  /* Check for intersection - Compute the signed distance (scaled by
     length of A or B for efficiency) of each end point with respect
     to the other line. The two line segments intersect only If the
     distance of LLa(1) and LLa(2) from line LLb are of opposite signs
     AND the distance of LLb(1) and LLb(2) from LLa are of opposite
     signs. */

  VCross3(W1,A,V1);         /* V1 = W1xA 
				||V1|| = ||A||*(distance of LLb[0] to LLa) */

  VCross3(W2,A,V2);         /* V2 = W2xA
				||V2|| = ||A||*(distance of LLb(2) to LLa) */


  dpv = VDot3(V1,V2);
  if (dpv > 10*MACH_EPS) return 0;

  /* If the end point of LLb is on or nearly on LLa, then we have to
     do extra checks */

  if (dpv > -10*MACH_EPS) {

    /* Do the real check to try and remove the ambiguity. Compute v1 =
       V1/(||W1||*||A||) and v2 = V2/(||W2||*||A||) and check the dot
       product. We computed dp = V1.V2, so to get v1.v2 we divide by
       A_Len2 = ||A||*||A|| times ||V1||*||V2|| */

    W1_len = VLen3(W1);
    W2_len = VLen3(W2);
    W1W2_len_calc = 1;

    dpv = dpv/(W1_len*W2_len*A_len2);

    if (dpv > 10*MACH_EPS) return 0;

    /* Just to be sure we will do some checks later */

    xtrachk = 1;
  }


  /* Since we are only interested in checking if U1.U2 < 0, we will not
     bother with negating W1 in the next two computations */

  VCross3(W1,B,U1);         /* U1 = -W1xB 
				||U1|| = ||B||*(distance of LLa(1) to LLb) */

  VCross3(W2,B,U2);         /* U2 = -W2xB
				||U2|| = ||B||*(distance of LLa(2) to LLb) */


  dpu = VDot3(U1,U2);
  if (dpu > 10*MACH_EPS) return 0;

  /* If the end point of LLa is on or nearly on LLb, then we have to do
     extra checks */

  if (dpu > -10*MACH_EPS) {

    /* Do the real check to try and remove the ambiguity. Compute u1 =
       U1/(||W1||*||B||) and u2 = U2/||W2||.||B|| and check the dot
       product. We computed dp = U1.U2, so to get u1.u2 we divide by
       B_Len2 = ||B||*||B|| times ||W1|*||W2||*/

    if (!W1W2_len_calc) { /* if we haven't calculated lengths of W1, W2 */
      W1_len = VLen3(W1);
      W2_len = VLen3(W2);
      W1W2_len_calc = 1;
    }

    dpu = dpu/(W1_len*W2_len*B_len2);

    if (dpu > 10*MACH_EPS) return 0;

    /* Just to be sure we will do some checks later */

    xtrachk = 1;
  }
      

  /* MAKE SURE THE LINES ARE COPLANAR - THIS IS AN EXPENSIVE CHECK SO
     LETS ELIMINATE THE OBVIOUS 2D CASE WHERE ALL THE Z-COORDINATES ARE
     THE SAME */

  coplanar = 0;

  if ((fabs(LLa[0][2]-LLa[1][2]) < MACH_EPS) && 
      (fabs(LLa[0][2]-LLb[0][2]) < MACH_EPS) && 
      (fabs(LLa[0][2]-LLb[1][2]) < MACH_EPS))

    coplanar = 1;

  else {    

    /* Have to do an explicit check for coplanarity

    The following checks if the distance of the first point of the
    segment LLa (or equivalently of vector A) to the plane formed by
    the segment LLb and the second point of vector A is ZERO. V2 is
    normal to this plane

    */

    double dpb, B_len;


    dpb = VDot3(B,V1);

    if (!W1W2_len_calc) {
      W1_len = VLen3(W1);
      W2_len = VLen3(W2);      
    }
    B_len = sqrt(B_len2);

    dpb = dpb/(B_len*W1_len);

    if (fabs(dpb) < 10*MACH_EPS) coplanar = 1;
				   
  }

  if (!coplanar) return 0;
      



  /* Compute real unsigned distances of LLa(1) and LLa(2) from LLb */

  A_len = sqrt(A_len2);

  s1 = VLen3(V1);
  s1 = s1/A_len;

  s2 = VLen3(V2);
  s2 = s2/A_len;


  /* Check for collinearity first */

  if (s1 <= ptol && s2 <= ptol) {

    /* Vectors are collinear  */

    /* Order the vertices from minimum to maximum distance w.r.t. some
       reference point on the line

       Index assignment LLa(1) = 1, LLa(2) = 2, LLb(1) = 3, LLb(2) = 4
    */

    for (i = 0; i < 4; i++)
      iop[i] = i;


    /* Distance of LLa[0] to itself */
    ldist[0] = 0.0;

    /* Signed Distance of LLa(2) to LLa(1) along A */
    ldist[1] = VLen3(A);

    /* Signed Distance of LLa[0] to LLb[0] - W1.A gives distance from
       LLa[0] to LLb[0] scaled by length of A (which is ldist[1]) */

    ldist[2] = VDot3(W1,A);
    ldist[2] = ldist[2]/ldist[1];

    /* Signed Distance of LLa[0] to LLb[1] - W3.A gives distance from LLa[0]
       to LLb[1] scaled by length of A (which is ldist[1]) */

    VDiff3(LLb[1],LLa[0],W3);

    ldist[3] = VDot3(W3,A);
    ldist[3] = ldist[3]/ldist[1];


    for (i = 0; i < 3; i++)
      for (j = i+1; j < 4; j++)
	if (ldist[i] > ldist[j]) {
	  itmp = iop[i];
	  iop[i] = iop[j];
	  iop[j] = itmp;

	  rtmp = ldist[i];
	  ldist[i] = ldist[j];
	  ldist[j] = rtmp;
	}


    /* Check for the collinear but disjoint case - if iop[0] and
       iop[1] belong to the same segment (<= 2 or > 2), and, iop[2]
       and iop[3] belong to the same segment (> 2 or <=2), then they
       are disjoint */

    if (((iop[0] <= 2) && (iop[1] <= 2) && (iop[2] > 2) && (iop[3] > 2)) ||
	((iop[0] > 2) && (iop[1] > 2) && (iop[2] <= 2) && (iop[3] <= 2)))
      return 0;


    /* Check the middle two points in the ordered set - They are the
       only ones that can represent the points of intersection of the
       coincident lines. This will catch coincident end points as
       well */

    for (i = 1; i <= 2; i++)
      if ((iop[i] <= 1)) { 	/* Point of line a */
	VCopy3(xpt[*nx],LLa[iop[i]]);
	*nx = *nx + 1;
      }
      else {  /* Point of line b */
	VCopy3(xpt[*nx],LLb[iop[i]-2]);
	*nx = *nx + 1;
      }

  }
  else if (dpv <= 10.0*MACH_EPS && dpu <= 10.0*MACH_EPS) {
    /* MOST RECENT CHANGE ^^^^^^^^^^ UNTESTED 08/30/08 */

    /* Find intersection point of _lines_ (not the segments)
       parametrically { check to see if the parameters are between 0
       and 1

       Formula is unfortunately unstable when (s1+s2) is close to 0.0

       HELPS THAT WE TAKE THE LONGER SEGMENT TO BE LINE b, BECAUSE THE
       DISTANCE OF THE END POINTS OF b TO LINE a, WILL BE LARGER THAN
       THE DISTANCE OF THE END POINTS OF a TO LINE b

    */


    sb =  s1/(s1+s2);

    if ((sb < -1e+01*MACH_EPS) || ((sb-1.0) > 1e+01*MACH_EPS)) {
     
      /* intersection is outside of line segment LLb */
      return 0;

    }


    /* Find the intersection point */

    for (k = 0; k < 3; k++) {
      xpt[0][k] = (1.0-sb)*LLb[0][k] + sb*LLb[1][k];
    }
    *nx = 1;


    if (sb < 0.0) {
      /* If we do accept a vertex that is below the bounds, Make sure
	 the point is within tolerance of the first vertex */
	
      dist2 = PP_Dist2(xpt[0],LLb[0]);
      if (dist2 > ptol2) {
	*nx = 0;
	return 0;
      }
      
    }
      
    if (sb > 1.0) {
      /* If we do accept a vertex that is above the bounds, Make sure
	 the point is within tolerance of the second vertex */
      
      dist2 = PP_Dist2(xpt[0],LLb[1]);
      if (dist2 > ptol2) {
	*nx = 0;
	return 0;
      }
      
    }


    if (xtrachk) {

      /* If there was ambiguity about where this is a true
	 intersection or not, do some extra checks */

      VDiff3(xpt[0],LLa[0],X1);
            
      dist = VDot3(X1,A);
      dist = dist/A_len;
                         /*  We have to divide by A_len because
                             to get the true distance of X1 from LLa(1)
                             we should have dotted X1 with a unit vector
                             along A */

      
      /* Compute parameter of intersection point on line segment LLa */

      sa = dist/A_len;


      if (sa < 0.0) {

	/* If we do accept a vertex that is below the bounds, make sure the
	   point is within tolerance of the first vertex */

	if (fabs(dist) > ptol) {
	  *nx = 0;
	  return 0;
	}
      }


      if (sa > 1.0) {
	/* If we do accept a vertex that is above the bounds, make sure the
	   point is within tolerance of the second vertex */

	dist2 = PP_Dist2(xpt[0],LLa[1]);
	if (dist2 > ptol2) {
	  *nx = 0;
	  return 0;
	}
      }

    }

  }
  
  return 1;
}
