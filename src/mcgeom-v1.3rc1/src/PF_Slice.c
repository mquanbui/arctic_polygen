#include <stdio.h>
#include <malloc.h>

#include "MCGeom.h"

#ifdef __cplusplus
extern "C" {
#endif

  /* Slice a *** CONVEX *** polygon by a line segment into two */
  /* The left part of the polygon is returned in pxyzl and the right
     part in pxyzr */
  /* Some aspects of this code assume that we are dealing with a
     planar polygon and that the line lies in the plane of the
     polygon. */
  
  int PF_Slice(int np, double (*pxyz)[3], double (*lxyz)[3], int *ie, 
	       double eps, 
	       int *npl, double (*pxyzl)[3], int *ptagl, 
	       int *npr, double (*pxyzr)[3], int *ptagr) {
    
    double dist2, maxdist2, dp, pcen[3], ixyz2[10][3], cpvec[3];
    double ldir[3], vec1[3], vec2[3], vec3[3], ixyz[10][3];
    int i, j, k, nint, nx, iedge[10], ie2[2], ii, jj, imax, jmax;

    
    /* Find the intersection points of this line with the polygon */

    X_LineFace(lxyz, np, pxyz, eps, &nx, ixyz, iedge);



    if (nx > 2) {

      /*
	If we got more than 2 intersection points it is most likely
	due to roundoff errors. Even though we try really hard,
	sometimes two nearly identical points (within a specified
	tolerance) are identified as being distinct. A less likely
	scenario is that the polygon is non-convex, but we do not
	expect it here.
	
	To deal with this, assuming a roundoff issue, we just pick the
	two intersection points that are the farthest apart. This
	assumes that the remaining points are really close to one of
	these two points.
      */

      imax = -1;
      jmax = -1;
      maxdist2 = -1.0e+16;
      
      for (ii = 0; ii < nx-1; ii++) {
	for (jj = ii+1; jj < nx; jj++) {
	  
	  dist2 = PP_Dist2(ixyz[ii],ixyz[jj]);
	  if (dist2 > maxdist2) {
	    imax = ii;
	    jmax = jj;
	    maxdist2 = dist2;
	  }
	}
      }
      
      if (imax != -1) {
	VCopy3(ixyz[0],ixyz[imax]);
	iedge[0] = iedge[imax];
      }
      if (jmax != -1) {	  
	VCopy3(ixyz[1],ixyz[jmax]);
	iedge[1] = iedge[jmax];
      }
      
      nx = 2;
    }



    if (nx == 2) { 
      /* First check if the two points are not identical */
      
      dist2 = PP_Dist2(ixyz[0],ixyz[1]);
      
      if (dist2 <  eps*eps) /* coincident intersection points */
	nx = 1;
      
    }


    /*** NOTE: EARLIER WE JUST RETURNED 0 if nx == 0. INTERFACE
	 RECONSTRUCTION ROUTINES WHICH USE THIS PROCEDURE HAVE NOT
	 BEEN TESTED WITH THE NEW CHANGES ****/



    if (nx == 0 || nx == 1) {

      PF_Center(np,pxyz,pcen);

      VDiff3(lxyz[1],lxyz[0],ldir);

      VDiff3(pcen,lxyz[0],vec2);

      VCross3(vec2,ldir,cpvec);

      if (cpvec[2] > 0) {
	/* polygon is to right of line */

	*npl = 0;

	*npr = np;
	for (j = 0; j < np; j++)
	  for (i = 0; i < 3; i++)
	    pxyzr[j][i] = pxyz[j][i];
      }
      else {
	/* polygon is to left of line */
	
	*npl = np;
	for (j = 0; j < np; j++)
	  for (i = 0; i < 3; i++)
	    pxyzl[j][i] = pxyz[j][i];
	
	*npr = 0;
      }

      for (i = 0; i < 3; i++)
	lxyz[0][i] = lxyz[1][i] = ixyz[0][i];

      ie[0] = iedge[0];
    }
    else {

      ie[0] = iedge[0];
      ie[1] = iedge[1];


      /* Regular intersection */


      /*
	Make sure the line segment passed to PF_Subdivide is oriented
	consistently with the line segment. We must also copy ixyz
	into ixyz2 before passing into PF_Subdivide because
	PF_Subdivide is expecting an array of dimensions (2,ndim)
      */

      VDiff3(ixyz[1],ixyz[0],vec1);
      VDiff3(ixyz[1],lxyz[0],ldir);

      dp = VDot3(vec1,ldir);

      if (dp <= 0) {

	/* Switch the first and second intersection points around */

	VCopy3(ixyz2[0],ixyz[1]);
	VCopy3(ixyz2[1],ixyz[0]);

	ie2[0] = ie[1];
	ie2[1] = ie[0];

      }
      else {

	VCopy3(ixyz2[0],ixyz[0]);
	VCopy3(ixyz2[1],ixyz[1]);

	ie2[0] = ie[0];
	ie2[1] = ie[1];

      }

      /* Subdivide the polygon */

      PF_Subdivide(np, pxyz, 2, ixyz2, ie2, eps, npl, pxyzl, ptagl,
		   npr, pxyzr, ptagr);


      /* Prepare info about subdividing line for sending out */
      
      VCopy3(lxyz[0],ixyz2[0]);
      VCopy3(lxyz[1],ixyz2[1]);

      ie[0] = ie2[0];
      ie[1] = ie2[1];

    }


    
    return 1;
  }


#ifdef __cplusplus
}
#endif
