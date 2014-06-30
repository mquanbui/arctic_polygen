#include <stdio.h>
#include <malloc.h>

#include "MCGeom.h"

#ifdef __cplusplus
extern "C" {
#endif

  /* Find points at which a line segment intersects a polygon */

  /* Some aspects of this code assume that we are dealing with a
     planar polygon in the XY plane and that the line lies in the
     plane of the polygon. */

  /* Procedure also returns the edge indices where the intersection
     was found. If an intersection occurs at a vertex, then the
     procedure picks the edge for which this vertex is the first
     vertex */

  /* We are not checking coplanarity here */
  
  int X_LineFace(double (*lxyz)[3], int np, double (*pxyz)[3], double ptol, 
		 int *nx, double (*ixyz)[3], int *ie) {

    int MAXINTPTS=10;
    double dist2, dp, lbb[2][3], pcen[3], plxyz[2][3];
    double intxyz[4][3];
    int status, nint, ip, i, j, k, nxv, skip;


    /* First check if __both__ end points of the intersecting line
       segment lie on an edge of the polygon, in which case, the end
       points are returned as the intersection points. If only one end
       point lies on a polygon edge, then proceed as usual */

    *nx = 0;
    for (i = 0; i < np; i++) {

      VCopy3(plxyz[0],pxyz[i]);
      VCopy3(plxyz[1],pxyz[(i+1)%np]);

      status = P_OnLineSeg(lxyz[0],plxyz,ptol);

      if (status) {

	status = P_OnLineSeg(lxyz[1],plxyz,ptol);

	if (status) {
	  
	  *nx = 2;
	  VCopy3(ixyz[0],lxyz[0]);
	  VCopy3(ixyz[1],lxyz[1]);
	 
	  return 1;
	}

      }
    }


    /* Also check explicitly if the vertices of the polygon are on
       the line */

    nxv = 0;
    *nx = 0;
    for (i = 0; i < np; i++) {

      status = P_OnLineSeg(pxyz[i],lxyz,ptol);

      if (status) {

	VCopy3(ixyz[*nx],pxyz[i]);

	ie[*nx] = i;

	nxv++;
	(*nx)++;

      }
      
    }
      
    
    for (i = 0; i < np; i++) {

      /* Skip this edge if one of its vertices is already detected as
	 being on the intersecting line */

      skip = 0;
      for (j = 0; j < nxv; j++) {
	if (i == ie[j] || (i+1)%np == ie[j]) {
	  skip = 1;
	  break;
	}
      }
      

      if (skip) continue;


      VCopy3(plxyz[0],pxyz[i]);
      VCopy3(plxyz[1],pxyz[(i+1)%np]);

      X_LineLine(lxyz, plxyz, ptol, &nint, intxyz);
      
      if (!nint) continue;
      

      for (j = 0; j < nint; j++) {

	if (*nx < MAXINTPTS) {

	  for (k = 0; k < 3; k++)
	    ixyz[*nx][k] = intxyz[j][k];
	  ie[*nx] = i;

	  (*nx)++;

	}

      }

    }
    

    return 1;
  }


#ifdef __cplusplus
}
#endif
