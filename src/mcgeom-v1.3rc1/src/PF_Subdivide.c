#include <stdio.h>
#include <malloc.h>

#include "MCGeom.h"

#ifdef __cplusplus
extern "C" {
#endif

  /* Cut a polygon along a piecewise linear curve (PLC) into two pieces */
  /* The left part of the polygon is returned in pxyzl and the right
     part in pxyzr */
  /* Assumptions: 
     1) Planar polygon
     2) Piecewise linear curve lies in the plane of the polygon. 
     3) End points of piecewise linear curve lie on boundary of the polygon
     4) If polygon is non-convex or the interior points of PLC lie on
        the polygon boundary, the polygon may be split into multiple
        pieces. However, we are going to return 2 polygons only - some
        of them may be multiply connected or degenerate polygons.
*/
  
  int PF_Subdivide(int np, double (*pxyz)[3], int nlp, double (*lxyz)[3], 
		   int *ie, double eps, 
		   int *npl, double (*pxyzl)[3], int *ptagl,
		   int *npr, double (*pxyzr)[3], int *ptagr) {
    
    double dist2, dp, lbb[2][3], pcen[3], plbb[2][3], plxyz[2][3];
    double dist2_0, dist2_1;
    int disjoint, found, i, ip, j, k, k1, k2, kk1, kk2, rev = 0;
    static int first = 1;

    if (ie[0] < -1 || ie[0] > np-1)
      fprintf(stderr,"PF_Subdivide: Improper initialization of edge indices, ie\n");
    if (ie[1] < -1 || ie[1] > np-1)
      fprintf(stderr,"PF_Subdivide: Improper initialization of edge indices, ie\n");
    
    /* Find the edges on which the end points of the PLC lie */
 
    if (ie[0] == -1) {

      found = 0;
      for (ip = 0; ip < np; ip++) {
      
	for (i = 0; i < 3; i++) 
	  plxyz[0][i] = pxyz[ip][i];
	for (i = 0; i < 3; i++) 
	  plxyz[1][i] = pxyz[(ip+1)%np][i];
      
	if (P_OnLineSeg(lxyz[0],plxyz,eps)) {
	          
	  found = 1;

	  /* point is on this polygon edge */

	  ie[0] = ip;
      
	  /* Check if the point is coincident with the second vertex
	     of the edge in which case we have to pick the next edge
	     as the one that it is on. Simplifies the code to build
	     the sub-polygons later on */

	  if (PP_Dist2(lxyz[0],plxyz[1]) < eps*eps) 
	    ie[0] = (ip+1)%np;

	  break;
	}
      }

    }

    if (ie[1] == -1) {

      for (ip = 0; ip < np; ip++) {
      
	for (i = 0; i < 3; i++) 
	  plxyz[0][i] = pxyz[ip][i];
	for (i = 0; i < 3; i++) 
	  plxyz[1][i] = pxyz[(ip+1)%np][i];
      
	if (P_OnLineSeg(lxyz[nlp-1],plxyz,eps)) {
	          
	  found = 1;

	  /* point is on this polygon edge */

	  ie[1] = ip;
      
	  /* Check if the point is coincident with the second vertex
	     of the edge in which case we have to pick the next edge
	     as the one that it is on. Simplifies the code to build
	     the sub-polygons later on */

	  if (PP_Dist2(lxyz[nlp-1],plxyz[1]) < eps*eps)
	    ie[1] = (ip+1)%np;

	  break;
	}
      }

    }

    if (ie[0] == -1 || ie[1] == -1) {
      if (first) {
	first = 0;
	fprintf(stderr,
		"Could not locate one of the PLC end points on polygon\n");
      }
      return 0;
    }


    /* Now we have to build the polygons. Eventually, we have to break
       into multiple polygons if the interior points of the PLC lie on
       the boundary of the polygon */
    
    *npl = *npr = 0;

    kk1 = ie[0]; kk2 = ie[1];

    
    if (kk1 == kk2) {

      /* We have to do a special geometric evaluation, since we can't
	 tell which side is which for the interface */
      
      /* Distance from point 0 to lxyz[0] */

      dist2_0 = PP_Dist2(pxyz[kk1],lxyz[0]);

      /* Distance from point 0 to lxyz[nlp-1] */

      dist2_1 = PP_Dist2(pxyz[kk1],lxyz[nlp-1]);

      /* Is PLC in the same direction as the polygon points or
	 in the reverse direction - Will help us determine uniquely
	 which polygon is on the left side */

      rev =  (dist2_0 > dist2_1) ? 1 : 0;

    }

    
    /* Polygon to the left of the PLC */

    /* First add the polygon points */

    i = kk2+1;
    if (kk1 == kk2) {
      if (rev == 0) { 
	/* All polygon points will be in the left sub-polygon */
	while (i <= (kk1+np)) {
	  VCopy3(pxyzl[*npl],pxyz[i%np]);
	  ptagl[*npl] = 0;                 /* point of original polygon */
	  (*npl)++;
	  i++;
	}
      }  
      else {
	/* No polygon points will be in the right sub-polygon */
      }
    }
    else if (kk1 < kk2) {
      while (i <= (kk1+np)) {
	VCopy3(pxyzl[*npl],pxyz[i%np]);
	ptagl[*npl] = 0;                   /* point of original polygon */
	(*npl)++;
	i++;
      }
    }
    else {
      while (i <= kk1) {
	VCopy3(pxyzl[*npl],pxyz[i%np]);
	ptagl[*npl] = 0;                   /* point of original polygon */
	(*npl)++;
	i++;
      }
    }

    if ((*npl) == 0) { /* No polygon points were included */
      if (nlp == 2) {
	/* We must be grazing an edge. So ignore this "degenerate" polygon */
	*npl = 0;
      }
      else {
	/* we add all the points of the PLC */
	for (i = 0; i < nlp; i++) {
	  VCopy3(pxyzl[*npl],lxyz[i]);
	  ptagl[*npl] = 1;                 /* New point */
	  (*npl)++;
	}
      }
    }
    else {
      /* Then add the PLC points (for the start and end points
	 make sure that they don't coincide with existing polygon
	 points before adding them) */

      if (PP_Dist2(lxyz[0],pxyzl[*npl-1]) > eps*eps) {
	VCopy3(pxyzl[*npl],lxyz[0]);
	ptagl[*npl] = 1;                 /* New point */
	(*npl)++;
      }

      for (i = 1; i < nlp-1; i++) {
	VCopy3(pxyzl[*npl],lxyz[i]);
	ptagl[*npl] = 1;                 /* New point */
	(*npl)++;
      }

      if (PP_Dist2(lxyz[nlp-1],pxyzl[0]) > eps*eps) {
	VCopy3(pxyzl[*npl],lxyz[nlp-1]);
	ptagl[*npl] = 1;                 /* New point */
	(*npl)++;
      }
    }




    /* Polygon to the right of the PLC */

    /* First add polygon points */

    i = kk1+1;
    if (kk1 == kk2) {
      if (rev == 0) {
	/* No polygon points will be included in right sub-polygon */
      }
      else {
	/* All polygon points will be included in right sub-polygon */
	while (i <= kk2+np) {
	  VCopy3(pxyzr[*npr],pxyz[i%np]);
	  ptagr[*npr] = 0;                   /* point of original polygon */
	  (*npr)++;
	  i++;
	}
       }
    }
    else if (kk1 < kk2) {
      while (i <= kk2) {
	VCopy3(pxyzr[*npr],pxyz[i%np]);
	ptagr[*npr] = 0;                   /* point of original polygon */
	(*npr)++;
	i++;
      }
    }
    else {
      while (i <= kk2+np) {
	VCopy3(pxyzr[*npr],pxyz[i%np]);
	ptagr[*npr] = 0;                   /* point of original polygon */
	(*npr)++;
	i++;
      }
    }


    if ((*npr) == 0) { /* No polygon points were included */
      if (nlp == 2) {
	/* We must be grazing an edge. So ignore this "degenerate" polygon */
	*npr = 0;
      }
      else {
	/* Add all the points of the PLC in reverse order */
	for (i = nlp-1; i >= 0; i--) {
	  VCopy3(pxyzr[*npr],lxyz[i]);
	  ptagr[*npr] = 1;                  /* New point */
	  (*npr)++;
	}
      }
    }
    else {
      /* Then add the PLC points in reverse order (for the start
	 and end points make sure that they don't coincide with
	 existing polygon points before adding them) */

      if (PP_Dist2(lxyz[nlp-1],pxyzr[*npr-1]) > eps*eps) {
	VCopy3(pxyzr[*npr],lxyz[nlp-1]);
	ptagr[*npr] = 1;                   /* New point */
	(*npr)++;
      }

      for (i = nlp-2; i >= 1; i--) {
	VCopy3(pxyzr[*npr],lxyz[i]);
	ptagr[*npr] = 1;                   /* New point */
	(*npr)++;
      }
      
      if (PP_Dist2(lxyz[0],pxyzr[0]) > eps*eps) {
	VCopy3(pxyzr[*npr],lxyz[0]);
	ptagr[*npr] = 1;                   /* New point */
	(*npr)++;
      }
    }

    return 1;
  }


#ifdef __cplusplus
}
#endif
