#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "MCGeom.h"


/* Procedure to compute minimum and maximum dihedral angles of general
   polyhedron */

/* If the polyhedron is of a specific type, say a tet, pyramid, prism
   or a hex, then an ordering of vertices is implied and the
   coordinates are to be provide according to that ordering. The only
   other argument required in these cases is the output "cosangs"
   argument. One can call 

   PR_MinMaxAngles(rxyz, n, &minang, &maxang, 0, 0, 0); 

   */

/* For a general polyhedron, connectivity info has to be provided to
   the routine so as to determine the various faces. 'nf' is the
   number of faces of the polyhedron, the array 'nfv' has the number
   of vertices of each face and rfverts has the indices of the face
   vertices. */

void PR_MinMaxAngles(double (*rxyz)[3], int n, double *minang, double *maxang, 
		     int **rfverts, int *nfv, int nf) {
  int i, j, k, l=0, found, nedges, iv0, iv1, iv2;
  int feverts[2], (*everts)[2];
  double vec01[3], vec02[3], len1_sqr, len2_sqr, mincos, maxcos;
  double dp, n1[3], n2[3];

  mincos = 99.0;
  maxcos = -99.0;

#ifdef DEBUG
  if (n < 4)
    fprintf(stderr,"PR_CosAngles: Input is not a polyhedron\n");
#endif

  /* ncos needs to be specified only in the case of  general polygons */

  switch (n) {
  case 4:
    Tet_MinMaxAngles(rxyz, minang, maxang);
    break;
  case 5:
    break;
  case 6:
    break;
  case 8:
    Hex_MinMaxAngles(rxyz, minang, maxang);
    break;
  default: /* General polyhedron */
    everts = (int (*)[2]) malloc(10*nf*sizeof(int [2]));
    nedges = 0;
    for (i = 0; i < nf; i++) {
      for (j = 0; j < nfv[i]; j++) {
	feverts[0] = rfverts[i][j];
	feverts[1] = rfverts[i][(j+1)%nfv[i]];

	/* Check if edge has been processed already */
	if (i != 0) {
	  found = 0;
	  for (k = 0; k < nedges; k++) {
	    if (((everts[k][0] == feverts[0]) && (everts[k][1] == feverts[1])) ||
		((everts[k][0] == feverts[1]) && (everts[k][1] == feverts[0]))) {
	      found = 1;
	      break;
	    }
	  }
	  if (found)
	    continue;
	}

	everts[nedges][0] = feverts[0];
	everts[nedges][1] = feverts[1];


	/* Pick one vertex adjacent to the edge to form a triangle or plane */
	iv0 = rfverts[i][(j+1)%nfv[i]];
	iv1 = rfverts[i][(j+2)%nfv[i]];
	iv2 = rfverts[i][j];
	VDiff3(rxyz[iv1], rxyz[iv0], vec01);
	VDiff3(rxyz[iv2], rxyz[iv0], vec02);
	VCross3(vec01,vec02,n1);
	len1_sqr = VLenSqr3(n1);

	/* Find another face using this edge */
	found = 0; k = 0;
	while (!found && k < nf) {
	  if (k == i)
	    continue;
	  for (l = 0; l < nfv[k]; l++) {
	    if (((rfverts[k][l] == feverts[0]) && 
		 (rfverts[k][(l+1)%nfv[k]] == feverts[1])) || 
		((rfverts[k][l] == feverts[1]) && 
		 (rfverts[k][(l+1)%nfv[k]] == feverts[0]))) {
	      found = 1;
	      break;
	    }
	  }
	  k++;
	}
	
	if (!found) {
	  fprintf(stderr,"Cannot find another face using edge\n");
	  nedges++;
	  continue;
	}


	/* Get normal for the adjacent face */
	iv0 = rfverts[k][(l+1)%nfv[k]];
	iv1 = rfverts[k][(l+2)%nfv[k]];
	iv2 = rfverts[k][l];
	VDiff3(rxyz[iv1], rxyz[iv0], vec01);
	VDiff3(rxyz[iv2], rxyz[iv0], vec02);
	VCross3(vec01,vec02,n2);
	len2_sqr = VLenSqr3(n2);


	dp = VDot3(n1,n2);
	dp = sqrt((dp*dp)/(len1_sqr*len2_sqr));
	if (dp < mincos)
	  mincos = dp;
	if (dp > maxcos)
	  maxcos = dp;
	nedges++;
      }
    }
    free(everts);

    *minang = acos(maxcos);
    *maxang = acos(mincos);
  }
   
}
