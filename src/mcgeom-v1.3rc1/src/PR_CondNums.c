#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "MCGeom.h"

#ifdef __cplusplus
extern "C" {
#endif

  /* Assumption: Each vertex in the region has a maximum of 10 edges
     of that region coming into it */

void PR_CondNums(double (*rxyz)[3], int n, double *condnums, int (*reverts)[2], int ne) {
  int i, j, iv0, iv1, iv2, iv3, gen;
  double vec01[3], vec02[3], vec03[3], cpvec[3], a, b, vol6;
  static int first = 1, maxn = 0, *nnbr, (*vnbr)[10];

  static int tetidx[4][3] = {{1,2,3},{2,0,3},{0,1,3},{0,2,1}};
  static int prsmidx[6][3] = {{1,2,3},{2,0,4},{0,1,5},{5,4,0},{3,5,1},{4,3,2}};
  static int hexidx[8][3] = {{1,3,4},{2,0,5},{3,1,6},{0,2,7},{7,5,0},{4,6,1},
			     {5,7,2},{6,4,3}};

  /* n - Number of vertices, ne - Number of Edges
     Case 1: n == 4 ---------> Tet
     Case 2: n == 6 && ne == 9 or not specified (0) ---------> Tri Prism 
     case 4: n == 8 && ne == 12 or not specified (0) ---------> Hex
  */

  if (n == 4 || 
      (n == 6 && (ne == 0 || ne == 9)) || 
      (n == 8 && (ne == 0 || ne == 12)))
    gen = 0;
  else {
    gen = 1;
      
    if (first) {
      first = 0;
      vnbr = (int (*)[10]) malloc(n*sizeof(int [10]));
      nnbr = (int *) malloc(n*sizeof(int));
      maxn = n;
    }
    else {
      if (n > maxn) {
	vnbr = (int (*)[10]) realloc(vnbr,2*n*sizeof(int [10]));
	nnbr = (int *) realloc(nnbr,2*n*sizeof(int));
	maxn = 2*n;
      }
    }

    /* Collect neighbor vertex data */
    for (i = 0; i < ne; i++) {
      iv0 = reverts[i][0];
      iv1 = reverts[i][1];
      
      vnbr[iv0][nnbr[iv0]] = iv1;
      nnbr[iv0]++;
      
      vnbr[iv1][nnbr[iv1]] = iv0;
      nnbr[iv1]++;
      
      if (nnbr[iv0] == 10 || nnbr[iv1] == 10) {
	fprintf(stderr,"Limit of 10 neighbors reached\n");
      }
    }
  
    for (i = 0; i < n; i++)
      if (nnbr[i] != 3) {
	fprintf(stderr,"Vertex is not trivalent\n");
	return;
      }
  }


  for (i = 0; i < n; i++) {
    if (gen) {
      iv0 = i; iv1 = vnbr[i][0]; iv2 = vnbr[i][1]; iv3 = vnbr[i][2];
    }
    else {
      switch (n) {
      case 4:
	iv0 = i; iv1 = tetidx[i][0]; iv2 = tetidx[i][1]; iv3 = tetidx[i][2];
	break;
      case 6:
	iv0 = i; iv1 = prsmidx[i][0]; iv2 = prsmidx[i][1]; iv3 = prsmidx[i][2];
	break;
      case 8:
	iv0 = i; iv1 = hexidx[i][0]; iv2 = hexidx[i][1]; iv3 = hexidx[i][2];
	break;
      }
    }
      
    VDiff3(rxyz[iv1],rxyz[iv0],vec01);
    VDiff3(rxyz[iv2],rxyz[iv0],vec02);
    VDiff3(rxyz[iv3],rxyz[iv0],vec03);
    
    a = VLenSqr3(vec01) + VLenSqr3(vec02) + VLenSqr3(vec03);
    VCross3(vec02,vec03,cpvec);
    b = VLenSqr3(cpvec);
    VCross3(vec03,vec01,cpvec);
    b += VLenSqr3(cpvec);
    VCross3(vec01,vec02,cpvec);
    b += VLenSqr3(cpvec);
    vol6 = VDot3(cpvec,vec03);
    
    condnums[i] = sqrt(a*b)/vol6;
  }

}

#ifdef __cplusplus
}
#endif
