#include "MCGeom.h"

/* Centroid of polygonal face - calculated by subdividing polygon into
   triangles and doing an area weighted sum of the triangle
   centroids */

void PF_Centroid(int n, double (*xyz)[3], double *pcen) {
  int i, j;
  double tcen[200][3], tarea[200];
  double vec0[3], vec1[3], cpvec[3], areasum=0.0;
  
  for (j = 0; j < 3; j++)
    pcen[j] = 0.0;
    
  for (i = 1; i < n-1; i++) {

    /* centroid of triangle (0,i,i+1) */

    for (j = 0; j < 3; j++)
      tcen[i][j] = (xyz[0][j]+xyz[i][j]+xyz[i+1][j])/3.0;

    /* area of triangle */

    VDiff3(xyz[i],xyz[0],vec0);
    VDiff3(xyz[i+1],xyz[0],vec1);
    VCross3(vec0,vec1,cpvec);
    tarea[i] = VLen3(cpvec);

    areasum += tarea[i];
  }
    
  for (i = 1; i < n-1; i++) { 
    for (j = 0; j < 3; j++)
      pcen[j] += tarea[i]*tcen[i][j];
  }

  for (j = 0; j < 3; j++)
    pcen[j] /= areasum;
}

  
