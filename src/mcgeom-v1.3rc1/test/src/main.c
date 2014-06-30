#include <stdlib.h>
#include <stdio.h>
#include "MCGeom.h"

int main() {

  double vxyz[3], pxyz[25][3], eps=1.0e-12;
  int i, j, c, np, done, flag, status;


  done = 0;
  while (!done) {

    fprintf(stderr,"Number of vertices in polygon: \n");
    status = fscanf(stdin,"%d",&np);
    if (status != 1) {
      done = 1;
      continue;
    }

    fprintf(stderr,"Polygon vertex coordinates:\n");
    for (i = 0; i < np; i++)
      fscanf(stdin,"%lf %lf %lf",&(pxyz[i][0]),&(pxyz[i][1]),&(pxyz[i][2]));
      

    int done1 = 0, status;
    while (!done1) {
      fprintf(stderr,"\n");
      fprintf(stderr,"Point coordinates:\n");
      status = fscanf(stdin,"%lf %lf %lf",&(vxyz[0]),&(vxyz[1]),&(vxyz[2]));
      if (status != 3) { 
	done1 = 1;
	continue;
      }

      flag = P_InPolyFace2D(vxyz,np,pxyz,1.0e-08,1);

      /*
      c = 0;
      for (i = 0, j = np-1; i < np; j = i++) {
        if ((((pxyz[i][1]<=vxyz[1]) && (vxyz[1]<pxyz[j][1])) ||
             ((pxyz[j][1]<=vxyz[1]) && (vxyz[1]<pxyz[i][1]))) &&
            (vxyz[0] < (pxyz[j][0] - pxyz[i][0]) * (vxyz[1] - pxyz[i][1]) / (pxyz[j][1] - pxyz[i][1]) +
pxyz[i][0]))

          c = !c;
      }
 
      flag = c;
      */

      if (flag)
	fprintf(stderr,"In !!\n\n");
      else
	fprintf(stderr,"Out 8-(\n\n");
    }

  }
}
