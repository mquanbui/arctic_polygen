#include <stdio.h>
#include <malloc.h>

#include "MCGeom.h"

#ifdef __cplusplus
extern "C" {
#endif

  /* Find intersection polygon of two CONVEX polygons */

  /* Doesn't do anything fancy - just does half-space intersections of
     one polygon with the lines of the other */

  int X_FaceFace(int np1, double (*pxyz1)[3], int np2, double (*pxyz2)[3],
		 double ptol, int *nip, double (*ipxyz)[3], int *iptag) {

    double xmin1, xmin2, xmax1, xmax2, ymin1, ymin2, ymax1, ymax2;
    double xmin, xmax, ymin, ymax;
    double plen, pcen[3], plxyz[2][3], vdir[3];
    double pxyzl[1000][3], pxyzr[1000][3];
    int i, j, k, npl, npr, ptagl[1000], ptagr[1000], ie[1000];

    
    for (i = 0; i < 1000; i++) ie[i] = -1;

    
    xmin1 = ymin1 =  1.0e+20; xmax1 = ymax1 = -1.0e+20;
    for (i = 0; i < np1; i++) {
      if (xmin1 > pxyz1[i][0])	xmin1 = pxyz1[i][0];
      if (ymin1 > pxyz1[i][1])	ymin1 = pxyz1[i][1];
      if (xmax1 < pxyz1[i][0])	xmax1 = pxyz1[i][0];
      if (ymax1 < pxyz1[i][1])	ymax1 = pxyz1[i][1];
    }
    
    xmin2 = ymin2 =  1.0e+20; xmax2 = ymax2 = -1.0e+20;
    for (i = 0; i < np2; i++) {
      if (xmin2 > pxyz2[i][0])	xmin2 = pxyz2[i][0];
      if (ymin2 > pxyz2[i][1])	ymin2 = pxyz2[i][1];
      if (xmax2 < pxyz2[i][0])	xmax2 = pxyz2[i][0];
      if (ymax2 < pxyz2[i][1])	ymax2 = pxyz2[i][1];
    }
    
    if (xmin1 > xmax2 || xmin2 > xmax1 || ymin1 > ymax2 || ymin2 > ymax1) {
      *nip = 0;
      return 0;
    }
    
    xmin = (xmin1 < xmin2) ? xmin1 : xmin2;
    xmax = (xmax1 > xmax2) ? xmax1 : xmax2;
    ymin = (ymin1 < ymin2) ? ymin1 : ymin2;
    ymax = (ymax1 > ymax2) ? ymax1 : ymax2;
    
    plen = 1.5*sqrt((xmax-xmin)*(xmax-xmin) + (ymax-ymin)*(ymax-ymin));

    
    *nip = np2;
    for (i = 0; i < np2; i++) {
      VCopy3(ipxyz[i],pxyz2[i]);
      iptag[i] = 0;
    }

    for (i = 0; i < np1; i++) {

      for (k = 0; k < 3; k++)
	pcen[k] = 0.5*(pxyz1[i][k]+pxyz1[(i+1)%np1][k]);

      VDiff3(pxyz1[(i+1)%np1],pxyz1[i],vdir);
      VNormalize3(vdir);

      for (k = 0; k < 3; k++) {
	plxyz[0][k] = pcen[k] - plen*vdir[k];
	plxyz[1][k] = pcen[k] + plen*vdir[k];
      }


      PF_Slice(*nip, ipxyz, plxyz, ie, ptol, &npl, pxyzl, ptagl, &npr, pxyzr, 
	       ptagr);
      
      if (!npl) {
	*nip = 0;
	return 0;
      }

      *nip = npl;
      for (j = 0; j < npl; j++) {
	VCopy3(ipxyz[j],pxyzl[j]);
	iptag[j] = iptag[j] | ptagl[j];
      }

    }
    

    return 1;
  }


#ifdef __cplusplus
}
#endif
