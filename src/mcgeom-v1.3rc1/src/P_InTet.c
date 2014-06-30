#include <stdio.h>
#include <stdlib.h>
#include "MCGeom.h"

/* Check if point is in on the boundary of the tet. The assumption is
   that the face formed by the first three vertices points towards the
   fourth vertex. Tolerance and flag are not used presently */

int P_InTet(double *pxyz, double (*txyz)[3], double tol, int flag) {
  int i;
  double xyz[4][3], vol;
  

  VCopy3(xyz[3],pxyz);


  /* Tet 1,2,3,P or in C notation, Tet 0,1,2,P */

  for (i = 0; i < 3; i++)
    VCopy3(xyz[i],txyz[i]);
  
  vol = Tet_Volume(xyz);

  if (vol < 0)
    return 0;



  /* Tet 1,3,4,P or in C notation, Tet 0,2,3,P */

  VCopy3(xyz[1],txyz[2]);
  VCopy3(xyz[2],txyz[3]);

  vol = Tet_Volume(xyz);

  if (vol < 0)
    return 0;



  /* Tet 1,4,2,P or in C notation, Tet 0,3,1,P */

  VCopy3(xyz[1],txyz[3]);
  VCopy3(xyz[2],txyz[1]);

  vol = Tet_Volume(xyz);

  if (vol < 0)
    return 0;



  /* Tet 2,4,3,P or in C notation, Tet 1,3,2,P */

  VCopy3(xyz[0],txyz[1]);
  VCopy3(xyz[1],txyz[3]);
  VCopy3(xyz[2],txyz[2]);

  vol = Tet_Volume(xyz);

  if (vol < 0)
    return 0;


  return 1;

} 
