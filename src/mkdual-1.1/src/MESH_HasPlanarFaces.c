#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "MSTK.h"
#include "MCGeom.h"

/* Check if all faces of a mesh are planar - return 1 if true */
/* Also return the maximum angle deviation (in radians) between
   normals at the corners of a face */

int MESH_HasPlanarFaces(Mesh_ptr mesh, double *maxdev) {
  int idx, i, j, k, nfv, status=0;
  double fxyz[MAXPV2][3], txyz[3][3], normal0[3], normal[3], dp, mindp;
  MFace_ptr mf;

  mindp = 0.0;
  idx = 0;
  while ((mf = MESH_Next_Face(mesh,&idx))) {

    if (MF_Num_Edges(mf) == 3) continue;

    MF_Coords(mf,&nfv,fxyz);

    for (i = 0; i < nfv; i++) {
      for (j = 0; j < 3; j++)
        for (k = 0; k < 3; k++)
          txyz[j][k] = fxyz[((i+nfv-1)+j)%nfv][k];
      
      if (i == 0)
        Tri_UnitNormal(txyz,normal0);
      else {
        Tri_UnitNormal(txyz,normal);
        
        dp = VDot3(normal0,normal);

        if (dp < mindp)
          mindp = dp;
      }
    }

  }

  status = (fabs(1.0-mindp) < 1.0e-6);
  *maxdev = acos(mindp);
  return status;
  
}
