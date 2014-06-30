#include <stdio.h>
#include <math.h>
#include "MCGeom.h"


double Hex_Volume(double (*rxyz)[3]) {
  int i, j, k, ind;
  double fxyz[6][3], zxyz[3], txyz[4][3], vol;

  /* Templates for forming the tets - for a given tet, i, of the hex,
     the tet vertices are rxyz[tet_ev_index[i][0]],
     rxyz[tet_ev_index[i][1]], fxyz[tet_f_index[i]] and zxyz */

  static int tet_ev_index[24][2] = {{0,1},{1,0},{1,2},{2,1},
				    {2,3},{3,2},{3,0},{0,3},
				    {4,5},{5,4},{5,6},{6,5},
				    {6,7},{7,6},{7,4},{4,7},
				    {0,4},{4,0},{1,5},{5,1},
				    {2,6},{6,2},{3,7},{7,3}};
  static int tet_f_index[24] = {0,2,0,3,0,4,0,5,
				2,1,3,1,4,1,5,1,
				2,5,3,2,4,3,5,4};

  /* Compute the volume of a hex by decomposing it into 24 tets */

  /* Compute face centers */

  /* Face 0 - bottom face */

  for (k = 0; k < 3; k++)
    fxyz[0][k] = (rxyz[0][k]+rxyz[1][k]+rxyz[2][k]+rxyz[3][k])/4.0;

  /* Face 1 - top face */

  for (k = 0; k < 3; k++)
    fxyz[1][k] = (rxyz[4][k]+rxyz[5][k]+rxyz[6][k]+rxyz[7][k])/4.0;

  /* Lateral faces */

  for (k = 0; k < 3; k++)
    fxyz[2][k] = (rxyz[0][k]+rxyz[1][k]+rxyz[5][k]+rxyz[4][k])/4.0;
    
  for (k = 0; k < 3; k++)
    fxyz[3][k] = (rxyz[1][k]+rxyz[2][k]+rxyz[6][k]+rxyz[5][k])/4.0;
    
  for (k = 0; k < 3; k++)
    fxyz[4][k] = (rxyz[2][k]+rxyz[3][k]+rxyz[7][k]+rxyz[6][k])/4.0;
    
  for (k = 0; k < 3; k++)
    fxyz[5][k] = (rxyz[3][k]+rxyz[0][k]+rxyz[4][k]+rxyz[7][k])/4.0;
    

  /* hex center - is it better to average the face centers??? */

  for (k = 0; k < 3; k++)
    zxyz[k] = (rxyz[0][k]+rxyz[1][k]+rxyz[2][k]+rxyz[3][k]+
	       rxyz[4][k]+rxyz[5][k]+rxyz[6][k]+rxyz[7][k])/8.0;


  vol = 0.0;

  /* Compute volumes of tets forming the hex */

  /* vertex 4 of tets stays the same */

  for (k = 0; k < 3; k++)
    txyz[3][k] = zxyz[k];


  for (i = 0; i < 24; i++) {

    /* vertex 1 of tet */

    ind = tet_ev_index[i][0];
    for (k = 0; k < 3; k++) 
      txyz[0][k] = rxyz[ind][k];

    /* vertex 2 of tet */

    ind = tet_ev_index[i][1];
    for (k = 0; k < 3; k++)
      txyz[1][k] = rxyz[ind][k];

    /* vertex 3 of tet - face center of hex */

    ind = tet_f_index[i];
    for (k = 0; k < 3; k++)
      txyz[2][k] = fxyz[ind][k];

    /* Tet volume */

    vol += Tet_Volume(txyz);

  }

  return vol;
}
