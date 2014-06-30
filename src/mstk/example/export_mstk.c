#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "MSTK.h"


/* Example program to read a mesh and print out some information about it */


int main(int argc, char *argv[]) {
  int i, k, len, suff, idx, idx2, ok, edir, nv, ne, natt, ival, ncomp;
  double xyz[3], rval, *rval_arr;
  void *pval;
  char mname[256], filename[256], attname[256];
  Mesh_ptr mesh;
  MVertex_ptr v;
  MEdge_ptr e;
  MFace_ptr f;
  MAttrib_ptr attrib;
  GEntity_ptr gent;
  List_ptr fedges;
  MType attentdim;
  MAttType atttype;

  
  fprintf(stderr,"\n");

  if (argc == 1) {
    fprintf(stderr,"Usage: %s testfile.mstk",argv[0]);
    exit(-1);
  }

  /* Initialize MSTK - Always do this even if it does not seem to
     matter in this version of MSTK */

  MSTK_Init();

  /* Load the mesh */

  strcpy(mname,argv[1]);
  len = strlen(mname);
  suff = 0;
  if (len > 5) {
    k = len-5;
    while (!suff && k > 1) {
      if (strncmp(&(mname[k]),".mstk",5) == 0) {
        suff = 1;
      }
      else
        k--;
    }
  }
  if (suff)
    mname[len-5] = '\0';

  strcpy(filename,mname);
  strcat(filename,".gmv");

  mesh = MESH_New(UNKNOWN_REP);
  ok = MESH_ImportFromGMV(mesh,filename,NULL);
  if (!ok) {
    fprintf(stderr,"Cannot open input file %s\n\n\n",filename);
    exit(-1);
  }



  /* Write out a copy of the mesh in the F1 format */

  strcpy(filename,mname);
  strcat(filename,"-copy.mstk");
  MESH_WriteToFile(mesh,filename,F1,NULL);

  
  /* Deleting of mesh is not necessary if this the end of the program */

  MESH_Delete(mesh);

  return 1;
}
