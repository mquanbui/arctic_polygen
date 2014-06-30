#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "MSTK.h"
#include "mkdual.h"

int use_centroids_only=0;
int kill_small_edges=1;
double small_edge_tol=1.0e-06;

int main(int argc, char **argv) {
  int len, ok, opts[5]={0,0,0,0,0}, native=1, outopts[10]={0,0,0,0,0,0,0,0,0,0};
  int natt=-1, i;
  char **attnames=NULL, format[32], ext[8], ans[256];
  char infname[256], outfname[256], gmvfname[256], mname[256], rtype_str[8];
  Mesh_ptr mesh1, mesh2;
  

  switch(argc) {
  case 2:
    strcpy(infname,argv[1]);
    break;
  case 3: case 4: case 5:
    strcpy(infname,argv[1]);
    for (i = 2; i < argc; i++) {
      if (strncmp(argv[i],"use_centroids=",14) == 0) {
        sscanf(argv[i]+14,"%s",ans);
        if (strncmp(ans,"y",1) == 0)
          use_centroids_only = 1;
      }
      if (strncmp(argv[i],"kill_small_edges=",17) == 0) {
        sscanf(argv[i]+17,"%s",ans);
        if (strncmp(ans,"y",1) == 0)
          kill_small_edges = 1;
        else if (strncmp(ans,"n",1) == 0)
          kill_small_edges = 0;
      }
      if (strncmp(argv[i],"small_edge_tol=",15) == 0) {
        sscanf(argv[i]+15,"%lf",&small_edge_tol);
      }
    }
    break;
  default:
    fprintf(stderr,"Usage: %s infname.ext <use_centroid=y|n> <kill_small_edges=y|n> <small_edge_tol=...\n",argv[0]);
    exit(-1);
  }



  MSTK_Init();
    

  len = strlen(infname);
  if (len > 5 && (strncmp(&(infname[len-5]),".mstk",5) == 0)) {
    strcpy(mname,infname);
    mname[len-5] = '\0';
    strcpy(format,"mstk");
    strcpy(ext,"mstk");
    native = 1;
  }
  else if (len > 4 && (strncmp(&(infname[len-4]),".gmv",4) == 0)) {
    strcpy(mname,infname);
    mname[len-4] = '\0';
    strcpy(format,"gmv");
    strcpy(ext,"gmv");
    native = 0;
  }
  else if (len > 4 && (strncmp(&(infname[len-4]),".exo",4) == 0)) {
    strcpy(mname,infname);
    mname[len-4] = '\0';
    strcpy(format,"exodusii");
    strcpy(ext,"exo");
    fprintf(stderr,"Format %s, Ext %s\n",format,ext);
    native = 0;
  }
  else if (len > 4 && (strncmp(&(infname[len-4]),".",1) == 0)) {
    fprintf(stderr,"Unknown input file extension\n");
    exit(-1);
  }
  else {
    strcpy(mname,infname);
    strcpy(format,"mstk");
    strcpy(ext,"mstk");
    strcat(infname,".mstk");
    native = 1;
  }
    
  strcpy(outfname,mname);
  strcat(outfname,"-dual.");
  strcat(outfname,ext);
    

  if (native) {
    mesh1 = MESH_New(F1);
    ok = MESH_InitFromFile(mesh1,infname,NULL);
    if (!ok) {
      fprintf(stderr,"Cannot open input file %s\n",infname);
      exit(-1);
    }
    
  }
  else {

    mesh1 = MESH_New(F1);
    ok = MESH_ImportFromFile(mesh1,infname,format,NULL,NULL);
    if (!ok) {
      fprintf(stderr,"Cannot open input file %s\n",infname);
      exit(-1);
    }

  }


  MESH_BuildClassfn(mesh1);


  strcpy(rtype_str,MESH_RepType_Str(mesh1));
  if (rtype_str[0] != 'F') {
    fprintf(stderr,"\n WARNING: Input mesh is not full representation\n");
    fprintf(stderr,"Voronoi mesh can be represented only by full representations\n");
    fprintf(stderr,"Choosing representation type F1\n\n");
    
    mesh2 = MESH_New(F1);
  } 
  else
    mesh2 = MESH_New(MESH_RepType(mesh1));



  /************* MAIN CALL TO CREATE DUAL MESH ***************************/



  if (MESH_Num_Regions(mesh1)) 
    ok = MESH_MakeDual3(mesh1,mesh2,opts);
  else
    ok = MESH_MakeDual2(mesh1,mesh2,opts);


  /***********************************************************************/



  /**************** VERIFY THAT THE MESH IS OK ***************************/

  if (!MESH_CheckTopo(mesh2))
    fprintf(stderr,"Dual mesh is not valid\n");


  double maxdev;
  if (!MESH_HasPlanarFaces(mesh2,&maxdev)) {
    double PI = 3.1415926;
    fprintf(stderr,"***WARNING*** Mesh has curved faces\n");
    fprintf(stderr,"Maximum deviation of face normals is %-lf degrees\n", maxdev*180/PI);
  }
  else
    fprintf(stderr,"Mesh has only planar faces\n");

  /***********************************************************************/


  if (native) {
    MESH_WriteToFile(mesh2,outfname,F1,NULL);

    strcpy(gmvfname,mname);
    strcat(gmvfname,"-dual.gmv");
    MESH_ExportToGMV(mesh2,gmvfname,0,NULL,NULL,NULL);
  }
  else {

    MESH_ExportToFile(mesh2,outfname,format,natt,(const char **)attnames,outopts,NULL);

  }

  MESH_Delete(mesh1);
  MESH_Delete(mesh2);

  exit(0);
}


  
