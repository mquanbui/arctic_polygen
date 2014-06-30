#include <stdio.h>
#include <stdlib.h>
#include "MSTK.h"
#include "MCGeom.h"


extern int use_centroids_only;
extern int kill_small_edges;
extern double small_edge_tol;

/* Convert a triangular mesh to a polygonal mesh */

int MESH_MakeDual2(Mesh_ptr dmesh, Mesh_ptr vmesh, int *opts) {
  MFace_ptr f0_dm, fcur_dm, fnext_dm, f_dm, ef;
  MFace_ptr *ef_vm, *vf_vm, *rfarray, f_vm;
  MEdge_ptr e0_dm, ecur_dm, enext_dm, e_dm, e_vm;
  MVertex_ptr v_dm, *fv_vm, v_vm, vf_dm;
  MAttrib_ptr vmvatt, vmfatt, vmnfatt, dmentatt;
  List_ptr efaces_dm, fedges_dm, vedges_dm;
  List_ptr vfaces_dm;
  List_ptr fedges_vm, efaces_vm, vedges_vm;
  List_ptr vedges_sub;
  int i, j, k;
  int del, found, done, nfv_vm, nfv_dm, markid, fedir, ival;
  int idx1, idx2, idx3, idx4;
  int nef_vm, nvf_vm, nfe_dm, evidx, dir, outerdone, topocase;
  int nnewf, *rfdirs, nrf_vm, firstwarn1=1, firstwarn2=1;
  int maxnfv_vm=100; /* max number of vertices in a voronoi face */
  int maxnef_vm=10;  /* this is the max number of voronoi faces
		        associated with a Delaunay edge */
  int maxnvf_vm=10;  /* this is the max number of special boundary 
			faces associated with a boundary vertex */
  double exyz[2][3], fxyz[MAXPV2][3], cxyz[3], rval;
  double evec[3], vec1[3], vec2[3], normal[3], dp, cen[3], vxyz[3];
  void *pval;

  
  del = opts[0];  /* whether mesh is delaunay or not */


  fv_vm = (MVertex_ptr *) malloc(maxnfv_vm*sizeof(MVertex_ptr));

  vmvatt = MAttrib_New(dmesh,"VMVATT",POINTER,MALLTYPE);
  vmfatt = MAttrib_New(dmesh,"VMFATT",POINTER,MALLTYPE);
  vmnfatt = MAttrib_New(dmesh,"VMNFATT",INT,MALLTYPE);
  dmentatt = MAttrib_New(vmesh,"DMENTATT",POINTER,MFACE);


  markid = MSTK_GetMarker();


  idx1 = 0;
  while ((v_dm = MESH_Next_Vertex(dmesh,&idx1))) {

    vf_vm = (MFace_ptr *) malloc(maxnvf_vm*sizeof(MFace_ptr));

    nvf_vm = 0;    /* Number of dual faces for this Delaunay vertex */

    vedges_dm = MV_Edges(v_dm);
    List_Mark(vedges_dm,markid);
   
    MEnt_Set_AttVal(v_dm,vmfatt,0,0.0,vf_vm);
    MEnt_Set_AttVal(v_dm,vmnfatt,0,0.0,NULL);

    /* Do different things based on what type of geometric model
       entity the edge is classified on */

    switch (MV_GEntDim(v_dm)) {
    case 0: case 1:

      idx2 = 0;
      while ((e0_dm = List_Next_Entry(vedges_dm,&idx2))) {

	if (ME_GEntDim(e0_dm) == 2) continue;

	/* e0_dm is on a model edge. Walk through tris from e0_dm to
	   the next mesh edge on a model edge */

	ecur_dm = e0_dm;

	efaces_dm = ME_Faces(ecur_dm);
	idx3 = 0;
	while ((fcur_dm = List_Next_Entry(efaces_dm,&idx3))) {
	  if (ME_Vertex(ecur_dm,0) == v_dm) {
	    if (MF_EdgeDir(fcur_dm,ecur_dm))
	      break;
	  }
	  else {
	    if (!MF_EdgeDir(fcur_dm,ecur_dm))
	      break;
	  }	      
	}
	List_Delete(efaces_dm);


	/* 
	   If this is a boundary edge which does not have a face on
	   the side we are interested in walking through, then skip
	   it. We will walk TO this edge from another boundary edge
	   connected to the vertex.
	*/

	if (!fcur_dm) continue;

	nfv_vm = 0;

	if (MV_GEntDim(v_dm) == 0) {

	  /* First vertex of Voronoi face is the Delaunay vertex
	     itself */
	  
	  fv_vm[nfv_vm] = MV_New(vmesh);
	  MV_Coords(v_dm,vxyz);
	  MV_Set_Coords(fv_vm[nfv_vm],vxyz);
	  MV_Set_GEntDim(fv_vm[nfv_vm],MV_GEntDim(v_dm));
	  MV_Set_GEntID(fv_vm[nfv_vm],MV_GEntID(v_dm));

	  nfv_vm++;
	}


	/* First vertex (if v_dm is on model edge) or second vertex
	   (if v_dm is on model vertex) is midpoint of mesh edge */

	MEnt_Get_AttVal(ecur_dm,vmvatt,&ival,&rval,&pval);

	/* Is there a vertex already associated with the midpoint of
	   this edge? */

	fv_vm[nfv_vm] = pval;

	if (!fv_vm[nfv_vm]) {

	  MV_Coords(ME_Vertex(ecur_dm,0),exyz[0]);
	  MV_Coords(ME_Vertex(ecur_dm,1),exyz[1]);

	  cxyz[0] = (exyz[0][0]+exyz[1][0])/2.0;
	  cxyz[1] = (exyz[0][1]+exyz[1][1])/2.0;
	  cxyz[2] = (exyz[0][2]+exyz[1][2])/2.0;

	  fv_vm[nfv_vm] = MV_New(vmesh);
	  MV_Set_Coords(fv_vm[nfv_vm],cxyz);
	  MV_Set_GEntDim(fv_vm[nfv_vm],ME_GEntDim(ecur_dm));
	  MV_Set_GEntID(fv_vm[nfv_vm],ME_GEntID(ecur_dm));

	  MEnt_Set_AttVal(ecur_dm,vmvatt,0,0.0,fv_vm[nfv_vm]);
	}

	nfv_vm++;




	/* Now traverse the faces in ccw order around the vertex until
	   we come to the next edge of the vertex classified on a
	   model edge */

	done = 0;
	while (!done) { 
	  
	  
	  /* Is there already a Voronoi vertex associated with the
	     circumcenter of this tri? */

	  MEnt_Get_AttVal(fcur_dm,vmvatt,&ival,&rval,&pval);	  
	  
	  fv_vm[nfv_vm] = (MVertex_ptr) pval;	  
	  
	  if (!fv_vm[nfv_vm]) {
	    
	    MF_Coords(fcur_dm,&nfv_dm,fxyz); 

	    if (use_centroids_only) {
	      for (i = 0; i < 3; i++)
		cxyz[i] = (fxyz[0][i]+fxyz[1][i]+fxyz[2][i])/3.0;
	    }
	    else {
	      Tri_CircumCen(fxyz,cxyz);
	      
	      if (!P_InPolyFace3D(cxyz,3,fxyz,0.0,0)) {
		
		/* Circumcenter not in tri - ok if this tri has an edge on
		   the boundary but otherwise the mesh is not Delaunay */
		
		int inttri=1;  /* assume interior tri */
		MEdge_ptr edum;
		
		fedges_dm = MF_Edges(fcur_dm,1,0);
		for (i = 0; i < 3; i++) {
		  edum = List_Entry(fedges_dm,i);
		  if (ME_GEntDim(edum) == 1) {
		    inttri = 0;
		    break;
		  }
		}
		List_Delete(fedges_dm);
		
		if (!inttri || !del) {
		  if (firstwarn1) {
		    fprintf(stderr,"Circumcenter not inside tri\n");
		    fprintf(stderr,"Using geometric center of tri\n");
		    fprintf(stderr,"Expect non-Voronoi mesh\n");
		    firstwarn1 = 0;
		  }
		  
		  for (i = 0; i < 3; i++)
		    cxyz[i] = (fxyz[0][i]+fxyz[1][i]+fxyz[2][i])/3.0;
		}
	      }
	    }

	    /* Make vertex at circumcenter of tri or geometric center
	       of tri if the cirumcenter is outside the tri */	  
	    
	    fv_vm[nfv_vm] = MV_New(vmesh);        
	    MV_Set_Coords(fv_vm[nfv_vm],cxyz);
	    MV_Set_GEntDim(fv_vm[nfv_vm],2);
	    MV_Set_GEntID(fv_vm[nfv_vm],MR_GEntID(fcur_dm));
	    
	    MEnt_Set_AttVal(fcur_dm,vmvatt,0,0.0,fv_vm[nfv_vm]);

	  }
	  nfv_vm++;
	    
	  /* Tri edges in ccw order starting from v_dm */
	  fedges_dm = MF_Edges(fcur_dm,1,v_dm);	 
	  ecur_dm = List_Entry(fedges_dm,2);  /* Last edge of triangle */
	  List_Delete(fedges_dm);
	  
	  if (ME_GEntDim(ecur_dm) == 1) 
	    done = 1;
	  else {
	    /* Get tri "in front of the edge" */
	    
	    efaces_dm = ME_Faces(ecur_dm);
	    idx3 = 0;
	    while ((fcur_dm = List_Next_Entry(efaces_dm,&idx3))) {
	      if (ME_Vertex(ecur_dm,0) == v_dm) {
		if (MF_EdgeDir(fcur_dm,ecur_dm))
		  break;
	      }
	      else {
		if (!MF_EdgeDir(fcur_dm,ecur_dm))
		  break;
	      }
	    }
	    List_Delete(efaces_dm);
	    
	    if (!fcur_dm) {
	      fprintf(stderr,"Trouble!! Could not find next face\n");
	      exit(-1);
	    }
	  }

	} /* while (!done) */


	/* Last vertex is the midpoint of the last edge */

	MEnt_Get_AttVal(ecur_dm,vmvatt,&ival,&rval,&pval);
	
	fv_vm[nfv_vm] = pval;

	if (!fv_vm[nfv_vm]) {

	  MV_Coords(ME_Vertex(ecur_dm,0),exyz[0]);
	  MV_Coords(ME_Vertex(ecur_dm,1),exyz[1]);

	  cxyz[0] = (exyz[0][0]+exyz[1][0])/2.0;
	  cxyz[1] = (exyz[0][1]+exyz[1][1])/2.0;
	  cxyz[2] = (exyz[0][2]+exyz[1][2])/2.0;

	  fv_vm[nfv_vm] = MV_New(vmesh);
	  MV_Set_Coords(fv_vm[nfv_vm],cxyz);
	  MV_Set_GEntDim(fv_vm[nfv_vm],ME_GEntDim(ecur_dm));
	  MV_Set_GEntID(fv_vm[nfv_vm],ME_GEntID(ecur_dm));

	  MEnt_Set_AttVal(ecur_dm,vmvatt,0,0.0,fv_vm[nfv_vm]);
	}

	nfv_vm++;

		
	
	/* Create Voronoi face corresponding to Delaunay vertex */
	
	vf_vm[nvf_vm] = MF_New(vmesh);
	MF_Set_Vertices(vf_vm[nvf_vm], nfv_vm, fv_vm);
	
	MF_Set_GEntDim(vf_vm[nvf_vm], 2);
	MF_Set_GEntID(vf_vm[nvf_vm], MR_GEntID(fcur_dm));

	/* Attach primary mesh entity to dual mesh face - DO WE NEED THIS?? */
	MEnt_Set_AttVal(vf_vm[nvf_vm],dmentatt,0,0.0,v_dm); 
	
	nef_vm++;

      } /* while ((f0_dm = List_Next_Entry(efaces,&idx2))) */
      
      MEnt_Set_AttVal(v_dm,vmfatt,0,0.0,vf_vm);
      MEnt_Set_AttVal(v_dm,vmnfatt,nef_vm,0.0,NULL);

      break;


    case 2: 


      e0_dm = List_Entry(vedges_dm,0);
      ecur_dm = e0_dm;

      efaces_dm = ME_Faces(ecur_dm);
      idx3 = 0;
      while ((fcur_dm = List_Next_Entry(efaces_dm,&idx3))) {
	if (ME_Vertex(ecur_dm,0) == v_dm) {
	  if (MF_EdgeDir(fcur_dm,ecur_dm))
	    break;
	}
	else {
	  if (!MF_EdgeDir(fcur_dm,ecur_dm))
	    break;
	}	      
      }
      List_Delete(efaces_dm);



      /* Traverse the faces in ccw order around the vertex until we
	 come to the next edge of the vertex classified on a model
	 edge */

      nfv_vm = 0;
      done = 0;
      while (!done) { 
	  
	/* Is there already a Voronoi vertex associated with the
	   circumcenter of this tri? */

	MEnt_Get_AttVal(fcur_dm,vmvatt,&ival,&rval,&pval);	  
	  
	fv_vm[nfv_vm] = (MVertex_ptr) pval;	  
	
	if (!fv_vm[nfv_vm]) {
	  
	  MF_Coords(fcur_dm,&nfv_dm,fxyz); 
	  
	  if (use_centroids_only) {
	    for (i = 0; i < 3; i++)
	      cxyz[i] = (fxyz[0][i]+fxyz[1][i]+fxyz[2][i])/3.0;
	  }
	  else {
	    Tri_CircumCen(fxyz,cxyz);
	  
	    if (!P_InPolyFace3D(cxyz,3,fxyz,0.0,0)) {
	    
	      /* Circumcenter not in tri - ok if this tri has an edge on
		 the boundary but otherwise the mesh is not Delaunay */
	    
	      int inttri=1;  /* assume interior tri */
	      MEdge_ptr edum;
	    
	      fedges_dm = MF_Edges(fcur_dm,1,0);
	      for (i = 0; i < 3; i++) {
		edum = List_Entry(fedges_dm,i);
		if (ME_GEntDim(edum) == 1) {
		  inttri = 0;
		  break;
		}
	      }
	      List_Delete(fedges_dm);
	    
	      if (!inttri || !del) {
		if (firstwarn1) {
		  fprintf(stderr,"Circumcenter not inside tri\n");
		  fprintf(stderr,"Using geometric center of tri\n");
		  fprintf(stderr,"Expect non-planar faces\n");
		  firstwarn1 = 0;
		}
	      
		for (i = 0; i < 3; i++)
		  cxyz[i] = (fxyz[0][i]+fxyz[1][i]+fxyz[2][i])/3.0;
	      }
	    }
	  }

	  /* Make vertex at circumcenter of tet or geometric center
	     of tri if the cirumcenter is outside the tri */	  
	  
	  fv_vm[nfv_vm] = MV_New(vmesh);        
	  MV_Set_Coords(fv_vm[nfv_vm],cxyz);
	  MV_Set_GEntDim(fv_vm[nfv_vm],2);
	  MV_Set_GEntID(fv_vm[nfv_vm],MR_GEntID(fcur_dm));
	  
	  MEnt_Set_AttVal(fcur_dm,vmvatt,0,0.0,fv_vm[nfv_vm]);
	  
	}
	nfv_vm++;
	
	
	fedges_dm = MF_Edges(fcur_dm,1,v_dm);	  
	ecur_dm = List_Entry(fedges_dm,2);
	List_Delete(fedges_dm);
	
	if (ecur_dm == e0_dm) {
	  done = 1;
	}
	else {
	  /* Get tri "in front of the edge" */
	  
	  efaces_dm = ME_Faces(ecur_dm);
	  idx3 = 0;
	  while ((fcur_dm = List_Next_Entry(efaces_dm,&idx3))) {
	    if (ME_Vertex(ecur_dm,0) == v_dm) {
	      if (MF_EdgeDir(fcur_dm,ecur_dm))
		break;
	    }
	    else {
	      if (!MF_EdgeDir(fcur_dm,ecur_dm))
		break;
	    }
	  }
	  List_Delete(efaces_dm);
	  
	  if (!fcur_dm) {
	    fprintf(stderr,"Trouble!! Could not find next face\n");
	    exit(-1);
	  }
	}

      } /* while (!done) */


      /* Create Voronoi face corresponding to Delaunay vertex */
	
      vf_vm[nvf_vm] = MF_New(vmesh);
      MF_Set_Vertices(vf_vm[nvf_vm], nfv_vm, fv_vm);
      
      MF_Set_GEntDim(vf_vm[nvf_vm], 2);
      MF_Set_GEntID(vf_vm[nvf_vm], MR_GEntID(fcur_dm));
      
      /* Attach primary mesh entity to dual mesh face - DO WE NEED THIS?? */
      MEnt_Set_AttVal(vf_vm[nvf_vm],dmentatt,0,0.0,v_dm); 
	
      nvf_vm++;
      
      MEnt_Set_AttVal(v_dm,vmfatt,0,0.0,vf_vm);
      MEnt_Set_AttVal(v_dm,vmnfatt,nvf_vm,0.0,NULL);

      break;
    }

    List_Unmark(vedges_dm,markid);
    List_Delete(vedges_dm);
    
  }



  
  /*** Collapse out degenerate edges ***/

  /*** HAVE TO INFORM THE PRIMARY MESH ENTITIES THAT SOME ENTITIES 
       GOT DELETED ***/


/*   idx1 = 0; */
/*   while ((e_vm = MESH_Next_Edge(vmesh,&idx1))) { */
/*     MVertex_ptr ev[2]; */
/*     MEdge_ptr vedge; */
/*     double len2; */
/*     int gdim[2], gid[2]; */

/*     ev[0] = ME_Vertex(e_vm,0); */
/*     ev[1] = ME_Vertex(e_vm,1); */

/*     MV_Coords(ev[0],exyz[0]); */
/*     MV_Coords(ev[1],exyz[1]); */

/*     len2 = PP_Dist2(exyz[0],exyz[1]); */

/*     if (len2 > 0.0) continue; */

/*     /\* Degenerate edge which we have to collapse *\/ */

/*     gdim[0] = MV_GEntDim(ev[0]); */
/*     gid[0] = MV_GEntID(ev[0]); */
/*     gdim[1] = MV_GEntDim(ev[1]); */
/*     gid[1] = MV_GEntID(ev[1]); */

/*     if (gdim[0] == gdim[1] && gid[0] != gid[1]) { */
/*       fprintf(stderr,"Cannot collapse out zero length edge due to topological considerations\n"); */
/*       continue; */
/*     } */


/*     /\* for simplicity of coding we want to always collapse ev[0] to */
/*        ev[1] (i.e. replace ev[0] with ev[1]) *\/ */
/*     /\* for topological validity, we always want to collapse the vertex */
/*        classified on the higher dimension entity to the vertex */
/*        classified on the lower dimension entity *\/ */
/*     /\* So if gdim[1] > gdim[0], switch the two vertices around *\/ */

/*     if (gdim[1] > gdim[0]) { */
/*       MVertex_ptr vtmp; */
/*       vtmp = ev[0]; */
/*       ev[0] = ev[1]; */
/*       ev[1] = vtmp; */
/*     } */

/*     vedges_vm = MV_Edges(ev[0]); */
/*     idx2 = 0; */
/*     while ((vedge = List_Next_Entry(vedges_vm,&idx2))) { */
/*       ME_Replace_Vertex(vedge,ev[0],ev[1]); */
/*     } */
/*     List_Delete(vedges_vm); */


/*     efaces_vm = ME_Faces(e_vm); */

/*     idx2 = 0; */
/*     while ((ef = List_Next_Entry(efaces_vm,&idx2))) { */
/*       List_ptr fedges; */
/*       int nfe, found; */
/*       MEdge_ptr fedge, olde[3], nue[2]; */

/*       fedges = MF_Edges(ef,1,0); */
/*       nfe = List_Num_Entries(fedges); */

/*       if (nfe == 2) { */
/* 	MEntity_ptr ment; */

/* 	/\* Face has become degenerate - delete it. Also tell primary mesh */
/* 	 that this face has been deleted *\/ */

/* 	MEnt_Get_AttVal(ef,dmentatt,&ival,&rval,&ment); */

/* 	MEnt_Get_AttVal(ment,vmnfatt,&nef_vm,&rval,&pval); */
/* 	MEnt_Get_AttVal(ment,vmfatt,&ival,&rval,&pval); */
/* 	ef_vm = (MFace_ptr *) pval; */

/* 	found = 0; */
/* 	for (i = 0; i < nef_vm; i++) { */
/* 	  if (ef_vm[i] == ef) { */
/* 	    found = 1; */
/* 	    for (j = i; j < nef_vm-1; j++) */
/* 	      ef_vm[i] = ef_vm[i+1]; */
/* 	    nef_vm--; */
/* 	    MEnt_Set_AttVal(ment,vmnfatt,nef_vm,0.0,NULL); */
/* 	    break; */
/* 	  } */
/* 	} */
/* 	if (!found) { */
/* 	  fprintf(stderr,"Dual mesh face not attached to correct primary mesh entity\n"); */
/* 	} */

/* 	MF_Delete(ef,0); */
/* 	List_Delete(fedges); */
/* 	continue; */
/*       } */

/*       found = 0; */
/*       idx3 = 0; */
/*       for (i = 0; i < nfe; i++) { */
/* 	fedge = List_Entry(fedges,i); */
/* 	if (fedge == e_vm) { */

/* 	  found = 1; */
/* 	  olde[0] = List_Entry(fedges,(i-1+nfe)%nfe); */
/* 	  olde[1] = e_vm; */
/* 	  olde[2] = List_Entry(fedges,(i+1)%nfe); */
/* 	  nue[0] = olde[0]; */
/* 	  nue[1] = olde[2]; */

/* 	  MF_Replace_Edges(ef,3,olde,2,nue); */

/* 	  break; */
/* 	} */
/*       } */

/*       List_Delete(fedges); */

/*     } /\* while (ef ....) *\/ */

/*     ME_Delete(e_vm,0); */
/*     if (ev[0] != ev[1]) */
/*       MV_Delete(ev[0],0); */

/*     List_Delete(efaces_vm); */

/*   } /\* while (e_vm ....) *\/ */



  free(fv_vm);



  /* Collapse out degenerate edges (and any degenerate faces and
     regions that result) */

  if (kill_small_edges) {
    idx1 = 0;
    while ((e_vm = MESH_Next_Edge(vmesh,&idx1))) {
      MVertex_ptr ev[2];
      double len2;

      ev[0] = ME_Vertex(e_vm,0);
      ev[1] = ME_Vertex(e_vm,1);
      
      MV_Coords(ev[0],exyz[0]);
      MV_Coords(ev[1],exyz[1]);
      
      len2 = PP_Dist2(exyz[0],exyz[1]);
      
      if (len2 < small_edge_tol*small_edge_tol) {
	ME_Collapse(e_vm,ev[1],1);
      }

    } /* while (e_vm ....) */
  }
  
  return 1;
  
}
