#include <stdio.h>
#include <stdlib.h>
#include "MSTK.h"
#include "MCGeom.h"

/* Convert a Tetrahedral mesh to a Polyhedral mesh */


extern int use_centroids_only;
extern int kill_small_edges;
extern double small_edge_tol;

int MESH_MakeDual3(Mesh_ptr dmesh, Mesh_ptr vmesh, int *opts) {
  MRegion_ptr rcur_dm, r_vm;
  MFace_ptr f0_dm, fcur_dm, fnext_dm, f_dm, ef;
  MFace_ptr *ef_vm, *vf_vm, *rfarray, f_vm;
  MEdge_ptr e0_dm, ecur_dm, enext_dm, e_dm, e_vm;
  MVertex_ptr v_dm, *fv_vm, v_vm, vf_dm;
  MAttrib_ptr vmvatt, vmfatt, vmnfatt, dmentatt;
  List_ptr rfaces_dm, efaces_dm, fedges_dm, vedges_dm;
  List_ptr redges_dm, vfaces_dm, vregions_dm;
  List_ptr fedges_vm, efaces_vm, vedges_vm;
  List_ptr vedges_sub, vregions_sub;
  int i, j, ir, it, converged;
  int del, found, done, nfv_vm, nfv_dm, nrv_dm, markid, fedir, ival;
  int idx1, idx2, idx3, idx4;
  int nef_vm, nvf_vm, nfe_dm, evidx, dir, outerdone, topocase;
  int nnewf, *rfdirs, nrf_vm, firstwarn1=1, firstwarn2=1;
  int ngreg, gregions[10], greg0;
  int maxnfv_vm=100; /* max number of vertices in a voronoi face */
  int maxnef_vm=10;  /* this is the max number of voronoi faces
		        associated with a Delaunay edge */
  int maxnvf_vm=10;  /* this is the max number of special boundary 
			faces associated with a boundary vertex */
  int maxnrf_vm=MAXPF3;
  double exyz[2][3], fxyz[MAXPV2][3], rxyz[4][3], cxyz[3], rval;
  double evec[3], vec1[3], vec2[3], normal[3], avenormal[3], dp, cen[3];
  double cxyz1[3], cxyz2[3], dist2;
  void *pval;

  
  del = opts[0];  /* whether mesh is delaunay or not */


  fv_vm = (MVertex_ptr *) malloc(maxnfv_vm*sizeof(MVertex_ptr));


  vmvatt = MAttrib_New(dmesh,"VMVATT",POINTER,MALLTYPE);
  vmfatt = MAttrib_New(dmesh,"VMFATT",POINTER,MALLTYPE);
  vmnfatt = MAttrib_New(dmesh,"VMNFATT",INT,MALLTYPE);
  dmentatt = MAttrib_New(vmesh,"DMENTATT",POINTER,MFACE);


  markid = MSTK_GetMarker();


  /* Number of Voronoi faces spawned by Delaunay edge */

  idx1 = 0;
  while ((e_dm = MESH_Next_Edge(dmesh,&idx1))) {

    ef_vm = (MFace_ptr *) malloc(maxnef_vm*sizeof(MFace_ptr));

    nef_vm = 0;    /* Number of dual faces for this Delaunay edge */

    efaces_dm = ME_Faces(e_dm);
    List_Mark(efaces_dm,markid);
   
    MEnt_Set_AttVal(e_dm,vmfatt,0,0.0,ef_vm);
    MEnt_Set_AttVal(e_dm,vmnfatt,0,0.0,NULL);

    /* Do different things based on what type of geometric model
       entity the edge is classified on */

    switch (ME_GEntDim(e_dm)) {
    case 1:

      idx2 = 0;
      while ((f0_dm = List_Next_Entry(efaces_dm,&idx2))) {

	if (MF_GEntDim(f0_dm) == 3) continue; /* internal face */

	/* f0_dm is on a model face. Walk through tets from f0_dm to
	   the next mesh face on a model face */

	fcur_dm = f0_dm;
	fnext_dm = 0;


	fedir = MF_EdgeDir(fcur_dm,e_dm);
	rcur_dm = MF_Region(fcur_dm,fedir);

	/* If this is a boundary face which does not have a region on
	   the side we are interested in walking through, then skip it
	   For example, in a 2-manifold model, a mesh edge on a model
	   face will have two mesh faces classified on the model
	   face. One of these will have a region on the appropriate
	   side while the other will not. Then we walk through the
	   mesh regions starting from the former to the latter. */

	if (!rcur_dm) continue;


	nfv_vm = 0;


	/* First vertex is midpoint of mesh edge */

	MEnt_Get_AttVal(e_dm,vmvatt,&ival,&rval,&pval);

	/* Is there a vertex already associated with the midpoint of
	   this edge? */

	fv_vm[nfv_vm] = pval;

	if (!fv_vm[nfv_vm]) {

	  MV_Coords(ME_Vertex(e_dm,0),exyz[0]);
	  MV_Coords(ME_Vertex(e_dm,1),exyz[1]);

	  cxyz[0] = (exyz[0][0]+exyz[1][0])/2.0;
	  cxyz[1] = (exyz[0][1]+exyz[1][1])/2.0;
	  cxyz[2] = (exyz[0][2]+exyz[1][2])/2.0;

	  fv_vm[nfv_vm] = MV_New(vmesh);
	  MV_Set_Coords(fv_vm[nfv_vm],cxyz);
	  MV_Set_GEntDim(fv_vm[nfv_vm],1);
	  MV_Set_GEntID(fv_vm[nfv_vm],ME_GEntID(e_dm));

	  MEnt_Set_AttVal(e_dm,vmvatt,0,0.0,fv_vm[nfv_vm]);
	}
	
	nfv_vm++;


	/* Second vertex is circumcenter of starting face */

	MEnt_Get_AttVal(fcur_dm,vmvatt,&ival,&rval,&pval);
	
	/* Is there a vertex already associated with the circumcenter
	   of this triangle? */

	fv_vm[nfv_vm] = pval;

	if (!fv_vm[nfv_vm]) {

	  MF_Coords(fcur_dm,&nfv_dm,fxyz);

	  if (use_centroids_only) {
	    for (i = 0; i < 3; i++)
	      cxyz[i] = (fxyz[0][i]+fxyz[1][i]+fxyz[2][i])/3.0;
	  }
	  else {
	    Tri_CircumCen(fxyz,cxyz);

	    if (!P_InPolyFace3D(cxyz,3,fxyz,0.0,0)) {
	    
	      /* Just use the geometric center */

	      if (firstwarn2) {
		fprintf(stderr,"Circumcenter not in triangle!\n");
		fprintf(stderr,"Using geometric center instead\n");
		fprintf(stderr,"Expect non-planar faces\n");
		firstwarn2=0;
	      }

	      for (i = 0; i < 3; i++)
		cxyz1[i] = (fxyz[0][i]+fxyz[1][i]+fxyz[2][i])/3.0;

	      /* Try to find a point on the line joining the circumcenter
		 and the centroid that is closest to the circumcenter but
		 still inside the tetrahedron */

	      /* converged = 0; it = 0; */
	      converged = 1; it = 0; cxyz[0]=cxyz1[0];cxyz[1]=cxyz1[1];cxyz[2]=cxyz1[2];
	      while (!converged) {
		for (i = 0; i < 3; i++)
		  cxyz2[i] = (cxyz[i]+cxyz1[i])/2.0;
		if (P_InPolyFace3D(cxyz2,3,fxyz,0.0,0)) {
		  for (i = 0; i < 3; i++) cxyz1[i] = cxyz2[i];
		}
		else {
		  for (i = 0; i < 3; i++) cxyz[i] = cxyz2[i];
		}
		dist2 = (cxyz[0]-cxyz1[0])*(cxyz[0]-cxyz1[0]) +
		  (cxyz[1]-cxyz1[1])*(cxyz[1]-cxyz1[1]) +
		  (cxyz[2]-cxyz1[2])*(cxyz[2]-cxyz1[2]);

		if (dist2 < 1.0e-16) converged = 1;
		it++;
		if (it > 10) break;
	      }
	  }
	  }

	  /* Make vertex at circumcenter of triangle or at the
	     geometric center of the triangle if the circumcenter is
	     outside the triangle */

	  fv_vm[nfv_vm] = MV_New(vmesh);
	  MV_Set_Coords(fv_vm[nfv_vm],cxyz);
	  MV_Set_GEntDim(fv_vm[nfv_vm],2);
	  MV_Set_GEntID(fv_vm[nfv_vm],MF_GEntID(fcur_dm));

	  MEnt_Set_AttVal(fcur_dm,vmvatt,0,0.0,fv_vm[nfv_vm]);
	}

	nfv_vm++;


	done = 0;
	while (!done) { 
	  
	  /* Get tet "in front of the face" when the face is oriented to 
	     use the edge in a positive sense */
	  
	  fedir = MF_EdgeDir(fcur_dm,e_dm);
	  rcur_dm = MF_Region(fcur_dm,fedir);
	  
	  MEnt_Get_AttVal(rcur_dm,vmvatt,&ival,&rval,&pval);
	  
	  /* Is there already a Voronoi vertex associated with the
	     circumcenter of this tet? */
	  
	  fv_vm[nfv_vm] = (MVertex_ptr) pval;
	  
	  
	  if (!fv_vm[nfv_vm]) {
	    
	    MR_Coords(rcur_dm,&nrv_dm,rxyz);       /* nrv_dm = 4 */
	    
	    if (use_centroids_only) {
	      for (i = 0; i < 3; i++)
		cxyz[i] = (rxyz[0][i]+rxyz[1][i]+rxyz[2][i]+rxyz[3][i])/4.0;
	    }
	    else {
	      Tet_CircumCen(rxyz,cxyz);
	    
	      if (!P_InTet(cxyz,rxyz,0.0,0)) {
	      
		/* Circumcenter not in tet - ok if this tet has a face on
		   the boundary but otherwise the mesh is not Delaunay */
	      
		int inttet=1;  /* assume interior tet */
		MFace_ptr fdum;
	      
		rfaces_dm = MR_Faces(rcur_dm);
		for (i = 0; i < 4; i++) {
		  fdum = List_Entry(rfaces_dm,i);
		  if (MF_GEntDim(fdum) == 2) {
		    inttet = 0;
		    break;
		  }
		}
		List_Delete(rfaces_dm);
	      
		if (!inttet || !del) {
		  if (firstwarn1) {
		    fprintf(stderr,"Circumcenter not inside tet\n");
		    fprintf(stderr,"Using geometric center of tet\n");
		    fprintf(stderr,"Expect non-planar faces\n");
		    firstwarn1 = 0;
		  }

		  for (i = 0; i < 3; i++)
		    cxyz1[i] = (rxyz[0][i]+rxyz[1][i]+rxyz[2][i]+rxyz[3][i])/4.0;

		  /* converged = 0; it = 0; */
		  converged = 1; it = 0; cxyz[0]=cxyz1[0];cxyz[1]=cxyz1[1];cxyz[2]=cxyz1[2];
		  while (!converged) {
		    for (i = 0; i < 3; i++)
		      cxyz2[i] = (cxyz[i]+cxyz1[i])/2.0;
		    if (P_InTet(cxyz2,rxyz,0.0,0)) {
		      for (i = 0; i < 3; i++) cxyz1[i] = cxyz2[i];
		    }
		    else {
		      for (i = 0; i < 3; i++) cxyz[i] = cxyz2[i];
		    }
		    dist2 = (cxyz[0]-cxyz1[0])*(cxyz[0]-cxyz1[0]) +
		      (cxyz[1]-cxyz1[1])*(cxyz[1]-cxyz1[1]) +
		      (cxyz[2]-cxyz1[2])*(cxyz[2]-cxyz1[2]);
		  
		    if (dist2 < 1.0e-16) converged = 1;
		    it++;
		    if (it > 10) break;
		  }
		}
	      }
	    }

	    /* Make vertex at circumcenter of tet or geometric center
	       of tet if the cirumcenter is outside the tet */	  
	    
	    fv_vm[nfv_vm] = MV_New(vmesh);        
	    MV_Set_Coords(fv_vm[nfv_vm],cxyz);
	    MV_Set_GEntDim(fv_vm[nfv_vm],3);
	    MV_Set_GEntID(fv_vm[nfv_vm],MR_GEntID(rcur_dm));
	    
	    MEnt_Set_AttVal(rcur_dm,vmvatt,0,0.0,fv_vm[nfv_vm]);

	  }
	  nfv_vm++;
	    

	  
	  rfaces_dm = MR_Faces(rcur_dm);
	  
	  found = 0; i = 0;
	  while (i < 4 && !found) {
	    f_dm = List_Entry(rfaces_dm,i);
	    if (f_dm != fcur_dm && MEnt_IsMarked(f_dm,markid)) { 
	      fnext_dm = f_dm;
	      found = 1;
	    }
	    i++;
	  }
	  
	  List_Delete(rfaces_dm);
	  
	  if (MF_GEntDim(fnext_dm) == 2) done = 1;
	  fcur_dm = fnext_dm;
	}


	/* Last vertex is vertex associated with circumcenter of
	   terminating face */
	
	MEnt_Get_AttVal(fcur_dm,vmvatt,&ival,&rval,&pval);
	
	/* Is there a vertex already associated with the circumcenter
	   of this triangle? */

	fv_vm[nfv_vm] = pval;

	if (!fv_vm[nfv_vm]) {

	  MF_Coords(fcur_dm,&nfv_dm,fxyz);

	  if (use_centroids_only) {
	    for (i = 0; i < 3; i++)
	      cxyz[i] = (fxyz[0][i]+fxyz[1][i]+fxyz[2][i])/3.0;	    		  }
	  else {
	    Tri_CircumCen(fxyz,cxyz);

	    if (!P_InPolyFace3D(cxyz,3,fxyz,0.0,0)) {
	    
	      /* Just use the geometric center */

	      if (firstwarn2) {
		fprintf(stderr,"Circumcenter not in triangle!\n");
		fprintf(stderr,"Using geometric center instead\n");
		fprintf(stderr,"Expect non-planar faces\n");
		firstwarn2=0;
	      }

	      for (i = 0; i < 3; i++)
		cxyz1[i] = (fxyz[0][i]+fxyz[1][i]+fxyz[2][i])/3.0;	    	    

	      /* converged = 0; it = 0; */
	      converged = 1; it = 0; cxyz[0]=cxyz1[0];cxyz[1]=cxyz1[1];cxyz[2]=cxyz1[2];
	      while (!converged) {
		for (i = 0; i < 3; i++)
		  cxyz2[i] = (cxyz[i]+cxyz1[i])/2.0;
		if (P_InPolyFace3D(cxyz2,3,fxyz,0.0,0)) {
		  for (i = 0; i < 3; i++) cxyz1[i] = cxyz2[i];
		}
		else {
		  for (i = 0; i < 3; i++) cxyz[i] = cxyz2[i];
		}
		dist2 = (cxyz[0]-cxyz1[0])*(cxyz[0]-cxyz1[0]) +
		  (cxyz[1]-cxyz1[1])*(cxyz[1]-cxyz1[1]) +
		  (cxyz[2]-cxyz1[2])*(cxyz[2]-cxyz1[2]);

		if (dist2 < 1.0e-16) converged = 1;
		it++;
		if (it > 10) break;
	      }

	    }
	  }

	  /* Make vertex at circumcenter of triangle or at the
	     geometric center of the triangle if the circumcenter is
	     outside the triangle */

	  fv_vm[nfv_vm] = MV_New(vmesh);
	  MV_Set_Coords(fv_vm[nfv_vm],cxyz);
	  MV_Set_GEntDim(fv_vm[nfv_vm],2);
	  MV_Set_GEntID(fv_vm[nfv_vm],MF_GEntID(fcur_dm));

	  MEnt_Set_AttVal(fcur_dm,vmvatt,0,0.0,fv_vm[nfv_vm]);
	}

	nfv_vm++;
	
	
	
	/* Add modified Voronoi face corresponding to Delaunay edge */
	
	ef_vm[nef_vm] = MF_New(vmesh);
	MF_Set_Vertices(ef_vm[nef_vm], nfv_vm, fv_vm);
	
	MF_Set_GEntDim(ef_vm[nef_vm], MREGION);
	MF_Set_GEntID(ef_vm[nef_vm], MR_GEntID(rcur_dm));

	/* Attach primary mesh entity to dual mesh face */
	MEnt_Set_AttVal(ef_vm[nef_vm],dmentatt,0,0.0,e_dm);
	
	nef_vm++;

      } /* while ((f0_dm = List_Next_Entry(efaces,&idx2))) */
      
      /*      MEnt_Set_AttVal(e_dm,vmfatt,0,0.0,ef_vm);	*/
      MEnt_Set_AttVal(e_dm,vmnfatt,nef_vm,0.0,NULL);

      break;



    case 2: 

      /* Just like for mesh edges on a model edge, except that we
	 don't have to introduce the midpoint of the edge as a
	 vertex */

      idx2 = 0;
      while ((f0_dm = List_Next_Entry(efaces_dm,&idx2))) {

	if (MF_GEntDim(f0_dm) == 3) continue; /* internal face */

	/* f0_dm is on a model face. Walk through tets from f0_dm to
	   the next mesh face on a model face */

	fcur_dm = f0_dm;
	fnext_dm = 0;


	fedir = MF_EdgeDir(fcur_dm,e_dm);
	rcur_dm = MF_Region(fcur_dm,fedir);

	/* If this is a boundary face which does not have a region on
	   the side we are interested in walking through, then skip it
	   For example, in a 2-manifold model, a mesh edge on a model
	   face will have two mesh faces classified on the model
	   face. One of these will have a region on the appropriate
	   side while the other will not. Then we walk through the
	   mesh regions starting from the former to the latter. */

	if (!rcur_dm) continue;

	
	nfv_vm = 0;

	/* First vertex is circumcenter of starting face */

	MEnt_Get_AttVal(fcur_dm,vmvatt,&ival,&rval,&pval);
	
	/* Is there a vertex already associated with the circumcenter
	   of this triangle? */

	fv_vm[nfv_vm] = pval;

	if (!fv_vm[nfv_vm]) {

	  MF_Coords(fcur_dm,&nfv_dm,fxyz);

	  if (use_centroids_only) {
	    for (i = 0; i < 3; i++)
	      cxyz1[i] = (fxyz[0][i]+fxyz[1][i]+fxyz[2][i])/3.0;	    		  }
	  else {
	    Tri_CircumCen(fxyz,cxyz);

	    if (!P_InPolyFace3D(cxyz,3,fxyz,0.0,0)) {
	    
	      /* Just use the geometric center */

	      if (firstwarn2) {
		fprintf(stderr,"Circumcenter not in triangle!\n");
		fprintf(stderr,"Using geometric center instead\n");
		fprintf(stderr,"Expect non-planar faces\n");
		firstwarn2=0;
	      }

	      for (i = 0; i < 3; i++)
		cxyz1[i] = (fxyz[0][i]+fxyz[1][i]+fxyz[2][i])/3.0;	    	    

	      /*	    converged = 0; it = 0; */
	      converged = 1; it = 0; cxyz[0]=cxyz1[0];cxyz[1]=cxyz1[1];cxyz[2]=cxyz1[2];
	      while (!converged) {
		for (i = 0; i < 3; i++)
		  cxyz2[i] = (cxyz[i]+cxyz1[i])/2.0;
		if (P_InPolyFace3D(cxyz2,3,fxyz,0.0,0)) {
		  for (i = 0; i < 3; i++) cxyz1[i] = cxyz2[i];
		}
		else {
		  for (i = 0; i < 3; i++) cxyz[i] = cxyz2[i];
		}
		dist2 = (cxyz[0]-cxyz1[0])*(cxyz[0]-cxyz1[0]) +
		  (cxyz[1]-cxyz1[1])*(cxyz[1]-cxyz1[1]) +
		  (cxyz[2]-cxyz1[2])*(cxyz[2]-cxyz1[2]);

		if (dist2 < 1.0e-16) converged = 1;
		it++;
		if (it > 10) break;
	      }

	    }
	  }

	  /* Make vertex at circumcenter of triangle or at the
	     geometric center of the triangle if the circumcenter is
	     outside the triangle */

	  fv_vm[nfv_vm] = MV_New(vmesh);
	  MV_Set_Coords(fv_vm[nfv_vm],cxyz);
	  MV_Set_GEntDim(fv_vm[nfv_vm],2);
	  MV_Set_GEntID(fv_vm[nfv_vm],MF_GEntID(fcur_dm));

	  MEnt_Set_AttVal(fcur_dm,vmvatt,0,0.0,fv_vm[nfv_vm]);
	}

	nfv_vm++;


	done = 0;
	while (!done) { 
	  
	  /* Get tet "in front of the face" when the face is oriented to 
	     use the edge in a positive sense */
	  
	  fedir = MF_EdgeDir(fcur_dm,e_dm);
	  rcur_dm = MF_Region(fcur_dm,fedir);
	  
	  MEnt_Get_AttVal(rcur_dm,vmvatt,&ival,&rval,&pval);
	  
	  /* Is there already a Voronoi vertex associated with the
	     circumcenter of this tet? */
	  
	  fv_vm[nfv_vm] = (MVertex_ptr) pval;
	  
	  
	  if (!fv_vm[nfv_vm]) {
	    
	    MR_Coords(rcur_dm,&nrv_dm,rxyz);       /* nrv_dm = 4 */
	    
	    if (use_centroids_only) {
	      for (i = 0; i < 3; i++)
		cxyz[i] = (rxyz[0][i]+rxyz[1][i]+rxyz[2][i]+rxyz[3][i])/4.0;
	    }
	    else {
	      Tet_CircumCen(rxyz,cxyz);
	    
	      if (!P_InTet(cxyz,rxyz,0.0,0)) {
	      
		/* Circumcenter not in tet - ok if this tet has a face on
		   the boundary but otherwise the mesh is not Delaunay */
	      
		int inttet=1;  /* assume interior tet */
		MFace_ptr fdum;
	      
		rfaces_dm = MR_Faces(rcur_dm);
		for (i = 0; i < 4; i++) {
		  fdum = List_Entry(rfaces_dm,i);
		  if (MF_GEntDim(fdum) == 2) {
		    inttet = 0;
		    break;
		  }
		}
		List_Delete(rfaces_dm);
	      
		if (!inttet || !del) {
		  if (firstwarn1) {
		    fprintf(stderr,"Circumcenter not inside tet\n");
		    fprintf(stderr,"Using geometric center of tet\n");
		    fprintf(stderr,"Expect non-planar faces\n");
		    firstwarn1=0;
		  }
	      
		  for (i = 0; i < 3; i++)
		    cxyz1[i] = (rxyz[0][i]+rxyz[1][i]+rxyz[2][i]+rxyz[3][i])/4.0;

		  /* converged = 0; it = 0; */
		  converged = 1; it = 0; cxyz[0]=cxyz1[0];cxyz[1]=cxyz1[1];cxyz[2]=cxyz1[2];
		  while (!converged) {
		    for (i = 0; i < 3; i++)
		      cxyz2[i] = (cxyz[i]+cxyz1[i])/2.0;
		    if (P_InTet(cxyz2,rxyz,0.0,0)) {
		      for (i = 0; i < 3; i++) cxyz1[i] = cxyz2[i];
		    }
		    else {
		      for (i = 0; i < 3; i++) cxyz[i] = cxyz2[i];
		    }
		    dist2 = (cxyz[0]-cxyz1[0])*(cxyz[0]-cxyz1[0]) +
		      (cxyz[1]-cxyz1[1])*(cxyz[1]-cxyz1[1]) +
		      (cxyz[2]-cxyz1[2])*(cxyz[2]-cxyz1[2]);
		  
		    if (dist2 < 1.0e-16) converged = 1;
		    it++;
		    if (it > 10) break;
		  }
		}
	      }
	    }


	    /* Make vertex at circumcenter of tet or geometric center
	       of tet if the cirumcenter is outside the tet */	  
	    
	    fv_vm[nfv_vm] = MV_New(vmesh);        
	    MV_Set_Coords(fv_vm[nfv_vm],cxyz);
	    MV_Set_GEntDim(fv_vm[nfv_vm],3);
	    MV_Set_GEntID(fv_vm[nfv_vm],MR_GEntID(rcur_dm));
	    
	    MEnt_Set_AttVal(rcur_dm,vmvatt,0,0.0,fv_vm[nfv_vm]);
	    
	  }
	  nfv_vm++;
	  

	  rfaces_dm = MR_Faces(rcur_dm);
	  
	  found = 0; i = 0;
	  while (i < 4 && !found) {
	    f_dm = List_Entry(rfaces_dm,i);
	    if (f_dm != fcur_dm && MEnt_IsMarked(f_dm,markid)) { 
	      fnext_dm = f_dm;
	      found = 1;
	    }
	    i++;
	  }
	  
	  List_Delete(rfaces_dm);
	  
	  if (MF_GEntDim(fnext_dm) == 2) done = 1;
	  fcur_dm = fnext_dm;
	}


	/* Last vertex is vertex associated with circumcenter of
	   terminating face */
	
	MEnt_Get_AttVal(fcur_dm,vmvatt,&ival,&rval,&pval);
	
	/* Is there a vertex already associated with the circumcenter
	   of this triangle? */

	fv_vm[nfv_vm] = pval;

	if (!fv_vm[nfv_vm]) {

	  MF_Coords(fcur_dm,&nfv_dm,fxyz);

	  if (use_centroids_only) {
	    for (i = 0; i < 3; i++)
	      cxyz[i] = (rxyz[0][i]+rxyz[1][i]+rxyz[2][i]+rxyz[3][i])/4.0;
	  }
	  else {
	    Tri_CircumCen(fxyz,cxyz);

	    if (!P_InPolyFace3D(cxyz,3,fxyz,0.0,0)) {
	    
	      /* Just use the geometric center */

	      if (firstwarn2) {
		fprintf(stderr,"Circumcenter not in triangle!\n");
		fprintf(stderr,"Using geometric center instead\n");
		fprintf(stderr,"Expect non-planar faces\n");
		firstwarn2=0;
	      }

	      for (i = 0; i < 3; i++)
		cxyz1[i] = (fxyz[0][i]+fxyz[1][i]+fxyz[2][i])/3.0;	    	    

	      /* converged = 0; it = 0; */
	      converged = 1; it = 0; cxyz[0]=cxyz1[0];cxyz[1]=cxyz1[1];cxyz[2]=cxyz1[2];

	      while (!converged) {
		for (i = 0; i < 3; i++)
		  cxyz2[i] = (cxyz[i]+cxyz1[i])/2.0;
		if (P_InPolyFace3D(cxyz2,3,fxyz,0.0,0)) {
		  for (i = 0; i < 3; i++) cxyz1[i] = cxyz2[i];
		}
		else {
		  for (i = 0; i < 3; i++) cxyz[i] = cxyz2[i];
		}
		dist2 = (cxyz[0]-cxyz1[0])*(cxyz[0]-cxyz1[0]) +
		  (cxyz[1]-cxyz1[1])*(cxyz[1]-cxyz1[1]) +
		  (cxyz[2]-cxyz1[2])*(cxyz[2]-cxyz1[2]);

		if (dist2 < 1.0e-16) converged = 1;
		it++;
		if (it > 10) break;
	      }

	    }
	  }

	  /* Make vertex at circumcenter of triangle or at the
	     geometric center of the triangle if the circumcenter is
	     outside the triangle */

	  fv_vm[nfv_vm] = MV_New(vmesh);
	  MV_Set_Coords(fv_vm[nfv_vm],cxyz);
	  MV_Set_GEntDim(fv_vm[nfv_vm],2);
	  MV_Set_GEntID(fv_vm[nfv_vm],MF_GEntID(fcur_dm));

	  MEnt_Set_AttVal(fcur_dm,vmvatt,0,0.0,fv_vm[nfv_vm]);
	}

	nfv_vm++;
	
	
	
	/* Add modified Voronoi face corresponding to Delaunay edge */
	
	ef_vm[nef_vm] = MF_New(vmesh);
	MF_Set_Vertices(ef_vm[nef_vm], nfv_vm, fv_vm);
	
	MF_Set_GEntDim(ef_vm[nef_vm], MREGION);
	MF_Set_GEntID(ef_vm[nef_vm], MR_GEntID(rcur_dm));

	/* Attach primary mesh entity to dual mesh face */
	MEnt_Set_AttVal(ef_vm[nef_vm],dmentatt,0,0.0,e_dm);
	
	nef_vm++;

      } /* while ((f0_dm = List_Next_Entry(efaces,&idx2))) */

      /* MEnt_Set_AttVal(e_dm,vmfatt,0,0.0,ef_vm);	*/
      MEnt_Set_AttVal(e_dm,vmnfatt,nef_vm,0.0,NULL);

      break;



    case 3:      

      f0_dm = List_Entry(efaces_dm,0);
      fcur_dm = f0_dm;
      fnext_dm = 0;

      nfv_vm = 0;
      done = 0;
      while (!done) { 

	/* Get tet "in front of the face" when the face is oriented to 
	   use the edge in a positive sense */

	fedir = MF_EdgeDir(fcur_dm,e_dm);
	rcur_dm = MF_Region(fcur_dm,fedir);

	MEnt_Get_AttVal(rcur_dm,vmvatt,&ival,&rval,&pval);

	/* Is there already a Voronoi vertex associated with this tet? */

	fv_vm[nfv_vm] = (MVertex_ptr) pval;

	if (!fv_vm[nfv_vm]) {

	  MR_Coords(rcur_dm,&nrv_dm,rxyz);       /* nrv_dm = 4 */

	  if (use_centroids_only) {
	    for (i = 0; i < 3; i++)
	      cxyz[i] = (rxyz[0][i]+rxyz[1][i]+rxyz[2][i]+rxyz[3][i])/4.0;
	  }
	  else {
	    Tet_CircumCen(rxyz,cxyz);

	    if (!P_InTet(cxyz,rxyz,0.0,0)) {

	      /* Circumcenter not in tet - ok if this tet has a face on
		 the boundary but otherwise the mesh is not Delaunay */

	      int inttet=1;  /* assume interior tet */
	      MFace_ptr fdum;

	      rfaces_dm = MR_Faces(rcur_dm);
	      for (i = 0; i < 4; i++) {
		fdum = List_Entry(rfaces_dm,i);
		if (MF_GEntDim(fdum) == 2) {
		  inttet = 0;
		  break;
		}
	      }
	      List_Delete(rfaces_dm);

	      if (!inttet || !del) {
		if (firstwarn1) {
		  fprintf(stderr,"Circumcenter not inside tet\n");
		  fprintf(stderr,"Using geometric center of tet\n");
		  fprintf(stderr,"Expect non-planar faces\n");
		  firstwarn1=0;
		}

		for (i = 0; i < 3; i++)
		  cxyz1[i] = (rxyz[0][i]+rxyz[1][i]+rxyz[2][i]+rxyz[3][i])/4.0;

		/* converged = 0; it = 0; */
		converged = 1; it = 0; cxyz[0]=cxyz1[0];cxyz[1]=cxyz1[1];cxyz[2]=cxyz1[2];
		while (!converged) {
		  for (i = 0; i < 3; i++)
		    cxyz2[i] = (cxyz[i]+cxyz1[i])/2.0;
		  if (P_InTet(cxyz2,rxyz,0.0,0)) {
		    for (i = 0; i < 3; i++) cxyz1[i] = cxyz2[i];
		  }
		  else {
		    for (i = 0; i < 3; i++) cxyz[i] = cxyz2[i];
		  }
		  dist2 = (cxyz[0]-cxyz1[0])*(cxyz[0]-cxyz1[0]) +
		    (cxyz[1]-cxyz1[1])*(cxyz[1]-cxyz1[1]) +
		    (cxyz[2]-cxyz1[2])*(cxyz[2]-cxyz1[2]);
		
		  if (dist2 < 1.0e-16) converged = 1;
		  it++;
		  if (it > 10) break;
		}
	      }
	    }
	  }

	  /* Make vertex at circumcenter of tet or geometric center
	     of tet if the cirumcenter is outside the tet */	  
	    
	  fv_vm[nfv_vm] = MV_New(vmesh);        
	  MV_Set_Coords(fv_vm[nfv_vm],cxyz);
	  MV_Set_GEntDim(fv_vm[nfv_vm],3);
	  MV_Set_GEntID(fv_vm[nfv_vm],MR_GEntID(rcur_dm));

	  MEnt_Set_AttVal(rcur_dm,vmvatt,0,0.0,fv_vm[nfv_vm]);

	}
	nfv_vm++;


	rfaces_dm = MR_Faces(rcur_dm);

	found = 0; i = 0;
	while (i < 4 && !found) {
	  f_dm = List_Entry(rfaces_dm,i);
	  if (f_dm != fcur_dm && MEnt_IsMarked(f_dm,markid)) { 
	    fnext_dm = f_dm;
	    found = 1;
	  }
	  i++;
	}

	List_Delete(rfaces_dm);
	  
	if (fnext_dm == f0_dm) done = 1;
	fcur_dm = fnext_dm;
      }
	    
	
      /* Add Voronoi face corresponding to Delaunay edge */
	
      ef_vm[nef_vm] = MF_New(vmesh);
      MF_Set_Vertices(ef_vm[nef_vm], nfv_vm, fv_vm);

      MF_Set_GEntDim(ef_vm[nef_vm], MREGION);
      MF_Set_GEntID(ef_vm[nef_vm], MR_GEntID(rcur_dm));

      /* Attach primary mesh entity to dual mesh face */
      MEnt_Set_AttVal(ef_vm[nef_vm],dmentatt,0,0.0,e_dm);

      nef_vm++;
      MEnt_Set_AttVal(e_dm,vmnfatt,nef_vm,0.0,NULL); /* nef_vm will be 1 */

      break;
    }

    List_Unmark(efaces_dm,markid);
    List_Delete(efaces_dm);

  }



  /* Now create special capping faces on the boundaries of the
     Delaunay mesh which will cap the normally semi-infinite Voronoi
     cells on the boundary */
  
  
  
  idx1 = 0;
  while ((v_dm = MESH_Next_Vertex(dmesh,&idx1))) {
    
    
    nvf_vm = 0;                   /* Number of dual faces for Delaunay vertex */
    
    vedges_dm = MV_Edges(v_dm);
    
    
    switch (MV_GEntDim(v_dm)) {
      
    case 3: /* Interior vertex - no special faces to be created */
      break;
    case 2: /* Vertex on model face */
      
      /* Make a face connecting the circumcenters of all the connected
	 faces that are on the model face */
      
      vf_vm = (MEdge_ptr *) malloc(maxnvf_vm*sizeof(MEdge_ptr));
      
      idx2 = 0;
      while ((e0_dm = List_Next_Entry(vedges_dm,&idx2))) {

	if (ME_GEntDim(e0_dm) == 3) continue; /* internal edge */

	/* e0_dm is on the model face. Starting from e0, walk through
	   the faces from 'vfaces' that are on this model face and
	   connect their circumcenters to form a boundary face */

	ecur_dm = e0_dm;
	enext_dm = 0;
	fcur_dm = 0;

	efaces_dm = ME_Faces(ecur_dm);
	found = 0;
	idx3 = 0;
	while ((ef = List_Next_Entry(efaces_dm,&idx3))) {
	  if (MF_GEntDim(ef) == 2 && ef != fcur_dm) {
	      
	    /* Pick the correct face so that we loop around the
	       vertex consistent with the face normals */
	    
	    fedir = MF_EdgeDir(ef,e0_dm);
	    if (ME_Vertex(e0_dm,!fedir) == v_dm) {		
	      found = 1;
	      fcur_dm = ef;
	      break;
	    }
	  }
	}
	List_Delete(efaces_dm);	  	
	
	nfv_vm = 0;
	done = 0;
	while (!done) {
	  	  
	  MEnt_Get_AttVal(fcur_dm,vmvatt,&ival,&rval,&pval);
	  
	  /* Is there already a Voronoi vertex associated with this
	     boundary face ? If not create it */
	  
	  fv_vm[nfv_vm] = (MVertex_ptr) pval;
	  
	  if (!fv_vm[nfv_vm]) {
	    
	    MF_Coords(fcur_dm,&nfv_dm,fxyz);  /* nfv_dm = 3 */

	    if (use_centroids_only) {
	      for (i = 0; i < 3; i++)
		cxyz[i] = (fxyz[0][i]+fxyz[1][i]+fxyz[2][i])/3.0;
	    }
	    else {

	      Tri_CircumCen(fxyz,cxyz);
	    
	      if (!P_InPolyFace3D(cxyz,3,fxyz,0.0,0)) {
	      
		/* Just use the geometric center */
	     
		if (firstwarn2) {
		  fprintf(stderr,"Circumcenter not in triangle!\n");
		  fprintf(stderr,"Using geometric center instead\n");
		  fprintf(stderr,"Expect non-planar faces\n");
		  firstwarn2 = 0;
		}
	      
		for (i = 0; i < 3; i++)
		  cxyz1[i] = (fxyz[0][i]+fxyz[1][i]+fxyz[2][i])/3.0;

		/* converged = 0; it = 0; */
		converged = 1; it = 0; cxyz[0]=cxyz1[0];cxyz[1]=cxyz1[1];cxyz[2]=cxyz1[2];

		while (!converged) {
		  for (i = 0; i < 3; i++)
		    cxyz2[i] = (cxyz[i]+cxyz1[i])/2.0;
		  if (P_InPolyFace3D(cxyz2,3,fxyz,0.0,0)) {
		    for (i = 0; i < 3; i++) cxyz1[i] = cxyz2[i];
		  }
		  else {
		    for (i = 0; i < 3; i++) cxyz[i] = cxyz2[i];
		  }
		  dist2 = (cxyz[0]-cxyz1[0])*(cxyz[0]-cxyz1[0]) +
		    (cxyz[1]-cxyz1[1])*(cxyz[1]-cxyz1[1]) +
		    (cxyz[2]-cxyz1[2])*(cxyz[2]-cxyz1[2]);
		
		  if (dist2 < 1.0e-16) converged = 1;
		  it++;
		  if (it > 10) break;
		}
	      
	      }
	    }
	    
	    /* Make vertex at circumcenter of triangle or at the
	       geometric center of the triangle if the circumcenter is
	       outside the triangle */
	    
	    fv_vm[nfv_vm] = MV_New(vmesh);
	    MV_Set_Coords(fv_vm[nfv_vm],cxyz);
	    MV_Set_GEntDim(fv_vm[nfv_vm],2);
	    MV_Set_GEntID(fv_vm[nfv_vm],MF_GEntID(fcur_dm));
	    
	    MEnt_Set_AttVal(fcur_dm,vmvatt,0,0.0,fv_vm[nfv_vm]);
	    
	  }
	  nfv_vm++;
	  
	  
	  /* Get the edges of the faces in direction 'dir', such that
	     the first edge in the list will be ecur_dm and the last
	     edge will be the other edge of the face connected to the
	     vertex v_dm. This depends on how the ecur_dm is used in the
	     face '(edir') and whether v_dm is the first or second
	     vertex of the edge */

	  fedir = MF_EdgeDir(fcur_dm,ecur_dm);
	  evidx = (ME_Vertex(ecur_dm,0) == v_dm) ? 0 : 1;
	  dir = evidx^fedir;
	  
	  fedges_dm = MF_Edges(fcur_dm,dir,v_dm);
	  nfe_dm = List_Num_Entries(fedges_dm);
	  
	  enext_dm = List_Entry(fedges_dm,nfe_dm-1);
	  
	  List_Delete(fedges_dm);
	  
	  if (enext_dm == e0_dm)
	    done = 1;
	  else {
	    ecur_dm = enext_dm;

	    efaces_dm = ME_Faces(ecur_dm);
	    evidx = (ME_Vertex(ecur_dm,0) == v_dm) ? 0 : 1;
	    found = 0;
	    idx3 = 0;
	    while ((ef = List_Next_Entry(efaces_dm,&idx3))) {
	      if (MF_GEntDim(ef) == 2 && ef != fcur_dm) {
		found = 1;
		fcur_dm = ef;
		break;
	      }
	    }
	    List_Delete(efaces_dm);
	  }
	}  /* while (!done) */
	
	if (done) break;
      } /* while (e0_dm = ....) */
      
      
      vf_vm[nvf_vm] = MF_New(vmesh);
      MF_Set_Vertices(vf_vm[nvf_vm], nfv_vm, fv_vm);
      
      MF_Set_GEntDim(vf_vm[nvf_vm], MFACE);
      MF_Set_GEntID(vf_vm[nvf_vm], MV_GEntID(v_dm));
      nvf_vm++;
      
      MEnt_Set_AttVal(v_dm,vmfatt,0,0.0,vf_vm);
      MEnt_Set_AttVal(v_dm,vmnfatt,nvf_vm,0.0,NULL); /* only 1 face */

      /* Attach primary mesh entity to dual mesh face */
      MEnt_Set_AttVal(vf_vm[0],dmentatt,0,0.0,v_dm);
      
      break;

      
    case 1: case 0: /* Vertex on model edge or model vertex */
      
      vf_vm = (MEdge_ptr *) malloc(maxnvf_vm*sizeof(MEdge_ptr));
      
      outerdone = 0;
      while (!outerdone) {
	
	nnewf = 0;

	idx2 = 0;
	while ((e0_dm = List_Next_Entry(vedges_dm,&idx2))) {
	  
	  if (ME_GEntDim(e0_dm) != 1) continue;  /* want edge on model edge */
	  
	  efaces_dm = ME_Faces(e0_dm);
	  idx3 = 0;
	  fcur_dm = 0;
	  while ((ef = List_Next_Entry(efaces_dm,&idx3))) {
	    if (MF_GEntDim(ef) == 2 && !MEnt_IsMarked(ef,markid)) {
	      
	      /* Pick the correct face so that we loop around the
		 vertex consistent with the face normals */
	      
	      fedir = MF_EdgeDir(ef,e0_dm);
	      if (ME_Vertex(e0_dm,!fedir) == v_dm) {		
		found = 1;
		fcur_dm = ef;
		break;
	      }
	      
	    }
	  }
	  List_Delete(efaces_dm);
	  
	  if (fcur_dm == 0) continue; /* Couldn't find connected face
					 that has not been processed */
	  
	  f0_dm = fcur_dm;
	  nfv_vm = 0;
	  
	  
	  /* Add vertex into dual mesh polygon */
	  
	  MEnt_Get_AttVal(v_dm,vmvatt,&ival,&rval,&pval);
	  fv_vm[nfv_vm] = pval;
	  
	  if (!fv_vm[nfv_vm]) {
	    fv_vm[nfv_vm] = MV_New(vmesh);
	    MV_Coords(v_dm,cxyz);
	    MV_Set_Coords(fv_vm[nfv_vm],cxyz);
	    MV_Set_GEntDim(fv_vm[nfv_vm],MV_GEntDim(v_dm));
	    MV_Set_GEntID(fv_vm[nfv_vm],MV_GEntID(v_dm));
	    MEnt_Set_AttVal(v_dm,vmvatt,0,0.0,fv_vm[nfv_vm]);
	  }
	  nfv_vm++;
	  
	  
	  /* Add mid-point of edge into dual mesh polygon */
	  
	  MEnt_Get_AttVal(e0_dm,vmvatt,&ival,&rval,&pval);
	  fv_vm[nfv_vm] = pval;
	  
	  if (!fv_vm[nfv_vm]) {
	    fv_vm[nfv_vm] = MV_New(vmesh);
	    MV_Coords(ME_Vertex(e0_dm,0),exyz[0]);
	    MV_Coords(ME_Vertex(e0_dm,1),exyz[1]);
	    for (i = 0; i < 3; i++) cxyz[i] = (exyz[0][i]+exyz[1][i])/2.0;
	    MV_Set_Coords(fv_vm[nfv_vm],cxyz);
	    MV_Set_GEntDim(fv_vm[nfv_vm],ME_GEntID(e0_dm));
	    MV_Set_GEntID(fv_vm[nfv_vm],ME_GEntID(e0_dm));
	    MEnt_Set_AttVal(e0_dm,vmvatt,0,0.0,fv_vm[nfv_vm]);
	  }
	  nfv_vm++;
	  
	  ecur_dm = e0_dm;
	  enext_dm = 0;
	  
	  done = 0;
	  while (!done) {
	    
	    MEnt_Mark(fcur_dm,markid);
	    
	    MEnt_Get_AttVal(fcur_dm,vmvatt,&ival,&rval,&pval);
	    fv_vm[nfv_vm] = pval;
	    
	    if (!fv_vm[nfv_vm]) {
	      
	      MF_Coords(fcur_dm,&nfv_dm,fxyz);
	      

	      if (use_centroids_only) {
		for (i = 0; i < 3; i++)
		  cxyz[i] = (rxyz[0][i]+rxyz[1][i]+rxyz[2][i]+rxyz[3][i])/4.0;
	      }
	      else {

		Tri_CircumCen(fxyz,cxyz);
	      
		if (!P_InPolyFace3D(cxyz,3,fxyz,0.0,0)) {
		
		  /* circumcenter not in triangle. Use geometric center */

		  if (firstwarn2) {
		    fprintf(stderr,"Circumcenter not in triangle!\n");
		    fprintf(stderr,"Using geometric center instead\n");
		    fprintf(stderr,"Expect non-planar faces\n");
		    firstwarn2=0;
		  }
		
		  for (i = 0; i < 3; i++)
		    cxyz1[i] = (fxyz[0][i]+fxyz[1][i]+fxyz[2][i])/3.0;
		
		  /* converged = 0; it = 0; */
		  converged = 1; it = 0; cxyz[0]=cxyz1[0];cxyz[1]=cxyz1[1];cxyz[2]=cxyz1[2];
		  while (!converged) {
		    for (i = 0; i < 3; i++)
		      cxyz2[i] = (cxyz[i]+cxyz1[i])/2.0;
		    if (P_InPolyFace3D(cxyz2,3,fxyz,0.0,0)) {
		      for (i = 0; i < 3; i++) cxyz1[i] = cxyz2[i];
		    }
		    else {
		      for (i = 0; i < 3; i++) cxyz[i] = cxyz2[i];
		    }
		    dist2 = (cxyz[0]-cxyz1[0])*(cxyz[0]-cxyz1[0]) +
		      (cxyz[1]-cxyz1[1])*(cxyz[1]-cxyz1[1]) +
		      (cxyz[2]-cxyz1[2])*(cxyz[2]-cxyz1[2]);
		  
		    if (dist2 < 1.0e-16) converged = 1;
		    it++;
		    if (it > 10) break;
		  }

		}
	      }
	      
	      /* Make vertex at circumcenter of triangle or at the
		 geometric center of ther triangle if the circumcenter
		 is outside the triangle */
	      
	      fv_vm[nfv_vm] = MV_New(vmesh);
	      MV_Set_Coords(fv_vm[nfv_vm],cxyz);
	      MV_Set_GEntDim(fv_vm[nfv_vm],MF_GEntDim(fcur_dm));
	      MV_Set_GEntID(fv_vm[nfv_vm],MF_GEntID(fcur_dm));
	      
	      
	      MEnt_Set_AttVal(fcur_dm,vmvatt,0,0.0,fv_vm[nfv_vm]);
	      
	    } /* if (!fv_vm[nfv_vm]) */
	    
	    nfv_vm++;
	    
	    
	    evidx = (ME_Vertex(ecur_dm,0) == v_dm) ? 0 : 1;
	    fedir = MF_EdgeDir(fcur_dm,ecur_dm);
	    dir = evidx^fedir;
	    
	    fedges_dm = MF_Edges(fcur_dm,dir,v_dm);
	    enext_dm = List_Entry(fedges_dm,2);
	    List_Delete(fedges_dm);
	    
	    
	    /* If we reached another edge on a model edge, we are done */
	    
	    if (ME_GEntDim(enext_dm) == 1) done = 1;
	    
	    ecur_dm = enext_dm;
	    
	    efaces_dm = ME_Faces(ecur_dm);
	    idx3 = 0;
	    while ((ef = List_Next_Entry(efaces_dm,&idx3))) {
	      if (MF_GEntDim(ef) == 2 && !MEnt_IsMarked(ef,markid)) {
		fcur_dm = ef;
		break;
	      }	    
	    }
	    List_Delete(efaces_dm);
	    
	  }

	  /* Add mid-point of edge into dual mesh polygon */
	  
	  MEnt_Get_AttVal(ecur_dm,vmvatt,&ival,&rval,&pval);
	  fv_vm[nfv_vm] = pval;
	  
	  if (!fv_vm[nfv_vm]) {
	    fv_vm[nfv_vm] = MV_New(vmesh);
	    MV_Coords(ME_Vertex(ecur_dm,0),exyz[0]);
	    MV_Coords(ME_Vertex(ecur_dm,1),exyz[1]);
	    for (i = 0; i < 3; i++) cxyz[i] = (exyz[0][i]+exyz[1][i])/2.0;
	    MV_Set_Coords(fv_vm[nfv_vm],cxyz);
	    MV_Set_GEntDim(fv_vm[nfv_vm],ME_GEntID(ecur_dm));
	    MV_Set_GEntID(fv_vm[nfv_vm],ME_GEntID(ecur_dm));
	    MEnt_Set_AttVal(ecur_dm,vmvatt,0,0.0,fv_vm[nfv_vm]);
	  }
	  nfv_vm++;
	  
	  
	  /* Make the polygonal face from the vertices */
	  
	  vf_vm[nvf_vm] = MF_New(vmesh);
	  MF_Set_Vertices(vf_vm[nvf_vm], nfv_vm, fv_vm);
	  
	  MF_Set_GEntDim(vf_vm[nvf_vm], MFACE);
	  MF_Set_GEntID(vf_vm[nvf_vm], MF_GEntID(f0_dm));
	
	  /* Attach primary mesh entity to dual mesh face */
	  MEnt_Set_AttVal(vf_vm[nvf_vm],dmentatt,0,0.0,v_dm);
  
	  nvf_vm++;
	  nnewf++;
	  
	} /* while (e0_dm ....) */

	if (nnewf == 0) outerdone = 1;
      } /* while (!outerdone) */

      MEnt_Set_AttVal(v_dm,vmfatt,0,0.0,vf_vm);	  
      MEnt_Set_AttVal(v_dm,vmnfatt,nvf_vm,0.0,NULL);


      vfaces_dm = MV_Faces(v_dm);
      List_Unmark(vfaces_dm,markid);
      List_Delete(vfaces_dm);

      break;
    }

    List_Delete(vedges_dm);

  } /* while ((v_dm = ....)) */




  
  /* /\*** Collapse out degenerate edges and faces ***\/ */

  /* /\*** HAVE TO INFORM THE PRIMARY MESH ENTITIES THAT SOME ENTITIES  */
  /*      GOT DELETED ***\/ */


  /* idx1 = 0; */
  /* while ((e_vm = MESH_Next_Edge(vmesh,&idx1))) { */
  /*   MVertex_ptr ev[2]; */
  /*   MEdge_ptr vedge; */
  /*   double len2; */
  /*   int gdim[2], gid[2]; */

  /*   ev[0] = ME_Vertex(e_vm,0); */
  /*   ev[1] = ME_Vertex(e_vm,1); */

  /*   MV_Coords(ev[0],exyz[0]); */
  /*   MV_Coords(ev[1],exyz[1]); */

  /*   len2 = PP_Dist2(exyz[0],exyz[1]); */

  /*   if (len2 > 0.0) continue; */

  /*   /\* Degenerate edge which we have to collapse *\/ */

  /*   gdim[0] = MV_GEntDim(ev[0]); */
  /*   gid[0] = MV_GEntID(ev[0]); */
  /*   gdim[1] = MV_GEntDim(ev[1]); */
  /*   gid[1] = MV_GEntID(ev[1]); */

  /*   if (gdim[0] == gdim[1] && gid[0] != gid[1]) { */
  /*     fprintf(stderr,"Cannot collapse out zero length edge due to topological considerations\n"); */
  /*     continue; */
  /*   } */


  /*   /\* for simplicity of coding we want to always collapse ev[0] to */
  /*      ev[1] (i.e. replace ev[0] with ev[1]) *\/ */
  /*   /\* for topological validity, we always want to collapse the vertex */
  /*      classified on the higher dimension entity to the vertex */
  /*      classified on the lower dimension entity *\/ */
  /*   /\* So if gdim[1] > gdim[0], switch the two vertices around *\/ */

  /*   if (gdim[1] > gdim[0]) { */
  /*     MVertex_ptr vtmp; */
  /*     vtmp = ev[0]; */
  /*     ev[0] = ev[1]; */
  /*     ev[1] = vtmp; */
  /*   } */

  /*   vedges_vm = MV_Edges(ev[0]); */
  /*   idx2 = 0; */
  /*   while ((vedge = List_Next_Entry(vedges_vm,&idx2))) { */
  /*     ME_Replace_Vertex(vedge,ev[0],ev[1]); */
  /*   } */
  /*   List_Delete(vedges_vm); */


  /*   efaces_vm = ME_Faces(e_vm); */

  /*   idx2 = 0; */
  /*   while ((ef = List_Next_Entry(efaces_vm,&idx2))) { */
  /*     List_ptr fedges; */
  /*     int nfe, found; */
  /*     MEdge_ptr fedge, olde[3], nue[2]; */

  /*     fedges = MF_Edges(ef,1,0); */
  /*     nfe = List_Num_Entries(fedges); */

  /*     if (nfe == 2) { */
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
  /*     } */

  /*     found = 0; */
  /*     idx3 = 0; */
  /*     for (i = 0; i < nfe; i++) { */
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
  /*     } */

  /*     List_Delete(fedges); */

  /*   } /\* while (ef ....) *\/ */

  /*   ME_Delete(e_vm,0); */
  /*   if (ev[0] != ev[1]) */
  /*     MV_Delete(ev[0],0); */

  /*   List_Delete(efaces_vm); */

  /* } /\* while (e_vm ....) *\/ */


  /* /\* If there are any faces with two long coincident edges because */
  /*    degenerate edges between them got deleted, delete them *\/ */

  /* idx1 = 0; */
  /* while ((f_vm = MESH_Next_Face(vmesh,&idx1))) { */
  /*   int nfe; */
  /*   List_ptr fedges; */
    
  /*   fedges = MF_Edges(f_vm,1,0); */
  /*   nfe = List_Num_Entries(fedges); */

  /*   if (nfe == 2) { */
  /*     MEntity_ptr ment; */
      
  /*     /\* Face has become degenerate *\/ */


  /*     /\* Tell primary mesh that this face has been deleted *\/ */

  /*     MEnt_Get_AttVal(f_vm,dmentatt,&ival,&rval,&ment); */

  /*     MEnt_Get_AttVal(ment,vmnfatt,&nef_vm,&rval,&pval); */
  /*     MEnt_Get_AttVal(ment,vmfatt,&ival,&rval,&pval); */
  /*     ef_vm = (MFace_ptr *) pval; */

  /*     found = 0; */
  /*     for (i = 0; i < nef_vm; i++) { */
  /* 	if (ef_vm[i] == f_vm) { */
  /* 	  found = 1; */
  /* 	  for (j = i; j < nef_vm-1; j++) */
  /* 	    ef_vm[i] = ef_vm[i+1]; */
  /* 	  nef_vm--; */
  /* 	  MEnt_Set_AttVal(ment,vmnfatt,nef_vm,0.0,NULL); */
  /* 	  break; */
  /* 	} */
  /*     } */
  /*     if (!found) { */
  /* 	fprintf(stderr,"Dual mesh face not attached to correct primary mesh entity\n"); */
  /*     } */

  /*     /\* Delete the face *\/ */

  /*     MF_Delete(f_vm,0); */

      
  /*     /\* Merge the edges *\/ */

  /*     MEs_Merge(List_Entry(fedges,0),List_Entry(fedges,1)); */

  /*     List_Delete(fedges); */
  /*   } */
  /* } */



      
  /**** Make Dual Mesh Regions ****/



  /* Now make polyhedral regions of the dual. Each interior vertex
     will have one corresponding polyhedral region in the dual. Each
     boundary vertex in a 2-manifold model will also have one polyhedral
     region in the dual. A boundary vertex on a non-manifold boundary
     may have more than one region in the dual. */

  rfarray = (MFace_ptr *) malloc(maxnrf_vm*sizeof(MFace_ptr));
  rfdirs = (int *) malloc(maxnrf_vm*sizeof(int));

  idx1 = 0;
  while ((v_dm = MESH_Next_Vertex(dmesh,&idx1))) {

    vedges_dm = MV_Edges(v_dm);

    
    switch (MV_GEntDim(v_dm)) {

    case 3: 
      /* Interior vertex */
      /* Make a dual mesh region by collecting the unique dual mesh
	 faces associated with the edges of this vertex */

      idx2 = 0;
      nrf_vm = 0;
      while ((e_dm = List_Next_Entry(vedges_dm,&idx2))) {

	MEnt_Get_AttVal(e_dm,vmnfatt,&nef_vm,&rval,&pval);
	if (nef_vm == 0)
	  continue;

	MEnt_Get_AttVal(e_dm,vmfatt,&ival,&rval,&pval);

	ef_vm = (MFace_ptr *) pval;
	f_vm = ef_vm[0];
	rfarray[nrf_vm] = f_vm;

	/* Because of the way the Voronoi face was created, its
	   "normal" will always be aligned with the direction of the
	   edge. This means that the Voronoi region associated with v0
	   will use the face in the positive sense and the Voronoi
	   region associated with v1 in the negative sense */

	dir = (ME_Vertex(e_dm,0) == v_dm);

	rfdirs[nrf_vm] = dir;
	nrf_vm++;
      }

      r_vm = MR_New(vmesh);
      MR_Set_Faces(r_vm,nrf_vm,rfarray,rfdirs);

      MR_Set_GEntID(r_vm,MV_GEntID(v_dm));

      break;

    case 0: case 1: case 2:

      /* First check if this is a 2-manifold case, simple non-manifold
	 case (model face with different model regions on either
	 side), or complex non-manifold case (model face with same
	 model region on either side) */

      vfaces_dm = MV_Faces(v_dm);
      topocase = 0; /* assume the default is 2-manifold case */
      idx2 = 0;
      while ((vf_dm = List_Next_Entry(vfaces_dm,&idx2))) {
	MRegion_ptr reg0, reg1;
	int greg0, greg1;

	if (MF_GEntDim(vf_dm) == 3) continue;

	reg0 = MF_Region(vf_dm,0);
	greg0 = reg0 ? MR_GEntID(reg0) : 0;
	reg1 = MF_Region(vf_dm,1);
	greg1 = reg1 ? MR_GEntID(reg1) : 0;

	if (greg0 && greg1) {
	  topocase = 1;  /* non-manifold */
	  if (reg0 == reg1) {
	    topocase = 2;
	    break;
	  }
	}
      }
      List_Delete(vfaces_dm);


      switch (topocase) {
      case 0: case 1:

	vregions_dm = MV_Regions(v_dm);
	idx2 = 0;
	ngreg = 0;
	while ((rcur_dm = List_Next_Entry(vregions_dm,&idx2))) {
	  
	  greg0 = MR_GEntID(rcur_dm);
	  found = 0;
	  for (i = 0; i < ngreg; i++) {
	    if (gregions[i] == greg0) {
	      found = 1;
	      break;
	    }
	  }
	  
	  if (!found) {
	    gregions[ngreg] = greg0;
	    ngreg++;
	  }
	}
	
	
	/* Collect dual faces associated with this vertex and its
	   edges that are in each model region and make mesh regions
	   out of them */
	
	for (ir = 0; ir < ngreg; ir++) {

	  nrf_vm = 0;
	  
	  /* collect primary mesh regions in this model region */
	  
	  vregions_sub = List_New(0);
	  idx3 = 0;
	  while ((rcur_dm = List_Next_Entry(vregions_dm,&idx3)))
	    if (MR_GEntID(rcur_dm) == gregions[ir])
	      List_Add(vregions_sub,rcur_dm);
	  
	  /* Make a unique list of edges belonging to these mesh regions */
	  
	  vedges_sub = List_New(0);
	  idx3 = 0;
	  while ((rcur_dm = List_Next_Entry(vregions_sub,&idx3))) {
	    redges_dm = MR_Edges(rcur_dm);
	    
	    idx4 = 0;
	    while ((e_dm = List_Next_Entry(redges_dm,&idx4))) {
	      if (!MEnt_IsMarked(e_dm,markid) &&
		  (ME_Vertex(e_dm,0) == v_dm || ME_Vertex(e_dm,1) == v_dm)) {
		List_Add(vedges_sub,e_dm);
		MEnt_Mark(e_dm,markid);
	      }
	    }
	    List_Delete(redges_dm);
	  }
	  List_Unmark(vedges_sub, markid);
	  List_Delete(vregions_sub);
	  
	  

	  /* First add dual mesh faces associated with internal edges
	     in this list */

	  idx3 = 0;
	  while ((e_dm = List_Next_Entry(vedges_sub,&idx3))) {
	    
	    if (ME_GEntDim(e_dm) != 3) continue;

	    MEnt_Get_AttVal(e_dm,vmnfatt,&nef_vm,&rval,&pval);
	    if (nef_vm == 0)
	      continue;
	    
	    MEnt_Get_AttVal(e_dm,vmfatt,&ival,&rval,&pval);
	    
	    ef_vm = (MFace_ptr *) pval;
	    rfarray[nrf_vm] = f_vm = ef_vm[0]; /* interior edge - only 1 assoc. face */
	    
	/* Because of the way the Voronoi face was created, its
	   "normal" will always be aligned with the direction of the
	   edge. This means that the Voronoi region associated with v0
	   will use the face in the positive sense and the Voronoi
	   region associated with v1 in the negative sense */

	    dir = (ME_Vertex(e_dm,0) == v_dm);
	    
	    rfdirs[nrf_vm] = dir;
	    nrf_vm++;
	    
	    /* Mark all of the face's edges - need in the next loop */
	    
	    fedges_vm = MF_Edges(f_vm,1,0);
	    List_Mark(fedges_vm,markid);
	    List_Delete(fedges_vm);
	  }



	  /* Then add dual mesh faces associated with boundary edges
	     in this list */
	  
	  idx3 = 0;
	  while ((e_dm = List_Next_Entry(vedges_sub,&idx3))) {
	    
	    if (ME_GEntDim(e_dm) == 3) continue;
	    
	    MEnt_Get_AttVal(e_dm,vmnfatt,&nef_vm,&rval,&pval);
	    MEnt_Get_AttVal(e_dm,vmfatt,&ival,&rval,&pval);
	    ef_vm = (MFace_ptr *) pval;
	    
	    /* edge may have multiple dual faces associated with it */
	    /* Add ALL of them - later we will eliminate those that */
	    /* are not connected to another face in the list along  */
	    /* at least one edge                                    */
	    
	    found = 0;
	    for (j = 0; j < nef_vm; j++) {
	      f_vm = ef_vm[j];
	    
	      rfarray[nrf_vm] = f_vm;
	    
	      /* Because of the way the Voronoi face was created, its
		 "normal" will always be aligned with the direction of
		 the edge. This means that the Voronoi region
		 associated with v0 will use the face in the positive
		 sense and the Voronoi region associated with v1 in
		 the negative sense */

	      dir = (ME_Vertex(e_dm,0) == v_dm);

	      rfdirs[nrf_vm] = dir;
	      nrf_vm++;
	    }
	  } /* while (e_dm ....) */
	  
	  List_Delete(vedges_sub);


	  /* Finally **ALL** the dual mesh faces associated with the vertex */

	  MEnt_Get_AttVal(v_dm,vmnfatt,&nvf_vm,&rval,&pval);
	  MEnt_Get_AttVal(v_dm,vmfatt,&ival,&rval,&pval);
	  vf_vm = (MFace_ptr *) pval;
	  
	  for (i = 0; i < nvf_vm; i++) {
	    
	    f_vm = vf_vm[i];

	    rfarray[nrf_vm] = f_vm;
	    
	    rfdirs[nrf_vm] = -1; /* Special value to indicate it is not computed */
	    nrf_vm++;
	    
	  } /* for (i = nvf_vm; .....) */
	  


	  /* Now eliminate the dual mesh faces which are not connected
	     along all edges to another dual mesh face from the
	     rfarray faces */

	  for (i = 0; i < nrf_vm; i++) {
	    MEnt_Mark(rfarray[i],markid);
	  }

	  i = 0;
	  while (i < nrf_vm) {

	    fedges_vm = MF_Edges(rfarray[i],1,0);
	    List_Unmark(fedges_vm,markid);

	    idx4 = 0;
	    while ((e_vm = List_Next_Entry(fedges_vm,&idx4))) {
	      found = 0;
	      for (j = 0; j < nrf_vm; j++) {
		if (i == j) continue;

		if (MF_UsesEntity(rfarray[j],e_vm,MEDGE)) {
		  found = 1;
		  break;
		}
	      }
	      if (!found)
		break;
	    }
	    List_Delete(fedges_vm);

	    if (found) 
	      i++;
	    else {
	      
	      /* face whose edge is not connected to any other face in
		 rfarray. This implies that this face does not belong
		 to the mesh region being built. Happens only at
		 multimaterial vertices */
	      
	      for (j = i+1; j < nrf_vm; j++) {
		rfarray[j-1] = rfarray[j];
		rfdirs[j-1] = rfdirs[j];
	      }
	      nrf_vm--;

	      /* we won't increment i because rfarray[i] now contains
		 a face that has not yet been checked */
	    }
	  }

	  for (i = 0; i < nrf_vm; i++) {
	    MEnt_Unmark(rfarray[i],markid);
	  }

	  
	  /* Calculate point inside the new region */
	  
	  cxyz[0] = cxyz[1] = cxyz[2] = 0.0;
	  for (i = 0; i < nrf_vm; i++) {
	    MF_Coords(rfarray[i],&nfv_vm,fxyz);
	    PF_Center(nfv_vm,fxyz,cen);
	    for (j = 0; j < 3; j++) cxyz[j] += cen[j];
	  }
	  for (j = 0; j < 3; j++) cxyz[j] /= nrf_vm;
	  
	  
	  /* Now compute the direction in which the faces connected to
	     the primary mesh vertex should be used by comparing the
	     vector from the center of the region to the center of the
	     face and the normal of the face */
	  
	  for (i = 0; i < nrf_vm; i++) {
	    f_vm = rfarray[nrf_vm-i-1];
	    if (rfdirs[nrf_vm-i-1] >= 0) continue;
	    
	    /* Calculate the center of the face and also its normal */
	    
	    MF_Coords(f_vm,&nfv_vm,fxyz);
	    avenormal[0] = avenormal[1] = avenormal[2] = 0.0;
	    for (j = 0; j < nfv_vm; j++) {
	      VDiff3(fxyz[(j+2)%nfv_vm],fxyz[(j+1)%nfv_vm],vec1);
	      VDiff3(fxyz[j],fxyz[(j+1)%nfv_vm],vec2);
	      VCross3(vec1,vec2,normal);	  
	      VSum3(avenormal,normal,avenormal);
	    }
	    
	    PF_Center(nfv_vm,fxyz,cen);
	    
	    VDiff3(cen,cxyz,vec1);
	    dp = VDot3(vec1,avenormal);
	    
	    rfdirs[nrf_vm-i-1] = (dp > 0);
	  }
	  
	  /* Find an interior dual face associated with this vertex */
	  /* The dual mesh region associated with the vertex will be
	     classified on the same model region as the dual mesh
	     face */
	  
	  greg0 = 0;
	  for (i = 0; i < nrf_vm; i++) {
	    f_vm = rfarray[i];
	    if (!greg0 && MF_GEntDim(f_vm) == 3) {
	      greg0 = MF_GEntID(f_vm);
	      break;
	    }
	  }

	  r_vm = MR_New(vmesh);
	  MR_Set_Faces(r_vm,nrf_vm,rfarray,rfdirs);	  
	  
	  if (!greg0) {
	    fprintf(stderr,"Cannot find classification for dual mesh region\n");
	  }

	  MR_Set_GEntID(r_vm,greg0);

	} /* for (ir = 0; ir < ngreg; ir++) */
	
	break;

	case 2:

	  fprintf(stderr,"Non-manifold case with surface embedded in region!!\n");
	  fprintf(stderr,"NOT IMPLEMENTED YET! NEED GEOMETRIC CHECKS\n");

	  break;
      } /* switch (topocase) */ 

      List_Delete(vedges_dm);

    } /* switch (MV_GEntDim(v_dm)) */

  } /* while (v_dm = MESH_Next_Vertex(dmesh,&idx1)) */

  free(rfarray);
  free(rfdirs);

  free(fv_vm);
  free(ef_vm);
  free(vf_vm);
  


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
    }

  } /* while (e_vm ....) */

  
  return 1;
  
}
