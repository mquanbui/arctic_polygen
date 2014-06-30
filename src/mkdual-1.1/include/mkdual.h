#ifndef __MAKE_DUAL
#define __MAKE_DUAL
#endif


#ifdef __cplusplus
extern "C" {
#endif

#include "MSTK.h"



/* Convert a simplicial mesh into a dual mesh truncated by the
   boundaries of the original mesh.

   Triangular mesh  ---> Polygonal mesh
   Tetrahedral mesh ---> Polyhedral mesh
   Delaunay mesh    ---> "Voronoi" mesh

   Each tetrahedron of the primary mesh gives rise to one node of the
   dual mesh. If the circumcenter of the tetrahedron is inside the
   tetrahedron, the circumcenter is used as the node location; if not,
   then a point inside the tet on the line connecting the circumcenter
   and the centroid is used.

   As a result of the above algorithm for locating the dual nodes
   general tetrahedral meshes will result in polyhedra with curved
   faces. Even duals of Delaunay meshes will have curved faces if
   their circumcenters are not inside the tetrahedra (this will
   typically occur at the boundary in Delaunay meshes) */

  int MESH_MakeDual(Mesh_ptr primesh, Mesh_ptr dualmesh, int *opts);
  int MESH_MakeDual2(Mesh_ptr primesh, Mesh_ptr dualmesh, int *opts);
  int MESH_MakeDual3(Mesh_ptr primesh, Mesh_ptr dualmesh, int *opts);


#ifdef __cplusplus
}
#endif
