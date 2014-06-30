#ifndef _H_DEL2VOR
#define _H_DEL2VOR

#ifdef __cplusplus
extern "C" {
#endif

#include "MSTK.h"


int convert_del2vor(Mesh_ptr delaunay_mesh, Mesh_ptr voronoi_mesh, int *options);


#ifdef __cplusplus
}
#endif

#endif

