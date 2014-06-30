#include "MCGeom.h"

#ifdef __cplusplus
extern "C" {
#endif

  double PP_Dist2(double *xyz1, double *xyz2) {

    return ((xyz1[0]-xyz2[0])*(xyz1[0]-xyz2[0]) +
	    (xyz1[1]-xyz2[1])*(xyz1[1]-xyz2[1]) +
	    (xyz1[2]-xyz2[2])*(xyz1[2]-xyz2[2]));
  }

#ifdef __cplusplus
}
#endif
