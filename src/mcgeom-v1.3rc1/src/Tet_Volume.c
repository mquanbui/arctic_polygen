#include <stdio.h>
#include <math.h>
#include "MCGeom.h"


double Tet_Volume(double (*rxyz)[3]) {
  double vec01[3], vec02[3], vec03[3], n0[3], vol;

  /* Assume that vertices of tet are ordered so that they obey the
     right-hand rule. In other words, the normal of the face formed by
     the first three vertices points towards the fourth vertex */


  VDiff3(rxyz[1],rxyz[0],vec01);
  VDiff3(rxyz[2],rxyz[0],vec02);
  VDiff3(rxyz[3],rxyz[0],vec03);
  VCross3(vec01,vec02,n0);
  
  vol = VDot3(n0,vec03)/6.0;
  return vol;

}
