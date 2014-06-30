#ifndef _H_MCGeom
#define _H_MCGeom

#ifdef __cplusplus
extern "C" {
#endif 

#include <math.h>

#ifdef M_PI
#define PI M_PI
#elif M_PIl
#define PI M_PIl
#else
#define PI 3.14159265358979323846
#endif


  /* VECTOR FUNCTIONS */

  void VNeg3(double *a);
  void VDiff3(double *a, double *b, double *c);
  void VSum3(double *a, double *b, double *c);
  void VScale3(double *a, double s);
  double VLen3(double *a);
  double VLenSqr3(double *a);
  void VNormalize3(double *a);
  double VDot3(double *a, double *b);
  void VCross3(double *a, double *b, double *x);
  void VCopy3(double *destination, double *source);

  /* MATRIX FUNCTIONS */
  /* General Matrices - matrices must be declared as double ** */
  void MDiff(double **a1, double **a2, int nr, int nc, double **a3);
  void MSum(double **a1, double **a2, int nr, int nc, double **a3);
  void MMult(double **a1, int nr1, int nc1, double **a2, int nr2, int nc2, double **a3);
  double MDet(double **a, int n);
  void MTransp(double **a, int nr, int nc, double **atrans);
  void MInv(double **a, int n, double **ainv);

  /* 3x3 or smaller matrices - matrices must be declared as double [3][3] */
  /* Even if you are working with a 2x3 matrix or a 3x1 vector declare
     it as [3][3] and fill in the appropriate entries */
  void MDiff3(double a1[3][3], double a2[3][3], int nr, int nc, double a3[3][3]);
  void MSum3(double a1[3][3], double a2[3][3], int nr, int nc, double a3[3][3]);
  void MMult3(double a1[3][3], int nr1, int nc1, double a2[3][3], int nr2, int nc2, double a3[3][3]);
  void MVMult3(double mat[3][3], int nr1, int nc1, double vec[], int nr2, double resvec[]);
  double MDet3(double a[3][3], int n);
  void MTransp3(double a[3][3], int nr, int nc, double atrans[3][3]);
  void MInv3(double a[3][3], int n, double ainv[3][3]);


  /* Computational Geometry Stuff */

  /* Point functions */

  double PP_Dist2(double *xyz1, double *xyz2);

  /* Given a point, find the closest point on a line segment */
  double P_ClosestPntOnLine(double *pxyz, double *lxyz0, double *lxyz1, 
			    double eps, double *cxyz);
  double P_MinDist2Quad(double vxyz[], double fxyz[][3], double eps);
  double P_MinDist2Tri(double vxyz[], double fxyz[][3], double eps);

  /* Find projection of point to an infinite line passing through
     lxyz[0] and lxyz[1] */
  void P_Prj2Line(double pxyz[], double lxyz[][3], double prxyz[], double *distSQR);

  void P_Prj2TriPlane(double pxyz[3], double (*fxyz)[3], double prxyz[], double *dist);

  void P_ReflectGen(double p[], double v[], double n[], double rp[]);
  void P_RotateGen(double p[], double theta, double rp[], double ra[], 
		   double p1[]);
  int P_InPolyFace2D(double *vxyz, int np, double (*pxyz)[3], double tol,
		     int flag);
  int P_InPolyFace3D(double *vxyz, int np, double (*pxyz)[3], double tol,
		     int flag);
  int P_OnLineSeg(double *vxyz, double (*lxyz)[3], double tol);
  int P_InTet(double *pxyz, double (*txyz)[3], double tol, int flag);
  int P_InPolyhedron(double *xyz, int nf, int *nfp, double (**pxyz)[3], 
		     double tol, int flag);

  
  /* Intersection Functions */

  int  X_LineLine(double (*La)[3], double (*Lb)[3], double eps, int *nx, 
		  double (*xpt)[3]);
  int  X_LineFace(double (*lxyz)[3], int np, double (*pxyz)[3], double eps,
		  int *nx, double (*ixyz)[3], int *ie);
  int  X_LineCircle(double (*La)[3], double X0[3], double R, double eps, 
		    int *nx, double (*xpt)[3]);

  /* convex polygons only */
  int  X_FaceFace(int np1, double (*pxyz1)[3], int np2, double (*pxyz2)[3],
		  double ptol, int *nip, double (*ipxyz)[3], int *iptag);

  int PF_Subdivide(int np, double (*pxyz)[3], int nlp, double (*lxyz)[3],
                   int *ie, double eps, int *npl, double (*pxyzl)[3], 
		   int *ptagl, int *npr, double (*pxyzr)[3], int *ptagr);

  int PF_Slice(int np, double (*pxyz)[3], double (*lxyz)[3], int *ie,
	       double eps, 
	       int *npl, double (*pxyzl)[3], int *ptagl, 
	       int *npr, double (*pxyzr)[3], int *ptagr);

  /* Polygonal Face Functions */

  void PF_Center(int n, double (*fxyz)[3], double *cen);
  void PF_Centroid(int n, double (*fxyz)[3], double *cen);
  double PF_Area(int n, double (*fxyz)[3]);
  void PF_Angles(double (*fxyz)[3], int n, double *angles);
  void PF_CosAngles(double (*fxyz)[3], int n, double *cosangles);
  void PF_MinMaxCosAngles(double (*fxyz)[3], int n, double *mincos, double *maxcos);
  void PF_MinMaxAngles(double (*fxyz)[3], int n, double *minang, double *maxang);
  void PF_CondNums(double (*fxyz)[3], int n, double *condnums);
  void PF_CondNum_i(double (*fxyz)[3], int n, int i, double *condnum);
  void PF_MinMaxCondNums(double (*fxyz)[3], int n, double *mincn, double *maxcn);
  int  PF_LSSphere(double (*fxyz)[3], int n, double *pcen, double *rad);


  /* Polyhedral region functions */

  double PR_Volume(double (*rxyz)[3], int n, int **rfverts, int *nfv, int nf, int *star_shaped);
  void PR_Angles(double (*rxyz)[3], int n, double *angles, int **rfverts, int *nfv, int nf, int *ncos);
  void PR_CosAngles(double (*rxyz)[3], int n, double *cosangles, int **rfverts, int *nfv, int nf, int *ncos);
  void PR_MinMaxCosAngles(double (*rxyz)[3], int n, double *mincos, double *maxcos, int **rfverts, int *nfv, int nf);
  void PR_MinMaxAngles(double (*rxyz)[3], int n, double *minang, double *maxang, int **rfverts, int *nfv, int nf);
  void PR_CondNums(double (*rxyz)[3], int n, double *condnums, int (*reverts)[2], int ne);
  void PR_MinMaxCondNums(double (*rxyz)[3], int n, double *mincn, double *maxcn, int (*reverts)[2], int ne);


  /* Triangle Functions */

  double TriCotAngle(double txyz[][3], int idx);
  void Tri_CotAngles(double txyz[][3], double cotangs[]);
  double Tri_Area(double (*xyz)[3]);
  void   Tri_CircumCen(double (*xyz)[3], double *cen);
  void   Tri_Normal(double (*xyz)[3], double *normal);
  void   Tri_UnitNormal(double (*xyz)[3], double *unormal);
  
  /* Subdivision of quad into triangles based on dihedral angles */

  void Quad_BestTriangles(double fxyz[][3], int t1ind[], double t1xyz[][3], int t2ind[], double t2xyz[][3]);

  /* Tet Functions */

  double Tet_Volume(double (*rxyz)[3]);
  void Tet_Angles(double (*rxyz)[3], double *angles);
  void Tet_CosAngles(double (*rxyz)[3], double *cosangs);
  void Tet_MinMaxCosAngles(double (*rxyz)[3], double *mincos, double *maxcos);
  void Tet_MinMaxAngles(double (*rxyz)[3], double *minang, double *maxang);
  void   Tet_CircumCen(double (*xyz)[3], double *cen);

  /* Hex Functions */

  double Hex_Volume(double (*rxyz)[3]);
  void Hex_Angles(double (*rxyz)[3], double *angles);
  void Hex_CosAngles(double (*rxyz)[3], double *cosangs);
  void Hex_MinMaxCosAngles(double (*rxyz)[3], double *mincos, double *maxcos);
  void Hex_MinMaxAngles(double (*rxyz)[3], double *minang, double *maxang);



  /* Condition Number Shape Measures for corners */

  double Corner_CondNum2(double xyz0[3], double xyz1[3], double xyz2[3], double refnormal[3]);
  double Corner_CondNum2(double xyz0[3], double xyz1[3], double xyz2[3], double refnormal[3]);
  double Corner_CondNum2_Mod(double xyz0[3], double xyz1[3], double xyz2[3], double refnormal[3]);
  double Corner_CondNum2_Mod(double xyz0[3], double xyz1[3], double xyz2[3], double refnormal[3]);
  void CornerNormal2(double xyz0[3], double xyz1[3], double xyz2[3], double *normal);
  double Corner_CondNum3(double xyz0[3], double xyz1[3], double xyz2[3], double xyz3[3]);
  double Corner_CondNum3_Mod(double xyz0[3], double xyz1[3], double xyz2[3], double xyz3[3]);



#ifdef __cplusplus
}
#endif

#endif
