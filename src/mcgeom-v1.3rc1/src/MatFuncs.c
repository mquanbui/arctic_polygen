#include <stdlib.h>
#include <stdio.h>

#include "MCGeom.h"


void MDiff(double **a1, double **a2, int nr, int nc, double **a3) {
  int i, j;

  for (i = 0; i < nr; i++) {
    for (j = 0; j < nc; j++) {
      a3[i][j] = a1[i][j]-a2[i][j];
    }
  }
}

void MSum(double **a1, double **a2, int nr, int nc, double **a3) {
  int i, j;

  for (i = 0; i < nr; i++)
    for (j = 0; j < nc; j++) {
      a3[i][j] = a1[i][j]+a2[i][j];
    }
}

void MMult(double **a1, int nr1, int nc1, double **a2, int nr2, int nc2, double **a3) {
  int i, j, k;
  double sum;

  if (nc1 != nr2) {
#ifdef DEBUG
    fprintf(stderr,"MMult: Incompatible matrices\n");
#endif
    return;
  }

  for (i = 0; i < nr1; i++) {
    for (j = 0; j < nc2; j++) {
      a3[i][j] = 0;
      for (k = 0; k < nc1; k++)
	a3[i][j] += a1[i][k]*a2[k][j];
    }
  }  
}

double MDet(double **a, int n) {
  int i, j, k, k1, k2, k3, sgn;
  double det=0.0, detsub;
  static double **asub;
  static int maxn=0;

  /* Matrix has to be square */

  if (n == 1) {
    return a[0][0];
  }
  else if (n == 2) {
    det = a[0][0]*a[1][1]-a[1][0]*a[0][1];
    return det;
  }
  else if (n == 3) {
    det = (a[0][0]*(a[1][1]*a[2][2]-a[1][2]*a[2][1]) -
	   a[0][1]*(a[1][0]*a[2][2]-a[1][2]*a[2][0]) +
	   a[0][2]*(a[1][0]*a[2][1]-a[1][1]*a[2][0]));
    return det;
  }

  /* Expand the sub matrix if necessary */
  if (maxn < n) {
    if (maxn == 0) {
      asub = (double **) malloc(n*sizeof(double *));
      for (i = 0; i < n; i++)
	asub[i] = (double *) malloc(n*sizeof(double));
    }
    else {
      asub = (double **) realloc(asub,n*sizeof(double *));
      for (i = 0; i < maxn; i++)
	asub[i] = (double *) realloc(asub[i],n*sizeof(double *));
      for (i = maxn; i < n; i++)
	asub[i] = (double *) malloc(n*sizeof(double *));
    }
    maxn = n;
  }

  sgn = 1;
  for (j = 0; j < n; j++) {
    if (a[0][j] == 0.0)
      continue;

    /* Form sub-matrix by skipping 1st row and j'th column */
    for (k1 = 0; k1 < n-1; k1++) {
      for (k2 = 0, k3 = 0; k2 < n; k2++) {
	if (k2 == j)
	  continue;
	asub[k1][k3] = a[1+k1][k2];
	k3++;
      }
    }

    detsub = MDet(asub,n-1);

    if (detsub == 0.0)
      continue;

    det += sgn*a[0][j]*detsub;
    sgn = -sgn;
  }

  return det;
}


void MTransp(double **a, int nr, int nc, double **atrans) {
  int i, j;

  for (i = 0; i < nr; i++)
    for (j = 0; j < nc; j++)
      atrans[j][i] = a[i][j];
}

/* Code from DSP Design Performance page by Dr. Jeffrey Tafts
   http://www.nauticom.net/www/jdtaft/FortranMatrix.htm */

void MInv(double **a, int n, double **ainv) {
  double alpha, beta;
  int i, j, k, n2;
  static double **D;
  static double maxn=0;

  /* Expand the temporary matrix if necessary */
  if (maxn < n) {
    if (maxn == 0) {
      D = (double **) malloc(n*sizeof(double *));
      for (i = 0; i < n; i++)
	D[i] = (double *) malloc(2*n*sizeof(double));
    }
    else {
      D = (double **) realloc(D,n*sizeof(double *));
      for (i = 0; i < 2*maxn; i++)
	D[i] = (double *) realloc(D[i],2*n*sizeof(double *));
      for (i = 2*maxn; i < 2*n; i++)
	D[i] = (double *) malloc(2*n*sizeof(double *));
    }
    maxn = n;
  }

  /* initialize the reduction matrix */
  n2 = 2*n;
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      D[i][j] = a[i][j];
      D[i][n+j] = 0.0;
    }
    D[i][n+i] = 1.0;
  }

  /*  do the reduction  */
  for (i = 0; i < n; i++) {
    alpha = D[i][i];
    if (alpha == 0.0) {
      fprintf(stderr,"MInverse: Singular Matrix\n");
      return;
    }
    
    for (j = 0; j < n2; j++)
      D[i][j] = D[i][j]/alpha;

    for (k = 0; k < n; k++) {
      if ((k-i)== 0) 
	continue;

      beta = D[k][i];
      for (j = 0; j < n2; j++)
	D[k][j] = D[k][j] - beta*D[i][j];
    }
  }

  /* copy result into output matrix */
  for (i =0; i < n; i++)
    for (j = 0; j < n; j++)
      ainv[i][j] = D[i][j+n];

}



/* Specialized code for 3x3 or smaller matrices */

void MDiff3(double a1[3][3], double a2[3][3], int nr, int nc, double a3[3][3]) {
  int i, j;

  for (i = 0; i < nr; i++) {
    for (j = 0; j < nc; j++) {
      a3[i][j] = a1[i][j]-a2[i][j];
    }
  }
}

void MSum3(double a1[3][3], double a2[3][3], int nr, int nc, double a3[3][3]) {
  int i, j;

  for (i = 0; i < nr; i++)
    for (j = 0; j < nc; j++) {
      a3[i][j] = a1[i][j]+a2[i][j];
    }
}

void MMult3(double a1[3][3], int nr1, int nc1, double a2[3][3], int nr2, int nc2, double a3[3][3]) {
  int i, j, k;
  double sum;

  if (nc1 != nr2) {
#ifdef DEBUG
    fprintf(stderr,"MMult: Incompatible matrices\n");
#endif
    return;
  }

  for (i = 0; i < nr1; i++) {
    for (j = 0; j < nc2; j++) {
      a3[i][j] = 0;
      for (k = 0; k < nc1; k++)
	a3[i][j] += a1[i][k]*a2[k][j];
    }
  }  
}

void MVMult3(double mat[3][3], int nr1, int nc1, double vec[3], int nr2, double resvec[3]) {
  int i, k;
  double sum;

  if (nc1 != nr2) {
#ifdef DEBUG
    fprintf(stderr,"MMult: Incompatible matrices/vectors\n");
#endif
    return;
  }

  for (i = 0; i < nr1; i++) {
    resvec[i] = 0;
    for (k = 0; k < nc1; k++)
      resvec[i] += mat[i][k]*vec[k];
  }  
}

double MDet3(double a[3][3], int n) {
  int i, j, k, k1, k2, k3, sgn;
  double det=0.0, detsub;
  double asub[3][3];
  static int maxn=0;

  /* Matrix has to be square */

  if (n == 1) {
    return a[0][0];
  }
  else if (n == 2) {
    det = a[0][0]*a[1][1]-a[1][0]*a[0][1];
    return det;
  }
  else if (n == 3) {
    det = (a[0][0]*(a[1][1]*a[2][2]-a[1][2]*a[2][1]) -
	   a[0][1]*(a[1][0]*a[2][2]-a[1][2]*a[2][0]) +
	   a[0][2]*(a[1][0]*a[2][1]-a[1][1]*a[2][0]));
    return det;
  }
  else
    return 0;

}


void MTransp3(double a[3][3], int nr, int nc, double atrans[3][3]) {
  int i, j;

  for (i = 0; i < nr; i++)
    for (j = 0; j < nc; j++)
      atrans[j][i] = a[i][j];
}

/* Code from DSP Design Performance page by Dr. Jeffrey Tafts
   http://www.nauticom.net/www/jdtaft/FortranMatrix.htm */

void MInv3(double a[3][3], int n, double ainv[3][3]) {
  double alpha, beta;
  int i, j, k, n2;
  static double D[3][6];
  static double maxn=0;


  /* initialize the reduction matrix */
  n2 = 2*n;
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      D[i][j] = a[i][j];
      D[i][n+j] = 0.0;
    }
    D[i][n+i] = 1.0;
  }

  /*  do the reduction  */
  for (i = 0; i < n; i++) {
    alpha = D[i][i];
    if (alpha == 0.0) {
      fprintf(stderr,"MInverse: Singular Matrix\n");
      return;
    }
    
    for (j = 0; j < n2; j++)
      D[i][j] = D[i][j]/alpha;

    for (k = 0; k < n; k++) {
      if ((k-i)== 0) 
	continue;

      beta = D[k][i];
      for (j = 0; j < n2; j++)
	D[k][j] = D[k][j] - beta*D[i][j];
    }
  }

  /* copy result into output matrix */
  for (i =0; i < n; i++)
    for (j = 0; j < n; j++)
      ainv[i][j] = D[i][j+n];

}

