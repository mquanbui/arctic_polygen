#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "MCGeom.h"

double func_ls(double (*pxyz)[3], int nv, double *cen);
double deriv_ls(double (*pxyz)[3], int nv, double *cen, int i);

int PF_LSSphere(double (*pxyz)[3], int nv, double *pcen, double *rad) {
  int i, j, niter, nlarge, nconst, converged, done, cgflag, k1d, cutback, nflat;
  double cen[4], ocen[4], df[4], odf[4], dfdf, odfodf=0.0, dnorm2, dnorm, B, dir[4];
  double cb0, fb0, ca, cb, cbp, cbm, fa, fb, fbp, fbm, dc, mindc, mult;
  double snorm, snorm2, eps=1.0e-06, mcen[3], vec0[3], vec1[3], normal[3], mlt;


  /* cen[0-2] is sphere center, cen[3] is the radius */

  
  /* Assume the polygon is a plane and take a large radius and distant
     center to start with */

  mcen[0] = mcen[1] = mcen[2] = 0.0; 

  for (i = 0; i < nv; i++)
    for (j = 0; j < 3; j++)
      mcen[j] += (1.0/nv)*pxyz[i][j];

  VDiff3(pxyz[2],pxyz[1],vec0);
  VDiff3(pxyz[0],pxyz[1],vec1);
  VCross3(vec0,vec1,normal);
  
  /* We can vary mlt to choose various points along a local normal
     going throught the mean of the polygon's corners */
  /* right now we just choose the mean of all the points */

  mlt = 0.0;  
  for (j = 0; j < 3; j++)
    cen[j] = mcen[j] + mlt*normal[j];


  cen[3] = sqrt((pxyz[0][0]-cen[0])*(pxyz[0][0]-cen[0]) +
               (pxyz[0][1]-cen[1])*(pxyz[0][1]-cen[1]) +
               (pxyz[0][2]-cen[2])*(pxyz[0][2]-cen[2]));

  for (i = 0; i < 4; i++)
    ocen[i] = cen[i];


  converged = 0;
  nconst = 0;
  nlarge = 0;
  niter = 1;
  while (!converged) {
    for (i = 0, dnorm2 = 0; i < 4; i++) {
      df[i] = -deriv_ls(pxyz,nv,cen,i);
      dnorm2 += df[i]*df[i];
    }
    if (dnorm2 < eps*eps*eps*eps) {
      converged = 1;
      continue;
    }

    /* Second part of search dir for NL conjugate gradient procedure
       Exceptions:
       1. For the first iter
       2. We want to restart the process after some number of iters
    */

    dfdf = dnorm2;

    cgflag = 1;
    if (niter > 1 && niter%100 != 0 && cgflag == 1) {
      B = dfdf/odfodf;
      for (i = 0; i < 4; i++) {
        df[i] += B*odf[i];
      }
    }

    for (i = 0, dnorm2 = 0; i < 4; i++)
      dnorm2 += df[i]*df[i];

    dnorm = sqrt(dnorm2);
    for (i = 0; i < 4; i++)
      dir[i] = df[i]/dnorm;

    /**********************************************************/
    /* 1-d minimization along search dir */
    /**********************************************************/

    ca = 0.0;
    dc = 10*eps;
    mindc = 0.1*eps;

    fa = func_ls(pxyz,nv,cen);

    cb0 = ca;     /* old cb */
    fb0 = fa;

    cb = ca;
    cb = ca + dc;   /* current cb */

    done = 0;  mult = 2.0;
    cutback = 0; nflat = 0;
    k1d = 0;
    while (!done) {

      for (i = 0; i < 4; i++)
        cen[i] = ocen[i] + cb*dir[i];

      fb = func_ls(pxyz,nv,cen);

      if (k1d > 100) {
        done = 1;
        continue;
      }

      /* If the function value is increasing, scale back the step size
         and try to approach the minimum with a smaller step size */

      if (fb-fb0 >= 0.0) {
        if (fb-fb0 > eps) {
          if (fabs(cb-cb0) < eps || fabs((cb-cb0)/cb0) < eps) {
            cb = cb0; /* revert back to last step and quit */
            done = 1;
          }
          else {
            /* Restart from last step and take a smaller stride */
            cutback = 1;
            dc = 0.5*dc;
            if (dc < mindc) {
              cb = cb0;
            done = 1;
            }
            else
              cb = cb0 + dc;
          }
          nflat = 0;
        }
        else { /* function is flat or increasing very slowly */
          if (nflat > 5) {
            cb = cb0;
            done = 1;
          }
          else {
            nflat++;
            cb = cb + dc;
          }
        }
      }
      else {  /* Function value has decreased */

        /* Make sure we are still on the left side of the minimum */
        cbp = cb + eps;
        for (i = 0; i < 4; i++)
          cen[i] = ocen[i] + cbp*dir[i];
        fbp = func_ls(pxyz,nv,cen);

        cbm = cb - eps; /* cb will never be 0.0 at this stage */
        for (i = 0; i < 4; i++)
          cen[i] = ocen[i] + cbm*dir[i];
        fbm = func_ls(pxyz,nv,cen);

        if (fbm > fb && fb > fbp) {
          /* We are indeed on the left side of the minimum */
          if (!cutback)
            dc = 5*dc;
          cb0 = cb;
          fb0 = fb;
          cb = cb + dc;
        }
        else {
          /* We overshot to the right side of the minimum. Cutback! */
          /* Restart from last step and take a smaller stride */
          cutback = 1;
          dc = 0.5*dc;
          if (dc < mindc) {
            cb = cb0;
            done = 1;
          }
          else
            cb = cb0 + dc;
        }
        nflat = 0;
      }
      k1d++;
    }

    snorm2 = ((ocen[0]-cen[0])*(ocen[0]-cen[0]) + 
	      (ocen[1]-cen[1])*(ocen[1]-cen[1]) +
	      (ocen[2]-cen[2])*(ocen[2]-cen[2]) + 
	      (ocen[3]-cen[3])*(ocen[3]-cen[3]));
    snorm = sqrt(snorm2);
 
    fb = func_ls(pxyz,nv,cen);

    niter++;
    if (niter > 1000) {
      /* it is EXTREMELY likely that this is some oscillation about the min */
      converged = 1;
      break;
    }

    /* Store old gradient and old conjugate gradient dir */
    for (j = 0; j < 4; j++) odf[j] = df[j];
    odfodf = dfdf;

    /* If the node movement is small then add to the count of
       how many times the movement has consecutively been small */

    if (snorm < eps)
      nconst++;


    /* Update cur coordinates to the newly calculated coordinates */
    for (j = 0; j < 4; j++) ocen[j] = cen[j];

    /* If there we previously small node movements but now the node
     movement is large again, then reset the counter. */

    if (nconst > 0 && snorm >= eps) {
      if (nlarge < 5) {
        nconst = 0;
        nlarge++;
      }
      else {
        /* solution is probably oscillating at a low value */
        converged = 1;
        continue;
      }
    }

    /* If there were five consecutive steps where the node movement
       was very small then consider the optimization to be converged */

    if (nconst >= 5) {
      converged = 1;
      continue;
    }

  }


  for (i = 0; i < 3; i++)
    pcen[i] = cen[i];
  *rad = cen[3];

  return 1;
}

double func_ls(double (*pxyz)[3], int nv, double *cen) {
  double val = 0.0, r_i, r, pcen[3] = {0.0,0.0,0.0}, ppt[3], sumwt;
  int i, j, k;
  static int maxnv = 0;
  static double *wt = NULL;

  if (nv > maxnv) {
    maxnv = nv;
    if (wt == NULL) 
      wt = (double *) malloc(maxnv*sizeof(double));
    else
      wt = (double *) realloc(wt,maxnv*sizeof(double));
  }

  r = cen[3];

  for (i = 0; i < nv; i++) {
    r_i = sqrt((pxyz[i][0]-cen[0])*(pxyz[i][0]-cen[0]) + 
               (pxyz[i][1]-cen[1])*(pxyz[i][1]-cen[1]) +
               (pxyz[i][2]-cen[2])*(pxyz[i][2]-cen[2]));
    val += (r_i - r)*(r_i - r);
  }

  
  /* return val; */


  /* Add distance of polygon "center" from the sphere */
  for (i = 0; i < nv; i++) {
    wt[i] = (1.0/nv);
    for (j = 0; j < 3; j++)
      pcen[j] += wt[i]*pxyz[i][j];
  }
  for (j = 0; j < 3; j++)
    pcen[j] /= nv;

  r_i = sqrt((pcen[0]-cen[0])*(pcen[0]-cen[0]) + 
	     (pcen[1]-cen[1])*(pcen[1]-cen[1]) +
	     (pcen[2]-cen[2])*(pcen[2]-cen[2]));
  val += (r_i - r)*(r_i - r);

  return val;
  


  /* Pick nv other random points inside the polygon and compute their 
     distance from the sphere */
  for (k = 0; k < nv; k++) {
    srand(k+1);
    for (i = 0, sumwt = 0; i < nv; i++) {
      wt[i] = rand();
      sumwt += wt[i];
    }
    for (i = 0; i < nv; i++)
      wt[i] /= sumwt;
      
    ppt[0] = ppt[1] = ppt[2] = 0.0;
    for (i = 0; i < nv; i++)
      for (j = 0; j < 3; j++)
	ppt[j] += wt[i]*pxyz[i][j];

    r_i = sqrt((ppt[0]-cen[0])*(ppt[0]-cen[0]) + 
	       (ppt[1]-cen[1])*(ppt[1]-cen[1]) +
	       (ppt[2]-cen[2])*(ppt[2]-cen[2]));
    val += (r_i - r)*(r_i - r);

  }
    

  return val;
  
}

double deriv_ls(double (*pxyz)[3], int nv, double *cen, int i) {
  double dx = 1.0e-06;
  double f0, f1, df, val=0.0;

  f0 = func_ls(pxyz,nv,cen);
  cen[i] += dx;
  f1 = func_ls(pxyz,nv,cen);
  cen[i] -= dx;

  df = f1-f0;
  val = df/dx;

  return val;
}
