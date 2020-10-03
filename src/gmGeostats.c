/*
 * SPDX-FileCopyrightText: 2020 Helmholtz-Zentrum Dresden-Rossendorf
 *  <support@boogaart.de>
 *
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
// attention: comment this if not compiling
#include <stdio.h>
#include <Rinternals.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#define inR   // attention: this must be uncommented if not compiling

#ifdef inR 
#include <R.h>
#include <Rmath.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h> 
#endif

#define maxIntervals 1000
short int binBuf[maxIntervals];
double doubleBuf[maxIntervals];
/* int intBuf[maxIntervals];*/

/* massive use of functions to verify: min, max, abs, sqrt, cos, sin */
/* particularly needed: to build the meta-function that calls fbandXXXX 
 * functions, from line 328 */

typedef void (*vgramDensityFunctionPtr)(int d,double *,double *); 
typedef double (*vgramFunctionPtr)(double,const double *); 
typedef void (*bandSimFunctionPtr)(int m,const double *,double *,double,const double *); 

/* invBitExp2
 inverts the sequence of the bits. 
 this is used for the generation of almost equally space directions in 2D
 
 Lantuejoul (2002), page 194
 */
double invBitExp2(int i) {
  int bit = 1;
  int inv = 0;
  while(i) {
    inv <<=1;
    inv |= (i&1);
    i>>=1;
    bit<<=1;
  }
  return ((double)inv)/bit;
}


/* invBitExp2
 inverts the sequence of the digits in b-adict representation. 
 this is used for the generation of almost equally space directions in 3D
 
 Lantuejoul (2002), page 194
 */
double invBitExp(int i,int b) {
  int bit = 1;
  int inv = 0;
  while(i) {
    inv *=b;
    inv += (i%b);
    i/=b;
    bit*=b;
  }
  return ((double)inv)/bit;
}

/*
 Gaussian covariance function
 */

double cGauss(double h,const double *extra) {
  return exp(-(h*h));
}

/*
 Spherical covariance function
 */
double cSph(double h,const double *extra) {
  return( h<1 ? 1-1.5*h+0.5*(h*h*h): 0 );
}

/*
 Exponential covariance function
 */
double cExp(double h,const double *extra) {
  return exp(-h);
  
}

/*
 Switcher for covariance functions 
 */
vgramFunctionPtr cgramFunctions[]={
  cGauss,cSph,cExp
};


static R_NativePrimitiveArgType CcalcCgram_t[] = {
  /* dimX  LDX   X       dimY    LDY    Y       dimC   C     Nugget  nCgr  typeCgr  A       Sill    moreC   ijEq*/
  INTSXP,INTSXP,REALSXP,INTSXP,INTSXP,REALSXP,INTSXP,REALSXP,REALSXP,INTSXP,REALSXP,REALSXP,REALSXP,REALSXP,INTSXP
};

void CcalcCgram(
    const int *dimX, 
    const int *LDX,
    const double *X,
    const int *dimY,
    const int *LDY,
    const double *Y,
    const int *dimC,
    double *C,
    const double *Nugget, /* d x d*/
const int *nCgrams,   /* 1 */
const int *typeCgram, /* length=nCgrams */
const double *A,  /* nCgrams x m x m sqrt inverse Matrices */
const double *Sill,   /* nCgrams x d x d*/
const double *moreCgramData, /* n x ?*/
const int *ijEqual
) {
  int d=dimC[0];
  int nX=dimC[1];
  int nY=dimC[3];
  int m =dimX[1];
  int ldx = *LDX;
  int ldy = *LDY;
  if( dimC[2]!=d )
    error("CcalcVgram: Expected covariance dimensions not compatible");
  if( dimX[0]!=nX )
    error("CcalcVgram: Output does not fit input size for X");
  if( dimY[0]!=nY )
    error("CcalcVgram: Output does not fit input size for Y");
  if( dimY[1]!=m )
    error("CcalcVgram: Column dimensions of X and Y do not fit");
  if( m<1 || m>3)
    error("Can not handel spatial dimensions outside 1-3");
  int outBufSize=d*d*nX*nY;
  int i,j,k,ev,lx,ly,s;
  double delta[3];
  double v[3];
  double h2,h,val;
  for(i=0;i<outBufSize;i++)
    C[i]=0.0;
  if( *ijEqual ) {
    if( nX!=nY )
      error("CcalcVgram: ijEqual and rows of X and Y don't fit");
    for(lx=0;lx<nX;lx++) {
      for(i=0;i<d;i++)
        for(j=0;j<d;j++) 
          C[i+d*(lx+nX*(j+d*lx))]=Nugget[i+d*j];
    }
  }
  for(s=0;s<*nCgrams;s++) { // structure s
    for(lx=0;lx<nX;lx++) // sample index on dataset X
      for(ly=0;ly<nY;ly++) { // sample index on dataset Y
        for(j=0;j<m;j++) // geographic coordinate index
          delta[j]=Y[ly+ldy*j]-X[lx+ldx*j]; // compute spatial lags between the selected locations
        h2=0;
        for(j=0;j<m;j++) { 
          v[j]=0;
          for(k=0;k<m;k++) {
            v[j]+=A[s+*nCgrams * (j+m*k)]*delta[k];
          }
          h2+=v[j]*v[j];
        }
        h=sqrt(h2);
        val=(*(cgramFunctions[typeCgram[s]]))(h,moreCgramData+s);
        for(i=0;i<d;i++)
          for(j=0;j<d;j++)
            C[i+d*(lx+nX*(j+d*ly))] += val*Sill[s+*nCgrams*(i+d*j)];
      }
  }
}



/* BEGIN deprecated function ? */
/*
 vsdfGauss (vector for spatial density function )
 generates a vector of independent random normals in a vector
 */
void vsdfGauss(int d,double *extra,double *omega) {
  int i;
  for(i=0;i<d;i++)
    omega[i]=norm_rand();
}
/* END deprecated function ? */


/*
 fbandGauss
 
 generates a sinus function with random phase and random frequence on
 a band for a gaussian covariance structure.
 
 This is a mixture of turning bands and spectral simulation.
 
 */
void fbandGauss(int n, /* number of locations */
const double *projs, /* projected locations */ 
double *band, /* output */ 
double range, /* projected range */
const double *extra /* extra parameter (unused), for consistency with other covariances */
){
  /* Extract a freq. from the 1D Gaussian density along the unitdirection 
   * where projs were calculated. Evaluate a wave of random phase at the
   * projs points. Return that wave */	
  int i;
  double phase,amp,omega,d1;   
  omega = norm_rand()  * M_SQRT2 / range; 
  phase = unif_rand() * M_2PI;
  amp  = M_SQRT2; /* Lantuejoul (2002), page 191 */
  for(i=0;i<n;i++) {
    d1 = phase + projs[i]*omega; // projs, *projs==projs[0], projs[i]==*(projs+i) 
    d1 = sin(d1);
    band[i] = amp*d1;
  }
}


// IS THIS RIGHT??
void fbandSph(int n, /* number of locations */
  const double *projs, /* projected locations */ 
  double *band, /* output */ 
  double range, /* projected range */
  const double *extra /* extra parameter (unused), for consistency with other covariances */
){
  
  int i,j,nIntervals;
  double t,x0,x1,effrange;
  /* Lantuejoul (2002), page 197 */ 					
  /* Find the smallest proj; select a uniform point "x0" left from it at 
   * most "range" away. Domain = (x0, max(projs)) */
  x0 = projs[0];  /* does min exist?? */
  x1 = projs[0];
  effrange = range;
  // effrange = range *2.0; // this was an attempt to get reasonable range fitting
  for(i=1;i<n;i++) {
    if( projs[i]>x1 )
      x1=projs[i];
    else if( projs[i]<x0 )
      x0=projs[i];
  }  
  x0 += -unif_rand() * effrange;
  nIntervals = (int) ceil((x1-x0)/effrange);
  if( nIntervals > maxIntervals )
    error("fbandSph: Exceeded maxIntervals");
  for(i=0;i<nIntervals;i++) 
    binBuf[i] = unif_rand()<0.5?1:-1;
  for(i=0;i<n;i++) {
    t=(projs[i]-x0)/effrange;
    j=(unsigned int)floor(t);
    //band[i]=binBuf[j]*(t-j-0.5)*M_SQRT_3;
    band[i]=binBuf[j]*(t-j-0.5)*2.0*M_SQRT_3;
  }
}

int bsearchDouble(double x,int n,double *s) {
  int j0=0;
  int j1=n-1;
  int j;
  while(j1-j0>1) {
    /* We know s[j0] <= x < s[j1] */
    j = (j0+j1)/2;
    if( x < s[j] )
      j1=j;
    else
      j0=j;
  }
  return(j0);
}


void fbandExp(int n, /* number of locations */
    const double *projs, /* projected locations */ 
    double *band, /* output */ 
    double range, /* projected range */
    const double *extra /* extra parameter (unused), for consistency with other covariances */
){
  int i,j,ns;
  double d1,x0,x1,sign,effrange;
  /* Lantuejoul (2002), page 196 */ 
  /* Find the smallest proj; select a random exponential point "x0" left 
   * from it with lambda "2range" away. Domain = (x0, max(projs)) */
  sign = unif_rand()>0.5? 1 : -1; /* start + or - randomly */
  effrange = range;  /*  ATTENTION: 3*range, but range is inverted ??? */
  x0 = projs[0];  /* does min exist?? */
  x1 = projs[0];
  for(i=1;i<n;i++) {
    if( projs[i]>x1 )
      x1=projs[i];
    else if( projs[i]<x0 )
      x0=projs[i];
  }  
  x0 -= 2*effrange*exp_rand();
  
  /* Partition the domain with a Poisson point process of lambda=2range. */
  ns = 0;
  doubleBuf[0] = x0;
  while( doubleBuf[ns]<x1 ){
    if( ns>=maxIntervals )
      error("fbandExp: too small range; merge with nugget?");
    doubleBuf[ns+1] = doubleBuf[ns] + 2*effrange*exp_rand();
    ns ++;
  }
  /* Assign values*/
  for(i=0;i<n;i++){
    j=bsearchDouble(projs[i],ns,doubleBuf);
    d1 = (doubleBuf[j+1] + doubleBuf[j])/2; /* midpoint */
  d1 = projs[i] - d1;
  band[i] = d1>0 ? sign : -sign; /* -1 if projs[i]<midpoint; +1 otherwise */
  }  
}


bandSimFunctionPtr bandSim[]={
  fbandGauss,fbandSph,fbandExp
};

void getUnitvec(
    int dimX, /* m = 2 or 3 */
  int ip, /* number of the band being simulated */
  double *unitvec /* out: m x 1*/ 
) {
  /* weak discrepancy sequence of pseudorandom directions in 2D or 3D:
   * Lantuejoul (2002), page 194, after Freulon (1992) */			  
  int i;
  double d1,d2,d3;
  if(dimX>3)
    error("no expression for unit vectors in dimension larger than 3");
  if( dimX==3) {
    d1=invBitExp2(ip)*M_2_PI;
    d2=invBitExp(ip,3);
    d3 = sqrt(1-d2*d2); 
    unitvec[2] = d2;
    unitvec[0] = cos(d1)*d3; 
    unitvec[1] = sin(d1)*d3;
  } else if( dimX==2 ) {
    d1=invBitExp2(ip);
    unitvec[0] = cos(d1*M_PI); 
    unitvec[1] = sin(d1*M_PI);   
  } else if( dimX== 1) {
    unitvec[0]=1;
  }
}



static R_NativePrimitiveArgType CMVTurningBands_t[] = {  /* INTSXP,REALSXP */
 /* dimX,    X,  dimZ    Z     nBands sqrtNug nCgram  typeCgr   A    sqrtSill moreCgr */
 INTSXP,REALSXP,INTSXP,REALSXP,INTSXP,REALSXP,INTSXP,INTSXP,REALSXP,REALSXP,REALSXP   
 };
  

void CMVTurningBands(
    const int *dimX, /* n,m */
  const double *X,
  const int *dimZ, /*d,n,nsim */
  double *Z, /* Output Simulation transposed*/
  const int *nBands,
  const double *sqrtNugget, /* d x d*/
  const int *nCgrams,   /* 1 */
  const int *typeCgram, /* length=nCgrams */
  const double *A,  /* nCgrams x m x m sqrt inverse Matrices */
  const double *sqrtSill,   /* nCgrams x d x d*/
  const double *moreCgramData /* n x ?*/ 
) {
  const int maxCgramType=2;
  const int nsim=dimZ[2];
  int i,j,k,s,ss,ev;
  double d1,d2,d3;
  const double sqrtNBands=sqrt((double) *nBands);
  double phase,amp;
  const int n=dimX[0];
  const int m=dimX[1];
  const int d=dimZ[0];
  double projs[n];
  double band[n];
  double v[3];
  double omega[3];
  int sim;
  Rprintf("Starting calculations\n");
  if( m<1 || m>3 )
    error("CMVTurningBands: illegal X column dimension");
  if( dimZ[1]!=n )
    error("CMVTurningBands: Z and X do not fit in dimension");
#ifdef inR 
  GetRNGstate();
#endif
  /*setting Z to 0*/
  for(sim=0;sim<nsim;sim++) { 
    for(i=0;i<n;i++)
      for(j=0;j<d;j++)
        Z[d*i+j]=0.0;
    for(s=0;s<*nBands;s++){/* band */
  getUnitvec(3, s+1, &(omega[0])); /* obtain a direction; always in 3D, in order for the spherical variogram to be correct */
  //getUnitvec(m, s+1, omega); /* obtain a direction */
  for(ss=0;ss< *nCgrams;ss++) { /* variogram structure */
  if( typeCgram[ss]<0 || typeCgram[ss] > maxCgramType )
    error("CMVTurningBands: Unknown variogram type");
  /* project all data onto the direction */
  for(i=0;i<n;i++){ /* location */
  for(j=0;j<m;j++) {
    v[j]=0;
    for(k=0;k<m;k++)
      v[j]+=A[ss+ *nCgrams *(j+m*k)]*X[i+n*k];
  }
  projs[i]=0;
    for(j=0;j<m;j++){ /* spatial dimension */
  projs[i]+=omega[j]*v[j];
    }
  }
  /* for each eigenvalue, ... */
  for(ev=0;ev<d;ev++) { /* eigenvector */
  /* ... obtain a curve at all proj points following the covariance model */  
  (*bandSim[typeCgram[ss]])(n,projs,band,1.0,moreCgramData+ss);  
    /* this function takes the projs and returns on band the  the curve */
    /* ... multiply the eigenvector by the curve, and accumulate */
#ifdef _OPENMP
#pragma omp parallel for		             \
    if(!omp_in_parallel()&&0)		        \
      num_threads(omp_get_num_procs())	\
      default(shared) private(i,j,d2) 
    for(i=0;i<n;i++){ /* location */
    for(j=0;j<d;j++){ /* variable */
    d2 = sqrtSill[ss + *nCgrams *(j+d*ev)]*band[i]; 
      Z[d*i+j]+=d2;  
    }	  
    }
#else
    for(i=0;i<n;i++){ /* location */
    for(j=0;j<d;j++){ /* variable */
    d2 = sqrtSill[ss + *nCgrams *(j+d*ev)]*band[i]; 
      Z[d*i+j]+=d2;  
    }	  
    }
#endif
  }
  }
    } 
    /* Rescale*/
    for(i=0;i<n;i++)
      for(j=0;j<d;j++)
        Z[d*i+j] /= sqrtNBands; // Lantuejoul 2002 p.193
    /* Nugget */
    for(i=0;i<n;i++)
      for(j=0;j<d;j++) {
        d1 = norm_rand();
        for(k=0;k<d;k++)
          Z[d*i+k] += d1*sqrtNugget[k+d*j]; /*check*/
      }
      Z+=d*n;
  }
#ifdef inR 
  PutRNGstate();
#endif
}




static R_NativePrimitiveArgType CCondSim_t[] = {  /* INTSXP,REALSXP */
 /* dimZin  Zin   Cinv   dimX,    X,     dimZ    Z     nBands sqrtNugget nugget nCgrams typeCgr   A  sqrtSill sill   moreCgr   cbuf   dbuf */
   INTSXP,REALSXP,REALSXP,INTSXP,REALSXP,INTSXP,REALSXP,INTSXP,REALSXP,REALSXP,INTSXP,INTSXP,REALSXP,REALSXP,REALSXP,REALSXP,REALSXP,REALSXP
};
 

void CCondSim(
    const int *dimZin, /* IN: d, nin */
    const double *Zin, /*IN: nin x d Randomfield data to condition to */
    const double *Cinv,/*IN: (nin * d) x (nin x d) inverse of Covariance*/
    const int *dimX, /* IN: n, m */
    const double *X, /* IN: All Lokations, first nin conditioning */ 
    const int *dimZ, /* IN: d, n,nsim */
    double *Z, /* OUT: t() Output Simulation */
    const int *nBands, /* IN: Desired number of Bands*/
    const double *sqrtNugget, /* IN: d x d */
    const double *nugget, /* IN: dxd */
    const int *nCgrams,   /* IN: number of variograms */
    const int *typeCgram, /* IN: type of each variogram,length=nCgrams */
    /* 0=Gauss, 1=Spherical, 2=Exponential */
    const double *A,  /* IN: Anisotropy matrices, nCgrams x m x m inverse Matrices */
    const double *sqrtSill,   /* IN: nCgrams x d x d*/
    const double *sill,       /* IN: nCgrams x d x d*/
    const double *moreCgramData, /* nGrams x 1 Extraparamter*/ 
    double *cbuf,  /* BUF: Buffer of length d*d*nin */
    double *dbuf  /* BUF: Buffer of length d*nin*nsim */
) {
  const int maxCgramType=2;
  const int n=dimX[0];
  const int m=dimX[1];
  const int d=dimZin[0];
  const int nin = dimZin[1];
  const int nsim= dimZ[2];
  const int nd=n*d;
  const char No='N';
  const char Transposed='T';
  const int dmnin=d*nin;
  const int oneI=1;
  const int zeroI = 0;
  const double zero=0.0;
  const double one=1.0;
  const double minus1=-1.0;
  int i,j,k,s,ss,ev,l;
  int sim,shift;
  double d1,d2,d3,cv;
  const int dimXin[2] = {nin,m}; 
  const int dimXout[2] = {1,m};
  const int dimCbuf[4] = {d,nin,d,1};
  // Unconditional Simulation
  Rprintf("starting unconditional simulation (%d)\n",nsim);
  CMVTurningBands(dimX,
                  X,
                  dimZ,
                  Z,
                  nBands,
                  sqrtNugget,
                  nCgrams,
                  typeCgram,
                  A,
                  sqrtSill,
                  moreCgramData
  );	  
  Rprintf("unconditional simulation done (%d)\n",nsim);
  Rprintf("starting conditioning by dual kriging\n");
  /* Kriging from simulated using Cinv */
  // Create differenes of obs and sim
#ifdef _OPENMP
#pragma omp parallel			               \
  if(!omp_in_parallel())			           \
    num_threads(omp_get_num_procs())		\
    default(shared) private(i,j,sim,shift,k)
    {
#pragma omp parallel for	      
      for(int sim=0;sim<nsim;sim++) {
        int shift=nd*sim;
        for(int i=0;i<nin;i++)
          for(int j=0;j<d;j++){
            int k=d*i+j;
            Z[k+shift]-=Zin[k];
          }
      }
    }
#else
  for(int sim=0;sim<nsim;sim++) {
    int shift=nd*sim;
    for(int i=0;i<nin;i++)
      for(int j=0;j<d;j++){
        int k=d*i+j;
        Z[k+shift]-=Zin[k];
      }
  }
#endif
  /* Z[hinten] -= \hat{Z[hinten]}(Z[vorn]) = cov(hinten,vorn)%*% Cinv %*% Z[vorn] */
  /* dbuf = Cinv %*% Z[vorn] 
   Cinv in R ^ d*nin x d*nin
   Z[vorn] in R^d*nin
   
   */
  // Dual Kriging preparation
#ifdef _OPENMP
#pragma omp parallel for			           \
  if(!omp_in_parallel()&&0)			        \
    num_threads(omp_get_num_procs())		\
    default(shared) private(sim)
    for(sim=0;sim<nsim;sim++) {
      F77_NAME(dgemv)(&No,
               &dmnin,
               &dmnin,
               &one,
               Cinv,
               &dmnin,
               Z+nd*sim,
               &oneI,
               &zero,
               dbuf+dmnin*sim,
               &oneI);
    }
#else
    for(sim=0;sim<nsim;sim++) {
      F77_NAME(dgemv)(&No,
               &dmnin,
               &dmnin,
               &one,
               Cinv,
               &dmnin,
               Z+nd*sim,
               &oneI,
               &zero,
               dbuf+dmnin*sim,
               &oneI);
    }
#endif
    // /* points in*/
    //for(i=0;i<nin;i++)
    //for(j=0;j<m;j++)
    //Xin[m*i+j] = X[m*i+j];
    // /* points out*/
    //for(i=nin;i<n;i++)
    //for(j=0;j<m;j++)
    //Xout[m*(i-nin)+j] = X[m*i+j];
    /* fill cbuf*/
    // ******** Iterate over points
    for(i=nin;i<n;i++) {
      CcalcCgram(
        dimXin,
        &n,
        X,
        dimXout,
        &n,
        X+i,
        dimCbuf,
        cbuf,
        nugget, /* d x d*/
    nCgrams,   /* 1 */
    typeCgram, /* length=nCgrams */
    A,  /* nCgrams x m x m sqrt inverse Matrices */
    sill,   /* nCgrams x d x d*/
    moreCgramData,
    &zeroI
      );
      
      //for(l=0;l<nin;l++) { /* points in*/
      //for(j=0;j<d;j++)  /* koor point in */
      //for(k=0;k<d;k++) /* koor point out */
      //cbuf[j+d*(k+d*l)]=0.0; /* geeignet d*d*nin, nugget no */
      //for(s=0;s<*nCgrams;s++) { /* structure */
      //d1=0;
      //for(j=0;j<m;j++) /* spatial dimension */
      //for(k=0;k<m;k++) /* spatial dimension */
      //d1+=A[s + *nCgrams *(j+m*k)]*omega[j]*omega[k]; /* check */
      //cv=(*cgramFunctions[typeCgram[s]])(sqrt(d1),moreCgramData+s);
      //for(j=0;j<d;j++)  /* koor point in */
      //for(k=0;k<d;k++) /* koor point out */
      //cbuf[j+d*(k+d*l)]+=cv*sill[s+*nCgrams*(j+d*k)];
      //}
      //}
      
      /* Z[,i] +=  t(cbuf(d*nin,d)) %*% dbuf(d*nin) */
#ifdef _OPENMP
#pragma omp parallel for			               \
      if(!omp_in_parallel()&&0)			        \
        num_threads(omp_get_num_procs())		\
        default(shared) private(sim)
        for(sim=0;sim<nsim;sim++) {
          F77_NAME(dgemv)(&Transposed,
                   &dmnin,
                   &d,
                   &minus1,
                   cbuf,
                   &dmnin,
                   dbuf+dmnin*sim,
                   &oneI,
                   &one,
                   Z+i*d+nd*sim,
                   &oneI);
        }
#else
        for(sim=0;sim<nsim;sim++) {
          F77_NAME(dgemv)(&Transposed,
                   &dmnin,
                   &d,
                   &minus1,
                   cbuf,
                   &dmnin,
                   dbuf+dmnin*sim,
                   &oneI,
                   &one,
                   Z+i*d+nd*sim,
                   &oneI);
        }
#endif
    }
}



//anaV(v1,m,mx,i*h,dimY,y,sigma0,sigma1);
extern void anaV(double *v,   // velocity of a datum
                 const int m,         // number of variables
                 const double *x,     // location of the datum
                 const double t,      // time moment
                 const int *dimY,     // dimension of the data nodes
                 const double *y,     // data nodes
                 const double *wY,    // weights of the data nodes
                 const double sigma0, // parameter sigma0
                 const double sigma1  // parameter sigma1
) {
  size_t nY=dimY[1];  // number of data nodes
  double sigmat=sigma0+t*(sigma1-sigma0);  // deviation at this moment
  double sigmaD=(sigma1-sigma0); 
  for(int i=0;i<m;i++){ // initialize velocity
    v[i]=0;
  }
  double ws=0; // sum of weights
  for(int j=0;j<nY;j++) { // loop on number of data nodes
    // double d=0; 
    double dz[m];    
    double s2=0;
    for(int i=0;i<m;i++) {
      const double ddz=(x[i]-(1-t)*y[m*j+i])/sigmat;
      dz[i]=ddz;
      s2+=ddz*ddz;
    }
    double w=exp(-s2/2.0)*wY[j]; // weighting
    ws+=w;
    for(int i=0;i<m;i++)
      v[i]+=w*(sigmaD*dz[i]-y[m*j+i]);
  }
  for(int i=0;i<m;i++){
    v[i]/=ws;
  } 
}



static R_NativePrimitiveArgType anaForwardC_t[] = {  /* INTSXP,REALSXP */
    /* dimX,    x,     dimY    y     wY    stepsp   sigma0p    sigma1p */
     INTSXP,REALSXP,INTSXP,REALSXP,REALSXP,INTSXP, REALSXP,  REALSXP
};


extern void anaForwardC(const int *dimX,
                        double *x,
                        const int *dimY,
                        const double *y,
                        const double *wY,    // weights of the data nodes
                        const int *stepsp,
                        const double *sigma0p,
                        const double *sigma1p
) {
  
  const double sigma0=*sigma0p;
  const double sigma1=*sigma1p;
  const size_t steps=*stepsp;
  const size_t m=dimX[0];
  const size_t nx=dimX[1];
  // const size_t ny=dimY[1];
  const double h=((double)1.0)/steps;
  if( dimY[0]!=m )
    error("anaForwardC: x and y have different number of variables / rows");
#ifdef _OPENMP
  #pragma omp parallel for 
  for(size_t i=0;i<nx;i++) {
    double v1[m];
    double v2[m];
    double xx[m];
    double *mx=x+m*i;
    for(size_t s=0;s<steps;s++) {
      // v1
      anaV(v1,m,mx,s*h,dimY,y,wY,sigma0,sigma1);
      // xx=x+h*v1
      for(int j=0;j<m;j++)
        xx[j]=mx[j]+h*v1[j];
      // v2
      anaV(v2,m,xx,(s+1)*h,dimY,y,wY,sigma0,sigma1);
      // x+=0.5*h*(v1+v2)
      for(int j=0;j<m;j++)
        mx[j]+=0.5*h*(v1[j]+v2[j]);
    }
  }
#else
  for(size_t i=0;i<nx;i++) {
    double v1[m];
    double v2[m];
    double xx[m];
    double *mx=x+m*i;
    for(size_t s=0;s<steps;s++) {
      // v1
      anaV(v1,m,mx,s*h,dimY,y,wY,sigma0,sigma1);
      // xx=x+h*v1
      for(int j=0;j<m;j++)
        xx[j]=mx[j]+h*v1[j];
      // v2
      anaV(v2,m,xx,(s+1)*h,dimY,y,wY,sigma0,sigma1);
      // x+=0.5*h*(v1+v2)
      for(int j=0;j<m;j++)
        mx[j]+=0.5*h*(v1[j]+v2[j]);
    }
  }
#endif
}


static R_NativePrimitiveArgType anaBackwardC_t[] = {  /* INTSXP,REALSXP */
    /* dimX,    x,     dimY    y     wY    stepsp   sigma0p    sigma1p */
    INTSXP,REALSXP,INTSXP,REALSXP,REALSXP,INTSXP, REALSXP,  REALSXP
};

extern void anaBackwardC(const int *dimX,
                         double *x,
                         const int *dimY,
                         const double *y,
                         const double *wY,    // weights of the data nodes
                         const int *stepsp,
                         const double *sigma0p,
                         const double *sigma1p
) {
  
  const double sigma0=*sigma0p;
  const double sigma1=*sigma1p;
  const size_t steps=*stepsp;
  const size_t m=dimX[0];
  const size_t nx=dimX[1];
  // const size_t ny=dimY[1];
  const double h=((double)1.0)/steps;
  if( dimY[0]!=m )
    error("anaBackwardC: x and y have different number of variables / rows");
#ifdef _OPENMP
#pragma omp parallel for 
  for(size_t i=0;i<nx;i++) {
    double v1[m];
    double v2[m];
    double xx[m];
    double *mx=x+m*i;
    for(size_t s=0;s<steps;s++) {
      // v1
      anaV(v1,m,mx,1-s*h,dimY,y,wY,sigma0,sigma1);
      // xx=x-h*v1
      for(int j=0;j<m;j++)
        xx[j]=mx[j]-h*v1[j];
      // v2
      anaV(v2,m,xx,1-(s+1)*h,dimY,y,wY,sigma0,sigma1);
      // x+=0.5*h*(v1+v2)
      for(int j=0;j<m;j++)
        mx[j]-=0.5*h*(v1[j]+v2[j]);
    }
  }
#else
  for(size_t i=0;i<nx;i++) {
    double v1[m];
    double v2[m];
    double xx[m];
    double *mx=x+m*i;
    for(size_t s=0;s<steps;s++) {
      // v1
      anaV(v1,m,mx,1-s*h,dimY,y,wY,sigma0,sigma1);
      // xx=x-h*v1
      for(int j=0;j<m;j++)
        xx[j]=mx[j]-h*v1[j];
      // v2
      anaV(v2,m,xx,1-(s+1)*h,dimY,y,wY,sigma0,sigma1);
      // x+=0.5*h*(v1+v2)
      for(int j=0;j<m;j++)
        mx[j]-=0.5*h*(v1[j]+v2[j]);
    }
  }
#endif
}





static R_CMethodDef cMethods[] = {
  {"CcalcCgram", (DL_FUNC) &CcalcCgram, 15, CcalcCgram_t},
  {"CMVTurningBands", (DL_FUNC) &CMVTurningBands, 11, CMVTurningBands_t},
  {"CCondSim", (DL_FUNC) & CCondSim, 18,  CCondSim_t},
  {"anaForwardC", (DL_FUNC) &anaForwardC, 8, anaForwardC_t},
  {"anaBackwardC", (DL_FUNC) &anaBackwardC, 8, anaBackwardC_t},
  {NULL, NULL, 0}
};


void R_init_compositions(DllInfo *info)
{
  R_registerRoutines(info, cMethods, NULL, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
  R_forceSymbols(info, TRUE);
}

