/* Copyright 1998, Numerical Algorithms Group Ltd, Oxford, UK.

   FORTRAN LIBRARY
   ***** Marks 15, 16, 17 and 18. *****
   PC - Microsoft C and Microsoft Powerstation Fortran 4.0 or
        Digital Visual Fortran version.

   Header file to enable the NAG library of Fortran routines to be called from
   within a C program.

   Author: Mike Dewar, University of Bath

   Version 2.1 by Malcolm Cohen, Numerical Algorithms Group Ltd., Oxford.
   Mark 15 version compiled by Richard Hann and Ian Hounam,
   Numerical Algorithms Group Ltd., Oxford, February 1992.
   Mark 16 version compiled by Ian Hounam,
   Numerical Algorithms Group Ltd., Oxford, August 1994.
   Mark 17 version compiled by Ian Hounam,
   Numerical Algorithms Group Ltd., Oxford, February 1996.
   Mark 18 version compiled by Ian Hounam,
   Numerical Algorithms Group Ltd., Oxford, March 1998.

*/

#ifndef NAG_FTN_INCLUDED
#define NAG_FTN_INCLUDED
typedef struct { double re,im; } Complex;
#define CONST const

#ifdef __cplusplus
extern "C" {
#endif

extern void __stdcall A00AAF(
void
);

extern void __stdcall A02AAF(
  CONST double *xxr,
  CONST double *xxi,
  double *yr,
  double *yi
);

extern double __stdcall A02ABF(
  CONST double *xxr,
  CONST double *xxi
);

extern void __stdcall A02ACF(
  CONST double *xxr,
  CONST double *xxi,
  CONST double *yyr,
  CONST double *yyi,
  double *zr,
  double *zi
);

extern void __stdcall C02AEF(
  double a[],
  int *n,
  double rez[],
  double imz[],
  double *tol,
  int *ifail
);

extern void __stdcall C02AFF(
  CONST double a[] /* 2 dimension */,
  CONST int *n,
  CONST int *scale,
  double z[] /* 2 dimension */,
  double work[],
  int *ifail
);

extern void __stdcall C02AGF(
  CONST double a[],
  CONST int *n,
  CONST int *scale,
  double z[] /* 2 dimension */,
  double work[],
  int *ifail
);

extern void __stdcall C02AHF(
  CONST double *ar,
  CONST double *ai,
  CONST double *br,
  CONST double *bi,
  CONST double *cr,
  CONST double *ci,
  double zsm[],
  double zlg[],
  int *ifail
);

extern void __stdcall C02AJF(
  CONST double *a,
  CONST double *b,
  CONST double *c,
  double zsm[],
  double zlg[],
  int *ifail
);

extern void __stdcall C05ADF(
  CONST double *a,
  CONST double *b,
  CONST double *eps,
  CONST double *eta,
  double (__stdcall *f)(double *),
  double *x,
  int *ifail
);

extern void __stdcall C05AGF(
  double *x,
  CONST double *hh,
  CONST double *eps,
  CONST double *eta,
  double (__stdcall *f)(double *),
  double *a,
  double *b,
  int *ifail
);

extern void __stdcall C05AJF(
  double *x,
  CONST double *eps,
  CONST double *eta,
  double (__stdcall *f)(double *),
  CONST int *nfmax,
  int *ifail
);

extern void __stdcall C05AVF(
  double *x,
  double *fx,
  double *h,
  CONST double *boundl,
  CONST double *boundu,
  double *a,
  double c[],
  int *ind,
  int *ifail
);

extern void __stdcall C05AXF(
  double *x,
  CONST double *fx,
  CONST double *tol,
  CONST int *ir,
  CONST double *scale,
  double c[],
  int *ind,
  int *ifail
);

extern void __stdcall C05AZF(
  double *x,
  double *y,
  double *fx,
  CONST double *tolx,
  CONST int *ir,
  double c[],
  int *ind,
  int *ifail
);

extern void __stdcall C05NBF(
  void (__stdcall *fcn)(int *, double[], double[], int *),
  CONST int *n,
  double x[],
  double fvec[],
  CONST double *tol,
  double wa[],
  CONST int *lwa,
  int *ifail
);

extern void __stdcall C05NCF(
  void (__stdcall *fcn)(int *, double[], double[], int *),
  CONST int *n,
  double x[],
  double fvec[],
  CONST double *xtol,
  CONST int *maxfev,
  CONST int *ml,
  CONST int *mu,
  CONST double *epsfcn,
  double diag[],
  CONST int *mode,
  CONST double *factor,
  CONST int *nprint,
  int *nfev,
  double fjac[] /* 2 dimension */,
  CONST int *ldfjac,
  double r[],
  CONST int *lr,
  double qtf[],
  double w[] /* 2 dimension */,
  int *ifail
);

extern void __stdcall C05NDF(
  int *irevcm,
  CONST int *n,
  double x[],
  double fvec[],
  CONST double *xtol,
  CONST int *ml,
  CONST int *mu,
  CONST double *epsfcn,
  double diag[],
  CONST int *mode,
  CONST double *factor,
  double fjac[] /* 2 dimension */,
  CONST int *ldfjac,
  double r[],
  CONST int *lr,
  double qtf[],
  double w[] /* 2 dimension */,
  int *ifail
);

extern void __stdcall C05PBF(
  void (__stdcall *fcn)(int *, double[], double[], double[], int *, int *),
  CONST int *n,
  double x[],
  double fvec[],
  double fjac[] /* 2 dimension */,
  CONST int *ldfjac,
  CONST double *tol,
  double wa[],
  CONST int *lwa,
  int *ifail
);

extern void __stdcall C05PCF(
  void (__stdcall *fcn)(int *, double[], double[], double[], int *, int *),
  CONST int *n,
  double x[],
  double fvec[],
  double fjac[] /* 2 dimension */,
  CONST int *ldfjac,
  CONST double *xtol,
  CONST int *maxfev,
  double diag[],
  CONST int *mode,
  CONST double *factor,
  CONST int *nprint,
  int *nfev,
  int *njev,
  double r[],
  CONST int *lr,
  double qtf[],
  double w[] /* 2 dimension */,
  int *ifail
);

extern void __stdcall C05PDF(
  int *irevcm,
  CONST int *n,
  double x[],
  double fvec[],
  double fjac[] /* 2 dimension */,
  CONST int *ldfjac,
  CONST double *xtol,
  double diag[],
  CONST int *mode,
  CONST double *factor,
  double r[],
  CONST int *lr,
  double qtf[],
  double w[] /* 2 dimension */,
  int *ifail
);

extern void __stdcall C05ZAF(
  CONST int *m,
  CONST int *n,
  CONST double x[],
  CONST double fvec[],
  CONST double fjac[] /* 2 dimension */,
  CONST int *ldfjac,
  double xp[],
  CONST double fvecp[],
  CONST int *mode,
  double err[]
);

extern void __stdcall C06BAF(
  CONST double *seqn,
  int *ncall,
  double *result,
  double *abserr,
  double work[],
  CONST int *iwork,
  int *ifail
);

extern double __stdcall C06DBF(
  CONST double *x,
  CONST double c[],
  CONST int *n,
  CONST int *s
);

extern void __stdcall C06EAF(
  double x[],
  CONST int *pts,
  int *ifail
);

extern void __stdcall C06EBF(
  double x[],
  CONST int *pts,
  int *ifail
);

extern void __stdcall C06ECF(
  double x[],
  double y[],
  CONST int *pts,
  int *ifail
);

extern void __stdcall C06EKF(
  CONST int *job,
  double x[],
  double y[],
  CONST int *n,
  int *ifail
);

extern void __stdcall C06FAF(
  double x[],
  CONST int *pts,
  double work[],
  int *ifail
);

extern void __stdcall C06FBF(
  double x[],
  CONST int *pts,
  double work[],
  int *ifail
);

extern void __stdcall C06FCF(
  double x[],
  double y[],
  CONST int *pts,
  double work[],
  int *ifail
);

extern void __stdcall C06FFF(
  CONST int *ndim,
  CONST int *l,
  CONST int nd[],
  CONST int *n,
  double x[],
  double y[],
  double work[],
  CONST int *lwork,
  int *ifail
);

extern void __stdcall C06FJF(
  CONST int *ndim,
  CONST int nd[],
  CONST int *n,
  double x[],
  double y[],
  double work[],
  CONST int *lwork,
  int *ifail
);

extern void __stdcall C06FKF(
  CONST int *job,
  double x[],
  double y[],
  CONST int *n,
  double work[],
  int *ifail
);

extern void __stdcall C06FPF(
  CONST int *m,
  CONST int *n,
  double x[],
  CONST char *init, CONST int init_len,
  double trig[],
  double work[],
  int *ifail
);

extern void __stdcall C06FQF(
  CONST int *m,
  CONST int *n,
  double x[],
  CONST char *init, CONST int init_len,
  double trig[],
  double work[],
  int *ifail
);

extern void __stdcall C06FRF(
  CONST int *m,
  CONST int *n,
  double x[],
  double y[],
  CONST char *init, CONST int init_len,
  double trig[],
  double work[],
  int *ifail
);

extern void __stdcall C06FUF(
  CONST int *m,
  CONST int *n,
  double x[],
  double y[],
  CONST char *init, CONST int init_len,
  double trigm[],
  double trign[],
  double work[],
  int *ifail
);

extern void __stdcall C06FXF(
  CONST int *n1,
  CONST int *n2,
  CONST int *n3,
  double x[],
  double y[],
  CONST char *init, CONST int init_len,
  double trign1[],
  double trign2[],
  double trign3[],
  double work[],
  int *ifail
);

extern void __stdcall C06GBF(
  double x[],
  CONST int *pts,
  int *ifail
);

extern void __stdcall C06GCF(
  double y[],
  CONST int *pts,
  int *ifail
);

extern void __stdcall C06GQF(
  CONST int *m,
  CONST int *n,
  double x[],
  int *ifail
);

extern void __stdcall C06GSF(
  CONST int *m,
  CONST int *n,
  CONST double x[],
  double u[],
  double v[],
  int *ifail
);

extern void __stdcall C06HAF(
  CONST int *m,
  CONST int *n,
  double x[] /* 2 dimension */,
  CONST char *init, CONST int init_len,
  double trig[],
  double work[] /* 2 dimension */,
  int *ifail
);

extern void __stdcall C06HBF(
  CONST int *m,
  CONST int *n,
  double x[] /* 2 dimension */,
  CONST char *init, CONST int init_len,
  double trig[],
  double work[] /* 2 dimension */,
  int *ifail
);

extern void __stdcall C06HCF(
  CONST char *direct, CONST int direct_len,
  CONST int *m,
  CONST int *n,
  double x[] /* 2 dimension */,
  CONST char *init, CONST int init_len,
  double trig[],
  double work[] /* 2 dimension */,
  int *ifail
);

extern void __stdcall C06HDF(
  CONST char *direct, CONST int direct_len,
  CONST int *m,
  CONST int *n,
  double x[] /* 2 dimension */,
  CONST char *init, CONST int init_len,
  double trig[],
  double work[] /* 2 dimension */,
  int *ifail
);

extern void __stdcall C06LAF(
  void (__stdcall *fun)(double *, double *, double *, double *),
  CONST int *n,
  CONST double t[],
  double valinv[],
  double errest[],
  CONST double *relerr,
  CONST double *alphab,
  CONST double *tfac,
  CONST int *mxterm,
  int *nterms,
  int *na,
  double *alow,
  double *ahigh,
  int *nfeval,
  double work[],
  int *ifail
);

extern void __stdcall C06LBF(
  void (__stdcall *f)(Complex *,Complex *),
  CONST double *sigma0,
  double *sigma,
  double *b,
  CONST double *epstol,
  CONST int *mmax,
  int *m,
  double acoef[],
  double errvec[],
  int *ifail
);

extern void __stdcall C06LCF(
  CONST double *t,
  CONST double *sigma,
  CONST double *b,
  CONST int *m,
  CONST double acoef[],
  CONST double errvec[],
  double *finv,
  int *ifail
);

extern double __stdcall D01AHF(
  CONST double *a,
  CONST double *b,
  CONST double *epr,
  int *npts,
  double *relerr,
  double (__stdcall *f)(double *),
  CONST int *nl,
  int *ifail
);

extern void __stdcall D01AJF(
  double (__stdcall *f)(double *),
  CONST double *a,
  CONST double *b,
  CONST double *epsabs,
  CONST double *epsrel,
  double *result,
  double *abserr,
  double work[],
  CONST int *lwork,
  int iwork[],
  CONST int *liwork,
  int *ifail
);

extern void __stdcall D01AKF(
  double (__stdcall *f)(double *),
  CONST double *a,
  CONST double *b,
  CONST double *epsabs,
  CONST double *epsrel,
  double *result,
  double *abserr,
  double work[],
  CONST int *lwork,
  int iwork[],
  CONST int *liwork,
  int *ifail
);

extern void __stdcall D01ALF(
  double (__stdcall *f)(double *),
  CONST double *a,
  CONST double *b,
  CONST int *npts,
  CONST double points[],
  CONST double *epsabs,
  CONST double *epsrel,
  double *result,
  double *abserr,
  double work[],
  CONST int *lwork,
  int iwork[],
  CONST int *liwork,
  int *ifail
);

extern void __stdcall D01AMF(
  double (__stdcall *f)(double *),
  CONST double *bound,
  CONST int *inf,
  CONST double *epsabs,
  CONST double *epsrel,
  double *result,
  double *abserr,
  double work[],
  CONST int *lwork,
  int iwork[],
  CONST int *liwork,
  int *ifail
);

extern void __stdcall D01ANF(
  double (__stdcall *f)(double *),
  CONST double *a,
  CONST double *b,
  CONST double *omega,
  CONST int *key,
  CONST double *epsabs,
  CONST double *epsrel,
  double *result,
  double *abserr,
  double work[],
  CONST int *lwork,
  int iwork[],
  CONST int *liwork,
  int *ifail
);

extern void __stdcall D01APF(
  double (__stdcall *f)(double *),
  CONST double *a,
  CONST double *b,
  CONST double *alfa,
  CONST double *beta,
  CONST int *key,
  CONST double *epsabs,
  CONST double *epsrel,
  double *result,
  double *abserr,
  double work[],
  CONST int *lwork,
  int iwork[],
  CONST int *liwork,
  int *ifail
);

extern void __stdcall D01AQF(
  double (__stdcall *f)(double *),
  CONST double *a,
  CONST double *b,
  CONST double *c,
  CONST double *epsabs,
  CONST double *epsrel,
  double *result,
  double *abserr,
  double work[],
  CONST int *lwork,
  int iwork[],
  CONST int *liwork,
  int *ifail
);

extern void __stdcall D01ARF(
  CONST double *a,
  CONST double *b,
  double (__stdcall *f)(double *),
  CONST double *relacc,
  CONST double *absacc,
  CONST int *maxrul,
  CONST int *iparm,
  double *acc,
  double *ans,
  int *n,
  double alpha[],
  int *ifail
);

extern void __stdcall D01ASF(
  double (__stdcall *f)(double *),
  CONST double *a,
  CONST double *omega,
  CONST int *key,
  CONST double *epsabs,
  double *result,
  double *abserr,
  CONST int *limlst,
  int *lst,
  double erlst[],
  double rslst[],
  int ierlst[],
  double work[],
  CONST int *lwork,
  int iwork[],
  CONST int *liwork,
  int *ifail
);

extern void __stdcall D01ATF(
  void (__stdcall *f)(double[], double[], int *),
  CONST double *a,
  CONST double *b,
  CONST double *epsabs,
  CONST double *epsrel,
  double *result,
  double *abserr,
  double work[],
  CONST int *lwork,
  int iwork[],
  CONST int *liwork,
  int *ifail
);

extern void __stdcall D01AUF(
  void (__stdcall *f)(double[], double[], int *),
  CONST double *a,
  CONST double *b,
  CONST int *key,
  CONST double *epsabs,
  CONST double *epsrel,
  double *result,
  double *abserr,
  double work[],
  CONST int *lwork,
  int iwork[],
  CONST int *liwork,
  int *ifail
);

extern double __stdcall D01BAF(
  void (__stdcall *wtfun)(double *, double *, int *, int *, double[], double[], int *),
  CONST double *a,
  CONST double *b,
  CONST int *npts,
  double (__stdcall *fun)(double *),
  int *ifail
);

extern void __stdcall D01BAW(
  CONST double *a,
  CONST double *b,
  CONST int *itype,
  CONST int *npts,
  double weight[],
  double abscis[],
  int *ifail
);

extern void __stdcall D01BAX(
  CONST double *a,
  CONST double *b,
  CONST int *itype,
  CONST int *npts,
  double weight[],
  double abscis[],
  int *ifail
);

extern void __stdcall D01BAY(
  CONST double *a,
  CONST double *b,
  CONST int *itype,
  CONST int *npts,
  double weight[],
  double abscis[],
  int *ifail
);

extern void __stdcall D01BAZ(
  CONST double *a,
  CONST double *b,
  CONST int *itype,
  CONST int *npts,
  double weight[],
  double abscis[],
  int *ifail
);

extern void __stdcall D01BBF(
  void (__stdcall *wtfun)(double *, double *, int *, int *, double[], double[], int *),
  CONST double *a,
  CONST double *b,
  CONST int *itype,
  CONST int *npts,
  CONST double weight[],
  CONST double abscis[],
  int *ifail
);

extern void __stdcall D01BCF(
  CONST int *itype,
  CONST double *aa,
  CONST double *bb,
  CONST double *cc,
  CONST double *dd,
  CONST int *npnts,
  double weight[],
  double abscis[],
  int *ifail
);

extern void __stdcall D01BDF(
  double (__stdcall *f)(double *),
  CONST double *a,
  CONST double *b,
  CONST double *epsabs,
  CONST double *epsrel,
  double *result,
  double *abserr
);

extern void __stdcall D01DAF(
  CONST double *ya,
  CONST double *yb,
  double (__stdcall *phi1)(double *),
  double (__stdcall *phi2)(double *),
  double (__stdcall *f)(double *, double *),
  CONST double *absacc,
  double *ans,
  int *npts,
  int *ifail
);

extern void __stdcall D01EAF(
  CONST int *ndim,
  CONST double a[],
  CONST double b[],
  int *mincls,
  CONST int *maxcls,
  CONST int *nfun,
  void (__stdcall *funsub)(int *, double[], int *, double[]),
  CONST double *absreq,
  CONST double *relreq,
  CONST int *lenwrk,
  double work[],
  double finest[],
  double absest[],
  int *ifail
);

extern double __stdcall D01FBF(
  CONST int *ndim,
  CONST int nptvec[],
  CONST int *lwa,
  CONST double weight[],
  CONST double abscis[],
  double (__stdcall *fun)(int *, double[]),
  int *ifail
);

extern void __stdcall D01FCF(
  CONST int *ndim,
  CONST double a[],
  CONST double b[],
  int *minpts,
  CONST int *maxpts,
  double (__stdcall *functn)(int *, double[]),
  CONST double *eps,
  double *acc,
  CONST int *lenwrk,
  double wrkstr[],
  double *finval,
  int *ifail
);

extern void __stdcall D01FDF(
  CONST int *n,
  double (__stdcall *f)(int *, double[]),
  CONST double *sigma,
  void (__stdcall *region)(int *, double[], int *, double *, double *),
  CONST int *limit,
  CONST double *r0,
  CONST double *u,
  double *result,
  int *npts,
  int *ifail
);

extern void __stdcall D01FDV(
  CONST int *ndim,
  CONST double x[],
  CONST int *j,
  CONST double *c,
  CONST double *d
);

extern void __stdcall D01GAF(
  CONST double x[],
  CONST double y[],
  CONST int *n,
  double *ans,
  double *er,
  int *ifail
);

extern void __stdcall D01GBF(
  CONST int *numvar,
  CONST double a[],
  CONST double b[],
  int *minpts,
  CONST int *maxpts,
  double (__stdcall *functn)(int *, double[]),
  CONST double *releps,
  double *relerr,
  CONST int *lenwrk,
  double wrkstr[],
  double *finval,
  int *ifail
);

extern void __stdcall D01GCF(
  CONST int *n,
  double (__stdcall *f)(int *, double[]),
  void (__stdcall *region)(int *, double[], int *, double *, double *),
  CONST int *npts,
  double vk[],
  CONST int *nrand,
  CONST int *itrans,
  double *res,
  double *err,
  int *ifail
);

extern void __stdcall D01GDF(
  CONST int *n,
  void (__stdcall *vecfun)(int *, double[], double[], int *),
  void (__stdcall *vecreg)(int *, double[], int *, double[], double[], int *),
  CONST int *npts,
  double vk[],
  CONST int *nrand,
  CONST int *itrans,
  double *res,
  double *err,
  int *ifail
);

extern void __stdcall D01GYF(
  CONST int *n,
  CONST int *npts,
  double vk[],
  int *ifail
);

extern void __stdcall D01GZF(
  CONST int *n,
  CONST int *np1,
  CONST int *np2,
  double vk[],
  int *ifail
);

extern void __stdcall D01JAF(
  double (__stdcall *f)(int *, double[]),
  CONST int *n,
  CONST double *radius,
  CONST double *epsa,
  CONST double *epsr,
  CONST int *method,
  CONST int *icoord,
  double *result,
  double *esterr,
  int *evals,
  int *ifail
);

extern void __stdcall D01PAF(
  CONST int *numvar,
  double vertex[] /* 2 dimension */,
  CONST int *iv1,
  CONST int *iv2,
  double (__stdcall *intgnd)(int *, double[]),
  int *minord,
  CONST int *maxord,
  double intvls[],
  double *esterr,
  int *ifail
);

extern void __stdcall D02AGF(
  double *h,
  CONST double error[],
  CONST double parerr[],
  double param[],
  double c[] /* 2 dimension */,
  CONST int *n,
  CONST int *n1,
  CONST int *m1,
  void (__stdcall *aux)(double[], double[], double *, double[]),
  void (__stdcall *bcaux)(double[], double[], double[]),
  void (__stdcall *raaux)(double *, double *, double *, double[]),
  void (__stdcall *prsol)(double[], double *, int *, double[]),
  double mat[] /* 2 dimension */,
  CONST double copy[] /* 2 dimension */,
  double wspace[] /* 2 dimension */,
  double wspac1[],
  CONST double wspac2[],
  int *ifail
);

extern void __stdcall D02BAF(
  double *x,
  CONST double *xend,
  CONST int *n,
  double y[],
  double *tol,
  void (__stdcall *fcn)(double *, double[], double[]),
  double w[] /* 2 dimension */,
  int *ifail
);

extern void __stdcall D02BBF(
  double *x,
  CONST double *xend,
  CONST int *n,
  double y[],
  double *tol,
  CONST int *irelab,
  void (__stdcall *fcn)(double *, double[], double[]),
  void (__stdcall *output)(double *, double[]),
  double w[] /* 2 dimension */,
  int *ifail
);

extern void __stdcall D02BDF(
  double *x,
  CONST double *xend,
  CONST int *n,
  double y[],
  CONST double *tol,
  CONST int *irelab,
  void (__stdcall *fcn)(double *, double[], double[]),
  double *stiff,
  CONST double *ynorm,
  double w[] /* 2 dimension */,
  CONST int *iw,
  CONST int *m,
  void (__stdcall *output)(double *, double[], double[], double *),
  int *ifail
);

extern void __stdcall D02BGF(
  double *x,
  CONST double *xend,
  CONST int *n,
  double y[],
  double *tol,
  CONST double *hmax,
  CONST int *m,
  CONST double *val,
  void (__stdcall *fcn)(double *, double[], double[]),
  double w[] /* 2 dimension */,
  int *ifail
);

extern void __stdcall D02BHF(
  double *x,
  CONST double *xend,
  CONST int *n,
  double y[],
  double *tol,
  CONST int *irelab,
  CONST double *hmax,
  void (__stdcall *fcn)(double *, double[], double[]),
  double (__stdcall *g)(double *, double[]),
  double w[] /* 2 dimension */,
  int *ifail
);

extern void __stdcall D02BJF(
  double *t,
  CONST double *tend,
  CONST int *neq,
  double y[],
  void (__stdcall *fcn)(double *, double[], double[]),
  CONST double *tol,
  CONST char *relabs, CONST int relabs_len,
  void (__stdcall *output)(double *, double[]),
  double (__stdcall *g)(double *, double[]),
  double rwork[],
  int *ifail
);

extern void __stdcall D02BJX(
  CONST double *rdum1,
  CONST double rdum2[]
);

extern void __stdcall D02CAF(
  double *t,
  CONST double *tend,
  CONST int *neq,
  double y[],
  CONST double *tol,
  void (__stdcall *fcn)(double *, double[], double[]),
  double rwork[],
  int *ifail
);

extern void __stdcall D02CBF(
  double *t,
  CONST double *tend,
  CONST int *neq,
  double y[],
  CONST double *tol,
  CONST int *irelab,
  void (__stdcall *fcn)(double *, double[], double[]),
  void (__stdcall *output)(double *, double[]),
  double rwork[],
  int *ifail
);

extern void __stdcall D02CGF(
  double *t,
  CONST double *tend,
  CONST int *neq,
  double y[],
  CONST double *tol,
  CONST double *hmax,
  CONST int *m,
  CONST double *val,
  void (__stdcall *fcn)(double *, double[], double[]),
  double rwork[],
  int *ifail
);

extern void __stdcall D02CHF(
  double *t,
  CONST double *tend,
  CONST int *neq,
  double y[],
  CONST double *tol,
  CONST int *irelab,
  CONST double *hmax,
  void (__stdcall *fcn)(double *, double[], double[]),
  double (__stdcall *g)(double *, double[]),
  double rwork[],
  int *ifail
);

extern void __stdcall D02CJF(
  double *t,
  CONST double *tend,
  CONST int *neq,
  double y[],
  void (__stdcall *fcn)(double *, double[], double[]),
  CONST double *tol,
  CONST char *relabs, CONST int relabs_len,
  void (__stdcall *output)(double *, double[]),
  double (__stdcall *g)(double *, double[]),
  double rwork[],
  int *ifail
);

extern double __stdcall D02CJW(
  CONST double *rdum1,
  CONST double *rdum2
);

extern void __stdcall D02CJX(
  CONST double *rdum1,
  CONST double *rdum2
);

extern void __stdcall D02EAF(
  double *x,
  double *xend,
  CONST int *n,
  double y[],
  double *tol,
  void (__stdcall *fcn)(double[], double[], double[]),
  double w[],
  CONST int *iw,
  int *ifail
);

extern void __stdcall D02EAZ(
  CONST double *x,
  CONST double y[],
  CONST double pw[] /* 2 dimension */
);

extern void __stdcall D02EBF(
  double *x,
  double *xend,
  CONST int *n,
  double y[],
  double *tol,
  CONST int *irelab,
  void (__stdcall *fcn)(double[], double[], double[]),
  CONST int *mped,
  void (__stdcall *pederv)(double[], double[], double[]),
  void (__stdcall *output)(double *, double[]),
  double w[],
  CONST int *iw,
  int *ifail
);

extern void __stdcall D02EGF(
  double *x,
  double *xend,
  CONST int *n,
  double y[],
  double *tol,
  double *hmax,
  CONST int *m,
  CONST double *val,
  void (__stdcall *fcn)(double[], double[], double[]),
  double w[],
  CONST int *iw,
  int *ifail
);

extern void __stdcall D02EHF(
  double *x,
  double *xend,
  CONST int *n,
  double y[],
  double *tol,
  CONST int *irelab,
  double *hmax,
  void (__stdcall *fcn)(double[], double[], double[]),
  CONST int *mped,
  void (__stdcall *pederv)(double[], double[], double[]),
  double (__stdcall *g)(double *, double[]),
  double w[],
  CONST int *iw,
  int *ifail
);

extern void __stdcall D02EJF(
  double *x,
  double *xend,
  CONST int *n,
  double y[],
  void (__stdcall *fcn)(double[], double[], double[]),
  void (__stdcall *pederv)(double *, double[], double[]),
  double *tol,
  CONST char *relabs, CONST int relabs_len,
  void (__stdcall *output)(double *, double[]),
  double (__stdcall *g)(double *, double[]),
  double w[],
  CONST int *iw,
  int *ifail
);

extern double __stdcall D02EJW(
  CONST double *rdum1,
  CONST double *rdum2
);

extern void __stdcall D02EJX(
  CONST double *rdum1,
  CONST double *rdum2
);

extern void __stdcall D02EJY(
  CONST double *rdum1,
  CONST double rdum2[],
  CONST double rdum3[]
);

extern void __stdcall D02GAF(
  CONST double u[] /* 2 dimension */,
  CONST double v[] /* 2 dimension */,
  CONST int *n,
  CONST double *a,
  CONST double *b,
  CONST double *tol,
  void (__stdcall *fcn)(double *, double[], double[]),
  CONST int *mnp,
  double x[],
  double y[] /* 2 dimension */,
  int *np,
  double w[],
  CONST int *lw,
  int iw[],
  CONST int *liw,
  int *ifail
);

extern void __stdcall D02GAX(
  CONST double *eps,
  CONST double y[],
  CONST double z[],
  CONST double a[],
  CONST int *m
);

extern void __stdcall D02GAZ(
  CONST double *x,
  CONST double *eps,
  CONST double y[],
  CONST double f[],
  CONST int *m
);

extern void __stdcall D02GBF(
  CONST double *a,
  CONST double *b,
  CONST int *n,
  CONST double *tol,
  void (__stdcall *fcnf)(double *, double[]),
  void (__stdcall *fcng)(double *, double[]),
  double c[] /* 2 dimension */,
  double d[] /* 2 dimension */,
  double gam[],
  CONST int *mnp,
  double x[],
  double y[] /* 2 dimension */,
  int *np,
  double w[],
  CONST int *lw,
  int iw[],
  CONST int *liw,
  int *ifail
);

extern void __stdcall D02HAF(
  double a[] /* 2 dimension */,
  CONST double b[] /* 2 dimension */,
  CONST int *n,
  CONST double *x,
  CONST double *x1,
  CONST double *tol,
  void (__stdcall *fcn)(double *, double[], double[]),
  double soln[] /* 2 dimension */,
  CONST int *m1,
  double w[] /* 2 dimension */,
  CONST int *iw,
  int *ifail
);

extern void __stdcall D02HBF(
  double p[],
  CONST int *n1,
  CONST double pe[],
  CONST double e[],
  CONST int *n,
  double soln[] /* 2 dimension */,
  CONST int *m1,
  void (__stdcall *fcn)(double *, double[], double[], double[]),
  void (__stdcall *bc)(double[], double[], double[]),
  void (__stdcall *range)(double[], double[], double[]),
  double w[] /* 2 dimension */,
  CONST int *iw,
  int *ifail
);

extern void __stdcall D02HBX(
  CONST int *istate,
  CONST int *iflag,
  CONST int *ifail1,
  CONST double p[],
  CONST int *m,
  CONST double f[],
  CONST double *pnorm,
  CONST double *pnorm1,
  CONST double *eps,
  CONST double d[]
);

extern int __stdcall D02HBY(
  CONST double p[],
  CONST int *m
);

extern void __stdcall D02HBZ(
  CONST double e[],
  CONST int *q,
  CONST double p[],
  CONST int *m
);

extern void __stdcall D02JAF(
  CONST int *n,
  double (__stdcall *cf)(int *, double *),
  void (__stdcall *bc)(int *, int *, double *),
  CONST double *x0,
  CONST double *x1,
  CONST int *k1,
  CONST int *kp,
  double c[],
  double w[],
  CONST int *lw,
  int iw[],
  int *ifail
);

extern void __stdcall D02JBF(
  CONST int *n,
  double (__stdcall *cf)(int *, int *, double *),
  void (__stdcall *bc)(int *, int *, double *),
  CONST double *x0,
  CONST double *x1,
  CONST int *k1,
  CONST int *kp,
  double c[] /* 2 dimension */,
  CONST int *ic,
  double w[],
  CONST int *lw,
  int iw[],
  CONST int *liw,
  int *ifail
);

extern void __stdcall D02KAF(
  CONST double *xl,
  CONST double *xr,
  void (__stdcall *coeffn)(double *, double *, double *, double *, double *, int *),
  double bcond[] /* 2 dimension */,
  CONST int *k,
  CONST double *tol,
  double *elam,
  double *delam,
  void (__stdcall *monit)(int *, int *, double *, double[]),
  int *ifail
);

extern void __stdcall D02KAY(
  CONST int *nit,
  CONST int *iflag,
  CONST double *elam,
  CONST double finfo[]
);

extern void __stdcall D02KDF(
  CONST double xpoint[],
  CONST int *nxp,
  void (__stdcall *coeffn)(double *, double *, double *, double *, double *, int *),
  void (__stdcall *bdyval)(double *, double *, double *, double[], double[]),
  CONST int *k,
  CONST double *tol,
  double *elam,
  double *delam,
  double hmax[] /* 2 dimension */,
  int *maxit,
  CONST int *maxfun,
  void (__stdcall *monit)(int *, int *, double *, double[]),
  int *ifail
);

extern void __stdcall D02KEF(
  CONST double xpoint[],
  CONST int *nxp,
  int *ic1,
  void (__stdcall *coeffn)(double *, double *, double *, double *, double *, int *),
  void (__stdcall *bdyval)(double *, double *, double *, double[], double[]),
  CONST int *k,
  CONST double *tol,
  double *elam,
  double *delam,
  double hmax[] /* 2 dimension */,
  int *maxit,
  CONST int *maxfun,
  void (__stdcall *monit)(int *, int *, double *, double[]),
  void (__stdcall *report)(double *, double[], int *),
  int *ifail
);

extern void __stdcall D02LAF(
  void (__stdcall *f)(int *, double *, double[], double[]),
  CONST int *neq,
  double *t,
  CONST double *tend,
  double y[],
  double yp[],
  double ydp[],
  double rwork[],
  CONST int *lrwork,
  int *ifail
);

extern void __stdcall D02LXF(
  CONST int *neq,
  CONST double *h,
  CONST double *tol,
  CONST double thres[],
  CONST double thresp[],
  CONST int *maxstp,
  int *start,
  CONST int *onestp,
  CONST int *high,
  double rwork[],
  CONST int *lrwork,
  int *ifail
);

extern void __stdcall D02LYF(
  CONST int *neq,
  double *hnext,
  double *hused,
  double *hstart,
  int *nsucc,
  int *nfail,
  int *natt,
  double thres[],
  double thresp[],
  CONST double rwork[],
  CONST int *lrwork,
  int *ifail
);

extern void __stdcall D02LZF(
  CONST int *neq,
  CONST double *t,
  CONST double y[],
  CONST double yp[],
  CONST int *nwant,
  CONST double *twant,
  double ywant[],
  double ypwant[],
  CONST double rwork[],
  CONST int *lrwork,
  int *ifail
);

extern void __stdcall D02MVF(
  CONST int *neqmax,
  CONST int *ny2dim,
  CONST int *maxord,
  double _const[],
  CONST double *tcrit,
  CONST double *hmin,
  CONST double *hmax,
  CONST double *h0,
  CONST int *maxstp,
  CONST int *mxhnil,
  CONST char *norm, CONST int norm_len,
  double rwork[],
  int *ifail
);

extern void __stdcall D02MZF(
  CONST double *tsol,
  double sol[],
  CONST int *m,
  CONST int *neqmax,
  CONST int *neq,
  CONST double ysave[] /* 2 dimension */,
  CONST int *ny2dim,
  CONST double rwork[],
  int *ifail
);

extern void __stdcall D02NBF(
  CONST int *neq,
  CONST int *neqmax,
  double *t,
  double *tout,
  double y[],
  double ydoti[],
  double rwork[],
  CONST double rtol[],
  CONST double atol[],
  CONST int *itol,
  int inform[],
  void (__stdcall *fcn)(int *, double[], double[], double[], int *),
  double ysave[] /* 2 dimension */,
  CONST int *ny2dim,
  void (__stdcall *jac)(int *, double[], double[], double[], double[], double[]),
  double wkjac[],
  CONST int *nwkjac,
  void (__stdcall *monitr)(int *, int *, double[], double[], double[], double[], 
                 double[], double[], double[], double[], int *, int *, 
                 double[], double[], int *),
  CONST int *itask,
  CONST int *itrace,
  int *ifail
);

extern void __stdcall D02NBY(
  CONST int *neq,
  CONST int *neqmax,
  CONST double *t,
  CONST double *hlast,
  CONST double *hnext,
  CONST double y[],
  CONST double ydot[],
  CONST double ysave[] /* 2 dimension */,
  CONST double r[],
  CONST double acor[] /* 2 dimension */,
  CONST int *imon,
  CONST int *inln,
  CONST double *hmin,
  CONST double *hmax,
  CONST int *nqu
);

extern void __stdcall D02NBZ(
  CONST int *neq,
  CONST double *t,
  CONST double y[],
  CONST double *h,
  CONST double *d,
  CONST double p[]
);

extern void __stdcall D02NCF(
  CONST int *neq,
  CONST int *neqmax,
  double *t,
  double *tout,
  double y[],
  double ydoti[],
  double rwork[],
  CONST double rtol[],
  CONST double atol[],
  CONST int *itol,
  int inform[],
  void (__stdcall *fcn)(int *, double[], double[], double[], int *),
  double ysave[] /* 2 dimension */,
  CONST int *ny2dim,
  void (__stdcall *jac)(int *, double[], double[], double[], double[], int *, int *, 
              double[]),
  double wkjac[],
  CONST int *nwkjac,
  int jacpvt[],
  CONST int *njcpvt,
  void (__stdcall *monitr)(int *, int *, double[], double[], double[], double[], 
                 double[], double[], double[], double[], int *, int *, 
                 double[], double[], int *),
  CONST int *itask,
  CONST int *itrace,
  int *ifail
);

extern void __stdcall D02NCZ(
  CONST int *neq,
  CONST double *t,
  CONST double y[],
  CONST double *h,
  CONST double *d,
  CONST int *ml,
  CONST int *mu,
  CONST double p[]
);

extern void __stdcall D02NDF(
  CONST int *neq,
  CONST int *neqmax,
  double *t,
  double *tout,
  double y[],
  double ydoti[],
  double rwork[],
  CONST double rtol[],
  CONST double atol[],
  CONST int *itol,
  int inform[],
  void (__stdcall *fcn)(int *, double[], double[], double[], int *),
  double ysave[] /* 2 dimension */,
  CONST int *ny2dim,
  void (__stdcall *jac)(int *, double[], double[], double[], double[], int *, double[]),
  double wkjac[],
  CONST int *nwkjac,
  int jacpvt[],
  CONST int *njcpvt,
  void (__stdcall *monitr)(int *, int *, double[], double[], double[], double[], 
                 double[], double[], double[], double[], int *, int *, 
                 double[], double[], int *),
  CONST int *itask,
  CONST int *itrace,
  int *ifail
);

extern void __stdcall D02NDZ(
  CONST int *neq,
  CONST double *t,
  CONST double y[],
  CONST double *h,
  CONST double *d,
  CONST int *j,
  CONST double p[]
);

extern void __stdcall D02NGF(
  CONST int *neq,
  CONST int *neqmax,
  double *t,
  double *tout,
  double y[],
  double ydoti[],
  double rwork[],
  CONST double rtol[],
  CONST double atol[],
  CONST int *itol,
  int inform[],
  void (__stdcall *resid)(int *, double[], double[], double[], double[], int *),
  double ysave[] /* 2 dimension */,
  CONST int *ny2dim,
  void (__stdcall *jac)(int *, double[], double[], double[], double[], double[], 
              double[]),
  double wkjac[],
  CONST int *nwkjac,
  void (__stdcall *monitr)(int *, int *, double[], double[], double[], double[], 
                 double[], double[], double[], double[], int *, int *, 
                 double[], double[], int *),
  int lderiv[],
  CONST int *itask,
  CONST int *itrace,
  int *ifail
);

extern void __stdcall D02NGZ(
  CONST int *neq,
  CONST double *t,
  CONST double y[],
  CONST double ydot[],
  CONST double *h,
  CONST double *d,
  CONST double p[]
);

extern void __stdcall D02NHF(
  CONST int *neq,
  CONST int *neqmax,
  double *t,
  double *tout,
  double y[],
  double ydoti[],
  double rwork[],
  CONST double rtol[],
  CONST double atol[],
  CONST int *itol,
  int inform[],
  void (__stdcall *resid)(int *, double[], double[], double[], double[], int *),
  double ysave[] /* 2 dimension */,
  CONST int *ny2dim,
  void (__stdcall *jac)(int *, double[], double[], double[], double[], double[], int *, 
              int *, double[]),
  double wkjac[],
  CONST int *nwkjac,
  int jacpvt[],
  CONST int *njcpvt,
  void (__stdcall *monitr)(int *, int *, double[], double[], double[], double[], 
                 double[], double[], double[], double[], int *, int *, 
                 double[], double[], int *),
  int lderiv[],
  CONST int *itask,
  CONST int *itrace,
  int *ifail
);

extern void __stdcall D02NHZ(
  CONST int *neq,
  CONST double *t,
  CONST double y[],
  CONST double ydot[],
  CONST double *h,
  CONST double *d,
  CONST int *ml,
  CONST int *mu,
  CONST double p[]
);

extern void __stdcall D02NJF(
  CONST int *neq,
  CONST int *neqmax,
  double *t,
  double *tout,
  double y[],
  double ydoti[],
  double rwork[],
  CONST double rtol[],
  CONST double atol[],
  CONST int *itol,
  int inform[],
  void (__stdcall *resid)(int *, double[], double[], double[], double[], int *),
  double ysave[] /* 2 dimension */,
  CONST int *ny2dim,
  void (__stdcall *jac)(int *, double[], double[], double[], double[], double[], int *, 
              double[]),
  double wkjac[],
  CONST int *nwkjac,
  int jacpvt[],
  CONST int *njcpvt,
  void (__stdcall *monitr)(int *, int *, double[], double[], double[], double[], 
                 double[], double[], double[], double[], int *, int *, 
                 double[], double[], int *),
  int lderiv[],
  CONST int *itask,
  CONST int *itrace,
  int *ifail
);

extern void __stdcall D02NJZ(
  CONST int *neq,
  CONST double *t,
  CONST double y[],
  CONST double ydot[],
  CONST double *h,
  CONST double *d,
  CONST int *j,
  CONST double p[]
);

extern void __stdcall D02NMF(
  CONST int *neq,
  CONST int *neqmax,
  double *t,
  double *tout,
  double y[],
  double ydoti[],
  double rwork[],
  CONST double rtol[],
  CONST double atol[],
  CONST int *itol,
  int inform[],
  double ysave[] /* 2 dimension */,
  CONST int *ny2dim,
  double wkjac[],
  CONST int *nwkjac,
  int jacpvt[],
  CONST int *njcpvt,
  int *imon,
  int *inln,
  int *ires,
  int *irevcm,
  CONST int *itask,
  CONST int *jtrace,
  int *ifail
);

extern void __stdcall D02NNF(
  CONST int *neq,
  CONST int *neqmax,
  double *t,
  double *tout,
  double y[],
  double ydoti[],
  double rwork[],
  CONST double rtol[],
  CONST double atol[],
  CONST int *itol,
  int inform[],
  double ysave[] /* 2 dimension */,
  CONST int *ny2dim,
  double wkjac[],
  CONST int *nwkjac,
  int jacpvt[],
  CONST int *njcpvt,
  int *imon,
  int *inln,
  int *ires,
  int *irevcm,
  int lderiv[],
  CONST int *itask,
  CONST int *jtrace,
  int *ifail
);

extern void __stdcall D02NRF(
  int *j,
  int *iplace,
  CONST int inform[]
);

extern void __stdcall D02NSF(
  CONST int *neq,
  CONST int *neqmax,
  CONST char *jceval, CONST int jceval_len,
  CONST int *nwkjac,
  double rwork[],
  int *ifail
);

extern void __stdcall D02NTF(
  CONST int *neq,
  CONST int *neqmax,
  CONST char *jceval, CONST int jceval_len,
  CONST int *ml,
  CONST int *mu,
  CONST int *nwkjac,
  CONST int *njcpvt,
  double rwork[],
  int *ifail
);

extern void __stdcall D02NUF(
  CONST int *neq,
  CONST int *neqmax,
  CONST char *jceval, CONST int jceval_len,
  CONST int *nwkjac,
  CONST int ia[],
  CONST int *nia,
  CONST int ja[],
  CONST int *nja,
  int jacpvt[],
  CONST int *njcpvt,
  CONST double *sens,
  CONST double *u,
  CONST double *eta,
  CONST int *lblock,
  CONST int *isplit,
  double rwork[],
  int *ifail
);

extern void __stdcall D02NVF(
  CONST int *neqmax,
  CONST int *ny2dim,
  CONST int *maxord,
  CONST char *method, CONST int method_len,
  CONST int *petzld,
  double _const[],
  CONST double *tcrit,
  CONST double *hmin,
  CONST double *hmax,
  CONST double *h0,
  CONST int *maxstp,
  CONST int *mxhnil,
  CONST char *norm, CONST int norm_len,
  double rwork[],
  int *ifail
);

extern void __stdcall D02NWF(
  CONST int *neqmax,
  CONST int *ny2dim,
  CONST int *maxord,
  double _const[],
  CONST double *tcrit,
  CONST double *hmin,
  CONST double *hmax,
  CONST double *h0,
  CONST int *maxstp,
  CONST int *mxhnil,
  CONST char *norm, CONST int norm_len,
  double rwork[],
  int *ifail
);

extern void __stdcall D02NXF(
  CONST int *icall,
  int *liwreq,
  int *liwusd,
  int *lrwreq,
  int *lrwusd,
  int *nlu,
  int *nnz,
  int *ngp,
  int *isplit,
  int *igrow,
  CONST int *lblock,
  int *nblock,
  CONST int inform[]
);

extern void __stdcall D02NYF(
  CONST int *neq,
  CONST int *neqmax,
  double *hu,
  double *h,
  double *tcur,
  double *tolsf,
  CONST double rwork[],
  int *nst,
  int *nre,
  int *nje,
  int *nqu,
  int *nq,
  int *niter,
  int *imxer,
  int algequ[],
  CONST int inform[],
  int *ifail
);

extern void __stdcall D02NZF(
  CONST int *neqmax,
  CONST double *tcrit,
  CONST double *h,
  CONST double *hmin,
  CONST double *hmax,
  CONST int *maxstp,
  CONST int *maxhnl,
  double rwork[],
  int *ifail
);

extern void __stdcall D02PAF(
  double *x,
  CONST double *xend,
  CONST int *n,
  double y[],
  double cin[],
  CONST double *tol,
  void (__stdcall *fcn)(double *, double[], double[]),
  double comm[],
  double _const[],
  double _cout[],
  double w[] /* 2 dimension */,
  CONST int *iw,
  CONST int *iw1,
  int *ifail
);

extern void __stdcall D02PCF(
  void (__stdcall *f)(double *, double[], double[]),
  CONST double *twant,
  double *tgot,
  double ygot[],
  double ypgot[],
  double ymax[],
  double work[],
  int *ifail
);

extern void __stdcall D02PDF(
  void (__stdcall *f)(double *, double[], double[]),
  double *tnow,
  double ynow[],
  double ypnow[],
  double work[],
  int *ifail
);

extern void __stdcall D02PVF(
  CONST int *neq,
  CONST double *tstart,
  CONST double ystart[],
  CONST double *tend,
  CONST double *tol,
  CONST double thres[],
  CONST int *method,
  CONST char *task, CONST int task_len,
  CONST int *errass,
  CONST double *hstart,
  double work[],
  CONST int *lenwrk,
  int *ifail
);

extern void __stdcall D02PWF(
  CONST double *tendnu,
  int *ifail
);

extern void __stdcall D02PXF(
  CONST double *twant,
  CONST char *reqest, CONST int reqest_len,
  CONST int *nwant,
  double ywant[],
  double ypwant[],
  void (__stdcall *f)(double *, double[], double[]),
  double work[],
  double wrkint[],
  CONST int *lenint,
  int *ifail
);

extern void __stdcall D02PYF(
  int *totfcn,
  int *stpcst,
  double *waste,
  int *stpsok,
  double *hnext,
  int *ifail
);

extern void __stdcall D02PZF(
  double rmserr[],
  double *errmax,
  double *terrmx,
  CONST double work[],
  int *ifail
);

extern void __stdcall D02QDF(
  double *x,
  double *xend,
  CONST int *n,
  double y[],
  double cin[],
  CONST double rtol[],
  CONST double atol[],
  void (__stdcall *fcn)(double[], double[], double[]),
  double comm[],
  double _const[],
  double _cout[],
  CONST char *jacstr, CONST int jacstr_len,
  CONST int mbands[],
  void (__stdcall *pederv)(double *, double[], double[]),
  double w[],
  CONST int *iw,
  int iwk[],
  int *ifail
);

extern void __stdcall D02QFF(
  void (__stdcall *f)(int *, double *, double[], double[]),
  CONST int *neqf,
  double *t,
  double y[],
  CONST double *tout,
  double (__stdcall *g)(int *, double *, double[], double[], int *),
  CONST int *neqg,
  int *root,
  double rwork[],
  CONST int *lrwork,
  int iwork[],
  CONST int *liwork,
  int *ifail
);

extern double __stdcall D02QFZ(
  CONST int *neqf,
  CONST double *t,
  CONST double y[],
  CONST double yp[],
  CONST int *k
);

extern void __stdcall D02QGF(
  CONST int *neqf,
  double *t,
  double y[],
  CONST double *tout,
  CONST int *neqg,
  int *root,
  int *irevcm,
  double *trvcm,
  int *yrvcm,
  int *yprvcm,
  CONST double *grvcm,
  int *kgrvcm,
  double rwork[],
  CONST int *lrwork,
  int iwork[],
  CONST int *liwork,
  int *ifail
);

extern void __stdcall D02QQF(
  CONST double comm[],
  CONST double chk[],
  CONST int *n,
  double w[],
  CONST int *iw,
  int *ifail
);

extern void __stdcall D02QWF(
  char *statef, CONST int statef_len,
  CONST int *neqf,
  CONST int *vectol,
  CONST double atol[],
  CONST int *latol,
  CONST double rtol[],
  CONST int *lrtol,
  CONST int *onestp,
  CONST int *crit,
  CONST double *tcrit,
  CONST double *hmax,
  CONST int *maxstp,
  CONST int *neqg,
  int *alterg,
  CONST int *sophst,
  double rwork[],
  CONST int *lrwork,
  int iwork[],
  CONST int *liwork,
  int *ifail
);

extern void __stdcall D02QXF(
  CONST int *neqf,
  double yp[],
  double *tcurr,
  double *hlast,
  double *hnext,
  int *odlast,
  int *odnext,
  int *nsucc,
  int *nfail,
  double *tolfac,
  int *badcmp,
  CONST double rwork[],
  CONST int *lrwork,
  CONST int iwork[],
  CONST int *liwork,
  int *ifail
);

extern void __stdcall D02QYF(
  CONST int *neqg,
  int *index,
  int *type,
  int events[],
  double resids[],
  CONST double rwork[],
  CONST int *lrwork,
  CONST int iwork[],
  CONST int *liwork,
  int *ifail
);

extern void __stdcall D02QZF(
  CONST int *neqf,
  CONST double *twant,
  CONST int *nwant,
  double ywant[],
  double ypwant[],
  CONST double rwork[],
  CONST int *lrwork,
  CONST int iwork[],
  CONST int *liwork,
  int *ifail
);

extern void __stdcall D02RAF(
  CONST int *m,
  CONST int *nmax,
  int *n,
  CONST int *numbeg,
  CONST int *nummix,
  CONST double *tol,
  CONST int *init,
  double x[],
  double y[] /* 2 dimension */,
  CONST int *iy,
  double abt[],
  void (__stdcall *fcn)(double *, double *, double[], double[], int *),
  void (__stdcall *g)(double *, double[], double[], double[], int *),
  CONST int *ijac,
  void (__stdcall *jacobf)(double *, double *, double[], double[], int *),
  void (__stdcall *jacobg)(double *, double[], double[], double[], double[], int *),
  double *deleps,
  void (__stdcall *jaceps)(double *, double *, double[], double[], int *),
  void (__stdcall *jacgep)(double *, double[], double[], double[], int *),
  double work[],
  CONST int *lwork,
  int iwork[],
  CONST int *liwork,
  int *ifail
);

extern void __stdcall D02SAF(
  double p[],
  CONST int *m,
  CONST int *n,
  CONST int *n1,
  CONST double pe[],
  double pf[],
  CONST double e[],
  double dp[],
  int *npoint,
  double wp[] /* 2 dimension */,
  CONST int *iwp,
  CONST int *icount,
  void (__stdcall *range)(double[], int *, double[], int *),
  void (__stdcall *bc)(double[], double[], double[], int *, int *),
  void (__stdcall *fcn)(double *, double[], double[], int *, double[], int *, int *),
  void (__stdcall *eqn)(double[], int *, double[], int *),
  int (__stdcall *constr)(double[], int *),
  CONST double *ymax,
  void (__stdcall *monit)(int *, int *, int *, double[], int *, double[], double *, 
                double *, double *, double[]),
  void (__stdcall *prsol)(double *, double[], int *),
  double w[] /* 2 dimension */,
  CONST int *iw1,
  CONST int *iw2,
  int *ifail
);

extern void __stdcall D02SAS(
  CONST int *istate,
  CONST int *iflag,
  CONST int *ifail1,
  CONST double p[],
  CONST int *m,
  CONST double f[],
  CONST double *pnorm,
  CONST double *pnorm1,
  CONST double *eps,
  CONST double d[]
);

extern void __stdcall D02TGF(
  CONST int *n,
  CONST int m[],
  CONST int l[],
  CONST double *x0,
  CONST double *x1,
  CONST int *k1,
  CONST int *kp,
  double c[] /* 2 dimension */,
  CONST int *ic,
  void (__stdcall *coeff)(double *, int *, double[], int *, int *, double *),
  void (__stdcall *bdyc)(double *, int *, int *, double[], int *, int *, double *),
  double w[],
  CONST int *lw,
  int iw[],
  CONST int *liw,
  int *ifail
);

extern void __stdcall D02TKF(
  void (__stdcall *ffun)(),
  void (__stdcall *fjac)(),
  void (__stdcall *gafun)(),
  void (__stdcall *gbfun)(),
  void (__stdcall *gajac)(),
  void (__stdcall *gbjac)(),
  void (__stdcall *guess)(),
  double work[],
  int iwork[],
  int *ifail
);

extern void __stdcall D02TVF(
  CONST int *neq,
  CONST int m[],
  CONST int *nlbc,
  CONST int *nrbc,
  CONST int *ncol,
  CONST double tols[],
  CONST int *mxmesh,
  CONST int *nmesh,
  CONST double mesh[],
  CONST int ipmesh[],
  double rwork[],
  CONST int *lrwork,
  int iwork[],
  CONST int *liwork,
  int *ifail
);

extern void __stdcall D02TXF(
  CONST int *mxmesh,
  CONST int *nmesh,
  CONST double mesh[],
  CONST int ipmesh[],
  double rwork[],
  int iwork[],
  int *ifail
);

extern void __stdcall D02TYF(
  CONST double *x,
  double y[] /* 2 dimension */,
  CONST int *neq,
  CONST int *mmax,
  double rwork[],
  CONST int iwork[],
  int *ifail
);

extern void __stdcall D02TZF(
  CONST int *mxmesh,
  int *nmesh,
  double mesh[],
  int ipmesh[],
  double *ermx,
  int *iermx,
  int *ijermx,
  CONST double rwork[],
  CONST int iwork[],
  int *ifail
);

extern void __stdcall D02XAF(
  CONST double *xsol,
  CONST double *x,
  CONST double _cout[],
  CONST int *n,
  CONST double y[],
  CONST double w[] /* 2 dimension */,
  CONST int *iw,
  double sol[],
  int *ifail
);

extern void __stdcall D02XBF(
  CONST double *xsol,
  CONST double *x,
  CONST double _cout[],
  CONST int *n,
  CONST double y[],
  CONST double w[] /* 2 dimension */,
  CONST int *iw,
  CONST int *m,
  double *sol,
  int *ifail
);

extern void __stdcall D02XJF(
  CONST double *xsol,
  double sol[],
  CONST int *m,
  CONST double w[] /* 2 dimension */,
  CONST int *neqmax,
  CONST int *iw,
  CONST int *neq,
  CONST double *x,
  CONST int *nq,
  CONST double *hu,
  CONST double *h,
  int *ifail
);

extern void __stdcall D02XKF(
  CONST double *xsol,
  double sol[],
  CONST int *m,
  CONST double w[] /* 2 dimension */,
  CONST int *neqmax,
  CONST int *iw,
  CONST double w2[],
  CONST int *neq,
  CONST double *x,
  CONST int *nq,
  CONST double *hu,
  CONST double *h,
  int *ifail
);

extern void __stdcall D02YAF(
  CONST double *x,
  CONST double *h,
  CONST int *n,
  double y[],
  void (__stdcall *fcn)(double *, double[], double[]),
  double w[] /* 2 dimension */,
  CONST int *iw1,
  CONST int *iw2
);

extern double __stdcall D02ZAF(
  CONST int *n,
  CONST double v[],
  CONST double w[],
  int *ifail
);

extern void __stdcall D03EAF(
  CONST int *stage1,
  CONST int *ext,
  CONST int *dorm,
  CONST int *n,
  CONST double *p,
  CONST double *q,
  CONST double x[],
  CONST double y[],
  CONST int *n1p1,
  double phi[],
  double phid[],
  double *alpha,
  double c[] /* 2 dimension */,
  CONST int *ic,
  CONST int *np4,
  int icint[],
  CONST int *np1,
  int *ifail
);

extern void __stdcall D03EBF(
  CONST int *n1,
  CONST int *n2,
  CONST int *n1m,
  CONST double a[] /* 2 dimension */,
  CONST double b[] /* 2 dimension */,
  CONST double c[] /* 2 dimension */,
  CONST double d[] /* 2 dimension */,
  CONST double e[] /* 2 dimension */,
  CONST double q[] /* 2 dimension */,
  double t[] /* 2 dimension */,
  CONST double *aparam,
  CONST int *itmax,
  int *itcoun,
  int *itused,
  CONST int *ndir,
  CONST int *ixn,
  CONST int *iyn,
  CONST double *conres,
  CONST double *conchn,
  double resids[],
  double chngs[],
  double wrksp1[] /* 2 dimension */,
  double wrksp2[] /* 2 dimension */,
  double wrksp3[] /* 2 dimension */,
  int *ifail
);

extern void __stdcall D03ECF(
  CONST int *n1,
  CONST int *n2,
  CONST int *n3,
  CONST int *n1m,
  CONST int *n2m,
  CONST double a[] /* 3 dimension */,
  CONST double b[] /* 3 dimension */,
  CONST double c[] /* 3 dimension */,
  CONST double d[] /* 3 dimension */,
  CONST double e[] /* 3 dimension */,
  CONST double f[] /* 3 dimension */,
  CONST double g[] /* 3 dimension */,
  CONST double q[] /* 3 dimension */,
  double t[] /* 3 dimension */,
  CONST double *aparam,
  CONST int *itmax,
  int *itcoun,
  int *itused,
  CONST int *ndir,
  CONST int *ixn,
  CONST int *iyn,
  CONST int *izn,
  CONST double *conres,
  CONST double *conchn,
  double resids[],
  double chngs[],
  double wrksp1[] /* 3 dimension */,
  double wrksp2[] /* 3 dimension */,
  double wrksp3[] /* 3 dimension */,
  double wrksp4[] /* 3 dimension */,
  int *ifail
);

extern void __stdcall D03EDF(
  CONST int *ngx,
  CONST int *ngy,
  CONST int *lda,
  double a[] /* 2 dimension */,
  double rhs[],
  double ub[],
  CONST int *maxit,
  CONST double *acc,
  double us[],
  double u[],
  CONST int *iout,
  int *numit,
  int *ifail
);

extern void __stdcall D03EEF(
  CONST double *xmin,
  CONST double *xmax,
  CONST double *ymin,
  CONST double *ymax,
  void (__stdcall *pdef)(double *, double *, double *, double *, double *, double *, 
               double *, double *, double *),
  void (__stdcall *bndy)(double *, double *, double *, double *, double *, int *),
  CONST int *ngx,
  CONST int *ngy,
  CONST int *lda,
  double a[] /* 2 dimension */,
  double rhs[],
  CONST char *scheme, CONST int scheme_len,
  int *ifail
);

extern void __stdcall D03FAF(
  CONST double *xs,
  CONST double *xf,
  CONST int *l,
  CONST int *lbdcnd,
  CONST double bdxs[] /* 2 dimension */,
  CONST double bdxf[] /* 2 dimension */,
  CONST double *ys,
  CONST double *yf,
  CONST int *m,
  CONST int *mbdcnd,
  CONST double bdys[] /* 2 dimension */,
  CONST double bdyf[] /* 2 dimension */,
  CONST double *zs,
  CONST double *zf,
  CONST int *n,
  CONST int *nbdcnd,
  CONST double bdzs[] /* 2 dimension */,
  CONST double bdzf[] /* 2 dimension */,
  CONST double *lambda,
  CONST int *ldimf,
  CONST int *mdimf,
  double f[] /* 3 dimension */,
  double *pertrb,
  double w[],
  CONST int *lwrk,
  int *ifail
);

extern void __stdcall D03MAF(
  CONST double *h,
  CONST int *m,
  CONST int *n,
  CONST int *nb,
  int *npts,
  double places[] /* 2 dimension */,
  int index[] /* 2 dimension */,
  CONST int *idim,
  int (__stdcall *in)(double *, double *),
  double dist[] /* 2 dimension */,
  CONST int *ld,
  int *ifail
);

extern void __stdcall D03PAF(
  CONST int *m,
  CONST double *a,
  CONST double *b,
  double *ts,
  double *tout,
  double u[],
  CONST int *npts,
  CONST double *acc,
  double work[],
  CONST int *iwk,
  int *ind,
  int *ifail
);

extern void __stdcall D03PAZ(
  CONST int *npde,
  CONST double x[],
  CONST int *npts,
  CONST double *ts,
  CONST double *tlast,
  CONST double u[] /* 2 dimension */,
  CONST int *iu,
  CONST double *tout,
  CONST double *dt
);

extern void __stdcall D03PBF(
  CONST int *npde,
  CONST int *m,
  void (__stdcall *pdef)(int *, double *, double *, double[], double[], double[], 
               double[], double[]),
  void (__stdcall *bndy)(int *, double *, double[], int *, double[], double[], double[]),
  CONST double *a,
  CONST double *b,
  double *ts,
  double *tout,
  double u[] /* 2 dimension */,
  CONST int *npts,
  CONST int *imesh,
  double x[],
  CONST double *acc,
  double work[],
  CONST int *iwk,
  int *ind,
  int *ifail
);

extern void __stdcall D03PCF(
  CONST int *npde,
  CONST int *m,
  double *ts,
  double *tout,
  void (__stdcall *pdedef)(int *, double *, double *, double[], double[], double[], 
                 double[], double[], int *),
  void (__stdcall *bndary)(int *, double *, double[], double[], int *, double[], 
                 double[], int *),
  double u[] /* 2 dimension */,
  CONST int *npts,
  double x[],
  CONST double *acc,
  double w[],
  CONST int *nw,
  int iw[],
  CONST int *niw,
  CONST int *itask,
  CONST int *itrace,
  int *ind,
  int *ifail
);

extern void __stdcall D03PCK(
  CONST int *npde,
  CONST double *t,
  CONST int *nv,
  CONST double v[],
  CONST double vdot[],
  CONST int *nxi,
  CONST double xi[],
  CONST double u[] /* 2 dimension */,
  CONST double ux[] /* 2 dimension */,
  CONST double ri[] /* 2 dimension */,
  CONST double uti[] /* 2 dimension */,
  CONST double utxi[] /* 2 dimension */,
  CONST double vres[],
  CONST int *ires
);

extern void __stdcall D03PDF(
  CONST int *npde,
  CONST int *m,
  double *ts,
  double *tout,
  void (__stdcall *pdedef)(int *, double *, double[], int *, double[], double[], 
                 double[], double[], double[], int *),
  void (__stdcall *bndary)(int *, double *, double[], double[], int *, double[], 
                 double[], int *),
  double u[] /* 2 dimension */,
  CONST int *nbkpts,
  CONST double xbkpts[],
  CONST int *npoly,
  CONST int *npts,
  double x[],
  void (__stdcall *uinit)(int *, int *, double[], double[]),
  CONST double *acc,
  double w[],
  CONST int *nw,
  int iw[],
  CONST int *niw,
  CONST int *itask,
  CONST int *itrace,
  int *ind,
  int *ifail
);

extern void __stdcall D03PEF(
  CONST int *npde,
  double *ts,
  double *tout,
  void (__stdcall *pdedef)(),
  void (__stdcall *bndary)(),
  double u[] /* 2 dimension */,
  CONST int *npts,
  double x[],
  CONST int *nleft,
  CONST double *acc,
  double w[],
  CONST int *nw,
  int iw[],
  CONST int *niw,
  CONST int *itask,
  CONST int *itrace,
  int *ind,
  int *ifail
);

extern void __stdcall D03PEK(
  CONST int *npde,
  CONST double *t,
  CONST int *nv,
  CONST double v[],
  CONST double vdot[],
  CONST int *nxi,
  CONST double xi[],
  CONST double u[] /* 2 dimension */,
  CONST double ux[] /* 2 dimension */,
  CONST double uti[] /* 2 dimension */,
  CONST double vres[],
  CONST int *ires
);

extern void __stdcall D03PEL(
  CONST double *time,
  CONST int *nip,
  CONST int *npde,
  CONST double x[],
  CONST double u[] /* 2 dimension */,
  CONST double fmon[]
);

extern void __stdcall D03PFF(
  CONST int *npde,
  double *ts,
  double *tout,
  void (__stdcall *pdedef)(),
  void (__stdcall *numflx)(),
  void (__stdcall *bndary)(),
  double u[] /* 2 dimension */,
  CONST int *npts,
  double x[],
  CONST double acc[],
  CONST double *tsmax,
  double w[],
  CONST int *nw,
  int iw[],
  CONST int *niw,
  CONST int *itask,
  CONST int *itrace,
  int *ind,
  int *ifail
);

extern void __stdcall D03PFP(
  CONST int *npde,
  CONST double *t,
  CONST double *x,
  CONST double u[],
  CONST double ux[],
  double p[] /* 2 dimension */,
  double c[],
  double d[],
  double s[],
  CONST int *ires
);

extern void __stdcall D03PGF(
  CONST int *npde,
  CONST int *m,
  void (__stdcall *pdef)(int *, double[], double *, double[], double[], double[], 
               double[], double[]),
  void (__stdcall *bndy)(int *, double *, double[], int *, double[], double[], double[]),
  double *ts,
  double *tout,
  double u[] /* 2 dimension */,
  CONST int *iu,
  CONST int *npts,
  double x[],
  CONST double *relerr,
  CONST double *abserr,
  CONST int *inorm,
  void (__stdcall *montr)(int *, double[], int *, double *, double *, double[], int *, 
                double *, double *),
  CONST int *imon,
  int *iband,
  double work[],
  CONST int *iwk,
  int *ind,
  int *ifail
);

extern void __stdcall D03PHF(
  CONST int *npde,
  CONST int *m,
  double *ts,
  double *tout,
  void (__stdcall *pdedef)(int *, double *, double *, double[], double[], int *, 
                 double[], double[], double[], double[], double[], int *),
  void (__stdcall *bndary)(int *, double *, double[], double[], int *, double[], 
                 double[], int *, double[], double[], int *),
  double u[],
  CONST int *npts,
  double x[],
  CONST int *ncode,
  void (__stdcall *odedef)(int *, double *, int *, double[], double[], int *, double[], 
                 double[], double[], double[], double[], double[], double[], 
                 int *),
  CONST int *nxi,
  CONST double xi[],
  CONST int *neqn,
  CONST double rtol[],
  CONST double atol[],
  CONST int *itol,
  CONST char *norm, CONST int norm_len,
  CONST char *laopt, CONST int laopt_len,
  double algopt[],
  double w[],
  CONST int *nw,
  int iw[],
  CONST int *niw,
  CONST int *itask,
  CONST int *itrace,
  int *ind,
  int *ifail
);

extern void __stdcall D03PJF(
  CONST int *npde,
  CONST int *m,
  double *ts,
  double *tout,
  void (__stdcall *pdedef)(int *, double *, double[], int *, double[], double[], int *, 
                 double[], double[], double[], double[], double[], int *),
  void (__stdcall *bndary)(int *, double *, double[], double[], int *, double[], 
                 double[], int *, double[], double[], int *),
  double u[],
  CONST int *nbkpts,
  CONST double xbkpts[],
  CONST int *npoly,
  CONST int *npts,
  double x[],
  CONST int *ncode,
  void (__stdcall *odedef)(int *, double *, int *, double[], double[], int *, double[], 
                 double[], double[], double[], double[], double[], double[], 
                 int *),
  CONST int *nxi,
  CONST double xi[],
  CONST int *neqn,
  void (__stdcall *uvinit)(int *, int *, double[], double[], int *, double[]),
  CONST double rtol[],
  CONST double atol[],
  CONST int *itol,
  CONST char *norm, CONST int norm_len,
  CONST char *laopt, CONST int laopt_len,
  double algopt[],
  double w[],
  CONST int *nw,
  int iw[],
  CONST int *niw,
  CONST int *itask,
  CONST int *itrace,
  int *ind,
  int *ifail
);

extern void __stdcall D03PKF(
  CONST int *npde,
  CONST double *ts,
  CONST double *tout,
  void (__stdcall *pdedef)(),
  void (__stdcall *bndary)(),
  CONST double u[],
  CONST int *npts,
  CONST double x[],
  CONST int *nleft,
  CONST int *ncode,
  void (__stdcall *odedef)(),
  CONST int *nxi,
  CONST double xi[],
  CONST int *neqn,
  CONST double rtol[],
  CONST double atol[],
  CONST int *itol,
  CONST char *norm, CONST int norm_len,
  CONST char *laopt, CONST int laopt_len,
  double algopt[],
  CONST double w[],
  CONST int *nw,
  CONST int iw[],
  CONST int *niw,
  CONST int *itask,
  CONST int *itrace,
  CONST int *ind,
  int *ifail
);

extern void __stdcall D03PLF(
  CONST int *npde,
  double *ts,
  double *tout,
  void (__stdcall *pdedef)(),
  void (__stdcall *numflx)(),
  void (__stdcall *bndary)(),
  double u[],
  CONST int *npts,
  double x[],
  CONST int *ncode,
  void (__stdcall *odedef)(),
  CONST int *nxi,
  CONST double xi[],
  CONST int *neqn,
  CONST double rtol[],
  CONST double atol[],
  CONST int *itol,
  CONST char *norm, CONST int norm_len,
  CONST char *laopt, CONST int laopt_len,
  double algopt[],
  double w[],
  CONST int *nw,
  int iw[],
  CONST int *niw,
  CONST int *itask,
  CONST int *itrace,
  int *ind,
  int *ifail
);

extern void __stdcall D03PLP(
  CONST int *npde,
  CONST double *t,
  CONST double *x,
  CONST double u[],
  CONST double ux[],
  CONST int *nv,
  CONST double v[],
  CONST double vdot[],
  double p[] /* 2 dimension */,
  double c[],
  double d[],
  double s[],
  CONST int *ires
);

extern void __stdcall D03PPF(
  CONST int *npde,
  CONST int *m,
  CONST double *ts,
  CONST double *tout,
  void (__stdcall *pdedef)(),
  void (__stdcall *bndary)(),
  void (__stdcall *uvinit)(),
  CONST double u[],
  CONST int *npts,
  CONST double x[],
  CONST int *ncode,
  void (__stdcall *odedef)(),
  CONST int *nxi,
  CONST double xi[],
  CONST int *neqn,
  CONST double rtol[],
  CONST double atol[],
  CONST int *itol,
  CONST char *norm, CONST int norm_len,
  CONST char *laopt, CONST int laopt_len,
  double algopt[],
  CONST int *remesh,
  CONST int *nxfix,
  double xfix[],
  CONST int *nrmesh,
  CONST double *dxmesh,
  CONST double *trmesh,
  CONST int *ipminf,
  CONST double *xratio,
  CONST double *_const,
  void (__stdcall *monffd)(),
  CONST double w[],
  CONST int *nw,
  int iw[],
  CONST int *niw,
  CONST int *itask,
  CONST int *itrace,
  CONST int *ind,
  int *ifail
);

extern void __stdcall D03PRF(
  CONST int *npde,
  CONST double *ts,
  CONST double *tout,
  void (__stdcall *pdedef)(),
  void (__stdcall *bndary)(),
  void (__stdcall *uvinit)(),
  CONST double u[],
  CONST int *npts,
  CONST double x[],
  CONST int *nleft,
  CONST int *ncode,
  void (__stdcall *odedef)(),
  CONST int *nxi,
  CONST double xi[],
  CONST int *neqn,
  CONST double rtol[],
  CONST double atol[],
  CONST int *itol,
  CONST char *norm, CONST int norm_len,
  CONST char *laopt, CONST int laopt_len,
  double algopt[],
  CONST int *remesh,
  CONST int *nxfix,
  double xfix[],
  CONST int *nrmesh,
  CONST double *dxmesh,
  CONST double *trmesh,
  CONST int *ipminf,
  CONST double *xratio,
  CONST double *_const,
  void (__stdcall *monfkb)(),
  CONST double w[],
  CONST int *nw,
  int iw[],
  CONST int *niw,
  CONST int *itask,
  CONST int *itrace,
  CONST int *ind,
  int *ifail
);

extern void __stdcall D03PSF(
  CONST int *npde,
  double *ts,
  double *tout,
  void (__stdcall *pdedef)(),
  void (__stdcall *numflx)(),
  void (__stdcall *bndary)(),
  void (__stdcall *uvinit)(),
  double u[],
  CONST int *npts,
  double x[],
  CONST int *ncode,
  void (__stdcall *odedef)(),
  CONST int *nxi,
  CONST double xi[],
  CONST int *neqn,
  CONST double rtol[],
  CONST double atol[],
  CONST int *itol,
  CONST char *norm, CONST int norm_len,
  CONST char *laopt, CONST int laopt_len,
  double algopt[],
  int *remesh,
  CONST int *nxfix,
  double xfix[],
  CONST int *nrmesh,
  CONST double *dxmesh,
  CONST double *trmesh,
  CONST int *ipminf,
  CONST double *xratio,
  CONST double *_const,
  void (__stdcall *monitf)(),
  double w[],
  CONST int *nw,
  int iw[],
  CONST int *niw,
  CONST int *itask,
  CONST int *itrace,
  int *ind,
  int *ifail
);

extern void __stdcall D03PUF(
  CONST double uleft[],
  CONST double uright[],
  CONST double *gamma,
  double flux[],
  int *ifail
);

extern void __stdcall D03PVF(
  CONST double uleft[],
  CONST double uright[],
  CONST double *gamma,
  CONST char *path, CONST int path_len,
  double flux[],
  int *ifail
);

extern void __stdcall D03PWF(
  CONST double uleft[],
  CONST double uright[],
  CONST double *gamma,
  double flux[],
  int *ifail
);

extern void __stdcall D03PXF(
  CONST double uleft[],
  CONST double uright[],
  CONST double *gamma,
  double *tol,
  int *niter,
  double flux[],
  int *ifail
);

extern void __stdcall D03PYF(
  CONST int *npde,
  CONST double u[],
  CONST int *nbkpts,
  CONST double xbkpts[],
  CONST int *npoly,
  CONST int *npts,
  CONST double xp[],
  CONST int *intpts,
  CONST int *itype,
  double uout[] /* 3 dimension */,
  double w[],
  CONST int *nw,
  int *ifail
);

extern void __stdcall D03PZF(
  CONST int *npde,
  CONST int *m,
  CONST double u[] /* 2 dimension */,
  CONST int *npts,
  CONST double x[],
  CONST double xp[],
  CONST int *intpts,
  CONST int *itype,
  double uout[] /* 3 dimension */,
  int *ifail
);

extern void __stdcall D03RAF(
  CONST int *npde,
  double *ts,
  CONST double *tout,
  double dt[],
  double *xmin,
  double *xmax,
  double *ymin,
  double *ymax,
  CONST int *nx,
  CONST int *ny,
  CONST double *tols,
  CONST double *tolt,
  void (__stdcall *pdedef)(),
  void (__stdcall *bndary)(),
  void (__stdcall *pdeiv)(),
  void (__stdcall *monitr)(),
  int opti[],
  CONST double optr[] /* 2 dimension */,
  double rwk[],
  CONST int *lenrwk,
  int iwk[],
  CONST int *leniwk,
  int lwk[],
  CONST int *lenlwk,
  CONST int *itrace,
  int *ind,
  int *ifail
);

extern void __stdcall D03RBF(
  CONST int *npde,
  double *ts,
  CONST double *tout,
  double dt[],
  CONST double *tols,
  CONST double *tolt,
  void (__stdcall *inidom)(),
  void (__stdcall *pdedef)(),
  void (__stdcall *bndary)(),
  void (__stdcall *pdeiv)(),
  void (__stdcall *monitr)(),
  int opti[],
  CONST double optr[] /* 2 dimension */,
  double rwk[],
  CONST int *lenrwk,
  int iwk[],
  CONST int *leniwk,
  int lwk[],
  CONST int *lenlwk,
  CONST int *itrace,
  int *ind,
  int *ifail
);

extern void __stdcall D03RYF(
  CONST int *nx,
  CONST int *ny,
  CONST int *npts,
  CONST int *nrows,
  CONST int *nbnds,
  CONST int *nbpts,
  CONST int lrow[],
  CONST int irow[],
  CONST int icol[],
  CONST int llbnd[],
  CONST int ilbnd[],
  CONST int lbnd[],
  int iwk[],
  CONST int *leniwk,
  char *pgrid, CONST int pgrid_len,
  int *ifail
);

extern void __stdcall D03RZF(
  CONST int *level,
  CONST int *nlev,
  CONST double *xmin,
  CONST double *ymin,
  CONST double *dxb,
  CONST double *dyb,
  CONST int lgrid[],
  CONST int istruc[],
  int *npts,
  double x[],
  double y[],
  CONST int *lenxy,
  int *ifail
);

extern void __stdcall D03UAF(
  CONST int *n1,
  CONST int *n2,
  CONST int *n1m,
  CONST double a[] /* 2 dimension */,
  CONST double b[] /* 2 dimension */,
  CONST double c[] /* 2 dimension */,
  CONST double d[] /* 2 dimension */,
  CONST double e[] /* 2 dimension */,
  CONST double *aparam,
  CONST int *it,
  double r[] /* 2 dimension */,
  double wrksp1[] /* 2 dimension */,
  double wrksp2[] /* 2 dimension */,
  int *ifail
);

extern void __stdcall D03UBF(
  CONST int *n1,
  CONST int *n2,
  CONST int *n3,
  CONST int *n1m,
  CONST int *n2m,
  CONST double a[] /* 3 dimension */,
  CONST double b[] /* 3 dimension */,
  CONST double c[] /* 3 dimension */,
  CONST double d[] /* 3 dimension */,
  CONST double e[] /* 3 dimension */,
  CONST double f[] /* 3 dimension */,
  CONST double g[] /* 3 dimension */,
  CONST double *aparam,
  CONST int *it,
  double r[] /* 3 dimension */,
  double wrksp1[] /* 3 dimension */,
  double wrksp2[] /* 3 dimension */,
  double wrksp3[] /* 3 dimension */,
  int *ifail
);

extern void __stdcall D04AAF(
  CONST double *xval,
  CONST int *nder,
  CONST double *hbase,
  double der[],
  double erest[],
  double (__stdcall *fun)(double *),
  int *ifail
);

extern void __stdcall D05AAF(
  CONST double *lambda,
  CONST double *a,
  CONST double *b,
  double (__stdcall *k1)(double *, double *),
  double (__stdcall *k2)(double *, double *),
  double (__stdcall *g)(double *),
  double f[],
  double c[],
  CONST int *n,
  CONST int *ind,
  double w1[] /* 2 dimension */,
  double w2[] /* 2 dimension */,
  double wd[],
  CONST int *nmax,
  CONST int *mn,
  int *ifail
);

extern void __stdcall D05ABF(
  double (__stdcall *k)(double *, double *),
  double (__stdcall *g)(double *),
  CONST double *lambda,
  CONST double *a,
  CONST double *b,
  CONST int *odorev,
  CONST int *ev,
  CONST int *n,
  double cm[] /* 2 dimension */,
  double f1[] /* 2 dimension */,
  double wk[] /* 2 dimension */,
  CONST int *nmax,
  CONST int *nt2p1,
  double f[],
  double c[],
  int *ifail
);

extern void __stdcall D05BAF(
  double (__stdcall *ck)(double *),
  double (__stdcall *cg)(double *, double[]),
  double (__stdcall *cf)(double *),
  CONST char *method, CONST int method_len,
  CONST int *iorder,
  CONST double *alim,
  CONST double *tlim,
  double yn[],
  double errest[],
  CONST int *nout,
  CONST double *tol,
  CONST double *thresh,
  double work[],
  CONST int *iwk,
  int *ifail
);

extern void __stdcall D05BDF(
  double (__stdcall *ck)(),
  double (__stdcall *cf)(),
  double (__stdcall *cg)(),
  CONST char *initwt, CONST int initwt_len,
  CONST int *iorder,
  CONST double *tlim,
  CONST double *tolnl,
  CONST int *nmesh,
  double yn[],
  double work[],
  CONST int *lwk,
  int nct[],
  int *ifail
);

extern void __stdcall D05BEF(
  double (__stdcall *ck)(),
  double (__stdcall *cf)(),
  double (__stdcall *cg)(),
  CONST char *initwt, CONST int initwt_len,
  CONST int *iorder,
  CONST double *tlim,
  CONST double *tolnl,
  CONST int *nmesh,
  double yn[],
  double work[],
  CONST int *lwk,
  int nct[],
  int *ifail
);

extern void __stdcall D05BWF(
  CONST char *method, CONST int method_len,
  CONST int *iorder,
  double omega[],
  CONST int *nomg,
  int *lensw,
  double sw[] /* 2 dimension */,
  CONST int *ldsw,
  CONST int *nwt,
  int *ifail
);

extern void __stdcall D05BYF(
  CONST int *iorder,
  CONST int *iq,
  CONST int *lenfw,
  double wt[],
  double sw[] /* 2 dimension */,
  CONST int *ldsw,
  double work[],
  CONST int *lwk,
  int *ifail
);

extern void __stdcall E01AAF(
  double a[],
  double b[],
  double c[],
  CONST int *n1,
  CONST int *n2,
  CONST int *n,
  CONST double *x
);

extern void __stdcall E01ABF(
  CONST int *n,
  CONST double *p,
  double a[],
  double g[],
  CONST int *n1,
  CONST int *n2,
  int *ifail
);

extern void __stdcall E01AEF(
  CONST int *m,
  CONST double *xmin,
  CONST double *xmax,
  CONST double x[],
  CONST double y[],
  CONST int ip[],
  CONST int *n,
  CONST int *itmin,
  CONST int *itmax,
  double a[],
  double wrk[],
  CONST int *lwrk,
  int iwrk[],
  CONST int *liwrk,
  int *ifail
);

extern void __stdcall E01BAF(
  CONST int *m,
  CONST double x[],
  CONST double y[],
  double k[],
  double c[],
  CONST int *lck,
  double wrk[],
  CONST int *lwrk,
  int *ifail
);

extern void __stdcall E01BEF(
  CONST int *n,
  CONST double x[],
  CONST double f[],
  double d[],
  int *ifail
);

extern void __stdcall E01BFF(
  CONST int *n,
  CONST double x[],
  CONST double f[],
  CONST double d[],
  CONST int *m,
  CONST double px[],
  double pf[],
  int *ifail
);

extern void __stdcall E01BGF(
  CONST int *n,
  CONST double x[],
  CONST double f[],
  CONST double d[],
  CONST int *m,
  CONST double px[],
  double pf[],
  double pd[],
  int *ifail
);

extern void __stdcall E01BHF(
  CONST int *n,
  CONST double x[],
  CONST double f[],
  CONST double d[],
  CONST double *a,
  CONST double *b,
  double *pint,
  int *ifail
);

extern void __stdcall E01DAF(
  CONST int *mx,
  CONST int *my,
  CONST double x[],
  CONST double y[],
  CONST double f[],
  int *px,
  int *py,
  double lamda[],
  double mu[],
  double c[],
  double wrk[],
  int *ifail
);

extern void __stdcall E01RAF(
  CONST int *n,
  CONST double x[],
  CONST double f[],
  int *m,
  double a[],
  double u[],
  int iw[],
  int *ifail
);

extern void __stdcall E01RBF(
  CONST int *m,
  CONST double a[],
  CONST double u[],
  CONST double *x,
  double *f,
  int *ifail
);

extern void __stdcall E01SAF(
  CONST int *m,
  CONST double x[],
  CONST double y[],
  CONST double f[],
  int triang[],
  double grads[] /* 2 dimension */,
  int *ifail
);

extern void __stdcall E01SBF(
  CONST int *m,
  CONST double x[],
  CONST double y[],
  CONST double f[],
  CONST int triang[],
  CONST double grads[] /* 2 dimension */,
  CONST double *px,
  CONST double *py,
  double *pf,
  int *ifail
);

extern void __stdcall E01SEF(
  CONST int *m,
  CONST double x[],
  CONST double y[],
  CONST double f[],
  double *rnw,
  double *rnq,
  CONST int *nw,
  CONST int *nq,
  double fnodes[],
  int *minnq,
  double wrk[],
  int *ifail
);

extern void __stdcall E01SFF(
  CONST int *m,
  CONST double x[],
  CONST double y[],
  CONST double f[],
  CONST double *rnw,
  CONST double fnodes[] /* 2 dimension */,
  CONST double *px,
  CONST double *py,
  double *pf,
  int *ifail
);

extern void __stdcall E01SGF(
  CONST int *m,
  CONST double x[],
  CONST double y[],
  CONST double f[],
  CONST int *nw,
  CONST int *nq,
  int iq[],
  CONST int *liq,
  double rq[],
  CONST int *lrq,
  int *ifail
);

extern void __stdcall E01SHF(
  CONST int *m,
  CONST double x[],
  CONST double y[],
  CONST double f[],
  CONST int iq[],
  CONST int *liq,
  CONST double rq[],
  CONST int *lrq,
  CONST int *n,
  CONST double u[],
  CONST double v[],
  double q[],
  double qx[],
  double qy[],
  int *ifail
);

extern void __stdcall E01TGF(
  CONST int *m,
  CONST double x[],
  CONST double y[],
  CONST double z[],
  CONST double f[],
  CONST int *nw,
  CONST int *nq,
  int iq[],
  CONST int *liq,
  double rq[],
  CONST int *lrq,
  int *ifail
);

extern void __stdcall E01THF(
  CONST int *m,
  CONST double x[],
  CONST double y[],
  CONST double z[],
  CONST double f[],
  CONST int iq[],
  CONST int *liq,
  CONST double rq[],
  CONST int *lrq,
  CONST int *n,
  CONST double u[],
  CONST double v[],
  CONST double w[],
  double q[],
  double qx[],
  double qy[],
  double qz[],
  int *ifail
);

extern void __stdcall E02ACF(
  CONST double x[],
  CONST double y[],
  CONST int *n,
  double aa[],
  CONST int *m1,
  double *ref
);

extern void __stdcall E02ADF(
  CONST int *m,
  CONST int *kplus1,
  CONST int *nrows,
  CONST double x[],
  CONST double y[],
  CONST double w[],
  double work1[] /* 2 dimension */,
  double work2[] /* 2 dimension */,
  double a[] /* 2 dimension */,
  double s[],
  int *ifail
);

extern void __stdcall E02AEF(
  CONST int *nplus1,
  CONST double a[],
  CONST double *xcap,
  double *p,
  int *ifail
);

extern void __stdcall E02AFF(
  CONST int *nplus1,
  CONST double f[],
  double a[],
  int *ifail
);

extern void __stdcall E02AGF(
  CONST int *m,
  CONST int *kplus1,
  CONST int *nrows,
  CONST double *xmin,
  CONST double *xmax,
  CONST double x[],
  CONST double y[],
  CONST double w[],
  CONST int *mf,
  CONST double xf[],
  CONST double yf[],
  CONST int *lyf,
  CONST int ip[],
  double a[] /* 2 dimension */,
  double s[],
  int *np1,
  double wrk[],
  CONST int *lwrk,
  int iwrk[],
  CONST int *liwrk,
  int *ifail
);

extern void __stdcall E02AHF(
  CONST int *np1,
  CONST double *xmin,
  CONST double *xmax,
  CONST double a[],
  CONST int *ia1,
  CONST int *la,
  double *patm1,
  double adif[],
  CONST int *iadif1,
  CONST int *ladif,
  int *ifail
);

extern void __stdcall E02AJF(
  CONST int *np1,
  CONST double *xmin,
  CONST double *xmax,
  CONST double a[],
  CONST int *ia1,
  CONST int *la,
  CONST double *qatm1,
  double ain[],
  CONST int *iaint1,
  CONST int *laint,
  int *ifail
);

extern void __stdcall E02AKF(
  CONST int *np1,
  CONST double *xmin,
  CONST double *xmax,
  CONST double a[],
  CONST int *ia1,
  CONST int *la,
  CONST double *x,
  double *result,
  int *ifail
);

extern void __stdcall E02BAF(
  CONST int *m,
  CONST int *ncap7,
  CONST double x[],
  CONST double y[],
  CONST double w[],
  double k[],
  double work1[],
  double work2[] /* 2 dimension */,
  double c[],
  double *ss,
  int *ifail
);

extern void __stdcall E02BBF(
  CONST int *ncap7,
  CONST double k[],
  CONST double c[],
  CONST double *x,
  double *s,
  int *ifail
);

extern void __stdcall E02BCF(
  CONST int *ncap7,
  CONST double k[],
  CONST double c[],
  CONST double *x,
  CONST int *left,
  double s[],
  int *ifail
);

extern void __stdcall E02BDF(
  CONST int *ncap7,
  CONST double k[],
  CONST double c[],
  double *defint,
  int *ifail
);

extern void __stdcall E02BEF(
  CONST char *start, CONST int start_len,
  CONST int *m,
  CONST double x[],
  CONST double y[],
  CONST double w[],
  CONST double *s,
  CONST int *nest,
  int *n,
  double k[],
  double c[],
  double *fp,
  double wrk[],
  CONST int *lwrk,
  int iwrk[],
  int *ifail
);

extern void __stdcall E02CAF(
  CONST int m[],
  CONST int *n,
  CONST int *k,
  CONST int *l,
  CONST double x[],
  CONST double y[],
  CONST double f[],
  CONST double w[],
  CONST int *nx,
  double a[],
  CONST int *na,
  CONST double xmin[],
  CONST double xmax[],
  CONST double nux[],
  CONST int *inuxp1,
  CONST double nuy[],
  CONST int *inuyp1,
  double work[],
  CONST int *nwork,
  int *ifail
);

extern void __stdcall E02CBF(
  CONST int *mfirst,
  CONST int *mlast,
  CONST int *k,
  CONST int *l,
  CONST double x[],
  CONST double *xmin,
  CONST double *xmax,
  CONST double *y,
  CONST double *ymin,
  CONST double *ymax,
  double ff[],
  CONST double a[],
  CONST int *na,
  double work[],
  CONST int *nwork,
  int *ifail
);

extern void __stdcall E02DAF(
  CONST int *m,
  CONST int *px,
  CONST int *py,
  CONST double x[],
  CONST double y[],
  CONST double f[],
  CONST double w[],
  double lamda[],
  double mu[],
  CONST int point[],
  CONST int *npoint,
  double dl[],
  double c[],
  CONST int *nc,
  double ws[],
  CONST int *nws,
  CONST double *eps,
  double *sigma,
  int *rank,
  int *ifail
);

extern void __stdcall E02DBF(
  CONST int *m,
  CONST int *px,
  CONST int *py,
  CONST double x[],
  CONST double y[],
  double ff[],
  CONST double lamda[],
  CONST double mu[],
  CONST int point[],
  CONST int *npoint,
  CONST double c[],
  CONST int *nc,
  int *ifail
);

extern void __stdcall E02DCF(
  CONST char *start, CONST int start_len,
  CONST int *mx,
  CONST double x[],
  CONST int *my,
  CONST double y[],
  CONST double f[],
  CONST double *s,
  CONST int *nxest,
  CONST int *nyest,
  int *nx,
  double lamda[],
  int *ny,
  double mu[],
  double c[],
  double *fp,
  double wrk[],
  CONST int *lwrk,
  int iwrk[],
  CONST int *liwrk,
  int *ifail
);

extern void __stdcall E02DDF(
  CONST char *start, CONST int start_len,
  CONST int *m,
  double x[],
  double y[],
  CONST double f[],
  CONST double w[],
  CONST double *s,
  CONST int *nxest,
  CONST int *nyest,
  int *nx,
  double lamda[],
  int *ny,
  double mu[],
  double c[],
  double *fp,
  int *rank,
  double wrk[],
  CONST int *lwrk,
  int iwrk[],
  CONST int *liwrk,
  int *ifail
);

extern void __stdcall E02DEF(
  CONST int *m,
  CONST int *px,
  CONST int *py,
  CONST double x[],
  CONST double y[],
  CONST double lamda[],
  CONST double mu[],
  CONST double c[],
  double ff[],
  double wrk[],
  int iwrk[],
  int *ifail
);

extern void __stdcall E02DFF(
  CONST int *mx,
  CONST int *my,
  CONST int *px,
  CONST int *py,
  CONST double x[],
  CONST double y[],
  CONST double lamda[],
  CONST double mu[],
  CONST double c[],
  double ff[],
  double wrk[],
  CONST int *lwrk,
  int iwrk[],
  CONST int *liwrk,
  int *ifail
);

extern void __stdcall E02GAF(
  CONST int *m,
  double a[] /* 2 dimension */,
  CONST int *la,
  double b[],
  CONST int *nplus2,
  CONST double *tol,
  double x[],
  double *resid,
  int *irank,
  int *iter,
  int iwork[],
  int *ifail
);

extern void __stdcall E02GBF(
  CONST int *m,
  CONST int *n,
  CONST int *mpl1,
  double e[] /* 2 dimension */,
  CONST int *ier,
  CONST double f[],
  double x[],
  CONST int *mxs1,
  void (__stdcall *monit)(int *, double[], int *, int *, double *),
  CONST int *iprint,
  int *k,
  double *el1n,
  int indx[],
  double w[],
  CONST int *iw,
  int *ifail
);

extern void __stdcall E02GCF(
  CONST int *m,
  CONST int *n,
  CONST int *mdim,
  CONST int *ndim,
  double a[] /* 2 dimension */,
  double b[],
  CONST double *tol1,
  double *reler,
  double x[],
  double *resmax,
  int *irank,
  int *iter,
  int *ifail
);

extern void __stdcall E02RAF(
  CONST int *ia,
  CONST int *ib,
  CONST double c[],
  CONST int *ic,
  double a[],
  double b[],
  double w[],
  CONST int *jw,
  int *ifail
);

extern void __stdcall E02RBF(
  CONST double a[],
  CONST int *ia,
  CONST double b[],
  CONST int *ib,
  CONST double *x,
  double *ans,
  int *ifail
);

extern void __stdcall E02ZAF(
  CONST int *px,
  CONST int *py,
  CONST double lamda[],
  CONST double mu[],
  CONST int *m,
  CONST double x[],
  CONST double y[],
  int point[],
  CONST int *npoint,
  int adres[],
  CONST int *nadres,
  int *ifail
);

extern void __stdcall E04ABF(
  void (__stdcall *fun)(double *, double *),
  double *eps,
  double *t,
  double *a,
  double *b,
  int *maxcal,
  double *x,
  double *f,
  int *ifail
);

extern void __stdcall E04BBF(
  void (__stdcall *fun)(double *, double *, double *),
  double *eps,
  double *t,
  double *a,
  double *b,
  int *maxcal,
  double *x,
  double *f,
  double *g,
  int *ifail
);

extern void __stdcall E04CCF(
  CONST int *n,
  double x[],
  double *fmin,
  CONST double *eps,
  CONST int *n1,
  double pdstar[],
  double pstar[],
  double pbar[],
  double step[],
  double y[],
  double p[] /* 2 dimension */,
  void (__stdcall *funct)(int *, double[], double *),
  void (__stdcall *monit)(double *, double *, double[], int *, int *, int *),
  CONST int *maxit,
  int *ifail
);

extern void __stdcall E04DGF(
  CONST int *n,
  void (__stdcall *fungrd)(int *, int *, double[], double *, double[], int *, int[], 
                 double[]),
  int *iter,
  double *objf,
  double objgrd[],
  double x[],
  int iwork[],
  double work[],
  CONST int iuser[],
  CONST double user[],
  int *ifail
);

extern void __stdcall E04DJF(
  CONST int *ioptns,
  int *inform
);

extern void __stdcall E04DKF(
  CONST char *string, int string_len
);

extern void __stdcall E04FCF(
  CONST int *m,
  CONST int *n,
  void (__stdcall *lsfn1)(int *, int *, int *, double[], double[], int[], int *, 
                double[], int *),
  void (__stdcall *lsmon)(int *, int *, double[], double[], double[], int *, double[], 
                int *, int *, int *, int[], int *, double[], int *),
  CONST int *iprint,
  CONST int *maxcal,
  CONST double *eta,
  CONST double *xtol,
  CONST double *stepmx,
  double x[],
  double *fsumsq,
  double fvec[],
  double fjac[] /* 2 dimension */,
  CONST int *lj,
  double s[],
  double vt[] /* 2 dimension */,
  CONST int *lvt,
  int *niter,
  int *nftotl,
  CONST int iw[],
  CONST int *liw,
  double w[],
  CONST int *lw,
  int *ifail
);

extern void __stdcall E04FCV(
  int *m,
  int *n,
  void (__stdcall *lsfjs)(int *, int *, int *, void(*), double[], double[], double[], 
                int *, int[], int *, double[], int *),
  void (__stdcall *lsf)(),
  CONST double *eps,
  CONST double *t,
  CONST double *eta,
  CONST double *sftbnd,
  CONST double *xlamda,
  CONST double p[],
  CONST double *gtp,
  double x[],
  double *f,
  double *alpha,
  CONST double fjac[] /* 2 dimension */,
  CONST int *lj,
  double fvec[],
  CONST double g[],
  int *nftotl,
  int *iflag,
  CONST int iw[],
  CONST int *liw,
  double w[],
  CONST int *lw
);

extern void __stdcall E04FDF(
  CONST int *m,
  CONST int *n,
  double x[],
  double *fsumsq,
  CONST int iw[],
  CONST int *liw,
  double w[],
  CONST int *lw,
  int *ifail
);

extern void __stdcall E04FYF(
  CONST int *m,
  CONST int *n,
  void (__stdcall *lsfun1)(),
  double x[],
  double *fsumsq,
  double w[],
  CONST int *lw,
  CONST int iuser[],
  CONST double user[],
  int *ifail
);

extern void __stdcall E04GBF(
  int *m,
  int *n,
  void (__stdcall *lsqlin)(int *, int *, void(*), void(*), double *, double *, double *, 
                 double *, double *, double[], double *, double[], double *, 
                 double *, double[], int *, double[], double[], int *, int *, 
                 int[], int *, double[], int *),
  void (__stdcall *lsfjc)(int *, int *, int *, double[], double[], double[], int *, 
                int[], int *, double[], int *),
  void (__stdcall *lsmon)(int *, int *, double[], double[], double[], int *, double[], 
                int *, int *, int *, int[], int *, double[], int *),
  CONST int *iprint,
  CONST int *maxcal,
  CONST double *eta,
  CONST double *xtol,
  CONST double *stepmx,
  double x[],
  double *fsumsq,
  double fvec[],
  double fjac[] /* 2 dimension */,
  int *lj,
  double s[],
  double vt[] /* 2 dimension */,
  CONST int *lvt,
  int *niter,
  int *nftotl,
  CONST int iw[],
  CONST int *liw,
  double w[],
  CONST int *lw,
  int *ifail
);

extern void __stdcall E04GCF(
  CONST int *m,
  CONST int *n,
  double x[],
  double *fsumsq,
  CONST int iw[],
  CONST int *liw,
  double w[],
  CONST int *lw,
  int *ifail
);

extern void __stdcall E04GDF(
  CONST int *m,
  CONST int *n,
  void (__stdcall *lsfjc)(int *, int *, int *, double[], double[], double[], int *, 
                int[], int *, double[], int *),
  void (__stdcall *lsmon)(int *, int *, double[], double[], double[], int *, double[], 
                int *, int *, int *, int[], int *, double[], int *),
  CONST int *iprint,
  CONST int *maxfun,
  CONST double *eta,
  CONST double *xtol,
  CONST double *stepmx,
  double x[],
  double *fsumsq,
  double fvec[],
  double fjac[] /* 2 dimension */,
  CONST int *lj,
  double s[],
  double vt[] /* 2 dimension */,
  CONST int *lvt,
  int *niter,
  int *nftotl,
  CONST int iw[],
  CONST int *liw,
  double w[],
  CONST int *lw,
  int *ifail
);

extern void __stdcall E04GEF(
  CONST int *m,
  CONST int *n,
  double x[],
  double *fsumsq,
  CONST int iw[],
  CONST int *liw,
  double w[],
  CONST int *lw,
  int *ifail
);

extern void __stdcall E04GYF(
  CONST int *m,
  CONST int *n,
  void (__stdcall *lsfun2)(int *, int *, double[], double[], double[], int *, 
                        int[], double[]),
  double x[],
  double *fsumsq,
  double w[],
  CONST int *lw,
  CONST int iuser[],
  CONST double user[],
  int *ifail
);

extern void __stdcall E04GZF(
  CONST int *m,
  CONST int *n,
  void (__stdcall *lsfun2)(int *, int *, double[], double[], double[], int *, 
                        int[], double[]),
  double x[],
  double *fsumsq,
  double w[],
  CONST int *lw,
  CONST int iuser[],
  CONST double user[],
  int *ifail
);

extern void __stdcall E04HBF(
  int *n,
  void (__stdcall *sfun)(int *, int *, double[], double *, double[], int[], int *, 
               double[], int *),
  double x[],
  int *nf,
  double delta[],
  double hesl[],
  CONST int *lh,
  double hesd[],
  double *f,
  double g[],
  CONST int iw[],
  CONST int *liw,
  double w[],
  CONST int *lw,
  int *ifail
);

extern void __stdcall E04HCF(
  int *n,
  void (__stdcall *sfun)(int *, int *, double[], double *, double[], int[], int *, 
               double[], int *),
  double x[],
  double *f,
  double g[],
  CONST int iw[],
  CONST int *liw,
  double w[],
  CONST int *lw,
  int *ifail
);

extern void __stdcall E04HDF(
  int *n,
  void (__stdcall *sfun)(int *, int *, double[], double *, double[], int[], int *, 
               double[], int *),
  void (__stdcall *shess)(int *, int *, double[], double[], int *, double[], int[], 
                int *, double[], int *),
  double x[],
  double g[],
  double hesl[],
  int *lh,
  double hesd[],
  CONST int iw[],
  CONST int *liw,
  double w[],
  CONST int *lw,
  int *ifail
);

extern void __stdcall E04HEF(
  CONST int *m,
  CONST int *n,
  void (__stdcall *lsfjc)(int *, int *, int *, double[], double[], double[], int *, 
                int[], int *, double[], int *),
  void (__stdcall *lshes)(int *, int *, int *, double[], double[], double[], int *, 
                int[], int *, double[], int *),
  void (__stdcall *lsmon)(int *, int *, double[], double[], double[], int *, double[], 
                int *, int *, int *, int[], int *, double[], int *),
  CONST int *iprint,
  CONST int *maxcal,
  CONST double *eta,
  CONST double *xtol,
  CONST double *stepmx,
  double x[],
  double *fsumsq,
  double fvec[],
  double fjac[] /* 2 dimension */,
  CONST int *lj,
  double s[],
  double vt[] /* 2 dimension */,
  CONST int *lvt,
  int *niter,
  int *nftotl,
  CONST int iw[],
  CONST int *liw,
  double w[],
  CONST int *lw,
  int *ifail
);

extern void __stdcall E04HEV(
  int *m,
  int *n,
  void (__stdcall *lsfjs)(int *, int *, int *, void(*), double[], double[], double[], 
                int *, int[], int *, double[], int *),
  void (__stdcall *lsf)(),
  CONST double *eps,
  CONST double *t,
  CONST double *eta,
  CONST double *sftbnd,
  CONST double *xlamda,
  CONST double p[],
  CONST double *gtp,
  double x[],
  double *f,
  double *alpha,
  double fjac[] /* 2 dimension */,
  CONST int *lj,
  double fvec[],
  double g[],
  int *nftotl,
  int *iflag,
  CONST int iw[],
  CONST int *liw,
  double w[],
  CONST int *lw
);

extern void __stdcall E04HFF(
  CONST int *m,
  CONST int *n,
  double x[],
  double *fsumsq,
  CONST int iw[],
  CONST int *liw,
  double w[],
  CONST int *lw,
  int *ifail
);

extern void __stdcall E04HYF(
  CONST int *m,
  CONST int *n,
  void (__stdcall *lsfun2)(int *, int *, double[], double[], double[], int *, 
                        int[], double[]),
  void (__stdcall *lshes2)(int *, int *, double[], double[], double[], int *, 
                        int[], double[]),
  double x[],
  double *fsumsq,
  double w[],
  CONST int *lw,
  CONST int iuser[],
  CONST double user[],
  int *ifail
);

extern void __stdcall E04JAF(
  CONST int *n,
  CONST int *ibound,
  double bl[],
  double bu[],
  double x[],
  double *f,
  int iw[],
  CONST int *liw,
  double w[],
  CONST int *lw,
  int *ifail
);

extern void __stdcall E04JBF(
  int *n,
  void (__stdcall *sfun)(int *, int *, double[], double *, double[], int[], int *, 
               double[], int *),
  void (__stdcall *monit)(int *, double[], double *, double[], int[], double *, double *, 
                int *, int *, int *, int[], int *, double[], int *),
  CONST int *iprint,
  CONST int *locsch,
  CONST int *intyp,
  void (__stdcall *minlin)(int *, void(*), double *, double *, double *, double *, 
                 double *, double[], double *, double[], double *, double *, 
                 double[], int *, int *, int[], int *, double[], int *),
  CONST int *maxfun,
  CONST double *eta,
  CONST double *xtol,
  CONST double *stepmx,
  CONST double *fest,
  CONST double delta[],
  CONST int *ibound,
  double bl[],
  double bu[],
  double x[],
  double hesl[],
  CONST int *lh,
  double hesd[],
  int istate[],
  double *f,
  double g[],
  int iw[],
  CONST int *liw,
  double w[],
  CONST int *lw,
  int *ifail
);

extern void __stdcall E04JBQ(
  int *n,
  void (__stdcall *sfun)(int *, int *, double[], double *, double[], int[], int *, 
               double[], int *),
  CONST double *eps,
  CONST double *t,
  CONST double *eta,
  CONST double *sftbnd,
  CONST double *xlamda,
  CONST double p[],
  CONST double *gtp,
  double x[],
  double *f,
  double *alpha,
  CONST double g[],
  int *nftotl,
  int *iflag,
  CONST int iw[],
  CONST int *liw,
  double w[],
  CONST int *lw
);

extern void __stdcall E04JYF(
  CONST int *n,
  CONST int *ibound,
  void (__stdcall *funct1)(int *, double[], double *, int[], double[]),
  double bl[],
  double bu[],
  double x[],
  double *f,
  int iw[],
  CONST int *liw,
  double w[],
  CONST int *lw,
  CONST int iuser[],
  CONST double user[],
  int *ifail
);

extern void __stdcall E04KAF(
  CONST int *n,
  CONST int *ibound,
  double bl[],
  double bu[],
  double x[],
  double *f,
  double g[],
  int iw[],
  CONST int *liw,
  double w[],
  CONST int *lw,
  int *ifail
);

extern void __stdcall E04KBF(
  int *n,
  void (__stdcall *sfun)(int *, int *, double[], double *, double[], int[], int *, 
               double[], int *),
  void (__stdcall *monit)(int *, double[], double *, double[], int[], double *, double *, 
                int *, int *, int *, int[], int *, double[], int *),
  CONST int *iprint,
  CONST int *locsch,
  CONST int *intyp,
  void (__stdcall *minlin)(int *, void(*), double *, double *, double *, double *, 
                 double *, double[], double *, double[], double *, double *, 
                 double[], int *, int *, int[], int *, double[], int *),
  CONST int *maxfun,
  CONST double *eta,
  CONST double *xtol,
  CONST double *stepmx,
  CONST double *fest,
  CONST int *ibound,
  double bl[],
  double bu[],
  double x[],
  double hesl[],
  CONST int *lh,
  double hesd[],
  int istate[],
  double *f,
  double g[],
  int iw[],
  CONST int *liw,
  double w[],
  CONST int *lw,
  int *ifail
);

extern void __stdcall E04KCF(
  CONST int *n,
  CONST int *ibound,
  double bl[],
  double bu[],
  double x[],
  double *f,
  double g[],
  int iw[],
  CONST int *liw,
  double w[],
  CONST int *lw,
  int *ifail
);

extern void __stdcall E04KDF(
  CONST int *n,
  void (__stdcall *sfun)(int *, int *, double[], double *, double[], int[], int *, 
               double[], int *),
  void (__stdcall *monit)(int *, double[], double *, double[], int[], double *, double *, 
                int *, int *, int *, int[], int *, double[], int *),
  CONST int *iprint,
  CONST int *maxfun,
  CONST double *eta,
  CONST double *xtol,
  CONST double *delta,
  CONST double *stepmx,
  CONST int *ibound,
  double bl[],
  double bu[],
  double x[],
  double hesl[],
  CONST int *lh,
  double hesd[],
  int istate[],
  double *f,
  double g[],
  int iw[],
  CONST int *liw,
  double w[],
  CONST int *lw,
  int *ifail
);

extern void __stdcall E04KYF(
  CONST int *n,
  CONST int *ibound,
  void (__stdcall *funct2)(int *, double[], double *, double[], int[], double[]),
  double bl[],
  double bu[],
  double x[],
  double *f,
  double g[],
  int iw[],
  CONST int *liw,
  double w[],
  CONST int *lw,
  CONST int iuser[],
  CONST double user[],
  int *ifail
);

extern void __stdcall E04KZF(
  CONST int *n,
  CONST int *ibound,
  void (__stdcall *funct2)(int *, double[], double *, double[], int[], double[]),
  double bl[],
  double bu[],
  double x[],
  double *f,
  double g[],
  int iw[],
  CONST int *liw,
  double w[],
  CONST int *lw,
  CONST int iuser[],
  CONST double user[],
  int *ifail
);

extern void __stdcall E04LAF(
  CONST int *n,
  CONST int *ibound,
  double bl[],
  double bu[],
  double x[],
  double *f,
  double g[],
  int iw[],
  CONST int *liw,
  double w[],
  CONST int *lw,
  int *ifail
);

extern void __stdcall E04LBF(
  CONST int *n,
  void (__stdcall *sfun)(int *, int *, double[], double *, double[], int[], int *, 
               double[], int *),
  void (__stdcall *shess)(int *, int *, double[], double[], int *, double[], int[], 
                int *, double[], int *),
  void (__stdcall *monit)(int *, double[], double *, double[], int[], double *, double *, 
                int *, int *, int *, int[], int *, double[], int *),
  CONST int *iprint,
  CONST int *maxfun,
  CONST double *eta,
  CONST double *xtol,
  CONST double *stepmx,
  CONST int *ibound,
  double bl[],
  double bu[],
  double x[],
  double hesl[],
  CONST int *lh,
  double hesd[],
  int istate[],
  double *f,
  double g[],
  int iw[],
  CONST int *liw,
  double w[],
  CONST int *lw,
  int *ifail
);

extern void __stdcall E04LBS(
  int *n,
  void (__stdcall *sfun)(int *, int *, double[], double *, double[], int[], int *, 
               double[], int *),
  CONST double *eps,
  CONST double *t,
  CONST double *eta,
  CONST double *sftbnd,
  CONST double *xlamda,
  CONST double p[],
  CONST double *gtp,
  double x[],
  double *f,
  double *alpha,
  double g[],
  int *nftotl,
  int *iflag,
  CONST int iw[],
  CONST int *liw,
  double w[],
  CONST int *lw
);

extern void __stdcall E04LYF(
  CONST int *n,
  CONST int *ibound,
  void (__stdcall *funct2)(int *, double[], double *, double[], int[], double[]),
  void (__stdcall *hess2)(int *, double[], double[], int *, double[], int[], 
                       double[]),
  double bl[],
  double bu[],
  double x[],
  double *f,
  double g[],
  int iw[],
  CONST int *liw,
  double w[],
  CONST int *lw,
  CONST int iuser[],
  CONST double user[],
  int *ifail
);

extern void __stdcall E04MBF(
  CONST int *itmax,
  CONST int *msglvl,
  CONST int *n,
  CONST int *nclin,
  CONST int *nctotl,
  CONST int *nrowa,
  CONST double a[] /* 2 dimension */,
  CONST double bl[],
  CONST double bu[],
  CONST double cvec[],
  CONST int *linobj,
  double x[],
  int istate[],
  double *objlp,
  double clamda[],
  int iw[],
  CONST int *leniw,
  double w[],
  CONST int *lenw,
  int *ifail
);

extern void __stdcall E04MFF(
  CONST int *n,
  CONST int *nclin,
  CONST double a[] /* 2 dimension */,
  CONST int *lda,
  CONST double bl[],
  CONST double bu[],
  CONST double cvec[],
  int istate[],
  double x[],
  int *iter,
  double *obj,
  double ax[],
  double clamda[],
  int iw[],
  CONST int *leniw,
  double w[],
  CONST int *lenw,
  int *ifail
);

extern void __stdcall E04MGF(
  CONST int *ioptns,
  int *inform
);

extern void __stdcall E04MHF(
  CONST char *string, int string_len
);

extern void __stdcall E04MZF(
  CONST int *infile,
  CONST int *maxn,
  CONST int *maxm,
  CONST int *maxnnz,
  CONST double *xbldef,
  CONST double *xbudef,
  CONST int *mpslst,
  int *n,
  int *m,
  int *nnz,
  int *iobj,
  int *ncolh,
  double a[],
  int ha[],
  int ka[],
  double bl[],
  double bu[],
  char *start, CONST int start_len,
  char *names, CONST int names_len,
  int *nname,
  char *crname, CONST int crname_len,
  double xs[],
  int istate[],
  int *ifail
);

extern void __stdcall E04NAF(
  CONST int *itmax,
  CONST int *msglvl,
  CONST int *n,
  CONST int *nclin,
  CONST int *nctotl,
  CONST int *nrowa,
  CONST int *nrowh,
  CONST int *ncolh,
  CONST double *bigbnd,
  CONST double a[] /* 2 dimension */,
  CONST double bl[],
  CONST double bu[],
  CONST double cvec[],
  CONST double featol[],
  CONST double hess[] /* 2 dimension */,
  void (__stdcall *qphess)(int *, int *, int *, int *, double[], double[], double[]),
  CONST int *cold,
  CONST int *lp,
  CONST int *orthog,
  double x[],
  int istate[],
  int *iter,
  double *obj,
  double clamda[],
  int iw[],
  CONST int *leniw,
  double w[],
  CONST int *lenw,
  int *ifail
);

extern void __stdcall E04NAN(
  CONST int *n,
  CONST int *nrowh,
  CONST int *ncolh,
  CONST int *jthcol,
  CONST double hess[] /* 2 dimension */,
  CONST double x[],
  CONST double hx[]
);

extern void __stdcall E04NCF(
  CONST int *mm,
  CONST int *n,
  CONST int *nclin,
  CONST int *lda,
  CONST int *ldr,
  CONST double a[] /* 2 dimension */,
  CONST double bl[],
  CONST double bu[],
  CONST double cvec[],
  int istate[],
  int kx[],
  double x[],
  double r[] /* 2 dimension */,
  double b[],
  int *iter,
  double *obj,
  double clamda[],
  int iw[],
  CONST int *leniw,
  double w[],
  CONST int *lenw,
  int *ifail
);

extern void __stdcall E04NDF(
  CONST int *ioptns,
  int *inform
);

extern void __stdcall E04NEF(
  CONST char *string, int string_len
);

extern void __stdcall E04NFF(
  CONST int *n,
  CONST int *nclin,
  CONST double a[] /* 2 dimension */,
  CONST int *lda,
  CONST double bl[],
  CONST double bu[],
  CONST double cvec[],
  CONST double h[] /* 2 dimension */,
  CONST int *ldh,
  void (__stdcall *qphess)(int *, int *, double[], int *, double[], double[]),
  int istate[],
  double x[],
  int *iter,
  double *obj,
  double ax[],
  double clamda[],
  int iw[],
  CONST int *leniw,
  double w[],
  CONST int *lenw,
  int *ifail
);

extern void __stdcall E04NGF(
  CONST int *ioptns,
  int *inform
);

extern void __stdcall E04NHF(
  CONST char *string, int string_len
);

extern void __stdcall E04NKF(
  CONST int *n,
  CONST int *m,
  CONST int *nnz,
  CONST int *iobj,
  CONST int *ncolh,
  void (__stdcall *qphx)(),
  double a[],
  CONST int ha[],
  CONST int ka[],
  double bl[],
  double bu[],
  CONST char *start, CONST int start_len,
  char *names, CONST int names_len,
  CONST int *nname,
  CONST char *crname, CONST int crname_len,
  int *ns,
  double xs[],
  int istate[],
  int *miniz,
  int *minz,
  int *ninf,
  double *sinf,
  double *obj,
  double clamda[],
  int iz[],
  CONST int *leniz,
  double z[],
  CONST int *lenz,
  int *ifail
);

extern void __stdcall E04NLF(
  CONST int *ioptns,
  int *inform
);

extern void __stdcall E04NMF(
  CONST char *string, int string_len
);

extern void __stdcall E04UCF(
  CONST int *n,
  CONST int *nclin,
  CONST int *ncnln,
  CONST int *lda,
  CONST int *ldcju,
  CONST int *ldr,
  double a[] /* 2 dimension */,
  CONST double bl[],
  CONST double bu[],
  void (__stdcall *confun)(int *, int *, int *, int *, int[], double[], double[], 
                 double[], int *, int[], double[]),
  void (__stdcall *objfun)(int *, int *, double[], double *, double[], int *, int[], 
                 double[]),
  int *iter,
  int istate[],
  double c[],
  double cjacu[] /* 2 dimension */,
  double clamda[],
  double *objf,
  double gradu[],
  double r[] /* 2 dimension */,
  double x[],
  int iw[],
  CONST int *leniw,
  double w[],
  CONST int *lenw,
  CONST int iuser[],
  CONST double user[],
  int *ifail
);

extern void __stdcall E04UDF(
  CONST int *ioptns,
  int *inform
);

extern void __stdcall E04UDM(
  CONST int *mode,
  CONST int *ncnln,
  CONST int *n,
  CONST int *nrowj,
  CONST int needc[],
  CONST double x[],
  CONST double c[],
  CONST double cjac[] /* 2 dimension */,
  CONST int *nstate,
  CONST int iuser[],
  CONST double user[]
);

extern void __stdcall E04UEF(
  CONST char *string, int string_len
);

extern void __stdcall E04UFF(
  int *irevcm,
  CONST int *n,
  CONST int *nclin,
  CONST int *ncnln,
  CONST int *lda,
  CONST int *ldcju,
  CONST int *ldr,
  double a[] /* 2 dimension */,
  CONST double bl[],
  CONST double bu[],
  int *iter,
  int istate[],
  double c[],
  double cjacu[] /* 2 dimension */,
  double clamda[],
  double *objf,
  double gradu[],
  double r[] /* 2 dimension */,
  double x[],
  int needc[],
  int iw[],
  CONST int *leniw,
  double w[],
  CONST int *lenw,
  int *ifail
);

extern void __stdcall E04UNF(
  CONST int *m,
  CONST int *n,
  CONST int *nclin,
  CONST int *ncnln,
  CONST int *lda,
  CONST int *ldcju,
  CONST int *ldfju,
  CONST int *ldr,
  double a[] /* 2 dimension */,
  CONST double bl[],
  CONST double bu[],
  CONST double y[],
  void (__stdcall *confun)(int *, int *, int *, int *, int[], double[], double[], 
                 double[], int *, int[], double[]),
  void (__stdcall *objfun)(int *, int *, int *, int *, double[], double[], double[], 
                 int *, int[], double[]),
  int *iter,
  int istate[],
  double c[],
  double cjacu[] /* 2 dimension */,
  double f[],
  double fjacu[] /* 2 dimension */,
  double clamda[],
  double *objf,
  double r[] /* 2 dimension */,
  double x[],
  int iw[],
  CONST int *leniw,
  double w[],
  CONST int *lenw,
  CONST int iuser[],
  CONST double user[],
  int *ifail
);

extern void __stdcall E04UPF(
  CONST int *m,
  CONST int *n,
  CONST int *nclin,
  CONST int *ncnln,
  CONST int *lda,
  CONST int *ldcju,
  CONST int *ldfju,
  CONST int *ldr,
  double a[] /* 2 dimension */,
  CONST double bl[],
  CONST double bu[],
  void (__stdcall *confun)(int *, int *, int *, int *, int[], double[], double[], 
                 double[], int *, int[], double[]),
  void (__stdcall *objfun)(int *, int *, int *, int *, double[], double[], double[], 
                 int *, int[], double[]),
  int *iter,
  int istate[],
  double c[],
  double cjacu[] /* 2 dimension */,
  double f[],
  double fjacu[] /* 2 dimension */,
  double clamda[],
  double *objf,
  double r[] /* 2 dimension */,
  double x[],
  int iw[],
  CONST int *leniw,
  double w[],
  CONST int *lenw,
  CONST int iuser[],
  CONST double user[],
  int *ifail
);

extern void __stdcall E04UQF(
  CONST int *ioptns,
  int *inform
);

extern void __stdcall E04URF(
  CONST char *string, int string_len
);

extern void __stdcall E04VCF(
  CONST int *itmax,
  CONST int *msglvl,
  CONST int *n,
  CONST int *nclin,
  CONST int *ncnln,
  CONST int *nctotl,
  CONST int *nrowa,
  CONST int *nrowj,
  CONST int *nrowr,
  CONST double *bigbnd,
  CONST double *epsaf,
  CONST double *eta,
  CONST double *ftol,
  double a[] /* 2 dimension */,
  CONST double bl[],
  CONST double bu[],
  CONST double featol[],
  void (__stdcall *confun)(int *, int *, int *, int *, double[], double[], double[], 
                 int *),
  void (__stdcall *objfun)(int *, int *, double[], double *, double[], int *),
  CONST int *cold,
  CONST int *fealin,
  CONST int *orthog,
  double x[],
  int istate[],
  double r[] /* 2 dimension */,
  int *iter,
  double c[],
  double cjac[] /* 2 dimension */,
  CONST double *objf,
  double objgrd[],
  double clamda[],
  int iw[],
  CONST int *leniw,
  double w[],
  CONST int *lenw,
  int *ifail
);

extern void __stdcall E04VDF(
  CONST int *itmax,
  CONST int *msglvl,
  CONST int *n,
  CONST int *nclin,
  CONST int *ncnln,
  CONST int *nctotl,
  CONST int *nrowa,
  CONST int *nrowj,
  CONST double *ctol,
  CONST double *ftol,
  double a[] /* 2 dimension */,
  CONST double bl[],
  CONST double bu[],
  void (__stdcall *confun)(int *, int *, int *, int *, double[], double[], double[], 
                 int *),
  void (__stdcall *objfun)(int *, int *, double[], double *, double[], int *),
  double x[],
  int istate[],
  double c[],
  double cjac[] /* 2 dimension */,
  CONST double *objf,
  double objgrd[],
  double clamda[],
  int iw[],
  CONST int *leniw,
  double w[],
  CONST int *lenw,
  int *ifail
);

extern void __stdcall E04VDM(
  CONST int *mode,
  CONST int *ncnln,
  CONST int *n,
  CONST int *nrowj,
  CONST double x[],
  CONST double c[],
  CONST double cjac[] /* 2 dimension */,
  CONST int *nstate
);

extern void __stdcall E04XAF(
  CONST int *msglvl,
  CONST int *n,
  CONST double *epsrf,
  double x[],
  CONST int *mode,
  void (__stdcall *objfun)(int *, int *, double[], double *, double[], int *, int[], 
                 double[]),
  CONST int *lhes,
  double hforw[],
  CONST double *fx,
  double grad[],
  double hcntrl[],
  double hesian[] /* 2 dimension */,
  int *iwarn,
  double work[],
  CONST int iuser[],
  CONST double user[],
  int info[],
  int *ifail
);

extern void __stdcall E04YAF(
  int *m,
  int *n,
  void (__stdcall *lsqfun)(int *, int *, int *, double[], double[], double[], int *, 
                 int[], int *, double[], int *),
  double x[],
  double fvec[],
  double fjac[] /* 2 dimension */,
  int *lj,
  CONST int iw[],
  CONST int *liw,
  double w[],
  CONST int *lw,
  int *ifail
);

extern void __stdcall E04YBF(
  int *m,
  int *n,
  void (__stdcall *lsqfun)(int *, int *, int *, double[], double[], double[], int *, 
                 int[], int *, double[], int *),
  void (__stdcall *lsqhes)(int *, int *, int *, double[], double[], double[], int *, 
                 int[], int *, double[], int *),
  double x[],
  double fvec[],
  double fjac[] /* 2 dimension */,
  int *lj,
  double b[],
  int *lb,
  CONST int iw[],
  CONST int *liw,
  double w[],
  CONST int *lw,
  int *ifail
);

extern void __stdcall E04YCF(
  CONST int *job,
  CONST int *m,
  CONST int *n,
  CONST double *fsumsq,
  CONST double s[],
  double v[] /* 2 dimension */,
  CONST int *lv,
  double cj[],
  double work[],
  int *ifail
);

extern void __stdcall E04ZCF(
  CONST int *n,
  CONST int *m,
  CONST int *la,
  void (__stdcall *con)(int *, int *, int *, int *, double[], double[], double[], int *),
  void (__stdcall *fun)(int *, int *, double[], double *, double[], int *),
  CONST double c[],
  CONST double a[] /* 2 dimension */,
  CONST double *f,
  CONST double g[],
  CONST double x[],
  double w[],
  CONST int *lw,
  int *ifail
);

extern void __stdcall F01AAF(
  double a[] /* 2 dimension */,
  CONST int *lda,
  CONST int *n,
  double x[] /* 2 dimension */,
  CONST int *ldx,
  double p[],
  int *ifail
);

extern void __stdcall F01ABF(
  double a[] /* 2 dimension */,
  CONST int *ia,
  CONST int *n,
  double b[] /* 2 dimension */,
  CONST int *ib,
  double z[],
  int *ifail
);

extern void __stdcall F01ACF(
  CONST int *n,
  CONST double *eps,
  double a[] /* 2 dimension */,
  CONST int *ia,
  double b[] /* 2 dimension */,
  CONST int *ib,
  double z[],
  int *l,
  int *ifail
);

extern void __stdcall F01ADF(
  CONST int *n,
  double a[] /* 2 dimension */,
  CONST int *ia,
  int *ifail
);

extern void __stdcall F01AEF(
  CONST int *n,
  double a[] /* 2 dimension */,
  CONST int *ia,
  double b[] /* 2 dimension */,
  CONST int *ib,
  double dl[],
  int *ifail
);

extern void __stdcall F01AFF(
  CONST int *n,
  CONST int *im1,
  CONST int *im2,
  double b[] /* 2 dimension */,
  CONST int *ib,
  double dl[],
  double z[] /* 2 dimension */,
  CONST int *iz
);

extern void __stdcall F01AGF(
  CONST int *n,
  CONST double *atol,
  double a[] /* 2 dimension */,
  CONST int *ia,
  double d[],
  double e[],
  double e2[]
);

extern void __stdcall F01AHF(
  CONST int *n,
  CONST int *im1,
  CONST int *im2,
  CONST double a[] /* 2 dimension */,
  CONST int *ia,
  CONST double e[],
  double z[] /* 2 dimension */,
  CONST int *iz
);

extern void __stdcall F01AJF(
  CONST int *n,
  CONST double *atol,
  CONST double a[] /* 2 dimension */,
  CONST int *ia,
  double d[],
  double e[],
  double z[] /* 2 dimension */,
  CONST int *iz
);

extern void __stdcall F01AKF(
  CONST int *n,
  CONST int *k,
  CONST int *l,
  double a[] /* 2 dimension */,
  CONST int *ia,
  int intger[]
);

extern void __stdcall F01ALF(
  CONST int *kl,
  CONST int *l,
  CONST int *ir,
  CONST double a[] /* 2 dimension */,
  CONST int *ia,
  CONST int intger[],
  double z[] /* 2 dimension */,
  CONST int *iz,
  CONST int *n
);

extern void __stdcall F01AMF(
  CONST int *n,
  CONST int *k,
  CONST int *l,
  double ar[] /* 2 dimension */,
  CONST int *iar,
  double ai[] /* 2 dimension */,
  CONST int *iai,
  int intger[]
);

extern void __stdcall F01ANF(
  CONST int *kl,
  CONST int *l,
  CONST int *ir,
  double ar[] /* 2 dimension */,
  CONST int *iar,
  double ai[] /* 2 dimension */,
  CONST int *iai,
  CONST int intger[],
  double zr[] /* 2 dimension */,
  CONST int *izr,
  double zi[] /* 2 dimension */,
  CONST int *izi,
  CONST int *n
);

extern void __stdcall F01APF(
  CONST int *n,
  CONST int *low,
  CONST int *iupp,
  CONST int intger[],
  CONST double h[] /* 2 dimension */,
  CONST int *ih,
  double v[] /* 2 dimension */,
  CONST int *iv
);

extern void __stdcall F01ATF(
  CONST int *n,
  CONST int *ib,
  double a[] /* 2 dimension */,
  CONST int *ia,
  int *low,
  int *lhi,
  double d[]
);

extern void __stdcall F01AUF(
  CONST int *n,
  CONST int *low,
  CONST int *lhi,
  CONST int *m,
  CONST double d[],
  double z[] /* 2 dimension */,
  CONST int *iz
);

extern void __stdcall F01AVF(
  CONST int *n,
  CONST int *ib,
  double ar[] /* 2 dimension */,
  CONST int *iar,
  double ai[] /* 2 dimension */,
  CONST int *iai,
  int *low,
  int *lhi,
  double d[]
);

extern void __stdcall F01AWF(
  CONST int *n,
  CONST int *low,
  CONST int *lhi,
  CONST int *m,
  CONST double d[],
  double zr[] /* 2 dimension */,
  CONST int *izr,
  double zi[] /* 2 dimension */,
  CONST int *izi
);

extern void __stdcall F01AXF(
  CONST int *m,
  CONST int *n,
  double qr[] /* 2 dimension */,
  CONST int *iqr,
  double alpha[],
  int ipivot[],
  double y[],
  double e[],
  int *ifail
);

extern void __stdcall F01AYF(
  CONST int *n,
  CONST double *tol,
  double a[],
  CONST int *ia,
  double d[],
  double e[],
  double e2[]
);

extern void __stdcall F01AZF(
  CONST int *n,
  CONST int *m1,
  CONST int *m2,
  CONST double a[],
  CONST int *ia,
  double z[] /* 2 dimension */,
  CONST int *izi
);

extern void __stdcall F01BCF(
  CONST int *n,
  CONST double *tol,
  double z[] /* 2 dimension */,
  CONST int *iz,
  double w[] /* 2 dimension */,
  CONST int *iw,
  double d[],
  double e[],
  double c[],
  double s[]
);

extern void __stdcall F01BDF(
  CONST int *n,
  double a[] /* 2 dimension */,
  CONST int *ia,
  double b[] /* 2 dimension */,
  CONST int *ib,
  double dl[],
  int *ifail
);

extern void __stdcall F01BEF(
  CONST int *n,
  CONST int *im1,
  CONST int *im2,
  CONST double b[] /* 2 dimension */,
  CONST int *ib,
  CONST double dl[],
  double z[] /* 2 dimension */,
  CONST int *iz
);

extern void __stdcall F01BKF(
  CONST int *m,
  CONST int *n,
  CONST double *t,
  double a[] /* 2 dimension */,
  CONST int *ia,
  double d[],
  double aijmax[],
  int inc[],
  int *irank,
  int *ifail
);

extern void __stdcall F01BLF(
  CONST int *m,
  CONST int *n,
  CONST double *t,
  double a[] /* 2 dimension */,
  CONST int *ia,
  double aijmax[],
  int *irank,
  int inc[],
  double d[],
  double u[] /* 2 dimension */,
  CONST int *iu,
  double du[],
  int *ifail
);

extern void __stdcall F01BNF(
  CONST int *n,
  Complex a[] /* 2 dimension */,
  CONST int *ia,
  double p[],
  int *ifail
);

extern void __stdcall F01BPF(
  CONST int *n,
  Complex a[] /* 2 dimension */,
  CONST int *ia,
  CONST Complex v[],
  int *ifail
);

extern void __stdcall F01BQF(
  CONST int *n,
  double *eps,
  double rl[],
  CONST int *irl,
  double d[],
  int *ifail
);

extern void __stdcall F01BRF(
  CONST int *n,
  CONST int *nz,
  double a[],
  CONST int *licn,
  int irn[],
  CONST int *lirn,
  int icn[],
  CONST double *u,
  int ikeep[] /* 2 dimension */,
  int iw[] /* 2 dimension */,
  double w[],
  CONST int *lblock,
  CONST int *grow,
  CONST int abort[],
  int idisp[],
  int *ifail
);

extern void __stdcall F01BSF(
  CONST int *n,
  CONST int *nz,
  double a[],
  CONST int *licn,
  CONST int ivect[],
  CONST int jvect[],
  int icn[],
  CONST int ikeep[] /* 2 dimension */,
  int iw[] /* 2 dimension */,
  double w[],
  CONST int *grow,
  CONST double *eps,
  double *rmin,
  CONST int *abort,
  CONST int idisp[],
  int *ifail
);

extern void __stdcall F01BTF(
  CONST int *n,
  double a[] /* 2 dimension */,
  CONST int *ia,
  double p[],
  double *dp,
  int *ifail
);

extern void __stdcall F01BUF(
  CONST int *n,
  CONST int *m1,
  CONST int *k,
  double a[] /* 2 dimension */,
  CONST int *ia,
  double w[],
  int *ifail
);

extern void __stdcall F01BVF(
  CONST int *n,
  CONST int *ma1,
  CONST int *mb1,
  CONST int *m3,
  CONST int *k,
  double a[] /* 2 dimension */,
  CONST int *ia,
  double b[] /* 2 dimension */,
  CONST int *ib,
  double c[] /* 2 dimension */,
  CONST int *ic,
  double w[],
  int *ifail
);

extern void __stdcall F01BWF(
  CONST int *n,
  CONST int *m1,
  double a[] /* 2 dimension */,
  CONST int *ia,
  double d[],
  double e[]
);

extern void __stdcall F01BXF(
  CONST int *n,
  double a[] /* 2 dimension */,
  CONST int *ia,
  double p[],
  int *ifail
);

extern void __stdcall F01CKF(
  double a[] /* 2 dimension */,
  double b[] /* 2 dimension */,
  double c[] /* 2 dimension */,
  CONST int *n,
  CONST int *p,
  CONST int *m,
  double z[],
  CONST int *iz,
  CONST int *opt,
  int *ifail
);

extern void __stdcall F01CLF(
  double a[] /* 2 dimension */,
  CONST double b[] /* 2 dimension */,
  CONST double c[] /* 2 dimension */,
  CONST int *n,
  CONST int *p,
  CONST int *m,
  int *ifail
);

extern void __stdcall F01CRF(
  double a[],
  CONST int *m,
  CONST int *n,
  CONST int *mn,
  int move[],
  CONST int *iwrk,
  int *ifail
);

extern void __stdcall F01CTF(
  CONST char *transa, CONST int transa_len,
  CONST char *transb, CONST int transb_len,
  CONST int *m,
  CONST int *n,
  CONST double *alpha,
  CONST double a[] /* 2 dimension */,
  CONST int *lda,
  CONST double *beta,
  CONST double b[] /* 2 dimension */,
  CONST int *ldb,
  double c[] /* 2 dimension */,
  CONST int *ldc,
  int *ifail
);

extern void __stdcall F01CWF(
  CONST char *transa, CONST int transa_len,
  CONST char *transb, CONST int transb_len,
  CONST int *m,
  CONST int *n,
  CONST Complex *alpha,
  CONST Complex a[] /* 2 dimension */,
  CONST int *lda,
  CONST Complex *beta,
  CONST Complex b[] /* 2 dimension */,
  CONST int *ldb,
  Complex c[] /* 2 dimension */,
  CONST int *ldc,
  int *ifail
);

extern void __stdcall F01LBF(
  CONST int *n,
  CONST int *m1,
  CONST int *m2,
  double a[] /* 2 dimension */,
  CONST int *ia,
  double al[] /* 2 dimension */,
  CONST int *il,
  int in[],
  int *iv,
  int *ifail
);

extern void __stdcall F01LEF(
  CONST int *n,
  double a[],
  CONST double *lambda,
  double b[],
  double c[],
  CONST double *tol,
  double d[],
  int in[],
  int *ifail
);

extern void __stdcall F01LHF(
  CONST int *n,
  CONST int *nbloks,
  CONST int blkstr[] /* 2 dimension */,
  double a[],
  CONST int *lena,
  int pivot[],
  double *tol,
  int *index,
  int *ifail
);

extern void __stdcall F01LZF(
  CONST int *n,
  CONST double a[] /* 2 dimension */,
  CONST int *nra,
  double c[] /* 2 dimension */,
  CONST int *nrc,
  CONST int *wantb,
  double b[],
  CONST int *wantq,
  CONST int *wanty,
  double y[] /* 2 dimension */,
  CONST int *nry,
  CONST int *ly,
  CONST int *wantz,
  double z[] /* 2 dimension */,
  CONST int *nrz,
  CONST int *ncz,
  double d[],
  double e[],
  double work1[],
  double work2[],
  int *ifail
);

extern void __stdcall F01MAF(
  CONST int *n,
  CONST int *nz,
  double a[],
  CONST int *la,
  int ini[],
  CONST int *lini,
  int inj[],
  double *droptl,
  CONST double *densw,
  double w[] /* 2 dimension */,
  int ik[] /* 2 dimension */,
  int iw[] /* 2 dimension */,
  CONST int abort[],
  int inform[],
  int *iflag
);

extern void __stdcall F01MCF(
  CONST int *n,
  CONST double a[],
  CONST int *lal,
  int nrow[],
  double l[],
  double d[],
  int *ifail
);

extern void __stdcall F01NAF(
  CONST int *n,
  CONST int *ml,
  CONST int *mu,
  Complex a[] /* 2 dimension */,
  CONST int *nra,
  CONST double *tol,
  int in[],
  double scale[],
  int *ifail
);

extern void __stdcall F01QAF(
  CONST int *m,
  CONST int *n,
  CONST double a[] /* 2 dimension */,
  CONST int *nra,
  double c[] /* 2 dimension */,
  CONST int *nrc,
  double z[],
  int *ifail
);

extern void __stdcall F01QBF(
  CONST int *m,
  CONST int *n,
  CONST double a[] /* 2 dimension */,
  CONST int *nra,
  double c[] /* 2 dimension */,
  CONST int *nrc,
  double work[],
  int *ifail
);

extern void __stdcall F01QCF(
  CONST int *m,
  CONST int *n,
  double a[] /* 2 dimension */,
  CONST int *lda,
  double zeta[],
  int *ifail
);

extern void __stdcall F01QDF(
  CONST char *trans, CONST int trans_len,
  CONST char *wheret, CONST int wheret_len,
  CONST int *m,
  CONST int *n,
  double a[] /* 2 dimension */,
  CONST int *lda,
  CONST double zeta[],
  CONST int *ncolb,
  double b[] /* 2 dimension */,
  CONST int *ldb,
  double work[],
  int *ifail
);

extern void __stdcall F01QEF(
  CONST char *wheret, CONST int wheret_len,
  CONST int *m,
  CONST int *n,
  CONST int *ncolq,
  double a[] /* 2 dimension */,
  CONST int *lda,
  CONST double zeta[],
  double work[],
  int *ifail
);

extern void __stdcall F01QFF(
  CONST char *pivot, CONST int pivot_len,
  CONST int *m,
  CONST int *n,
  double a[] /* 2 dimension */,
  CONST int *lda,
  double zeta[],
  int perm[],
  double work[],
  int *ifail
);

extern void __stdcall F01QGF(
  CONST int *m,
  CONST int *n,
  double a[] /* 2 dimension */,
  CONST int *lda,
  double zeta[],
  int *ifail
);

extern void __stdcall F01QJF(
  CONST int *m,
  CONST int *n,
  double a[] /* 2 dimension */,
  CONST int *lda,
  double zeta[],
  int *ifail
);

extern void __stdcall F01QKF(
  CONST char *wheret, CONST int wheret_len,
  CONST int *m,
  CONST int *n,
  CONST int *nrowp,
  double a[] /* 2 dimension */,
  CONST int *lda,
  CONST double zeta[],
  double work[],
  int *ifail
);

extern void __stdcall F01RCF(
  CONST int *m,
  CONST int *n,
  Complex a[] /* 2 dimension */,
  CONST int *lda,
  Complex theta[],
  int *ifail
);

extern void __stdcall F01RDF(
  CONST char *trans, CONST int trans_len,
  CONST char *wheret, CONST int wheret_len,
  CONST int *m,
  CONST int *n,
  Complex a[] /* 2 dimension */,
  CONST int *lda,
  CONST Complex theta[],
  CONST int *ncolb,
  Complex b[] /* 2 dimension */,
  CONST int *ldb,
  Complex work[],
  int *ifail
);

extern void __stdcall F01REF(
  CONST char *wheret, CONST int wheret_len,
  CONST int *m,
  CONST int *n,
  CONST int *ncolq,
  Complex a[] /* 2 dimension */,
  CONST int *lda,
  CONST Complex theta[],
  Complex work[],
  int *ifail
);

extern void __stdcall F01RFF(
  CONST char *pivot, CONST int pivot_len,
  CONST int *m,
  CONST int *n,
  Complex a[] /* 2 dimension */,
  CONST int *lda,
  Complex theta[],
  int perm[],
  double work[],
  int *ifail
);

extern void __stdcall F01RGF(
  CONST int *m,
  CONST int *n,
  Complex a[] /* 2 dimension */,
  CONST int *lda,
  Complex theta[],
  int *ifail
);

extern void __stdcall F01RJF(
  CONST int *m,
  CONST int *n,
  Complex a[] /* 2 dimension */,
  CONST int *lda,
  Complex theta[],
  int *ifail
);

extern void __stdcall F01RKF(
  CONST char *wheret, CONST int wheret_len,
  CONST int *m,
  CONST int *n,
  CONST int *nrowp,
  Complex a[] /* 2 dimension */,
  CONST int *lda,
  CONST Complex theta[],
  Complex work[],
  int *ifail
);

extern void __stdcall F01ZAF(
  CONST char *job, CONST int job_len,
  CONST char *uplo, CONST int uplo_len,
  CONST char *diag, CONST int diag_len,
  CONST int *n,
  double a[] /* 2 dimension */,
  CONST int *lda,
  double b[],
  int *ifail
);

extern void __stdcall F01ZBF(
  CONST char *job, CONST int job_len,
  CONST char *uplo, CONST int uplo_len,
  CONST char *diag, CONST int diag_len,
  CONST int *n,
  Complex a[] /* 2 dimension */,
  CONST int *lda,
  Complex b[],
  int *ifail
);

extern void __stdcall F01ZCF(
  CONST char *job, CONST int job_len,
  CONST int *m,
  CONST int *n,
  CONST int *kl,
  CONST int *ku,
  double a[] /* 2 dimension */,
  CONST int *lda,
  double b[] /* 2 dimension */,
  CONST int *ldb,
  int *ifail
);

extern void __stdcall F01ZDF(
  CONST char *job, CONST int job_len,
  CONST int *m,
  CONST int *n,
  CONST int *kl,
  CONST int *ku,
  Complex a[] /* 2 dimension */,
  CONST int *lda,
  Complex b[] /* 2 dimension */,
  CONST int *ldb,
  int *ifail
);

extern void __stdcall F02AAF(
  double a[] /* 2 dimension */,
  CONST int *ia,
  CONST int *n,
  double r[],
  double e[],
  int *ifail
);

extern void __stdcall F02ABF(
  CONST double a[] /* 2 dimension */,
  CONST int *ia,
  CONST int *n,
  double r[],
  double v[] /* 2 dimension */,
  CONST int *iv,
  double e[],
  int *ifail
);

extern void __stdcall F02ADF(
  double a[] /* 2 dimension */,
  CONST int *ia,
  double b[] /* 2 dimension */,
  CONST int *ib,
  CONST int *n,
  double r[],
  double de[],
  int *ifail
);

extern void __stdcall F02AEF(
  double a[] /* 2 dimension */,
  CONST int *ia,
  double b[] /* 2 dimension */,
  CONST int *ib,
  CONST int *n,
  double r[],
  double v[] /* 2 dimension */,
  CONST int *iv,
  double dl[],
  double e[],
  int *ifail
);

extern void __stdcall F02AFF(
  double a[] /* 2 dimension */,
  CONST int *ia,
  CONST int *n,
  double rr[],
  double ri[],
  int intger[],
  int *lfail
);

extern void __stdcall F02AGF(
  double a[] /* 2 dimension */,
  CONST int *ia,
  CONST int *n,
  double rr[],
  double ri[],
  double vr[] /* 2 dimension */,
  CONST int *ivr,
  double vi[] /* 2 dimension */,
  CONST int *ivi,
  int intger[],
  int *ifail
);

extern void __stdcall F02AJF(
  double ar[] /* 2 dimension */,
  CONST int *iar,
  double ai[] /* 2 dimension */,
  CONST int *iai,
  CONST int *n,
  double rr[],
  double ri[],
  int intger[],
  int *lfail
);

extern void __stdcall F02AKF(
  double ar[] /* 2 dimension */,
  CONST int *iar,
  double ai[] /* 2 dimension */,
  CONST int *iai,
  CONST int *n,
  double wr[],
  double wi[],
  double vr[] /* 2 dimension */,
  CONST int *ivr,
  double vi[] /* 2 dimension */,
  CONST int *ivi,
  int intger[],
  int *ifail
);

extern void __stdcall F02AMF(
  CONST int *n,
  CONST double *eps,
  double d[],
  double e[],
  double z[] /* 2 dimension */,
  CONST int *iz,
  int *ifail
);

extern void __stdcall F02ANF(
  CONST int *nn,
  CONST double *acc,
  double hr[] /* 2 dimension */,
  CONST int *ihr,
  double hi[] /* 2 dimension */,
  CONST int *ihi,
  double wr[],
  double wi[],
  int *lfail
);

extern void __stdcall F02APF(
  CONST int *nn,
  CONST double *acc,
  double h[] /* 2 dimension */,
  CONST int *ih,
  double wr[],
  double wi[],
  int icnt[],
  int *lfail
);

extern void __stdcall F02AQF(
  CONST int *n,
  CONST int *low,
  CONST int *upp,
  CONST double *machep,
  double h[] /* 2 dimension */,
  CONST int *ih,
  double vecs[] /* 2 dimension */,
  CONST int *ivecs,
  double wr[],
  double wi[],
  int cnt[],
  int *ifail
);

extern void __stdcall F02ARF(
  CONST int *n,
  CONST int *low,
  CONST int *iupp,
  CONST double *acheps,
  CONST int intger[],
  double hr[] /* 2 dimension */,
  CONST int *ihr,
  double hi[] /* 2 dimension */,
  CONST int *ihi,
  double wr[],
  double wi[],
  double vr[] /* 2 dimension */,
  CONST int *ivr,
  double vi[] /* 2 dimension */,
  CONST int *ivi,
  int *ifail
);

extern void __stdcall F02AVF(
  CONST int *n,
  CONST double *acheps,
  double d[],
  double e[],
  int *ifail
);

extern void __stdcall F02AWF(
  double ar[] /* 2 dimension */,
  CONST int *iar,
  double ai[] /* 2 dimension */,
  CONST int *iai,
  CONST int *n,
  double wr[],
  double wk1[],
  double wk2[],
  double wk3[],
  int *ifail
);

extern void __stdcall F02AXF(
  CONST double ar[] /* 2 dimension */,
  CONST int *iar,
  CONST double ai[] /* 2 dimension */,
  CONST int *iai,
  CONST int *n,
  double wr[],
  double vr[] /* 2 dimension */,
  CONST int *ivr,
  double vi[] /* 2 dimension */,
  CONST int *ivi,
  double wk1[],
  double wk2[],
  double wk3[],
  int *ifail
);

extern void __stdcall F02AYF(
  CONST int *n,
  CONST double *eps,
  double d[],
  double e[],
  double z[] /* 2 dimension */,
  CONST int *iz,
  double w[] /* 2 dimension */,
  CONST int *iw,
  int *ifail
);

extern void __stdcall F02BBF(
  double a[] /* 2 dimension */,
  CONST int *ia,
  CONST int *n,
  CONST double *alb,
  CONST double *ub,
  CONST int *m,
  int *mm,
  double r[],
  double v[] /* 2 dimension */,
  CONST int *iv,
  double d[],
  double e[],
  double e2[],
  double x[] /* 2 dimension */,
  double g[],
  int c[],
  int icount[],
  int *ifail
);

extern void __stdcall F02BCF(
  double a[] /* 2 dimension */,
  CONST int *ia,
  CONST int *n,
  CONST double *rlb,
  CONST double *rub,
  CONST int *m,
  int *mm,
  double rr[],
  double ri[],
  double vr[] /* 2 dimension */,
  CONST int *ivr,
  double vi[] /* 2 dimension */,
  CONST int *ivi,
  int intger[],
  int icnt[],
  int c[],
  double b[] /* 2 dimension */,
  CONST int *ib,
  double u[],
  double v[],
  int *lfail
);

extern void __stdcall F02BDF(
  double ar[] /* 2 dimension */,
  CONST int *iar,
  double ai[] /* 2 dimension */,
  CONST int *iai,
  CONST int *n,
  CONST double *rlb,
  CONST double *rub,
  CONST int *m,
  int *mm,
  double rr[],
  double ri[],
  double vr[] /* 2 dimension */,
  CONST int *ivr,
  double vi[] /* 2 dimension */,
  CONST int *ivi,
  int intger[],
  int c[],
  double br[] /* 2 dimension */,
  CONST int *ibr,
  double bi[] /* 2 dimension */,
  CONST int *ibi,
  double u[],
  double v[],
  int *lfail
);

extern void __stdcall F02BEF(
  CONST int *n,
  CONST double c[],
  CONST double *alb,
  CONST double *ub,
  CONST double *acheps,
  CONST double *eps,
  CONST double b[],
  double beta[],
  CONST int *m,
  int *mm,
  double root[],
  double vec[] /* 2 dimension */,
  CONST int *ivec,
  int icount[],
  double x[] /* 2 dimension */,
  int log[],
  int *ifail
);

extern void __stdcall F02BFF(
  CONST double c[],
  double b[],
  double beta[],
  CONST int *n,
  CONST int *m1,
  CONST int *m2,
  CONST int *mm12,
  double *eps1,
  CONST double *relfeh,
  double *eps2,
  int *iz,
  double x[],
  double wu[]
);

extern void __stdcall F02BJF(
  CONST int *n,
  double a[] /* 2 dimension */,
  CONST int *ia,
  double b[] /* 2 dimension */,
  CONST int *ib,
  CONST double *eps1,
  double alfr[],
  double alfi[],
  double beta[],
  CONST int *matz,
  double z[] /* 2 dimension */,
  CONST int *iz,
  int iter[],
  int *ifail
);

extern void __stdcall F02BKF(
  CONST int *n,
  CONST int *mm,
  CONST double a[] /* 2 dimension */,
  CONST int *ia,
  CONST double wi[],
  CONST int c[],
  double wr[],
  double z[] /* 2 dimension */,
  CONST int *iz,
  double b[] /* 2 dimension */,
  CONST int *ib,
  double u[],
  double v[],
  int *ifail
);

extern void __stdcall F02BLF(
  CONST int *n,
  CONST int *mm,
  CONST double ar[] /* 2 dimension */,
  CONST int *iar,
  CONST double ai[] /* 2 dimension */,
  CONST int *iai,
  CONST double wi[],
  CONST int c[],
  double wr[],
  double zr[] /* 2 dimension */,
  CONST int *izr,
  double zi[] /* 2 dimension */,
  CONST int *izi,
  double br[] /* 2 dimension */,
  CONST int *ibr,
  double bi[] /* 2 dimension */,
  CONST int *ibi,
  double u[],
  double v[],
  int *ifail
);

extern void __stdcall F02EAF(
  CONST char *job, CONST int job_len,
  CONST int *n,
  double a[] /* 2 dimension */,
  CONST int *lda,
  double wr[],
  double wi[],
  double z[] /* 2 dimension */,
  CONST int *ldz,
  double work[],
  CONST int *lwork,
  int *ifail
);

extern void __stdcall F02EBF(
  CONST char *job, CONST int job_len,
  CONST int *n,
  double a[] /* 2 dimension */,
  CONST int *lda,
  double wr[],
  double wi[],
  double vr[] /* 2 dimension */,
  CONST int *ldvr,
  double vi[] /* 2 dimension */,
  CONST int *ldvi,
  double work[],
  CONST int *lwork,
  int *ifail
);

extern void __stdcall F02ECF(
  CONST char *crit, CONST int crit_len,
  CONST int *n,
  double a[] /* 2 dimension */,
  CONST int *lda,
  CONST double *wl,
  CONST double *wu,
  CONST int *mest,
  int *m,
  double wr[],
  double wi[],
  double vr[] /* 2 dimension */,
  CONST int *ldvr,
  double vi[] /* 2 dimension */,
  CONST int *ldvi,
  double work[],
  CONST int *lwork,
  int iwork[],
  int bwork[],
  int *ifail
);

extern void __stdcall F02FAF(
  CONST char *job, CONST int job_len,
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  double a[] /* 2 dimension */,
  CONST int *lda,
  double w[],
  double work[],
  CONST int *lwork,
  int *ifail
);

extern void __stdcall F02FCF(
  CONST char *job, CONST int job_len,
  CONST char *range, CONST int range_len,
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  double a[] /* 2 dimension */,
  CONST int *lda,
  CONST double *wl,
  CONST double *wu,
  CONST int *il,
  CONST int *iu,
  CONST int *mest,
  int *m,
  double w[],
  double z[] /* 2 dimension */,
  CONST int *ldz,
  double work[],
  CONST int *lwork,
  int iwork[],
  int *ifail
);

extern void __stdcall F02FDF(
  CONST int *itype,
  CONST char *job, CONST int job_len,
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  double a[] /* 2 dimension */,
  CONST int *lda,
  double b[] /* 2 dimension */,
  CONST int *ldb,
  double w[],
  double work[],
  CONST int *lwork,
  int *ifail
);

extern void __stdcall F02FHF(
  CONST int *n,
  CONST int *ma,
  double a[] /* 2 dimension */,
  CONST int *nra,
  CONST int *mb,
  double b[] /* 2 dimension */,
  CONST int *nrb,
  double d[],
  double work[],
  CONST int *lwork,
  int *ifail
);

extern void __stdcall F02FJF(
  CONST int *n,
  int *em,
  CONST int *p,
  int *km,
  CONST double *eps,
  double (__stdcall *ip)(int *, int *, double[], double[], double[], int *, int[], int *),
  void (__stdcall *op)(int *, int *, double[], double[], double[], int *, int[], int *),
  void (__stdcall *inf)(int *, int *, int *, int *, int *, double[], double[]),
  CONST int *novecs,
  double x[] /* 2 dimension */,
  CONST int *mn,
  double d[],
  double work[],
  CONST int *lwork,
  CONST double rwork[],
  CONST int *lrwork,
  CONST int iwork[],
  CONST int *liwork,
  int *ifail
);

extern void __stdcall F02FJZ(
  CONST int *istate,
  CONST int *nextit,
  CONST int *nevals,
  CONST int *nevecs,
  CONST int *k,
  CONST double f[],
  CONST double d[]
);

extern void __stdcall F02GAF(
  CONST char *job, CONST int job_len,
  CONST int *n,
  Complex a[] /* 2 dimension */,
  CONST int *lda,
  Complex w[],
  Complex z[] /* 2 dimension */,
  CONST int *ldz,
  double rwork[],
  Complex work[],
  CONST int *lwork,
  int *ifail
);

extern void __stdcall F02GBF(
  CONST char *job, CONST int job_len,
  CONST int *n,
  Complex a[] /* 2 dimension */,
  CONST int *lda,
  Complex w[],
  Complex v[] /* 2 dimension */,
  CONST int *ldv,
  double rwork[],
  Complex work[],
  CONST int *lwork,
  int *ifail
);

extern void __stdcall F02GCF(
  CONST char *crit, CONST int crit_len,
  CONST int *n,
  Complex a[] /* 2 dimension */,
  CONST int *lda,
  CONST double *wl,
  CONST double *wu,
  CONST int *mest,
  int *m,
  Complex w[],
  Complex v[] /* 2 dimension */,
  CONST int *ldv,
  Complex work[],
  CONST int *lwork,
  double rwork[],
  int iwork[],
  int bwork[],
  int *ifail
);

extern void __stdcall F02GJF(
  CONST int *n,
  double ar[] /* 2 dimension */,
  CONST int *iar,
  double ai[] /* 2 dimension */,
  CONST int *iai,
  double br[] /* 2 dimension */,
  CONST int *ibr,
  double bi[] /* 2 dimension */,
  CONST int *ibi,
  CONST double *eps1,
  double alfr[],
  double alfi[],
  double beta[],
  CONST int *matz,
  double zr[] /* 2 dimension */,
  CONST int *izr,
  double zi[] /* 2 dimension */,
  CONST int *izi,
  int iter[],
  int *ifail
);

extern void __stdcall F02HAF(
  CONST char *job, CONST int job_len,
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  Complex a[] /* 2 dimension */,
  CONST int *lda,
  double w[],
  double rwork[],
  Complex work[],
  CONST int *lwork,
  int *ifail
);

extern void __stdcall F02HCF(
  CONST char *job, CONST int job_len,
  CONST char *range, CONST int range_len,
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  Complex a[] /* 2 dimension */,
  CONST int *lda,
  CONST double *wl,
  CONST double *wu,
  CONST int *il,
  CONST int *iu,
  CONST int *mest,
  int *m,
  double w[],
  Complex z[] /* 2 dimension */,
  CONST int *ldz,
  Complex work[],
  CONST int *lwork,
  double rwork[],
  int iwork[],
  int *ifail
);

extern void __stdcall F02HDF(
  CONST int *itype,
  CONST char *job, CONST int job_len,
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  Complex a[] /* 2 dimension */,
  CONST int *lda,
  Complex b[] /* 2 dimension */,
  CONST int *ldb,
  double w[],
  double rwork[],
  Complex work[],
  CONST int *lwork,
  int *ifail
);

extern void __stdcall F02SDF(
  CONST int *n,
  CONST int *ma1,
  CONST int *mb1,
  double a[] /* 2 dimension */,
  CONST int *ia,
  double b[] /* 2 dimension */,
  CONST int *ib,
  CONST int *sym,
  CONST double *relep,
  CONST double *rmu,
  double vec[],
  double d[],
  int _int[],
  double work[],
  CONST int *lwork,
  int *ifail
);

extern void __stdcall F02SWF(
  CONST int *n,
  double a[] /* 2 dimension */,
  CONST int *lda,
  double d[],
  double e[],
  CONST int *ncoly,
  double y[] /* 2 dimension */,
  CONST int *ldy,
  CONST int *wantq,
  double q[] /* 2 dimension */,
  CONST int *ldq,
  int *ifail
);

extern void __stdcall F02SXF(
  CONST int *n,
  double a[] /* 2 dimension */,
  CONST int *lda,
  CONST int *ncoly,
  double y[] /* 2 dimension */,
  CONST int *ldy,
  double work[],
  int *ifail
);

extern void __stdcall F02SYF(
  CONST int *n,
  double d[],
  double e[],
  CONST int *ncolb,
  double b[] /* 2 dimension */,
  CONST int *ldb,
  CONST int *nrowy,
  double y[] /* 2 dimension */,
  CONST int *ldy,
  CONST int *ncolz,
  double z[] /* 2 dimension */,
  CONST int *ldz,
  double work[],
  int *ifail
);

extern void __stdcall F02SZF(
  CONST int *n,
  CONST double d[],
  CONST double e[],
  double sv[],
  CONST int *wantb,
  double b[],
  CONST int *wanty,
  double y[] /* 2 dimension */,
  CONST int *nry,
  CONST int *ly,
  CONST int *wantz,
  double z[] /* 2 dimension */,
  CONST int *nrz,
  CONST int *ncz,
  double work1[],
  double work2[],
  double work3[],
  int *ifail
);

extern void __stdcall F02UWF(
  CONST int *n,
  Complex a[] /* 2 dimension */,
  CONST int *lda,
  double d[],
  double e[],
  CONST int *ncoly,
  Complex y[] /* 2 dimension */,
  CONST int *ldy,
  CONST int *wantq,
  Complex q[] /* 2 dimension */,
  CONST int *ldq,
  Complex work[],
  int *ifail
);

extern void __stdcall F02UXF(
  CONST int *n,
  Complex a[] /* 2 dimension */,
  CONST int *lda,
  CONST int *ncoly,
  Complex y[] /* 2 dimension */,
  CONST int *ldy,
  double rwork[],
  Complex cwork[],
  int *ifail
);

extern void __stdcall F02UYF(
  CONST int *n,
  double d[],
  double e[],
  CONST int *ncolb,
  Complex b[] /* 2 dimension */,
  CONST int *ldb,
  CONST int *nrowy,
  Complex y[] /* 2 dimension */,
  CONST int *ldy,
  CONST int *ncolz,
  Complex z[] /* 2 dimension */,
  CONST int *ldz,
  double work[],
  int *ifail
);

extern void __stdcall F02WAF(
  CONST int *m,
  CONST int *n,
  double a[] /* 2 dimension */,
  CONST int *nra,
  CONST int *wantb,
  double b[],
  double sv[],
  double work[],
  CONST int *lwork,
  int *ifail
);

extern void __stdcall F02WBF(
  CONST int *m,
  CONST int *n,
  double a[] /* 2 dimension */,
  CONST int *nra,
  CONST int *wantb,
  double b[],
  double sv[],
  double work[],
  CONST int *lwork,
  int *ifail
);

extern void __stdcall F02WCF(
  CONST int *m,
  CONST int *n,
  CONST int *minmn,
  CONST double a[] /* 2 dimension */,
  CONST int *nra,
  double q[] /* 2 dimension */,
  CONST int *nrq,
  double sv[],
  double pt[] /* 2 dimension */,
  CONST int *nrpt,
  double work[],
  CONST int *lwork,
  int *ifail
);

extern void __stdcall F02WDF(
  CONST int *m,
  CONST int *n,
  double a[] /* 2 dimension */,
  CONST int *nra,
  CONST int *wantb,
  double b[],
  CONST double *tol,
  int *svd,
  int *irank,
  double z[],
  double sv[],
  CONST int *wantr,
  double r[] /* 2 dimension */,
  CONST int *nrr,
  CONST int *wantpt,
  double pt[] /* 2 dimension */,
  CONST int *nrpt,
  double work[],
  CONST int *lwork,
  int *ifail
);

extern void __stdcall F02WEF(
  CONST int *m,
  CONST int *n,
  double a[] /* 2 dimension */,
  CONST int *lda,
  CONST int *ncolb,
  double b[] /* 2 dimension */,
  CONST int *ldb,
  CONST int *wantq,
  double q[] /* 2 dimension */,
  CONST int *ldq,
  double sv[],
  CONST int *wantp,
  double pt[] /* 2 dimension */,
  CONST int *ldpt,
  double work[],
  int *ifail
);

extern void __stdcall F02WUF(
  CONST int *n,
  double a[] /* 2 dimension */,
  CONST int *lda,
  CONST int *ncolb,
  double b[] /* 2 dimension */,
  CONST int *ldb,
  CONST int *wantq,
  double q[] /* 2 dimension */,
  CONST int *ldq,
  double sv[],
  CONST int *wantp,
  double work[],
  int *ifail
);

extern void __stdcall F02XEF(
  CONST int *m,
  CONST int *n,
  Complex a[] /* 2 dimension */,
  CONST int *lda,
  CONST int *ncolb,
  Complex b[] /* 2 dimension */,
  CONST int *ldb,
  CONST int *wantq,
  Complex q[] /* 2 dimension */,
  CONST int *ldq,
  double sv[],
  CONST int *wantp,
  Complex ph[] /* 2 dimension */,
  CONST int *ldph,
  double rwork[],
  Complex cwork[],
  int *ifail
);

extern void __stdcall F02XUF(
  CONST int *n,
  Complex a[] /* 2 dimension */,
  CONST int *lda,
  CONST int *ncolb,
  Complex b[] /* 2 dimension */,
  CONST int *ldb,
  CONST int *wantq,
  Complex q[] /* 2 dimension */,
  CONST int *ldq,
  double sv[],
  CONST int *wantp,
  double rwork[],
  Complex cwork[],
  int *ifail
);

extern void __stdcall F03AAF(
  double a[] /* 2 dimension */,
  CONST int *ia,
  CONST int *n,
  double *det,
  double wkspce[],
  int *ifail
);

extern void __stdcall F03ABF(
  double a[] /* 2 dimension */,
  CONST int *ia,
  CONST int *n,
  double *det,
  double wkspce[],
  int *ifail
);

extern void __stdcall F03ACF(
  CONST double a[] /* 2 dimension */,
  CONST int *ia,
  CONST int *n,
  CONST int *m,
  double *det,
  double rl[] /* 2 dimension */,
  CONST int *il,
  CONST int *m1,
  int *lfail
);

extern void __stdcall F03ADF(
  Complex a[] /* 2 dimension */,
  CONST int *ia,
  CONST int *n,
  double *detr,
  double *deti,
  double wkspce[],
  int *ifail
);

extern void __stdcall F03AEF(
  CONST int *n,
  double a[] /* 2 dimension */,
  CONST int *ia,
  double p[],
  double *d1,
  int *id,
  int *ifail
);

extern void __stdcall F03AFF(
  CONST int *n,
  CONST double *eps,
  double a[] /* 2 dimension */,
  CONST int *ia,
  double *d1,
  int *id,
  double p[],
  int *ifail
);

extern void __stdcall F03AGF(
  CONST int *n,
  CONST int *m,
  CONST double a[] /* 2 dimension */,
  CONST int *ia,
  double rl[] /* 2 dimension */,
  CONST int *il,
  CONST int *mm,
  double *d1,
  int *id,
  int *lfail
);

extern void __stdcall F03AHF(
  CONST int *n,
  Complex a[] /* 2 dimension */,
  CONST int *ia,
  double *detr,
  double *deti,
  int *id,
  double rint[],
  int *ifail
);

extern void __stdcall F03AMF(
  CONST int *n,
  CONST int *ten,
  CONST double p[],
  double *detm,
  int *idete
);

extern void __stdcall F04AAF(
  double a[] /* 2 dimension */,
  CONST int *ia,
  CONST double b[] /* 2 dimension */,
  CONST int *ib,
  CONST int *n,
  CONST int *m,
  double c[] /* 2 dimension */,
  CONST int *ic,
  double wkspce[],
  int *ifail
);

extern void __stdcall F04ABF(
  double a[] /* 2 dimension */,
  CONST int *ia,
  CONST double b[] /* 2 dimension */,
  CONST int *ib,
  CONST int *n,
  CONST int *m,
  double c[] /* 2 dimension */,
  CONST int *ic,
  double wkspce[],
  double bb[] /* 2 dimension */,
  CONST int *ibb,
  int *ifail
);

extern void __stdcall F04ACF(
  CONST double a[] /* 2 dimension */,
  CONST int *ia,
  CONST double b[] /* 2 dimension */,
  CONST int *ib,
  CONST int *n,
  CONST int *m,
  CONST int *lr,
  double c[] /* 2 dimension */,
  CONST int *ic,
  double rl[] /* 2 dimension */,
  CONST int *il,
  CONST int *m1,
  int *lfail
);

extern void __stdcall F04ADF(
  Complex a[] /* 2 dimension */,
  CONST int *ia,
  CONST Complex b[] /* 2 dimension */,
  CONST int *ib,
  CONST int *n,
  CONST int *m,
  Complex c[] /* 2 dimension */,
  CONST int *ic,
  double wkspce[],
  int *ifail
);

extern void __stdcall F04AEF(
  CONST double a[] /* 2 dimension */,
  CONST int *ia,
  CONST double b[] /* 2 dimension */,
  CONST int *ib,
  CONST int *n,
  CONST int *m,
  double c[] /* 2 dimension */,
  CONST int *ic,
  double wkspce[],
  double aa[] /* 2 dimension */,
  CONST int *iaa,
  double bb[] /* 2 dimension */,
  CONST int *ibb,
  int *ifail
);

extern void __stdcall F04AFF(
  CONST int *n,
  CONST int *ir,
  double a[] /* 2 dimension */,
  CONST int *ia,
  double p[],
  CONST double b[] /* 2 dimension */,
  CONST int *ib,
  CONST double *eps,
  double x[] /* 2 dimension */,
  CONST int *ix,
  double bb[] /* 2 dimension */,
  CONST int *ibb,
  int *l,
  int *ifail
);

extern void __stdcall F04AGF(
  CONST int *n,
  CONST int *ir,
  double a[] /* 2 dimension */,
  CONST int *ia,
  double p[],
  CONST double b[] /* 2 dimension */,
  CONST int *ib,
  double x[] /* 2 dimension */,
  CONST int *ix
);

extern void __stdcall F04AHF(
  CONST int *n,
  CONST int *ir,
  CONST double a[] /* 2 dimension */,
  CONST int *ia,
  CONST double aa[] /* 2 dimension */,
  CONST int *iaa,
  CONST double p[],
  CONST double b[] /* 2 dimension */,
  CONST int *ib,
  CONST double *eps,
  double x[] /* 2 dimension */,
  CONST int *ix,
  double bb[] /* 2 dimension */,
  CONST int *ibb,
  int *l,
  int *ifail
);

extern void __stdcall F04AJF(
  CONST int *n,
  CONST int *ir,
  CONST double a[] /* 2 dimension */,
  CONST int *ia,
  CONST double p[],
  double b[] /* 2 dimension */,
  CONST int *ib
);

extern void __stdcall F04AKF(
  CONST int *n,
  CONST int *ir,
  CONST Complex a[] /* 2 dimension */,
  CONST int *ia,
  CONST double rint[],
  Complex b[] /* 2 dimension */,
  CONST int *ib
);

extern void __stdcall F04ALF(
  CONST int *n,
  CONST int *m,
  CONST int *ir,
  CONST double rl[] /* 2 dimension */,
  CONST int *il,
  CONST int *m1,
  CONST double b[] /* 2 dimension */,
  CONST int *ib,
  double x[] /* 2 dimension */,
  CONST int *ix
);

extern void __stdcall F04AMF(
  CONST double a[] /* 2 dimension */,
  CONST int *ia,
  double x[] /* 2 dimension */,
  CONST int *ix,
  CONST double b[] /* 2 dimension */,
  CONST int *ib,
  CONST int *m,
  CONST int *n,
  CONST int *ip,
  CONST double *eta,
  double qr[] /* 2 dimension */,
  CONST int *iqr,
  double alpha[],
  double e[],
  double y[],
  double z[],
  double r[],
  int ipivot[],
  int *ifail
);

extern void __stdcall F04ANF(
  CONST int *m,
  CONST int *n,
  double qr[] /* 2 dimension */,
  CONST int *iqr,
  double alpha[],
  CONST int ipivot[],
  double r[],
  double y[],
  double z[]
);

extern void __stdcall F04AQF(
  CONST int *n,
  CONST int *m,
  CONST double rl[],
  CONST double d[],
  CONST double b[],
  double x[]
);

extern void __stdcall F04ARF(
  double a[] /* 2 dimension */,
  CONST int *ia,
  CONST double b[],
  CONST int *n,
  double c[],
  double wks[],
  int *ifail
);

extern void __stdcall F04ASF(
  double a[] /* 2 dimension */,
  CONST int *ia,
  CONST double b[],
  CONST int *n,
  double c[],
  double wk1[],
  double wk2[],
  int *ifail
);

extern void __stdcall F04ATF(
  CONST double a[] /* 2 dimension */,
  CONST int *ia,
  CONST double b[],
  CONST int *n,
  double c[],
  double aa[] /* 2 dimension */,
  CONST int *iaa,
  double wks1[],
  double wks2[],
  int *ifail
);

extern void __stdcall F04AUF(
  CONST int *m,
  CONST int *n,
  CONST double a[] /* 2 dimension */,
  CONST int *ia,
  CONST int *ip,
  CONST double d[],
  CONST int inc[],
  CONST int *ir,
  double b[] /* 2 dimension */,
  CONST int *ib,
  double x[] /* 2 dimension */,
  CONST int *ix,
  double u[] /* 2 dimension */,
  CONST int *iu,
  double du[],
  int *ifail
);

extern void __stdcall F04AWF(
  CONST int *n,
  CONST int *ir,
  Complex a[] /* 2 dimension */,
  CONST int *ia,
  double p[],
  CONST Complex b[] /* 2 dimension */,
  CONST int *ib,
  Complex x[] /* 2 dimension */,
  CONST int *ix
);

extern void __stdcall F04AXF(
  CONST int *n,
  CONST double a[],
  CONST int *licn,
  CONST int icn[],
  int ikeep[] /* 2 dimension */,
  double rhs[],
  double w[],
  CONST int *mtype,
  CONST int idisp[],
  double *resid
);

extern void __stdcall F04AYF(
  CONST int *n,
  CONST int *ir,
  CONST double a[] /* 2 dimension */,
  CONST int *ia,
  CONST double p[],
  double b[] /* 2 dimension */,
  CONST int *ib,
  int *ifail
);

extern void __stdcall F04AZF(
  CONST int *n,
  CONST int *ir,
  CONST double a[] /* 2 dimension */,
  CONST int *ia,
  CONST double p[],
  double b[] /* 2 dimension */,
  CONST int *ib,
  int *ifail
);

extern void __stdcall F04EAF(
  CONST int *n,
  double d[],
  double du[],
  double dl[],
  double b[],
  int *ifail
);

extern void __stdcall F04FAF(
  CONST int *job,
  CONST int *n,
  double d[],
  double e[],
  double b[],
  int *ifail
);

extern void __stdcall F04FEF(
  CONST int *n,
  CONST double t[],
  double x[],
  CONST int *wantp,
  double p[],
  CONST int *wantv,
  double v[],
  double *vlast,
  double work[],
  int *ifail
);

extern void __stdcall F04FFF(
  CONST int *n,
  CONST double t[],
  CONST double b[],
  double x[],
  CONST int *wantp,
  double p[],
  double work[],
  int *ifail
);

extern void __stdcall F04JAF(
  CONST int *m,
  CONST int *n,
  double a[] /* 2 dimension */,
  CONST int *nra,
  double b[],
  CONST double *tol,
  double *sigma,
  int *irank,
  double work[],
  CONST int *lwork,
  int *ifail
);

extern void __stdcall F04JDF(
  CONST int *m,
  CONST int *n,
  double a[] /* 2 dimension */,
  CONST int *nra,
  double b[],
  CONST double *tol,
  double *sigma,
  int *irank,
  double work[],
  CONST int *lwork,
  int *ifail
);

extern void __stdcall F04JGF(
  CONST int *m,
  CONST int *n,
  double a[] /* 2 dimension */,
  CONST int *nra,
  double b[],
  CONST double *tol,
  int *svd,
  double *sigma,
  int *irank,
  double work[],
  CONST int *lwork,
  int *ifail
);

extern void __stdcall F04JLF(
  CONST int *m,
  CONST int *n,
  CONST int *p,
  double a[] /* 2 dimension */,
  CONST int *lda,
  double b[] /* 2 dimension */,
  CONST int *ldb,
  double d[],
  double x[],
  double y[],
  double work[],
  CONST int *lwork,
  int *ifail
);

extern void __stdcall F04JMF(
  CONST int *m,
  CONST int *n,
  CONST int *p,
  double a[] /* 2 dimension */,
  CONST int *lda,
  double b[] /* 2 dimension */,
  CONST int *ldb,
  double c[],
  double d[],
  double x[],
  double work[],
  CONST int *lwork,
  int *ifail
);

extern void __stdcall F04KLF(
  CONST int *m,
  CONST int *n,
  CONST int *p,
  Complex a[] /* 2 dimension */,
  CONST int *lda,
  Complex b[] /* 2 dimension */,
  CONST int *ldb,
  Complex d[],
  Complex x[],
  Complex y[],
  Complex work[],
  CONST int *lwork,
  int *ifail
);

extern void __stdcall F04KMF(
  CONST int *m,
  CONST int *n,
  CONST int *p,
  Complex a[] /* 2 dimension */,
  CONST int *lda,
  Complex b[] /* 2 dimension */,
  CONST int *ldb,
  Complex c[],
  Complex d[],
  Complex x[],
  Complex work[],
  CONST int *lwork,
  int *ifail
);

extern void __stdcall F04LDF(
  CONST int *n,
  CONST int *m1,
  CONST int *m2,
  CONST int *ir,
  CONST double a[] /* 2 dimension */,
  CONST int *ia,
  CONST double al[] /* 2 dimension */,
  CONST int *il,
  CONST int in[],
  double b[] /* 2 dimension */,
  CONST int *ib,
  int *ifail
);

extern void __stdcall F04LEF(
  CONST int *job,
  CONST int *n,
  CONST double a[],
  CONST double b[],
  CONST double c[],
  CONST double d[],
  CONST int in[],
  double y[],
  double *tol,
  int *ifail
);

extern void __stdcall F04LHF(
  CONST char *trans, CONST int trans_len,
  CONST int *n,
  CONST int *nbloks,
  CONST int blkstr[] /* 2 dimension */,
  CONST double a[],
  CONST int *lena,
  CONST int pivot[],
  double b[] /* 2 dimension */,
  CONST int *ldb,
  CONST int *ir,
  int *ifail
);

extern void __stdcall F04MAF(
  CONST int *n,
  CONST int *nz,
  CONST double a[],
  CONST int *la,
  CONST int ini[],
  CONST int *lini,
  CONST int inj[],
  double b[],
  double eps[],
  int kmax[],
  double w[] /* 2 dimension */,
  double w1[] /* 2 dimension */,
  CONST int ik[] /* 2 dimension */,
  CONST int inform[],
  int *iflag
);

extern void __stdcall F04MBF(
  CONST int *n,
  CONST double b[],
  double x[],
  void (__stdcall *aprod)(int *, int *, double[], double[], double[], int *, int[], 
                int *),
  void (__stdcall *msolve)(int *, int *, double[], double[], double[], int *, int[], 
                 int *),
  CONST int *precon,
  CONST double *shift,
  CONST double *rtol,
  CONST int *itnlim,
  CONST int *msglvl,
  int *itn,
  double *anorm,
  double *acond,
  double *rnorm,
  double *xnorm,
  double work[] /* 2 dimension */,
  CONST double rwork[],
  CONST int *lrwork,
  CONST int iwork[],
  CONST int *liwork,
  int *inform,
  int *ifail
);

extern void __stdcall F04MCF(
  CONST int *n,
  CONST double l[],
  CONST int *ll,
  CONST double d[],
  CONST int nrow[],
  CONST int *p,
  CONST double b[] /* 2 dimension */,
  CONST int *nrb,
  CONST int *iselct,
  double x[] /* 2 dimension */,
  CONST int *nrx,
  int *ifail
);

extern void __stdcall F04MEF(
  CONST int *n,
  CONST double t[],
  double x[],
  double *v,
  double work[],
  int *ifail
);

extern void __stdcall F04MFF(
  CONST int *n,
  CONST double t[],
  CONST double b[],
  double x[],
  double *p,
  double work[],
  int *ifail
);

extern void __stdcall F04NAF(
  CONST int *job,
  CONST int *n,
  CONST int *ml,
  CONST int *mu,
  CONST Complex a[] /* 2 dimension */,
  CONST int *nra,
  CONST int in[],
  Complex b[],
  double *tol,
  int *ifail
);

extern void __stdcall F04QAF(
  CONST int *m,
  CONST int *n,
  double b[],
  double x[],
  double se[],
  void (__stdcall *aprod)(int *, int *, int *, double[], double[], double[], int *, 
                int[], int *),
  CONST double *damp,
  CONST double *atol,
  CONST double *btol,
  CONST double *conlim,
  int *itnlim,
  CONST int *msglvl,
  int *itn,
  double *anorm,
  double *acond,
  double *rnorm,
  double *arnorm,
  double *xnorm,
  double work[] /* 2 dimension */,
  CONST double rwork[],
  CONST int *lrwork,
  CONST int iwork[],
  CONST int *liwork,
  int *inform,
  int *ifail
);

extern void __stdcall F04YAF(
  CONST int *job,
  CONST int *p,
  CONST double *sigma,
  double a[] /* 2 dimension */,
  CONST int *nra,
  CONST int *svd,
  CONST int *irank,
  CONST double sv[],
  double cj[],
  double work[],
  int *ifail
);

extern void __stdcall F04YCF(
  int *icase,
  CONST int *n,
  double x[],
  double *estnrm,
  double work[],
  int iwork[],
  int *ifail
);

extern void __stdcall F04ZCF(
  int *icase,
  CONST int *n,
  Complex x[],
  double *estnrm,
  Complex work[],
  int *ifail
);

extern void __stdcall F05AAF(
  double a[] /* 2 dimension */,
  CONST int *ia,
  CONST int *m,
  CONST int *n1,
  CONST int *n2,
  double s[],
  double *cc,
  int *icol,
  int *ifail
);

extern void __stdcall F06AAF(
  double *a,
  double *b,
  double *c,
  double *s
);

extern void __stdcall F06BAF(
  double *a,
  double *b,
  double *c,
  double *s
);

extern void __stdcall F06BCF(
  CONST double *t,
  double *c,
  double *s
);

extern void __stdcall F06BEF(
  CONST char *job, CONST int job_len,
  double *x,
  double *y,
  double *z,
  double *c,
  double *s
);

extern void __stdcall F06BHF(
  double *x,
  double *y,
  double *z,
  CONST double *c,
  CONST double *s
);

extern double __stdcall F06BLF(
  CONST double *a,
  CONST double *b,
  int *fail
);

extern double __stdcall F06BMF(
  CONST double *scale,
  CONST double *ssq
);

extern double __stdcall F06BNF(
  CONST double *a,
  CONST double *b
);

extern double __stdcall F06BPF(
  CONST double *a,
  CONST double *b,
  CONST double *c
);

extern void __stdcall F06CAF(
  Complex *a,
  Complex *b,
  double *c,
  Complex *s
);

extern void __stdcall F06CBF(
  Complex *a,
  Complex *b,
  Complex *c,
  double *s
);

extern void __stdcall F06CCF(
  CONST Complex *t,
  double *c,
  Complex *s
);

extern void __stdcall F06CDF(
  Complex *t,
  Complex *c,
  double *s
);

extern void __stdcall F06CHF(
  Complex *x,
  Complex *y,
  Complex *z,
  CONST double *c,
  CONST Complex *s
);

extern void __stdcall F06CLF(
  Complex * F06CLF_,
  CONST Complex *a,
  CONST Complex *b,
  int *fail
);

extern void __stdcall F06DBF(
  CONST int *n,
  CONST int *_const,
  int x[],
  CONST int *incx
);

extern void __stdcall F06DFF(
  CONST int *n,
  CONST int x[],
  CONST int *incx,
  int y[],
  CONST int *incy
);

extern double __stdcall F06EAF(
  CONST int *n,
  CONST double x[],
  CONST int *incx,
  CONST double y[],
  CONST int *incy
);

extern void __stdcall F06ECF(
  CONST int *n,
  CONST double *alpha,
  CONST double x[],
  CONST int *incx,
  double y[],
  CONST int *incy
);

extern void __stdcall F06EDF(
  CONST int *n,
  CONST double *alpha,
  double x[],
  CONST int *incx
);

extern void __stdcall F06EFF(
  CONST int *n,
  CONST double x[],
  CONST int *incx,
  double y[],
  CONST int *incy
);

extern void __stdcall F06EGF(
  CONST int *n,
  double x[],
  CONST int *incx,
  double y[],
  CONST int *incy
);

extern double __stdcall F06EJF(
  CONST int *n,
  CONST double x[],
  CONST int *incx
);

extern double __stdcall F06EKF(
  CONST int *n,
  CONST double x[],
  CONST int *incx
);

extern void __stdcall F06EPF(
  CONST int *n,
  double x[],
  CONST int *incx,
  double y[],
  CONST int *incy,
  CONST double *c,
  CONST double *s
);

extern double __stdcall F06ERF(
  CONST int *nz,
  CONST double x[],
  CONST int indx[],
  CONST double y[]
);

extern void __stdcall F06ETF(
  CONST int *nz,
  CONST double *a,
  CONST double x[],
  CONST int indx[],
  double y[]
);

extern void __stdcall F06EUF(
  CONST int *nz,
  CONST double y[],
  double x[],
  CONST int indx[]
);

extern void __stdcall F06EVF(
  CONST int *nz,
  double y[],
  double x[],
  CONST int indx[]
);

extern void __stdcall F06EWF(
  CONST int *nz,
  CONST double x[],
  CONST int indx[],
  double y[]
);

extern void __stdcall F06EXF(
  CONST int *nz,
  double x[],
  CONST int indx[],
  double y[],
  CONST double *c,
  CONST double *s
);

extern double __stdcall F06FAF(
  CONST int *n,
  CONST int *j,
  CONST double *tolx,
  CONST double x[],
  CONST int *incx,
  CONST double *toly,
  CONST double y[],
  CONST int *incy
);

extern void __stdcall F06FBF(
  CONST int *n,
  CONST double *_const,
  double x[],
  CONST int *incx
);

extern void __stdcall F06FCF(
  CONST int *n,
  CONST double d[],
  CONST int *incd,
  double x[],
  CONST int *incx
);

extern void __stdcall F06FDF(
  CONST int *n,
  CONST double *alpha,
  CONST double x[],
  CONST int *incx,
  double y[],
  CONST int *incy
);

extern void __stdcall F06FGF(
  CONST int *n,
  double x[],
  CONST int *incx
);

extern void __stdcall F06FJF(
  CONST int *n,
  CONST double x[],
  CONST int *incx,
  double *scale,
  double *sumsq
);

extern double __stdcall F06FKF(
  CONST int *n,
  CONST double w[],
  CONST int *incw,
  CONST double x[],
  CONST int *incx
);

extern void __stdcall F06FLF(
  CONST int *n,
  CONST double x[],
  CONST int *incx,
  double *xmax,
  double *xmin
);

extern void __stdcall F06FPF(
  CONST int *n,
  double x[],
  CONST int *incx,
  double y[],
  CONST int *incy,
  CONST double *c,
  CONST double *s
);

extern void __stdcall F06FQF(
  CONST char *pivot, CONST int pivot_len,
  CONST char *direct, CONST int direct_len,
  CONST int *n,
  double *alpha,
  double x[],
  CONST int *incx,
  double c[],
  double s[]
);

extern void __stdcall F06FRF(
  CONST int *n,
  double *alpha,
  double x[],
  CONST int *incx,
  CONST double *tol,
  double *zeta
);

extern void __stdcall F06FSF(
  CONST int *n,
  double *alpha,
  double x[],
  CONST int *incx,
  CONST double *tol,
  double *z1
);

extern void __stdcall F06FTF(
  CONST int *n,
  double *delta,
  double y[],
  CONST int *incy,
  CONST double *zeta,
  CONST double z[],
  CONST int *incz
);

extern void __stdcall F06FUF(
  CONST int *n,
  CONST double z[],
  CONST int *incz,
  CONST double *z1,
  double *alpha,
  double x[],
  CONST int *incx
);

extern void __stdcall F06GAF(
  Complex * F06GAF_,
  CONST int *n,
  CONST Complex x[],
  CONST int *incx,
  CONST Complex y[],
  CONST int *incy
);

extern void __stdcall F06GBF(
  Complex * F06GBF_,
  CONST int *n,
  CONST Complex x[],
  CONST int *incx,
  CONST Complex y[],
  CONST int *incy
);

extern void __stdcall F06GCF(
  CONST int *n,
  CONST Complex *alpha,
  CONST Complex x[],
  CONST int *incx,
  Complex y[],
  CONST int *incy
);

extern void __stdcall F06GDF(
  CONST int *n,
  CONST Complex *alpha,
  Complex x[],
  CONST int *incx
);

extern void __stdcall F06GFF(
  CONST int *n,
  CONST Complex x[],
  CONST int *incx,
  Complex y[],
  CONST int *incy
);

extern void __stdcall F06GGF(
  CONST int *n,
  Complex x[],
  CONST int *incx,
  Complex y[],
  CONST int *incy
);

extern void __stdcall F06GRF(
  Complex * F06GRF_,
  CONST int *nz,
  CONST Complex x[],
  CONST int indx[],
  CONST Complex y[]
);

extern void __stdcall F06GSF(
  Complex * F06GSF_,
  CONST int *nz,
  CONST Complex x[],
  CONST int indx[],
  CONST Complex y[]
);

extern void __stdcall F06GTF(
  CONST int *nz,
  CONST Complex *a,
  CONST Complex x[],
  CONST int indx[],
  Complex y[]
);

extern void __stdcall F06GUF(
  CONST int *nz,
  CONST Complex y[],
  Complex x[],
  CONST int indx[]
);

extern void __stdcall F06GVF(
  CONST int *nz,
  Complex y[],
  Complex x[],
  CONST int indx[]
);

extern void __stdcall F06GWF(
  CONST int *nz,
  CONST Complex x[],
  CONST int indx[],
  Complex y[]
);

extern void __stdcall F06HBF(
  CONST int *n,
  CONST Complex *_const,
  Complex x[],
  CONST int *incx
);

extern void __stdcall F06HCF(
  CONST int *n,
  CONST Complex d[],
  CONST int *incd,
  Complex x[],
  CONST int *incx
);

extern void __stdcall F06HDF(
  CONST int *n,
  CONST Complex *alpha,
  CONST Complex x[],
  CONST int *incx,
  Complex y[],
  CONST int *incy
);

extern void __stdcall F06HGF(
  CONST int *n,
  Complex x[],
  CONST int *incx
);

extern void __stdcall F06HPF(
  CONST int *n,
  Complex x[],
  CONST int *incx,
  Complex y[],
  CONST int *incy,
  CONST Complex *c,
  CONST Complex *s
);

extern void __stdcall F06HQF(
  CONST char *pivot, CONST int pivot_len,
  CONST char *direct, CONST int direct_len,
  CONST int *n,
  Complex *alpha,
  Complex x[],
  CONST int *incx,
  double c[],
  Complex s[]
);

extern void __stdcall zdotu(
  Complex * zdotu_,
  CONST int *n,
  CONST Complex x[],
  CONST int *incx,
  CONST Complex y[],
  CONST int *incy
);

extern void __stdcall zdotc(
  Complex * zdotc_,
  CONST int *n,
  CONST Complex x[],
  CONST int *incx,
  CONST Complex y[],
  CONST int *incy
);

extern void __stdcall zaxpy(
  CONST int *n,
  CONST Complex *alpha,
  CONST Complex x[],
  CONST int *incx,
  Complex y[],
  CONST int *incy
);

extern void __stdcall zscal(
  CONST int *n,
  CONST Complex *alpha,
  Complex x[],
  CONST int *incx
);

extern void __stdcall zcopy(
  CONST int *n,
  CONST Complex x[],
  CONST int *incx,
  Complex y[],
  CONST int *incy
);

extern void __stdcall zswap(
  CONST int *n,
  Complex x[],
  CONST int *incx,
  Complex y[],
  CONST int *incy
);

extern void __stdcall zdotui(
  Complex * zdotui_,
  CONST int *nz,
  CONST Complex x[],
  CONST int indx[],
  CONST Complex y[]
);

extern void __stdcall zdotci(
  Complex * zdotci_,
  CONST int *nz,
  CONST Complex x[],
  CONST int indx[],
  CONST Complex y[]
);

extern void __stdcall zaxpyi(
  CONST int *nz,
  CONST Complex *a,
  CONST Complex x[],
  CONST int indx[],
  Complex y[]
);

extern void __stdcall zgthr(
  CONST int *nz,
  CONST Complex y[],
  Complex x[],
  CONST int indx[]
);

extern void __stdcall zgthrz(
  CONST int *nz,
  Complex y[],
  Complex x[],
  CONST int indx[]
);

extern void __stdcall zsctr(
  CONST int *nz,
  CONST Complex x[],
  CONST int indx[],
  Complex y[]
);

extern void __stdcall F06HRF(
  CONST int *n,
  Complex *alpha,
  Complex x[],
  CONST int *incx,
  CONST double *tol,
  Complex *theta
);

extern void __stdcall F06HTF(
  CONST int *n,
  Complex *delta,
  Complex y[],
  CONST int *incy,
  CONST Complex *theta,
  CONST Complex z[],
  CONST int *incz
);

extern void __stdcall F06JDF(
  CONST int *n,
  CONST double *alpha,
  Complex x[],
  CONST int *incx
);

extern double __stdcall F06JJF(
  CONST int *n,
  CONST Complex x[],
  CONST int *incx
);

extern double __stdcall F06JKF(
  CONST int *n,
  CONST Complex x[],
  CONST int *incx
);

extern int __stdcall F06JLF(
  CONST int *n,
  CONST double x[],
  CONST int *incx
);

extern int __stdcall F06JMF(
  CONST int *n,
  CONST Complex x[],
  CONST int *incx
);

extern void __stdcall F06KCF(
  CONST int *n,
  CONST double d[],
  CONST int *incd,
  Complex x[],
  CONST int *incx
);

extern void __stdcall F06KDF(
  CONST int *n,
  CONST double *alpha,
  CONST Complex x[],
  CONST int *incx,
  Complex y[],
  CONST int *incy
);

extern void __stdcall F06KFF(
  CONST int *n,
  CONST double x[],
  CONST int *incx,
  Complex y[],
  CONST int *incy
);

extern void __stdcall F06KJF(
  CONST int *n,
  CONST Complex x[],
  CONST int *incx,
  double *scale,
  double *sumsq
);

extern int __stdcall F06KLF(
  CONST int *n,
  CONST double x[],
  CONST int *incx,
  CONST double *tol
);

extern void __stdcall F06KPF(
  CONST int *n,
  Complex x[],
  CONST int *incx,
  Complex y[],
  CONST int *incy,
  CONST double *c,
  CONST double *s
);

extern void __stdcall F06PAF(
  CONST char *trans, CONST int trans_len,
  CONST int *m,
  CONST int *n,
  CONST double *alpha,
  CONST double a[] /* 2 dimension */,
  CONST int *lda,
  CONST double x[],
  CONST int *incx,
  CONST double *beta,
  double y[],
  CONST int *incy
);

extern void __stdcall F06PBF(
  CONST char *trans, CONST int trans_len,
  CONST int *m,
  CONST int *n,
  CONST int *kl,
  CONST int *ku,
  CONST double *alpha,
  CONST double a[] /* 2 dimension */,
  CONST int *lda,
  CONST double x[],
  CONST int *incx,
  CONST double *beta,
  double y[],
  CONST int *incy
);

extern void __stdcall F06PCF(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST double *alpha,
  CONST double a[] /* 2 dimension */,
  CONST int *lda,
  CONST double x[],
  CONST int *incx,
  CONST double *beta,
  double y[],
  CONST int *incy
);

extern void __stdcall F06PDF(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST int *k,
  CONST double *alpha,
  CONST double a[] /* 2 dimension */,
  CONST int *lda,
  CONST double x[],
  CONST int *incx,
  CONST double *beta,
  double y[],
  CONST int *incy
);

extern void __stdcall F06PEF(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST double *alpha,
  CONST double ap[],
  CONST double x[],
  CONST int *incx,
  CONST double *beta,
  double y[],
  CONST int *incy
);

extern void __stdcall F06PFF(
  CONST char *uplo, CONST int uplo_len,
  CONST char *trans, CONST int trans_len,
  CONST char *diag, CONST int diag_len,
  CONST int *n,
  CONST double a[] /* 2 dimension */,
  CONST int *lda,
  double x[],
  CONST int *incx
);

extern void __stdcall F06PGF(
  CONST char *uplo, CONST int uplo_len,
  CONST char *trans, CONST int trans_len,
  CONST char *diag, CONST int diag_len,
  CONST int *n,
  CONST int *k,
  CONST double a[] /* 2 dimension */,
  CONST int *lda,
  double x[],
  CONST int *incx
);

extern void __stdcall F06PHF(
  CONST char *uplo, CONST int uplo_len,
  CONST char *trans, CONST int trans_len,
  CONST char *diag, CONST int diag_len,
  CONST int *n,
  CONST double ap[],
  double x[],
  CONST int *incx
);

extern void __stdcall F06PJF(
  CONST char *uplo, CONST int uplo_len,
  CONST char *trans, CONST int trans_len,
  CONST char *diag, CONST int diag_len,
  CONST int *n,
  CONST double a[] /* 2 dimension */,
  CONST int *lda,
  double x[],
  CONST int *incx
);

extern void __stdcall F06PKF(
  CONST char *uplo, CONST int uplo_len,
  CONST char *trans, CONST int trans_len,
  CONST char *diag, CONST int diag_len,
  CONST int *n,
  CONST int *k,
  CONST double a[] /* 2 dimension */,
  CONST int *lda,
  double x[],
  CONST int *incx
);

extern void __stdcall F06PLF(
  CONST char *uplo, CONST int uplo_len,
  CONST char *trans, CONST int trans_len,
  CONST char *diag, CONST int diag_len,
  CONST int *n,
  CONST double ap[],
  double x[],
  CONST int *incx
);

extern void __stdcall F06PMF(
  CONST int *m,
  CONST int *n,
  CONST double *alpha,
  CONST double x[],
  CONST int *incx,
  CONST double y[],
  CONST int *incy,
  double a[] /* 2 dimension */,
  CONST int *lda
);

extern void __stdcall F06PPF(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST double *alpha,
  CONST double x[],
  CONST int *incx,
  double a[] /* 2 dimension */,
  CONST int *lda
);

extern void __stdcall F06PQF(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST double *alpha,
  CONST double x[],
  CONST int *incx,
  double ap[]
);

extern void __stdcall F06PRF(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST double *alpha,
  CONST double x[],
  CONST int *incx,
  CONST double y[],
  CONST int *incy,
  double a[] /* 2 dimension */,
  CONST int *lda
);

extern void __stdcall F06PSF(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST double *alpha,
  CONST double x[],
  CONST int *incx,
  CONST double y[],
  CONST int *incy,
  double ap[]
);

extern void __stdcall F06QFF(
  CONST char *matrix, CONST int matrix_len,
  CONST int *m,
  CONST int *n,
  CONST double a[] /* 2 dimension */,
  CONST int *lda,
  double b[] /* 2 dimension */,
  CONST int *ldb
);

extern void __stdcall zdscal(
  CONST int *n,
  CONST double *alpha,
  Complex x[],
  CONST int *incx
);

extern double __stdcall dznrm2(
  CONST int *n,
  CONST Complex x[],
  CONST int *incx
);

extern double __stdcall dzasum(
  CONST int *n,
  CONST Complex x[],
  CONST int *incx
);

extern int __stdcall idamax(
  CONST int *n,
  CONST double x[],
  CONST int *incx
);

extern int __stdcall izamax(
  CONST int *n,
  CONST Complex x[],
  CONST int *incx
);

extern void __stdcall dgemv(
  CONST char *trans, CONST int trans_len,
  CONST int *m,
  CONST int *n,
  CONST double *alpha,
  CONST double a[] /* 2 dimension */,
  CONST int *lda,
  CONST double x[],
  CONST int *incx,
  CONST double *beta,
  double y[],
  CONST int *incy
);

extern void __stdcall dgbmv(
  CONST char *trans, CONST int trans_len,
  CONST int *m,
  CONST int *n,
  CONST int *kl,
  CONST int *ku,
  CONST double *alpha,
  CONST double a[] /* 2 dimension */,
  CONST int *lda,
  CONST double x[],
  CONST int *incx,
  CONST double *beta,
  double y[],
  CONST int *incy
);

extern void __stdcall dsymv(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST double *alpha,
  CONST double a[] /* 2 dimension */,
  CONST int *lda,
  CONST double x[],
  CONST int *incx,
  CONST double *beta,
  double y[],
  CONST int *incy
);

extern void __stdcall dsbmv(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST int *k,
  CONST double *alpha,
  CONST double a[] /* 2 dimension */,
  CONST int *lda,
  CONST double x[],
  CONST int *incx,
  CONST double *beta,
  double y[],
  CONST int *incy
);

extern void __stdcall dspmv(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST double *alpha,
  CONST double ap[],
  CONST double x[],
  CONST int *incx,
  CONST double *beta,
  double y[],
  CONST int *incy
);

extern void __stdcall dtrmv(
  CONST char *uplo, CONST int uplo_len,
  CONST char *trans, CONST int trans_len,
  CONST char *diag, CONST int diag_len,
  CONST int *n,
  CONST double a[] /* 2 dimension */,
  CONST int *lda,
  double x[],
  CONST int *incx
);

extern void __stdcall dtbmv(
  CONST char *uplo, CONST int uplo_len,
  CONST char *trans, CONST int trans_len,
  CONST char *diag, CONST int diag_len,
  CONST int *n,
  CONST int *k,
  CONST double a[] /* 2 dimension */,
  CONST int *lda,
  double x[],
  CONST int *incx
);

extern void __stdcall dtpmv(
  CONST char *uplo, CONST int uplo_len,
  CONST char *trans, CONST int trans_len,
  CONST char *diag, CONST int diag_len,
  CONST int *n,
  CONST double ap[],
  double x[],
  CONST int *incx
);

extern void __stdcall dtrsv(
  CONST char *uplo, CONST int uplo_len,
  CONST char *trans, CONST int trans_len,
  CONST char *diag, CONST int diag_len,
  CONST int *n,
  CONST double a[] /* 2 dimension */,
  CONST int *lda,
  double x[],
  CONST int *incx
);

extern void __stdcall dtbsv(
  CONST char *uplo, CONST int uplo_len,
  CONST char *trans, CONST int trans_len,
  CONST char *diag, CONST int diag_len,
  CONST int *n,
  CONST int *k,
  CONST double a[] /* 2 dimension */,
  CONST int *lda,
  double x[],
  CONST int *incx
);

extern void __stdcall dtpsv(
  CONST char *uplo, CONST int uplo_len,
  CONST char *trans, CONST int trans_len,
  CONST char *diag, CONST int diag_len,
  CONST int *n,
  CONST double ap[],
  double x[],
  CONST int *incx
);

extern void __stdcall dger(
  CONST int *m,
  CONST int *n,
  CONST double *alpha,
  CONST double x[],
  CONST int *incx,
  CONST double y[],
  CONST int *incy,
  double a[] /* 2 dimension */,
  CONST int *lda
);

extern void __stdcall dsyr(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST double *alpha,
  CONST double x[],
  CONST int *incx,
  double a[] /* 2 dimension */,
  CONST int *lda
);

extern void __stdcall dspr(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST double *alpha,
  CONST double x[],
  CONST int *incx,
  double ap[]
);

extern void __stdcall dsyr2(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST double *alpha,
  CONST double x[],
  CONST int *incx,
  CONST double y[],
  CONST int *incy,
  double a[] /* 2 dimension */,
  CONST int *lda
);

extern void __stdcall dspr2(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST double *alpha,
  CONST double x[],
  CONST int *incx,
  CONST double y[],
  CONST int *incy,
  double ap[]
);

extern double __stdcall F06QGF(
  CONST char *norm, CONST int norm_len,
  CONST char *matrix, CONST int matrix_len,
  CONST int *m,
  CONST int *n,
  CONST double a[] /* 2 dimension */,
  CONST int *lda
);

extern void __stdcall F06QHF(
  CONST char *matrix, CONST int matrix_len,
  CONST int *m,
  CONST int *n,
  CONST double *_const,
  CONST double *diag,
  double a[] /* 2 dimension */,
  CONST int *lda
);

extern void __stdcall F06QJF(
  CONST char *side, CONST int side_len,
  CONST char *trans, CONST int trans_len,
  CONST int *n,
  CONST int perm[],
  CONST int *k,
  double b[] /* 2 dimension */,
  CONST int *ldb
);

extern void __stdcall F06QKF(
  CONST char *side, CONST int side_len,
  CONST char *trans, CONST int trans_len,
  CONST int *n,
  CONST double perm[],
  CONST int *k,
  double b[] /* 2 dimension */,
  CONST int *ldb
);

extern void __stdcall F06QMF(
  CONST char *uplo, CONST int uplo_len,
  CONST char *pivot, CONST int pivot_len,
  CONST char *direct, CONST int direct_len,
  CONST int *n,
  CONST int *k1,
  CONST int *k2,
  CONST double c[],
  CONST double s[],
  double a[] /* 2 dimension */,
  CONST int *lda
);

extern void __stdcall F06QPF(
  CONST int *n,
  CONST double *alpha,
  double x[],
  CONST int *incx,
  CONST double y[],
  CONST int *incy,
  double a[] /* 2 dimension */,
  CONST int *lda,
  double c[],
  double s[]
);

extern void __stdcall F06QQF(
  CONST int *n,
  CONST double *alpha,
  double x[],
  CONST int *incx,
  double a[] /* 2 dimension */,
  CONST int *lda,
  double c[],
  double s[]
);

extern void __stdcall F06QRF(
  CONST char *side, CONST int side_len,
  CONST int *n,
  CONST int *k1,
  CONST int *k2,
  double c[],
  double s[],
  double a[] /* 2 dimension */,
  CONST int *lda
);

extern void __stdcall F06QSF(
  CONST char *side, CONST int side_len,
  CONST int *n,
  CONST int *k1,
  CONST int *k2,
  double c[],
  double s[],
  double a[] /* 2 dimension */,
  CONST int *lda
);

extern void __stdcall F06QTF(
  CONST char *side, CONST int side_len,
  CONST int *n,
  CONST int *k1,
  CONST int *k2,
  double c[],
  double s[],
  double a[] /* 2 dimension */,
  CONST int *lda
);

extern void __stdcall F06QVF(
  CONST char *side, CONST int side_len,
  CONST int *n,
  CONST int *k1,
  CONST int *k2,
  CONST double c[],
  double s[],
  double a[] /* 2 dimension */,
  CONST int *lda
);

extern void __stdcall F06QWF(
  CONST char *side, CONST int side_len,
  CONST int *n,
  CONST int *k1,
  CONST int *k2,
  CONST double c[],
  double s[],
  double a[] /* 2 dimension */,
  CONST int *lda
);

extern void __stdcall F06QXF(
  CONST char *side, CONST int side_len,
  CONST char *pivot, CONST int pivot_len,
  CONST char *direct, CONST int direct_len,
  CONST int *m,
  CONST int *n,
  CONST int *k1,
  CONST int *k2,
  CONST double c[],
  CONST double s[],
  double a[] /* 2 dimension */,
  CONST int *lda
);

extern double __stdcall F06RAF(
  CONST char *norm, CONST int norm_len,
  CONST int *m,
  CONST int *n,
  CONST double a[] /* 2 dimension */,
  CONST int *lda,
  double work[]
);

extern double __stdcall F06RBF(
  CONST char *norm, CONST int norm_len,
  CONST int *n,
  CONST int *kl,
  CONST int *ku,
  CONST double ab[] /* 2 dimension */,
  CONST int *ldab,
  double work[]
);

extern double __stdcall F06RCF(
  CONST char *norm, CONST int norm_len,
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST double a[] /* 2 dimension */,
  CONST int *lda,
  double work[]
);

extern double __stdcall F06RDF(
  CONST char *norm, CONST int norm_len,
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST double ap[],
  double work[]
);

extern double __stdcall F06REF(
  CONST char *norm, CONST int norm_len,
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST int *k,
  CONST double ab[] /* 2 dimension */,
  CONST int *ldab,
  double work[]
);

extern double __stdcall F06RJF(
  CONST char *norm, CONST int norm_len,
  CONST char *uplo, CONST int uplo_len,
  CONST char *diag, CONST int diag_len,
  CONST int *m,
  CONST int *n,
  CONST double a[] /* 2 dimension */,
  CONST int *lda,
  double work[]
);

extern double __stdcall F06RKF(
  CONST char *norm, CONST int norm_len,
  CONST char *uplo, CONST int uplo_len,
  CONST char *diag, CONST int diag_len,
  CONST int *n,
  CONST double ap[],
  double work[]
);

extern double __stdcall F06RLF(
  CONST char *norm, CONST int norm_len,
  CONST char *uplo, CONST int uplo_len,
  CONST char *diag, CONST int diag_len,
  CONST int *n,
  CONST int *k,
  CONST double ab[] /* 2 dimension */,
  CONST int *ldab,
  double work[]
);

extern double __stdcall F06RMF(
  CONST char *norm, CONST int norm_len,
  CONST int *n,
  CONST double a[] /* 2 dimension */,
  CONST int *lda,
  double work[]
);

extern void __stdcall F06SAF(
  CONST char *trans, CONST int trans_len,
  CONST int *m,
  CONST int *n,
  CONST Complex *alpha,
  CONST Complex a[] /* 2 dimension */,
  CONST int *lda,
  CONST Complex x[],
  CONST int *incx,
  CONST Complex *beta,
  Complex y[],
  CONST int *incy
);

extern void __stdcall F06SBF(
  CONST char *trans, CONST int trans_len,
  CONST int *m,
  CONST int *n,
  CONST int *kl,
  CONST int *ku,
  CONST Complex *alpha,
  CONST Complex a[] /* 2 dimension */,
  CONST int *lda,
  CONST Complex x[],
  CONST int *incx,
  CONST Complex *beta,
  Complex y[],
  CONST int *incy
);

extern void __stdcall F06SCF(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST Complex *alpha,
  CONST Complex a[] /* 2 dimension */,
  CONST int *lda,
  CONST Complex x[],
  CONST int *incx,
  CONST Complex *beta,
  Complex y[],
  CONST int *incy
);

extern void __stdcall F06SDF(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST int *k,
  CONST Complex *alpha,
  CONST Complex a[] /* 2 dimension */,
  CONST int *lda,
  CONST Complex x[],
  CONST int *incx,
  CONST Complex *beta,
  Complex y[],
  CONST int *incy
);

extern void __stdcall F06SEF(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST Complex *alpha,
  CONST Complex ap[],
  CONST Complex x[],
  CONST int *incx,
  CONST Complex *beta,
  Complex y[],
  CONST int *incy
);

extern void __stdcall F06SFF(
  CONST char *uplo, CONST int uplo_len,
  CONST char *trans, CONST int trans_len,
  CONST char *diag, CONST int diag_len,
  CONST int *n,
  CONST Complex a[] /* 2 dimension */,
  CONST int *lda,
  Complex x[],
  CONST int *incx
);

extern void __stdcall zgemv(
  CONST char *trans, CONST int trans_len,
  CONST int *m,
  CONST int *n,
  CONST Complex *alpha,
  CONST Complex a[] /* 2 dimension */,
  CONST int *lda,
  CONST Complex x[],
  CONST int *incx,
  CONST Complex *beta,
  Complex y[],
  CONST int *incy
);

extern void __stdcall zgbmv(
  CONST char *trans, CONST int trans_len,
  CONST int *m,
  CONST int *n,
  CONST int *kl,
  CONST int *ku,
  CONST Complex *alpha,
  CONST Complex a[] /* 2 dimension */,
  CONST int *lda,
  CONST Complex x[],
  CONST int *incx,
  CONST Complex *beta,
  Complex y[],
  CONST int *incy
);

extern void __stdcall zhemv(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST Complex *alpha,
  CONST Complex a[] /* 2 dimension */,
  CONST int *lda,
  CONST Complex x[],
  CONST int *incx,
  CONST Complex *beta,
  Complex y[],
  CONST int *incy
);

extern void __stdcall zhbmv(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST int *k,
  CONST Complex *alpha,
  CONST Complex a[] /* 2 dimension */,
  CONST int *lda,
  CONST Complex x[],
  CONST int *incx,
  CONST Complex *beta,
  Complex y[],
  CONST int *incy
);

extern void __stdcall zhpmv(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST Complex *alpha,
  CONST Complex ap[],
  CONST Complex x[],
  CONST int *incx,
  CONST Complex *beta,
  Complex y[],
  CONST int *incy
);

extern void __stdcall ztrmv(
  CONST char *uplo, CONST int uplo_len,
  CONST char *trans, CONST int trans_len,
  CONST char *diag, CONST int diag_len,
  CONST int *n,
  CONST Complex a[] /* 2 dimension */,
  CONST int *lda,
  Complex x[],
  CONST int *incx
);

extern void __stdcall F06SGF(
  CONST char *uplo, CONST int uplo_len,
  CONST char *trans, CONST int trans_len,
  CONST char *diag, CONST int diag_len,
  CONST int *n,
  CONST int *k,
  CONST Complex a[] /* 2 dimension */,
  CONST int *lda,
  Complex x[],
  CONST int *incx
);

extern void __stdcall F06SHF(
  CONST char *uplo, CONST int uplo_len,
  CONST char *trans, CONST int trans_len,
  CONST char *diag, CONST int diag_len,
  CONST int *n,
  CONST Complex ap[],
  Complex x[],
  CONST int *incx
);

extern void __stdcall F06SJF(
  CONST char *uplo, CONST int uplo_len,
  CONST char *trans, CONST int trans_len,
  CONST char *diag, CONST int diag_len,
  CONST int *n,
  CONST Complex a[] /* 2 dimension */,
  CONST int *lda,
  Complex x[],
  CONST int *incx
);

extern void __stdcall F06SKF(
  CONST char *uplo, CONST int uplo_len,
  CONST char *trans, CONST int trans_len,
  CONST char *diag, CONST int diag_len,
  CONST int *n,
  CONST int *k,
  CONST Complex a[] /* 2 dimension */,
  CONST int *lda,
  Complex x[],
  CONST int *incx
);

extern void __stdcall F06SLF(
  CONST char *uplo, CONST int uplo_len,
  CONST char *trans, CONST int trans_len,
  CONST char *diag, CONST int diag_len,
  CONST int *n,
  CONST Complex ap[],
  Complex x[],
  CONST int *incx
);

extern void __stdcall F06SMF(
  CONST int *m,
  CONST int *n,
  CONST Complex *alpha,
  CONST Complex x[],
  CONST int *incx,
  CONST Complex y[],
  CONST int *incy,
  Complex a[] /* 2 dimension */,
  CONST int *lda
);

extern void __stdcall F06SNF(
  CONST int *m,
  CONST int *n,
  CONST Complex *alpha,
  CONST Complex x[],
  CONST int *incx,
  CONST Complex y[],
  CONST int *incy,
  Complex a[] /* 2 dimension */,
  CONST int *lda
);

extern void __stdcall F06SPF(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST double *alpha,
  CONST Complex x[],
  CONST int *incx,
  Complex a[] /* 2 dimension */,
  CONST int *lda
);

extern void __stdcall F06SQF(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST double *alpha,
  CONST Complex x[],
  CONST int *incx,
  Complex ap[]
);

extern void __stdcall F06SRF(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST Complex *alpha,
  CONST Complex x[],
  CONST int *incx,
  CONST Complex y[],
  CONST int *incy,
  Complex a[] /* 2 dimension */,
  CONST int *lda
);

extern void __stdcall F06SSF(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST Complex *alpha,
  CONST Complex x[],
  CONST int *incx,
  CONST Complex y[],
  CONST int *incy,
  Complex ap[]
);

extern void __stdcall F06TFF(
  CONST char *matrix, CONST int matrix_len,
  CONST int *m,
  CONST int *n,
  CONST Complex a[] /* 2 dimension */,
  CONST int *lda,
  Complex b[] /* 2 dimension */,
  CONST int *ldb
);

extern void __stdcall F06THF(
  CONST char *matrix, CONST int matrix_len,
  CONST int *m,
  CONST int *n,
  CONST Complex *_const,
  CONST Complex *diag,
  Complex a[] /* 2 dimension */,
  CONST int *lda
);

extern void __stdcall F06TMF(
  CONST char *uplo, CONST int uplo_len,
  CONST char *pivot, CONST int pivot_len,
  CONST char *direct, CONST int direct_len,
  CONST int *n,
  CONST int *k1,
  CONST int *k2,
  CONST double c[],
  CONST Complex s[],
  Complex a[] /* 2 dimension */,
  CONST int *lda
);

extern void __stdcall F06TPF(
  CONST int *n,
  CONST Complex *alpha,
  Complex x[],
  CONST int *incx,
  CONST Complex y[],
  CONST int *incy,
  Complex a[] /* 2 dimension */,
  CONST int *lda,
  double c[],
  Complex s[]
);

extern void __stdcall F06TQF(
  CONST int *n,
  CONST Complex *alpha,
  Complex x[],
  CONST int *incx,
  Complex a[] /* 2 dimension */,
  CONST int *lda,
  double c[],
  Complex s[]
);

extern void __stdcall F06TRF(
  CONST char *side, CONST int side_len,
  CONST int *n,
  CONST int *k1,
  CONST int *k2,
  Complex c[],
  double s[],
  Complex a[] /* 2 dimension */,
  CONST int *lda
);

extern void __stdcall F06TSF(
  CONST char *side, CONST int side_len,
  CONST int *n,
  CONST int *k1,
  CONST int *k2,
  double c[],
  Complex s[],
  Complex a[] /* 2 dimension */,
  CONST int *lda
);

extern void __stdcall F06TTF(
  CONST char *side, CONST int side_len,
  CONST int *n,
  CONST int *k1,
  CONST int *k2,
  double c[],
  Complex s[],
  Complex a[] /* 2 dimension */,
  CONST int *lda
);

extern void __stdcall F06TVF(
  CONST char *side, CONST int side_len,
  CONST int *n,
  CONST int *k1,
  CONST int *k2,
  CONST Complex c[],
  double s[],
  Complex a[] /* 2 dimension */,
  CONST int *lda
);

extern void __stdcall F06TWF(
  CONST char *side, CONST int side_len,
  CONST int *n,
  CONST int *k1,
  CONST int *k2,
  CONST double c[],
  Complex s[],
  Complex a[] /* 2 dimension */,
  CONST int *lda
);

extern void __stdcall F06TXF(
  CONST char *side, CONST int side_len,
  CONST char *pivot, CONST int pivot_len,
  CONST char *direct, CONST int direct_len,
  CONST int *m,
  CONST int *n,
  CONST int *k1,
  CONST int *k2,
  CONST double c[],
  CONST Complex s[],
  Complex a[] /* 2 dimension */,
  CONST int *lda
);

extern void __stdcall F06TYF(
  CONST char *side, CONST int side_len,
  CONST char *pivot, CONST int pivot_len,
  CONST char *direct, CONST int direct_len,
  CONST int *m,
  CONST int *n,
  CONST int *k1,
  CONST int *k2,
  CONST Complex c[],
  CONST double s[],
  Complex a[] /* 2 dimension */,
  CONST int *lda
);

extern double __stdcall F06UAF(
  CONST char *norm, CONST int norm_len,
  CONST int *m,
  CONST int *n,
  CONST Complex a[] /* 2 dimension */,
  CONST int *lda,
  double work[]
);

extern double __stdcall F06UBF(
  CONST char *norm, CONST int norm_len,
  CONST int *n,
  CONST int *kl,
  CONST int *ku,
  CONST Complex ab[] /* 2 dimension */,
  CONST int *ldab,
  double work[]
);

extern double __stdcall F06UCF(
  CONST char *norm, CONST int norm_len,
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST Complex a[] /* 2 dimension */,
  CONST int *lda,
  double work[]
);

extern double __stdcall F06UDF(
  CONST char *norm, CONST int norm_len,
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST Complex ap[],
  double work[]
);

extern double __stdcall F06UEF(
  CONST char *norm, CONST int norm_len,
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST int *k,
  CONST Complex ab[] /* 2 dimension */,
  CONST int *ldab,
  double work[]
);

extern double __stdcall F06UFF(
  CONST char *norm, CONST int norm_len,
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST Complex a[] /* 2 dimension */,
  CONST int *lda,
  double work[]
);

extern double __stdcall F06UGF(
  CONST char *norm, CONST int norm_len,
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST Complex ap[],
  double work[]
);

extern void __stdcall ztbmv(
  CONST char *uplo, CONST int uplo_len,
  CONST char *trans, CONST int trans_len,
  CONST char *diag, CONST int diag_len,
  CONST int *n,
  CONST int *k,
  CONST Complex a[] /* 2 dimension */,
  CONST int *lda,
  Complex x[],
  CONST int *incx
);

extern void __stdcall ztpmv(
  CONST char *uplo, CONST int uplo_len,
  CONST char *trans, CONST int trans_len,
  CONST char *diag, CONST int diag_len,
  CONST int *n,
  CONST Complex ap[],
  Complex x[],
  CONST int *incx
);

extern void __stdcall ztrsv(
  CONST char *uplo, CONST int uplo_len,
  CONST char *trans, CONST int trans_len,
  CONST char *diag, CONST int diag_len,
  CONST int *n,
  CONST Complex a[] /* 2 dimension */,
  CONST int *lda,
  Complex x[],
  CONST int *incx
);

extern void __stdcall ztbsv(
  CONST char *uplo, CONST int uplo_len,
  CONST char *trans, CONST int trans_len,
  CONST char *diag, CONST int diag_len,
  CONST int *n,
  CONST int *k,
  CONST Complex a[] /* 2 dimension */,
  CONST int *lda,
  Complex x[],
  CONST int *incx
);

extern void __stdcall ztpsv(
  CONST char *uplo, CONST int uplo_len,
  CONST char *trans, CONST int trans_len,
  CONST char *diag, CONST int diag_len,
  CONST int *n,
  CONST Complex ap[],
  Complex x[],
  CONST int *incx
);

extern void __stdcall zgeru(
  CONST int *m,
  CONST int *n,
  CONST Complex *alpha,
  CONST Complex x[],
  CONST int *incx,
  CONST Complex y[],
  CONST int *incy,
  Complex a[] /* 2 dimension */,
  CONST int *lda
);

extern void __stdcall zgerc(
  CONST int *m,
  CONST int *n,
  CONST Complex *alpha,
  CONST Complex x[],
  CONST int *incx,
  CONST Complex y[],
  CONST int *incy,
  Complex a[] /* 2 dimension */,
  CONST int *lda
);

extern void __stdcall zher(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST double *alpha,
  CONST Complex x[],
  CONST int *incx,
  Complex a[] /* 2 dimension */,
  CONST int *lda
);

extern void __stdcall zhpr(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST double *alpha,
  CONST Complex x[],
  CONST int *incx,
  Complex ap[]
);

extern void __stdcall zher2(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST Complex *alpha,
  CONST Complex x[],
  CONST int *incx,
  CONST Complex y[],
  CONST int *incy,
  Complex a[] /* 2 dimension */,
  CONST int *lda
);

extern void __stdcall zhpr2(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST Complex *alpha,
  CONST Complex x[],
  CONST int *incx,
  CONST Complex y[],
  CONST int *incy,
  Complex ap[]
);

extern double __stdcall F06UHF(
  CONST char *norm, CONST int norm_len,
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST int *k,
  CONST Complex ab[] /* 2 dimension */,
  CONST int *ldab,
  double work[]
);

extern double __stdcall F06UJF(
  CONST char *norm, CONST int norm_len,
  CONST char *uplo, CONST int uplo_len,
  CONST char *diag, CONST int diag_len,
  CONST int *m,
  CONST int *n,
  CONST Complex a[] /* 2 dimension */,
  CONST int *lda,
  double work[]
);

extern double __stdcall F06UKF(
  CONST char *norm, CONST int norm_len,
  CONST char *uplo, CONST int uplo_len,
  CONST char *diag, CONST int diag_len,
  CONST int *n,
  CONST Complex ap[],
  double work[]
);

extern double __stdcall F06ULF(
  CONST char *norm, CONST int norm_len,
  CONST char *uplo, CONST int uplo_len,
  CONST char *diag, CONST int diag_len,
  CONST int *n,
  CONST int *k,
  CONST Complex ab[] /* 2 dimension */,
  CONST int *ldab,
  double work[]
);

extern double __stdcall F06UMF(
  CONST char *norm, CONST int norm_len,
  CONST int *n,
  CONST Complex a[] /* 2 dimension */,
  CONST int *lda,
  double work[]
);

extern double __stdcall F06VGF(
  CONST char *norm, CONST int norm_len,
  CONST char *matrix, CONST int matrix_len,
  CONST int *m,
  CONST int *n,
  CONST Complex a[] /* 2 dimension */,
  CONST int *lda
);

extern void __stdcall F06VJF(
  CONST char *side, CONST int side_len,
  CONST char *trans, CONST int trans_len,
  CONST int *n,
  CONST int perm[],
  CONST int *k,
  Complex b[] /* 2 dimension */,
  CONST int *ldb
);

extern void __stdcall F06VKF(
  CONST char *side, CONST int side_len,
  CONST char *trans, CONST int trans_len,
  CONST int *n,
  CONST double perm[],
  CONST int *k,
  Complex b[] /* 2 dimension */,
  CONST int *ldb
);

extern void __stdcall F06VXF(
  CONST char *side, CONST int side_len,
  CONST char *pivot, CONST int pivot_len,
  CONST char *direct, CONST int direct_len,
  CONST int *m,
  CONST int *n,
  CONST int *k1,
  CONST int *k2,
  CONST double c[],
  CONST double s[],
  Complex a[] /* 2 dimension */,
  CONST int *lda
);

extern void __stdcall F06YAF(
  CONST char *transa, CONST int transa_len,
  CONST char *transb, CONST int transb_len,
  CONST int *m,
  CONST int *n,
  CONST int *k,
  CONST double *alpha,
  CONST double a[] /* 2 dimension */,
  CONST int *lda,
  CONST double b[] /* 2 dimension */,
  CONST int *ldb,
  CONST double *beta,
  double c[] /* 2 dimension */,
  CONST int *ldc
);

extern void __stdcall F06YCF(
  CONST char *side, CONST int side_len,
  CONST char *uplo, CONST int uplo_len,
  CONST int *m,
  CONST int *n,
  CONST double *alpha,
  CONST double a[] /* 2 dimension */,
  CONST int *lda,
  CONST double b[] /* 2 dimension */,
  CONST int *ldb,
  CONST double *beta,
  double c[] /* 2 dimension */,
  CONST int *ldc
);

extern void __stdcall F06YFF(
  CONST char *side, CONST int side_len,
  CONST char *uplo, CONST int uplo_len,
  CONST char *transa, CONST int transa_len,
  CONST char *diag, CONST int diag_len,
  CONST int *m,
  CONST int *n,
  CONST double *alpha,
  CONST double a[] /* 2 dimension */,
  CONST int *lda,
  double b[] /* 2 dimension */,
  CONST int *ldb
);

extern void __stdcall F06YJF(
  CONST char *side, CONST int side_len,
  CONST char *uplo, CONST int uplo_len,
  CONST char *transa, CONST int transa_len,
  CONST char *diag, CONST int diag_len,
  CONST int *m,
  CONST int *n,
  CONST double *alpha,
  CONST double a[] /* 2 dimension */,
  CONST int *lda,
  double b[] /* 2 dimension */,
  CONST int *ldb
);

extern void __stdcall F06YPF(
  CONST char *uplo, CONST int uplo_len,
  CONST char *trans, CONST int trans_len,
  CONST int *n,
  CONST int *k,
  CONST double *alpha,
  CONST double a[] /* 2 dimension */,
  CONST int *lda,
  CONST double *beta,
  double c[] /* 2 dimension */,
  CONST int *ldc
);

extern void __stdcall F06YRF(
  CONST char *uplo, CONST int uplo_len,
  CONST char *trans, CONST int trans_len,
  CONST int *n,
  CONST int *k,
  CONST double *alpha,
  CONST double a[] /* 2 dimension */,
  CONST int *lda,
  CONST double b[] /* 2 dimension */,
  CONST int *ldb,
  CONST double *beta,
  double c[] /* 2 dimension */,
  CONST int *ldc
);

extern void __stdcall F06ZAF(
  CONST char *transa, CONST int transa_len,
  CONST char *transb, CONST int transb_len,
  CONST int *m,
  CONST int *n,
  CONST int *k,
  CONST Complex *alpha,
  CONST Complex a[] /* 2 dimension */,
  CONST int *lda,
  CONST Complex b[] /* 2 dimension */,
  CONST int *ldb,
  CONST Complex *beta,
  Complex c[] /* 2 dimension */,
  CONST int *ldc
);

extern void __stdcall F06ZCF(
  CONST char *side, CONST int side_len,
  CONST char *uplo, CONST int uplo_len,
  CONST int *m,
  CONST int *n,
  CONST Complex *alpha,
  CONST Complex a[] /* 2 dimension */,
  CONST int *lda,
  CONST Complex b[] /* 2 dimension */,
  CONST int *ldb,
  CONST Complex *beta,
  Complex c[] /* 2 dimension */,
  CONST int *ldc
);

extern void __stdcall F06ZFF(
  CONST char *side, CONST int side_len,
  CONST char *uplo, CONST int uplo_len,
  CONST char *transa, CONST int transa_len,
  CONST char *diag, CONST int diag_len,
  CONST int *m,
  CONST int *n,
  CONST Complex *alpha,
  CONST Complex a[] /* 2 dimension */,
  CONST int *lda,
  Complex b[] /* 2 dimension */,
  CONST int *ldb
);

extern void __stdcall F06ZJF(
  CONST char *side, CONST int side_len,
  CONST char *uplo, CONST int uplo_len,
  CONST char *transa, CONST int transa_len,
  CONST char *diag, CONST int diag_len,
  CONST int *m,
  CONST int *n,
  CONST Complex *alpha,
  CONST Complex a[] /* 2 dimension */,
  CONST int *lda,
  Complex b[] /* 2 dimension */,
  CONST int *ldb
);

extern void __stdcall F06ZPF(
  CONST char *uplo, CONST int uplo_len,
  CONST char *trans, CONST int trans_len,
  CONST int *n,
  CONST int *k,
  CONST double *alpha,
  CONST Complex a[] /* 2 dimension */,
  CONST int *lda,
  CONST double *beta,
  Complex c[] /* 2 dimension */,
  CONST int *ldc
);

extern void __stdcall F06ZRF(
  CONST char *uplo, CONST int uplo_len,
  CONST char *trans, CONST int trans_len,
  CONST int *n,
  CONST int *k,
  CONST Complex *alpha,
  CONST Complex a[] /* 2 dimension */,
  CONST int *lda,
  CONST Complex b[] /* 2 dimension */,
  CONST int *ldb,
  CONST double *beta,
  Complex c[] /* 2 dimension */,
  CONST int *ldc
);

extern void __stdcall F06ZTF(
  CONST char *side, CONST int side_len,
  CONST char *uplo, CONST int uplo_len,
  CONST int *m,
  CONST int *n,
  CONST Complex *alpha,
  CONST Complex a[] /* 2 dimension */,
  CONST int *lda,
  CONST Complex b[] /* 2 dimension */,
  CONST int *ldb,
  CONST Complex *beta,
  Complex c[] /* 2 dimension */,
  CONST int *ldc
);

extern void __stdcall F06ZUF(
  CONST char *uplo, CONST int uplo_len,
  CONST char *trans, CONST int trans_len,
  CONST int *n,
  CONST int *k,
  CONST Complex *alpha,
  CONST Complex a[] /* 2 dimension */,
  CONST int *lda,
  CONST Complex *beta,
  Complex c[] /* 2 dimension */,
  CONST int *ldc
);

extern void __stdcall F06ZWF(
  CONST char *uplo, CONST int uplo_len,
  CONST char *trans, CONST int trans_len,
  CONST int *n,
  CONST int *k,
  CONST Complex *alpha,
  CONST Complex a[] /* 2 dimension */,
  CONST int *lda,
  CONST Complex b[] /* 2 dimension */,
  CONST int *ldb,
  CONST Complex *beta,
  Complex c[] /* 2 dimension */,
  CONST int *ldc
);

extern void __stdcall F07ADF(
  CONST int *m,
  CONST int *n,
  double a[] /* 2 dimension */,
  CONST int *lda,
  int ipiv[],
  int *info
);

extern void __stdcall dgemm(
  CONST char *transa, CONST int transa_len,
  CONST char *transb, CONST int transb_len,
  CONST int *m,
  CONST int *n,
  CONST int *k,
  CONST double *alpha,
  CONST double a[] /* 2 dimension */,
  CONST int *lda,
  CONST double b[] /* 2 dimension */,
  CONST int *ldb,
  CONST double *beta,
  double c[] /* 2 dimension */,
  CONST int *ldc
);

extern void __stdcall dsymm(
  CONST char *side, CONST int side_len,
  CONST char *uplo, CONST int uplo_len,
  CONST int *m,
  CONST int *n,
  CONST double *alpha,
  CONST double a[] /* 2 dimension */,
  CONST int *lda,
  CONST double b[] /* 2 dimension */,
  CONST int *ldb,
  CONST double *beta,
  double c[] /* 2 dimension */,
  CONST int *ldc
);

extern void __stdcall dtrmm(
  CONST char *side, CONST int side_len,
  CONST char *uplo, CONST int uplo_len,
  CONST char *transa, CONST int transa_len,
  CONST char *diag, CONST int diag_len,
  CONST int *m,
  CONST int *n,
  CONST double *alpha,
  CONST double a[] /* 2 dimension */,
  CONST int *lda,
  double b[] /* 2 dimension */,
  CONST int *ldb
);

extern void __stdcall dtrsm(
  CONST char *side, CONST int side_len,
  CONST char *uplo, CONST int uplo_len,
  CONST char *transa, CONST int transa_len,
  CONST char *diag, CONST int diag_len,
  CONST int *m,
  CONST int *n,
  CONST double *alpha,
  CONST double a[] /* 2 dimension */,
  CONST int *lda,
  double b[] /* 2 dimension */,
  CONST int *ldb
);

extern void __stdcall dsyrk(
  CONST char *uplo, CONST int uplo_len,
  CONST char *trans, CONST int trans_len,
  CONST int *n,
  CONST int *k,
  CONST double *alpha,
  CONST double a[] /* 2 dimension */,
  CONST int *lda,
  CONST double *beta,
  double c[] /* 2 dimension */,
  CONST int *ldc
);

extern void __stdcall dsyr2k(
  CONST char *uplo, CONST int uplo_len,
  CONST char *trans, CONST int trans_len,
  CONST int *n,
  CONST int *k,
  CONST double *alpha,
  CONST double a[] /* 2 dimension */,
  CONST int *lda,
  CONST double b[] /* 2 dimension */,
  CONST int *ldb,
  CONST double *beta,
  double c[] /* 2 dimension */,
  CONST int *ldc
);

extern void __stdcall zgemm(
  CONST char *transa, CONST int transa_len,
  CONST char *transb, CONST int transb_len,
  CONST int *m,
  CONST int *n,
  CONST int *k,
  CONST Complex *alpha,
  CONST Complex a[] /* 2 dimension */,
  CONST int *lda,
  CONST Complex b[] /* 2 dimension */,
  CONST int *ldb,
  CONST Complex *beta,
  Complex c[] /* 2 dimension */,
  CONST int *ldc
);

extern void __stdcall zhemm(
  CONST char *side, CONST int side_len,
  CONST char *uplo, CONST int uplo_len,
  CONST int *m,
  CONST int *n,
  CONST Complex *alpha,
  CONST Complex a[] /* 2 dimension */,
  CONST int *lda,
  CONST Complex b[] /* 2 dimension */,
  CONST int *ldb,
  CONST Complex *beta,
  Complex c[] /* 2 dimension */,
  CONST int *ldc
);

extern void __stdcall ztrmm(
  CONST char *side, CONST int side_len,
  CONST char *uplo, CONST int uplo_len,
  CONST char *transa, CONST int transa_len,
  CONST char *diag, CONST int diag_len,
  CONST int *m,
  CONST int *n,
  CONST Complex *alpha,
  CONST Complex a[] /* 2 dimension */,
  CONST int *lda,
  Complex b[] /* 2 dimension */,
  CONST int *ldb
);

extern void __stdcall ztrsm(
  CONST char *side, CONST int side_len,
  CONST char *uplo, CONST int uplo_len,
  CONST char *transa, CONST int transa_len,
  CONST char *diag, CONST int diag_len,
  CONST int *m,
  CONST int *n,
  CONST Complex *alpha,
  CONST Complex a[] /* 2 dimension */,
  CONST int *lda,
  Complex b[] /* 2 dimension */,
  CONST int *ldb
);

extern void __stdcall zherk(
  CONST char *uplo, CONST int uplo_len,
  CONST char *trans, CONST int trans_len,
  CONST int *n,
  CONST int *k,
  CONST double *alpha,
  CONST Complex a[] /* 2 dimension */,
  CONST int *lda,
  CONST double *beta,
  Complex c[] /* 2 dimension */,
  CONST int *ldc
);

extern void __stdcall zher2k(
  CONST char *uplo, CONST int uplo_len,
  CONST char *trans, CONST int trans_len,
  CONST int *n,
  CONST int *k,
  CONST Complex *alpha,
  CONST Complex a[] /* 2 dimension */,
  CONST int *lda,
  CONST Complex b[] /* 2 dimension */,
  CONST int *ldb,
  CONST double *beta,
  Complex c[] /* 2 dimension */,
  CONST int *ldc
);

extern void __stdcall zsymm(
  CONST char *side, CONST int side_len,
  CONST char *uplo, CONST int uplo_len,
  CONST int *m,
  CONST int *n,
  CONST Complex *alpha,
  CONST Complex a[] /* 2 dimension */,
  CONST int *lda,
  CONST Complex b[] /* 2 dimension */,
  CONST int *ldb,
  CONST Complex *beta,
  Complex c[] /* 2 dimension */,
  CONST int *ldc
);

extern void __stdcall zsyrk(
  CONST char *uplo, CONST int uplo_len,
  CONST char *trans, CONST int trans_len,
  CONST int *n,
  CONST int *k,
  CONST Complex *alpha,
  CONST Complex a[] /* 2 dimension */,
  CONST int *lda,
  CONST Complex *beta,
  Complex c[] /* 2 dimension */,
  CONST int *ldc
);

extern void __stdcall zsyr2k(
  CONST char *uplo, CONST int uplo_len,
  CONST char *trans, CONST int trans_len,
  CONST int *n,
  CONST int *k,
  CONST Complex *alpha,
  CONST Complex a[] /* 2 dimension */,
  CONST int *lda,
  CONST Complex b[] /* 2 dimension */,
  CONST int *ldb,
  CONST Complex *beta,
  Complex c[] /* 2 dimension */,
  CONST int *ldc
);

extern void __stdcall dgetrf(
  CONST int *m,
  CONST int *n,
  double a[] /* 2 dimension */,
  CONST int *lda,
  int ipiv[],
  int *info
);

extern void __stdcall F07AEF(
  CONST char *trans, CONST int trans_len,
  CONST int *n,
  CONST int *nrhs,
  CONST double a[] /* 2 dimension */,
  CONST int *lda,
  CONST int ipiv[],
  double b[] /* 2 dimension */,
  CONST int *ldb,
  int *info
);

extern void __stdcall F07AGF(
  CONST char *norm, CONST int norm_len,
  CONST int *n,
  CONST double a[] /* 2 dimension */,
  CONST int *lda,
  CONST double *anorm,
  double *rcond,
  double work[],
  int iwork[],
  int *info
);

extern void __stdcall F07AHF(
  CONST char *trans, CONST int trans_len,
  CONST int *n,
  CONST int *nrhs,
  CONST double a[] /* 2 dimension */,
  CONST int *lda,
  CONST double af[] /* 2 dimension */,
  CONST int *ldaf,
  CONST int ipiv[],
  CONST double b[] /* 2 dimension */,
  CONST int *ldb,
  double x[] /* 2 dimension */,
  CONST int *ldx,
  double ferr[],
  double berr[],
  double work[],
  int iwork[],
  int *info
);

extern void __stdcall F07AJF(
  CONST int *n,
  double a[] /* 2 dimension */,
  CONST int *lda,
  CONST int ipiv[],
  double work[],
  CONST int *lwork,
  int *info
);

extern void __stdcall F07ARF(
  CONST int *m,
  CONST int *n,
  Complex a[] /* 2 dimension */,
  CONST int *lda,
  int ipiv[],
  int *info
);

extern void __stdcall F07ASF(
  CONST char *trans, CONST int trans_len,
  CONST int *n,
  CONST int *nrhs,
  CONST Complex a[] /* 2 dimension */,
  CONST int *lda,
  CONST int ipiv[],
  Complex b[] /* 2 dimension */,
  CONST int *ldb,
  int *info
);

extern void __stdcall F07AUF(
  CONST char *norm, CONST int norm_len,
  CONST int *n,
  CONST Complex a[] /* 2 dimension */,
  CONST int *lda,
  CONST double *anorm,
  double *rcond,
  Complex work[],
  double rwork[],
  int *info
);

extern void __stdcall F07AVF(
  CONST char *trans, CONST int trans_len,
  CONST int *n,
  CONST int *nrhs,
  CONST Complex a[] /* 2 dimension */,
  CONST int *lda,
  CONST Complex af[] /* 2 dimension */,
  CONST int *ldaf,
  CONST int ipiv[],
  CONST Complex b[] /* 2 dimension */,
  CONST int *ldb,
  Complex x[] /* 2 dimension */,
  CONST int *ldx,
  double ferr[],
  double berr[],
  Complex work[],
  double rwork[],
  int *info
);

extern void __stdcall F07AWF(
  CONST int *n,
  Complex a[] /* 2 dimension */,
  CONST int *lda,
  CONST int ipiv[],
  Complex work[],
  CONST int *lwork,
  int *info
);

extern void __stdcall F07BDF(
  CONST int *m,
  CONST int *n,
  CONST int *kl,
  CONST int *ku,
  double ab[] /* 2 dimension */,
  CONST int *ldab,
  int ipiv[],
  int *info
);

extern void __stdcall F07BEF(
  CONST char *trans, CONST int trans_len,
  CONST int *n,
  CONST int *kl,
  CONST int *ku,
  CONST int *nrhs,
  CONST double ab[] /* 2 dimension */,
  CONST int *ldab,
  CONST int ipiv[],
  double b[] /* 2 dimension */,
  CONST int *ldb,
  int *info
);

extern void __stdcall F07BGF(
  CONST char *norm, CONST int norm_len,
  CONST int *n,
  CONST int *kl,
  CONST int *ku,
  CONST double ab[] /* 2 dimension */,
  CONST int *ldab,
  CONST int ipiv[],
  CONST double *anorm,
  double *rcond,
  double work[],
  int iwork[],
  int *info
);

extern void __stdcall F07BHF(
  CONST char *trans, CONST int trans_len,
  CONST int *n,
  CONST int *kl,
  CONST int *ku,
  CONST int *nrhs,
  CONST double ab[] /* 2 dimension */,
  CONST int *ldab,
  CONST double afb[] /* 2 dimension */,
  CONST int *ldafb,
  CONST int ipiv[],
  CONST double b[] /* 2 dimension */,
  CONST int *ldb,
  double x[] /* 2 dimension */,
  CONST int *ldx,
  double ferr[],
  double berr[],
  double work[],
  int iwork[],
  int *info
);

extern void __stdcall F07BRF(
  CONST int *m,
  CONST int *n,
  CONST int *kl,
  CONST int *ku,
  Complex ab[] /* 2 dimension */,
  CONST int *ldab,
  int ipiv[],
  int *info
);

extern void __stdcall F07BSF(
  CONST char *trans, CONST int trans_len,
  CONST int *n,
  CONST int *kl,
  CONST int *ku,
  CONST int *nrhs,
  CONST Complex ab[] /* 2 dimension */,
  CONST int *ldab,
  CONST int ipiv[],
  Complex b[] /* 2 dimension */,
  CONST int *ldb,
  int *info
);

extern void __stdcall F07BUF(
  CONST char *norm, CONST int norm_len,
  CONST int *n,
  CONST int *kl,
  CONST int *ku,
  CONST Complex ab[] /* 2 dimension */,
  CONST int *ldab,
  CONST int ipiv[],
  CONST double *anorm,
  double *rcond,
  Complex work[],
  double rwork[],
  int *info
);

extern void __stdcall F07BVF(
  CONST char *trans, CONST int trans_len,
  CONST int *n,
  CONST int *kl,
  CONST int *ku,
  CONST int *nrhs,
  CONST Complex ab[] /* 2 dimension */,
  CONST int *ldab,
  CONST Complex afb[] /* 2 dimension */,
  CONST int *ldafb,
  CONST int ipiv[],
  CONST Complex b[] /* 2 dimension */,
  CONST int *ldb,
  Complex x[] /* 2 dimension */,
  CONST int *ldx,
  double ferr[],
  double berr[],
  Complex work[],
  double rwork[],
  int *info
);

extern void __stdcall F07FDF(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  double a[] /* 2 dimension */,
  CONST int *lda,
  int *info
);

extern void __stdcall dgetrs(
  CONST char *trans, CONST int trans_len,
  CONST int *n,
  CONST int *nrhs,
  CONST double a[] /* 2 dimension */,
  CONST int *lda,
  CONST int ipiv[],
  double b[] /* 2 dimension */,
  CONST int *ldb,
  int *info
);

extern void __stdcall dgecon(
  CONST char *norm, CONST int norm_len,
  CONST int *n,
  CONST double a[] /* 2 dimension */,
  CONST int *lda,
  CONST double *anorm,
  double *rcond,
  double work[],
  int iwork[],
  int *info
);

extern void __stdcall dgerfs(
  CONST char *trans, CONST int trans_len,
  CONST int *n,
  CONST int *nrhs,
  CONST double a[] /* 2 dimension */,
  CONST int *lda,
  CONST double af[] /* 2 dimension */,
  CONST int *ldaf,
  CONST int ipiv[],
  CONST double b[] /* 2 dimension */,
  CONST int *ldb,
  double x[] /* 2 dimension */,
  CONST int *ldx,
  double ferr[],
  double berr[],
  double work[],
  int iwork[],
  int *info
);

extern void __stdcall dgetri(
  CONST int *n,
  double a[] /* 2 dimension */,
  CONST int *lda,
  CONST int ipiv[],
  double work[],
  CONST int *lwork,
  int *info
);

extern void __stdcall zgetrf(
  CONST int *m,
  CONST int *n,
  Complex a[] /* 2 dimension */,
  CONST int *lda,
  int ipiv[],
  int *info
);

extern void __stdcall zgetrs(
  CONST char *trans, CONST int trans_len,
  CONST int *n,
  CONST int *nrhs,
  CONST Complex a[] /* 2 dimension */,
  CONST int *lda,
  CONST int ipiv[],
  Complex b[] /* 2 dimension */,
  CONST int *ldb,
  int *info
);

extern void __stdcall zgecon(
  CONST char *norm, CONST int norm_len,
  CONST int *n,
  CONST Complex a[] /* 2 dimension */,
  CONST int *lda,
  CONST double *anorm,
  double *rcond,
  Complex work[],
  double rwork[],
  int *info
);

extern void __stdcall zgerfs(
  CONST char *trans, CONST int trans_len,
  CONST int *n,
  CONST int *nrhs,
  CONST Complex a[] /* 2 dimension */,
  CONST int *lda,
  CONST Complex af[] /* 2 dimension */,
  CONST int *ldaf,
  CONST int ipiv[],
  CONST Complex b[] /* 2 dimension */,
  CONST int *ldb,
  Complex x[] /* 2 dimension */,
  CONST int *ldx,
  double ferr[],
  double berr[],
  Complex work[],
  double rwork[],
  int *info
);

extern void __stdcall zgetri(
  CONST int *n,
  Complex a[] /* 2 dimension */,
  CONST int *lda,
  CONST int ipiv[],
  Complex work[],
  CONST int *lwork,
  int *info
);

extern void __stdcall dgbtrf(
  CONST int *m,
  CONST int *n,
  CONST int *kl,
  CONST int *ku,
  double ab[] /* 2 dimension */,
  CONST int *ldab,
  int ipiv[],
  int *info
);

extern void __stdcall dgbtrs(
  CONST char *trans, CONST int trans_len,
  CONST int *n,
  CONST int *kl,
  CONST int *ku,
  CONST int *nrhs,
  CONST double ab[] /* 2 dimension */,
  CONST int *ldab,
  CONST int ipiv[],
  double b[] /* 2 dimension */,
  CONST int *ldb,
  int *info
);

extern void __stdcall dgbcon(
  CONST char *norm, CONST int norm_len,
  CONST int *n,
  CONST int *kl,
  CONST int *ku,
  CONST double ab[] /* 2 dimension */,
  CONST int *ldab,
  CONST int ipiv[],
  CONST double *anorm,
  double *rcond,
  double work[],
  int iwork[],
  int *info
);

extern void __stdcall dgbrfs(
  CONST char *trans, CONST int trans_len,
  CONST int *n,
  CONST int *kl,
  CONST int *ku,
  CONST int *nrhs,
  CONST double ab[] /* 2 dimension */,
  CONST int *ldab,
  CONST double afb[] /* 2 dimension */,
  CONST int *ldafb,
  CONST int ipiv[],
  CONST double b[] /* 2 dimension */,
  CONST int *ldb,
  double x[] /* 2 dimension */,
  CONST int *ldx,
  double ferr[],
  double berr[],
  double work[],
  int iwork[],
  int *info
);

extern void __stdcall zgbtrf(
  CONST int *m,
  CONST int *n,
  CONST int *kl,
  CONST int *ku,
  Complex ab[] /* 2 dimension */,
  CONST int *ldab,
  int ipiv[],
  int *info
);

extern void __stdcall zgbtrs(
  CONST char *trans, CONST int trans_len,
  CONST int *n,
  CONST int *kl,
  CONST int *ku,
  CONST int *nrhs,
  CONST Complex ab[] /* 2 dimension */,
  CONST int *ldab,
  CONST int ipiv[],
  Complex b[] /* 2 dimension */,
  CONST int *ldb,
  int *info
);

extern void __stdcall zgbcon(
  CONST char *norm, CONST int norm_len,
  CONST int *n,
  CONST int *kl,
  CONST int *ku,
  CONST Complex ab[] /* 2 dimension */,
  CONST int *ldab,
  CONST int ipiv[],
  CONST double *anorm,
  double *rcond,
  Complex work[],
  double rwork[],
  int *info
);

extern void __stdcall zgbrfs(
  CONST char *trans, CONST int trans_len,
  CONST int *n,
  CONST int *kl,
  CONST int *ku,
  CONST int *nrhs,
  CONST Complex ab[] /* 2 dimension */,
  CONST int *ldab,
  CONST Complex afb[] /* 2 dimension */,
  CONST int *ldafb,
  CONST int ipiv[],
  CONST Complex b[] /* 2 dimension */,
  CONST int *ldb,
  Complex x[] /* 2 dimension */,
  CONST int *ldx,
  double ferr[],
  double berr[],
  Complex work[],
  double rwork[],
  int *info
);

extern void __stdcall dpotrf(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  double a[] /* 2 dimension */,
  CONST int *lda,
  int *info
);

extern void __stdcall F07FEF(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST int *nrhs,
  CONST double a[] /* 2 dimension */,
  CONST int *lda,
  double b[] /* 2 dimension */,
  CONST int *ldb,
  int *info
);

extern void __stdcall F07FGF(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST double a[] /* 2 dimension */,
  CONST int *lda,
  CONST double *anorm,
  double *rcond,
  double work[],
  int iwork[],
  int *info
);

extern void __stdcall F07FHF(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST int *nrhs,
  CONST double a[] /* 2 dimension */,
  CONST int *lda,
  CONST double af[] /* 2 dimension */,
  CONST int *ldaf,
  CONST double b[] /* 2 dimension */,
  CONST int *ldb,
  double x[] /* 2 dimension */,
  CONST int *ldx,
  double ferr[],
  double berr[],
  double work[],
  int iwork[],
  int *info
);

extern void __stdcall F07FJF(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  double a[] /* 2 dimension */,
  CONST int *lda,
  int *info
);

extern void __stdcall F07FRF(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  Complex a[] /* 2 dimension */,
  CONST int *lda,
  int *info
);

extern void __stdcall F07FSF(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST int *nrhs,
  CONST Complex a[] /* 2 dimension */,
  CONST int *lda,
  Complex b[] /* 2 dimension */,
  CONST int *ldb,
  int *info
);

extern void __stdcall F07FUF(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST Complex a[] /* 2 dimension */,
  CONST int *lda,
  CONST double *anorm,
  double *rcond,
  Complex work[],
  double rwork[],
  int *info
);

extern void __stdcall F07FVF(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST int *nrhs,
  CONST Complex a[] /* 2 dimension */,
  CONST int *lda,
  CONST Complex af[] /* 2 dimension */,
  CONST int *ldaf,
  CONST Complex b[] /* 2 dimension */,
  CONST int *ldb,
  Complex x[] /* 2 dimension */,
  CONST int *ldx,
  double ferr[],
  double berr[],
  Complex work[],
  double rwork[],
  int *info
);

extern void __stdcall F07FWF(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  Complex a[] /* 2 dimension */,
  CONST int *lda,
  int *info
);

extern void __stdcall F07GDF(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  double ap[],
  int *info
);

extern void __stdcall F07GEF(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST int *nrhs,
  CONST double ap[],
  double b[] /* 2 dimension */,
  CONST int *ldb,
  int *info
);

extern void __stdcall F07GGF(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST double ap[],
  CONST double *anorm,
  double *rcond,
  double work[],
  int iwork[],
  int *info
);

extern void __stdcall F07GHF(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST int *nrhs,
  CONST double ap[],
  CONST double afp[],
  CONST double b[] /* 2 dimension */,
  CONST int *ldb,
  double x[] /* 2 dimension */,
  CONST int *ldx,
  double ferr[],
  double berr[],
  double work[],
  int iwork[],
  int *info
);

extern void __stdcall F07GJF(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  double ap[],
  int *info
);

extern void __stdcall F07GRF(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  Complex ap[],
  int *info
);

extern void __stdcall F07GSF(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST int *nrhs,
  CONST Complex ap[],
  Complex b[] /* 2 dimension */,
  CONST int *ldb,
  int *info
);

extern void __stdcall F07GUF(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST Complex ap[],
  CONST double *anorm,
  double *rcond,
  Complex work[],
  double rwork[],
  int *info
);

extern void __stdcall F07GVF(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST int *nrhs,
  CONST Complex ap[],
  CONST Complex afp[],
  CONST Complex b[] /* 2 dimension */,
  CONST int *ldb,
  Complex x[] /* 2 dimension */,
  CONST int *ldx,
  double ferr[],
  double berr[],
  Complex work[],
  double rwork[],
  int *info
);

extern void __stdcall F07GWF(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  Complex ap[],
  int *info
);

extern void __stdcall F07HDF(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST int *kd,
  double ab[] /* 2 dimension */,
  CONST int *ldab,
  int *info
);

extern void __stdcall F07HEF(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST int *kd,
  CONST int *nrhs,
  CONST double ab[] /* 2 dimension */,
  CONST int *ldab,
  double b[] /* 2 dimension */,
  CONST int *ldb,
  int *info
);

extern void __stdcall F07HGF(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST int *kd,
  CONST double ab[] /* 2 dimension */,
  CONST int *ldab,
  CONST double *anorm,
  double *rcond,
  double work[],
  int iwork[],
  int *info
);

extern void __stdcall F07HHF(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST int *kd,
  CONST int *nrhs,
  CONST double ab[] /* 2 dimension */,
  CONST int *ldab,
  CONST double afb[] /* 2 dimension */,
  CONST int *ldafb,
  CONST double b[] /* 2 dimension */,
  CONST int *ldb,
  double x[] /* 2 dimension */,
  CONST int *ldx,
  double ferr[],
  double berr[],
  double work[],
  int iwork[],
  int *info
);

extern void __stdcall dpotrs(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST int *nrhs,
  CONST double a[] /* 2 dimension */,
  CONST int *lda,
  double b[] /* 2 dimension */,
  CONST int *ldb,
  int *info
);

extern void __stdcall dpocon(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST double a[] /* 2 dimension */,
  CONST int *lda,
  CONST double *anorm,
  double *rcond,
  double work[],
  int iwork[],
  int *info
);

extern void __stdcall dporfs(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST int *nrhs,
  CONST double a[] /* 2 dimension */,
  CONST int *lda,
  CONST double af[] /* 2 dimension */,
  CONST int *ldaf,
  CONST double b[] /* 2 dimension */,
  CONST int *ldb,
  double x[] /* 2 dimension */,
  CONST int *ldx,
  double ferr[],
  double berr[],
  double work[],
  int iwork[],
  int *info
);

extern void __stdcall dpotri(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  double a[] /* 2 dimension */,
  CONST int *lda,
  int *info
);

extern void __stdcall zpotrf(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  Complex a[] /* 2 dimension */,
  CONST int *lda,
  int *info
);

extern void __stdcall zpotrs(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST int *nrhs,
  CONST Complex a[] /* 2 dimension */,
  CONST int *lda,
  Complex b[] /* 2 dimension */,
  CONST int *ldb,
  int *info
);

extern void __stdcall zpocon(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST Complex a[] /* 2 dimension */,
  CONST int *lda,
  CONST double *anorm,
  double *rcond,
  Complex work[],
  double rwork[],
  int *info
);

extern void __stdcall zporfs(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST int *nrhs,
  CONST Complex a[] /* 2 dimension */,
  CONST int *lda,
  CONST Complex af[] /* 2 dimension */,
  CONST int *ldaf,
  CONST Complex b[] /* 2 dimension */,
  CONST int *ldb,
  Complex x[] /* 2 dimension */,
  CONST int *ldx,
  double ferr[],
  double berr[],
  Complex work[],
  double rwork[],
  int *info
);

extern void __stdcall zpotri(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  Complex a[] /* 2 dimension */,
  CONST int *lda,
  int *info
);

extern void __stdcall dpptrf(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  double ap[],
  int *info
);

extern void __stdcall dpptrs(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST int *nrhs,
  CONST double ap[],
  double b[] /* 2 dimension */,
  CONST int *ldb,
  int *info
);

extern void __stdcall dppcon(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST double ap[],
  CONST double *anorm,
  double *rcond,
  double work[],
  int iwork[],
  int *info
);

extern void __stdcall dpprfs(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST int *nrhs,
  CONST double ap[],
  CONST double afp[],
  CONST double b[] /* 2 dimension */,
  CONST int *ldb,
  double x[] /* 2 dimension */,
  CONST int *ldx,
  double ferr[],
  double berr[],
  double work[],
  int iwork[],
  int *info
);

extern void __stdcall dpptri(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  double ap[],
  int *info
);

extern void __stdcall zpptrf(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  Complex ap[],
  int *info
);

extern void __stdcall zpptrs(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST int *nrhs,
  CONST Complex ap[],
  Complex b[] /* 2 dimension */,
  CONST int *ldb,
  int *info
);

extern void __stdcall zppcon(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST Complex ap[],
  CONST double *anorm,
  double *rcond,
  Complex work[],
  double rwork[],
  int *info
);

extern void __stdcall zpprfs(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST int *nrhs,
  CONST Complex ap[],
  CONST Complex afp[],
  CONST Complex b[] /* 2 dimension */,
  CONST int *ldb,
  Complex x[] /* 2 dimension */,
  CONST int *ldx,
  double ferr[],
  double berr[],
  Complex work[],
  double rwork[],
  int *info
);

extern void __stdcall zpptri(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  Complex ap[],
  int *info
);

extern void __stdcall dpbtrf(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST int *kd,
  double ab[] /* 2 dimension */,
  CONST int *ldab,
  int *info
);

extern void __stdcall dpbtrs(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST int *kd,
  CONST int *nrhs,
  CONST double ab[] /* 2 dimension */,
  CONST int *ldab,
  double b[] /* 2 dimension */,
  CONST int *ldb,
  int *info
);

extern void __stdcall dpbcon(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST int *kd,
  CONST double ab[] /* 2 dimension */,
  CONST int *ldab,
  CONST double *anorm,
  double *rcond,
  double work[],
  int iwork[],
  int *info
);

extern void __stdcall dpbrfs(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST int *kd,
  CONST int *nrhs,
  CONST double ab[] /* 2 dimension */,
  CONST int *ldab,
  CONST double afb[] /* 2 dimension */,
  CONST int *ldafb,
  CONST double b[] /* 2 dimension */,
  CONST int *ldb,
  double x[] /* 2 dimension */,
  CONST int *ldx,
  double ferr[],
  double berr[],
  double work[],
  int iwork[],
  int *info
);

extern void __stdcall F07HRF(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST int *kd,
  Complex ab[] /* 2 dimension */,
  CONST int *ldab,
  int *info
);

extern void __stdcall F07HSF(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST int *kd,
  CONST int *nrhs,
  CONST Complex ab[] /* 2 dimension */,
  CONST int *ldab,
  Complex b[] /* 2 dimension */,
  CONST int *ldb,
  int *info
);

extern void __stdcall F07HUF(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST int *kd,
  CONST Complex ab[] /* 2 dimension */,
  CONST int *ldab,
  CONST double *anorm,
  double *rcond,
  Complex work[],
  double rwork[],
  int *info
);

extern void __stdcall F07HVF(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST int *kd,
  CONST int *nrhs,
  CONST Complex ab[] /* 2 dimension */,
  CONST int *ldab,
  CONST Complex afb[] /* 2 dimension */,
  CONST int *ldafb,
  CONST Complex b[] /* 2 dimension */,
  CONST int *ldb,
  Complex x[] /* 2 dimension */,
  CONST int *ldx,
  double ferr[],
  double berr[],
  Complex work[],
  double rwork[],
  int *info
);

extern void __stdcall F07MDF(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  double a[] /* 2 dimension */,
  CONST int *lda,
  int ipiv[],
  double work[],
  CONST int *lwork,
  int *info
);

extern void __stdcall F07MEF(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST int *nrhs,
  CONST double a[] /* 2 dimension */,
  CONST int *lda,
  CONST int ipiv[],
  double b[] /* 2 dimension */,
  CONST int *ldb,
  int *info
);

extern void __stdcall F07MGF(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST double a[] /* 2 dimension */,
  CONST int *lda,
  CONST int ipiv[],
  CONST double *anorm,
  double *rcond,
  double work[],
  int iwork[],
  int *info
);

extern void __stdcall F07MHF(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST int *nrhs,
  CONST double a[] /* 2 dimension */,
  CONST int *lda,
  CONST double af[] /* 2 dimension */,
  CONST int *ldaf,
  CONST int ipiv[],
  CONST double b[] /* 2 dimension */,
  CONST int *ldb,
  double x[] /* 2 dimension */,
  CONST int *ldx,
  double ferr[],
  double berr[],
  double work[],
  int iwork[],
  int *info
);

extern void __stdcall F07MJF(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  double a[] /* 2 dimension */,
  CONST int *lda,
  CONST int ipiv[],
  double work[],
  int *info
);

extern void __stdcall F07MRF(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  Complex a[] /* 2 dimension */,
  CONST int *lda,
  int ipiv[],
  Complex work[],
  CONST int *lwork,
  int *info
);

extern void __stdcall F07MSF(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST int *nrhs,
  CONST Complex a[] /* 2 dimension */,
  CONST int *lda,
  CONST int ipiv[],
  Complex b[] /* 2 dimension */,
  CONST int *ldb,
  int *info
);

extern void __stdcall F07MUF(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST Complex a[] /* 2 dimension */,
  CONST int *lda,
  CONST int ipiv[],
  CONST double *anorm,
  double *rcond,
  Complex work[],
  int *info
);

extern void __stdcall F07MVF(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST int *nrhs,
  CONST Complex a[] /* 2 dimension */,
  CONST int *lda,
  CONST Complex af[] /* 2 dimension */,
  CONST int *ldaf,
  CONST int ipiv[],
  CONST Complex b[] /* 2 dimension */,
  CONST int *ldb,
  Complex x[] /* 2 dimension */,
  CONST int *ldx,
  double ferr[],
  double berr[],
  Complex work[],
  double rwork[],
  int *info
);

extern void __stdcall F07MWF(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  Complex a[] /* 2 dimension */,
  CONST int *lda,
  CONST int ipiv[],
  Complex work[],
  int *info
);

extern void __stdcall F07NRF(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  Complex a[] /* 2 dimension */,
  CONST int *lda,
  int ipiv[],
  Complex work[],
  CONST int *lwork,
  int *info
);

extern void __stdcall F07NSF(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST int *nrhs,
  CONST Complex a[] /* 2 dimension */,
  CONST int *lda,
  CONST int ipiv[],
  Complex b[] /* 2 dimension */,
  CONST int *ldb,
  int *info
);

extern void __stdcall F07NUF(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST Complex a[] /* 2 dimension */,
  CONST int *lda,
  CONST int ipiv[],
  CONST double *anorm,
  double *rcond,
  Complex work[],
  int *info
);

extern void __stdcall zpbtrf(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST int *kd,
  Complex ab[] /* 2 dimension */,
  CONST int *ldab,
  int *info
);

extern void __stdcall zpbtrs(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST int *kd,
  CONST int *nrhs,
  CONST Complex ab[] /* 2 dimension */,
  CONST int *ldab,
  Complex b[] /* 2 dimension */,
  CONST int *ldb,
  int *info
);

extern void __stdcall zpbcon(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST int *kd,
  CONST Complex ab[] /* 2 dimension */,
  CONST int *ldab,
  CONST double *anorm,
  double *rcond,
  Complex work[],
  double rwork[],
  int *info
);

extern void __stdcall zpbrfs(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST int *kd,
  CONST int *nrhs,
  CONST Complex ab[] /* 2 dimension */,
  CONST int *ldab,
  CONST Complex afb[] /* 2 dimension */,
  CONST int *ldafb,
  CONST Complex b[] /* 2 dimension */,
  CONST int *ldb,
  Complex x[] /* 2 dimension */,
  CONST int *ldx,
  double ferr[],
  double berr[],
  Complex work[],
  double rwork[],
  int *info
);

extern void __stdcall dsytrf(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  double a[] /* 2 dimension */,
  CONST int *lda,
  int ipiv[],
  double work[],
  CONST int *lwork,
  int *info
);

extern void __stdcall dsytrs(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST int *nrhs,
  CONST double a[] /* 2 dimension */,
  CONST int *lda,
  CONST int ipiv[],
  double b[] /* 2 dimension */,
  CONST int *ldb,
  int *info
);

extern void __stdcall dsycon(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST double a[] /* 2 dimension */,
  CONST int *lda,
  CONST int ipiv[],
  CONST double *anorm,
  double *rcond,
  double work[],
  int iwork[],
  int *info
);

extern void __stdcall dsyrfs(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST int *nrhs,
  CONST double a[] /* 2 dimension */,
  CONST int *lda,
  CONST double af[] /* 2 dimension */,
  CONST int *ldaf,
  CONST int ipiv[],
  CONST double b[] /* 2 dimension */,
  CONST int *ldb,
  double x[] /* 2 dimension */,
  CONST int *ldx,
  double ferr[],
  double berr[],
  double work[],
  int iwork[],
  int *info
);

extern void __stdcall dsytri(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  double a[] /* 2 dimension */,
  CONST int *lda,
  CONST int ipiv[],
  double work[],
  int *info
);

extern void __stdcall zhetrf(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  Complex a[] /* 2 dimension */,
  CONST int *lda,
  int ipiv[],
  Complex work[],
  CONST int *lwork,
  int *info
);

extern void __stdcall zhetrs(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST int *nrhs,
  CONST Complex a[] /* 2 dimension */,
  CONST int *lda,
  CONST int ipiv[],
  Complex b[] /* 2 dimension */,
  CONST int *ldb,
  int *info
);

extern void __stdcall zhecon(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST Complex a[] /* 2 dimension */,
  CONST int *lda,
  CONST int ipiv[],
  CONST double *anorm,
  double *rcond,
  Complex work[],
  int *info
);

extern void __stdcall zherfs(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST int *nrhs,
  CONST Complex a[] /* 2 dimension */,
  CONST int *lda,
  CONST Complex af[] /* 2 dimension */,
  CONST int *ldaf,
  CONST int ipiv[],
  CONST Complex b[] /* 2 dimension */,
  CONST int *ldb,
  Complex x[] /* 2 dimension */,
  CONST int *ldx,
  double ferr[],
  double berr[],
  Complex work[],
  double rwork[],
  int *info
);

extern void __stdcall zhetri(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  Complex a[] /* 2 dimension */,
  CONST int *lda,
  CONST int ipiv[],
  Complex work[],
  int *info
);

extern void __stdcall zsytrf(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  Complex a[] /* 2 dimension */,
  CONST int *lda,
  int ipiv[],
  Complex work[],
  CONST int *lwork,
  int *info
);

extern void __stdcall zsytrs(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST int *nrhs,
  CONST Complex a[] /* 2 dimension */,
  CONST int *lda,
  CONST int ipiv[],
  Complex b[] /* 2 dimension */,
  CONST int *ldb,
  int *info
);

extern void __stdcall zsycon(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST Complex a[] /* 2 dimension */,
  CONST int *lda,
  CONST int ipiv[],
  CONST double *anorm,
  double *rcond,
  Complex work[],
  int *info
);

extern void __stdcall F07NVF(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST int *nrhs,
  CONST Complex a[] /* 2 dimension */,
  CONST int *lda,
  CONST Complex af[] /* 2 dimension */,
  CONST int *ldaf,
  CONST int ipiv[],
  CONST Complex b[] /* 2 dimension */,
  CONST int *ldb,
  Complex x[] /* 2 dimension */,
  CONST int *ldx,
  double ferr[],
  double berr[],
  Complex work[],
  double rwork[],
  int *info
);

extern void __stdcall F07NWF(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  Complex a[] /* 2 dimension */,
  CONST int *lda,
  CONST int ipiv[],
  Complex work[],
  int *info
);

extern void __stdcall F07PDF(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  double ap[],
  int ipiv[],
  int *info
);

extern void __stdcall F07PEF(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST int *nrhs,
  CONST double ap[],
  CONST int ipiv[],
  double b[] /* 2 dimension */,
  CONST int *ldb,
  int *info
);

extern void __stdcall F07PGF(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST double ap[],
  CONST int ipiv[],
  CONST double *anorm,
  double *rcond,
  double work[],
  int iwork[],
  int *info
);

extern void __stdcall F07PHF(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST int *nrhs,
  CONST double ap[],
  CONST double afp[],
  CONST int ipiv[],
  CONST double b[] /* 2 dimension */,
  CONST int *ldb,
  double x[] /* 2 dimension */,
  CONST int *ldx,
  double ferr[],
  double berr[],
  double work[],
  int iwork[],
  int *info
);

extern void __stdcall F07PJF(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  double ap[],
  CONST int ipiv[],
  double work[],
  int *info
);

extern void __stdcall F07PRF(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  Complex ap[],
  int ipiv[],
  int *info
);

extern void __stdcall F07PSF(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST int *nrhs,
  CONST Complex ap[],
  CONST int ipiv[],
  Complex b[] /* 2 dimension */,
  CONST int *ldb,
  int *info
);

extern void __stdcall F07PUF(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST Complex ap[],
  CONST int ipiv[],
  CONST double *anorm,
  double *rcond,
  Complex work[],
  int *info
);

extern void __stdcall F07PVF(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST int *nrhs,
  CONST Complex ap[],
  CONST Complex afp[],
  CONST int ipiv[],
  CONST Complex b[] /* 2 dimension */,
  CONST int *ldb,
  Complex x[] /* 2 dimension */,
  CONST int *ldx,
  double ferr[],
  double berr[],
  Complex work[],
  double rwork[],
  int *info
);

extern void __stdcall F07PWF(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  Complex ap[],
  CONST int ipiv[],
  Complex work[],
  int *info
);

extern void __stdcall F07QRF(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  Complex ap[],
  int ipiv[],
  int *info
);

extern void __stdcall F07QSF(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST int *nrhs,
  CONST Complex ap[],
  CONST int ipiv[],
  Complex b[] /* 2 dimension */,
  CONST int *ldb,
  int *info
);

extern void __stdcall F07QUF(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST Complex ap[],
  CONST int ipiv[],
  CONST double *anorm,
  double *rcond,
  Complex work[],
  int *info
);

extern void __stdcall F07QVF(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST int *nrhs,
  CONST Complex ap[],
  CONST Complex afp[],
  CONST int ipiv[],
  CONST Complex b[] /* 2 dimension */,
  CONST int *ldb,
  Complex x[] /* 2 dimension */,
  CONST int *ldx,
  double ferr[],
  double berr[],
  Complex work[],
  double rwork[],
  int *info
);

extern void __stdcall F07QWF(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  Complex ap[],
  CONST int ipiv[],
  Complex work[],
  int *info
);

extern void __stdcall F07TEF(
  CONST char *uplo, CONST int uplo_len,
  CONST char *trans, CONST int trans_len,
  CONST char *diag, CONST int diag_len,
  CONST int *n,
  CONST int *nrhs,
  CONST double a[] /* 2 dimension */,
  CONST int *lda,
  double b[] /* 2 dimension */,
  CONST int *ldb,
  int *info
);

extern void __stdcall F07TGF(
  CONST char *norm, CONST int norm_len,
  CONST char *uplo, CONST int uplo_len,
  CONST char *diag, CONST int diag_len,
  CONST int *n,
  CONST double a[] /* 2 dimension */,
  CONST int *lda,
  double *rcond,
  double work[],
  int iwork[],
  int *info
);

extern void __stdcall F07THF(
  CONST char *uplo, CONST int uplo_len,
  CONST char *trans, CONST int trans_len,
  CONST char *diag, CONST int diag_len,
  CONST int *n,
  CONST int *nrhs,
  CONST double a[] /* 2 dimension */,
  CONST int *lda,
  CONST double b[] /* 2 dimension */,
  CONST int *ldb,
  CONST double x[] /* 2 dimension */,
  CONST int *ldx,
  double ferr[],
  double berr[],
  double work[],
  int iwork[],
  int *info
);

extern void __stdcall F07TJF(
  CONST char *uplo, CONST int uplo_len,
  CONST char *diag, CONST int diag_len,
  CONST int *n,
  double a[] /* 2 dimension */,
  CONST int *lda,
  int *info
);

extern void __stdcall F07TSF(
  CONST char *uplo, CONST int uplo_len,
  CONST char *trans, CONST int trans_len,
  CONST char *diag, CONST int diag_len,
  CONST int *n,
  CONST int *nrhs,
  CONST Complex a[] /* 2 dimension */,
  CONST int *lda,
  Complex b[] /* 2 dimension */,
  CONST int *ldb,
  int *info
);

extern void __stdcall F07TUF(
  CONST char *norm, CONST int norm_len,
  CONST char *uplo, CONST int uplo_len,
  CONST char *diag, CONST int diag_len,
  CONST int *n,
  CONST Complex a[] /* 2 dimension */,
  CONST int *lda,
  double *rcond,
  Complex work[],
  double rwork[],
  int *info
);

extern void __stdcall F07TVF(
  CONST char *uplo, CONST int uplo_len,
  CONST char *trans, CONST int trans_len,
  CONST char *diag, CONST int diag_len,
  CONST int *n,
  CONST int *nrhs,
  CONST Complex a[] /* 2 dimension */,
  CONST int *lda,
  CONST Complex b[] /* 2 dimension */,
  CONST int *ldb,
  CONST Complex x[] /* 2 dimension */,
  CONST int *ldx,
  double ferr[],
  double berr[],
  Complex work[],
  double rwork[],
  int *info
);

extern void __stdcall zsyrfs(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST int *nrhs,
  CONST Complex a[] /* 2 dimension */,
  CONST int *lda,
  CONST Complex af[] /* 2 dimension */,
  CONST int *ldaf,
  CONST int ipiv[],
  CONST Complex b[] /* 2 dimension */,
  CONST int *ldb,
  Complex x[] /* 2 dimension */,
  CONST int *ldx,
  double ferr[],
  double berr[],
  Complex work[],
  double rwork[],
  int *info
);

extern void __stdcall zsytri(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  Complex a[] /* 2 dimension */,
  CONST int *lda,
  CONST int ipiv[],
  Complex work[],
  int *info
);

extern void __stdcall dsptrf(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  double ap[],
  int ipiv[],
  int *info
);

extern void __stdcall dsptrs(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST int *nrhs,
  CONST double ap[],
  CONST int ipiv[],
  double b[] /* 2 dimension */,
  CONST int *ldb,
  int *info
);

extern void __stdcall dspcon(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST double ap[],
  CONST int ipiv[],
  CONST double *anorm,
  double *rcond,
  double work[],
  int iwork[],
  int *info
);

extern void __stdcall dsprfs(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST int *nrhs,
  CONST double ap[],
  CONST double afp[],
  CONST int ipiv[],
  CONST double b[] /* 2 dimension */,
  CONST int *ldb,
  double x[] /* 2 dimension */,
  CONST int *ldx,
  double ferr[],
  double berr[],
  double work[],
  int iwork[],
  int *info
);

extern void __stdcall dsptri(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  double ap[],
  CONST int ipiv[],
  double work[],
  int *info
);

extern void __stdcall zhptrf(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  Complex ap[],
  int ipiv[],
  int *info
);

extern void __stdcall zhptrs(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST int *nrhs,
  CONST Complex ap[],
  CONST int ipiv[],
  Complex b[] /* 2 dimension */,
  CONST int *ldb,
  int *info
);

extern void __stdcall zhpcon(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST Complex ap[],
  CONST int ipiv[],
  CONST double *anorm,
  double *rcond,
  Complex work[],
  int *info
);

extern void __stdcall zhprfs(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST int *nrhs,
  CONST Complex ap[],
  CONST Complex afp[],
  CONST int ipiv[],
  CONST Complex b[] /* 2 dimension */,
  CONST int *ldb,
  Complex x[] /* 2 dimension */,
  CONST int *ldx,
  double ferr[],
  double berr[],
  Complex work[],
  double rwork[],
  int *info
);

extern void __stdcall zhptri(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  Complex ap[],
  CONST int ipiv[],
  Complex work[],
  int *info
);

extern void __stdcall zsptrf(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  Complex ap[],
  int ipiv[],
  int *info
);

extern void __stdcall zsptrs(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST int *nrhs,
  CONST Complex ap[],
  CONST int ipiv[],
  Complex b[] /* 2 dimension */,
  CONST int *ldb,
  int *info
);

extern void __stdcall zspcon(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST Complex ap[],
  CONST int ipiv[],
  CONST double *anorm,
  double *rcond,
  Complex work[],
  int *info
);

extern void __stdcall zsprfs(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST int *nrhs,
  CONST Complex ap[],
  CONST Complex afp[],
  CONST int ipiv[],
  CONST Complex b[] /* 2 dimension */,
  CONST int *ldb,
  Complex x[] /* 2 dimension */,
  CONST int *ldx,
  double ferr[],
  double berr[],
  Complex work[],
  double rwork[],
  int *info
);

extern void __stdcall zsptri(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  Complex ap[],
  CONST int ipiv[],
  Complex work[],
  int *info
);

extern void __stdcall dtrtrs(
  CONST char *uplo, CONST int uplo_len,
  CONST char *trans, CONST int trans_len,
  CONST char *diag, CONST int diag_len,
  CONST int *n,
  CONST int *nrhs,
  CONST double a[] /* 2 dimension */,
  CONST int *lda,
  double b[] /* 2 dimension */,
  CONST int *ldb,
  int *info
);

extern void __stdcall dtrcon(
  CONST char *norm, CONST int norm_len,
  CONST char *uplo, CONST int uplo_len,
  CONST char *diag, CONST int diag_len,
  CONST int *n,
  CONST double a[] /* 2 dimension */,
  CONST int *lda,
  double *rcond,
  double work[],
  int iwork[],
  int *info
);

extern void __stdcall dtrrfs(
  CONST char *uplo, CONST int uplo_len,
  CONST char *trans, CONST int trans_len,
  CONST char *diag, CONST int diag_len,
  CONST int *n,
  CONST int *nrhs,
  CONST double a[] /* 2 dimension */,
  CONST int *lda,
  CONST double b[] /* 2 dimension */,
  CONST int *ldb,
  CONST double x[] /* 2 dimension */,
  CONST int *ldx,
  double ferr[],
  double berr[],
  double work[],
  int iwork[],
  int *info
);

extern void __stdcall dtrtri(
  CONST char *uplo, CONST int uplo_len,
  CONST char *diag, CONST int diag_len,
  CONST int *n,
  double a[] /* 2 dimension */,
  CONST int *lda,
  int *info
);

extern void __stdcall ztrtrs(
  CONST char *uplo, CONST int uplo_len,
  CONST char *trans, CONST int trans_len,
  CONST char *diag, CONST int diag_len,
  CONST int *n,
  CONST int *nrhs,
  CONST Complex a[] /* 2 dimension */,
  CONST int *lda,
  Complex b[] /* 2 dimension */,
  CONST int *ldb,
  int *info
);

extern void __stdcall ztrcon(
  CONST char *norm, CONST int norm_len,
  CONST char *uplo, CONST int uplo_len,
  CONST char *diag, CONST int diag_len,
  CONST int *n,
  CONST Complex a[] /* 2 dimension */,
  CONST int *lda,
  double *rcond,
  Complex work[],
  double rwork[],
  int *info
);

extern void __stdcall ztrrfs(
  CONST char *uplo, CONST int uplo_len,
  CONST char *trans, CONST int trans_len,
  CONST char *diag, CONST int diag_len,
  CONST int *n,
  CONST int *nrhs,
  CONST Complex a[] /* 2 dimension */,
  CONST int *lda,
  CONST Complex b[] /* 2 dimension */,
  CONST int *ldb,
  CONST Complex x[] /* 2 dimension */,
  CONST int *ldx,
  double ferr[],
  double berr[],
  Complex work[],
  double rwork[],
  int *info
);

extern void __stdcall F07TWF(
  CONST char *uplo, CONST int uplo_len,
  CONST char *diag, CONST int diag_len,
  CONST int *n,
  Complex a[] /* 2 dimension */,
  CONST int *lda,
  int *info
);

extern void __stdcall F07UEF(
  CONST char *uplo, CONST int uplo_len,
  CONST char *trans, CONST int trans_len,
  CONST char *diag, CONST int diag_len,
  CONST int *n,
  CONST int *nrhs,
  CONST double ap[],
  double b[] /* 2 dimension */,
  CONST int *ldb,
  int *info
);

extern void __stdcall F07UGF(
  CONST char *norm, CONST int norm_len,
  CONST char *uplo, CONST int uplo_len,
  CONST char *diag, CONST int diag_len,
  CONST int *n,
  CONST double ap[],
  double *rcond,
  double work[],
  int iwork[],
  int *info
);

extern void __stdcall F07UHF(
  CONST char *uplo, CONST int uplo_len,
  CONST char *trans, CONST int trans_len,
  CONST char *diag, CONST int diag_len,
  CONST int *n,
  CONST int *nrhs,
  CONST double ap[],
  CONST double b[] /* 2 dimension */,
  CONST int *ldb,
  CONST double x[] /* 2 dimension */,
  CONST int *ldx,
  double ferr[],
  double berr[],
  double work[],
  int iwork[],
  int *info
);

extern void __stdcall F07UJF(
  CONST char *uplo, CONST int uplo_len,
  CONST char *diag, CONST int diag_len,
  CONST int *n,
  double ap[],
  int *info
);

extern void __stdcall F07USF(
  CONST char *uplo, CONST int uplo_len,
  CONST char *trans, CONST int trans_len,
  CONST char *diag, CONST int diag_len,
  CONST int *n,
  CONST int *nrhs,
  CONST Complex ap[],
  Complex b[] /* 2 dimension */,
  CONST int *ldb,
  int *info
);

extern void __stdcall F07UUF(
  CONST char *norm, CONST int norm_len,
  CONST char *uplo, CONST int uplo_len,
  CONST char *diag, CONST int diag_len,
  CONST int *n,
  CONST Complex ap[],
  double *rcond,
  Complex work[],
  double rwork[],
  int *info
);

extern void __stdcall F07UVF(
  CONST char *uplo, CONST int uplo_len,
  CONST char *trans, CONST int trans_len,
  CONST char *diag, CONST int diag_len,
  CONST int *n,
  CONST int *nrhs,
  CONST Complex ap[],
  CONST Complex b[] /* 2 dimension */,
  CONST int *ldb,
  CONST Complex x[] /* 2 dimension */,
  CONST int *ldx,
  double ferr[],
  double berr[],
  Complex work[],
  double rwork[],
  int *info
);

extern void __stdcall F07UWF(
  CONST char *uplo, CONST int uplo_len,
  CONST char *diag, CONST int diag_len,
  CONST int *n,
  Complex ap[],
  int *info
);

extern void __stdcall F07VEF(
  CONST char *uplo, CONST int uplo_len,
  CONST char *trans, CONST int trans_len,
  CONST char *diag, CONST int diag_len,
  CONST int *n,
  CONST int *kd,
  CONST int *nrhs,
  CONST double ab[] /* 2 dimension */,
  CONST int *ldab,
  double b[] /* 2 dimension */,
  CONST int *ldb,
  int *info
);

extern void __stdcall F07VGF(
  CONST char *norm, CONST int norm_len,
  CONST char *uplo, CONST int uplo_len,
  CONST char *diag, CONST int diag_len,
  CONST int *n,
  CONST int *kd,
  CONST double ab[] /* 2 dimension */,
  CONST int *ldab,
  double *rcond,
  double work[],
  int iwork[],
  int *info
);

extern void __stdcall F07VHF(
  CONST char *uplo, CONST int uplo_len,
  CONST char *trans, CONST int trans_len,
  CONST char *diag, CONST int diag_len,
  CONST int *n,
  CONST int *kd,
  CONST int *nrhs,
  CONST double ab[] /* 2 dimension */,
  CONST int *ldab,
  CONST double b[] /* 2 dimension */,
  CONST int *ldb,
  CONST double x[] /* 2 dimension */,
  CONST int *ldx,
  double ferr[],
  double berr[],
  double work[],
  int iwork[],
  int *info
);

extern void __stdcall F07VSF(
  CONST char *uplo, CONST int uplo_len,
  CONST char *trans, CONST int trans_len,
  CONST char *diag, CONST int diag_len,
  CONST int *n,
  CONST int *kd,
  CONST int *nrhs,
  CONST Complex ab[] /* 2 dimension */,
  CONST int *ldab,
  Complex b[] /* 2 dimension */,
  CONST int *ldb,
  int *info
);

extern void __stdcall F07VUF(
  CONST char *norm, CONST int norm_len,
  CONST char *uplo, CONST int uplo_len,
  CONST char *diag, CONST int diag_len,
  CONST int *n,
  CONST int *kd,
  CONST Complex ab[] /* 2 dimension */,
  CONST int *ldab,
  double *rcond,
  Complex work[],
  double rwork[],
  int *info
);

extern void __stdcall F07VVF(
  CONST char *uplo, CONST int uplo_len,
  CONST char *trans, CONST int trans_len,
  CONST char *diag, CONST int diag_len,
  CONST int *n,
  CONST int *kd,
  CONST int *nrhs,
  CONST Complex ab[] /* 2 dimension */,
  CONST int *ldab,
  CONST Complex b[] /* 2 dimension */,
  CONST int *ldb,
  CONST Complex x[] /* 2 dimension */,
  CONST int *ldx,
  double ferr[],
  double berr[],
  Complex work[],
  double rwork[],
  int *info
);

extern void __stdcall F08AEF(
  CONST int *m,
  CONST int *n,
  double a[] /* 2 dimension */,
  CONST int *lda,
  double tau[],
  double work[],
  CONST int *lwork,
  int *info
);

extern void __stdcall F08AFF(
  CONST int *m,
  CONST int *n,
  CONST int *k,
  double a[] /* 2 dimension */,
  CONST int *lda,
  CONST double tau[],
  double work[],
  CONST int *lwork,
  int *info
);

extern void __stdcall F08AGF(
  CONST char *side, CONST int side_len,
  CONST char *trans, CONST int trans_len,
  CONST int *m,
  CONST int *n,
  CONST int *k,
  double a[] /* 2 dimension */,
  CONST int *lda,
  CONST double tau[],
  double c[] /* 2 dimension */,
  CONST int *ldc,
  double work[],
  CONST int *lwork,
  int *info
);

extern void __stdcall F08AHF(
  CONST int *m,
  CONST int *n,
  double a[] /* 2 dimension */,
  CONST int *lda,
  double tau[],
  double work[],
  CONST int *lwork,
  int *info
);

extern void __stdcall F08AJF(
  CONST int *m,
  CONST int *n,
  CONST int *k,
  double a[] /* 2 dimension */,
  CONST int *lda,
  CONST double tau[],
  double work[],
  CONST int *lwork,
  int *info
);

extern void __stdcall F08AKF(
  CONST char *side, CONST int side_len,
  CONST char *trans, CONST int trans_len,
  CONST int *m,
  CONST int *n,
  CONST int *k,
  double a[] /* 2 dimension */,
  CONST int *lda,
  CONST double tau[],
  double c[] /* 2 dimension */,
  CONST int *ldc,
  double work[],
  CONST int *lwork,
  int *info
);

extern void __stdcall F08ASF(
  CONST int *m,
  CONST int *n,
  Complex a[] /* 2 dimension */,
  CONST int *lda,
  Complex tau[],
  Complex work[],
  CONST int *lwork,
  int *info
);

extern void __stdcall F08ATF(
  CONST int *m,
  CONST int *n,
  CONST int *k,
  Complex a[] /* 2 dimension */,
  CONST int *lda,
  CONST Complex tau[],
  Complex work[],
  CONST int *lwork,
  int *info
);

extern void __stdcall F08AUF(
  CONST char *side, CONST int side_len,
  CONST char *trans, CONST int trans_len,
  CONST int *m,
  CONST int *n,
  CONST int *k,
  Complex a[] /* 2 dimension */,
  CONST int *lda,
  CONST Complex tau[],
  Complex c[] /* 2 dimension */,
  CONST int *ldc,
  Complex work[],
  CONST int *lwork,
  int *info
);

extern void __stdcall F08AVF(
  CONST int *m,
  CONST int *n,
  Complex a[] /* 2 dimension */,
  CONST int *lda,
  Complex tau[],
  Complex work[],
  CONST int *lwork,
  int *info
);

extern void __stdcall F08AWF(
  CONST int *m,
  CONST int *n,
  CONST int *k,
  Complex a[] /* 2 dimension */,
  CONST int *lda,
  CONST Complex tau[],
  Complex work[],
  CONST int *lwork,
  int *info
);

extern void __stdcall F08AXF(
  CONST char *side, CONST int side_len,
  CONST char *trans, CONST int trans_len,
  CONST int *m,
  CONST int *n,
  CONST int *k,
  Complex a[] /* 2 dimension */,
  CONST int *lda,
  CONST Complex tau[],
  Complex c[] /* 2 dimension */,
  CONST int *ldc,
  Complex work[],
  CONST int *lwork,
  int *info
);

extern void __stdcall F08BEF(
  CONST int *m,
  CONST int *n,
  double a[] /* 2 dimension */,
  CONST int *lda,
  int jpvt[],
  double tau[],
  double work[],
  int *info
);

extern void __stdcall F08BSF(
  CONST int *m,
  CONST int *n,
  Complex a[] /* 2 dimension */,
  CONST int *lda,
  int jpvt[],
  Complex tau[],
  Complex work[],
  double rwork[],
  int *info
);

extern void __stdcall F08FEF(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  double a[] /* 2 dimension */,
  CONST int *lda,
  double d[],
  double e[],
  double tau[],
  double work[],
  CONST int *lwork,
  int *info
);

extern void __stdcall F08FFF(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  double a[] /* 2 dimension */,
  CONST int *lda,
  CONST double tau[],
  double work[],
  CONST int *lwork,
  int *info
);

extern void __stdcall F08FGF(
  CONST char *side, CONST int side_len,
  CONST char *uplo, CONST int uplo_len,
  CONST char *trans, CONST int trans_len,
  CONST int *m,
  CONST int *n,
  double a[] /* 2 dimension */,
  CONST int *lda,
  CONST double tau[],
  double c[] /* 2 dimension */,
  CONST int *ldc,
  double work[],
  CONST int *lwork,
  int *info
);

extern void __stdcall F08FSF(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  Complex a[] /* 2 dimension */,
  CONST int *lda,
  double d[],
  double e[],
  Complex tau[],
  Complex work[],
  CONST int *lwork,
  int *info
);

extern void __stdcall F08FTF(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  Complex a[] /* 2 dimension */,
  CONST int *lda,
  CONST Complex tau[],
  Complex work[],
  CONST int *lwork,
  int *info
);

extern void __stdcall F08FUF(
  CONST char *side, CONST int side_len,
  CONST char *uplo, CONST int uplo_len,
  CONST char *trans, CONST int trans_len,
  CONST int *m,
  CONST int *n,
  Complex a[] /* 2 dimension */,
  CONST int *lda,
  CONST Complex tau[],
  Complex c[] /* 2 dimension */,
  CONST int *ldc,
  Complex work[],
  CONST int *lwork,
  int *info
);

extern void __stdcall F08GEF(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  double ap[],
  double d[],
  double e[],
  double tau[],
  int *info
);

extern void __stdcall F08GFF(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST double ap[],
  CONST double tau[],
  double q[] /* 2 dimension */,
  CONST int *ldq,
  double work[],
  int *info
);

extern void __stdcall F08GGF(
  CONST char *side, CONST int side_len,
  CONST char *uplo, CONST int uplo_len,
  CONST char *trans, CONST int trans_len,
  CONST int *m,
  CONST int *n,
  double ap[],
  CONST double tau[],
  double c[] /* 2 dimension */,
  CONST int *ldc,
  double work[],
  int *info
);

extern void __stdcall F08GSF(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  Complex ap[],
  double d[],
  double e[],
  Complex tau[],
  int *info
);

extern void __stdcall F08GTF(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST Complex ap[],
  CONST Complex tau[],
  Complex q[] /* 2 dimension */,
  CONST int *ldq,
  Complex work[],
  int *info
);

extern void __stdcall F08GUF(
  CONST char *side, CONST int side_len,
  CONST char *uplo, CONST int uplo_len,
  CONST char *trans, CONST int trans_len,
  CONST int *m,
  CONST int *n,
  Complex ap[],
  CONST Complex tau[],
  Complex c[] /* 2 dimension */,
  CONST int *ldc,
  Complex work[],
  int *info
);

extern void __stdcall F08HEF(
  CONST char *vect, CONST int vect_len,
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST int *kd,
  double ab[] /* 2 dimension */,
  CONST int *ldab,
  double d[],
  double e[],
  double q[] /* 2 dimension */,
  CONST int *ldq,
  double work[],
  int *info
);

extern void __stdcall F08HSF(
  CONST char *vect, CONST int vect_len,
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST int *kd,
  Complex ab[] /* 2 dimension */,
  CONST int *ldab,
  double d[],
  double e[],
  Complex q[] /* 2 dimension */,
  CONST int *ldq,
  Complex work[],
  int *info
);

extern void __stdcall F08JEF(
  CONST char *compz, CONST int compz_len,
  CONST int *n,
  double d[],
  double e[],
  double z[] /* 2 dimension */,
  CONST int *ldz,
  double work[],
  int *info
);

extern void __stdcall F08JFF(
  CONST int *n,
  double d[],
  double e[],
  int *info
);

extern void __stdcall F08JGF(
  CONST char *compz, CONST int compz_len,
  CONST int *n,
  double d[],
  double e[],
  double z[] /* 2 dimension */,
  CONST int *ldz,
  double work[],
  int *info
);

extern void __stdcall F08JJF(
  CONST char *range, CONST int range_len,
  CONST char *order, CONST int order_len,
  CONST int *n,
  CONST double *vl,
  CONST double *vu,
  CONST int *il,
  CONST int *iu,
  CONST double *abstol,
  CONST double d[],
  CONST double e[],
  int *m,
  int *nsplit,
  double w[],
  int iblock[],
  int isplit[],
  double work[],
  int iwork[],
  int *info
);

extern void __stdcall F08JKF(
  CONST int *n,
  CONST double d[],
  CONST double e[],
  CONST int *m,
  CONST double w[],
  CONST int iblock[],
  CONST int isplit[],
  double z[] /* 2 dimension */,
  CONST int *ldz,
  double work[],
  int iwork[],
  int ifail[],
  int *info
);

extern void __stdcall F08JSF(
  CONST char *compz, CONST int compz_len,
  CONST int *n,
  double d[],
  double e[],
  Complex z[] /* 2 dimension */,
  CONST int *ldz,
  double work[],
  int *info
);

extern void __stdcall F08JUF(
  CONST char *compz, CONST int compz_len,
  CONST int *n,
  double d[],
  double e[],
  Complex z[] /* 2 dimension */,
  CONST int *ldz,
  double work[],
  int *info
);

extern void __stdcall F08JXF(
  CONST int *n,
  CONST double d[],
  CONST double e[],
  CONST int *m,
  CONST double w[],
  CONST int iblock[],
  CONST int isplit[],
  Complex z[] /* 2 dimension */,
  CONST int *ldz,
  double work[],
  int iwork[],
  int ifail[],
  int *info
);

extern void __stdcall F08KEF(
  CONST int *m,
  CONST int *n,
  double a[] /* 2 dimension */,
  CONST int *lda,
  double d[],
  double e[],
  double tauq[],
  double taup[],
  double work[],
  CONST int *lwork,
  int *info
);

extern void __stdcall F08KFF(
  CONST char *vect, CONST int vect_len,
  CONST int *m,
  CONST int *n,
  CONST int *k,
  double a[] /* 2 dimension */,
  CONST int *lda,
  CONST double tau[],
  double work[],
  CONST int *lwork,
  int *info
);

extern void __stdcall F08KGF(
  CONST char *vect, CONST int vect_len,
  CONST char *side, CONST int side_len,
  CONST char *trans, CONST int trans_len,
  CONST int *m,
  CONST int *n,
  CONST int *k,
  CONST double a[] /* 2 dimension */,
  CONST int *lda,
  CONST double tau[],
  double c[] /* 2 dimension */,
  CONST int *ldc,
  double work[],
  CONST int *lwork,
  int *info
);

extern void __stdcall F08KSF(
  CONST int *m,
  CONST int *n,
  Complex a[] /* 2 dimension */,
  CONST int *lda,
  double d[],
  double e[],
  Complex tauq[],
  Complex taup[],
  Complex work[],
  CONST int *lwork,
  int *info
);

extern void __stdcall F08KTF(
  CONST char *vect, CONST int vect_len,
  CONST int *m,
  CONST int *n,
  CONST int *k,
  Complex a[] /* 2 dimension */,
  CONST int *lda,
  CONST Complex tau[],
  Complex work[],
  CONST int *lwork,
  int *info
);

extern void __stdcall F08KUF(
  CONST char *vect, CONST int vect_len,
  CONST char *side, CONST int side_len,
  CONST char *trans, CONST int trans_len,
  CONST int *m,
  CONST int *n,
  CONST int *k,
  CONST Complex a[] /* 2 dimension */,
  CONST int *lda,
  CONST Complex tau[],
  Complex c[] /* 2 dimension */,
  CONST int *ldc,
  Complex work[],
  CONST int *lwork,
  int *info
);

extern void __stdcall F08MEF(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST int *ncvt,
  CONST int *nru,
  CONST int *ncc,
  double d[],
  double e[],
  double vt[] /* 2 dimension */,
  CONST int *ldvt,
  double u[] /* 2 dimension */,
  CONST int *ldu,
  double c[] /* 2 dimension */,
  CONST int *ldc,
  double work[],
  int *info
);

extern void __stdcall F08MSF(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST int *ncvt,
  CONST int *nru,
  CONST int *ncc,
  double d[],
  double e[],
  Complex vt[] /* 2 dimension */,
  CONST int *ldvt,
  Complex u[] /* 2 dimension */,
  CONST int *ldu,
  Complex c[] /* 2 dimension */,
  CONST int *ldc,
  double work[],
  int *info
);

extern void __stdcall F08NEF(
  CONST int *n,
  CONST int *ilo,
  CONST int *ihi,
  double a[] /* 2 dimension */,
  CONST int *lda,
  double tau[],
  double work[],
  CONST int *lwork,
  int *info
);

extern void __stdcall F08NFF(
  CONST int *n,
  CONST int *ilo,
  CONST int *ihi,
  double a[] /* 2 dimension */,
  CONST int *lda,
  CONST double tau[],
  double work[],
  CONST int *lwork,
  int *info
);

extern void __stdcall F08NGF(
  CONST char *side, CONST int side_len,
  CONST char *trans, CONST int trans_len,
  CONST int *m,
  CONST int *n,
  CONST int *ilo,
  CONST int *ihi,
  CONST double a[] /* 2 dimension */,
  CONST int *lda,
  CONST double tau[],
  double c[] /* 2 dimension */,
  CONST int *ldc,
  double work[],
  CONST int *lwork,
  int *info
);

extern void __stdcall F08NHF(
  CONST char *job, CONST int job_len,
  CONST int *n,
  double a[] /* 2 dimension */,
  CONST int *lda,
  int *ilo,
  int *ihi,
  double scale[],
  int *info
);

extern void __stdcall F08NJF(
  CONST char *job, CONST int job_len,
  CONST char *side, CONST int side_len,
  CONST int *n,
  CONST int *ilo,
  CONST int *ihi,
  CONST double scale[],
  CONST int *m,
  double v[] /* 2 dimension */,
  CONST int *ldv,
  int *info
);

extern void __stdcall F08NSF(
  CONST int *n,
  CONST int *ilo,
  CONST int *ihi,
  Complex a[] /* 2 dimension */,
  CONST int *lda,
  Complex tau[],
  Complex work[],
  CONST int *lwork,
  int *info
);

extern void __stdcall F08NTF(
  CONST int *n,
  CONST int *ilo,
  CONST int *ihi,
  Complex a[] /* 2 dimension */,
  CONST int *lda,
  CONST Complex tau[],
  Complex work[],
  CONST int *lwork,
  int *info
);

extern void __stdcall F08NUF(
  CONST char *side, CONST int side_len,
  CONST char *trans, CONST int trans_len,
  CONST int *m,
  CONST int *n,
  CONST int *ilo,
  CONST int *ihi,
  CONST Complex a[] /* 2 dimension */,
  CONST int *lda,
  CONST Complex tau[],
  Complex c[] /* 2 dimension */,
  CONST int *ldc,
  Complex work[],
  CONST int *lwork,
  int *info
);

extern void __stdcall F08NVF(
  CONST char *job, CONST int job_len,
  CONST int *n,
  Complex a[] /* 2 dimension */,
  CONST int *lda,
  int *ilo,
  int *ihi,
  double scale[],
  int *info
);

extern void __stdcall F08NWF(
  CONST char *job, CONST int job_len,
  CONST char *side, CONST int side_len,
  CONST int *n,
  CONST int *ilo,
  CONST int *ihi,
  CONST double scale[],
  CONST int *m,
  Complex v[] /* 2 dimension */,
  CONST int *ldv,
  int *info
);

extern void __stdcall F08PEF(
  CONST char *job, CONST int job_len,
  CONST char *compz, CONST int compz_len,
  CONST int *n,
  CONST int *ilo,
  CONST int *ihi,
  double h[] /* 2 dimension */,
  CONST int *ldh,
  double wr[],
  double wi[],
  double z[] /* 2 dimension */,
  CONST int *ldz,
  double work[],
  CONST int *lwork,
  int *info
);

extern void __stdcall F08PKF(
  CONST char *job, CONST int job_len,
  CONST char *eigsrc, CONST int eigsrc_len,
  CONST char *initv, CONST int initv_len,
  int select[],
  CONST int *n,
  CONST double h[] /* 2 dimension */,
  CONST int *ldh,
  double wr[],
  CONST double wi[],
  double vl[] /* 2 dimension */,
  CONST int *ldvl,
  double vr[] /* 2 dimension */,
  CONST int *ldvr,
  CONST int *mm,
  int *m,
  double work[],
  int ifaill[],
  int ifailr[],
  int *info
);

extern void __stdcall F08PSF(
  CONST char *job, CONST int job_len,
  CONST char *compz, CONST int compz_len,
  CONST int *n,
  CONST int *ilo,
  CONST int *ihi,
  Complex h[] /* 2 dimension */,
  CONST int *ldh,
  Complex w[],
  Complex z[] /* 2 dimension */,
  CONST int *ldz,
  Complex work[],
  CONST int *lwork,
  int *info
);

extern void __stdcall F08PXF(
  CONST char *job, CONST int job_len,
  CONST char *eigsrc, CONST int eigsrc_len,
  CONST char *initv, CONST int initv_len,
  CONST int select[],
  CONST int *n,
  CONST Complex h[] /* 2 dimension */,
  CONST int *ldh,
  Complex w[],
  Complex vl[] /* 2 dimension */,
  CONST int *ldvl,
  Complex vr[] /* 2 dimension */,
  CONST int *ldvr,
  CONST int *mm,
  int *m,
  Complex work[],
  double rwork[],
  int ifaill[],
  int ifailr[],
  int *info
);

extern void __stdcall F08QFF(
  CONST char *compq, CONST int compq_len,
  CONST int *n,
  double t[] /* 2 dimension */,
  CONST int *ldt,
  double q[] /* 2 dimension */,
  CONST int *ldq,
  int *ifst,
  int *ilst,
  double work[],
  int *info
);

extern void __stdcall F08QGF(
  CONST char *job, CONST int job_len,
  CONST char *compq, CONST int compq_len,
  CONST int select[],
  CONST int *n,
  double t[] /* 2 dimension */,
  CONST int *ldt,
  CONST double q[] /* 2 dimension */,
  CONST int *ldq,
  double wr[],
  double wi[],
  int *m,
  double *s,
  double *sep,
  double work[],
  CONST int *lwork,
  int iwork[],
  CONST int *liwork,
  int *info
);

extern void __stdcall F08QHF(
  CONST char *trana, CONST int trana_len,
  CONST char *tranb, CONST int tranb_len,
  CONST int *isgn,
  CONST int *m,
  CONST int *n,
  CONST double a[] /* 2 dimension */,
  CONST int *lda,
  CONST double b[] /* 2 dimension */,
  CONST int *ldb,
  double c[] /* 2 dimension */,
  CONST int *ldc,
  double *scale,
  int *info
);

extern void __stdcall F08QKF(
  CONST char *job, CONST int job_len,
  CONST char *howmny, CONST int howmny_len,
  int select[],
  CONST int *n,
  CONST double t[] /* 2 dimension */,
  CONST int *ldt,
  double vl[] /* 2 dimension */,
  CONST int *ldvl,
  double vr[] /* 2 dimension */,
  CONST int *ldvr,
  CONST int *mm,
  int *m,
  double work[],
  int *info
);

extern void __stdcall F08QLF(
  CONST char *job, CONST int job_len,
  CONST char *howmny, CONST int howmny_len,
  CONST int select[],
  CONST int *n,
  CONST double t[] /* 2 dimension */,
  CONST int *ldt,
  CONST double vl[] /* 2 dimension */,
  CONST int *ldvl,
  CONST double vr[] /* 2 dimension */,
  CONST int *ldvr,
  double s[],
  double sep[],
  CONST int *mm,
  int *m,
  double work[] /* 2 dimension */,
  CONST int *ldwork,
  int iwork[],
  int *info
);

extern void __stdcall F08QTF(
  CONST char *compq, CONST int compq_len,
  CONST int *n,
  Complex t[] /* 2 dimension */,
  CONST int *ldt,
  Complex q[] /* 2 dimension */,
  CONST int *ldq,
  CONST int *ifst,
  CONST int *ilst,
  int *info
);

extern void __stdcall F08QUF(
  CONST char *job, CONST int job_len,
  CONST char *compq, CONST int compq_len,
  CONST int select[],
  CONST int *n,
  Complex t[] /* 2 dimension */,
  CONST int *ldt,
  Complex q[] /* 2 dimension */,
  CONST int *ldq,
  Complex w[],
  int *m,
  double *s,
  double *sep,
  Complex work[],
  CONST int *lwork,
  int *info
);

extern void __stdcall F08QVF(
  CONST char *trana, CONST int trana_len,
  CONST char *tranb, CONST int tranb_len,
  CONST int *isgn,
  CONST int *m,
  CONST int *n,
  CONST Complex a[] /* 2 dimension */,
  CONST int *lda,
  CONST Complex b[] /* 2 dimension */,
  CONST int *ldb,
  Complex c[] /* 2 dimension */,
  CONST int *ldc,
  double *scale,
  int *info
);

extern void __stdcall F08QXF(
  CONST char *job, CONST int job_len,
  CONST char *howmny, CONST int howmny_len,
  CONST int select[],
  CONST int *n,
  Complex t[] /* 2 dimension */,
  CONST int *ldt,
  Complex vl[] /* 2 dimension */,
  CONST int *ldvl,
  Complex vr[] /* 2 dimension */,
  CONST int *ldvr,
  CONST int *mm,
  int *m,
  Complex work[],
  double rwork[],
  int *info
);

extern void __stdcall F08QYF(
  CONST char *job, CONST int job_len,
  CONST char *howmny, CONST int howmny_len,
  CONST int select[],
  CONST int *n,
  CONST Complex t[] /* 2 dimension */,
  CONST int *ldt,
  CONST Complex vl[] /* 2 dimension */,
  CONST int *ldvl,
  CONST Complex vr[] /* 2 dimension */,
  CONST int *ldvr,
  double s[],
  double sep[],
  CONST int *mm,
  int *m,
  Complex work[] /* 2 dimension */,
  CONST int *ldwork,
  double rwork[],
  int *info
);

extern void __stdcall F08SEF(
  CONST int *itype,
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  double a[] /* 2 dimension */,
  CONST int *lda,
  CONST double b[] /* 2 dimension */,
  CONST int *ldb,
  int *info
);

extern void __stdcall F08SSF(
  CONST int *itype,
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  Complex a[] /* 2 dimension */,
  CONST int *lda,
  Complex b[] /* 2 dimension */,
  CONST int *ldb,
  int *info
);

extern void __stdcall F08TEF(
  CONST int *itype,
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  double ap[],
  CONST double bp[],
  int *info
);

extern void __stdcall F08TSF(
  CONST int *itype,
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  Complex ap[],
  CONST Complex bp[],
  int *info
);

extern void __stdcall F11BAF(
  CONST char *method, CONST int method_len,
  CONST char *precon, CONST int precon_len,
  CONST char *norm, CONST int norm_len,
  CONST char *weight, CONST int weight_len,
  CONST int *iterm,
  CONST int *n,
  CONST int *m,
  CONST double *tol,
  CONST int *maxitn,
  CONST double *anorm,
  CONST double *sigmax,
  CONST int *monit,
  int *lwreq,
  int *ifail
);

extern void __stdcall F11BBF(
  int *irevcm,
  double u[],
  double v[],
  double work[],
  CONST int *lwork,
  int *ifail
);

extern void __stdcall F11BCF(
  int *itn,
  double *stplhs,
  double *stprhs,
  double *anorm,
  double *sigmax,
  int *ifail
);

extern void __stdcall F11DAF(
  CONST int *n,
  CONST int *nnz,
  double a[],
  CONST int *la,
  int irow[],
  int icol[],
  CONST int *lfill,
  CONST double *dtol,
  CONST char *pstrat, CONST int pstrat_len,
  CONST char *milu, CONST int milu_len,
  int ipivp[],
  int ipivq[],
  int istr[],
  int idiag[],
  int *nnzc,
  int *npivm,
  int iwork[],
  CONST int *liwork,
  int *ifail
);

extern void __stdcall F11DBF(
  CONST char *trans, CONST int trans_len,
  CONST int *n,
  CONST double a[],
  CONST int *la,
  CONST int irow[],
  CONST int icol[],
  int ipivp[],
  int ipivq[],
  CONST int istr[],
  CONST int idiag[],
  CONST char *check, CONST int check_len,
  CONST double y[],
  double x[],
  int *ifail
);

extern void __stdcall F11DCF(
  CONST char *method, CONST int method_len,
  CONST int *n,
  CONST int *nnz,
  CONST double a[],
  CONST int *la,
  CONST int irow[],
  CONST int icol[],
  int ipivp[],
  int ipivq[],
  CONST int istr[],
  CONST int idiag[],
  CONST double b[],
  CONST int *m,
  CONST double *tol,
  CONST int *maxitn,
  double x[],
  double *rnorm,
  int *itn,
  double work[],
  CONST int *lwork,
  int *ifail
);

extern void __stdcall F11DDF(
  CONST char *trans, CONST int trans_len,
  CONST int *n,
  CONST int *nnz,
  CONST double a[],
  CONST int irow[],
  CONST int icol[],
  CONST double rdiag[],
  CONST double *omega,
  CONST char *check, CONST int check_len,
  CONST double y[],
  double x[],
  int iwork[],
  int *ifail
);

extern void __stdcall F11DEF(
  CONST char *method, CONST int method_len,
  CONST char *precon, CONST int precon_len,
  CONST int *n,
  CONST int *nnz,
  CONST double a[],
  CONST int irow[],
  CONST int icol[],
  CONST double *omega,
  CONST double b[],
  CONST int *m,
  CONST double *tol,
  CONST int *maxitn,
  double x[],
  double *rnorm,
  int *itn,
  double work[],
  CONST int *lwork,
  int iwork[],
  int *ifail
);

extern void __stdcall F11GAF(
  CONST char *method, CONST int method_len,
  CONST char *precon, CONST int precon_len,
  CONST char *sigcmp, CONST int sigcmp_len,
  CONST char *norm, CONST int norm_len,
  CONST char *weight, CONST int weight_len,
  CONST int *iterm,
  CONST int *n,
  CONST double *tol,
  CONST int *maxitn,
  CONST double *anorm,
  CONST double *sigmax,
  CONST double *sigtol,
  CONST int *maxits,
  CONST int *monit,
  int *lwreq,
  int *ifail
);

extern void __stdcall F11GBF(
  int *irevcm,
  double u[],
  double v[],
  double work[],
  CONST int *lwork,
  int *ifail
);

extern void __stdcall F11GCF(
  int *itn,
  double *stplhs,
  double *stprhs,
  double *anorm,
  double *sigmax,
  int *its,
  double *sigerr,
  int *ifail
);

extern void __stdcall F11JAF(
  CONST int *n,
  CONST int *nnz,
  double a[],
  CONST int *la,
  int irow[],
  int icol[],
  CONST int *lfill,
  CONST double *dtol,
  CONST char *mic, CONST int mic_len,
  CONST double *dscale,
  CONST char *pstrat, CONST int pstrat_len,
  int ipiv[],
  int istr[],
  int *nnzc,
  int *npivm,
  int iwork[],
  CONST int *liwork,
  int *ifail
);

extern void __stdcall F11JBF(
  CONST int *n,
  CONST double a[],
  CONST int *la,
  CONST int irow[],
  CONST int icol[],
  int ipiv[],
  int istr[],
  CONST char *check, CONST int check_len,
  CONST double y[],
  double x[],
  int *ifail
);

extern void __stdcall F11JCF(
  CONST char *method, CONST int method_len,
  CONST int *n,
  CONST int *nnz,
  CONST double a[],
  CONST int *la,
  CONST int irow[],
  CONST int icol[],
  int ipiv[],
  int istr[],
  CONST double b[],
  CONST double *tol,
  CONST int *maxitn,
  double x[],
  double *rnorm,
  int *itn,
  double work[],
  CONST int *lwork,
  int *ifail
);

extern void __stdcall F11JDF(
  CONST int *n,
  CONST int *nnz,
  CONST double a[],
  CONST int irow[],
  CONST int icol[],
  CONST double rdiag[],
  CONST double *omega,
  CONST char *check, CONST int check_len,
  CONST double y[],
  double x[],
  int iwork[],
  int *ifail
);

extern void __stdcall F11JEF(
  CONST char *method, CONST int method_len,
  CONST char *precon, CONST int precon_len,
  CONST int *n,
  CONST int *nnz,
  CONST double a[],
  CONST int irow[],
  CONST int icol[],
  CONST double *omega,
  CONST double b[],
  CONST double *tol,
  CONST int *maxitn,
  double x[],
  double *rnorm,
  int *itn,
  double work[],
  CONST int *lwork,
  int iwork[],
  int *ifail
);

extern void __stdcall F11XEF(
  CONST int *n,
  CONST int *nnz,
  CONST double a[],
  CONST int irow[],
  CONST int icol[],
  CONST char *check, CONST int check_len,
  CONST double x[],
  double y[],
  int *ifail
);

extern void __stdcall F11XAF(
  CONST char *trans, CONST int trans_len,
  CONST int *n,
  CONST int *nnz,
  CONST double a[],
  CONST int irow[],
  CONST int icol[],
  CONST char *check, CONST int check_len,
  CONST double x[],
  double y[],
  int *ifail
);

extern void __stdcall F11ZAF(
  CONST int *n,
  int *nnz,
  double a[],
  int irow[],
  int icol[],
  CONST char *dup, CONST int dup_len,
  CONST char *zero, CONST int zero_len,
  int istr[],
  int iwork[],
  int *ifail
);

extern void __stdcall F11ZBF(
  CONST int *n,
  int *nnz,
  double a[],
  int irow[],
  int icol[],
  CONST char *dup, CONST int dup_len,
  CONST char *zero, CONST int zero_len,
  int istr[],
  int iwork[],
  int *ifail
);

extern double __stdcall ddot(
  CONST int *n,
  CONST double x[],
  CONST int *incx,
  CONST double y[],
  CONST int *incy
);

extern void __stdcall drotg(
  double *a,
  double *b,
  double *c,
  double *s
);

extern void __stdcall daxpy(
  CONST int *n,
  CONST double *alpha,
  CONST double x[],
  CONST int *incx,
  double y[],
  CONST int *incy
);

extern void __stdcall dscal(
  CONST int *n,
  CONST double *alpha,
  double x[],
  CONST int *incx
);

extern void __stdcall dcopy(
  CONST int *n,
  CONST double x[],
  CONST int *incx,
  double y[],
  CONST int *incy
);

extern void __stdcall dswap(
  CONST int *n,
  double x[],
  CONST int *incx,
  double y[],
  CONST int *incy
);

extern double __stdcall dnrm2(
  CONST int *n,
  CONST double x[],
  CONST int *incx
);

extern double __stdcall dasum(
  CONST int *n,
  CONST double x[],
  CONST int *incx
);

extern void __stdcall drot(
  CONST int *n,
  double x[],
  CONST int *incx,
  double y[],
  CONST int *incy,
  CONST double *c,
  CONST double *s
);

extern double __stdcall ddoti(
  CONST int *nz,
  CONST double x[],
  CONST int indx[],
  CONST double y[]
);

extern void __stdcall daxpyi(
  CONST int *nz,
  CONST double *a,
  CONST double x[],
  CONST int indx[],
  double y[]
);

extern void __stdcall dgthr(
  CONST int *nz,
  CONST double y[],
  double x[],
  CONST int indx[]
);

extern void __stdcall dgthrz(
  CONST int *nz,
  double y[],
  double x[],
  CONST int indx[]
);

extern void __stdcall dsctr(
  CONST int *nz,
  CONST double x[],
  CONST int indx[],
  double y[]
);

extern void __stdcall droti(
  CONST int *nz,
  double x[],
  CONST int indx[],
  double y[],
  CONST double *c,
  CONST double *s
);

extern void __stdcall ztrtri(
  CONST char *uplo, CONST int uplo_len,
  CONST char *diag, CONST int diag_len,
  CONST int *n,
  Complex a[] /* 2 dimension */,
  CONST int *lda,
  int *info
);

extern void __stdcall dtptrs(
  CONST char *uplo, CONST int uplo_len,
  CONST char *trans, CONST int trans_len,
  CONST char *diag, CONST int diag_len,
  CONST int *n,
  CONST int *nrhs,
  CONST double ap[],
  double b[] /* 2 dimension */,
  CONST int *ldb,
  int *info
);

extern void __stdcall dtpcon(
  CONST char *norm, CONST int norm_len,
  CONST char *uplo, CONST int uplo_len,
  CONST char *diag, CONST int diag_len,
  CONST int *n,
  CONST double ap[],
  double *rcond,
  double work[],
  int iwork[],
  int *info
);

extern void __stdcall dtprfs(
  CONST char *uplo, CONST int uplo_len,
  CONST char *trans, CONST int trans_len,
  CONST char *diag, CONST int diag_len,
  CONST int *n,
  CONST int *nrhs,
  CONST double ap[],
  CONST double b[] /* 2 dimension */,
  CONST int *ldb,
  CONST double x[] /* 2 dimension */,
  CONST int *ldx,
  double ferr[],
  double berr[],
  double work[],
  int iwork[],
  int *info
);

extern void __stdcall dtptri(
  CONST char *uplo, CONST int uplo_len,
  CONST char *diag, CONST int diag_len,
  CONST int *n,
  double ap[],
  int *info
);

extern void __stdcall ztptrs(
  CONST char *uplo, CONST int uplo_len,
  CONST char *trans, CONST int trans_len,
  CONST char *diag, CONST int diag_len,
  CONST int *n,
  CONST int *nrhs,
  CONST Complex ap[],
  Complex b[] /* 2 dimension */,
  CONST int *ldb,
  int *info
);

extern void __stdcall ztpcon(
  CONST char *norm, CONST int norm_len,
  CONST char *uplo, CONST int uplo_len,
  CONST char *diag, CONST int diag_len,
  CONST int *n,
  CONST Complex ap[],
  double *rcond,
  Complex work[],
  double rwork[],
  int *info
);

extern void __stdcall ztprfs(
  CONST char *uplo, CONST int uplo_len,
  CONST char *trans, CONST int trans_len,
  CONST char *diag, CONST int diag_len,
  CONST int *n,
  CONST int *nrhs,
  CONST Complex ap[],
  CONST Complex b[] /* 2 dimension */,
  CONST int *ldb,
  CONST Complex x[] /* 2 dimension */,
  CONST int *ldx,
  double ferr[],
  double berr[],
  Complex work[],
  double rwork[],
  int *info
);

extern void __stdcall ztptri(
  CONST char *uplo, CONST int uplo_len,
  CONST char *diag, CONST int diag_len,
  CONST int *n,
  Complex ap[],
  int *info
);

extern void __stdcall dtbtrs(
  CONST char *uplo, CONST int uplo_len,
  CONST char *trans, CONST int trans_len,
  CONST char *diag, CONST int diag_len,
  CONST int *n,
  CONST int *kd,
  CONST int *nrhs,
  CONST double ab[] /* 2 dimension */,
  CONST int *ldab,
  double b[] /* 2 dimension */,
  CONST int *ldb,
  int *info
);

extern void __stdcall dtbcon(
  CONST char *norm, CONST int norm_len,
  CONST char *uplo, CONST int uplo_len,
  CONST char *diag, CONST int diag_len,
  CONST int *n,
  CONST int *kd,
  CONST double ab[] /* 2 dimension */,
  CONST int *ldab,
  double *rcond,
  double work[],
  int iwork[],
  int *info
);

extern void __stdcall dtbrfs(
  CONST char *uplo, CONST int uplo_len,
  CONST char *trans, CONST int trans_len,
  CONST char *diag, CONST int diag_len,
  CONST int *n,
  CONST int *kd,
  CONST int *nrhs,
  CONST double ab[] /* 2 dimension */,
  CONST int *ldab,
  CONST double b[] /* 2 dimension */,
  CONST int *ldb,
  CONST double x[] /* 2 dimension */,
  CONST int *ldx,
  double ferr[],
  double berr[],
  double work[],
  int iwork[],
  int *info
);

extern void __stdcall ztbtrs(
  CONST char *uplo, CONST int uplo_len,
  CONST char *trans, CONST int trans_len,
  CONST char *diag, CONST int diag_len,
  CONST int *n,
  CONST int *kd,
  CONST int *nrhs,
  CONST Complex ab[] /* 2 dimension */,
  CONST int *ldab,
  Complex b[] /* 2 dimension */,
  CONST int *ldb,
  int *info
);

extern void __stdcall ztbcon(
  CONST char *norm, CONST int norm_len,
  CONST char *uplo, CONST int uplo_len,
  CONST char *diag, CONST int diag_len,
  CONST int *n,
  CONST int *kd,
  CONST Complex ab[] /* 2 dimension */,
  CONST int *ldab,
  double *rcond,
  Complex work[],
  double rwork[],
  int *info
);

extern void __stdcall ztbrfs(
  CONST char *uplo, CONST int uplo_len,
  CONST char *trans, CONST int trans_len,
  CONST char *diag, CONST int diag_len,
  CONST int *n,
  CONST int *kd,
  CONST int *nrhs,
  CONST Complex ab[] /* 2 dimension */,
  CONST int *ldab,
  CONST Complex b[] /* 2 dimension */,
  CONST int *ldb,
  CONST Complex x[] /* 2 dimension */,
  CONST int *ldx,
  double ferr[],
  double berr[],
  Complex work[],
  double rwork[],
  int *info
);

extern void __stdcall dgeqrf(
  CONST int *m,
  CONST int *n,
  double a[] /* 2 dimension */,
  CONST int *lda,
  double tau[],
  double work[],
  CONST int *lwork,
  int *info
);

extern void __stdcall dorgqr(
  CONST int *m,
  CONST int *n,
  CONST int *k,
  double a[] /* 2 dimension */,
  CONST int *lda,
  CONST double tau[],
  double work[],
  CONST int *lwork,
  int *info
);

extern void __stdcall dormqr(
  CONST char *side, CONST int side_len,
  CONST char *trans, CONST int trans_len,
  CONST int *m,
  CONST int *n,
  CONST int *k,
  double a[] /* 2 dimension */,
  CONST int *lda,
  CONST double tau[],
  double c[] /* 2 dimension */,
  CONST int *ldc,
  double work[],
  CONST int *lwork,
  int *info
);

extern void __stdcall dgelqf(
  CONST int *m,
  CONST int *n,
  double a[] /* 2 dimension */,
  CONST int *lda,
  double tau[],
  double work[],
  CONST int *lwork,
  int *info
);

extern void __stdcall dorglq(
  CONST int *m,
  CONST int *n,
  CONST int *k,
  double a[] /* 2 dimension */,
  CONST int *lda,
  CONST double tau[],
  double work[],
  CONST int *lwork,
  int *info
);

extern void __stdcall dormlq(
  CONST char *side, CONST int side_len,
  CONST char *trans, CONST int trans_len,
  CONST int *m,
  CONST int *n,
  CONST int *k,
  double a[] /* 2 dimension */,
  CONST int *lda,
  CONST double tau[],
  double c[] /* 2 dimension */,
  CONST int *ldc,
  double work[],
  CONST int *lwork,
  int *info
);

extern void __stdcall zgeqrf(
  CONST int *m,
  CONST int *n,
  Complex a[] /* 2 dimension */,
  CONST int *lda,
  Complex tau[],
  Complex work[],
  CONST int *lwork,
  int *info
);

extern void __stdcall zungqr(
  CONST int *m,
  CONST int *n,
  CONST int *k,
  Complex a[] /* 2 dimension */,
  CONST int *lda,
  CONST Complex tau[],
  Complex work[],
  CONST int *lwork,
  int *info
);

extern void __stdcall zunmqr(
  CONST char *side, CONST int side_len,
  CONST char *trans, CONST int trans_len,
  CONST int *m,
  CONST int *n,
  CONST int *k,
  Complex a[] /* 2 dimension */,
  CONST int *lda,
  CONST Complex tau[],
  Complex c[] /* 2 dimension */,
  CONST int *ldc,
  Complex work[],
  CONST int *lwork,
  int *info
);

extern void __stdcall zgelqf(
  CONST int *m,
  CONST int *n,
  Complex a[] /* 2 dimension */,
  CONST int *lda,
  Complex tau[],
  Complex work[],
  CONST int *lwork,
  int *info
);

extern void __stdcall zunglq(
  CONST int *m,
  CONST int *n,
  CONST int *k,
  Complex a[] /* 2 dimension */,
  CONST int *lda,
  CONST Complex tau[],
  Complex work[],
  CONST int *lwork,
  int *info
);

extern void __stdcall zunmlq(
  CONST char *side, CONST int side_len,
  CONST char *trans, CONST int trans_len,
  CONST int *m,
  CONST int *n,
  CONST int *k,
  Complex a[] /* 2 dimension */,
  CONST int *lda,
  CONST Complex tau[],
  Complex c[] /* 2 dimension */,
  CONST int *ldc,
  Complex work[],
  CONST int *lwork,
  int *info
);

extern void __stdcall dgeqpf(
  CONST int *m,
  CONST int *n,
  double a[] /* 2 dimension */,
  CONST int *lda,
  int jpvt[],
  double tau[],
  double work[],
  int *info
);

extern void __stdcall zgeqpf(
  CONST int *m,
  CONST int *n,
  Complex a[] /* 2 dimension */,
  CONST int *lda,
  int jpvt[],
  Complex tau[],
  Complex work[],
  double rwork[],
  int *info
);

extern void __stdcall dsytrd(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  double a[] /* 2 dimension */,
  CONST int *lda,
  double d[],
  double e[],
  double tau[],
  double work[],
  CONST int *lwork,
  int *info
);

extern void __stdcall dorgtr(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  double a[] /* 2 dimension */,
  CONST int *lda,
  CONST double tau[],
  double work[],
  CONST int *lwork,
  int *info
);

extern void __stdcall dormtr(
  CONST char *side, CONST int side_len,
  CONST char *uplo, CONST int uplo_len,
  CONST char *trans, CONST int trans_len,
  CONST int *m,
  CONST int *n,
  double a[] /* 2 dimension */,
  CONST int *lda,
  CONST double tau[],
  double c[] /* 2 dimension */,
  CONST int *ldc,
  double work[],
  CONST int *lwork,
  int *info
);

extern void __stdcall zhetrd(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  Complex a[] /* 2 dimension */,
  CONST int *lda,
  double d[],
  double e[],
  Complex tau[],
  Complex work[],
  CONST int *lwork,
  int *info
);

extern void __stdcall zungtr(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  Complex a[] /* 2 dimension */,
  CONST int *lda,
  CONST Complex tau[],
  Complex work[],
  CONST int *lwork,
  int *info
);

extern void __stdcall zunmtr(
  CONST char *side, CONST int side_len,
  CONST char *uplo, CONST int uplo_len,
  CONST char *trans, CONST int trans_len,
  CONST int *m,
  CONST int *n,
  Complex a[] /* 2 dimension */,
  CONST int *lda,
  CONST Complex tau[],
  Complex c[] /* 2 dimension */,
  CONST int *ldc,
  Complex work[],
  CONST int *lwork,
  int *info
);

extern void __stdcall dsptrd(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  double ap[],
  double d[],
  double e[],
  double tau[],
  int *info
);

extern void __stdcall dopgtr(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST double ap[],
  CONST double tau[],
  double q[] /* 2 dimension */,
  CONST int *ldq,
  double work[],
  int *info
);

extern void __stdcall dopmtr(
  CONST char *side, CONST int side_len,
  CONST char *uplo, CONST int uplo_len,
  CONST char *trans, CONST int trans_len,
  CONST int *m,
  CONST int *n,
  double ap[],
  CONST double tau[],
  double c[] /* 2 dimension */,
  CONST int *ldc,
  double work[],
  int *info
);

extern void __stdcall zhptrd(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  Complex ap[],
  double d[],
  double e[],
  Complex tau[],
  int *info
);

extern void __stdcall zupgtr(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST Complex ap[],
  CONST Complex tau[],
  Complex q[] /* 2 dimension */,
  CONST int *ldq,
  Complex work[],
  int *info
);

extern void __stdcall zupmtr(
  CONST char *side, CONST int side_len,
  CONST char *uplo, CONST int uplo_len,
  CONST char *trans, CONST int trans_len,
  CONST int *m,
  CONST int *n,
  Complex ap[],
  CONST Complex tau[],
  Complex c[] /* 2 dimension */,
  CONST int *ldc,
  Complex work[],
  int *info
);

extern void __stdcall dsbtrd(
  CONST char *vect, CONST int vect_len,
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST int *kd,
  double ab[] /* 2 dimension */,
  CONST int *ldab,
  double d[],
  double e[],
  double q[] /* 2 dimension */,
  CONST int *ldq,
  double work[],
  int *info
);

extern void __stdcall zhbtrd(
  CONST char *vect, CONST int vect_len,
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST int *kd,
  Complex ab[] /* 2 dimension */,
  CONST int *ldab,
  double d[],
  double e[],
  Complex q[] /* 2 dimension */,
  CONST int *ldq,
  Complex work[],
  int *info
);

extern void __stdcall dsteqr(
  CONST char *compz, CONST int compz_len,
  CONST int *n,
  double d[],
  double e[],
  double z[] /* 2 dimension */,
  CONST int *ldz,
  double work[],
  int *info
);

extern void __stdcall dsterf(
  CONST int *n,
  double d[],
  double e[],
  int *info
);

extern void __stdcall dpteqr(
  CONST char *compz, CONST int compz_len,
  CONST int *n,
  double d[],
  double e[],
  double z[] /* 2 dimension */,
  CONST int *ldz,
  double work[],
  int *info
);

extern void __stdcall dstebz(
  CONST char *range, CONST int range_len,
  CONST char *order, CONST int order_len,
  CONST int *n,
  CONST double *vl,
  CONST double *vu,
  CONST int *il,
  CONST int *iu,
  CONST double *abstol,
  CONST double d[],
  CONST double e[],
  int *m,
  int *nsplit,
  double w[],
  int iblock[],
  int isplit[],
  double work[],
  int iwork[],
  int *info
);

extern void __stdcall dstein(
  CONST int *n,
  CONST double d[],
  CONST double e[],
  CONST int *m,
  CONST double w[],
  CONST int iblock[],
  CONST int isplit[],
  double z[] /* 2 dimension */,
  CONST int *ldz,
  double work[],
  int iwork[],
  int ifail[],
  int *info
);

extern void __stdcall zsteqr(
  CONST char *compz, CONST int compz_len,
  CONST int *n,
  double d[],
  double e[],
  Complex z[] /* 2 dimension */,
  CONST int *ldz,
  double work[],
  int *info
);

extern void __stdcall zpteqr(
  CONST char *compz, CONST int compz_len,
  CONST int *n,
  double d[],
  double e[],
  Complex z[] /* 2 dimension */,
  CONST int *ldz,
  double work[],
  int *info
);

extern void __stdcall zstein(
  CONST int *n,
  CONST double d[],
  CONST double e[],
  CONST int *m,
  CONST double w[],
  CONST int iblock[],
  CONST int isplit[],
  Complex z[] /* 2 dimension */,
  CONST int *ldz,
  double work[],
  int iwork[],
  int ifail[],
  int *info
);

extern void __stdcall dgebrd(
  CONST int *m,
  CONST int *n,
  double a[] /* 2 dimension */,
  CONST int *lda,
  double d[],
  double e[],
  double tauq[],
  double taup[],
  double work[],
  CONST int *lwork,
  int *info
);

extern void __stdcall dorgbr(
  CONST char *vect, CONST int vect_len,
  CONST int *m,
  CONST int *n,
  CONST int *k,
  double a[] /* 2 dimension */,
  CONST int *lda,
  CONST double tau[],
  double work[],
  CONST int *lwork,
  int *info
);

extern void __stdcall dormbr(
  CONST char *vect, CONST int vect_len,
  CONST char *side, CONST int side_len,
  CONST char *trans, CONST int trans_len,
  CONST int *m,
  CONST int *n,
  CONST int *k,
  CONST double a[] /* 2 dimension */,
  CONST int *lda,
  CONST double tau[],
  double c[] /* 2 dimension */,
  CONST int *ldc,
  double work[],
  CONST int *lwork,
  int *info
);

extern void __stdcall zgebrd(
  CONST int *m,
  CONST int *n,
  Complex a[] /* 2 dimension */,
  CONST int *lda,
  double d[],
  double e[],
  Complex tauq[],
  Complex taup[],
  Complex work[],
  CONST int *lwork,
  int *info
);

extern void __stdcall zungbr(
  CONST char *vect, CONST int vect_len,
  CONST int *m,
  CONST int *n,
  CONST int *k,
  Complex a[] /* 2 dimension */,
  CONST int *lda,
  CONST Complex tau[],
  Complex work[],
  CONST int *lwork,
  int *info
);

extern void __stdcall zunmbr(
  CONST char *vect, CONST int vect_len,
  CONST char *side, CONST int side_len,
  CONST char *trans, CONST int trans_len,
  CONST int *m,
  CONST int *n,
  CONST int *k,
  CONST Complex a[] /* 2 dimension */,
  CONST int *lda,
  CONST Complex tau[],
  Complex c[] /* 2 dimension */,
  CONST int *ldc,
  Complex work[],
  CONST int *lwork,
  int *info
);

extern void __stdcall dbdsqr(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST int *ncvt,
  CONST int *nru,
  CONST int *ncc,
  double d[],
  double e[],
  double vt[] /* 2 dimension */,
  CONST int *ldvt,
  double u[] /* 2 dimension */,
  CONST int *ldu,
  double c[] /* 2 dimension */,
  CONST int *ldc,
  double work[],
  int *info
);

extern void __stdcall zbdsqr(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST int *ncvt,
  CONST int *nru,
  CONST int *ncc,
  double d[],
  double e[],
  Complex vt[] /* 2 dimension */,
  CONST int *ldvt,
  Complex u[] /* 2 dimension */,
  CONST int *ldu,
  Complex c[] /* 2 dimension */,
  CONST int *ldc,
  double work[],
  int *info
);

extern void __stdcall dgehrd(
  CONST int *n,
  CONST int *ilo,
  CONST int *ihi,
  double a[] /* 2 dimension */,
  CONST int *lda,
  double tau[],
  double work[],
  CONST int *lwork,
  int *info
);

extern void __stdcall dorghr(
  CONST int *n,
  CONST int *ilo,
  CONST int *ihi,
  double a[] /* 2 dimension */,
  CONST int *lda,
  CONST double tau[],
  double work[],
  CONST int *lwork,
  int *info
);

extern void __stdcall dormhr(
  CONST char *side, CONST int side_len,
  CONST char *trans, CONST int trans_len,
  CONST int *m,
  CONST int *n,
  CONST int *ilo,
  CONST int *ihi,
  CONST double a[] /* 2 dimension */,
  CONST int *lda,
  CONST double tau[],
  double c[] /* 2 dimension */,
  CONST int *ldc,
  double work[],
  CONST int *lwork,
  int *info
);

extern void __stdcall dgebal(
  CONST char *job, CONST int job_len,
  CONST int *n,
  double a[] /* 2 dimension */,
  CONST int *lda,
  int *ilo,
  int *ihi,
  double scale[],
  int *info
);

extern void __stdcall dgebak(
  CONST char *job, CONST int job_len,
  CONST char *side, CONST int side_len,
  CONST int *n,
  CONST int *ilo,
  CONST int *ihi,
  CONST double scale[],
  CONST int *m,
  double v[] /* 2 dimension */,
  CONST int *ldv,
  int *info
);

extern void __stdcall zgehrd(
  CONST int *n,
  CONST int *ilo,
  CONST int *ihi,
  Complex a[] /* 2 dimension */,
  CONST int *lda,
  Complex tau[],
  Complex work[],
  CONST int *lwork,
  int *info
);

extern void __stdcall zunghr(
  CONST int *n,
  CONST int *ilo,
  CONST int *ihi,
  Complex a[] /* 2 dimension */,
  CONST int *lda,
  CONST Complex tau[],
  Complex work[],
  CONST int *lwork,
  int *info
);

extern void __stdcall zunmhr(
  CONST char *side, CONST int side_len,
  CONST char *trans, CONST int trans_len,
  CONST int *m,
  CONST int *n,
  CONST int *ilo,
  CONST int *ihi,
  CONST Complex a[] /* 2 dimension */,
  CONST int *lda,
  CONST Complex tau[],
  Complex c[] /* 2 dimension */,
  CONST int *ldc,
  Complex work[],
  CONST int *lwork,
  int *info
);

extern void __stdcall zgebal(
  CONST char *job, CONST int job_len,
  CONST int *n,
  Complex a[] /* 2 dimension */,
  CONST int *lda,
  int *ilo,
  int *ihi,
  double scale[],
  int *info
);

extern void __stdcall zgebak(
  CONST char *job, CONST int job_len,
  CONST char *side, CONST int side_len,
  CONST int *n,
  CONST int *ilo,
  CONST int *ihi,
  CONST double scale[],
  CONST int *m,
  Complex v[] /* 2 dimension */,
  CONST int *ldv,
  int *info
);

extern void __stdcall dhseqr(
  CONST char *job, CONST int job_len,
  CONST char *compz, CONST int compz_len,
  CONST int *n,
  CONST int *ilo,
  CONST int *ihi,
  double h[] /* 2 dimension */,
  CONST int *ldh,
  double wr[],
  double wi[],
  double z[] /* 2 dimension */,
  CONST int *ldz,
  double work[],
  CONST int *lwork,
  int *info
);

extern void __stdcall dhsein(
  CONST char *job, CONST int job_len,
  CONST char *eigsrc, CONST int eigsrc_len,
  CONST char *initv, CONST int initv_len,
  int select[],
  CONST int *n,
  CONST double h[] /* 2 dimension */,
  CONST int *ldh,
  double wr[],
  CONST double wi[],
  double vl[] /* 2 dimension */,
  CONST int *ldvl,
  double vr[] /* 2 dimension */,
  CONST int *ldvr,
  CONST int *mm,
  int *m,
  double work[],
  int ifaill[],
  int ifailr[],
  int *info
);

extern void __stdcall zhseqr(
  CONST char *job, CONST int job_len,
  CONST char *compz, CONST int compz_len,
  CONST int *n,
  CONST int *ilo,
  CONST int *ihi,
  Complex h[] /* 2 dimension */,
  CONST int *ldh,
  Complex w[],
  Complex z[] /* 2 dimension */,
  CONST int *ldz,
  Complex work[],
  CONST int *lwork,
  int *info
);

extern void __stdcall zhsein(
  CONST char *job, CONST int job_len,
  CONST char *eigsrc, CONST int eigsrc_len,
  CONST char *initv, CONST int initv_len,
  CONST int select[],
  CONST int *n,
  CONST Complex h[] /* 2 dimension */,
  CONST int *ldh,
  Complex w[],
  Complex vl[] /* 2 dimension */,
  CONST int *ldvl,
  Complex vr[] /* 2 dimension */,
  CONST int *ldvr,
  CONST int *mm,
  int *m,
  Complex work[],
  double rwork[],
  int ifaill[],
  int ifailr[],
  int *info
);

extern void __stdcall dtrexc(
  CONST char *compq, CONST int compq_len,
  CONST int *n,
  double t[] /* 2 dimension */,
  CONST int *ldt,
  double q[] /* 2 dimension */,
  CONST int *ldq,
  int *ifst,
  int *ilst,
  double work[],
  int *info
);

extern void __stdcall dtrsen(
  CONST char *job, CONST int job_len,
  CONST char *compq, CONST int compq_len,
  CONST int select[],
  CONST int *n,
  double t[] /* 2 dimension */,
  CONST int *ldt,
  CONST double q[] /* 2 dimension */,
  CONST int *ldq,
  double wr[],
  double wi[],
  int *m,
  double *s,
  double *sep,
  double work[],
  CONST int *lwork,
  int iwork[],
  CONST int *liwork,
  int *info
);

extern void __stdcall dtrsyl(
  CONST char *trana, CONST int trana_len,
  CONST char *tranb, CONST int tranb_len,
  CONST int *isgn,
  CONST int *m,
  CONST int *n,
  CONST double a[] /* 2 dimension */,
  CONST int *lda,
  CONST double b[] /* 2 dimension */,
  CONST int *ldb,
  double c[] /* 2 dimension */,
  CONST int *ldc,
  double *scale,
  int *info
);

extern void __stdcall dtrevc(
  CONST char *job, CONST int job_len,
  CONST char *howmny, CONST int howmny_len,
  int select[],
  CONST int *n,
  CONST double t[] /* 2 dimension */,
  CONST int *ldt,
  double vl[] /* 2 dimension */,
  CONST int *ldvl,
  double vr[] /* 2 dimension */,
  CONST int *ldvr,
  CONST int *mm,
  int *m,
  double work[],
  int *info
);

extern void __stdcall dtrsna(
  CONST char *job, CONST int job_len,
  CONST char *howmny, CONST int howmny_len,
  CONST int select[],
  CONST int *n,
  CONST double t[] /* 2 dimension */,
  CONST int *ldt,
  CONST double vl[] /* 2 dimension */,
  CONST int *ldvl,
  CONST double vr[] /* 2 dimension */,
  CONST int *ldvr,
  double s[],
  double sep[],
  CONST int *mm,
  int *m,
  double work[] /* 2 dimension */,
  CONST int *ldwork,
  int iwork[],
  int *info
);

extern void __stdcall ztrexc(
  CONST char *compq, CONST int compq_len,
  CONST int *n,
  Complex t[] /* 2 dimension */,
  CONST int *ldt,
  Complex q[] /* 2 dimension */,
  CONST int *ldq,
  CONST int *ifst,
  CONST int *ilst,
  int *info
);

extern void __stdcall ztrsen(
  CONST char *job, CONST int job_len,
  CONST char *compq, CONST int compq_len,
  CONST int select[],
  CONST int *n,
  Complex t[] /* 2 dimension */,
  CONST int *ldt,
  Complex q[] /* 2 dimension */,
  CONST int *ldq,
  Complex w[],
  int *m,
  double *s,
  double *sep,
  Complex work[],
  CONST int *lwork,
  int *info
);

extern void __stdcall ztrsyl(
  CONST char *trana, CONST int trana_len,
  CONST char *tranb, CONST int tranb_len,
  CONST int *isgn,
  CONST int *m,
  CONST int *n,
  CONST Complex a[] /* 2 dimension */,
  CONST int *lda,
  CONST Complex b[] /* 2 dimension */,
  CONST int *ldb,
  Complex c[] /* 2 dimension */,
  CONST int *ldc,
  double *scale,
  int *info
);

extern void __stdcall ztrevc(
  CONST char *job, CONST int job_len,
  CONST char *howmny, CONST int howmny_len,
  CONST int select[],
  CONST int *n,
  Complex t[] /* 2 dimension */,
  CONST int *ldt,
  Complex vl[] /* 2 dimension */,
  CONST int *ldvl,
  Complex vr[] /* 2 dimension */,
  CONST int *ldvr,
  CONST int *mm,
  int *m,
  Complex work[],
  double rwork[],
  int *info
);

extern void __stdcall ztrsna(
  CONST char *job, CONST int job_len,
  CONST char *howmny, CONST int howmny_len,
  CONST int select[],
  CONST int *n,
  CONST Complex t[] /* 2 dimension */,
  CONST int *ldt,
  CONST Complex vl[] /* 2 dimension */,
  CONST int *ldvl,
  CONST Complex vr[] /* 2 dimension */,
  CONST int *ldvr,
  double s[],
  double sep[],
  CONST int *mm,
  int *m,
  Complex work[] /* 2 dimension */,
  CONST int *ldwork,
  double rwork[],
  int *info
);

extern void __stdcall dsygst(
  CONST int *itype,
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  double a[] /* 2 dimension */,
  CONST int *lda,
  CONST double b[] /* 2 dimension */,
  CONST int *ldb,
  int *info
);

extern void __stdcall zhegst(
  CONST int *itype,
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  Complex a[] /* 2 dimension */,
  CONST int *lda,
  Complex b[] /* 2 dimension */,
  CONST int *ldb,
  int *info
);

extern void __stdcall dspgst(
  CONST int *itype,
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  double ap[],
  CONST double bp[],
  int *info
);

extern void __stdcall zhpgst(
  CONST int *itype,
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  Complex ap[],
  CONST Complex bp[],
  int *info
);

extern void __stdcall G01AAF(
  CONST int *n,
  CONST double x[],
  int *iwt,
  double wt[],
  double *xmean,
  double *s2,
  double *s3,
  double *s4,
  double *xmin,
  double *xmax,
  double *wsum,
  int *ifail
);

extern void __stdcall G01ABF(
  CONST int *n,
  CONST double x1[],
  CONST double x2[],
  int *iwt,
  double wt[],
  double res[],
  int *ifail
);

extern void __stdcall G01ADF(
  CONST int *k,
  CONST double x[],
  CONST int ifreq[],
  double *xmean,
  double *s2,
  double *s3,
  double *s4,
  int *n,
  int *ifail
);

extern void __stdcall G01AEF(
  CONST int *n,
  CONST int *k2,
  CONST double x[],
  CONST int *iclass,
  double cint[],
  int ifreq[],
  double *xmin,
  double *xmax,
  int *ifail
);

extern void __stdcall G01AFF(
  CONST int *inob,
  CONST int *ipred,
  CONST int *m,
  CONST int *n,
  int nobs[] /* 2 dimension */,
  int *num,
  double pred[] /* 2 dimension */,
  double *chis,
  double p[],
  int *npos,
  int *ndf,
  int *m1,
  int *n1,
  int *ifail
);

extern void __stdcall G01AGF(
  CONST double x[],
  double y[],
  CONST int *nobs,
  int isort[],
  CONST int *nstepx,
  CONST int *nstepy,
  int *ifail
);

extern void __stdcall G01AHF(
  CONST double x[],
  CONST int *nobs,
  CONST int *nstepx,
  CONST int *nstepy,
  CONST int *istand,
  int iwork[],
  double work[],
  CONST int *lwork,
  double xsort[],
  double *xbar,
  double *xstd,
  int *ifail
);

extern void __stdcall G01AJF(
  CONST double x[],
  CONST int *n,
  int *nstepx,
  int *nstepy,
  CONST int *itype,
  int *ispace,
  double *xmin,
  double *xmax,
  double *xstep,
  int *n1,
  int *multy,
  int *ifail
);

extern void __stdcall G01ALF(
  CONST int *n,
  CONST double x[],
  int iwrk[],
  double res[],
  int *ifail
);

extern void __stdcall G01ARF(
  CONST char *range, CONST int range_len,
  CONST char *prt, CONST int prt_len,
  CONST int *n,
  CONST double y[],
  CONST int *nstepx,
  CONST int *nstepy,
  double *unit,
  char *plot /* 2 dimension */, CONST int plot_len,
  CONST int *ldp,
  int *lines,
  double sorty[],
  int iwork[],
  int *ifail
);

extern void __stdcall G01ASF(
  CONST char *prt, CONST int prt_len,
  CONST int *m,
  CONST int n[],
  CONST double x[] /* 2 dimension */,
  CONST int *ldx,
  CONST int *nstepx,
  CONST int *nstepy,
  char *plot /* 2 dimension */, CONST int plot_len,
  CONST int *ldp,
  double work[],
  int iwork[],
  int *ifail
);

extern double __stdcall G01BAF(
  CONST int *idf,
  CONST double *t,
  int *ifail
);

extern double __stdcall G01BBF(
  CONST int *i1,
  CONST int *i2,
  CONST double *a,
  int *ifail
);

extern double __stdcall G01BCF(
  CONST double *x,
  CONST int *n,
  int *ifail
);

extern double __stdcall G01BDF(
  CONST double *x,
  CONST double *a,
  CONST double *b,
  int *ifail
);

extern void __stdcall G01BJF(
  CONST int *n,
  CONST double *p,
  CONST int *k,
  double *plek,
  double *pgtk,
  double *peqk,
  int *ifail
);

extern void __stdcall G01BKF(
  CONST double *rlamda,
  CONST int *k,
  double *plek,
  double *pgtk,
  double *peqk,
  int *ifail
);

extern void __stdcall G01BLF(
  CONST int *n,
  CONST int *l,
  CONST int *m,
  CONST int *k,
  double *plek,
  double *pgtk,
  double *peqk,
  int *ifail
);

extern double __stdcall G01CAF(
  CONST double *p,
  CONST int *n,
  int *ifail
);

extern double __stdcall G01CBF(
  CONST double *p,
  CONST int *m,
  CONST int *n,
  int *ifail
);

extern double __stdcall G01CCF(
  CONST double *p,
  CONST int *n,
  int *ifail
);

extern double __stdcall G01CDF(
  CONST double *p,
  CONST double *a,
  CONST double *b,
  int *ifail
);

extern double __stdcall G01CEF(
  CONST double *p,
  int *ifail
);

extern void __stdcall G01DAF(
  CONST int *n,
  double pp[],
  CONST double *etol,
  double *errest,
  double work[],
  CONST int *iw,
  int *ifail
);

extern void __stdcall G01DBF(
  CONST int *n,
  double pp[],
  int *ifail
);

extern void __stdcall G01DCF(
  CONST int *n,
  CONST double *ex1,
  CONST double *ex2,
  CONST double *summ2,
  double vapvec[],
  int *ifail
);

extern void __stdcall G01DDF(
  CONST double x[],
  CONST int *n,
  CONST int *calwts,
  double a[],
  double *w,
  double *pw,
  int *ifail
);

extern void __stdcall G01DHF(
  CONST char *scores, CONST int scores_len,
  CONST char *ties, CONST int ties_len,
  CONST int *n,
  CONST double x[],
  double r[],
  int iwrk[],
  int *ifail
);

extern double __stdcall G01EAF(
  CONST char *tail, CONST int tail_len,
  CONST double *x,
  int *ifail
);

extern double __stdcall G01EBF(
  CONST char *tail, CONST int tail_len,
  CONST double *t,
  CONST double *df,
  int *ifail
);

extern double __stdcall G01ECF(
  CONST char *tail, CONST int tail_len,
  CONST double *x,
  CONST double *df,
  int *ifail
);

extern double __stdcall G01EDF(
  CONST char *tail, CONST int tail_len,
  CONST double *f,
  CONST double *df1,
  CONST double *df2,
  int *ifail
);

extern void __stdcall G01EEF(
  CONST double *x,
  CONST double *a,
  CONST double *b,
  CONST double *tol,
  double *p,
  double *q,
  double *pdf,
  int *ifail
);

extern double __stdcall G01EFF(
  CONST char *tail, CONST int tail_len,
  CONST double *g,
  CONST double *a,
  CONST double *b,
  int *ifail
);

extern double __stdcall G01EMF(
  CONST double *q,
  CONST double *v,
  CONST int *ir,
  int *ifail
);

extern void __stdcall G01EPF(
  CONST int *n,
  CONST int *ip,
  CONST double *d,
  double *pdl,
  double *pdu,
  double work[],
  int *ifail
);

extern double __stdcall G01ERF(
  CONST double *t,
  CONST double *vk,
  int *ifail
);

extern double __stdcall G01EYF(
  CONST int *n,
  CONST double *d,
  int *ifail
);

extern double __stdcall G01EZF(
  CONST int *n1,
  CONST int *n2,
  CONST double *d,
  int *ifail
);

extern double __stdcall G01FAF(
  CONST char *tail, CONST int tail_len,
  CONST double *p,
  int *ifail
);

extern double __stdcall G01FBF(
  CONST char *tail, CONST int tail_len,
  CONST double *p,
  CONST double *df,
  int *ifail
);

extern double __stdcall G01FCF(
  CONST double *p,
  CONST double *df,
  int *ifail
);

extern double __stdcall G01FDF(
  CONST double *p,
  CONST double *df1,
  CONST double *df2,
  int *ifail
);

extern double __stdcall G01FEF(
  CONST double *p,
  CONST double *a,
  CONST double *b,
  CONST double *tol,
  int *ifail
);

extern double __stdcall G01FFF(
  CONST double *p,
  CONST double *a,
  CONST double *b,
  CONST double *tol,
  int *ifail
);

extern double __stdcall G01FMF(
  CONST double *p,
  CONST double *v,
  CONST int *ir,
  int *ifail
);

extern double __stdcall G01GBF(
  CONST double *t,
  CONST double *df,
  CONST double *delta,
  CONST double *tol,
  CONST int *maxit,
  int *ifail
);

extern double __stdcall G01GCF(
  CONST double *x,
  CONST double *df,
  CONST double *rlamda,
  CONST double *tol,
  CONST int *maxit,
  int *ifail
);

extern double __stdcall G01GDF(
  CONST double *f,
  CONST double *df1,
  CONST double *df2,
  CONST double *rlamda,
  CONST double *tol,
  CONST int *maxit,
  int *ifail
);

extern double __stdcall G01GEF(
  CONST double *x,
  CONST double *a,
  CONST double *b,
  CONST double *rlamda,
  CONST double *tol,
  CONST int *maxit,
  int *ifail
);

extern double __stdcall G01HAF(
  CONST double *x,
  CONST double *y,
  CONST double *rho,
  int *ifail
);

extern double __stdcall G01HBF(
  CONST char *tail, CONST int tail_len,
  CONST int *n,
  CONST double a[],
  CONST double b[],
  CONST double xmu[],
  CONST double sig[] /* 2 dimension */,
  CONST int *ldsig,
  CONST double *tol,
  double wk[],
  CONST int *lwk,
  int *ifail
);

extern void __stdcall G01JCF(
  CONST double a[],
  CONST int mult[],
  CONST double rlamda[],
  CONST int *n,
  CONST double *c,
  double *p,
  double *pdf,
  CONST double *tol,
  CONST int *maxit,
  double wrk[],
  int *ifail
);

extern void __stdcall G01JDF(
  CONST char *method, CONST int method_len,
  CONST int *n,
  CONST double rlam[],
  CONST double *d,
  CONST double *c,
  double *prob,
  double work[],
  int *ifail
);

extern double __stdcall G01MBF(
  CONST double *x
);

extern void __stdcall G01NAF(
  CONST char *mom, CONST int mom_len,
  CONST char *mean, CONST int mean_len,
  CONST int *n,
  CONST double a[] /* 2 dimension */,
  CONST int *lda,
  CONST double emu[],
  CONST double sigma[] /* 2 dimension */,
  CONST int *ldsig,
  CONST int *l,
  double rkum[],
  double rmom[],
  double wk[],
  int *ifail
);

extern void __stdcall G01NBF(
  CONST char *case1, CONST int case1_len,
  CONST char *mean, CONST int mean_len,
  CONST int *n,
  CONST double a[] /* 2 dimension */,
  CONST int *lda,
  CONST double b[] /* 2 dimension */,
  CONST int *ldb,
  CONST double c[] /* 2 dimension */,
  CONST int *ldc,
  CONST double ela[],
  CONST double emu[],
  CONST double sigma[] /* 2 dimension */,
  CONST int *ldsig,
  CONST int *l1,
  CONST int *l2,
  int *lmax,
  double rmom[],
  double *abserr,
  CONST double *eps,
  double wk[],
  int *ifail
);

extern void __stdcall G02BAF(
  CONST int *n,
  CONST int *m,
  CONST double x[] /* 2 dimension */,
  CONST int *ix,
  double xbar[],
  double std[],
  double ssp[] /* 2 dimension */,
  CONST int *issp,
  double r[] /* 2 dimension */,
  CONST int *ir,
  int *ifail
);

extern void __stdcall G02BBF(
  CONST int *n,
  CONST int *m,
  CONST double x[] /* 2 dimension */,
  CONST int *ix,
  int miss[],
  double xmiss[],
  double xbar[],
  double std[],
  double ssp[] /* 2 dimension */,
  CONST int *issp,
  double r[] /* 2 dimension */,
  CONST int *ir,
  int *ncases,
  int *ifail
);

extern void __stdcall G02BCF(
  CONST int *n,
  CONST int *m,
  CONST double x[] /* 2 dimension */,
  CONST int *ix,
  CONST int miss[],
  CONST double xmiss[],
  double xbar[],
  double std[],
  double ssp[] /* 2 dimension */,
  CONST int *issp,
  double r[] /* 2 dimension */,
  CONST int *ir,
  int *ncases,
  double count[] /* 2 dimension */,
  CONST int *ic,
  int *ifail
);

extern void __stdcall G02BDF(
  CONST int *n,
  CONST int *m,
  CONST double x[] /* 2 dimension */,
  CONST int *ix,
  double xbar[],
  double std[],
  double sspz[] /* 2 dimension */,
  CONST int *isspz,
  double rz[] /* 2 dimension */,
  CONST int *irz,
  int *ifail
);

extern void __stdcall G02BEF(
  CONST int *n,
  CONST int *m,
  CONST double x[] /* 2 dimension */,
  CONST int *ix,
  int miss[],
  double xmiss[],
  double xbar[],
  double std[],
  double sspz[] /* 2 dimension */,
  CONST int *isspz,
  double rz[] /* 2 dimension */,
  CONST int *irz,
  int *ncases,
  int *ifail
);

extern void __stdcall G02BFF(
  CONST int *n,
  CONST int *m,
  CONST double x[] /* 2 dimension */,
  CONST int *ix,
  CONST int miss[],
  CONST double xmiss[],
  double xbar[],
  double std[],
  double sspz[] /* 2 dimension */,
  CONST int *isspz,
  double rz[] /* 2 dimension */,
  CONST int *irz,
  int *ncases,
  double count[] /* 2 dimension */,
  CONST int *ic,
  int *ifail
);

extern void __stdcall G02BGF(
  CONST int *n,
  CONST int *m,
  CONST double x[] /* 2 dimension */,
  CONST int *ix,
  CONST int *nvars,
  CONST int kvar[],
  double xbar[],
  double std[],
  double ssp[] /* 2 dimension */,
  CONST int *issp,
  double r[] /* 2 dimension */,
  CONST int *ir,
  int *ifail
);

extern void __stdcall G02BHF(
  CONST int *n,
  CONST int *m,
  CONST double x[] /* 2 dimension */,
  CONST int *ix,
  int miss[],
  double xmiss[],
  CONST int *mistyp,
  CONST int *nvars,
  CONST int kvar[],
  double xbar[],
  double std[],
  double ssp[] /* 2 dimension */,
  CONST int *issp,
  double r[] /* 2 dimension */,
  CONST int *ir,
  int *ncases,
  int *ifail
);

extern void __stdcall G02BJF(
  CONST int *n,
  CONST int *m,
  CONST double x[] /* 2 dimension */,
  CONST int *ix,
  CONST int miss[],
  CONST double xmiss[],
  CONST int *nvars,
  CONST int kvar[],
  double xbar[],
  double std[],
  double ssp[] /* 2 dimension */,
  CONST int *issp,
  double r[] /* 2 dimension */,
  CONST int *ir,
  int *ncases,
  double count[] /* 2 dimension */,
  CONST int *ic,
  int *ifail
);

extern void __stdcall G02BKF(
  CONST int *n,
  CONST int *m,
  CONST double x[] /* 2 dimension */,
  CONST int *ix,
  CONST int *nvars,
  CONST int kvar[],
  double xbar[],
  double std[],
  double sspz[] /* 2 dimension */,
  CONST int *isspz,
  double rz[] /* 2 dimension */,
  CONST int *irz,
  int *ifail
);

extern void __stdcall G02BLF(
  CONST int *n,
  CONST int *m,
  CONST double x[] /* 2 dimension */,
  CONST int *ix,
  int miss[],
  double xmiss[],
  CONST int *mistyp,
  CONST int *nvars,
  CONST int kvar[],
  double xbar[],
  double std[],
  double sspz[] /* 2 dimension */,
  CONST int *isspz,
  double rz[] /* 2 dimension */,
  CONST int *irz,
  int *ncases,
  int *ifail
);

extern void __stdcall G02BMF(
  CONST int *n,
  CONST int *m,
  CONST double x[] /* 2 dimension */,
  CONST int *ix,
  CONST int miss[],
  CONST double xmiss[],
  CONST int *nvars,
  CONST int kvar[],
  double xbar[],
  double std[],
  double sspz[] /* 2 dimension */,
  CONST int *isspz,
  double rz[] /* 2 dimension */,
  CONST int *irz,
  int *ncases,
  double count[] /* 2 dimension */,
  CONST int *ic,
  int *ifail
);

extern void __stdcall G02BNF(
  CONST int *n,
  CONST int *m,
  double x[] /* 2 dimension */,
  CONST int *ix,
  CONST int *itype,
  double rr[] /* 2 dimension */,
  CONST int *irr,
  int kworka[],
  int kworkb[],
  double work1[],
  double work2[],
  int *ifail
);

extern void __stdcall G02BPF(
  CONST int *n,
  CONST int *m,
  double x[] /* 2 dimension */,
  CONST int *ix,
  int miss[],
  double xmiss[],
  CONST int *itype,
  double rr[] /* 2 dimension */,
  CONST int *irr,
  int *ncases,
  int incase[],
  int kworka[],
  int kworkb[],
  int kworkc[],
  double work1[],
  double work2[],
  int *ifail
);

extern void __stdcall G02BQF(
  CONST int *n,
  CONST int *m,
  CONST double x[] /* 2 dimension */,
  CONST int *ix,
  CONST int *itype,
  double rr[] /* 2 dimension */,
  CONST int *irr,
  int kworka[],
  int kworkb[],
  double work1[],
  double work2[],
  int *ifail
);

extern void __stdcall G02BRF(
  CONST int *n,
  CONST int *m,
  CONST double x[] /* 2 dimension */,
  CONST int *ix,
  int miss[],
  double xmiss[],
  CONST int *itype,
  double rr[] /* 2 dimension */,
  CONST int *irr,
  int *ncases,
  int incase[],
  int kworka[],
  int kworkb[],
  int kworkc[],
  double work1[],
  double work2[],
  int *ifail
);

extern void __stdcall G02BSF(
  CONST int *n,
  CONST int *m,
  CONST double x[] /* 2 dimension */,
  CONST int *ix,
  CONST int miss[],
  CONST double xmiss[],
  CONST int *itype,
  double rr[] /* 2 dimension */,
  CONST int *irr,
  int *ncases,
  double count[] /* 2 dimension */,
  CONST int *ic,
  int kworka[],
  int kworkb[],
  int kworkc[],
  int kworkd[],
  double work1[],
  double work2[],
  int *ifail
);

extern void __stdcall G02BTF(
  CONST char *mean, CONST int mean_len,
  CONST int *m,
  CONST double *wt,
  CONST double x[],
  CONST int *incx,
  double *sw,
  double xbar[],
  double c[],
  int *ifail
);

extern void __stdcall G02BUF(
  CONST char *mean, CONST int mean_len,
  CONST char *weight, CONST int weight_len,
  CONST int *n,
  CONST int *m,
  CONST double x[] /* 2 dimension */,
  CONST int *ldx,
  CONST double wt[],
  double *sw,
  double wmean[],
  double c[],
  int *ifail
);

extern void __stdcall G02BWF(
  CONST int *m,
  double r[],
  int *ifail
);

extern void __stdcall G02BXF(
  CONST char *weight, CONST int weight_len,
  CONST int *n,
  CONST int *m,
  CONST double x[] /* 2 dimension */,
  CONST int *ldx,
  CONST double wt[],
  double xbar[],
  double std[],
  double v[] /* 2 dimension */,
  CONST int *ldv,
  double r[] /* 2 dimension */,
  int *ifail
);

extern void __stdcall G02BYF(
  CONST int *m,
  CONST int *ny,
  CONST int *nx,
  CONST int isz[],
  CONST double r[] /* 2 dimension */,
  CONST int *ldr,
  double p[] /* 2 dimension */,
  CONST int *ldp,
  double wk[],
  int *ifail
);

extern void __stdcall G02CAF(
  CONST int *n,
  CONST double x[],
  CONST double y[],
  double result[],
  int *ifail
);

extern void __stdcall G02CBF(
  CONST int *n,
  CONST double x[],
  CONST double y[],
  double result[],
  int *ifail
);

extern void __stdcall G02CCF(
  CONST int *n,
  CONST double x[],
  CONST double y[],
  CONST double *xmiss,
  CONST double *ymiss,
  double result[],
  int *ifail
);

extern void __stdcall G02CDF(
  CONST int *n,
  CONST double x[],
  CONST double y[],
  CONST double *xmiss,
  CONST double *ymiss,
  double result[],
  int *ifail
);

extern void __stdcall G02CEF(
  CONST int *nvars,
  CONST double xbar[],
  CONST double std[],
  CONST double ssp[] /* 2 dimension */,
  CONST int *issp,
  CONST double r[] /* 2 dimension */,
  CONST int *ir,
  CONST int *nvars2,
  CONST int korder[],
  double xbar2[],
  double std2[],
  double ssp2[] /* 2 dimension */,
  CONST int *issp2,
  double r2[] /* 2 dimension */,
  CONST int *ir2,
  int *ifail
);

extern void __stdcall G02CFF(
  CONST int *nvars,
  CONST int korder[],
  double xbar[],
  double std[],
  double ssp[] /* 2 dimension */,
  CONST int *issp,
  double r[] /* 2 dimension */,
  CONST int *ir,
  int kwork[],
  int *ifail
);

extern void __stdcall G02CGF(
  CONST int *ncases,
  CONST int *nvars,
  CONST int *nind,
  CONST double xbar[],
  CONST double ssp[] /* 2 dimension */,
  CONST int *issp,
  double r[] /* 2 dimension */,
  CONST int *ir,
  double result[],
  double coeff[] /* 2 dimension */,
  CONST int *icoeff,
  double _const[],
  double rinv[] /* 2 dimension */,
  CONST int *irinv,
  double c[] /* 2 dimension */,
  CONST int *ic,
  double wkz[] /* 2 dimension */,
  CONST int *iwkz,
  int *ifail
);

extern void __stdcall G02CHF(
  CONST int *ncases,
  CONST int *nvars,
  CONST int *nind,
  CONST double sspz[] /* 2 dimension */,
  CONST int *isspz,
  double rz[] /* 2 dimension */,
  CONST int *irz,
  double result[],
  double coeff[] /* 2 dimension */,
  CONST int *icoeff,
  double rzinv[] /* 2 dimension */,
  CONST int *irzinv,
  double cz[] /* 2 dimension */,
  CONST int *icz,
  double wkz[] /* 2 dimension */,
  CONST int *iwkz,
  int *ifail
);

extern void __stdcall G02CJF(
  CONST double x[] /* 2 dimension */,
  CONST int *ix,
  CONST double y[] /* 2 dimension */,
  CONST int *iy,
  CONST int *n,
  CONST int *m,
  CONST int *ir,
  double theta[] /* 2 dimension */,
  CONST int *it,
  double sigsq[],
  double c[] /* 2 dimension */,
  CONST int *ic,
  int ipiv[],
  double wk1[] /* 2 dimension */,
  double wk2[],
  int *ifail
);

extern void __stdcall G02DAF(
  CONST char *mean, CONST int mean_len,
  CONST char *weight, CONST int weight_len,
  CONST int *n,
  CONST double x[] /* 2 dimension */,
  CONST int *ldx,
  CONST int *m,
  CONST int isx[],
  CONST int *ip,
  CONST double y[],
  CONST double wt[],
  double *rss,
  int *idf,
  double b[],
  double se[],
  double cov[],
  double res[],
  double h[],
  double q[] /* 2 dimension */,
  CONST int *ldq,
  int *svd,
  int *irank,
  double p[],
  CONST double *tol,
  double wk[],
  int *ifail
);

extern void __stdcall G02DCF(
  CONST char *update, CONST int update_len,
  CONST char *mean, CONST int mean_len,
  CONST char *weight, CONST int weight_len,
  CONST int *m,
  CONST int isx[],
  double q[] /* 2 dimension */,
  CONST int *ldq,
  CONST int *ip,
  CONST double x[],
  CONST int *ix,
  CONST double *y,
  CONST double *wt,
  double *rss,
  double wk[],
  int *ifail
);

extern void __stdcall G02DDF(
  CONST int *n,
  CONST int *ip,
  CONST double q[] /* 2 dimension */,
  CONST int *ldq,
  double *rss,
  int *idf,
  double b[],
  double se[],
  double cov[],
  int *svd,
  int *irank,
  double p[],
  CONST double *tol,
  double wk[],
  int *ifail
);

extern void __stdcall G02DEF(
  CONST char *weight, CONST int weight_len,
  CONST int *n,
  CONST int *ip,
  double q[] /* 2 dimension */,
  CONST int *ldq,
  double p[],
  CONST double wt[],
  CONST double x[],
  double *rss,
  CONST double *tol,
  int *ifail
);

extern void __stdcall G02DFF(
  CONST int *ip,
  double q[] /* 2 dimension */,
  CONST int *ldq,
  CONST int *indx,
  double *rss,
  double wk[],
  int *ifail
);

extern void __stdcall G02DGF(
  CONST char *weight, CONST int weight_len,
  CONST int *n,
  CONST double wt[],
  double *rss,
  CONST int *ip,
  CONST int *irank,
  double cov[],
  double q[] /* 2 dimension */,
  CONST int *ldq,
  CONST int *svd,
  CONST double p[],
  CONST double y[],
  double b[],
  double se[],
  double res[],
  CONST double wk[],
  int *ifail
);

extern void __stdcall G02DKF(
  CONST int *ip,
  CONST int *iconst,
  CONST double p[],
  CONST double c[] /* 2 dimension */,
  CONST int *ldc,
  double b[],
  CONST double *rss,
  CONST int *idf,
  double se[],
  double cov[],
  double wk[],
  int *ifail
);

extern void __stdcall G02DNF(
  CONST int *ip,
  CONST int *irank,
  CONST double b[],
  CONST double cov[],
  CONST double p[],
  CONST double f[],
  int *est,
  double *stat,
  double *sestat,
  double *t,
  CONST double *tol,
  double wk[],
  int *ifail
);

extern void __stdcall G02EAF(
  CONST char *mean, CONST int mean_len,
  CONST char *weight, CONST int weight_len,
  CONST int *n,
  CONST int *m,
  CONST double x[] /* 2 dimension */,
  CONST int *ldx,
  CONST char *name, CONST int name_len,
  CONST int isx[],
  CONST double y[],
  CONST double wt[],
  int *nmod,
  char *model /* 2 dimension */, CONST int model_len,
  CONST int *ldm,
  double rss[],
  int nterms[],
  int mrank[],
  double wk[],
  int *ifail
);

extern void __stdcall G02ECF(
  CONST char *mean, CONST int mean_len,
  CONST int *n,
  CONST double *sigsq,
  CONST double *tss,
  CONST int *nmod,
  CONST int nterms[],
  CONST double rss[],
  double rsq[],
  double cp[],
  int *ifail
);

extern void __stdcall G02EEF(
  int *istep,
  CONST char *mean, CONST int mean_len,
  CONST char *weight, CONST int weight_len,
  CONST int *n,
  CONST int *m,
  CONST double x[] /* 2 dimension */,
  CONST int *ldx,
  CONST char *name, CONST int name_len,
  CONST int isx[],
  CONST int *maxip,
  CONST double y[],
  CONST double wt[],
  CONST double *fin,
  int *addvar,
  char *newvar, CONST int newvar_len,
  double *chrss,
  double *f,
  char *model, CONST int model_len,
  int *nterm,
  double *rss,
  int *idf,
  int *ifr,
  char *free, CONST int free_len,
  double exss[],
  double q[] /* 2 dimension */,
  CONST int *ldq,
  double p[],
  double wk[],
  int *ifail
);

extern void __stdcall G02FAF(
  CONST int *n,
  CONST int *ip,
  CONST int *nres,
  CONST double res[],
  CONST double h[],
  CONST double *rms,
  double sres[] /* 2 dimension */,
  CONST int *lds,
  int *ifail
);

extern void __stdcall G02FCF(
  CONST int *n,
  CONST int *ip,
  CONST double res[],
  double *d,
  double *pdl,
  double *pdu,
  double work[],
  int *ifail
);

extern void __stdcall G02GAF(
  CONST char *link, CONST int link_len,
  CONST char *mean, CONST int mean_len,
  CONST char *offset, CONST int offset_len,
  CONST char *weight, CONST int weight_len,
  CONST int *n,
  CONST double x[] /* 2 dimension */,
  CONST int *ldx,
  CONST int *m,
  CONST int isx[],
  CONST int *ip,
  CONST double y[],
  CONST double wt[],
  double *s,
  double *a,
  double *rss,
  int *idf,
  double b[],
  int *irank,
  double se[],
  double cov[],
  double v[] /* 2 dimension */,
  CONST int *ldv,
  CONST double *tol,
  CONST int *maxit,
  CONST int *iprint,
  CONST double *eps,
  double wk[],
  int *ifail
);

extern void __stdcall G02GBF(
  CONST char *link, CONST int link_len,
  CONST char *mean, CONST int mean_len,
  CONST char *offset, CONST int offset_len,
  CONST char *weight, CONST int weight_len,
  CONST int *n,
  CONST double x[] /* 2 dimension */,
  CONST int *ldx,
  CONST int *m,
  CONST int isx[],
  CONST int *ip,
  CONST double y[],
  double t[],
  CONST double wt[],
  double *dev,
  int *idf,
  double b[],
  int *irank,
  double se[],
  double cov[],
  double v[] /* 2 dimension */,
  CONST int *ldv,
  CONST double *tol,
  CONST int *maxit,
  CONST int *iprint,
  CONST double *eps,
  double wk[],
  int *ifail
);

extern void __stdcall G02GCF(
  CONST char *link, CONST int link_len,
  CONST char *mean, CONST int mean_len,
  CONST char *offset, CONST int offset_len,
  CONST char *weight, CONST int weight_len,
  CONST int *n,
  CONST double x[] /* 2 dimension */,
  CONST int *ldx,
  CONST int *m,
  CONST int isx[],
  CONST int *ip,
  CONST double y[],
  CONST double wt[],
  double *a,
  double *dev,
  int *idf,
  double b[],
  int *irank,
  double se[],
  double cov[],
  double v[] /* 2 dimension */,
  CONST int *ldv,
  CONST double *tol,
  CONST int *maxit,
  CONST int *iprint,
  CONST double *eps,
  double wk[],
  int *ifail
);

extern void __stdcall G02GDF(
  CONST char *link, CONST int link_len,
  CONST char *mean, CONST int mean_len,
  CONST char *offset, CONST int offset_len,
  CONST char *weight, CONST int weight_len,
  CONST int *n,
  CONST double x[] /* 2 dimension */,
  CONST int *ldx,
  CONST int *m,
  CONST int isx[],
  CONST int *ip,
  CONST double y[],
  CONST double wt[],
  double *s,
  double *a,
  double *dev,
  int *idf,
  double b[],
  int *irank,
  double se[],
  double cov[],
  double v[] /* 2 dimension */,
  CONST int *ldv,
  CONST double *tol,
  CONST int *maxit,
  CONST int *iprint,
  CONST double *eps,
  double wk[],
  int *ifail
);

extern void __stdcall G02GKF(
  CONST int *ip,
  CONST int *iconst,
  CONST double v[] /* 2 dimension */,
  CONST int *ldv,
  CONST double c[] /* 2 dimension */,
  CONST int *ldc,
  double b[],
  CONST double *s,
  double se[],
  double cov[],
  double wk[],
  int *ifail
);

extern void __stdcall G02GNF(
  CONST int *ip,
  CONST int *irank,
  CONST double b[],
  CONST double cov[],
  CONST double v[] /* 2 dimension */,
  CONST int *ldv,
  CONST double f[],
  int *est,
  double *stat,
  double *sestat,
  double *z,
  CONST double *tol,
  double wk[],
  int *ifail
);

extern void __stdcall G02HAF(
  CONST int *indw,
  CONST int *ipsi,
  CONST int *isigma,
  CONST int *indc,
  CONST int *n,
  CONST int *m,
  double x[] /* 2 dimension */,
  CONST int *ix,
  double y[],
  CONST double *cpsi,
  CONST double *h1,
  CONST double *h2,
  CONST double *h3,
  CONST double *cucv,
  CONST double *dchi,
  double theta[],
  double *sigma,
  double c[] /* 2 dimension */,
  CONST int *ic,
  double rs[],
  double wgt[],
  CONST double *tol,
  CONST int *maxit,
  CONST int *nitmon,
  double work[],
  int *ifail
);

extern void __stdcall G02HBF(
  double (__stdcall *ucv)(double *),
  CONST int *n,
  CONST int *m,
  CONST double x[] /* 2 dimension */,
  CONST int *ix,
  double a[],
  double z[],
  CONST double *bl,
  CONST double *bd,
  CONST double *tol,
  CONST int *maxit,
  CONST int *nitmon,
  int *nit,
  double wk[],
  int *ifail
);

extern void __stdcall G02HDF(
  double (__stdcall *chi)(double *),
  double (__stdcall *psi)(double *),
  CONST double *psip0,
  CONST double *beta,
  CONST int *indw,
  CONST int *isigma,
  CONST int *n,
  CONST int *m,
  double x[] /* 2 dimension */,
  CONST int *ix,
  double y[],
  double wgt[],
  double theta[],
  int *k,
  double *sigma,
  double rs[],
  CONST double *tol,
  CONST double *eps,
  CONST int *maxit,
  CONST int *nitmon,
  int *nit,
  double wk[],
  int *ifail
);

extern double __stdcall G02HDZ(
  CONST double *t
);

extern void __stdcall G02HFF(
  double (__stdcall *psi)(double *),
  double (__stdcall *psp)(double *),
  CONST int *indw,
  CONST int *indc,
  CONST double *sigma,
  CONST int *n,
  CONST int *m,
  CONST double x[] /* 2 dimension */,
  CONST int *ix,
  CONST double rs[],
  CONST double wgt[],
  double c[] /* 2 dimension */,
  CONST int *ic,
  double wk[],
  int *ifail
);

extern void __stdcall G02HKF(
  CONST int *n,
  CONST int *m,
  CONST double x[] /* 2 dimension */,
  CONST int *ldx,
  CONST double *eps,
  double cov[],
  double theta[],
  CONST int *maxit,
  CONST int *nitmon,
  CONST double *tol,
  int *nit,
  double wk[],
  int *ifail
);

extern void __stdcall G02HLF(
  void (__stdcall *ucv)(double *, double[], double *, double *, double *, double *),
  CONST double userp[],
  CONST int *indm,
  CONST int *n,
  CONST int *m,
  CONST double x[] /* 2 dimension */,
  CONST int *ldx,
  double cov[],
  double a[],
  double wt[],
  double theta[],
  CONST double *bl,
  CONST double *bd,
  CONST int *maxit,
  CONST int *nitmon,
  CONST double *tol,
  int *nit,
  double wk[],
  int *ifail
);

extern void __stdcall G02HMF(
  void (__stdcall *ucv)(double *, double[], double *, double *),
  CONST double userp[],
  CONST int *indm,
  CONST int *n,
  CONST int *m,
  CONST double x[] /* 2 dimension */,
  CONST int *ldx,
  double cov[],
  double a[],
  double wt[],
  double theta[],
  CONST double *bl,
  CONST double *bd,
  CONST int *maxit,
  CONST int *nitmon,
  CONST double *tol,
  int *nit,
  double wk[],
  int *ifail
);

extern void __stdcall G03AAF(
  CONST char *matrix, CONST int matrix_len,
  CONST char *std, CONST int std_len,
  CONST char *weight, CONST int weight_len,
  CONST int *n,
  CONST int *m,
  CONST double x[] /* 2 dimension */,
  CONST int *ldx,
  CONST int isx[],
  double s[],
  CONST double wt[],
  CONST int *nvar,
  double e[] /* 2 dimension */,
  CONST int *lde,
  double p[] /* 2 dimension */,
  CONST int *ldp,
  double v[] /* 2 dimension */,
  CONST int *ldv,
  double wk[],
  int *ifail
);

extern void __stdcall G03ACF(
  CONST char *weight, CONST int weight_len,
  CONST int *n,
  CONST int *m,
  CONST double x[] /* 2 dimension */,
  CONST int *ldx,
  CONST int isx[],
  CONST int *nx,
  CONST int ing[],
  CONST int *ng,
  CONST double wt[],
  int nig[],
  double cvm[] /* 2 dimension */,
  CONST int *ldcvm,
  double e[] /* 2 dimension */,
  CONST int *lde,
  int *ncv,
  double cvx[] /* 2 dimension */,
  CONST int *ldcvx,
  CONST double *tol,
  int *irankx,
  double wk[],
  CONST int *iwk,
  int *ifail
);

extern void __stdcall G03ADF(
  CONST char *weight, CONST int weight_len,
  CONST int *n,
  CONST int *m,
  CONST double z[] /* 2 dimension */,
  CONST int *ldz,
  CONST int isz[],
  int *nx,
  int *ny,
  CONST double wt[],
  double e[] /* 2 dimension */,
  CONST int *lde,
  int *ncv,
  double cvx[] /* 2 dimension */,
  CONST int *ldcvx,
  CONST int *mcv,
  double cvy[] /* 2 dimension */,
  CONST int *ldcvy,
  CONST double *tol,
  double wk[],
  CONST int *iwk,
  int *ifail
);

extern void __stdcall G03BAF(
  CONST char *stand, CONST int stand_len,
  CONST double *g,
  CONST int *nvar,
  CONST int *k,
  double fl[] /* 2 dimension */,
  CONST int *ldf,
  double flr[] /* 2 dimension */,
  double r[] /* 2 dimension */,
  CONST int *ldr,
  CONST double *acc,
  CONST int *maxit,
  int *iter,
  double wk[],
  int *ifail
);

extern void __stdcall G03BCF(
  CONST char *stand, CONST int stand_len,
  CONST char *pscale, CONST int pscale_len,
  CONST int *n,
  CONST int *m,
  double x[] /* 2 dimension */,
  CONST int *ldx,
  double y[] /* 2 dimension */,
  CONST int *ldy,
  double yhat[] /* 2 dimension */,
  double r[] /* 2 dimension */,
  CONST int *ldr,
  double *alpha,
  double *rss,
  double res[],
  double wk[],
  int *ifail
);

extern void __stdcall G03CAF(
  CONST char *matrix, CONST int matrix_len,
  CONST char *weight, CONST int weight_len,
  CONST int *n,
  CONST int *m,
  CONST double x[] /* 2 dimension */,
  CONST int *ldx,
  CONST int *nvar,
  CONST int isx[],
  CONST int *nfac,
  CONST double wt[],
  double e[],
  double stat[],
  double com[],
  double psi[],
  double res[],
  double fl[] /* 2 dimension */,
  CONST int *ldfl,
  CONST int iop[],
  int iwk[],
  double wk[],
  CONST int *lwk,
  int *ifail
);

extern void __stdcall G03CCF(
  CONST char *method, CONST int method_len,
  CONST char *rotate, CONST int rotate_len,
  CONST int *nvar,
  CONST int *nfac,
  CONST double fl[] /* 2 dimension */,
  CONST int *ldfl,
  CONST double psi[],
  CONST double e[],
  CONST double r[] /* 2 dimension */,
  CONST int *ldr,
  double fs[] /* 2 dimension */,
  CONST int *ldfs,
  double wk[],
  int *ifail
);

extern void __stdcall G03DAF(
  CONST char *weight, CONST int weight_len,
  CONST int *n,
  CONST int *m,
  CONST double x[] /* 2 dimension */,
  CONST int *ldx,
  CONST int isx[],
  CONST int *nvar,
  CONST int ing[],
  CONST int *ng,
  CONST double wt[],
  int nig[],
  double gmean[] /* 2 dimension */,
  CONST int *ldg,
  double det[],
  double gc[],
  double *stat,
  double *df,
  double *sig,
  double wk[],
  int iwk[],
  int *ifail
);

extern void __stdcall G03DBF(
  CONST char *equal, CONST int equal_len,
  CONST char *mode, CONST int mode_len,
  CONST int *nvar,
  CONST int *ng,
  CONST double gmean[] /* 2 dimension */,
  CONST int *ldg,
  CONST double gc[],
  CONST int *nobs,
  CONST int *m,
  CONST int isx[],
  CONST double x[] /* 2 dimension */,
  CONST int *ldx,
  double d[] /* 2 dimension */,
  CONST int *ldd,
  double wk[],
  int *ifail
);

extern void __stdcall G03DCF(
  CONST char *type, CONST int type_len,
  CONST char *equal, CONST int equal_len,
  CONST char *priors, CONST int priors_len,
  CONST int *nvar,
  CONST int *ng,
  CONST int nig[],
  CONST double gmean[] /* 2 dimension */,
  CONST int *ldg,
  CONST double gc[],
  CONST double det[],
  CONST int *nobs,
  CONST int *m,
  CONST int isx[],
  CONST double x[] /* 2 dimension */,
  CONST int *ldx,
  double prior[],
  double p[] /* 2 dimension */,
  CONST int *ldp,
  int iag[],
  CONST int *atiq,
  double ati[] /* 2 dimension */,
  double wk[],
  int *ifail
);

extern void __stdcall G03EAF(
  CONST char *update, CONST int update_len,
  CONST char *dist, CONST int dist_len,
  CONST char *scale, CONST int scale_len,
  CONST int *n,
  CONST int *m,
  CONST double x[] /* 2 dimension */,
  CONST int *ldx,
  CONST int isx[],
  double s[],
  double d[],
  int *ifail
);

extern void __stdcall G03ECF(
  CONST int *method,
  CONST int *n,
  double d[],
  int ilc[],
  int iuc[],
  double cd[],
  int iord[],
  double dord[],
  int iwk[],
  int *ifail
);

extern void __stdcall G03EFF(
  CONST char *weight, CONST int weight_len,
  CONST int *n,
  CONST int *m,
  CONST double x[] /* 2 dimension */,
  CONST int *ldx,
  CONST int isx[],
  CONST int *nvar,
  CONST int *k,
  double cmeans[] /* 2 dimension */,
  CONST int *ldc,
  CONST double wt[],
  int inc[],
  int nic[],
  double css[],
  double csw[],
  CONST int *maxit,
  int iwk[],
  double wk[],
  int *ifail
);

extern void __stdcall G03EHF(
  CONST char *orient, CONST int orient_len,
  CONST int *n,
  CONST double dord[],
  CONST double *dmin,
  CONST double *dstep,
  CONST int *nsym,
  char *c, CONST int c_len,
  CONST int *lenc,
  int *ifail
);

extern void __stdcall G03EJF(
  CONST int *n,
  CONST double cd[],
  CONST int iord[],
  CONST double dord[],
  int *k,
  double *dlevel,
  int ic[],
  int *ifail
);

extern void __stdcall G03FAF(
  CONST char *roots, CONST int roots_len,
  CONST int *n,
  CONST double d[],
  CONST int *ndim,
  double x[] /* 2 dimension */,
  CONST int *ldx,
  double eval[],
  double wk[],
  int iwk[],
  int *ifail
);

extern void __stdcall G03FCF(
  CONST char *type, CONST int type_len,
  CONST int *n,
  CONST int *ndim,
  CONST double d[],
  double x[] /* 2 dimension */,
  CONST int *ldx,
  double *stress,
  double dfit[],
  CONST int *iter,
  CONST int *iopt,
  double wk[],
  int iwk[],
  int *ifail
);

extern void __stdcall G03ZAF(
  CONST int *n,
  CONST int *m,
  CONST double x[] /* 2 dimension */,
  CONST int *ldx,
  CONST int *nvar,
  CONST int isx[],
  CONST double s[],
  CONST double e[],
  double z[] /* 2 dimension */,
  CONST int *ldz,
  int *ifail
);

extern void __stdcall G04ADF(
  CONST double data[],
  double var[],
  double amr[],
  double amc[],
  double amt[],
  CONST int lcode[] /* 2 dimension */,
  CONST int *ia,
  CONST int *n,
  CONST int *nn
);

extern void __stdcall G04AEF(
  CONST double y[],
  CONST int *n,
  CONST int *k,
  CONST int nobs[],
  double gbar[],
  double *gm,
  double ss[],
  int idf[],
  double *f,
  double *fp,
  int *ifail
);

extern void __stdcall G04AFF(
  CONST double y[] /* 3 dimension */,
  CONST int *iy1,
  CONST int *iy2,
  CONST int *m,
  CONST int *nr,
  CONST int *nc,
  double row[],
  double col[],
  double cell[] /* 2 dimension */,
  CONST int *icell,
  double *gm,
  double ss[],
  int idf[],
  double f[],
  double fp[],
  int *ifail
);

extern void __stdcall G04AGF(
  CONST double y[],
  CONST int *n,
  CONST int *k,
  CONST int lsub[],
  CONST int nobs[],
  CONST int *l,
  int ngp[],
  double gbar[],
  double sgbar[],
  double *gm,
  double ss[],
  int idf[],
  double f[],
  double fp[],
  int *ifail
);

extern void __stdcall G04BBF(
  CONST int *n,
  CONST double y[],
  CONST int *iblock,
  CONST int *nt,
  CONST int it[],
  double *gmean,
  double bmean[],
  double tmean[],
  double table[] /* 2 dimension */,
  CONST int *ldt,
  double c[] /* 2 dimension */,
  CONST int *ldc,
  int irep[],
  double r[],
  double ef[],
  CONST double *tol,
  CONST int *irdf,
  double wk[],
  int *ifail
);

extern void __stdcall G04BCF(
  CONST int *nrep,
  CONST int *nrow,
  CONST int *ncol,
  CONST double y[],
  CONST int *nt,
  CONST int it[],
  double *gmean,
  double tmean[],
  double table[] /* 2 dimension */,
  CONST int *ldt,
  double c[] /* 2 dimension */,
  CONST int *ldc,
  int irep[],
  double rpmean[],
  double rmean[],
  double cmean[],
  double r[],
  double ef[],
  CONST double *tol,
  CONST int *irdf,
  double wk[],
  int *ifail
);

extern void __stdcall G04CAF(
  CONST int *n,
  CONST double y[],
  CONST int *nfac,
  CONST int lfac[],
  CONST int *nblock,
  CONST int *inter,
  CONST int *irdf,
  CONST int *mterm,
  double table[] /* 2 dimension */,
  int *itotal,
  double tmean[],
  CONST int *maxt,
  double e[],
  int imean[],
  double semean[],
  double bmean[],
  double r[],
  int iwk[],
  int *ifail
);

extern void __stdcall G04DAF(
  CONST int *nt,
  CONST double tmean[],
  CONST int irep[],
  CONST double *rms,
  CONST double *rdf,
  CONST int *nc,
  CONST double ct[] /* 2 dimension */,
  CONST int *ldct,
  double est[],
  double table[] /* 2 dimension */,
  CONST int *ldt,
  CONST double *tol,
  CONST int *usetx,
  CONST double tx[],
  int *ifail
);

extern void __stdcall G04DBF(
  CONST char *type, CONST int type_len,
  CONST int *nt,
  CONST double tmean[],
  CONST double *rdf,
  CONST double c[] /* 2 dimension */,
  CONST int *ldc,
  CONST double *clevel,
  double cil[],
  double ciu[],
  int isig[],
  int *ifail
);

extern void __stdcall G04EAF(
  CONST char *type, CONST int type_len,
  CONST int *n,
  CONST int *levels,
  CONST int ifact[],
  double x[] /* 2 dimension */,
  CONST int *ldx,
  CONST double v[],
  double rep[],
  int *ifail
);

extern double __stdcall G05CAF(
  CONST double *x
);

extern void __stdcall G05CBF(
  CONST int *i
);

extern void __stdcall G05CCF(
void
);

extern void __stdcall G05CFF(
  int ia[],
  CONST int *ni,
  double xa[],
  CONST int *nx,
  int *ifail
);

extern void __stdcall G05CGF(
  CONST int ia[],
  CONST int *ni,
  CONST double xa[],
  CONST int *nx,
  int *ifail
);

extern double __stdcall G05DAF(
  CONST double *a,
  CONST double *b
);

extern double __stdcall G05DBF(
  CONST double *a
);

extern double __stdcall G05DCF(
  CONST double *a,
  CONST double *b
);

extern double __stdcall G05DDF(
  CONST double *a,
  CONST double *b
);

extern double __stdcall G05DEF(
  CONST double *a,
  CONST double *b
);

extern double __stdcall G05DFF(
  CONST double *a,
  CONST double *b
);

extern double __stdcall G05DGF(
  CONST double *a,
  CONST double *b,
  int *ifail
);

extern double __stdcall G05DHF(
  CONST int *n,
  int *ifail
);

extern double __stdcall G05DJF(
  CONST int *n,
  int *ifail
);

extern double __stdcall G05DKF(
  CONST int *m,
  CONST int *n,
  int *ifail
);

extern double __stdcall G05DLF(
  CONST double *g,
  CONST double *h,
  int *ifail
);

extern double __stdcall G05DMF(
  CONST double *g,
  CONST double *h,
  int *ifail
);

extern double __stdcall G05DPF(
  CONST double *a,
  CONST double *b,
  int *ifail
);

extern int __stdcall G05DRF(
  CONST double *alamda,
  int *ifail
);

extern int __stdcall G05DYF(
  CONST int *m,
  CONST int *n
);

extern int __stdcall G05DZF(
  CONST double *p
);

extern void __stdcall G05EAF(
  CONST double a[],
  CONST int *n,
  CONST double c[] /* 2 dimension */,
  CONST int *ic,
  CONST double *eps,
  double r[],
  CONST int *nr,
  int *ifail
);

extern void __stdcall G05EBF(
  CONST int *m,
  CONST int *n,
  double r[],
  CONST int *nr,
  int *ifail
);

extern void __stdcall G05ECF(
  CONST double *t,
  double r[],
  CONST int *nr,
  int *ifail
);

extern void __stdcall G05EDF(
  CONST int *n,
  CONST double *p,
  double r[],
  CONST int *nr,
  int *ifail
);

extern void __stdcall G05EEF(
  CONST int *n,
  CONST double *p,
  double r[],
  CONST int *nr,
  int *ifail
);

extern void __stdcall G05EFF(
  CONST int *l,
  CONST int *m,
  CONST int *n,
  double r[],
  CONST int *nr,
  int *ifail
);

extern void __stdcall G05EGF(
  CONST double *e,
  CONST double a[],
  CONST int *na,
  CONST double b[],
  CONST int *nb,
  double r[],
  CONST int *nr,
  double *var,
  int *ifail
);

extern void __stdcall G05EHF(
  int index[],
  CONST int *n,
  int *ifail
);

extern void __stdcall G05EJF(
  CONST int ia[],
  CONST int *n,
  int iz[],
  CONST int *m,
  int *ifail
);

extern double __stdcall G05EWF(
  double r[],
  CONST int *nr,
  int *ifail
);

extern void __stdcall G05EXF(
  CONST double p[],
  CONST int *np,
  CONST int *ip,
  CONST int *lp,
  double r[],
  CONST int *nr,
  int *ifail
);

extern int __stdcall G05EYF(
  CONST double r[],
  CONST int *nr
);

extern void __stdcall G05EZF(
  double z[],
  CONST int *n,
  CONST double r[],
  CONST int *nr,
  int *ifail
);

extern void __stdcall G05FAF(
  CONST double *a,
  CONST double *b,
  CONST int *n,
  double x[]
);

extern void __stdcall G05FBF(
  CONST double *a,
  CONST int *n,
  double x[]
);

extern void __stdcall G05FDF(
  CONST double *a,
  CONST double *b,
  CONST int *n,
  double x[]
);

extern void __stdcall G05FEF(
  CONST double *a,
  CONST double *b,
  CONST int *n,
  double x[],
  int *ifail
);

extern void __stdcall G05FFF(
  CONST double *a,
  CONST double *b,
  CONST int *n,
  double x[],
  int *ifail
);

extern void __stdcall G05FSF(
  CONST double *vk,
  CONST int *n,
  double t[],
  int *ifail
);

extern void __stdcall G05GAF(
  CONST char *side, CONST int side_len,
  CONST char *init, CONST int init_len,
  CONST int *m,
  CONST int *n,
  double a[] /* 2 dimension */,
  CONST int *lda,
  double wk[],
  int *ifail
);

extern void __stdcall G05GBF(
  CONST int *n,
  CONST double d[],
  double c[] /* 2 dimension */,
  CONST int *ldc,
  CONST double *eps,
  double wk[],
  int *ifail
);

extern void __stdcall G05HDF(
  CONST char *mode, CONST int mode_len,
  CONST int *k,
  CONST int *ip,
  CONST int *iq,
  CONST char *mean, CONST int mean_len,
  CONST double par[],
  CONST int *lpar,
  double qq[] /* 2 dimension */,
  CONST int *ik,
  CONST int *n,
  double w[] /* 2 dimension */,
  double ref[],
  CONST int *lref,
  int iwork[],
  CONST int *liwork,
  int *ifail
);

extern void __stdcall G07AAF(
  CONST int *n,
  CONST int *k,
  CONST double *clevel,
  double *pl,
  double *pu,
  int *ifail
);

extern void __stdcall G07ABF(
  CONST int *n,
  CONST double *xmean,
  CONST double *clevel,
  double *tl,
  double *tu,
  int *ifail
);

extern void __stdcall G07BBF(
  CONST char *method, CONST int method_len,
  CONST int *n,
  CONST double x[],
  CONST double xc[],
  CONST int ic[],
  double *xmu,
  double *xsig,
  CONST double *tol,
  CONST int *maxit,
  double *sexmu,
  double *sexsig,
  double *corr,
  double *dev,
  int nobs[],
  int *nit,
  double wk[],
  int *ifail
);

extern void __stdcall G07BEF(
  CONST char *cens, CONST int cens_len,
  CONST int *n,
  CONST double x[],
  CONST int ic[],
  double *beta,
  double *gamma,
  CONST double *tol,
  CONST int *maxit,
  double *sebeta,
  double *segam,
  double *corr,
  double *dev,
  int *nit,
  double wk[],
  int *ifail
);

extern void __stdcall G07CAF(
  CONST char *tail, CONST int tail_len,
  CONST char *equal, CONST int equal_len,
  CONST int *nx,
  CONST int *ny,
  CONST double *xmean,
  CONST double *ymean,
  CONST double *xstd,
  CONST double *ystd,
  CONST double *clevel,
  double *t,
  double *df,
  double *prob,
  double *dl,
  double *du,
  int *ifail
);

extern void __stdcall G07DAF(
  CONST int *n,
  CONST double x[],
  double y[],
  double *xme,
  double *xmd,
  double *xsd,
  int *ifail
);

extern void __stdcall G07DBF(
  CONST int *isigma,
  CONST int *n,
  CONST double x[],
  CONST int *ipsi,
  CONST double *c,
  CONST double *h1,
  CONST double *h2,
  CONST double *h3,
  CONST double *dchi,
  double *theta,
  double *sigma,
  CONST int *maxit,
  CONST double *tol,
  double rs[],
  int *nit,
  double wrk[],
  int *ifail
);

extern void __stdcall G07DCF(
  double (__stdcall *chi)(double *),
  double (__stdcall *psi)(double *),
  CONST int *isigma,
  CONST int *n,
  CONST double x[],
  CONST double *beta,
  double *theta,
  double *sigma,
  CONST int *maxit,
  CONST double *tol,
  double rs[],
  int *nit,
  double wrk[],
  int *ifail
);

extern void __stdcall G07DDF(
  CONST int *n,
  CONST double x[],
  CONST double *alpha,
  double *tmean,
  double *wmean,
  double *tvar,
  double *wvar,
  int *k,
  double sx[],
  int *ifail
);

extern void __stdcall G07EAF(
  CONST char *method, CONST int method_len,
  CONST int *n,
  CONST double x[],
  CONST double *clevel,
  double *theta,
  double *thetal,
  double *thetau,
  double *estcl,
  double *wlower,
  double *wupper,
  double wrk[],
  int iwrk[],
  int *ifail
);

extern void __stdcall G07EBF(
  CONST char *method, CONST int method_len,
  CONST int *n,
  CONST double x[],
  CONST int *m,
  CONST double y[],
  CONST double *clevel,
  double *theta,
  double *thetal,
  double *thetau,
  double *estcl,
  double *ulower,
  double *uupper,
  double wrk[],
  int iwrk[],
  int *ifail
);

extern void __stdcall G08AAF(
  CONST double x[],
  CONST double y[],
  CONST int *n,
  int *is,
  int *n1,
  double *p,
  int *ifail
);

extern void __stdcall G08ABF(
  CONST double x[],
  CONST double y[],
  CONST int *n,
  double w1[],
  double w2[],
  double *w,
  int *n1,
  double *p,
  int *ifail
);

extern void __stdcall G08ACF(
  CONST double x[],
  CONST int *n,
  CONST int *n1,
  double w1[],
  int *i1,
  int *i2,
  double *p,
  int *ifail
);

extern void __stdcall G08ADF(
  CONST double x[],
  CONST int *n,
  CONST int *n1,
  double w1[],
  double *u,
  double *p,
  int *ifail
);

extern void __stdcall G08AEF(
  CONST double x[] /* 2 dimension */,
  CONST int *ix,
  CONST int *k,
  CONST int *n,
  double w1[],
  double w2[],
  double *fr,
  double *p,
  int *ifail
);

extern void __stdcall G08AFF(
  CONST double x[],
  CONST int *lx,
  CONST int l[],
  CONST int *k,
  double w1[],
  double *h,
  double *p,
  int *ifail
);

extern void __stdcall G08AGF(
  CONST int *n,
  CONST double x[],
  CONST double *xme,
  CONST char *tail, CONST int tail_len,
  CONST char *zeros, CONST int zeros_len,
  double *rs,
  double *rsnor,
  double *p,
  int *nz1,
  double wrk[],
  int *ifail
);

extern void __stdcall G08AHF(
  CONST int *n1,
  CONST double x[],
  CONST int *n2,
  CONST double y[],
  CONST char *tail, CONST int tail_len,
  double *u,
  double *unor,
  double *p,
  int *ties,
  double ranks[],
  double wrk[],
  int *ifail
);

extern void __stdcall G08AJF(
  CONST int *n1,
  CONST int *n2,
  CONST char *tail, CONST int tail_len,
  CONST double *u,
  double *p,
  double wrk[],
  CONST int *lwrk,
  int *ifail
);

extern void __stdcall G08AKF(
  CONST int *n1,
  CONST int *n2,
  CONST char *tail, CONST int tail_len,
  CONST double ranks[],
  CONST double *u,
  double *p,
  double wrk[],
  CONST int *lwrk,
  int iwrk[],
  int *ifail
);

extern void __stdcall G08ALF(
  CONST int *n,
  CONST int *k,
  CONST double x[] /* 2 dimension */,
  CONST int *ldx,
  double *q,
  double *prob,
  int *ifail
);

extern void __stdcall G08BAF(
  CONST double x[],
  CONST int *n,
  CONST int *n1,
  double r[],
  CONST int *itest,
  double *w,
  double *v,
  double *pw,
  double *pv,
  int *ifail
);

extern void __stdcall G08CAF(
  CONST int *n,
  CONST double x[],
  CONST int *null,
  CONST int *np,
  double p[],
  CONST int *nest,
  CONST int *ntype,
  double *d,
  double *prob,
  double s[],
  int ind[],
  int *ifail
);

extern void __stdcall G08CBF(
  CONST int *n,
  CONST double x[],
  CONST char *dist, CONST int dist_len,
  double par[],
  CONST char *estima, CONST int estima_len,
  CONST int *ntype,
  double *d,
  double *z,
  double *p,
  double sx[],
  int *ifail
);

extern void __stdcall G08CCF(
  CONST int *n,
  CONST double x[],
  double (__stdcall *cdf)(double[]),
  CONST int *ntype,
  double *d,
  double *z,
  double *p,
  double sx[],
  int *ifail
);

extern void __stdcall G08CDF(
  CONST int *n1,
  CONST double x[],
  CONST int *n2,
  CONST double y[],
  CONST int *ntype,
  double *d,
  double *z,
  double *p,
  double sx[],
  double sy[],
  int *ifail
);

extern void __stdcall G08CGF(
  CONST int *k2,
  CONST int ifreq[],
  CONST double cint[],
  CONST char *dist, CONST int dist_len,
  CONST double par[],
  CONST int *iparam,
  CONST double prob[],
  double *chisq,
  double *p,
  int *ndf,
  double eval[],
  double chisqi[],
  int *ifail
);

extern void __stdcall G08DAF(
  CONST double x[] /* 2 dimension */,
  CONST int *ix,
  CONST int *k,
  CONST int *n,
  double rnk[] /* 2 dimension */,
  double *w,
  double *p,
  int *ifail
);

extern void __stdcall G08EAF(
  CONST char *cl, CONST int cl_len,
  CONST int *n,
  CONST double x[],
  CONST int *m,
  CONST int *maxr,
  int *nruns,
  int ncount[],
  double ex[],
  double c[] /* 2 dimension */,
  CONST int *ldc,
  double *chi,
  double *df,
  double *prob,
  double wrk[],
  CONST int *lwrk,
  int *ifail
);

extern void __stdcall G08EBF(
  CONST char *cl, CONST int cl_len,
  CONST int *n,
  CONST double x[],
  CONST int *msize,
  CONST int *lag,
  int ncount[] /* 2 dimension */,
  CONST int *ldc,
  double *ex,
  double *chi,
  double *df,
  double *p,
  double wrk[],
  int *ifail
);

extern void __stdcall G08ECF(
  CONST char *cl, CONST int cl_len,
  CONST int *n,
  CONST double x[],
  CONST int *msize,
  int ncount[] /* 3 dimension */,
  CONST int *ldc,
  double *ex,
  double *chi,
  double *df,
  double *p,
  int *ifail
);

extern void __stdcall G08EDF(
  CONST char *cl, CONST int cl_len,
  CONST int *n,
  CONST double x[],
  CONST int *m,
  CONST int *k,
  CONST double *rl,
  CONST double *ru,
  CONST double *til,
  int *ngaps,
  int ncount[],
  double ex[],
  double *chi,
  double *df,
  double *prob,
  int *ifail
);

extern void __stdcall G08RAF(
  CONST int *ns,
  int nv[],
  CONST int *nsum,
  CONST double y[],
  CONST int *ip,
  CONST double x[] /* 2 dimension */,
  CONST int *nx,
  CONST int *idist,
  CONST int *nmax,
  CONST double *tol,
  double parvar[] /* 2 dimension */,
  CONST int *npvar,
  int irank[],
  double zin[],
  double eta[],
  double vapvec[],
  double parest[],
  double work[],
  CONST int *lwork,
  int iwa[],
  int *ifail
);

extern void __stdcall G08RBF(
  CONST int *ns,
  CONST int nv[],
  CONST int *nsum,
  CONST double y[],
  CONST int *ip,
  CONST double x[] /* 2 dimension */,
  CONST int *nx,
  CONST int icen[],
  CONST double *gamma,
  CONST int *nmax,
  CONST double *tol,
  double parvar[] /* 2 dimension */,
  CONST int *npvar,
  int irank[],
  double zin[],
  double eta[],
  double vapvec[],
  double parest[],
  double work[],
  CONST int *lwork,
  int iwa[],
  int *ifail
);

extern void __stdcall G10ABF(
  CONST char *mode, CONST int mode_len,
  CONST char *weight, CONST int weight_len,
  CONST int *n,
  CONST double x[],
  CONST double y[],
  CONST double wt[],
  CONST double *rho,
  double yhat[],
  double c[] /* 2 dimension */,
  CONST int *ldc,
  double *rss,
  double *df,
  double res[],
  double h[],
  double wk[],
  int *ifail
);

extern void __stdcall G10ACF(
  CONST char *method, CONST int method_len,
  CONST char *weight, CONST int weight_len,
  CONST int *n,
  CONST double x[],
  CONST double y[],
  CONST double wt[],
  double yhat[],
  double c[] /* 2 dimension */,
  CONST int *ldc,
  double *rss,
  double *rdf,
  double res[],
  double h[],
  double *crit,
  double *rho,
  CONST double *u,
  CONST double *tol,
  int *maxcal,
  double wk[],
  int *ifail
);

extern void __stdcall G10BAF(
  CONST int *n,
  CONST double x[],
  CONST double *window,
  CONST double *slo,
  CONST double *shi,
  CONST int *ns,
  double smooth[],
  double t[],
  CONST int *usefft,
  double fft[],
  int *ifail
);

extern void __stdcall G10CAF(
  CONST int *itype,
  CONST int *n,
  CONST double y[],
  double smooth[],
  double rough[],
  int *ifail
);

extern void __stdcall G10ZAF(
  CONST char *weight, CONST int weight_len,
  CONST int *n,
  CONST double x[],
  CONST double y[],
  CONST double wt[],
  int *nord,
  double xord[],
  double yord[],
  double wwt[],
  double *rss,
  int iwrk[],
  int *ifail
);

extern void __stdcall G11AAF(
  CONST int *nrow,
  CONST int *ncol,
  CONST int nobst[] /* 2 dimension */,
  CONST int *ldt,
  double expt[] /* 2 dimension */,
  double chist[] /* 2 dimension */,
  double *prob,
  double *chi,
  double *g,
  double *df,
  int *ifail
);

extern void __stdcall G11BAF(
  CONST char *stat, CONST int stat_len,
  CONST char *update, CONST int update_len,
  CONST char *weight, CONST int weight_len,
  CONST int *n,
  CONST int *nfac,
  CONST int isf[],
  CONST int lfac[],
  CONST int ifac[] /* 2 dimension */,
  CONST int *ldf,
  CONST double y[],
  CONST double wt[],
  double table[],
  CONST int *maxt,
  int *ncells,
  int *ndim,
  int idim[],
  int icount[],
  double auxt[],
  int iwk[],
  int *ifail
);

extern void __stdcall G11BBF(
  CONST char *type, CONST int type_len,
  CONST char *weight, CONST int weight_len,
  CONST int *n,
  CONST int *nfac,
  CONST int isf[],
  CONST int lfac[],
  CONST int ifac[] /* 2 dimension */,
  CONST int *ldf,
  CONST double *percnt,
  CONST double y[],
  CONST double wt[],
  double table[],
  CONST int *maxt,
  int *ncells,
  int *ndim,
  int idim[],
  int icount[],
  int iwk[],
  double wk[],
  int *ifail
);

extern void __stdcall G11BCF(
  CONST char *stat, CONST int stat_len,
  CONST double table[],
  CONST int *ncells,
  CONST int *ndim,
  CONST int idim[],
  CONST int isdim[],
  double stable[],
  CONST int *maxst,
  int *mcells,
  int *mdim,
  int mlevel[],
  double auxt[],
  int iwk[],
  double wk[],
  int *ifail
);

extern void __stdcall G11SAF(
  CONST int *n2,
  CONST int *n,
  CONST int *gprob,
  CONST int *s,
  int x[] /* 2 dimension */,
  CONST int *nrowxr,
  int rl[],
  double a[],
  double c[],
  CONST int *iprint,
  CONST double *cgetol,
  CONST int *maxit,
  CONST int *chisqr,
  CONST int *ishow,
  int *niter,
  double alpha[],
  double gamma[],
  double var[] /* 2 dimension */,
  CONST int *iaa,
  double g[],
  double expp[] /* 2 dimension */,
  CONST int *ia,
  double obs[] /* 2 dimension */,
  double p[],
  double y[],
  double xl[],
  int ob[],
  double *ll,
  double *chi,
  int *idf,
  double *siglev,
  double w[],
  CONST int *lw,
  int *ifail
);

extern void __stdcall G11SBF(
  CONST int *n2,
  CONST int *n,
  int *s,
  int x[] /* 2 dimension */,
  CONST int *nrx,
  int rl[],
  int *ifail
);

extern void __stdcall G12AAF(
  CONST int *n,
  CONST double t[],
  CONST int ic[],
  CONST char *freq, CONST int freq_len,
  CONST int ifreq[],
  int *nd,
  double tp[],
  double p[],
  double psig[],
  int iwk[],
  int *ifail
);

extern void __stdcall G12BAF(
  CONST char *offset, CONST int offset_len,
  CONST int *n,
  CONST int *m,
  CONST int *ns,
  CONST double z[] /* 2 dimension */,
  CONST int *ldz,
  CONST int isz[],
  CONST int *ip,
  CONST double t[],
  CONST int ic[],
  CONST double omega[],
  CONST int isi[],
  double *dev,
  double b[],
  double se[],
  double sc[],
  double cov[],
  double res[],
  int *nd,
  double tp[],
  double sur[] /* 2 dimension */,
  CONST int *ndmax,
  CONST double *tol,
  CONST int *maxit,
  CONST int *iprint,
  double wk[],
  int iwk[],
  int *ifail
);

extern void __stdcall G13AAF(
  CONST double x[],
  CONST int *nx,
  CONST int *nd,
  CONST int *nds,
  CONST int *ns,
  double xd[],
  int *nxd,
  int *ifail
);

extern void __stdcall G13ABF(
  CONST double x[],
  CONST int *nx,
  CONST int *nk,
  double *xm,
  double *xv,
  double r[],
  double *stat,
  int *ifail
);

extern void __stdcall G13ACF(
  CONST double r[],
  CONST int *nk,
  CONST int *nl,
  double p[],
  double v[],
  double ar[],
  int *nvl,
  int *ifail
);

extern void __stdcall G13ADF(
  CONST int mr[],
  CONST double r[],
  CONST int *nk,
  CONST double *xv,
  CONST int *npar,
  double wa[],
  CONST int *nwa,
  double par[],
  double *rv,
  int isf[],
  int *ifail
);

extern void __stdcall G13AEF(
  CONST int mr[],
  double par[],
  CONST int *npar,
  double *c,
  CONST int *kfc,
  CONST double x[],
  CONST int *nx,
  int icount[],
  double ex[],
  double exr[],
  double al[],
  CONST int *iex,
  double *s,
  double g[],
  CONST int *igh,
  double sd[],
  double h[] /* 2 dimension */,
  CONST int *ih,
  double st[],
  CONST int *ist,
  int *nst,
  void (__stdcall *piv)(int[], double[], int *, double *, int *, int[], double *, 
              double[], double[], int *, int *, int *, double[]),
  CONST int *kpiv,
  CONST int *nit,
  int *itc,
  double zsp[],
  CONST int *kzsp,
  int isf[],
  double wa[],
  CONST int *iwa,
  double hc[] /* 2 dimension */,
  int *ifail
);

extern void __stdcall G13AFF(
  CONST int mr[],
  double par[],
  CONST int *npar,
  double *c,
  CONST int *kfc,
  CONST double x[],
  CONST int *nx,
  double *s,
  int *ndf,
  double sd[],
  CONST int *nppc,
  double cm[] /* 2 dimension */,
  CONST int *icm,
  double st[],
  int *nst,
  CONST int *kpiv,
  CONST int *nit,
  int *itc,
  int isf[],
  double res[],
  CONST int *ires,
  int *nres,
  int *ifail
);

extern void __stdcall G13AGF(
  double st[],
  CONST int *nst,
  CONST int mr[],
  CONST double par[],
  CONST int *npar,
  CONST double *c,
  double anx[],
  CONST int *nuv,
  double anexr[],
  double wa[],
  CONST int *nwa,
  int *ifail
);

extern void __stdcall G13AHF(
  CONST double st[],
  CONST int *nst,
  CONST int mr[],
  CONST double par[],
  CONST int *npar,
  double *c,
  CONST double *rms,
  CONST int *nfv,
  double fva[],
  double fsd[],
  double wa[],
  CONST int *nwa,
  int *ifail
);

extern void __stdcall G13AJF(
  CONST int mr[],
  CONST double par[],
  CONST int *npar,
  CONST double *c,
  CONST int *kfc,
  CONST double x[],
  CONST int *nx,
  double *rms,
  double st[],
  CONST int *ist,
  int *nst,
  CONST int *nfv,
  double fva[],
  double fsd[],
  CONST int *ifv,
  int isf[],
  double w[],
  CONST int *iw,
  int *ifail
);

extern void __stdcall G13ASF(
  CONST int *n,
  CONST double v[],
  CONST int mr[],
  CONST int *m,
  CONST double par[],
  CONST int *npar,
  CONST int *ishow,
  double c[],
  double acfvar[] /* 2 dimension */,
  CONST int *im,
  double *sum2,
  int *idf,
  double *siglev,
  int intgr[],
  CONST int *lmax,
  double work[],
  CONST int *lwork,
  int *ifail
);

extern void __stdcall G13AUF(
  CONST int *n,
  CONST double z[],
  CONST int *m,
  CONST int *k,
  CONST char *rs, CONST int rs_len,
  double y[],
  double mean[],
  int *ifail
);

extern void __stdcall G13BAF(
  CONST double y[],
  CONST int *ny,
  CONST int mr[],
  CONST int *nmr,
  CONST double par[],
  CONST int *npar,
  double *cy,
  double wa[],
  CONST int *nwa,
  double b[],
  CONST int *nb,
  int *ifail
);

extern void __stdcall G13BBF(
  CONST double y[],
  CONST int *ny,
  CONST int mr[],
  CONST int *nmr,
  CONST double par[],
  CONST int *npar,
  double *cy,
  double wa[],
  CONST int *iwa,
  double b[],
  CONST int *nb,
  int *ifail
);

extern void __stdcall G13BCF(
  CONST double x[],
  CONST double y[],
  CONST int *nxy,
  CONST int *nl,
  double *s,
  double *r0,
  double r[],
  double *stat,
  int *ifail
);

extern void __stdcall G13BDF(
  CONST double *r0,
  CONST double r[],
  CONST int *nl,
  CONST int nna[],
  CONST double *s,
  CONST int *nwds,
  double wa[],
  CONST int *iwa,
  double wds[],
  int isf[],
  int *ifail
);

extern void __stdcall G13BEF(
  CONST int mr[],
  CONST int *nser,
  CONST int mt[] /* 2 dimension */,
  double para[],
  CONST int *npara,
  CONST int *kfc,
  CONST int *nxxy,
  double xxy[] /* 2 dimension */,
  CONST int *ixxy,
  CONST int *kef,
  CONST int *nit,
  CONST int *kzsp,
  double zsp[],
  int *itc,
  double sd[],
  double cm[] /* 2 dimension */,
  CONST int *icm,
  double *s,
  double *d,
  int *ndf,
  CONST int *kzef,
  double res[],
  double sttf[],
  CONST int *isttf,
  int *nsttf,
  double wa[],
  CONST int *iwa,
  int mwa[],
  CONST int *imwa,
  CONST int *kpriv,
  int *ifail
);

extern void __stdcall G13BGF(
  double sttf[],
  CONST int *nsttf,
  CONST int mr[],
  CONST int *nser,
  CONST int mt[] /* 2 dimension */,
  CONST double para[],
  CONST int *npara,
  CONST int *nnv,
  double xxyn[] /* 2 dimension */,
  CONST int *ixxyn,
  CONST int *kzef,
  double res[],
  double wa[],
  CONST int *iwa,
  int *ifail
);

extern void __stdcall G13BHF(
  double sttf[],
  CONST int *nsttf,
  int mr[],
  CONST int *nser,
  CONST int mt[] /* 2 dimension */,
  CONST double para[],
  CONST int *npara,
  CONST int *nfv,
  double xxyn[] /* 2 dimension */,
  CONST int *ixxyn,
  int mrx[] /* 2 dimension */,
  CONST double parx[] /* 2 dimension */,
  CONST int *iparx,
  CONST double rmsxy[],
  CONST int *kzef,
  double fva[],
  double fsd[],
  double wa[],
  CONST int *iwa,
  int *ifail
);

extern void __stdcall G13BJF(
  int mr[],
  CONST int *nser,
  CONST int mt[] /* 2 dimension */,
  double para[],
  CONST int *npara,
  CONST int *kfc,
  CONST int *nev,
  CONST int *nfv,
  double xxy[] /* 2 dimension */,
  CONST int *ixxy,
  CONST int *kzef,
  double rmsxy[],
  int mrx[] /* 2 dimension */,
  CONST double parx[] /* 2 dimension */,
  CONST int *iparx,
  double fva[],
  double fsd[],
  double sttf[],
  CONST int *isttf,
  int *nsttf,
  double wa[],
  CONST int *iwa,
  int mwa[],
  CONST int *imwa,
  int *ifail
);

extern void __stdcall G13CAF(
  CONST int *nx,
  CONST int *mtx,
  CONST double *px,
  CONST int *iw,
  CONST int *mw,
  CONST int *ic,
  CONST int *nc,
  double c[],
  CONST int *kc,
  CONST int *l,
  CONST int *lg,
  CONST int *nxg,
  double xg[],
  int *ng,
  double stats[],
  int *ifail
);

extern void __stdcall G13CBF(
  CONST int *nx,
  CONST int *mtx,
  CONST double *px,
  CONST int *mw,
  CONST double *pw,
  CONST int *l,
  CONST int *kc,
  CONST int *lg,
  double xg[],
  int *ng,
  double stats[],
  int *ifail
);

extern void __stdcall G13CCF(
  CONST int *nxy,
  CONST int *mtxy,
  CONST double *pxy,
  CONST int *iw,
  CONST int *mw,
  CONST int *is,
  CONST int *ic,
  CONST int *nc,
  double cxy[],
  double cyx[],
  CONST int *kc,
  CONST int *l,
  CONST int *nxyg,
  double xg[],
  double yg[],
  int *ng,
  int *ifail
);

extern void __stdcall G13CDF(
  CONST int *nxy,
  CONST int *mtxy,
  CONST double *pxy,
  CONST int *mw,
  CONST int *is,
  CONST double *pw,
  CONST int *l,
  CONST int *kc,
  double xg[],
  double yg[],
  int *ng,
  int *ifail
);

extern void __stdcall G13CEF(
  CONST double xg[],
  CONST double yg[],
  CONST double xyrg[],
  CONST double xyig[],
  CONST int *ng,
  CONST double stats[],
  double ca[],
  double calw[],
  double caup[],
  double *t,
  double sc[],
  double sclw[],
  double scup[],
  int *ifail
);

extern void __stdcall G13CFF(
  CONST double xg[],
  CONST double yg[],
  CONST double xyrg[],
  CONST double xyig[],
  CONST int *ng,
  CONST double stats[],
  double gn[],
  double gnlw[],
  double gnup[],
  double ph[],
  double phlw[],
  double phup[],
  int *ifail
);

extern void __stdcall G13CGF(
  CONST double xg[],
  CONST double yg[],
  CONST double xyrg[],
  CONST double xyig[],
  CONST int *ng,
  CONST double stats[],
  CONST int *l,
  CONST int *n,
  double er[],
  double *erlw,
  double *erup,
  double rf[],
  double *rfse,
  int *ifail
);

extern void __stdcall G13DAF(
  CONST double x[] /* 2 dimension */,
  CONST int *nxm,
  CONST int *nx,
  CONST int *nsm,
  CONST int *ns,
  CONST int *nl,
  CONST int *icr,
  double c0[] /* 2 dimension */,
  double c[] /* 3 dimension */,
  int *ifail
);

extern void __stdcall G13DBF(
  CONST double c0[] /* 2 dimension */,
  CONST double c[] /* 3 dimension */,
  CONST int *nsm,
  CONST int *ns,
  CONST int *nl,
  CONST int *nk,
  double p[],
  double *v0,
  double v[],
  double d[] /* 3 dimension */,
  double db[] /* 2 dimension */,
  double w[] /* 3 dimension */,
  double wb[] /* 3 dimension */,
  int *nvp,
  double wa[],
  CONST int *iwa,
  int *ifail
);

extern void __stdcall G13DCF(
  CONST int *k,
  CONST int *n,
  CONST int *p,
  CONST int *q,
  CONST int *mean,
  double x[],
  CONST int *n4,
  double qq[] /* 2 dimension */,
  CONST int *ik,
  CONST double w[] /* 2 dimension */,
  CONST int parhld[],
  int *conds,
  int *iprint,
  CONST double *cgetol,
  CONST int *maxcal,
  CONST int *ishow,
  int *niter,
  double *logl,
  double v[] /* 2 dimension */,
  double g[],
  double disp[] /* 2 dimension */,
  CONST int *idisp,
  double w2[],
  CONST int *lw,
  int iw[],
  CONST int *liw,
  int *ifail
);

extern void __stdcall G13DJF(
  CONST int *k,
  CONST int *n,
  CONST double z[] /* 2 dimension */,
  CONST int *ik,
  CONST char *tr, CONST int tr_len,
  CONST int id[],
  CONST double delta[] /* 2 dimension */,
  CONST int *ip,
  CONST int *iq,
  CONST char *mean, CONST int mean_len,
  CONST double par[],
  CONST int *lpar,
  double qq[] /* 2 dimension */,
  CONST double v[] /* 2 dimension */,
  CONST int *lmax,
  double predz[] /* 2 dimension */,
  double sefz[] /* 2 dimension */,
  double ref[],
  CONST int *lref,
  double work[],
  CONST int *lwork,
  int iwork[],
  CONST int *liwork,
  int *ifail
);

extern void __stdcall G13DKF(
  CONST int *k,
  CONST int *lmax,
  CONST int *m,
  int *mlast,
  CONST double z[] /* 2 dimension */,
  CONST int *ik,
  double ref[],
  CONST int *lref,
  double v[] /* 2 dimension */,
  double predz[] /* 2 dimension */,
  double sefz[] /* 2 dimension */,
  double work[],
  int *ifail
);

extern void __stdcall G13DLF(
  CONST int *k,
  CONST int *n,
  CONST double z[] /* 2 dimension */,
  CONST int *ik,
  CONST char *tr, CONST int tr_len,
  CONST int id[],
  CONST double delta[] /* 2 dimension */,
  double w[] /* 2 dimension */,
  int *nd,
  double work[],
  int *ifail
);

extern void __stdcall G13DMF(
  CONST char *matrix, CONST int matrix_len,
  CONST int *k,
  CONST int *n,
  CONST int *m,
  CONST double w[] /* 2 dimension */,
  CONST int *ik,
  double wmean[],
  double r0[] /* 2 dimension */,
  double r[] /* 3 dimension */,
  int *ifail
);

extern void __stdcall G13DNF(
  CONST int *k,
  CONST int *n,
  CONST int *m,
  CONST int *ik,
  CONST double r0[] /* 2 dimension */,
  CONST double r[] /* 3 dimension */,
  int *maxlag,
  double parlag[] /* 3 dimension */,
  double x[],
  double pvalue[],
  double work[],
  CONST int *lwork,
  int *ifail
);

extern void __stdcall G13DPF(
  CONST int *k,
  CONST int *n,
  CONST double z[] /* 2 dimension */,
  CONST int *ik,
  CONST int *m,
  int *maxlag,
  double parlag[] /* 3 dimension */,
  double se[] /* 3 dimension */,
  double qq[] /* 3 dimension */,
  double x[],
  double pvalue[],
  double loglhd[],
  double work[],
  CONST int *lwork,
  int iwork[],
  int *ifail
);

extern void __stdcall G13DSF(
  CONST int *k,
  CONST int *n,
  CONST double v[] /* 2 dimension */,
  CONST int *ik,
  CONST int *p,
  CONST int *q,
  CONST int *m,
  CONST double par[],
  CONST int parhld[],
  double qq[] /* 2 dimension */,
  CONST int *ishow,
  double r0[] /* 2 dimension */,
  double c[] /* 3 dimension */,
  double acfvar[] /* 2 dimension */,
  CONST int *im,
  double *chi,
  int *idf,
  double *siglev,
  int iw[],
  CONST int *liw,
  double work[],
  CONST int *lwork,
  int *ifail
);

extern void __stdcall G13DXF(
  CONST int *k,
  CONST int *ip,
  CONST double par[],
  double rr[],
  double ri[],
  double rmod[],
  double work[],
  int iwork[],
  int *ifail
);

extern void __stdcall G13EAF(
  CONST int *n,
  CONST int *m,
  CONST int *l,
  CONST double a[] /* 2 dimension */,
  CONST int *lds,
  CONST double b[] /* 2 dimension */,
  CONST int *stq,
  CONST double q[] /* 2 dimension */,
  CONST int *ldq,
  CONST double c[] /* 2 dimension */,
  CONST int *ldm,
  CONST double r[] /* 2 dimension */,
  double s[] /* 2 dimension */,
  double k[] /* 2 dimension */,
  double h[] /* 2 dimension */,
  CONST double *tol,
  CONST int iwk[],
  double wk[],
  int *ifail
);

extern void __stdcall G13EBF(
  CONST char *transf, CONST int transf_len,
  CONST int *n,
  CONST int *m,
  CONST int *l,
  double a[] /* 2 dimension */,
  CONST int *lds,
  double b[] /* 2 dimension */,
  CONST int *stq,
  CONST double q[] /* 2 dimension */,
  CONST int *ldq,
  double c[] /* 2 dimension */,
  CONST int *ldm,
  CONST double r[] /* 2 dimension */,
  double s[] /* 2 dimension */,
  double k[] /* 2 dimension */,
  double h[] /* 2 dimension */,
  double u[] /* 2 dimension */,
  CONST double *tol,
  CONST int iwk[],
  double wk[],
  int *ifail
);

extern void __stdcall H02BBF(
  int *itmax,
  CONST int *msglvl,
  CONST int *n,
  CONST int *m,
  CONST double a[] /* 2 dimension */,
  CONST int *lda,
  CONST double bl[],
  CONST double bu[],
  CONST int intvar[],
  CONST double cvec[],
  CONST int *maxnod,
  CONST int *intfst,
  CONST int *maxdpt,
  double *toliv,
  double *tolfes,
  double *bigbnd,
  double x[],
  double *objmip,
  int iwork[],
  CONST int *liwork,
  double rwork[],
  CONST int *lrwork,
  int *ifail
);

extern void __stdcall H02BFF(
  CONST int *infile,
  CONST int *maxn,
  CONST int *maxm,
  CONST char *optim, CONST int optim_len,
  CONST double *xbldef,
  CONST double *xbudef,
  CONST int *maxdpt,
  CONST int *msglvl,
  int *n,
  int *m,
  double x[],
  char *crname, CONST int crname_len,
  int iwork[],
  CONST int *liwork,
  double rwork[],
  CONST int *lrwork,
  int *ifail
);

extern void __stdcall H02BUF(
  CONST int *infile,
  CONST int *maxn,
  CONST int *maxm,
  CONST char *optim, CONST int optim_len,
  CONST double *xbldef,
  CONST double *xbudef,
  char *nmobj, CONST int nmobj_len,
  char *nmrhs, CONST int nmrhs_len,
  char *nmrng, CONST int nmrng_len,
  char *nmbnd, CONST int nmbnd_len,
  CONST int *mpslst,
  int *n,
  int *m,
  double a[] /* 2 dimension */,
  double bl[],
  double bu[],
  double cvec[],
  double x[],
  int intvar[],
  char *crname, CONST int crname_len,
  char *nmprob, CONST int nmprob_len,
  int iwork[],
  int *ifail
);

extern void __stdcall H02BVF(
  CONST int *n,
  CONST int *m,
  CONST double a[] /* 2 dimension */,
  CONST int *lda,
  CONST double bl[],
  CONST double bu[],
  CONST double x[],
  CONST double clamda[],
  CONST int istate[],
  CONST char *crname, CONST int crname_len,
  int *ifail
);

extern void __stdcall H02BZF(
  CONST int *n,
  CONST int *m,
  double bl[],
  double bu[],
  double clamda[],
  int istate[],
  CONST int iwork[],
  CONST int *liwork,
  CONST double rwork[],
  CONST int *lrwork,
  int *ifail
);

extern void __stdcall H03ABF(
  CONST int kost[] /* 2 dimension */,
  CONST int *mmm,
  CONST int *ma,
  CONST int *mb,
  CONST int *m,
  int k15[],
  CONST int *maxit,
  int k7[],
  int k9[],
  int *numit,
  int k6[],
  int k8[],
  int k11[],
  int k12[],
  double *z,
  int *ifail
);

extern void __stdcall H03ADF(
  CONST int *n,
  CONST int *ns,
  CONST int *ne,
  CONST int *direct,
  CONST int *nnz,
  CONST double d[],
  int irow[],
  CONST int icol[],
  double *splen,
  int path[],
  int iwork[],
  double work[],
  int *ifail
);

extern void __stdcall M01AJF(
  double a[],
  double w[],
  int ind[],
  int indw[],
  CONST int *n,
  CONST int *nw,
  int *ifail
);

extern void __stdcall M01AKF(
  double a[],
  double w[],
  int ind[],
  int indw[],
  CONST int *n,
  CONST int *nw,
  int *ifail
);

extern void __stdcall M01APF(
  double a[],
  CONST int *ii,
  CONST int *jj,
  int *ifail
);

extern void __stdcall M01CAF(
  double rv[],
  CONST int *m1,
  CONST int *m2,
  CONST char *order, CONST int order_len,
  int *ifail
);

extern void __stdcall M01CBF(
  int iv[],
  CONST int *m1,
  CONST int *m2,
  CONST char *order, CONST int order_len,
  int *ifail
);

extern void __stdcall M01CCF(
  char *ch, CONST int ch_len,
  CONST int *m1,
  CONST int *m2,
  CONST int *l1,
  CONST int *l2,
  CONST char *order, CONST int order_len,
  int *ifail
);

extern void __stdcall M01DAF(
  CONST double rv[],
  CONST int *m1,
  CONST int *m2,
  CONST char *order, CONST int order_len,
  int irank[],
  int *ifail
);

extern void __stdcall M01DBF(
  CONST int iv[],
  CONST int *m1,
  CONST int *m2,
  CONST char *order, CONST int order_len,
  int irank[],
  int *ifail
);

extern void __stdcall M01DCF(
  CONST char *ch, CONST int ch_len,
  CONST int *m1,
  CONST int *m2,
  CONST int *l1,
  CONST int *l2,
  CONST char *order, CONST int order_len,
  int irank[],
  int *ifail
);

extern void __stdcall M01DEF(
  CONST double rm[] /* 2 dimension */,
  CONST int *ldm,
  CONST int *m1,
  CONST int *m2,
  CONST int *n1,
  CONST int *n2,
  CONST char *order, CONST int order_len,
  int irank[],
  int *ifail
);

extern void __stdcall M01DFF(
  CONST int iv[] /* 2 dimension */,
  CONST int *ldm,
  CONST int *m1,
  CONST int *m2,
  CONST int *n1,
  CONST int *n2,
  CONST char *order, CONST int order_len,
  int irank[],
  int *ifail
);

extern void __stdcall M01DJF(
  CONST double rm[] /* 2 dimension */,
  CONST int *ldm,
  CONST int *m1,
  CONST int *m2,
  CONST int *n1,
  CONST int *n2,
  CONST char *order, CONST int order_len,
  int irank[],
  int *ifail
);

extern void __stdcall M01DKF(
  CONST int im[] /* 2 dimension */,
  CONST int *ldm,
  CONST int *m1,
  CONST int *m2,
  CONST int *n1,
  CONST int *n2,
  CONST char *order, CONST int order_len,
  int irank[],
  int *ifail
);

extern void __stdcall M01DZF(
  int (__stdcall *compar)(int *, int *),
  CONST int *m1,
  CONST int *m2,
  int irank[],
  int *ifail
);

extern void __stdcall M01EAF(
  double rv[],
  CONST int *m1,
  CONST int *m2,
  int irank[],
  int *ifail
);

extern void __stdcall M01EBF(
  int iv[],
  CONST int *m1,
  CONST int *m2,
  int irank[],
  int *ifail
);

extern void __stdcall M01ECF(
  char *ch, CONST int ch_len,
  CONST int *m1,
  CONST int *m2,
  int irank[],
  int *ifail
);

extern void __stdcall M01ZAF(
  int iperm[],
  CONST int *m1,
  CONST int *m2,
  int *ifail
);

extern void __stdcall M01ZBF(
  int iperm[],
  CONST int *m1,
  CONST int *m2,
  int *ifail
);

extern void __stdcall M01ZCF(
  int iperm[],
  CONST int *m1,
  CONST int *m2,
  int icycl[],
  int *ifail
);

extern int __stdcall P01ABF(
  CONST int *ifail,
  CONST int *ierror,
  CONST char *srname, CONST int srname_len,
  CONST int *nrec,
  CONST char *rec, int rec_len
);

extern int __stdcall P01ACF(
  CONST int *ifail,
  CONST int *ierror,
  CONST char *srname, CONST int srname_len,
  CONST char *varbnm, CONST int varbnm_len,
  CONST int *nrec,
  CONST char *rec, int rec_len
);

extern double __stdcall S01BAF(
  CONST double *x,
  int *ifail
);

extern void __stdcall S01EAF(
  Complex * S01EAF_,
  CONST Complex *z,
  int *ifail
);

extern double __stdcall S07AAF(
  CONST double *x,
  int *ifail
);

extern double __stdcall S09AAF(
  CONST double *x,
  int *ifail
);

extern double __stdcall S09ABF(
  CONST double *x,
  int *ifail
);

extern double __stdcall S10AAF(
  CONST double *x,
  int *ifail
);

extern double __stdcall S10ABF(
  CONST double *x,
  int *ifail
);

extern double __stdcall S10ACF(
  CONST double *x,
  int *ifail
);

extern double __stdcall S11AAF(
  CONST double *x,
  int *ifail
);

extern double __stdcall S11ABF(
  CONST double *x,
  int *ifail
);

extern double __stdcall S11ACF(
  CONST double *x,
  int *ifail
);

extern double __stdcall S13AAF(
  CONST double *x,
  int *ifail
);

extern double __stdcall S13ACF(
  CONST double *x,
  int *ifail
);

extern double __stdcall S13ADF(
  CONST double *x,
  int *ifail
);

extern double __stdcall S14AAF(
  CONST double *x,
  int *ifail
);

extern double __stdcall S14ABF(
  CONST double *x,
  int *ifail
);

extern double __stdcall S14ACF(
  CONST double *x,
  int *ifail
);

extern void __stdcall S14ADF(
  CONST double *x,
  CONST int *n,
  CONST int *m,
  double w[],
  int *ifail
);

extern void __stdcall S14BAF(
  CONST double *a,
  CONST double *x,
  CONST double *tol,
  double *p,
  double *q,
  int *ifail
);

extern double __stdcall S15ABF(
  CONST double *x,
  int *ifail
);

extern double __stdcall S15ACF(
  CONST double *x,
  int *ifail
);

extern double __stdcall S15ADF(
  CONST double *x,
  int *ifail
);

extern double __stdcall S15AEF(
  CONST double *x,
  int *ifail
);

extern double __stdcall S15AFF(
  CONST double *x,
  int *ifail
);

extern void __stdcall S15DDF(
  Complex * S15DDF_,
  CONST Complex *z,
  int *ifail
);

extern double __stdcall S17ACF(
  CONST double *x,
  int *ifail
);

extern double __stdcall S17ADF(
  CONST double *x,
  int *ifail
);

extern double __stdcall S17AEF(
  CONST double *x,
  int *ifail
);

extern double __stdcall S17AFF(
  CONST double *x,
  int *ifail
);

extern double __stdcall S17AGF(
  CONST double *x,
  int *ifail
);

extern double __stdcall S17AHF(
  CONST double *x,
  int *ifail
);

extern double __stdcall S17AJF(
  CONST double *x,
  int *ifail
);

extern double __stdcall S17AKF(
  CONST double *x,
  int *ifail
);

extern void __stdcall S17DCF(
  CONST double *fnu,
  CONST double z[],
  CONST int *n,
  CONST char *scale, CONST int scale_len,
  double cy[] /* 2 dimension */,
  int *nz,
  double cwrk[] /* 2 dimension */,
  int *ifail
);

extern void __stdcall S17DEF(
  CONST double *fnu,
  CONST double z[],
  CONST int *n,
  CONST char *scale, CONST int scale_len,
  double cy[] /* 2 dimension */,
  int *nz,
  int *ifail
);

extern void __stdcall S17DGF(
  CONST char *deriv, CONST int deriv_len,
  CONST double z[],
  CONST char *scale, CONST int scale_len,
  double ai[],
  int *nz,
  int *ifail
);

extern void __stdcall S17DHF(
  CONST char *deriv, CONST int deriv_len,
  CONST double z[],
  CONST char *scale, CONST int scale_len,
  double bi[],
  int *ifail
);

extern void __stdcall S17DLF(
  CONST int *m,
  CONST double *fnu,
  CONST double z[],
  CONST int *n,
  CONST char *scale, CONST int scale_len,
  double cy[] /* 2 dimension */,
  int *nz,
  int *ifail
);

extern double __stdcall S18ACF(
  CONST double *x,
  int *ifail
);

extern double __stdcall S18ADF(
  CONST double *x,
  int *ifail
);

extern double __stdcall S18AEF(
  CONST double *x,
  int *ifail
);

extern double __stdcall S18AFF(
  CONST double *x,
  int *ifail
);

extern double __stdcall S18CCF(
  CONST double *x,
  int *ifail
);

extern double __stdcall S18CDF(
  CONST double *x,
  int *ifail
);

extern double __stdcall S18CEF(
  CONST double *x,
  int *ifail
);

extern double __stdcall S18CFF(
  CONST double *x,
  int *ifail
);

extern void __stdcall S18DCF(
  CONST double *fnu,
  CONST double z[],
  CONST int *n,
  CONST char *scale, CONST int scale_len,
  double cy[] /* 2 dimension */,
  int *nz,
  int *ifail
);

extern void __stdcall S18DEF(
  CONST double *fnu,
  CONST double z[],
  CONST int *n,
  CONST char *scale, CONST int scale_len,
  double cy[] /* 2 dimension */,
  int *nz,
  int *ifail
);

extern double __stdcall S19AAF(
  CONST double *x,
  int *ifail
);

extern double __stdcall S19ABF(
  CONST double *x,
  int *ifail
);

extern double __stdcall S19ACF(
  CONST double *x,
  int *ifail
);

extern double __stdcall S19ADF(
  CONST double *x,
  int *ifail
);

extern double __stdcall S20ACF(
  CONST double *x,
  int *ifail
);

extern double __stdcall S20ADF(
  CONST double *x,
  int *ifail
);

extern double __stdcall S21BAF(
  CONST double *x,
  CONST double *y,
  int *ifail
);

extern double __stdcall S21BBF(
  CONST double *x,
  CONST double *y,
  CONST double *z,
  int *ifail
);

extern double __stdcall S21BCF(
  CONST double *x,
  CONST double *y,
  CONST double *z,
  int *ifail
);

extern double __stdcall S21BDF(
  CONST double *x,
  CONST double *y,
  CONST double *z,
  CONST double *r,
  int *ifail
);

extern void __stdcall S21CAF(
  CONST double *u,
  CONST double *m,
  double *sn,
  double *cn,
  double *dn,
  int *ifail
);

extern double __stdcall X01AAF(
  CONST double *x
);

extern double __stdcall X01ABF(
  CONST double *x
);

extern double __stdcall X02AAF(
  CONST double *x
);

extern double __stdcall X02ABF(
  CONST double *x
);

extern double __stdcall X02ACF(
  CONST double *x
);

extern double __stdcall X02AGF(
  CONST double *x
);

extern double __stdcall X02AHF(
  CONST double *x
);

extern double __stdcall X02AJF(
void
);

extern double __stdcall X02AKF(
void
);

extern double __stdcall X02ALF(
void
);

extern double __stdcall X02AMF(
void
);

extern double __stdcall X02ANF(
void
);

extern int __stdcall X02BBF(
  CONST double *x
);

extern int __stdcall X02BEF(
  CONST double *x
);

extern int __stdcall X02BHF(
void
);

extern int __stdcall X02BJF(
void
);

extern int __stdcall X02BKF(
void
);

extern int __stdcall X02BLF(
void
);

extern int __stdcall X02CAF(
  CONST double *x
);

extern int __stdcall X02DAF(
  CONST double *x
);

extern int __stdcall X02DJF(
void
);

extern void __stdcall X03AAF(
  CONST double a[],
  CONST int *isizea,
  CONST double b[],
  CONST int *isizeb,
  CONST int *n,
  CONST int *istepa,
  CONST int *istepb,
  CONST double *c1,
  CONST double *c2,
  double *d1,
  double *d2,
  CONST int *sw,
  int *ifail
);

extern void __stdcall X03ABF(
  CONST Complex a[],
  CONST int *isizea,
  CONST Complex b[],
  CONST int *isizeb,
  CONST int *n,
  CONST int *istepa,
  CONST int *istepb,
  CONST Complex *cx,
  Complex *dx,
  CONST int *sw,
  int *ifail
);

extern void __stdcall X04AAF(
  CONST int *i,
  int *nerr
);

extern void __stdcall X04ABF(
  CONST int *i,
  int *nadv
);

extern void __stdcall X04BAF(
  CONST int *nout,
  CONST char *rec, int rec_len
);

extern void __stdcall X04BBF(
  CONST int *nin,
  char *rec, CONST int rec_len,
  int *ifail
);

extern void __stdcall X04CAF(
  CONST char *matrix, CONST int matrix_len,
  CONST char *diag, CONST int diag_len,
  CONST int *m,
  CONST int *n,
  CONST double a[] /* 2 dimension */,
  CONST int *lda,
  CONST char *title, CONST int title_len,
  int *ifail
);

extern void __stdcall X04CBF(
  CONST char *matrix, CONST int matrix_len,
  CONST char *diag, CONST int diag_len,
  CONST int *m,
  CONST int *n,
  CONST double a[] /* 2 dimension */,
  CONST int *lda,
  CONST char *format, CONST int format_len,
  CONST char *title, CONST int title_len,
  CONST char *labrow, CONST int labrow_len,
  CONST char *rlabs, CONST int rlabs_len,
  CONST char *labcol, CONST int labcol_len,
  CONST char *clabs, CONST int clabs_len,
  CONST int *ncols,
  CONST int *indent,
  int *ifail
);

extern void __stdcall X04CCF(
  CONST char *uplo, CONST int uplo_len,
  CONST char *diag, CONST int diag_len,
  CONST int *n,
  CONST double a[],
  CONST char *title, CONST int title_len,
  int *ifail
);

extern void __stdcall X04CDF(
  CONST char *uplo, CONST int uplo_len,
  CONST char *diag, CONST int diag_len,
  CONST int *n,
  CONST double a[],
  CONST char *format, CONST int format_len,
  CONST char *title, CONST int title_len,
  CONST char *labrow, CONST int labrow_len,
  CONST char *rlabs, CONST int rlabs_len,
  CONST char *labcol, CONST int labcol_len,
  CONST char *clabs, CONST int clabs_len,
  CONST int *ncols,
  CONST int *indent,
  int *ifail
);

extern void __stdcall X04CEF(
  CONST int *m,
  CONST int *n,
  CONST int *kl,
  CONST int *ku,
  CONST double a[] /* 2 dimension */,
  CONST int *lda,
  CONST char *title, CONST int title_len,
  int *ifail
);

extern void __stdcall X04CFF(
  CONST int *m,
  CONST int *n,
  CONST int *kl,
  CONST int *ku,
  CONST double a[] /* 2 dimension */,
  CONST int *lda,
  CONST char *format, CONST int format_len,
  CONST char *title, CONST int title_len,
  CONST char *labrow, CONST int labrow_len,
  CONST char *rlabs, CONST int rlabs_len,
  CONST char *labcol, CONST int labcol_len,
  CONST char *clabs, CONST int clabs_len,
  CONST int *ncols,
  CONST int *indent,
  int *ifail
);

extern void __stdcall X04DAF(
  CONST char *matrix, CONST int matrix_len,
  CONST char *diag, CONST int diag_len,
  CONST int *m,
  CONST int *n,
  CONST Complex a[] /* 2 dimension */,
  CONST int *lda,
  CONST char *title, CONST int title_len,
  int *ifail
);

extern void __stdcall X04DBF(
  CONST char *matrix, CONST int matrix_len,
  CONST char *diag, CONST int diag_len,
  CONST int *m,
  CONST int *n,
  CONST Complex a[] /* 2 dimension */,
  CONST int *lda,
  CONST char *usefrm, CONST int usefrm_len,
  CONST char *format, CONST int format_len,
  CONST char *title, CONST int title_len,
  CONST char *labrow, CONST int labrow_len,
  CONST char *rlabs, CONST int rlabs_len,
  CONST char *labcol, CONST int labcol_len,
  CONST char *clabs, CONST int clabs_len,
  CONST int *ncols,
  CONST int *indent,
  int *ifail
);

extern void __stdcall X04DCF(
  CONST char *uplo, CONST int uplo_len,
  CONST char *diag, CONST int diag_len,
  CONST int *n,
  CONST Complex a[],
  CONST char *title, CONST int title_len,
  int *ifail
);

extern void __stdcall X04DDF(
  CONST char *uplo, CONST int uplo_len,
  CONST char *diag, CONST int diag_len,
  CONST int *n,
  CONST Complex a[],
  CONST char *usefrm, CONST int usefrm_len,
  CONST char *format, CONST int format_len,
  CONST char *title, CONST int title_len,
  CONST char *labrow, CONST int labrow_len,
  CONST char *rlabs, CONST int rlabs_len,
  CONST char *labcol, CONST int labcol_len,
  CONST char *clabs, CONST int clabs_len,
  CONST int *ncols,
  CONST int *indent,
  int *ifail
);

extern void __stdcall X04DEF(
  CONST int *m,
  CONST int *n,
  CONST int *kl,
  CONST int *ku,
  CONST Complex a[] /* 2 dimension */,
  CONST int *lda,
  CONST char *title, CONST int title_len,
  int *ifail
);

extern void __stdcall X04DFF(
  CONST int *m,
  CONST int *n,
  CONST int *kl,
  CONST int *ku,
  CONST Complex a[] /* 2 dimension */,
  CONST int *lda,
  CONST char *usefrm, CONST int usefrm_len,
  CONST char *format, CONST int format_len,
  CONST char *title, CONST int title_len,
  CONST char *labrow, CONST int labrow_len,
  CONST char *rlabs, CONST int rlabs_len,
  CONST char *labcol, CONST int labcol_len,
  CONST char *clabs, CONST int clabs_len,
  CONST int *ncols,
  CONST int *indent,
  int *ifail
);

extern void __stdcall X04EAF(
  CONST char *matrix, CONST int matrix_len,
  CONST char *diag, CONST int diag_len,
  CONST int *m,
  CONST int *n,
  CONST int a[] /* 2 dimension */,
  CONST int *lda,
  CONST char *title, CONST int title_len,
  int *ifail
);

extern void __stdcall X04EBF(
  CONST char *matrix, CONST int matrix_len,
  CONST char *diag, CONST int diag_len,
  CONST int *m,
  CONST int *n,
  CONST int a[] /* 2 dimension */,
  CONST int *lda,
  CONST char *format, CONST int format_len,
  CONST char *title, CONST int title_len,
  CONST char *labrow, CONST int labrow_len,
  CONST char *rlabs, CONST int rlabs_len,
  CONST char *labcol, CONST int labcol_len,
  CONST char *clabs, CONST int clabs_len,
  CONST int *ncols,
  CONST int *indent,
  int *ifail
);

extern void __stdcall X05AAF(
  int itime[]
);

extern void __stdcall X05ABF(
  char *, CONST int _len,
  CONST int itime[]
);

extern int __stdcall X05ACF(
  CONST char *ctime1, CONST int ctime1_len,
  CONST char *ctime2, int ctime2_len
);

extern double __stdcall X05BAF(
void
);

extern void __stdcall Y07AAF(
  CONST int *nout
);

extern void __stdcall Y07ABF(
  CONST int *nout
);

extern void __stdcall Y07ACF(
  CONST int *nout
);

extern void __stdcall Y07ADF(
  CONST int *nin,
  CONST int *nout
);

extern void __stdcall Y07AEF(
  CONST int *nin,
  CONST int *nout
);

extern void __stdcall Y07BAF(
  CONST int *n,
  CONST double x[],
  CONST int *incx,
  CONST char *title, CONST int title_len,
  CONST int *nsig,
  CONST int *ncols,
  CONST int *nout
);

extern void __stdcall Y07BBF(
  CONST int *n,
  CONST double x[],
  CONST int *incx,
  CONST char *title, CONST int title_len,
  CONST int *nfigb,
  CONST int *nfiga,
  CONST int *ncols,
  CONST int *nout
);

extern void __stdcall Y07BDF(
  CONST char *job, CONST int job_len,
  CONST int *m,
  CONST int *n,
  CONST double a[] /* 2 dimension */,
  CONST int *lda,
  CONST char *title, CONST int title_len,
  CONST int *nsig,
  CONST int *ncols,
  CONST int *nout
);

extern void __stdcall Y07BEF(
  CONST char *job, CONST int job_len,
  CONST int *m,
  CONST int *n,
  CONST double a[] /* 2 dimension */,
  CONST int *lda,
  CONST char *title, CONST int title_len,
  CONST int *nfigb,
  CONST int *nfiga,
  CONST int *ncols,
  CONST int *nout
);

extern void __stdcall Y07BFF(
  CONST char *job, CONST int job_len,
  CONST int *n,
  CONST double a[],
  CONST char *title, CONST int title_len,
  CONST int *nsig,
  CONST int *ncols,
  CONST int *nout
);

extern void __stdcall Y07BGF(
  CONST char *job, CONST int job_len,
  CONST int *n,
  CONST double a[],
  CONST char *title, CONST int title_len,
  CONST int *nfigb,
  CONST int *nfiga,
  CONST int *ncols,
  CONST int *nout
);

extern void __stdcall Y07CAF(
  CONST int *n,
  CONST Complex x[],
  CONST int *incx,
  CONST char *title, CONST int title_len,
  CONST int *nsig,
  CONST int *ncols,
  CONST int *nout
);

extern void __stdcall Y07CBF(
  CONST int *n,
  CONST Complex x[],
  CONST int *incx,
  CONST char *title, CONST int title_len,
  CONST int *nfigb,
  CONST int *nfiga,
  CONST int *ncols,
  CONST int *nout
);

extern void __stdcall Y07CDF(
  CONST char *job, CONST int job_len,
  CONST int *m,
  CONST int *n,
  CONST Complex a[] /* 2 dimension */,
  CONST int *lda,
  CONST char *title, CONST int title_len,
  CONST int *nsig,
  CONST int *ncols,
  CONST int *nout
);

extern void __stdcall Y07CEF(
  CONST char *job, CONST int job_len,
  CONST int *m,
  CONST int *n,
  CONST Complex a[] /* 2 dimension */,
  CONST int *lda,
  CONST char *title, CONST int title_len,
  CONST int *nfigb,
  CONST int *nfiga,
  CONST int *ncols,
  CONST int *nout
);

extern void __stdcall Y07CFF(
  CONST char *job, CONST int job_len,
  CONST int *n,
  CONST Complex a[],
  CONST char *title, CONST int title_len,
  CONST int *nsig,
  CONST int *ncols,
  CONST int *nout
);

extern void __stdcall Y07CGF(
  CONST char *job, CONST int job_len,
  CONST int *n,
  CONST Complex a[],
  CONST char *title, CONST int title_len,
  CONST int *nfigb,
  CONST int *nfiga,
  CONST int *ncols,
  CONST int *nout
);

extern void __stdcall Y07DAF(
  CONST int *n,
  CONST int k[],
  CONST int *inck,
  CONST char *title, CONST int title_len,
  CONST int *nfig,
  CONST int *ncols,
  CONST int *nout
);

extern void __stdcall Y90AAF(
  CONST char *nag, CONST int nag_len,
  CONST int *tnag
);

extern void __stdcall Y90ABF(
  CONST char *prog, CONST int prog_len,
  CONST int *tprog
);

extern void __stdcall Y90ACF(
  CONST char *task, CONST int task_len,
  CONST int *ttask
);

extern void __stdcall Y90ADF(
  CONST char *nag, CONST int nag_len,
  int *xnag
);

extern void __stdcall Y90AEF(
  CONST char *prog, CONST int prog_len,
  int *xprog
);

extern void __stdcall Y90AFF(
  CONST char *task, CONST int task_len,
  int *xtask
);

extern void __stdcall Y90CAF(
  CONST char *matrix, CONST int matrix_len,
  CONST int *dtype,
  CONST int *ttype,
  CONST int *m,
  CONST int *n,
  CONST int *kl,
  CONST int *ku,
  Complex a[] /* 2 dimension */,
  CONST int *ia,
  Complex d[],
  CONST Complex diag[],
  CONST Complex odiag[],
  CONST double *cond,
  CONST Complex *scale,
  Complex *detman,
  int *detexp,
  CONST int *dist,
  int seed[]
);

extern void __stdcall Y90CBF(
  CONST int *m,
  CONST int *n,
  CONST int *kl,
  CONST int *ku,
  int seed[],
  CONST Complex d[],
  Complex a[] /* 2 dimension */,
  CONST int *lda
);

extern void __stdcall Y90CCF(
  CONST int *type,
  CONST int *dtype,
  CONST int *m,
  CONST int *n,
  Complex a[] /* 2 dimension */,
  CONST int *ia,
  Complex d[],
  CONST Complex diag[],
  CONST double *cond,
  CONST Complex *scale,
  Complex *detman,
  int *detexp,
  CONST int *dist,
  int seed[],
  CONST char *pvrow, CONST int pvrow_len,
  CONST char *pvcol, CONST int pvcol_len,
  CONST int ipvrow[],
  CONST int ipvcol[],
  Complex u[] /* 2 dimension */,
  CONST int *iu,
  Complex v[] /* 2 dimension */,
  CONST int *iv,
  Complex work1[] /* 2 dimension */,
  CONST int *iwork1,
  Complex work2[] /* 2 dimension */,
  CONST int *iwork2
);

extern void __stdcall Y90CDF(
  CONST char *side, CONST int side_len,
  CONST char *init, CONST int init_len,
  CONST int *m,
  CONST int *n,
  Complex a[] /* 2 dimension */,
  CONST int *lda,
  Complex x[],
  Complex d[],
  int seed[]
);

extern void __stdcall Y90CEF(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST int *kb,
  int seed[],
  CONST Complex d[],
  Complex a[] /* 2 dimension */,
  CONST int *lda
);

extern void __stdcall Y90CFF(
  CONST int *ttype,
  CONST char *uplo, CONST int uplo_len,
  CONST int *m,
  CONST int *n,
  Complex a[] /* 2 dimension */,
  CONST int *ia,
  CONST Complex diag[],
  CONST Complex odiag[],
  Complex *detman,
  int *detexp,
  CONST int *dist,
  int seed[]
);

extern void __stdcall Y90CGF(
  CONST char *det, CONST int det_len,
  CONST int *vtype,
  CONST int *n,
  Complex v[],
  CONST Complex vbound[],
  CONST double *cond,
  CONST Complex *scale,
  Complex *detman,
  int *detexp,
  CONST int *dist,
  int seed[]
);

extern void __stdcall Y90CJF(
  CONST int *n,
  CONST int *k,
  CONST int *l,
  CONST int *dtype,
  CONST double *cond,
  CONST Complex eig[],
  Complex a[] /* 2 dimension */,
  CONST int *lda,
  Complex vec[] /* 2 dimension */,
  CONST int *ldvec,
  Complex veci[] /* 2 dimension */,
  CONST int *ldveci,
  int seed[],
  double x[],
  Complex cwk1[] /* 2 dimension */,
  CONST int *ldcwk1
);

extern void __stdcall Y90CPF(
  CONST char *uplo, CONST int uplo_len,
  CONST char *diag, CONST int diag_len,
  CONST int *ttype,
  CONST int *n,
  int *kd,
  Complex a[] /* 2 dimension */,
  CONST int *lda,
  Complex ab[] /* 2 dimension */,
  CONST int *ldab,
  Complex d[],
  CONST double *cond,
  CONST Complex *scale,
  Complex *detman,
  int *detexp,
  CONST int *dist,
  int seed[]
);

extern void __stdcall Y90DBF(
  CONST char *conv, CONST int conv_len,
  CONST char *matrix, CONST int matrix_len,
  CONST int *m,
  CONST int *n,
  CONST int *kl,
  CONST int *ku,
  CONST Complex a[] /* 2 dimension */,
  CONST int *ia,
  Complex b[] /* 2 dimension */,
  CONST int *ib
);

extern void __stdcall Y90DCF(
  CONST char *matrix, CONST int matrix_len,
  CONST int *m,
  CONST int *mn,
  CONST int *n,
  CONST int *kl,
  CONST int *ku,
  CONST Complex a[] /* 2 dimension */,
  CONST int *lda,
  CONST Complex b[] /* 2 dimension */,
  CONST int *ldb,
  Complex c[] /* 2 dimension */,
  CONST int *ldc
);

extern double __stdcall Y90DDF(
  CONST char *norm, CONST int norm_len,
  CONST char *matrix, CONST int matrix_len,
  CONST int *m,
  CONST int *n,
  CONST int *kl,
  CONST int *ku,
  CONST Complex a[] /* 2 dimension */,
  CONST int *lda
);

extern void __stdcall Y90DEF(
  CONST int *m,
  CONST int *mn,
  CONST int *n,
  CONST int *kl,
  CONST int *ku,
  CONST Complex l[] /* 2 dimension */,
  CONST int *il,
  CONST Complex u[] /* 2 dimension */,
  CONST int *iu,
  Complex a[] /* 2 dimension */,
  CONST int *ia
);

extern void __stdcall Y90DFF(
  CONST char *trans, CONST int trans_len,
  CONST int *m,
  CONST int *n,
  CONST int *kl,
  CONST int *ku,
  CONST Complex a[] /* 2 dimension */,
  CONST int *ia,
  Complex b[] /* 2 dimension */,
  CONST int *ib
);

extern void __stdcall Y90DGF(
  CONST char *type, CONST int type_len,
  CONST int *nlose,
  Complex a[] /* 2 dimension */,
  CONST int *lda,
  CONST int *m,
  CONST int *n
);

extern void __stdcall Y90DHF(
  CONST char *matrix, CONST int matrix_len,
  CONST int *m,
  CONST int *n,
  CONST Complex a[] /* 2 dimension */,
  CONST int *nda,
  Complex b[] /* 2 dimension */,
  CONST int *ndb
);

extern void __stdcall Y90DKF(
  CONST char *matrix, CONST int matrix_len,
  Complex a[] /* 2 dimension */,
  CONST int *nda,
  CONST int *n
);

extern void __stdcall Y90DMF(
  Complex *vman,
  int *vexp,
  CONST int *scale
);

extern void __stdcall Y90DNF(
  CONST char *matra, CONST int matra_len,
  CONST char *matrb, CONST int matrb_len,
  CONST int *m,
  CONST int *n,
  CONST int *mn,
  CONST Complex *alpha,
  CONST Complex a[] /* 2 dimension */,
  CONST int *ia,
  CONST Complex b[] /* 2 dimension */,
  CONST int *ib,
  CONST Complex *beta,
  Complex c[] /* 2 dimension */,
  CONST int *ic
);

extern void __stdcall Y90EBF(
  Complex * Y90EBF_,
  CONST int *idist,
  int seed[]
);

extern double __stdcall Y90GAF(
  CONST int *k
);

extern void __stdcall Y90HAF(
  CONST int *inag,
  CONST char *nag, CONST int nag_len,
  CONST int *tnag,
  CONST int *wnag,
  CONST int *xnag
);

extern void __stdcall Y90HBF(
  CONST char *prog, CONST int prog_len,
  CONST int *tprog,
  CONST int *xprog
);

extern void __stdcall Y90HCF(
  CONST int *itask,
  CONST char *task, CONST int task_len,
  CONST int *ttask,
  CONST int *wtask,
  CONST int *xtask,
  CONST char *nag, CONST int nag_len,
  CONST int *nnag
);

extern void __stdcall Y90HDF(
  CONST int *itest,
  CONST int *test,
  CONST int *warn,
  CONST char *ttest, int ttest_len
);

extern void __stdcall Y90PAF(
  CONST char *line, int line_len
);

extern void __stdcall Y90PBF(
  CONST char *matrix, CONST int matrix_len,
  CONST char *line, CONST int line_len,
  CONST int *m,
  CONST int *n,
  CONST Complex a[] /* 2 dimension */,
  CONST int *ia
);

extern void __stdcall Y90PCF(
  CONST char *matrix, CONST int matrix_len,
  CONST char *line, CONST int line_len,
  CONST char *mname1, CONST int mname1_len,
  CONST char *mname2, CONST int mname2_len,
  CONST int *m,
  CONST int *n,
  CONST Complex a[] /* 2 dimension */,
  CONST int *ia,
  CONST Complex b[] /* 2 dimension */,
  CONST int *ib
);

extern void __stdcall Y90PDF(
  CONST char *matrix, CONST int matrix_len,
  CONST char *line, CONST int line_len,
  CONST int *m,
  CONST int *n,
  CONST double a[] /* 2 dimension */,
  CONST int *lda
);

extern void __stdcall Y90PEF(
  CONST char *matrix, CONST int matrix_len,
  CONST char *line, CONST int line_len,
  CONST char *mname1, CONST int mname1_len,
  CONST char *mname2, CONST int mname2_len,
  CONST int *m,
  CONST int *n,
  CONST double a[] /* 2 dimension */,
  CONST int *ia,
  CONST double b[] /* 2 dimension */,
  CONST int *ib
);

extern void __stdcall Y90PFF(
  CONST char *line, CONST int line_len,
  CONST Complex *cvar
);

extern void __stdcall Y90PGF(
  CONST char *line, CONST int line_len,
  CONST char *varch, int varch_len
);

extern void __stdcall Y90PHF(
  CONST char *line, CONST int line_len,
  CONST int *ivar
);

extern void __stdcall Y90PJF(
  CONST char *line, CONST int line_len,
  CONST int *lvar
);

extern void __stdcall Y90PKF(
  CONST char *line, CONST int line_len,
  CONST double *varr
);

extern void __stdcall Y90PLF(
  CONST char *line, CONST int line_len,
  CONST int *nvecc,
  CONST Complex vecc[],
  CONST int *ivecc
);

extern void __stdcall Y90PMF(
  CONST char *line, CONST int line_len,
  CONST char *vname1, CONST int vname1_len,
  CONST char *vname2, CONST int vname2_len,
  CONST int *nvecc,
  CONST Complex vecc1[],
  CONST int *ivecc1,
  CONST Complex vecc2[],
  CONST int *ivecc2
);

extern void __stdcall Y90PNF(
  CONST char *line, CONST int line_len,
  CONST int *nveci,
  CONST int veci[],
  CONST int *iveci
);

extern void __stdcall Y90PPF(
  CONST char *line, CONST int line_len,
  CONST char *vname1, CONST int vname1_len,
  CONST char *vname2, CONST int vname2_len,
  CONST int *nveci,
  CONST int veci1[],
  CONST int *iveci1,
  CONST int veci2[],
  CONST int *iveci2
);

extern void __stdcall Y90PRF(
  CONST char *line, CONST int line_len,
  CONST int *nvecr,
  CONST double vecr[],
  CONST int *ivecr
);

extern void __stdcall Y90PSF(
  CONST char *line, CONST int line_len,
  CONST char *vname1, CONST int vname1_len,
  CONST char *vname2, CONST int vname2_len,
  CONST int *nvecr,
  CONST double vecr1[],
  CONST int *ivecr1,
  CONST double vecr2[],
  CONST int *ivecr2
);

extern void __stdcall Y90PTF(
  CONST char *line, CONST int line_len,
  CONST int *nvecl,
  CONST int vecl[],
  CONST int *ivecl
);

extern void __stdcall Y90PUF(
  CONST char *line, CONST int line_len,
  CONST char *vname1, CONST int vname1_len,
  CONST char *vname2, CONST int vname2_len,
  CONST int *nvecl,
  CONST int vecl1[],
  CONST int *ivecl1,
  CONST int vecl2[],
  CONST int *ivecl2
);

extern void __stdcall Y90PVF(
  CONST char *matrix, CONST int matrix_len,
  CONST char *line, CONST int line_len,
  CONST int *m,
  CONST int *n,
  CONST int a[] /* 2 dimension */,
  CONST int *ia
);

extern void __stdcall Y90PWF(
  CONST char *matrix, CONST int matrix_len,
  CONST char *line, CONST int line_len,
  CONST char *mname1, CONST int mname1_len,
  CONST char *mname2, CONST int mname2_len,
  CONST int *m,
  CONST int *n,
  CONST int a[] /* 2 dimension */,
  CONST int *ia,
  CONST int b[] /* 2 dimension */,
  CONST int *ib
);

extern void __stdcall Y90PXF(
  CONST char *line, CONST int line_len,
  CONST int *nvecch,
  CONST char *vecch, CONST int vecch_len,
  CONST int *ivecch
);

extern void __stdcall Y90PYF(
  CONST char *line, CONST int line_len,
  CONST char *vname1, CONST int vname1_len,
  CONST char *vname2, CONST int vname2_len,
  CONST int *nvecx,
  CONST char *vecx1, CONST int vecx1_len,
  CONST int *ivecx1,
  CONST char *vecx2, CONST int vecx2_len,
  CONST int *ivecx2
);

extern void __stdcall Y90RAF(
  CONST char *matrix, CONST int matrix_len,
  CONST int *dtype,
  CONST int *ttype,
  CONST int *m,
  CONST int *n,
  CONST int *kl,
  CONST int *ku,
  double a[] /* 2 dimension */,
  CONST int *ia,
  double d[],
  CONST double diag[],
  CONST double odiag[],
  CONST double *cond,
  CONST double *scale,
  double *detman,
  int *detexp,
  CONST int *dist,
  int seed[]
);

extern void __stdcall Y90RBF(
  CONST int *m,
  CONST int *n,
  CONST int *kl,
  CONST int *ku,
  int seed[],
  CONST double d[],
  double a[] /* 2 dimension */,
  CONST int *lda
);

extern void __stdcall Y90RCF(
  CONST int *type,
  CONST int *dtype,
  CONST int *m,
  CONST int *n,
  double a[] /* 2 dimension */,
  CONST int *ia,
  double d[],
  CONST double diag[],
  CONST double *cond,
  CONST double *scale,
  double *detman,
  int *detexp,
  CONST int *dist,
  int seed[],
  CONST char *pvrow, CONST int pvrow_len,
  CONST char *pvcol, CONST int pvcol_len,
  CONST int ipvrow[],
  CONST int ipvcol[],
  double u[] /* 2 dimension */,
  CONST int *iu,
  double v[] /* 2 dimension */,
  CONST int *iv,
  double work1[] /* 2 dimension */,
  CONST int *iwork1,
  double work2[] /* 2 dimension */,
  CONST int *iwork2
);

extern void __stdcall Y90RDF(
  CONST char *side, CONST int side_len,
  CONST char *init, CONST int init_len,
  CONST int *m,
  CONST int *n,
  double a[] /* 2 dimension */,
  CONST int *lda,
  double x[],
  double d[],
  int seed[]
);

extern void __stdcall Y90REF(
  CONST char *uplo, CONST int uplo_len,
  CONST int *n,
  CONST int *kb,
  int seed[],
  CONST double d[],
  double a[] /* 2 dimension */,
  CONST int *lda
);

extern void __stdcall Y90RFF(
  CONST int *ttype,
  CONST char *uplo, CONST int uplo_len,
  CONST int *m,
  CONST int *n,
  double a[] /* 2 dimension */,
  CONST int *ia,
  CONST double diag[],
  CONST double odiag[],
  double *detman,
  int *detexp,
  CONST int *dist,
  int seed[]
);

extern void __stdcall Y90RGF(
  CONST char *det, CONST int det_len,
  CONST int *vtype,
  CONST int *n,
  double v[],
  CONST double vbound[],
  CONST double *cond,
  CONST double *scale,
  double *detman,
  int *detexp,
  CONST int *dist,
  int seed[]
);

extern void __stdcall Y90RHF(
  CONST int *dtype,
  CONST int *type,
  CONST int *n,
  CONST int *kl,
  int nrow[],
  double a[],
  double d[],
  CONST double diag[],
  CONST double odiag[],
  CONST double *cond,
  CONST double *scale,
  double *detman,
  int *detexp,
  CONST int *dist,
  int seed[],
  int irow[],
  double work1[],
  double work2[] /* 2 dimension */,
  CONST int *iwork2
);

extern void __stdcall Y90RJF(
  CONST int *n,
  CONST int *k,
  CONST int *l,
  CONST int *dtype,
  CONST double *condv,
  CONST double diag[],
  CONST double rrx[],
  CONST double rix[],
  CONST int *ncj,
  CONST int icj[],
  double a[] /* 2 dimension */,
  CONST int *ia,
  double vrx[] /* 2 dimension */,
  CONST int *ivrx,
  double vrix[] /* 2 dimension */,
  CONST int *ivrix,
  int seed[],
  double x[],
  double wk1[] /* 2 dimension */,
  CONST int *iwk1,
  double wk2[] /* 2 dimension */,
  CONST int *iwk2,
  double wk3[] /* 2 dimension */,
  CONST int *iwk3
);

extern void __stdcall Y90RKF(
  CONST int *type,
  CONST int *dtype,
  CONST int *n,
  int *nblock,
  int block[] /* 2 dimension */,
  double a[] /* 2 dimension */,
  CONST int *ia,
  double d[],
  CONST double diag[],
  CONST double *cond,
  CONST double *scale,
  double *detman,
  int *detexp,
  CONST int *dist,
  int seed[],
  CONST int iwk1[],
  CONST int iwk2[],
  double rwk1[] /* 2 dimension */,
  CONST int *irwk1,
  double rwk2[] /* 2 dimension */,
  CONST int *irwk2,
  double rwk3[] /* 2 dimension */,
  CONST int *irwk3,
  double rwk4[] /* 2 dimension */,
  CONST int *irwk4
);

extern void __stdcall Y90RLF(
  CONST int *n,
  int *nblock,
  int block[] /* 2 dimension */,
  int seed[]
);

extern void __stdcall Y90RMF(
  CONST char *det, CONST int det_len,
  CONST char *perm, CONST int perm_len,
  CONST char *bstruc, CONST int bstruc_len,
  CONST char *bdiag, CONST int bdiag_len,
  CONST int *type,
  CONST int *m,
  CONST int *n,
  CONST int *dist,
  CONST int *dtype,
  CONST double *cond,
  CONST double *scale,
  double *detman,
  int *detexp,
  double d[],
  CONST double diag[],
  double a[],
  int *na,
  int icol[],
  int irow[],
  double v[],
  double w[],
  int *nblock,
  int block[],
  int seed[]
);

extern void __stdcall Y90RNF(
  CONST int *n,
  CONST double a[] /* 2 dimension */,
  CONST int *ia,
  int *nnvstr,
  int nvstr[],
  int vstr[],
  int stack[],
  int path[],
  int vertex[],
  int number[],
  int lowlnk[],
  int nedge[],
  int loedge[]
);

extern void __stdcall Y90RPF(
  CONST char *uplo, CONST int uplo_len,
  CONST char *diag, CONST int diag_len,
  CONST int *ttype,
  CONST int *n,
  int *kd,
  double a[] /* 2 dimension */,
  CONST int *lda,
  double ab[] /* 2 dimension */,
  CONST int *ldab,
  double d[],
  CONST double *cond,
  CONST double *scale,
  double *detman,
  int *detexp,
  CONST int *dist,
  int seed[]
);

extern void __stdcall Y90RQF(
  CONST int *job,
  CONST int *n,
  double d[],
  double *t0,
  double t[],
  double q[] /* 2 dimension */,
  CONST int *ldq,
  double *err,
  double work[],
  int *ifail
);

extern void __stdcall Y90SAF(
  CONST char *matrix, CONST int matrix_len,
  CONST int *m,
  CONST int *n,
  CONST int *kl,
  CONST int *ku,
  double a[] /* 2 dimension */,
  CONST int *ia
);

extern void __stdcall Y90SBF(
  CONST char *conv, CONST int conv_len,
  CONST char *matrix, CONST int matrix_len,
  CONST int *m,
  CONST int *n,
  CONST int *kl,
  CONST int *ku,
  CONST double a[] /* 2 dimension */,
  CONST int *ia,
  double b[] /* 2 dimension */,
  CONST int *ib
);

extern void __stdcall Y90SCF(
  CONST char *matrix, CONST int matrix_len,
  CONST int *m,
  CONST int *mn,
  CONST int *n,
  CONST int *kl,
  CONST int *ku,
  CONST double a[] /* 2 dimension */,
  CONST int *lda,
  CONST double b[] /* 2 dimension */,
  CONST int *ldb,
  double c[] /* 2 dimension */,
  CONST int *ldc
);

extern double __stdcall Y90SDF(
  CONST char *norm, CONST int norm_len,
  CONST char *matrix, CONST int matrix_len,
  CONST int *m,
  CONST int *n,
  CONST int *kl,
  CONST int *ku,
  CONST double a[] /* 2 dimension */,
  CONST int *lda
);

extern void __stdcall Y90SEF(
  CONST int *m,
  CONST int *mn,
  CONST int *n,
  CONST int *kl,
  CONST int *ku,
  CONST double l[] /* 2 dimension */,
  CONST int *il,
  CONST double u[] /* 2 dimension */,
  CONST int *iu,
  double a[] /* 2 dimension */,
  CONST int *ia
);

extern void __stdcall Y90SFF(
  CONST int *m,
  CONST int *n,
  CONST int *kl,
  CONST int *ku,
  CONST double a[] /* 2 dimension */,
  CONST int *ia,
  double b[] /* 2 dimension */,
  CONST int *ib
);

extern void __stdcall Y90SGF(
  CONST char *type, CONST int type_len,
  CONST int *nlose,
  double a[] /* 2 dimension */,
  CONST int *lda,
  CONST int *m,
  CONST int *n
);

extern void __stdcall Y90SHF(
  CONST char *matrix, CONST int matrix_len,
  CONST int *m,
  CONST int *n,
  CONST double a[] /* 2 dimension */,
  CONST int *nda,
  double b[] /* 2 dimension */,
  CONST int *ndb
);

extern void __stdcall Y90SKF(
  CONST char *matrix, CONST int matrix_len,
  double a[] /* 2 dimension */,
  CONST int *nda,
  CONST int *n
);

extern void __stdcall Y90SLF(
  CONST char *matrix, CONST int matrix_len,
  CONST char *prec, CONST int prec_len,
  CONST int *n,
  CONST int *m,
  CONST double a[] /* 2 dimension */,
  CONST int *ia,
  CONST double b[] /* 2 dimension */,
  CONST int *ib,
  CONST double x[] /* 2 dimension */,
  CONST int *ix,
  double r[] /* 2 dimension */,
  CONST int *ir
);

extern void __stdcall Y90SMF(
  double *vman,
  int *vexp,
  CONST int *scale
);

extern void __stdcall Y90SNF(
  CONST char *matra, CONST int matra_len,
  CONST char *matrb, CONST int matrb_len,
  CONST int *m,
  CONST int *n,
  CONST int *mn,
  CONST double *alpha,
  CONST double a[] /* 2 dimension */,
  CONST int *ia,
  CONST double b[] /* 2 dimension */,
  CONST int *ib,
  CONST double *beta,
  double c[] /* 2 dimension */,
  CONST int *ic
);

extern void __stdcall Y90SPF(
  CONST int *select,
  CONST int *m,
  CONST int *n,
  CONST int irow[],
  CONST double l[],
  CONST double d[],
  CONST double a[] /* 2 dimension */,
  CONST int *ia,
  double b[] /* 2 dimension */,
  CONST int *ib,
  double work[] /* 2 dimension */,
  CONST int *iwork
);

extern void __stdcall Y90SQF(
  CONST int *n,
  CONST int irow[],
  CONST double l[],
  CONST double d[],
  double a[]
);

extern double __stdcall Y90SRF(
  CONST char *norm, CONST int norm_len,
  CONST int *n,
  CONST int irow[],
  double a[],
  double work[]
);

extern double __stdcall Y90TAF(
  int seed[]
);

extern double __stdcall Y90TBF(
  CONST int *dist,
  int seed[]
);

extern void __stdcall Y90TCF(
  CONST char *order, CONST int order_len,
  double v[],
  CONST int *n
);

extern void __stdcall Y90TDF(
  CONST int *n,
  CONST int *k,
  CONST int *l,
  CONST double d[],
  CONST int intger[],
  CONST double a[] /* 2 dimension */,
  CONST int *ia,
  double w[] /* 2 dimension */,
  CONST int *iw
);

extern void __stdcall Y90TEF(
  CONST int *n,
  CONST int *k,
  CONST int *l,
  CONST double d[],
  CONST int intger[],
  CONST double a[] /* 2 dimension */,
  CONST int *ia,
  double w1[] /* 2 dimension */,
  CONST int *iw1,
  double w2[] /* 2 dimension */,
  CONST int *iw2
);

extern void __stdcall Y90TFF(
  CONST int *n,
  CONST double rr[],
  CONST double ri[],
  int iord1[],
  int iord2[]
);

extern int __stdcall Y90WAF(
  CONST char *ca, CONST int ca_len,
  CONST char *cb, int cb_len
);

extern void __stdcall Y90ZFF(
  CONST char *prtype, CONST int prtype_len,
  CONST char *matype, CONST int matype_len,
  CONST char *select, CONST int select_len,
  int *n,
  double *x
);

#ifdef __cplusplus
}
#endif
#endif
