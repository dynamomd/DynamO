/*   quartic.c
    version 54

    a test of relative accuracy of various quartic routines.

    use:
      quartic [-a] [-c n] [-q n] [-d n]
      (no parameters)  : a loop through 10,000 prechosen quartic coefficients
      -a               : improve cubic roots (default: no iteration)
      -c  n            : keep reading n cubic roots from standard input
                         (if n=0 read coefficients), terminate with 'q'.
      -d  n            : set debug value to 'n' (default: 5)
      -q  n            : keep reading n quartic roots from standard input
                         (if n=0 read coefficients), terminate with 'q'.

    debug <= 1 : print cases used

    method-
     http://linus.socs.uts.edu.au/~don/pubs/solving.html

     4 Apr 2004 bringing yacfraid notation into line with solving.html
     4 Apr 2004 fixed bug in final root calculation in yacfraid
     8 Mar 2004 corrected modified yacfraid algorithm
     8 Mar 2004 modified yacfraid algorithm
    31 Jan 2004 printing results when number of roots vary
    30 Jan 2004 printing error of chosen cubic if debug<1
    27 Jan 2004 seeking worst coefficients for each algorithm
    27 Jan 2004 choosing best route for k in chris
    26 Jan 2004 conforming variables to solving.html
    26 Jan 2004 fixed bug in yacfraid for multiplicity 3
    25 Jan 2004 speeding up chris
    24 Jan 2004 fixed chris error (A <-> B)
    23 Jan 2004 use debug<1 for diagnostics
    21 Jan 2004 solve cubic in neumark, yacfraid and chris if d=0
    20 Jan 2004 3 roots to cubic in chris if e1=0
    19 Jan 2004 debugging Christianson routine
    14 Jan 2004 adding Christianson, cut determinant in quadratic call
     5 Jan 2004 improving cubic roots by Newton-Raphson iteration
    23 Dec 2003 seeking snag in yacfraid
    21 Dec 2003 putting diagnostic printout in yacfraid
    18 Dec 2003 allowing input of equation coefficients
    18 Dec 2003 recording cubic root giving least worst error
    17 Dec 2003 recording cubic root giving the most quartic roots
    16 Dec 2003 using the cubic root giving the most quartic roots
    16 Dec 2003 trying consistency of 3 cubic roots where available
    15 Dec 2003 trying all 3 cubic roots where available
    13 Dec 2003 trying to fix Neumark
    13 Dec 2003 cleaning up diagnostic format
    12 Dec 2003 initialising n,m,po3 in cubic
    12 Dec 2003 allow cubic to return 3 zero roots if p=q=r=0
    10 Dec 2003 added setargs
     2 Dec 2003 finding worst cases
     2 Dec 2003 negating j if p>0 in cubic
     1 Dec 2003 changing v in cubic from (sinsqk > doub0) to (sinsqk >= doub0)
     1 Dec 2003 test changing v in cubic from po3sq+po3sq to doub2*po3sq
    30 Nov 2003 counting cases in all solving routines
    29 Nov 2003 testing wsq >= doub0
    29 Nov 2003 mult by doub2 for v in cubic
    29 Nov 2003 cut po3cu from cubic
    29 Nov 2003 better quadratic
    29 Nov 2003 count agreements
    17 Nov 2003 count special cases
    16 Nov 2003 option of loop or read coefficients from input
    15 Nov 2003 fixed cubic() bug
    11 Nov 2003 added Brown and Yacoub & Fraidenraich's algorithms
    21 Jan 1989 quartic selecting Descartes, Ferrari, or Neumark algorithms
    16 Jul 1981 Don Herbison-Evans

"Solving Quartics and Cubics for Graphics", D. Herbison-Evans,
Graphics Gems V (ed.: A. Paeth) Academic Press (Chesnut Hill), pp. 3-15 (1995).

"Solving Quartics and Cubics for Graphics", D. Herbison-Evans,
Research Report CS-86-56, Department of Computer Science, University of Waterloo
(1986)

"Caterpillars and the Inaccurate Solution of Cubic and Quartic Equations",
D. Herbison-Evans, Australian Computer Science Communications,
Vol. 5, No. 1, pp. 80-90 (1983)

  subroutines:
     setcns      - set constants
     setargs     - determine which equations to test
     looptest    - test 10000 coefficient combinations
     cases       - print statistics
     docoeffs    - read coefficients of a particular equation to solve
     cubictest   - test cubic solutions
     quartictest - test quartic solutions
     compare     - check accuracy of results
     errors      - find errors in a set of roots
     acos3       - find arcos(x/3)
     curoot      - find cube root
     quadratic   - solve a quadratic
     cubic       - solve a cubic
     cubnewton   - improve cubic roots by iteration
     quanewton   - improve quatric roots by iteration
     quartic     - solve a quartic
     descartes   - use Descartes' algorithm to solve a quartic
     ferrari     - use Ferrari's algorithm to solve a quartic
     neumark     - use Neumark's algorithm to solve a quartic
     yacfraid    - use Yacoub & Fraidenraich's algorithm to solve a quartic
     chris       - use Christianson's algorithm to solve a quartic
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <iostream>
#include <limits>

#define TRUE 1
#define FALSE 0
#define NCASES 40

double a, b, c, d; /* quartic coefficients */
double cc[4], cd[4], cf[4], cn[4], cy[4];
double d3o8, d3o256; /* double precision constants */
double doub0;
double doub1, doub2;
double doub3, doub4;
double doub5, doub6;
double doub8, doub9, doub12;
double doub16, doub24;
double doub27, doub64;
double doubmax; /* approx square root of max double number */
double doubmin; /* smallest double number */
double doubtol; /* tolerance of double numbers */
double errq;
double inv2, inv3, inv4;
double inv8, inv16, inv32, inv64, inv128;
double maxc, maxd, maxf, maxn, maxy;
double p, q, r;    /* cubic coefficients */
double qrts[4][3]; /* quartic roots for each cubic root */
double rt3;        /* square root of 3 */
double rtsc[4], rtsd[4], rtsf[4], rtsn[4], rtsq[4],
    rtsy[4]; /* roots given by each algorithm */
double rterc[4], rterd[4], rterf[4], rtern[4], rterq[4],
    rtery[4]; /* errors of roots */
double worstc, worstd, worstf, worstn,
    worsty;       /* worst error for each algotithm */
double worst3[3]; /* worst error for a given cubic root */
double x0, x1, x2, x3;

int agr, dis, zer;
int agreecd, agreecf, agreecn, agreecy;
int agreedf, agreedn, agreedy;
int agreefn, agreefy, agreeny;
int debug = 10000; /* <1 for lots of diagnostics, >5 for none */
int docubic, doquartic;
int iterate;
int j3;
int n1, n2, n3, n4[3];
int n, nc, nd, nf, nn, ny;
int nrtsc, nrtsd, nrtsf, nrtsn, nrtsq, nrtsy;
int nt3, nt2, nt2sq, nden, nm4, ndet1, ndet2, ndet3, ndet4, ndet5;
int nqud[NCASES];
int ncub[NCASES];
int nchr[NCASES];
int ndes[NCASES];
int nfer[NCASES];
int nneu[NCASES];
int nyac[NCASES];
int nzc, nzd, nzf, nzn, nzy;
int zerodf, zerodn, zerody, zerofn, zerofy, zerony;
int zerodfn, zerodfy, zerodny, zerofny;
int tot;

double acos3(double);
int cubic(double p, double q, double r, double v3[4]);
double curoot(double x);
int quadratic(double b, double c, double rts[4]);
void cubnewton(double p, double q, double r, int n3, double v3[4]);
int descartes(double a, double b, double c, double d, double rts[4]);
int chris(double a, double b, double c, double d, double rts[4]);
int ferrari(double a, double b, double c, double d, double rts[4]);
int neumark(double a, double b, double c, double d, double rts[4]);
int yacfraid(double a, double b, double c, double d, double rts[4]);

void compare(), cases(), cubictest(), quartictest();

double errors(double a, double b, double c, double d, double rts[4],
              double rterr[4], int nrts);
void quanewton(double p, double q, double r, double s, int n4, double v4[4]);

void setcns()
/*
    set up constants.

    called by main, quartictest.
*/
{
  int j;

  doub0 = (double)0;
  doub1 = (double)1;
  doub2 = (double)2;
  doub3 = (double)3;
  doub4 = (double)4;
  doub5 = (double)5;
  doub6 = (double)6;
  doub8 = (double)8;
  doub9 = (double)9;
  doub12 = (double)12;
  doub16 = (double)16;
  doub24 = (double)24;
  doub27 = (double)27;
  doub64 = (double)64;
  inv2 = doub1 / doub2;
  inv3 = doub1 / doub3;
  inv4 = doub1 / doub4;
  inv8 = doub1 / doub8;
  inv16 = doub1 / doub16;
  inv32 = doub1 / (double)32;
  inv64 = doub1 / doub64;
  inv128 = doub1 / (double)128;
  d3o8 = doub3 / doub8;
  d3o256 = doub3 / (double)256;
  rt3 = sqrt(doub3);

  doubtol = doub1;
  for (j = 1; doub1 + doubtol > doub1; ++j) {
    doubtol *= inv2;
  }
  doubtol = sqrt(doubtol);

  doubmin = inv2;
  for (j = 1; j <= 100; ++j) {
    doubmin = doubmin * doubmin;
    if ((doubmin * doubmin) <= (doubmin * doubmin * inv2))
      break;
  }
  // doubmax=doub1/sqrt(doub2*doubmin);
  doubmax = sqrt(std::numeric_limits<double>::max());

} /* setcns */
/***********************************/

void setargs(int argc, char *argv[])
/*
  determine whether to do looptest or just individual equation

  10 Dec 2003 Don Herbison-Evans

  called by main
*/
{
  int j;

  docubic = FALSE;
  doquartic = FALSE;
  iterate = FALSE;
  for (j = 1; j < argc; ++j) {
    if (argv[j][0] == '-') {
      if (argv[j][1] == 'a') {
        iterate = TRUE;
        printf("iterate cubics\n");
      } else if (argv[j][1] == 'c') {
        docubic = TRUE;
        sscanf(argv[++j], "%d", &n);
        printf("do cubic %d\n", n);
      } else if (argv[j][1] == 'q') {
        doquartic = TRUE;
        sscanf(argv[++j], "%d", &n);
        printf("do quartic %d\n", n);
      } else if (argv[j][1] == 'd') {
        sscanf(argv[++j], "%d", &debug);
        printf("debug %d\n", debug);
      }
    }
  }
} /* setargs */
/********************************************************/

void looptest()
/*
  loop through a range of coefficient values
  comparing the algorithms

  called by main.
  calls compare.
*/
{
  double ten4, ten8;
  double v[10];
  int j, k, l, m;

  printf("quartic loop test\n");
  nc = 0;
  nd = 0;
  nf = 0;
  nn = 0;
  ny = 0;
  maxc = doub0;
  maxd = doub0;
  maxf = doub0;
  maxn = doub0;
  maxy = doub0;
  tot = 0;
  zer = 0;
  agr = 0;
  dis = 0;
  agreecd = 0;
  agreecf = 0;
  agreecn = 0;
  agreecy = 0;
  agreedf = 0;
  agreedn = 0;
  agreedy = 0;
  agreefn = 0;
  agreefy = 0;
  agreeny = 0;
  for (j = 0; j < NCASES; ++j) {
    nqud[j] = 0;
    ncub[j] = 0;
    ndes[j] = 0;
    nfer[j] = 0;
    nneu[j] = 0;
    nyac[j] = 0;
  }
  ten4 = (double)10000;
  ten8 = ten4 * ten4;
  j = 0;
  for (a = ten8; a > inv2 / ten8; a /= ten4) {
    v[j] = a;
    ++j;
  }
  for (j = 0; j < 5; ++j)
    v[j + 5] = -v[j];
  for (j = 0; j < 10; ++j) {
    a = v[j];
    for (k = 0; k < 10; ++k) {
      b = v[k];
      for (l = 0; l < 10; ++l) {
        c = v[l];
        for (m = 0; m < 10; ++m) {
          d = v[m];
          compare();
        } /* for m */
      } /* for l */
    } /* for k */
  } /* for j */
  printf("total cases: %d\n", tot);
  printf("number of real roots: five agree %d, disagree %d\n", agr, dis);
  printf("two agree: cd %d, cf %d, cn %d, cy %d\n", agreecd, agreecf, agreecn,
         agreecy);
  printf("           df %d, dn %d, dy %d, fn %d, fy %d, ny %d\n", agreedf,
         agreedn, agreedy, agreefn, agreefy, agreeny);
  printf("5 agree on no real roots: %d\n", zer);
  printf("no real roots: chris %d, desc %d, ferr %d, neum %d, yacf %d\n", nzc,
         nzd, nzf, nzn, nzy);
  printf("Christianson : best %d, worst %g %g %g %g, error %g\n", nc, cc[0],
         cc[1], cc[2], cc[3], maxc);
  printf("Descartes : best %d, worst %g %g %g %g, error %g\n", nd, cd[0], cd[1],
         cd[2], cd[3], maxd);
  printf("Ferrari : best %d, worst %g %g %g %g, error %g\n", nf, cf[0], cf[1],
         cf[2], cf[3], maxf);
  printf("Neumark : best %d, worst %g %g %g %g, error %g\n", nn, cn[0], cn[1],
         cn[2], cn[3], maxn);
  printf("Yacoub  : best %d, worst %g %g %g %g, error %g\n", ny, cy[0], cy[1],
         cy[2], cy[3], maxy);
  cases();
} /* looptest */
/********************************************/

void cases()
/*
  print statistics

   2 Dec 2003 Don Herbison-Evans

  called by looptest, quartictest.
*/
{
  int j;

  printf("quadratic cases:\n");
  for (j = 0; j < NCASES; ++j) {
    printf("  %2d %d,", j, nqud[j]);
    if (j % 5 == 4)
      printf("\n");
  }
  printf("\n");

  printf("cubic cases:\n");
  for (j = 0; j < NCASES; ++j) {
    printf("  %2d %d,", j, ncub[j]);
    if (j % 5 == 4)
      printf("\n");
  }
  printf("\n");

  printf("descartes cases:\n");
  for (j = 0; j < NCASES; ++j) {
    printf("  %2d %d,", j, ndes[j]);
    if (j % 5 == 4)
      printf("\n");
  }
  printf("\n");

  printf("ferrari cases:\n");
  for (j = 0; j < NCASES; ++j) {
    printf("  %2d %d,", j, nfer[j]);
    if (j % 5 == 4)
      printf("\n");
  }
  printf("\n");

  printf("neumark cases:\n");
  for (j = 0; j < NCASES; ++j) {
    printf("  %2d %d,", j, nneu[j]);
    if (j % 5 == 4)
      printf("\n");
  }
  printf("\n");

  printf("yacfraid cases:\n");
  for (j = 0; j < NCASES; ++j) {
    printf("  %2d %d,", j, nyac[j]);
    if (j % 5 == 4)
      printf("\n");
  }
  printf("\n");

  printf("chris cases:\n");
  for (j = 0; j < NCASES; ++j) {
    printf("  %2d %d,", j, nchr[j]);
    if (j % 5 == 4)
      printf("\n");
  }
  printf("\n");

} /* cases */
/************************************/

void docoeffs()
/*
   read a set of coefficients from standard input
   and compare results of the different algorithms.

   16 Nov 2003 Don Herbison-Evans

   called by main.
   calls cubictest, quartictest.
*/
{
  // float fa,fb,fc,fd;
  float fp, fq, fr;

  printf("docoeffs %d\n", n);
  x0 = doub0;
  x1 = doub0;
  x2 = doub0;
  x3 = doub0;
  {
    if (n == 0) {
      x0 = doub1;
      if (docubic == TRUE) {
        std::cin >> fp >> fq >> fr; // sscanf(buf,"%g %g %g",&fp,&fq,&fr);
        p = (double)fp;
        q = (double)fq;
        r = (double)fr;
      } else {
        // sscanf(buf,"%g %g %g %g",&fa,&fb,&fc,&fd);
        std::cin >> a >> b >> c >> d;
        // sscanf(buf,"%g %g %g %g",&a,&b,&c,&d);
        // printf("\n%g %g %g %g\n",fa,fb,fc,fd);
        // a = (double)fa;
        // b = (double)fb;
        // c = (double)fc;
        // d = (double)fd;
        printf("\n%g %g %g %g\n", a, b, c, d);
      }
    } else if (n == 1) {
      // sscanf(buf,"%d\n",&ix0);
      // x0 = (double)ix0;
      std::cin >> x0;
    } else if (n == 2) {
      // sscanf(buf,"%d %d\n",&ix0,&ix1);
      // x0 = (double)ix0;
      // x1 = (double)ix1;
      std::cin >> x0 >> x1;
    } else if (n == 3) {
      // sscanf(buf,"%d %d %d\n",&ix0,&ix1,&ix2);
      // x0 = (double)ix0;
      // x1 = (double)ix1;
      // x2 = (double)ix2;
      std::cin >> x0 >> x1 >> x2;
    } else if (n == 4) {
      // sscanf(buf,"%d %d %d %d\n",&ix0,&ix1,&ix2,&ix3);
      // x0 = (double)ix0;
      // x1 = (double)ix1;
      // x2 = (double)ix2;
      // x3 = (double)ix3;
      std::cin >> x0 >> x1 >> x2 >> x3;
    }
  }
  if ((fabs(x0) + fabs(x1) + fabs(x2) + fabs(x3)) != doub0) {
    if (docubic == TRUE)
      cubictest();
    if (doquartic == TRUE)
      quartictest();
  }
} /* docoeffs */
/*****************************************/

void cubictest()
/*
  called by docoeffs.
*/
{
  int j;
  int nrtsc;
  double rtsc[4];

  if (n == 1) {
    p = -x0;
    q = doub1;
    r = -x0;
  } else if (n == 3) {
    p = -(x0 + x1 + x2);
    q = x0 * x1 + x0 * x2 + x1 * x2;
    r = -x0 * x1 * x2;
  }
  printf("\ncubic test: %g %g %g\n", x0, x1, x2);
  printf(" x^3 + %gx^2 + %gx + %g\n", p, q, r);
  nrtsc = cubic(p, q, r, rtsc);
  printf("%d roots\n", nrtsc);
  for (j = 0; j < nrtsc; ++j)
    printf("%g  ", rtsc[j]);
  printf("\n");
} /* cubictest */
/******************************************/

void quartictest()
/*
  called by docoeffs.
  calls descartes, ferrari, neumark, yacfraid,
        errors, cases, setcns.
*/
{
  int j;

  printf("quartictest %d\n", n);
  if (n == 4) {
    a = -(x0 + x1 + x2 + x3);
    b = x0 * x1 + x0 * x2 + x0 * x3 + x1 * x2 + x1 * x3 + x2 * x3;
    c = -(x0 * x1 * x2 + x0 * x1 * x3 + x0 * x2 * x3 + x1 * x2 * x3);
    d = x0 * x1 * x2 * x3;
  } else if (n == 2) {
    a = -(x0 + x1);
    b = x0 * x1;
    c = a;
    d = b;
  }
  printf("\nquartic test: %g %g %g %g\n", x0, x1, x2, x3);
  printf("x^4 + %gx^3 + %gx^2 + %gx + %g\n", a, b, c, d);

  setcns();
  printf("\nDescartes\n");
  nrtsd = descartes(a, b, c, d, rtsd);
  worstd = errors(a, b, c, d, rtsd, rterd, nrtsd);
  if (nrtsd > 0)
    for (j = 0; j < nrtsd; ++j)
      printf(" %g %g\n", rtsd[j], rterd[j]);
  else
    printf("no real roots found\n");

  setcns();
  printf("\nFerrari\n");
  nrtsf = ferrari(a, b, c, d, rtsf);
  worstf = errors(a, b, c, d, rtsf, rterf, nrtsf);
  if (nrtsf > 0)
    for (j = 0; j < nrtsf; ++j)
      printf(" %g %g\n", rtsf[j], rterf[j]);
  else
    printf("no real roots found\n");

  setcns();
  printf("\nNeumark\n");
  nrtsn = neumark(a, b, c, d, rtsn);
  worstn = errors(a, b, c, d, rtsn, rtern, nrtsn);
  if (nrtsn > 0)
    for (j = 0; j < nrtsn; ++j)
      printf(" %g %g\n", rtsn[j], rtern[j]);
  else
    printf("no real roots found\n");

  setcns();
  printf("\nYacoub and Fraidenraich\n");
  nrtsy = yacfraid(a, b, c, d, rtsy);
  worsty = errors(a, b, c, d, rtsy, rtery, nrtsy);
  if (nrtsy > 0)
    for (j = 0; j < nrtsy; ++j)
      printf(" %g %g\n", rtsy[j], rtery[j]);
  else
    printf("no real roots found\n");

  setcns();
  printf("\nChristianson\n");
  nrtsc = chris(a, b, c, d, rtsc);
  worstc = errors(a, b, c, d, rtsc, rterc, nrtsc);
  if (nrtsc > 0)
    for (j = 0; j < nrtsc; ++j)
      printf(" %g %g\n", rtsc[j], rterc[j]);
  else
    printf("no real roots found\n");
  printf("\n x^4 + %gx^3 + %gx^2 + %gx + %g\n", a, b, c, d);
  if (debug < 2)
    cases();
} /* quartictest */
/***************************************/

void compare()
/*
   called by looptest.
   calls errors, descartes, ferrari, neumark, yacfraid, chris.
*/
{
  ++tot;
  nrtsd = descartes(a, b, c, d, rtsd);
  if (iterate == TRUE)
    quanewton(a, b, c, d, nrtsd, rtsd);
  worstd = errors(a, b, c, d, rtsd, rterd, nrtsd);
  nrtsf = ferrari(a, b, c, d, rtsf);
  if (iterate == TRUE)
    quanewton(a, b, c, d, nrtsf, rtsf);
  worstf = errors(a, b, c, d, rtsf, rterf, nrtsf);
  nrtsn = neumark(a, b, c, d, rtsn);
  if (iterate == TRUE)
    quanewton(a, b, c, d, nrtsn, rtsn);
  worstn = errors(a, b, c, d, rtsn, rtern, nrtsn);
  nrtsy = yacfraid(a, b, c, d, rtsy);
  if (iterate == TRUE)
    quanewton(a, b, c, d, nrtsy, rtsy);
  worsty = errors(a, b, c, d, rtsy, rtery, nrtsy);
  nrtsc = chris(a, b, c, d, rtsc);
  if (iterate == TRUE)
    quanewton(a, b, c, d, nrtsc, rtsc);
  worstc = errors(a, b, c, d, rtsc, rterc, nrtsc);

  if (nrtsd == 0)
    ++nzd;
  if (nrtsf == 0)
    ++nzf;
  if (nrtsn == 0)
    ++nzn;
  if (nrtsy == 0)
    ++nzy;
  if (nrtsc == 0)
    ++nzc;

  if (nrtsc == nrtsd)
    ++agreecd;
  if (nrtsc == nrtsf)
    ++agreecf;
  if (nrtsc == nrtsn)
    ++agreecn;
  if (nrtsc == nrtsy)
    ++agreecy;
  if (nrtsd == nrtsf)
    ++agreedf;
  if (nrtsd == nrtsn)
    ++agreedn;
  if (nrtsd == nrtsy)
    ++agreedy;
  if (nrtsf == nrtsn)
    ++agreefn;
  if (nrtsf == nrtsy)
    ++agreefy;
  if (nrtsn == nrtsy)
    ++agreeny;

  if ((nrtsc == nrtsd) && (nrtsc == nrtsf) && (nrtsc == nrtsn) &&
      (nrtsc == nrtsy)) {
    if (nrtsf == 0)
      ++zer;
    ++agr;
  } else
    ++dis;

  if ((nrtsc == nrtsd) && (nrtsc == nrtsf) && (nrtsc == nrtsn) &&
      (nrtsc == nrtsy) && (nrtsc != 0)) {
    if ((worstd < worstf) && (worstd < worstn) && (worstd < worsty) &&
        (worstd < worstc))
      ++nd;

    if ((worstf < worstd) && (worstf < worstn) && (worstf < worsty) &&
        (worstf < worstc))
      ++nf;

    if ((worstn < worstd) && (worstn < worstf) && (worstn < worsty) &&
        (worstn < worstc))
      ++nn;

    if ((worsty < worstd) && (worsty < worstf) && (worsty < worstn) &&
        (worsty < worstc))
      ++ny;

    if ((worstc < worstd) && (worstc < worstf) && (worstc < worstn) &&
        (worstc < worsty))
      ++nc;

    if (maxd < worstd) {
      maxd = worstd;
      cd[0] = a;
      cd[1] = b;
      cd[2] = c;
      cd[3] = d;
    }
    if (maxf < worstf) {
      maxf = worstf;
      cf[0] = a;
      cf[1] = b;
      cf[2] = c;
      cf[3] = d;
    }
    if (maxn < worstn) {
      maxn = worstn;
      cn[0] = a;
      cn[1] = b;
      cn[2] = c;
      cn[3] = d;
    }
    if (maxy < worsty) {
      maxy = worsty;
      cy[0] = a;
      cy[1] = b;
      cy[2] = c;
      cy[3] = d;
    }
    if (maxc < worstc) {
      maxc = worstc;
      cc[0] = a;
      cc[1] = b;
      cc[2] = c;
      cc[3] = d;
    }

  } /* nrts agree */
} /* compare */
/****************************************************/

double errors(double a, double b, double c, double d, double rts[4],
              double rterr[4], int nrts)
/*
  find the errors

  called by quartictest, docoeff, compare,
            chris, descartes, ferrari, neumark, yacfraid.
*/
{
  int k;
  double deriv, test, worst;

  worst = doubmax;
  if (nrts > 0) {
    worst = doub0;
    for (k = 0; k < nrts; ++k) {
      test = (((rts[k] + a) * rts[k] + b) * rts[k] + c) * rts[k] + d;
      if (test == doub0)
        rterr[k] = doub0;
      else {
        deriv =
            ((doub4 * rts[k] + doub3 * a) * rts[k] + doub2 * b) * rts[k] + c;
        if (deriv != doub0)
          rterr[k] = fabs(test / deriv);
        else {
          deriv = (doub12 * rts[k] + doub6 * a) * rts[k] + doub2 * b;
          if (deriv != doub0)
            rterr[k] = sqrt(fabs(test / deriv));
          else {
            deriv = doub24 * rts[k] + doub6 * a;
            if (deriv != doub0)
              rterr[k] = curoot(fabs(test / deriv));
            else
              rterr[k] = sqrt(sqrt(fabs(test) / doub24));
          }
        }
      }
      if (rts[k] != doub0)
        rterr[k] /= rts[k];
      if (rterr[k] < doub0)
        rterr[k] = -rterr[k];
      if (rterr[k] > worst)
        worst = rterr[k];
    }
  }
  return (worst);
} /* errors */
/**********************************************/

double acos3(double x)
/*
    find cos(acos(x)/3)

    16 Jul 1981   Don Herbison-Evans

    called by cubic .
*/
{
  double value;

  value = cos(acos(x) * inv3);
  return (value);
} /* acos3 */
/***************************************/

double curoot(double x)
/*
    find cube root of x.

    30 Jan 1989   Don Herbison-Evans

    called by cubic .
*/
{
  double value;
  double absx;
  int neg;

  neg = 0;
  absx = x;
  if (x < doub0) {
    absx = -x;
    neg = 1;
  }
  if (absx != doub0)
    value = exp(log(absx) * inv3);
  else
    value = doub0;
  if (neg == 1)
    value = -value;
  return (value);
} /* curoot */
/****************************************************/

int quadratic(double b, double c, double rts[4])
/*
    solve the quadratic equation -

        x**2 + b*x + c = 0

    14 Jan 2004   cut determinant in quadratic call
    29 Nov 2003   improved
    16 Jul 1981   Don Herbison-Evans

    called by  cubic,quartic,chris,descartes,ferrari,neumark.
*/
{
  int nquad;
  double dis, rtdis;

  dis = b * b - doub4 * c;
  rts[0] = doub0;
  rts[1] = doub0;
  if (b == doub0) {
    if (c == doub0) {
      nquad = 2;
      ++nqud[0];
    } else {
      if (c < doub0) {
        nquad = 2;
        rts[0] = sqrt(-c);
        rts[1] = -rts[0];
        ++nqud[1];
      } else {
        nquad = 0;
        ++nqud[2];
      }
    }
  } else if (c == doub0) {
    nquad = 2;
    rts[0] = -b;
    ++nqud[3];
  } else if (dis >= doub0) {
    nquad = 2;
    rtdis = sqrt(dis);
    if (b > doub0) {
      rts[0] = (-b - rtdis) * inv2;
      ++nqud[4];
    } else {
      rts[0] = (-b + rtdis) * inv2;
      ++nqud[5];
    }
    if (rts[0] == doub0) {
      rts[1] = -b;
      ++nqud[6];
    } else {
      rts[1] = c / rts[0];
      ++nqud[7];
    }
  } else {
    nquad = 0;
    ++nqud[8];
  }
  if (debug < 1) {
    printf("quad  b %g   c %g  dis %g\n", b, c, dis);
    printf("      %d %g %g\n", nquad, rts[0], rts[1]);
  }
  return (nquad);
} /* quadratic */
/**************************************************/

int cubic(double p, double q, double r, double v3[4])
/*
  find the real roots of the cubic -
      x**3 + p*x**2 + q*x + r = 0

    12 Dec 2003 initialising n,m,po3
    12 Dec 2003 allow return of 3 zero roots if p=q=r=0
     2 Dec 2003 negating j if p>0
     1 Dec 2003 changing v from (sinsqk > doub0) to (sinsqk >= doub0)
     1 Dec 2003 test changing v from po3sq+po3sq to doub2*po3sq
    16 Jul 1981 Don Herbison-Evans

  input parameters -
    p,q,r - coeffs of cubic equation.

  output-
    the number of real roots
    v3 - the roots.

  global constants -
    rt3 - sqrt(3)
    inv3 - 1/3
    doubmax - square root of largest number held by machine

    method -
    see D.E. Littlewood, "A University Algebra" pp.173 - 6

    15 Nov 2003 output 3 real roots: Don Herbison-Evans
       Apr 1981 initial version: Charles Prineas

    called by  cubictest,quartic,chris,yacfraid,neumark,descartes,ferrari.
    calls      quadratic,acos3,curoot,cubnewton.
*/
{
  // printf("pqr: %g %g %g\n",p,q,r);
  int j, n3;
  double po3, po3sq, qo3, po3q;
  double uo3, u2o3, uo3sq4, uo3cu4;
  double v, vsq, wsq;
  double m1, m2, mcube;
  double muo3, s, scube, t, cosk, rt3sink, sinsqk;

  m1 = doub0;
  m2 = doub0;
  po3 = doub0;
  v = doub0;
  uo3 = doub0;
  cosk = doub0;
  if (r == doub0) {
    ++ncub[0];
    // printf("cubic root case 0\n");
    n3 = quadratic(p, q, v3);
    v3[n3++] = doub0;
    goto done;
  }
  if ((p == doub0) && (q == doub0)) {
    ++ncub[1];
    // printf("cubic root case 1\n");
    v3[0] = curoot(-r);
    v3[1] = v3[0];
    v3[2] = v3[0];
    n3 = 3;
    goto done;
  }
  n3 = 1;
  if ((p > doubmax) || (p < -doubmax)) {
    v3[0] = -p;
    ++ncub[2];
    // printf("cubic root case 2\n");
    goto done;
  }
  if ((q > doubmax) || (q < -doubmax)) {
    if (q > doub0) {
      v3[0] = -r / q;
      ++ncub[3];
      // printf("cubic root case 3\n");
      goto done;
    } else if (q < doub0) {
      v3[0] = -sqrt(-q);
      ++ncub[4];
      // printf("cubic root case 4\n");
      goto done;
    } else {
      v3[0] = doub0;
      ++ncub[5];
      // printf("cubic root case 5\n");
      goto done;
    }
  } else if ((r > doubmax) || (r < -doubmax)) {
    v3[0] = -curoot(r);
    ++ncub[6];
    // printf("cubic root case 6\n");
    goto done;
  } else {
    po3 = p * inv3;
    po3q = po3 * q;
    po3sq = po3 * po3;
    if (po3sq > doubmax) {
      v3[0] = -p;
      ++ncub[7];
      // printf("cubic root case 7\n");
      goto done;
    } else {
      v = r + po3 * (po3sq + po3sq - q);
      if ((v > doubmax) || (v < -doubmax)) {
        v3[0] = -p;
        ++ncub[8];
        // printf("cubic root case 8\n");
        goto done;
      } else {
        vsq = v * v;
        qo3 = q * inv3;
        uo3 = qo3 - po3sq;
        u2o3 = uo3 + uo3;
        if ((u2o3 > doubmax) || (u2o3 < -doubmax)) {
          if (p == doub0) {
            if (q > doub0) {
              v3[0] = -r / q;
              ++ncub[9];
              // printf("cubic root case 9\n");
              goto done;
            } else if (q < doub0) {
              v3[0] = -sqrt(-q);
              ++ncub[10];
              // printf("cubic root case 10\n");
              goto done;
            } else {
              v3[0] = doub0;
              ++ncub[11];
              // printf("cubic root case 11\n");
              goto done;
            }
          } else {
            v3[0] = -q / p;
            ++ncub[12];
            // printf("cubic root case 12\n");
            goto done;
          }
        }
        uo3sq4 = u2o3 * u2o3;
        if (uo3sq4 > doubmax) {
          if (p == doub0) {
            if (q > doub0) {
              v3[0] = -r / q;
              ++ncub[13];
              // printf("cubic root case 13\n");
              goto done;
            } else if (q < doub0) {
              v3[0] = -sqrt(-q);
              ++ncub[14];
              // printf("cubic root case 14\n");
              goto done;
            } else {
              v3[0] = doub0;
              ++ncub[15];
              // printf("cubic root case 15\n");
              goto done;
            }
          } else {
            v3[0] = -q / p;
            ++ncub[16];
            // printf("cubic root case 16\n");
            goto done;
          }
        }
        uo3cu4 = uo3sq4 * uo3;
        wsq = uo3cu4 + vsq;
        if (wsq > doub0) {
          /*
              cubic has one real root -
          */
          if (v <= doub0) {
            mcube = (-v + sqrt(wsq)) * inv2;
            ++ncub[17];
            // printf("cubic root case 17\n");
          } else {
            mcube = (-v - sqrt(wsq)) * inv2;
            ++ncub[18];
            // printf("cubic root case 18\n");
          }
          m1 = curoot(mcube);
          if (m1 != doub0) {
            m2 = -uo3 / m1;
            ++ncub[19];
            // printf("cubic root case 19\n");
          } else {
            m2 = doub0;
            ++ncub[20];
            // printf("cubic root case 20\n");
          }
          v3[0] = m1 + m2 - po3;
        } else {
          /*
              cubic has three real roots -
          */
          if (uo3 < doub0) {
            muo3 = -uo3;
            if (muo3 > doub0) {
              s = sqrt(muo3);
              ++ncub[21];
              // printf("cubic root case 21\n");
              if (p > doub0) {
                s = -s;
                ++ncub[22];
                // printf("cubic root case 22\n");
              }
            } else {
              s = doub0;
              ++ncub[23];
              // printf("cubic root case 23\n");
            }
            scube = s * muo3;
            if (scube == doub0) {
              v3[0] = m1 + m2 - po3;
              n3 = 1;
              ++ncub[24];
              // printf("cubic root case 24\n");
            } else {
              t = -v / (scube + scube);
              cosk = acos3(t);
              v3[0] = (s + s) * cosk - po3;
              n3 = 1;
              sinsqk = doub1 - cosk * cosk;
              if (sinsqk >= doub0) {
                rt3sink = rt3 * sqrt(sinsqk);
                v3[1] = s * (-cosk + rt3sink) - po3;
                v3[2] = s * (-cosk - rt3sink) - po3;
                n3 = 3;
                ++ncub[25];
                // printf("cubic root case 25\n");
              } else {
                ++ncub[26];
                // printf("cubic root case 26\n");
              }
            }
          } else
          /*
              cubic has multiple root -
          */
          {
            ++ncub[27];
            // printf("cubic root case 27\n");
            v3[0] = curoot(v) - po3;
            v3[1] = v3[0];
            v3[2] = v3[0];
            n3 = 3;
          }
        }
      }
    }
  }
done:
  if (debug < 1) {
    // printf("cubic %d %g %g %g\n",n3,p,q,r);
    for (j = 0; j < n3; ++j)
      printf("   %d %13g %13g\n", j, v3[j],
             r + v3[j] * (q + v3[j] * (p + v3[j])));
    printf("v %g,  uo3 %g,  m1 %g,   m2 %g,  po3 %g, cosk %g\n", v, uo3, m1, m2,
           po3, cosk);
    for (j = 0; j < 28; ++j) {
      printf("  %d", ncub[j]);
      if ((j % 10) == 9)
        printf("\n");
    }
    printf("\n");
  }
  if (iterate == TRUE)
    cubnewton(p, q, r, n3, v3);
  return (n3);
} /* cubic */
/****************************************************/

void quanewton(double p, double q, double r, double s, int n4, double v4[4])
/*
   improve roots of quartic by Newton-Raphson iteration

   6 Feb 2007  Arno Formella

*/
{
  int j, k;
  double corr, deriv, err, root;

  if (debug < 2)
    printf("quanewtona %d %g\n", n4, v4[0]);
  for (j = 0; j < n4; ++j) {
    for (k = 0; k < 4; ++k) {
      root = v4[j];
      err = (((root + p) * root + q) * root + r) * root + s;
      deriv = ((doub4 * root + doub3 * p) * root + doub2 * q) * root + r;
      if (deriv != doub0)
        corr = err / deriv;
      else
        corr = doub0;
      v4[j] -= corr;
      if (debug < 1)
        printf("quanewtonb %d %d %g %g %g %g %g\n", j, k, root, err, deriv,
               corr, v4[j]);
    }
  }
} /* quanewton */
/****************************************************/

void cubnewton(double p, double q, double r, int n3, double v3[4])
/*
   improve roots of cubic by Newton-Raphson iteration

   5 Jan 2004  Don Herbison-Evans

   called by cubic.
*/
{
  int j, k;
  double corr, deriv, err, root;

  if (debug < 2)
    printf("cubnewtona %d %g\n", n3, v3[0]);
  for (j = 0; j < n3; ++j) {
    for (k = 0; k < 4; ++k) {
      root = v3[j];
      err = ((root + p) * root + q) * root + r;
      deriv = (doub3 * root + doub2 * p) * root + q;
      if (deriv != doub0)
        corr = err / deriv;
      else
        corr = doub0;
      v3[j] -= corr;
      if (debug < 1)
        printf("cubnewtonb %d %d %g %g %g %g %g\n", j, k, root, err, deriv,
               corr, v3[j]);
    }
  }
} /* cubnewton */
/****************************************************/

int quartic(double a, double b, double c, double d, double rts[4])
/*
  Solve quartic equation using either
  quadratic, Ferrari's or Neumark's algorithm.

  called by
  calls  descartes, ferrari, neumark, yacfraid.

    15 Dec 2003  added yacfraid
    10 Dec 2003  added descartes with neg coeffs
    21 Jan 1989  Don Herbison-Evans
*/
{
  int j, k, nq, nr;
  double roots[4];

  if (fabs(a) > doubmax)
    nr = yacfraid(a, b, c, d, rts);
  else if ((a == doub0) && (c == doub0)) {
    nq = quadratic(b, d, roots);
    nr = 0;
    for (j = 0; j < nq; ++j) {
      if (roots[0] >= doub0) {
        rts[0] = sqrt(roots[0]);
        rts[1] = -rts[0];
        nr = 2;
      }
      if (roots[1] >= doub0) {
        rts[nr] = sqrt(roots[1]);
        rts[nr + 1] = -rts[nr];
        nr += 2;
      }
    }
  } else {
    k = 0;
    if (a < doub0)
      k += 2;
    if (b < doub0)
      k += 1;
    if (c < doub0)
      k += 8;
    if (d < doub0)
      k += 4;
    switch (k) {
    case 0:
      nr = neumark(a, b, c, d, rts);
      break;
    case 1:
      nr = neumark(a, b, c, d, rts);
      break;
    case 2:
      nr = neumark(a, b, c, d, rts);
      break;
    case 3:
      nr = ferrari(a, b, c, d, rts);
      break;
    case 4:
      nr = neumark(a, b, c, d, rts);
      break;
    case 5:
      nr = descartes(a, b, c, d, rts);
      break;
    case 6:
      nr = neumark(a, b, c, d, rts);
      break;
    case 7:
      nr = neumark(a, b, c, d, rts);
      break;
    case 8:
      nr = neumark(a, b, c, d, rts);
      break;
    case 9:
      nr = ferrari(a, b, c, d, rts);
      break;
    case 10:
      nr = neumark(a, b, c, d, rts);
      break;
    case 11:
      nr = neumark(a, b, c, d, rts);
      break;
    case 12:
      nr = neumark(a, b, c, d, rts);
      break;
    case 13:
      nr = neumark(a, b, c, d, rts);
      break;
    case 14:
      nr = neumark(a, b, c, d, rts);
      break;
    case 15:
      nr = descartes(-a, b, -c, d, rts);
      break;
    }
    if (k == 15)
      for (j = 0; j < nr; ++j)
        rts[j] = -rts[j];
  }
  return (nr);
} /* quartic */
/*****************************************/

int descartes(double a, double b, double c, double d, double rts[4])
/*
  Solve quartic equation using
  Descartes-Euler-Cardano algorithm

  called by quartic, compare, quartictest.

  Strong, T. "Elemementary and Higher Algebra"
     Pratt and Oakley, p. 469 (1859)

  16 Jul 1981  Don Herbison-Evans
*/
{
  int j;
  double v1[4], v2[4], v3[4];
  double k, y;
  double p, q, r;
  double e0, e1, e2;
  double g, h;
  double asq;
  double ainv4;
  double e1invk;

  asq = a * a;
  e2 = b - asq * d3o8;
  e1 = c + a * (asq * inv8 - b * inv2);
  e0 = d + asq * (b * inv16 - asq * d3o256) - a * c * inv4;

  p = doub2 * e2;
  q = e2 * e2 - doub4 * e0;
  r = -e1 * e1;

  n3 = cubic(p, q, r, v3);
  for (j3 = 0; j3 < n3; ++j3) {
    y = v3[j3];
    if (y <= doub0) {
      n4[j3] = 0;
      ++ndes[0];
    } /* y<0 */
    else {
      k = sqrt(y);
      ainv4 = a * inv4;
      e1invk = e1 / k;
      g = (y + e2 + e1invk) * inv2;
      h = (y + e2 - e1invk) * inv2;
      n1 = quadratic(-k, g, v1);
      n2 = quadratic(k, h, v2);
      qrts[0][j3] = v1[0] - ainv4;
      qrts[1][j3] = v1[1] - ainv4;
      qrts[n1][j3] = v2[0] - ainv4;
      qrts[n1 + 1][j3] = v2[1] - ainv4;
      n4[j3] = n1 + n2;
      ++ndes[1];
    } /* y>=0 */
    for (j = 0; j < n4[j3]; ++j)
      rts[j] = qrts[j][j3];
    worst3[j3] = errors(a, b, c, d, rts, rterd, n4[j3]);
  } /* j3 loop */
  j3 = 0;
  if (n3 != 1) {
    if ((n4[0] == n4[1]) && (n4[1] == n4[2]))
      ++ndes[NCASES - 2];
    else
      ++ndes[NCASES - 1];
    if ((n4[1] > n4[j3]) || ((worst3[1] < worst3[j3]) && (n4[1] == n4[j3])))
      j3 = 1;
    if ((n4[2] > n4[j3]) || ((worst3[2] < worst3[j3]) && (n4[2] == n4[j3])))
      j3 = 2;
  }
  for (j = 0; j < n4[j3]; ++j)
    rts[j] = qrts[j][j3];
  if (debug < 1)
    printf("descartes chose cubic %d %g %g\n\n", j3, v3[j3], worst3[j3]);
  ++ndes[30 + n4[j3]];
  ++ndes[35 + j3];
  return (n4[j3]);
} /* descartes */
/****************************************************/

int ferrari(double a, double b, double c, double d, double rts[4])
/*
    solve the quartic equation -

    x**4 + a*x**3 + b*x**2 + c*x + d = 0

    called by quartic, compare, quartictest.
    calls     cubic, quadratic.

    input -
    a,b,c,e - coeffs of equation.

    output -
    n4 - number of real roots.
    rts - array of root values.

    method :  Ferrari - Lagrange
    Theory of Equations, H.W. Turnbull p. 140 (1947)

    16 Jul 1981 Don Herbison-Evans

    calls  cubic, quadratic
*/
{
  int j;
  double asqinv4;
  double ainv2;
  double d4;
  double yinv2;
  double v1[4], v2[4], v3[4];
  double p, q, r;
  double y;
  double e, f, esq, fsq, ef;
  double g, gg, h, hh;

  ainv2 = a * inv2;
  asqinv4 = ainv2 * ainv2;
  d4 = d * doub4;

  p = b;
  q = a * c - d4;
  r = (asqinv4 - b) * d4 + c * c;
  n3 = cubic(p, q, r, v3);
  for (j3 = 0; j3 < n3; ++j3) {
    y = v3[j3];
    yinv2 = y * inv2;
    esq = asqinv4 - b - y;
    fsq = yinv2 * yinv2 - d;
    if ((esq < doub0) && (fsq < doub0)) {
      n4[j3] = 0;
      ++nfer[0];
    } else {
      ef = -(inv4 * a * y + inv2 * c);
      if (((a > doub0) && (y > doub0) && (c > doub0)) ||
          ((a > doub0) && (y < doub0) && (c < doub0)) ||
          ((a < doub0) && (y > doub0) && (c < doub0)) ||
          ((a < doub0) && (y < doub0) && (c > doub0)) || (a == doub0) ||
          (y == doub0) || (c == doub0))
      /* use ef - */
      {
        if ((b < doub0) && (y < doub0)) {
          e = sqrt(esq);
          f = ef / e;
          ++nfer[1];
        } else if (d < doub0) {
          f = sqrt(fsq);
          e = ef / f;
          ++nfer[2];
        } else {
          if (esq > doub0) {
            e = sqrt(esq);
            ++nfer[3];
          } else {
            e = doub0;
            ++nfer[4];
          }
          if (fsq > doub0) {
            f = sqrt(fsq);
            ++nfer[5];
          } else {
            f = doub0;
            ++nfer[6];
          }
          if (ef < doub0) {
            f = -f;
            ++nfer[7];
          }
        }
      } else
      /* use esq and fsq - */
      {
        if (esq > doub0) {
          e = sqrt(esq);
          ++nfer[8];
        } else {
          e = doub0;
          ++nfer[9];
        }
        if (fsq > doub0) {
          f = sqrt(fsq);
          ++nfer[10];
        } else {
          f = doub0;
          ++nfer[11];
        }
        if (ef < doub0) {
          f = -f;
          ++nfer[12];
        }
      }
      /* note that e >= doub0 */
      g = ainv2 - e;
      gg = ainv2 + e;
      if (((b > doub0) && (y > doub0)) || ((b < doub0) && (y < doub0))) {
        if (((a > doub0) && (e > doub0)) || ((a < doub0) && (e < doub0))) {
          g = (b + y) / gg;
          ++nfer[13];
        } else if (((a > doub0) && (e < doub0)) ||
                   ((a < doub0) && (e > doub0))) {
          gg = (b + y) / g;
          ++nfer[14];
        } else
          ++nfer[15];
      }
      hh = -yinv2 + f;
      h = -yinv2 - f;
      if (((f > doub0) && (y < doub0)) || ((f < doub0) && (y > doub0))) {
        h = d / hh;
        ++nfer[16];
      } else if (((f < doub0) && (y < doub0)) || ((f > doub0) && (y > doub0))) {
        hh = d / h;
        ++nfer[17];
      } else
        ++nfer[18];

      n1 = quadratic(gg, hh, v1);
      n2 = quadratic(g, h, v2);
      n4[j3] = n1 + n2;
      qrts[0][j3] = v1[0];
      qrts[1][j3] = v1[1];
      qrts[n1 + 0][j3] = v2[0];
      qrts[n1 + 1][j3] = v2[1];
    }
    for (j = 0; j < n4[j3]; ++j)
      rts[j] = qrts[j][j3];
    worst3[j3] = errors(a, b, c, d, rts, rterf, n4[j3]);
  } /* j3 loop */
  j3 = 0;
  if (n3 != 1) {
    if ((n4[0] == n4[1]) && (n4[1] == n4[2]))
      ++nfer[NCASES - 2];
    else
      ++nfer[NCASES - 1];
    if ((n4[1] > n4[j3]) || ((worst3[1] < worst3[j3]) && (n4[1] == n4[j3])))
      j3 = 1;
    if ((n4[2] > n4[j3]) || ((worst3[2] < worst3[j3]) && (n4[2] == n4[j3])))
      j3 = 2;
  }
  for (j = 0; j < n4[j3]; ++j)
    rts[j] = qrts[j][j3];
  if (debug < 1)
    printf("ferrari chose cubic %d %g %g\n\n", j3, v3[j3], worst3[j3]);
  ++nfer[30 + n4[j3]];
  ++nfer[35 + j3];
  return (n4[j3]);
} /* ferrari */
/*****************************************/

int neumark(double a, double b, double c, double d, double rts[4])
/*
    solve the quartic equation -

    x**4 + a*x**3 + b*x**2 + c*x + d = 0

    called by quartic, compare, quartictest.
    calls     cubic, quadratic.

    input parameters -
    a,b,c,e - coeffs of equation.

    output parameters -
    n4 - number of real roots.
    rts - array of root values.

    method -  S. Neumark
       "Solution of Cubic and Quartic Equations" - Pergamon 1965

     1 Dec 1985   translated to C with help of Shawn Neely
    16 Jul 1981   Don Herbison-Evans

*/
{
  int j;
  double y, g, gg, h, hh, gdis, gdisrt, hdis, hdisrt, g1, g2, h1, h2;
  double bmy, gerr, herr, y4, bmysq;
  double v1[4], v2[4], v3[4];
  double asq;
  double d4;
  double p, q, r;
  double hmax, gmax;

  if (d == doub0) {
    n3 = 0;
    n4[0] = cubic(a, b, c, rts);
    for (j = 0; j < n4[0]; ++j)
      qrts[j][0] = rts[j];
    qrts[n4[0]++][0] = doub0;
    goto done;
  }
  asq = a * a;
  d4 = d * doub4;
  p = -b * doub2;
  q = b * b + a * c - d4;
  r = (c - a * b) * c + asq * d;
  if (debug < 3)
    printf("neumarka %g %g %g %g,  %g %g %g\n", a, b, c, d, p, q, r);
  n3 = cubic(p, q, r, v3);
  // printf("cubic %d %g\n",n3,v3[0]);
  for (j3 = 0; j3 < n3; ++j3) {
    y = v3[j3];
    bmy = b - y;
    y4 = y * doub4;
    bmysq = bmy * bmy;
    gdis = asq - y4;
    hdis = bmysq - d4;
    if (debug < 3)
      printf("neumarkb %g %g\n", gdis, hdis);
    if ((gdis < doub0) || (hdis < doub0)) {
      n4[j3] = 0;
      ++nneu[0];
    } else {
      g1 = a * inv2;
      h1 = bmy * inv2;
      gerr = asq + y4;
      herr = hdis;
      if (d > doub0) {
        herr = bmysq + d4;
        ++nneu[1];
      }
      if ((y < doub0) || (herr * gdis > gerr * hdis)) {
        gdisrt = sqrt(gdis);
        g2 = gdisrt * inv2;
        if (gdisrt != doub0) {
          h2 = (a * h1 - c) / gdisrt;
          ++nneu[2];
        } else {
          h2 = doub0;
          ++nneu[3];
        }
      } else {
        hdisrt = sqrt(hdis);
        h2 = hdisrt * inv2;
        if (hdisrt != doub0) {
          g2 = (a * h1 - c) / hdisrt;
          ++nneu[4];
        } else {
          g2 = doub0;
          ++nneu[5];
        }
      }
      /*
          note that in the following, the tests ensure non-zero
          denominators -
      */
      h = h1 - h2;
      hh = h1 + h2;
      hmax = hh;
      if (hmax < doub0) {
        hmax = -hmax;
        ++nneu[6];
      }
      if (hmax < h) {
        hmax = h;
        ++nneu[7];
      }
      if (hmax < -h) {
        hmax = -h;
        ++nneu[8];
      }
      if ((h1 > doub0) && (h2 > doub0)) {
        h = d / hh;
        ++nneu[9];
      }
      if ((h1 < doub0) && (h2 < doub0)) {
        h = d / hh;
        ++nneu[10];
      }
      if ((h1 > doub0) && (h2 < doub0)) {
        hh = d / h;
        ++nneu[11];
      }
      if ((h1 < doub0) && (h2 > doub0)) {
        hh = d / h;
        ++nneu[12];
      }
      if (h > hmax) {
        h = hmax;
        ++nneu[13];
      }
      if (h < -hmax) {
        h = -hmax;
        ++nneu[14];
      }
      if (hh > hmax) {
        hh = hmax;
        ++nneu[15];
      }
      if (hh < -hmax) {
        hh = -hmax;
        ++nneu[16];
      }

      g = g1 - g2;
      gg = g1 + g2;
      gmax = gg;
      if (gmax < doub0) {
        gmax = -gmax;
        ++nneu[17];
      }
      if (gmax < g) {
        gmax = g;
        ++nneu[18];
      }
      if (gmax < -g) {
        gmax = -g;
        ++nneu[19];
      }
      if ((g1 > doub0) && (g2 > doub0)) {
        g = y / gg;
        ++nneu[20];
      }
      if ((g1 < doub0) && (g2 < doub0)) {
        g = y / gg;
        ++nneu[21];
      }
      if ((g1 > doub0) && (g2 < doub0)) {
        gg = y / g;
        ++nneu[22];
      }
      if ((g1 < doub0) && (g2 > doub0)) {
        gg = y / g;
        ++nneu[23];
      }
      if (g > gmax) {
        g = gmax;
        ++nneu[24];
      }
      if (g < -gmax) {
        g = -gmax;
        ++nneu[25];
      }
      if (gg > gmax) {
        gg = gmax;
        ++nneu[26];
      }
      if (gg < -gmax) {
        gg = -gmax;
        ++nneu[27];
      }

      n1 = quadratic(gg, hh, v1);
      n2 = quadratic(g, h, v2);
      // printf("quadratic roots 1: %d %g %g\n",n1,v1[0],v1[2]);
      // printf("quadratic roots 2: %d %g %g\n",n2,v2[0],v2[2]);
      n4[j3] = n1 + n2;
      qrts[0][j3] = v1[0];
      qrts[1][j3] = v1[1];
      qrts[n1 + 0][j3] = v2[0];
      qrts[n1 + 1][j3] = v2[1];
    }
    for (j = 0; j < n4[j3]; ++j)
      rts[j] = qrts[j][j3];
    worst3[j3] = errors(a, b, c, d, rts, rtern, n4[j3]);
  } /* j3 loop */
done:
  j3 = 0;
  if (n3 > 1) {
    if ((n4[0] == n4[1]) && (n4[1] == n4[2]))
      ++nneu[NCASES - 2];
    else
      ++nneu[NCASES - 1];
    if ((n4[1] > n4[j3]) || ((worst3[1] < worst3[j3]) && (n4[1] == n4[j3])))
      j3 = 1;
    if ((n4[2] > n4[j3]) || ((worst3[2] < worst3[j3]) && (n4[2] == n4[j3])))
      j3 = 2;
  }
  for (j = 0; j < n4[j3]; ++j)
    rts[j] = qrts[j][j3];
  if (debug < 1)
    printf("neumark chose cubic %d %g %g\n\n", j3, v3[j3], worst3[j3]);
  ++nneu[30 + n4[j3]];
  ++nneu[35 + j3];
  return (n4[j3]);
} /* neumark */
/****************************************************/

int yacfraid(double a, double b, double c, double d, double rts[4])
/*
    solve the quartic equation -

    x**4 + a*x**3 + b*x**2 + c*x + d = 0

    called by quartic, compare, quartictest.
    calls     cubic, quadratic.

    input parameters -
    a,b,c,e - coeffs of equation.

    output parameters -
    n4 - number of real roots.
    rts - array of root values.

    method -
       K.S. Brown
       Reducing Quartics to Cubics,
       http://www.seanet.com/~ksbrown/kmath296.htm (1967)

       Michael Daoud Yacoub & Gustavo Fraidenraich
       "A new simple solution of the general quartic equation"
       Revised 16 Feb 2004

    14 Nov 2003 Don Herbison-Evans
*/
{
  int j;
  double y;
  double v3[4];
  double asq, acu;
  double b4;
  double det0, det1, det2, det3;
  double det0rt, det1rt, det2rt, det3rt;
  double e, f, g = 0, h = 0, k;
  double fsq, gsq, hsq, invk;
  double P, Q, R, U;

  if (d == doub0) {
    n3 = 0;
    n4[0] = cubic(a, b, c, rts);
    for (j = 0; j < n4[0]; ++j)
      qrts[j][0] = rts[j];
    qrts[n4[0]++][0] = doub0;
    goto done;
  }
  asq = a * a;
  acu = a * asq;
  b4 = b * doub4;
  n3 = 0;

  P = asq * b - b4 * b + doub2 * a * c + doub16 * d;
  Q = asq * c - b4 * c + doub8 * a * d;
  R = asq * d - c * c;
  U = acu - b4 * a + doub8 * c;
  n4[0] = 0;
  if (U == doub0) {
    if (P == doub0) {
      det0 = doub3 * asq - doub8 * b;
      if (det0 < doub0) {
        ++nyac[0];
        goto done;
      }
      det0rt = sqrt(det0);
      qrts[0][0] = (-a + det0rt) * inv4;
      qrts[1][0] = qrts[0][0];
      qrts[2][0] = (-a - det0rt) * inv4;
      qrts[3][0] = qrts[2][0];
      ++nyac[1];
      n4[0] = 4;
      goto done;
    } /* P=0 */
    else {
      det1 = asq * asq - doub8 * asq * b + doub16 * b * b - doub64 * d;
      if (det1 < doub0) {
        ++nyac[2];
        goto done;
        ;
      }
      n4[0] = 0;
      det1rt = sqrt(det1);
      det2 = doub3 * asq - doub8 * b + doub2 * det1rt;
      if (det2 >= doub0) {
        det2rt = sqrt(det2);
        qrts[0][0] = (-a + det2rt) * inv4;
        qrts[1][0] = (-a - det2rt) * inv4;
        n4[0] = 2;
        ++nyac[3];
      }
      det3 = doub3 * asq - doub8 * b - doub2 * det1rt;
      if (det3 >= doub0) {
        det3rt = sqrt(det3);
        qrts[n4[0]++][0] = (-a + det3rt) * inv4;
        qrts[n4[0]++][0] = (-a - det3rt) * inv4;
        ++nyac[5];
      }
      if (n4[0] == 0)
        ++nyac[6];
      goto done;
    } /* P<>0 */
  }

  n3 = cubic(P / U, Q / U, R / U, v3);
  for (j3 = 0; j3 < n3; ++j3) {
    y = v3[j3];
    j = 0;
    k = a + doub4 * y;
    if (k == doub0) {
      ++nyac[9];
      goto donej3;
    }
    invk = doub1 / k;
    e = (acu - doub4 * c - doub2 * a * b + (doub6 * asq - doub16 * b) * y) *
        invk;
    fsq = (acu + doub8 * c - doub4 * a * b) * invk;
    if (fsq < doub0) {
      ++nyac[10];
      goto donej3;
    }
    f = sqrt(fsq);
    gsq = doub2 * (e + f * k);
    hsq = doub2 * (e - f * k);
    if (gsq >= doub0) {
      ++nyac[11];
      g = sqrt(gsq);
      qrts[j++][j3] = (-a - f - g) * inv4;
      qrts[j++][j3] = (-a - f + g) * inv4;
    }
    if (hsq >= doub0) {
      ++nyac[12];
      h = sqrt(hsq);
      qrts[j++][j3] = (-a + f - h) * inv4;
      qrts[j++][j3] = (-a + f + h) * inv4;
    }
    if (debug < 1) {
      printf("j3 %d y %g k %g fsq %g gsq %g hsq %g\n", j3, y, k, fsq, gsq, hsq);
      printf("e %g f %g g %g h %g\n", e, f, g, h);
    }
  donej3:
    n4[j3] = j;
    for (j = 0; j < n4[j3]; ++j)
      rts[j] = qrts[j][j3];
    worst3[j3] = errors(a, b, c, d, rts, rtery, n4[j3]);
  } /* j3 loop */
done:
  j3 = 0;
  if (n3 > 1) {
    if ((n4[0] == n4[1]) && (n4[1] == n4[2]))
      ++nyac[NCASES - 2];
    else {
      ++nyac[NCASES - 1];
      if ((n4[0] != n4[1]) && (n4[0] != n4[2]) && (n4[1] != n4[2]))
        printf("yace %d %d %d %g %g %g %g\n", n4[0], n4[1], n4[2], a, b, c, d);
    }
    if ((n4[1] > n4[j3]) || ((worst3[1] < worst3[j3]) && (n4[1] == n4[j3])))
      j3 = 1;
    if ((n4[2] > n4[j3]) || ((worst3[2] < worst3[j3]) && (n4[2] == n4[j3])))
      j3 = 2;
  }
  for (j = 0; j < n4[j3]; ++j)
    rts[j] = qrts[j][j3];
  ++nyac[30 + n4[j3]];
  ++nyac[35 + j3];
  if (debug < 1)
    printf("yacfraid chose cubic %d %g  %g\n\n", j3, v3[j3], worst3[j3]);
  return (n4[j3]);
} /* yacfraid */
/*****************************************/

int chris(double a, double b, double c, double d, double rts[4])
/*
  Solve quartic equation using
  Christianson's algorithm

  called by compare, quartictest.
  calls errors, quadratic, cubic.

  Bruce Christianson, Solving Quartics Using Palindromes,
     Mathematical Gazette, Vol. 75, pp. 327-328 (1991)

  14 Jan 2004  Don Herbison-Evans
*/
{
  int j;
  int n0, n1, n2;
  double v0[4], v1[4], v2[4], v3[4];
  double y, ysq, ycu;
  double p, q, r;
  double asq, acu, aqu, ao4;
  double e0, e1, e2;
  double k, ksq, kcu, kqu, kquinv, k1, k2;
  double g, h, g1, Z1, Z2;

  if (d == doub0) {
    n3 = 0;
    n4[0] = cubic(a, b, c, rts);
    for (j = 0; j < n4[0]; ++j)
      qrts[j][0] = rts[j];
    qrts[n4[0]++][0] = doub0;
    goto done;
  }
  asq = a * a;
  acu = asq * a;
  aqu = acu * a;
  ao4 = a * inv4;
  e2 = b - d3o8 * asq;
  e1 = c - inv2 * b * a + inv8 * acu;
  e0 = d - inv4 * c * a + inv16 * b * asq - d3o256 * acu * a;
  if (debug < 1)
    printf("chrisa e0 %g   e1 %g   e2 %g\n", e0, e1, e2);
  if (e1 == doub0) {
    n3 = 0;
    n4[0] = 0;
    ++nchr[0];
    n3 = 3;
    v3[2] = doub0;
    v3[1] = doub0;
    v3[0] = -inv8 * (doub16 * e0 - doub4 * e2 * e2);
  } else {
    p = (inv2 * b * asq - inv2 * b * b - inv2 * c * a + doub2 * d -
         doub3 * inv32 * aqu) /
        e1;
    q = doub3 * inv16 * asq - inv2 * b;
    r = inv16 * b * a - inv8 * c - inv64 * acu;
    if (debug < 1)
      printf("chrisb %g %g %g\n", p, q, r);
    ++nchr[1];
    n3 = cubic(p, q, r, v3);
  }
  for (j3 = 0; j3 < n3; ++j3) {
    y = v3[j3];
    n4[j3] = 0;
    ysq = y * y;
    ycu = y * ysq;
    if ((y == doub0) ||
        ((y < doub0) && (e1 <= doub0) && (e0 >= doub0) && (e2 >= doub0))) {
      kqu = y * ycu + e2 * ysq + e1 * y + e0;
      if (kqu <= doub0) {
        ksq = doub0;
        ++nchr[2];
        goto donej3;
      }
      ksq = sqrt(kqu);
    } else {
      ksq = ysq + inv2 * e2 + inv4 * e1 / y;
      if (ksq <= doub0) {
        kqu = doub0;
        ++nchr[3];
        goto donej3;
      }
      kqu = ksq * ksq;
    }
    k = sqrt(ksq);
    kcu = k * ksq;
    kquinv = doub1 / kqu;
    g = doub4 * y * kcu;
    h = (doub6 * ysq + e2) * ksq;
    if (debug < 1) {
      k2 = sqrt(sqrt(fabs(y * ycu + e2 * ysq + e1 * y + e0)));
      if (y == doub0)
        k1 = doub0;
      else
        k1 = sqrt(fabs(ysq + inv2 * e2 + inv4 * e1 / y));
      g1 = (doub4 * ycu + doub2 * y * e2 + e1) * k;
      // printf("chrisc %g %g %g %g %g\n",k,k1,y,g,g1,h);
    }
    n0 = quadratic(g * kquinv, h * kquinv - doub2, v0);
    if (n0 < 1) {
      ++nchr[4];
      goto donej3;
    }
    Z1 = v0[0];
    Z2 = v0[1];
    n1 = quadratic(-Z1, doub1, v1);
    if (n1 > 0) {
      ++nchr[5];
      n4[j3] = n1;
      qrts[0][j3] = y + k * v1[0] - ao4;
      qrts[1][j3] = y + k * v1[1] - ao4;
    }
    n2 = quadratic(-Z2, doub1, v2);
    if (n2 > 0) {
      ++nchr[6];
      n4[j3] += n2;
      qrts[n1][j3] = y + k * v2[0] - ao4;
      qrts[n1 + 1][j3] = y + k * v2[1] - ao4;
    }
  donej3:
    if (debug < 1)
      printf("chrisd %g %g\n", ksq, kqu);
    for (j = 0; j < n4[j3]; ++j)
      rts[j] = qrts[j][j3];
    worst3[j3] = errors(a, b, c, d, rts, rterc, n4[j3]);
  } /* j3 loop */
done:
  j3 = 0;
  if (n3 > 1) {
    if ((n4[0] == n4[1]) && (n4[1] == n4[2]))
      ++nchr[NCASES - 2];
    else {
      ++nchr[NCASES - 1];
      if (debug < 1)
        printf("chrise %d %d %d %g %g %g %g\n", n4[0], n4[1], n4[2], a, b, c,
               d);
    }
    if ((n4[1] > n4[j3]) || ((worst3[1] < worst3[j3]) && (n4[1] == n4[j3])))
      j3 = 1;
    if ((n4[2] > n4[j3]) || ((worst3[2] < worst3[j3]) && (n4[2] == n4[j3])))
      j3 = 2;
  }
  for (j = 0; j < n4[j3]; ++j)
    rts[j] = qrts[j][j3];
  if (debug < 1)
    printf("chris chose cubic %d %g %g\n\n", j3, v3[j3], worst3[j3]);
  ++nchr[30 + n4[j3]];
  ++nchr[35 + j3];
  return (n4[j3]);
} /* chris */

/***********************************/
