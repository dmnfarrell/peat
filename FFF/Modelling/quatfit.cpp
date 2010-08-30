/*
 #
 # FFF - Flexible Force Field
 # Copyright (C) 2010 Jens Erik Nielsen
 #
 # This program is free software: you can redistribute it and/or modify
 # it under the terms of the GNU General Public License as published by
 # the Free Software Foundation, either version 3 of the License, or
 # (at your option) any later version.
 #
 # This program is distributed in the hope that it will be useful,
 # but WITHOUT ANY WARRANTY; without even the implied warranty of
 # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 # GNU General Public License for more details.
 #
 # You should have received a copy of the GNU General Public License
 # along with this program.  If not, see <http://www.gnu.org/licenses/>.
 #
 # Contact information: 
 # Email: Jens.Nielsen_at_gmail.com
 # Normal mail:
 # Jens Nielsen
 # SBBS, Conway Institute
 # University College Dublin
 # Dublin 4, Ireland
 
  The program to superimpose atoms of two molecules by quaternion method

 David J. Heisterberg
 The Ohio Supercomputer Center
 1224 Kinnear Rd.
 Columbus, OH  43212-1163
 (614)292-6036
 djh@osc.edu    djh@ohstpy.bitnet    ohstpy::djh

 Translated to C from fitest.f program and interfaced with Xmol program
 by Jan Labanowski,  jkl@osc.edu   jkl@ohstpy.bitnet   ohstpy::jkl

 To complile:
   cc -o quatfit quatfit.c -lm

 Copyright: Ohio Supercomputer Center, David J. Heisterberg, 1990.
 The program can be copied and distributed freely, provided that
 this copyright in not removed. You may acknowledge the use of the
 program in published material as:
 David J. Heisterberg, 1990, unpublished results.

*/

#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include <string>
#include <stdlib.h>
#include "quatfit.h"

/* options
  -r refmol    reference molecule xmol file. If this option is not given, the
               information is read from standard input.

  -f fitmol    fitted molecule input xmol file If this option is not given, the
               information is read from standard input.

  -p pairs     file with the list of fitted atom pairs and weights. If this
               option is not specified, pairs and weights are taken from
               stdin. If file name "none" is used (i.e. -p none), atoms of the
               fitted molecule are fitted to atoms of the reference
               molecule with the same number, and all weights are assumed 1.0.
               If molecules do not have the same number of atoms, the
               smaller number of atoms is fitted.

  -o outmol    output file for transformed fitted molecule. If this option
               is not given, the information is written to standard output.

  -s statfile  file with fit statistics (distances between fitted atoms, etc).
               If this option is not given, the information is written to
               standard output.

  If any files are read from stdin, the order is: refmol, fitmol, pairs.
  If any files are written to stdout, the order is: outmol, statfile.

The file formats are:
 The refmol, fitmol and outmol files are in the XYZ format used with xmol
 program from Minnesota Supercomputer Institute. The format is:
   1st line: number of atoms
   2nd line: title
   3rd and next lines have the format depending on the kind of information:
       T  X  Y  Z                (total of 4 columns)
       T  X  Y  Z  C             (total of 5 columns)
       T  X  Y  Z  Mx My Mz      (total of 7 columns)
       T  X  Y  Z  C Mx My Mz    (total of 8 columns)
     where T is atom type (usually elemnt symbol), X, Y, Z are cartesian
     coordinates of an atom in Angstroms, C is atomic charge, and Mx, My, Mz
     are normal modes.

 The pairs file format is:
   1st line: number of pairs
   2nd and next lines:
       Ar   Af    W
     where Ar is the atom number of the reference molecule, Af is the atom
     number of fitted molecule, w is the statistical weight. Weights W are
     related to the square of expected deviation "sigma" between the reference
     and fitted molecule atoms and allow to make fit of some atom pairs more
     tight. W is proportional to 1/sigma^2. The larger the weight, the more
     tight will be the resulting fit for the given pair.

 The statfile lists results of the fit with explanation.

*/

/*====================================================================
CENTER
 center or translate a molecule.
 n - number of atoms
 x - on input  - original xyz coordinates of a molecule
     on output - moved xyz coordinates (see io for modes).

 w - if io=1, weights of atoms
     if io=2 or 3, unused

 io - 1 weighted geometric center of the molecule will be at (0,0,0)
      2 molecule will be moved by a vector -o (i.e., components of a vector o
        will be subtracted from atom coordinates).
      3 molecule will be moved by a vector +o (i.e., components of a vector o
        will be added atom coordinates).

 o - if io=1, output, center of original coordinates
     if io=2, input, vector o will be subtracted from atomic coordinates
     if io=3, input, vector o will be added to atomic coordinates

=====================================================================*/

void center(int n, double x[4][MAXPOINTS],double w[MAXPOINTS], int io, double o[4]) {
  double wnorm, modif;
  int i;

  if (io == 2) {
    modif = -1.0;
  }
  else if (io == 3) {
    modif = 1.0;
  }
  else {
    modif = -1.0;
    o[1] = 0.0;
    o[2] = 0.0;
    o[3] = 0.0;
    wnorm = 0.0;
    for (i = 1; i <= n; i++) {
      o[1] = o[1] + x[1][i] * sqrt(w[i]);
      o[2] = o[2] + x[2][i] * sqrt(w[i]);
      o[3] = o[3] + x[3][i] * sqrt(w[i]);
      wnorm = wnorm + sqrt(w[i]);
    }
    o[1] = o[1] / wnorm;
    o[2] = o[2] / wnorm;
    o[3] = o[3] / wnorm;
  }
  
  
  for (i = 1; i <= n; i++) {
    x[1][i] = x[1][i] + modif*o[1];
    x[2][i] = x[2][i] + modif*o[2];
    x[3][i] = x[3][i] + modif*o[3];
  }

}

/*================================
 ROTMOL
 rotate a molecule
 n - number of atoms
 x - input coordinates
 y - rotated coordinates y = u * x
 u - left rotation matrix
==================================*/

void rotmol (int n, double x[4][MAXPOINTS], double y[4][MAXPOINTS], double u[4][4]) {
  double yx, yy, yz;
  int i;

  for (i = 1; i <= n; i++) {
    yx = u[1][1] * x[1][i] + u[1][2] * x[2][i] + u[1][3] * x[3][i];
    yy = u[2][1] * x[1][i] + u[2][2] * x[2][i] + u[2][3] * x[3][i];
    yz = u[3][1] * x[1][i] + u[3][2] * x[2][i] + u[3][3] * x[3][i];

    y[1][i] = yx;
    y[2][i] = yy;
    y[3][i] = yz;
  }
}

/*=======================================================
 JACOBI
 Jacobi diagonalizer with sorted output. It is only good for 4x4 matrices.
 (was too lazy to do pointers...)
 a - input: matrix to diagonalize
 v - output: eigenvectors
 d - output: eigenvalues
 nrot - input: maximum number of sweeps
=========================================================*/

void jacobi (double a[4][4], double d[4], double v[4][4], int nrot) {
  double onorm, dnorm;
  double b, dma, q, t, c, s;
  double atemp, vtemp, dtemp;
  int i, j, k, l;

  for (j = 0; j <= 3; j++) {
    for (i = 0; i <= 3; i++) {
      v[i][j] = 0.0;
    }
    v[j][j] = 1.0;
    d[j] = a[j][j];
  }

  for (l = 1; l <= nrot; l++) {
    dnorm = 0.0;
    onorm = 0.0;
    for (j = 0; j <= 3; j++) {
      dnorm = dnorm + fabs(d[j]);
      for (i = 0; i <= j - 1; i++) {
	onorm = onorm + fabs(a[i][j]);
      }
    }
    if((onorm/dnorm) <= 1.0e-12) goto Exit_now;
    for (j = 1; j <= 3; j++) {
      for (i = 0; i <= j - 1; i++) {
	b = a[i][j];
	if(fabs(b) > 0.0) {
	  dma = d[j] - d[i];
	  if((fabs(dma) + fabs(b)) <=  fabs(dma)) {
	    t = b / dma;
	  }
	  else {
	    q = 0.5 * dma / b;
	    t = 1.0/(fabs(q) + sqrt(1.0+q*q));
	    if(q < 0.0) {
	      t = -t;
	    }
	  }
	  c = 1.0/sqrt(t * t + 1.0);
	  s = t * c;
	  a[i][j] = 0.0;
	  for (k = 0; k <= i-1; k++) {
	    atemp = c * a[k][i] - s * a[k][j];
	    a[k][j] = s * a[k][i] + c * a[k][j];
	    a[k][i] = atemp;
	  }
	  for (k = i+1; k <= j-1; k++) {
	    atemp = c * a[i][k] - s * a[k][j];
	    a[k][j] = s * a[i][k] + c * a[k][j];
	    a[i][k] = atemp;
	  }
	  for (k = j+1; k <= 3; k++) {
	    atemp = c * a[i][k] - s * a[j][k];
           a[j][k] = s * a[i][k] + c * a[j][k];
           a[i][k] = atemp;
	  }
	  for (k = 0; k <= 3; k++) {
	    vtemp = c * v[k][i] - s * v[k][j];
	    v[k][j] = s * v[k][i] + c * v[k][j];
	    v[k][i] = vtemp;
	  }
         dtemp = c*c*d[i] + s*s*d[j] - 2.0*c*s*b;
         d[j] = s*s*d[i] + c*c*d[j] +  2.0*c*s*b;
         d[i] = dtemp;
	}  /* end if */
      } /* end for i */
    } /* end for j */
  } /* end for l */

 Exit_now:

  nrot = l;

  for (j = 0; j <= 2; j++) {
    k = j;
    dtemp = d[k];
    for (i = j+1; i <= 3; i++) {
      if(d[i] < dtemp) {
	k = i;
	dtemp = d[k];
      }
    }

    if(k > j) {
      d[k] = d[j];
      d[j] = dtemp;
      for (i = 0; i <= 3; i++) {
	dtemp = v[i][k];
	v[i][k] = v[i][j];
	v[i][j] = dtemp;
      }
    }
  }
}



/*==========================================
 Q2MAT
 Generate a left rotation matrix from a normalized quaternion

 INPUT
   q      - normalized quaternion

 OUTPUT
   u      - the rotation matrix
===========================================*/

void q2mat (double q[4], double u[4][4]) {
  u[1][1] = q[0]*q[0] + q[1]*q[1] - q[2]*q[2] - q[3]*q[3];
  u[2][1] = 2.0 * (q[1] * q[2] - q[0] * q[3]);
  u[3][1] = 2.0 * (q[1] * q[3] + q[0] * q[2]);

  u[1][2] = 2.0 * (q[2] * q[1] + q[0] * q[3]);
  u[2][2] = q[0]*q[0] - q[1]*q[1] + q[2]*q[2] - q[3]*q[3];
  u[3][2] = 2.0 * (q[2] * q[3] - q[0] * q[1]);

  u[1][3] = 2.0 *(q[3] * q[1] - q[0] * q[2]);
  u[2][3] = 2.0 * (q[3] * q[2] + q[0] * q[1]);
  u[3][3] = q[0]*q[0] - q[1]*q[1] - q[2]*q[2] + q[3]*q[3];
}


/*==========================================
 QTRFIT
 Find the quaternion, q,[and left rotation matrix, u] that minimizes

   |qTXq - Y| ^ 2  [|uX - Y| ^ 2]

 This is equivalent to maximizing Re (qTXTqY).

 This is equivalent to finding the largest eigenvalue and corresponding
 eigenvector of the matrix

 [A2   AUx  AUy  AUz ]
 [AUx  Ux2  UxUy UzUx]
 [AUy  UxUy Uy2  UyUz]
 [AUz  UzUx UyUz Uz2 ]

 where

   A2   = Xx Yx + Xy Yy + Xz Yz
   Ux2  = Xx Yx - Xy Yy - Xz Yz
   Uy2  = Xy Yy - Xz Yz - Xx Yx
   Uz2  = Xz Yz - Xx Yx - Xy Yy
   AUx  = Xz Yy - Xy Yz
   AUy  = Xx Yz - Xz Yx
   AUz  = Xy Yx - Xx Yy
   UxUy = Xx Yy + Xy Yx
   UyUz = Xy Yz + Xz Yy
   UzUx = Xz Yx + Xx Yz

 The left rotation matrix, u, is obtained from q by

   u = qT1q

 INPUT
   n      - number of points
   x      - fitted molecule coordinates
   y      - reference molecule coordinates
   w      - weights

 OUTPUT
   q      - the best-fit quaternion
   u      - the best-fit left rotation matrix
   nr     - max number of jacobi sweeps

=====================================*/

void qtrfit (int n, double x[4][MAXPOINTS], double y[4][MAXPOINTS], double w[MAXPOINTS], double q[4], double u[4][4], int nr) {
  double xxyx, xxyy, xxyz;
  double xyyx, xyyy, xyyz;
  double xzyx, xzyy, xzyz;
  double c[4][4], v[4][4];
  double d[4];
  int i, j;


  /* generate the upper triangle of the quadratic form matrix */

  xxyx = 0.0;
  xxyy = 0.0;
  xxyz = 0.0;
  xyyx = 0.0;
  xyyy = 0.0;
  xyyz = 0.0;
  xzyx = 0.0;
  xzyy = 0.0;
  xzyz = 0.0;

  for (i = 1; i <= n; i++) {
    xxyx = xxyx + x[1][i] * y[1][i] * w[i];
    xxyy = xxyy + x[1][i] * y[2][i] * w[i];
    xxyz = xxyz + x[1][i] * y[3][i] * w[i];
    xyyx = xyyx + x[2][i] * y[1][i] * w[i];
    xyyy = xyyy + x[2][i] * y[2][i] * w[i];
    xyyz = xyyz + x[2][i] * y[3][i] * w[i];
    xzyx = xzyx + x[3][i] * y[1][i] * w[i];
    xzyy = xzyy + x[3][i] * y[2][i] * w[i];
    xzyz = xzyz + x[3][i] * y[3][i] * w[i];
  }

  for(i = 0; i <= 3; i++) {
    for(j = 0; j <= 3; j++) {
      c[i][j] = 0.0;
    }
  }

  c[0][0] = xxyx + xyyy + xzyz;

  c[0][1] = xzyy - xyyz;
  c[1][1] = xxyx - xyyy - xzyz;

  c[0][2] = xxyz - xzyx;
  c[1][2] = xxyy + xyyx;
  c[2][2] = xyyy - xzyz - xxyx;

  c[0][3] = xyyx - xxyy;
  c[1][3] = xzyx + xxyz;
  c[2][3] = xyyz + xzyy;
  c[3][3] = xzyz - xxyx - xyyy;

  /* diagonalize c */

  jacobi (c, d, v, nr);

  /* extract the desired quaternion */

  q[0] = v[0][3];
  q[1] = v[1][3];
  q[2] = v[2][3];
  q[3] = v[3][3];

  /* generate the rotation matrix */

  q2mat (q, u);

}

/*============================================
  scan_line - scans line of XYZ file of Xmol and parses it into individual
  components.
  ==============================================*/

int scan_line(FILE *file, int *n_fields, char *symb, double *x, double *y, double *z, double *charge, double *mx, double *my, double *mz) {
  int n, i;
  char line[MAXLINELEN];

  *x = 0.0;
  *y = 0.0;
  *z = 0.0;
  *charge = 0.0;
  *mx = 0.0;
  *my = 0.0;
  *mz = 0.0;

  /* get a line of text from the file, exit EOF reached */
  if(fgets(line, MAXLINELEN, file) == NULL) {
    return((int)-1);
  }

  /* find number of fields in the line */
  i = 1;
  n = 0;
  while (line[i] != '\0') {
    if((!isspace(line[i-1])) && isspace(line[i])) {
      n++;
    }
    i++;
  }

  /* check if n_fields is set already */
  if(*n_fields == 0) {
    *n_fields = n;
  }

  /* check if n_fields is correct */
  if(*n_fields != n) {
    return((int)-2);
  }

  switch (n) {

  case 4: if(sscanf(line, "%s%le%le%le",  symb, x, y, z) != 4) {
    return((int)-3);
  }
  else {
    return((int)0);
  }

  case 5: if(sscanf(line, "%s%le%le%le%le",  symb, x, y, z, charge) != 5) {
    return((int)-3);
  }
  else {
    return((int)0);
  }

  case 7: if(sscanf(line, "%s%le%le%le%le%le%le",  symb, x, y, z, mx, my, mz) != 7) {
    return((int)-3);
  }
  else {
    return((int)0);
  }

  case 8: if(sscanf(line, "%s%le%le%le%le%le%le%le",
		    symb, x, y, z, charge, mx, my, mz) != 8) {
    return((int)-3);
  }
  else {
    return((int)0);
  }

  default: return((int)-3);

  } /* end switch */
}

/*=============================================
 write one line of XYZ xmol coordinate file.
n_fields - number of columns in the XYZ file of Xmol.
               4     At.Symb  X  Y  Z
               5     At.Symb  X  Y  Z Charge
               7     At.Symb  X  Y  Z  X-amplitude  Y-amplitude  Z-amplitude
               8     At.Symb  X  Y  Z  Charge X-ampl.  Y-ampl.  Z-ampl.
===============================================*/

int write_line(FILE *file,int  n_fields,char * symb, double x, double y, double z, double charge, double mx, double my, double mz) {

  switch (n_fields) {

  case 4: fprintf(file, "%10s%11.6f%11.6f%11.6f\n", symb, x, y, z);
    break;

  case 5: fprintf(file, "%10s%11.6f%11.6f%11.6f%11.6f\n", symb, x, y, z, charge);
    break;

  case 7: fprintf(file, "%10s%11.6f%11.6f%11.6f%11.6f%11.6f%11.6f\n",
                symb, x, y, z, mx, my, mz);
  break;

  case 8: fprintf(file, "%10s%11.6f%11.6f%11.6f%11.6f%11.6f%11.6f%11.6f\n",
		  symb, x, y, z, charge, mx, my, mz);
  break;

  default: return((int)-3);

  } /* end switch */

  return((int)0);

}

/*=======================================================*/


/* FITEST
 rigid fit test driver
 reads in data, fits, and writes out
*/

//  int old_main(argc, argv)
//  int argc;
//  char *argv[];
//  {
//   int n_fields_r;                /* no of fields in xmol file for ref. molec */
//   int n_fields_f;                /* no of fields in xmol file for fit. molec */
//   int nat_r;                     /* number of all atoms in reference molecule */
//   int nat_f;                     /* number of all atoms in fitted molecule */
//   char title_r[MAXLINELEN];      /* title line in reference mol. xmol file */
//   char title_f[MAXLINELEN];      /* title line in fitted mol. xmol file */
//   char symb_r[MAXPOINTS][10];    /* atom type symbols for reference molecule */
//   char symb_f[MAXPOINTS][10];    /* atom type symbols for fitted molecule */
//   char line[MAXLINELEN];         /* scratch space for line */
//   double xyz_r[4][MAXPOINTS];    /* coordinates for reference molecule */
//   double xyz_f[4][MAXPOINTS];    /* coordinates for fitted molecule */
//   double charge_r[MAXPOINTS];    /* charges for reference molecule */
//   double charge_f[MAXPOINTS];    /* charges for fitted molecule */
//   double modes_r[4][MAXPOINTS];  /* normal modes for reference */
//   double modes_f[4][MAXPOINTS];  /* normal modes for fitted */
//   int npairs;                    /* no of fitted atom pairs */
//   int atoms_r[MAXPOINTS];        /* atoms of ref. molecule to be superimposed */
//   int atoms_f[MAXPOINTS];        /* atoms of fit. molecule to be superimposed */
//   double ref_xyz[4][MAXPOINTS];  /* ref. molecule atom coordinates to fit */
//   double fit_xyz[4][MAXPOINTS];  /* fit. molecule atom coordinates to fit */
//   double weight[MAXPOINTS];      /* fitted atom pair weights */
//   double ref_center[4];          /* center of ref. molecule fitted atoms */
//   double fit_center[4];          /* center of ref. molecule fitted atoms */
//   double q[4];                   /* quaternion */
//   double u[4][4];                /* left rotation matrix for coordinates */
//   int i, j, ch;                  /* aux variables */
//   double s, d, wd, rms, wnorm, dotm;   /* aux variables */
//   FILE *refmol;                  /* file variable for reference mol xmol file */
//   FILE *fitmol;                  /* file variable for reference mol xmol file */
//   FILE *pairinp;                 /* file with fitted atom pairs and weights */
//   FILE *outfile;                 /* rotated and translated fit. mol. xmol file*/
//   FILE *statfile;                /* file with fit goodness values */
//   int max_sweeps;                /* max number of iterations in jacobi */
//   int opt;                       /* option letter */
//   int read_pairs;                /* 1 - read pairs, 0 - make pairs, weights */
//   static char usage[]=
//    "Usage : quatfit [-r ref] [-f fit] [-p pairs] [-o out] [-s stat]\n";
//   extern char *optarg;           /* option argument from getopt */

//   /* set defaults */
//   max_sweeps = 30;
//   read_pairs = 1;
//   refmol = stdin;
//   fitmol = stdin;
//   pairinp = stdin;
//   outfile = stdout;
//   statfile = stdout;

//   while ((opt = getopt(argc, argv, OPTIONS)) != EOF) {
//     switch (opt) {
//       case 'r':
//         if((refmol = fopen(optarg, "r")) == NULL)  {
//           fprintf(stderr,"Error: Could not find ref-mol-file: %s\n", optarg);
//           return(1);
//           }
//         break;

//       case 'f':
//         if((fitmol = fopen(optarg, "r")) == NULL)  {
//           fprintf(stderr,"Error: Could not find fit-mol-file: %s\n", optarg);
//           return(2);
//           }
//         break;

//       case 'p':
//         if(strcmp(optarg,"none") == 0) { /* if argument is "none" */
//           read_pairs = 0;
//           break;
//           }
//         if((pairinp = fopen(optarg, "r")) == NULL)  {
//           fprintf(stderr,
//                   "Error: Could not find pairs-wieight-file: %s\n", optarg);
//           return(3);
//           }
//         break;

//       case 's':
//         if((statfile = fopen(optarg, "w")) == NULL) {
//           fprintf(stderr, "Error: Could not open out-file %s\n", optarg);
//           return(4);
//           }
//         break;

//       case 'o':
//         if((outfile = fopen(optarg, "w")) == NULL)  {
//           fprintf(stderr,
//                   "Error: Could not open fitted-mol-out-file: %s\n", optarg);
//           return(5);
//           }
//         break;

//       case '?':
//         fprintf(stderr,"Error: %s\n", usage);
//         return(6);
//       }  /* end switch */
//     } /* end while */


//   /* Now read in the ref molecule */
//   while (fgets(line, MAXLINELEN, refmol) != NULL) { /* scan until n_at found */
//     i = 0;
//     while ((ch = line[i++]) != '\0') {
//       if(!isspace(ch)) {
//         i = -1;
//         break;
//         }
//       }
//     if(i == -1) {
//       break;
//       }
//     }

//   if(i != -1) { /* if white space only to the end */
//    fprintf(stderr, "Error: Error in line 1 of refmol file\n");
//    return(7);
//    }

//   if(sscanf(line, "%d\n", &nat_r) != 1) {  /* n  atoms */
//     fprintf(stderr, "Error: Error in line 1 of refmol file.\n");
//     return(8);
//     }

//   if(nat_r >= MAXPOINTS) {
//     fprintf(stderr,
//     "Error: Molecule too big. Recompile program with larger MAXPOINTS\n");
//     return(8);
//     }

//   if(fgets(title_r,  MAXLINELEN, refmol) == NULL) {  /* title line */
//     fprintf(stderr, "Error: Error in line 2 of refmol file.\n");
//     return(9);
//     }

//   /* read coordinates of ref molecule */
//   n_fields_r = 0;
//   for (i = 1; i <=  nat_r; i++) {
//     if(scan_line(refmol, &n_fields_r, &symb_r[i][0], &xyz_r[1][i], &xyz_r[2][i],
//                  &xyz_r[3][i], &charge_r[i], &modes_r[1][i], &modes_r[2][i],
//                  &modes_r[3][i]) != 0) {
//       fprintf(stderr,"Error: Error in line %d of refmol file.\n", i+2);
//       return(10);
//       }
//     }

//   /* Now read in the fitted molecule */
//   while (fgets(line, MAXLINELEN, fitmol) != NULL) { /* scan until n_at found */
//     i = 0;
//     while ((ch = line[i++]) != '\0') {
//       if(!isspace(ch)) {
//         i = -1;
//         break;
//         }
//       }
//     if(i == -1) {
//       break;
//       }
//     }

//   if(i != -1) { /* if white space only to the end */
//    fprintf(stderr, "Error: Error in line 1 of fitmol file\n");
//    return(11);
//    }

//   if(sscanf(line, "%d\n", &nat_f) != 1) {
//     fprintf(stderr, "Error: Error in line 1 of fitmol file.\n");
//     return(12);
//     }

//   if(nat_f >= MAXPOINTS) {
//     fprintf(stderr,
//     "Error: Molecule too big. Recompile program with larger MAXPOINTS\n");
//     return(13);
//     }

//   if(fgets(title_f,  MAXLINELEN, fitmol) == NULL) {
//     fprintf(stderr, "Error: Error in line 2 of fitmol file\n");
//     return(13);
//     }

//   /* read coordinates of fitted molecule */
//   n_fields_f = 0;
//   for (i = 1; i <=  nat_f; i++) {
//     if(scan_line(fitmol, &n_fields_f, &symb_f[i][0], &xyz_f[1][i], &xyz_f[2][i],
//                  &xyz_f[3][i], &charge_f[i], &modes_f[1][i], &modes_f[2][i],
//                  &modes_f[3][i]) != 0) {
//       fprintf(stderr,"Error: Error in line %d of fitmol file.\n", i+2);
//       return(14);
//       }
//     }

//   /* Read or set pairs and weights */
//   if(read_pairs == 1) {
//     while (fgets(line, MAXLINELEN, pairinp) != NULL) { /* skip white space */
//       i = 0;
//       while ((ch = line[i++]) != '\0') {
//         if(!isspace(ch)) {
//           i = -1;
//           break;
//           }
//         }
//       if(i == -1) {
//         break;
//         }
//       }

//     if(i != -1) { /* if white space only to the end */
//       fprintf(stderr, "Error: Error in line 1 of pairs and weights file\n");
//       return(15);
//       }

//     if(sscanf(line, "%d", &npairs) != 1) {
//       fprintf(stderr,
//               "Error: Error in line 1 of file with pairs and weights\n");
//       return(16);
//       }

//     if(npairs < 2) {
//       fprintf(stderr,
//               "Error: Cannot fit a single atom. Need at least 2\n");
//       return(16);
//       }

//     for (i = 1; i <= npairs; i++) {
//       if(fscanf(pairinp,"%d%d%le",&atoms_r[i], &atoms_f[i], &weight[i]) != 3) {
//         fprintf(stderr,
//                 "Error: Error in line %d of pairs and weights file.\n", i+2);
//         return(17);
//         }
//       if((atoms_r[i] < 1) || (atoms_r[i] > nat_r) || (atoms_f[i] < 1) ||
//          (atoms_f[i] > nat_f) ) {
//         fprintf(stderr,
//         "Error: Error in line %d of pairs and weights (number out of range)\n",
//         i+2);
//         return(18);
//         }
//       }
//     }
//   else {  /* initialize pairs to consectutive numbers */
//     npairs = nat_r;
//     if(nat_r > nat_f) {
//       npairs = nat_f;
//       }

//     for(i = 1; i <= npairs; i++) {
//       atoms_r[i] = i;
//       atoms_f[i] = i;
//       weight[i] = 1.0;
//       }
//     }

//   /* extract fitted atoms to tables */
//   for (i = 1; i <= npairs; i++) {
//     for (j = 1; j <= 3; j++) {
//       ref_xyz[j][i] = xyz_r[j][atoms_r[i]];
//       fit_xyz[j][i] = xyz_f[j][atoms_f[i]];
//       }
//     }

//   /* ===  Atom coordinates are fit in both modes === */
//   /* center ref molecule fitted atoms around (0,0,0) */
//   center (npairs, ref_xyz, weight, 1, ref_center);

//   /* center fitted molecule fitted atoms around (0,0,0) */
//   center (npairs, fit_xyz, weight, 1, fit_center);

//   /* fit specified atoms of fit_molecule to those of ref_molecule */
//   qtrfit(npairs, fit_xyz, ref_xyz, weight, q, u, max_sweeps);

//   /* subtract coordinates of the center of fitted atoms of the fitted molecule
//      from all atom coordinates of the fitted molecule (note that weight is
//      a dummy parameter) */
//   center(nat_f, xyz_f, weight, 2, fit_center);

//   /* rotate the fitted molecule by the rotation matrix u */
//   rotmol(nat_f, xyz_f, xyz_f, u);
//   /* same with set of fitted atoms of the fitted molecule */
//   rotmol(npairs, fit_xyz, fit_xyz, u);

//   /* if modes given in fitted molecule, rotate the modes too */
//   if(n_fields_f > 4) {
//     rotmol(nat_f, modes_f, modes_f, u);

//     /* calculate dot product of reference and fitted molecule modes */
//     if(n_fields_r > 4) {
//       dotm = 0.0;
//       for (i = 1; i <= npairs; i++) {
//         for (j = 1; j <= 3; j++) {
//           dotm += modes_r[j][atoms_r[i]]*modes_f[j][atoms_f[i]];
//         }
//       }
//     }
//   }

//   /* translate atoms of the fitted molecule to the center
//        of fitted atoms of the reference molecule */
//   center(nat_f, xyz_f, weight, 3, ref_center);
//   /* same with set of fitted atoms of the fitted molecule */
//   center(npairs, fit_xyz, weight, 3, ref_center);
//   /* translate fitted atoms of reference molecule to their orig. location */
//   center(npairs, ref_xyz, weight, 3, ref_center);

//   /* write modified XYZ file for fitted molecule */
//   fprintf(outfile,"%d\n%s", nat_f, title_f);
//   for (i = 1; i <= nat_f; i++) {
//     write_line(outfile, n_fields_f, symb_f[i], xyz_f[1][i], xyz_f[2][i],
//                xyz_f[3][i], charge_f[i], modes_f[1][i], modes_f[2][i],
//                modes_f[3][i]);
//     }

//   fflush(outfile);
//   /* find distances between fitted and reference atoms and print them in
//      out file */

//   fprintf(statfile,
//           "\nDistances and weighted distances between fitted atoms\n");
//   fprintf(statfile,"Ref.At. Fit.At.  Distance  Dist*sqrt(weight)  weight\n");
//   rms = 0.0;
//   wnorm = 0.0;
//   for (i = 1; i <= npairs; i++) {
//     d = 0.0;
//     for (j = 1; j <= 3; j++) {
//       s = ref_xyz[j][i] - fit_xyz[j][i];
//       d += s*s;
//       }
//     d = sqrt(d);
//     s = sqrt(weight[i]);
//     wd = s*d;
//     rms += wd*wd;
//     wnorm += s;
//     fprintf(statfile, "  %3d    %3d  %11.6f   %11.6f    %11.6f\n",
//     atoms_r[i], atoms_f[i], d, wd, weight[i]);
//     }

//   rms = sqrt(rms/wnorm);
//   fprintf(statfile, "\n\nWeighted root mean square=%10.6f\n\n", rms);
//   fprintf(statfile, "\n\nCenter of reference molecule fitted atoms\n");
//   fprintf(statfile, "Xc = %11.6f Yc = %11.6f Zc = %11.6f\n",
//           ref_center[1], ref_center[2], ref_center[3]);

//   fprintf(statfile, "\n\nCenter of fitted molecule fitted atoms\n");
//   fprintf(statfile, "Xc = %11.6f Yc = %11.6f Zc = %11.6f\n",
//           fit_center[1], fit_center[2], fit_center[3]);

//   fprintf(statfile,"\n\nLeft rotation matrix\n");
//   for (i = 1; i <= 3; i++) {
//     fprintf(statfile, " %11.6f  %11.6f  %11.6f\n",
//             u[1][i], u[2][i], u[3][i]);
//     }

//   if((n_fields_f > 4) && (n_fields_r > 4)) {
//     fprintf(statfile,
//     "\nDot product of normal modes on fitted atom pairs =%11.6f\n", dotm);
//     }

//   fclose(outfile);
//   fclose(statfile);
//   fclose(pairinp);
//   fclose(fitmol);
//   fclose(refmol);

//   return(0);

//   }



