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
 */
#ifndef QUATFIT_H
#define QUATFIT_H
//
// Header files
//
#include <vector>
#include <stdio.h>
#include <math.h>
//
// Local header files
//
#include "fff.h"
#include "atom_class.h"
#include "vector3d.h"
//
// MAXPOINTS is the max number of points that can be transformed in a quatfit operation
// 
#define MAXPOINTS     4000
#define MAXLINELEN    250
#define OPTIONS       "r:f:p:o:s:"

void center (int n, double x[4][MAXPOINTS],double w[MAXPOINTS], int io, double o[4]);
void rotmol (int n, double x[4][MAXPOINTS], double y[4][MAXPOINTS], double u[4][4]);
void jacobi (double a[4][4], double d[4], double v[4][4], int nrot);
void q2mat (double q[4], double u[4][4]);
void qtrfit (int n, double x[4][MAXPOINTS], double y[4][MAXPOINTS], double w[MAXPOINTS], double q[4], double u[4][4], int nr);
int scan_line(FILE *file, int *n_fields, char *symb, double *x, double *y, double *z, double *charge, double *mx, double *my, double *mz);
int write_line(FILE *file,int  n_fields,char * symb, double x, double y, double z, double charge, double mx, double my, double mz);

class quatfit_class {
 public:
  vector<double> tvector;
  vector<double> rotmatrix;
  quatfit_class();
  //~quatfit_class();
  void fit(int numpoints,vector<atom_class> refpoints,vector<atom_class> fitpoints);
  void transform_all(FFF &Protein);

  void transform(int numpoints,vector<atom_class> & fitcoords);
  void rotate_around_vector(vector3d V,double angle);
  //
  //
  //
  double fit_center[4];
  double ref_center[4];
  //
  // Flags
  //
  bool center_on_transform;
  //
  // Private functions
  //
private:
  void _transform(vector<atom_class> & fitcoords);
};

//
// Convert degrees to radians
//
inline float deg2rad(double angle) {
  return angle*atan(1.0)/45.0;
}


#endif
