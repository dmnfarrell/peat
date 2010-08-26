#ifndef ALGEBRA
#define ALGEBRA
/*
 #
 # UFFBAPS - Unified Force Field for Binding And Protein Stability
 # Copyright (C) 2010 Jens Erik Nielsen & Chresten Soendergaard
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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string>
#include <vector>

using namespace std;

struct Vector {
  
  float x,y,z;

  Vector();
  
  Vector(float xi, float yi, float zi);


  void normalise();

  Vector operator-(Vector other);
  Vector operator+(Vector other);
  void operator=(Vector other);
  Vector operator*(float constant);

  float operator*(Vector other); //dot product
  Vector operator^(Vector other); //cross product

  float length();
  Vector get_perpendicular();
  void rotate(Vector axis, float angle);
  void rotate(float angle, Vector p1,Vector p2);

  void display(string name);
  void test();

};


struct Matrix4x4 {


  float a[4][4];

  Matrix4x4(float a11, float a12, float a13, float a14,
	    float a21, float a22, float a23, float a24,
	    float a31, float a32, float a33, float a34,
	    float a41, float a42, float a43, float a44);
  




};

Vector mult4(Vector v, Matrix4x4 m);


struct Sigmoidal {

  Sigmoidal();
  float value(float desolvation_factor);
  float calculate(float x);

  vector<float> precalcs;
  float offset, factor, step;
  int index;

};

#endif
