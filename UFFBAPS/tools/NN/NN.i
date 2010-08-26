%module NN

%{
#include "NN.h"

%}



%include "std_vector.i"
%include "std_string.i"
// Instantiate templates used by example
namespace std {
  %template(IntVector) vector<int>;
  %template(DoubleVector) vector<double>;
  %template(FloatVector) vector<float>;
  %template(FloatVectorVector) vector<vector<float> >;
}

%include "NN.h"

