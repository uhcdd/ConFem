%module LinAlg2
%{
#define SWIG_FILE_WITH_INIT
#include "LinAlg2.h"
extern int n;
%}

%include "numpy.i"

%init %{
import_array();
%}

%apply (double* IN_ARRAY1    , int DIM1) {(double* a,  int dim_a),
					                      (double* b,  int dim_b),
					                      (double* r,  int dim_r),
					                      (double* le, int dim_le),
					                      (double* ri, int dim_ri),
										  (double* KVecU, int dim_KVecU),
										  (double* KVecL, int dim_KVecL),
										  (double* ub, int dim_ub),
										  (double* uub, int dim_uub),
										  (double* pp, int dim_pp),
										  (double* pp0, int dim_pp0)};
//										  (double* A, int dim_A),
//										  (double *DataOut, int dim_DataOut)};
%apply (int* IN_ARRAY1    , int DIM1)    {(int* b2, int dim_b2),
										  (int* Skyline, int dim_Skyline)};
%apply (long long* IN_ARRAY1, int DIM1)  {(long long *b1, int dim_b1),
										  (long long *SDiag, int dim_SDiag)};
//%apply (double* IN_PLACE_ARRAY1, int DIM1) {(double* r, int dim_r)}; 

%include "LinAlg2.h"

