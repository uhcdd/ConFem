%module ConFemElemC
%{
	/* the resulting C file should be built as a python extension */
	#define SWIG_FILE_WITH_INIT
	/*  Includes the header in the wrapper code */
	#include "ConFemElemC.h"
%}

/*  include the numpy typemaps */
%include "numpy.i"
/*  need this for correct module initialization */
%init %{
	import_array();
%}

%apply (double* IN_ARRAY1 , int DIM1) {(double* det_, int dim_det_),
									   (double* Data, int dim_D),
									   (double* a, int dim_a),
									   (double* EdgeDir, int dim_Ed),
									   (double* sig, int dim_sig),
									   (double* rvec, int dim_rvec)};

%apply (double* IN_ARRAY2 , int DIM1, int DIM2) {(double* BC,  int dim_BC0, int dim_BC1),
												 (double* XX,  int dim_X0,  int dim_X1),
												 (double* Vn,  int dim_Vn0, int dim_Vn1),
												 (double* gg0, int dim_g00, int dim_g01),
												 (double* gg1, int dim_g10, int dim_g11),
												 (double* BB,  int dim_BB0, int dim_BB1),
												 (double* TD,  int dim_TD0, int dim_TD1),
												 (double* GeomK, int dim_GeomK0, int dim_GeomK1),
												 (double* kmat,int dim_kmat0, int dim_kmat1),
												 (double* CC,  int dim_CC0, int dim_CC1),
												 (double* X,   int dim_X0, int dim_X1)};

/*  Parse the header file to generate wrappers */
%include "ConFemElemC.h"

