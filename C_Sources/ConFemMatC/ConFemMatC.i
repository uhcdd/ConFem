%module ConFemMatC
%{
	/* the resulting C file should be built as a python extension */
	#define SWIG_FILE_WITH_INIT
	/*  Includes the header in the wrapper code */
	#include "ConFemMatC.h"
	/* ??? */
/*	extern int ElemDim; */
/*	extern double Emod, nu; */
/*	extern bool ElemPlSt; */
%}

/*  include the numpy typemaps */
%include "numpy.i"
/*  need this for correct module initialization */
%init %{
	import_array();
%}

%apply (double* IN_ARRAY1 , int DIM1) {(double* Eps, int dim_E),
                                       (double* ElemStateVarN, int dim_StN),
									   (double* ElemStateVar, int dim_St),
									   (double* ElemDataP, int dim_DaP),
									   (double* MatM, int dim_Ma),
									   (double* sig, int dim_si),
									   (double* ww, int dim_ww),
									   (double* EpsR, int dim_EpsR),
									   (double* Dps, int dim_D),
									   (double* sigR, int dim_sigR),
									   (double* CR, int dim_CR),
									   (double* DataOut, int dim_DataOut),
									   (double* la, int dim_la),
									   (double* vv, int dim_vv)};

%apply (double* IN_ARRAY2 , int DIM1, int DIM2) {(double* MatM2,  int dim_Ma21, int dim_Ma22)};

/*  Parse the header file to generate wrappers */
%include "ConFemMatC.h"

