// LinAlg.cpp : Defines the entry point for the console application.
//


#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <fstream>
#include <iostream>
#include <iomanip>

#include "LinAlg2.h"
using namespace std;


//*******************************************************************
//
//  triangularization with gauss
//  Sparse Matrix A,B
//  unsymmetrisch
//  long n        dimension of matrix    
//  double *a     pointer to start address of stiffness vector, upper right part including diagonal, 
//                                        here columns are stored in a sequence
//  double *b     pointer to start address of stiffness vector, lower left part without diagonal
//                                        here rows are stored in a sequence 
//  double *r     pointer to start address of right side vector, overwritten with solution
//  int *b1       pointer to start address of vector with indices for diagonalelements  +1
//  int *b2       pointer to vector with skyline heights of a resp. bandwidths of b
//
//*******************************************************************
void sim0_so(double *r, int dim_r, double *a, int dim_a, double *b, int dim_b, long long *b1, int dim_b1, int *b2, int dim_b2, int n)
{
   int i,k,j,band,diff,end;
   long long atop, diag;
   // modification of right side
   for( k=0 ; k<n-1 ; k++ ) {
      for( i=k+1 ; i<n ; i++ ) {
		 band = *(b2+i);
		 diff = i-k;
	     if( diff < band ) {
//			atop = *(b1+i)+i-k-1;
         	atop = *(b1+i)+i-k;    // uhc
            *(r+i) = *(r+i) - *(b+atop) * *(r+k);
		 }
	  }
   }
    // backsubstitution
   for( k=n-1 ; k>=0 ; k-- ) {
//      diag = *(b1+k)-1;
      diag = *(b1+k);     // uhc
      *(r+k) = *(r+k) / *(a+diag);
	  band = *(b2+k);
	  end = k - band;
      for( i=k-1,j=1 ; i>end ; i--,j++ )
         *(r+i) = *(r+i) - *(a+diag+j) * *(r+k); 
   }
//    cout << r[0] << endl;
}  
void sim0_lu( double *a, int dim_a, double* b, int dim_b, long long *b1, int dim_b1, int* b2, int dim_b2, int n)
{
   int i,k,band;
   long long j, diag, radr, diff, atop;
   register int jj, *pb2;
   long long *pb1;
   double pivot;
   double tol = 1.e-9; //GRNULL;
   register double value;
   
   for( k=0 ; k<n-1 ; k++ ) {
//    diag = *(b1+k)-1;
      diag = *(b1+k);  // uhc
//       if( fabs( *(a+diag) ) < tol ) { printf("solve::sim0_lu-1"); exit(0); }
//	  cout << k << "     " << *(a+diag) << endl;
	  pivot = 1./(*(a+diag)) ;
	  pb1 = b1 + k;
	  pb2 = b2 + k;

      for( j=k+1 ; j<n ; j++ ) {  
	     band = *(b2+j);
		 diff = j-k;
	     if( diff < band ) {
//          diag = *(b1+j)-1;
            diag = *(b1+j);  // uhc
			radr = diag+diff;
	        *(b+radr) = *(b+radr) * pivot; 
      }  }
 	  // modification lower left part p
      for( i=1 ; i<n-k ; i++ ) {
	     if( i <   *(pb2+i) ) { 
//			radr = *(pb1+i)+i-1;
			radr = *(pb1+i)+i;   // uhc
			value= *(b+radr);
//			for( j=*(pb1+i),jj=i-1; j<radr; j++,jj-- ) { 
			for( j=*(pb1+i)+1,jj=i-1; j<radr; j++,jj-- ) {   // uhc 
			   if( jj   < *(pb2+jj) ) {    
//			      atop  = *(pb1+jj)+jj-1;
			      atop  = *(pb1+jj)+jj; // uhc
	              *(b+j)= *(b+j) - value * *(a+atop);
      }  }  }  }
	  // modification of upper right part a 
      for( i=1 ; i<n-k ; i++ ) {
	     if( i <   *(pb2+i) ) {   
//			radr = *(pb1+i)+i-1; 
			radr = *(pb1+i)+i;  // uhc 
			value= *(a+radr);
//			for( j=*(pb1+i)-1,jj=i; j<radr; j++,jj-- ) {     
			for( j=*(pb1+i),jj=i; j<radr; j++,jj-- ) {   // uhc     
			   if( jj   < *(pb2+jj) ) {                     
//			      atop  = *(pb1+jj)+jj-1;                   
			      atop  = *(pb1+jj)+jj;  // uhc                   
	              *(a+j)= *(a+j) - *(b+atop) * value;
	  }  }  }  }
//    diag  = *(b1+n-1)-1;
      diag  = *(b1+n-1);
	  value = *(a+diag);
//       if( fabs(value) < tol ) { printf("solve::sim0_lu-2"); exit(0); }
   }
}
void sim0_mmul( double *le, int dim_le, double *a, int dim_a, double* b, int dim_b, double *ri, int dim_ri, long long *b1, int dim_b1, int* b2, int dim_b2, int n)
{
   int i,k,band,diff;
   long long diag;
   // loop over all rows of matrix and vector
   for( k=0 ; k<n ; k++ ) 
//	k=2;
   {
       // lower left part
	   *(le+k) = 0.;
	   diag = *(b1+k);
	   band = *(b2+k);
	   for( i=0 ; i<band; i++ ) {
		   diff = k-i;
		   *(le+k) = *(le+k) + *(b+diag+i) * *(ri+k-i); 
//                cout << "XX" << endl;
//             cout << i  <<"  " << band<< "  " << diff << "  " << *(b+diag+i) << "  " << *(ri+k-i) << " :  " << *(le+k) << endl; 
       }
       // upper right part
       for( i=k; i<n ; i++ ) {
	     diag = *(b1+i);
		 band = *(b2+i);
		 diff = i-k; 
	     if( diff < band ) {
            *(le+k) = *(le+k) + *(a+diag+diff) * *(ri+i);
 //               cout << "YY" << endl;
 //             cout << i  <<"  " << band<< "  " << diff << endl; 
		 }
	  }
   }
}

//*******************************************************************
//
//  Matrix vector triangularization with gauss
//  Sparse Matrix A
//  symmetrisch
//  long n        dimension of matrix    
//  double *a     pointer to start address of stiffness vector, upper right part including diagonal, 
//                                        here columns are stored in a sequence
//  double *r     pointer to start address of right side vector
//  int *b1       pointer to start address of vector with indices for diagonalelements  +1
//  int *b2       pointer to vector with skyline heights of a resp. bandwidths of b
//
//*******************************************************************
void sim1_so(double *r, int dim_r, double *a, int dim_a, long long *b1, int dim_b1, int *b2, int dim_b2, int n)
{
   int i,k,j,band,diff,end;
   long long diag, atop;
   double pivot; 
   // modification of right side 
   for( k=0 ; k<n-1 ; k++ ) { 
//    diag = *(b1+k)-1;
      diag = *(b1+k);                      // uhc 
	  pivot = 1./(*(a+diag)) ;             // ! division by zero not catched
      for( i=k+1 ; i<n ; i++ ) {
         band = *(b2+i); 
         diff = i-k; 
         if( diff < band ) { 
//			atop = *(b1+i)+i-k-1;
			atop = *(b1+i)+i-k;            // uhc
			*(r+i) = *(r+i) - *(a+atop) * *(r+k) * pivot;
   }  }  }  
   // backsubstitution
   for( k=n-1 ; k>=0 ; k-- ) {
//      diag = *(b1+k)-1;
        diag = *(b1+k);                    // uhc
		*(r+k) = *(r+k) / *(a+diag);
		band = *(b2+k);
		end = k - band;
		for( i=k-1,j=1 ; i>end ; i--,j++ )
			*(r+i) = *(r+i) - *(a+diag+j) * *(r+k);
   }
}
void sim1_lu( double *a, int dim_a, long long *b1, int dim_b1, int* b2, int dim_b2, int n)
{
// int i,j,k,band,diag,diff,radr;
   int k;
   long long diag, i, j, radr, jj, atop;
   register int *pb2;
   long long *pb1;
   double pivot;
   double tol = 1.e-9;
   register double value;
   for( k=0 ; k<n-1 ; k++ )
   {
//    diag = *(b1+k)-1;
      diag = *(b1+k);                              // uhc
      pivot = 1./(*(a+diag)) ;                     // ! not catched for division by zero
      pb1 = b1 + k;
      pb2 = b2 + k;
	  // modification of upper right part 
      for( i=1 ; i<n-k ; i++ ) {
         if( i < *(pb2+i) ) {
//			radr = *(pb1+i)+i-1; 
			radr = *(pb1+i)+i;                     // uhc 
            value= *(a+radr) * pivot;
//			for( j=*(pb1+i)-1,jj=i; j<radr; j++,jj-- ) {     
			for( j=*(pb1+i),  jj=i; j<radr; j++,jj-- ) {   // uhc      
               if( jj   < *(pb2+jj) ) {
//			      atop  = *(pb1+jj)+jj-1;                   
			      atop  = *(pb1+jj)+jj;            // uhc                   
                  *(a+j)= *(a+j) - *(a+atop) * value;
      }  }  }  }
   }
}
//void sim1_mmul( double *le, int dim_le, double *a, int dim_a, double* b, int dim_b, double *ri, int dim_ri, long long *b1, int dim_b1, int* b2, int dim_b2, int n)
void   sim1_mmul( double *le, int dim_le, double *a, int dim_a,                       double *ri, int dim_ri, long long *b1, int dim_b1, int* b2, int dim_b2, int n)
{
   int i,k,band,diff;
   long long diag;
   // loop over all rows of matrix and vector
   for( k=0 ; k<n ; k++ ) 
   {
       // lower left part of matrix multiplication whereby using matrix upper right values
	   *(le+k) = 0.;
	   diag = *(b1+k);
	   band = *(b2+k);
//	   for( i=0 ; i<band; i++ ) {
	   for( i=1 ; i<band; i++ ) {
		   diff = k-i;
//		   *(le+k) = *(le+k) + *(b+diag+i) * *(ri+k-i); 
		   *(le+k) = *(le+k) + *(a+diag+i) * *(ri+k-i); 
       }
       // upper right part
       for( i=k; i<n ; i++ ) {
	     diag = *(b1+i);
		 band = *(b2+i);
		 diff = i-k; 
	     if( diff < band ) {
            *(le+k) = *(le+k) + *(a+diag+diff) * *(ri+i);
		 }
	  }
   }
}

// luPAR from LUStandAlone
void sim1_luP(double *a, int dim_a, long long *b1, int dim_b1, int* b2, int dim_b2, int n, int NK)
{
	int i, k, *UpRow;
	double pivot, *AR;
	UpRow = (int *)malloc((n) * sizeof(int));
	AR = (double *)malloc((n) * sizeof(double));									// holds rows / columns of previous LU-decomposition step
	for (i = 0; i < n; i++) *(AR + i) = 0.;
	for (i = 1; i < n; i++) *(UpRow + i) = i - *(b2 + i) + 1;    					// Upper row (with lowest index) of respective column

	for (k = 1; k < n; k++) {														// loop over submatrices
		pivot = 1. / (*(a + *(b1 + k - 1)));
		for (int ii = k; ii < n; ii++)												// loop over columns
		{
			if (*(UpRow + ii) < k) {												// see Gauss.py											
				*(AR + n - ii) = *(a + *(b1 + ii) + ii - k + 1); 					// a[k-1,ii] --> a[ii,k] --> diag[ii] + ii - (k-1)
			}
		}
		for (int ii = k; ii < n; ii++)	 											// loop over columns
		{
			if (*(UpRow + ii) < k)													// current row k is below skyline in column ii / has larger index
			{
				double value = -pivot * *(AR + n - ii);
				double *aa = a + *(b1 + ii) - k;									// base (botom of skyline) for a-addresses in column i
				double *ra = AR + n - k - ii;
				for (int jj = k; jj < (ii + 1); jj++) 								// loop over rows
				{
					*(aa + jj) = *(aa + jj) + value * *(ra + jj);
				}
			}
		}
	}
}

//*******************************************************************
//
//  from CaeFem3
//
//****************
void BoundFinishLC( int N, int k, long long *SDiag, int dim_SDiag, int *Skyline, int dim_Skyline, double *KVecU, int dim_KVecU, int SymSys, double *KVecL, int dim_KVecL)
{
	int j, diff;
	long long jj, ind;
	jj = *(SDiag+k);							//	jj = SDiag[k]
//	modify values connected with system matrix column
	for (j=0; j<*(Skyline + k); j++) {
		*(KVecU + jj + j) = 0.;					//	for j in range(Skyline[k]) : KVecU[jj + j] = 0.		   # upper right part of system matrix
	}
	if (SymSys == 0) {							//	if not SymSys:
		for (j=0; j<*(Skyline + k); j++) {
			*(KVecL + jj + j) = 0.;				//	for j in range(Skyline[k]) : KVecL[jj + j] = 0.	   # lower left part of system matrix
		}
	}
//	modify values connected with system matrix row-- "diff" within loops is not constant due to banded vector storage
	for (j = k; j<N; j++) {						//	for j in range(k, N) :
		diff = j - k;							//	diff = j - k
		if ( diff<*(Skyline+j) ) {
			ind = *(SDiag + j);
			*(KVecU + ind + diff) = 0.;			//	if diff<Skyline[j] : KVecU[SDiag[j] + diff] = 0.
		}
	}
	if (SymSys == 0) {							//	if not SymSys : # unsymmetric system
		for (j = k; j<N; j++) {					//	for j in range(k, N) :
			diff = j - k;						//	diff = j - k
			if (diff < *(Skyline + j)) {
				ind = *(SDiag + j);
				*(KVecL+ind+diff) = 0.;			//	if diff<Skyline[j] : KVecL[SDiag[j] + diff] = 0.
			}
		}
	}
	*(KVecU + jj) = 1.; 						//	KVecU[jj] = 1.
}
void RightHandSideLC(int N, int k, double LF, double LT, double *ub, int dim_ub, double *uub, int dim_uub, long long *SDiag, int dim_SDiag, int *Skyline, int dim_Skyline, \
	                 double *pp, int dim_pp, double *pp0, int dim_pp0, double *KVecU, int dim_KVecU, int SymSys, double *KVecL, int dim_KVecL)
//def RightHandSideL(      N,     k,        LF,        LT,         ub,                     uub,                   SDiag,                     Skyline,                                                                            KVecU,                    SymSys, KVecL) :
{
	int j, diff;
	long long jj, ind;
	double val, valT;
	val = LF * *(ub+k) - *(uub+k);				//	val = LF * ub[k] - uub[k]
	valT= LT * *(ub+k);							//	valT = LT * ub[k]
	jj = *(SDiag+k);							//	jj = SDiag[k]
//	 modify values connected with system matrix column
	for (j = 0; j < *(Skyline + k); j++) {		//	for j in range(Skyline[k]) :
		*(pp+k-j) = *(pp+k-j) - *(KVecU+jj+j) * val; //	pp[k - j] = pp[k - j] - KVecU[jj + j] * val
		*(pp0+k-j)= *(pp0+k-j)- *(KVecU+jj+j) * valT; //pp0[k - j] = pp0[k - j] - KVecU[jj + j] * valT
	}
//	modify values connected with system matrix row-- "diff" within loops unfortunately is not constant due to banded vector storage
	if (SymSys == 0) {							//			if not SymSys:										  	# unsymmetric system
		for (j = k; j < N; j++) {				//				for j in range(k, N) :
			diff = j - k;						//					diff = j - k
			if (diff < *(Skyline + j) ) {		//					if diff<Skyline[j] :
				ind = *(SDiag + j);
				*(pp+j) = *(pp+j) - *(KVecL+ind+diff) * val; //	pp[j] = pp[j] - KVecL[SDiag[j] + diff] * val
				*(pp0+j)= *(pp0+j)- *(KVecL+ind+diff) * valT;//pp0[j] = pp0[j] - KVecL[SDiag[j] + diff] * valT
			}
		}
	}
	else {										//					else : # use KVecU for rows instead KVecL in case of symmetric system
		for (j=k; j<N; j++) {					// for j in range(k, N):
			diff = j - k;						//	diff = j - k
			if (diff < *(Skyline+j)) {			//								if diff<Skyline[j] :
				ind = *(SDiag + j);
				*(pp+j) = *(pp+j) - *(KVecU+ind+diff)* val;	// pp[j] = pp[j] - KVecU[SDiag[j] + diff] * val
				*(pp0+j)= *(pp0+j)- *(KVecU+ind+diff)* valT;// pp0[j] = pp0[j] - KVecU[SDiag[j] + diff] * valT
			}
		}
	}
}


void RightHandSideLC_(int N, int k, long long *SDiag, int dim_SDiag, int *Skyline, int dim_Skyline, \
	double *pp, int dim_pp, double *pp0, int dim_pp0, double *KVecU, int dim_KVecU, int SymSys, double *KVecL, int dim_KVecL, double val, double valT)
	//def RightHandSideL(      N,     k,        LF,        LT,         ub,                     uub,                   SDiag,                     Skyline,                                                                            KVecU,                    SymSys, KVecL) :
{
	int j, diff;
	long long jj, ind;
//	double val, valT;
//	val = LF * *(ub + k) - *(uub + k);				//	val = LF * ub[k] - uub[k]
//	valT = LT * *(ub + k);							//	valT = LT * ub[k]
	jj = *(SDiag + k);							//	jj = SDiag[k]
												//	 modify values connected with system matrix column
	for (j = 0; j < *(Skyline + k); j++) {		//	for j in range(Skyline[k]) :
		*(pp + k - j) = *(pp + k - j) - *(KVecU + jj + j) * val; //	pp[k - j] = pp[k - j] - KVecU[jj + j] * val
		*(pp0 + k - j) = *(pp0 + k - j) - *(KVecU + jj + j) * valT; //pp0[k - j] = pp0[k - j] - KVecU[jj + j] * valT
	}
	//	modify values connected with system matrix row-- "diff" within loops unfortunately is not constant due to banded vector storage
	if (SymSys == 0) {							//			if not SymSys:										  	# unsymmetric system
		for (j = k; j < N; j++) {				//				for j in range(k, N) :
			diff = j - k;						//					diff = j - k
			if (diff < *(Skyline + j)) {		//					if diff<Skyline[j] :
				ind = *(SDiag + j);
				*(pp + j) = *(pp + j) - *(KVecL + ind + diff) * val; //	pp[j] = pp[j] - KVecL[SDiag[j] + diff] * val
				*(pp0 + j) = *(pp0 + j) - *(KVecL + ind + diff) * valT;//pp0[j] = pp0[j] - KVecL[SDiag[j] + diff] * valT
			}
		}
	}
	else {										//					else : # use KVecU for rows instead KVecL in case of symmetric system
		for (j = k; j<N; j++) {					// for j in range(k, N):
			diff = j - k;						//	diff = j - k
			if (diff < *(Skyline + j)) {			//								if diff<Skyline[j] :
				ind = *(SDiag + j);
				*(pp + j) = *(pp + j) - *(KVecU + ind + diff)* val;	// pp[j] = pp[j] - KVecU[SDiag[j] + diff] * val
				*(pp0 + j) = *(pp0 + j) - *(KVecU + ind + diff)* valT;// pp0[j] = pp0[j] - KVecU[SDiag[j] + diff] * valT
			}
		}
	}
}


/* already available in ConFemMatCx64
double EigenJacobiSymC(double *A, int NN, double *DataOut, int dim_DataOut) // from applied numerical methods 4.8, DataOut for control, N=3 used for the following
{
	const int itmax = 10, N = 3;					//	itmax = 10 # 10
	const double Eps1 = 1.0e-10, Eps2 = 1.0e-10, Eps3 = 1.0e-7; //
	int i, j, iter, k, ind;
	double T[3][3], AIK[3], EIGEN[3], sigma1, sigma2, OffDsq, S, Q, siQ, P, alpha, CSA,SNA,HOLDKI, val; // T = zeros((N, N), dtype = float), AIK = zeros((N), dtype = float), EIGEN = zeros((N), dtype = float),sigma1, OffDsq = 0., 0.
	sigma1 = 0.;
	OffDsq = 0.;
	for (i = 0; i < N; i++) {					//		for i in range(N) :
		for (j = 0; j < N; j++) {
			T[i][j] = 0.;						// 		T[i, i] = 1.0
		}
		T[i][i] = 1.;
		sigma1 += *(A+N*i+i) * *(A+N*i+i);		// 		sigma1 = sigma1 + A[i, i] * *2
		for (j = i + 1; j > N; j++) {			//		for j in range(i + 1, N) :
			OffDsq += *(A +N*i+j) * *(A+N*i+j); //		OffDsq = OffDsq + A[i, j] * *2      # considers symmetry
		}
	}
	S = 2.0*OffDsq + sigma1;					//		S = 2.0*OffDsq + sigma1                                     # this should not change during transformations
	if (sqrt(S) < Eps1) {						// 		if sqrt(S)<Eps1:
		*(A + 0) = 1.;							//		eV = zeros((N), dtype = float); eV[0] = 1.
		*(A + 1) = 0.;
		*(A + 2) = 0.;
		return(0.);								// 		return 0., eV							// 
	}
	// iteration loop
	for (iter = 0; iter < itmax; iter++) {		//		for iter in range(itmax) :
		for (i=0; i<N-1; i++) {					//for i in range(N - 1) :
			for (j=i+1; j<N; j++) {				//			for j in range(i + 1, N) :
				if (*(A + i * N + j) <= Eps2) { break; } //if abs(A[i, j]) <= Eps2 : break
				Q = abs(*(A+j*N+j) - *(A+i*N+i));//Q = abs(A[j, j] - A[i, i])                          # diagonal i, j difference
				// compute sine and cosine of rotation angle
				if (Q >= Eps1) {				//     if Q >= Eps1:                                     # see carnahan : applied numerical methods, p. 251
					siQ = Q/(*(A+i*N+i) - *(A+j*N+j));				// siQ = Q / (A[i, i] - A[j, j])
					P = 2. * *(A+i*N+j) * siQ;	//					P = 2.*A[i, j] * siQ
					alpha = P / Q;				// alpha = P / Q
					CSA = sqrt( 0.5*(1. + 1./sqrt(1.+    alpha*alpha)) );//CSA = sqrt(0.5*(1. + 1. / sqrt(1. + alpha * *2))) # this is obviously larger than sqrt(1 / 2)
					if (P > 0) { SNA = -1. / (2.*CSA* sqrt(1. + 1. / (alpha*alpha))); }			// SNA = sign(P) / (2.*CSA*sqrt(1. + 1. / alpha * *2))
					else       { SNA =  1. / (2.*CSA* sqrt(1. + 1. / (alpha*alpha))); }
				}
				else {							// else:
					CSA = 1. / sqrt(2.);		// CSA = 1. / sqrt(2.)
					SNA = CSA;					// SNA = CSA
				}
				// update columns i and j of T--> T_old * T
				for (k = 0; k < N; k++) {			//for k in range(N) :
					HOLDKI = T[k][i];			// HOLDKI = T[k, i]
					T[k][i] = HOLDKI * CSA + T[k][j] * SNA; //	T[k, i] = HOLDKI * CSA + T[k, j] * SNA
					T[k][j] = HOLDKI * SNA - T[k][j] * CSA; //	T[k, j] = HOLDKI * SNA - T[k, j] * CSA
				}
				// rows i and j--> transpose(T) * A
				for ( k=i; k<N; k++) {			//	for k in range(i, N) :
					if (k <= j) {					//if k <= j :
						AIK[k] = *(A+i*N+k);	// AIK[k] = A[i, k]
						*(A+i*N+k) = AIK[k] * CSA + *(A+k*N+j) * SNA; //	A[i, k] = AIK[k] * CSA + A[k, j] * SNA
						if (k == j)	{			//	if k == j : 
							*(A+j*N+k) = AIK[k] * SNA - *(A+j*N+k) * CSA;  // A[j, k] = AIK[k] * SNA - A[j, k] * CSA
						}
					}
					else {						//else :
						HOLDKI = *(A+i*N+k);	// HOLDKI = A[i, k]
						*(A+i*N+k) = HOLDKI * CSA + *(A+j*N+k) * SNA;  //A[i, k] = HOLDKI * CSA + A[j, k] * SNA
						*(A+j*N+k) = HOLDKI * SNA - *(A+j*N+k) * CSA;  //	A[j, k] = HOLDKI * SNA - A[j, k] * CSA
					}
				}
				AIK[j] = SNA * AIK[i] - CSA * AIK[j];  // AIK[j] = SNA * AIK[i] - CSA * AIK[j]
				// columns i and j--> A * T
				for (k = 0; k < (j + 1); k++) { //for k in range(j + 1) :
					if (k > i) { *(A+k*N+j) = SNA * AIK[k] - CSA * *(A+k*N+j); } //if k > i:  A[k, j] = SNA * AIK[k] - CSA * A[k, j]
					else {						//else :
						HOLDKI = *(A+k*N+i);	// HOLDKI = A[k, i]
						*(A+k*N+i) = HOLDKI * CSA + *(A+k*N+j) * SNA;  // A[k, i] = HOLDKI * CSA + A[k, j] * SNA
						*(A+k*N+j) = HOLDKI * SNA - *(A+k*N+j) * CSA;  // A[k, j] = HOLDKI * SNA - A[k, j] * CSA
					}
				}
			}
		}
		// find sigma2 for transformed A and test for convergence
		sigma2 = 0.;							// sigma2 = 0.
		for (i=0; i<N; i++) {					//for i in range(N) :
			EIGEN[i] = *(A + i * N + i);		// EIGEN[i] = A[i, i]
			sigma2 += EIGEN[i] * EIGEN[i];      // sigma2 += EIGEN[i] * *2
		}
		if (((1. - sigma1 / sigma2) < Eps3) && (abs(sigma1 - S) < Eps3)) { //:  // if (1. - sigma1 / sigma2) < Eps3 and abs(sigma1 - S)<Eps3 :  # convergence
			val = -1.0e3;						// val, ind = -1.0e3, -1
			ind = -1;
			for (i = 0; i < N; i++) { // for i in range(N) : # find largest eigenvalue
				if (EIGEN[i] > val) {			//if EIGEN[i] > val:
					ind = i;                    // ind = i
					val = EIGEN[i];             // val = EIGEN[i] 
				}
			}
			for (i=0; i<N; i++) { *(A+i) = T[i][ind]; }  //	return T[:, ind]
			return(EIGEN[ind]);				// return EIGEN[ind]
		}
		sigma1 = sigma2;
	}
	return(-1.);
}
*/