//

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <fstream>
#include <iostream>
#include <iomanip>

#include "ConFemMatC.h"
#include "dsyevj3.c"
#include "dsyevv3.c"
#include "dsyevc3.c"

using namespace std;
const double ZeroD = 1.e-9;
const double pi = 3.14159265359;
const double TOL = -1.e-15;
const double kapZero = 1.0e-9;
//const double ff = 1.e7;      // used as scaling factor for strains in eighC, be careful with this: has considerable impact in case of two same principal values in connection with rc 114, see also PrincipalStrains.mws
const double ff = 1.e9; //1.e8;      // used as scaling factor for strains in eighC, be careful with this: has considerable impact in case of two same principal values in connection with rc 114, see also PrincipalStrains.mws
// following for eigJacobiSym
const double Eps1 = 1.0e-10;
const double Eps2 = 1.0e-10;
const double Eps3 = 1.0e-7;
// Microplane
const double I21Points[21][3] = { 1., 0., 0.,\
								 0., 1., 0.,\
								 0., 0., 1.,\
								 0.707106781187, 0.707106781187, 0.,\
								 0.707106781187, -0.707106781187, 0.,\
								 0.707106781187, 0., 0.707106781187,\
								 0.707106781187, 0., -0.707106781187,\
								 0., 0.707106781187, 0.707106781187,\
								 0., 0.707106781187, -0.707106781187,\
								 0.387907304067, 0.387907304067, 0.836095596749,\
								 0.387907304067, 0.387907304067, -0.836095596749,\
								 0.387907304067, -0.387907304067, 0.836095596749,\
								 0.387907304067, -0.387907304067, -0.836095596749,\
								 0.387907304067, 0.836095596749, 0.387907304067,\
								 0.387907304067, 0.836095596749, -0.387907304067,\
								 0.387907304067, -0.836095596749, 0.387907304067,\
								 0.387907304067, -0.836095596749, -0.387907304067,\
								 0.836095596749, 0.387907304067, 0.387907304067,\
								 0.836095596749, 0.387907304067, -0.387907304067,\
								 0.836095596749, -0.387907304067, 0.387907304067,\
								 0.836095596749, -0.387907304067, -0.387907304067 };
const double I21Weights[21] = { 0.0265214244093, 0.0265214244093, 0.0265214244093,\
	0.0199301476312, 0.0199301476312, 0.0199301476312, 0.0199301476312, 0.0199301476312, 0.0199301476312,\
	0.0250712367487, 0.0250712367487, 0.0250712367487, 0.0250712367487, 0.0250712367487, 0.0250712367487, 0.0250712367487, 0.0250712367487, 0.0250712367487, 0.0250712367487, 0.0250712367487, 0.0250712367487 };

void PrinCLT_1(double vx, double vy, double vxy, double *pEps)
{
	double a, b, c, si;
    a = 0.5*(vx-vy);
    b = sqrt(a*a+vxy*vxy);
    if( b>ZeroD ) {
        c = 0.5*(vx+vy);
        *(pEps+0) = c + b;
        *(pEps+2) = c - b;
		if( vxy<0. ) { si = -1.; }
		else         { si = 1.; }
        *(pEps+1) = 0.5*si*acos(a/b);
	}
    else {
		if( vx>vy ) {*(pEps+0)= vx; *(pEps+1)= 0.;     *(pEps+2)= vy;}
		else        {*(pEps+0)= vy; *(pEps+1)= 0.5*pi; *(pEps+2)= vx;}
	}
}

double PrinCLT_2(double vx, double vy, double vxy)
{
	double a, b, c;
    a = 0.5*(vx-vy);
    b = sqrt(a*a+vxy*vxy);
    c = 0.5*(vx+vy);
	return c + b;
}

/* IsoDamage; return codes 100 - 110: no exception;  111-199: exception */

/* eigenvalues and eigenvectors of strain */
int eighC( double I1, double J2S, double *Eps, double *la, double *vv)  // replaced by eighC1
{
	double I2, I3, aa, bb, xx, phi, n0, n1, n2, na, deno, lambda, *vvp, p; //, ff;
	int i;
//	ff = 1.e3;          // for upscaling of small strains for eigenvector computation
	if (fabs(J2S)<ZeroD) {         // nearly three equal principal values
		p = I1/3.;
		 *la    = p;
		*(la+1) = p;
		*(la+2) = p;
		*vv     = 1.;
		*(vv+1) = 0.;
		*(vv+2) = 0.;
		*(vv+3) = 0.;
		*(vv+4) = 1.;
		*(vv+5) = 0.;
		*(vv+6) = 0.;
		*(vv+7) = 0.;
		*(vv+8) = 1.;
		return 101;
	}
	I2 = Eps[0]*Eps[1]+Eps[1]*Eps[2]+Eps[2]*Eps[0]-Eps[3]*Eps[3]-Eps[4]*Eps[4]-Eps[5]*Eps[5];
	I3 = Eps[0]*Eps[1]*Eps[2] - Eps[0]*Eps[3]*Eps[3] - Eps[1]*Eps[4]*Eps[4] - Eps[2]*Eps[5]*Eps[5] +2*Eps[3]*Eps[4]*Eps[5];
	bb = I1*I1-3*I2;
	if (bb<0.) { return 111; }
	aa = sqrt(bb);
	xx = (2.*I1*I1*I1-9*I1*I2+27*I3)/(2.*aa*aa*aa);
	if (xx>1.+ ZeroD) { return 112; }
	if (xx>1.) { xx=1.; }
	phi = acos(xx)/3.;
	 *la    = I1/3.+2./3.*aa*cos(phi);					// eigenvalues - principal stresses
	*(la+1) = I1/3.+2./3.*aa*cos(phi-2.*pi/3.);
	*(la+2) = I1/3.+2./3.*aa*cos(phi-4.*pi/3.);
	if ( (*la<(*(la+1)-ZeroD)) || (*(la+1)<(*(la+2)-ZeroD)) ) { return 113; } 
	n0 = 1.;
	for (i=0;i<6;i++) { Eps[i] = ff*Eps[i]; }
	for (i=0;i<3;i++) {								// eigenvectors
		vvp = vv+3*i;								// pointer for eigenvector in matrix of eigenvectors
		lambda = ff*(*(la+i));
		deno = Eps[5]*Eps[2]-Eps[5]*lambda-Eps[3]*Eps[4];
		if ( fabs(deno)<ZeroD ) { 
			*(vvp)   = 0.;
			*(vvp+1) = 0.;
			*(vvp+2) = 0.;
			*(vvp+i) = 1.;   // it has not been checked whether this leads to consistent results under all conditions
			continue;
		}
		n1 = -(Eps[0]*Eps[2]-Eps[0]*lambda-lambda*Eps[2]+lambda*lambda-Eps[4]*Eps[4])/deno;
		n2 = (-Eps[5]*Eps[4]+Eps[3]*Eps[0]-Eps[3]*lambda)                            /deno;
		na = 1./sqrt(n0*n0+n1*n1+n2*n2);
		*(vvp)   = n0*na;
		*(vvp+1) = n1*na;
		*(vvp+2) = n2*na;
	}	
	return 100;
}
int eighC2( double I1, double J2S, double *Eps, double *la, double *vv)
{
	double I2, I3, aa, bb, xx, phi, n0, n1, n2, na, lambda, *vvp, p, deno0, deno1, deno2; //, ff;
	int i;
	if (fabs(J2S)<ZeroD) {         // nearly three equal principal values
		p = I1/3.;
		*la     = p;
		*(la+1) = p;
		*(la+2) = p;
		*vv     = 1.;
		*(vv+1) = 0.;
		*(vv+2) = 0.;
		*(vv+3) = 0.;
		*(vv+4) = 1.;
		*(vv+5) = 0.;
		*(vv+6) = 0.;
		*(vv+7) = 0.;
		*(vv+8) = 1.;
		return 101;
	}
	I2 = Eps[0]*Eps[1]+Eps[1]*Eps[2]+Eps[2]*Eps[0]-Eps[3]*Eps[3]-Eps[4]*Eps[4]-Eps[5]*Eps[5];
	I3 = Eps[0]*Eps[1]*Eps[2] - Eps[0]*Eps[3]*Eps[3] - Eps[1]*Eps[4]*Eps[4] - Eps[2]*Eps[5]*Eps[5] +2*Eps[3]*Eps[4]*Eps[5];
	bb = I1*I1-3.*I2;
	if (bb<0.) { return 111; }
	aa = sqrt(bb);
	xx = (2.*I1*I1*I1-9.*I1*I2+27*I3)/(2.*aa*aa*aa);
	if (xx>1.+ ZeroD) { return 112; }
	if (xx>1.) { xx=1.; }
	phi = acos(xx)/3.;
	 *la    = I1/3.+2./3.*aa*cos(phi);					// eigenvalues - principal stresses
	*(la+1) = I1/3.+2./3.*aa*cos(phi-2.*pi/3.);
	*(la+2) = I1/3.+2./3.*aa*cos(phi-4.*pi/3.);
//	if ( (*la<(*(la+1)-ZeroD)) || (*(la+1)<(*(la+2)-ZeroD)) ) { return 113; }  / should always be fulfilled 
	for (i=0;i<6;i++) { Eps[i] = ff*Eps[i]; }
	for (i=0;i<2;i++) {								// first two eigenvectors
		vvp = vv+3*i;								// pointer for eigenvector in matrix of eigenvectors
		lambda = ff*(*(la+i));
		deno0 = -Eps[3]*Eps[3]+Eps[1]*Eps[2]-Eps[1]*lambda-lambda*Eps[2]+lambda*lambda;
		deno1 = -Eps[4]*Eps[4]+Eps[2]*Eps[0]-Eps[2]*lambda-lambda*Eps[0]+lambda*lambda;
		deno2 = -Eps[5]*Eps[5]+Eps[1]*Eps[0]-Eps[0]*lambda-lambda*Eps[1]+lambda*lambda;
		if ( fabs(deno0)>ZeroD ) { 
			n0 =   1.;
			n1 = ( Eps[4]*Eps[3]-Eps[2]*Eps[5]+lambda*Eps[5])/deno0;
			n2 =-( Eps[1]*Eps[4]-Eps[3]*Eps[5]-lambda*Eps[4])/deno0;
		}
		else if ( fabs(deno1)>ZeroD ) {
			n0 = ( Eps[3]*Eps[4]+lambda*Eps[5]-Eps[2]*Eps[5])/deno1;
			n1 =   1.;
			n2 =-( Eps[0]*Eps[3]-lambda*Eps[3]-Eps[5]*Eps[4])/deno1;
		}
		else if ( fabs(deno2)>ZeroD ) {
			n0 =-(-Eps[5]*Eps[3]+Eps[1]*Eps[4]-lambda*Eps[4])/deno2;
			n1 =-( Eps[0]*Eps[3]-lambda*Eps[3]-Eps[5]*Eps[4])/deno2;
			n2 = 1.;
		}
		else {
/*
			*(vv+9) = *(la);
			*(vv+10) = *(la+1);
			*(vv+11) = *(la+2);
			*(vv+12) = lambda;
			*(vv+13) = deno0;
			*(vv+14) = deno1;
			*(vv+15) = deno2;
			*(vv+16) = deno_[0]; //Eps[1];
			*(vv+17) = deno_[1]; //Eps[2];
*/
			return 114; 
		}
		na = 1./sqrt(n0*n0+n1*n1+n2*n2);
		*(vvp)   = n0*na;
		*(vvp+1) = n1*na;
		*(vvp+2) = n2*na;
	}	
	*(vv+6) =  (*(vv+1))*(*(vv+5)) - (*(vv+2))*(*(vv+4));           // third eigenvector from dyadic product 
	*(vv+7) =  (*(vv+2))*(*(vv+3)) - (*(vv+0))*(*(vv+5));
	*(vv+8) =  (*(vv+0))*(*(vv+4)) - (*(vv+1))*(*(vv+3));
/*
*(vv+9) = *(la);
*(vv+10) = *(la+1);
*(vv+11) = *(la+2);
*(vv+12) = lambda;
*(vv+13) = deno0;
*(vv+14) = deno1;
*(vv+15) = deno2;
*(vv+16) = bb;
*(vv+17) = xx;
*/
	return 100;
}

int eigJacobiSymWrapper(double *Eps, int dim_Eps, double *la, int dim_la, double *vv, int dim_vv)
{
	int rc;
	rc = eigJacobiSym( Eps, la, vv);
	return rc;
}
int eigJacobiSym(double *Eps, double *la, double *vv)
{
	//	see starting lines for Eps1, Eps2, Eps3 = 1.0e-10, 1.0e-10, 1.0e-7
	int i,j,k, iter, itmax, N, ind;
	double sigma1, sigma2, OffDsq, S, T[3][3], AIK[3], EIGEN[3], A[3][3], Q, Q_, siQ, P, siP, alpha, CSA, SNA, HOLDKI, val;
	itmax = 10;
	N = 3;
	sigma1 = 0.;
	OffDsq = 0.;
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			T[i][j] = 0.;
		}
		T[i][i] = 1.;
		AIK[i]  = 0.;
		EIGEN[i]= 0.;
	}
	A[0][0] = *(Eps + 0);
	A[1][1] = *(Eps + 1);
	A[2][2] = *(Eps + 2);
	A[1][2] = *(Eps + 3);											// engineering strains of voigt notation not considered here
	A[0][2] = *(Eps + 4);
	A[0][1] = *(Eps + 5);
	for (i=0; i<N; i++) {
		sigma1 = sigma1 + A[i][i]*A[i][i];
		for (j =i+1; j<N; j++) {
			OffDsq = OffDsq + A[i][j]*A[i][j];
		}
	}
	S = 2.*OffDsq + sigma1;											// square sum of elements of A
	if (sqrt(S) < Eps1) { 											// numerically only zeros in A
		*(la) = 0.;
		*(vv + 0) = 1.;
		*(vv + 1) = 0.;
		*(vv + 2) = 0.;
		return 100;
	}
	// iteration loop
	for (iter = 0; iter < itmax; iter++) {
		for (i=0; i<N-1; i++) {
			for (j=i +1; j<N; j++) {
				// commented out as error uhc 190612. If, e.g. A[0][1] is zero, and A[0][2] not the loop is left to early
				// if (abs(A[i][j]) <= Eps2) break;					//  actual out of diagonal element numerically zero
				//
				Q_ = A[i][i] - A[j][j];
				Q = abs(Q_);										// diagonal i, j difference
				//	 compute sine and cosine of rotation angle
				if (Q > Eps1) {
					if (Q_> 0.) siQ = 1.;
					else        siQ = -1.;
					P = 2.*A[i][j] * siQ;
					if      (P > 0.)    siP =  1.;					// P!=0 as A[i][j] !=0
					else                siP = -1.;
					alpha = P / Q;
					CSA = sqrt(0.5*(1. + 1. / sqrt(1. + alpha * alpha))); // this is obviously larger than sqrt(1 / 2)
					SNA = siP / (2.*CSA*sqrt(1. + 1. /(alpha * alpha)));
				}
				else {
					CSA = 1. / sqrt(2.);
					SNA = CSA;
				}
				for (k = 0; k < N; k++) {
					HOLDKI = T[k][i];
					T[k][i] = HOLDKI * CSA + T[k][j] * SNA;
					T[k][j] = HOLDKI * SNA - T[k][j] * CSA;
				}
				// rows i and j--> transpose(T) * A
				for (k = i; k < N; k++) {
					if (k <= j) {
						AIK[k] = A[i][k];
						A[i][k] = AIK[k] * CSA + A[k][j] * SNA;
						if (k == j) { A[j][k] = AIK[k] * SNA - A[j][k] * CSA; }
					}
					else {
						HOLDKI = A[i][k];
						A[i][k] = HOLDKI * CSA + A[j][k] * SNA;
						A[j][k] = HOLDKI * SNA - A[j][k] * CSA;
					}
				}
				AIK[j] = SNA * AIK[i] - CSA * AIK[j];
				// columns i and j--> A * T
				for (k = 0; k < (j + 1); k++) {   //in range(j + 1) :
					if (k > i) { A[k][j] = SNA * AIK[k] - CSA * A[k][j]; }
					else {
						HOLDKI = A[k][i];
						A[k][i] = HOLDKI * CSA + A[k][j] * SNA;
						A[k][j] = HOLDKI * SNA - A[k][j] * CSA;
					}
				}
//              A[i, j] = 0.;
			}
		}

		// find sigma2 for transformed A  
		sigma2 = 0.;
		for (i = 0; i < N; i++) {
			EIGEN[i] = A[i][i];
			sigma2 = sigma2 + EIGEN[i] * EIGEN[i];
		}
		// and test for convergence
		if ( ((1. - sigma1 / sigma2) < Eps3) & (abs(sigma1 - S) < Eps3) ) {
			val = -1.0e6;									// might destroy slution in extreme situations
			ind = -1;
			for (i = 0; i < N; i++)							// find largest eigenvalue
				if (EIGEN[i] > val) {
					ind = i;
					val = EIGEN[i];
				}
			*(la + 0) = val;
			*(vv + 0) = T[0][ind];
			*(vv + 1) = T[1][ind];
			*(vv + 2) = T[2][ind];
//			*(la + 1) = CSA;							// for control purposes
//			*(la + 2) = SNA;
			return 100;
		}
		sigma1 = sigma2;
	}
	//	 no convergence
	return 114;
}

int EquivStrain1C( int CalcType, double I1, double J2, double J2S, double *Eps, double *EpsD, double *kap_, double *nd, double *Hd, double cc0, double cc1, double cc2, double cc3,
						double *laMax_, double *eVec, double *data) // data to bring data out for control
{
	double la[3], xxx, eins[6], nd_[6], f1, f2, f3, laMax, kap, vv[9];
	int rc;
	double A[3][3], Q[3][3], w[3];
	rc = eigJacobiSym( Eps, la, vv);
	//	rc = dsyevj3(A, Q, w);
	A[0][0] = *(Eps + 0);
	A[1][1] = *(Eps + 1);
	A[2][2] = *(Eps + 2);
	A[1][2] = *(Eps + 3);											// engineering strains of voigt notation not considered here
	A[0][2] = *(Eps + 4);
	A[0][1] = *(Eps + 5);
//	rc = dsyevv3(A, Q, w, la, vv);
//	rc = dsyevv3(A, la, vv);
	if (rc>110) {
/*
		*(data+0) = *(vv+9);
		*(data+1) = *(vv+10);
		*(data+2) = *(vv+11);
		*(data+3) = *(vv+12);
		*(data+4) = *(vv+13);
		*(data+5) = *(vv+14);
		*(data+6) = *(vv+15);
		*(data+7) = *(vv+16);
		*(data+8) = *(vv+17);
*/
		return rc; 
	}   
	laMax = la[0]; //  ila = 0  # largest principal value   \  if la[1]>la[0]: ila=1 \  if la[2]>la[1] and la[2]>la[0]: ila=2 \  laMax = la[ila]
    xxx = 0.25*pow((cc1*J2S+cc2*laMax+cc3*I1),2.0) + cc0*J2;
	kap = 0.5*(cc1*J2S+cc2*laMax+cc3*I1)+(sqrt(xxx)); 
	*kap_ = kap;
	*laMax_ = laMax;
	*(eVec + 0) = *(vv + 0);
	*(eVec + 1) = *(vv + 1);
	*(eVec + 2) = *(vv + 2);
	if (CalcType == 2) {
		if (J2S > ZeroD) {
			eins[0] = 1.;
			eins[1] = 1.;
			eins[2] = 1.;
			eins[3] = 0.;
			eins[4] = 0.;
			eins[5] = 0.;
			nd_[0] = (*(vv + 0)) * (*(vv + 0));
			nd_[1] = (*(vv + 1)) * (*(vv + 1));
			nd_[2] = (*(vv + 2)) * (*(vv + 2));
			nd_[3] = (*(vv + 1)) * (*(vv + 2));
			nd_[4] = (*(vv + 0)) * (*(vv + 2));
			nd_[5] = (*(vv + 0)) * (*(vv + 1));
			f1 = cc0 + cc1 * (kap) / (2.0 * J2S);
			f2 = cc2 * (kap);
			f3 = cc3 * (kap);
			*(nd + 0) = f1 * (*(EpsD + 0)) + f2 * nd_[0] + f3 * eins[0]; // nd  = (self.cc0+self.cc1*kap_/(2.0*J2S))*EpsD + self.cc2*kap_*nd_ +self.cc3*kap_*eins # gradient of damage function
			*(nd + 1) = f1 * (*(EpsD + 1)) + f2 * nd_[1] + f3 * eins[1];
			*(nd + 2) = f1 * (*(EpsD + 2)) + f2 * nd_[2] + f3 * eins[2];
			*(nd + 3) = f1 * (*(EpsD + 3)) + f2 * nd_[3] + f3 * eins[3];
			*(nd + 4) = f1 * (*(EpsD + 4)) + f2 * nd_[4] + f3 * eins[4];
			*(nd + 5) = f1 * (*(EpsD + 5)) + f2 * nd_[5] + f3 * eins[5];
			*Hd = -(cc1 * J2S + cc2 * laMax + cc3 * I1 - 2 * (kap)); // Hd = -( self.cc1*J2S + self.cc2*laMax + self.cc3*I1 -2*kap_ ) # H_d: dF/dkappa local 
		}
		else {
			*nd = 0.;
			*(nd + 1) = 0.;
			*(nd + 2) = 0.;
			*(nd + 3) = 0.;
			*(nd + 4) = 0.;
			*(nd + 5) = 0.;
			*Hd = 1.;   //
		}
	}

//	*(data + 0) = la[0];
//	*(data + 1) = vv[0];
//	*(data + 2) = vv[1];
//	*(data + 3) = vv[2];
//	*(data + 4) = la[0];
//	*(data + 5) = vv[0];
//	*(data + 6) = vv[1];
//	*(data + 7) = vv[2];
//	*(data + 4) = w[0];
//	*(data + 5) = w[1];
//	*(data + 6) = w[2];
//	*(data + 7) = rc;

	return rc;  
}
double ViscExten3DC1( double Dt, double eta, double *Dps, double *ElemStateVar, double *ElemStateVarN, int sI, double *Veps)  // def ViscExten3D(self, Dt, eta, Dps, Elem, ipI, sI): 
{
	double zz, NormDps, normVepsOld, dotDV, *ElSVar, *ElSVarN, VepsOld[6];
	int i;
	ElSVar = ElemStateVar+sI;  // VepsOld = array([Elem.StateVar[ipI,sI],Elem.StateVar[ipI,sI+1],Elem.StateVar[ipI,sI+2],Elem.StateVar[ipI,sI+3],Elem.StateVar[ipI,sI+4],Elem.StateVar[ipI,sI+5]]) # strain rate of previous time step 
	ElSVarN= ElemStateVarN+sI; // 
	for (i=0;i<6;i++) { VepsOld[i] = *(ElSVar+i); }
	normVepsOld = VepsOld[0]*VepsOld[0] + VepsOld[1]*VepsOld[1] + VepsOld[2]*VepsOld[2] + VepsOld[3]*VepsOld[3] + VepsOld[4]*VepsOld[4] + VepsOld[5]*VepsOld[5]; 
	dotDV =       VepsOld[0]*(*(Dps+0)) + VepsOld[1]*(*(Dps+1)) + VepsOld[2]*(*(Dps+2)) + VepsOld[3]*(*(Dps+3)) + VepsOld[4]*(*(Dps+4)) + VepsOld[5]*(*(Dps+5));
	NormDps =     (*(Dps+0))*(*(Dps+0)) + (*(Dps+1))*(*(Dps+1)) + (*(Dps+2))*(*(Dps+2)) + (*(Dps+3))*(*(Dps+3)) + (*(Dps+4))*(*(Dps+4)) + (*(Dps+5))*(*(Dps+5));
	normVepsOld = sqrt(normVepsOld);
	NormDps = sqrt(NormDps);
	if ((normVepsOld<ZeroD) || (dotDV<0.)) {                //   if norm(VepsOld)<ZeroD or dot(Dps,VepsOld)<0.: 
		if ( Dt>ZeroD)                 { for (i=0;i<6;i++) { VepsOld[i] = (*(Dps+i))/Dt; } }            //             if Dt>ZeroD: VepsOld = Dps/Dt
		else                           { for (i=0;i<6;i++) { VepsOld[i] =    0.;         } }            //             else:        VepsOld = zeros((6), dtype=float)
	}
	if ((Dt>ZeroD) && (NormDps>ZeroD)) { for (i=0;i<6;i++) { *(Veps+i) = 2.*(*(Dps+i))/Dt - VepsOld[i]; }	} //   if Dt>ZeroD and norm(Dps)>ZeroD: Veps = 2.*Dps/Dt - VepsOld             # actual strain rate
	else {                               for (i=0;i<6;i++) { *(Veps+i) =                    VepsOld[i]; }	} //         else:                              Veps = VepsOld
	if (Dt>ZeroD)                      { zz=2.*eta/Dt; }  //  if Dt>ZeroD: zz = 2.*self.eta/Dt 
	else                               { zz=0.;}   // else:        zz = 0. 
	for (i=0;i<6;i++) { *(ElSVarN+i) = *(Veps+i); } //         Elem.StateVarN[ipI,sI]   = Veps[0] / Elem.StateVarN[ipI,sI+1] = Veps[1] / Elem.StateVarN[ipI,sI+2] = Veps[2] / Elem.StateVarN[ipI,sI+3] = Veps[3] / Elem.StateVarN[ipI,sI+4] = Veps[4] / Elem.StateVarN[ipI,sI+5] = Veps[5]
	return zz;  // return zz, Veps
}


int IsoDamC1( int CalcType, int ElemDim, bool ElemPlSt, bool PrinStrains, double ElemLch, double *ElemStateVar, int dim_St, double *ElemStateVarN, int dim_StN,
					double *Eps, int dim_E, double *sig, int dim_si, double *MatM2, int dim_Ma21, int dim_Ma22, int LiTy,
					double cc0, double cc1, double cc2, double cc3, int RType, double *EpsR, int dim_EpsR, double kapStrength, double ElemCrBws, double gam2, double kapUlt,
					double edt, double ed, double gd, double nu, double Emod, double *Dps, int dim_D, double eta, double RegPar, int ElemScaleType,
					double *sigR, int dim_sigR, double *CR, int dim_CR, double Dt, double sVsTol, double *DataOut, int dim_DataOut) 
{
//	double I1, pp, J2, J2S, J3, xi, Eps_[6], Dps__[6], EpsD[6], kap_, nd[6], Hd, kap, dkk, beta, hdI, alpha, xxx,xx2, sig0[6], sigV[6], CD[6][6], zz, Veps[6], ccc, fact, data[10], D, kapOld, zzz;
	double I1, pp, J2, J2S, J3, xi, Eps_[6], Dps__[6], EpsD[6], kap_, nd[6], Hd, kap, dkk, beta, hdI, alpha, xxx, xx2, sig0[6], sigV[6], CD[6][6], zz, Veps[6], ccc, fact, data[1], D, kapOld, zzz;
	double sigD[6]; // , tt[3];
	double xxnu, x1nu, x2nu, ff0, ff1, ff5, svs, svs_, LaMax_, eVec[3];
	double laMax;
	double* MatM;
	int rc, i, j;
	MatM = MatM2;
	D      = *(ElemStateVar+0);
	kapOld = *(ElemStateVar+1);
	if (ElemDim==1) { 
		Eps_[0]  =        *Eps; 
		Eps_[1]  = -nu * (*Eps); 
		Eps_[2]  = -nu * (*Eps);
		Eps_[3]  =        0.;
		Eps_[4]  =        0.;
		Eps_[5]  =        0.;
		if (eta > 0.) {
			Dps__[0] = *Dps;
			Dps__[1] = -nu * (*Dps);
			Dps__[2] = -nu * (*Dps);
			Dps__[3] = 0.;
			Dps__[4] = 0.;
			Dps__[5] = 0.;
		}
	}
	else if (ElemDim==2) {  // elif Elem.dim==2:
		if (ElemPlSt) { 
			alpha = nu/(1-nu);
			Eps_[0]  = *Eps; 
			Eps_[1]  = *(Eps+1); 
			Eps_[2]  = -alpha*((*Eps)+(*(Eps+1)));
			Eps_[3]  = 0.;
			Eps_[4]  = 0.;
			Eps_[5]  = 0.5*(*(Eps+2));        // tensorial strain in voigt notation
			if (eta > 0.) {
				Dps__[0] = *(Dps);
				Dps__[1] = *(Dps + 1);
				Dps__[2] = -alpha * ((*Dps) + (*(Dps + 1)));
				Dps__[3] = 0.;
				Dps__[4] = 0.;
				Dps__[5] = *(Dps + 2);
			}
		}
		else { 
			Eps_[0]  = *(Eps); 
			Eps_[1]  = *(Eps+1); 
			Eps_[2]  = 0.;
			Eps_[3]  = 0.;
			Eps_[4]  = 0.;                        // plane strain
			Eps_[5]  = 0.5*(*(Eps+2));            // tensorial strain in voigt notation
			if (eta > 0.) {
				Dps__[0] = *(Dps);
				Dps__[1] = *(Dps + 1);
				Dps__[2] = 0.;
				Dps__[3] = 0.;
				Dps__[4] = 0.;
				Dps__[5] = *(Dps + 2);                  // engineering strain voigt notation
			}
		}
	}
	else if (ElemDim==3) { 
		Eps_[0]  = *(Eps); 
		Eps_[1]  = *(Eps+1); 
		Eps_[2]  = *(Eps+2);
		Eps_[3]  = 0.5*(*(Eps+3));
		Eps_[4]  = 0.5*(*(Eps+4));
		Eps_[5]  = 0.5*(*(Eps+5));
		if (eta > 0.) {
			Dps__[0] = *(Dps);
			Dps__[1] = *(Dps + 1);
			Dps__[2] = *(Dps + 2);
			Dps__[3] = *(Dps + 3);
			Dps__[4] = *(Dps + 4);
			Dps__[5] = *(Dps + 5);
		}
	}
	else if (ElemDim == 4) {               // axisymmetric CP4
		Eps_[0] = *(Eps);
		Eps_[1] = *(Eps + 1);
		Eps_[2] = *(Eps + 2);
		Eps_[3] = 0.;
		Eps_[4] = 0.;
		Eps_[5] = 0.5 * (*(Eps + 3));
		if (eta > 0.) {
			Dps__[0] = *(Dps);
			Dps__[1] = *(Dps + 1);
			Dps__[2] = *(Dps + 2);
			Dps__[3] = 0.;
			Dps__[4] = 0.;
			Dps__[5] = *(Dps + 3);
		}
	}
	else if (ElemDim==21) { // continuum based shell  elif Elem.dim==21:                          # continuum based shell
		xxx = -nu/(1-nu)*( (*(Eps)) + (*(Eps+1)) );  
		Eps_[0]  = *(Eps);  
		Eps_[1]  = *(Eps+1); 
		Eps_[2]  =   xxx;
		Eps_[3]  = 0.5*(*(Eps+3));
		Eps_[4]  = 0.5*(*(Eps+4));
		Eps_[5]  = 0.5*(*(Eps+5));
		if (eta > 0.) {
			xxx = -nu / (1 - nu) * ((*(Dps)) + (*(Dps + 1)));  // xxx = -self.nu/(1-self.nu)*(Dps_[0]+Dps_[1])
			Dps__[0] = *(Dps);                            // Dps__ = array([Dps_[0],Dps_[1],xxx,Dps_[3],Dps_[4],Dps_[5]])                            # triaxial strain Voigt notation
			Dps__[1] = *(Dps + 1);
			Dps__[2] = xxx;
			Dps__[3] = *(Dps + 3);
			Dps__[4] = *(Dps + 4);
			Dps__[5] = *(Dps + 5);
		}
	}
	else { return 121;}   // else: raise NameError("ConFemMaterials::IsoDamage.Sig: element type not implemented for this material")

	I1 = *Eps_+*(Eps_+1)+*(Eps_+2);   // 1st invariant strain tensor    I1=Eps[0,0]+Eps[1,1]+Eps[2,2]               # 1st invariant strain tensor
	pp = I1/3.;                     // volumetric stress             pp=I1/3.                                    # volumetric stress
	EpsD[0] =  *Eps_   - pp;      
	EpsD[1] = *(Eps_+1)- pp;
	EpsD[2] = *(Eps_+2)- pp;
	EpsD[3] = *(Eps_+3);
	EpsD[4] = *(Eps_+4);
	EpsD[5] = *(Eps_+5);
    J2 =0.5*(EpsD[0]*EpsD[0]+EpsD[1]*EpsD[1]+EpsD[2]*EpsD[2])+EpsD[3]*EpsD[3]+EpsD[4]*EpsD[4]+EpsD[5]*EpsD[5]; //# 2nd invariant of deviator
    J2S=sqrt(J2);
    J3=EpsD[0]*EpsD[1]*EpsD[2]-EpsD[0]*EpsD[4]*EpsD[4]-EpsD[3]*EpsD[3]*EpsD[2]+2.*EpsD[3]*EpsD[5]*EpsD[4] -EpsD[5]*EpsD[5]*EpsD[1]; // J3=EpsD[0]*EpsD[1]*EpsD[2]-EpsD[0]*EpsD[4]**2-EpsD[3]**2*EpsD[2]+2 *EpsD[3]*EpsD[5]*EpsD[4] -EpsD[5]**2*EpsD[1]# 3rd invariant of strain deviator
//  xi seems not to be used																																	///	if (J2S>ZeroD) {
//		xi = 0.5*(3.*sqrt(3.)/2.*J3/sqrt(J2*J2*J2) + 1.); // if J2S>ZeroD: xi = 0.5 * (StrNum0*J3/sqrt(J2**3) + 1.)    # tension indicator
//	}
//	else {
//		xi = -1.;							// else:         xi = -1.
//	}
	if (LiTy==1) {     //  if self.LiTy==1: kap_, nd, Hd = self.EquivStrain1( I1, J2, J2S, Eps, EpsD, ipI) / else: raise NameError("ConFemMaterials::IsoDamage.Sig: unknown type of limit function") 
		rc = EquivStrain1C( CalcType, I1, J2, J2S, Eps_, EpsD, &kap_, nd, &Hd, cc0, cc1, cc2, cc3, &LaMax_, eVec, data);          // note! Eps- might have been modified during principal value computation
		// solution for largest principal strain direction (--> nd ) not unique for, e.g., unaxial compression. Might lead to differences regarding python version.

//		*(DataOut + 0) = data[0];
//		*(DataOut + 1) = data[1];
//		*(DataOut + 2) = data[2];
//		*(DataOut + 3) = data[3];
//		*(DataOut + 4) = data[4];
//		*(DataOut + 5) = data[5];
//		*(DataOut + 6) = data[6];
//		*(DataOut + 7) = data[7];

		if (rc>110) { return rc; }
	}
	else { return 120; }
	/* Regularization */
	if (RType==1) {							//  if self.RType==1: /  if Elem.dim==  1: kap = EpsR[1] / elif Elem.dim==2: kap = EpsR[3] / elif Elem.dim==3: kap = EpsR[6] /  dkk = 1.
		if      (ElemDim==1 ) { kap = *(EpsR+1); }
		else if (ElemDim==2 ) { kap = *(EpsR+3); }
		else if (ElemDim==3 ) { kap = *(EpsR+6); }
		dkk = 1.;
	}
	else if (RType==2) {					//  elif self.RType==2:              
		if (kap_>kapStrength) {				//  if kap_> self.kapStrength:
			beta = ElemCrBws;				//  beta = Elem.CrBwS                   # element specific scaling factor
            if (ElemScaleType==1) {
				kap = beta*(kap_-kapStrength) + (1.-beta)/gam2 * (1-exp(-gam2*(kap_-kapStrength))) + kapStrength; // # scaled damage strain kap = beta*(kap_-self.kapStrength) + (1.-beta)/self.gam2 * (1-exp(-self.gam2*(kap_-self.kapStrength))) + self.kapStrength
				dkk = beta                    + (1.-beta)      *    exp(-gam2*(kap_-kapStrength)); // # scaling factor for tangential material stiffness dkk = beta                         + (1.-beta)           *    exp(-self.gam2*(kap_-self.kapStrength))
			}
			else {
				kap = (1.-beta)*kapStrength * log((kap_ - beta*kapStrength) / ((1.-beta)*kapStrength)) + kapStrength;      //(1 - beta)*self.kapStrength * log((kap_ - beta * self.kapStrength) / ((1 - beta)*self.kapStrength)) + self.kapStrength
				dkk = (1.-beta)*kapStrength /     (kap_ - beta*kapStrength);              //(1 - beta)*self.kapStrength / (kap_ - beta * self.kapStrength)
			}
		}
		else {								//   else:
			kap = kap_;						// kap = kap_
			dkk = 1.;						// dkk = 1.
		}
	}
	else {									// else:
		kap = kap_;							// kap = kap_ 
		dkk = 1.;							//  dkk = 1.  
	}
	xxx = Emod/(1+nu)/(1-2.*nu);
	xx2 = 2. * xxx;
	sig0[0] = xxx*((1 - nu)*Eps_[0] + nu*Eps_[1] + nu*Eps_[2]);    // sig0  = dot(self.C3_,Eps__)
	sig0[1] = xxx*(nu*Eps_[0] + (1 - nu)*Eps_[1] + nu*Eps_[2]);
	sig0[2] = xxx*(nu*Eps_[0] + nu*Eps_[1] + (1 - nu)*Eps_[2]);            
	sig0[3] = xx2*((0.5 - nu)*Eps_[3]);					// Eps_ is tensor component, must swirch to eng notation
	sig0[4] = xx2*((0.5 - nu)*Eps_[4]);
	sig0[5] = xx2*((0.5 - nu)*Eps_[5]);
	/*  */
	if (kap>kapUlt) { kap=kapUlt; } //  if kap>self.kapUlt: kap=self.kapUlt         # Damage should not be zero to avoid numerical singularity. This constrains D to dDestr 
//	if ((kap>edt) && (kap>kapOld+kapZero) && (J2S>=ZeroD)) {  // if kap>self.edt and kap>=(kapOld+ZeroD) and J2S>ZeroD: # case loading with nonzero strain deviator
	if ((kap>edt) && (kap>=kapOld) && (J2S >= ZeroD)) {  // if kap>self.edt and kap>=(kapOld+ZeroD) and J2S>ZeroD: # case loading with nonzero strain deviator
		D = 1.0 - exp(-pow(( kap-edt)/ed ,gd)); // D = 1.0 - exp(-pow(( kap-self.edt)/self.ed ,self.gd))  # scalar damage
		if (CalcType == 2) {
			hdI = pow((kap - edt) / ed, gd) * gd / (kap - edt) * exp(-pow((kap - edt) / ed, gd)) * dkk; 
			if (RType == 1) {  
				for (i = 0; i < 6; i++) {
					for (j = 0; j < 6; j++) {
						CD[i][j] = 0.;
			}	}	}
			else {								//    else:             CD = hdI/Hd * outer(sig0,nd)
				zzz = hdI / Hd;
				for (i = 0; i < 6; i++) {
					for (j = 0; j < 6; j++) {
						CD[i][j] = zzz * sig0[i] * nd[j];   // !!!! not symmetric
			}	}	}
		}
	}
	else {									//     else:
		if (CalcType == 2) {
			for (i = 0; i < 6; i++) {
				for (j = 0; j < 6; j++) {
					CD[i][j] = 0.;				//      CD = zeros((6,6),dtype=float)           # case of unloading or zero strain deviator
				}
				nd[i] = 0.;						// nd = zeros((6))
			}
			Hd = 1.;							//                 Hd = 1
			hdI = 0.;							//          hdI = 0
		}
	}
	*(ElemStateVarN+0) = D;					// store damage of actual iteration  // Elem.StateVarN[ipI,0] = D                   # store damage of actual iteration
	*(ElemStateVarN+1) = kap;				// store equivalent damage strain of actual iteration // Elem.StateVarN[ipI,1] = kap                 # 

	for (i = 0; i < 6; i++) { sigD[i] = (1 - D) * sig0[i]; }
	if (RType == 3) {
		// Rankine data
		rc = Rankine(PrinStrains, ElemDim, Eps_, sigD, eVec, &laMax);  // laMax may also be stresses, check whether this confused with LaMax_
		if (rc > 110) { return rc; }
	}
	zz = 0.;
	if ( eta>0.) {							// if self.eta>0.:                             # viscous regularization
		zz=ViscExten3DC1( Dt, eta, Dps__, ElemStateVar, ElemStateVarN, 2, Veps);  // def ViscExten3D(self, Dt, eta, Dps, Elem, ipI, sI):
		for (i=0;i<6;i++) { sigV[i] = eta*Veps[i]; }
		svs = 0.;								// svs: to monitor viscous stress related to isodam stress
		for (i = 0; i < 6; i++) {
			if (abs(*(sigD + i)) > sVsTol) { svs_ = sigV[i] / (*(sigD + i)); }
			else                           { svs_ = 0.; }
			if (abs(svs_) > abs(svs)) { svs = svs_; }
				*(DataOut + 0) = svs;								// output used in ConFem
		}
		for (i=0;i<6;i++) { sigD[i] = sigD[i] + sigV[i]; }
	}

	if (RType == 3) {
		// move Rankine data to state variables with viscous contributions for further processing regarding SDA
		rc = RankineUpd(sigD, eVec, ElemStateVarN, &laMax, 3);
	}
	/* crack energy increment  */
	if ((kap>kapStrength) && (kap>kapOld)) { ElemStateVarN[8] = *(sig+0)*Dps__[0]+*(sig+1)*Dps__[1]+*(sig+2)*Dps__[2]+*(sig+3)*Dps__[3]+*(sig+4)*Dps__[4]+*(sig+5)*Dps__[5];} // if kap>self.kapStrength and kap>kapOld: Elem.StateVarN[ipI,8] = dot(sig,Dps__) # Crack energy increment for step
	else                                   { ElemStateVarN[8] = 0.; }   // else: Elem.StateVarN[ipI,8] = 0.
	/* Specify general case for particular element dimensions */
	xxx = (1.-D)*xxx ; //Emod/(1+nu)/(1-2.*nu);   
	xxnu = xxx*(    nu);
	x1nu = xxx*(1. -nu);
	x2nu = xxx*(0.5-nu);
	if (ElemDim==1) {                       // if Elem.dim==1:
		if (RType==1) {                     // if self.RType==1:
			ccc = 0.5*RegPar*RegPar;		//              ccc = 0.5*self.RegPar**2  /               fact = 0.5*self.Emod / nd[0]=nd[0]-self.nu*nd[1]-self.nu*nd[2]
            fact = 0.5*Emod;
            nd[0]=nd[0]-nu*nd[1]-nu*nd[2];
			*(sig) = sigD[0]; // (1 - D)* sig0[0] + eta * Veps[0];                 // return  [sig[0],ccc*Eps_[1]*fact] 
			*(sig+1)  = ccc*Eps[1]*fact;
			*(MatM+0) = xxx*(1.-nu)+zz; 	*(MatM+1) = 0.;              //   return array([ [ CC[0,0] , 0 ] , [ 0 , ccc*fact ] ])
			*(MatM+2) = 0.;      			*(MatM+3) = ccc*fact;
            *(sigR+0) = 0.;                 //  return  [0,(EpsR[1]-kap_)*fact]
			*(sigR+1) = (EpsR[1]-kap_)*fact;
			*(CR+0)   = 0.;                 *(CR+1) = -Emod*Eps[0]*hdI;      //  CR = array([[0,-self.Emod*Eps_[0]*hdI],[-fact*nd[0]/Hd,fact]])
			*(CR+2)   = -fact*nd[0]/Hd;     *(CR+3) =  fact;               //  return [sig[0],ccc*Eps_[1]*fact], array([[CC[0,0],0],[0,ccc*fact]]), [0,(EpsR[1]-kap_)*fact], CR, [Eps_[0], sig[0]] 
		}
		else {								// else: / return [sig[0],0], array([[CC[0,0],0],[0,0]]), [Eps_[0], sig[0]]
			*(sig) = sigD[0]; // (1 - D)* sig0[0] + eta * Veps[0];
			*(sig+1) = 0.;
			*(MatM+0) = xxx*(1.-nu)-CD[0][0]+zz; *(MatM+1) = 0.; 
			*(MatM+2) = 0.;                      *(MatM+3) = 0.; 
		}
	}
	else if (ElemDim==2) {
		if (RType==1) {
			ccc = 0.5*RegPar*RegPar;		//              ccc = 0.5*self.RegPar**2  /               fact = 0.5*self.Emod 
            fact = 0.5*Emod;
			if (ElemPlSt) { 				return 122;			}
			else { 				return 123;			}
		}
		else {
			if (CalcType == 2) {
				if (ElemPlSt) {
					ff0 = (xxnu - CD[2][0]) / (x1nu - CD[2][2] + zz);
					ff1 = (xxnu - CD[2][1]) / (x1nu - CD[2][2] + zz);
					ff5 = (-CD[2][5]) / (x1nu - CD[2][2] + zz);
					*(MatM + 0) = (x1nu - CD[0][0] + zz) - (xxnu - CD[0][2]) * ff0; *(MatM + 1) = (xxnu - CD[0][1]) - (xxnu - CD[0][2]) * ff1; *(MatM + 2) = (-CD[0][5]) - (xxnu - CD[0][2]) * ff5;
					*(MatM + 3) = (xxnu - CD[1][0]) - (xxnu - CD[1][2]) * ff0; *(MatM + 4) = (x1nu - CD[1][1] + zz) - (xxnu - CD[1][2]) * ff1; *(MatM + 5) = (-CD[1][5]) - (xxnu - CD[1][2]) * ff5;
					*(MatM + 6) = (-CD[5][0]) - (-CD[5][2]) * ff0; *(MatM + 7) = (-CD[5][1]) - (-CD[5][2]) * ff1; *(MatM + 8) = (x2nu - CD[5][5] + zz) - (-CD[5][2]) * ff5;
				}
				else {
					*(MatM + 0) = (x1nu - CD[0][0] + zz);                     *(MatM + 1) = (xxnu - CD[0][1]);                     *(MatM + 2) = (-CD[0][5]);
					*(MatM + 3) = (xxnu - CD[1][0]);                     *(MatM + 4) = (x1nu - CD[1][1] + zz);                     *(MatM + 5) = (-CD[1][5]);
					*(MatM + 6) = (-CD[5][0]), * (MatM + 7) = (-CD[5][1]);                     *(MatM + 8) = (x2nu - CD[5][5] + zz);
				}
			}
			*(sig + 0) = sigD[0];
			*(sig + 1) = sigD[1];
			*(sig + 2) = sigD[5];         // voigt notation broken down to 2D
		}
	}
	else if (ElemDim == 3) {
		if (RType == 1) { return 124; }
		else {
			*(MatM + 0) = x1nu - CD[0][0] + zz; *(MatM + 1) = xxnu - CD[0][1];    *(MatM + 2) = xxnu - CD[0][2];    *(MatM + 3) = -CD[0][3];    *(MatM + 4) = -CD[0][4];    *(MatM + 5) = -CD[0][5];
			*(MatM + 6) = xxnu - CD[1][0];    *(MatM + 7) = x1nu - CD[1][1] + zz; *(MatM + 8) = xxnu - CD[1][2];    *(MatM + 9) = -CD[1][3];    *(MatM + 10) = -CD[1][4];    *(MatM + 11) = -CD[1][5];
			*(MatM + 12) = xxnu - CD[2][0];    *(MatM + 13) = xxnu - CD[2][1];    *(MatM + 14) = x1nu - CD[2][2] + zz; *(MatM + 15) = -CD[2][3];    *(MatM + 16) = -CD[2][4];    *(MatM + 17) = -CD[2][5];
			*(MatM + 18) = -CD[3][0];    *(MatM + 19) = -CD[3][1];    *(MatM + 20) = -CD[3][2];    *(MatM + 21) = x2nu - CD[3][3] + zz; *(MatM + 22) = -CD[3][4];    *(MatM + 23) = -CD[3][5];
			*(MatM + 24) = -CD[4][0];    *(MatM + 25) = -CD[4][1];    *(MatM + 26) = -CD[4][2];    *(MatM + 27) = -CD[4][3];    *(MatM + 28) = x2nu - CD[4][4] + zz; *(MatM + 29) = -CD[4][5];
			*(MatM + 30) = -CD[5][0];    *(MatM + 31) = -CD[5][1];    *(MatM + 32) = -CD[5][2];    *(MatM + 33) = -CD[5][3];    *(MatM + 34) = -CD[5][4];   *(MatM + 35) = x2nu - CD[5][5] + zz;
			for (i = 0; i < 6; i++) { *(sig + i) = sigD[i]; }
		}
	}
	else if (ElemDim == 4) {					// axisymmetric CP4
		if (RType == 1) {
			return 125;
		}
		else {
			*(MatM + 0)  = x1nu - CD[0][0] + zz; *(MatM + 1)  = xxnu - CD[0][1];      *(MatM + 2)  = xxnu - CD[0][2];        *(MatM + 3)  =      - CD[0][5];
			*(MatM + 4)  = xxnu - CD[1][0];      *(MatM + 5)  = x1nu - CD[1][1] + zz; *(MatM + 6)  = xxnu - CD[1][2];        *(MatM + 7)  =      - CD[1][5];
			*(MatM + 8)  = xxnu - CD[2][0];      *(MatM + 9)  = xxnu - CD[2][1];      *(MatM + 10) = x1nu - CD[2][2] + zz;   *(MatM + 11) =      - CD[2][5];
			*(MatM + 12) =      - CD[5][0];      *(MatM + 13) =       -CD[5][1];      *(MatM + 14) =      - CD[5][2];        *(MatM + 15) = x2nu - CD[5][5] + zz;
			*(sig + 0) = sigD[0];
			*(sig + 1) = sigD[1];
			*(sig + 2) = sigD[2];
			*(sig + 3) = sigD[5];         // voigt notation broken down to 2D
		}
	}
	else if (ElemDim == 21) {
		*(MatM + 0) = x1nu - CD[0][0] + zz; *(MatM + 1) = xxnu - CD[0][1];    *(MatM + 2) = xxnu - CD[0][2];    *(MatM + 3) = -CD[0][3];    *(MatM + 4) = -CD[0][4];    *(MatM + 5) = -CD[0][5];
		*(MatM + 6) = xxnu - CD[1][0];    *(MatM + 7) = x1nu - CD[1][1] + zz; *(MatM + 8) = xxnu - CD[1][2];    *(MatM + 9) = -CD[1][3];    *(MatM + 10) = -CD[1][4];    *(MatM + 11) = -CD[1][5];
		*(MatM + 12) = 0.;               *(MatM + 13) = 0.;               *(MatM + 14) = 0.;               *(MatM + 15) = 0.;               *(MatM + 16) = 0.;               *(MatM + 17) = 0.;
		*(MatM + 18) = -CD[3][0];    *(MatM + 19) = -CD[3][1];    *(MatM + 20) = -CD[3][2];    *(MatM + 21) = x2nu - CD[3][3] + zz; *(MatM + 22) = -CD[3][4];    *(MatM + 23) = -CD[3][5];
		*(MatM + 24) = -CD[4][0];    *(MatM + 25) = -CD[4][1];    *(MatM + 26) = -CD[4][2];    *(MatM + 27) = -CD[4][3];    *(MatM + 28) = x2nu - CD[4][4] + zz; *(MatM + 29) = -CD[4][5];
		*(MatM + 30) = -CD[5][0];    *(MatM + 31) = -CD[5][1];    *(MatM + 32) = -CD[5][2];    *(MatM + 33) = -CD[5][3];    *(MatM + 34) = -CD[5][4];   *(MatM + 35) = x2nu - CD[5][5] + zz;
		for (i = 0; i < 6; i++) { *(sig + i) = sigD[i]; }//  (1 - D)* sig0[i]; }
	}
	else { return 125; }       // else: raise NameError ("ConFemMaterials::Isodam.Sig: not implemented for this element type") 

	return 110;
}

int Rankine( bool PrinStrains, int ElemDim, double *Eps_, double *sig0, double *eVec, double *laMax)
{
	int rc;
	double la[3]; // , vv[3];
	if (PrinStrains) {
//		rc = eigJacobiSym(Eps_, la, vv);  // 	# principal strains(EvFlag = True), stresses   		laMax = zeros((1), dtype = float) 		eVec = zeros((3), dtype = float)
		rc = eigJacobiSym(Eps_, la, eVec);  // 	# principal strains(EvFlag = True), stresses   		laMax = zeros((1), dtype = float) 		eVec = zeros((3), dtype = float)
	}
	else {
//		rc = eigJacobiSym(sig0, la, vv);  // 	# principal strains(EvFlag = True), stresses   		laMax = zeros((1), dtype = float) 		eVec = zeros((3), dtype = float)
		rc = eigJacobiSym(sig0, la, eVec);  // 	# principal strains(EvFlag = True), stresses   		laMax = zeros((1), dtype = float) 		eVec = zeros((3), dtype = float)
	}
	if (rc>110) { return rc; }
	*laMax = la[0];	//  ila = 0  # largest principal value   \  if la[1]>la[0]: ila=1 \  if la[2]>la[1] and la[2]>la[0]: ila=2 \  laMax = la[ila]
//	*(eVec+0) = *(vv + 0);
//	*(eVec+1) = *(vv + 1);
//	*(eVec+2) = *(vv + 2);
	if ((ElemDim == 2) && (fabs(eVec[2]) > ZeroD)) {   //		if Elem.dim == 2 and abs(eVec[2])>ZeroD:  laMax = 0.    # lateral extension to plane, no lateral cracks allowed
		*laMax = 0.;
	}
	if (eVec[0] + eVec[1] + eVec[2] < 0.) { // # directions indicate the same plane (approximately for different IPs) but might have opposite directions. Has to be avoided 
		*(eVec+0) = -*(eVec+0);
		*(eVec+1) = -*(eVec+1);
		*(eVec+2) = -*(eVec+2);
	}
	return 100;
}
int RankineUpd( double *sig0, double *eVec, double *ElemStateVarN, double *laMax, int off)
{
	double tt[3];
	tt[0] = *(sig0+0) * *(eVec+0) + *(sig0+5) * *(eVec+1) + *(sig0+4) * *(eVec+2);
	tt[1] = *(sig0+5) * *(eVec+0) + *(sig0+1) * *(eVec+1) + *(sig0+3) * *(eVec+2);
	tt[2] = *(sig0+4) * *(eVec+0) + *(sig0+3) * *(eVec+1) + *(sig0+2) * *(eVec+2);
	//
	if (sqrt(tt[0] * tt[0] + tt[1] * tt[1] + tt[2] * tt[2]) < ZeroD) *laMax = 0.;
	*(ElemStateVarN + off + 6) = *laMax;
	*(ElemStateVarN + off + 7) = *(eVec+0);
	*(ElemStateVarN + off + 8) = *(eVec+1);
	*(ElemStateVarN + off + 9) = *(eVec+2);
	// tractions - after viscous contribution
	*(ElemStateVarN + off + 10) = tt[0];  //sig0[0] * eVec[0] + sig0[5] * eVec[1] + sig0[4] * eVec[2];  // tt[0];
	*(ElemStateVarN + off + 11) = tt[1];  // sig0[5] * eVec[0] + sig0[1] * eVec[1] + sig0[3] * eVec[2];  // tt[1];
	*(ElemStateVarN + off + 12) = tt[2];  // sig0[4] * eVec[0] + sig0[3] * eVec[1] + sig0[2] * eVec[2];  // tt[2];
	return 100;
}

inline int State1C( double ww1_1, double MatML[][3], double* epsL, double* sigL, int i, double* ElemStateVarN,
	         double selffct, double ds1, double etaH, double epsu, double epst, double nu, double xi, double Emod, double epst2, double eps_c11_k, double etaB, double deps_c11_k, double* ww) {			//						def State1C(i, ww) : # single crack - loading
	;
	//	_MatML = zeros((3, 3), dtype = float)
	//	_sigL = zeros((3), dtype = float)
	MatML[0][0] = -selffct / ds1 * (-etaH * epsu + 1 + etaH * epst);
	MatML[0][1] = -selffct / ds1 * (etaH * epst * nu - etaH * epsu * nu - etaH * epst * xi * nu + nu + etaH * epsu * xi * nu - xi * nu);
	MatML[1][0] = -Emod / ds1 * (-nu * etaH * epsu * epst + nu * etaH * epst2 + nu * epst);
	MatML[1][1] = -Emod / ds1 * (-xi * epsu + etaH * epsu * xi * epst - etaH * epst2 * xi + etaH * epst2 - etaH * epsu * epst + epst);
	sigL[0] = -selffct / ds1 * (xi * etaH * eps_c11_k * epsu + epsL[0] - xi * epsu + etaH * epsu * xi * nu * epsL[1] + etaH * epst * epsL[0] + etaH * epst * nu * epsL[1] + xi * etaB * deps_c11_k * epsu - etaH * epsu * nu * epsL[1] - xi * etaB * deps_c11_k * epst - etaH * epsu * epsL[0] + nu * epsL[1] - etaH * epst * xi * nu * epsL[1] - xi * nu * epsL[1] - xi * etaH * eps_c11_k * epst);
	sigL[1] = -Emod / ds1 * (epsL[1] * epst + nu * epsL[0] * epst - nu * epsu * xi * epst - epsL[1] * xi * epsu + epsL[1] * etaH * epst2 - epsL[1] * etaH * epsu * epst + nu * xi * etaH * eps_c11_k * epsu * epst + nu * etaH * epst2 * epsL[0] + nu * xi * etaB * deps_c11_k * epsu * epst - nu * xi * etaB * deps_c11_k * epst2 - nu * etaH * epsu * epsL[0] * epst - nu * xi * etaH * eps_c11_k * epst2 + epsL[1] * etaH * epsu * xi * epst - epsL[1] * etaH * epst2 * xi);
	*(ElemStateVarN + i + 2 ) = ww1_1;				//	Elem.StateVarN[ipI, i + 2] = ww
	*(ElemStateVarN + i + 3 ) = sigL[i];		//	Elem.StateVarN[ipI, i + 3] = _sigL[i]
	ww[0] = ww1_1;
	return 1;//						return 1, ww, _sigL, _MatML
}
inline int State2C( double ww1_2, double MatML[][3], double* epsL, double* sigL,
	         double selffct, double ds2, double tP1, double selfbw_, double wP1, double etaH, double nu, double xi, double Emod, double epst, double eps_c11_k, double etaB, double deps_c11_k, double* ww ) {				//						def State2C(ww) : # single crack - unloading
	;
	// _MatML = zeros((3, 3), dtype = float)
	//_sigL = zeros((3), dtype = float)
	MatML[0][0] = -selffct / ds2 * (-tP1 * selfbw_ - selffct * wP1 * etaH);
	MatML[0][1] = -selffct / ds2 * (-selffct * wP1 * etaH * nu + selffct * wP1 * etaH * xi * nu + selfbw_ * tP1 * xi * nu - nu * tP1 * selfbw_);
	MatML[1][0] = -Emod / ds2 * (-nu * selfbw_ * tP1 * epst - nu * selffct * wP1 * etaH * epst);
	MatML[1][1] = -Emod / ds2 * (-xi * selffct * wP1 - tP1 * selfbw_ * epst + selfbw_ * tP1 * xi * epst + selffct * wP1 * etaH * xi * epst - selffct * wP1 * etaH * epst);
	sigL[0] = -selffct / ds2 * (-selfbw_ * epsL[0] * tP1 + xi * selffct * wP1 * etaH * eps_c11_k + xi * selffct * wP1 * etaB * deps_c11_k - selffct * wP1 * etaH * nu * epsL[1] - selffct * wP1 * etaH * epsL[0] + selffct * wP1 * etaH * xi * nu * epsL[1] + selfbw_ * tP1 * xi * nu * epsL[1] - tP1 * selfbw_ * nu * epsL[1]);
	sigL[1] = -Emod / ds2 * (-nu * selfbw_ * epsL[0] * tP1 * epst + nu * xi * selffct * wP1 * etaH * eps_c11_k * epst + nu * xi * selffct * wP1 * etaB * deps_c11_k * epst - nu * selffct * wP1 * etaH * epsL[0] * epst - epsL[1] * selffct * xi * wP1 - epsL[1] * tP1 * selfbw_ * epst + epsL[1] * selfbw_ * tP1 * xi * epst + epsL[1] * selffct * wP1 * etaH * xi * epst - epsL[1] * selffct * wP1 * etaH * epst);
	ww[0] = ww1_2;
	return 2;				//						return 2, ww, _sigL, _MatML
}
inline int State3C( double MatML[][3], double* epsL, double* sigL, 
	        double Emod, double nu, double* ww) {			//						def State3C() : # single crack - closure
	//	_MatML = zeros((3, 3), dtype = float)
	//	_sigL = zeros((3), dtype = float)
	double xxx = Emod / (1 - nu * nu);
	MatML[0][0] = xxx;
	MatML[0][1] = xxx * nu;
	MatML[1][0] = xxx * nu;
	MatML[1][1] = xxx;
	sigL[0] = MatML[0][0] * epsL[0] + MatML[0][1] * epsL[1];
	sigL[1] = MatML[1][0] * epsL[0] + MatML[1][1] * epsL[1];
	ww[0] = 0.;
	return 3;		//						return 3, 0., _sigL, _MatML
}
inline int State4C( double ww1_4, double MatML[][3], double* epsL, double* sigL, double* ElemStateVarN,
			 double selffct, double ds4,  double xi, double etaH, double eps_c11_k, double etaB, double deps_c11_k, double nu, double rho, double Emod, double epst, double* ww) { //	def State4C(ww) : # single cracking break through after softening
	//	_MatML = zeros((3, 3), dtype = float)
	//	_sigL = zeros((3), dtype = float)
	sigL[0]     = -selffct / ds4 * (xi * etaH * eps_c11_k + xi * etaB * deps_c11_k - etaH * nu * epsL[1] - epsL[0] * etaH + etaH * xi * nu * epsL[1] - xi * rho);
	MatML[0][0] = selffct * etaH / ds4;
	MatML[0][1] = selffct * (etaH * nu - etaH * xi * nu) / ds4;
	sigL[1]     = -Emod / ds4 * (nu * xi * etaH * eps_c11_k * epst + nu * xi * etaB * deps_c11_k * epst - nu * epsL[0] * etaH * epst - epsL[1] * xi - epsL[1] * etaH * epst + epsL[1] * etaH * xi * epst - nu * xi * rho * epst);
	MatML[1][0] = selffct * nu * etaH / ds4;
	MatML[1][1] = Emod * (xi + etaH * epst - etaH * xi * epst) / ds4;
	*(ElemStateVarN + 2) = ww1_4; // Elem.StateVarN[ipI, 2] = ww
	ww[0] = ww1_4;
	return 4;	//						return  4, ww, _sigL, _MatML
}
//						# dual cracking - assumes that first crack is in state 4
inline int State1C2(double ww2_1, double MatML[][3], double* epsL, double* sigL, double* ElemStateVarN, 
			 double selffct, double dd4, double xi, double etaH, double eps_c11_k, double etaB, double deps_c11_k, double rho, double dd1, double epsu, double eps_c21_k, double epst, double deps_c21_k, double* ww) {	//						def State1C2(ww) : # dual cracking - loading - 1st crack with full crack
	;
	//	_MatML = zeros((3, 3), dtype = float)
	//	_sigL = zeros((3), dtype = float)
	sigL[0] = selffct / dd4 * (xi * etaH * eps_c11_k + xi * etaB * deps_c11_k - etaH * epsL[0] - xi * rho);
	sigL[1] = -selffct / dd1 * (xi * epsu - xi * etaH * eps_c21_k * epsu + xi * etaH * eps_c21_k * epst - xi * etaB * deps_c21_k * epsu + xi * etaB * deps_c21_k * epst - epsL[1] - etaH * epst * epsL[1] + etaH * epsu * epsL[1]);
	MatML[0][0] = -selffct / dd4 * etaH; //             # !!!!!dd4 presumably negative
	MatML[1][1] = -selffct / dd1 * (-1 - etaH * epst + etaH * epsu);
	*(ElemStateVarN + 7) = ww2_1;
	*(ElemStateVarN + 8) = sigL[1];
	ww[1] = ww2_1;
	return 1;		//						return 1, ww, _sigL, _MatML
}
/*
*/
inline int State2C2(double ww2_2, double MatML[][3], double* epsL, double* sigL, 
	         double selffct, double dd4, double xi, double etaH, double eps_c11_k, double etaB, double deps_c11_k, double rho, double dd2, double selfbw_, double tP2, double wP2, double* ww) {	//						def State2C2(ww) : # dual cracking - unloading
	;
	//	_MatML = zeros((3, 3), dtype = float)
	//	_sigL = zeros((3), dtype = float)
	sigL[0] = selffct / dd4 * (xi * etaH * eps_c11_k + xi * etaB * deps_c11_k - etaH * epsL[0] - xi * rho);
	sigL[1] = selffct / dd2 * (-selfbw_ * tP2 * epsL[1] + xi * selffct * wP2 * etaH * eps_c11_k + xi * selffct * wP2 * etaB * deps_c11_k - selffct * wP2 * etaH * epsL[1]);
	MatML[0][0] = -selffct / dd4 * etaH;
	MatML[1][1] = selffct / dd2 * (-tP2 * selfbw_ - selffct * wP2 * etaH);
	ww[1] = ww2_2;
	return 2;			//						return 2, ww, _sigL, _MatML
}
inline int State3C2( double MatML[][3], double* epsL, double* sigL,
			  double selffct, double dd4, double xi, double etaH, double eps_c11_k, double etaB, double deps_c11_k, double rho, double Emod, double* ww) {	//						def State3C2() : # dual cracking - crack closure
	//	_MatML = zeros((3, 3), dtype = float)
	//	_sigL = zeros((3), dtype = float)
	sigL[0] = selffct / dd4 * (xi * etaH * eps_c11_k + xi * etaB * deps_c11_k - etaH * epsL[0] - xi * rho);
	sigL[1] = MatML[1][1] * epsL[1];
	MatML[0][0] = -selffct / dd4 * etaH;
	MatML[1][1] = Emod;
	ww[1] = 0.;
	return 3;		//						return 3, 0., _sigL, _MatML
}
inline int State4C2(double ww2_4, double MatML[][3], double* epsL, double* sigL, double* ElemStateVarN, 
			 double selffct, double dd4, double xi, double etaH, double eps_c11_k, double etaB, double deps_c11_k, double rho, double eps_c21_k, double deps_c21_k, double* ww) {		//	def State4C2(ww) : # dual cracking - full crack
	//	_MatML = zeros((3, 3), dtype = float)
	//	_sigL = zeros((3), dtype = float)
	sigL[0] = selffct / dd4 * (xi * etaH * eps_c11_k + xi * etaB * deps_c11_k - etaH * epsL[0] - xi * rho);
	sigL[1] = selffct / dd4 * (xi * etaH * eps_c21_k + xi * etaB * deps_c21_k - etaH * epsL[1] - xi * rho);
	MatML[0][0] = -selffct / dd4 * etaH;
	MatML[1][1] = MatML[0][0];
	*(ElemStateVarN + 7) = ww2_4;	//						Elem.StateVarN[ipI, 7] = ww
	ww[1] = ww2_4;
	return 4;			//						return  4, ww, _sigL, _MatML
}

int ElasticLTC2(int ElemDim, bool ElemPlSt, double *ElemStateVar, int dim_St, double *ElemStateVarN, int dim_StN,
	double *Eps, int dim_E, double *sig, int dim_si, double *MatM, int dim_Ma, double* ww, int dim_ww, double selfnu, double selfEmod, double *Dps, int dim_D, double Dt,
	double selfrho, double selfwcr_, double selfbw_, double selfepsct, double ElemLch_, double selfCrackVisc, double selffct, double selfepscu,
	double *DataOut, int dim_DataOut)
{
	double Eps_[6]; 
	int i; 
	double nu, Emod, pEps[3], phe, rho, epsL[3], epsLP[3], wP1, tP1, wP2, tP2, wcr, xi, etaB, etaH, epsu, epst, nu2,
			epst2, ds1, ds2, ds4, eps_c11_k, deps_c11_k, ww1_1, ww1_2, ww1_4, dd1, dd2, dd4, eps_c21_k, deps_c21_k, ww2_1, ww2_2, ww2_4, cc, pSig[3];
//	double nori[2], nact1[2], nact2[2], phe_;
	int State, State2, StateP1, StateP2, State_, State2_;
	double Trans[3][3];
	double MatML[3][3] = { {0.,0.,0.},{0.,0.,0.},{0.,0.,0.} };
	double sigL[3] = { 0.,0.,0. };
	if (ElemDim == 2) {  		// if Elem.dim == 2:
		if (!ElemPlSt) return 10;		//	if not Elem.PlSt : raise NameError("ConFemMaterials::ElasticLT.sig: ElasticLT not yet defined for plane strain")
		Eps_[0] = *Eps;		//	Eps_ = Eps
		Eps_[1] = *(Eps + 1);
		Eps_[2] = *(Eps + 2);
	}
	else if (ElemDim == 21) {  // //	elif Elem.dim == 21 :
		Eps_[0] = *Eps;		//		Eps_ = array([Eps[0], Eps[1], Eps[5]])
		Eps_[1] = *(Eps + 1);
		Eps_[2] = *(Eps + 5);
	}
	else { return 20; }			// 	raise NameError("CaeFemMaterials::ElasticLT.Sig: element type not implemented for this material")
	nu = selfnu;				//	nu = self.nu                                            # Poisson's ratio
	Emod = selfEmod;			//	Emod = self.Emod
	PrinCLT_1(Eps_[0], Eps_[1], 0.5*Eps_[2], pEps); //	pep, phe, pep_ = PrinCLT_(Eps_[0], Eps_[1], 0.5 * Eps_[2]) # principal strains, largest value, corr.direction, lower value
	phe = pEps[1];
	ElemStateVarN[9]  = pEps[0]; //	Elem.StateVarN[ipI, 9] = pep                             # 1st larger principal strain
	ElemStateVarN[10] = pEps[2]; //	Elem.StateVarN[ipI, 10] = pep_                            # 2nd principal strain
	State = int(*ElemStateVarN); //	State = int(Elem.StateVarN[ipI, 0])                      # current state of actual time step
//	if (State == 0) { ElemStateVarN[5] = pEps[1]; }//	if State == 0: Elem.StateVarN[ipI, 5] = phe              # save for updated potential 1st crack direction
	// # cracked state
	if (State > 0) { //	if State > 0:
		rho = selfrho;		//		rho = self.rho # 0.005 #0.02 # 0.01 #0.                # related residual strength with full crack
		//		TOL = -1.e-15                                       #
		State2 = int(*(ElemStateVarN + 6)); //		State2 = int(Elem.StateVarN[ipI, 6])                 # current state 2 of actual time step
		if (State2 > 0) { nu = 0.; }//		if State2 > 0: nu = 0.
		// # transformation to local system
//		Trans = array([[cos(phe)**2, sin(phe)**2, cos(phe) * sin(phe)],[sin(phe)**2, cos(phe)**2, -cos(phe) * sin(phe)],[-2 * cos(phe) * sin(phe), 2 * cos(phe) * sin(phe), cos(phe)**2 - sin(phe)**2]] )# transformation matrix for strains from global to local
		Trans[0][0] = pow(cos(phe), 2.0); // transformation matrix for strains from global to local
		Trans[0][1] = pow(sin(phe), 2.0);
		Trans[0][2] = cos(phe) * sin(phe);
		Trans[1][0] = Trans[0][1];			//  pow(sin(phe),2.0);
		Trans[1][1] = Trans[0][0];			//  pow(cos(phe),2.0);
		Trans[1][2] = -Trans[0][2];			// -cos(phe)*sin(phe);
		Trans[2][0] = -2 * cos(phe) * sin(phe);
		Trans[2][1] = -Trans[2][0];			// -2*cos(phe)*sin(phe);
		Trans[2][2] = pow(cos(phe), 2.0) - pow(sin(phe), 2.0);
		for (i = 0; i < 3; i++) { epsL[i] = Trans[i][0] * Eps_[0] + Trans[i][1] * Eps_[1] + Trans[i][2] * Eps_[2]; } // local strain //		epsL = dot(Trans, Eps_) # local strain
		//		# state of previous time step
		epsLP[0] = *(ElemStateVar + 9);		//		epsLP = array([Elem.StateVar[ipI, 9], Elem.StateVar[ipI, 10], 0.]) # local strain previous step
		epsLP[1] = *(ElemStateVar + 10);
		StateP1 = int(*(ElemStateVar + 0));	//		StateP1 = Elem.StateVar[ipI, 0]                      # final state of last time step 1st crack
		wP1 = *(ElemStateVar + 2);			//		wP1 = Elem.StateVar[ipI, 2]                          # largest crack ever reached of last time step first crack
		tP1 = *(ElemStateVar + 3);			//		tP1 = Elem.StateVar[ipI, 3]                          # corresponding crack traction
		StateP2 = int(*(ElemStateVar + 6));	//		StateP2 = Elem.StateVar[ipI, 6]                      # final state of last time step 2nd crack
		wP2 = *(ElemStateVar + 7);			//		wP2 = Elem.StateVar[ipI, 7]                          # largest crack ever reached of last time step first crack
		tP2 = *(ElemStateVar + 8);			//		tP2 = Elem.StateVar[ipI, 8]                          # corresponding crack traction
		//		# check for largest principal stress direction relative to initial 1st crack direction
//		nori[0] = cos(*(ElemStateVarN + 5));//		nori = array([cos(Elem.StateVarN[ipI, 5]), sin(Elem.StateVarN[ipI, 5])]) # initial 1st crack direction
//		nori[1] = sin(*(ElemStateVarN + 5));
//		nact1[0] = cos(phe);				//		nact1 = array([cos(phe), sin(phe)])                  # current direction of largest principal stress /strain?
//		nact1[1] = sin(phe);
//		nact2[0] = cos(phe + 0.5*pi);		//		nact2 = array([cos(phe + 0.5 * pi), sin(phe + 0.5 * pi)])    # current direction of minor principal stress / strain?
//		nact2[1] = sin(phe + 0.5*pi);
//		if ( abs(nori[0]*nact1[0]+nori[1]*nact1[1]) < abs(nori[0]*nact2[0]+nori[1]*nact2[1]) ) { //		if abs(dot(nori, nact1)) < abs(dot(nori, nact2)) : # larger principal stress direction nact 1 far away from initial crack direction-- switch to nact2
//			if (nori[0] * nact2[0] + nori[1] * nact2[1] > 0.) { phe_ = phe + 0.5 * pi; } 		 //			if dot(nori, nact2) > 0.: phe_ = phe + 0.5 * pi        # same orientation of nori and nact2
//			else                                              { phe_ = phe - 0.5 * pi; }		//			else:                  phe_ = phe - 0.5 * pi        # opposite orientation of nori and nact
//			phe = phe_;																			//			phe = phe_
//		}
		//		# auxiliary values
		wcr = selfwcr_ - rho * (selfwcr_ - selfbw_ * selfepsct); //		wcr = self.wcr_ - rho * (self.wcr_ - self.bw_ * self.epsct)                                # effective critical crack width
		xi = selfbw_ / ElemLch_;			//		xi = self.bw_ / Elem.Lch_                           # related crack band width
		etaB = selfCrackVisc / selffct;		//		etaB = self.CrackVisc / self.fct                      # related crack viscosity
		etaH = 2. * etaB / Dt;				//		etaH = 2. * etaB / Dt                                   # related crack viscosity 2
		epsu = selfepscu;					//		epsu = self.epscu                                   # failure strain
		epst = selfepsct;					//		epst = self.epsct                                   # tensile strength strain
		nu2 = nu * nu;						//		nu2 = nu * nu
		epst2 = epst * epst;				//		epst2 = epst * epst
		//		# auxiliary values - 1st crack
		ds1 = xi * epsu + nu2 * epst - etaH * epst2 + etaH * epsu * epst - etaH * epsu * xi * epst + etaH * epst2 * xi - xi * nu2 * epst + etaH * epsu * xi * nu2 * epst + etaH * epst2 * nu2 - etaH * epsu * nu2 * epst - etaH * epst2 * xi * nu2 - epst;
		ds2 = selffct * wP1 * xi + tP1 * selfbw_ * epst - tP1 * selfbw_ * nu2 * epst - selfbw_ * tP1 * xi * epst + selfbw_ * tP1 * xi * nu2 * epst - selffct * wP1 * etaH * xi * epst + selffct * wP1 * etaH * xi * nu2 * epst + selffct * wP1 * etaH * epst - selffct * wP1 * etaH * nu2 * epst;
		ds4 = xi + etaH * epst - nu2 * etaH * epst - xi * etaH * epst + nu2 * xi * etaH * epst;
		eps_c11_k  = *(ElemStateVar + 11) / selfbw_;	//		eps_c11_k = Elem.StateVar[ipI, 11] / self.bw_          # 1st crack strain of previous time step for viscous contributions
		deps_c11_k = *(ElemStateVar + 13) / selfbw_;	//		deps_c11_k = Elem.StateVar[ipI, 13] / self.bw_          # 1st crack strain velocity of previous time step for viscous contributions
		//		# predictors crack width 1st crack
		ww1_1 = selfbw_ / ds1 * (-epsu * epst - epsL[0] * epst + epsL[0] * epsu - nu2 * etaH * eps_c11_k * epsu * epst - xi * etaH * eps_c11_k * epsu * epst - xi * etaB * deps_c11_k * epsu * epst + nu2 * xi * etaH * eps_c11_k * epsu * epst - epsu * xi * nu2 * epst + nu2 * xi * etaB * deps_c11_k * epsu * epst - nu * epsL[1] * epst + xi * nu * epsL[1] * epst + nu * epsL[1] * epsu + nu2 * etaH * eps_c11_k * epst2 - etaH * eps_c11_k * epst2 - nu * epsL[1] * xi * epsu - etaB * deps_c11_k * epst2 + nu2 * etaB * deps_c11_k * epst2 + xi * etaH * eps_c11_k * epst2 + xi * etaB * deps_c11_k * epst2 - nu2 * xi * etaB * deps_c11_k * epst2 - nu2 * xi * etaH * eps_c11_k * epst2 + epsu * nu2 * epst + xi * epsu * epst - nu2 * etaB * deps_c11_k * epsu * epst + etaH * eps_c11_k * epsu * epst + etaB * deps_c11_k * epsu * epst);
		if (abs(ds2) > ZeroD) { ww1_2 = selfbw_ * wP1 * selffct / ds2 * (eps_c11_k * etaH * epst + deps_c11_k * etaB * epst - nu2 * etaB * deps_c11_k * epst - nu2 * etaH * eps_c11_k * epst + nu * epsL[1] + nu2 * xi * etaB * deps_c11_k * epst + nu2 * xi * etaH * eps_c11_k * epst - nu * epsL[1] * xi - xi * etaH * eps_c11_k * epst - xi * etaB * deps_c11_k * epst + epsL[0]);}
		else				  {	ww1_2 = 0.; }			//		else:              ww1_2 = 0.
		ww1_4 = selfbw_ / ds4 * (eps_c11_k * etaH * epst + deps_c11_k * etaB * epst - nu2 * etaH * eps_c11_k * epst - nu2 * etaB * deps_c11_k * epst + nu * epsL[1] - xi * etaH * eps_c11_k * epst - xi * etaB * deps_c11_k * epst + nu2 * xi * etaH * eps_c11_k * epst + nu2 * xi * etaB * deps_c11_k * epst - nu * epsL[1] * xi + epsL[0] + rho * epst * (nu2 - 1 + xi - nu2 * xi));
		//		#
		State_  = -1,						//		State_ = None                                       # initial value current states
		State2_ = 0;						//		State2_ = 0
		// # actual state of cracking depending on previous state and crack width predictors
		if ( StateP1 >= 4 )  {				// if StateP1 >= 4:                                      # state 4 or 5 - maybe single or dual cracking as direction 1 decoupled from 2 for the following(nu = 0)
			if (ww1_4 > 0. ) {				//if ww1_4 > 0:                                     # state 4 open crack loading or unloading
				State_ = State4C(ww1_4, MatML, epsL, sigL, ElemStateVarN, selffct,ds4,xi,etaH,eps_c11_k,etaB,deps_c11_k,nu,rho,Emod,epst, ww); // State_, ww1, sigL, MatML = State4C(ww1_4)          #
				if ( State2 > 0 )  { //				if State2 > 0:                                # dual cracking crack 2
					//	# auxiliary values - 2nd crack
					dd1 = -xi * epsu - etaH * epsu * epst + epst + etaH * epst2 - etaH * epst2 * xi + etaH * epsu * xi * epst;
					dd2 = -tP2 * selfbw_ * epst + selfbw_ * tP2 * xi * epst - selffct * wP2 * xi + selffct * wP2 * etaH * xi * epst - selffct * wP2 * etaH * epst;
					dd4 = -xi + etaH * epst * xi - etaH * epst;
					eps_c21_k = *(ElemStateVar + 12) / selfbw_;        //  # 2nd crack strain of previous time step for viscous contributions
					deps_c21_k = *(ElemStateVar + 14) / selfbw_;       //   # 2nd crack strain velocity of previous time step for viscous contributions
					//					# predictors for crack width 2nd crack
					ww2_1 = -selfbw_ / dd1 * (-epsu * epst + etaH * eps_c21_k * epsu * epst - etaH * eps_c21_k * epst2 + etaB * deps_c21_k * epsu * epst - etaB * deps_c21_k * epst2 + epsu * xi * epst - xi * etaH * eps_c21_k * epsu * epst + xi * etaH * eps_c21_k * epst2 - xi * etaB * deps_c21_k * epsu * epst + xi * etaB * deps_c21_k * epst2 - epst * epsL[1] + epsu * epsL[1]);
					if (abs(dd2) > ZeroD) {
						ww2_2 = selfbw_ * wP2 * selffct / dd2 * (-etaH * eps_c21_k * epst - etaB * deps_c21_k * epst + xi * etaH * eps_c21_k * epst + xi * etaB * deps_c21_k * epst - epsL[1]);
					}
					else {
						ww2_2 = 0.;
					}
					ww2_4 = selfbw_ / dd4 * (-etaH * eps_c21_k * epst - etaB * deps_c21_k * epst + xi * etaH * eps_c21_k * epst + xi * etaB * deps_c21_k * epst - epsL[1]);
					//					#
					if ( StateP2 >= 4 ) {		//	if StateP2 >= 4 : # state 4 or 5
						if (ww2_4 > 0) { 
							State2_ = State4C2(ww2_4, MatML, epsL, sigL, ElemStateVarN, selffct, dd4, xi, etaH, eps_c11_k, etaB, deps_c11_k, rho, eps_c21_k, deps_c21_k, ww); 
						}						//	if   ww2_4 > 0:   State2_, ww2, sigL, MatML = State4C2(ww2_4)         # state 4 open crack loading or unloading
						else {
							State2_ = State3C2( MatML, epsL, sigL, selffct, dd4, xi, etaH, eps_c11_k, etaB, deps_c11_k, rho, Emod, ww);
							State2_ = 5;
						} 						//	else:             State2_, ww2, sigL, MatML = State3C2(); State2_ = 5  # state 5 crack closure
					}
					else {						//					else:                                   # crack 2 state 1 or 2 or 3
						if ( ww2_1 > wcr ) { State2_ = State4C2( ww2_4, MatML, epsL, sigL, ElemStateVarN, selffct, dd4, xi, etaH, eps_c11_k, etaB, deps_c11_k, rho, eps_c21_k, deps_c21_k, ww); }//	if   ww2_1 > wcr: State2_, ww2, sigL, MatML = State4C2(ww2_4)         # new state 4
						else if ( ( epsL[1] > epsLP[1] + TOL) && ( ww2_1 > wP2 ) ) { State2_ =	State1C2( ww2_1, MatML, epsL, sigL, ElemStateVarN, selffct, dd4, xi, etaH, eps_c11_k, etaB, deps_c11_k, rho, dd1,
							        epsu,  eps_c21_k, epst, deps_c21_k, ww ); } //	elif epsL[1] > epsLP[1] + TOL and ww2_1 > wP2:  State2_, ww2, sigL, MatML = State1C2(ww2_1)     # state 1 loading
						else if ( ww2_2 > 0 ) { State2_ = State2C2( ww2_2, MatML, epsL, sigL, selffct, dd4, xi, etaH, eps_c11_k, etaB, deps_c11_k, rho, dd2, selfbw_, tP2, wP2, ww ); }	//	elif ww2_2 > 0: State2_, ww2, sigL, MatML = State2C2(ww2_2)         # state 2 unloading
						else if ( ww2_2 <= 0 ) { State2_ = State3C2( MatML, epsL, sigL, selffct, dd4, xi, etaH, eps_c11_k, etaB, deps_c11_k, rho, Emod, ww );  }	//	elif ww2_2 <= 0: State2_, ww2, sigL, MatML = State3C2()              # state 3 crack closure
						else { return 30 ; } //						else: raise NameError("ConFemMaterials::ElasticLT.sig: crack exception 2")
					}
				}
				else { State2_ = 0; }
			}
			else {
				State_ = State3C( MatML, epsL, sigL, Emod, nu, ww );		//:  State_, ww1, sigL, MatML = State3C(); 
				State_ = 5;		//State_ = 5    # state 5 - crack closure
			}
		}
		else {		//		else:
			if  (ww1_1 > wcr ) { State_ = State4C( ww1_4, MatML, epsL, sigL, ElemStateVarN, selffct, ds4, xi, etaH, eps_c11_k, etaB, deps_c11_k, nu, rho, Emod, epst, ww); } //	if ww1_1 > wcr:  State_, ww1, sigL, MatML = State4C(ww1_4)        # starting state 4 loading beyond critical crack width
			else if ( (epsL[0] > epsLP[0] + TOL) && (ww1_1 > wP1) ) { State_ = State1C( ww1_1,  MatML, epsL, sigL, 0, ElemStateVarN, selffct, ds1, etaH, epsu, epst, nu, xi, Emod, epst2, eps_c11_k, etaB, deps_c11_k, ww); }//	elif epsL[0] > epsLP[0] + TOL and ww1_1 > wP1: State_, ww1, sigL, MatML = State1C(0, ww1_1)        # state 1 loading
			else if ( ww1_2 > 0. ) { State_ = State2C( ww1_2, MatML, epsL, sigL, selffct, ds2, tP1, selfbw_, wP1, etaH, nu, xi, Emod, epst, eps_c11_k, etaB, deps_c11_k, ww); }//	elif ww1_2 > 0: State_, ww1, sigL, MatML = State2C(ww1_2)        # state 2 unloading
			else if ( ww1_2 <= 0.) { State_ = State3C( MatML, epsL, sigL, Emod, nu, ww); }//	elif ww1_2 <= 0: State_, ww1, sigL, MatML = State3C() # state 3 crack closure
			else { return 40 ; }//			else :
		}
		//		# finish
		*(ElemStateVarN + 0) = State_;		//		Elem.StateVarN[ipI, 0] = State_
		*(ElemStateVarN + 6) = State2_;		//		Elem.StateVarN[ipI, 6] = State2_
		*(ElemStateVarN + 1) = max(sigL[0], sigL[1]); //		Elem.StateVarN[ipI, 1] = max(sigL[0], sigL[1])       # presumably relevant for dual cracking only
		*(ElemStateVarN + 11) = ww[0];		//		Elem.StateVarN[ipI, 11] = ww1                        # current crack width 1
		*(ElemStateVarN + 12) = ww[1]; //		Elem.StateVarN[ipI, 12] = ww2                        # current crack width 2
		*(ElemStateVarN + 13) = 2. * (ww[0] - *(ElemStateVar + 11)) / Dt - *(ElemStateVar + 13); //Elem.StateVarN[ipI, 13] = 2. * (ww1 - Elem.StateVar[ipI, 11]) / Dt - Elem.StateVar[ipI, 13]  # current crack width velocity 1
		*(ElemStateVarN + 14) = 2. * (ww[1] - *(ElemStateVar + 12)) / Dt - *(ElemStateVar + 14); //Elem.StateVarN[ipI, 14] = 2. * (ww2 - Elem.StateVar[ipI, 12]) / Dt - Elem.StateVar[ipI, 14]  # current crack width velocity 2
		//		#
		for (i = 0; i < 3; i++) { sig[i] = Trans[0][i] * sigL[0] + Trans[1][i] * sigL[1] + Trans[2][i] * sigL[2]; } // global stress	//		sig = dot(Trans.transpose(), sigL)                   # transformation of stress into global system
		*(MatM + 0) = Trans[0][0] * (MatML[0][0] * Trans[0][0] + MatML[0][1] * Trans[1][0]) + Trans[1][0] * (MatML[1][0] * Trans[0][0] + MatML[1][1] * Trans[1][0]); // 00	//		MatM = dot(Trans.transpose(), dot(MatML, Trans))     # transformation of tangent stiffness into global system
		*(MatM + 1) = Trans[0][0] * (MatML[0][0] * Trans[0][1] + MatML[0][1] * Trans[1][1]) + Trans[1][0] * (MatML[1][0] * Trans[0][1] + MatML[1][1] * Trans[1][1]); // 01
		*(MatM + 2) = Trans[0][0] * (MatML[0][0] * Trans[0][2] + MatML[0][1] * Trans[1][2]) + Trans[1][0] * (MatML[1][0] * Trans[0][2] + MatML[1][1] * Trans[1][2]); // 02
		*(MatM + 3) = Trans[0][1] * (MatML[0][0] * Trans[0][0] + MatML[0][1] * Trans[1][0]) + Trans[1][1] * (MatML[1][0] * Trans[0][0] + MatML[1][1] * Trans[1][0]); // 10
		*(MatM + 4) = Trans[0][1] * (MatML[0][0] * Trans[0][1] + MatML[0][1] * Trans[1][1]) + Trans[1][1] * (MatML[1][0] * Trans[0][1] + MatML[1][1] * Trans[1][1]); // 11
		*(MatM + 5) = Trans[0][1] * (MatML[0][0] * Trans[0][2] + MatML[0][1] * Trans[1][2]) + Trans[1][1] * (MatML[1][0] * Trans[0][2] + MatML[1][1] * Trans[1][2]); // 12
		*(MatM + 6) = Trans[0][2] * (MatML[0][0] * Trans[0][0] + MatML[0][1] * Trans[1][0]) + Trans[1][2] * (MatML[1][0] * Trans[0][0] + MatML[1][1] * Trans[1][0]); // 20
		*(MatM + 7) = Trans[0][2] * (MatML[0][0] * Trans[0][1] + MatML[0][1] * Trans[1][1]) + Trans[1][2] * (MatML[1][0] * Trans[0][1] + MatML[1][1] * Trans[1][1]); // 21
		*(MatM + 8) = Trans[0][2] * (MatML[0][0] * Trans[0][2] + MatML[0][1] * Trans[1][2]) + Trans[1][2] * (MatML[1][0] * Trans[0][2] + MatML[1][1] * Trans[1][2]); // 22
		if (ElemDim == 21 ) {			//		if Elem.dim == 21:
			cc = 0.5 * Emod / (2. + 0. * nu);	//                        # zero poisson's ratio assumed for out of plane shear
			sig[3] = cc * Eps[3];				//			sig_ = array([sig[0], sig[1], 0., cc * Eps[3], cc * Eps[4], sig[2]])
			sig[4] = cc * Eps[4];
			sig[5] = sig[2];
			sig[2] = 0.;
//			MatM_ = array([ [MatM[0, 0], MatM[0, 1], 0., 0., 0., MatM[0, 2]],
//							[MatM[1, 0], MatM[1, 1], 0., 0., 0., MatM[1, 2]],
//							[0.,         0.,         0., 0., 0., 0.],
//							[0.,         0.,         0., cc, 0., 0.],
//							[0.,         0.,         0., 0., cc, 0.],
//							[MatM[2, 0], MatM[2, 1], 0., 0., 0., MatM[2, 2]]] )
//			*(MatM + 0)  = *(MatM + 0);			// 0,0
//			*(MatM + 1)  = *(MatM + 1);			// 0,1
			*(MatM + 11) = *(MatM + 5);			// 1,5 <- 1,2 
			*(MatM + 5)  = *(MatM + 2);			// 0,5 <- 0,2
			*(MatM + 30) = *(MatM + 6); ;		// 5,0 <- 2,0
			*(MatM + 31) = *(MatM + 7);			// 5,1 <- 2,1
			*(MatM + 6)  = *(MatM + 3);			// 1,0
			*(MatM + 7)  = *(MatM + 4);			// 1,1
			*(MatM + 35) = *(MatM + 8);			// 5,5 <- 2,2
			//
			*(MatM + 2)  = 0.;				    // 0,2 
			*(MatM + 3)  = 0.;					// 0,3
			*(MatM + 4)  = 0.;					// 0,4
			*(MatM + 8)  = 0.;					// 1,2 
			*(MatM + 9)  = 0.;					// 1,3
			*(MatM + 10) = 0.;					// 1,4
			//
			for (i = 12; i < 30; i++) *(MatM + i) = 0.;
			*(MatM + 21) = cc;					// 3,3
			*(MatM + 28) = cc;					// 4,4
			*(MatM + 32) = 0.;					// 5,2 
			*(MatM + 33) = 0.;					// 5,3
			*(MatM + 34) = 0.;					// 5,4
		}
	}
	//	# uncracked state
	else {		//	else:
		cc = Emod / (1 - nu * nu);			//			xxx = Emod / (1 - nu * *2)
		if (ElemDim == 2) {						// isotropic linear elastic plane stress
			*(MatM + 0) = cc;					//			MatM = array([[xxx, xxx * nu, 0], [xxx * nu, xxx, 0], [0, 0, 0.5 * xxx * (1 - nu)]] ) # isotropic linear elastic plane stress
			*(MatM + 1) = cc * nu;
			*(MatM + 2) = 0.;
			*(MatM + 3) = cc * nu;
			*(MatM + 4) = cc;
			*(MatM + 5) = 0.;
			*(MatM + 6) = 0.;
			*(MatM + 7) = 0.;
			*(MatM + 8) = 0.5 * cc * (1 - nu);
			sig[0] = *(MatM + 0) * Eps_[0] + *(MatM + 1) * Eps_[1];			//			sig = dot(MatM, Eps_) # stress
			sig[1] = *(MatM + 3) * Eps_[0] + *(MatM + 4) * Eps_[1];
			sig[2] = *(MatM + 8) * Eps_[2];
			PrinCLT_1(sig[0], sig[1], sig[2], pSig);	//		pig, _, _ = PrinCLT_(sig[0], sig[1], sig[2])       # principal stresses, larger value, corr.direction, lower value
		}
		else if ( ElemDim == 21 ) { //	shell	elif Elem.dim == 21:
//			MatM_ = array([[Emod / (1 - nu * *2), nu * Emod / (1 - nu * *2), 0., 0., 0., 0.],
//			[nu * Emod / (1 - nu * *2), Emod / (1 - nu * *2), 0., 0., 0., 0.],
//			[0., 0., 0., 0., 0., 0.],
//			[0., 0., 0., Emod / (2 + 2 * nu), 0., 0.],
//			[0., 0., 0., 0., Emod / (2 + 2 * nu), 0],
//			[0., 0., 0., 0., 0., Emod / (2 + 2 * nu)]] )
			*(MatM + 0) = cc;
			*(MatM + 1) = nu * cc;
			*(MatM + 2) = 0.; *(MatM + 3) = 0.; *(MatM + 4) = 0.; *(MatM + 5) = 0.;
			*(MatM + 6) = nu * cc;
			*(MatM + 7) = cc;
			*(MatM + 8) = 0.; *(MatM + 9) = 0.; *(MatM + 10) = 0.; *(MatM + 11) = 0.;
			for (i = 12; i < 32; i++) *(MatM + i) = 0.;
			*(MatM + 21) = Emod / (2. + 2. * nu);
			*(MatM + 28) = Emod / (2. + 2. * nu);
			*(MatM + 35) = Emod / (2. + 2. * nu);
			//			sig_ = dot(MatM_, Eps)
			//			sig = array([sig_[0], sig_[1], sig_[5]])          #
			sig[0] = *(MatM + 0) * Eps[0] + *(MatM + 1) * Eps[1];
			sig[1] = *(MatM + 6) * Eps[0] + *(MatM + 7) * Eps[1];
			sig[2] = 0.;
			sig[3] = *(MatM + 21) * Eps[3];
			sig[4] = *(MatM + 28) * Eps[4];
			sig[5] = *(MatM + 35) * Eps[5];
			PrinCLT_1(sig[0], sig[1], sig[5], pSig);	//		pig, _, _ = PrinCLT_(sig[0], sig[1], sig[2])       # principal stresses, larger value, corr.direction, lower value
		}
		else { return 50 ; }
		*(ElemStateVarN + 1) = pSig[0];				//		Elem.StateVarN[ipI, 1] = pig                         # larger 1st principal stress for dual cracking control
//		if Elem.dim == 2:
//			return sig, MatM, [Eps_[0], Eps_[1], 0., Eps_[2], sig[0], sig[1], 0., sig[2], ww1, ww2]
//		elif Elem.dim == 21 :
//			return sig_, MatM_, [sig_[0], sig_[1], sig_[2], sig_[3], sig_[4], sig_[5], Eps_[0], Eps_[1], Eps_[2], sig[0], sig[1], sig[2], ww1, ww2]
	}
return 110;
}

int IsoDamUniaxC1(int ElemDim, bool ElemPlSt, double ElemLch, double *ElemStateVar, int dim_St, double *ElemStateVarN, int dim_StN,
	double *Eps, int dim_E, double *sig, int dim_si, double *MatM, int dim_Ma, int LiTy,
	double alpha, double cc1, double cc2, double cc3, int RType, double *EpsR, int dim_EpsR, double kapStrength, double ElemCrBws, double gam2, double kapUlt,
	double edt, double ed, double gd, double nu, double Emod, double *Dps, int dim_D, double eta, double RegPar,
	double *sigR, int dim_sigR, double *CR, int dim_CR, double Dt, double sVsTol, double *DataOut, int dim_DataOut)
{
	double D, kapOld, eps, Dps__[6], kap, dkk, hdI, xxx, sigL, CC, zz, sigV, Veps[6];
	int i;
	/*
	if CalcType == 0: return[], [], []
	*/
	D      = *(ElemStateVar + 0);	//	D = Elem.StateVar[ipI, 0]                        # D of last converged load / time increment
	kapOld = *(ElemStateVar + 1);	//		kapOld = Elem.StateVar[ipI, 1]                   # kappa of last converged load / time increment in loading branch
	if (ElemDim == 1) {				// 		if Elem.dim == 1:
			eps = *Eps;				// eps = Eps_[0]
			Dps__[0] = *Dps; 		//Dps__ = array([Dps_[0], -nu * Dps_[0], -nu * Dps_[0], 0, 0, 0]) # voigt notation
			Dps__[1] = 0.;
			Dps__[2] = 0.;
			Dps__[3] = 0.;
			Dps__[4] = 0.;
			Dps__[5] = 0.;
	}
	else { return 121; }   // 	else: raise NameError("ConFemMaterials::IsoDamageUniax:MatC: elem dim > 1")
	if (eps < 0.0) {
		kap = -eps;
	}
	else {
		kap = alpha * eps;
	}
	dkk = 1.;
	if (kap>kapUlt) { kap = kapUlt; } 	//  if kap>self.kapUlt: kap=self.kapUlt         # Damage should not be zero to avoid numerical singularity. This constrains D to dDestr 
	if ((kap > edt) && (kap > kapOld + kapZero) && (abs(eps) / sqrt(3.)) >= ZeroD) {  //	   if kap>self.edt and kap>=(kapOld+kapZero) and (abs(eps) / sqrt(3.))>ZeroD:  # case loading -- strain deviator condition for compatibility with IsoDamage
		xxx = pow((kap - edt) / ed, gd);
		D   = 1.0 - exp(-xxx);  
		hdI = xxx * gd*exp(-xxx) / (kap - edt) * dkk; // # dD / dkappa
	}
	else {
		hdI = 0;
	}
	sigL = (1. - D)*Emod*eps;
	CC   = (1. - D - hdI * kap)*Emod;
	*(ElemStateVarN + 0) = D; // store damage of actual iteration  // Elem.StateVarN[ipI,0] = D                   # store damage of actual iteration
	*(ElemStateVarN + 1) = kap; // store equivalent damage strain of actual iteration // Elem.StateVarN[ipI,1] = kap                 # 

	if (eta>0.) {  // if self.eta>0.:                             # viscous regularization
		zz = ViscExten3DC1(Dt, eta, Dps__, ElemStateVar, ElemStateVarN, 2, Veps);  // def ViscExten3D(self, Dt, eta, Dps, Elem, ipI, sI):
		sigV = eta * Veps[0]; 
	}
	else {					// else:
		zz = 0.;			// zz = 0.
		for (i = 0; i<6; i++) { Veps[i] = 0.; }  //  Veps = zeros((6),dtype=float) 
		sigV = 0.; 
	}
	sigL = sigL + sigV;
	*(sig)     = sigL;
	*(sig + 1) = 0.;
	*(MatM + 0) = CC + zz; *(MatM + 1) = 0.;
	*(MatM + 2) = 0.;      *(MatM + 3) = 0.;

	return 110;
}

// microplane
void DamFunc1(double kap0V, double alphaV, double betaV, double kapOld, double eta, double *kap, double *dd, double *Dp)
{
	double kap_, dd_, Dp_;
	if (kap0V > kapOld) { kap_ = kap0V; }		//	kap = max(self.kap0V, kapOld)
	else                { kap_ = kapOld; }
	if (eta > kap_)     { kap_ = eta; }			//		if eta > kap: kap = eta
	dd_ = 1. - kap0V / kap_ * (1. + alphaV*(exp(betaV*(kap0V - kap_)) - 1.)); //			dd = 1. - self.kap0V / kap * (1. + self.alphaV*(exp(self.betaV*(self.kap0V - kap)) - 1.))
	if (dd_ > ZeroD) {
		Dp_ = kap0V / (kap_*kap_) * (1 + alphaV * (exp(betaV*(kap0V - kap_)) - 1.)) + kap0V / kap_ * alphaV*betaV*exp(betaV*(kap0V - kap_));
//if dd > ZeroD: Dp = self.kap0V / kap * *2 * (1 + self.alphaV*(exp(self.betaV*(self.kap0V - kap)) - 1.)) + self.kap0V / kap * self.alphaV*self.betaV*exp(self.betaV*(self.kap0V - kap))
	}
	else { Dp_ = 0.;  }						//			else : Dp = 0.
	*kap = kap_; *dd = dd_; *Dp = Dp_;		//				return kap, dd, Dp
}

void DevStiffness( double m00, double m01, double m02, double m11, double m12, double m22, 
					double n02, double n12, double n22, double tbt, double tbt2, double obt, double obt2,
					double *nn, double *DVD)
{
	*(DVD + 0 * 6 + 0) = m00 * n02*tbt2 + 2.*m01* *(nn + 0) * tbt* *(nn + 1) * obt + 2.*m02* *(nn + 0) * obt* *(nn + 2) * tbt + m11 * n12*obt2 + 2.*m12* *(nn + 1) * obt2* *(nn + 2) + m22 * n22*obt2;
	*(DVD + 0 * 6 + 1) = m00 * n02*tbt*obt + m01 * *(nn + 0) * tbt2* *(nn + 1) + m02 * *(nn + 0) * obt* *(nn + 2) * tbt + m01 * *(nn + 0) * obt2* *(nn + 1) + m11 * n12*obt*tbt + m12 * *(nn + 1) * obt2* *(nn + 2) + m02 * *(nn + 0) * obt2* *(nn + 2) + m12 * *(nn + 1) * obt* *(nn + 2) * tbt + m22 * n22*obt2;
	*(DVD + 0 * 6 + 2) = m00 * n02*tbt*obt + m01 * *(nn + 0) * tbt* *(nn + 1) * obt + m02 * *(nn + 0) * tbt2* *(nn + 2) + m01 * *(nn + 0) * obt2* *(nn + 1) + m11 * n12*obt2 + m12 * *(nn + 1) * obt* *(nn + 2) * tbt + m02 * *(nn + 0) * obt2* *(nn + 2) + m12 * *(nn + 1) * obt2* *(nn + 2) + m22 * n22*obt*tbt;
	*(DVD + 0 * 6 + 3) = -0.50*m01* *(nn + 0) * tbt* *(nn + 2) - 0.50*m02* *(nn + 0) * tbt* *(nn + 1) - 0.50*m11* *(nn + 2) *  *(nn + 1) * obt - 0.50*m12*n12*obt - 0.50*m12*n22*obt - 0.50*m22* *(nn + 2) * obt* *(nn + 1);
	*(DVD + 0 * 6 + 4) = -0.50*m00* *(nn + 0) * tbt* *(nn + 2) - 0.50*m02*n02*tbt - 0.50*m01* *(nn + 2) *  *(nn + 1) * obt - 0.50*m12* *(nn + 0) * obt* *(nn + 1) - 0.50*m02*n22*obt - 0.50*m22* *(nn + 0) * obt* *(nn + 2);
	*(DVD + 0 * 6 + 5) = -0.50*m00* *(nn + 0) * tbt* *(nn + 1) - 0.50*m01*n12*obt - 0.50*m02* *(nn + 2) *  *(nn + 1) * obt - 0.50*m01*n02*tbt - 0.50*m11* *(nn + 0) * obt* *(nn + 1) - 0.50*m12* *(nn + 0) * obt* *(nn + 2);
	*(DVD + 1 * 6 + 1) = m00 * n02*obt2 + 2.*m01* *(nn + 0) * tbt* *(nn + 1) * obt + 2.*m02* *(nn + 0) * obt2* *(nn + 2) + m11 * n12*tbt2 + 2.*m12* *(nn + 1) * obt* *(nn + 2) * tbt + m22 * n22*obt2;
	*(DVD + 1 * 6 + 2) = m00 * n02*obt2 + m01 * *(nn + 0) * obt2* *(nn + 1) + m02 * *(nn + 0) * obt* *(nn + 2) * tbt + m01 * *(nn + 0) * tbt* *(nn + 1) * obt + m11 * n12*obt*tbt + m12 * *(nn + 1) * tbt2* *(nn + 2) + m02 * *(nn + 0) * obt2* *(nn + 2) + m12 * *(nn + 1) * obt2* *(nn + 2) + m22 * n22*obt*tbt;
	*(DVD + 1 * 6 + 3) = -0.50*m01* *(nn + 0) * obt* *(nn + 2) - 0.50*m02* *(nn + 0) * obt* *(nn + 1) - 0.50*m11* *(nn + 1) * tbt* *(nn + 2) - 0.50*m12*n12*tbt - 0.50*m12*n22*obt - 0.50*m22* *(nn + 2) * obt* *(nn + 1);
	*(DVD + 1 * 6 + 4) = -0.50*m00* *(nn + 0) * obt* *(nn + 2) - 0.50*m02*n02*obt - 0.50*m01* *(nn + 2) *  *(nn + 1) * tbt - 0.50*m12* *(nn + 1) * tbt* *(nn + 0) - 0.50*m02*n22*obt - 0.50*m22* *(nn + 0) * obt* *(nn + 2);
	*(DVD + 1 * 6 + 5) = -0.50*m00* *(nn + 0) * obt* *(nn + 1) - 0.50*m01*n12*tbt - 0.50*m02* *(nn + 2) *  *(nn + 1) * obt - 0.50*m01*n02*obt - 0.50*m11* *(nn + 1) * tbt* *(nn + 0) - 0.50*m12* *(nn + 0) * obt* *(nn + 2);
	*(DVD + 2 * 6 + 2) = m00 * n02*obt2 + 2.*m01* *(nn + 0) * obt2* *(nn + 1) + 2.*m02* *(nn + 0) * obt* *(nn + 2) * tbt + m11 * n12*obt2 + 2.*m12* *(nn + 1) * obt* *(nn + 2) * tbt + m22 * n22*tbt2;
	*(DVD + 2 * 6 + 3) = -0.50*m01* *(nn + 0) * obt* *(nn + 2) - 0.50*m02* *(nn + 0) * obt* *(nn + 1) - 0.50*m11* *(nn + 2) *  *(nn + 1) * obt - 0.50*m12*n12*obt - 0.50*m12*n22*tbt - 0.50*m22* *(nn + 2) * tbt* *(nn + 1);
	*(DVD + 2 * 6 + 4) = -0.50*m00* *(nn + 0) * obt* *(nn + 2) - 0.50*m02*n02*obt - 0.50*m01* *(nn + 2) *  *(nn + 1) * obt - 0.50*m12* *(nn + 0) * obt* *(nn + 1) - 0.50*m02*n22*tbt - 0.50*m22* *(nn + 0) *  *(nn + 2) * tbt;
	*(DVD + 2 * 6 + 5) = -0.50*m00* *(nn + 0) * obt* *(nn + 1) - 0.50*m01*n12*obt - 0.50*m02* *(nn + 1) *  *(nn + 2) * tbt - 0.50*m01*n02*obt - 0.50*m11* *(nn + 0) * obt* *(nn + 1) - 0.50*m12* *(nn + 0) *  *(nn + 2) * tbt;
	*(DVD + 3 * 6 + 3) = 0.25*m11*n22 + 0.50*m12* *(nn + 1) *  *(nn + 2) + 0.25*m22*n12;
	*(DVD + 3 * 6 + 4) = 0.25*m01*n22 + 0.25*m12* *(nn + 0) *  *(nn + 2) + 0.25*m02* *(nn + 1) *  *(nn + 2) + 0.25*m22* *(nn + 0) *  *(nn + 1);
	*(DVD + 3 * 6 + 5) = 0.25*m01* *(nn + 1) *  *(nn + 2) + 0.25*m02*n12 + 0.25*m11* *(nn + 0) *  *(nn + 2) + 0.25*m12* *(nn + 0) *  *(nn + 1);
	*(DVD + 4 * 6 + 4) = 0.25*m00*n22 + 0.50*m02* *(nn + 0) *  *(nn + 2) + 0.25*m22*n02;
	*(DVD + 4 * 6 + 5) = 0.25*m00* *(nn + 1) *  *(nn + 2) + 0.25*m02* *(nn + 0) *  *(nn + 1) + 0.25*m01* *(nn + 0) *  *(nn + 2) + 0.25*m12*n02;
	*(DVD + 5 * 6 + 5) = 0.25*m00*n12 + 0.50*m01* *(nn + 0) *  *(nn + 1) + 0.25*m11*n02;	DVD[0, 0] = m00 * n02*tbt2 + 2.*m01*nn[0] * tbt*nn[1] * obt + 2.*m02*nn[0] * obt*nn[2] * tbt + m11 * n12*obt2 + 2.*m12*nn[1] * obt2*nn[2] + m22 * n22*obt2;
	*(DVD + 1 * 6 + 0) = *(DVD + 0 * 6 + 1);
	*(DVD + 2 * 6 + 0) = *(DVD + 0 * 6 + 2);
	*(DVD + 2 * 6 + 1) = *(DVD + 1 * 6 + 2);
	*(DVD + 3 * 6 + 0) = *(DVD + 0 * 6 + 3);
	*(DVD + 3 * 6 + 1) = *(DVD + 1 * 6 + 3);
	*(DVD + 3 * 6 + 2) = *(DVD + 2 * 6 + 3);
	*(DVD + 4 * 6 + 0) = *(DVD + 0 * 6 + 4);
	*(DVD + 4 * 6 + 1) = *(DVD + 1 * 6 + 4);
	*(DVD + 4 * 6 + 2) = *(DVD + 2 * 6 + 4);
	*(DVD + 4 * 6 + 3) = *(DVD + 3 * 6 + 4);
	*(DVD + 5 * 6 + 0) = *(DVD + 0 * 6 + 5);
	*(DVD + 5 * 6 + 1) = *(DVD + 1 * 6 + 5);
	*(DVD + 5 * 6 + 2) = *(DVD + 2 * 6 + 5);
	*(DVD + 5 * 6 + 3) = *(DVD + 3 * 6 + 5);
	*(DVD + 5 * 6 + 4) = *(DVD + 4 * 6 + 5);
}
int MicroPlaneC1(int ElemDim, bool ElemPlSt, bool PrinStrains, double ElemLch, double *ElemStateVar, int dim_St, double *ElemStateVarN, int dim_StN,
		double *Eps, int dim_E, double *sig, int dim_si, double *MatM, int dim_Ma, int type, double E_V, double E_D,
		double kV0, double kV1, double kV2, double kap0V, double alphaV, double betaV, int RType, double ElemCrBws, double gam2, double kapUlt, int ElemScaleType,
		double eps_ct, double e0, double ed, double gd, double nu, double Emod, double *Dps, int dim_D, double etaV,
		double Dt, double sVsTol, double *DataOut, int dim_DataOut, int nState, int nInt, int PlStressI, double PlStressL, int iS)
{
	const double obt = 1. / 3., obs = 1. / 6., obn = 1. / 9., tbt = -2. / 3.;
	const double obt2 = obt * obt;
	const double tbt2 = tbt * tbt;
	const double VV[3][3] = {obt, 0, 0, 0, obt, 0, 0, 0, obt};		// VV = array([[obt, 0, 0], [0, obt, 0], [0, 0, obt]]) # projection tensor for volumetric part
	const double VVV[6][6] = {obn,obn,obn,0.,0.,0., obn,obn,obn,0.,0.,0., obn,obn,obn,0.,0.,0., 0.,0.,0.,0.,0.,0., 0.,0.,0.,0.,0.,0., 0.,0.,0.,0.,0.,0.}; //	VVV = zeros((6, 6), dtype = float) VVV[0, 0] = obn VVV[0, 1] = obn VVV[0, 2] = obn VVV[1, 1] = obn VVV[1, 2] = obn VVV[2, 2] = obn 	VVV[1, 0], VVV[2, 0], VVV[2, 1] = VVV[0, 1], VVV[0, 2], VVV[1, 2]
	double ep2, epsT[3][3], Dps__[6], epsVol, dd_iso, I1, J2, kapOld, nn[3], epsDD[3], I1mp, J2mp, eta, eta_, kap, dd, Dp, beta, xx, dkk, sigVV, sigDD[3], zz, Veps[6], sigV[6], cD;
	double sigT[3][3] = { 0.,0.,0., 0.,0.,0., 0.,0.,0.};
	double s00, s01, s02, m00, m11, m22, m01, m02, m12, n02, n12, n22, thet, psi, alph, bet, gam, delt, xyz;
	double D0[3][3], D1[3][3], D2[3][3];
	double DDD[6][6], MMatS[6][6], MMatT[6][6], DV1[6][6];
	double laMax, eVec[3], Eps_[6], sig_[6];
	int ns, ii, i, j, k, rc;
	bool PlStFlag;
	rc = 100;										// normal termination
	ep2 = *(ElemStateVarN + nState*nInt + 4);						// 	ep2 = Elem.StateVarN[ipI, self.nState*self.nInt + 4]  # value of last iteration taken, not from last converged step
	if (ElemDim == 2) {												// 		if Elem.dim == 2:
		if (ElemPlSt) {												// 			if Elem.PlSt :
			epsT[0][0] = *(Eps + 0);								// epsT = array([[Eps[0], 0.5*Eps[2], 0.], [0.5*Eps[2], Eps[1], 0.], [0., 0., ep2]])
			epsT[0][1] = *(Eps + 2)*0.5;
			epsT[0][2] = 0.;
			epsT[1][0] = *(Eps + 2)*0.5;
			epsT[1][1] = *(Eps + 1);
			epsT[1][2] = 0.;
			epsT[2][0] = 0.;
			epsT[2][1] = 0.;
			epsT[2][2] = ep2;
			Dps__[0] = *(Dps + 0);									// Dps__ = array([Dps[0], Dps[1], 0., 0., 0., Dps[2]]) # plane stress--> strain Voigt notation --  used for viscous regularization only with diagonal stiffness
			Dps__[1] = *(Dps + 1);
			Dps__[2] = 0.;
			Dps__[3] = 0.;
			Dps__[4] = 0.;
			Dps__[5] = *(Dps + 2);
		}
		else {  // 	else :
			epsT[0][0] = *(Eps + 0);								//  epsT = array([[Eps[0], 0.5*Eps[2], 0.], [0.5*Eps[2], Eps[1], 0.], [0., 0., 0.]])        #
			epsT[0][1] = *(Eps + 2)*0.5;
			epsT[0][2] = 0.;
			epsT[1][0] = *(Eps + 2)*0.5;
			epsT[1][1] = *(Eps + 1);
			epsT[1][2] = 0.;
			epsT[2][0] = 0.;
			epsT[2][1] = 0.;
			epsT[2][2] = 0.;
			Dps__[0] = *(Dps + 0);									//	Dps__ = array([Dps[0], Dps[1], 0., 0., 0., Dps[2]])  # plane strain--> strain Voigt notation
			Dps__[1] = *(Eps + 1);
			Dps__[2] = 0.;
			Dps__[3] = 0.;
			Dps__[4] = 0.;
			Dps__[5] = *(Dps + 2);
		}
	}
	else if (ElemDim==3) {											// 	elif Elem.dim == 3:
		epsT[0][0] = *(Eps + 0);									//epsT = array([[Eps[0], 0.5*Eps[5], 0.5*Eps[4]], [0.5*Eps[5], Eps[1], 0.5*Eps[3]], [0.5*Eps[4], 0.5*Eps[3], Eps[2]]]) # strain tensor arrangement
		epsT[0][1] = *(Eps + 5)*0.5;
		epsT[0][2] = *(Eps + 4)*0.5;
		epsT[1][0] = *(Eps + 5)*0.5;
		epsT[1][1] = *(Eps + 1);
		epsT[1][2] = *(Eps + 3)*0.5;
		epsT[2][0] = *(Eps + 4)*0.5;
		epsT[2][1] = *(Eps + 3)*0.5;
		epsT[2][2] = *(Eps + 2);
		Dps__[0] = *(Dps + 0);										//Dps__ = array([Dps[0], Dps[1], Dps[2], Dps[3], Dps[4], Dps[5]])                # triaxial strain increment Voigt notation
		Dps__[1] = *(Dps + 1);
		Dps__[2] = *(Dps + 2);
		Dps__[3] = *(Dps + 3);
		Dps__[4] = *(Dps + 4);
		Dps__[5] = *(Dps + 5);
	}
	else if (ElemDim == 21) {										// 	elif Elem.dim == 21:                                                                  # continuum based shell
		epsT[0][0] = *(Eps + 0);									//	epsT = array([[Eps[0], 0.5*Eps[5], 0.5*Eps[4]], [0.5*Eps[5], Eps[1], 0.5*Eps[3]], [0.5*Eps[4], 0.5*Eps[3], ep2]]) # strain tensor arrangement
		epsT[0][1] = *(Eps + 5)*0.5;
		epsT[0][2] = *(Eps + 4)*0.5;
		epsT[1][0] = *(Eps + 5)*0.5;
		epsT[1][1] = *(Eps + 1);
		epsT[1][2] = *(Eps + 3)*0.5;
		epsT[2][0] = *(Eps + 4)*0.5;
		epsT[2][1] = *(Eps + 3)*0.5;
		epsT[2][2] = ep2;
		Dps__[0] = *(Dps + 0);										// Dps__ = array([Dps[0], Dps[1], 0., Dps[3], Dps[4], Dps[5]])                # triaxial strain increment for continuum bases shell Voigt notation
		Dps__[1] = *(Dps + 1);
		Dps__[2] = 0.;
		Dps__[3] = *(Dps + 3);
		Dps__[4] = *(Dps + 4);
		Dps__[5] = *(Dps + 5);
	}
	else { return 200; }											// else: raise NameError("ConFemMaterials::Microplane.Sig: not implemented")
//	epsVol = VV[0][0] * epsT[0][0] + VV[1][1] * epsT[1][1] + VV[2][2] * epsT[2][2];
	ns = nState;													//ns = self.nState
	PlStFlag = 0;													//	PlStFlag = False
	if ((ElemDim == 2 and ElemPlSt) or (ElemDim == 21)) { PlStFlag = 1; } // if (Elem.dim == 2 and Elem.PlSt) or Elem.dim == 21: PlStFlag = True
	// eventual plane stress loop
	for (ii = 0; ii < PlStressI; ii++) {							// 		for ii in range(self.PlStressI) : # for eventual iteration of plane stress
		epsVol = VV[0][0] * epsT[0][0] + VV[1][1] * epsT[1][1] + VV[2][2] * epsT[2][2]; // moved this from above due to updating of epsT[2,2]
		dd_iso = 0.;  // dd_iso, I1, J2 = 0., 0., 0.
		I1 = 0.;
		J2 = 0.;
		for (i = 0; i < 3; i++) { for (j = 0; j < 3; j++) { sigT[i][j] = 0.; } }
		for (i = 0; i < 6; i++) { for (j = 0; j < 6; j++) { MMatS[i][j] = 0.; } }
		for (i = 0; i < 6; i++) { for (j = 0; j < 6; j++) { MMatT[i][j] = 0.; } }
		for (i = 0; i < 6; i++) { for (j = 0; j < 6; j++) { DV1[i][j] = 0.; } }
		// microplane loop
		for (i = 0; i < nInt; i++) {								// 		for i in range(self.nInt) :
			kapOld = *(ElemStateVar + ns * i);						// kapOld = Elem.StateVar[ipI, ns*i]
			nn[0] = I21Points[i][0];								// nn = array([I21Points[i, 0], I21Points[i, 1], I21Points[i, 2]])
			nn[1] = I21Points[i][1];
			nn[2] = I21Points[i][2];
			// V - D - Split projection tensors
			D0[0][0] = nn[0] - nn[0] * obt;  D0[0][1] = 0.5*nn[1]; 		      D0[0][2] = 0.5*nn[2];				// D0 = array([[nn[0] - nn[0] * obt, 0.5*nn[1], 0.5*nn[2]],
			D0[1][0] = 0.5*nn[1];            D0[1][1] = -nn[0] * obt;        D0[1][2] = 0.;					//             [0.5*nn[1], -nn[0] * obt, 0.],
			D0[2][0] = 0.5*nn[2];            D0[2][1] = 0.;                   D0[2][2] = -nn[0] * obt;			//             [0.5*nn[2], 0., -nn[0] * obt]])
			D1[0][0] = -nn[1] * obt;		 D1[0][1] = 0.5*nn[0];            D1[0][2] = 0.;					// D1 = array([[-nn[1] * obt, 0.5*nn[0], 0.],
			D1[1][0] = 0.5*nn[0];			 D1[1][1] = nn[1] - nn[1] * obt;  D1[1][2] = 0.5*nn[2];				// [0.5*nn[0], nn[1] - nn[1] * obt, 0.5*nn[2]],
			D1[2][0] = 0.;					 D1[2][1] = 0.5*nn[2];			  D1[2][2] = -nn[1] * obt;			// [0., 0.5*nn[2], -nn[1] * obt]])
			D2[0][0] = -nn[2] * obt;		 D2[0][1] = 0.;                   D2[0][2] = 0.5*nn[0];				// D2 = array([[-nn[2] * obt, 0., 0.5*nn[0]],
			D2[1][0] = 0.;					 D2[1][1] = -nn[2] * obt;         D2[1][2] = 0.5*nn[1];				// [0., -nn[2] * obt, 0.5*nn[1]],
			D2[2][0] = 0.5*nn[0];		     D2[2][1] = 0.5*nn[1];			  D2[2][2] = nn[2] - nn[2] * obt;	// [0.5*nn[0], 0.5*nn[1], nn[2] - nn[2] * obt]])
			// vector deviator strains
//			epsDD = array([D0[0, 0] * epsT[0, 0] + D0[0, 1] * epsT[0, 1] + D0[0, 2] * epsT[0, 2] + D0[1, 0] * epsT[1, 0] + D0[1, 1] * epsT[1, 1] + D0[1, 2] * epsT[1, 2] + D0[2, 0] * epsT[2, 0] + D0[2, 1] * epsT[2, 1] + D0[2, 2] * epsT[2, 2], \
//				D1[0, 0] * epsT[0, 0] + D1[0, 1] * epsT[0, 1] + D1[0, 2] * epsT[0, 2] + D1[1, 0] * epsT[1, 0] + D1[1, 1] * epsT[1, 1] + D1[1, 2] * epsT[1, 2] + D1[2, 0] * epsT[2, 0] + D1[2, 1] * epsT[2, 1] + D1[2, 2] * epsT[2, 2], \
//				D2[0, 0] * epsT[0, 0] + D2[0, 1] * epsT[0, 1] + D2[0, 2] * epsT[0, 2] + D2[1, 0] * epsT[1, 0] + D2[1, 1] * epsT[1, 1] + D2[1, 2] * epsT[1, 2] + D2[2, 0] * epsT[2, 0] + D2[2, 1] * epsT[2, 1] + D2[2, 2] * epsT[2, 2]])
			epsDD[0] = D0[0][0] * epsT[0][0] + D0[0][1] * epsT[0][1] + D0[0][2] * epsT[0][2] + D0[1][0] * epsT[1][0] + D0[1][1] * epsT[1][1] + D0[1][2] * epsT[1][2] + D0[2][0] * epsT[2][0] + D0[2][1] * epsT[2][1] + D0[2][2] * epsT[2][2];
			epsDD[1] = D1[0][0] * epsT[0][0] + D1[0][1] * epsT[0][1] + D1[0][2] * epsT[0][2] + D1[1][0] * epsT[1][0] + D1[1][1] * epsT[1][1] + D1[1][2] * epsT[1][2] + D1[2][0] * epsT[2][0] + D1[2][1] * epsT[2][1] + D1[2][2] * epsT[2][2];
			epsDD[2] = D2[0][0] * epsT[0][0] + D2[0][1] * epsT[0][1] + D2[0][2] * epsT[0][2] + D2[1][0] * epsT[1][0] + D2[1][1] * epsT[1][1] + D2[1][2] * epsT[1][2] + D2[2][0] * epsT[2][0] + D2[2][1] * epsT[2][1] + D2[2][2] * epsT[2][2];
			// microplane strain invariants, state variable 
			I1mp = 3.*epsVol;
			J2mp = 3. / 2.*(epsDD[0] * epsDD[0] + epsDD[1] * epsDD[1] + epsDD[2] * epsDD[2]);
			eta_ = kV0 * I1mp + sqrt((kV1*I1mp)*(kV1*I1mp) + kV2 * J2mp);
			// damage functons
			if (type == 1) {										//			if   self.type == 1:
				DamFunc1(kap0V, alphaV, betaV, kapOld, eta, &kap, &dd, &Dp);   // kap, dd, Dp = self.DamFunc1(kapOld, eta_)
			}
			else {													// elif self.type == 2 :
				if (RType == 2) {									// if self.RType == 2 : # crack band regularization 
					if (eta_ > eps_ct) {							//	if eta_>self.eps_ct:                                                    # limit strain exceeded
						beta = ElemCrBws;							//beta = Elem.CrBwS
						if (ElemScaleType == 1) {
							xx = eta_ - eps_ct;							//xx = eta_ - self.eps_ct
							eta = beta * xx + (1. - beta) / gam2 * (1. - exp(-gam2 * xx)) + eps_ct;	//eta = beta * xx + (1. - beta) / (self.gam2) * (1. - exp(-self.gam2*xx)) + self.eps_ct
							dkk = beta + (1. - beta) * exp(-gam2 * xx);	//dkk = beta + (1. - beta)             *     exp(-self.gam2*xx)
						}
						else {
							eta = (1. - beta)*eps_ct * log((eta_ - beta * eps_ct) / ((1 - beta)*eps_ct)) + eps_ct; 	//eta = (1 - beta)*self.eps_ct * log((eta_ - beta * self.eps_ct) / ((1 - beta)*self.eps_ct)) + self.eps_ct
							dkk = (1. - beta)*eps_ct / (eta_ - beta * eps_ct);	//dkk = (1 - beta)*self.eps_ct / (eta_ - beta * self.eps_ct)
						}
					}
					else {											// 		else:              # below limit strain
						eta = eta_;
						dkk = 1.;
					}
				}
				else {												//	else:              # no regularization
					eta = eta_;
					dkk = 1.;
				}
				if (eta > kapUlt) { eta = kapUlt; }					//  if eta>self.kapUlt: eta = self.kapUlt         # to have some residual stiffness
				// damage function
				kap = eta;											// 		kap = max(self.e0, kapOld, eta)
				if (kap < kapOld) { kap = kapOld; }
				if (kap <= e0) {									// 				if kap <= self.e0:
					dd = 0.;
					Dp = 0.;
				}
				else {
					dd = 1. - exp(-pow(((kap - e0) / ed), gd));		//dd = 1. - exp(-((kap - self.e0) / self.ed)**self.gd)
					Dp = (1. - dd) * gd / pow(ed, gd) * pow((kap - e0), (gd - 1.));	// Dp = (1. - dd) * self.gd / self.ed**self.gd * (kap - self.e0)**(self.gd - 1.) ???????????????????
					Dp = Dp * dkk;									// Dp = Dp * dkk                                                                 # regularization scaling of derivative dD / dkap tangential stiffness
				}
			}
			// stresses
			sigVV = (1. - dd)*E_V*epsVol; 							// sigVV = (1. - dd)*    self.E_V*epsVol                                             # volumetric stress(scalar)
			sigDD[0] = (1. - dd)*E_D*epsDD[0];						// sigDD = (1. - dd)*dot(self.E_DM, epsDD)                                            # deviatoric stress(vector)
			sigDD[1] = (1. - dd)*E_D*epsDD[1];
			sigDD[2] = (1. - dd)*E_D*epsDD[2];
			// sig = sig + 6.*I21Weights[i] * (sigVV*VV + sigDD[0] * D0 + sigDD[1] * D1 + sigDD[2] * D2)   # stress tensor integration # factor 6. calibrates the weighting factors
			sigT[0][0] = sigT[0][0] + +6.*I21Weights[i] * (sigVV*VV[0][0] + sigDD[0] * D0[0][0] + sigDD[1] * D1[0][0] + sigDD[2] * D2[0][0]);
			sigT[0][1] = sigT[0][1] + +6.*I21Weights[i] * (sigVV*VV[0][1] + sigDD[0] * D0[0][1] + sigDD[1] * D1[0][1] + sigDD[2] * D2[0][1]);
			sigT[0][2] = sigT[0][2] + +6.*I21Weights[i] * (sigVV*VV[0][2] + sigDD[0] * D0[0][2] + sigDD[1] * D1[0][2] + sigDD[2] * D2[0][2]);
			sigT[1][1] = sigT[1][1] + +6.*I21Weights[i] * (sigVV*VV[1][1] + sigDD[0] * D0[1][1] + sigDD[1] * D1[1][1] + sigDD[2] * D2[1][1]);
			sigT[1][2] = sigT[1][2] + +6.*I21Weights[i] * (sigVV*VV[1][2] + sigDD[0] * D0[1][2] + sigDD[1] * D1[1][2] + sigDD[2] * D2[1][2]);
			sigT[2][2] = sigT[2][2] + +6.*I21Weights[i] * (sigVV*VV[2][2] + sigDD[0] * D0[2][2] + sigDD[1] * D1[2][2] + sigDD[2] * D2[2][2]);
			sigT[1][0] = sigT[0][1];
			sigT[2][0] = sigT[0][2];
			sigT[2][1] = sigT[1][2];
			// secant material stiffness
			n02 = nn[0] * nn[0];									//n02 = nn[0]*nn[0]
			n12 = nn[1] * nn[1];									//n12 = nn[1] * *2
			n22 = nn[2] * nn[2];									//n22 = nn[2] * *2
			s00 = 1.; s01 = 1.; s02 = 1.;							//s00, s01, s02 = 1., 1., 1.
			m00 = 1.; m11 = 1.; m22 = 1.;							//m00, m11, m22 = 1., 1., 1.
			m01 = 0.; m02 = 0.; m12 = 0.;							//m01, m02, m12 = 0., 0., 0.
			DevStiffness(m00, m01, m02, m11, m12, m22, n02, n12, n22, tbt, tbt2, obt, obt2,	nn, &DDD[0][0]); 			//DDD = DevStiffness()
			// MMatS = MMatS + 6.*I21Weights[i] * ((1. - dd)*self.E_V*VVV + (1. - dd)*self.E_D*DDD) # presumably symmetric, i.e.only upper right part builded
			for (j = 0; j < 6; j++) {
				for (k = 0; k < 6; k++) {
					MMatS[j][k]= MMatS[j][k]+6.*I21Weights[i] * ((1. - dd)*E_V*VVV[j][k] + (1. - dd)*E_D*DDD[j][k]);
				}
			}
			// corrector for tangential material stiffness
			if ((kap > kapOld) && (dd > ZeroD)) {       			// if kap>kapOld and dd>ZeroD:
				s00 = epsDD[0]; s01 = epsDD[1]; s02 = epsDD[2];		//s00, s01, s02 = epsDD[0], epsDD[1], epsDD[2]
				m00 = s00 * s00; m11 = s01 * s01; m22 = s02 * s02;	//m00, m11, m22 = s00 * *2, s01**2, s02**2
				m01 = s00 * s01; m02 = s00 * s02; m12 = s01 * s02;	//m01, m02, m12 = s00 * s01, s00*s02, s01*s02
				DevStiffness(m00, m01, m02, m11, m12, m22, n02, n12, n22, tbt, tbt2, obt, obt2, nn, &DDD[0][0]); 			//DDD = DevStiffness()
				DV1[0][0] = -obt * s00*nn[0] * tbt - obt2 * s01*nn[1] - obt2 * s02*nn[2];
				DV1[0][1] = DV1[0][0];
				DV1[0][2] = DV1[0][0];
				DV1[1][0] = -obt2 * s00*nn[0] - obt * s01*nn[1] * tbt - obt2 * s02*nn[2];
				DV1[1][1] = DV1[1][0];
				DV1[1][2] = DV1[1][0];
				DV1[2][0] = -obt2 * s00*nn[0] - obt2 * s01*nn[1] - obt * s02*nn[2] * tbt;
				DV1[2][1] = DV1[2][0];
				DV1[2][2] = DV1[2][0];
				DV1[3][0] = obs * s01*nn[2] + obs * s02*nn[1];
				DV1[3][1] = DV1[3][0];
				DV1[3][2] = DV1[3][0];
				DV1[4][0] = obs * s00*nn[2] + obs * s02*nn[0];
				DV1[4][1] = DV1[4][0];
				DV1[4][2] = DV1[4][0];
				DV1[5][0] = obs * s00*nn[1] + obs * s01*nn[0];
				DV1[5][1] = DV1[5][0];
				DV1[5][2] = DV1[5][0];
				xx = sqrt(kV1*kV1 * I1mp*I1mp + kV2 * J2mp);		//xx = sqrt(self.kV1**2 * I1mp**2 + self.kV2 * J2mp)
				thet = 9.*kV1*kV1 / xx;								//thet = 9.*self.kV1**2 / xx
				psi  = 1.5*kV2 / xx;								//psi = 1.5*self.kV2 / xx
				alph = (3.*kV0 + thet * epsVol)*E_V*epsVol;			//alph = (3.*self.kV0 + thet * epsVol)*self.E_V*epsVol
				bet  = (3.*kV0 + thet * epsVol)*E_D;				//bet = (3.*self.kV0 + thet * epsVol)*self.E_D
				gam  = psi * E_V*epsVol;							//gam = psi * self.E_V*epsVol
				delt = psi * E_D;									//delt = psi * self.E_D
				for (j = 0; j < 6; j++) {
					for (k = 0; k < 6; k++) {
						MMatT[j][k] = MMatT[j][k] - 6.*I21Weights[i] * Dp*(alph * VVV[j][k] + bet * DV1[j][k] + gam * DV1[k][j] + delt * DDD[j][k]); //MMatT = MMatT - 6.*I21Weights[i] * Dp*(alph * VVV + bet * DV1 + gam * transpose(DV1) + delt * DDD)
					}
				}
			}
			dd_iso = dd_iso+2.*I21Weights[i] * dd;					//dd_iso += 2.*I21Weights[i] * dd
			I1 = I1 + 2.*I21Weights[i] * I1mp;						//I1 += 2.*I21Weights[i] * I1mp
			J2 = J2 + 2.*I21Weights[i] * J2mp;						//J2 += 2.*I21Weights[i] * J2mp
			// microplane state variable updates
//
//			int ll_ = 7;
//			*(DataOut + i *ll_ + 0) = i;
//			*(DataOut + i *ll_ + 1) = dd;
//			*(DataOut + i *ll_ + 2) = sigVV;
//			*(DataOut + i *ll_ + 3) = sigDD[0];
//			*(DataOut + i *ll_ + 4) = sigDD[1];
//			*(DataOut + i *ll_ + 5) = sigDD[2];
//			*(DataOut + i * ll_ + 6) = sigT[0][0];
//
			*(ElemStateVarN + ns * i) = kap;						//Elem.StateVarN[ipI, ns*i] = kap
		}
		ep2   = -(MMatS[2][0] * epsT[0][0] + MMatS[2][1] * epsT[1][1] + 2.0*MMatS[2][5] * epsT[0][1]) / MMatS[2][2];
		if (abs(epsT[2][2]) > ZeroD*1.e-3)  { xyz = abs((epsT[2][2] - ep2) / epsT[2][2]); }  //if abs(epsT[2][2])>ZeroD*1.e-3: xyz = abs((epsT[2, 2] - ep2) / epsT[2, 2])
		else								{ xyz = 0.; }			//else:                          xyz = 0.
		epsT[2][2] = ep2; 											//		if PlStFlag : epsT[2, 2] = ep2
		if ((!PlStFlag) || (xyz < PlStressL)) { break; }			//	if not PlStFlag or xyz<self.PlStressL : break
//		epsT[2][2] = ep2; 											//		if PlStFlag : epsT[2, 2] = ep2
	}
	if (PlStFlag && ii >= PlStressI) { rc = 110; } //if ii >= self.PlStressI - 1 : raise NameError("ConFemMaterials::Microplane.Sig: plane stress iteration", ii, sig[2, 2])
	*(ElemStateVarN + nState * nInt)     = dd_iso;					//Elem.StateVarN[ipI, self.nState*self.nInt] = dd_iso
	*(ElemStateVarN + nState * nInt + 1) = I1;						//Elem.StateVarN[ipI, self.nState*self.nInt + 1] = I1
	*(ElemStateVarN + nState * nInt + 2) = J2;						//Elem.StateVarN[ipI, self.nState*self.nInt + 2] = J2
	*(ElemStateVarN + nState * nInt + 4) = ep2;						//Elem.StateVarN[ipI, self.nState*self.nInt + 4] = ep2
	// Rankine data
	if (RType == 3) {
		Eps_[0] = epsT[0][0];
		Eps_[1] = epsT[1][1];
		Eps_[2] = epsT[2][2];
		Eps_[3] = epsT[1][2];
		Eps_[4] = epsT[0][2];
		Eps_[5] = epsT[0][1];
		sig_[0] = sigT[0][0];
		sig_[1] = sigT[1][1];
		sig_[2] = sigT[2][2];
		sig_[3] = sigT[1][2];
		sig_[4] = sigT[0][2];
		sig_[5] = sigT[0][1];
		rc = Rankine(PrinStrains, ElemDim, Eps_, sig_, eVec, &laMax);  // laMax may also be stresses, check whether this confused with LaMax_
		if (rc > 110) { return rc; }
	}
	// viscous contribution
	if (etaV>0.) {													// if self.eta>0.:                             # viscous regularization
		zz = ViscExten3DC1(Dt, etaV, Dps__, ElemStateVar, ElemStateVarN, iS, Veps);  // def ViscExten3D(self, Dt, etaV, Dps, Elem, ipI, sI):
		for (i = 0; i<6; i++) { sigV[i]     = etaV * Veps[i]; }
		for (i = 0; i<6; i++) { MMatS[i][i] = MMatS[i][i] + zz; }
	}
	else {															// else:
		zz = 0.;													// zz = 0.
		for (i = 0; i<6; i++) { sigV[i] = 0.; }
	}
	// move Rankine data to state variables with viscous contributions for further processing regarding SDA
	if (RType == 3) {
		for (i = 0; i<6; i++) { sig_[i] = sig_[i] + sigV[i]; }
		rc = RankineUpd(sig_, eVec, ElemStateVarN, &laMax, iS);
	}
	// collect
	for (j = 0; j < 6; j++) {										//		MatM = MMatS + MMatT
		for (k = 0; k < 6; k++) {
			MMatS[j][k] = MMatS[j][k] + MMatT[j][k];
		}
	}
	if (ElemDim == 2) {												//		if Elem.dim == 2 : # 2D
		if (ElemPlSt) {												// 			if Elem.PlSt:                                                                   # plane stress
			cD = 1. / MMatS[2][2];									// cD = 1. / MatM[2, 2]
			*(MatM + 0) = MMatS[0][0] - MMatS[0][2]*MMatS[2][0]*cD;  *(MatM + 1) = MMatS[0][1] - MMatS[0][2]*MMatS[2][1]*cD;  *(MatM + 2) = MMatS[0][5] - MMatS[0][2]*MMatS[2][5]*cD;  //MatM_ = array([[MatM[0, 0] - MatM[0, 2] * MatM[2, 0] * cD, MatM[0, 1] - MatM[0, 2] * MatM[2, 1] * cD, MatM[0, 5] - MatM[0, 2] * MatM[2, 5] * cD],
			*(MatM + 3) = MMatS[1][0] - MMatS[1][2]*MMatS[2][0]*cD;  *(MatM + 4) = MMatS[1][1] - MMatS[1][2]*MMatS[2][1]*cD;  *(MatM + 5) = MMatS[1][5] - MMatS[1][2]*MMatS[2][5]*cD;  //[MatM[1, 0] - MatM[1, 2] * MatM[2, 0] * cD, MatM[1, 1] - MatM[1, 2] * MatM[2, 1] * cD, MatM[1, 5] - MatM[1, 2] * MatM[2, 5] * cD],
			*(MatM + 6) = MMatS[5][0] - MMatS[5][2]*MMatS[2][0]*cD;  *(MatM + 7) = MMatS[5][1] - MMatS[5][2]*MMatS[2][1]*cD;  *(MatM + 8) = MMatS[5][5] - MMatS[5][2]*MMatS[2][5]*cD;  //[MatM[5, 0] - MatM[5, 2] * MatM[2, 0] * cD, MatM[5, 1] - MatM[5, 2] * MatM[2, 1] * cD, MatM[5, 5] - MatM[5, 2] * MatM[2, 5] * cD]])
		}
		else {														// 	else:                # plane strain
			*(MatM + 0) = MMatS[0][0];  *(MatM + 1) = MMatS[0][1];  *(MatM + 2) = MMatS[0][5]; //MatM_=array([[MatM[0,0],MatM[0,1],MatM[0,5]],[MatM[1,0],MatM[1,1],MatM[1,5]],[MatM[5,0],MatM[5,1],MatM[5,5]]])
			*(MatM + 3) = MMatS[1][0];  *(MatM + 4) = MMatS[1][1];  *(MatM + 5) = MMatS[1][5];
			*(MatM + 6) = MMatS[5][0];  *(MatM + 7) = MMatS[5][1];  *(MatM + 8) = MMatS[5][5];
		}
		*(sig + 0) = sigT[0][0] + sigV[0];							// return array([sig[0, 0], sig[1, 1], sig[0, 1]]), MatM_, [Eps[0], Eps[1], Eps[2], sig[0, 0], sig[1, 1], sig[0, 1]]
		*(sig + 1) = sigT[1][1] + sigV[1];
		*(sig + 2) = sigT[0][1] + sigV[5];							// voigt notation broken down to 2D
	}
	else if (ElemDim == 3) {										// elif Elem.dim == 3:                                                                   # 3D
		*(MatM + 0)  = MMatS[0][0]; *(MatM + 1)  = MMatS[0][1]; *(MatM + 2)  = MMatS[0][2]; *(MatM + 3)  = MMatS[0][3]; *(MatM + 4)  = MMatS[0][4]; *(MatM + 5)  = MMatS[0][5];
		*(MatM + 6)  = MMatS[1][0]; *(MatM + 7)  = MMatS[1][1]; *(MatM + 8)  = MMatS[1][2]; *(MatM + 9)  = MMatS[1][3]; *(MatM + 10) = MMatS[1][4]; *(MatM + 11) = MMatS[1][5];//
		*(MatM + 12) = MMatS[2][0]; *(MatM + 13) = MMatS[2][1]; *(MatM + 14) = MMatS[2][2]; *(MatM + 15) = MMatS[2][3]; *(MatM + 16) = MMatS[2][4]; *(MatM + 17) = MMatS[2][5];
		*(MatM + 18) = MMatS[3][0]; *(MatM + 19) = MMatS[3][1]; *(MatM + 20) = MMatS[3][2]; *(MatM + 21) = MMatS[3][3]; *(MatM + 22) = MMatS[3][4]; *(MatM + 23) = MMatS[3][5];
		*(MatM + 24) = MMatS[4][0]; *(MatM + 25) = MMatS[4][1]; *(MatM + 26) = MMatS[4][2]; *(MatM + 27) = MMatS[4][3]; *(MatM + 28) = MMatS[4][4]; *(MatM + 29) = MMatS[4][5];
		*(MatM + 30) = MMatS[5][0]; *(MatM + 31) = MMatS[5][1]; *(MatM + 32) = MMatS[5][2]; *(MatM + 33) = MMatS[5][3]; *(MatM + 34) = MMatS[5][4]; *(MatM + 35) = MMatS[5][5];
		*(sig + 0) = sigT[0][0] + sigV[0];							// return array([sig[0,0],sig[1,1],sig[2,2],sig[1,2],sig[0,2],sig[0,1]]),MatM,[sig[0,0],sig[1,1],sig[2,2],sig[1,2],sig[0,2],sig[0,1]]
		*(sig + 1) = sigT[1][1] + sigV[1];
		*(sig + 2) = sigT[2][2] + sigV[2];
		*(sig + 3) = sigT[1][2] + sigV[3];
		*(sig + 4) = sigT[0][2] + sigV[4];
		*(sig + 5) = sigT[0][1] + sigV[5];
	}
	else if (ElemDim == 21) {										// 	elif Elem.dim == 21:                                                                  # Continuum based shell
		if (abs(MMatS[2][2]) < ZeroD) { return 220; }
		cD = 1. / MMatS[2][2];										// 		cD = 1. / MatM[2, 2]
		*(MatM + 0)  = MMatS[0][0] - MMatS[0][2]*MMatS[2][0]*cD; *(MatM + 1)  = MMatS[0][1] - MMatS[0][2]*MMatS[2][1]*cD; *(MatM + 2)  = 0.; *(MatM + 3)  = MMatS[0][3] - MMatS[0][2]*MMatS[2][3]*cD; *(MatM + 4)  = MMatS[0][4] - MMatS[0][2]*MMatS[2][4]*cD; *(MatM + 5)  = MMatS[0][5] - MMatS[0][2]*MMatS[2][5]*cD;
		*(MatM + 6)  = MMatS[1][0] - MMatS[1][2]*MMatS[2][0]*cD; *(MatM + 7)  = MMatS[1][1] - MMatS[1][2]*MMatS[2][1]*cD; *(MatM + 8)  = 0.; *(MatM + 9)  = MMatS[1][3] - MMatS[1][2]*MMatS[2][3]*cD; *(MatM + 10) = MMatS[1][4] - MMatS[1][2]*MMatS[2][4]*cD; *(MatM + 11) = MMatS[1][5] - MMatS[1][2]*MMatS[2][5]*cD;
		*(MatM + 12) = 0.;                                       *(MatM + 13) = 0.;                                       *(MatM + 14) = 0.; *(MatM + 15) = 0.;                                       *(MatM + 16) = 0.;                                       *(MatM + 17) = 0.;
		*(MatM + 18) = MMatS[3][0] - MMatS[3][2]*MMatS[2][0]*cD; *(MatM + 19) = MMatS[3][1] - MMatS[3][2]*MMatS[2][1]*cD; *(MatM + 20) = 0.; *(MatM + 21) = MMatS[3][3] - MMatS[3][2]*MMatS[2][3]*cD; *(MatM + 22) = MMatS[3][4] - MMatS[3][2]*MMatS[2][4]*cD; *(MatM + 23) = MMatS[3][5] - MMatS[3][2]*MMatS[2][5]*cD;
		*(MatM + 24) = MMatS[4][0] - MMatS[4][2]*MMatS[2][0]*cD; *(MatM + 25) = MMatS[4][1] - MMatS[4][2]*MMatS[2][1]*cD; *(MatM + 26) = 0.; *(MatM + 27) = MMatS[4][3] - MMatS[4][2]*MMatS[2][3]*cD; *(MatM + 28) = MMatS[4][4] - MMatS[4][2]*MMatS[2][4]*cD; *(MatM + 29) = MMatS[4][5] - MMatS[4][2]*MMatS[2][5]*cD;
		*(MatM + 30) = MMatS[5][0] - MMatS[5][2]*MMatS[2][0]*cD; *(MatM + 31) = MMatS[5][1] - MMatS[5][2]*MMatS[2][1]*cD; *(MatM + 32) = 0.; *(MatM + 33) = MMatS[5][3] - MMatS[5][2]*MMatS[2][3]*cD; *(MatM + 34) = MMatS[5][4] - MMatS[5][2]*MMatS[2][4]*cD; *(MatM + 35) = MMatS[5][5] - MMatS[5][2]*MMatS[2][5]*cD;
		*(sig + 0) = sigT[0][0] + sigV[0];							// return array([sig[0,0],sig[1,1],sig[2,2],sig[1,2],sig[0,2],sig[0,1]]),MatM,[sig[0,0],sig[1,1],sig[2,2],sig[1,2],sig[0,2],sig[0,1]]
		*(sig + 1) = sigT[1][1] + sigV[1];
		*(sig + 2) = sigT[2][2] + sigV[2];
		*(sig + 3) = sigT[1][2] + sigV[3];
		*(sig + 4) = sigT[0][2] + sigV[4];
		*(sig + 5) = sigT[0][1] + sigV[5];
	}
	else { return 225; }											// else: raise NameError("ConFemMaterials::Microplane.Sig: not implemented for this element type")
	return rc;
}

int MicroPlaneC2(int ElemDim, bool ElemPlSt, bool PrinStrains, double ElemLch, double *ElemStateVar, int dim_St, double *ElemStateVarN, int dim_StN,
	double *Eps, int dim_E, double *sig, int dim_si, double *MatM, int dim_Ma, int type, double E_V, double E_D,
	double kV0, double kV1, double kV2, double kap0V, double alphaV, double betaV, int RType, double ElemCrBws, double gam2, double kapUlt, int ElemScaleType,
	double eps_ct, double e0, double ed, double gd, double nu, double Emod, double *Dps, int dim_D, double etaV,
	double Dt, double sVsTol, double *DataOut, int dim_DataOut, int nState, int nInt, int PlStressI, double PlStressL, int iS)
{
	int rc, ns, ii, i, j, k;
	const double ob3 = 1. / 3., ob6 = 1. / 6., ob9 = 1. / 9., tb3 = 2. / 3., ob2 = 0.5; 
	const double fb9 = 4. / 9., tb9 = 2. / 9., ob4 = 1. / 4.; 
	const double VVV[6][6] = { ob9,ob9,ob9,0.,0.,0., ob9,ob9,ob9,0.,0.,0., ob9,ob9,ob9,0.,0.,0., 0.,0.,0.,0.,0.,0., 0.,0.,0.,0.,0.,0., 0.,0.,0.,0.,0.,0. }; //	VVV = VolStiffness() VVV = zeros((6, 6), dtype = float) VVV[0, 0] = obn VVV[0, 1] = obn VVV[0, 2] = obn VVV[1, 1] = obn VVV[1, 2] = obn VVV[2, 2] = obn 	VVV[1, 0], VVV[2, 0], VVV[2, 1] = VVV[0, 1], VVV[0, 2], VVV[1, 2]
	double ep2, dd_iso, I1, J2, kapOld, eps[6], nn[3], MMatS[6][6], MMatT[6][6], epsVol, D0[6], D1[6], D2[6], epsDD[3], I1mp, J2mp, eta_, beta, xx, eta, dkk, kap, dd, Dp, sigVV, sigDD[3], Dps__[6], sigV[6];
	double DevStiff[6][6], CCStiff[6][6], VCStiff[6][6], q, Q, A, B, sigV0, c[6], xyz, sig_[6], cD, zz, Veps[6];
	bool PlStFlag;
	rc = 100;										// normal termination
	ep2 = *(ElemStateVarN + nState * nInt + 4);						// 	ep2 = Elem.StateVarN[ipI, self.nState*self.nInt + 4]  # value of last iteration taken, not from last converged step
	if (ElemDim == 2) {												// 		if Elem.dim == 2:
		if (ElemPlSt) {												// 			if Elem.PlSt :
			eps[0] = *(Eps + 0);									//			if Elem.PlSt:  Eps = array([Eps_[0], Eps_[1], ep2, 0., 0., Eps_[2]])
			eps[1] = *(Eps + 1);
			eps[2] = ep2;
			eps[3] = 0.;
			eps[4] = 0.;
			eps[5] = *(Eps + 2);
			Dps__[0] = *(Dps + 0);									// Dps__ = array([Dps[0], Dps[1], 0., 0., 0., Dps[2]]) # plane stress--> strain Voigt notation --  used for viscous regularization only with diagonal stiffness
			Dps__[1] = *(Dps + 1);
			Dps__[2] = 0.;
			Dps__[3] = 0.;
			Dps__[4] = 0.;
			Dps__[5] = *(Dps + 2);
		}
		else {  // 	else :	else:          Eps = array([Eps_[0], Eps_[1], 0., 0., 0., Eps_[2]])
			eps[0] = *(Eps + 0);
			eps[1] = *(Eps + 1);
			eps[2] = 0.;
			eps[3] = 0.;
			eps[4] = 0.;
			eps[5] = *(Eps + 2);
			Dps__[0] = *(Dps + 0);									//	Dps__ = array([Dps[0], Dps[1], 0., 0., 0., Dps[2]])  # plane strain--> strain Voigt notation
			Dps__[1] = *(Eps + 1);
			Dps__[2] = 0.;
			Dps__[3] = 0.;
			Dps__[4] = 0.;
			Dps__[5] = *(Dps + 2);
		}
	}
	else if (ElemDim == 4) {
		eps[0] = *(Eps + 0);
		eps[1] = *(Eps + 1);
		eps[2] = *(Eps + 2);
		eps[3] = 0.;
		eps[4] = 0.;
		eps[5] = *(Eps + 3);
		Dps__[0] = *(Dps + 0);
		Dps__[1] = *(Dps + 1);
		Dps__[2] = *(Dps + 2);
		Dps__[3] = 0.;
		Dps__[4] = 0.;
		Dps__[5] = *(Dps + 3);
	}
	else if (ElemDim == 3) {											// 	elif Elem.dim == 3:		elif Elem.dim == 3:  Eps = array([Eps_[0], Eps_[1], Eps_[2], Eps_[3], Eps_[4], Eps_[5]])
		eps[0] = *(Eps + 0);
		eps[1] = *(Eps + 1);
		eps[2] = *(Eps + 2);
		eps[3] = *(Eps + 3);
		eps[4] = *(Eps + 4);
		eps[5] = *(Eps + 5);
		Dps__[0] = *(Dps + 0);										//Dps__ = array([Dps[0], Dps[1], Dps[2], Dps[3], Dps[4], Dps[5]])                # triaxial strain increment Voigt notation
		Dps__[1] = *(Dps + 1);
		Dps__[2] = *(Dps + 2);
		Dps__[3] = *(Dps + 3);
		Dps__[4] = *(Dps + 4);
		Dps__[5] = *(Dps + 5);
	}
	else if (ElemDim == 21) {										// 	elif Elem.dim == 21: Eps = array([Eps_[0], Eps_[1], 0., Eps_[3], Eps_[4], Eps_[5]])                      # continuum based shell
		eps[0] = *(Eps + 0);
		eps[1] = *(Eps + 1);
		eps[2] = ep2;
		eps[3] = *(Eps + 3);
		eps[4] = *(Eps + 4);
		eps[5] = *(Eps + 5);
		Dps__[0] = *(Dps + 0);										// Dps__ = array([Dps[0], Dps[1], 0., Dps[3], Dps[4], Dps[5]])                # triaxial strain increment for continuum bases shell Voigt notation
		Dps__[1] = *(Dps + 1);
		Dps__[2] = 0.;
		Dps__[3] = *(Dps + 3);
		Dps__[4] = *(Dps + 4);
		Dps__[5] = *(Dps + 5);
	}
	else { return 200; }											// else: raise NameError("ConFemMaterials::Microplane.Sig: not implemented")
//	epsVol = ob3* (eps[0] +eps[1] + eps[2]);						//	epsVol = ob3 * (Eps[0] + Eps[1] + Eps[2])
	ns = nState;													//ns = self.nState
	PlStFlag = 0;													//	PlStFlag = False
	if ((ElemDim == 2 and ElemPlSt) or (ElemDim == 21)) { PlStFlag = 1; } // if (Elem.dim == 2 and Elem.PlSt) or Elem.dim == 21: PlStFlag = True
																		  // eventual plane stress loop
	for (ii = 0; ii < PlStressI; ii++) {							// 		for ii in range(self.PlStressI) : # for eventual iteration of plane stress
		epsVol = ob3 * (eps[0] + eps[1] + eps[2]);						//	epsVol = ob3 * (Eps[0] + Eps[1] + Eps[2])
		dd_iso = 0.;												// dd_iso, I1, J2 = 0., 0., 0.
		I1 = 0.;
		J2 = 0.;
		for (i = 0; i < 6; i++)                           { sig_[i] = 0.; } 
		for (i = 0; i < 6; i++) { for (j = 0; j < 6; j++) { MMatS[i][j] = 0.; } }
		for (i = 0; i < 6; i++) { for (j = 0; j < 6; j++) { MMatT[i][j] = 0.; } }
		// microplane loop
		for (i = 0; i < nInt; i++) {								// 		for i in range(self.nInt) :
			kapOld = *(ElemStateVar + ns * i);						// kapOld = Elem.StateVar[ipI, ns*i]
			nn[0] = I21Points[i][0];								// nn = array([I21Points[i, 0], I21Points[i, 1], I21Points[i, 2]])
			nn[1] = I21Points[i][1];
			nn[2] = I21Points[i][2];
			// V - D - Split projection tensors
			D0[0] =  tb3 * nn[0];									// 			D0 = array([tb3*nn[0], -ob3 * nn[0], -ob3 * nn[0], 0, ob2*nn[2], ob2*nn[1]])
			D0[1] = -ob3 * nn[0];
			D0[2] = -ob3 * nn[0];
			D0[3] =  0.;
			D0[4] =  ob2 * nn[2];
			D0[5] =  ob2 * nn[1];
			D1[0] = -ob3 * nn[1];									// 			D1 = array([-ob3 * nn[1], tb3*nn[1], -ob3 * nn[1], ob2*nn[2], 0, ob2*nn[0]])
			D1[1] =  tb3 * nn[1];
			D1[2] = -ob3 * nn[1];
			D1[3] =  ob2 * nn[2];
			D1[4] =  0.;
			D1[5] =  ob2 * nn[0];
			D2[0] = -ob3 * nn[2];									// 		D2 = array([-ob3 * nn[2], -ob3 * nn[2], tb3*nn[2], ob2*nn[1], ob2*nn[0], 0])
			D2[1] = -ob3 * nn[2];
			D2[2] =  tb3 * nn[2];
			D2[3] =  ob2 * nn[1];
			D2[4] =  ob2 * nn[0];
			D2[5] = 0.;
			// vector deviator strains
			epsDD[0] = D0[0] * eps[0] + D0[1] * eps[1] + D0[2] * eps[2] + D0[3] * eps[3] + D0[4] * eps[4] + D0[5] * eps[5]; // 	epsDD = array([dot(D0, Eps), dot(D1, Eps), dot(D2, Eps)])
			epsDD[1] = D1[0] * eps[0] + D1[1] * eps[1] + D1[2] * eps[2] + D1[3] * eps[3] + D1[4] * eps[4] + D1[5] * eps[5];
			epsDD[2] = D2[0] * eps[0] + D2[1] * eps[1] + D2[2] * eps[2] + D2[3] * eps[3] + D2[4] * eps[4] + D2[5] * eps[5];
			// microplane strain invariants, state variable
			I1mp = 3.*epsVol;										// I1mp = 3.*epsVol
			J2mp = 1.5* (epsDD[0] * epsDD[0] + epsDD[1] * epsDD[1] + epsDD[2] * epsDD[2]);	//J2mp = 3. / 2.*dot(epsDD, epsDD)
			eta_ = kV0 * I1mp + sqrt((kV1*I1mp)*(kV1*I1mp) + kV2 * J2mp);		//eta_ = self.kV0*I1mp + sqrt((self.kV1*I1mp)**2 + self.kV2*J2mp)
			// equivalent strain
			if (RType == 2) {									// if self.RType == 2 : # crack band regularization 
				if (eta_ > eps_ct) {							//	if eta_>self.eps_ct:                                                    # limit strain exceeded
					beta = ElemCrBws;							//beta = Elem.CrBwS
					if (ElemScaleType == 1) {
						xx = eta_ - eps_ct;							//xx = eta_ - self.eps_ct
						eta = beta * xx + (1. - beta) / gam2 * (1. - exp(-gam2 * xx)) + eps_ct;	//eta = beta * xx + (1. - beta) / (self.gam2) * (1. - exp(-self.gam2*xx)) + self.eps_ct
						dkk = beta + (1. - beta) * exp(-gam2 * xx);	//dkk = beta + (1. - beta)             *     exp(-self.gam2*xx)
					}
					else {
						eta = (1. - beta)*eps_ct * log((eta_ - beta * eps_ct) / ((1 - beta)*eps_ct)) + eps_ct; 	//eta = (1 - beta)*self.eps_ct * log((eta_ - beta * self.eps_ct) / ((1 - beta)*self.eps_ct)) + self.eps_ct
						dkk = (1. - beta)*eps_ct / (eta_ - beta * eps_ct);	//dkk = (1 - beta)*self.eps_ct / (eta_ - beta * self.eps_ct)
					}
				}
				else {											// 		else:              # below limit strain
					eta = eta_;
					dkk = 1.;
				}
			}
			else {												//	else:              # no regularization
				eta = eta_;
				dkk = 1.;
			}
			if (eta > kapUlt) { eta = kapUlt; }					//  if eta>self.kapUlt: eta = self.kapUlt         # to have some residual stiffness
			// damage function
			kap = eta;											// 		kap = max(self.e0, kapOld, eta)
			if (kap < kapOld) { kap = kapOld; }
			if (kap <= e0) {									// 				if kap <= self.e0:
				dd = 0.;
				Dp = 0.;
			}
			else {
				dd = 1. - exp(-pow(((kap - e0) / ed), gd));		//dd = 1. - exp(-((kap - self.e0) / self.ed)**self.gd)
				Dp = (1. - dd) * gd / pow(ed, gd) * pow((kap - e0), (gd - 1.));	// Dp = (1. - dd) * self.gd / self.ed**self.gd * (kap - self.e0)**(self.gd - 1.) ???????????????????
				Dp = Dp * dkk;									// Dp = Dp * dkk                                                                 # regularization scaling of derivative dD / dkap tangential stiffness
			}
			// stresses 
			sigVV = (1. - dd)*E_V*epsVol;						//sigVV = (1. - dd)*    self.E_V*epsVol                     # volumetric stress(scalar)
			sigDD[0] = (1. - dd)*E_D*epsDD[0];					//sigDD = (1. - dd)*dot(self.E_DM, epsDD)                    # deviatoric stress(vector)
			sigDD[1] = (1. - dd)*E_D*epsDD[1];
			sigDD[2] = (1. - dd)*E_D*epsDD[2];
			sigV[0] = ob3 * sigVV;								//sigV = array([ob3*sigVV, ob3*sigVV, ob3*sigVV, 0., 0., 0.])
			sigV[1] = ob3 * sigVV;
			sigV[2] = ob3 * sigVV;
			sigV[3] = 0.;
			sigV[4] = 0.;
			sigV[5] = 0.;
			for (j=0; j<6; j++)  	//sig = sig + 6.*I21Weights[i] * (sigV + sigDD[0] * D0 + sigDD[1] * D1 + sigDD[2] * D2) # stress tensor integration # factor 6. calibrates the weighting factors
				{
					sig_[j] = sig_[j] + 6.*I21Weights[i] * ( sigV[j] + sigDD[0] * D0[j] + sigDD[1] * D1[j] + sigDD[2] * D2[j] );
				}
			// secant stiffness
			DevStiff[0][0] =  fb9 * nn[0]*nn[0] + ob9 * nn[1]*nn[1] + ob9 * nn[2]*nn[2];  // 			DDD = DevStiffness()
			DevStiff[0][1] = -tb9 * nn[0]*nn[0] - tb9 * nn[1]*nn[1] + ob9 * nn[2]*nn[2];
			DevStiff[0][2] = -tb9 * nn[0]*nn[0] + ob9 * nn[1]*nn[1] - tb9 * nn[2]*nn[2];
			DevStiff[0][3] = -ob3 * nn[2] * nn[1];
			DevStiff[0][4] =  ob6 * nn[0] * nn[2];
			DevStiff[0][5] =  ob6 * nn[0] * nn[1];
			DevStiff[1][0] = DevStiff[0][1]; // -tb9 * nn[0] * *2 - tb9 * nn[1] * *2 + ob9 * nn[2] * *2
			DevStiff[1][1] =  ob9 * nn[0]*nn[0] + fb9 * nn[1]*nn[1] + ob9 * nn[2]*nn[2];
			DevStiff[1][2] =  ob9 * nn[0]*nn[0] - tb9 * nn[1]*nn[1] - tb9 * nn[2]*nn[2];
			DevStiff[1][3] =  ob6 * nn[2] * nn[1];
			DevStiff[1][4] = -ob3 * nn[0] * nn[2];
			DevStiff[1][5] =  ob6 * nn[0] * nn[1];
			DevStiff[2][0] = DevStiff[0][2]; // -tb9 * nn[0] * *2 + ob9 * nn[1] * *2 - tb9 * nn[2] * *2 
			DevStiff[2][1] = DevStiff[1][2]; // ob9*nn[0] * *2 - tb9 * nn[1] * *2 - tb9 * nn[2] * *2 
			DevStiff[2][2] =  ob9 * nn[0]*nn[0] + ob9 * nn[1]*nn[1] + fb9 * nn[2]*nn[2];
			DevStiff[2][3] =  ob6 * nn[2] * nn[1];
			DevStiff[2][4] =  ob6 * nn[0] * nn[2];
			DevStiff[2][5] = -ob3 * nn[0] * nn[1];
			DevStiff[3][0] = DevStiff[0][3];  // -ob3 * nn[2] * nn[1]
			DevStiff[3][1] = DevStiff[1][3];  // ob6*nn[2] * nn[1]
			DevStiff[3][2] = DevStiff[2][3];  // ob6*nn[2] * nn[1] 
			DevStiff[3][3] = ob4 * nn[2]*nn[2] + ob4 * nn[1]*nn[1];
			DevStiff[3][4] = ob4 * nn[0] * nn[1];
			DevStiff[3][5] = ob4 * nn[0] * nn[2];
			DevStiff[4][0] = DevStiff[0][4]; // ob6*nn[0] * nn[2]
			DevStiff[4][1] = DevStiff[1][4]; // -ob3 * nn[0] * nn[2]
			DevStiff[4][2] = DevStiff[2][4]; // ob6*nn[0] * nn[2]
			DevStiff[4][3] = DevStiff[3][4]; // ob4*nn[0] * nn[1]
			DevStiff[4][4] = ob4 * nn[2]*nn[2] + ob4 * nn[0]*nn[0];
			DevStiff[4][5] = ob4 * nn[2] * nn[1];
			DevStiff[5][0] = DevStiff[0][5]; // ob6*nn[0] * nn[1]
			DevStiff[5][1] = DevStiff[1][5]; // ob6*nn[0] * nn[1]
			DevStiff[5][2] = DevStiff[2][5]; // -ob3 * nn[0] * nn[1]
			DevStiff[5][3] = DevStiff[3][5]; // ob4*nn[0] * nn[2]
			DevStiff[5][4] = DevStiff[4][5]; // ob4*nn[2] * nn[1]
			DevStiff[5][5] = ob4 * nn[1]*nn[1] + ob4 * nn[0]*nn[0];
			for (j = 0; j < 6; j++) {
				for (k = 0; k < 6; k++) {
					MMatS[j][k] = MMatS[j][k] + 6.*I21Weights[i] * ((1. - dd)*E_V*VVV[j][k] + (1. - dd)*E_D*DevStiff[j][k]);// 	MMatS = MMatS + 6.*I21Weights[i] * ((1. - dd)*self.E_V*VVV + (1. - dd)*self.E_D*DDD) # presumably symmetric, i.e.only upper right part builded
				}
			}
			// corrector for tangential material stiffness 
			if ((kap > kapOld) && (dd > ZeroD)) {       			// if kap>kapOld and dd>ZeroD:
				q = 3.*kV1*epsVol;
				Q = sqrt(q*q + 1.5*kV2*(epsDD[0] * epsDD[0] + epsDD[1] * epsDD[1] + epsDD[2] * epsDD[2]));//Q = sqrt((3.*self.kV1*epsVol)**2 + 1.5*self.kV2*(epsDD[0] * epsDD[0] + epsDD[1] * epsDD[1] + epsDD[2] * epsDD[2]))
				Q = 1. / Q;
				A = Dp * (3.*kV0 + Q * 9.*kV1*kV1 * epsVol);	//A = Dp * (3.*self.kV0 + Q * (3.*self.kV1)**2 * epsVol)
				B = Dp * (1.5*Q*kV2);							//B = Dp * (1.5*Q*self.kV2)
				sigV0 = E_V * epsVol;
				c[0] =  tb3 * epsDD[0] * nn[0] - ob3 * epsDD[1] * nn[1] - ob3 * epsDD[2] * nn[2]; //c0 = tb3 * epsDD[0] * nn[0] - ob3 * epsDD[1] * nn[1] - ob3 * epsDD[2] * nn[2]
				c[1] = -ob3 * epsDD[0] * nn[0] + tb3 * epsDD[1] * nn[1] - ob3 * epsDD[2] * nn[2];//	c1 = -ob3 * epsDD[0] * nn[0] + tb3 * epsDD[1] * nn[1] - ob3 * epsDD[2] * nn[2]
				c[2] = -ob3 * epsDD[0] * nn[0] - ob3 * epsDD[1] * nn[1] + tb3 * epsDD[2] * nn[2];//	c2 = -ob3 * epsDD[0] * nn[0] - ob3 * epsDD[1] * nn[1] + tb3 * epsDD[2] * nn[2]
				c[3] =  ob2 * epsDD[1] * nn[2] + ob2 * epsDD[2] * nn[1];//	c3 = ob2 * epsDD[1] * nn[2] + ob2 * epsDD[2] * nn[1]
				c[4] =  ob2 * epsDD[0] * nn[2] + ob2 * epsDD[2] * nn[0]; //	c4 = ob2 * epsDD[0] * nn[2] + ob2 * epsDD[2] * nn[0]
				c[5] =  ob2 * epsDD[0] * nn[1] + ob2 * epsDD[1] * nn[0]; //	c5 = ob2 * epsDD[0] * nn[1] + ob2 * epsDD[1] * nn[0]
				for ( j=0; j<6; j++ ) {					//					CC, VC = CVStiffness()
					for ( k=0; k<6; k++ )
						CCStiff[j][k] = c[j] * c[k];
				}
				for (j = 0; j <3; j++ ) {
					for (k = 0; k<6; k++)
						VCStiff[j][k] = ob3 * c[k]; //	array([[ob3*c0, ob3*c1, ob3*c2, ob3*c3, ob3*c4, ob3*c5],[ob3*c0, ob3*c1, ob3*c2, ob3*c3, ob3*c4, ob3*c5],[ob3*c0, ob3*c1, ob3*c2, ob3*c3, ob3*c4, ob3*c5],
				}
				for (j = 3; j < 6; j++) {
					for (k = 0; k<6; k++)
						VCStiff[j][k] = 0.;//	[0, 0, 0, 0, 0, 0],						[0, 0, 0, 0, 0, 0],						[0, 0, 0, 0, 0, 0]])
				}
				for (j = 0; j < 6;  j++) {
					for (k = 0; k<6; k++)
						MMatT[j][k] = MMatT[j][k] - 6.*I21Weights[i] * (sigV0*(A*VVV[j][k] + B * VCStiff[j][k]) + E_D * (A*VCStiff[k][j] + B * CCStiff[j][k])); //MMatT = MMatT - 6.*I21Weights[i] * (sigV0*(A*VVV + B * VC) + self.E_D*(A*transpose(VC) + B * CC))
				}
			}
			dd_iso = dd_iso + 2.*I21Weights[i] * dd;					//dd_iso += 2.*I21Weights[i] * dd
			I1 = I1 + 2.*I21Weights[i] * I1mp;						//I1 += 2.*I21Weights[i] * I1mp
			J2 = J2 + 2.*I21Weights[i] * J2mp;						//J2 += 2.*I21Weights[i] * J2mp
			// microplane state variable updates
			*(ElemStateVarN + ns * i) = kap;						//Elem.StateVarN[ipI, ns*i] = kap
		}
		ep2 = -(MMatS[2][0]*eps[0] + MMatS[2][1]*eps[1] + MMatS[2][5]*eps[5]) / MMatS[2][2];//ep2 = -(MMatS[2, 0] * Eps[0] + MMatS[2, 1] * Eps[1] + MMatS[2, 5] * Eps[5]) / MMatS[2, 2] # lateral normal strain in case of zero lateral stress
		if (abs(eps[2]) > ZeroD*1.e-3) { xyz = abs((eps[2] - ep2) / eps[2]); } //			if abs(Eps[2])>ZeroD*1.e-3: xyz = abs((Eps[2] - ep2) / Eps[2])
		else						   { xyz = 0.; }								//else:                          xyz = 0.
		eps[2] = ep2; 												//		if PlStFlag: Eps[2] = ep2
		if ((!PlStFlag) || (xyz < PlStressL)) { break; }			//	if not PlStFlag or xyz<self.PlStressL : break
	}
	if (PlStFlag && ii >= PlStressI) { rc = 110; } //if ii >= self.PlStressI - 1 : raise NameError("ConFemMaterials::Microplane.Sig: plane stress iteration", ii, sig[2, 2])
	*(ElemStateVarN + nState * nInt) = dd_iso;					//Elem.StateVarN[ipI, self.nState*self.nInt] = dd_iso
	*(ElemStateVarN + nState * nInt + 1) = I1;						//Elem.StateVarN[ipI, self.nState*self.nInt + 1] = I1
	*(ElemStateVarN + nState * nInt + 2) = J2;						//Elem.StateVarN[ipI, self.nState*self.nInt + 2] = J2
	*(ElemStateVarN + nState * nInt + 4) = ep2;						//Elem.StateVarN[ipI, self.nState*self.nInt + 4] = ep2
//
//	*(DataOut + 0) = dd_iso;
//	*(DataOut + 1) = I1;
//	*(DataOut + 2) = J2;
//	*(DataOut + 3) = ep2;
//
	// viscous contribution
	if (etaV > 0.) {													// if self.eta>0.:                             # viscous regularization
		zz = ViscExten3DC1(Dt, etaV, Dps__, ElemStateVar, ElemStateVarN, iS, Veps);  // def ViscExten3D(self, Dt, etaV, Dps, Elem, ipI, sI):
		for (i = 0; i < 6; i++) { sig_[i] = sig_[i] + etaV * Veps[i]; }
		for (i = 0; i < 6; i++) { MMatS[i][i] = MMatS[i][i] + zz; }
	}
	// collect
	for (j = 0; j < 6; j++) {										//		MatM = MMatS + MMatT
		for (k = 0; k < 6; k++) {
			MMatS[j][k] = MMatS[j][k] + MMatT[j][k];
		}
	}
	if (ElemDim == 2) {												//		if Elem.dim == 2 : # 2D
		if (ElemPlSt) {												// 			if Elem.PlSt:                                                                   # plane stress
			cD = 1. / MMatS[2][2];									// cD = 1. / MatM[2, 2]
			*(MatM + 0) = MMatS[0][0] - MMatS[0][2] * MMatS[2][0] * cD;  *(MatM + 1) = MMatS[0][1] - MMatS[0][2] * MMatS[2][1] * cD;  *(MatM + 2) = MMatS[0][5] - MMatS[0][2] * MMatS[2][5] * cD;  //MatM_ = array([[MatM[0, 0] - MatM[0, 2] * MatM[2, 0] * cD, MatM[0, 1] - MatM[0, 2] * MatM[2, 1] * cD, MatM[0, 5] - MatM[0, 2] * MatM[2, 5] * cD],
			*(MatM + 3) = MMatS[1][0] - MMatS[1][2] * MMatS[2][0] * cD;  *(MatM + 4) = MMatS[1][1] - MMatS[1][2] * MMatS[2][1] * cD;  *(MatM + 5) = MMatS[1][5] - MMatS[1][2] * MMatS[2][5] * cD;  //[MatM[1, 0] - MatM[1, 2] * MatM[2, 0] * cD, MatM[1, 1] - MatM[1, 2] * MatM[2, 1] * cD, MatM[1, 5] - MatM[1, 2] * MatM[2, 5] * cD],
			*(MatM + 6) = MMatS[5][0] - MMatS[5][2] * MMatS[2][0] * cD;  *(MatM + 7) = MMatS[5][1] - MMatS[5][2] * MMatS[2][1] * cD;  *(MatM + 8) = MMatS[5][5] - MMatS[5][2] * MMatS[2][5] * cD;  //[MatM[5, 0] - MatM[5, 2] * MatM[2, 0] * cD, MatM[5, 1] - MatM[5, 2] * MatM[2, 1] * cD, MatM[5, 5] - MatM[5, 2] * MatM[2, 5] * cD]])
		}
		else {														// 	else:                # plane strain
			*(MatM + 0) = MMatS[0][0];  *(MatM + 1) = MMatS[0][1];  *(MatM + 2) = MMatS[0][5]; //MatM_=array([[MatM[0,0],MatM[0,1],MatM[0,5]],[MatM[1,0],MatM[1,1],MatM[1,5]],[MatM[5,0],MatM[5,1],MatM[5,5]]])
			*(MatM + 3) = MMatS[1][0];  *(MatM + 4) = MMatS[1][1];  *(MatM + 5) = MMatS[1][5];
			*(MatM + 6) = MMatS[5][0];  *(MatM + 7) = MMatS[5][1];  *(MatM + 8) = MMatS[5][5];
		}
		*(sig + 0) = sig_[0]; //  + sigV[0];							// return array([sig[0, 0], sig[1, 1], sig[0, 1]]), MatM_, [Eps[0], Eps[1], Eps[2], sig[0, 0], sig[1, 1], sig[0, 1]]
		*(sig + 1) = sig_[1]; //  + sigV[1];
		*(sig + 2) = sig_[5]; //  [1] + sigV[5];							// voigt notation broken down to 2D
	}
	else if (ElemDim == 4) {
		*(MatM + 0) = MMatS[0][0];  *(MatM + 1) = MMatS[0][1];  *(MatM + 2) = MMatS[0][2];   *(MatM + 3) = MMatS[0][5];
		*(MatM + 4) = MMatS[1][0];  *(MatM + 5) = MMatS[1][1];  *(MatM + 6) = MMatS[1][2];   *(MatM + 7) = MMatS[1][5];
		*(MatM + 8) = MMatS[2][0];  *(MatM + 9) = MMatS[2][1];  *(MatM +10) = MMatS[2][2];   *(MatM +11) = MMatS[2][5];
		*(MatM +12) = MMatS[5][0];  *(MatM +13) = MMatS[5][1];  *(MatM +14) = MMatS[5][2];   *(MatM +15) = MMatS[5][5];
		*(sig + 0) = sig_[0]; //  + sigV[0];
		*(sig + 1) = sig_[1]; //  + sigV[1];
		*(sig + 2) = sig_[2];
		*(sig + 3) = sig_[5]; //  [1] + sigV[5];							// voigt notation broken down to 2D
	}
	else if (ElemDim == 3) {										// elif Elem.dim == 3:                                                                   # 3D
		*(MatM + 0) = MMatS[0][0]; *(MatM + 1) = MMatS[0][1]; *(MatM + 2) = MMatS[0][2]; *(MatM + 3) = MMatS[0][3]; *(MatM + 4) = MMatS[0][4]; *(MatM + 5) = MMatS[0][5];
		*(MatM + 6) = MMatS[1][0]; *(MatM + 7) = MMatS[1][1]; *(MatM + 8) = MMatS[1][2]; *(MatM + 9) = MMatS[1][3]; *(MatM + 10) = MMatS[1][4]; *(MatM + 11) = MMatS[1][5];//
		*(MatM + 12) = MMatS[2][0]; *(MatM + 13) = MMatS[2][1]; *(MatM + 14) = MMatS[2][2]; *(MatM + 15) = MMatS[2][3]; *(MatM + 16) = MMatS[2][4]; *(MatM + 17) = MMatS[2][5];
		*(MatM + 18) = MMatS[3][0]; *(MatM + 19) = MMatS[3][1]; *(MatM + 20) = MMatS[3][2]; *(MatM + 21) = MMatS[3][3]; *(MatM + 22) = MMatS[3][4]; *(MatM + 23) = MMatS[3][5];
		*(MatM + 24) = MMatS[4][0]; *(MatM + 25) = MMatS[4][1]; *(MatM + 26) = MMatS[4][2]; *(MatM + 27) = MMatS[4][3]; *(MatM + 28) = MMatS[4][4]; *(MatM + 29) = MMatS[4][5];
		*(MatM + 30) = MMatS[5][0]; *(MatM + 31) = MMatS[5][1]; *(MatM + 32) = MMatS[5][2]; *(MatM + 33) = MMatS[5][3]; *(MatM + 34) = MMatS[5][4]; *(MatM + 35) = MMatS[5][5];
		*(sig + 0) = sig_[0];// +sigV[0];							// return array([sig[0,0],sig[1,1],sig[2,2],sig[1,2],sig[0,2],sig[0,1]]),MatM,[sig[0,0],sig[1,1],sig[2,2],sig[1,2],sig[0,2],sig[0,1]]
		*(sig + 1) = sig_[1];// +sigV[1];
		*(sig + 2) = sig_[2];// +sigV[2];
		*(sig + 3) = sig_[3];// +sigV[3];
		*(sig + 4) = sig_[4];// +sigV[4];
		*(sig + 5) = sig_[5];// +sigV[5];
	}
	else if (ElemDim == 21) {										// 	elif Elem.dim == 21:                                                                  # Continuum based shell
		if (abs(MMatS[2][2]) < ZeroD) { return 220; }
		cD = 1. / MMatS[2][2];										// 		cD = 1. / MatM[2, 2]
		*(MatM + 0) = MMatS[0][0] - MMatS[0][2] * MMatS[2][0] * cD; *(MatM + 1) = MMatS[0][1] - MMatS[0][2] * MMatS[2][1] * cD; *(MatM + 2) = 0.; *(MatM + 3) = MMatS[0][3] - MMatS[0][2] * MMatS[2][3] * cD; *(MatM + 4) = MMatS[0][4] - MMatS[0][2] * MMatS[2][4] * cD; *(MatM + 5) = MMatS[0][5] - MMatS[0][2] * MMatS[2][5] * cD;
		*(MatM + 6) = MMatS[1][0] - MMatS[1][2] * MMatS[2][0] * cD; *(MatM + 7) = MMatS[1][1] - MMatS[1][2] * MMatS[2][1] * cD; *(MatM + 8) = 0.; *(MatM + 9) = MMatS[1][3] - MMatS[1][2] * MMatS[2][3] * cD; *(MatM + 10) = MMatS[1][4] - MMatS[1][2] * MMatS[2][4] * cD; *(MatM + 11) = MMatS[1][5] - MMatS[1][2] * MMatS[2][5] * cD;
		*(MatM + 12) = 0.;                                       *(MatM + 13) = 0.;                                       *(MatM + 14) = 0.; *(MatM + 15) = 0.;                                       *(MatM + 16) = 0.;                                       *(MatM + 17) = 0.;
		*(MatM + 18) = MMatS[3][0] - MMatS[3][2] * MMatS[2][0] * cD; *(MatM + 19) = MMatS[3][1] - MMatS[3][2] * MMatS[2][1] * cD; *(MatM + 20) = 0.; *(MatM + 21) = MMatS[3][3] - MMatS[3][2] * MMatS[2][3] * cD; *(MatM + 22) = MMatS[3][4] - MMatS[3][2] * MMatS[2][4] * cD; *(MatM + 23) = MMatS[3][5] - MMatS[3][2] * MMatS[2][5] * cD;
		*(MatM + 24) = MMatS[4][0] - MMatS[4][2] * MMatS[2][0] * cD; *(MatM + 25) = MMatS[4][1] - MMatS[4][2] * MMatS[2][1] * cD; *(MatM + 26) = 0.; *(MatM + 27) = MMatS[4][3] - MMatS[4][2] * MMatS[2][3] * cD; *(MatM + 28) = MMatS[4][4] - MMatS[4][2] * MMatS[2][4] * cD; *(MatM + 29) = MMatS[4][5] - MMatS[4][2] * MMatS[2][5] * cD;
		*(MatM + 30) = MMatS[5][0] - MMatS[5][2] * MMatS[2][0] * cD; *(MatM + 31) = MMatS[5][1] - MMatS[5][2] * MMatS[2][1] * cD; *(MatM + 32) = 0.; *(MatM + 33) = MMatS[5][3] - MMatS[5][2] * MMatS[2][3] * cD; *(MatM + 34) = MMatS[5][4] - MMatS[5][2] * MMatS[2][4] * cD; *(MatM + 35) = MMatS[5][5] - MMatS[5][2] * MMatS[2][5] * cD;
		*(sig + 0) = sig_[0]; //+sigV[0];							// return array([sig[0,0],sig[1,1],sig[2,2],sig[1,2],sig[0,2],sig[0,1]]),MatM,[sig[0,0],sig[1,1],sig[2,2],sig[1,2],sig[0,2],sig[0,1]]
		*(sig + 1) = sig_[1]; //+sigV[1];
		*(sig + 2) = sig_[2]; //+sigV[2];
		*(sig + 3) = sig_[3]; //+sigV[3];
		*(sig + 4) = sig_[4]; //+sigV[4];
		*(sig + 5) = sig_[5]; //+sigV[5];
	}
	else { return 225; }											// else: raise NameError("ConFemMaterials::Microplane.Sig: not implemented for this element type")
	return rc;
}


