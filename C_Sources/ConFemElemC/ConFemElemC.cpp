//

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <fstream>
#include <iostream>
#include <iomanip>

#include "ConFemElemC.h"
using namespace std;
const double ZeroD = 1.e-9;

int CPS4FormBC(double X0, double Y0, double X1, double Y1, double X2, double Y2, double X3, double Y3, double r, double s, double t, double *BC, int dim_BC0, int dim_BC1,
	double *det_, int dim_det_, double* Data, int dim_D)
{
	double h00, h01, h10, h11, h20, h21, h30, h31, JAC00, JAC01, JAC10, JAC11, det, deti, a1, a2, b1, b2, B00, B10, B20, B30, B01, B11, B21, B31;
	h00 = ( 1 + s)*0.25;
	h01 = ( 1 + r)*0.25;
	h10 =-( 1 + s)*0.25;
	h11 = ( 1 - r)*0.25;
	h20 = (-1 + s)*0.25;
	h21 = (-1 + r)*0.25;
	h30 = ( 1 - s)*0.25;
	h31 =-( 1 + r)*0.25;
	JAC00 = h00 * X0 + h10 * X1 + h20 * X2 + h30 * X3;
	JAC01 = h00 * Y0 + h10 * Y1 + h20 * Y2 + h30 * Y3;
	JAC10 = h01 * X0 + h11 * X1 + h21 * X2 + h31 * X3;
	JAC11 = h01 * Y0 + h11 * Y1 + h21 * Y2 + h31 * Y3;
	det = JAC00 * JAC11 - JAC01 * JAC10;
	if (fabs(det)<ZeroD) { return 11; }
	deti = 1. / det;
	a1 = Y0*h01 + Y1*h11 + Y2*h21 + Y3*h31;
	a2 = Y0*h00 + Y1*h10 + Y2*h20 + Y3*h30;
	b1 = X0*h01 + X1*h11 + X2*h21 + X3*h31;
	b2 = X0*h00 + X1*h10 + X2*h20 + X3*h30;
	B00 = deti * (h00 * a1 - h01 * a2);
	B10 = deti * (h10 * a1 - h11 * a2);
	B20 = deti * (h20 * a1 - h21 * a2);
	B30 = deti * (h30 * a1 - h31 * a2);
	B01 = deti * (-h00 * b1 + h01 * b2);
	B11 = deti * (-h10 * b1 + h11 * b2);
	B21 = deti * (-h20 * b1 + h21 * b2);
	B31 = deti * (-h30 * b1 + h31 * b2);
/*
	*(Data + 0) = B00;
	*(Data + 1) = B10;
	*(Data + 2) = B20;
	*(Data + 3) = B30;
	*(Data + 4) = B01;
	*(Data + 5) = B11;
	*(Data + 6) = B21;
	*(Data + 7) = B31;
*/
	*det_ = det;
	*(BC + 0)     = B00;
	*(BC + 0 + 2) = B10;
	*(BC + 0 + 4) = B20;
	*(BC + 0 + 6) = B30;
	*(BC + 8 + 1) = B01;
	*(BC + 8 + 3) = B11;
	*(BC + 8 + 5) = B21;
	*(BC + 8 + 7) = B31;
	*(BC + 16+ 0) = B01;
	*(BC + 16+ 1) = B00;
	*(BC + 16+ 2) = B11;
	*(BC + 16+ 3) = B10;
	*(BC + 16+ 4) = B21;
	*(BC + 16+ 5) = B20;
	*(BC + 16+ 6) = B31;
	*(BC + 16+ 7) = B30;
//	   [[B00, 0,   B10, 0,   B20, 0,   B30, 0],
//		[0,   B01, 0,   B11, 0,   B21, 0,   B31],
//		[B01, B00, B11, B10, B21, B20, B31, B30]])
	return 10;
}



/*
def Formb(self, r, s, t) :
	return B00, B10, B20, B30, B01, B11, B21, B31, det

def FormB(self, r, s, t, val) :
	B00, B10, B20, B30, B01, B11, B21, B31, det = self.Formb(r, s, t)
	B = array([[B00, 0, B10, 0, B20, 0, B30, 0],
	[0, B01, 0, B11, 0, B21, 0, B31],
	[B01, B00, B11, B10, B21, B20, B31, B30]])
	return B, det

def DetJ(self, r, s, t) :
	h00 = -(1 - s)*0.25
	h01 = -(1 - r)*0.25
	h10 = (1 - s)*0.25
	h11 = -(1 + r)*0.25
	h20 = (1 + s)*0.25
	h21 = (1 + r)*0.25
	h30 = -(1 + s)*0.25
	h31 = (1 - r)*0.25
	JAC00 = h00 * self.X0 + h10 * self.X1 + h20 * self.X2 + h30 * self.X3
	JAC01 = h00 * self.Y0 + h10 * self.Y1 + h20 * self.Y2 + h30 * self.Y3
	JAC10 = h01 * self.X0 + h11 * self.X1 + h21 * self.X2 + h31 * self.X3
	JAC11 = h01 * self.Y0 + h11 * self.Y1 + h21 * self.Y2 + h31 * self.Y3
	det = JAC00 * JAC11 - JAC01 * JAC10
	return det



*/



int SH4BasicsC( double r, double s, double t, double *N, double *br, double *bs, double *JJ, double *JI, double *vv, double *XX, double *a, double *Vn, double *EdgeDir)
{
	double detJ, ll, x0, x1, x2, xx, t2;
	int k, k4, k3;
	*(N+0)  =  (1-r) *(1-s)*0.25;  // N =  array([(1-r)*(1-s)*0.25, (1+r)*(1-s)*0.25, (1+r)*(1+s)*0.25, (1-r)*(1+s)*0.25]) 
	*(N+1)  =  (1+r) *(1-s)*0.25;
	*(N+2)  =  (1+r) *(1+s)*0.25;
	*(N+3)  =  (1-r) *(1+s)*0.25;
	*(br+0) =  (-1+s)*0.25;          // br = array([(-1+s)*0.25,      ( 1-s)*0.25,     ( 1+s)*0.25,      -( 1+s)*0.25]) 
	*(br+1) =  ( 1-s)*0.25;
	*(br+2) =  ( 1+s)*0.25;
	*(br+3) = -( 1+s)*0.25;
	*(bs+0) =  (-1+r)*0.25;         // bs = array([(-1+r)*0.25,     -( 1+r)*0.25,     ( 1+r)*0.25,       ( 1-r)*0.25]) 
	*(bs+1) = -( 1+r)*0.25;
	*(bs+2) =  ( 1+r)*0.25;
	*(bs+3) =  ( 1-r)*0.25;

	for (k=0;k<3;k++) {
		t2 = t/2.;
		k4 = 4*k;
		k3 = 3*k;
		*(JJ+k3)=*(br)*((*(XX+k4))+t2*(*(a))*(*(Vn+k4)))+*(br+1)*((*(XX+k4+1))+t2*(*(a+1))*(*(Vn+k4+1)))+*(br+2)*((*(XX+k4+2))+t2*(*(a+2))*(*(Vn+k4+2)))+*(br+3)*((*(XX+k4+3))+t2*(*(a+3))*(*(Vn+k4+3))); 
		*(JJ+k3+1)=*(bs)*((*(XX+k4))+t2*(*(a))*(*(Vn+k4)))+*(bs+1)*((*(XX+k4+1))+t2*(*(a+1))*(*(Vn+k4+1)))+*(bs+2)*((*(XX+k4+2))+t2*(*(a+2))*(*(Vn+k4+2)))+*(bs+3)*((*(XX+k4+3))+t2*(*(a+3))*(*(Vn+k4+3)));
		*(JJ+k3+2)=(*(N))*(0.5*(*(a))*(*(Vn+k4)))+(*(N+1))*(0.5*(*(a+1))*(*(Vn+k4+1)))+(*(N+2))*(0.5*(*(a+2))*(*(Vn+k4+2)))+(*(N+3))*(0.5*(*(a+3))*(*(Vn+k4+3)));
//		printf("%4i %6.3f %6.3f %6.3f    %6.3f\n",k,*(N+0),*(a+0),*(Vn+4*k+0),*(JJ+3*k+2));
	}   // JJ = zeros((3,3),dtype=float)
// for k in xrange(3): JJ[k,0]    =br[0]  *(self.XX[k,0]+t/2.*self.a[0]*self.Vn[k,0])+br[1]*(self.XX[k,1]+t/2.*self.a[1]*self.Vn[k,1])+br[2]*(self.XX[k,2]+t/2.*self.a[2]*self.Vn[k,2])+br[3]*(self.XX[k,3]+t/2.*self.a[3]*self.Vn[k,3])
// for k in xrange(3): JJ[k,1]=bs[0]*(self.XX[k,0]+t/2.*self.a[0]*self.Vn[k,0])+bs[1]*(self.XX[k,1]+t/2.*self.a[1]*self.Vn[k,1])+bs[2]*(self.XX[k,2]+t/2.*self.a[2]*self.Vn[k,2])+bs[3]*(self.XX[k,3]+t/2.*self.a[3]*self.Vn[k,3])
// for k in xrange(3): JJ[k,2]=N[0]*(1/2.*self.a[0]*self.Vn[k,0])+N[1]*(1/2.*self.a[1]*self.Vn[k,1])+N[2]*(1/2.*self.a[2]*self.Vn[k,2])+N[3]*(1/2.*self.a[3]*self.Vn[k,3])

	detJ = (*(JJ))*(*(JJ+3+1))*(*(JJ+6+2))-(*(JJ))*(*(JJ+3+2))*(*(JJ+6+1))-(*(JJ+3))*(*(JJ+1))*(*(JJ+6+2))+(*(JJ+3))*(*(JJ+2))*(*(JJ+6+1))+(*(JJ+6))*(*(JJ+3+2))*(*(JJ+1))-(*(JJ+6))*(*(JJ+3+1))*(*(JJ+2));
	detJ = 1./detJ;
	if (fabs(detJ)<ZeroD) {return 11;}
	(*(JI))     = detJ*((*(JJ+3+1))*(*(JJ+6+2))-(*(JJ+3+2))*(*(JJ+6+1)));  //         JI = inv(JJ)
	(*(JI+1))   =-detJ*((*(JJ+1))*(*(JJ+6+2))-(*(JJ+2))*(*(JJ+6+1)));
	(*(JI+2))   = detJ*((*(JJ+1))*(*(JJ+3+2))-(*(JJ+2))*(*(JJ+3+1)));
	(*(JI+3))   =-detJ*((*(JJ+3))*(*(JJ+6+2))-(*(JJ+3+2))*(*(JJ+6))); 
	(*(JI+3+1)) = detJ*((*(JJ))*(*(JJ+6+2))-(*(JJ+2))*(*(JJ+6)));
	(*(JI+3+2)) =-detJ*((*(JJ))*(*(JJ+3+2))-(*(JJ+2))*(*(JJ+3)));
	(*(JI+6))   = detJ*((*(JJ+3))*(*(JJ+6+1))-(*(JJ+3+1))*(*(JJ+6)));
	(*(JI+6+1)) =-detJ*((*(JJ))*(*(JJ+6+1))-(*(JJ+1))*(*(JJ+6)));
	(*(JI+6+2)) = detJ*((*(JJ))*(*(JJ+3+1))-(*(JJ+1))*(*(JJ+3)));
	ll = sqrt((*(JJ+2))*(*(JJ+2))+(*(JJ+3+2))*(*(JJ+3+2))+(*(JJ+6+2))*(*(JJ+6+2)));   // ll=sqrt(JJ[0,2]**2+JJ[1,2]**2+JJ[2,2]**2)   # normal of local coordinate system, 3rd column
	if (fabs(ll)<ZeroD) {return 12;}
	ll = 1./ll;

	(*(vv))     = 0.;							// vv[0,0]     vv = array([[0.,0.,JJ[0,2]/ll],[0.,0.,JJ[1,2]/ll],[0.,0.,JJ[2,2]/ll]]) 
	(*(vv+1))   = 0.;							// vv[0,1]
	(*(vv+2))   = (*(JJ+2))*ll;					// vv[0,2]
	(*(vv+3))   = 0.;							// vv[1,0]
	(*(vv+3+1)) = 0.;							// vv[1,1]
	(*(vv+3+2)) = (*(JJ+3+2))*ll;				// vv[1,2]
	(*(vv+6))   = 0.;							// vv[2,0]
	(*(vv+6+1)) = 0.;							// vv[2,1]
	(*(vv+6+2)) = (*(JJ+6+2))*ll;				// vv[2,2]
	if (0) {             //        LoC = True         if LoC:                                     # local right handed coordinate system with 1st direction / column aligned to element edge
		x0 = (*(EdgeDir+1))*(*(vv+6+2))-(*(EdgeDir+2))*(*(vv+3+2));   //  x0 = self.EdgeDir[1]*vv[2,2]-self.EdgeDir[2]*vv[1,2]
        x1 = (*(EdgeDir+2))*(*(vv+2))-(*(EdgeDir))*(*(vv+6+2));   //  x1 = self.EdgeDir[2]*vv[0,2]-self.EdgeDir[0]*vv[2,2]
        x2 = (*(EdgeDir))*(*(vv+3+2))-(*(EdgeDir+1))*(*(vv+2));   //  x2 = self.EdgeDir[0]*vv[1,2]-self.EdgeDir[1]*vv[0,2]
        xx = sqrt(x0*x0+x1*x1+x2*x2);                 //  xx = sqrt(x0**2+x1**2+x2**2)
		if (xx<ZeroD) {return 13;}
		xx = 1./xx;
        (*(vv+1))   = -x0*xx;                             //  vv[0,1] = -x0/xx                        # 2nd column, approx perp. to element edge, reversed in sign to preserve right handedness
        (*(vv+3+1)) = -x1*xx;                             //  vv[1,1] = -x1/xx
        (*(vv+6+1)) = -x2*xx;                             //  vv[2,1] = -x2/xx 
        x0 = (*(vv+3+2))*(*(vv+6+1))-(*(vv+6+2))*(*(vv+3+1));         //  x0 = vv[1,2]*vv[2,1]-vv[2,2]*vv[1,1]
        x1 = (*(vv+6+2))*(*(vv+1))-(*(vv+2))*(*(vv+6+1));         //  x1 = vv[2,2]*vv[0,1]-vv[0,2]*vv[2,1]
        x2 = (*(vv+2))*(*(vv+3+1))-(*(vv+3+2))*(*(vv+1));         //  x2 = vv[0,2]*vv[1,1]-vv[1,2]*vv[0,1]
        xx = sqrt(x0*x0+x1*x1+x2*x2);                 //  xx = sqrt(x0**2+x1**2+x2**2)
		if (xx<ZeroD) {return 14;}
		xx = 1./xx;
        (*(vv))    = -x0*xx;                             //  vv[0,0] = -x0/xx                        # 1st column, approx aligned to element edge, sign reversal of 2nd column is implicitely corrected
        (*(vv+3))  = -x1*xx;                             //  vv[1,0] = -x1/xx
        (*(vv+6))  = -x2*xx;                             //  vv[2,0] = -x2/xx
	}
	else {
		if(abs((*(vv+3+2)))<0.99) {							 //            if abs(vv[1,2])<0.99:                       # local coordinate system V_1 from cross product of V_n and e_y ( V_1 in e_x - e_z plane) if V_n is not to much aligned to e_y  
			xx = sqrt((*(vv+6+2))*(*(vv+6+2))+(*(vv+2))*(*(vv+2)));	 //           ll = sqrt(vv[2,2]**2+vv[0,2]**2)        # length of V_1
			(*(vv)) = (*(vv+6+2))/xx;						 //                vv[0,0] = vv[2,2]/ll                    # V_1[0] normalized;  V1[1] = 0
			(*(vv+6)) =-(*(vv+2))/xx;						 //                vv[2,0] =-vv[0,2]/ll                    # V_1[2] normalized
			(*(vv+1))  = (*(vv+3+2))*(*(vv+6));					 //                vv[0,1] = vv[1,2]*vv[2,0]               # as V_n and V_1 are orthogonal and both have unit length V_2 also should have unit length
			(*(vv+3+1)) = (*(vv+6+2))*(*(vv))-(*(vv+2))*(*(vv+6));	 //                vv[1,1] = vv[2,2]*vv[0,0]-vv[0,2]*vv[2,0]
			(*(vv+6+1)) =-(*(vv+3+2))*(*(vv)); //                vv[2,1] =-vv[1,2]*vv[0,0]
		}
		else { //            else:                                       # local coordinate system V_1 from cross product of V_n and e_x ( V_1 in e_y - e_z plane)
			xx = sqrt((*(vv+6+2))*(*(vv+6+2))+(*(vv+3+2))*(*(vv+3+2)));  //                        # length of V_1
			(*(vv+3)) =-(*(vv+6+2))/xx;						 //                vv[1,0] =-vv[2,2]/ll                    # V_1[0] normalized;  V1[0] = 0
			(*(vv+6)) = (*(vv+3+2))/xx;						 //                vv[2,0] = vv[1,2]/ll                    # V_1[2] normalized
			(*(vv+1))  = (*(vv+3+2))*(*(vv+6))-(*(vv+6+2))*(*(vv+3));	 //                vv[0,1] = vv[1,2]*vv[2,0]-vv[2,2]*vv[1,0] # as V_n and V_1 are orthogonal and both have unit length V_2 also should have unit length
			(*(vv+3+1)) =-(*(vv+2))*(*(vv+6));					 //                vv[1,1] =-vv[0,2]*vv[2,0]
			(*(vv+6+1)) = (*(vv+2))*(*(vv+3));					 //                vv[2,1] = vv[0,2]*vv[1,0]
		}
	}
	return 1;
}
int SH4FormBC( double r, double s, double t, double *BB, int dim_BB0, int dim_BB1, double *Data, int dim_D, double *XX, int dim_X0, int dim_X1, double *a, int dim_a, double *Vn, int dim_Vn0,int dim_Vn1,
				double *EdgeDir, int dim_Ed, double *gg0, int dim_g00, int dim_g01, double *gg1, int dim_g10, int dim_g11, double *TD, int dim_TD0, int dim_TD1)  //  def FormB(self, r, s, t_, NLg):   t = t_
{
	int rc, k;
	double N[4], br[4], bs[4], JJ[9], JI[9], vv[9], HH[3][2][4], det; // HH = zeros((3,2,4),dtype=float)
	double J01D, J11D, J21D, J01B, J11B, J21B, J00A, J10A, J20A, J00C, J10C, J20C, J02D, J12D, J22D, J02B, J12B, J22B, J02A, J12A, J22A, J02C, J12C, J22C, t25;
	double td[3][3];

	rc = SH4BasicsC( r, s, t, N, br, bs, JJ, JI, vv, XX, a, Vn, EdgeDir);  // N, br, bs, JJ, JI, vv = self.Basics( r, s, t)
//	for (k=0;k<12;k++) { *(Data+k) = JJ[k]; }
	if (rc>10) { return rc; }
/* 
	JJ[0,0] -> J[0], JJ[0,1] -> J[1], JJ[0,2] -> J[2] 
	JJ[1,0] -> J[3], JJ[1,1] -> J[4], JJ[1,2] -> J[5]
	JJ[2,0] -> J[6], JJ[2,1] -> J[7], JJ[2,2] -> J[8] 
*/
	det = JJ[0]*JJ[4]*JJ[8]-JJ[0]*JJ[5]*JJ[7]-JJ[3]*JJ[1]*JJ[8]+JJ[3]*JJ[2]*JJ[7]+JJ[6]*JJ[1]*JJ[5]-JJ[6]*JJ[2]*JJ[4]; // det = JJ[0,0]*JJ[1,1]*JJ[2,2]-JJ[0,0]*JJ[1,2]*JJ[2,1]-JJ[1,0]*JJ[0,1]*JJ[2,2]+JJ[1,0]*JJ[0,2]*JJ[2,1]+JJ[2,0]*JJ[0,1]*JJ[1,2]-JJ[2,0]*JJ[0,2]*JJ[1,1] 
	*(Data+0)=det;
	for (k=0;k<4;k++) {             // for k in xrange(4):
		HH[0][0][k]=JJ[0]*(*(gg0+k))+JJ[3]*(*(gg0+4+k))+JJ[6]*(*(gg0+8+k));   //  HH[0,0,k]=JJ[0,0]*self.gg[0,0,k]+JJ[1,0]*self.gg[0,1,k]+JJ[2,0]*self.gg[0,2,k]
        HH[0][1][k]=JJ[0]*(*(gg1+k))+JJ[3]*(*(gg1+4+k))+JJ[6]*(*(gg1+8+k));   //  HH[0,1,k]=JJ[0,0]*self.gg[1,0,k]+JJ[1,0]*self.gg[1,1,k]+JJ[2,0]*self.gg[1,2,k]
        HH[1][0][k]=JJ[1]*(*(gg0+k))+JJ[4]*(*(gg0+4+k))+JJ[7]*(*(gg0+8+k));   //  HH[1,0,k]=JJ[0,1]*self.gg[0,0,k]+JJ[1,1]*self.gg[0,1,k]+JJ[2,1]*self.gg[0,2,k]
        HH[1][1][k]=JJ[1]*(*(gg1+k))+JJ[4]*(*(gg1+4+k))+JJ[7]*(*(gg1+8+k));   //  HH[1,1,k]=JJ[0,1]*self.gg[1,0,k]+JJ[1,1]*self.gg[1,1,k]+JJ[2,1]*self.gg[1,2,k]
        HH[2][0][k]=JJ[2]*(*(gg0+k))+JJ[5]*(*(gg0+4+k))+JJ[8]*(*(gg0+8+k));   //  HH[2,0,k]=JJ[0,2]*self.gg[0,0,k]+JJ[1,2]*self.gg[0,1,k]+JJ[2,2]*self.gg[0,2,k]  
		HH[2][1][k]=JJ[2]*(*(gg1+k))+JJ[5]*(*(gg1+4+k))+JJ[8]*(*(gg1+8+k));   //  HH[2,1,k]=JJ[0,2]*self.gg[1,0,k]+JJ[1,2]*self.gg[1,1,k]+JJ[2,2]*self.gg[1,2,k] 
//printf("%4i %6.3f %6.3f %6.3f\n",k,(*(gg0+0*4+k)),(*(gg0+1*4+k)),(*(gg0+2*4+k)) );

	}
	for (k=0;k<4;k++) {             // BB = zeros((6,20),dtype=float) for k in xrange(4):
        *(BB  +k*5)     = JJ[0]*br[k];                               //        BB[0,k*5+0]=JJ[0,0]*br[k]
        *(BB  +k*5+1)   = JJ[3]*br[k];                               //        BB[0,k*5+1]=JJ[1,0]*br[k]
        *(BB  +k*5+2)   = JJ[6]*br[k];                               //        BB[0,k*5+2]=JJ[2,0]*br[k]
        *(BB  +k*5+3)   =               t*br[k]*HH[0][0][k];         //        BB[0,k*5+3]=               t*br[k]*HH[0,0,k]
        *(BB  +k*5+4)   =               t*br[k]*HH[0][1][k];         //        BB[0,k*5+4]=               t*br[k]*HH[0,1,k]
        *(BB+20 +k*5)   = JJ[1]*bs[k];                               //        BB[1,k*5+0]=JJ[0,1]*bs[k]
        *(BB+20 +k*5+1) = JJ[4]*bs[k];                               //        BB[1,k*5+1]=JJ[1,1]*bs[k]
        *(BB+20 +k*5+2) = JJ[7]*bs[k];                               //        BB[1,k*5+2]=JJ[2,1]*bs[k]
        *(BB+20 +k*5+3) =               t*bs[k]*HH[1][0][k];         //        BB[1,k*5+3]=               t*bs[k]*HH[1,0,k]
        *(BB+20 +k*5+4) =               t*bs[k]*HH[1][1][k];         //        BB[1,k*5+4]=               t*bs[k]*HH[1,1,k]
        *(BB+40 +k*5+3) = N[k]*HH[2][0][k];                          //        BB[2,k*5+3]=N[k]*HH[2,0,k]
        *(BB+40 +k*5+4) = N[k]*HH[2][1][k];                          //        BB[2,k*5+4]=N[k]*HH[2,1,k]
        *(BB+100+k*5)   = JJ[0]*bs[k]  +JJ[1]*br[k];                 //        BB[5,k*5+0]=JJ[0,0]*bs[k]  +JJ[0,1]*br[k]        
        *(BB+100+k*5+1) = JJ[3]*bs[k]  +JJ[4]*br[k];                 //        BB[5,k*5+1]=JJ[1,0]*bs[k]  +JJ[1,1]*br[k]
        *(BB+100+k*5+2) = JJ[6]*bs[k]  +JJ[7]*br[k];                 //        BB[5,k*5+2]=JJ[2,0]*bs[k]  +JJ[2,1]*br[k]
        *(BB+100+k*5+3) =     t*(bs[k]*HH[0][0][k]+br[k]*HH[1][0][k]);//       BB[5,k*5+3]=     t*(bs[k]*HH[0,0,k]+br[k]*HH[1,0,k])
        *(BB+100+k*5+4) =     t*(bs[k]*HH[0][1][k]+br[k]*HH[1][1][k]);//       BB[5,k*5+4]=     t*(bs[k]*HH[0,1,k]+br[k]*HH[1,1,k])
	}
	// flag = True    #  True --> Assumed-Natural-Strain-Method to avoid transverse shear locking          if not flag:   .... else:
	t25 = 0.25*t;
    J01D =-0.5*(*(XX+1)  -*(XX+2))  -t25*(a[1]*(*(Vn+1))  -a[2]*(*(Vn+2))); //        J01D=-0.5*(self.XX[0,1]-self.XX[0,2])-0.25*t*(self.a[1]*self.Vn[0,1]-self.a[2]*self.Vn[0,2])
    J11D =-0.5*(*(XX+4+1)-*(XX+4+2))-t25*(a[1]*(*(Vn+4+1))-a[2]*(*(Vn+4+2))); //        J11D=-0.5*(self.XX[1,1]-self.XX[1,2])-0.25*t*(self.a[1]*self.Vn[1,1]-self.a[2]*self.Vn[1,2])
    J21D =-0.5*(*(XX+8+1)-*(XX+8+2))-t25*(a[1]*(*(Vn+8+1))-a[2]*(*(Vn+8+2))); //        J21D=-0.5*(self.XX[2,1]-self.XX[2,2])-0.25*t*(self.a[1]*self.Vn[2,1]-self.a[2]*self.Vn[2,2])
    J01B =-0.5*(*(XX)    -*(XX+3))  -t25*(a[0]*(*(Vn))    -a[3]*(*(Vn+3))); //        J01B=-0.5*(self.XX[0,0]-self.XX[0,3])-0.25*t*(self.a[0]*self.Vn[0,0]-self.a[3]*self.Vn[0,3])
    J11B =-0.5*(*(XX+4)  -*(XX+4+3))-t25*(a[0]*(*(Vn+4))  -a[3]*(*(Vn+4+3))); //        J11B=-0.5*(self.XX[1,0]-self.XX[1,3])-0.25*t*(self.a[0]*self.Vn[1,0]-self.a[3]*self.Vn[1,3])
    J21B =-0.5*(*(XX+8)  -*(XX+8+3))-t25*(a[0]*(*(Vn+8))  -a[3]*(*(Vn+8+3))); //        J21B=-0.5*(self.XX[2,0]-self.XX[2,3])-0.25*t*(self.a[0]*self.Vn[2,0]-self.a[3]*self.Vn[2,3])       
    J00A = 0.5*(*(XX+2)  -*(XX+3))  +t25*(a[2]*(*(Vn+2))  -a[3]*(*(Vn+3))); //        J00A= 0.5*(self.XX[0,2]-self.XX[0,3])+0.25*t*(self.a[2]*self.Vn[0,2]-self.a[3]*self.Vn[0,3])
    J10A = 0.5*(*(XX+4+2)-*(XX+4+3))+t25*(a[2]*(*(Vn+4+2))-a[3]*(*(Vn+4+3))); //        J10A= 0.5*(self.XX[1,2]-self.XX[1,3])+0.25*t*(self.a[2]*self.Vn[1,2]-self.a[3]*self.Vn[1,3])
    J20A = 0.5*(*(XX+8+2)-*(XX+8+3))+t25*(a[2]*(*(Vn+8+2))-a[3]*(*(Vn+8+3))); //        J20A= 0.5*(self.XX[2,2]-self.XX[2,3])+0.25*t*(self.a[2]*self.Vn[2,2]-self.a[3]*self.Vn[2,3])
    J00C =-0.5*(*(XX)    -*(XX+1))  -t25*(a[0]*(*(Vn))    -a[1]*(*(Vn+1))); //        J00C=-0.5*(self.XX[0,0]-self.XX[0,1])-0.25*t*(self.a[0]*self.Vn[0,0]-self.a[1]*self.Vn[0,1])
    J10C =-0.5*(*(XX+4)  -*(XX+4+1))-t25*(a[0]*(*(Vn+4))  -a[1]*(*(Vn+4+1))); //        J10C=-0.5*(self.XX[1,0]-self.XX[1,1])-0.25*t*(self.a[0]*self.Vn[1,0]-self.a[1]*self.Vn[1,1])
    J20C =-0.5*(*(XX+8)  -*(XX+8+1))-t25*(a[0]*(*(Vn+8))  -a[1]*(*(Vn+8+1))); //        J20C=-0.5*(self.XX[2,0]-self.XX[2,1])-0.25*t*(self.a[0]*self.Vn[2,0]-self.a[1]*self.Vn[2,1])
    J02D = 0.25*(a[1]*(*(Vn+1))  +a[2]*(*(Vn+2)));                         //        J02D=0.25*(self.a[1]*self.Vn[0,1]+self.a[2]*self.Vn[0,2])
    J12D = 0.25*(a[1]*(*(Vn+4+1))+a[2]*(*(Vn+4+2)));                         //        J12D=0.25*(self.a[1]*self.Vn[1,1]+self.a[2]*self.Vn[1,2])
    J22D = 0.25*(a[1]*(*(Vn+8+1))+a[2]*(*(Vn+8+2)));                         //        J22D=0.25*(self.a[1]*self.Vn[2,1]+self.a[2]*self.Vn[2,2])
    J02B = 0.25*(a[0]*(*(Vn))    +a[3]*(*(Vn+3)));                         //        J02B=0.25*(self.a[0]*self.Vn[0,0]+self.a[3]*self.Vn[0,3])
    J12B = 0.25*(a[0]*(*(Vn+4))  +a[3]*(*(Vn+4+3)));                         //        J12B=0.25*(self.a[0]*self.Vn[1,0]+self.a[3]*self.Vn[1,3])
    J22B = 0.25*(a[0]*(*(Vn+8))  +a[3]*(*(Vn+8+3)));                         //        J22B=0.25*(self.a[0]*self.Vn[2,0]+self.a[3]*self.Vn[2,3])
    J02A = 0.25*(a[2]*(*(Vn+2))  +a[3]*(*(Vn+3)));                         //        J02A=0.25*(self.a[2]*self.Vn[0,2]+self.a[3]*self.Vn[0,3])
    J12A = 0.25*(a[2]*(*(Vn+4+2))+a[3]*(*(Vn+4+3)));                         //        J12A=0.25*(self.a[2]*self.Vn[1,2]+self.a[3]*self.Vn[1,3])
    J22A = 0.25*(a[2]*(*(Vn+8+2))+a[3]*(*(Vn+8+3)));                         //        J22A=0.25*(self.a[2]*self.Vn[2,2]+self.a[3]*self.Vn[2,3])
    J02C = 0.25*(a[0]*(*(Vn))    +a[1]*(*(Vn+1)));                         //        J02C=0.25*(self.a[0]*self.Vn[0,0]+self.a[1]*self.Vn[0,1])
    J12C = 0.25*(a[0]*(*(Vn+4))  +a[1]*(*(Vn+4+1)));                         //        J12C=0.25*(self.a[0]*self.Vn[1,0]+self.a[1]*self.Vn[1,1])
    J22C = 0.25*(a[0]*(*(Vn+8))  +a[1]*(*(Vn+8+1)));                         //        J22C=0.25*(self.a[0]*self.Vn[2,0]+self.a[1]*self.Vn[2,1])

	*(BB+60)   =-0.25*(1-r)*J02B;
    *(BB+60+1) =-0.25*(1-r)*J12B;
    *(BB+60+2) =-0.25*(1-r)*J22B;
    *(BB+60+3) = 0.25*(1-r)*((J01B*(*(gg0))+J11B*(*(gg0+4))+J21B*(*(gg0+8)))-t*(J02B*(*(gg0))+J12B*(*(gg0+4))+J22B*(*(gg0+8))));
    *(BB+60+4) = 0.25*(1-r)*((J01B*(*(gg1))+J11B*(*(gg1+4))+J21B*(*(gg1+8)))-t*(J02B*(*(gg1))+J12B*(*(gg1+4))+J22B*(*(gg1+8))));
    *(BB+60+5) =-0.25*(1+r)*J02D; 
    *(BB+60+6) =-0.25*(1+r)*J12D;
    *(BB+60+7) =-0.25*(1+r)*J22D;
    *(BB+60+8) = 0.25*(1+r)*((J01D*(*(gg0+1))+J11D*(*(gg0+4+1))+J21D*(*(gg0+8+1)))-t*(J02D*(*(gg0+1))+J12D*(*(gg0+4+1))+J22D*(*(gg0+8+1))));
    *(BB+60+9) = 0.25*(1+r)*((J01D*(*(gg1+1))+J11D*(*(gg1+4+1))+J21D*(*(gg1+8+1)))-t*(J02D*(*(gg1+1))+J12D*(*(gg1+4+1))+J22D*(*(gg1+8+1))));
    *(BB+60+10)= 0.25*(1+r)*J02D; 
    *(BB+60+11)= 0.25*(1+r)*J12D;
    *(BB+60+12)= 0.25*(1+r)*J22D;
    *(BB+60+13)= 0.25*(1+r)*((J01D*(*(gg0+2))+J11D*(*(gg0+4+2))+J21D*(*(gg0+8+2)))+t*(J02D*(*(gg0+2))+J12D*(*(gg0+4+2))+J22D*(*(gg0+8+2))));
    *(BB+60+14)= 0.25*(1+r)*((J01D*(*(gg1+2))+J11D*(*(gg1+4+2))+J21D*(*(gg1+8+2)))+t*(J02D*(*(gg1+2))+J12D*(*(gg1+4+2))+J22D*(*(gg1+8+2))));
    *(BB+60+15)= 0.25*(1-r)*J02B;
    *(BB+60+16)= 0.25*(1-r)*J12B;
    *(BB+60+17)= 0.25*(1-r)*J22B;
    *(BB+60+18)= 0.25*(1-r)*((J01B*(*(gg0+3))+J11B*(*(gg0+4+3))+J21B*(*(gg0+8+3)))+t*(J02B*(*(gg0+3))+J12B*(*(gg0+4+3))+J22B*(*(gg0+8+3))));
    *(BB+60+19)= 0.25*(1-r)*((J01B*(*(gg1+3))+J11B*(*(gg1+4+3))+J21B*(*(gg1+8+3)))+t*(J02B*(*(gg1+3))+J12B*(*(gg1+4+3))+J22B*(*(gg1+8+3))));
    
	*(BB+80)   =-0.25*(1-s)*J02C;
    *(BB+80+1) =-0.25*(1-s)*J12C;
    *(BB+80+2) =-0.25*(1-s)*J22C;
    *(BB+80+3) = 0.25*(1-s)*((J00C*(*(gg0))+J10C*(*(gg0+4))+J20C*(*(gg0+8)))-t*(J02C*(*(gg0))+J12C*(*(gg0+4))+J22C*(*(gg0+8))));
    *(BB+80+4) = 0.25*(1-s)*((J00C*(*(gg1))+J10C*(*(gg1+4))+J20C*(*(gg1+8)))-t*(J02C*(*(gg1))+J12C*(*(gg1+4))+J22C*(*(gg1+8))));
    *(BB+80+5) = 0.25*(1-s)*J02C; 
    *(BB+80+6) = 0.25*(1-s)*J12C;
    *(BB+80+7) = 0.25*(1-s)*J22C;
    *(BB+80+8) = 0.25*(1-s)*((J00C*(*(gg0+1))+J10C*(*(gg0+4+1))+J20C*(*(gg0+8+1)))+t*(J02C*(*(gg0+1))+J12C*(*(gg0+4+1))+J22C*(*(gg0+8+1))));
    *(BB+80+9) = 0.25*(1-s)*((J00C*(*(gg1+1))+J10C*(*(gg1+4+1))+J20C*(*(gg1+8+1)))+t*(J02C*(*(gg1+1))+J12C*(*(gg1+4+1))+J22C*(*(gg1+8+1))));
    *(BB+80+10)= 0.25*(1+s)*J02A; 
    *(BB+80+11)= 0.25*(1+s)*J12A;
    *(BB+80+12)= 0.25*(1+s)*J22A;
    *(BB+80+13)= 0.25*(1+s)*((J00A*(*(gg0+2))+J10A*(*(gg0+4+2))+J20A*(*(gg0+8+2)))+t*(J02A*(*(gg0+2))+J12A*(*(gg0+4+2))+J22A*(*(gg0+8+2))));
    *(BB+80+14)= 0.25*(1+s)*((J00A*(*(gg1+2))+J10A*(*(gg1+4+2))+J20A*(*(gg1+8+2)))+t*(J02A*(*(gg1+2))+J12A*(*(gg1+4+2))+J22A*(*(gg1+8+2))));
    *(BB+80+15)=-0.25*(1+s)*J02A; 
    *(BB+80+16)=-0.25*(1+s)*J12A;
    *(BB+80+17)=-0.25*(1+s)*J22A;
    *(BB+80+18)= 0.25*(1+s)*((J00A*(*(gg0+3))+J10A*(*(gg0+4+3))+J20A*(*(gg0+8+3)))-t*(J02A*(*(gg0+3))+J12A*(*(gg0+4+3))+J22A*(*(gg0+8+3))));
    *(BB+80+19)= 0.25*(1+s)*((J00A*(*(gg1+3))+J10A*(*(gg1+4+3))+J20A*(*(gg1+8+3)))-t*(J02A*(*(gg1+3))+J12A*(*(gg1+4+3))+J22A*(*(gg1+8+3))));

//	for (k=0;k<120;k++) { *(Data+k+1) = *(BB+k); }
/*	JI[0,0] -> JI[0], JI[0,1] -> JI[1], JI[0,2] -> JI[2] 
	JI[1,0] -> JI[3], JI[1,1] -> JI[4], JI[1,2] -> JI[5]
	JI[2,0] -> JI[6], JI[2,1] -> JI[7], JI[2,2] -> JI[8]
*/
	td[0][0] = JI[0]*(*(vv))  +JI[1]*(*(vv+3))  +JI[2]*(*(vv+6));
	td[0][1] = JI[0]*(*(vv+1))+JI[1]*(*(vv+3+1))+JI[2]*(*(vv+6+1));
	td[0][2] = JI[0]*(*(vv+2))+JI[1]*(*(vv+3+2))+JI[2]*(*(vv+6+2));
    td[1][0] = JI[3]*(*(vv))  +JI[4]*(*(vv+3))  +JI[5]*(*(vv+6));
	td[1][1] = JI[3]*(*(vv+1))+JI[4]*(*(vv+3+1))+JI[5]*(*(vv+6+1));
	td[1][2] = JI[3]*(*(vv+2))+JI[4]*(*(vv+3+2))+JI[5]*(*(vv+6+2));
	td[2][0] = JI[6]*(*(vv))  +JI[7]*(*(vv+3))  +JI[8]*(*(vv+6));
	td[2][1] = JI[6]*(*(vv+1))+JI[7]*(*(vv+3+1))+JI[8]*(*(vv+6+1));
	td[2][2] = JI[6]*(*(vv+2))+JI[7]*(*(vv+3+2))+JI[8]*(*(vv+6+2));

    *(TD+ 0) =   td[0][0]*td[0][0];
	*(TD+ 0+1) = td[1][0]*td[1][0];
	*(TD+ 0+2) = td[2][0]*td[2][0];
	*(TD+ 0+3) = td[1][0]*td[2][0];
	*(TD+ 0+4) = td[0][0]*td[2][0];
	*(TD+ 0+5) = td[0][0]*td[1][0];

	*(TD+ 6) =   td[0][1]*td[0][1];
	*(TD+ 6+1) = td[1][1]*td[1][1];
	*(TD+ 6+2) = td[2][1]*td[2][1];
	*(TD+ 6+3) = td[1][1]*td[2][1];
	*(TD+ 6+4) = td[0][1]*td[2][1];
	*(TD+ 6+5) = td[0][1]*td[1][1];

    *(TD+12) =   td[0][2]*td[0][2];
	*(TD+12+1) = td[1][2]*td[1][2];
	*(TD+12+2) = td[2][2]*td[2][2];
	*(TD+12+3) = td[1][2]*td[2][2];
	*(TD+12+4) = td[0][2]*td[2][2];
	*(TD+12+5) = td[0][2]*td[1][2];

	*(TD+18) =   2.*td[0][1]*td[0][2];
	*(TD+18+1) = 2.*td[1][1]*td[1][2];
	*(TD+18+2) = 2.*td[2][1]*td[2][2];
	*(TD+18+3) = td[1][1]*td[2][2]+td[2][1]*td[1][2];
	*(TD+18+4) = td[0][1]*td[2][2]+td[0][2]*td[2][1];
	*(TD+18+5) = td[0][1]*td[1][2]+td[1][1]*td[0][2];

	*(TD+24) =   2.*td[0][0]*td[0][2];
	*(TD+24+1) = 2.*td[1][0]*td[1][2];
	*(TD+24+2) = 2.*td[2][0]*td[2][2];
	*(TD+24+3) = td[1][0]*td[2][2]+td[2][0]*td[1][2];
	*(TD+24+4) = td[0][0]*td[2][2]+td[0][2]*td[2][0];
	*(TD+24+5) = td[0][0]*td[1][2]+td[1][0]*td[0][2];

	*(TD+30) =   2.*td[0][0]*td[0][1];
	*(TD+30+1) = 2.*td[1][0]*td[1][1];
	*(TD+30+2) = 2.*td[2][0]*td[2][1];
	*(TD+30+3) = td[1][0]*td[2][1]+td[2][0]*td[1][1];
	*(TD+30+4) = td[0][0]*td[2][1]+td[0][1]*td[2][0];
	*(TD+30+5) = td[0][0]*td[1][1]+td[1][0]*td[0][1];

//	for (i=0;i<dim_Vn0;i++) { 
//		for (j=0;j<dim_Vn1;j++) {
//			k = i*dim_Vn1+j;
//			*(BB+k) = *(Vn+k);
//	}	}
	return rc;
}
/*
                BB[3,0] =-0.25*(1-r)*J02B
                BB[3,1] =-0.25*(1-r)*J12B
                BB[3,2] =-0.25*(1-r)*J22B
                BB[3,3] = 0.25*(1-r)*((J01B*self.gg[0,0,0]+J11B*self.gg[0,1,0]+J21B*self.gg[0,2,0])-t*(J02B*self.gg[0,0,0]+J12B*self.gg[0,1,0]+J22B*self.gg[0,2,0]))
                BB[3,4] = 0.25*(1-r)*((J01B*self.gg[1,0,0]+J11B*self.gg[1,1,0]+J21B*self.gg[1,2,0])-t*(J02B*self.gg[1,0,0]+J12B*self.gg[1,1,0]+J22B*self.gg[1,2,0]))
                BB[3,5] =-0.25*(1+r)*J02D 
                BB[3,6] =-0.25*(1+r)*J12D
                BB[3,7] =-0.25*(1+r)*J22D
                BB[3,8] = 0.25*(1+r)*((J01D*self.gg[0,0,1]+J11D*self.gg[0,1,1]+J21D*self.gg[0,2,1])-t*(J02D*self.gg[0,0,1]+J12D*self.gg[0,1,1]+J22D*self.gg[0,2,1]))
                BB[3,9]= 0.25*(1+r)*((J01D*self.gg[1,0,1]+J11D*self.gg[1,1,1]+J21D*self.gg[1,2,1])-t*(J02D*self.gg[1,0,1]+J12D*self.gg[1,1,1]+J22D*self.gg[1,2,1]))
                BB[3,10]= 0.25*(1+r)*J02D 
                BB[3,11]= 0.25*(1+r)*J12D
                BB[3,12]= 0.25*(1+r)*J22D
                BB[3,13]= 0.25*(1+r)*((J01D*self.gg[0,0,2]+J11D*self.gg[0,1,2]+J21D*self.gg[0,2,2])+t*(J02D*self.gg[0,0,2]+J12D*self.gg[0,1,2]+J22D*self.gg[0,2,2]))
                BB[3,14]= 0.25*(1+r)*((J01D*self.gg[1,0,2]+J11D*self.gg[1,1,2]+J21D*self.gg[1,2,2])+t*(J02D*self.gg[1,0,2]+J12D*self.gg[1,1,2]+J22D*self.gg[1,2,2]))
                BB[3,15]= 0.25*(1-r)*J02B
                BB[3,16]= 0.25*(1-r)*J12B
                BB[3,17]= 0.25*(1-r)*J22B
                BB[3,18]= 0.25*(1-r)*((J01B*self.gg[0,0,3]+J11B*self.gg[0,1,3]+J21B*self.gg[0,2,3])+t*(J02B*self.gg[0,0,3]+J12B*self.gg[0,1,3]+J22B*self.gg[0,2,3]))
                BB[3,19]= 0.25*(1-r)*((J01B*self.gg[1,0,3]+J11B*self.gg[1,1,3]+J21B*self.gg[1,2,3])+t*(J02B*self.gg[1,0,3]+J12B*self.gg[1,1,3]+J22B*self.gg[1,2,3]))
                BB[4,0] =-0.25*(1-s)*J02C
                BB[4,1] =-0.25*(1-s)*J12C
                BB[4,2] =-0.25*(1-s)*J22C
                BB[4,3] = 0.25*(1-s)*((J00C*self.gg[0,0,0]+J10C*self.gg[0,1,0]+J20C*self.gg[0,2,0])-t*(J02C*self.gg[0,0,0]+J12C*self.gg[0,1,0]+J22C*self.gg[0,2,0]))
                BB[4,4] = 0.25*(1-s)*((J00C*self.gg[1,0,0]+J10C*self.gg[1,1,0]+J20C*self.gg[1,2,0])-t*(J02C*self.gg[1,0,0]+J12C*self.gg[1,1,0]+J22C*self.gg[1,2,0]))
                BB[4,5] = 0.25*(1-s)*J02C 
                BB[4,6] = 0.25*(1-s)*J12C
                BB[4,7] = 0.25*(1-s)*J22C
                BB[4,8] = 0.25*(1-s)*((J00C*self.gg[0,0,1]+J10C*self.gg[0,1,1]+J20C*self.gg[0,2,1])+t*(J02C*self.gg[0,0,1]+J12C*self.gg[0,1,1]+J22C*self.gg[0,2,1]))
                BB[4,9] = 0.25*(1-s)*((J00C*self.gg[1,0,1]+J10C*self.gg[1,1,1]+J20C*self.gg[1,2,1])+t*(J02C*self.gg[1,0,1]+J12C*self.gg[1,1,1]+J22C*self.gg[1,2,1]))
                BB[4,10]= 0.25*(1+s)*J02A 
                BB[4,11]= 0.25*(1+s)*J12A 
                BB[4,12]= 0.25*(1+s)*J22A
                BB[4,13]= 0.25*(1+s)*((J00A*self.gg[0,0,2]+J10A*self.gg[0,1,2]+J20A*self.gg[0,2,2])+t*(J02A*self.gg[0,0,2]+J12A*self.gg[0,1,2]+J22A*self.gg[0,2,2]))
                BB[4,14]= 0.25*(1+s)*((J00A*self.gg[1,0,2]+J10A*self.gg[1,1,2]+J20A*self.gg[1,2,2])+t*(J02A*self.gg[1,0,2]+J12A*self.gg[1,1,2]+J22A*self.gg[1,2,2]))
                BB[4,15]=-0.25*(1+s)*J02A 
                BB[4,16]=-0.25*(1+s)*J12A
                BB[4,17]=-0.25*(1+s)*J22A
                BB[4,18]= 0.25*(1+s)*((J00A*self.gg[0,0,3]+J10A*self.gg[0,1,3]+J20A*self.gg[0,2,3])-t*(J02A*self.gg[0,0,3]+J12A*self.gg[0,1,3]+J22A*self.gg[0,2,3]))
                BB[4,19]= 0.25*(1+s)*((J00A*self.gg[1,0,3]+J10A*self.gg[1,1,3]+J20A*self.gg[1,2,3])-t*(J02A*self.gg[1,0,3]+J12A*self.gg[1,1,3]+J22A*self.gg[1,2,3]))

        td=array([[JI[0,0]*vv[0,0]+JI[0,1]*vv[1,0]+JI[0,2]*vv[2,0],JI[0,0]*vv[0,1]+JI[0,1]*vv[1,1]+JI[0,2]*vv[2,1],JI[0,0]*vv[0,2]+JI[0,1]*vv[1,2]+JI[0,2]*vv[2,2]],
                  [JI[1,0]*vv[0,0]+JI[1,1]*vv[1,0]+JI[1,2]*vv[2,0],JI[1,0]*vv[0,1]+JI[1,1]*vv[1,1]+JI[1,2]*vv[2,1],JI[1,0]*vv[0,2]+JI[1,1]*vv[1,2]+JI[1,2]*vv[2,2]],
                  [JI[2,0]*vv[0,0]+JI[2,1]*vv[1,0]+JI[2,2]*vv[2,0],JI[2,0]*vv[0,1]+JI[2,1]*vv[1,1]+JI[2,2]*vv[2,1],JI[2,0]*vv[0,2]+JI[2,1]*vv[1,2]+JI[2,2]*vv[2,2]]])
        TD=array([[td[0,0]**2,         td[1,0]**2,       td[2,0]**2,     td[1,0]*td[2,0],                td[0,0]*td[2,0],                td[0,0]*td[1,0]],
                  [td[0,1]**2,         td[1,1]**2,       td[2,1]**2,     td[1,1]*td[2,1],                td[0,1]*td[2,1],                td[0,1]*td[1,1]],
                  [td[0,2]**2,         td[1,2]**2,       td[2,2]**2,     td[1,2]*td[2,2],                td[0,2]*td[2,2],                td[0,2]*td[1,2]],
                  [2*td[0,1]*td[0,2],2*td[1,1]*td[1,2],2*td[2,1]*td[2,2],td[1,1]*td[2,2]+td[2,1]*td[1,2],td[0,1]*td[2,2]+td[0,2]*td[2,1],td[0,1]*td[1,2]+td[1,1]*td[0,2]],
                  [2*td[0,0]*td[0,2],2*td[1,0]*td[1,2],2*td[2,0]*td[2,2],td[1,0]*td[2,2]+td[2,0]*td[1,2],td[0,0]*td[2,2]+td[0,2]*td[2,0],td[0,0]*td[1,2]+td[1,0]*td[0,2]],
                  [2*td[0,0]*td[0,1],2*td[1,0]*td[1,1],2*td[2,0]*td[2,1],td[1,0]*td[2,1]+td[2,0]*td[1,1],td[0,0]*td[2,1]+td[0,1]*td[2,0],td[0,0]*td[1,1]+td[1,0]*td[0,1]]])

	return BB, det, TD

*/

int SH4FormGeomC( double r, double s, double t, double *GeomK, int dim_GeomK0, int dim_GeomK1, double *Data, int dim_D, double *sig, int dim_sig,
				double *gg0, int dim_g00, int dim_g01, double *gg1, int dim_g10, int dim_g11)
{
	int rc, i, j, k, bi, bj, n2, bi0, bi1, bi2, bi3, bi4, bj0, bj1, bj2, bj3, bj4;
	double N[4], br[4], bs[4], t2, S11, S22, S33, S23si, S23sj, S13rj, S13ri, S12, ggg00, ggg11, ggg01, ggg10;

	N[+0] =  (1-r) *(1-s)*0.25;  // N =  array([(1-r)*(1-s)*0.25, (1+r)*(1-s)*0.25, (1+r)*(1+s)*0.25, (1-r)*(1+s)*0.25]) 
	N[+1] =  (1+r) *(1-s)*0.25;
	N[+2] =  (1+r) *(1+s)*0.25;
	N[+3] =  (1-r) *(1+s)*0.25;
	br[0] =  (-1+s)*0.25;          // br = array([(-1+s)*0.25,      ( 1-s)*0.25,     ( 1+s)*0.25,      -( 1+s)*0.25]) 
	br[1] =  ( 1-s)*0.25;
	br[2] =  ( 1+s)*0.25;
	br[3] = -( 1+s)*0.25;
	bs[0] =  (-1+r)*0.25;         // bs = array([(-1+r)*0.25,     -( 1+r)*0.25,     ( 1+r)*0.25,       ( 1-r)*0.25]) 
	bs[1] = -( 1+r)*0.25;
	bs[2] =  ( 1+r)*0.25;
	bs[3] =  ( 1-r)*0.25;

	t2 = t*t;							//  t2 = t*t
	n2 = 20;							// number of nodes times number of dofs 
	for (i=0;i<4;i++) {					//        for i in xrange(4):
		bi  = 5*i;
		bi0 = n2*  bi;
		bi1 = n2*(bi+1);
		bi2 = n2*(bi+2);
		bi3 = n2*(bi+3);
		bi4 = n2*(bi+4);
		for (j=0;j<4;j++) {				//  for j in xrange(4):
			bj  = 5*j;
			bj0 = bj;
			bj1 = bj+1;
			bj2 = bj+2;
			bj3 = bj+3;
			bj4 = bj+4;
			S11   = sig[0] * br[i]*br[j];    // S11   = sig[0] * br[i]*br[j]              #  s_11 sequence according to voigt notation
            S22   = sig[1] * bs[i]*bs[j];    // S22   = sig[1] * bs[i]*bs[j]              #  s_22
            S33   = sig[2] * N[i] *N[j];    // S33   = sig[2] * N[i] *N[j]               #  s_33
            S23si = sig[3] * bs[i]*N[j];    // S23si = sig[3] * bs[i]*N[j]               #  s_23
            S23sj = sig[3] * N[i] *bs[j];    // S23sj = sig[3] * N[i] *bs[j]              #  s_23
            S13ri = sig[4] * br[i]*N[j];    // S13ri = sig[4] * br[i]*N[j]               #  s_13
            S13rj = sig[4] * N[i] *br[j];    // S13rj = sig[4] * N[i] *br[j]              #  s_13
            S12   = sig[5] *(br[i]*bs[j]+bs[i]*br[j]);    // S12   = sig[5] *(br[i]*bs[j]+bs[i]*br[j]) #  s_12 
//	self.gg[0,0,0] -> (*(gg0+0+0))
//	self.gg[0,0,1] -> (*(gg0+0+1))
//	self.gg[0,0,2] -> (*(gg0+0+2))
//	self.gg[0,0,3] -> (*(gg0+0+3))
//	self.gg[0,1,0] -> (*(gg0+4+0))
//	self.gg[0,1,1] -> (*(gg0+4+1))
//	self.gg[0,1,2] -> (*(gg0+4+2))
//	self.gg[0,1,3] -> (*(gg0+4+3))
//	self.gg[0,2,0] -> (*(gg0+8+0))
//	self.gg[0,2,1] -> (*(gg0+8+1))
//	self.gg[0,2,2] -> (*(gg0+8+2))
//	self.gg[0,2,3] -> (*(gg0+8+3))
//  self.gg[0,0,i] -> (*(gg0+0+i))
//  self.gg[0,0,j] -> (*(gg0+0+j))
//	self.gg[0,1,i] -> (*(gg0+4+i))
//	self.gg[0,1,j] -> (*(gg0+4+j))
//	self.gg[0,2,i] -> (*(gg0+8+i))
//	self.gg[0,2,j] -> (*(gg0+8+j))
            ggg00 = ((*(gg0+0+i))*(*(gg0+0+j)) + (*(gg0+4+i))*(*(gg0+4+j)) + (*(gg0+8+i))*(*(gg0+8+j)));    //ggg00 = (self.gg[0,0,i]*self.gg[0,0,j] + self.gg[0,1,i]*self.gg[0,1,j] + self.gg[0,2,i]*self.gg[0,2,j])
            ggg11 = ((*(gg1+0+i))*(*(gg1+0+j)) + (*(gg1+4+i))*(*(gg1+4+j)) + (*(gg1+8+i))*(*(gg1+8+j)));    //ggg11 = (self.gg[1,0,i]*self.gg[1,0,j] + self.gg[1,1,i]*self.gg[1,1,j] + self.gg[1,2,i]*self.gg[1,2,j])
            ggg01 = ((*(gg0+0+i))*(*(gg1+0+j)) + (*(gg0+4+i))*(*(gg1+4+j)) + (*(gg0+8+i))*(*(gg1+8+j)));    //ggg01 = (self.gg[0,0,i]*self.gg[1,0,j] + self.gg[0,1,i]*self.gg[1,1,j] + self.gg[0,2,i]*self.gg[1,2,j])
            ggg10 = ((*(gg1+0+i))*(*(gg0+0+j)) + (*(gg1+4+i))*(*(gg0+4+j)) + (*(gg1+8+i))*(*(gg0+8+j)));    //ggg10 = (self.gg[1,0,i]*self.gg[0,0,j] + self.gg[1,1,i]*self.gg[0,1,j] + self.gg[1,2,i]*self.gg[0,2,j])

            *(GeomK+bi0+bj0) = S11+S22+S12; //        GeomK[bi,  bj]   = S11+S22+S12
            *(GeomK+bi0+bj1) = 0.; //        GeomK[bi,  bj+1] = 0. 
            *(GeomK+bi0+bj2) = 0.; //        GeomK[bi,  bj+2] = 0.
            *(GeomK+bi0+bj3) = (S11*t + S22*t + S12*t)*(*(gg0+0+j)) + (S13ri + S23si)*(*(gg0+0+j)); //        GeomK[bi,  bj+3] = (S11*t + S22*t + S12*t)*self.gg[0,0,j] + (S13ri + S23si)*self.gg[0,0,j]
            *(GeomK+bi0+bj4) = (S11*t + S22*t + S12*t)*(*(gg1+0+j)) + (S13ri + S23si)*(*(gg1+0+j)); //        GeomK[bi,  bj+4] = (S11*t + S22*t + S12*t)*self.gg[1,0,j] + (S13ri + S23si)*self.gg[1,0,j]

            *(GeomK+bi1+bj0) = 0.; //        GeomK[bi+1,bj]   = 0.
			*(GeomK+bi1+bj1) = S11+S22+S12; //        GeomK[bi+1,bj+1] = S11+S22+S12
            *(GeomK+bi1+bj2) = 0.; //        GeomK[bi+1,bj+2] = 0.
            *(GeomK+bi1+bj3) = (S11*t + S22*t + S12*t)*(*(gg0+4+j)) + (S13ri + S23si)*(*(gg0+4+j)); //        GeomK[bi+1,bj+3] = (S11*t + S22*t + S12*t)*self.gg[0,1,j] + (S13ri + S23si)*self.gg[0,1,j]
            *(GeomK+bi1+bj4) = (S11*t + S22*t + S12*t)*(*(gg1+4+j)) + (S13ri + S23si)*(*(gg1+4+j)); //        GeomK[bi+1,bj+4] = (S11*t + S22*t + S12*t)*self.gg[1,1,j] + (S13ri + S23si)*self.gg[1,1,j]

            *(GeomK+bi2+bj0) = 0.; //        GeomK[bi+2,bj]   = 0.
            *(GeomK+bi2+bj1) = 0.; //        GeomK[bi+2,bj+1] = 0.
			*(GeomK+bi2+bj2) = S11+S22+S12; //        GeomK[bi+2,bj+2] = S11+S22+S12
            *(GeomK+bi2+bj3) = (S11*t + S22*t + S12*t)*(*(gg0+8+j)) + (S13ri + S23si)*(*(gg0+8+j)); //        GeomK[bi+2,bj+3] = (S11*t + S22*t + S12*t)*self.gg[0,2,j] + (S13ri + S23si)*self.gg[0,2,j]
            *(GeomK+bi2+bj4) = (S11*t + S22*t + S12*t)*(*(gg1+8+j)) + (S13ri + S23si)*(*(gg1+8+j)); //        GeomK[bi+2,bj+4] = (S11*t + S22*t + S12*t)*self.gg[1,2,j] + (S13ri + S23si)*self.gg[1,2,j]

            *(GeomK+bi3+bj0) = (S11*t + S22*t + S12*t)*(*(gg0+0+i)) + (S13rj + S23sj)*(*(gg0+0+i)); //        GeomK[bi+3,bj]   = (S11*t + S22*t + S12*t)*self.gg[0,0,i] + (S13rj + S23sj)*self.gg[0,0,i]
            *(GeomK+bi3+bj1) = (S11*t + S22*t + S12*t)*(*(gg0+4+i)) + (S13rj + S23sj)*(*(gg0+4+i)); //        GeomK[bi+3,bj+1] = (S11*t + S22*t + S12*t)*self.gg[0,1,i] + (S13rj + S23sj)*self.gg[0,1,i]
            *(GeomK+bi3+bj2) = (S11*t + S22*t + S12*t)*(*(gg0+8+i)) + (S13rj + S23sj)*(*(gg0+8+i)); //        GeomK[bi+3,bj+2] = (S11*t + S22*t + S12*t)*self.gg[0,2,i] + (S13rj + S23sj)*self.gg[0,2,i]
			*(GeomK+bi3+bj3) = (S11*t2 + S22*t2 + S12*t2 + S33)*ggg00 + (S13ri*t + S13rj*t + S23si*t + S23sj*t)*ggg00; //        GeomK[bi+3,bj+3] = (S11*t2 + S22*t2 + S12*t2 + S33)*ggg00 + (S13ri*t + S13rj*t + S23si*t + S23sj*t)*ggg00
            *(GeomK+bi3+bj4) = (S11*t2 + S22*t2 + S12*t2 + S33)*ggg01 + (S13ri*t + S13rj*t + S23si*t + S23sj*t)*ggg01; //        GeomK[bi+3,bj+4] = (S11*t2 + S22*t2 + S12*t2 + S33)*ggg01 + (S13ri*t + S13rj*t + S23si*t + S23sj*t)*ggg01

            *(GeomK+bi4+bj0) = (S11*t + S22*t + S12*t)*(*(gg1+0+i)) + (S13rj + S23sj)*(*(gg1+0+i)); //        GeomK[bi+4,bj]   = (S11*t + S22*t + S12*t)*self.gg[1,0,i] + (S13rj + S23sj)*self.gg[1,0,i]
            *(GeomK+bi4+bj1) = (S11*t + S22*t + S12*t)*(*(gg1+4+i)) + (S13rj + S23sj)*(*(gg1+4+i)); //        GeomK[bi+4,bj+1] = (S11*t + S22*t + S12*t)*self.gg[1,1,i] + (S13rj + S23sj)*self.gg[1,1,i]
            *(GeomK+bi4+bj2) = (S11*t + S22*t + S12*t)*(*(gg1+8+i)) + (S13rj + S23sj)*(*(gg1+8+i)); //        GeomK[bi+4,bj+2] = (S11*t + S22*t + S12*t)*self.gg[1,2,i] + (S13rj + S23sj)*self.gg[1,2,i]
            *(GeomK+bi4+bj3) = (S11*t2 + S22*t2 + S12*t2 + S33)*ggg10 + (S13ri*t + S13rj*t + S23si*t + S23sj*t)*ggg10; //        GeomK[bi+4,bj+3] = (S11*t2 + S22*t2 + S12*t2 + S33)*ggg10 + (S13ri*t + S13rj*t + S23si*t + S23sj*t)*ggg10
			*(GeomK+bi4+bj4) = (S11*t2 + S22*t2 + S12*t2 + S33)*ggg11 + (S13ri*t + S13rj*t + S23si*t + S23sj*t)*ggg11; //        GeomK[bi+4,bj+4] = (S11*t2 + S22*t2 + S12*t2 + S33)*ggg11 + (S13ri*t + S13rj*t + S23si*t + S23sj*t)*ggg11
		}	
	}
//	*(Data+0) = sig[0];
//	*(Data+1) = sig[1];
//	*(Data+2) = sig[2];
//	*(Data+3) = sig[3];
//	*(Data+4) = sig[4];
//	*(Data+5) = sig[5];
//	*(Data+6) = S11;
//	*(Data+7) = S22;
//	*(Data+8) = S12;
	rc = 0;
	return rc;
}

int BTransSig( double *rvec, int dim_rvec, double fact, double* BB, int dim_BB0, int dim_BB1, double *sig, int dim_sig)
{
	int i, j;
	if (dim_rvec != dim_BB1) {
		return 10;
	}
	if (dim_sig != dim_BB0) {
		return 20;
	}
	for (i = 0; i < dim_BB1; i++) {
		for (j = 0; j < dim_BB0; j++) {
			*(rvec+i) = *(rvec+i) + fact * *(BB + j*dim_BB1 + i) * *(sig+j);
		}
	}
	return 0;
}

int BTxCxB( double *kmat, int dim_kmat0, int dim_kmat1, double fact, double* BB, int dim_BB0, int dim_BB1, double *CC, int dim_CC0, int dim_CC1, double* Data, int dim_D,
			double *X, int dim_X0, int dim_X1)
{
	int i, j, k, k_;
	double XX;
	if (dim_kmat0 != dim_BB1) {
		return 10;
	}
	if (dim_CC0 != dim_BB0) {
		return 20;
	}

//	for (k = 0; k < dim_X0; k++) {
//		for (j = 0; j < dim_X1; j++) {
//			for (k_ = 0; k_ < dim_BB0; k_++) {
//				*(X + k*dim_X1 + j) = *(X + k*dim_X1 + j) + *(CC + k*dim_CC1+ k_) * *(BB + k_*dim_BB1 + j);
//			}
//		}
//	}

	for (j = 0; j < dim_kmat1; j++) {
		for (k = 0; k < dim_BB0; k++) {
			XX = 0.;
			for (k_ = 0; k_ < dim_BB0; k_++) {
				XX = XX + *(CC + k*dim_CC1 + k_) * *(BB + k_*dim_BB1 + j);
			}
			XX = fact * XX;
			for (i = 0; i < dim_kmat0; i++) {
//				*(kmat + i*dim_kmat1 + j) = *(kmat + i*dim_kmat1 + j) + fact * *(BB + k*dim_BB1 + i) * *(X + k*dim_X1 + j);
				*(kmat + i*dim_kmat1 + j) = *(kmat + i*dim_kmat1 + j) + *(BB + k*dim_BB1 + i) * XX; // *(X + k * dim_X1 + j);
			}
		}
	}
	return 0;
}