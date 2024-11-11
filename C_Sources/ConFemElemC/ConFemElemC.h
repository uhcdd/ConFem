int CPS4FormBC(double X0, double Y0, double X1, double Y1, double X2, double Y2, double X3, double Y3, double r, double s, double t, double *BC, int dim_BC0, int dim_BC1,
	            double *det_, int dim_det_, double* Data, int dim_D);

int SH4BasicsC( double r, double s, double t, double *N, double *br, double *bs, double *JJ, double *JI, double *vv, double *XX, double *a, double *Vn, double *EdgeDir);
int SH4FormBC( double r, double s, double t, double *BB, int dim_BB0, int dim_BB1, double *Data, int dim_D, double *XX, int dim_X0, int dim_X1, double *a, int dim_a, double *Vn, int dim_Vn0,int dim_Vn1,
				double *EdgeDir, int dim_Ed, double *gg0, int dim_g00, int dim_g01, double *gg1, int dim_g10, int dim_g11, double *TD, int dim_TD0, int dim_TD1);
int SH4FormGeomC( double r, double s, double t, double *GeomK, int dim_GeomK0, int dim_GeomK1, double *Data, int dim_D, double *sig, int dim_sig,
				double *gg0, int dim_g00, int dim_g01, double *gg1, int dim_g10, int dim_g11);

int BTransSig(double* rvec, int dim_rvec, double fact, double* BB, int dim_BB0, int dim_BB1, double* sig, int dim_sig);
int BTxCxB(double* kmat, int dim_kmat0, int dim_kmat1, double fact, double* BB, int dim_BB0, int dim_BB1, double* CC, int dim_CC0, int dim_CC1, double* Data, int dim_D,
				double* X, int dim_X0, int dim_X1);

