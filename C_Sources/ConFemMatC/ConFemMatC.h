
void PrinCLT_1(double vx, double vy, double vxy, double *pEps);
double PrinCLT_2(double vx, double vy, double vxy);

int ElasticLTC2(int ElemDim, bool ElemPlSt, double* ElemStateVar, int dim_St, double* ElemStateVarN, int dim_StN,
	double* Eps, int dim_E, double* sig, int dim_si, double* MatM, int dim_Ma, double *ww, int dim_ww, double nu, double Emod, double* Dps, int dim_D, double Dt, 
	double rho, double selfwcr_, double selfbw_, double selfepsct, double ElemLch_, double selfCrackVisc, double selffct, double selfepscu,
	double* DataOut, int dim_DataOut);

/* IsoDamage */

int eighC( double I1, double J2, double *Eps, double *la, double *vv);
int EquivStrain1C( int CalcType, double I1, double J2, double J2S, double *Eps, double *EpsD, double *kap_, double *nd, double *Hd, double cc0, double cc1, double cc2, double cc3,
						double *LaMax, double *eVec, double *data);
int eighC2( double I1, double J2D, double *Eps, double *la, double *vv);
int eigJacobiSym(double *Eps, double *la, double *vv);
double ViscExten3DC1( double Dt, double eta, double *Dps, double *ElemStateVar, double *ElemStateVarN, int sI, double *Veps);

int eigJacobiSymWrapper(double *Eps, int dim_E, double *la, int dim_la, double *vv, int dim_vv);
int IsoDamC1( int CalcType, int ElemDim, bool ElemPlSt, bool PrinStrains, double ElemLch, double *ElemStateVar, int dim_St, double *ElemStateVarN, int dim_StN,
					double *Eps, int dim_E, double *sig, int dim_si, double *MatM2, int dim_Ma21, int dim_Ma22, int LiTy,
					double cc0, double cc1, double cc2, double cc3, int RType, double *EpsR, int dim_EpsR, double kapStrength, double ElemCrBws, double gam2, double kapUlt,
					double edt, double ed, double gd, double nu, double Emod, double *Dps, int dim_D, double eta, double RegPar, int ElemScaleType,
					double *sigR, int dim_sigR, double *CR, int dim_CR, double Dt, double sVsTol, double *DataOut, int dim_DataOut);

int IsoDamUniaxC1(int ElemDim, bool ElemPlSt,               double ElemLch, double *ElemStateVar, int dim_St, double *ElemStateVarN, int dim_StN,
					double *Eps, int dim_E, double *sig, int dim_si, double *MatM, int dim_Ma, int LiTy,
					double alpha,double cc1, double cc2, double cc3, int RType, double *EpsR, int dim_EpsR, double kapStrength, double ElemCrBws, double gam2, double kapUlt,
					double edt, double ed, double gd, double nu, double Emod, double *Dps, int dim_D, double eta, double RegPar,
					double *sigR, int dim_sigR, double *CR, int dim_CR, double Dt, double sVsTol, double *DataOut, int dim_DataOut);

/* SDA */
int Rankine(bool PrinStrains, int ElemDim, double *Eps_, double *sig0, double *eVec, double *laMax);
int RankineUpd(double *sig0, double *eVec, double *ElemStateVarN, double *laMax, int off);

/* ElasticLT SDA */

/* Microplane */
void DamFunc1(double kap0V, double alphaV, double betaV, double kapOld, double eta, double *kap, double *dd, double *Dp);
void DevStiffness(double m00, double m01, double m02, double m11, double m12, double m22,
	double n02, double n12, double n22, double tbt, double tbt2, double obt, double obt2,
	double *nn, double *DVD);

int MicroPlaneC1(int ElemDim, bool ElemPlSt, bool PrinStrains, double ElemLch, double *ElemStateVar, int dim_St, double *ElemStateVarN, int dim_StN,
	double *Eps, int dim_E, double *sig, int dim_si, double *MatM, int dim_Ma, int type, double E_V, double E_DM,
	double KV0, double kV1, double kV2, double kap0V, double alphaV, double betaV, int RType, double ElemCrBws, double gam2, double kapUlt, int ElemScaleType,
	double eps_ct, double e0, double ed, double gd, double nu, double Emod, double *Dps, int dim_D, double etaV,
	double Dt, double sVsTol, double *DataOut, int dim_DataOut, int nState, int nInt, int PlStressI, double PlStressL, int iS);

/* modified version - should be more efficient with same results */
int MicroPlaneC2(int ElemDim, bool ElemPlSt, bool PrinStrains, double ElemLch, double *ElemStateVar, int dim_St, double *ElemStateVarN, int dim_StN,
	double *Eps, int dim_E, double *sig, int dim_si, double *MatM, int dim_Ma, int type, double E_V, double E_DM,
	double KV0, double kV1, double kV2, double kap0V, double alphaV, double betaV, int RType, double ElemCrBws, double gam2, double kapUlt, int ElemScaleType,
	double eps_ct, double e0, double ed, double gd, double nu, double Emod, double *Dps, int dim_D, double etaV,
	double Dt, double sVsTol, double *DataOut, int dim_DataOut, int nState, int nInt, int PlStressI, double PlStressL, int iS);

