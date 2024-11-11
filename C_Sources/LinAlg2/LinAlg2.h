

void sim0_so(double *r, int dim_r,      double *a, int dim_a, double *b, int dim_b,                         long long *b1, int dim_b1, int *b2, int dim_b2, int n);
void sim0_lu(                           double *a, int dim_a, double *b, int dim_b,                         long long *b1, int dim_b1, int *b2, int dim_b2, int n);
void sim0_mmul(double *le, int dim_le,  double *a, int dim_a, double* b, int dim_b, double *ri, int dim_ri, long long *b1, int dim_b1, int *b2, int dim_b2, int n);
void sim1_so(double *r, int dim_r,      double *a, int dim_a,                                               long long *b1, int dim_b1, int *b2, int dim_b2, int n);
void sim1_lu(                           double *a, int dim_a,                                               long long *b1, int dim_b1, int *b2, int dim_b2, int n);
void sim1_mmul(double *le, int dim_le,  double *a, int dim_a,                       double *ri, int dim_ri, long long *b1, int dim_b1, int *b2, int dim_b2, int n);
void sim1_luP(                          double *a, int dim_a,                                               long long *b1, int dim_b1, int* b2, int dim_b2, int n, int NK);
//void sim1_mmul(double *le, int dim_le,  double *a, int dim_a, double* b, int dim_b, double *ri, int dim_ri, int* b1, int dim_b1, int* b2, int dim_b2, int n);
void BoundFinishLC(int N, int k, long long *SDiag, int dim_SDiag, int *Skyline, int dim_Skyline, double *KVecU, int dim_KVecU, int SymSys, double *KVecL, int dim_KVecL);
void RightHandSideLC(int N, int k, double LF, double LT, double *ub, int dim_ub, double *uub, int dim_uub, long long *SDiag, int dim_SDiag, int *Skyline, int dim_Skyline, \
	double *pp, int dim_pp, double *pp0, int dim_pp0, double *KVecU, int dim_KVecU, int SymSys, double *KVecL, int dim_KVecL);
void RightHandSideLC_(int N, int k,  long long *SDiag, int dim_SDiag, int *Skyline, int dim_Skyline, \
	double *pp, int dim_pp, double *pp0, int dim_pp0, double *KVecU, int dim_KVecU, int SymSys, double *KVecL, int dim_KVecL, double val, double valT);
//double EigenJacobiSymC(double *A, int N, double *DataOut, int dim_DataOut);

