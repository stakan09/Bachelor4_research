extern char* DEtype;

typedef struct LUs_data{
     double* L1;
     double* L2;
     double* L3;
     double* U1;
     double* U2;
     double* y;
}LUmat;

void make_mat(char DE,int bunkatu,LUmat* mat_data);

void __LUbunkai0(int N,double* a,double* b,double* c,LUmat* mat_data);
void __LUbunkai1(int N,double* a,double* b,double* c,LUmat* mat_data);
void __LUbunkai2(int N,double* a,double* b,double* c,LUmat* mat_data);
void LU(char DE,int N,double* a,double* b,double* c,LUmat* mat_data);


void __solveD(int N,LUmat* mat_data,double* u_old,double* u_new);
void __solveN(int N,LUmat* mat_data,double* u_old,double* u_new);
void __solveS(int N,LUmat* mat_data,double* u_old,double* u_new);

void SOLVE(char DE,int N,LUmat* mat_data,double* u_old,double* u_new);

void destroymat(LUmat* mat_data);


