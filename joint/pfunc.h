#define PARAMS(prefix) prefix##_SCN, prefix##_Phe, prefix##_S2O3


#define PARAMS2(prefix,pf) pf(prefix##_SCN) pf(prefix##_Phe) pf(prefix##_S2O3)


enum{
     PARAMS(X),
     PARAMS(S),
     SHURUI
};

#define P_LIST(pf) \
     PARAMS2(alpha,pf)\
     PARAMS2(X_max,pf)\
     PARAMS2(d_SCN,pf)\
     PARAMS2(d_Phe,pf)\
     PARAMS2(d_S2O3,pf)\
     PARAMS2(c,pf)\
     PARAMS2(beta,pf)\
     PARAMS2(gamma,pf)\
     PARAMS2(A,pf)\
     PARAMS2(B,pf)\
     PARAMS2(P,pf)\
     pf(Pin)

#define P_ENUM_DEF(x) x,
enum{
     P_LIST(P_ENUM_DEF)
     /*
     PARAMS2(alpha,P), 
     PARAMS(X_max), 
     PARAMS(d_SCN),
     PARAMS(d_Phe),
     PARAMS(d_S2O3),
     PARAMS(c),
     PARAMS(P),
     */
     PARAM_N
};

enum{
     PARAMS(initX),
     PARAMS(initS),
};

typedef void(*FUNCPTR)(double*,double*,double*);
typedef double(*FUNCPTR0)(double*,double*);


#ifdef __cplusplus
extern "C"{
#endif

void model1(double* x,double* p,double* dx);
void model2(double* x,double* p,double* dx);
void model3(double* x,double* p,double* dx);
void model4(double* x,double* p,double* dx);
void model5(double* x,double* p,double* dx);
void model6(double* x,double* p,double* dx);
void model7(double* x,double* p,double* dx);
void rk(double dt,double* x,double* p,FUNCPTR func);
double target_func3sp(double*x,double*p);
double target_funcSCN(double*x,double*p);
double target_funcPhe(double*x,double*p);
double target_funcS2O3(double*x,double*p);
double zeroscore(double*x,double*p);
void gradient(double* x,double* p,const int* swch,double* grad,FUNCPTR0 f1,FUNCPTR0 f2,FUNCPTR0 f3,FUNCPTR0 f4);
double norm(double* x,int SIZE);
void controlCG(double* x,double* p,int* swch,FUNCPTR0 f1,FUNCPTR0 f2,FUNCPTR0 f3,FUNCPTR0 f4);
void show_param(const double*p);
void golden_separate(double*x,double* p1,double* p2,double* p3,const int* swch,FUNCPTR0 f1,FUNCPTR0 f2,FUNCPTR0 f3);
extern double speeddata3sp[9][3];
extern double speeddata3spnew[12][3];
extern double speeddataSCN[9][3];
extern double speeddataPhe[9][3];
extern double speeddataS2O3[10][3];
extern double Xdata3sp[9][3];
extern double Xdata3spnew[12][3];
extern double XdataPhe[9][3];
extern double XdataSCN[9][3];
extern double XdataS2O3[10][3];

extern double Lalp;
extern int flag;
#ifdef __cplusplus
}
#endif
