#define Num_V 1
#define datasize 9
#define end_time 10.
#define PLOTSIZE 100

enum {
     x0,
     a,
     k,
     c,
     d,
     e,
     PARAMSIZE
} PARAM;


extern double Vardata[datasize*Num_V];
extern double timedata[datasize];
extern double Vardataplot[PLOTSIZE*Num_V];
extern double timedataplot[PLOTSIZE];
extern int TOO_MUCH;
extern double lambdaL1;
extern double lambdaL2;

typedef void(*ODE)(double*,double*,double*);
typedef double(*ENTIREFUNC)(double*,double*,double*);


void rk(double* __dt,double* x,double* p,ODE func);
void rkf(double* __dt,double* x,double* p,ODE f);
void testfunc2(double* x,double* p,double* dx);
void gradient(double* p,double* d,double* b,double* grad,int* swch,ENTIREFUNC f);
double norm(int SIZE,double* x);
void vsub(int N,double* u,double* v,double* result);
double dist(int N,double* u,double* v);
double dist1(int N,double* u,double* trueval);
void ODErecord(double* p,double* result,double* timedata,int data_s,ODE f);
double ODEcost(double* p);
double wrap1(double* p,double* d,double* b);
double wrap2(double* p,double* d,double* b);
double shrink(double x,double a);
void controlCG(double* p,double* d,double* b,int* swch,ENTIREFUNC f);
void SplitBregman(double* p,int*swch,ENTIREFUNC f);
void replace(int* lft,int* block,int right,int *rnd,double piv,double* data);
void quicksort(int st,int end,double* data);

