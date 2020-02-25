#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <signal.h>
#include <assert.h>
#include <float.h>
#include "pfunc.h"

#ifdef DBG
#define __dbg do{fprintf(stderr,"%s %d %s\n",__FILE__,__LINE__,__PRETTY_FUNCTION__);}while(0)
#else
#define __dbg do{}while(0)
#endif
#ifdef DBG
#define __dbgf(...) do{fprintf(stderr,"%s %d %s " ,__FILE__,__LINE__,__PRETTY_FUNCTION__);fprintf(stderr,__VA_ARGS__);}while(0)
#else
#define __dbgf(...) do{}while(0)
#endif

#define P_NAME_DEF(x) #x,
static const char* PARAM_NAME[PARAM_N+1]={
     P_LIST(P_NAME_DEF)
          "PARAM_N"
};

static double input_HRTper3sp[9]={
     1./24,
     1./18,
     1./12,
     1./8,
     1./10,
     1./24,
     1./24,
     1./24,
     1./24};
static double input_HRTperPhe[9]={
     1./24,
     1./12,
     1./8,
     1./8,
     1./8,
     1./8,
     1./8,
     1./8,
     1./8};
static double input_HRTperSCN[9]={
     1./24,
     1./18,
     1./12,
     1./8,
     1./6,
     1./6,
     1./6,
     1./6,
     0
};
static double input_HRTperS2O3[10]={
     1./24,
     1./12,
     1./8,
     1./8,
     1./8,
     1./8,
     1./8,
     1./8,
     1./8,
     1./8
};

double Lalp = 1;
int flag = 0;
static double source_data3sp[9][3]={{225.5, 97.8, 232.95}, {280.7, 119.6, 310.2}, {418.1, 189.9, 451.1},
     {601.4, 251.9, 634.5}, {459.9, 225.5, 500.5}, {210.3, 83.45, 219.2}, {382.5, 36.6, 67.5},
     {59.7, 251.0, 73.2}, {75.9, 37.5, 452.4}};

static double source_data3spnew[9][3]={{228.65, 124.133, 243.95}, {285.956,170.667, 323.493}, {448.12,247.08,487.12},
     {660.47,384,745.455}, {498.32,307.2,584.504}, {230,128,248.48}, {422.933, 42.667,71.81},
     {56.67,256,76.89}, {75.23,42.667, 494.06}};


static double source_dataPhe[9][3]={{0,128,0},{0,256,0},{0,384,0},{0,512,0},{0,768,0},{0,1536,0},{0,2048,0},{0,3072,0},{0,3737,0}};

static double source_dataSCN[9][3]={{199.6,0,0},{303.7,0,0},{449,0,0},{682,0,0},{894,0,0},{1359,0,0},{1125,0,0},{881.6,0,0},{0,0,0}};

static double source_dataS2O3[10][3]={{0,0,250.1},{0,0,490.6},{0,0,745.05},{0,0,1487.4},{0,0,3009.75},{0,0,6042.9},{0,0,7499.4},{0,0,12087.6},{0,0,18118.8},{0,0,9269.55}};

static double init3sp[6]={1,1,1,230,128,252};
static double initSCN[6]={1,0,0,230,0,0};
static double initPhe[6]={0,1,0,0,128,0};
static double initS2O3[6]={0,0,1,0,0,252};
static double init[6]={0};

double speeddataSCN[9][3]={0};
double speeddataS2O3[10][3]={0};
double speeddataPhe[9][3]={0};
double speeddata3sp[9][3]={0};
double speeddata3spnew[12][3]={0};
double Xdata3sp[9][3]={0};
double Xdata3spnew[12][3]={0};
double XdataPhe[9][3]={0};
double XdataSCN[9][3]={0};
double XdataS2O3[10][3]={0};

void model1(double* x,double* p,double* dx){
     dx[X_SCN]  = p[alpha_SCN] *x[S_SCN] *(p[X_max_SCN]  - x[X_SCN] )*x[X_SCN]  - (p[d_SCN_SCN]           + p[d_SCN_S2O3]*x[X_S2O3] + p[d_SCN_Phe] *x[X_Phe])*x[X_SCN] ;
     dx[X_Phe]  = p[alpha_Phe] *x[S_Phe] *(p[X_max_Phe]  - x[X_Phe] )*x[X_Phe]  - (p[d_Phe_SCN]*x[X_SCN]  + p[d_Phe_S2O3]*x[X_S2O3] + p[d_Phe_Phe]          )*x[X_Phe] ;
     dx[X_S2O3] = p[alpha_S2O3]*x[S_S2O3]*(p[X_max_S2O3] - x[X_S2O3])*x[X_S2O3] - (p[d_S2O3_SCN]*x[X_SCN] + p[d_S2O3_S2O3]          + p[d_S2O3_Phe]*x[X_Phe])*x[X_S2O3];
     dx[S_SCN]  = -p[c_SCN] *x[X_SCN] *x[S_SCN]  + p[Pin]*(p[P_SCN]*init[initS_SCN]  - x[S_SCN]) ;
     dx[S_Phe]  = -p[c_Phe] *x[X_Phe] *x[S_Phe]  + p[Pin]*(p[P_Phe]*init[initS_Phe]  - x[S_Phe]) ;
     dx[S_S2O3] = -p[c_S2O3]*x[X_S2O3]*x[S_S2O3] + p[Pin]*(p[P_S2O3]*init[initS_S2O3] - x[S_S2O3]);
} 
void model2(double* x,double* p,double* dx){
     dx[X_SCN]  = p[alpha_SCN] *pow(x[S_SCN],1)/(pow(p[beta_SCN],1)+pow(x[S_SCN],1)) *(p[X_max_SCN]  - x[X_SCN] )*x[X_SCN]  - (p[d_SCN_SCN]           + p[d_SCN_S2O3]*x[X_S2O3] + p[d_SCN_Phe] *x[X_Phe])*x[X_SCN] ;
     dx[X_Phe]  = p[alpha_Phe] *pow(x[S_Phe],1)/(pow(p[beta_Phe],1)+pow(x[S_Phe],1))*(p[X_max_Phe]  - x[X_Phe] )*x[X_Phe]  - (p[d_Phe_SCN]*x[X_SCN]  + p[d_Phe_S2O3]*x[X_S2O3] + p[d_Phe_Phe]          )*x[X_Phe] ;
     dx[X_S2O3] = p[alpha_S2O3]*pow(x[S_S2O3],1)/(pow(p[beta_S2O3],1)+pow(x[S_S2O3],1))*(p[X_max_S2O3] - x[X_S2O3])*x[X_S2O3] - (p[d_S2O3_SCN]*x[X_SCN] + p[d_S2O3_S2O3]          + p[d_S2O3_Phe]*x[X_Phe])*x[X_S2O3];
     dx[S_SCN]  = -p[c_SCN] *pow(x[X_SCN],1) /(pow(p[gamma_SCN],1)+pow(x[X_SCN],1))*x[S_SCN]  + p[Pin]*(p[P_SCN]*init[initS_SCN]  - x[S_SCN]) ;
     dx[S_Phe]  = -p[c_Phe] *pow(x[X_Phe],1) /(pow(p[gamma_Phe],1)+pow(x[X_Phe],1))*x[S_Phe]  + p[Pin]*(p[P_Phe]*init[initS_Phe]  - x[S_Phe]) ;
     dx[S_S2O3] = -p[c_S2O3]*pow(x[X_S2O3],1)/(pow(p[gamma_S2O3],1)+pow(x[X_S2O3],1))*x[S_S2O3] + p[Pin]*(p[P_S2O3]*init[initS_S2O3] - x[S_S2O3]);
} 
void model3(double* x,double* p,double* dx){
     dx[X_SCN]  = p[alpha_SCN] *pow(x[S_SCN],1)/(pow(p[beta_SCN],1)+pow(x[S_SCN],1))*(p[X_max_SCN]  - x[X_SCN] )*x[X_SCN]  - (p[d_SCN_Phe]*x[X_Phe]/*/(1.+x[X_Phe])*/  + p[d_SCN_S2O3]*x[X_S2O3]/*/(1.+x[X_S2O3])*/ + p[d_SCN_SCN] + 
               p[A_SCN]*pow(x[S_SCN],3)/(pow(p[B_SCN],3)+pow(x[S_SCN],3))/*+p[A_Phe]*pow(x[S_Phe],3)/(pow(p[B_Phe],3)+pow(x[S_Phe],3))+p[A_S2O3]*pow(x[S_S2O3],3)/(pow(p[B_S2O3],3)+pow(x[S_S2O3],3))*/)*x[X_SCN] ;
     dx[X_Phe]  = p[alpha_Phe] *pow(x[S_Phe],1)/(pow(p[beta_Phe],1)+pow(x[S_Phe],1))*(p[X_max_Phe]  - x[X_Phe] )*x[X_Phe]  - (p[d_Phe_SCN]*x[X_SCN]/*/(1.+x[X_SCN])*/  + p[d_Phe_S2O3]*x[X_S2O3]/*/(1.+x[X_S2O3])*/ + p[d_Phe_Phe] + 
               p[A_Phe]*pow(x[S_Phe],3)/(pow(p[B_Phe],3)+pow(x[S_Phe],3))/*+p[A_SCN]*pow(x[S_SCN],3)/(pow(p[B_SCN],3)+pow(x[S_SCN],3))+p[A_S2O3]*pow(x[S_S2O3],3)/(pow(p[B_S2O3],3)+pow(x[S_S2O3],3))*/)*x[X_Phe] ;
     dx[X_S2O3] = p[alpha_S2O3] *pow(x[S_S2O3],1)/(pow(p[beta_S2O3],1)+pow(x[S_S2O3],1))*(p[X_max_S2O3]  - x[X_S2O3] )*x[X_S2O3]  - (p[d_S2O3_SCN]*x[X_SCN]/*/(1.+x[X_SCN])*/  + p[d_S2O3_Phe]*x[X_Phe]/*/(1.+x[X_Phe])*/ + p[d_S2O3_S2O3] + 
               p[A_S2O3]*pow(x[S_S2O3],3)/(pow(p[B_S2O3],3)+pow(x[S_S2O3],3))/*+p[A_SCN]*pow(x[S_SCN],3)/(pow(p[B_SCN],3)+pow(x[S_SCN],3))+p[A_Phe]*pow(x[S_Phe],3)/(pow(p[B_Phe],3)+pow(x[S_Phe],3))*/)*x[X_S2O3];
     dx[S_SCN]  = -p[c_SCN] *pow(x[X_SCN],1) /(pow(p[gamma_SCN],1)+pow(x[X_SCN],1))*x[S_SCN]  + p[Pin]*(p[P_SCN]*init[initS_SCN]  - x[S_SCN]) ;
     dx[S_Phe]  = -p[c_Phe] *pow(x[X_Phe],1) /(pow(p[gamma_Phe],1)+pow(x[X_Phe],1))*x[S_Phe]  + p[Pin]*(p[P_Phe]*init[initS_Phe]  - x[S_Phe]) ;
     dx[S_S2O3] = -p[c_S2O3]*pow(x[X_S2O3],1)/(pow(p[gamma_S2O3],1)+pow(x[X_S2O3],1))*x[S_S2O3] + p[Pin]*(p[P_S2O3]*init[initS_S2O3] - x[S_S2O3]);
}
void model4(double* x,double* p,double* dx){
     dx[X_SCN]  = p[alpha_SCN] *pow(x[S_SCN],1)/(pow(p[beta_SCN],1)+pow(x[S_SCN],1))*(p[X_max_SCN]  - x[X_SCN] - p[d_SCN_Phe]*x[X_Phe]  - p[d_SCN_S2O3]*x[X_S2O3])*x[X_SCN] - p[A_SCN]*pow(x[S_SCN],3)/(pow(p[B_SCN],3)+pow(x[S_SCN],3))*x[X_SCN] - p[d_SCN_SCN]*x[X_SCN] ;
     dx[X_Phe]  = p[alpha_Phe] *pow(x[S_Phe],1)/(pow(p[beta_Phe],1)+pow(x[S_Phe],1))*(p[X_max_Phe]  - x[X_Phe] - p[d_Phe_SCN]*x[X_SCN]  - p[d_Phe_S2O3]*x[X_S2O3])*x[X_Phe] - p[A_Phe]*pow(x[S_Phe],3)/(pow(p[B_Phe],3)+pow(x[S_Phe],3))*x[X_Phe] - p[d_Phe_Phe]*x[X_Phe];
     dx[X_S2O3] = p[alpha_S2O3] *pow(x[S_S2O3],1)/(pow(p[beta_S2O3],1)+pow(x[S_S2O3],1))*(p[X_max_S2O3]  - x[X_S2O3] - p[d_S2O3_SCN]*x[X_SCN]  - p[d_S2O3_Phe]*x[X_Phe])*x[X_S2O3] - p[A_S2O3]*pow(x[S_S2O3],3)/(pow(p[B_S2O3],3)+pow(x[S_S2O3],3))*x[X_S2O3] - p[d_S2O3_S2O3]*x[X_S2O3] ;
     dx[S_SCN]  = -p[c_SCN] *pow(x[X_SCN],1) /(pow(p[gamma_SCN],1)+pow(x[X_SCN],1))*x[S_SCN]  + p[Pin]*(p[P_SCN]*init[initS_SCN]  - x[S_SCN]) ;
     dx[S_Phe]  = -p[c_Phe] *pow(x[X_Phe],1) /(pow(p[gamma_Phe],1)+pow(x[X_Phe],1))*x[S_Phe]  + p[Pin]*(p[P_Phe]*init[initS_Phe]  - x[S_Phe]) ;
     dx[S_S2O3] = -p[c_S2O3]*pow(x[X_S2O3],1)/(pow(p[gamma_S2O3],1)+pow(x[X_S2O3],1))*x[S_S2O3] + p[Pin]*(p[P_S2O3]*init[initS_S2O3] - x[S_S2O3]);
}
void model5(double* x,double* p,double* dx){
     dx[X_SCN]  = p[alpha_SCN] *(pow(x[S_SCN],1)/(pow(p[beta_SCN],1)+pow(x[S_SCN],1))*p[X_max_SCN]  - x[X_SCN] )*x[X_SCN]  - (p[d_SCN_Phe]*x[X_Phe]  + p[d_SCN_S2O3]*x[X_S2O3] + p[d_SCN_SCN] + p[A_SCN]*pow(x[S_SCN],3)/(pow(p[B_SCN],3)+pow(x[S_SCN],3)))*x[X_SCN] ;
     dx[X_Phe]  = p[alpha_Phe] *(pow(x[S_Phe],1)/(pow(p[beta_Phe],1)+pow(x[S_Phe],1))*p[X_max_Phe]  - x[X_Phe] )*x[X_Phe]  - (p[d_Phe_SCN]*x[X_SCN]  + p[d_Phe_S2O3]*x[X_S2O3] + p[d_Phe_Phe] + p[A_Phe]*pow(x[S_Phe],3)/(pow(p[B_Phe],3)+pow(x[S_Phe],3)))*x[X_Phe] ;
     dx[X_S2O3] = p[alpha_S2O3] *(pow(x[S_S2O3],1)/(pow(p[beta_S2O3],1)+pow(x[S_S2O3],1))*p[X_max_S2O3]  - x[X_S2O3] )*x[X_S2O3]  - (p[d_S2O3_SCN]*x[X_SCN]  + p[d_S2O3_Phe]*x[X_Phe] + p[d_S2O3_S2O3] + p[A_S2O3]*pow(x[S_S2O3],3)/(pow(p[B_S2O3],3)+pow(x[S_S2O3],3)))*x[X_S2O3] ;
     dx[S_SCN]  = -p[c_SCN] *pow(x[X_SCN],1) /(pow(p[gamma_SCN],1)+pow(x[X_SCN],1))*x[S_SCN]  + p[Pin]*(p[P_SCN]*init[initS_SCN]  - x[S_SCN]) ;
     dx[S_Phe]  = -p[c_Phe] *pow(x[X_Phe],1) /(pow(p[gamma_Phe],1)+pow(x[X_Phe],1))*x[S_Phe]  + p[Pin]*(p[P_Phe]*init[initS_Phe]  - x[S_Phe]) ;
     dx[S_S2O3] = -p[c_S2O3]*pow(x[X_S2O3],1)/(pow(p[gamma_S2O3],1)+pow(x[X_S2O3],1))*x[S_S2O3] + p[Pin]*(p[P_S2O3]*init[initS_S2O3] - x[S_S2O3]);
}

void rk(double dt,double* x,double* p,FUNCPTR func){
     double dth=1./2*dt;
     double tempx[SHURUI]={},x1[SHURUI]={},x2[SHURUI]={},x3[SHURUI]={},x4[SHURUI]={};
     for(int i=0;i<SHURUI;i++)tempx[i]=x[i];
     func(tempx,p,x1);
     for(int i=0;i<SHURUI;i++) tempx[i]=x[i]+dth*x1[i];
     func(tempx,p,x2);
     for(int i=0;i<SHURUI;i++) tempx[i]=x[i]+dth*x2[i];
     func(tempx,p,x3);
     for(int i=0;i<SHURUI;i++) tempx[i]=x[i]+dt*x3[i];
     func(tempx,p,x4);
     for(int i=0;i<SHURUI;i++) x[i]=x[i]+dt*(x1[i]+2*x2[i]+2*x3[i]+x4[i])/6.0;
}

void rkf(double *dt, double *x0, double *x1,double *p,FUNCPTR f){
     double k1[SHURUI]={}, k2[SHURUI]={}, k3[SHURUI]={}, k4[SHURUI]={}, k5[SHURUI]={}, k6[SHURUI]={}, tmp1[SHURUI]={};
     int loop = 1;
     while(loop){
          f(x0, p, k1);
          for(int i = 0; i < SHURUI; ++i)
               tmp1[i] = x0[i] + 1./4**dt*k1[i];
          f(tmp1, p , k2);
          for(int i = 0; i < SHURUI; ++i)
               tmp1[i] = x0[i] + *dt*(3./32*k1[i] + 9./32*k2[i]);
          f(tmp1, p, k3);
          for(int i = 0; i < SHURUI; ++i)
               tmp1[i] = x0[i] + *dt*(1932./2197*k1[i] - 7200./2197*k2[i] + 7296./2197*k3[i]);
          f(tmp1, p, k4);
          for(int i = 0; i < SHURUI; ++i)
               tmp1[i] = x0[i] + *dt*(439./216*k1[i] - 8*k2[i] + 3680./513*k3[i] - 845./4104*k4[i]);
          f(tmp1, p, k5);
          for(int i = 0; i < SHURUI; ++i)
               tmp1[i] = x0[i] + *dt*(-8./27*k1[i] + 2*k2[i] - 3544./2565*k3[i] + 1859./4104*k4[i] - 11./40*k5[i]);
          f(tmp1, p, k6);
          for(int i = 0; i < SHURUI; ++i){
               tmp1[i] = *dt*(1./360*k1[i] - 128./4275*k3[i] - 2197./75240*k4[i] + 1./50*k5[i] + 2./55*k6[i]);
               //     tmp1[i] = x0[i] + dt*(16./135*k1[i] + 6656./12825*k3[i] + 28561./56430*k4[i] - 9./50*k5[i] + 2./55*k6[i]);
               //     tmp2[i] = x0[i] + dt*(25./216*k1[i] + 1408./2565*k3[i] + 2197./4104*k4[i] - 1./5*k5[i]);
          }
          double err = 0;
          for(int i = 0; i < SHURUI; ++i)err += tmp1[i]*tmp1[i];
          err = sqrt(err);
          const double tol = 1e-10;
          if(err >tol){
               *dt *= 0.5;
               //fprintf(stderr, "tol = %e, err = %e, restart, new dt = %f\n", tol, err, *dt);
          }else *dt *= 1.1, loop = 0;

     }
     for(int i = 0; i < SHURUI; ++i)
          x1[i] = x0[i] + *dt*(16./135*k1[i] + 6656./12825*k3[i] + 28561./56430*k4[i] - 9./50*k5[i] + 2./55*k6[i]);
}
static double input_concentrationS2O3[11] = {0,1,1,1,8./4,8./2,8./1,8./0.8,8./0.5,8./0.3,8./0.25};
static double input_concentrationSCN[10] = {0,1,1,1,1,1,6./4,6./2,6./3,0};//SCN
static double input_concentrationPhe[10] = {0,1,1,1,8./6,8./4,8./2,8./1.5,8./1,8./0.75};//PHE
static double input_concentration3sp[3] = {1,2,1./3};
static double input_concentration3spnew[12] = {1,2,6,12,18,24,30,40,50,60,80,90};

double target_func3sp(double*x,double*p){
     int condition = 9,count=1;
     double difference=0;
     double dt=0.001;
     double gosa = 0;
     double gosaold = 0;
     double eps=1e-2*dt;
     double speed[3]={};
     double xtemp[SHURUI]={};
     __dbg;
     

     for(int i=0;i<SHURUI;i++) init[i]=init3sp[i];
     for(int i=0;i<SHURUI;i++) x[i]=xtemp[i]=init[i];

     for(int j=0;j<condition;j++){
          
        /*  p[Pin]=1./24;
           p[P_SCN]=input_concentration3spnew[j];
               p[P_Phe]=input_concentration3spnew[j];
               p[P_S2O3]=input_concentration3spnew[j];
*/
          p[Pin]=input_HRTper3sp[j];
          if(j<6){
               p[P_SCN]=input_concentration3sp[0];
               p[P_Phe]=input_concentration3sp[0];
               p[P_S2O3]=input_concentration3sp[0];
          }else if(j==6){
               p[P_SCN]=input_concentration3sp[1];
               p[P_Phe]=input_concentration3sp[2];
               p[P_S2O3]=input_concentration3sp[2];
          }else if(j==7){
               p[P_SCN]=input_concentration3sp[2];
               p[P_Phe]=input_concentration3sp[1];
               p[P_S2O3]=input_concentration3sp[2];
          }else if(j==8){
               p[P_SCN]=input_concentration3sp[2];
               p[P_Phe]=input_concentration3sp[2];
               p[P_S2O3]=input_concentration3sp[1];
          }
          
          __dbg;
          do{
               __dbg;
               count+=1;
               rk(dt,xtemp,p,model3);
               __dbg;
               gosa = 0;
               for(int i=0;i<SHURUI;i++) gosa+=fabs(xtemp[i]-x[i])/fmax(fabs(x[i]),fabs(xtemp[i]));
               __dbgf("gosa:%.15e\n",gosa);
               for(int i=0;i<SHURUI;i++) x[i]=xtemp[i];
               __dbg;
               if(count==1e10){
                    flag=1;
                    break;
               }
               /*   if(fabs(gosaold-gosa)/fmax(gosaold,gosa)<1e-10){
                    break;
                    }*/
               gosaold = gosa;
          }while(gosa>eps);
          __dbg;
          for(int i=0;i<3;i++) speed[i]=speeddata3sp[j][i]=(init[initS_SCN+i]*p[P_SCN+i]-x[S_SCN+i])*24.*p[Pin];
          for(int i=0;i<3;i++) Xdata3sp[j][i]=x[i];
         difference = difference
               /*+pow(fabs(speed[0]-source_data3spnew[j][0]),2)
                 +pow(fabs(speed[1]-source_data3spnew[j][1]),2)
                 +pow(fabs(speed[2]-source_data3spnew[j][2]),2)*/
               +pow(100*fabs(speed[0]-source_data3sp[j][0])/source_data3sp[j][0],2)
               +pow(100*fabs(speed[1]-source_data3sp[j][1])/source_data3sp[j][1],2)
               +pow(100*fabs(speed[2]-source_data3sp[j][2])/source_data3sp[j][2],2)
               ;

     }

     __dbg;
     if(flag==1){
          difference=1e15;
     }

     difference+=Lalp*pow(norm(p,PARAM_N-4),2);
     return difference;
}

double target_funcSCN(double*x,double*p){
     int condition = 7,count=1;
     double difference=0;
     double dt=0.001;
     double gosa = 0;
     double gosaold = 0;
     double eps=1e-2*dt;
     double speed[3]={};
     double xtemp[SHURUI]={};
     __dbg;
     for(int i=0;i<SHURUI;i++) init[i]=initSCN[i];
     for(int i=0;i<SHURUI;i++) x[i]=xtemp[i]=init[i];
     for(int j=0;j<condition;j++){



          p[Pin]=input_HRTperSCN[j];
          p[P_SCN]=input_concentrationSCN[j+1];
          p[P_Phe]=input_concentrationSCN[0];
          p[P_S2O3]=input_concentrationSCN[0];
          __dbg;
          do{
               __dbg;
               count+=1;
               rk(dt,xtemp,p,model3);
               __dbg;
               gosa = 0;
               //   for(int i=0;i<4;i++) gosa+=fabs(xtemp[i]-x[i])/fmax(fabs(x[i]),fabs(xtemp[i]));
               __dbgf("gosa:%.15e\n",gosa);
               gosa=fabs(xtemp[0]-x[0])/fmax(fabs(x[0]),fabs(xtemp[0]))+fabs(xtemp[3]-x[3])/fmax(fabs(x[3]),fabs(xtemp[3]));
               for(int i=0;i<SHURUI;i++) x[i]=xtemp[i];
               __dbg;
               if(count==1e7){
                    flag=1;
                    break;
               }
               /* if(fabs(gosaold-gosa)/fmax(gosaold,gosa)<1e-10){
                  break;
                  }*/
               //gosaold = gosa;
          }while(gosa>eps);
          __dbg;
          speed[0]=speeddataSCN[j][0]=(init[initS_SCN]*p[P_SCN]-x[S_SCN])*24.*p[Pin];
          for(int i=0;i<3;i++) XdataSCN[j][i]=x[i];
          difference = difference
               //+pow(fabs(speed[0]-source_dataSCN[j][0]),2)
               +pow(100*fabs(speed[0]-source_dataSCN[j][0])/source_dataSCN[j][0],2)
               ;
     }

     if(flag==1){
          difference=1e15;
          flag=0;
     }

     __dbg;
     difference+=Lalp*pow(norm(p,PARAM_N-4),2);
     return difference;
}

double zeroscore(double* x,double* p){
     double noscore=0.0;
     return noscore;
}

double target_funcPhe(double*x,double*p){
     int condition = 9,count=1;
     double difference=0;
     double dt=0.001;
     double gosa = 0;
     double gosaold = 0;
     double eps=1e-2*dt;
     double speed[3]={};
     double xtemp[SHURUI]={};
     __dbg;
     for(int i=0;i<SHURUI;i++) init[i]=initPhe[i];
     for(int i=0;i<SHURUI;i++) x[i]=xtemp[i]=init[i];
     for(int j=0;j<condition;j++){
          p[Pin]=input_HRTperPhe[j];
          p[P_SCN]=input_concentrationPhe[0];
          p[P_Phe]=input_concentrationPhe[j+1];
          p[P_S2O3]=input_concentrationPhe[0];

          __dbg;
          do{
               __dbg;
               count+=1;
               rk(dt,xtemp,p,model3);
               __dbg;
               gosa = 0;
               //for(int i=0;i<3;i++) gosa+=fabs(xtemp[i]-x[i])/fmax(fabs(x[i]),fabs(xtemp[i]));
               gosa=fabs(xtemp[1]-x[1])/fmax(fabs(x[1]),fabs(xtemp[1]))+fabs(xtemp[4]-x[4])/fmax(fabs(x[4]),fabs(xtemp[4]));
               __dbgf("gosa:%.15e\n",gosa);
               for(int i=0;i<SHURUI;i++) x[i]=xtemp[i];
               __dbg;
               if(count==1e10){
                    flag=1;
                    break;
               }
               /* if(fabs(gosaold-gosa)/fmax(gosaold,gosa)<1e-10){
                  break;
                  }*/
               gosaold=gosa;
          }while(gosa>eps);
          __dbg;
          for(int i=0;i<3;i++) XdataPhe[j][i]=x[i];
          speed[1]=speeddataPhe[j][1]=(init[initS_Phe]*p[P_Phe]-x[S_Phe])*24.*p[Pin];
          difference = difference
               +pow(100*fabs(speed[1]-source_dataPhe[j][1])/source_dataPhe[j][1],2)
               // +pow(fabs(speed[1]-source_dataPhe[j][1]),2)
               ;
     }

     __dbg;
     if(flag==1){
          difference=1e15;
          flag=0;
     }
     difference+=Lalp*pow(norm(p,PARAM_N-4),2);
     return difference;
}

double target_funcS2O3(double*x,double*p){
     int condition = 10,count=1;
     double difference=0;
     double dt=0.001;
     double gosa = 0;
     double gosaold = 0;
     double eps=1e-2*dt;
     double speed[3]={};
     double xtemp[SHURUI]={};
     // char filename[128];
     // FILE *fp;
     __dbg;
     for(int i=0;i<SHURUI;i++) init[i]=initS2O3[i];
     for(int i=0;i<SHURUI;i++) x[i]=xtemp[i]=init[i];

     for(int j=0;j<condition;j++){

          p[Pin]=input_HRTperS2O3[j];
          p[P_SCN]=input_concentrationS2O3[0];
          p[P_Phe]=input_concentrationS2O3[0];
          p[P_S2O3]=input_concentrationS2O3[j+1];
          __dbg;
          //sprintf(filename,"filename%d",j);
          //fp=fopen(filename,"w");  
          do{  

               //    fprintf(fp,"%lf %lf %lf %lf %lf %lf\n",x[0],x[1],x[2],x[3],x[4],x[5]);
               __dbg;
               count+=1;
               rk(dt,xtemp,p,model3);
               __dbg;
               gosa = 0;
               //for(int i=0;i<3;i++) gosa+=fabs(xtemp[i]-x[i])/fmax(fabs(x[i]),fabs(xtemp[i]));
               gosa=fabs(xtemp[2]-x[2])/fmax(fabs(x[2]),fabs(xtemp[2]))+fabs(xtemp[5]-x[5])/fmax(fabs(x[5]),fabs(xtemp[5]));
               __dbgf("gosa:%.15e\n",gosa);
               for(int i=0;i<SHURUI;i++) x[i]=xtemp[i];
               __dbg;
               if(count==1e10){
                    flag=1;
                    break;
               }
               /* if(fabs(gosaold-gosa)/fmax(gosaold,gosa)<1e-10){
                  break;
                  }*/
               gosaold = gosa;
          }while(gosa>eps);
          // fclose(fp);
          __dbg;
          speed[2]=speeddataS2O3[j][2]=(init[initS_S2O3]*p[P_S2O3]-x[S_S2O3])*24.*p[Pin];
          for(int i=0;i<3;i++) XdataS2O3[j][i]=x[i];
          difference = difference
               //  +pow(fabs(speed[2]-source_dataS2O3[j][2]),2)
               +pow(100*fabs(speed[2]-source_dataS2O3[j][2])/source_dataS2O3[j][2],2)
               ;
     }

     if(flag==1){
          difference=1e15;
          flag=0;
     }

     __dbg;
     difference+=Lalp*pow(norm(p,PARAM_N-4),2);
     return difference;
}



void gradient(double* x,double* p,const int* swch,double* grad,FUNCPTR0 f1,FUNCPTR0 f2,FUNCPTR0 f3,FUNCPTR0 f4){
     double delta=1e-10;
     double rel_tol=1E-7;
     double abs_tol=1E-7;
     double meps=DBL_EPSILON;
     double msqrteps=sqrt(PARAM_N*meps);
     double backward=0,forward=0;
     double temppara1[PARAM_N],temppara2[PARAM_N];
     for(int i=0;i<PARAM_N;i++){ 
          temppara1[i]=p[i];
          temppara2[i]=p[i];
     }
     for(int i=0;i<PARAM_N;i++){
          double reps=delta;//fmax(msqrteps,fabs(p[i])*rel_tol + abs_tol);
          if(swch[i]==0){//kaeru wo 0

               temppara1[i]=temppara1[i]-fabs(temppara1[i])*reps;
               backward = f1(x,temppara1)+f2(x,temppara1)+f3(x,temppara1)+f4(x,temppara1);
               temppara2[i]=temppara2[i]+fabs(temppara2[i])*reps;
               forward  = f1(x,temppara2)+f2(x,temppara2)+f3(x,temppara2)+f4(x,temppara2);
               grad[i]=(forward-backward)/(2.*delta*p[i]);
          }else{
               grad[i]=0;
          }
          temppara1[i]=p[i];
          temppara2[i]=p[i];
     }
}
double norm(double* x,int SIZE){
     double sum=0.0;
     for(int i=0;i<SIZE;i++) sum+=pow(x[i],2);
     return sqrt(sum);
}

void show_param(const double*p){
     FILE*fp;
     fp=fopen("pdata.dat","w");
#define TO_STDOUT_FP(fp,...)\
     printf(__VA_ARGS__);fprintf(fp,__VA_ARGS__)

     for(int i=0;i<PARAM_N-4;i++){
          TO_STDOUT_FP(fp,"%s %lf\n",PARAM_NAME[i],p[i]);
          //fprintf(fp,"%s %lf\n",PARAM_AME[i],i,p[i]);
     }
     fclose(fp);
}

static double* param_ref=0;
//static int ctrl_while_done=0;
static int called_sig=0;

void show_result_on_signal(int sig){
     called_sig+=1;
     if(called_sig>1){
          abort();
     }
     /*
        while(ctrl_while_done==0){
        }
        show_param(param_ref);
        */
}
void controlCG(double* x,double* p,int* swch,FUNCPTR0 f1,FUNCPTR0 f2,FUNCPTR0 f3,FUNCPTR0 f4){
     signal(SIGINT,show_result_on_signal);
     param_ref=p;
     int score=0;
     int k=0;
     double alpha=1e-10;
     double alphasub[10]={0};
     double grad[PARAM_N];for(int i=0;i<PARAM_N;++i)grad[i]=0.0;
     double d[PARAM_N]={0};
     double goodp[PARAM_N]={0};
     double goodd[PARAM_N]={0};
     double golden[PARAM_N]={0};
     // double otamesi[10][PARAM_N]={0};
     double value,lastvalue=DBL_MAX;
     int index=0;
     double gradp=0,gradc=0;
     //ctrl_while_done=0;
     called_sig=0;
     for(int i=1;i<10;i++) alphasub[i]=i/10.;
     alphasub[0]=1e-2;
     while(1){
          for(int i=0;i<PARAM_N;i++){
               goodp[i]=p[i];
               goodd[i]=d[i];
          }
          lastvalue=f1(x,goodp)+f2(x,goodp)+f3(x,goodp)+f4(x,goodp);
          if(score>=6) alpha=alpha*10;
          score+=1;
          k=k+1;
          gradient(x,p,swch,grad,f1,f2,f3,f4);
          for(int i=0;i<PARAM_N-4;i++) fprintf(stderr,"%d %lf\n",i,grad[i]);
          gradc=pow(norm(grad,PARAM_N),2);
          fprintf(stderr," %d times : cost is %15.15e, gradnorm is %lf\n",k,lastvalue,sqrt(gradc));
          if(!(sqrt(gradc)>=0.1)){
               break;
          }
          //for(int j=0;j<10;j++){
          //     alpha=alphasub[j];
          //      for(int i=0;i<PARAM_N;i++) p[i]=goodp[i];
          if(k==1){
               for(int i=0;i<PARAM_N;i++){
                    d[i]=-grad[i];
                    p[i]+=alpha*d[i];
               }
          }else{
               for(int i=0;i<PARAM_N;i++){
                    d[i]=-grad[i]+gradc/gradp*d[i];
                    p[i]=fmax(p[i]+alpha*d[i],0.);
               }
          }
          //     for(int i=0;i<PARAM_N;i++) otamesi[j][i]=p[i];
          // }
          /* for(int i=0;i<10;i++){
             for(int j=0;j<PARAM_N;j++) p[j]=otamesi[i][j];
             const double tmpcv=f(x,p);
             fprintf(stderr,"COST:%d == %g\n",i,tmpcv); 
             if(!(value<tmpcv)){
             value=tmpcv;
             index=i;
             }
             }
             for(int i=0;i<PARAM_N;i++) p[i]=otamesi[index][i];
             for(int i=0;i<PARAM_N;i++) fprintf(stderr,"p[%d] = %lf\n",i,p[i]);*/
          value=f1(x,p)+f2(x,p)+f3(x,p)+f4(x,p);
          if(fabs(lastvalue-value)<1e-4) break;
          if(!(lastvalue>value)){
               fprintf(stderr,"over\n");
               score=0;
               /*      golden_separate(x,goodp,p,golden,swch,f1,f2,f3);
                       value = f1(x,golden)+f2(x,golden)+f3(x,golden);
                       if(lastvalue<value){*/
               for(int i=0;i<PARAM_N;i++){
                    p[i]=goodp[i];
                    d[i]=goodd[i];
               }
               /*             fprintf(stderr,"good\n");
                              }else{
                              for(int i=0;i<PARAM_N;i++) p[i]=golden[i];
                              fprintf(stderr,"golden\n");
                              }*/
               if(alpha<1e-15) break;
               alpha=alpha/10;
               continue;       
          }
          //     goto LOOP_END;
          gradp=gradc;
          if(called_sig==1){
               break;
          }
     }


     //LOOP_END:;
     //ctrl_while_done=1;
}

void golden_separate(double*x,double* p1,double* p2,double* p3,const int* swch,FUNCPTR0 f1,FUNCPTR0 f2,FUNCPTR0 f3){
     double tempp1[PARAM_N],tempp2[PARAM_N],tempp3[PARAM_N],tempp4[PARAM_N];
     double fp3=100,fp4=0;
     double L=0,a=0,b=0,c=0;
     double phi=(1.+sqrt(5.0))/2.0;
     for(int i=0;i<PARAM_N;i++){
          tempp1[i]=p1[i];
          tempp2[i]=p2[i];
     }
     for(int i=0;i<PARAM_N;i++){
          if(swch==0){
               while(fabs(fp3-fp4)<=1e-4){
                    L=tempp1[i]-tempp2[i];
                    b=L/(1+phi);
                    a=phi*b;
                    c=phi*a;

                    tempp3[i]=tempp1[i]+a;
                    tempp4[i]=tempp1[i]+a+c;

                    fp3=f1(x,tempp3)+f2(x,tempp3)+f3(x,tempp3);
                    fp4=f1(x,tempp4)+f2(x,tempp4)+f3(x,tempp4);

                    if(fp3>fp4){
                         tempp1[i]=tempp3[i];
                    }else{
                         tempp2[i]=tempp4[i];
                    }
               }
               tempp1[i]=p1[i];
               tempp2[i]=p2[i];
               p3[i]=tempp3[i];
          }
     }
}


/*void jacobi(int N,double *x,double *p,FUNCPTR f){
  double* jacobi;
  double dx=0.001;
  double 
  jacobi=(double *)malloc(N*N);
  for(int i=0;i<N;i++){
  for(int j=0;j<N;j++){
  jacobi[j+i*N]=f(*/
