#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define datasize 20
#define PARAMSIZE 10

extern int TOO_MUCH;
typedef double(*FUNCPTR)(double*,double*,double*);
void gradient(double* p,double* d,double* b,double* grad,FUNCPTR f);
double norm(int SIZE,double* x);
void vsub(int N,double* u,double* v,double* result);
double soutai(int N,double* u,double* v);
double testfunc(double* p);
double wrap1(double* p,double* d,double* b);
double shrink(double x,double a);
void controlCG(double* p,double* d,double* b,FUNCPTR f);
void SplitBregman(double* p,FUNCPTR f);
void testfunc1(double *p);
void testfunc2(double *p);
