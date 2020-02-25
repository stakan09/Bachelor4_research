#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <random>
#include <iostream>
#include <float.h>
#include "bregsplit.h"

int TOO_MUCH=0;
static double lambda=1;

void gradient(double* p,double* d,double* b,double* grad,FUNCPTR f){
     double delta=1e-8;
     double rel_tol=1e-8;
     double abs_tol=1e-8;
     double backward=0,forward=0;
     double temppara1[PARAMSIZE],temppara2[PARAMSIZE];
     for(int i=0;i<PARAMSIZE;i++){ 
          temppara1[i]=p[i];
          temppara2[i]=p[i];
     }
     for(int i=0;i<PARAMSIZE;i++){
          double eps=fabs(p[i])*rel_tol + abs_tol;

          temppara1[i]=temppara1[i]-eps;
          backward = f(temppara1,d,b);
          temppara2[i]=temppara2[i]+eps;
          forward  = f(temppara2,d,b);
          grad[i]=(forward-backward)/(2.*eps);
          temppara1[i]=p[i];
          temppara2[i]=p[i];
     }
}

double norm(int SIZE,double* x){
     double sum=0.0;
     for(int i=0;i<SIZE;i++) sum+=pow(x[i],2);
     return sqrt(sum);
}

void vsub(int N,double* u,double* v,double* result){
     for(int i=0;i<N;i++){
          result[i]= u[i]-v[i];
     }
}

double dist(int N,double* u,double* v){
     double vec[N]={0};
     vsub(N,u,v,vec);
     return norm(N,vec);
}

double source[datasize][2]=
{
  {-1.464493, 2.144741},
{-1.454372, 2.115197},
{-0.195140, 0.038080},
{-1.915903, 3.670685},
{-0.596408, 0.355702},
{1.645432, 2.707447},
{-0.116991, 0.013687},
{-1.702300, 2.897825},
{0.279389, 0.078058},
{0.540925, 0.292600},
{-1.642187, 2.696779},
{0.224716, 0.050497},
{1.158608, 1.342372},
{-1.113465, 1.239805},
{-0.325326, 0.105837},
{-1.000888, 1.001777},
{-0.832541, 0.693125},
{1.212945, 1.471236},
{-0.101625, 0.010328},
{-0.920242, 0.846845}
};

double testfunc(double* p){
     double value[datasize]={0};
     double cost=0;
     double x=0;
     for(int j=0;j<datasize;j++){
     double pwx=1.0;
          x=source[j][0];
          for(int i=0;i<PARAMSIZE;i++,pwx*=x) value[j] += pwx*p[i];
          //printf("cccc:%g %g %g\n",source[j][0],source[j][1],value[j]);
          cost += fabs(value[j]-source[j][1])/fabs(source[j][1]);
     }
     return cost;
}
void testfunc1(double *p){
     FILE*fp;
     fp=fopen("test.txt","w");
     double dx=1e-3;
     double x=0;
     for(int j=0;j<2*1e3;j++){
          x=dx*j-1;
     double y=0;
          for(int i=0;i<PARAMSIZE;i++) y+=pow(x,i)*p[i];
          fprintf(fp,"%lf %lf\n",x,y);
     }
     fclose(fp);
}

double wrap1(double* p,double* d,double* b){
     double gsub[PARAMSIZE]={0};
     vsub(PARAMSIZE,d,p,gsub);
     vsub(PARAMSIZE,gsub,b,gsub);
     return lambda/2*pow(norm(PARAMSIZE,gsub),2)+testfunc(p);
}


double shrink(double x,double a){
     return x/fabs(x)*fmax(fabs(x)-a,0);
}

void controlCG(double* p,double* d,double* b,FUNCPTR f){
     int score=0;
     int k=0;
     double alpha=1e-2;
     double grad[PARAMSIZE];
     for(int i=0;i<PARAMSIZE;++i)grad[i]=0.0;
     double dvec[PARAMSIZE]={0};
     double goodp[PARAMSIZE]={0};
     double gooddvec[PARAMSIZE]={0};
     double value,lastvalue=DBL_MAX;
     double gradp=0,gradc=0;
     while(1){
          for(int i=0;i<PARAMSIZE;i++){
               goodp[i]=p[i];
               gooddvec[i]=dvec[i];
          }
          lastvalue=f(p,d,b);
          if(score>=6) alpha=alpha*10;
          score+=1;
          k=k+1;
          gradient(p,d,b,grad,f);
          gradc=pow(norm(PARAMSIZE,grad),2);
          //fprintf(stderr," %d times : cost is %15.15e, gradnorm is %lf\n",k,lastvalue,sqrt(gradc));
          //for(int i=0;i<PARAMSIZE;i++) fprintf(stderr," grad[%d] = %lf\n",i,grad[i]);
          if(sqrt(gradc)<10){
               break;
          }
          for(int i=0;i<PARAMSIZE;i++){
               dvec[i]=-grad[i];
               p[i]+=alpha*dvec[i];
          }
          value=f(p,d,b);
          if(fabs(lastvalue)<fabs(value)){
          //     fprintf(stderr,"over\n");
               score=0;
               for(int i=0;i<PARAMSIZE;i++){
                    p[i]=goodp[i];
                    dvec[i]=gooddvec[i];
               }
               if(alpha<1e-15) break;
               alpha=alpha/10;
               continue;       
          }
          gradp=gradc;
          if(k>1e5) {
               TOO_MUCH=1;
               break;
          }
     }
}



void SplitBregman(double* p,FUNCPTR f){
     double tempp[PARAMSIZE]={0};
     double tempd[PARAMSIZE]={0};
     double d[PARAMSIZE]={0};
     double tempb[PARAMSIZE]={0};
     double b[PARAMSIZE]={0};
     double grad[PARAMSIZE]={0};
     double tol=1e-6;
     double usub[PARAMSIZE]={0};
     do{  
          for(int i=0;i<PARAMSIZE;i++) tempp[i]=p[i];
         //fprintf(stderr,"CG\n");
          controlCG(p,d,b,f);
         //fprintf(stderr,"d,b\n");
          for(int i=0;i<PARAMSIZE;i++){
               d[i]=shrink(tempp[i]+b[i],1/lambda);
          }

          for(int i=0;i<PARAMSIZE;i++){
               b[i]+=p[i]-d[i];
          }
     }while(dist(PARAMSIZE,p,tempp)>tol);
}


int main(int argc, const char **argv){
     if(argc < 3){
          fprintf(stderr, "Usage: ./plot seed\n");
          exit(1);
     }
     double p[PARAMSIZE]={0};
     int count=atoi(argv[1]);
     int flag=atoi(argv[2]);
     std::mt19937_64 engine;
     std::uniform_real_distribution<double> range0(-1,1);
     std::uniform_int_distribution<int> range1(-1000,1000);
     engine.seed(count);
     for(int i=0;i<PARAMSIZE;i++){
          //p[i]=range0(engine);
          //p[i]=i==2?1.0:0.0;
          p[i]=(double)range1(engine)/1000.;
          printf("p[%d] = %lf\n",i,p[i]);
     }
    // printf("before %lf\n",testfunc(p));
     SplitBregman(p,wrap1);
     if(TOO_MUCH==0){
          printf("%d %lf\n",count,testfunc(p));
          for(int i=0;i<PARAMSIZE;i++) printf("p[%d] = %lf\n",i,p[i]);
     }
     TOO_MUCH=0;
     if(flag==0) testfunc1(p);
     return 0;
}

