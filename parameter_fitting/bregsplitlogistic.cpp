#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <random>
#include <iostream>
#include <float.h>
#include <assert.h>
#include "bregsplit.h"
//#include <boost/numeric/odeint.hpp>
#define bugsearch(name) for(int z=0;z<Num_V;z++) fprintf(stderr, "%d:%lf\n",__LINE__,name[z])
#define bsx16 bugsearch(x1);bugsearch(x2);bugsearch(x3);bugsearch(x4);bugsearch(x5);bugsearch(x6)

#define check(st,end) fprintf(stderr,"%d\n",__LINE__);for(int i=st;i<end;i++) fprintf(stderr,"data[%d] == %lf\n",i,data[i])
int TOO_MUCH=0;
double lambdaL1=1e-3;
double lambdaL2=1e-5;


///runge kutta 4ji///
void rk(double *__dt,double* x,double* p,ODE func){
     double dt=*__dt;
     double dth=1./2*dt;
     double tempx[Num_V]={},x1[Num_V]={},x2[Num_V]={},x3[Num_V]={},x4[Num_V]={};
     for(int i=0;i<Num_V;i++)tempx[i]=x[i];
     func(tempx,p,x1);
     for(int i=0;i<Num_V;i++) tempx[i]=x[i]+dth*x1[i];
     func(tempx,p,x2);
     for(int i=0;i<Num_V;i++) tempx[i]=x[i]+dth*x2[i];
     func(tempx,p,x3);
     for(int i=0;i<Num_V;i++) tempx[i]=x[i]+dt*x3[i];
     func(tempx,p,x4);
     for(int i=0;i<Num_V;i++) x[i]=x[i]+dt*(x1[i]+2*x2[i]+2*x3[i]+x4[i])/6.0;
}


///runge kutta fehlberg///
void rkf(double* __dt,double* x,double* p,ODE f){
     double dt=*__dt;
     double x1[Num_V]={},x2[Num_V]={},x3[Num_V]={},x4[Num_V]={},x5[Num_V]={},x6[Num_V]={},tempx[Num_V]={};
     double eps=1E-5;
     double error = 0.0;
     double s = 0.0;
     double tol=1e-8;
     while(1){          

          error=0.0;
          for(int i=0;i<Num_V;i++) tempx[i]=x[i];

          f(tempx,p,x1);

          for(int i=0;i<Num_V;i++) tempx[i]=x[i] + 1./4*dt*x1[i];

          f(tempx,p,x2);

          for(int i=0;i<Num_V;i++) tempx[i]=x[i] + 3./32*dt*x1[i] + 9./32*dt*x2[i];

          f(tempx,p,x3);

          for(int i=0;i<Num_V;i++) tempx[i]=x[i] + 1932./2197*dt*x1[i] -7200./2197*dt*x2[i] + 7296./2197*dt*x3[i];

          f(tempx,p,x4);

          for(int i=0;i<Num_V;i++) tempx[i]=x[i] + 439./216*dt*x1[i] -8*dt*x2[i] + 3680./513*dt*x3[i] -845./4104*dt*x4[i];

          f(tempx,p,x5);

          for(int i=0;i<Num_V;i++) tempx[i]=x[i] -8./27*dt*x1[i] + 2.*dt*x2[i] -3544./2565*dt*x3[i] + 1859./4104*dt*x4[i] -11./40*dt*x5[i];

          f(tempx,p,x6);

          for(int i=0;i<Num_V;i++) error += 1./(fabs(x[i])+tol)*dt*fabs((16./135-25./216)*x1[i] + (6656./12825-1408./2565)*x3[i] + (28561./56430-2197./4104)*x4[i] - (9./50-1./5)*x5[i] + 2./55*x6[i]);
          //for(int i=0;i<Num_V;i++) fprintf(stderr,"%d %lf %lf %lf %lf %lf %lf\n",i,x1[i],x2[i],x3[i],x4[i],x5[i],x6[i]);
          s=fmin(fmax(pow(eps*dt/(2.*(error+tol)),1.0/4.0),1e-1),10);
          // fprintf(stderr,"S:%.10E\n",s);
          // fprintf(stderr,"DT:%.10E\n",dt);
          //fprintf(stderr,"ERROR:%.10E\n",error);
          assert(!(isnan(error)));
          if(s>1){
               for(int i=0;i<Num_V;i++) x[i] += dt*(25./216*x1[i] + 1408./2565*x3[i] + 2197./4104*x4[i] -1./5*x5[i]);
               dt=s*dt;
               break;
          }else{
               dt=s*dt;
          }
          /*if(dt<1e-50) {
            TOO_MUCH=2;
            return;
            }*/

     }
     *__dt=dt;
}



///tokitai ODE///
void testfunc2(double* x,double* p,double* dx){
     //for(int i=0;i<Num_V;i++) dx[i]=p[1]+p[2]*x[i]+p[3]*pow(x[i],2)+p[4]*pow(x[i],3)+p[5]*pow(x[i],4);  
     //for(int i=0;i<Num_V;i++) dx[i]=p[1]*x[i]*(p[2]-x[i]) + p[3] + p[4]*pow(x[i] + p[5]*pow(x[i],4);
     for(int i=0;i<Num_V;i++) dx[i]=p[1]*x[i]*(p[2]-x[i]) + p[3] + p[4]*x[i]/(1+x[i]) + p[5]*pow(x[i],2)/(1+pow(x[i],2)) ;
     for(int i=0;i<Num_V;i++) if(isfinite(x[i])==0) exit(1);
}



///gradient keisan///
void gradient(double* p,double* d,double* b,double* grad,int *swch,ENTIREFUNC f){
     double rel_tol=1e-8;
     double abs_tol=1e-8;
     double backward=0,forward=0;
     double temppara1[PARAMSIZE]={},temppara2[PARAMSIZE]={};
     for(int i=0;i<PARAMSIZE;i++) temppara1[i]=temppara2[i]=p[i];
     for(int i=0;i<PARAMSIZE;i++){
          if(swch[i]==0){
               double eps=fabs(p[i])*rel_tol + abs_tol;
               //fprintf(stderr,"eps == %10E\n",eps);
               temppara1[i]=p[i]-eps;
               /*fprintf(stderr,"grad[%d] == %10E\n",i,grad[i]);
                 fprintf(stderr,"p[%d] == %10E\n",i,p[i]);*/

               backward = f(temppara1,d,b);
               //fprintf(stderr,"backward\n");
               temppara2[i]=p[i]+eps;
               /* fprintf(stderr,"grad[%d] == %10E\n",i,grad[i]);
                  fprintf(stderr,"p[%d] == %10E\n",i,p[i]);*/

               forward  = f(temppara2,d,b);
               /* fprintf(stderr,"forward\n");
                  fprintf(stderr,"grad[%d] == %.10E\n",i,grad[i]);
                  fprintf(stderr,"p[%d] == %.10E\n",i,p[i]);
                  fprintf(stderr,"result\n");*/
               grad[i]=(forward-backward)/(2.*eps);
               //fprintf(stderr,"grad[%d] == %.10E\n",i,grad[i]);
               //fprintf(stderr,"p[%d] == %.10E\n",i,p[i]);
               temppara1[i]=p[i];
               temppara2[i]=p[i];
          }else{
               grad[i]=0;
          }
     }
}



///norm///
double norm(int SIZE,double* x){
     double sum=0.0;
     for(int i=0;i<SIZE;i++) sum+=pow(x[i],2);
     return sqrt(sum);
}



///vector hikizan///
void vsub(int N,double* u,double* v,double* result){
     for(int i=0;i<N;i++){
          result[i]= u[i]-v[i];
     }
}



///kyori keisan///
double dist(int N,double* u,double* v){
     std::vector<double> vec(N);
     vsub(N,u,v,vec.data());
     return norm(N,vec.data());
}



double dist1(int N,double* u,double* trueval){
     double cost=0;
     std::vector<double> vec(N);
     vsub(N,u,trueval,vec.data());
     for(int i=0;i<N;i++) cost+=fabs(vec[i]-trueval[i])/fabs(trueval[i]);
     return sqrt(cost);
}

double Vardata[datasize*Num_V]=
{};// line;Num_V,column:datasize

double timedata[datasize]=
{};

double timedataplot[PLOTSIZE]=
{}; 
double Vardataplot[PLOTSIZE*Num_V]=
{};



///ODE keisan///
void ODErecord(double* p,double* result,double* time,int data_s,ODE f){
     double x[Num_V]={};
     double xprev[Num_V]={};
     double t=0;
     double tprev=0;
     double dtprev=0;
     double dt=1e-4;
     int i=0;
     int j=0;
     for(i=0;i<Num_V;i++) x[i]=xprev[i]=p[0];
     //size_t steps = integrate( harmonic_oscillator, x ,0.0,4.0,0.00
     do{  
          for(i=0;i<Num_V;i++){
               xprev[i]=x[i];
          }
          dtprev=dt;
          tprev=t;
          rk(&dt,x,p,f);
          t+=dt;

          if(j<data_s && t>time[j] && fabs(time[j]-t)>1e-8){
               //fprintf(stderr,"aaaaaaaa\n");
               dt=time[j]-tprev;
               //fprintf(stderr,"%lf\n",dt);
               rk(&dt,xprev,p,f);
               //for(int i=0;i<Num_V;i++) fprintf(stderr,"%lf\n",xprev[i]);
               t=time[j];
               for(i=0;i<Num_V;i++) result[j+i*data_s]=x[i]=xprev[i];
               j+=1;
               dt=1e-4;
          }else if(j<data_s && t>time[j] && !(fabs(time[j]-t)>1e-8)){
               //fprintf(stderr,"bbbbbbbb\n");
               t=time[j];
               for(i=0;i<Num_V;i++) result[j+i*data_s]=xprev[i]=x[i];
               j+=1;
          }
                 //fprintf(stderr,"time == %lf, dt == %10E\n",t,dt);
     }while(t<end_time);
}


///ODE no moto data tono sa wo keisan///
double ODEcost(double* p){
     double result[datasize*Num_V]={};
     double cost[Num_V]={};
     ODErecord(p,result,timedata,datasize,testfunc2);
     for(int i=0;i<Num_V;i++) cost[i]=dist(datasize,&result[i*datasize],&Vardata[i*datasize]);
     return pow(norm(Num_V,cost),2);
}


///wrapper, cost zentai wo keisan,L1 saitekika you///
double wrap1(double* p,double* d,double* b){
     double gsub[PARAMSIZE]={0};
     vsub(PARAMSIZE,d,p,gsub);
     vsub(PARAMSIZE,gsub,b,gsub);
     //fprintf(stderr,"hidari == %.10E\n",lambdaL1/2*pow(norm(PARAMSIZE,gsub),2));
     //fprintf(stderr,"migi== %.10E\n",ODEcost(p));
     return lambdaL1/2.*pow(norm(PARAMSIZE,gsub),2)+ODEcost(p);
}


///wrapper, cost zentai wo keisan, L2 norm saitekika you///
double wrap2(double* p,double* d,double* b){
     return ODEcost(p)+lambdaL2*pow(norm(PARAMSIZE,p),2);
}


///L1 norm you no shrink///
double shrink(double x,double a){
     return x/fabs(x)*fmax(fabs(x)-a,0);
}


///L2 norm saitekika////
void controlCG(double* p,double* d,double* b,int* swch,ENTIREFUNC f){
     int score=0;
     int k=0;
     double alpha=1e-5;
     double grad[PARAMSIZE]={};
     double dvec[PARAMSIZE]={0};
     double goodp[PARAMSIZE]={0};
     double gooddvec[PARAMSIZE]={0};
     double value,lastvalue=DBL_MAX;
     double gradc=0,paradist=0;
     while(1){
          if(score>3){
               alpha=fmin(10*alpha,1e-1);
          }

         /* for(int i=0;i<PARAMSIZE;i++)  fprintf(stderr,"p[%d] = %lf\n",i,p[i]);
          fprintf(stderr,"\n");
          *///fprintf(stderr,"%d\n",score);
          lastvalue=f(p,d,b);
          score+=1;
          k=k+1;
          //fprintf(stderr,"====grad====\n");
          gradient(p,d,b,grad,swch,f);
          //fprintf(stderr,"grad fin\n");
          if(TOO_MUCH==2 || TOO_MUCH==3){
               alpha*=0.9;
               // fprintf(stderr,"TOO MUCH\n");
               for(int i=0;i<PARAMSIZE;i++){
                   if(i<c){
                    p[i]=fmax(p[i]+alpha*dvec[i],0);
                  }else{
                    p[i]=p[i]+alpha*dvec[i];
                    }
               }
               TOO_MUCH=0;
               continue;
          }

          // fprintf(stderr,"calculate gradc\n");
          gradc=norm(PARAMSIZE,grad);
          if(gradc<1e-3){
               break;
          }
          for(int i=0;i<PARAMSIZE;i++){
               goodp[i]=p[i];
               gooddvec[i]=dvec[i];
          }
          
          for(int i=0;i<PARAMSIZE;i++){
               dvec[i]=-grad[i];
               if(i<c){
                   p[i]=fmax(p[i]+alpha*dvec[i],0);
                }else{
               p[i]=p[i]+alpha*dvec[i];
                }
          }
          //fprintf(stderr,"calculate value\n");
          paradist=dist(PARAMSIZE,p,goodp);
          if(paradist<1e-7)break;
          value=f(p,d,b);
          //fprintf(stderr,"last=%lf current=%lf\n",lastvalue,value);
          //if(fabs(value-lastvalue)/fmax(lastvalue,value)<1e-7) break;
          if(!(fabs(lastvalue)>fabs(value)) || isfinite(value)==0 || TOO_MUCH==2){
               //fprintf(stderr,"===over===\n");
               score=0;
               while(1){
                    for(int i=0;i<PARAMSIZE;i++){
                         p[i]=goodp[i];
                         dvec[i]=gooddvec[i];
                    }
                    if(alpha<1e-20) {
                         TOO_MUCH=1;
                         return;
                    }
                    alpha=alpha/10;
                    for(int i=0;i<PARAMSIZE;i++){
                         dvec[i]=-grad[i];
                                   if(i<c){
                                      p[i]=fmax(p[i]+alpha*dvec[i],0);
                               }else{
                         p[i]=p[i]+alpha*dvec[i];
                             }
                    }
                    value=f(p,d,b);
                    if(fabs(lastvalue)>fabs(value) && !(isfinite(value)==0)) {
                         break;
                         TOO_MUCH=0;
                    }
               }
               //fprintf(stderr,"==========\n");
          }
          /* if(k>1e5) {
             TOO_MUCH=1;
             break;
             }*/
     }
}


///L1 norm saitekika you no SplitBregman hou///
void SplitBregman(double* p,int *swch,ENTIREFUNC f){
     char filename[64]={};
     FILE *fp;
     double value=0,lastvalue=0;
     double tempp[PARAMSIZE]={0};
     double d[PARAMSIZE]={0};
     double b[PARAMSIZE]={0};
     double tol=1e-7;
     int count = 0;
     /*   gradient(p,d,b,grad,swch,f);
          sprintf(filename,"gradb_%f.dat",log10(lambdaL1));
          fp=fopen(filename,"w");
          for(int i=0;i<PARAMSIZE;i++) fprintf(fp,"grad[%d] == %lf\n",i,grad[i]);
          fclose(fp);
          */
     do{  
          for(int i=0;i<c;i++) p[i]=fmax(p[i],0);
          count+=1;
          //fprintf(stderr," %d time\n",count);
          for(int i=0;i<PARAMSIZE;i++) tempp[i]=p[i];
          //fprintf(stderr,"========CG=========\n");
          controlCG(p,d,b,swch,f);
          if(TOO_MUCH==1) exit(1);
          //fprintf(stderr,"======d and b======\n");
          for(int i=0;i<PARAMSIZE;i++){
               if(swch[i]==0)d[i]=shrink(p[i]+b[i],1./lambdaL1);
          }
          for(int i=0;i<PARAMSIZE;i++){
               if(swch[i]==0)b[i]+=p[i]-d[i];
          }
          value=f(p,d,b);
     }while(dist(PARAMSIZE,p,tempp)/norm(PARAMSIZE,p)>tol);
         /*gradient(p,d,b,grad,swch,f);
       sprintf(filename,"grada_%f.dat",log10(lambdaL1));
       fp=fopen(filename,"w");
       for(int i=0;i<PARAMSIZE;i++) fprintf(fp,"grad[%d] == %lf\n",i,grad[i]);
       fclose(fp);
       */
     //for(int i=0;i<PARAMSIZE;i++) gradient(p,d,b,grad,f); kesukoukawo gradientsanshousite ireru
}


void replace(int* lft,int* block,int right,int *rnd,double piv,double* data){
     double rep=0;
     int left = *lft;
     int rend=*rnd;
     if(data[left] >= piv && data[right] < piv && right > left){
          rep = data[left];
          data[left] = data[right];
          data[right] = rep;
          *block = left;
          left+=1;
          rend-=1;
     }
     *lft = left;
     *rnd = rend;
}



void quicksort(int st,int end,double* data){
     int left,right,middle;
     right=left=0;
     int ok=0;
     int block=-100;
     middle=st;
     int rend=end;
     double piv=data[middle];
     while(1){
          for(left=st;left<end;left++){
               for(right=rend-1;right>=st;right--){
                    replace(&left,&block,right,&rend,piv,data);
               }
          }
          if(block==-100){
               middle+=1;
               if(middle==end) break;
               piv=data[middle];
               continue;
          }
          break;
     }
     for(int i=st;i<end-1;i++){
          if(data[i]>data[i+1]) ok+=1;
     }
     if(ok>0 && end>middle && middle>=0){ 
          quicksort(st,block+1,data);
          quicksort(block+1,end,data);
     }
     ok=0;
     for(int i=st;i<end-1;i++){
          if(data[i]>data[i+1]) ok+=1;
     }
     if(ok>0 && end>middle && middle>=0){
          quicksort(st,end,data);
     }

}


