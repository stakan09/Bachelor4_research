#include <stdlib.h>
#include "SOLVE.h"

///type wake///
char DEtypedesu[3]={'D','N','S'};
char* DEtype=DEtypedesu;


///matrix tukuru///
void make_mat(char DE,int bunkatu,LUmat* mat_data){
     if(DE==DEtype[2]){
          mat_data->L1=(double*)malloc((bunkatu-1)*sizeof(double));
          mat_data->L2=(double*)malloc((bunkatu-1)*sizeof(double));
          mat_data->L3=(double*)malloc(bunkatu*sizeof(double));
          mat_data->U1=(double*)malloc((bunkatu-2)*sizeof(double));
          mat_data->U2=(double*)malloc((bunkatu-1)*sizeof(double));
          mat_data->y=(double*)malloc(bunkatu*sizeof(double));
     }
     else{
          mat_data->L1=(double*)malloc(bunkatu*sizeof(double));
          mat_data->L2=(double*)malloc(bunkatu*sizeof(double));
          mat_data->L3=NULL;
          mat_data->U1=(double*)malloc(bunkatu*sizeof(double));
          mat_data->U2=NULL;
          mat_data->y=(double*)malloc(bunkatu*sizeof(double));
     }

}


///dirchlet boundary condition no LU bunkai///
void __LUbunkai0(int N,double* a,double* b,double* c,LUmat* mat_data){
     double* L1=mat_data->L1;
     double* L2=mat_data->L2;
     double* L3=mat_data->L3;
     double* U1=mat_data->U1;
     double* U2=mat_data->U2;
     double* y=mat_data->y;

     int i;
     L1[1]=a[1];
     U1[1]=c[1]/L1[1];
     for(i=2;i<=N-1;i++){
          L2[i]=b[i];
          L1[i]=a[i]-L2[i]*U1[i-1];
          U1[i]=c[i]/L1[i];
     }
}



///neumann boundary condition no LU bunkai///
void __LUbunkai1(int N,double* a,double* b,double* c,LUmat* mat_data){
     double* L1=mat_data->L1;
     double* L2=mat_data->L2;
     double* L3=mat_data->L3;
     double* U1=mat_data->U1;
     double* U2=mat_data->U2;
     double* y=mat_data->y;

     int i;
     L1[0]=a[0];
     U1[0]=c[0]/L1[0];
     for(i=1;i<=N-1;i++){
          L2[i]=b[i];
          L1[i]=a[i]-L2[i]*U1[i-1];
          U1[i]=c[i]/L1[i];
     }
}    



///Periodic boundary condition no LU bunkai///
void __LUbunkai2(int N,double* a,double* b,double* c,LUmat* mat_data){
     double* L1=mat_data->L1;
     double* L2=mat_data->L2;
     double* L3=mat_data->L3;
     double* U1=mat_data->U1;
     double* U2=mat_data->U2;
     double* y=mat_data->y;

     int i,j;
     L1[0]=a[0];
     U2[0]=b[0]/L1[0];
     U1[0]=c[0]/L1[0];
     for(i=1;i<=N-3;i++){
          L2[i]=b[i];
          L1[i]=a[i]-L2[i]*U1[i-1];
          U1[i]=c[i]/L1[i];
          U2[i]=-L2[i]*U2[i-1]/L1[i];
     }
     // i=N-2 no toki
     L2[N-2]=b[N-2];
     L1[N-2]=a[N-2]-L2[N-2]*U1[N-3];
     U2[N-2]=(c[N-2]-L2[N-2]*U2[N-3])/L1[N-2];

     L3[0]=c[N-1];
     for(i=1;i<=N-3;i++){
          L3[i]=-L3[i-1]*U1[i-1];
     }
     L3[N-2]=b[N-1]-L3[N-3]*U1[N-3];
     double sum=0.0;
     for(i=0;i<=N-2;i++){//siguma no siki
          sum = sum + L3[i]*U2[i];
     }
     L3[N-1]=a[N-1]-sum;
}


///LU bunkai///
void LU(char DE,int N,double *a,double* b,double* c,LUmat* mat_data){
     if(DE==DEtype[0]){
          __LUbunkai0(N,a,b,c,mat_data);
     }else if(DE==DEtype[1]){
          __LUbunkai1(N,a,b,c,mat_data);
     }else if(DE==DEtype[2]){
          __LUbunkai2(N,a,b,c,mat_data);
     }
}



///dirichlet boundary condition no solver///
void __solveD(int N,LUmat* mat_data,double* u_old,double* u_new){
     double* L1=mat_data->L1;
     double* L2=mat_data->L2;
     double* L3=mat_data->L3;
     double* U1=mat_data->U1;
     double* U2=mat_data->U2;
     double* y=mat_data->y;


     int i;
     y[1]=u_old[1]/L1[1];
     for(i=2;i<=N-1;i++){
          y[i]=(u_old[i]-L2[i]*y[i-1])/L1[i];
     }
     u_new[N-1]=y[N-1];
     for(i=N-2;i>=1;i--){
          u_new[i]=y[i]-U1[i]*u_new[i+1];
     }
}




///neumann boundary condition no solver///
void __solveN(int N,LUmat* mat_data,double* u_old,double* u_new){
     double* L1=mat_data->L1;
     double* L2=mat_data->L2;
     double* L3=mat_data->L3;
     double* U1=mat_data->U1;
     double* U2=mat_data->U2;
     double* y=mat_data->y;


     int i;
     y[0]=u_old[0]/L1[0];
     for(i=1;i<=N-1;i++){
          y[i]=(u_old[i]-L2[i]*y[i-1])/L1[i];
     }
     u_new[N-1]=y[N-1];
     for(i=N-2;i>=0;i--){
          u_new[i]=y[i]-U1[i]*u_new[i+1];
     }
}



///Periodic boundary condition no solver///
void __solveS(int N,LUmat* mat_data,double* u_old,double* u_new){
     double* L1=mat_data->L1;
     double* L2=mat_data->L2;
     double* L3=mat_data->L3;
     double* U1=mat_data->U1;
     double* U2=mat_data->U2;
     double* y=mat_data->y;



     int i,j;
     y[0]=u_old[0]/L1[0];//Ly=u_old
     for(i=1;i<=N-2;i++){
          y[i]=(u_old[i]-L2[i]*y[i-1])/L1[i];
     }
     double sum = 0.0;
     for(i=0;i<=N-2;i++){//sigma tukuru
          sum = sum + L3[i]*y[i];
     }
     y[N-1]=(u_old[N-1]-sum)/L3[N-1];

     u_new[N-1]=y[N-1];//Uu_new=y
     u_new[N-2]=y[N-2]-U2[N-2]*u_new[N-1];
     for(i=N-3;i>=0;i--){
          u_new[i]=y[i]-U1[i]*u_new[i+1]-U2[i]*u_new[N-1];
     }
     for(i=0;i<=N-1;i++){
          u_old[i]=u_new[i];
     }
}


///SOLVER wakeru////
void SOLVE(char DE,int N,LUmat* mat_data,double* u_old,double* u_new){
     if(DE==DEtype[0]){
         __solveD(N,mat_data,u_old,u_new);
     }
     else if(DE==DEtype[1]){
          __solveN(N,mat_data,u_old,u_new);
     }
     else if(DE==DEtype[2]){
         __solveS(N,mat_data,u_old,u_new);
     }
}


///matrix kaihou///
void destroymat(LUmat* mat_data) {

     if (mat_data->L1) {
          free(mat_data->L1);
          mat_data->L1 = NULL;
     }

     if (mat_data->L2) {
          free(mat_data->L2);
          mat_data->L2 = NULL;
     }

     if (mat_data->L3) {
          free(mat_data->L3);
          mat_data->L3 = NULL;
     }

     if (mat_data->U1) {
          free(mat_data->U1);
          mat_data->U1 = NULL;
     }

     if (mat_data->U2) {
          free(mat_data->U2);
          mat_data->U2 = NULL;
     }

     if (mat_data->y) {
          free(mat_data->y);
          mat_data->y = NULL;
     }
}

