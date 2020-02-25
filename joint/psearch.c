#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "pfunc.h"


int main(){
     int swch[PARAM_N];
     double x[SHURUI] = {}, p[PARAM_N] = {};
     char dataname[128];
     FILE *fp1;
     // FILE *fp2;
     // fp2=fopen("LdataPHE","w");
     // do{
     //    count+=1;
     //    fprintf(stderr,"Count:%d\n",count);
     //Lalp=1./N;
     //////////////////
     x[X_SCN]  =-1;
     x[X_Phe]  =-1;
     x[X_S2O3] =-1;
     x[S_SCN]  =-1;
     x[S_Phe]  =-1;
     x[S_S2O3] =-1;
     //////////////////
    /* p[alpha_SCN]   = 72;//50.064424; 
     p[alpha_Phe]   = 78;//4.966141;
     p[alpha_S2O3]  = 15;
     p[X_max_SCN]   = 1;
     p[X_max_Phe]   = 1;
     p[X_max_S2O3]  = 1;
     p[d_SCN_SCN]   = 0.009;
     p[d_Phe_Phe]   = 0.009;
     p[d_S2O3_S2O3] = 0.009;
     p[d_SCN_Phe]   = 0.537594;
     p[d_SCN_S2O3]  = 0.353003;
     p[d_Phe_SCN]   = 0.704144;
     p[d_Phe_S2O3]  = 0.679520;
     p[d_S2O3_SCN]  = 0.138756;
     p[d_S2O3_Phe]  = 0.885283;
     p[c_SCN]       = 48; 
     p[c_Phe]       = 28; 
     p[c_S2O3]      = 43;
     p[beta_SCN]    = 80;
     p[beta_Phe]    = 63 ;
     p[beta_S2O3]   = 23;
     p[gamma_SCN]   = 76;
     p[gamma_Phe]   = 100;
     p[gamma_S2O3]  = 35;
     p[A_SCN]       = 30; 
     p[A_Phe]       = 40;
     p[A_S2O3]      = 26;
     p[B_SCN]       = 93 ;
     p[B_Phe]       = 83;
     p[B_S2O3]      = 27;*/
     for(int i=0;i<PARAM_N-4;i++) scanf("%lf",&p[i]);


     /////////////////////
     p[Pin]=-1;
     p[P_SCN]=-1;
     p[P_Phe]=-1;
     p[P_S2O3]=-1;

     //////////////SWITCHER//////////////////
     swch[alpha_SCN]   = 1;
     swch[alpha_Phe]   = 1;
     swch[alpha_S2O3]  = 1;
     swch[X_max_SCN]   = 1;
     swch[X_max_Phe]   = 1;
     swch[X_max_S2O3]  = 1;
     swch[d_SCN_SCN]   = 0;
     swch[d_Phe_Phe]   = 1;
     swch[d_S2O3_S2O3] = 1;
     swch[d_SCN_Phe]   = 1;
     swch[d_SCN_S2O3]  = 1;
     swch[d_Phe_SCN]   = 1;
     swch[d_Phe_S2O3]  = 1;
     swch[d_S2O3_SCN]  = 1;
     swch[d_S2O3_Phe]  = 1;
     swch[c_SCN]       = 0;
     swch[c_Phe]       = 1;
     swch[c_S2O3]      = 1;
     swch[beta_SCN]    = 0;
     swch[beta_Phe]    = 1;
     swch[beta_S2O3]   = 1;
     swch[gamma_SCN]   = 0;
     swch[gamma_Phe]   = 1;
     swch[gamma_S2O3]  = 1;
     swch[Pin]         = 1;
     swch[P_SCN]       = 1;
     swch[P_Phe]       = 1;
     swch[P_S2O3]      = 1;
     swch[A_SCN]       = 0;
     swch[A_Phe]       = 1;
     swch[A_S2O3]      = 1;
     swch[B_SCN]       = 0 ;
     swch[B_Phe]       = 1;
     swch[B_S2O3]      = 1;

     ////////////////////////////////////////
     controlCG(x,p,swch,target_func3sp,zeroscore,zeroscore,zeroscore);
     printf("final cost : 3sp=%.15e\n , SCN=%.15e\n , Phe=%.15e\n , S2O3=%.15e\n ",target_func3sp(x,p),target_funcSCN(x,p),target_funcPhe(x,p),target_funcS2O3(x,p));
     show_param(p);
     //fprintf(fp2,"%lf %.15e\n",Lalp,target_func(x,p));
    target_func3sp(x,p);
     target_funcSCN(x,p);
     target_funcPhe(x,p);
     target_funcS2O3(x,p);

     sprintf(dataname,"optdataS3sp");
     fp1=fopen(dataname,"w");
     for(int i=0;i<9;i++){
          fprintf(fp1,"%d %lf %lf %lf\n",i,speeddata3sp[i][0],speeddata3sp[i][1],speeddata3sp[i][2]);
          fprintf(fp1,"%d %lf %lf %lf\n",i+1,speeddata3sp[i][0],speeddata3sp[i][1],speeddata3sp[i][2]);
     }
     fclose(fp1);
     sprintf(dataname,"optdataSSCN");
     fp1=fopen(dataname,"w");
     for(int i=0;i<7;i++){
          fprintf(fp1,"%d %lf %lf %lf\n",i,speeddataSCN[i][0],speeddataSCN[i][1],speeddataSCN[i][2]);
          fprintf(fp1,"%d %lf %lf %lf\n",i+1,speeddataSCN[i][0],speeddataSCN[i][1],speeddataSCN[i][2]);
     }
     fclose(fp1);
     sprintf(dataname,"optdataSPhe");
     fp1=fopen(dataname,"w");
     for(int i=0;i<9;i++){
          fprintf(fp1,"%d %lf %lf %lf\n",i,speeddataPhe[i][0],speeddataPhe[i][1],speeddataPhe[i][2]);
          fprintf(fp1,"%d %lf %lf %lf\n",i+1,speeddataPhe[i][0],speeddataPhe[i][1],speeddataPhe[i][2]);
     }
     fclose(fp1);
     sprintf(dataname,"optdataSS2O3");
     fp1=fopen(dataname,"w");
     for(int i=0;i<10;i++){
          fprintf(fp1,"%d %lf %lf %lf\n",i,speeddataS2O3[i][0],speeddataS2O3[i][1],speeddataS2O3[i][2]);
          fprintf(fp1,"%d %lf %lf %lf\n",i+1,speeddataS2O3[i][0],speeddataS2O3[i][1],speeddataS2O3[i][2]);
     }
     fclose(fp1);

     sprintf(dataname,"optdataX3sp");
     fp1=fopen(dataname,"w");
     for(int i=0;i<9;i++) fprintf(fp1,"%d %lf %lf %lf\n",i,Xdata3sp[i][0],Xdata3sp[i][1],Xdata3sp[i][2]);
     fclose(fp1);
     sprintf(dataname,"optdataXSCN");
     fp1=fopen(dataname,"w");
     for(int i=0;i<7;i++) fprintf(fp1,"%d %lf %lf %lf\n",i,XdataSCN[i][0],XdataSCN[i][1],XdataSCN[i][2]);
     fclose(fp1);
     sprintf(dataname,"optdataXPhe");
     fp1=fopen(dataname,"w");
     for(int i=0;i<9;i++) fprintf(fp1,"%d %lf %lf %lf\n",i,XdataPhe[i][0],XdataPhe[i][1],XdataPhe[i][2]);
     fclose(fp1);
     sprintf(dataname,"optdataXS2O3");
     fp1=fopen(dataname,"w");
     for(int i=0;i<10;i++) fprintf(fp1,"%d %lf %lf %lf\n",i,XdataS2O3[i][0],XdataS2O3[i][1],XdataS2O3[i][2]);
     fclose(fp1);




     // }while(count<N);
     //fclose(fp2);
     return 0;
}

