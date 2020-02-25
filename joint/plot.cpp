#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "pfunc.h"
#include <random>
#include <iostream>


int main(int argc, const char **argv){

     if(argc < 2){
          fprintf(stderr, "Usage: ./plot seed\n");
          exit(1);
     }

     double x[SHURUI],p[PARAM_N] = {};
     char filename[128]={0};
     int count=atoi(argv[1]);
     x[X_SCN]=-1;
     x[X_Phe]=-1;
     x[X_S2O3]=-1;
     x[S_SCN]=-1;
     x[S_Phe]=-1;
     x[S_S2O3]=-1;
     p[alpha_SCN] = 1;
     p[alpha_Phe] = 1;
     p[alpha_S2O3] = 1;
     p[X_max_SCN] = 1.000000 ;
     p[X_max_Phe] = 1.000000 ;
     p[X_max_S2O3] = 1.000000 ;
     p[d_SCN_SCN] = 0.009000 ;
     p[d_SCN_Phe] = 0.652937 ;
     p[d_SCN_S2O3] = 0.301683 ;
     p[d_Phe_SCN] = 1.173825 ;
     p[d_Phe_Phe] = 0.009000 ;
     p[d_Phe_S2O3] = 1.263849 ;
     p[d_S2O3_SCN] = 0.503777 ;
     p[d_S2O3_Phe] = 0.201693 ;
     p[d_S2O3_S2O3] = 0.009000 ;
     p[c_SCN] = 63.997726 ;
     p[c_Phe] = 20.141059 ;
     p[c_S2O3] = 45.441319 ;
     p[beta_SCN] = 45.995207 ;
     p[beta_Phe] = 9.997722 ;
     p[beta_S2O3] = 33.982746 ;
     p[gamma_SCN] = 4.002079 ;
     p[gamma_Phe] = 0.926137 ;
     p[gamma_S2O3] = 84.563511 ;
     p[A_SCN] = 44.937399 ;
     p[A_Phe] = 4.335593 ;
     p[A_S2O3] = 3.966940 ;
     p[B_SCN] = 75.995412 ;
     p[B_Phe] = 27.031486 ;
     p[B_S2O3] = 58.056881 ;
     /*     p[alpha_SCN] = 55.585437;
            p[alpha_Phe] = 71.156726;
            p[alpha_S2O3] = 0.305910;
            p[X_max_SCN] = 1.000000 ;
            p[X_max_Phe] = 1.000000 ;
            p[X_max_S2O3] = 1.000000 ;
            p[d_SCN_SCN] = 0.009000 ;
            p[d_SCN_Phe] = 0.008845;
            p[d_SCN_S2O3] = 0.006627;
            p[d_Phe_SCN] = 0.007094;
            p[d_Phe_Phe] = 0.009000;
            p[d_Phe_S2O3] = 0.007127;
            p[d_S2O3_SCN] = 0.000187;
            p[d_S2O3_Phe] = 0.000288;
            p[d_S2O3_S2O3] = 0.009000;
            p[c_SCN] = 36.803757;
            p[c_Phe] = 52.886314;
            p[c_S2O3] = 9.075364;
            p[beta_SCN] = 19.516453;
            p[beta_Phe] = 5.058998;
            p[beta_S2O3] = 630.980520;
            p[gamma_SCN] = 7.504556;
            p[gamma_Phe] = 2.109228;
            p[gamma_S2O3] = 0.986614 ;
            p[A_SCN] = 52.301366;
            p[A_Phe] = 63.860668;
            p[A_S2O3] = 845.096672 ;
            p[B_SCN] = 84.647485 ;
            p[B_Phe] = 23.881614;
            p[B_S2O3] = 7280.455332 ;
            */
     //  for(int i=0;i<PARAM_N-4;i++) scanf("%lf",&p[i]);
     p[Pin]=-1;
     p[P_SCN]=-1;
     p[P_Phe]=-1;
     p[P_S2O3]=-1;

     std::mt19937_64 engine;
     std::uniform_real_distribution<double> range0(0, 1);
     std::uniform_real_distribution<double> range1(0, 10000);
     std::uniform_real_distribution<double> range2(6000, 8000);
     std::uniform_real_distribution<double> range3(100, 1000);
     std::uniform_real_distribution<double> range4(0,0.1);
     std::uniform_real_distribution<double> range5(1000,2000);

     /*for(size_t i = 0; i < 256; ++i)
       std::cout << range1(engine) << "\n";*/

     double answer;
     //while(1){
     engine.seed(count);

     //     p[alpha_SCN]  = range4(engine);
     p[c_SCN]      = range0(engine);

     //     p[alpha_Phe]  = range4(engine);
     p[c_Phe]      = range0(engine);

     // p[alpha_S2O3] = range4(engine);
     p[c_S2O3]     = range0(engine);

     p[A_SCN]      = range0(engine);
     p[B_SCN]      = range1(engine);

     p[A_Phe]      = range0(engine);
     p[B_Phe]      = range1(engine);

     p[A_S2O3]     = range0(engine);
     p[B_S2O3]     = range1(engine);

     p[beta_SCN]   = range1(engine);
     p[beta_Phe]   = range1(engine);
     p[beta_S2O3]  = range1(engine);

     p[gamma_SCN]  = range0(engine);
     p[gamma_Phe]  = range0(engine);
     p[gamma_S2O3] = range0(engine);

     p[d_SCN_SCN]   = range0(engine);
     p[d_Phe_Phe]   = range0(engine);
     p[d_S2O3_S2O3] = range0(engine);

     p[d_SCN_Phe]  = range0(engine);
     p[d_SCN_S2O3] = range0(engine);
     p[d_Phe_SCN]  = range0(engine);
     p[d_Phe_S2O3] = range0(engine);
     p[d_S2O3_SCN] = range0(engine);
     p[d_S2O3_Phe] = range0(engine);

     answer=/*target_func3sp(x,p)+target_funcSCN(x,p)+target_funcPhe(x,p)+target_funcS2O3(x,p)*/target_funcSCN(x,p);
     FILE *fp;
     sprintf(filename,"plotdataS3sp");
     fp=fopen(filename,"w");
     for(int i=0;i<9;i++){
          fprintf(fp,"%d %lf %lf %lf\n",i,speeddata3sp[i][0],speeddata3sp[i][1],speeddata3sp[i][2]);
          fprintf(fp,"%d %lf %lf %lf\n",i+1,speeddata3sp[i][0],speeddata3sp[i][1],speeddata3sp[i][2]);
     }
     fclose(fp);
     sprintf(filename,"plotdataSSCN");
     fp=fopen(filename,"w");
     for(int i=0;i<7;i++){
          fprintf(fp,"%d %lf %lf %lf\n",i,speeddataSCN[i][0],speeddataSCN[i][1],speeddataSCN[i][2]);
          fprintf(fp,"%d %lf %lf %lf\n",i+1,speeddataSCN[i][0],speeddataSCN[i][1],speeddataSCN[i][2]);
     }
     fclose(fp);
     sprintf(filename,"plotdataSPhe");
     fp=fopen(filename,"w");
     for(int i=0;i<9;i++){
          fprintf(fp,"%d %lf %lf %lf\n",i,speeddataPhe[i][0],speeddataPhe[i][1],speeddataPhe[i][2]);
          fprintf(fp,"%d %lf %lf %lf\n",i+1,speeddataPhe[i][0],speeddataPhe[i][1],speeddataPhe[i][2]);
     }
     fclose(fp);
     sprintf(filename,"plotdataSS2O3");
     fp=fopen(filename,"w");
     for(int i=0;i<10;i++){
          fprintf(fp,"%d %lf %lf %lf\n",i,speeddataS2O3[i][0],speeddataS2O3[i][1],speeddataS2O3[i][2]);
          fprintf(fp,"%d %lf %lf %lf\n",i+1,speeddataS2O3[i][0],speeddataS2O3[i][1],speeddataS2O3[i][2]);
     }
     fclose(fp);
     sprintf(filename,"plotdataX3sp");
     fp=fopen(filename,"w");
     for(int i=0;i<9;i++) fprintf(fp,"%d %lf %lf %lf\n",i,Xdata3sp[i][0],Xdata3sp[i][1],Xdata3sp[i][2]);
     fclose(fp);
     sprintf(filename,"plotdataXSCN");
     fp=fopen(filename,"w");
     for(int i=0;i<7;i++) fprintf(fp,"%d %lf %lf %lf\n",i,XdataSCN[i][0],XdataSCN[i][1],XdataSCN[i][2]);
     fclose(fp);
     sprintf(filename,"plotdataXPhe");
     fp=fopen(filename,"w");
     for(int i=0;i<9;i++) fprintf(fp,"%d %lf %lf %lf\n",i,XdataPhe[i][0],XdataPhe[i][1],XdataPhe[i][2]);
     fclose(fp);
     sprintf(filename,"plotdataXS2O3");
     fp=fopen(filename,"w");
     for(int i=0;i<10;i++) fprintf(fp,"%d %lf %lf %lf\n",i,XdataS2O3[i][0],XdataS2O3[i][1],XdataS2O3[i][2]);
     fclose(fp);

     show_param(p);

     if(flag==0){
          fprintf(stderr,"%d %lf\n",count,answer);
     }
     flag=0;
     //if(answer<1e4) break;

     count+=1;
     //}
     return 0;
}


