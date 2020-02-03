#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <random>
#include <iostream>
#include <float.h>
#include <vector>
#include <assert.h>
#include "bregsplit.h"


int main(int argc, const char **argv){
     FILE* fp;
     if(argc < 2){
          fprintf(stderr, "Usage: ./logistic lambda \n");
          exit(1);
     }
     int seed = 0;
     double p[PARAMSIZE]={0};
     double p2[PARAMSIZE]={0};
     char filename[128]={};
     double data[datasize]={};
     double data2[datasize]={};
     double dataplot[PLOTSIZE]={};
     double data2plot[PLOTSIZE]={};
     double grad[PARAMSIZE]={};
     int swch[PARAMSIZE]={0};
     double* null1=0;
     double* null2=0;
    // lambdaL1 = lambdaL2 = atof(argv[1]);
     seed=atoi(argv[1]);   

     fprintf(stderr,"%d\n",seed);
     
     std::mt19937_64 engine(seed);
     std::uniform_real_distribution<double> rangex0(0,1);
     std::uniform_real_distribution<double> rangea(0,3);
     std::uniform_real_distribution<double> rangek(0,3);
     std::uniform_real_distribution<double> rangeb(-0.5,0.5);
     std::uniform_real_distribution<double> rangeC(-0.5,0.5);
     /*std::uniform_int_distribution<> rangeA(0,2000);
     std::uniform_int_distribution<> rangeB(2000,4000);
     std::uniform_int_distribution<> rangeC(4000,6000);
     std::uniform_int_distribution<> rangeD(6000,8000);
     std::uniform_int_distribution<> rangeE(8000,10000);
*/
     std::uniform_real_distribution<double> range3(-0.10,0.10);
    
     std::normal_distribution<double> nrx0(0.1,0.02);
     std::normal_distribution<double> nra(1.0,0.2);
     std::normal_distribution<double> nrk(2.0,0.4);
     std::normal_distribution<double> nro(0,0.5);
     
     for(int i=0;i<PLOTSIZE;i++) {
        timedataplot[i]=end_time/PLOTSIZE*i;
     }
     for(int i=0;i<datasize;i++){
          timedata[i]=0.5+i;
     }



     //for(int i=0;i<datasize;i++) fprintf(stderr,"%lf\n",timedata[i]);
     quicksort(0,datasize,timedata);
     //for(int i=0;i<datasize;i++) fprintf(stderr,"%lf\n",timedata[i]);
     
    p[x0]=p2[x0]=0.1;
     p[a]=p2[a]=1.0;
     p[k]=p2[k]=2.0;
     p[c]=p2[c]=0;
     p[d]=p2[d]=0;
     p[e]=p2[e]=0;
     /*p[0]=1.0;
     p[1]=0.0;
     p[2]=-2.0;
     p[3]=0.0;
     p[4]=0.0;
     p[5]=0.0;
*/
     //////////////////////////////////////
     swch[x0]=0;
     swch[a]=0;
     swch[k]=0;
     swch[c]=0;
     swch[d]=0;
     swch[e]=0;
     /////////////////////////////////////
     ODErecord(p,Vardata,timedata,datasize,testfunc2);
    /* sprintf(filename,"calc_%f.dat",log10(lambdaL1));
     fp=fopen(filename,"w");
     for(int j=0;j<Num_V;j++){
          for(int i=0;i<datasize;i++) fprintf(fp,"%lf %lf\n",timedata[i],Vardata[i+j*datasize]);
     }
     fclose(fp);
*/
    // for(int i=0;i<datasize;i++) fprintf(stderr,"%lf\n",Vardata[i]);
     
     ////gosa///
    /*for(int i=0;i<Num_V;i++){
          for(int j=0;j<datasize;j++) Vardata[j+i*datasize]+=range3(engine)*fabs(Vardata[j+i*datasize]);
     }*/
     sprintf(filename,"calcc_%f.dat",log10(lambdaL1));
     fp=fopen(filename,"w");
     for(int j=0;j<Num_V;j++){
          for(int i=0;i<datasize;i++) fprintf(fp,"%lf %lf\n",timedata[i],Vardata[i+j*datasize]);
     }
     fclose(fp);


     
          p[0]=p2[0]=nrx0(engine);
          p[1]=p2[1]=nra(engine);
          p[2]=p2[2]=nrk(engine);
          p[3]=p2[3]=nro(engine);
          p[4]=p2[4]=nro(engine);
          p[5]=p2[5]=nro(engine);
     /*sprintf(filename,"calcsubp_%f.dat",log10(lambdaL1));
     fp=fopen(filename,"w");
     ODErecord(p,dataplot,timedataplot,PLOTSIZE,testfunc2);
     for(int j=0;j<Num_V;j++){
          for(int i=0;i<PLOTSIZE;i++) fprintf(fp,"%lf %lf\n",timedataplot[i],dataplot[i+j*PLOTSIZE]);
     }*/
    /* sprintf(filename,"p1b_%d.dat",seed);
     fp=fopen(filename,"w");
          for(int i=0;i<PARAMSIZE;i++) fprintf(fp,"p[%d] == %lf\n",i,p[i]);
     fclose(fp);

    //fprintf(stderr,"start lambda==%f\n",log10(lambdaL1));
     //fprintf(stderr,"L1\n");
     SplitBregman(p,swch,wrap1);
     //sprintf(filename,"p1_%f.dat",log10(lambdaL1));
     sprintf(filename,"p1a_%d.dat",seed);
     fp=fopen(filename,"w");
          for(int i=0;i<PARAMSIZE;i++) fprintf(fp,"p[%d] == %lf\n",i,p[i]);
     fclose(fp);
     //sprintf(filename,"calc1p_%f.dat",log10(lambdaL1));
     sprintf(filename,"calc1p_%d.dat",seed);
     fp=fopen(filename,"w");
     ODErecord(p,dataplot,timedataplot,PLOTSIZE,testfunc2);
     for(int j=0;j<Num_V;j++){
          for(int i=0;i<PLOTSIZE;i++) fprintf(fp,"%lf %lf\n",timedataplot[i],dataplot[i+j*PLOTSIZE]);
     }
     fclose(fp);
     //sprintf(filename,"calc1_%f.dat",log10(lambdaL1));
     sprintf(filename,"calc1_%d.dat",seed);
     fp=fopen(filename,"w");
     ODErecord(p,data,timedata,datasize,testfunc2);
     for(int j=0;j<Num_V;j++){
          for(int i=0;i<datasize;i++) fprintf(fp,"%lf %lf\n",timedata[i],data[i+j*datasize]);
     }
     fclose(fp);
     */
    //fprintf(stderr,"L2\n");
     controlCG(p2,null1,null2,swch,wrap2);


     //sprintf(filename,"p2_%f.dat",log10(lambdaL2));
     sprintf(filename,"p2a_%d.dat",seed);
     fp=fopen(filename,"w");
     for(int i=0;i<PARAMSIZE;i++) fprintf(fp,"p2[%d] == %lf\n",i,p2[i]);
     fclose(fp);

     //sprintf(filename,"calc2p_%f.dat",log10(lambdaL2));
     sprintf(filename,"calc2p_%d.dat",seed);
     fp=fopen(filename,"w");
     ODErecord(p2,data2plot,timedataplot,PLOTSIZE,testfunc2);
     for(int j=0;j<Num_V;j++){
          for(int i=0;i<PLOTSIZE;i++) fprintf(fp,"%lf %lf\n",timedataplot[i],data2plot[i+j*PLOTSIZE]);
     }
     fclose(fp);
      //sprintf(filename,"calc2_%f.dat",log10(lambdaL2));
      sprintf(filename,"calc2_%d.dat",seed);
     fp=fopen(filename,"w");
     ODErecord(p2,data2,timedata,datasize,testfunc2);
     for(int j=0;j<Num_V;j++){
          for(int i=0;i<datasize;i++) fprintf(fp,"%lf %lf\n",timedata[i],data2[i+j*datasize]);
     }
     fclose(fp);
     

     //fprintf(stderr,"end lambda==%f\n",log10(lambdaL1));
     return 0;
}

