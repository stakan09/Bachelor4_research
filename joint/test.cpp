#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <random>
#include <iostream>
#include <float.h>


double func(double x){
     //return x+pow(x,4)/1000.;
     return x*x;
}
// (1 - x + x**2 - x**3 + x**4 - x**5)/(1 - x**2 + x**4 - x**6)
int main(){
     double x;
     std::mt19937_64 engine;
     std::uniform_real_distribution<double> range0(-2,2);
     engine.seed(1);
     FILE* fp;
     fp=fopen("testdata.txt","w");
     for(int i=0;i<20;i++){
          x=range0(engine);
          fprintf(fp,"{%lf, %lf},\n",x,func(x));
     }
     fclose(fp);
}
