#include <math.h>
#include <stdio.h>
void mk_color(int Nx,int Ny,double* u,double umin,double umax,unsigned char* image){
     int i,j;
     double color;
     for(i=0;i<=Nx;i++){
          for(j=0;j<=Ny;j++){
               color=(u[i+j*(Nx+1)]-umin)/(umax-umin);
               if(color < 0.0){
                    image[4*(i+j*(Nx+1))+0] = 0;
                    image[4*(i+j*(Nx+1))+1] = 0;
                    image[4*(i+j*(Nx+1))+2] = 255;
               }else if(color >= 1.0){
                    image[4*(i+j*(Nx+1))+0] = 255;
                    image[4*(i+j*(Nx+1))+1] = 0;
                    image[4*(i+j*(Nx+1))+2] = 0;
               }else if( (color >= 0.0) && (color < 0.250)){
                    image[4*(i+j*(Nx+1))+0] = 0;
                    image[4*(i+j*(Nx+1))+1] =round(255*4*color);
                    image[4*(i+j*(Nx+1))+2] = 255;
               }else if( (color >= 0.250) && (color < 0.50)){
                    image[4*(i+j*(Nx+1))+0] = 0;
                    image[4*(i+j*(Nx+1))+1] = 255;
                    image[4*(i+j*(Nx+1))+2] = round(255*(2-4*color));
               }else if( (color >= 0.50) && (color < 0.750)){
                    image[4*(i+j*(Nx+1))+0] = round(255*(4*color-2));
                    image[4*(i+j*(Nx+1))+1] = 255;
                    image[4*(i+j*(Nx+1))+2] = 0;
               }else if( (color >= 0.750) && (color < 1.0)){
                    image[4*(i+j*(Nx+1))+0] = 255;
                    image[4*(i+j*(Nx+1))+1] = round(255*4*(1-color));
                    image[4*(i+j*(Nx+1))+2] = 0;
               }
               image[4*(i+j*(Nx+1))+3] = 255;
          }
     }
}
