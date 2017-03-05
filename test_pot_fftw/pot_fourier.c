#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>

#define Nx 1024
#define Nv 1024
#define L 2.0
#define L_min -1.0
#define V 2.0
#define V_min -1.0
#define pi 3.141592654
#define G 6.67408E-11
#define FLOAT double
#define T 1
#define skip 1
#define deltat 0.1

FLOAT gauss(FLOAT pos, FLOAT vel, FLOAT amp, FLOAT sigma);
FLOAT * densidad(FLOAT *fase);
int ndx(int fila, int column);
FLOAT * potfourier(FLOAT *rho);

int main(){

  int i,j,k;
  FLOAT L_max = L_min+L;
  FLOAT V_max = V_min+V;
  FLOAT delx = L/(Nx);
  FLOAT delv = V/(Nv-1);

  FLOAT *phase;
  phase = malloc(sizeof(FLOAT)*Nx*Nv);
  FLOAT *dens;
  dens=malloc(sizeof(FLOAT)*Nx);
  FLOAT *pot;
  pot=malloc(sizeof(FLOAT)*Nx);

  for(i=0;i<Nv;i++){
    for(j=0;j<Nx;j++){
      phase[ndx(i,j)]=gauss(L_min+j*delx, V_min+i*delv, 4.0, 0.08);
    }
  }

  dens=densidad(phase);

  pot=potfourier(dens);
  for(i=0;i<Nx;i++){
    printf("%f \n", pot[i]);
  }
}

FLOAT gauss(FLOAT pos, FLOAT vel, FLOAT amp, FLOAT sigma){
  return amp*exp(-(pow(pos,2)+pow(vel,2))/sigma);
}
FLOAT * densidad(FLOAT *fase){
  FLOAT * rho;
  rho=malloc(sizeof(FLOAT)*Nx);
  int i, j;
  for(i=0;i<Nx;i++){
    for(j=0;j<Nv;j++){
      rho[i]+=fase[ndx(j,i)];
    }
  }
  return rho;
}
int ndx(int fila, int column){
  return fila*Nx+column;
}
FLOAT * potfourier(FLOAT *rho){
  int i;
  FLOAT ks;
  FLOAT * res;
  res=malloc(sizeof(FLOAT)*Nx);
  fftw_complex *rho_in, *rho_out, *rho_fin;
  fftw_plan rho_plan;
  rho_in=malloc(sizeof(fftw_complex)*Nx);
  rho_out=malloc(sizeof(fftw_complex)*Nx);
  rho_fin=malloc(sizeof(fftw_complex)*Nx);
  rho_plan = fftw_plan_dft_1d(Nx, rho_in, rho_out, FFTW_FORWARD, FFTW_ESTIMATE);

  for(i=0;i<Nx;i++){
    rho_in[i]=rho[i];
  }
  fftw_execute(rho_plan);
  fftw_destroy_plan(rho_plan);
  rho_plan = fftw_plan_dft_1d(Nx, rho_out, rho_fin, FFTW_BACKWARD, FFTW_ESTIMATE);
  for(i=1;i<Nx;i++){
    ks=2*pi/L*(FLOAT)i;
    rho_out[i]=rho_out[i]/(Nx*(-pow(ks,2)));
  }
  fftw_execute(rho_plan);

  for(i=0;i<Nx;i++){
    res[i]=creal(rho_fin[i]/Nx);
  }
  return res;
}
