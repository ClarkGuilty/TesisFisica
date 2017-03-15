#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>

#define Nx 1024
#define Nv 1024
#define FLOAT double
#define pi 3.141592654
#define L 4.
#define A 0.3

int i,j,k;
FLOAT delx=L/Nx;

FLOAT * potential(FLOAT *rho);
FLOAT * potfourier(FLOAT *rho);
void init_rho(FLOAT *rho);
FLOAT * anal();
FLOAT sinc(FLOAT x);

int main(){
  FLOAT * dens;
  FLOAT * rela;
  FLOAT * fourier;
  FLOAT * sol_an;
  dens=malloc(sizeof(FLOAT)*Nx);
  rela=malloc(sizeof(FLOAT)*Nx);
  fourier=malloc(sizeof(FLOAT)*Nx);
  sol_an=malloc(sizeof(FLOAT)*Nx);
  init_rho(dens);
  rela=potential(dens);
  fourier=potfourier(dens);
  sol_an=anal();

  for (i=0;i<Nx;i++){
    printf("%f %f %f %f %f \n", i*delx, dens[i], rela[i], fourier[i], sol_an[i]);
  }
  return 0;
}

void init_rho(FLOAT *rho){
  for(i=0;i<Nx;i++){
    rho[i]=A*sin(2*pi*i*delx/L);
  }
}
FLOAT * potential(FLOAT *rho){
  FLOAT *Va;
  FLOAT *V_temp;
  Va=malloc(sizeof(FLOAT)*Nx);
  V_temp=malloc(sizeof(FLOAT)*Nx);
  for(i=0;i<Nx;i++){
    V_temp[i]=0.0;
  }
  for(j=0;j<Nx*Nx/4;j++){
    for(i=1;i<Nx-1;i++){
      Va[i]=0.5*(V_temp[i-1]+V_temp[i+1] - 0.005*rho[i]*delx*delx);
    }
    Va[0]=0.5*(V_temp[Nx-1]+V_temp[1] - 0.005*rho[i]*delx*delx);
    Va[Nx-1]=0.5*(V_temp[Nx-2]+V_temp[0] - 0.005*rho[i]*delx*delx);
    for(i=0;i<Nx;i++){
      V_temp[i]=Va[i];
    }
  }
  return Va;
}
FLOAT * potfourier(FLOAT *rho){
  FLOAT Kx;
  FLOAT kx;
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
  rho_out[0]=0.0;
  for(i=1;i<Nx;i++){
    kx=1/L*(FLOAT)i;
    Kx=kx*sinc(0.5*kx*delx);
    rho_out[i]=rho_out[i]/(-pow(Kx,2));
  }
  fftw_execute(rho_plan);

  for(i=0;i<Nx;i++){
    res[i]=creal(rho_fin[i]/Nx);
  }
  return res;
}
FLOAT * anal(){
  FLOAT * sol;
  sol=malloc(sizeof(FLOAT)*Nx);
  for(i=0;i<Nx;i++){
    sol[i]=-A*sin(2*pi*i*delx/L);
  }
  return sol;
}
FLOAT sinc(FLOAT x){
  if (x==0){
    return 1.0;
  }
  return sin(x)/x;
}
