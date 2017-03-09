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
#define T 30
#define skip 1
#define deltat 0.1

FLOAT gauss(FLOAT pos, FLOAT vel, FLOAT amp, FLOAT sigma);
FLOAT jeans(FLOAT pos, FLOAT vel, FLOAT rho, FLOAT amp, FLOAT sig, int n);
FLOAT * densidad(FLOAT *fase);
FLOAT * potential(FLOAT *rho);
FLOAT * potfourier(FLOAT *rho);
FLOAT * acceleration(FLOAT *Va);
FLOAT * update(FLOAT * fase, FLOAT * azz);
int ndx(int fila, int column);
void printINFO(FLOAT * density, FILE * dens_file, FLOAT * azz, FILE * azz_file, FLOAT * potencial, FILE * pot_file, FLOAT * fase, FILE * fase_file);
void printCONS();
FLOAT sinc(FLOAT x);

int main(){

  printCONS();
  int i,j,k;
  FLOAT L_max = L_min+L;
  FLOAT V_max = V_min+V;
  FLOAT delx = L/(Nx);
  FLOAT delv = V/(Nv-1);

  FLOAT *phase;
  phase = malloc(sizeof(FLOAT)*Nx*Nv);
  FILE * phase_dat;
  phase_dat=fopen("phase_dat.txt", "w");
  FLOAT *phase_new;
  phase_new = malloc(sizeof(FLOAT)*Nx*Nv);
  FLOAT *dens;
  dens=malloc(sizeof(FLOAT)*Nx);
  FILE * dens_dat;
  dens_dat=fopen("dens_dat.txt", "w");
  FLOAT *acc;
  acc=malloc(sizeof(FLOAT)*Nx);
  FILE * acc_dat;
  acc_dat=fopen("acc_dat.txt", "w");
  FLOAT *pot;
  pot=malloc(sizeof(FLOAT)*Nx);
  FILE * pot_dat;
  pot_dat=fopen("pot_dat.txt", "w");
  FLOAT *pot_new;
  pot_new=malloc(sizeof(FLOAT)*Nx);


  for(i=0;i<Nv;i++){
    for(j=0;j<Nx;j++){
      phase[ndx(i,j)]=gauss(L_min+j*delx, V_min+i*delv, 4.0, 0.08);
      //phase[ndx(i,j)]=jeans(L_min+j*delx, V_min+i*delv, 5.0, 0.01, 0.5, 2);
      phase_new[ndx(i,j)]=phase[ndx(i,j)];

    }
  }

  for(k=0;k<T;k++){

    dens=densidad(phase);
    //pot_new=potential(dens);
    pot_new=potfourier(dens);
    acc=acceleration(pot_new);

    if(k%skip==0){
        printINFO(dens, dens_dat, acc, acc_dat, pot_new, pot_dat, phase, phase_dat);
    }

    phase_new=update(phase, acc);

    for(i=0;i<Nv;i++){
      for(j=0;j<Nx;j++){
        phase[ndx(i,j)]=phase_new[ndx(i,j)];
      }
    }
    for(i=1;i<Nx;i++){
      pot[i]=pot_new[i];
    }
  }
  return 0;
}

FLOAT gauss(FLOAT pos, FLOAT vel, FLOAT amp, FLOAT sigma){
  return amp*exp(-(pow(pos,2)+pow(vel,2))/sigma);
}
FLOAT jeans(FLOAT pos, FLOAT vel, FLOAT rho, FLOAT amp, FLOAT sig, int n){
  FLOAT k0 = 2.0*pi/L;
  FLOAT k = n*k0;
  return rho/(pow(2*pi*sig*sig,0.5))*exp(-pow(vel,2)/(2*sig*sig))*(1+amp*cos(k*pos));
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
FLOAT * potential(FLOAT *rho){

  FLOAT delx = L/(Nx-1);
  int i,j;
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
  int i;
  FLOAT Kx;
  FLOAT kx;
  FLOAT delx = L/(Nx);
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
    kx=1/L*(FLOAT)i;
    Kx=kx*sinc(0.5*kx*delx);
    rho_out[i]=rho_out[i]/(-pow(Kx,2)*Nx);
  }
  fftw_execute(rho_plan);

  for(i=0;i<Nx;i++){
    res[i]=creal(rho_fin[i]/(4*Nx));
  }
  return res;
}
FLOAT * acceleration(FLOAT *Va){
  int i;
  FLOAT delx = L/(Nx-1);
  FLOAT *aceleracion;
  aceleracion=malloc(sizeof(FLOAT)*Nx);
  aceleracion[0]=0; aceleracion[Nx-1]=0;
  for(i=1;i<Nx-1;i++){
    aceleracion[i]=-(Va[i+1]-Va[i-1])/(2.0*delx);
  }
  return aceleracion;
}
FLOAT * update(FLOAT * fase, FLOAT * azz){
  int i,j;
  int i_v_new, j_x_new;
  FLOAT x, v, x_new, v_new;
  FLOAT L_max = L_min+L;
  FLOAT V_max = V_min+V;
  FLOAT delx = L/(Nx-1);
  FLOAT delv = V/(Nv-1);

  FLOAT * phase_temp;
  phase_temp = malloc(sizeof(FLOAT)*Nv*Nv);
  for(i=0;i<Nv;i++){
    for(j=0;j<Nx;j++){
      phase_temp[ndx(i,j)]=0.0;
    }
  }

  for(i=0;i<Nv;i++){
    for(j=0;j<Nx;j++){
      v=V_min+i*delv;
      x=L_min+j*delx;
      v_new=v+deltat*azz[j];
      x_new=x+deltat*v_new;
      if(v_new >= V_min && v_new <= V_max){ // && x_new >= L_min && x_new <= L_max){
        i_v_new= (int) round((v_new-V_min)/delv);
        j_x_new= (int) round((x_new-L_min)/delx);
        if(x_new < L_min){
          j_x_new = (int) Nx+j_x_new-1;
        }
        else if(x_new > L_max){
          j_x_new = (int) j_x_new%Nx;
        }
        phase_temp[ndx(i_v_new, j_x_new)]=phase_temp[ndx(i_v_new, j_x_new)]+fase[ndx(i,j)];
      }
    }
  }
  return phase_temp;
}
int ndx(int fila, int column){
  return fila*Nx+column;
}
void printINFO(FLOAT * density, FILE * dens_file, FLOAT * azz, FILE * azz_file, FLOAT * potencial, FILE * pot_file, FLOAT * fase, FILE * fase_file){
  int i, j;
  for(i=0;i<Nv;i++){
    for(j=0;j<Nx;j++){
      fprintf(fase_file, "%lf ", fase[ndx(i,j)]);
    }
    fprintf(fase_file, "\n");
  }
  for(j=0;j<Nx;j++){
    fprintf(dens_file, "%lf \n", density[j]);
    fprintf(azz_file, "%lf \n", azz[j]);
    fprintf(pot_file, "%lf \n", potencial[j]);
  }
}
void printCONS(){
  FILE *CONS;
  CONS=fopen("Constantes.txt", "w");
  fprintf(CONS, " Nx= %d\n Nv= %d\n L= %lf\n L_min= %lf\n V= %lf\n V_min= %lf\n T= %d\n skip= %d\n deltat= %lf", Nx, Nv, L, L_min, V, V_min, T, skip, deltat);
}
FLOAT sinc(FLOAT x){
  if (x==0){
    return 1.0;
  }
  return sin(x)/x;
}
