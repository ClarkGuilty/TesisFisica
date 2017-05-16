#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>

//-------------------------Constantes-------------------------//
#define Nx 1024
#define Nv 1024

#define L 1.0 //1.0
#define L_min -0.5 // -0.5
#define V 2.0
#define V_min -1.0

#define pi 3.141592654
#define FLOAT float

#define T 20
#define skip 1
#define deltat 0.1

//-------------------------Variables globales-------------------------//
FLOAT delx=L/(Nx);
FLOAT delv=V/(Nv);
FLOAT L_max = L_min+L;
FLOAT V_max = V_min+V;

int i,j,k;
int i_v_new, j_x_new;
FLOAT x, v, x_new, v_new;

FILE *phase_rela_dat, *phase_four_dat, *dens_dat, *acc_dat, *pot_dat, *vels_dat;
FLOAT *phase, *phase_new, *dens, *acc, *pot, *pot_temp, *vels;

FLOAT Kx;
FLOAT kx;
fftw_complex *rho_out, *rho_in, *rho_fin;
fftw_plan rho_plan;

char *method;

//-------------------------Declaracion de funciones-------------------------//
void gauss(FLOAT *arreglo, FLOAT *arreglo_new, FLOAT amp, FLOAT sigma);
void bullet(FLOAT *arreglo, FLOAT *arreglo_new, FLOAT amp1, FLOAT sigma1, FLOAT x1, FLOAT amp2, FLOAT sigma2, FLOAT x2);
void jeans(FLOAT *arreglo, FLOAT *arreglo_new, FLOAT rho, FLOAT amp, FLOAT sig, int n);
void densidad(FLOAT *fase, FLOAT *rho);
void potential(FLOAT *rho, FLOAT *Va, FLOAT *V_temp);
void potfourier_real(FLOAT *rho, FLOAT *res);
void acceleration(FLOAT *Va, FLOAT *aceleracion);
void update(FLOAT * fase, FLOAT * azz, FLOAT * phase_temp);
int ndx(int fila, int column);
void printINFO(int indice, FLOAT * density, FILE * dens_file, FLOAT * azz, FILE * azz_file, FLOAT * potencial, FILE * pot_file, FLOAT * fase, FILE * fase_file, FLOAT * speed, FILE * speed_file);
void printCONS(char *state);
FLOAT sinc(FLOAT x);
void check(FLOAT *arreglo);
void check2(fftw_complex *arreglo);
void dens_vel(FLOAT *fase, FLOAT *rho_v);
void RELAX();
void FOURIER();
//-------------------------Main-------------------------//
int main(){
  dens_dat=fopen("dens_dat.txt", "w");
  acc_dat=fopen("acc_dat.txt", "w");
  pot_dat=fopen("pot_dat.txt", "w");
  vels_dat=fopen("vels_dat.txt", "w");

  phase = malloc(sizeof(FLOAT)*Nx*Nv);
  phase_new = malloc(sizeof(FLOAT)*Nx*Nv);
  dens=malloc(sizeof(FLOAT)*Nx);
  acc=malloc(sizeof(FLOAT)*Nx);
  pot=malloc(sizeof(FLOAT)*Nx);
  pot_temp=malloc(sizeof(FLOAT)*Nx);
  vels=malloc(sizeof(FLOAT)*Nv);
  check(phase); check(phase_new); check(dens); check(acc); check(pot); check(pot_temp); check(vels);

  //gauss(phase, phase_new, 3, 0.08);
  //bullet(phase, phase_new, 5, 0.04, -0.45, 5, 0.04, 0.45);
  //jeans(phase, phase_new, 0.25, 0.05, 0.5, 2);
  jeans(phase, phase_new, 10.0, 0.05, 0.5, 2);

  //RELAX();
  FOURIER();

  return 0;
}

//-------------------------Funciones-------------------------//
void gauss(FLOAT *arreglo, FLOAT *arreglo_new, FLOAT amp, FLOAT sigma){
  for(i=0;i<Nv;i++){
    for(j=0;j<Nx;j++){
      FLOAT pos=L_min+j*delx+0.5*delx;
      FLOAT vel=V_min+i*delv+0.5*delv;
      arreglo[ndx(i,j)]=amp*exp(-(pow(pos,2)+pow(vel,2))/sigma);
    }
  }
}
void bullet(FLOAT *arreglo, FLOAT *arreglo_new, FLOAT amp1, FLOAT sigma1, FLOAT x1, FLOAT amp2, FLOAT sigma2, FLOAT x2){
  for(i=0;i<Nv;i++){
    for(j=0;j<Nx;j++){
      FLOAT pos=L_min+j*delx;
      FLOAT vel=V_min+i*delv;
      arreglo[ndx(i,j)]=amp1*exp(-(pow(pos-x1,2)+pow(vel,2))/sigma1)+amp2*exp(-(pow(pos-x2,2)+pow(vel,2))/sigma2);
    }
  }
}
void jeans(FLOAT *arreglo, FLOAT *arreglo_new, FLOAT rho, FLOAT amp, FLOAT sig, int n){
  FLOAT k0 = 2.0*pi/L;
  FLOAT k = n*k0;
  for(i=0;i<Nv;i++){
    for(j=0;j<Nx;j++){
      FLOAT pos=L_min+j*delx;
      FLOAT vel=V_min+i*delv;
      arreglo[ndx(i,j)]=rho/(pow(2*pi*sig*sig,0.5))*exp(-pow(vel,2)/(2*sig*sig))*(1+amp*cos(k*pos));
    }
  }
}
void densidad(FLOAT *fase, FLOAT *rho){
  for(i=0;i<Nx;i++){
    rho[i]=0.0;
    for(j=0;j<Nv;j++){
      rho[i]+=fase[ndx(j,i)]*delv;
    }
  }
}
void potential(FLOAT *rho, FLOAT *Va, FLOAT *V_temp){
  for(i=0;i<Nx;i++){
    V_temp[i]=0.0;
  }
  for(j=0;j<Nx*Nx;j++){
    for(i=1;i<Nx-1;i++){
      Va[i]=0.5*(V_temp[i-1]+V_temp[i+1] - rho[i]*delx*delx);
    }
    Va[0]=0.5*(V_temp[Nx-1]+V_temp[1] - rho[i]*delx*delx);
    Va[Nx-1]=0.5*(V_temp[Nx-2]+V_temp[0] - rho[i]*delx*delx);
    for(i=0;i<Nx;i++){
      V_temp[i]=Va[i];
    }
  }
}
void potfourier_real(FLOAT *rho, FLOAT *res){

  rho_in=fftw_malloc(sizeof(fftw_complex)*Nx);
  rho_out=fftw_malloc(sizeof(fftw_complex)*Nx);
  rho_fin=fftw_malloc(sizeof(fftw_complex)*Nx);
  check2(rho_in); check2(rho_fin); check2(rho_out);

  for(i=0;i<Nx;i++){
    rho_in[i]=rho[i];
  }
  rho_plan = fftw_plan_dft_1d(Nx, rho_in, rho_out, 1, FFTW_ESTIMATE);
  fftw_execute(rho_plan);
  fftw_destroy_plan(rho_plan);

  rho_out[0]=0.0;
  for(i=1;i<Nx;i++){
    kx=2*pi/L*(FLOAT)i;
    Kx=kx*sinc(0.5*kx*delx);
    rho_out[i]=rho_out[i]/(-pow(Kx,2));
  }
  rho_plan = fftw_plan_dft_1d(Nx, rho_out, rho_fin, -1, FFTW_ESTIMATE);
  fftw_execute(rho_plan);
  fftw_destroy_plan(rho_plan);
  for(i=0;i<Nx;i++){
    res[i]=rho_fin[i]/(Nx);
  }
}
void acceleration(FLOAT *Va, FLOAT *aceleracion){
  for(i=0;i<Nx-1;i++){
    aceleracion[i]=-(Va[i+1]-Va[i])/delx;
  }
  aceleracion[Nx-1]=-(Va[0]-Va[Nx-1])/delx;
}
void update(FLOAT * fase, FLOAT * azz, FLOAT * fase_new){
  for(i=0;i<Nv;i++){
    for(j=0;j<Nx;j++){
      fase_new[ndx(i,j)]=0.0;
    }
  }

  for(i=0;i<Nv;i++){
    for(j=0;j<Nx;j++){
      v=V_min+i*delv+0.5*delv;
      x=L_min+j*delx+0.5*delx;
      v_new=v+deltat*azz[j];
      x_new=x+deltat*v_new;
      i_v_new= (int) ((v_new-V_min)/delv);
      j_x_new= (int) ((x_new-L_min)/delx);

      if(i_v_new >= 0 && i_v_new < Nv){
        if(j_x_new < 0){
          j_x_new = (int) (Nx+j_x_new-1);
        }
        else if(j_x_new >= Nx){
          j_x_new = (int) (j_x_new%Nx);
        }
        fase_new[ndx(i_v_new, j_x_new)]=fase_new[ndx(i_v_new, j_x_new)]+fase[ndx(i,j)];
      }
    }
  }
  for(i=0;i<Nv;i++){
    for(j=0;j<Nx;j++){
      fase[ndx(i,j)]=fase_new[ndx(i,j)];
    }
  }
}
int ndx(int fila, int column){
  return fila*Nx+column;
}
void printINFO(int indice, FLOAT * density, FILE * dens_file, FLOAT * azz, FILE * azz_file, FLOAT * potencial, FILE * pot_file, FLOAT * fase, FILE * fase_file, FLOAT * speed, FILE * speed_file){
  if (indice%skip==0){
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
    for(j=0;j<Nv;j++){
      fprintf(speed_file, "%lf \n", speed[j]);
    }
  }
}
void printCONS(char *state){
  FILE *CONS;
  CONS=fopen("Constantes.txt", "w");
  fprintf(CONS, " Nx= %d\n Nv= %d\n L= %lf\n L_min= %lf\n V= %lf\n V_min= %lf\n T= %d\n skip= %d\n deltat= %lf \n Metodo= %s", Nx, Nv, L, L_min, V, V_min, T, skip, deltat, state);
}
FLOAT sinc(FLOAT x){
  if (x==0){
    return 1.0;
  }
  return sin(x)/x;
}
void check(FLOAT *arreglo){
  if(!arreglo){
    printf("Un arreglo no se definio correctamente \n");
    exit(0);
  }
}
void check2(fftw_complex *arreglo){
  if(!arreglo){
    printf("Un arreglo tipo fftw_complex no se definio correctamente \n");
    exit(0);
  }
}
void dens_vel(FLOAT *fase, FLOAT *rho_v){
  for(j=0;j<Nv;j++){
    rho_v[j]=0.0;
    for(i=0;i<Nx;i++){
      rho_v[j]+=fase[ndx(j,i)]*delx;
    }
  }
}
void RELAX(){
  phase_rela_dat=fopen("phase_rela_dat.txt", "w");
  method="Relaxation";
  printCONS(method);
  for(k=0;k<T;k++){
    printf("Paso %d/%d \n", k+1, T);
    densidad(phase, dens);
    potential(dens, pot, pot_temp);
    acceleration(pot, acc);

    dens_vel(phase, vels);

    printINFO(k, dens, dens_dat, acc, acc_dat, pot, pot_dat, phase, phase_rela_dat, vels, vels_dat);

    update(phase, acc, phase_new);
  }
}
void FOURIER(){
  phase_four_dat=fopen("phase_four_dat.txt", "w");
  method="Fourier";
  printCONS(method);
  for(k=0;k<T;k++){
    printf("Paso %d/%d \n", k+1, T);
    densidad(phase, dens);
    potfourier_real(dens, pot);
    acceleration(pot, acc);

    dens_vel(phase, vels);

    printINFO(k, dens, dens_dat, acc, acc_dat, pot, pot_dat, phase, phase_four_dat, vels, vels_dat);

    update(phase, acc, phase_new);
  }
}
