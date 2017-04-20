#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>

#define Nx 1024
#define Ny 1024
#define Nz 1024
#define Nv 1024
#define Nu 1024
#define Nw 1024

#define L 2.0 //1.0
#define L_min -1.0 // -0.5
#define V 2.0
#define V_min -1.0

#define pi 3.141592654
#define G 6.67408E-11
#define FLOAT double

#define T 30
#define skip 3
#define deltat 0.1

FLOAT delx=L/(Nx);
FLOAT delv=V/(Nv);
FLOAT dely=L/(Ny);
FLOAT delu=V/(Nu);
FLOAT delz=L/(Nz);
FLOAT delw=V/(Nw);
FLOAT L_max = L_min+L;
FLOAT V_max = V_min+V;
int i,j,k;

void gauss(FLOAT *arreglo, FLOAT *arreglo_new, FLOAT amp, FLOAT sigma);
void jeans(FLOAT *arreglo, FLOAT *arreglo_new, FLOAT rho, FLOAT amp, FLOAT sig, int n);
void densidad(FLOAT *fase, FLOAT *rho);
void potential(FLOAT *rho, FLOAT *Va, FLOAT *V_temp);
void potfourier_real(FLOAT *rho, FLOAT *res);
void acceleration(FLOAT *Va, FLOAT *aceleracion);
void update(FLOAT * fase, FLOAT * azz, FLOAT * phase_temp);
int ndx(int fila, int column);
void printINFO(int indice, FLOAT * density, FILE * dens_file, FLOAT * azz, FILE * azz_file, FLOAT * potencial, FILE * pot_file, FLOAT * fase, FILE * fase_file, FLOAT *potencial_rela, FILE * pot_rela_file, FLOAT *azz_rela, FILE *azz_rela_file);
void printCONS();
void check(FLOAT *arreglo);
FLOAT sinc(FLOAT x);

int main(){
  FILE *phase_dat, *dens_dat, *acc_dat, *pot_dat, *acc_rela_dat, *pot_rela_dat;
  phase_dat=fopen("phase_dat.txt", "w");
  dens_dat=fopen("dens_dat.txt", "w");
  acc_dat=fopen("acc_dat.txt", "w");
  pot_dat=fopen("pot_dat.txt", "w");
  acc_rela_dat=fopen("acc_rela_dat.txt", "w");
  pot_rela_dat=fopen("pot_rela_dat.txt", "w");

  FLOAT *phase, *phase_new, *dens, *acc, *pot, *acc_rela, *pot_temp, *pot_rela;
  phase = malloc(sizeof(FLOAT)*Nx*Nv);
  phase_new = malloc(sizeof(FLOAT)*Nx*Nv);
  dens=malloc(sizeof(FLOAT)*Nx);
  acc=malloc(sizeof(FLOAT)*Nx);
  pot=malloc(sizeof(FLOAT)*Nx);
  acc_rela=malloc(sizeof(FLOAT)*Nx);
  pot_temp=malloc(sizeof(FLOAT)*Nx);
  pot_rela=malloc(sizeof(FLOAT)*Nx);
  check(phase); check(phase_new); check(dens); check(acc); check(pot);
  check(acc_rela); check(pot_temp); check(pot_rela);

  printCONS();

  gauss(phase, phase_new, 4.0, 0.08);
  //jeans(phase, phase_new, 5.0, 0.01, 0.5, 2);

  for(k=0;k<T;k++){

    printf("Paso %d/%d \n", k, T);

    densidad(phase, dens);

    //potential(dens, pot_rela, pot_temp);
    potfourier_real(dens, pot);

    //acceleration(pot_rela, acc_rela);
    acceleration(pot, acc);

    printINFO(k, dens, dens_dat, acc, acc_dat, pot, pot_dat, phase, phase_dat, pot_rela, pot_rela_dat, acc_rela, acc_rela_dat);

    update(phase, acc_rela, phase_new);
  }

  return 0;
}

void gauss(FLOAT *arreglo, FLOAT *arreglo_new, FLOAT amp, FLOAT sigma){
  for(i=0;i<Nv;i++){
    for(j=0;j<Nx;j++){
      FLOAT pos=L_min+j*delx;
      FLOAT vel=V_min+i*delv;
      arreglo[ndx(i,j)]=amp*exp(-(pow(pos,2)+pow(vel,2))/sigma);
      arreglo_new[ndx(i,j)]=arreglo[ndx(i,j)];
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
      arreglo_new[ndx(i,j)]=arreglo[ndx(i,j)];
    }
  }
}
void densidad(FLOAT *fase, FLOAT *rho){
  for(i=0;i<Nx;i++){
    rho[i]=0.0;
    for(j=0;j<Nv;j++){
      rho[i]+=fase[ndx(j,i)];
    }
  }
}
void potential(FLOAT *rho, FLOAT *Va, FLOAT *V_temp){
  for(i=0;i<Nx;i++){
    V_temp[i]=0.0;
  }
  for(j=0;j<Nx*Nx/4;j++){
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
  FLOAT Kx;
  FLOAT kx;
  int ncx=(Nx/2+1);
  FLOAT *rho_in, *rho_fin;
  fftw_complex *rho_out;
  fftw_plan rho_plan;

  rho_in=fftw_malloc(sizeof(FLOAT)*Nx);
  rho_out=fftw_malloc(sizeof(fftw_complex)*ncx);
  rho_fin=fftw_malloc(sizeof(FLOAT)*Nx);

  for(i=0;i<Nx;i++){
    rho_in[i]=rho[i];
  }
  rho_plan = fftw_plan_dft_r2c_1d(Nx, rho_in, rho_out, FFTW_ESTIMATE);
  fftw_execute(rho_plan);
  fftw_destroy_plan(rho_plan);

  rho_out[0]=0.0;
  for(i=1;i<Nx;i++){
    kx=2*pi/L*(FLOAT)i;
    Kx=kx*sinc(0.5*kx*delx);
    rho_out[i]=rho_out[i]/(-pow(Kx,2)*10);
  }
  rho_plan = fftw_plan_dft_c2r_1d(Nx, rho_out, rho_fin, FFTW_ESTIMATE);
  fftw_execute(rho_plan);
  fftw_destroy_plan(rho_plan);

  for(i=0;i<Nx;i++){
    res[i]=rho_fin[i]/(Nx);
  }
  fftw_free(rho_in); fftw_free(rho_out); fftw_free(rho_fin);
}
void acceleration(FLOAT *Va, FLOAT *aceleracion){
  for(i=1;i<Nx;i++){
    aceleracion[i]=-(Va[i]-Va[i-1])/delx;
  }
  aceleracion[0]=-(Va[0]-Va[Nx-1])/delx;
}
void update(FLOAT * fase, FLOAT * azz, FLOAT * fase_new){
  int i_v_new, j_x_new;
  FLOAT x, v, x_new, v_new;
  for(i=0;i<Nv;i++){
    for(j=0;j<Nx;j++){
      fase_new[ndx(i,j)]=0.0;
    }
  }

  for(i=0;i<Nv;i++){
    for(j=0;j<Nx;j++){
      v=V_min+i*delv;
      x=L_min+j*delx;
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
void printINFO(int indice, FLOAT * density, FILE * dens_file, FLOAT * azz, FILE * azz_file, FLOAT * potencial, FILE * pot_file, FLOAT * fase, FILE * fase_file, FLOAT *potencial_rela, FILE * pot_rela_file, FLOAT *azz_rela, FILE *azz_rela_file){
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
      fprintf(pot_rela_file, "%lf \n", potencial_rela[j]);
      fprintf(azz_rela_file, "%lf \n", azz_rela[j]);
    }
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
void check(FLOAT *arreglo){
  if(!arreglo){
    printf("Un arreglo no se definio correctamente \n");
    exit(0);
  }
}
