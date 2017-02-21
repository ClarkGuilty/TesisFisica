#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define Nx 1024
#define Nv 1024
#define L 1.0
#define L_min -0.5
#define V 2.0
#define V_min -1.0
#define pi 3.141592654
#define G 6.67408E-11
#define FLOAT double
#define T 1

FLOAT gauss(FLOAT pos, FLOAT vel, FLOAT amp, FLOAT sigma);
FLOAT jeans(FLOAT pos, FLOAT vel, FLOAT rho, FLOAT amp, FLOAT sig, int n);
FLOAT * densidad(FLOAT *fase);
FLOAT * potential(FLOAT *rho, FLOAT *V_prev);
FLOAT * acceleration(FLOAT *Va);
FLOAT * update(FLOAT * fase, FLOAT * azz, FLOAT deltat);
int ndx(int fila, int column);

int main(){

  int i,j,k;
  FLOAT L_max = L_min+L;
  FLOAT V_max = V_min+V;
  FLOAT delx = L/(Nx-1);
  FLOAT delv = V/(Nv-1);
  FLOAT delt = 0.1;

  FLOAT *phase;
  phase = malloc(sizeof(FLOAT)*Nx*Nv);

  FLOAT *dens;
  dens=malloc(sizeof(FLOAT)*Nx);
  FLOAT *acc;
  acc=malloc(sizeof(FLOAT)*Nx);
  FLOAT *pot;
  pot=malloc(sizeof(FLOAT)*Nx);
  FLOAT *pot_new;
  pot_new=malloc(sizeof(FLOAT)*Nx);

  for(i=0;i<Nv;i++){
    for(j=0;j<Nx;j++){
      //phase[ndx(i,j)]=gauss(L_min+j*delx, V_min+i*delv, 4.0, 0.08);
      phase[ndx(i,j)]=jeans(L_min+j*delx, V_min+i*delv, 5.0, 0.01, 0.5, 2);
    }
  }

  for(i=0;i<Nx;i++){
    pot[i]=0.0;
    pot_new[i]=0.0;
  }

  FLOAT *phase_new;
  phase_new = malloc(sizeof(FLOAT)*Nx*Nv);

  for(i=0; i<Nv; i++){
    for(j=0; j<Nx; j++){
      printf("%lf  ", phase[ndx(i,j)]);
    }
    printf("\n");
  }

  for(k=0;k<T;k++){

    dens=densidad(phase);
    pot_new=potential(dens, pot);
    acc=acceleration(pot_new);

    phase_new=update(phase, acc, delt);

    for(i=0;i<Nv;i++){
      for(j=0;j<Nx;j++){
        phase[ndx(i,j)]=phase_new[ndx(i,j)];
      }
    }
    for(i=1;i<Nx;i++){
      pot[i]=pot_new[i];
    }

    if(k%1==0){
      for(i=0; i<Nv; i++){
        for(j=0; j<Nx; j++){
          printf("%lf  ", phase_new[ndx(i,j)]);
        }
        printf("\n");
      }
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

FLOAT * potential(FLOAT *rho, FLOAT *V_prev){

  FLOAT delx = L/(Nx-1);
  int i,j;
  FLOAT *Va;
  FLOAT *V_temp;
  Va=malloc(sizeof(FLOAT)*Nx);
  V_temp=malloc(sizeof(FLOAT)*Nx);
  for(i=0;i<Nx;i++){
    V_temp[i]=V_prev[i];
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

FLOAT * update(FLOAT * fase, FLOAT * azz, FLOAT deltat){
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
