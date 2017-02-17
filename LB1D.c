#include <stdio.h>
#include <stdlib.h>
#include <math.h>

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

float gauss(float pos, float vel);
float * densidad(float **fase);
float * potential(float *rho, float *V_prev);
float * acceleration(float *Va);
float ** update(float ** fase, float * azz, float deltat);

int main(){

  int i,j,k;
  float L_max = L_min+L;
  float V_max = V_min+V;
  float delx = L/(Nx-1);
  float delv = V/(Nv-1);
  float delt = 0.01;


  float **phase;
  phase = malloc(sizeof(float *) * Nv);
  for(i=0;i<Nv;i++){
    phase[i]=malloc(sizeof(float)*Nx);
  }
  float *dens;
  dens=malloc(sizeof(float)*Nx);
  float *acc;
  acc=malloc(sizeof(float)*Nx);
  float *pot;
  pot=malloc(sizeof(float)*Nx);
  float *pot_new;
  pot_new=malloc(sizeof(float)*Nx);

  for(i=0;i<Nv;i++){
    for(j=0;j<Nx;j++){
      phase[i][j]=gauss(L_min+j*delx, V_min+i*delv);
    }
  }

  for(i=0;i<Nx;i++){
    pot[i]=0.0;
    pot_new[i]=0.0;
  }

  /*for(i=0;i<Nx;i++){
    printf("%f %f %f \n", dens[i], pot[i], acc[i]);
  }*/

  float ** phase_new;
  phase_new = malloc(sizeof(float *) * Nv);
  for(i=0;i<Nv;i++){
    phase_new[i]=malloc(sizeof(float) * Nx);
  }

  for(i=0; i<Nv; i++){
    for(j=0; j<Nx; j++){
      printf("%f  ", phase[i][j]);
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
        phase[i][j]=phase_new[i][j];
      }
    }
    for(i=0;i<Nx;i++){
      pot[i]=pot_new[i];
    }

    if(k%1==0){
      for(i=0; i<Nv; i++){
        for(j=0; j<Nx; j++){
          printf("%f  ", phase_new[i][j]);
        }
        printf("\n");
      }
    }
  }
  return 0;
}

float gauss(float pos, float vel){
  return 4*exp(-(pow(pos,2)+pow(vel,2))/0.08);
}

float * densidad(float **fase){
  float * rho;
  rho=malloc(sizeof(float)*Nx);
  int i, j;
  for(i=0;i<Nx;i++){
    for(j=0;j<Nv;j++){
      rho[i]+=fase[j][i];
    }
  }
  return rho;
}

float * potential(float *rho, float *V_prev){

  float delx = L/(Nx-1);
  int i,j;
  float *Va;
  float *V_temp;
  Va=malloc(sizeof(float)*Nx);
  V_temp=malloc(sizeof(float)*Nx);
  Va[0]=0; Va[Nx-1]=0;
  for(i=0;i<Nx;i++){
    V_temp[i]=V_prev[i];
  }
  for(j=0;j<2*Nx;j++){
    for(i=1;i<Nx-1;i++){
      Va[i]=0.5*(V_temp[i-1]+V_temp[i+1] - rho[i]*delx*delx);
    }
    for(i=0;i<Nx;i++){
      V_temp[i]=Va[i];
    }
  }
  return Va;
}

float * acceleration(float *Va){
  int i;
  float delx = L/(Nx-1);
  float *aceleracion;
  aceleracion=malloc(sizeof(float)*Nx);
  aceleracion[0]=0; aceleracion[Nx-1]=0;
  for(i=1;i<Nx-1;i++){
    aceleracion[i]=-(Va[i+1]-Va[i-1])/(2.0*delx);
  }
  return aceleracion;
}

float ** update(float ** fase, float * azz, float deltat){
  int i,j;
  int i_v_new, j_x_new;
  float x, v, x_new, v_new;
  float L_max = L_min+L;
  float V_max = V_min+V;
  float delx = L/(Nx-1);
  float delv = V/(Nv-1);
  float **phase_temp;
  phase_temp = malloc(sizeof(float *) * Nv);
  for(i=0;i<Nv;i++){
    phase_temp[i]=malloc(sizeof(float) * Nx);
  }
  for(i=0;i<Nv;i++){
    for(j=0;j<Nx;j++){
      phase_temp[i][j]=0.0;
    }
  }

  for(i=0;i<Nv;i++){
    for(j=0;j<Nx;j++){
      v=V_min+i*delv;
      x=L_min+j*delx;
      v_new=v+deltat*azz[j];
      x_new=x+deltat*v_new;
      if(v_new > V_min && v_new < V_max && x_new > L_min && x_new < L_max){
        i_v_new= (int) round((v_new-V_min)/delv);
        j_x_new= (int) round((x_new-L_min)/delx);
        phase_temp[i_v_new][j_x_new]=phase_temp[i_v_new][j_x_new]+fase[i][j];
      }
    }
  }
  return phase_temp;
}
