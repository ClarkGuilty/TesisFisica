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
#define T 200

FLOAT gauss(FLOAT pos, FLOAT vel);
FLOAT * densidad(FLOAT **fase);
FLOAT * potential(FLOAT *rho, FLOAT *V_prev);
FLOAT * acceleration(FLOAT *Va);
FLOAT ** update(FLOAT ** fase, FLOAT * azz, FLOAT deltat);

int main(){

  int i,j,k;
  FLOAT L_max = L_min+L;
  FLOAT V_max = V_min+V;
  FLOAT delx = L/(Nx-1);
  FLOAT delv = V/(Nv-1);
  FLOAT delt = 0.1;


  FLOAT **phase;
  phase = malloc(sizeof(FLOAT *) * Nv);
  for(i=0;i<Nv;i++){
    phase[i]=malloc(sizeof(FLOAT)*Nx);
  }
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

  FLOAT ** phase_new;
  phase_new = malloc(sizeof(FLOAT *) * Nv);
  for(i=0;i<Nv;i++){
    phase_new[i]=malloc(sizeof(FLOAT) * Nx);
  }

  for(i=0; i<Nv; i++){
    for(j=0; j<Nx; j++){
      printf("%lf  ", phase[i][j]);
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
    for(i=1;i<Nx;i++){
      pot[i]=pot_new[i];
    }

    if(k%7==0){
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

FLOAT gauss(FLOAT pos, FLOAT vel){
  return 4*exp(-(pow(pos,2)+pow(vel,2))/0.08);
}

FLOAT * densidad(FLOAT **fase){
  FLOAT * rho;
  rho=malloc(sizeof(FLOAT)*Nx);
  int i, j;
  for(i=0;i<Nx;i++){
    for(j=0;j<Nv;j++){
      rho[i]+=fase[j][i];
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
  Va[0]=0; Va[Nx-1]=0;
  for(i=0;i<Nx;i++){
    V_temp[i]=V_prev[i];
  }
  for(j=0;j<Nx*Nx/4;j++){
    for(i=1;i<Nx-1;i++){
      Va[i]=0.5*(V_temp[i-1]+V_temp[i+1] - 0.005*rho[i]*delx*delx);
    }
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

FLOAT ** update(FLOAT ** fase, FLOAT * azz, FLOAT deltat){
  int i,j;
  int i_v_new, j_x_new;
  FLOAT x, v, x_new, v_new;
  FLOAT L_max = L_min+L;
  FLOAT V_max = V_min+V;
  FLOAT delx = L/(Nx-1);
  FLOAT delv = V/(Nv-1);
  FLOAT **phase_temp;
  phase_temp = malloc(sizeof(FLOAT *) * Nv);
  for(i=0;i<Nv;i++){
    phase_temp[i]=malloc(sizeof(FLOAT) * Nx);
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
