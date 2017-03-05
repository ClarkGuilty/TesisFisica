#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>

#define L 2.0
#define pi 3.1415
#define N 1024

double function1(double t);

int main(){

  int i;
  double t=0.0;
  double delx=L/(N-1);
  double * x;
  x=malloc(sizeof(double)*N);
  for(i=0;i<N;i++){
    x[i]=function1(t);
    t+=delx;
  }
  fftw_complex *in, *out;
  fftw_plan my_plan;
  in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N);
  for(i=0;i<N;i++){
    in[i] = x[i];
  }
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N);
  my_plan = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(my_plan);
  fftw_destroy_plan(my_plan);
  for(i=0;i<N;i++){
    printf("%f %f %f \n", creal(in[i]), creal(out[i]), cimag(out[i]));
  }
  fftw_free(in);
  fftw_free(out);

  return 0;
}

double function1(double var){
	return sin(2*M_PI*var/L);
}
