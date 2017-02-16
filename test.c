#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(){
  float a=4.33;
  float b=3.7;
  printf("%f \n", a/b);
  printf("%d \n", (int) a/b);
  printf("%f \n", round(a/b));
  printf("%d \n", (int) round(a/b));
  return 0;
}
