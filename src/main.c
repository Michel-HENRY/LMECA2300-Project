#include "vector.h"
#include "kernel.h"

int main(){
  Kernel kernel = Lucy;
  double h = 1.1;
  int DIM = 2;
  Vector* v1 = Vector_new(DIM);
  Vector* v2 = Vector_new(DIM);

  double x1[2] = {0,0};
  double x2[2] = {1,0};

  Vector_initialise(v1,x1);
  Vector_initialise(v2,x2);

  Vector_print(v1);
  Vector_print(v2);


  double eval = eval_kernel(v1,v2,h,kernel);
  printf("eval = %f\n", eval);

  Vector* dW = grad_kernel(v1,v2,h,kernel);
  Vector_print(dW);

  Vector_free(v1);
  Vector_free(v2);

  exit(0);
}
