#include "kernel.h"

double eval_Cubic_kernel(double R, double h);
double eval_Lucy_kernel(double R, double h);
double derivative_kernel(Vector* v1, Vector* v2, double h, Kernel kernel, int axis);
double derivative_Cubic_kernel(double q, double h);
double derivative_Lucy_kernel(double q, double h);


double eval_kernel(Vector* v1, Vector* v2, double h, Kernel kernel) {
  double d = dist(v1,v2);
  double R;
  double W = 0;
  if(kernel == Cubic) {
      R = d/(2*h);
      W =  eval_Cubic_kernel(R, h);
  }
	else if (kernel == Lucy) {
		R = d/h;
		W =  eval_Lucy_kernel(R, h);
	}
  return W;
}

double eval_Cubic_kernel(double R, double h) {
  double alpha = 15.0/(7.0*M_PI*h*h);
  double eps = 1e-8;
  if(R >= -eps && R <= 1-eps) return alpha*(2.0/3.0 - pow(R,2) + pow(R,3)/2);
  else if(R > 1-eps && R <= 2+eps) return alpha*(pow(2-R, 3)/6);
  else return 0;
}

double eval_Lucy_kernel(double R, double h) {
 double alpha = 5.0 / (M_PI*pow(h, 2));
  double eps = 1e-8;
 if (R >= -eps && R <= 1+eps) return alpha*(1+3*R)*(1-pow(R,3));
 else return 0;
}

Vector* grad_kernel(Vector* v1, Vector* v2, double h, Kernel kernel) {
  if(equal(v1,v2)){
    return Vector_new(v1->DIM);
  }
  Vector* dW = Vector_new(v1->DIM);
  for(int i = 0; i < v1->DIM ; i++){
    dW->X[i] = derivative_kernel(v1, v2, h, kernel, i);
  }
  return dW;
}


double derivative_kernel(Vector* v1, Vector* v2, double h, Kernel kernel, int axis){
  double d = dist(v1,v2);
  double dWdr,drdx;
  if(kernel == Cubic){
    dWdr = derivative_Cubic_kernel(d/(2*h), h);
    drdx = (v1->X[axis] - v2->X[axis])/(d*2*h);
		return dWdr*drdx;
	}
  else if(kernel == Lucy){
    dWdr = derivative_Lucy_kernel(d/h, h);
    drdx = (v1->X[axis] - v2->X[axis])/(d*h);
		return dWdr*drdx;
	}
	return 0;
}

double derivative_Cubic_kernel(double R,double h){
	 double alpha = 15.0/(7.0*M_PI*pow(h, 2));
	 double dW = 0;
   double eps = 1e-8;
	 if (R >= -eps && R <= 1-eps) dW = -2*R + 1.5*pow(R, 2);
	 else if (R > 1-eps && R <= 2+eps) dW = -0.5*pow(2-R, 2);
	 return alpha*dW;
}

 double derivative_Lucy_kernel(double R, double h){
	 double alpha = 5.0/(M_PI*pow(h, 2));
	 double dW = 0;
   double eps = 1e-8;
	 if (R >= -eps && R <= 1+eps){
     dW = -12*R+24*pow(R,2)-12*pow(R,3);
   }
	 return alpha*dW;
 }

void kernel_validation(){
  Kernel kernel = Lucy;
  double h = 2.0;
  int DIM = 2;

  Vector* v1 = Vector_new(DIM);
  Vector* v2 = Vector_new(DIM);
  double x1[2] = {0,0};
  double x2[2] = {0,0};
  Vector_initialise(v1,x1);
  Vector_initialise(v2,x2);
  Vector_print(v1);
  Vector_print(v2);


  double eval = eval_kernel(v1,v2,h,kernel);
  Vector* dW = grad_kernel(v1,v2,h,kernel);
  printf("eval = %f\n", eval);
  Vector_print(dW);

  Vector_free(v1);
  Vector_free(v2);
  Vector_free(dW);

  exit(0);
}
