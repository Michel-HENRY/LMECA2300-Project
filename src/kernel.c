#include "kernel.h"

double eval_Cubic_kernel(double R, double h);
double eval_Lucy_kernel(double R, double h);
double derivative_kernel(Vector* v1, Vector* v2, double h, Kernel kernel, int axis);
double derivative_Cubic_kernel(double q, double h);
double derivative_Lucy_kernel(double q, double h);


double eval_kernel(Vector* v1, Vector* v2, double h, Kernel kernel) {
  double R = dist(v1,v2);
  double W = 0;
  if(kernel == Cubic) {
      W =  eval_Cubic_kernel(R, h);
  }
	else if (kernel == Lucy) {
		W =  eval_Lucy_kernel(R, h);
	}
  else if (kernel == Quartic){
    W = eval_Quartic_kernel(R,h);
  }
  return W;
}

static double eval_Cubic_kernel(double R, double h) {
  double alpha = 40.0/(7.0*M_PI*h*h);
  double q = R/h;
  double eps = 1e-8;
  if(q >= eps && q <= 1/2-eps) return alpha*(6*(pow(q,3) - pow(q,2)) + 1) ;
  else if(q > 1/2-eps && q <= 1-eps) return alpha*(2*pow(1-q,3));
  else return 0;
}

static double eval_Lucy_kernel(double R, double h) {
 double alpha = 5.0 / (M_PI*pow(h, 2));
 double q = R/h;
  double eps = 1e-8;
 if (q >= eps && q <= 1-eps) return alpha*(1+3*q)*pow(1-q,3);
 else return 0;
}

static double eval_Quartic_kernel(double R, double h){
  double alpha = 15.0/(7.0*M_PI*h*h);
  double q = 2*R/h;
  double eps = 1e-8;
  if(q >= eps && q <= 2-eps) return alpha*(2.0/3.0 - 9.0/8.0*q*q + 19/24.0 * pow(q,3) - 5/32.0 * pow(q,4));
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


static double derivative_kernel(Vector* v1, Vector* v2, double h, Kernel kernel, int axis){
  double R = dist(v1,v2);
  double dWdr,drdx = 1;
  double eta = 0.01*h*h;
  drdx = (v1->X[axis] - v2->X[axis])/(R+eta);
  if(kernel == Cubic){
    dWdr = derivative_Cubic_kernel(R, h);
		return dWdr*drdx;
	}
  else if(kernel == Lucy){
    dWdr = derivative_Lucy_kernel(R, h);
		return dWdr*drdx;
	}
  else if(kernel = Quartic){
    dWdr = derivative_Quartic_kernel(R,h);
    return dWdr*drdx;
  }
	return 0;
}

static double derivative_Cubic_kernel(double R,double h){
   double alpha = 40/(7*M_PI*h*h);
   double dW = 0;
   double q = R/h;
   double eps = 1e-8;
   if(q > -eps && q <= 0.5+eps) dW = 6*(3*q*q - 2*q);
   else if (q > 0.5+eps && q < 1+eps) dW = -6*pow(1-q,2);
   return alpha*dW/h;
}

static double derivative_Lucy_kernel(double R, double h){
	 double alpha = 5.0/(M_PI*pow(h, 2));
	 double dW = 0;
   double eps = 1e-8;
   double q = R/h;
	 if (q >= eps && q <= 1-eps){
     dW = -12*q*(1 - 2*q + q*q);
   }
	 return alpha*dW/h;
 }
static double derivative_Quartic_kernel(double R, double h){
  double alpha = 15.0/(7.0*M_PI*h*h);
  double dW = 0;
  double eps = 1e-8;
  double q = 2*R/h;
  if(q >= eps && q <= 2-eps) dW = q*(-18.0/8.0 + 57.0*q/24.0 - 20.0/32.0*q*q);
  return 2*alpha*dW/h;
}
void kernel_validation(){
  Kernel kernel = Cubic;
  double h = 1;
  int DIM = 2;

  Vector* v1 = Vector_new(DIM);
  Vector* v2 = Vector_new(DIM);
  double x1[2] = {0,0};
  double x2[2] = {0.5,0};
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
