#include "kernel.h"
double derivative_kernel(double d, double h, Kernel kernel);
double eval_Cubic_kernel(double R, double h);
double eval_Lucy_kernel(double R, double h);
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

// everything here should be double checked because there are still problems with signs
Vector* grad_kernel(Vector* v1, Vector* v2, double h, Kernel kernel) {
  if(equal(v1,v2)) return Vector_new(v1->DIM);

	double d = dist(v1,v2);
  Vector* dW = Vector_new(v1->DIM);
  double dWdr = derivative_kernel(d,h,kernel);
  for(int i = 0; i < v1->DIM ; i++){
    double abs_v1v2i = fabs(v1->X[i] - v2->X[i]);
    dW->X[i] = dWdr * abs_v1v2i;
  }

  return dW;
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


double derivative_kernel(double d, double h, Kernel kernel){
  double dWdr,drdx;
  if(kernel == Cubic){
    dWdr = derivative_Cubic_kernel(d/(2*h), h);
    drdx = 1/(2*h*d);
		return dWdr*drdx;
	}
  else if(kernel == Lucy){
    dWdr = derivative_Lucy_kernel(d/h, h);
    drdx = 1/(h*d);
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
     printf("R = %f\n", R);
     printf("dW = %f\n", dW);
     // dW += 3*pow(1-R,3);
     // dW += 3*pow(1-R,2)*(1+3*R);
   }
	 return alpha*dW;
 }
