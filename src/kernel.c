#include "kernel.h"

double eval_Cubic_kernel(double q, double h);
double eval_Lucy_kernel(double q, double h);
double eval_NewQuartic_kernel(double q, double h);
double eval_Quintic_kernel(double q, double h);

double eval_kernel(xy *p1, xy *p2, double kh, Kernel kernel) {
    double d = sqrt(squared(p1->x-p2->x) + squared(p1->y-p2->y));
    double h;
    if(kernel == Cubic) {
        h = kh/2;
        return eval_Cubic_kernel(d/h, h);
    }
	else if (kernel == Lucy) {
		h = kh;
		return eval_Lucy_kernel(d/h, h);
	}
	else if (kernel == NewQuartic) {
		h = kh/2;
		return eval_NewQuartic_kernel(d/h, h);
	}
	else if (kernel == Quintic) {
		h = kh/3;
		return eval_Quintic_kernel(d/h, h);
	}
}

double eval_Cubic_kernel(double q, double h) {
    double alpha = 15.0/(7.0*M_PI*pow(h, 2));
    if(q >= 0 && q <= 1) return alpha * (2.0/3 - q*q + q*q*q/2);
    else if(q > 1 && q <= 2) return alpha * (pow(2-q, 3)/6);
    else return 0;
}

double eval_Lucy_kernel(double q, double h) {
	double alpha = 5.0 / (M_PI*pow(h, 2));
	if (q >= 0 && q <= 1) return alpha * (1+3*q)*(1-q*q*q);
	else return 0;
}

double eval_NewQuartic_kernel(double q, double h) {
	double alpha = 15.0 / (7.0*M_PI*pow(h, 2));
	if (q >= 0 && q <= 2) return alpha * (2.0/3 - 9.0/8*q*q + 19.0/24* pow(q,3) -5.0/32*pow(q,4));
	else return 0;
}

double eval_Quintic_kernel(double q, double h) {
	double alpha = 7.0 / (478.0*M_PI*pow(h, 2));
	if (q >= 0 && q <= 1) return alpha * (pow(3-q,5)-6*pow(2-q,5)+15*pow(1-q,5));
	else if (q > 1 && q <= 2) return alpha *(pow(3-q, 5)-6*pow(2 - q, 5));
	else if (q > 2 && q <= 3) return alpha * pow(3 - q, 5);
	else return 0;
}

// double eval_Wendland_kernel(double q, double h) {
//     if(q <= 3) return (7/(144*M_PI*h*h)) * pow(2-2*q/3, 4) * (1+4*q/3);
//     else return 0;
// }

double derivative_Cubic_kernel(double q, double h);
double derivative_Lucy_kernel(double q, double h);
double derivative_NewQuartic_kernel(double q, double h);
double derivative_Quintic_kernel(double q, double h);

// everything here should be double checked because there are still problems with signs
xy* grad_kernel(xy* p1, xy* p2, double kh, Kernel kernel) {
     if(p1->x == p2->x && p1->y == p2->y) return xy_new(0,0);

	double d = sqrt(pow(p1->x - p2->x, 2) + pow(p1->y - p2->y, 2));
    double d_x = p1->x-p2->x;
    double d_y = p1->y-p2->y;

	double g;
	double h;
	if (kernel == Cubic) {
		h = kh / 2;
		g = derivative_Cubic_kernel(d / h, h);
	}
	else if (kernel == Lucy) {
		h = kh;
		g = derivative_Lucy_kernel(d / h, h);
	}
	else if (kernel == NewQuartic) {
		h = kh / 2;
		g = derivative_NewQuartic_kernel(d / h, h);
	}
	else if (kernel == Quintic) {
		h = kh / 3;
		g = derivative_Quintic_kernel(d / h, h);
	}

	double g_x = g*(d_x / (h*d));
	double g_y = g*(d_y / (h*d));

	return xy_new(g_x, g_y);
}

double derivative_kernel(double r, double h, Kernel kernel){
	if(kernel == Cubic){
		return derivative_Cubic_kernel(2*r/h, h/2);
	} else if(kernel == Lucy){
		return derivative_Lucy_kernel(r/h, h);
	} else if(kernel == NewQuartic){
		return derivative_NewQuartic_kernel(2*r/h, h/2);
	} else if(kernel == Quintic){
		return derivative_Quintic_kernel(3*r/h, h/3);
	}
	return 0;
}

double derivative_Cubic_kernel(double q,double h)
{
	 double alpha = 15.0/(7.0*M_PI*pow(h, 2));
	 double g;
	 if (q >= 0 && q <= 1)
		 g = -2 * q + 1.5*pow(q, 2);
	 else if (q > 1 && q <= 2)
		 g = -0.5*pow((2 - q), 2);
	 else
		 g = 0;
	 return alpha*g;
 }

 double derivative_Lucy_kernel(double q, double h)
 {
	 double alpha = 5.0/(M_PI*pow(h, 2));
	 double g;
	 if (q >= 0 && q <= 1)
		 g = 3 * (1 - pow(q, 2) - 4 * pow(q, 3));
	 else
		 g = 0;
	 return alpha*g;
 }

 double derivative_NewQuartic_kernel(double q, double h)
 {
	 double alpha = 15.0/(7.0*M_PI*pow(h,2));
	 double g;
	 if (q >= 0 && q <= 2)
		 g = (-9.0/4.0)*q+(19.0/8.0)*pow(q,2)-(5.0/8.0)*pow(q,3);
	 else
		 g = 0;
	 return alpha*g;
 }

 double derivative_Quintic_kernel(double q, double h)
 {
	 double alpha = 7.0 / (478.0*M_PI*pow(h, 2));
	 double g;
	 if (q >= 0 && q <= 1)
		 g = -5*pow(3-q,4)+30*pow(2-q,4)-75*pow(1-q,4);
	 else if (q > 1 && q <= 2)
		 g = -5 * pow(3 - q, 4) + 30 * pow(2 - q, 4);
	 else if (q > 2 && q <= 3)
		 g = -5 * pow(3 - q, 4);
	 else
		 g = 0;
	 return alpha*g;
 }
