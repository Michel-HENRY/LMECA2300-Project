#ifndef KERNEL_H
#define KERNEL_H
#include "utils.h"

typedef enum Kernel Kernel;

enum Kernel {Cubic,Lucy,NewQuartic,Quintic};

xy* grad_kernel(xy* p1, xy* p2, double kh, Kernel kernel);
double eval_kernel(xy *p1, xy *p2, double kh, Kernel kernel);
double derivative_kernel(double r, double h, Kernel kernel);

#endif
