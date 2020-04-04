#ifndef KERNEL_H
#define KERNEL_H

#include "vector.h"
#include <math.h>

#define M_PI 3.14159265358979323846

typedef enum Kernel Kernel;

enum Kernel {Cubic,Lucy};

double eval_kernel(Vector* v1,Vector* v2, double h, Kernel kernel);
Vector* grad_kernel(Vector* v1, Vector* v2, double h, Kernel kernel);

#endif
