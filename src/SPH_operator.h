// This file contains the following functions :
//  - get_div() : which compute the divergence of a function in 2D
//  - get_grad() : which compute the gradient of scalar function in 2D
//  - get_lapl() : which compute the laplacian of a function in 2D
#ifndef SPH_OPERATOR_H
#define SPH_OPERATOR_H

#include "vector.h"
#include "kernel.h"
#include "particle.h"

double div_u(Particle* pi, Kernel kernel);
double div_local(Vector* fi, Vector* fj, Vector* dWij, double mj);
Vector* grad_P(Particle* pi, Kernel kernel);
Vector* grad_local(double fi, double fj, Vector* dWij, double mj, double rhoi,double rhoj);
Vector* lapl_u(Particle* pi, Kernel kernel);
Vector* lapl_u_shao(Particle* pi, Kernel kernel);
Vector* lapl_local(Vector* fi, Vector* fj, Vector* xi, Vector* xj, Vector* dWij, double Vj, double h);
void kernel_validation();
#endif
