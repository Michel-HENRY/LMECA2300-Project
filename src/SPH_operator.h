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
Vector* grad_P(Particle* pi, Kernel kernel);
Vector* lapl_u(Particle* pi, Kernel kernel);
Vector* lapl_u_shao(Particle* pi, Kernel kernel);
void kernel_validation();
#endif
