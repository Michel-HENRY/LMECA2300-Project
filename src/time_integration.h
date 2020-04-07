#ifndef TIME_INTEGRATION_H
#define TIME_INTEGRATION_H

#include "vector.h"
#include "particle.h"
#include "SPH_operator.h"
#include "kernel.h"

double* rhs_mass_conservation(Particle** p, int n_p, Kernel kernel);
Vector** rhs_momentum_conservation(Particle** p, int n_p, Kernel kernel);

void time_integration(Particle** p, int n_p,Kernel kernel, double dt);
void time_integration_mass(Particle** p, int n_p,double* rhs_mass,double dt);
void time_integration_momentum(Particle** p, int n_p, Vector** rhs_momentum,double dt);
void time_integration_position(Particle** p, int n_p, double dt);

#endif
