#ifndef TIME_INTEGRATION_H
#define TIME_INTEGRATION_H

#include "vector.h"
#include "particle.h"
#include "SPH_operator.h"
#include "kernel.h"
#include "boundary.h"

double* rhs_mass_conservation(Particle** p, int n_p, Kernel kernel);
Vector** rhs_momentum_conservation(Particle** p, int n_p, Kernel kernel);
Vector** CSPM_rhs_momentum_conservation(Particle** p, int n_p, Kernel kernel);

void time_integration(Particle** p, int n_p,Kernel kernel, double dt, Edges* edges);
void time_integration_XSPH(Particle** p, int n_p, Kernel kernel, double dt, Edges* edges, double eta);
void time_integration_CSPM(Particle** p, int n_p, Kernel kernel, double dt, Edges* edges, double eta);
void time_integration_mass(Particle** p, int n_p,double* rhs_mass,double dt);
void time_integration_momentum(Particle** p, int n_p, Vector** rhs_momentum,double dt);
void time_integration_position(Particle** p, int n_p, double dt);

void update_pressureMod(Particle** p, int n_p, double rho_0, double g, double H);
void update_pressure(Particle** p, int n_p, double rho_0);
void XSPH_correction(Particle** p, int n_p, Kernel kernel, double eta);
void CSPM_density(Particle** p, int n_p, Kernel kernel);
Vector* CSPM_pressure(Particle* pi, Kernel kernel);
#endif
