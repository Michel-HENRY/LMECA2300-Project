#ifndef TIME_INTEGRATION_H
#define TIME_INTEGRATION_H

#include "vector.h"
#include "particle.h"
#include "SPH_operator.h"
#include "kernel.h"
#include "boundary.h"
#include "free_surface.h"

void print_rhoMax(Particle** p, int n_p);
void print_rhoMin(Particle** p, int n_p);
// For initialisation
void density(Particle** p, int n_p, Kernel kernel);

// Fot time integraiton
void time_integration(Particle** p, int n_p, Kernel kernel, double dt, Edges* edges);
void time_integration_XSPH(Particle** p, int n_p, Kernel kernel, double dt, Edges* edges, double eta);
void time_integration_CSPM(Particle** p, int n_p, Kernel kernel, double dt, Edges* edges);
void time_integration_mass(Particle** p, int n_p,double* rhs_mass,double dt);
void time_integration_momentum(Particle** p, int n_p, Vector** rhs_momentum,double dt);
void time_integration_position(Particle** p, int n_p, double dt);


double* rhs_mass_conservation(Particle** p, int n_p, Kernel kernel);
Vector** rhs_momentum_conservation(Particle** p, int n_p, Kernel kernel);
Vector** get_dP(Particle** p, int n_p, Kernel kernel);

// Correction Methods
void XSPH_correction(Particle** p, int n_p, Kernel kernel, double eta);
double* CSPM_div(Particle** p, int n_p, Kernel kernel);
void CSPM_density(Particle** p, int n_p, Kernel kernel);
void CSPM_pressure(Particle** p, int n_p, Kernel kernel, Vector** dP);

#endif
