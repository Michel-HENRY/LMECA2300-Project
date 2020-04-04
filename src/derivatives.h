#ifndef DERIVATIVES_H
#define DERIVATIVES_H

#include "utils.h"
#include "particle.h"
#include "kernel.h"
#include "SPH.h"
#include "consistency.h"

typedef double(*scalar_getter)(Particle* particle); // scalar getter function (e.g. pressure getter)
typedef xy*(*xy_getter)(Particle* particle); // xy getter function (e.g. velocity getter)

double compute_div(Particle * particle, xy_getter get, Kernel kernel, double kh);
void compute_grad(Particle * particle, scalar_getter get, Kernel kernel, double kh, xy* grad);
double compute_lapl(Particle *particle, scalar_getter get, Kernel kernel, double kh);

#endif
