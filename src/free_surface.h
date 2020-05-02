#ifndef FREE_SURFACE_H
#define FREE_SURFACE_H

#include "vector.h"
#include "particle.h"
#include "kernel.h"

// ADD Cs in particle
Vector** get_force_surface(Particle** p, int n_p, Kernel kernel);
static double get_curvature(Particle* p, double norm_n, Kernel kernel);
double get_lapl_Cs(Particle*p, Kernel kernel);
double get_lapl_CsKGC(Particle*p, Kernel kernel, Vector** dWi);
static Vector* get_n(Particle*p, Kernel kernel);
static double get_smooth_Cs(Particle* p, Kernel kernel);

#endif
