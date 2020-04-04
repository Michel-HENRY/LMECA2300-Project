#ifndef CONSISTENCY_H
#define CONSISTENCY_H
#include "particle.h"
#include "print_particules.h"
#include "kernel.h"
#include "derivatives.h"
#include "SPH.h"

void get_M0(Particle** p, int n_p, double kh, Kernel kernel);
void get_M1(Particle** p, int n_p, double kh, Kernel kernel);
double get_M0_local(Particle* pi, double kh, Kernel kernel);
double get_M1_local(Particle* pi, double kh, Kernel kernel);
void correct_grad(xy *current_grad, Particle *p, double kh, Kernel kernel);
void correct_grad_local(xy *current_grad, Particle *p, Particle *j, double kh, Kernel kernel);
void density_correction_MLS(Particle** p, int n_p, double kh, Kernel kernel);
double density_correction_MLS_local(Particle* pi, double kh, Kernel kernel);
double get_W_MLS(Particle* pi, Particle* pj, double kh, Kernel kernel, double* beta);
void get_A(Particle* pi, double kh, Kernel kernel, double A[3][3]);
void get_beta(double A[3][3], double* beta);
void Corrective_Smoothed_Particle_Method(Particle *p,Particle_derivatives *dp, double kh, Kernel kernel);
#endif
