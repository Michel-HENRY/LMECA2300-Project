#ifndef SPH_H
#define SPH_H
#include "particle.h"
#include "print_particules.h"
#include "kernel.h"
#include "derivatives.h"
#include "consistency.h"

typedef struct Boundary Boundary;
typedef struct Setup Setup;
typedef struct Residual Residual;
typedef void(*update_positions)(Grid* grid, Particle** particles, Particle_derivatives** particles_derivatives, Residual** residuals, int n_p, Setup* setup);
typedef enum Free_surface_detection Free_surface_detection;

//CSF       : if  ||n_i|| > threshold => particle i belongs to interface
//DIVERGENCE: if div(pos) < threshold => particle i belongs to interface
//NONE      : no free surface detection
enum Free_surface_detection {CSF,DIVERGENCE,NONE};

struct Setup{
	int itermax;
	double timestep;
	double kh;
	Verlet* verlet;
	Kernel kernel;
	Free_surface_detection free_surface_detection;  // Strategy to estimate if a particle should be considered on the free surface or not
	double interface_threshold; // Threshold for the detection of particles on the free surface
	double XSPH_epsilon; // Parameter between 0 and 1 multiplying the XSPH correction;0 if no correction wanted
	bool gravity;
};

struct Residual {
	double mass_eq;
	double momentum_x_eq;
	double momentum_y_eq;
};

struct Boundary{
  double xleft;
  double xright;
  double ybottom;
  double ytop;
	double CR;
  double CF;
};

// Setup* Setup_new(int iter, double timestep,double kh, Verlet* verlet, Kernel kernel,Free_surface_detection free_surface_detection,double interface_threshold, double XSPH_epsilon);
Setup* Setup_new(int iter, double timestep,double kh, Verlet* verlet, Kernel kernel,Free_surface_detection free_surface_detection,double interface_threshold, double XSPH_epsilon, bool gravity);
void Setup_free(Setup* setup);

Residual* Residual_new();
void free_Residuals(Residual** residuals,int N);

void simulate(Grid* grid, Particle** particles, Particle_derivatives** particles_derivatives, Residual** residuals, int n_p, update_positions update_positions, Setup* setup, Animation* animation);
void simulate_boundary(Grid* grid, Particle** particles, Particle_derivatives** particles_derivatives, Residual** residuals, int n_p, update_positions update_positions, Setup* setup, Animation* animation, Boundary* boundary);
void random_moves(Grid* grid, Particle** particles, Particle_derivatives** particles_derivatives, Residual** residuals, int n_p, Setup* setup);
void update_positions_Colagrossi(Grid* grid, Particle** particles, Particle_derivatives** particles_derivatives, Residual** residuals, int n_p, Setup* setup);
void update_positions_seminar_5(Grid* grid, Particle** particles, Particle_derivatives** particles_derivatives, Residual** residuals, int n_p, Setup* setup);
void update_positions_ellipse(Grid* grid, Particle** particles, Particle_derivatives** particles_derivatives, Residual** residuals, int n_p, Setup* setup);
void update_positions_test_static_bubble(Grid* grid, Particle** particles, Particle_derivatives** particles_derivatives, Residual** residuals, int n_p, Setup* setup);

void compute_Cs(Particle *particle, Kernel kernel, double kh);
void assemble_residual_NS(Particle* particle, Particle_derivatives* Particle_derivatives, Residual* residual, Setup* setup);
void time_integrate(Particle* particle, Residual* residual, double delta_t);
void time_integrate_CSPM(Particle* particle, Particle_derivatives* dpi, Residual* residual, Setup* setup);

double compute_curvature(Particle *particle, Setup *setup, double epsilon);
void compute_XSPH_correction(Particle *particle, Kernel kernel, double kh,double epsilon);

void assemble_residual_NS_test(Particle* particle, Particle_derivatives* particle_derivatives, Residual* residual, double radius_circle, Setup* setup);

double compute_admissible_dt(double safety_param, double h_p, double c_0, double rho_0, double mu, double sigma);

void compute_normal(Particle *particle, Particle_derivatives* particle_derivatives);

Boundary* Boundary_new(double xleft, double xright, double ybottom, double ytop, double CR, double CF);
void Boundary_free(Boundary* boundary);
void reflective_boundary(Particle** p, int n_p, Boundary* boundary,double Rp);

#endif
