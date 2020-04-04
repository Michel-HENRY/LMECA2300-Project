#include "print_particules.h"
#include "particle.h"
#include "SPH.h"
#include "derivatives.h"
#include <math.h>
#include "kernel.h"
#include "consistency.h"

//#include "crtdbg.h" // for memory leak detection; comment if you're on Linux

void script_circle_to_ellipse();
void dam_break();
void box();

int main() {
	// _CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF); // comment if on Linux
	 script_circle_to_ellipse();
	// dam_break();
	return EXIT_SUCCESS;
}
void dam_break(){
	// Parameters of the problem
	double R = 0.1;
	double l = 0.057; // particle distribution on [-l,l] x [-l,l]
	double L = 1; // size of the domain: [-L,L] x [-L,L]
	double H = 1;
	double dt = 1.0e-4; // physical time step
	double T = 0.2; // duration of simulation
	bool gravity = 1; // 1 if we consider the gravity

	// Physical parameters
	double rho_0 = 1000.0; // initial (physical) density of water at 20°C (in kg/m^3)
	double mu = 1.0016e-3; // dynamic viscosity of water at 20°C (in N.s/m^2)
	double gamma = 7.0; // typical value for liquid (dimensionless)
	double c_0 = 1.0;//1481; // sound speed in water at 20°C (in m/s)
	double sigma = 72.86e-3; // surface tension of water-air interface at 20°C (in N/m)


	// SPH parameters
	int n_per_dim = 51; // number of particles per dimension
	double kh = sqrt(21) * 2 * l / n_per_dim; // kernel width to ensure 21 particles in the neighborhood
	int n_iter = (int)(T / dt); // number of iterations to perform
	Kernel kernel = Cubic; // kernel choice
	double interface_threshold = 1.2;//1.5; // If ||n_i|| > threshold => particle i belongs to interface (first detection approach)
	Verlet *verlet = NULL; // don't use Verlet (for now)
	double XSPH_epsilon = 0.5;
	Free_surface_detection surface_detection = DIVERGENCE;
	double CR = 1.0;
	double CF = 0.0;

	printf("n_iter = %d\n", n_iter);

	// Animation parameter
	double T_anim = 10; // duration of animation
	double dt_anim = T_anim / n_iter; // time step of animation

	// Initialize particles on a square
	int n_p = squared(n_per_dim); // total number of particles
	double h = 2 * l / (n_per_dim - 1); // step between neighboring particles
	double m = rho_0 * h*h;
	Particle** particles = (Particle**)malloc(n_p * sizeof(Particle*));
	Particle_derivatives** particles_derivatives = malloc(n_p * sizeof(Particle_derivatives*));
	Residual** residuals = malloc(n_p * sizeof(Residual*));
	for (int i = 0; i < n_per_dim; i++) {
		for (int j = 0; j < n_per_dim; j++) {
			int index = i * n_per_dim + j;
			xy *pos = xy_new(-l + i * h, -l + j * h);
			xy *v = xy_new(0.0, 0.0); // initial velocity = 0
			// Insert initial condition on velocity here
			particles[index] = Particle_new(index, m, pos, v, rho_0, mu, c_0, gamma, sigma);
			particles_derivatives[index] = Particle_derivatives_new(index);
			residuals[index] = Residual_new();
		}
	}
	// Setup grid
	Grid *grid = Grid_new(-L, 3*L, -H, 2*H, kh);
	// Setup BOUNDARY
	double lb = 0.420;
	double hb = 0.440;
	double Rp = 0.001; //particle radius
	Boundary* boundary = Boundary_new(-l-Rp,lb-l+Rp,-l-Rp,hb-l+Rp,CR,CF);

	// Setup setup
	Setup *setup = Setup_new(n_iter, dt, kh, verlet, kernel, surface_detection, interface_threshold, XSPH_epsilon, gravity);
	// Setup animation
	Animation *animation = Animation_new(n_p, dt_anim, grid, 1);
	// Simulation
	simulate_boundary(grid, particles, particles_derivatives, residuals, n_p, update_positions_seminar_5, setup, animation, boundary);
	// Free memory
	Boundary_free(boundary);
	free_particles(particles, n_p);
	free_particles_derivatives(particles_derivatives, n_p);
	free_Residuals(residuals, n_p);
	Grid_free(grid);
	Setup_free(setup);
	Animation_free(animation);
}
// Evolution of a 2D circle with non-zero initial velocities (no surface tension force)
// Test case from "Simulating Free Surface Flows with SPH", Monaghan (1994)
void script_circle_to_ellipse() {

	// Parameters of the problem
	double l = 1.0; // radius of the circle
	double L = 2.0*l; // size of the domain: [-L,L] x [-L,L]
	double dt = 1.0e-5; // physical time step
	double T = 0.0076; // duration of simulation
	bool gravity = 0;

	// Physical parameters
	double rho_0 = 1000.0; // initial (physical) density of water at 20°C (in kg/m^3)
	double mu = 1.0016e-3; // dynamic viscosity of water at 20°C (in N.s/m^2)
	double gamma = 7.0; // typical value for liquid (dimensionless)
	double c_0 = 1400.0;//1481; // sound speed in water at 20°C (in m/s)
	double sigma = 0.0; // surface tension of water-air interface at 20°C (in N/m)

	// SPH parameters
	Verlet *verlet = NULL; // don't use Verlet (for now)
	Kernel kernel = Cubic; // kernel choice
	double interface_threshold = 1.5;//1000.0; // If ||n_i|| > threshold => particle i belongs to interface (first detection approach)
	double XSPH_epsilon = 0.5;
	Free_surface_detection surface_detection = DIVERGENCE;
	int N_c = 30; // number of circonferences on which points are placed
	int N_p = 6; // number of points on the first circonference (doubled for every circonference)
	int N_tot = 1; // total number of points
	for (int i = 1; i < N_c; i++) {
		N_tot += i * N_p;
	}
	printf("N_tot = %d \n", N_tot);
	int n_iter = (int)(T / dt); // number of iterations to perform
	double kh = 0.2*l;// is ideal to reach t = 0.0076; // kernel width


	// Animation parameter
	double T_anim = 0.1; // duration of animation
	double dt_anim = T_anim / n_iter; // time step of animation

	// Initialize particles in a circle
	double m = rho_0 * M_PI * l * l / N_tot; // mass of each particle

	Particle** particles = (Particle**)malloc(N_tot * sizeof(Particle*));
	Particle_derivatives** particles_derivatives = malloc(N_tot * sizeof(Particle_derivatives*));
	Residual** residuals = malloc(N_tot * sizeof(Residual*));

	// parameters defining the circle
	double b, delta_s, k, theta;
	delta_s = l / ((double)N_c - 1.0);
	theta = (2*M_PI) / ((double)N_p);

	int index = 0;
	for (int i = 0; i < N_c; i++) {
		b = i;
		if (b == 0) {
			xy *pos = xy_new(0.0, 0.0);
			xy *v = xy_new(0.0, 0.0);
			particles[index] = Particle_new(index, m, pos, v, rho_0, mu, c_0, gamma, sigma);
			particles_derivatives[index] = Particle_derivatives_new(index);
			residuals[index] = Residual_new();
			index++;
		}
		else {
			for (int j = 0; j < i*N_p; j++) {
				k = (double)j / b;
				xy *pos = xy_new(b*delta_s*cos(k*theta), b*delta_s*sin(k*theta));
				xy *v = xy_new(-100.0*pos->x, 100.0*pos->y);
				particles[index] = Particle_new(index, m, pos, v, rho_0, mu, c_0, gamma, sigma);
				particles_derivatives[index] = Particle_derivatives_new(index);
				residuals[index] = Residual_new();
				index++;
			}
		}
	}


	// Setup grid
	Grid *grid = Grid_new(-L, L, -L, L, kh);
	// Setup animation
	Animation *animation = Animation_new(N_tot, dt_anim, grid, 1);
	// Setup setup
	Setup *setup = Setup_new(n_iter, dt, kh, verlet, kernel, surface_detection, interface_threshold, XSPH_epsilon,gravity);

	// Start simulation
	simulate(grid, particles, particles_derivatives, residuals, N_tot, update_positions_ellipse, setup, animation);

	// Free stuff
	free_particles(particles, N_tot);
	free_particles_derivatives(particles_derivatives, N_tot);
	free_Residuals(residuals, N_tot);
	Grid_free(grid);
	Setup_free(setup);
	Animation_free(animation);
}
