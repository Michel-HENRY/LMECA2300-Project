#include "vector.h"
#include "kernel.h"
#include "particle.h"
#include "animation.h"


int main(){
  // List_validation();
  // Vector_validation();
  // Particle_validation();

  //particle domain
  double l = 1.0;
  // Grid definition
  // double L = 1.2;
  // double H = 1.2;

  double timeout = 2.5;

  // Parameters
  double rho = 1e3;
  double dynamic_viscosity = 1e-3;
  int n_p_dim = 100;
  int n_p = n_p_dim*n_p_dim;
  double h = l/(n_p_dim - 1); // step between neighboring particles
  double mass = rho * h*h;
  double R_p = l/(2*n_p_dim);

  // Create Particles
  Particle** particles = (Particle**) malloc(n_p*sizeof(Particle*));

  for(int i = 0; i < n_p_dim; i++){
    for(int j = 0; j < n_p_dim; j++){
      int index = i*n_p_dim + j;
      Parameters* param = Parameters_new(rho, mass, dynamic_viscosity, h);
      Vector* x = Vector_new(2);
      Vector* u = Vector_new(2);
      Vector* f = Vector_new(2);
      double P = l-j*h;
      double pos[2] = {i*h,j*h};

      Vector_initialise(x,pos);
      Fields* fields = Fields_new(x,u,f,P);
      particles[index] = Particle_new(param, fields);
    }
  }

  // // Create GRID
  double extra = 0.1;
  Grid* grid = Grid_new(-extra, l+extra, -extra, l+extra, h);
  double domain[4] = {0,l,0,l};
  Animation* animation = Animation_new(n_p, timeout, grid, R_p, domain);
  show(particles, animation, 0, false, false);
  // // Link grid and particles
  update_cells(grid, particles, n_p);
  update_neighbors(grid, particles, n_p, 0);



  // Free memory
  Particles_free(particles, n_p);
  printf("END FREE PARTICLES\n");
  Grid_free(grid);
  printf("END FREE GRID\n");
  Animation_free(animation);
  printf("END FREE ANIMATION\n");
  return EXIT_SUCCESS;
}
