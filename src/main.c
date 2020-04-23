#include "vector.h"
#include "kernel.h"
#include "particle.h"
#include "animation.h"
#include "time_integration.h"
#include "boundary.h"

#include <time.h> // To generate random number
#include <stdlib.h>

int dam_break();
int one_particle();

Particles** fluidProblem();

int main(){
  // dam_break();
  one_particle();
}
int one_particle(){
  double lx = 1;                         // Longueur du domaine de particule
  double ly = 1;                          // Hauteur du domaine de particle
  int n_p_dim = 10;

  // Parameters
  double rho_0 = 1e3;                     // Densité initiale
  double dynamic_viscosity = 0;        // Viscosité dynamique
  double g = 0.00;                        // Gravité
  int n_p_dim_x = n_p_dim*lx;             // Nombre de particule par dimension
  int n_p_dim_y = n_p_dim*ly;
  int n_p = n_p_dim_x*n_p_dim_y;          // Nombre de particule total
  double h = lx/25;                // step between neighboring particles
  double kh = sqrt(21)*lx/n_p_dim_x;      // Rayon du compact pour l'approximation
  double mass = rho_0 * h*h;              // Masse d'une particule, constant
  double Rp = h/2;                        // Rayon d'une particule
  double eta = 0.0;                       // XSPH parameter from 0 to 1
  double treshold = 20;                   // Critère pour la surface libre
  double tension = 72*1e-3;               // Tension de surface de l'eau
  double P0 = 1e5;                        // Pression atmosphérique

  // ------------------------------------------------------------------
  // ------------------------ SET Particles ---------------------------
  // ------------------------------------------------------------------
  Particle** particles = (Particle**) malloc(n_p*sizeof(Particle*));
  for(int i = 0; i < n_p_dim_x; i++){
    for(int j = 0; j < n_p_dim_y; j++){
      int index = i*n_p_dim_y + j;
      Parameters* param = Parameters_new(rho_0, mass, dynamic_viscosity, kh, Rp, tension, treshold,P0);
      Vector* x = Vector_new(2);
      Vector* u = Vector_new(2);
      Vector* f = Vector_new(2);

      float posx = (float)rand()/(float)(RAND_MAX/lx);
      float posy = (float)rand()/(float)(RAND_MAX/ly);

      double pos[2] = {posx,posy};
      double P = 0;

      float ux = 10*(float)rand()/(float)(RAND_MAX/lx);
      float uy = 10*(float)rand()/(float)(RAND_MAX/ly);

      u->X[0] = ux;
      u->X[1] = uy;

      Vector_initialise(x,pos);
      Fields* fields = Fields_new(x,u,f,P);
      particles[index] = Particle_new(param, fields);
    }
  }
  // ------------------------------------------------------------------
  // ------------------------ SET Edges -------------------------------
  // ------------------------------------------------------------------

  double L = 1;
  double H = 1;
  int n_e = 4;
  double CF = 0.0;
  double CR = 1.0;

  Vector** vertices = (Vector**) malloc(n_e*sizeof(vertices));
  for(int i = 0; i < n_e; i++){
    vertices[i] = Vector_new(2);
  }
  vertices[0]->X[0] = 0;                vertices[2]->X[0] = L;
  vertices[0]->X[1] = 0;                vertices[2]->X[1] = H;
  vertices[1]->X[0] = L;                vertices[3]->X[0] = 0;
  vertices[1]->X[1] = 0;                vertices[3]->X[1] = H;


  Vector** edge = (Vector**) malloc(n_e*2*sizeof(vertices));
  for(int i = 0;i < n_e; i++){
    edge[2*i] = vertices[i];
    edge[2*i+1] = vertices[(i+1)%n_e];
  }
  Edges* edges = Edges_new(n_e, edge, CR,CF);
  double domain[4] = {0,L,0,H};
  // ------------------------------------------------------------------
  // ------------------------ SET Grid --------------------------------
  // ------------------------------------------------------------------
  double extra = 0.0;
  extra = 0.5;
  Grid* grid = Grid_new(0-extra, L+extra, 0-extra, H+extra, kh);

  // ------------------------------------------------------------------
  // ------------------------ SET Animation ---------------------------
  // ------------------------------------------------------------------
  double timeout = 0.001;                 // Durée d'une frame
  Animation* animation = Animation_new(n_p, timeout, grid, Rp, domain);

  // ------------------------------------------------------------------
  // ------------------------ Start integration -----------------------
  // ------------------------------------------------------------------
  double t = 0;
  double tEnd = 1;
  double dt = 0.001;
  int iter_max = (int) (tEnd-t)/dt;
  int output = 1;
  printf("iter max = %d\n",iter_max);
  // // Temporal loop
  Kernel kernel = Cubic;
  int i = 0;
  while (t < tEnd){
    printf("-----------\t t/tEnd : %.3f/%.1f\t-----------\n", t,tEnd);
    if (i%output == 0)
      show(particles, animation, i, false, false);
    // update_cells(grid, particles, n_p);
    // update_neighbors(grid, particles, n_p, i);
    // update_pressureDam(particles, n_p, rho_0, g, H);
    // time_integration(particles, n_p, kernel, dt, edges);
    time_integration_CSPM(particles, n_p, kernel, dt, edges,eta);
    if (i%output == 0)
      show(particles, animation, i, false, false);
    printf("Time integration completed\n");

    i++;
    t += dt;
  }
  show(particles,animation, iter_max, true, false);



  // ------------------------------------------------------------------
  // ------------------------ FREE Memory -----------------------------
  // ------------------------------------------------------------------
  Particles_free(particles, n_p);
  printf("END FREE PARTICLES\n");
  Edges_free(edges);
  printf("END FREE EDGES\n");
  Grid_free(grid);
  printf("END FREE GRID\n");
  // Animation_free(animation);
  // printf("END FREE ANIMATION\n");
  return EXIT_SUCCESS;
}

int dam_break(){
  double lx = 1;                         // Longueur du domaine de particule
  double ly = 2;                          // Hauteur du domaine de particle
  int n_p_dim = 30;

  // Parameters
  double rho_0 = 1e3;                     // Densité initiale
  double dynamic_viscosity = 1e-6;        // Viscosité dynamique
  double g = 9.81;                        // Gravité
  int n_p_dim_x = n_p_dim*lx;             // Nombre de particule par dimension
  int n_p_dim_y = n_p_dim*ly;
  int n_p = n_p_dim_x*n_p_dim_y;          // Nombre de particule total
  double h = lx/n_p_dim_x;                // step between neighboring particles
  double kh = sqrt(21)*lx/n_p_dim_x;      // Rayon du compact pour l'approximation
  double mass = rho_0 * h*h;              // Masse d'une particule, constant
  double Rp = h/2;                        // Rayon d'une particule
  double eta = 0.0;                       // XSPH parameter from 0 to 1
  double treshold = 20;                   // Critère pour la surface libre
  double tension = 72*1e-3;               // Tension de surface de l'eau
  double P0 = 1e5;                        // Pression atmosphérique

  // ------------------------------------------------------------------
  // ------------------------ SET Particles ---------------------------
  // ------------------------------------------------------------------
  Particle** particles = (Particle**) malloc(n_p*sizeof(Particle*));
  for(int i = 0; i < n_p_dim_x; i++){
    for(int j = 0; j < n_p_dim_y; j++){
      int index = i*n_p_dim_y + j;
      Parameters* param = Parameters_new(rho_0, mass, dynamic_viscosity, kh, Rp, tension, treshold,P0);
      Vector* x = Vector_new(2);
      Vector* u = Vector_new(2);
      Vector* f = Vector_new(2);

      f->X[1] = -g;

      double pos[2] = {Rp + i*h ,Rp + j*h};
      double P = 0;
      // if(i == 0){
        // u->X[0] = -2*(1 - pos[1]*pos[1]);
      // }

      Vector_initialise(x,pos);
      Fields* fields = Fields_new(x,u,f,P);
      particles[index] = Particle_new(param, fields);
    }
  }
  // ------------------------------------------------------------------
  // ------------------------ SET Edges -------------------------------
  // ------------------------------------------------------------------

  double L = 4;
  double H = 4;
  int n_e = 4;
  double CF = 0.0;
  double CR = 1.0;

  Vector** vertices = (Vector**) malloc(n_e*sizeof(vertices));
  for(int i = 0; i < n_e; i++){
    vertices[i] = Vector_new(2);
  }
  vertices[0]->X[0] = 0;                vertices[2]->X[0] = L;
  vertices[0]->X[1] = 0;                vertices[2]->X[1] = H;
  vertices[1]->X[0] = L;                vertices[3]->X[0] = 0;
  vertices[1]->X[1] = 0;                vertices[3]->X[1] = H;


  Vector** edge = (Vector**) malloc(n_e*2*sizeof(vertices));
  for(int i = 0;i < n_e; i++){
    edge[2*i] = vertices[i];
    edge[2*i+1] = vertices[(i+1)%n_e];
  }
  Edges* edges = Edges_new(n_e, edge, CR,CF);
  double domain[4] = {0,L,0,H};
  // ------------------------------------------------------------------
  // ------------------------ SET Grid --------------------------------
  // ------------------------------------------------------------------
  double extra = 0.0;
  extra = 0.5;
  Grid* grid = Grid_new(0-extra, L+extra, 0-extra, H+extra, kh);

  // ------------------------------------------------------------------
  // ------------------------ SET Animation ---------------------------
  // ------------------------------------------------------------------
  double timeout = 0.001;                 // Durée d'une frame
  Animation* animation = Animation_new(n_p, timeout, grid, Rp, domain);

  // ------------------------------------------------------------------
  // ------------------------ Start integration -----------------------
  // ------------------------------------------------------------------
  double t = 0;
  double tEnd = 1;
  double dt = 0.0001;
  int iter_max = (int) (tEnd-t)/dt;
  int output = 1;
  printf("iter max = %d\n",iter_max);
  // // Temporal loop
  Kernel kernel = Cubic;
  int i = 0;
  while (t < tEnd){
    printf("-----------\t t/tEnd : %.3f/%.1f\t-----------\n", t,tEnd);
    if (i%output == 0)
      show(particles, animation, i, false, false);
    update_cells(grid, particles, n_p);
    update_neighbors(grid, particles, n_p, i);
    update_pressureDam(particles, n_p, rho_0, g, H);
    // time_integration(particles, n_p, kernel, dt, edges);
    time_integration_CSPM(particles, n_p, kernel, dt, edges,eta);
    if (i%output == 0)
      show(particles, animation, i, false, false);
    printf("Time integration completed\n");

    i++;
    t += dt;
  }
  show(particles,animation, iter_max, true, false);



  // ------------------------------------------------------------------
  // ------------------------ FREE Memory -----------------------------
  // ------------------------------------------------------------------
  Particles_free(particles, n_p);
  printf("END FREE PARTICLES\n");
  Edges_free(edges);
  printf("END FREE EDGES\n");
  Grid_free(grid);
  printf("END FREE GRID\n");
  // Animation_free(animation);
  // printf("END FREE ANIMATION\n");
  return EXIT_SUCCESS;
}
