#include "vector.h"
#include "kernel.h"
#include "particle.h"
#include "animation.h"
#include "time_integration.h"
#include "boundary.h"

void print_if_neighbors(Particle* pi){
  ListNode* node = pi->neighbors->head;
  int i = 0;
  while(node != NULL){
    printf("%d\t", i);
    node = node->next;
    i++;
  }
}
int main(){
  //particle domain
  double l = 1.0;
  double ly = l;
  // Grid definition

  double timeout = 0.001; //Pour accelerer la simu

  // Parameters
  double rho_0 = 1e3;
  double dynamic_viscosity = 1e-6;
  double g = 9.81;
  int n_p_dim = 75;
  int n_p = n_p_dim*n_p_dim;
  double h = l/n_p_dim; // step between neighboring particles
  double kh = sqrt(21)*l/n_p_dim;
  double mass = rho_0 * h*h;
  double Rp = h/2;
  double eta = 0.5;

  // ------------------------------------------------------------------
  // ------------------------ SET Particles ---------------------------
  // ------------------------------------------------------------------
  Particle** particles = (Particle**) malloc(n_p*sizeof(Particle*));
  for(int i = 0; i < n_p_dim; i++){
    for(int j = 0; j < n_p_dim; j++){
      int index = i*n_p_dim + j;
      Parameters* param = Parameters_new(rho_0, mass, dynamic_viscosity, kh, Rp);
      Vector* x = Vector_new(2);
      Vector* u = Vector_new(2);
      Vector* f = Vector_new(2);

      // f->X[1] = -g;

      double pos[2] = {Rp + i*h,Rp + j*h};
      double P = 0;
      // if(i == 0){
      //   u->X[0] = 2*(1 - pos[1]*pos[1]);
      // }

      Vector_initialise(x,pos);
      Fields* fields = Fields_new(x,u,f,P);
      particles[index] = Particle_new(param, fields);
    }
  }
  // ------------------------------------------------------------------
  // ------------------------ SET Edges -------------------------------
  // ------------------------------------------------------------------

  double L = 1.5;
  double H = 1;
  int n_e = 4;
  double CF = 0.5;
  double CR = 0.5;

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
  Animation* animation = Animation_new(n_p, timeout, grid, Rp, domain);

  // ------------------------------------------------------------------
  // ------------------------ Start integration -----------------------
  // ------------------------------------------------------------------
  double t = 0;
  double tEnd = 100;
  double dt = 0.1;
  int iter_max = (int) (tEnd-t)/dt;
  int output = 1;
  printf("iter max = %d\n",iter_max);
  // // Temporal loop
  Kernel kernel = Cubic;
  int i = 0;
  while (t < tEnd){
    printf("-----------\t t/tEnd : %.3f/%.1f\t-----------\n", t,tEnd);
    if (i%output == 0)
      show(particles, animation, i, false, true);
    update_cells(grid, particles, n_p);
    update_neighbors(grid, particles, n_p, i);
    update_pressureMod(particles, n_p, rho_0,g,ly);
    printf("P = %f\n",particles[0]->fields->P);
    // time_integration(particles, n_p, kernel, dt, edges);
    time_integration_CSPM(particles, n_p, kernel, dt, edges,eta);
    printf("Time integration completed\n");

    i++;
    t += dt;
  }
  show(particles,animation, iter_max, false, false);



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
