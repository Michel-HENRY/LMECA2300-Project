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

  double timeout = 1;

  // Parameters
  double rho = 1e3;
  double dynamic_viscosity = 1e-3;
  double g = 0;
  int n_p_dim = 10;
  int n_p = n_p_dim*n_p_dim;
  double h = l/n_p_dim; // step between neighboring particles
  double kh = sqrt(21) * l / n_p_dim;
  double mass = rho * h*h;
  double Rp = h/2;

  // ------------------------------------------------------------------
  // ------------------------ SET Particles ---------------------------
  // ------------------------------------------------------------------
  Particle** particles = (Particle**) malloc(n_p*sizeof(Particle*));

  for(int i = 0; i < n_p_dim; i++){
    for(int j = 0; j < n_p_dim; j++){
      int index = i*n_p_dim + j;
      Parameters* param = Parameters_new(rho, mass, dynamic_viscosity, kh, Rp);
      Vector* x = Vector_new(2);
      Vector* u = Vector_new(2);
      Vector* f = Vector_new(2);
      u->X[0] = 1.3/2;
      // u->X[1] = 0.7/2;
      double  P = rho*g*(l - j*h);
      // double P = 0;
      double pos[2] = {Rp + i*h,Rp + j*h};

      Vector_initialise(x,pos);
      Fields* fields = Fields_new(x,u,f,P);
      particles[index] = Particle_new(param, fields);
    }
  }
  // ------------------------------------------------------------------
  // ------------------------ SET Edges -------------------------------
  // ------------------------------------------------------------------

  double L = 1.1;
  double H = 1.1;
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
  extra = ly/2;
  Grid* grid = Grid_new(0-extra, L+extra, -extra, H+extra, kh);

  // ------------------------------------------------------------------
  // ------------------------ SET Animation ---------------------------
  // ------------------------------------------------------------------
  Animation* animation = Animation_new(n_p, timeout, grid, Rp, domain);

  // ------------------------------------------------------------------
  // ------------------------ Start integration -----------------------
  // ------------------------------------------------------------------
  double t = 0;
  double tEnd = 50;
  double dt = 0.05;
  int iter_max = (int) (tEnd-t)/dt;
  int output = 1;
  printf("iter max = %d\n",iter_max);
  // // Temporal loop
  Kernel kernel = Lucy;
  int i = 0;
  while (t < tEnd){
    printf("-----------\t t/tEnd : %.3f/%.1f\t-----------\n", t,tEnd);
    if (i%output == 0)
      show(particles, animation, i, false, true);
    update_cells(grid, particles, n_p);
    update_neighbors(grid, particles, n_p, i);
    // update_pressure(particles, n_p, rho, g, l); // Pressure dyn + P hydro
    time_integration(particles, n_p, kernel, dt, edges);
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
