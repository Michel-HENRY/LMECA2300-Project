#include "vector.h"
#include "kernel.h"
#include "particle.h"
#include "animation.h"
#include "time_integration.h"
#include "boundary.h"

#include <time.h> // To generate random number
#include <stdlib.h>

// Validation function
int dam_break();
int boundary_validation();
int SPH_operator_validation();

// Useful function to create the usual sructures
Particle** fluidProblem(Parameters* param, int n_p_dim_x, int n_p_dim_y, double g, double rho, double P, bool isUniform);
Edges* get_box(double L, double H, int n_e , double CF, double CR, double domain[4]);

int main(){
  // dam_break();
  boundary_validation();
  // SPH_operator_validation();
  // kernel_validation();
}
Particle** fluidProblem(Parameters* param, int n_p_dim_x, int n_p_dim_y, double g, double rho, double P, bool isUniform){
  int n_p = n_p_dim_x*n_p_dim_y;
  Particle** particles = (Particle**) malloc(n_p*sizeof(Particle*));
  for(int i = 0; i < n_p_dim_x; i++){
      for(int j = 0; j < n_p_dim_y; j++){
        int index = i*n_p_dim_y + j;

        Vector* x = Vector_new(2);
        Vector* u = Vector_new(2);
        Vector* f = Vector_new(2);

        f->X[1] = -g;

        double Rp = param->Rp;
        if (isUniform){
          x->X[0] = Rp*(2*i + 1);      x->X[1] = Rp*(2*j + 1);
        }
        else{
          double l = Rp*n_p_dim_x;
          x->X[0] = (float)rand()/(float)(RAND_MAX/l);
          x->X[1] = (float)rand()/(float)(RAND_MAX/l);
          u->X[0] = 10*(float)rand()/(float)(RAND_MAX/l);
          u->X[1 ]= 10*(float)rand()/(float)(RAND_MAX/l);
          //
          //     u->X[0] = ux;
          //     u->X[1] = uy;
        }
        Fields* fields = Fields_new(x,u,f,P,rho);
        particles[index] = Particle_new(param, fields);
      }
    }
    return particles;
}
Edges* get_box(double L, double H, int n_e , double CF, double CR, double domain[4]){
  // Create the vertices
  Vector** vertices = (Vector**) malloc(n_e*sizeof(vertices));
  for(int i = 0; i < n_e; i++){
    vertices[i] = Vector_new(2);
  }
  vertices[0]->X[0] = 0;                vertices[2]->X[0] = L;
  vertices[0]->X[1] = 0;                vertices[2]->X[1] = H;
  vertices[1]->X[0] = L;                vertices[3]->X[0] = 0;
  vertices[1]->X[1] = 0;                vertices[3]->X[1] = H;

  // Create the edges structure
  Vector** edge = (Vector**) malloc(n_e*2*sizeof(vertices));
  for(int i = 0;i < n_e; i++){
    edge[2*i] = copy(vertices[i]);
    edge[2*i+1] = copy(vertices[(i+1)%n_e]);
  }
  for(int i = 0; i < n_e ; i++){
    Vector_free(vertices[i]);
  }
  free(vertices);

  Edges* edges = Edges_new(n_e, edge, CR,CF);
  domain[0] = 0;
  domain[1] = L;
  domain[2] = 0;
  domain[3] = H;
  return edges;
}

int boundary_validation(){
  double lx = 1;                          // Longueur du domaine de particule
  double ly = 1;                          // Hauteur du domaine de particle
  int n_p_dim = 100;                       // Nombre moyen de particule par dimension

  // Parameters
  double rho_0 = 1e3;                     // Densité initiale
  double dynamic_viscosity = 0;           // Viscosité dynamique
  double g = 0.00;                        // Gravité
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
  double P0 = 0;                          // Pression atmosphérique

  // ------------------------------------------------------------------
  // ------------------------ SET Particles ---------------------------
  // ------------------------------------------------------------------
  Parameters* param = Parameters_new(mass, dynamic_viscosity, kh, Rp, tension, treshold, P0);
  Particle** particles = fluidProblem(param, n_p_dim_x, n_p_dim_y, g, rho_0, P0,false);

  // ------------------------------------------------------------------
  // ------------------------ SET Edges -------------------------------
  // ------------------------------------------------------------------

  double L = 1;
  double H = 1;
  int n_e = 4;
  double CF = 0.0;
  double CR = 1.0;

  double domain[4];
  Edges* edges = get_box(L,H,n_e,CF,CR,domain);
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
  double tEnd = 10;
  double dt = 0.001;
  int iter_max = (int) (tEnd-t)/dt;
  int output = 1;
  printf("iter max = %d\n",iter_max);
  // // Temporal loop
  Kernel kernel = Cubic;
  int i = 0;
  while (t < tEnd){
    printf("-----------\t t/tEnd : %.3f/%f\t-----------\n", t,tEnd);
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
  Animation_free(animation);
  printf("END FREE ANIMATION\n");
  return EXIT_SUCCESS;
}

int SPH_operator_validation(){
  double lx = 1;                          // Longueur du domaine de particule
  double ly = 1;                          // Hauteur du domaine de particle
  int n_p_dim = 25;                       // Nombre moyen de particule par dimension

  // Parameters
  double rho_0 = 1;                       // Densité initiale
  double dynamic_viscosity = 1;           // Viscosité dynamique
  double g = 0.00;                        // Gravité
  int n_p_dim_x = n_p_dim*lx;             // Nombre de particule par dimension
  int n_p_dim_y = n_p_dim*ly;
  int n_p = n_p_dim_x*n_p_dim_y;          // Nombre de particule total
  double h = lx/n_p_dim_x;                // step between neighboring particles
  double kh = sqrt(21)*lx/n_p_dim_x;      // Rayon du compact pour l'approximation
  double mass = rho_0*lx*ly/n_p;          // Masse d'une particule, constant
  double Rp = h/2;                        // Rayon d'une particule
  double eta = 0.0;                       // XSPH parameter from 0 to 1
  double treshold = 20;                   // Critère pour la surface libre
  double tension = 72*1e-3;               // Tension de surface de l'eau
  double P0 = 0;                          // Pression atmosphérique

  // ------------------------------------------------------------------
  // ------------------------ SET Particles ---------------------------
  // ------------------------------------------------------------------
  Parameters* param = Parameters_new(mass, dynamic_viscosity, kh, Rp, tension, treshold, P0);
  Particle** particles = fluidProblem(param, n_p_dim_x, n_p_dim_y, g, rho_0, P0,true);

  // Impose linear gradient pressure :
  for(int i = 0; i < n_p ; i ++){
    double x,y;
    x = particles[i]->fields->x->X[0];
    y = particles[i]->fields->x->X[1];
    // particles[i]->fields->P = 2*x + y;
    // particles[i]->fields->u->X[0] = x;
    // particles[i]->fields->u->X[1] = x*x;
  }

  // ------------------------------------------------------------------
  // ------------------------ SET Edges -------------------------------
  // ------------------------------------------------------------------

  double L = 1;
  double H = 1;
  int n_e = 4;
  double CF = 0.0;
  double CR = 1.0;

  double domain[4];
  Edges* edges = get_box(L,H,n_e,CF,CR,domain);
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

  // ------------------------------------------------------------------
  // -------------------------- Computation ---------------------------
  // ------------------------------------------------------------------

  // // Temporal loop
  Kernel kernel = Cubic;

  update_cells(grid, particles, n_p);
  update_neighbors(grid, particles, n_p, 0);
  double meanx = 0;
  double meany = 0;
  double max = 0;
  double min = 0;
  for(int i = 0; i < n_p ; i++){
    // Vector* gradP = grad_P(particles[i], kernel);         //Pas précis
    // Vector* gradP = CSPM_pressure(particles[i], kernel);  //Ordre 1 sur les noeuds intérieurs
    // Vector_print(gradP[i]);
    // Vector_free(gradP);


    // Vector* lapl_s = lapl_u_shao(particles[i],kernel);      //Pas tres precis : Ordre 0
    // Vector* lapl = lapl_u(particles[i], kernel);            //Pareil que shao : Ordre 0
    // printf("----------------\n");
    // Vector_print(lapl_s);
    // Vector_print(lapl);
    // Vector_free(lapl);
    // Vector_free(lapl_s);

    // double divergence = div_u(particles[i],kernel);
    // printf("div = %f\n",divergence);
  }

  show(particles,animation, (int) 1, true, false);



  // ------------------------------------------------------------------
  // ------------------------ FREE Memory -----------------------------
  // ------------------------------------------------------------------
  Particles_free(particles, n_p);
  printf("END FREE PARTICLES\n");
  Edges_free(edges);
  printf("END FREE EDGES\n");
  Grid_free(grid);
  printf("END FREE GRID\n");
  Animation_free(animation);
  printf("END FREE ANIMATION\n");
  return EXIT_SUCCESS;
}

int dam_break(){
  double lx = 1;                          // Longueur du domaine de particule
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
  double P0 = 0;                          // Pression atmosphérique

  // ------------------------------------------------------------------
  // ------------------------ SET Particles ---------------------------
  // ------------------------------------------------------------------
  // Particle** particles = (Particle**) malloc(n_p*sizeof(Particle*));
  Parameters* param = Parameters_new(mass, dynamic_viscosity, kh, Rp, tension, treshold,P0);
  Particle** particles = fluidProblem(param, n_p_dim_x, n_p_dim_y, g, rho_0, P0,true);

  // ------------------------------------------------------------------
  // ------------------------ SET Edges -------------------------------
  // ------------------------------------------------------------------

  double L = 4;
  double H = 4;
  int n_e = 4;
  double CF = 0.0;
  double CR = 1.0;

  double domain[4];
  Edges* edges = get_box(L,H,n_e,CF,CR,domain);
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
  Animation_free(animation);
  printf("END FREE ANIMATION\n");
  return EXIT_SUCCESS;
}
