#include "vector.h"
#include "kernel.h"
#include "particle.h"
#include "animation.h"
#include "time_integration.h"
#include "boundary.h"

#include <time.h> // To generate random number
#include <stdlib.h>

#define M_PI 3.14159265358979323846

// Validation function
static int dam_break();
static int boundary_validation();
static int SPH_operator_validation();
static int hydrostatic_eq();
static int waves();
static int vortices();

// Useful function to create the usual sructures
Particle** fluidProblem(Parameters* param, int n_p_dim_x, int n_p_dim_y, double g, double rho, double P, bool isUniform);
Edges* get_box(double L, double H, int n_e , double CF, double CR, double domain[4]);

int main(){
  dam_break();
  // boundary_validation();
  // SPH_operator_validation();
  // hydrostatic_eq();
  // waves();
}
Particle** fluidProblem(Parameters* param, int n_p_dim_x, int n_p_dim_y, double g, double rho, double P, bool isUniform){
  int n_p = n_p_dim_x*n_p_dim_y;
  Particle** particles = (Particle**) malloc(n_p*sizeof(Particle*));
  for(int j = 0; j< n_p_dim_y; j++){
      for(int i = 0; i < n_p_dim_x; i++){
        int index = j*n_p_dim_x + i;

        Vector* x = Vector_new(2);
        Vector* u = Vector_new(2);
        Vector* f = Vector_new(2);

        // f->X[1] = -g*rho;

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

static int boundary_validation(){
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
  Parameters* param = Parameters_new(mass, dynamic_viscosity, kh, Rp, tension, treshold, P0,g);
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

static int SPH_operator_validation(){
  double lx = 1;                          // Longueur du domaine de particule
  double ly = 1;                          // Hauteur du domaine de particle
  int n_p_dim = 51;                       // Nombre moyen de particule par dimension

  // Parameters
  double rho_0 = 1000;                       // Densité initiale
  double dynamic_viscosity = 1;           // Viscosité dynamique
  double g = 0.00;                        // Gravité
  int n_p_dim_x = n_p_dim;             // Nombre de particule par dimension
  int n_p_dim_y = n_p_dim*(ly/lx);
  int n_p = n_p_dim_x*n_p_dim_y;          // Nombre de particule total
  double h = lx/n_p_dim_x;                // step between neighboring particles
  double kh = 6*h;      // Rayon du compact pour l'approximation
  double mass = rho_0*lx*ly/n_p;          // Masse d'une particule, constant
  double Rp = h/2;                        // Rayon d'une particule
  double eta = 0.0;                       // XSPH parameter from 0 to 1
  double treshold = 20;                   // Critère pour la surface libre
  double tension = 72*1e-3;               // Tension de surface de l'eau
  double P0 = 0;                          // Pression atmosphérique

  // ------------------------------------------------------------------
  // ------------------------ SET Particles ---------------------------
  // ------------------------------------------------------------------
  Parameters* param = Parameters_new(mass, dynamic_viscosity, kh, Rp, tension, treshold, P0,g);
  Particle** particles = fluidProblem(param, n_p_dim_x, n_p_dim_y, g, rho_0, P0,true);

  // Impose linear gradient pressure :
  for(int i = 0; i<n_p ; i++){
    double x,y;
    x = particles[i]->fields->x->X[0];
    y = particles[i]->fields->x->X[1];
    particles[i]->fields->P = y;
    particles[i]->fields->u->X[0] = x;
    particles[i]->fields->u->X[1] = x;
  }

  // ------------------------------------------------------------------
  // ------------------------ SET Edges -------------------------------
  // ------------------------------------------------------------------

  double L = lx;
  double H = ly;
  int n_e = 4;
  double CF = 0.0;
  double CR = 1.0;

  double domain[4];
  Edges* edges = get_box(L,H,n_e,CF,CR,domain);
  // ------------------------------------------------------------------
  // ------------------------ SET Grid --------------------------------
  // ------------------------------------------------------------------
  double extra = 0.0;
  extra = lx/2;
  Grid* grid = Grid_new(0-extra, L+extra, 0-extra, H+extra, kh);

  // ------------------------------------------------------------------
  // -------------------------- Computation ---------------------------
  // ------------------------------------------------------------------

  // // Temporal loop
  Kernel kernel = Quartic;
  // Kernel kernel = Cubic;
  // Kernel kernel = Lucy;

  update_cells(grid, particles, n_p);
  update_neighbors(grid, particles, n_p, 0);
  CSPM_density(particles,n_p,kernel);
  Vector** dP = get_dP(particles, n_p, kernel);
  Vector** dP_CSPM = get_dP(particles, n_p, kernel);

  for(int i = 0; i < n_p ; i++){
    times_into(dP[i], 1/particles[i]->fields->rho);
    times_into(dP_CSPM[i], 1/particles[i]->fields->rho);
  }
  CSPM_pressure(particles,n_p,kernel, dP_CSPM);
  for(int i = 0; i < n_p; i++){
    printf("dPx[%d] = %f\t dPx_CSPM[%d] = %f\t\t",i,dP[i]->X[0], i, dP_CSPM[i]->X[0]);
    printf("dPy[%d] = %f\t dPy_CSPM[%d] = %f",i,dP[i]->X[1], i, dP_CSPM[i]->X[1]);
    // Vector* lapl_s = lapl_u_shao(particles[i],kernel);      //Pas tres precis : Ordre 0
    // Vector* lapl = lapl_u(particles[i], kernel);            //Pareil que shao : Ordre 0
    // printf("----------------\n");
    // Vector_print(lapl_s);
    // Vector_print(lapl);
    // Vector_free(lapl);
    // Vector_free(lapl_s);

    double divergence = div_u(particles[i],kernel);
    printf("\t\tdiv = %f\n",divergence);

    printf("\n");
  }

  printf("h = %f\t",h);
  printf("kh = %f\n", kh );

  // show(particles,animation, (int) 1, true, false);



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

static int dam_break(){
  double lx = 0.114;                          // Longueur du domaine de particule
  double ly = 0.228;                          // Hauteur du domaine de particle
  int n_p_dim = 36;

  // Parameters
  double rho_0 = 1e3;                     // Densité initiale
  double dynamic_viscosity = 1e-3;        // Viscosité dynamique
  double g = 9.81;                        // Gravité
  int n_p_dim_x = n_p_dim;                // Nombre de particule par dimension
  int n_p_dim_y = n_p_dim*(ly/lx);
  int n_p = n_p_dim_x*n_p_dim_y;          // Nombre de particule total
  double delta = lx/n_p_dim_x;                // step between neighboring particles
  double h = sqrt(21)*lx/n_p_dim_x;      // Rayon du compact pour l'approximation
  // double h = 4*delta;
  double mass = rho_0*delta*delta;                // Masse d'une particule, constant
  double Rp = delta/2;                        // Rayon d'une particule
  double eta = 0.00;                       // XSPH parameter from 0 to 1
  double treshold = 20;                   // Critère pour la surface libre
  double tension = 0.0;//7*1e-2;                     // Tension de surface de l'eau
  double P0 = 0;                          // Pression atmosphérique

  // ------------------------------------------------------------------
  // ------------------------ SET Particles ---------------------------
  // ------------------------------------------------------------------
  Parameters* param = Parameters_new(mass, dynamic_viscosity, h, Rp, tension, treshold,P0,g);
  Particle** particles = fluidProblem(param, n_p_dim_x, n_p_dim_y, g, rho_0, P0,true);
  set_artificialViscosity(param, 0.3, 0);

  //Set density
  double B = 0.85*1e5;
  double gamma = 7;
  set_density(particles, n_p_dim_x, n_p_dim_y, rho_0, B, gamma);

  // ------------------------------------------------------------------
  // ------------------------ SET Edges -------------------------------
  // ------------------------------------------------------------------

  double L = 0.440;
  double H = 0.440;
  int n_e = 4;
  double CF = 0.0;
  double CR = 1.0;

  double domain[4];
  Edges* edges = get_box(L,H,n_e,CF,CR,domain);
  // ------------------------------------------------------------------
  // ------------------------ SET Grid --------------------------------
  // ------------------------------------------------------------------
  double extra = 0;
  extra = 0.01;
  Grid* grid = Grid_new(0-extra, L+extra, 0-extra, H+extra, h);

  // ------------------------------------------------------------------
  // ------------------------ SET Animation ---------------------------
  // ------------------------------------------------------------------
  double timeout = 0.1;                 // Durée d'une frame
  Animation* animation = Animation_new(n_p, timeout, grid, Rp, domain);

  // ------------------------------------------------------------------
  // ------------------------ Start integration -----------------------
  // ------------------------------------------------------------------
  double t = 0;
  double tEnd = 10;
  double dt = 1e-5;
  int iter_max = (int) (tEnd-t)/dt;
  int output = 1;
  printf("iter max = %d\n",iter_max);
  // // Temporal loop
  Kernel kernel = Cubic;
  int i = 0;
  while (t < tEnd){
    printf("-----------\t t/tEnd : %.3f/%.1f\t-----------\n", t,tEnd);
    update_cells(grid, particles, n_p);
    update_neighbors(grid, particles, n_p, i);
    update_pressure(particles, n_p, rho_0, B, gamma);
    imposeFScondition(particles, n_p_dim_x, n_p_dim_y);
    time_integration_CSPM(particles, n_p, kernel, dt, edges,eta);
    show(particles, animation, i, false, false);

    printf("Time integration completed\n");
    i++;
    t += dt;
  }
  show(particles,animation, iter_max, true, true);



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

static int hydrostatic_eq(){
    double lx = 1;                          // Longueur du domaine de particule
    double ly = 1;                          // Hauteur du domaine de particle
    int n_p_dim = 25;

    // Parameters
    double rho_0 = 1e3;                     // Densité initiale
    double dynamic_viscosity = 0;        // Viscosité dynamique
    double g = 0.00;                        // Gravité
    int n_p_dim_x = n_p_dim;                // Nombre de particule par dimension
    int n_p_dim_y = n_p_dim*(ly/lx);
    int n_p = n_p_dim_x*n_p_dim_y;          // Nombre de particule total
    double delta = lx/n_p_dim_x;                // step between neighboring particles
    double h = 6*delta;      // Rayon du compact pour l'approximation
    double mass = rho_0*delta*delta;                // Masse d'une particule, constant
    double Rp = delta/2;                        // Rayon d'une particule
    double eta = 0.00;                       // XSPH parameter from 0 to 1
    double treshold = 20;                   // Critère pour la surface libre
    double tension = 0.0;//7*1e-2;                     // Tension de surface de l'eau
    double P0 = 0;                          // Pression atmosphérique
    double c = 1500;
    double gamma = 7;
    double B = c*c*rho_0/gamma;


    // ------------------------------------------------------------------
    // ------------------------ SET Particles ---------------------------
    // ------------------------------------------------------------------
    Parameters* param = Parameters_new(mass, dynamic_viscosity, h, Rp, tension, treshold,P0,g);
    Particle** particles = fluidProblem(param, n_p_dim_x, n_p_dim_y, g, rho_0, P0,true);
    set_density(particles,n_p_dim_x, n_p_dim_y, rho_0, B, gamma);
    set_artificialViscosity(param,0.015,0);

    // ------------------------------------------------------------------
    // ------------------------ SET Edges -------------------------------
    // ------------------------------------------------------------------

    double L = 1;
    double H = 1;
    int n_e = 4;
    double CF = 0.0;
    double CR = 0.0;

    double domain[4];
    Edges* edges = get_box(L,H,n_e,CF,CR,domain);
    // ------------------------------------------------------------------
    // ------------------------ SET Grid --------------------------------
    // ------------------------------------------------------------------
    double extra = 0;
    extra = 0.01;
    Grid* grid = Grid_new(0-extra, L+extra, 0-extra, H+extra, h);

    // ------------------------------------------------------------------
    // ------------------------ SET Animation ---------------------------
    // ------------------------------------------------------------------
    double timeout = 0.0001;                 // Durée d'une frame
    Animation* animation = Animation_new(n_p, timeout, grid, Rp, domain);

    // ------------------------------------------------------------------
    // ------------------------ Start integration -----------------------
    // ------------------------------------------------------------------
    double t = 0;
    double tEnd = 1;
    double dt = 1e-5;
    int iter_max = (int) (tEnd-t)/dt;
    int output = 1;
    printf("iter max = %d\n",iter_max);
    // // Temporal loop
    Kernel kernel = Cubic;
    int i = 0;
    while (t < tEnd){
      printf("-----------\t t/tEnd : %.3f/%.1f\t-----------\n", t,tEnd);
      printf("-------- Stabilité condition dt <= %e -------\n", delta*1e-1);
      update_cells(grid, particles, n_p);
      update_neighbors(grid, particles, n_p, i);
      update_pressure(particles,n_p,rho_0,B, gamma);
      // update_pressureEq(particles, n_p);
      // time_integration(particles, n_p, kernel, dt, edges);
      time_integration_CSPM(particles, n_p, kernel, dt, edges,0);
      show(particles, animation, i, false, false);

      printf("Time integration completed\n");
      i++;
      t += dt;
    }
    show(particles,animation, iter_max, false, true);



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

static int waves(){
  double lx = 2;                          // Longueur du domaine de particule
  double ly = 1.0/3;                          // Hauteur du domaine de particle
  int n_p_dim = 100;

  // Parameters
  double rho_0 = 1e3;                     // Densité initiale
  double dynamic_viscosity = 1e-3;        // Viscosité dynamique
  double g = 9.81;                        // Gravité
  int n_p_dim_x = n_p_dim;                // Nombre de particule par dimension
  int n_p_dim_y = n_p_dim*(ly/lx);
  int n_p = n_p_dim_x*n_p_dim_y;          // Nombre de particule total
  double delta = lx/n_p_dim_x;                // step between neighboring particles
  double h = sqrt(21)*lx/n_p_dim_x;      // Rayon du compact pour l'approximation
  double mass = rho_0*delta*delta;                // Masse d'une particule, constant
  double Rp = delta/2;                        // Rayon d'une particule
  double eta = 0.00;                       // XSPH parameter from 0 to 1
  double treshold = 20;                   // Critère pour la surface libre
  double tension = 7*1e-2;                     // Tension de surface de l'eau
  double P0 = 0;                          // Pression atmosphérique

  // ------------------------------------------------------------------
  // ------------------------ SET Particles ---------------------------
  // ------------------------------------------------------------------
  Parameters* param = Parameters_new(mass, dynamic_viscosity, h, Rp, tension, treshold,P0,g);
  Particle** particles = fluidProblem(param, n_p_dim_x, n_p_dim_y, g, rho_0, P0,true);
  // Apply perturbation
  double w = 10*M_PI*lx;                        // Pulsation de la perturbation si périodique
  double eps = lx/15;                           // Amplitude de la perturbation
  double sigma = 0.05;                          // Ecart-type de la perturbation
  for(int j = 0; j < n_p_dim_y; j++){
    for(int i = 0; i < n_p_dim_x; i++){
      int index = j*n_p_dim_x + i;
      Vector* x = particles[index]->fields->x;
      x->X[1] += (x->X[1]/ly) * eps * exp(-pow(x->X[0]-lx/2,2)/sigma);
    }
  }

  double B = 1e5;
  double gamma = 7;
  set_density(particles, n_p_dim_x, n_p_dim_y, rho_0, B, gamma);

  double a = 0.15;
  double b = 0;
  set_artificialViscosity(param,a,b);

  // ------------------------------------------------------------------
  // ------------------------ SET Edges -------------------------------
  // ------------------------------------------------------------------

  double L = 2;
  double H = 1;
  int n_e = 4;
  double CF = 0.0;
  double CR = 0.0;

  double domain[4];
  Edges* edges = get_box(L,H,n_e,CF,CR,domain);
  // ------------------------------------------------------------------
  // ------------------------ SET Grid --------------------------------
  // ------------------------------------------------------------------
  double extra = 0;
  extra = 0.01;
  Grid* grid = Grid_new(0-extra, L+extra, 0-extra, H+extra, h);

  // ------------------------------------------------------------------
  // ------------------------ SET Animation ---------------------------
  // ------------------------------------------------------------------
  double timeout = 0.00001;                 // Durée d'une frame
  Animation* animation = Animation_new(n_p, timeout, grid, Rp, domain);

  // ------------------------------------------------------------------
  // ------------------------ Start integration -----------------------
  // ------------------------------------------------------------------
  double t = 0;
  double tEnd = 10;
  double dt = 1e-4;
  int iter_max = (int) (tEnd-t)/dt;
  int output = 1;
  printf("iter max = %d\n",iter_max);
  // // Temporal loop
  Kernel kernel = Cubic;
  int i = 0;
  while (t < tEnd){
    printf("-----------\t t/tEnd : %.3f/%.1f\t-----------\n", t,tEnd);
    update_cells(grid, particles, n_p);
    update_neighbors(grid, particles, n_p, i);

    update_pressure(particles, n_p, rho_0, B,gamma);
    // imposeFScondition(particles, n_p_dim_x, n_p_dim_y);

    double* rhs_mass = rhs_mass_conservation(particles,n_p,kernel);
    time_integration_mass(particles,n_p,rhs_mass,dt);
    free(rhs_mass);
    if ((i+1)%30) CSPM_density(particles,n_p, kernel);

    print_rhoMax(particles,n_p);
    print_rhoMin(particles,n_p);

    Vector** rhs_momentum = CSPM_rhs_momentum_conservation(particles,n_p,kernel);
    time_integration_momentum(particles,n_p,rhs_momentum,dt);
    for(int i = 0; i < n_p; i++)
      Vector_free(rhs_momentum[i]);
    free(rhs_momentum);

    // XSPH_correction(p, n_p, kernel, eta);

    time_integration_position(particles,n_p,dt);
    reflective_boundary(particles, n_p, edges);

    show(particles, animation, i, false, false);

    printf("Time integration completed\n");
    i++;
    t += dt;
  }
  show(particles,animation, iter_max, false, true);



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
//
// static int vortices(){
//   double lx = 1;                          // Longueur du domaine de particule
//   double ly = 1;                          // Hauteur du domaine de particle
//   int n_p_dim = 50;
//
//   // Parameters
//   double rho_0 = 1e3;                     // Densité initiale
//   double dynamic_viscosity = 1e-3;        // Viscosité dynamique
//   double g = 0.0;                        // Gravité
//   int n_p_dim_x = n_p_dim;                // Nombre de particule par dimension
//   int n_p_dim_y = n_p_dim*(ly/lx);
//   int n_p = n_p_dim_x*n_p_dim_y;          // Nombre de particule total
//   double delta = lx/n_p_dim_x;                // step between neighboring particles
//   double h = sqrt(21)*lx/n_p_dim_x;      // Rayon du compact pour l'approximation
//   double mass = rho_0*delta*delta;                // Masse d'une particule, constant
//   double Rp = delta/2;                        // Rayon d'une particule
//   double eta = 0.00;                       // XSPH parameter from 0 to 1
//   double treshold = 20;                   // Critère pour la surface libre
//   double tension = 0.0;//7*1e-2;                     // Tension de surface de l'eau
//   double P0 = 0;                          // Pression atmosphérique
//
//   // ------------------------------------------------------------------
//   // ------------------------ SET Particles ---------------------------
//   // ------------------------------------------------------------------
//   Parameters* param = Parameters_new(mass, dynamic_viscosity, h, Rp, tension, treshold,P0,g);
//   Particle** particles = fluidProblem(param, n_p_dim_x, n_p_dim_y, g, rho_0, P0,true);
//
//   // Set initial conditions
//   Vector* O = Vector_new(2);
//   Vector_initialise(O,(double[2]){lx/2,ly/2});
//
//   double A = 1;
//   for(int j = 0; j < n_p_dim_y; j ++){
//     for(int i = 0; i < n_p_dim_x; i ++){
//       int index = j*n_p_dim_x + i;
//       Particle* pi = particles[index];
//       Vector* x = pi->fields->x;
//
//       double theta = atan2((x->X[1] - ly/2),(x->X[0] - lx/2));
//       // printf("theta = %f\n",theta);
//       double r = dist(O,x);
//       double u,v;
//       u = (1 - r/lx)* A * cos(theta);
//       v = (1 - r/lx)* A * sin(theta);
//
//       pi->fields->u->X[0] = u;
//       pi->fields->u->X[1] = v;
//
//       // Vector_print(pi->fields->u);
//     }
//   }
//
//
//   // ------------------------------------------------------------------
//   // ------------------------ SET Edges -------------------------------
//   // ------------------------------------------------------------------
//
//   double L = 1;
//   double H = 1;
//   int n_e = 4;
//   double CF = 0.0;
//   double CR = 1.0;
//
//   double domain[4];
//   Edges* edges = get_box(L,H,n_e,CF,CR,domain);
//   // ------------------------------------------------------------------
//   // ------------------------ SET Grid --------------------------------
//   // ------------------------------------------------------------------
//   double extra = 0;
//   extra = 0.01;
//   Grid* grid = Grid_new(0-extra, L+extra, 0-extra, H+extra, h);
//
//   // ------------------------------------------------------------------
//   // ------------------------ SET Animation ---------------------------
//   // ------------------------------------------------------------------
//   double timeout = 0.00001;                 // Durée d'une frame
//   Animation* animation = Animation_new(n_p, timeout, grid, Rp, domain);
//
//   // ------------------------------------------------------------------
//   // ------------------------ Start integration -----------------------
//   // ------------------------------------------------------------------
//   // ------------------------------------------------------------------
//   // ------------------------ Start integration -----------------------
//   // ------------------------------------------------------------------
//   double t = 0;
//   double tEnd = 10;
//   double dt = 1e-4;
//   int iter_max = (int) (tEnd-t)/dt;
//   int output = 1;
//   printf("iter max = %d\n",iter_max);
//   // // Temporal loop
//   Kernel kernel = Cubic;
//   int i = 0;
//   while (t < tEnd){
//     printf("-----------\t t/tEnd : %.3f/%.1f\t-----------\n", t,tEnd);
//     update_cells(grid, particles, n_p);
//     update_neighbors(grid, particles, n_p, i);
//
//     double* rhs_mass = rhs_mass_conservation(particles,n_p,kernel);
//     time_integration_mass(particles,n_p,rhs_mass,dt);
//     free(rhs_mass);
//     CSPM_density(particles,n_p, kernel);
//
//     print_rhoMax(particles,n_p);
//     print_rhoMin(particles,n_p);
//
//     update_pressure(particles, n_p, rho_0);
//
//     Vector** rhs_momentum = CSPM_rhs_momentum_conservation(particles,n_p,kernel);
//     time_integration_momentum(particles,n_p,rhs_momentum,dt);
//     for(int i = 0; i < n_p; i++)
//       Vector_free(rhs_momentum[i]);
//     free(rhs_momentum);
//
//     // XSPH_correction(p, n_p, kernel, eta);
//
//     time_integration_position(particles,n_p,dt);
//     reflective_boundary(particles, n_p, edges);
//
//     show(particles, animation, i, false, false);
//
//     printf("Time integration completed\n");
//     i++;
//     t += dt;
//   }
//   show(particles,animation, iter_max, false, true);
//
//
//
//   // ------------------------------------------------------------------
//   // ------------------------ FREE Memory -----------------------------
//   // ------------------------------------------------------------------
//   Particles_free(particles, n_p);
//   printf("END FREE PARTICLES\n");
//   Edges_free(edges);
//   printf("END FREE EDGES\n");
//   Grid_free(grid);
//   printf("END FREE GRID\n");
//   Animation_free(animation);
//   printf("END FREE ANIMATION\n");
//   return EXIT_SUCCESS;
// }
