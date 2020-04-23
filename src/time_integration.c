#include "time_integration.h"

// ------------------------------------------------------------------
// ------------------------ Print functions -------------------------
// ------------------------------------------------------------------

static void print_vMax(Particle** p, int n_p){
  double max_y = 0;

  for(int i = 0; i < n_p; i++){
    double v = fabs(p[i]->fields->u->X[1]);
    if(max_y < v){
			max_y = v;
		}
  }
  printf("Max Velocity = %f\n", max_y);
}
static void print_rhoMax(Particle** p, int n_p){
  double max_y = p[0]->fields->rho;

  for(int i = 1; i < n_p; i++){
    double v = p[i]->fields->rho;
    if(max_y < v){
			max_y = v;
		}
  }
  printf("Max Density = %f\n", max_y);
}
static void print_rhoMin(Particle** p, int n_p){
  double min_y = p[0]->fields->rho;

  for(int i = 1; i < n_p; i++){
    double v = p[i]->fields->rho;
    if(min_y > v){
			min_y = v;
		}
  }
  printf("Min Density = %f\n", min_y);
}
static void print_momentumMax(Vector** momentum,  int n_p){
  double max_y = 0;
  double max_x = 0;
  for(int i = 0; i < n_p; i++){
    double qm_Y = fabs(momentum[i]->X[1]);
    double qm_X = fabs(momentum[i]->X[0]);
    printf("momentum ? %f %f\n",momentum[i]->X[0], momentum[i]->X[1]);
    if (max_y < qm_Y){
      max_y = qm_Y;
    }
    if (max_x < qm_X){
      max_x = qm_X;
    }
  }
  printf("Max Momentum = %f, %f\n", max_x, max_y);
}

// ------------------------------------------------------------------
// ------------------------ Compute all the rhs ---------------------
// ------------------------------------------------------------------
void update_pressureEq(Particle** p, int n_p){
  for(int i = 0; i < n_p ; i++){
    Particle* pi = p[i];
    pi->fields->P = pi->param->P0;
  }
}
void update_pressure(Particle** p, int n_p, double rho_0, double g, double H){
  double gamma = 7;
  double c = 1500;
  double B = rho_0*c*c/gamma;
  for(int i = 0; i < n_p ; i++){
    Particle* pi = p[i];
    double rho = pi->fields->rho;
    double Pdyn = B*(pow(rho/rho_0,gamma) - 1);

    double y = pi->fields->x->X[1];
    double Phydro = rho*g*(H - y);

    pi->fields->P = (Pdyn + pi->param->P0) + Phydro;
  }
}
void update_pressureDam(Particle** p, int n_p, double rho_0, double g, double H){

  double B = 0.85*1e5;
  double gamma = 7;
  for(int i = 0; i < n_p ; i++){
    Particle* pi = p[i];
    double rho = pi->fields->rho;
    double Pdyn = B*(pow(rho/rho_0,gamma) - 1);

    double y = pi->fields->x->X[1];
    double Phydro = rho*g*(H - y);

    pi->fields->P = (Pdyn + pi->param->P0) + Phydro;
  }
}
void update_pressureMod(Particle** p, int n_p, double rho_0){
    double gamma = 7;
    double c = 1500;
    double B = rho_0*c*c/gamma;
    for(int i = 0; i < n_p ; i++){
      Particle* pi = p[i];
      double rho = pi->fields->rho;
      double Pdyn = B*(pow(rho/rho_0,gamma) - 1);

      pi->fields->P = Pdyn + pi->param->P0;
    }
}
void XSPH_correction(Particle** p, int n_p, Kernel kernel, double eta){
  for(int i = 0; i < n_p; i++){
    Particle* pi = p[i];
    Vector* ui = pi->fields->u;
    Vector* corr = Vector_new(ui->DIM);
    ListNode* node = pi->neighbors->head;
    while(node != NULL){
      Particle* pj = node->v;
      Vector* uj = pj->fields->u;
      double W = eval_kernel(pi->fields->x, pj->fields->x, pi->param->h,kernel);
      // printf("W = %f\n", W);
      double mj = pj->param->mass;
      double rhoi = pi->fields->rho;
      double rhoj = pj->fields->rho;

      for(int d = 0; d < ui->DIM; d++){
        corr->X[d] += (2*mj)/(rhoi + rhoj)*(uj->X[d] - ui->X[d])*W;
      }
      node = node->next;
    }
    // printf("corr = \n");
    // Vector_print(corr);
    for(int d = 0; d < ui->DIM; d++){
      ui->X[d] += eta*corr->X[d];
    }

  }
}
double* rhs_mass_conservation(Particle** p, int n_p, Kernel kernel){
  double* rhs = (double*) malloc(sizeof(double)*n_p);
  for(int i = 0; i < n_p; i ++){
    Particle* pi = p[i];
    double rho = pi->fields->rho;
    double divergence_u = div_u(pi,kernel);
    // printf("div u = %f\n", divergence_u);
    rhs[i] = -rho*divergence_u;
  }
  return rhs;
}

Vector** rhs_momentum_conservation(Particle** p, int n_p, Kernel kernel){
  Vector** rhs = (Vector**) malloc(sizeof(Vector*)*n_p);
  Vector** force_surface = get_force_surface(p, n_p, kernel);
  for(int i = 0; i < n_p; i++){
    Particle* pi = p[i];
    double rho = pi->fields->rho;
    double viscosity = pi->param->dynamic_viscosity/pi->fields->rho;
    //TODO : Fix gradient de pression
    // printf("P = %f\n", pi->fields->P);
    Vector* grad_Pressure = grad_P(pi,kernel);
    times_into(grad_Pressure, -1);
    // printf("Gradient de pression %i: \n",i);
    // Vector_print(grad_Pressure);
    Vector* laplacian_u = lapl_u_Brookshaw(pi,kernel);
    // printf("Laplacian %i :\n", i);
    // Vector_print(laplacian_u);
    times_into(laplacian_u, viscosity);
    Vector* forces = pi->fields->f;
    sum_into(forces, force_surface[i]);

    rhs[i] = Vector_new(2);
    sum_into(rhs[i],grad_Pressure);
    sum_into(rhs[i],laplacian_u);
    sum_into(rhs[i],forces);

   Vector_free(grad_Pressure);
   Vector_free(laplacian_u);
   Vector_free(force_surface[i]);
  }
  free(force_surface);
  return rhs;
}
Vector** CSPM_rhs_momentum_conservation(Particle** p, int n_p, Kernel kernel){
  Vector** rhs = (Vector**) malloc(sizeof(Vector*)*n_p);
  Vector** force_surface = get_force_surface(p, n_p, kernel);
  for(int i = 0; i < n_p; i++){
    Particle* pi = p[i];
    double rho = pi->fields->rho;
    double viscosity = pi->param->dynamic_viscosity/pi->fields->rho;
    Vector* grad_Pressure = CSPM_pressure(pi,kernel);
    times_into(grad_Pressure,-1);
    // printf("Gradient de pression %i: \n",i);
    // Vector_print(grad_Pressure);
    Vector* laplacian_u = lapl_u_Brookshaw(pi,kernel);
    // printf("Laplacian %i :\n", i);
    // Vector_print(laplacian_u);
    times_into(laplacian_u, viscosity);
    Vector* forces = pi->fields->f;
    sum_into(forces, force_surface[i]);

    rhs[i] = Vector_new(2);
    sum_into(rhs[i],grad_Pressure);
    sum_into(rhs[i],laplacian_u);
    sum_into(rhs[i],forces);

   Vector_free(grad_Pressure);
   Vector_free(laplacian_u);
   Vector_free(force_surface[i]);
  }
  free(force_surface);
  return rhs;
}


void CSPM_density(Particle** p, int n_p, Kernel kernel){
  double* rho_CSPM = (double*) malloc(sizeof(double)*n_p);
  for(int i =0; i < n_p; i++){
    Particle* pi = p[i];
    double num = 1e-8;
    double den = 1e-8;

    ListNode *current = pi->neighbors->head;
    while(current != NULL){
      Particle* pj = current->v;
      double mj = pj->param->mass;
      double Vj = mj/pj->fields->rho;
      double W = eval_kernel(pi->fields->x, pj->fields->x, pi->param->h, kernel);
      num += mj*W;
      den += W*Vj;
      current = current->next;
    }
    rho_CSPM[i] = num/den;
  }
  for(int i = 0; i< n_p;i++){
    p[i]->fields->rho = rho_CSPM[i];
  }
  free(rho_CSPM);
}

Vector* CSPM_pressure(Particle* pi, Kernel kernel){
  double rhoi = pi->fields->rho;
  Vector* grad_Pressure = grad_P(pi,kernel);
  double num_x = 1e-8;
  double num_y = 1e-8;
  double den_x = 1e-8;
  double den_y = 1e-8;
  ListNode* current = pi->neighbors->head;
  while(current != NULL){
    Particle* pj = current->v;
    double Vj = pj->param->mass/pj->fields->rho;
    Vector* dW = grad_kernel(pi->fields->x, pj->fields->x, pi->param->h, kernel);
    double xjxi = pj->fields->x->X[0] - pi->fields->x->X[0];
    double yjyi = pj->fields->x->X[1] - pi->fields->x->X[1];
    den_x += Vj * dW->X[0] * xjxi;
    den_y += Vj * dW->X[1] * yjyi;

    Vector_free(dW);
    current = current->next;
  }
  num_x = grad_Pressure->X[0]/rhoi;
  num_y = grad_Pressure->X[1]/rhoi;

  Vector* dP = Vector_new(grad_Pressure->DIM);
  dP->X[0] = num_x/den_x;
  dP->X[1] = num_y/den_y;

  Vector_free(grad_Pressure);
  return dP;
}



// ------------------------------------------------------------------
// ------------------- Solve the material derivatives ---------------
// ------------------------------------------------------------------

void time_integration_mass(Particle** p, int n_p,double* rhs_mass,double dt){
  for(int i = 0; i < n_p; i++){
    Particle* pi = p[i];
    double update = rhs_mass[i]*dt;
    pi->fields->rho += update;
  }
}
void time_integration_momentum(Particle** p, int n_p, Vector** rhs_momentum,double dt){
  for(int i = 0; i < n_p; i++){
    Particle* pi = p[i];
    Vector* update = times(rhs_momentum[i],dt);
    sum_into(pi->fields->u, update);
    Vector_free(update);
  }
}
void time_integration_position(Particle** p, int n_p, double dt){
  for(int i = 0; i < n_p; i++){
    Particle* pi = p[i];
    Vector* update = times(pi->fields->u, dt);
    sum_into(pi->fields->x, update);
    Vector_free(update);
  }
}

void time_integration(Particle** p, int n_p, Kernel kernel, double dt, Edges* edges){
  Vector** rhs_momentum = rhs_momentum_conservation(p,n_p,kernel);
  double* rhs_mass = rhs_mass_conservation(p,n_p,kernel);
  time_integration_mass(p,n_p,rhs_mass,dt);
  time_integration_momentum(p,n_p,rhs_momentum,dt);

  // Check boundary;
  time_integration_position(p,n_p,dt);
  reflective_boundary(p, n_p, edges);

  print_rhoMax(p,n_p);
  print_rhoMin(p,n_p);

  for(int i = 0; i < n_p; i++){
    Vector_free(rhs_momentum[i]);
  }
  free(rhs_momentum);
  free(rhs_mass);
}

void time_integration_CSPM(Particle** p, int n_p, Kernel kernel, double dt, Edges* edges, double eta){

  double* rhs_mass = rhs_mass_conservation(p,n_p,kernel);
  time_integration_mass(p,n_p,rhs_mass,dt);
  CSPM_density(p,n_p, kernel);

  Vector** rhs_momentum = CSPM_rhs_momentum_conservation(p,n_p,kernel);
  time_integration_momentum(p,n_p,rhs_momentum,dt);

  // XSPH_correction(p, n_p, kernel, eta);

  time_integration_position(p,n_p,dt);
  reflective_boundary(p, n_p, edges);

  print_rhoMax(p,n_p);
  print_rhoMin(p,n_p);

  for(int i = 0; i < n_p; i++){
    Vector_free(rhs_momentum[i]);
  }
  free(rhs_momentum);
  free(rhs_mass);
}

void time_integration_XSPH(Particle** p, int n_p, Kernel kernel, double dt, Edges* edges, double eta){
  Vector** rhs_momentum = rhs_momentum_conservation(p,n_p,kernel);
  double* rhs_mass = rhs_mass_conservation(p,n_p,kernel);
  time_integration_mass(p,n_p,rhs_mass,dt);
  time_integration_momentum(p,n_p,rhs_momentum,dt);

  XSPH_correction(p, n_p, kernel, eta);

  time_integration_position(p,n_p,dt);
  reflective_boundary(p, n_p, edges);

  for(int i = 0; i < n_p; i++){
    Vector_free(rhs_momentum[i]);
  }
  free(rhs_momentum);
  free(rhs_mass);
}

// int main(){
//   //particle domain
//   double l = 1;                         // Longueur du domaine de particule
//   double ly = l;                          // Hauteur du domaine de particle
//
//   // Parameters
//   double rho_0 = 1e3;                     // Densité initiale
//   double dynamic_viscosity = 1e-6;        // Viscosité dynamique
//   double g = 9.81;                        // Gravité
//   int n_p_dim = 50;                       // Nombre de particule par dimension
//   int n_p = n_p_dim*n_p_dim;              // Nombre de particule total
//   double h = l/n_p_dim;                   // step between neighboring particles
//   double kh = sqrt(21)*l/n_p_dim;         // Rayon du compact pour l'approximation
//   double mass = rho_0 * h*h;              // Masse d'une particule, constant
//   double Rp = h/2;                        // Rayon d'une particule
//   double eta = 0.5;                       // XSPH parameter from 0 to 1
//   double treshold = 20;                   // Critère pour la surface libre
//   double tension = 72*1e-3;               // Tension de surface de l'eau
//   double P0 = 0;                        // Pression atmosphérique
//
//   // ------------------------------------------------------------------
//   // ------------------------ SET Particles ---------------------------
//   // ------------------------------------------------------------------
//   Particle** particles = (Particle**) malloc(n_p*sizeof(Particle*));
//   for(int i = 0; i < n_p_dim; i++){
//     for(int j = 0; j < n_p_dim; j++){
//       int index = i*n_p_dim + j;
//       Parameters* param = Parameters_new(rho_0, mass, dynamic_viscosity, kh, Rp, tension, treshold,P0);
//       Vector* x = Vector_new(2);
//       Vector* u = Vector_new(2);
//       Vector* f = Vector_new(2);
//
//       // f->X[1] = -g;
//
//       double pos[2] = {Rp + i*h ,Rp + j*h};
//       double P = 0;
//       // if(i == 0){
//         // u->X[0] = -2*(1 - pos[1]*pos[1]);
//       // }
//
//       Vector_initialise(x,pos);
//       Fields* fields = Fields_new(x,u,f,P);
//       particles[index] = Particle_new(param, fields);
//     }
//   }
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
//   Vector** vertices = (Vector**) malloc(n_e*sizeof(vertices));
//   for(int i = 0; i < n_e; i++){
//     vertices[i] = Vector_new(2);
//   }
//   vertices[0]->X[0] = 0;                vertices[2]->X[0] = L;
//   vertices[0]->X[1] = 0;                vertices[2]->X[1] = H;
//   vertices[1]->X[0] = L;                vertices[3]->X[0] = 0;
//   vertices[1]->X[1] = 0;                vertices[3]->X[1] = H;
//
//
//   Vector** edge = (Vector**) malloc(n_e*2*sizeof(vertices));
//   for(int i = 0;i < n_e; i++){
//     edge[2*i] = vertices[i];
//     edge[2*i+1] = vertices[(i+1)%n_e];
//   }
//   Edges* edges = Edges_new(n_e, edge, CR,CF);
//   double domain[4] = {0,L,0,H};
//   // ------------------------------------------------------------------
//   // ------------------------ SET Grid --------------------------------
//   // ------------------------------------------------------------------
//   double extra = 0.0;
//   extra = 0.5;
//   Grid* grid = Grid_new(0-extra, L+extra, 0-extra, H+extra, kh);
//
//   // ------------------------------------------------------------------
//   // ------------------------ SET Animation ---------------------------
//   // ------------------------------------------------------------------
//   double timeout = 0.001;                 // Durée d'une frame
//   Animation* animation = Animation_new(n_p, timeout, grid, Rp, domain);
//
//   // ------------------------------------------------------------------
//   // ------------------------ Start integration -----------------------
//   // ------------------------------------------------------------------
//   double t = 0;
//   double tEnd = 5;
//   double dt = 0.0001;
//   int iter_max = (int) (tEnd-t)/dt;
//   int output = 1;
//   printf("iter max = %d\n",iter_max);
//   // // Temporal loop
//   Kernel kernel = Cubic;
//   int i = 0;
//   while (t < tEnd){
//     printf("-----------\t t/tEnd : %.3f/%.1f\t-----------\n", t,tEnd);
//
//     update_cells(grid, particles, n_p);
//     update_neighbors(grid, particles, n_p, i);
//     update_pressureMod(particles, n_p, rho_0,g, H,P0);
//     // update_pressure(particles, n_p, rho_0);
//     printf("P = %f\n",particles[0]->fields->P);
//     // time_integration(particles, n_p, kernel, dt, edges);
//     time_integration_CSPM(particles, n_p, kernel, dt, edges,eta);
//     if (i%output == 0)
//       show(particles, animation, i, false, false);
//     printf("Time integration completed\n");
//
//     i++;
//     t += dt;
//   }
//   show(particles,animation, iter_max, true, false);
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
//   // Animation_free(animation);
//   // printf("END FREE ANIMATION\n");
//   return EXIT_SUCCESS;
// }
