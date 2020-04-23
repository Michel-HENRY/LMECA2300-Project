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
static void XSPH_correction(Particle** p, int n_p, Kernel kernel, double eta){
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
static double* rhs_mass_conservation(Particle** p, int n_p, Kernel kernel){
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

static Vector** rhs_momentum_conservation(Particle** p, int n_p, Kernel kernel){
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
    Vector* laplacian_u = lapl_u(pi,kernel);
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
static Vector** CSPM_rhs_momentum_conservation(Particle** p, int n_p, Kernel kernel){
  Vector** rhs = (Vector**) malloc(sizeof(Vector*)*n_p);
  Vector** force_surface = get_force_surface(p, n_p, kernel);
  for(int i = 0; i < n_p; i++){
    Particle* pi = p[i];
    double rho = pi->fields->rho;
    double viscosity = pi->param->dynamic_viscosity/pi->fields->rho;

    // Pressure Gradient
    Vector* grad_Pressure = CSPM_pressure(pi,kernel);
    times_into(grad_Pressure,-1);
    // printf("Gradient de pression %i: \n",i);
    // Vector_print(grad_Pressure);

    // Viscosity forces
    Vector* laplacian_u = lapl_u(pi,kernel);
    times_into(laplacian_u, viscosity);
    // printf("Laplacian %i :\n", i);
    // Vector_print(laplacian_u);

    // Surfaces forces
    Vector* forces = pi->fields->f;
    sum_into(forces, force_surface[i]);

    // Artificial viscosity
    double a = 0.3;
    double b = 0;
    Vector* pij = get_Pi_ij(pi,a,b,kernel);
    times_into(pij, -1);

    rhs[i] = Vector_new(2);
    sum_into(rhs[i],grad_Pressure);
    sum_into(rhs[i],laplacian_u);
    sum_into(rhs[i],forces);
    sum_into(rhs[i],pij);

   Vector_free(grad_Pressure);
   Vector_free(laplacian_u);
   Vector_free(force_surface[i]);
   Vector_free(pij);
  }
  free(force_surface);
  return rhs;
}


static void CSPM_density(Particle** p, int n_p, Kernel kernel){
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

static Vector* CSPM_pressure(Particle* pi, Kernel kernel){
  double rhoi = pi->fields->rho;
  Vector* grad_Pressure = grad_P(pi,kernel);

  // static int counter = 0;
  // printf("counter = %d\n", counter);
  // if (counter%3000){
    double num_x = 1e-8;    double den_x = 1e-8;
    double num_y = 1e-8;    double den_y = 1e-8;
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

    grad_Pressure->X[0] = num_x/den_x;
    grad_Pressure->X[1] = num_y/den_y;
  // }
  //
  // counter++;
  return grad_Pressure;
}

static Vector* get_Pi_ij(Particle* pi, double a, double b, Kernel kernel){
  double h = pi->param->h;
  double c = 1500;
  double phi = 0.01*h;
  double rhoi = pi->fields->rho;
  Vector* xi = pi->fields->x;
  Vector* vi = pi->fields->u;

  Vector* Pi_i = Vector_new(xi->DIM);

  ListNode* current = pi->neighbors->head;
  while(current!=NULL){
    Particle* pj = current->v;
    Vector* vj = pj->fields->u;
    Vector* xj = pj->fields->x;

    Vector* vij = diff(vi,vj);
    Vector* xij = diff(xi,xj);
    Vector* dW  = grad_kernel(xi,xj,h,kernel);


    if(dot(vij,xij) < 0){
      double rhoj = pj->fields->rho;
      double rhoij = (rhoi + rhoj)/2;
      double X = h*dot(vij,xij)/(dist(xi,xj)*dist(xi,xj) + phi*phi);

      double pij = (-a*c*X + b*X*X)/rhoij;
      double mj = pj->param->mass;

      times_into(dW, mj*pij);
      sum_into(Pi_i, dW);
    }

    Vector_free(dW);
    Vector_free(vij);
    Vector_free(xij);

    current = current->next;
  }
  return Pi_i;
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

  XSPH_correction(p, n_p, kernel, eta);

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
