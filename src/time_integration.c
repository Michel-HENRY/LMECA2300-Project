#include "time_integration.h"

// ------------------------------------------------------------------
// ------------------------ Print functions -------------------------
// ------------------------------------------------------------------

void print_vMax(Particle** p, int n_p){
  double max_y = 0;

  for(int i = 0; i < n_p; i++){
    double v = fabs(p[i]->fields->u->X[1]);
    if(max_y < v){
			max_y = v;
		}
  }
  printf("Max Velocity = %f\n", max_y);
}
void print_rhoMax(Particle** p, int n_p){
  double max_y = p[0]->param->rho;

  for(int i = 1; i < n_p; i++){
    double v = p[i]->param->rho;
    if(max_y < v){
			max_y = v;
		}
  }
  printf("Max Density = %f\n", max_y);
}
void print_rhoMin(Particle** p, int n_p){
  double min_y = p[0]->param->rho;

  for(int i = 1; i < n_p; i++){
    double v = p[i]->param->rho;
    if(min_y > v){
			min_y = v;
		}
  }
  printf("Min Density = %f\n", min_y);
}
void print_momentumMax(Vector** momentum,  int n_p){
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

void update_pressureMod(Particle** p, int n_p, double rho_0, double g, double H){
  double gamma = 7;
  double c = 1500;
  double B = rho_0*c*c/gamma;
  for(int i = 0; i < n_p ; i++){
    Particle* pi = p[i];
    double rho = pi->param->rho;
    double Pdyn = B*(pow(rho/rho_0,gamma) - 1);

    double y = pi->fields->x->X[1];
    double Phydro = rho*g*(H - y);

    pi->fields->P = Pdyn + Phydro;
  }
}
void update_pressure(Particle** p, int n_p, double rho_0){
    double gamma = 7;
    double c = 1500;
    double B = rho_0*c*c/gamma;
    for(int i = 0; i < n_p ; i++){
      Particle* pi = p[i];
      double rho = pi->param->rho;
      double Pdyn = B*(pow(rho/rho_0,gamma) - 1);

      pi->fields->P = Pdyn;
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
      double rhoi = pi->param->rho;
      double rhoj = pj->param->rho;

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
    double rho = pi->param->rho;
    double divergence_u = div_u(pi,kernel);
    // printf("div u = %f\n", divergence_u);
    rhs[i] = -rho*divergence_u;
  }
  return rhs;
}

Vector** rhs_momentum_conservation(Particle** p, int n_p, Kernel kernel){
  Vector** rhs = (Vector**) malloc(sizeof(Vector*)*n_p);
  for(int i = 0; i < n_p; i++){
    Particle* pi = p[i];
    double rho = pi->param->rho;
    double viscosity = pi->param->dynamic_viscosity;
    //TODO : Fix gradient de pression
    // printf("P = %f\n", pi->fields->P);
    Vector* grad_Pressure = grad_P(pi,kernel);
    times_into(grad_Pressure, -1/rho);
    // printf("Gradient de pression %i: \n",i);
    // Vector_print(grad_Pressure);
    Vector* laplacian_u = lapl_u_Brookshaw(pi,kernel);
    // printf("Laplacian %i :\n", i);
    // Vector_print(laplacian_u);
    times_into(laplacian_u, viscosity);
    Vector* forces = pi->fields->f;
    times_into(forces, 1/rho);

    rhs[i] = Vector_new(2);
    sum_into(rhs[i],grad_Pressure);
    sum_into(rhs[i],laplacian_u);
    sum_into(rhs[i],forces);

   Vector_free(grad_Pressure);
   Vector_free(laplacian_u);
  }
  return rhs;
}
Vector** CSPM_rhs_momentum_conservation(Particle** p, int n_p, Kernel kernel){
  Vector** rhs = (Vector**) malloc(sizeof(Vector*)*n_p);
  for(int i = 0; i < n_p; i++){
    Particle* pi = p[i];
    double rho = pi->param->rho;
    double viscosity = pi->param->dynamic_viscosity;
    Vector* grad_Pressure = CSPM_pressure(pi,kernel);
    // printf("Gradient de pression %i: \n",i);
    // Vector_print(grad_Pressure);
    Vector* laplacian_u = lapl_u_Brookshaw(pi,kernel);
    // printf("Laplacian %i :\n", i);
    // Vector_print(laplacian_u);
    times_into(laplacian_u, viscosity);
    Vector* forces = pi->fields->f;
    times_into(forces, 1/rho);

    rhs[i] = Vector_new(2);
    // sum_into(rhs[i],grad_Pressure);
    sum_into(rhs[i],laplacian_u);
    sum_into(rhs[i],forces);

   Vector_free(grad_Pressure);
   Vector_free(laplacian_u);
  }
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
      double Vj = mj/pj->param->rho;
      double W = eval_kernel(pi->fields->x, pj->fields->x, pi->param->h, kernel);
      num += mj*W;
      den += W*Vj;
      current = current->next;
    }
    rho_CSPM[i] = num/den;
  }
  for(int i = 0; i< n_p;i++){
    p[i]->param->rho = rho_CSPM[i];
  }
  free(rho_CSPM);
}

Vector* CSPM_pressure(Particle* pi, Kernel kernel){
  double rhoi = pi->param->rho;
  Vector* grad_Pressure = grad_P(pi,kernel);
  double num_x = 1e-8;
  double num_y = 1e-8;
  double den_x = 1e-8;
  double den_y = 1e-8;
  ListNode* current = pi->neighbors->head;
  while(current != NULL){
    Particle* pj = current->v;
    double Vj = pj->param->mass/pj->param->rho;
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
    pi->param->rho += update;
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
