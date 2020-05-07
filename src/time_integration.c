#include "time_integration.h"

// ------------------------------------------------------------------
// ------------------------ Compute fields ---- ---------------------
// ------------------------------------------------------------------
void density(Particle** p, int n_p, Kernel kernel){
  for(int i = 0; i < n_p ; i++){
    Particle* pi = p[i];
    double rho = 0;
    ListNode* current = pi->neighbors->head;
    while (current != NULL) {
      Particle* pj = current->v;
      double W = eval_kernel(pi->fields->x,pj->fields->x, pi->param->h, kernel);
      rho += pj->param->mass*W;
      current = current->next;
    }
    pi->fields->rho = rho;
  }
}
void get_dP(Particle** p, int n_p, Kernel kernel){
  for(int i = 0; i < n_p; i++){
    Vector_free(p[i]->fields->dP);
    p[i]->fields->dP = grad_P(p[i], kernel);
  }
}


// ------------------------------------------------------------------
// ------------------------ Compute all the rhs ---------------------
// ------------------------------------------------------------------

// Mass term
double* rhs_mass_conservation(Particle** p, int n_p, Kernel kernel){
  double* rhs = (double*) malloc(sizeof(double)*n_p);
  double check = 0;
  for(int i = 0; i < n_p; i ++){
    Particle* pi = p[i];
    double rho = pi->fields->rho;
    double du = div_u(pi,kernel);
    // printf("div u = %f\n", du);
    rhs[i] = -rho*du;
    check += du;
  }
  printf("SUM du = %f\n", check);
  return rhs;
}
double* CSPM_rhs_mass_conservation(Particle** p, int n_p, Kernel kernel){
  double* rhs = (double*) malloc(sizeof(double)*n_p);
  double* du = CSPM_div(p,n_p,kernel);
  for(int i = 0; i < n_p; i ++){
    Particle* pi = p[i];
    double rho = pi->fields->rho;
    // printf("div u = %f\n", du);
    rhs[i] = -rho*du[i];
  }
  free(du);
  return rhs;
}
// Momentum term
Vector** rhs_momentum_conservation(Particle** p, int n_p, Kernel kernel){
  Vector** rhs = (Vector**) malloc(sizeof(Vector*)*n_p);
  get_dP(p,n_p,kernel);
  for(int i = 0; i < n_p; i++){
    Particle* pi = p[i];
    double rho = pi->fields->rho;

    // External forces
    Vector* f = pi->fields->f;
    // Pressure Gradient
    Vector* dP = pi->fields->dP;
    rhs[i] = Vector_new(2);
    for(int d = 0; d < 2; d++){
      rhs[i]->X[d] = (-dP->X[d] + f->X[d])/rho;
    }
  }
  return rhs;
}

// ------------------------------------------------------------------
// ---------------------- Available correction ----------------------
// ------------------------------------------------------------------
double* CSPM_div(Particle** p, int n_p, Kernel kernel){
  double* dU = (double*) malloc(sizeof(double)*n_p);
  for(int i = 0;i < n_p;++i){
    Particle* pi = p[i];
    Vector* Xi = pi->fields->x;
    double rhoi = pi->fields->rho;
    double denx = 0;
    double deny = 0;
    double numx = 0;
    double numy = 0;
    double xi = Xi->X[0];
    double yi = Xi->X[1];


    ListNode* current = pi->neighbors->head;
    while(current != NULL){
      Particle* pj = current->v;
      Vector* Xj = pj->fields->x;
      double mj = pj->param->mass;
      double rhoj = pj->fields->rho;
      double xj = Xj->X[0];
      double yj = Xj->X[1];
      double h = pj->param->h;

      Vector* dW = grad_kernel(Xi,Xj,h,kernel);
      double dWx = dW->X[0];
      double dWy = dW->X[1];

      denx += dWx*(xj-xi)*mj/rhoj;
      deny += dWy*(yj-yi)*mj/rhoj;

      numx += (pj->fields->u->X[0] - pi->fields->u->X[0])*dWx*mj/rhoj;
      numy += (pj->fields->u->X[1] - pi->fields->u->X[1])*dWy*mj/rhoj;

      Vector_free(dW);
      current = current->next;
    }
    double dudx = numx/denx;
    double dvdy = numy/deny;

    dU[i] = dudx+dvdy;
  }
  return dU;
}
void CSPM_density(Particle** p, int n_p, Kernel kernel){
  double* rho_CSPM = (double*) malloc(sizeof(double)*n_p);

  for(int i = 0; i < n_p; i++){
    Particle* pi = p[i];
    double num = 0.0;
    double den = 0.0;

    ListNode *current = pi->neighbors->head;
    if(current == NULL){
      printf("CSPM_density : No neighbors for particle %d!!!\n", i);
      rho_CSPM[i] = pi->fields->rho;
      break;
    }
    while(current != NULL){
      Particle* pj = current->v;
      double mj = pj->param->mass;
      double rhoj = pj->fields->rho;
      double W = eval_kernel(pi->fields->x, pj->fields->x, pi->param->h, kernel);
      num += mj*W;
      den += W*mj/rhoj;
      current = current->next;
    }
    rho_CSPM[i] = num/den;
  }
  for(int i = 0; i< n_p;i++){
    // printf("density before = %f\t\tdensity CSPM = %f\n",p[i]->fields->rho, rho_CSPM[i]);
    p[i]->fields->rho = rho_CSPM[i];
  }
  free(rho_CSPM);
}
void CSPM_pressure(Particle** p, int n_p, Kernel kernel, Vector** dP){
  Vector** dP_CSPM = (Vector**) malloc(sizeof(Vector*)*n_p);
  for(int i = 0;i < n_p;++i){
    Particle* pi = p[i];
    Vector* Xi = pi->fields->x;
    dP_CSPM[i] = Vector_new(Xi->DIM);
    double rhoi = pi->fields->rho;
    double denx = 0;
    double deny = 0;
    double numx = 0;
    double numy = 0;
    double xi = Xi->X[0];
    double yi = Xi->X[1];


    ListNode* current = pi->neighbors->head;
    if(current == NULL){
      printf("CSPM_pressure : No neighbors for particle %d!!!\n", i);
      dP_CSPM[i] = dP[i];
      break;
    }
    while(current != NULL){
      Particle* pj = current->v;
      Vector* Xj = pj->fields->x;
      double mj = pj->param->mass;
      double rhoj = pj->fields->rho;
      double xj = Xj->X[0];
      double yj = Xj->X[1];
      double h = pj->param->h;

      Vector* dW = grad_kernel(Xi,Xj,h,kernel);
      double dWx = dW->X[0];

      double dWy = dW->X[1];

      denx += dWx*(xj-xi)*mj/rhoj;
      deny += dWy*(yj-yi)*mj/rhoj;

      numx += (pj->fields->P/rhoj - pi->fields->P/rhoi)*dWx*mj/rhoj;
      numy += (pj->fields->P/rhoj - pi->fields->P/rhoi)*dWy*mj/rhoj;

      Vector_free(dW);
      current = current->next;
    }
    // dP_CSPM[i]->X[0] = dP[i]->X[0]/denx;
    // dP_CSPM[i]->X[1] = dP[i]->X[1]/deny;
    dP_CSPM[i]->X[0] = rhoi*numx/denx;
    dP_CSPM[i]->X[1] = rhoi*numy/deny;
  }
  for(int i = 0; i < n_p; i++){
    dP[i]->X[0] = dP_CSPM[i]->X[0];
    dP[i]->X[1] = dP_CSPM[i]->X[1];
    Vector_free(dP_CSPM[i]);
  }
  free(dP_CSPM);
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

// ------------------------------------------------------------------
// --------------------- Time integration scheme ---- ---------------
// ------------------------------------------------------------------
void time_integration(Particle** p, int n_p, Kernel kernel, double dt, Edges* edges){
  Vector** rhs_momentum = rhs_momentum_conservation(p,n_p,kernel);

  time_integration_momentum(p,n_p,rhs_momentum,dt);

  time_integration_position(p,n_p,dt);
  reflective_boundary(p, n_p, edges);

  print_rhoMax(p,n_p);
  print_rhoMin(p,n_p);

  print_Pmax(p,n_p);
  print_Pmin(p,n_p);

  print_vMax(p,n_p);

  for(int i = 0; i < n_p; i++){
    Vector_free(rhs_momentum[i]);
  }
  free(rhs_momentum);
}
void time_integration_XSPH(Particle** p, int n_p, Kernel kernel, double dt, Edges* edges, double eta){
  Vector** rhs_momentum = rhs_momentum_conservation(p,n_p,kernel);

  time_integration_momentum(p,n_p,rhs_momentum,dt);
  XSPH_correction(p,n_p,kernel,eta);

  time_integration_position(p,n_p,dt);
  reflective_boundary(p, n_p, edges);

  print_rhoMax(p,n_p);
  print_rhoMin(p,n_p);

  print_Pmax(p,n_p);
  print_Pmin(p,n_p);

  print_vMax(p,n_p);

  for(int i = 0; i < n_p; i++){
    Vector_free(rhs_momentum[i]);
  }
  free(rhs_momentum);
}


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
  double max_y = p[0]->fields->rho;

  for(int i = 1; i < n_p; i++){
    double v = p[i]->fields->rho;
    if(max_y < v){
      max_y = v;
    }
  }
  printf("Max Density = %f\n", max_y);
}
void print_rhoMin(Particle** p, int n_p){
  double min_y = p[0]->fields->rho;

  for(int i = 1; i < n_p; i++){
    double v = p[i]->fields->rho;
    if(min_y > v){
      min_y = v;
    }
  }
  printf("Min Density = %f\n", min_y);
}
void print_Pmax(Particle** p, int n_p){
  double max_y = p[0]->fields->P;

  for(int i = 1; i < n_p; i++){
    double v = p[i]->fields->P;
    if(max_y < v){
      max_y = v;
    }
  }
  printf("Max Pressure = %f\n", max_y);
}
void print_Pmin(Particle** p, int n_p){
  double min_y = p[0]->fields->P;

  for(int i = 1; i < n_p; i++){
    double v = p[i]->fields->P;
    if(min_y > v){
      min_y = v;
    }
  }
  printf("Min Pressure = %f\n", min_y);
}
void print_dPminmax(Vector** dP, int n_p){
  double mindPx = dP[0]->X[0];    double mindPy = dP[0]->X[1];
  double maxdPx = dP[0]->X[0];    double maxdPy = dP[0]->X[1];

  for(int i = 1 ; i  < n_p; i++){
    if(dP[i]->X[0] > maxdPx){
      maxdPx = dP[i]->X[0];
    }
    if(dP[i]->X[1] > maxdPy){
      maxdPy = dP[i]->X[1];
    }
    if(dP[i]->X[0] < mindPx){
      mindPx = dP[i]->X[0];
    }
    if(dP[i]->X[1] < mindPy){
      mindPy = dP[i]->X[1];
    }
  }
  printf("\tMin\t\tMax\ndPx :\t %f\t\t%f\ndPy :\t %f\t\t%f\n",mindPx, maxdPx, mindPy, maxdPy);
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
void print_dP(Particle**p, int n_p,Kernel kernel){
  for(int i = 0; i < n_p; i++){
    Vector* dP = p[i]->fields->dP;
    printf("dPx[%d] = %f\t",i,dP->X[0]/p[i]->fields->rho);
    printf("dPy[%d] = %f\n",i,dP->X[1]/p[i]->fields->rho);
  }
}
