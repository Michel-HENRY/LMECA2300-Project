#include "time_integration.h"




void update_pressure(Particle** p, int n_p, double rho_0, double g, double H){
  double gamma = 7;
  double c = 1500;
  double B = rho_0*c*c/gamma;
  for(int i = 0; i < n_p ; i++){
    Particle* pi = p[i];
    double rho = pi->param->rho;
    double Pdyn = B*(pow(rho/rho_0,gamma) - 1);

    double y = pi->fields->x->X[1];
    double Phydro = -rho*g*(H - y);

    pi->fields->P = Pdyn + Phydro;
  }
}
double* rhs_mass_conservation(Particle** p, int n_p, Kernel kernel){
  double* rhs = (double*) malloc(sizeof(double)*n_p);
  for(int i = 0; i < n_p; i ++){
    Particle* pi = p[i];
    double rho = pi->param->rho;
    double divergence_u = div_u(pi,kernel);

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

    Vector* grad_Pressure = grad_P(pi,kernel);
    // printf("Gradient de pression %i: \n",i);
    // Vector_print(grad_Pressure);
    times_into(grad_Pressure, -1/rho);
    Vector* laplacian_u = lapl_u(pi,kernel);
    times_into(laplacian_u, viscosity);
    Vector* forces = pi->fields->f;

    rhs[i] = Vector_new(2);
    sum_into(rhs[i],grad_Pressure);
    sum_into(rhs[i],laplacian_u);
    sum_into(rhs[i],forces);

   Vector_free(grad_Pressure);
   Vector_free(laplacian_u);
  }
  return rhs;
}

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
  // time_integration_mass(p,n_p,rhs_mass,dt);
  // time_integration_momentum(p,n_p,rhs_momentum,dt);

  // Check boundary;
  time_integration_position(p,n_p,dt);
  reflective_boundary(p, n_p, dt, edges);

  for(int i = 0; i < n_p; i++){
    Vector_free(rhs_momentum[i]);
  }
  free(rhs_momentum);
  free(rhs_mass);
}
