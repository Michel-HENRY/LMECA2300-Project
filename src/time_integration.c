#include "time_integration.h"

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
    times_into(grad_Pressure, -1/rho);
    Vector* laplacian_u = lapl_u(pi,kernel);
    times_into(laplacian_u, viscosity);
    Vector* forces = pi->fields->f;

    rhs[i] = 0
    sum_into(rhs,grad_Pressure);
    sum_into(rhs,laplacian_u);
    sum_into(rhs,forces);

   Vector_free(grad_Pressure);
   Vector_free(laplacian_u);
  }
  return rhs;
}


void time_integration(Particle** p, int n_p,,Kernel kernel, double dt){
  Vector** rhs_momentum = rhs_momentum_conservation(p,n_p,kernel);
  double* rhs_mass = rhs_mass_conservation(p,n_p,kernel);

  time_integration_mass(p,n_p,rhs_mass,dt);
  time_integration_momentum(p,n_p,rhs_momentum,dt);
  time_integration_position(p,n_p,dt);

  for(int i = 0; i < n_p; i++){
    Vector_free(rhs_momentum[i]);
  }
  free(rhs_momentum);
  free(rhs_mass);
}
void time_integration_mass(Particle** p, int n_p,double* rhs_mass,double dt){
  for(int i = 0; i < n_p; i++){
    Particle* pi = p[i];
    double update = rhs_mass[i]*dt;
    p->Parameters->rho += update;
  }
}
void time_integration_momentum(Particle** p, int n_p, Vector** rhs_momentum,double dt);
{
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
