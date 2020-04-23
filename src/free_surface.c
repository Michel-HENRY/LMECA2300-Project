#include "free_surface.h"

Vector** get_force_surface(Particle** p, int n_p, Kernel kernel){
  double* Cs = (double*) malloc(sizeof(double)*n_p);
  Vector** forces = (Vector**) malloc(sizeof(Vector*)*n_p);
  for(int i = 0;i < n_p; i++){
    Cs[i] = get_smooth_Cs(p[i], kernel);
    // printf("Cs = %f\n", Cs[i]);
  }
  for(int i = 0; i < n_p; i++){
    p[i]->fields->Cs = Cs[i];
  }
  free(Cs);
  for(int i = 0; i < n_p; i++){
    Vector* n = get_n(p[i], kernel);
    double norm_n = norm(n);
    forces[i] = Vector_new(n->DIM);
    double tension = p[i]->param->tension;
    double treshold = p[i]->param->treshold;
    if(norm_n > treshold){
      double kappa = get_curvature(p[i],norm_n,kernel);
      for(int d = 0; d < n->DIM; d++){
        forces[i]->X[d] = -(tension*kappa)*n->X[d];
      }
      // La pression a la surface libre est nulle !!!
      p[i]->fields->P = p[i]->param->P0;
    }
    Vector_free(n);
  }
  return forces;
}


static double get_curvature(Particle* p, double norm_n, Kernel kernel){
  double lapl = get_lapl_Cs(p, kernel);
  return(-lapl/norm_n);
}



static double get_lapl_Cs(Particle*pi, Kernel kernel){
  double lapl = 0;

  double h = pi->param->h;
  double fi = pi->fields->Cs;
  Vector* xi = pi->fields->x;

  ListNode* node = pi->neighbors->head;
  while(node != NULL){
    Particle* pj = node->v;
    double Vj = pj->param->mass/pj->fields->rho;
    double fj = pj->fields->Cs;
    Vector* xj = pj->fields->x;

    Vector* dW = grad_kernel(xi ,xj,h, kernel);
    Vector* xi_xj = diff(xi,xj);
    double dist_xixj = dist(xi,xj); // => Distance nulle

    double res = Vj*(fi-fj)*dot(xi_xj,dW)/(dist_xixj*dist_xixj);

    lapl += res;
    Vector_free(dW);
    Vector_free(xi_xj);

    node = node->next;
  }
  return 2*lapl;
}
static Vector* get_n(Particle*pi, Kernel kernel){
  Vector* grad = Vector_new(pi->fields->u->DIM);

  double fi = pi->fields->Cs;
  double rho_i = pi->fields->rho;
  double h = pi->param->h;
  Vector* xi = pi->fields->x;

  ListNode* node = pi->neighbors->head;
  while(node != NULL){

    // Order -1 method
    Particle* pj = node->v;
    double fj = pj->fields->Cs;
    double mj = pj->param->mass;
    double rho_j = pj->fields->rho;
    Vector* xj = pj->fields->x;

    Vector* dW = grad_kernel(xi ,xj,h, kernel);
    double a = fi/(rho_i*rho_i) + fj/(rho_j*rho_j);
    Vector* inner = times(dW,a*mj);
    sum_into(grad, inner);

    Vector_free(inner);
    Vector_free(dW);
    node = node->next;
  }
  times_into(grad,rho_i);
  return grad;
}
static double get_smooth_Cs(Particle* pi, Kernel kernel){
  double Cs = 0;
  ListNode* current = pi->neighbors->head;
  while(current != NULL){
    Particle* pj = current->v;
    double Vj = pj->param->mass/pj->fields->rho;
    // double Csj = pj->fields->Cs;     //TODO : Compare which one is better
    double Csj = 1;
    double W = eval_kernel(pi->fields->x, pj->fields->x, pi->param->h, kernel);
    Cs += Vj*Csj*W;

    current = current->next;
  }
  return Cs;
}
