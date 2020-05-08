#include "SPH_operator.h"

double div_u(Particle* pi, Kernel kernel){
  double du = 0;

  double rhoi = pi->fields->rho;
  double kh = pi->param->h;
  Vector* ui = pi->fields->u;
  Vector* xi = pi->fields->x;

  ListNode* node = pi->neighbors->head;
  while(node != NULL){
    Particle* pj  = node->v;
    double mj     = pj->param->mass;
    Vector* uj    = pj->fields->u;
    Vector* xj    = pj->fields->x;

    Vector* dW    = grad_kernel(xi,xj,kh,kernel);
    du           += div_local(ui,uj,dW,mj);

    Vector_free(dW);
    node = node->next;
  }
  return du/rhoi;
}

double div_local(Vector* fi, Vector* fj, Vector* dWij, double mj){
  Vector* fj_fi = diff(fj,fi);
  double res    = mj*dot(fj_fi,dWij);
  Vector_free(fj_fi);
  return res;
}

Vector* grad_P(Particle* pi, Kernel kernel){
  Vector* dP = Vector_new(pi->fields->u->DIM);

  double Pi       = pi->fields->P;
  double rhoi     = pi->fields->rho;
  double kh       = pi->param->h;
  Vector* xi      = pi->fields->x;

  ListNode* node = pi->neighbors->head;
  while(node != NULL){
    Particle* pj    = node->v;
    double Pj       = pj->fields->P;
    double mj       = pj->param->mass;
    double rhoj     = pj->fields->rho;
    Vector* xj      = pj->fields->x;

    Vector* dW = grad_kernel(xi ,xj,kh, kernel);
    Vector* inner = grad_local(Pi,Pj,dW, mj, rhoi,rhoj);
    sum_into(dP, inner);

    Vector_free(inner);
    Vector_free(dW);

    node = node->next;
  }
  times_into(dP, rhoi);
  return dP;
}
Vector* grad_local(double fi, double fj, Vector* dWij, double mj, double rhoi,double rhoj){
  Vector* res = times(dWij, mj * (fi/(rhoi*rhoi) + fj/(rhoj*rhoj)));
  return res;
}

Vector* lapl_book(Particle* pi, Kernel kernel){
  double h = pi->param->h;
  Vector* xi = pi->fields->x;
  Vector* ui = pi->fields->u;
  Vector* ddu = Vector_new(ui->DIM);

  ListNode* current = pi->neighbors->head;
  while(current != NULL){
    Particle* pj = current->v;
    Vector* xj = pj->fields->x;
    Vector* uj = pj->fields->u;

    double mj = pj->param->mass;
    double rhoj = pj->fields->rho;
    double eta = 0.01*h*h;

    Vector* dW = grad_kernel(xi,xj,h,kernel);
    Vector* xij = diff(xi,xj);
    double r = norm(xij);

    Vector* inner = diff(ui,uj);
    times_into(inner, dot(dW,xij)*mj/(rhoj*(r*r + eta*eta)));

    sum_into(ddu, inner);
    current = current->next;

    Vector_free(dW);
    Vector_free(xij);
    Vector_free(inner);
  }
  times_into(ddu, 2);
  return(ddu);
}

Vector* lapl_u(Particle* pi, Kernel kernel){
  Vector* lapl = Vector_new(pi->fields->u->DIM);

  double h = pi->param->h;
  Vector* fi = pi->fields->u;
  Vector* xi = pi->fields->x;

  ListNode* node = pi->neighbors->head;
  while(node != NULL){
    Particle* pj = node->v;
    double mj = pj->param->mass;
    double rhoj = pj->fields->rho;
    Vector* fj = pj->fields->u;
    Vector* xj = pj->fields->x;

    Vector* dW = grad_kernel(xi ,xj,h, kernel);
    Vector* inner = lapl_local(fi,fj,xi,xj,dW, mj/rhoj,h);
    sum_into(lapl,inner);
    Vector_free(dW);
    Vector_free(inner);

    node = node->next;
  }
  times_into(lapl,2*(fi->DIM + 2));
  return lapl;
}

Vector* lapl_local(Vector* fi, Vector* fj, Vector* xi, Vector* xj, Vector* dWij, double Vj, double h){
  Vector* fi_fj = diff(fi,fj);
  Vector* rij = diff(xi,xj);
  double norm_r = norm(rij);
  Vector* res = times(dWij, Vj*dot(fi_fj, rij)/(norm_r*norm_r + 0.01*h*h));
  Vector_free(fi_fj);
  Vector_free(rij);
  return res;
}
Vector* lapl_u_shao(Particle* pi, Kernel kernel){ ////// Il y a déjà le *mu dans l'expression !!!
  Vector* lapl = Vector_new(pi->fields->u->DIM);

  double h = pi->param->h;
  double mui = pi->param->dynamic_viscosity;
  double eta = 1e-6;
  double rhoi = pi->fields->rho;
  Vector* fi = pi->fields->u;
  Vector* xi = pi->fields->x;



  ListNode* node = pi->neighbors->head;
  while(node != NULL){
    Particle* pj = node->v;
    double mj = pj->param->mass;
    double rhoj = pj->fields->rho;
    double muj = pj->param->dynamic_viscosity;
    Vector* fj = pj->fields->u;
    Vector* xj = pj->fields->x;

    Vector* dW = grad_kernel(xi,xj,h,kernel);
    Vector* fi_fj = diff(fi,fj);
    Vector* xi_xj = diff(xi,xj);
    double dist_xixj = dist(xi,xj);
    double num = 4*mj*(mui+muj)*dot(xi_xj, dW);
    double den = (rhoi + rhoj)*(rhoi + rhoj)*(dist_xixj*dist_xixj + eta*eta);
    Vector* inner = times(fi_fj,num/den);

    sum_into(lapl,inner);

    Vector_free(dW);
    Vector_free(fi_fj);
    Vector_free(xi_xj);
    Vector_free(inner);

    node = node->next;
  }
  return lapl;
}
