#include "derivatives.h"

double compute_div(Particle* pi, xy_getter get, Kernel kernel, double kh) {
    double div = 0;
    xy *fi = get(pi);
    ListNode *node = pi->neighborhood->head;
    while(node != NULL) {
        Particle *pj = node->v;

        xy *grad_W = grad_kernel(pi->pos, pj->pos, kh, kernel);
        // correct_grad_local(grad_W, pi, pj, kh, kernel);

        xy *fj = get(pj);
        div += ((fj->x - fi->x) * grad_W->x)*pj->m;
        div += ((fj->y - fi->y) * grad_W->y)*pj->m;
        free(grad_W);
        node = node->next;
    }
    div /= pi->rho;
    return div;
}

void compute_grad(Particle* pi, scalar_getter get, Kernel kernel, double kh,xy* grad) {

	double gx = 0;
	double gy = 0;
  double fi = get(pi);
  ListNode *node = pi->neighborhood->head;
  //printf("Computing gradient of (%lf, %lf), fi = %lf\n", particle->pos->x, particle->pos->y, fi);
  while(node != NULL) {
      Particle *pj = node->v;
      xy *grad_W = grad_kernel(pi->pos, pj->pos, kh, kernel);

      double fj = get(pj);
      //printf("Position of pj: (%lf, %lf), fj = %lf\n", pj->pos->x, pj->pos->y, fj);
      gx += pi->rho * pj->m * (fi/squared(pi->rho) + fj/squared(pj->rho)) * grad_W->x; // sign is not the same as in the def...
      gy += pi->rho * pj->m * (fi/squared(pi->rho) + fj/squared(pj->rho)) * grad_W->y;
      free(grad_W);
      //printf("grad = (%lf, %lf), fj = %lf\n", grad->x, grad->y, fj);
      node = node->next;
  }
  // xy *cg = xy_new(gx,gy);
  grad->x = gx;
  grad->y = gy;
  // correct_grad(grad, pi, kh, kernel);

  // free(cg);
}

double compute_lapl(Particle *pi, scalar_getter get, Kernel kernel, double kh) {
    double lapl = 0;
    double fi = get(pi);
    ListNode *node = pi->neighborhood->head;
    while(node != NULL) {
        Particle *pj = node->v;
        xy *grad_W = grad_kernel(pi->pos, pj->pos, kh, kernel);
        // correct_grad_local(grad_W, pi, pj, kh, kernel);

        double fj = get(pj);
        double d2 = squared(pi->pos->x - pj->pos->x) + squared(pi->pos->y - pj->pos->y); // squared distance between particles
        // xy *DXij = xy_new((pi->pos->x - pj->pos->x) / d2, (pi->pos->y - pj->pos->y) / d2); // Delta X_{ij}
        xy *DXij = xy_new((pi->pos->x - pj->pos->x) / d2, (pi->pos->y - pj->pos->y) / d2); // WARNING
        if(d2 != 0) lapl += 2 * (pj->m/pj->rho) * (fi - fj) * (DXij->x * grad_W->x + DXij->y * grad_W->y);
        free(grad_W);
        free(DXij);
        node = node->next;
    }
    return lapl;
}
