#include "boundary.h"

Edges* Edges_new(int n_e, Vector** edge, double CR, double CF){
  Edges* edges = (Edges*) malloc(sizeof(Edges));
  edges->n_e = n_e;
  edges->edge = edge;
  edges->CR = CR;
  edges->CF = CF;
  set_normal(edges);

  return(edges);
}

void Edges_free(Edges* edges){
  for(int i = 0; i < edges->n_e; i ++){
    Vector_free(edges->edge[i]);
    Vector_free(edges->n[i]);
  }
  free(edges->edge);
  free(edges->n);
  free(edges);
}

void set_normal(Edges* edges){
  Vector** n = malloc(sizeof(Vector**)*edges->n_e);
  for(int i = 0; i < edges->n_e; i ++){

    Vector* e0 = edges->edge[2*i];
    Vector* e1 = edges->edge[2*i+1];

    double ni[2];
    ni[0] = e1->X[1] - e0->X[1]; //dy
    ni[1] = e0->X[0] - e1->X[0]; //-dx
    n[i] = Vector_new(2);
    Vector_initialise(n[i], ni);
    times_into(n[i], 1/(norm(n[i]))); // Normalisation

    // printf("normale [%d] = \t", i);
    // Vector_print(n[i]);
  }
  edges->n = n;
}

bool isInside(Vector* C1, Edges* edges){
  for(int i = 0; i < edges->n_e; i++){
    Vector* e0 = edges->edge[2*i];
    Vector* e1 = edges->edge[2*i+1];

    double xa,xb,xc,ya,yb,yc;

    xa = e0->X[0];    ya = e0->X[1];
    xb = e1->X[0];    yb = e1->X[1];
    xc = C1->X[0];    yc = C1->X[1];

    if((xa-xb)*(yc-ya) > (ya-yb)*(xc-xa)){
      return false;
    }
  }
  return true;
}

double* distEdge(Vector* C1, Edges* edges){
  double* d_e = (double*) malloc(sizeof(double) * edges->n_e);
  for(int i = 0; i < edges->n_e; ++i){
    Vector* e0 = edges->edge[2*i];
    Vector* e1 = edges->edge[2*i+1];

    double A,B,C;
    A = -(e1->X[1] - e0->X[1]);
    B = (e1->X[0] - e0->X[0]);
    C = e1->X[0]*(e1->X[1] - e0->X[1]) - e1->X[1]*(e1->X[0] - e0->X[0]);

    double res = A*C1->X[0] + B*C1->X[1] + C;
    res /= sqrt(A*A + B*B);

    d_e[i] = fabs(res);
  }
  return d_e;
}

int indexCPlane(double* dist_edges, int n_e){
  int min = 0;
  for(int i = 1; i < n_e; i++){
    if(dist_edges[i] < dist_edges[min]){
      min = i;
    }
  }
  return min;
}

void update_mass_center(Vector* C1, double Rp, Edges* edges, double d, int index){

  Vector* n = edges->n[index];
  Vector* t = Vector_new(C1->DIM);

  t->X[0] = -n->X[1]; t->X[1] = n->X[0];

  double C1n = dot(C1,n);
  double C1t = dot(C1,t);

  double CR = edges->CR;
  C1n -= (1+CR)*(d-Rp);

  for(int i = 0; i < C1->DIM; i++){
    C1->X[i] = -C1n * n->X[i] + C1t * t->X[i]; // On change le signe perpeniculaire à la frontière
  }
}

void update_velocity(Particle* p, Edges* edges, int index){

  Vector* n = edges->n[index];
  Vector* t = Vector_new(p->fields->u->DIM);
  t->X[0] = -n->X[1]; t->X[1] = n->X[0];

  Vector* u = p->fields->u;

  double CR = edges->CR;
  double CF = edges->CF;
  double un = dot(u,n)*CR;
  double ut = dot(u,t)*(1-CF);

  for(int i = 0; i < u->DIM; i++){
    u->X[i] = -un*n->X[i] + ut* t->X[i];
  }
}

void reflective_boundary(Particle** p, int n_p, double dt, Edges* edges){
  for(int i = 0; i < n_p; i++){
    Particle* pi = p[i];

    int counter = 0;
    while(true){
      Vector* C1 = times(pi->fields->u,dt); // Compute u*dt
      sum_into(C1,pi->fields->x); // compute x_prev + u*dt

      if(isInside(C1, edges)){
        // printf("Inside !\n");
        break;
      }
      double* dist_edges = distEdge(C1, edges);
      int index = indexCPlane(dist_edges, edges->n_e);
      // update mass center and correct velocity
      double Rp = pi->param->Rp;
      update_mass_center(C1,Rp,edges, dist_edges[index], index);
      update_velocity(pi,edges,index);

      Vector_free(C1);
      free(dist_edges);
      counter++;
      if(counter > 1000){break;}
    }
  }
}
