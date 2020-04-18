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

bool isInside(Vector* C1, Edges* edges, double Rp){
  for(int i = 0; i < edges->n_e; i++){
    Vector* e0 = edges->edge[2*i];
    Vector* e1 = edges->edge[2*i+1];

    double xa,xb,xc,ya,yb,yc;

    xa = e0->X[0];    ya = e0->X[1];
    xb = e1->X[0];    yb = e1->X[1];
    xc = C1->X[0]+Rp*edges->n[i]->X[0];    yc = C1->X[1]+Rp*edges->n[i]->X[1];

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
    A = e1->X[1] - e0->X[1];
    B = e0->X[0] - e1->X[0];
    C = e0->X[1]*e1->X[0] - e1->X[1]*e0->X[0];

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
  C1n -= (1+CR)*(Rp+d); // Pas vrmt d'explication pour le -Rp au lieu de Rp
  // C1n -= (1+CR)*d;
  for(int i = 0; i < C1->DIM; i++){
    C1->X[i] = C1n*n->X[i] + C1t*t->X[i]; // On change le signe perpeniculaire à la frontière
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
    // printf("u_corr[%d] = %f\n",i, u->X[i]);
  }
}

void reflective_boundary(Particle** p, int n_p, double dt, Edges* edges){
  for(int i = 0; i < n_p; i++){
    Particle* pi = p[i];
    int counter = 0;
    Vector* C1 = NULL;
    while(true){
      C1 = pi->fields->x;
      double Rp = pi->param->Rp;
      if(isInside(C1, edges, Rp)){
        // printf("Inside !\n");
        break;
      }
      double* dist_edges = distEdge(C1, edges);
      // for(int j = 0; j < edges->n_e; j++){
      //   // printf("d[%d] = %.3f\t", j, dist_edges[j]);
      // }
      int index = indexCPlane(dist_edges, edges->n_e);
      // update mass center and correct velocity

      update_mass_center(C1,Rp,edges, dist_edges[index], index);
      update_velocity(pi,edges,index);

      free(dist_edges);
      counter++;
      if(counter > 10){break;}
    }
  }
}
// YOU NEED TO COPY THIS INTO THE MAIN FUNCTION OR TO ADD ALL THE NECESSARY INCLUDES.
// void validation_boundary(){
//   //particle domain
//   double l = 1.0;
//   double ly = l;
//   // Grid definition
//
//   double timeout = 1;
//
//   // Parameters
//   double rho = 1e3;
//   double dynamic_viscosity = 1e-3;
//   double g = 0;
//   int n_p_dim = 1;
//   int n_p = n_p_dim*n_p_dim;
//   // double h = l/n_p_dim; // step between neighboring particles
//   double h = 0.05;
//   double kh = sqrt(21) * l / n_p_dim;
//   double mass = rho * h*h;
//   double Rp = h/2;
//
//   // ------------------------------------------------------------------
//   // ------------------------ SET Particles ---------------------------
//   // ------------------------------------------------------------------
//   Particle** particles = (Particle**) malloc(n_p*sizeof(Particle*));
//
//   for(int i = 0; i < n_p_dim; i++){
//     for(int j = 0; j < n_p_dim; j++){
//       int index = i*n_p_dim + j;
//       Parameters* param = Parameters_new(rho, mass, dynamic_viscosity, kh, Rp);
//       Vector* x = Vector_new(2);
//       Vector* u = Vector_new(2);
//       Vector* f = Vector_new(2);
//       u->X[0] = 1.3/2;
//       u->X[1] = 0.7/2;
//       double  P = rho*g*(l - j*h);
//       // double P = 0;
//       // double pos[2] = {Rp + i*h,Rp + j*h};
//       double pos[2] = {l/2,ly/2};
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
//   double L = 1.1;
//   double H = 1.1;
//   int n_e = 4;
//   double CF = 0.3;
//   double CR = 0.7;
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
//   extra = ly/2;
//   Grid* grid = Grid_new(0-extra, L+extra, -extra, H+extra, kh);
//
//   // ------------------------------------------------------------------
//   // ------------------------ SET Animation ---------------------------
//   // ------------------------------------------------------------------
//   Animation* animation = Animation_new(n_p, timeout, grid, Rp, domain);
//
//   // ------------------------------------------------------------------
//   // ------------------------ Start integration -----------------------
//   // ------------------------------------------------------------------
//   double t = 0;
//   double tEnd = 5;
//   double dt = 0.05;
//   int iter_max = (int) (tEnd-t)/dt;
//   int output = 1;
//   printf("iter max = %d\n",iter_max);
//   // // Temporal loop
//   Kernel kernel = Lucy;
//   int i = 0;
//   while (t < tEnd){
//     printf("-----------\t t/tEnd : %.3f/%.1f\t-----------\n", t,tEnd);
//     if (i%output == 0)
//       show(particles, animation, i, false, true);
//     update_cells(grid, particles, n_p);
//     update_neighbors(grid, particles, n_p, i);
//     // update_pressure(particles, n_p, rho, g, l); // Pressure dyn + P hydro
//     time_integration(particles, n_p, kernel, dt, edges);
//     printf("Time integration completed\n");
//
//     i++;
//     t += dt;
//   }
//   show(particles,animation, iter_max, false, false);
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
