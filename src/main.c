#include "vector.h"
#include "kernel.h"
#include "particle.h"
#include "animation.h"
#include "time_integration.h"
#include "boundary.h"

#include <time.h> // To generate random number
#include <stdlib.h>
#include <stdio.h>

Kernel kernel = Cubic;

#define M_PI 3.14159265358979323846
int ListToTxt (double* X, int n){
    FILE* fichier=NULL;
    fichier = fopen("results.csv","a");
    if (fichier != NULL)
        {
            for (int i=0;i<n;i++){
                fprintf(fichier,"%f,",X[i]);
            }
            fclose(fichier);
            return 1;
        }
    else{
      printf("Impossible d'ouvrir le fichier\n");
      return 0;
    }
}

void ParticleToTxt(Particle** p, int n_p, int iter){
  FILE* fichier = NULL;
  fichier = fopen("results.csv","w");
  if(fichier != NULL){
    for(int i = 0; i < n_p ; i++){
      fprintf(fichier, "%d,\t%f,\t %f,\t%f,\t%f\n", iter,p[i]->fields->x->X[0],p[i]->fields->x->X[1] ,p[i]->fields->u->X[0], p[i]->fields->u->X[1]);
    }
    fclose(fichier);
    printf("Ecriture réussie\n");
  }
  else{
    printf("Impoosible d'ouvrir le fichier\n");
  }
}
// Validation function

static int moving_circle();

// Useful function to create the usual sructures
Particle** fluidProblem(Parameters* param, int n_p_dim_x, int n_p_dim_y, double g, double rho, double P, bool isUniform);
Edges* get_box(double L, double H, int n_e , double CF, double CR, double domain[4]);

int main(){
  moving_circle();
}
Particle** fluidProblem(Parameters* param, int n_p_dim_x, int n_p_dim_y, double g, double rho, double P, bool isUniform){
  int n_p = n_p_dim_x*n_p_dim_y;
  Particle** particles = (Particle**) malloc(n_p*sizeof(Particle*));
  for(int j = 0; j< n_p_dim_y; j++){
      for(int i = 0; i < n_p_dim_x; i++){
        int index = j*n_p_dim_x + i;

        Vector* x = Vector_new(2);
        Vector* u = Vector_new(2);
        Vector* f = Vector_new(2);

        // f->X[1] = -g*rho;

        double Rp = param->Rp;
        if (isUniform){
          x->X[0] = Rp*(2*i + 1);      x->X[1] = Rp*(2*j + 1);
        }
        else{
          double l = Rp*n_p_dim_x;
          x->X[0] = (float)rand()/(float)(RAND_MAX/l);
          x->X[1] = (float)rand()/(float)(RAND_MAX/l);
          u->X[0] = 10*(float)rand()/(float)(RAND_MAX/l);
          u->X[1 ]= 10*(float)rand()/(float)(RAND_MAX/l);
          //
          //     u->X[0] = ux;
          //     u->X[1] = uy;
        }
        Fields* fields = Fields_new(x,u,f,P,rho);
        particles[index] = Particle_new(param, fields);
      }
    }
    return particles;
}
Edges* get_box(double L, double H, int n_e , double CF, double CR, double domain[4]){
  // Create the vertices
  Vector** vertices = (Vector**) malloc(n_e*sizeof(vertices));
  for(int i = 0; i < n_e; i++){
    vertices[i] = Vector_new(2);
  }
  vertices[0]->X[0] = 0;                vertices[2]->X[0] = L;
  vertices[0]->X[1] = 0;                vertices[2]->X[1] = H;
  vertices[1]->X[0] = L;                vertices[3]->X[0] = 0;
  vertices[1]->X[1] = 0;                vertices[3]->X[1] = H;

  // Create the edges structure
  Vector** edge = (Vector**) malloc(n_e*2*sizeof(vertices));
  for(int i = 0;i < n_e; i++){
    edge[2*i] = copy(vertices[i]);
    edge[2*i+1] = copy(vertices[(i+1)%n_e]);
  }
  for(int i = 0; i < n_e ; i++){
    Vector_free(vertices[i]);
  }
  free(vertices);

  Edges* edges = Edges_new(n_e, edge, CR,CF);
  domain[0] = 0;
  domain[1] = L;
  domain[2] = 0;
  domain[3] = H;
  return edges;
}

static int moving_circle(){
  double R = 1;                         // rayon du cercle de particules
  int npR = 21;                          // Nombre de particule sur un rayon
  int np1 = 6;                          // Nombre de particule sur la 1ere circonférence
  int np = 1;                           // Nombre de particules total
  for(int i = 0; i < npR ; i++){
    np += np1*i;
  }
  printf("Nombre total de particules = %d\n",np);

  // Parameters
  double rho0 = 1000;                           // Densité initiale
  double dynamic_viscosity = 0;                 // Viscosité dynamique
  double h = R/npR;                             // Distance entre 2 particules
  double kh = 4*h;                              // Rayon du compact pour l'approximation
  double mass = M_PI*R*R*rho0/(double)np;       // Masse d'une particule, constant
  double Rp = h/2;                              // Rayon d'une particule
  double eta = 0.50;                            // XSPH parameter from 0 to 1
  double B = 1e5;
  double gamma = 7;
  // Set the parameters structures
  Parameters* param = Parameters_new(mass,kh,Rp);

  // Domain
  int ne = 100;
  double CR = 0;
  double CF = 0;
  Edges* boundary = EdgesCircle(ne,R,CR,CF);

  // Particles
  Particle** p = (Particle**)malloc(np*sizeof(Particle*));
  Fields* fields0 = Fields_new(Vector_new(2), Vector_new(2), Vector_new(2), 0, 0);
  p[0] = Particle_new(param,fields0);
  int index = 1;
  for(int i = 1; i < npR; i++){
    double Ri = 2*Rp*i;
    double ni = i*np1;
    for(int j = 0; j < ni; j++){
      double t = 2*M_PI*j/ni;
      // printf("t = %f\n",t);
      Vector* x = Vector_new(2);
      x->X[0] = Ri*cos(t);    x->X[1] = Ri*sin(t);

      Vector* u = Vector_new(2);
      double P = 0;
      double rho = 0;

      Vector* f = Vector_new(2);

      Fields* fields = Fields_new(x,u,f,P,rho);
      p[index] = Particle_new(param,fields);
      index ++;
      // printf("index = %d\n", index);
    }
  }
  // Grid
  Grid* grid = Grid_new(-R,R,-R,R,kh);

  // Animation
  double timeout = 0.001;
  Animation* animation = Animation_new(np,timeout,grid,Rp, boundary->edge, ne);

  // Initialisation
  update_cells(grid,p,np);
  update_neighbors(grid,p,np,0);
  density(p,np,kernel);
  CSPM_density(p,np,kernel);

  double t = 0;
  double tEnd = 0.5;
  double dt = 1e-4;
  int i = 0;
  while(t < tEnd){
    printf("-----------\t i = %d \t t/tEnd : %.3f/%f\t-----------\n",i,t,tEnd);
    update_cells(grid,p,np);
    update_neighbors(grid,p,np,0);
    CSPM_density(p,np,kernel);
    state_pressure(p,np,rho0,B,gamma);
    time_integration_XSPH(p,np,kernel,dt,boundary,eta);
    ParticleToTxt(p,np,i);
    show(p,animation,i,false,false);
    i++;
    t += dt;
  }
  show(p, animation, 0, true,true);

  Particles_free(p,np);
  Edges_free(boundary);
  Grid_free(grid);
  Animation_free(animation);
  return EXIT_SUCCESS;

}
