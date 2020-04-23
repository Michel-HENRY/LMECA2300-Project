#ifndef BOUNDARY_H
#define BOUNDARY_H

#include "vector.h"
#include "particle.h"
#include <math.h>

// TODO : Utiliser des listes pour pop les frontières qui ont déjà été validées.
typedef struct Edges Edges;

struct Edges {
  int n_e;        // Taille de edges
  Vector** edge;  // list contenantles paires de noeuds dans le sens horaire [a,b,b,c,c,a]
  Vector** n;      // Normale nortante
  double CR;
  double CF;
};

Edges* Edges_new(int n_e, Vector** edge, double CR, double CF);
void Edges_free(Edges* edges);
void set_normal(Edges* edges);

bool isInside(Vector* C1, Edges* edges, double Rp);
bool CenterIsInside(Vector* C1, Vector* e0, Vector* e1);
double* distEdge(Vector* C1, Edges* edges);
int indexCPlane(double* dist_edges, int n_e);
void update_mass_center(Vector* C1, double Rp, Edges* edges, double d, int index);
void update_velocity(Particle* p, Edges* edges, int index);

void reflective_boundary(Particle** p, int n_p, Edges* edges);
#endif
