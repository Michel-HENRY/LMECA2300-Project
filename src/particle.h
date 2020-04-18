#ifndef PARTICLE_H
#define PARTICLE_H

#include "vector.h"
#include "list.h"
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <math.h>



typedef struct Particle Particle;
typedef struct Parameters Parameters;
typedef struct Fields Fields;
typedef struct Cell Cell;
typedef struct Grid Grid;

struct Particle{
  Parameters* param;
  Fields* fields;

  Cell* cell;
  List* neighbors;
  List* potential_neighbors; // If Verlet
};

struct Parameters{
  double rho;
  double mass;
  double dynamic_viscosity;
  double h;
  double Rp;
};

struct Fields{
  Vector* x;
  Vector* u;
  Vector* f;
  double P;
};

struct Cell {
	int i,j;
	List* neighboring_cells;
	List* particles;
	bool visited; // for improved algorithm
};
struct Grid {
	int nCellx;
	int nCelly;
	Cell*** cells;	// 2D array of cell pointers
	double h; // side length of a cell
	double left;	// x-coordinate of left side
	double right;	// x-coordinate of right side
	double top;		// y-coordinate of top side
	double bottom;	// y-coordinate of bottom side
};
void Particle_validation();
// -------------------------------------------------------------------
// --------------------------- Grid + Cell ---------------------------
// -------------------------------------------------------------------
Grid* Grid_new(double left, double right, double bottom, double top, double h);
void Grid_free(Grid* grid);
void Cell_free(Cell* cell);
void reset_grid(Grid* grid);
Cell* localize_particle(Grid *grid, Particle *p);
void update_cells(Grid* grid, Particle** particles, int n_p);
void update_cells(Grid* grid, Particle** particles, int n_p);
void add_neighbors_from_cell(Particle* p, Cell* cell , double h);
void add_neighbors_from_cells(Grid* grid, Particle* p);
void update_from_potential_neighbors(Particle** particles, int n_p, double h);
void update_neighbors(Grid* grid, Particle** particles, int n_p, int iter);

// -------------------------------------------------------------------
// --------------------------- Parameters + fields--------------------
// -------------------------------------------------------------------
Parameters* Parameters_new(double rho, double mass, double dynamic_viscosity, double h, double Rp);
void Parameters_free(Parameters* param);
Fields* Fields_new(Vector* x, Vector* u, Vector* f, double P);
void Fields_free(Fields* fields);

// -------------------------------------------------------------------
// ------------------------------ Particle ---------------------------
// -------------------------------------------------------------------
Particle* Particle_new(Parameters* param, Fields* fields);
void Particle_free(Particle* particle);
void Particles_free(Particle** particles, int n_p);
void reset_particles(Particle** particles, int N, int iter);

#endif
