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
  double mass;
  double dynamic_viscosity;
  double h;
  double Rp;
  double tension;
  double treshold;
  double P0;
  double g;
  double a,b;
};

struct Fields{
  Vector* x;
  Vector* u;
  Vector* f;
  double P;
  double Cs;
  double rho;
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
void display_neighbors(Particle** particles, int n_p);
// -------------------------------------------------------------------
// --------------------------- Grid + Cell ---------------------------
// -------------------------------------------------------------------
Grid* Grid_new(double left, double right, double bottom, double top, double h);
void Grid_free(Grid* grid);
static void Cell_free(Cell* cell);
static void reset_grid(Grid* grid);
static Cell* localize_particle(Grid *grid, Particle *p);
void update_cells(Grid* grid, Particle** particles, int n_p);
static void add_neighbors_from_cell(Particle* p, Cell* cell , double h);
static void add_neighbors_from_cells(Grid* grid, Particle* p);
static void update_from_potential_neighbors(Particle** particles, int n_p, double h);
void update_neighbors(Grid* grid, Particle** particles, int n_p, int iter);

// -------------------------------------------------------------------
// --------------------------- Parameters + fields--------------------
// -------------------------------------------------------------------
Parameters* Parameters_new(double mass, double dynamic_viscosity, double h, double Rp, double tension, double treshold, double P0, double g);
void Parameters_free(Parameters* param);
Fields* Fields_new(Vector* x, Vector* u, Vector* f, double P, double rho);
void Fields_free(Fields* fields);
void set_artificialViscosity(Parameters* param, double a, double b);

// -------------------------------------------------------------------
// ------------------------------ Particle ---------------------------
// -------------------------------------------------------------------
Particle* Particle_new(Parameters* param, Fields* fields);
void Particle_free(Particle* particle);
void Particles_free(Particle** particles, int n_p);
static void reset_particles(Particle** particles, int N, int iter);
int get_n_neighbors(Particle* p);
void set_density(Particle** particles, int nx, int ny, double rho0, double B, double gamma);
#endif
