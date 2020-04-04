#ifndef PARTICLE_H
#define PARTICLE_H

#include "vector.h"
#include "list.h"


typedef struct Particle Particle;
typedef struct Parameters Parameters;
typedef struct Physical_fields Physical_fields;
typedef struct Cell Cell;
typedef struct Grid Grid;

struct Particle{
  Parameters* param;
  Physical_fields* fields;
  List* neighbors;

  Cell* cell;
  List* potential_neighborhood; // If Verlet
};

struct Parameters{
  double rho;
  double mass;
  double dynamic_viscosity;
};

struct Physical_fields{
  Vector* x;
  Vector* v;
  double P;
  Vector* f;
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
#endif
