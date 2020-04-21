#include "particle.h"


void Particle_validation(){
  //particle domain
  double l = 1.0;
  // Grid definition
  double L = 1.5;
  double H = 1.5;

  // Parameters
  double rho = 1e3;
  double dynamic_viscosity = 1e-3;
  int n_p_dim = 10;
  int n_p = n_p_dim*n_p_dim;
  double h = l/(n_p_dim - 1); // step between neighboring particles
	double mass = rho * h*h;

  // Create Particles
  Particle** particles = (Particle**) malloc(n_p*sizeof(Particle*));

  for(int i = 0; i < n_p_dim; i++){
    for(int j = 0; j < n_p_dim; j++){
      int index = i*n_p_dim + j;
      Parameters* param = Parameters_new(rho, mass, dynamic_viscosity, h, h, 0,0,0);
      Vector* x = Vector_new(2);
      Vector* u = Vector_new(2);
      Vector* f = Vector_new(2);
      double P = 0;
      double pos[2] = {i*h,j*h};

      Vector_initialise(x,pos);
      Fields* fields = Fields_new(x,u,f,P);
      particles[index] = Particle_new(param, fields);
    }
  }
  // // Create GRID
  Grid* grid = Grid_new(0, L, 0, H, h);
  // // Link grid and particles
  update_cells(grid, particles, n_p);
	update_neighbors(grid, particles, n_p, 0);

  // Free memory
  Particles_free(particles, n_p);
  printf("END FREE PARTICLES\n");
  Grid_free(grid);
  printf("END FREE GRID\n");
}

// -------------------------------------------------------------------
// -------------------------------------------------------------------
// -------------------------------------------------------------------
// --------------------------- Grid + Cell ---------------------------
// -------------------------------------------------------------------
// -------------------------------------------------------------------
// -------------------------------------------------------------------
Grid* Grid_new(double left, double right, double bottom, double top, double h) {
	// Build the grid
	int nCellx = ceil((right-left) / h);
	int nCelly = ceil((top-bottom) / h);
	printf("Grid size: (%d,%d)\n", nCellx, nCelly);
  // Can be replace by a CELL_NEW FUNCION.
	Cell*** cells = (Cell***) malloc(nCellx * sizeof(Cell**));

	for(int i = 0; i < nCellx; i++) {
		cells[i] = (Cell**) malloc(nCelly * sizeof(Cell*));
		for(int j = 0; j < nCelly; j++)
			cells[i][j] = (Cell*) malloc(sizeof(Cell));
	}

	// Build links between cells
	for(int i = 0; i < nCellx; i++) {
		for(int j = 0; j < nCelly; j++) {
			cells[i][j]->i = i, cells[i][j]->j = j;
			cells[i][j]->particles = List_new();
			cells[i][j]->visited = false;
			// Assign neighbor cells
			cells[i][j]->neighboring_cells = List_new();
			for (int di = -1; di <= 1; di++) for (int dj = -1; dj <= 1; dj++)
				if ((di || dj) && i+di >= 0 && i+di < nCellx && j+dj >= 0 && j+dj < nCelly)
					List_append(cells[i][j]->neighboring_cells, cells[i+di][j+dj]);
		}
	}
	Grid* grid = (Grid*)malloc(sizeof(Grid));
	grid->left = left,	grid->right = left + nCellx*h; // not very elegant but ok...
	grid->bottom = bottom,	grid->top = bottom + nCelly*h;
	grid->nCellx = nCellx;
	grid->nCelly = nCelly;
	grid->h = h;
	grid->cells = cells;
	return grid;
}
void Grid_free(Grid* grid) {
	for(int i = 0; i < grid->nCellx; i++) {
		for(int j = 0; j < grid->nCelly; j++)
			Cell_free(grid->cells[i][j]);
		free(grid->cells[i]);
	}
	free(grid->cells);
	free(grid);
}
void Cell_free(Cell* cell) {
	List_free(cell->neighboring_cells, NULL);
	List_free(cell->particles, NULL);
	free(cell);
}
void reset_grid(Grid* grid){
	for(int i = 0; i < grid->nCellx; i++) {
		for(int j = 0; j < grid->nCelly; j++) {
			Cell* cell = grid->cells[i][j];
			List_free(cell->particles, NULL);
			cell->particles = List_new();
			cell->visited = false;
		}
	}
}
// ------------------------- localize particle -------------------------
Cell* localize_particle(Grid *grid, Particle *p){
	int i = floor((p->fields->x->X[0] - grid->left) / grid->h);
	int j = floor((p->fields->x->X[1] - grid->bottom) / grid->h);
	if(i < 0 || i >= grid->nCellx || j < 0 || j >= grid->nCelly) {
		fprintf(stderr, "ERROR: Particle is outside the grid :(\n");
		exit(0);
	}
	return grid->cells[i][j];
}
// --------------------------- update links  --------------------------
void update_cells(Grid* grid, Particle** particles, int n_p){
	// Clean the grid before update
	reset_grid(grid);
	for(int i = 0; i < n_p; i++){
		Cell* cell = localize_particle(grid, particles[i]);
		particles[i]->cell = cell; // link cell to particle
		List_append(cell->particles, particles[i]); // link particle to cell
	}
}
// ---------------------- update neighbourhood  ----------------------
void add_neighbors_from_cell(Particle* p, Cell* cell , double h){
	// Iterate over particles in cell
	ListNode *node = cell->particles->head;
	while (node != NULL) {
		Particle* q = (Particle*)node->v;
    double dist_pq = dist(p->fields->x,q->fields->x);
		if(dist_pq <= h && !equal(p->fields->x,q->fields->x)){
      List_append(p->neighbors, q);
    }
		node = node->next;
	}
}
void add_neighbors_from_cells(Grid* grid, Particle* p){
	add_neighbors_from_cell(p, p->cell, grid->h);
	ListNode *node = p->cell->neighboring_cells->head;
	while (node != NULL) {
		Cell* cell = (Cell*) node->v;
		add_neighbors_from_cell(p, cell, grid->h);
		node = node->next;
	}
}
// Among potential neighbors, filter the valid ones
void update_from_potential_neighbors(Particle** particles, int n_p, double h){
	for (int i = 0; i < n_p; i++) {
		Particle* p = particles[i];
		ListNode *node = p->potential_neighbors->head;
		while (node != NULL) {
			Particle* q = (Particle*)node->v;
			if(dist(p->fields->x,q->fields->x) <= h){
        List_append(p->neighbors, q);
      }
			node = node->next;
		}
	}
}
void update_neighbors(Grid* grid, Particle** particles, int n_p, int iter){
	// Clean the particles before update
	reset_particles(particles, n_p, iter);
	for (int i = 0 ; i < n_p; i++)
		add_neighbors_from_cells(grid, particles[i]);
}

// -------------------------------------------------------------------
// -------------------------------------------------------------------
// -------------------------------------------------------------------
// --------------------------- Parameters + fields--------------------
// -------------------------------------------------------------------
// -------------------------------------------------------------------
// -------------------------------------------------------------------

Parameters* Parameters_new(double rho, double mass, double dynamic_viscosity, double h, double Rp, double tension, double treshold, double P0){
  Parameters* param = (Parameters*) malloc(sizeof(Parameters));
  param->rho = rho;
  param->mass = mass;
  param->dynamic_viscosity = dynamic_viscosity;
  param->h = h;
  param->Rp = Rp;
  param->tension = tension;
  param->treshold = treshold;
  param->P0 = P0;
  return param;
}
void Parameters_free(Parameters* param){
  free(param);
}

Fields* Fields_new(Vector* x, Vector* u, Vector* f, double P){
  Fields* fields = (Fields*) malloc(sizeof(Fields));
  fields->x = x;
  fields->u = u;
  fields->f = f;
  fields->P = P;
  fields->Cs = 1;

  return fields;
}
void Fields_free(Fields* fields){
  Vector_free(fields->x);
  Vector_free(fields->u);
  Vector_free(fields->f);

  free(fields);
}

// -------------------------------------------------------------------
// -------------------------------------------------------------------
// -------------------------------------------------------------------
// ------------------------------ Particle ---------------------------
// -------------------------------------------------------------------
// -------------------------------------------------------------------
// -------------------------------------------------------------------

Particle* Particle_new(Parameters* param, Fields* fields){
  Particle* particle = (Particle*) malloc(sizeof(Particle));
  particle->param = param;
  particle->fields = fields;
  particle->cell = NULL;
  particle->neighbors = List_new();
  particle->potential_neighbors = List_new();

  return particle;
}
void Particle_free(Particle* particle){
  // Cell are freed in the GRID
  Parameters_free(particle->param);
  Fields_free(particle->fields);
  List_free(particle->neighbors,NULL);
  List_free(particle->potential_neighbors,NULL);

  free(particle);
}

void Particles_free(Particle** particles, int n_p){
  for(int i = 0; i < n_p; i++){
    Particle_free(particles[i]);
  }
  free(particles);
}
// Empty neighbors of each particle
void reset_particles(Particle** particles, int N, int iter){
	for (int i = 0; i < N; i++) {
		Particle* p = particles[i];
		List_free(p->neighbors, NULL);
		p->neighbors = List_new();
	}
}
