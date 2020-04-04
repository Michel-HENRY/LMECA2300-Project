#include "particle.h"

// Private functions
bool check_distance(xy *p, xy *q, double r);

Grid* Grid_new(double x1, double x2, double y1, double y2, double h) {
	// Build the grid
	int nCellx = ceil((x2-x1) / h);
	int nCelly = ceil((y2-y1) / h);
	printf("Grid size: (%d,%d)\n", nCellx, nCelly);
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
	grid->left = x1,	grid->right = x1 + nCellx*h; // not very elegant but ok...
	grid->bottom = y1,	grid->top = y1 + nCelly*h;
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

Particle* Particle_new(int index, double m, xy* pos, xy* v, double rho_0, double mu, double c_0, double gamma, double sigma) {
	Particle *particle = malloc(sizeof(Particle));
	particle->index = index;
	particle->m = m;
	particle->pos = pos;
	particle->rho = rho_0;
	particle->v = v;
	particle->P = 0.0; // assuming that the fluid is at rest (P is the dynamic pressure and not the absolute one!)

	particle->normal = xy_new(0.0,0.0);
	particle->XSPH_correction = xy_new(0.0,0.0);
	particle->on_free_surface = false;

	particle->param = malloc(sizeof(Physical_parameters));
	particle->param->rho_0 = rho_0;
	particle->param->dynamic_viscosity = mu;
	particle->param->gamma = gamma;
	particle->param->sound_speed = c_0;
	particle->param->sigma = sigma;

	particle->cell = NULL;
	particle->neighborhood = List_new();
	particle->potential_neighborhood = List_new();
	return particle;
}

void Particle_free(Particle* particle) {
	free(particle->pos);
	free(particle->v);
	free(particle->param);
	free(particle->XSPH_correction);
	List_free(particle->neighborhood, NULL);
	List_free(particle->potential_neighborhood, NULL);
	free(particle);
}

void free_particles(Particle** particles, int N) {
	for (int i = 0; i < N; i++)
		Particle_free(particles[i]);
	free(particles);
}

Particle_derivatives* Particle_derivatives_new(int index) {
	Particle_derivatives *particle_derivatives = malloc(sizeof(Particle_derivatives));
	particle_derivatives->index = index;
	particle_derivatives->div_v = 0;
	particle_derivatives->grad_P = xy_new(0,0);
	particle_derivatives->lapl_v = xy_new(0,0);
	particle_derivatives->grad_Cs = xy_new(0,0);
	particle_derivatives->lapl_Cs = 0;
	return particle_derivatives;
}

void Particle_derivatives_free(Particle_derivatives* particle_derivatives) {
	free(particle_derivatives->grad_P);
	free(particle_derivatives->lapl_v);
	free(particle_derivatives->grad_Cs);
	free(particle_derivatives);
}
void Particle_derivatives_reset(Particle_derivatives *particle_derivatives) {
	particle_derivatives->div_v = 0;
	xy_reset(particle_derivatives->grad_P);
	xy_reset(particle_derivatives->lapl_v);
	xy_reset(particle_derivatives->grad_Cs);
	particle_derivatives->lapl_Cs = 0;
}

void free_particles_derivatives(Particle_derivatives** particles_derivatives, int N) {
	for (int i = 0;i < N;i++)
		Particle_derivatives_free(particles_derivatives[i]);
	free(particles_derivatives);
}

double Particle_get_P(Particle *particle) {	return particle->P; }
xy * Particle_get_v(Particle *particle) { return particle->v; }
xy * Particle_get_pos(Particle *particle) { return particle->pos; }
double Particle_get_v_x(Particle *particle) { return particle->v->x; }
double Particle_get_v_y(Particle *particle) { return particle->v->y; }
double Particle_get_Cs(Particle *particle) { return particle->Cs; }
xy * Particle_get_normal(Particle *particle) { return particle->normal; }


Verlet* Verlet_new(double kh, double L, int T) {
	Verlet *v = malloc(sizeof *v);
	v->kh = kh;
	v->L = L;
	v->T = T;
	return v;
}
///////////////////////////update cells///////////////////////////
Cell* localize_particle(Grid *grid, Particle *p) {
	int i = floor((p->pos->x - grid->left) / grid->h);
	int j = floor((p->pos->y - grid->bottom) / grid->h);
	if(i < 0 || i >= grid->nCellx || j < 0 || j >= grid->nCelly) {
		fprintf(stderr, "ERROR: Particle is outside the grid :(\n");
		exit(0);
	}
	return grid->cells[i][j];
}

// Update links between cells and particles
void update_cells(Grid* grid, Particle** particles, int N) {
	// Clean the grid before update
	reset_grid(grid);
	for(int i = 0; i < N; i++){
		Cell* cell = localize_particle(grid, particles[i]);
		particles[i]->cell = cell; // link cell to particle
		List_append(cell->particles, particles[i]); // link particle to cell
	}
}

/////////////////////////update neighborhood////////////////////////
// Add to the neighbors of particle p all particles q in cell s.t. |p-q| <= r
void add_neighbors_from_cell(Particle* p, Cell* cell , double r) {
	// Iterate over particles in cell
	ListNode *node = cell->particles->head;
	while (node != NULL) {
		Particle* q = (Particle*)node->v;
		// if((p->index != q->index) && check_distance(p->pos, q->pos, r))
		if(check_distance(p->pos, q->pos, r))
				List_append(p->neighborhood, q);
		node = node->next;
	}
}

// Add to particle p all its neighbors (from 9 cells)
void add_neighbors_from_cells(Grid* grid, Particle* p) {
	add_neighbors_from_cell(p, p->cell, grid->h);
	ListNode *node = p->cell->neighboring_cells->head;
	while (node != NULL) {
		Cell* cell = (Cell*) node->v;
		add_neighbors_from_cell(p, cell, grid->h);
		node = node->next;
	}
}

// Among potential neighbors, filter the valid ones
void update_from_potential_neighbors(Particle** particles, int N, double r) {
	for (int i = 0; i < N; i++) {
		Particle* p = particles[i];
		ListNode *node = p->potential_neighborhood->head;
		while (node != NULL) {
			Particle* q = (Particle*)node->v;
			if(check_distance(p->pos, q->pos, r))
				List_append(p->neighborhood, q);
			node = node->next;
		}
	}
}


void update_neighborhoods(Grid* grid, Particle** particles, int N, int iter, Verlet* verlet) {
	// Clean the particles before update
	reset_particles(particles, N, iter, verlet);
	if(verlet==NULL) {
		//update_neighborhoods_improved(grid);
		for (int i = 0 ; i < N; i++)
			add_neighbors_from_cells(grid, particles[i]);
	} else {
		if (iter%verlet->T == 0) {
			// update_neighborhoods_improved(grid);
			for (int i = 0; i < N; i++) {
				add_neighbors_from_cells(grid, particles[i]); // TODO: shouldn't verlet->L be used here?
				// Swap nbh and potential nbh
				// TODO: this is very ugly
				List* l = particles[i]->potential_neighborhood;
				particles[i]->potential_neighborhood = particles[i]->neighborhood;
				particles[i]->neighborhood = l;
			}
		}
		update_from_potential_neighbors(particles, N, verlet->kh);
	}
}

//////////////////update neighborhood-IMPROVED algorithm///////////////

// Compute all neighbors inside a cell
void update_neighbors_1(Cell* cell,double r) {
	ListNode *node = cell->particles->head;
	while (node != NULL) {
		Particle* p = (Particle*)node->v;

		ListNode *node2 = node->next;
		while (node2 != NULL) {
			Particle* q = (Particle*)node2->v;
			if(check_distance(p->pos, q->pos, r)) {
				List_append(p->neighborhood, q);
				List_append(q->neighborhood, p);
			}
			node2 = node2->next;
		}
		node = node->next;
	}
}

//compute all neighbors between 2 cells
void neighbors_in_2_cells(Cell* cell1, Cell* cell2, double r) {
	ListNode *node1 = cell1->particles->head;
	while (node1 != NULL) {
		Particle* p1 = (Particle*)node1->v;

		ListNode *node2 = cell2->particles->head;
		while (node2 != NULL) {
			Particle* p2 = (Particle*)node2->v;
			if(check_distance(p1->pos, p2->pos, r)) {
				List_append(p1->neighborhood, p2);
				List_append(p2->neighborhood, p1);
			}
			node2 = node2->next;
		}

		node1 = node1->next;
	}
}

// Compute all neighbors from this cell's particles
//// Compute the neighbors of particules of cell with the neighboring cells
void update_neighbors_2(Cell* cell, double r) {
	ListNode *node = cell->neighboring_cells->head;
	while (node != NULL) {
		Cell* neighborCell = (Cell*)node->v;
		if (!neighborCell->visited) { // if the neighboring cell was not already visited
			neighbors_in_2_cells(cell, neighborCell, r);
		}
		node = node->next;
	}
}

void update_neighborhoods_improved(Grid* grid) {
	Cell*** cells = grid->cells;
	for(int i = 0 ;i < grid->nCellx; i++) {
		for(int j = 0; j < grid->nCelly; j++) {
			update_neighbors_1(cells[i][j], grid->h);
			update_neighbors_2(cells[i][j], grid->h);
			cells[i][j]->visited = true;
		}
	}
}


////////////////////////////////////////////////////

// Empty the list of particles inside each Cell
void reset_grid(Grid* grid) {
	for(int i = 0; i < grid->nCellx; i++) {
		for(int j = 0; j < grid->nCelly; j++) {
			Cell* cell = grid->cells[i][j];
			List_free(cell->particles, NULL);

			cell->particles = List_new();
			cell->visited = false;
		}
	}
}

// Empty neighborhood of each particle
void reset_particles(Particle** particles, int N, int iter, Verlet* verlet) {
	for (int i = 0; i < N; i++) {
		Particle* p = particles[i];
		List_free(p->neighborhood, NULL);
		p->neighborhood = List_new();
		// If in Verlet mode, empty potential nbh
		if(verlet != NULL && iter%verlet->T == 0) {
			List_free(p->potential_neighborhood, NULL);
			p->potential_neighborhood = List_new();
		}
	}
}

// Check if |p-q| <= r
bool check_distance(xy *p, xy *q, double r) {
	return squared(p->x - q->x) + squared(p->y - q->y) <= squared(r);
}

// Generate N particles randomly located on [-L,L] x [-L,L]
// Velocity, rho and e are zero.
Particle** build_particles(int N, double L) {
	Particle** particles = (Particle**)malloc(N * sizeof(Particle*));
	for (int i = 0; i < N; i++) {
		double x = rand_interval(-L, L);
		double y = rand_interval(-L, L);
		xy* pos = xy_new(x, y);
		xy* vel = xy_new(0, 0);
		particles[i] = Particle_new(i, 0, pos, vel, 0, 0, 0, 0, 0);
	}
	return particles;
}
