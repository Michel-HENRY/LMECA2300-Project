#include "SPH.h"

Setup* Setup_new(int iter, double timestep,double kh,Verlet* verlet,Kernel kernel, Free_surface_detection free_surface_detection, double interface_threshold,double XSPH_epsilon, bool gravity) {
	Setup* setup = (Setup*)malloc(sizeof(Setup));
	setup->itermax = iter;
	setup->timestep = timestep;
	setup->kh = kh;
	setup->verlet = verlet;
	setup->kernel = kernel;
	setup->free_surface_detection = free_surface_detection;
	setup->interface_threshold = interface_threshold;
	setup->XSPH_epsilon = XSPH_epsilon;
	setup->gravity = gravity;
	return setup;
}

void Setup_free(Setup* setup) {
	if(setup->verlet != NULL)
		free(setup->verlet);
	free(setup);
}

Residual* Residual_new(){
	Residual* residual = (Residual*)malloc(sizeof(Residual));
	residual->mass_eq = 0;
	residual->momentum_x_eq = 0;
	residual->momentum_y_eq = 0;
	return residual;
}

void free_Residuals(Residual** residuals, int N) {
	for (int i = 0;i < N;i++)
		free(residuals[i]);
	free(residuals);
}
double max_velocity(Particle** p, int n_p){
	double max = 0;
	for(int i = 0 ; i < n_p; i++){
		double v = fabs(p[i]->v->x);
		if(max < v){
			max = v;
		}
	}
	return max;
}
void simulate(Grid* grid, Particle** particles, Particle_derivatives** particles_derivatives, Residual** residuals, int n_p, update_positions update_positions, Setup* setup, Animation* animation) {
	double current_time = 0.0;
	int ii = 5;
	printf("%d\n", setup->itermax);
	for (int iter = 0; iter < setup->itermax; iter++) {
		printf("----------------------------------------------------- \n");
		printf("iter %d / %d @ t = %lf \n", iter, setup->itermax, current_time);
		update_cells(grid, particles, n_p);
		update_neighborhoods(grid, particles, n_p, iter, setup->verlet);
		if (animation != NULL)
			display_particles(particles, animation, false, iter);
		update_positions(grid, particles, particles_derivatives, residuals, n_p, setup);
		printf("velocity_max = %f\n", max_velocity(particles,n_p));
		if (iter%ii == 0){
			// density_correction_MLS(particles, n_p, setup->kh, setup->kernel);
		}
		get_M0(particles,n_p,setup->kh,setup->kernel);
		get_M1(particles,n_p,setup->kh,setup->kernel);
		current_time += setup->timestep;
	}
	update_cells(grid, particles, n_p);
	update_neighborhoods(grid, particles, n_p, 0, setup->verlet);
	if (animation != NULL)
		display_particles(particles, animation, true, setup->itermax);
}
void simulate_boundary(Grid* grid, Particle** particles, Particle_derivatives** particles_derivatives, Residual** residuals, int n_p, update_positions update_positions, Setup* setup, Animation* animation, Boundary* boundary){
	double current_time = 0.0;
	double Rp = 0.001;
	int ii = 5;
	double bounds[4] = {boundary->xleft,boundary->xright, boundary->ybottom, boundary->ytop};
	printf("%d\n", setup->itermax);
	for (int iter = 0; iter < setup->itermax; iter++) {
		printf("----------------------------------------------------- \n");
		printf("iter %d / %d @ t = %lf \n", iter, setup->itermax, current_time);
		update_cells(grid, particles, n_p);
		update_neighborhoods(grid, particles, n_p, iter, setup->verlet);
		if (animation != NULL)
			display_particles_boundary(particles, animation, false,iter,bounds);
		update_positions(grid, particles, particles_derivatives, residuals, n_p, setup);

		printf("velocity_max = %f\n", max_velocity(particles,n_p));
		reflective_boundary(particles,n_p,boundary,Rp);
		if (iter%ii == 0){
			// density_correction_MLS(particles, n_p, setup->kh, setup->kernel);
		}
		get_M0(particles,n_p,setup->kh,setup->kernel);
		get_M1(particles,n_p,setup->kh,setup->kernel);
		current_time += setup->timestep;
	}
	update_cells(grid, particles, n_p);
	update_neighborhoods(grid, particles, n_p, 0, setup->verlet);
	if (animation != NULL)
		display_particles(particles, animation, true,setup->itermax);
}


//move randomly each particles
void random_moves(Grid* grid, Particle** particles, Particle_derivatives** particles_derivatives, Residual** residuals, int n_p, Setup* setup) {
	double max_speed = 2;
	for (int i = 0; i < n_p; i++) {
		double angle = rand_interval(0, 2)*M_PI;
		double speed = rand_interval(0, max_speed);
		Particle *p = particles[i];
		p->v->x = speed * cos(angle);
		p->v->y = speed * sin(angle);
		p->pos->x += p->v->x * setup->timestep;
		p->pos->y += p->v->y * setup->timestep;

		double s = 2;
		//bouncing with the wall
		if (p->pos->x < grid->left)
			p->pos->x = grid->left + s;
		if (p->pos->x > grid->right)
			p->pos->x = grid->right - s;
		if (p->pos->y < grid->bottom)
			p->pos->y = grid->bottom + s;
		if (p->pos->y > grid->top)
			p->pos->y = grid->top - s;
	}
}
void update_positions_Colagrossi(Grid* grid, Particle** particles, Particle_derivatives** particles_derivatives, Residual** residuals, int n_p, Setup* setup){
	// Colagrossi method
	// Compute derivatives
	for(int i = 0; i < n_p; i++){
		Particle* pi = particles[i];
		Particle_derivatives* dpi = particles_derivatives[i];
		Kernel kernel = setup->kernel;
		double kh = setup->kh;
		dpi->div_v = compute_div(pi, Particle_get_v, kernel, kh);//Verifier
		dpi->lapl_v->x = compute_lapl(pi, Particle_get_v_x,kernel,kh);
		dpi->lapl_v->y = compute_lapl(pi, Particle_get_v_y,kernel,kh);
		compute_grad(pi, Particle_get_P, kernel, kh, particles_derivatives[i]->grad_P);
	}
	// Compute the residuals
	double g = -9.81;
	for(int i = 0; i < n_p; i ++){
		Particle* pi = particles[i];
		Residual* res = residuals[i];
		Particle_derivatives* dpi = particles_derivatives[i];
		double rho = pi->rho;
		double mu = pi->param->dynamic_viscosity;
		// Mass conservation
		res->mass_eq = rho*dpi->div_v;
		//Momentum equation (on dirait qu'il manque un terme)
		// Pressure gradient
		res->momentum_x_eq = (-1.0/rho)*dpi->grad_P->x;
		res->momentum_y_eq = (-1.0/rho)*dpi->grad_P->y;
		// Viscosity constraint
		res->momentum_x_eq += (mu/rho) * dpi->lapl_v->x;
		res->momentum_y_eq += (mu/rho) * dpi->lapl_v->y;

		if (setup->gravity == 1){
			res->momentum_y_eq -= g;
		}
	}
	// Time integration
	double dt = setup->timestep;
	for(int i = 0; i < n_p; i ++){
		Particle* pi = particles[i];
		Residual* res = residuals[i];
		// Explicit Euler
		pi->rho += dt * res->mass_eq;
		pi->v->x += dt * res->momentum_x_eq;
		pi->v->y += dt * res->momentum_y_eq;
		// Implicit Euler
		pi->pos->x += pi->v->x*dt;
		pi->pos->y += pi->v->y*dt;
	}
	// Constitutive equation for the pressure
	for(int i = 0; i < n_p; i ++){
		Particle* pi = particles[i];
		double rho = pi->rho;
		double c = pi->param->sound_speed;
		double rho_0 = pi->param->rho_0;
		double gamma = pi->param->gamma;

		double B = c*c*rho_0/gamma;

		pi->P = B*(pow(rho/rho_0,gamma)-1);
	}
}
void update_positions_seminar_5(Grid* grid, Particle** particles, Particle_derivatives** particles_derivatives, Residual** residuals, int n_p, Setup* setup) {

	// Compute Cs, the XSPH correction on the velocity, and the divergence of the positions
	for (int i = 0; i < n_p; i++) {
		compute_Cs(particles[i], setup->kernel, setup->kh);
		if (setup->XSPH_epsilon != 0.0) compute_XSPH_correction(particles[i], setup->kernel, setup->kh,setup->XSPH_epsilon);
	}

	// Compute derivatives and normal
	for (int i = 0; i < n_p; i++) {
		particles_derivatives[i]->div_v = compute_div(particles[i], Particle_get_v, setup->kernel, setup->kh);
		particles_derivatives[i]->lapl_v->x = compute_lapl(particles[i], Particle_get_v_x, setup->kernel, setup->kh);
		particles_derivatives[i]->lapl_v->y = compute_lapl(particles[i], Particle_get_v_y, setup->kernel, setup->kh);
		compute_grad(particles[i], Particle_get_P, setup->kernel, setup->kh, particles_derivatives[i]->grad_P);
		compute_grad(particles[i], Particle_get_Cs, setup->kernel, setup->kh, particles_derivatives[i]->grad_Cs);
		particles_derivatives[i]->lapl_Cs = compute_lapl(particles[i], Particle_get_Cs, setup->kernel, setup->kh);
		// assemble_residual_NS(particles[i], particles_derivatives[i], residuals[i], setup);
		compute_normal(particles[i], particles_derivatives[i]);
	}

	// Assemble residual and compute curvature
	for (int i = 0; i < n_p; i++) {
	    //particles[i]->kappa = 2.0*compute_div(particles[i], Particle_get_normal, setup->kernel, setup->kh);
	    assemble_residual_NS(particles[i], particles_derivatives[i], residuals[i], setup);
	}

	// Integrate (obtain new values, i.e. density, velocities, pressure and positions, at time t+1)
	for (int i = 0; i < n_p; i++)
		time_integrate(particles[i], residuals[i], setup->timestep);
		// time_integrate_CSPM(particles[i],particles_derivatives[i], residuals[i], setup);
}

void compute_Cs(Particle *particle, Kernel kernel, double kh) {
	particle->Cs = 0;
	Particle *pi = particle;
	ListNode *node = pi->neighborhood->head;
	while (node != NULL) {
		Particle *pj = node->v;
		particle->Cs += (pj->m / pj->rho) * eval_kernel(pi->pos, pj->pos, kh, kernel);
		node = node->next;
	}
	// printf("pos = (%lf, %lf), Cs = %lf\n", particle->pos->x, particle->pos->y, particle->Cs);
}

void compute_normal(Particle *particle, Particle_derivatives* particle_derivatives) {
	particle->normal = xy_new(0.0, 0.0);
	xy *n = particle_derivatives->grad_Cs; // surface normal inward
	double norm_n = norm(n); // norm of n
	particle->normal->x = n->x / norm_n;
	particle->normal->y = n->y / norm_n;
}

// Assemble the residual of the (incompressible) Navier-Stokes equations based on the derivatives available
void assemble_residual_NS(Particle* particle, Particle_derivatives* particle_derivatives, Residual* residual,Setup* setup) {
	double mu_i = particle->param->dynamic_viscosity;

	double rho_i = particle->rho;
	double div_vel_i = particle_derivatives->div_v;
	xy* grad_P = particle_derivatives->grad_P;
	xy* lapl_v = particle_derivatives->lapl_v;

	// Compute UNIT normal vector
	xy *n = particle_derivatives->grad_Cs; // surface normal inward
	double norm_n = norm(n);
	n->x /= norm_n, n->y /= norm_n;

	double lapl_Cs = particle_derivatives->lapl_Cs;
	// Choose between curvature estimated with Laplacian of colour field or with divergence of normal
	// 	double kappa = - lapl_Cs / norm_n; // curvature with Laplacian of colour field

	double fs_x = 0; double fs_y = 0;
	// Apply surface tension only on particles in the vicinity the interface
	bool criterion;
	// Identification based on the norm of the normal
	if (setup->free_surface_detection == CSF)
		criterion = norm_n > setup->interface_threshold;
	// Identification based on the divergence of the position vector
	else if (setup->free_surface_detection == DIVERGENCE)
		criterion = compute_div(particle, Particle_get_pos, setup->kernel, setup->kh) <= setup->interface_threshold;
	else
		criterion = false;
	if (criterion) {
		particle->on_free_surface = true;
	  double kappa = compute_curvature(particle, setup, 0.5);
		particle->P= 0;
	}
	else
		particle->on_free_surface = false;

	residual->mass_eq = -rho_i * div_vel_i;
	residual->momentum_x_eq = (-1.0/rho_i) * grad_P->x + (mu_i/rho_i) * lapl_v->x;
	residual->momentum_y_eq = (-1.0/rho_i) * grad_P->y + (mu_i/rho_i) * lapl_v->y;
	if (setup->gravity == 1){
		double g = 9.81;
		residual->momentum_y_eq -= g;
	}
}



// Time integrate the Navier-Stokes equations based on the residual already assembled
void time_integrate(Particle* particle, Residual* residual, double delta_t) {

	// Update density and velocity with an Euler explicit scheme (TODO: implement more accurate and more stable schemes)
	particle->rho += delta_t * residual->mass_eq;
	particle->v->x += delta_t * residual->momentum_x_eq;
	particle->v->y += delta_t * residual->momentum_y_eq;

	// Update position with an Euler Implicit scheme
	particle->pos->x += delta_t * particle->v->x - delta_t * particle->XSPH_correction->x;
	particle->pos->y += delta_t * particle->v->y - delta_t * particle->XSPH_correction->y;


	// Update pressure with Tait's equation of state
	double B = squared(particle->param->sound_speed) * particle->param->rho_0 / particle->param->gamma;
	// double B = 0.85*1e5;
	particle->P = B * (pow(particle->rho / particle->param->rho_0, particle->param->gamma) - 1);

}
void time_integrate_CSPM(Particle* particle, Particle_derivatives *dp, Residual* residual, Setup* setup) {

	// Update density and velocity with an Euler explicit scheme (TODO: implement more accurate and more stable schemes)
	double delta_t = setup->timestep;
	particle->rho += delta_t * residual->mass_eq;
	Corrective_Smoothed_Particle_Method(particle , dp, setup->kh, setup->kernel);
	particle->v->x += delta_t * residual->momentum_x_eq;
	particle->v->y += delta_t * residual->momentum_y_eq;

	// Update position with an Euler Implicit scheme
	particle->pos->x += delta_t * particle->v->x - delta_t * particle->XSPH_correction->x;
	particle->pos->y += delta_t * particle->v->y - delta_t * particle->XSPH_correction->y;


	// Update pressure with Tait's equation of state
	double B = squared(particle->param->sound_speed) * particle->param->rho_0 / particle->param->gamma;
	// double B = 0.85*1e5;
	particle->P = B * (pow(particle->rho / particle->param->rho_0, particle->param->gamma) - 1);

}


// Normal should be available everywhere!
/*
double compute_curvature(Particle *particle, Setup *setup, double epsilon) {
	double num = epsilon * 2 * compute_div(particle, Particle_get_normal, setup->kernel, setup->kh);
	// printf("%lf\n", compute_div(particle, Particle_get_normal, setup->kernel, setup->kh));
	double denom = 0;
	Particle *pi = particle;
	ListNode *node = pi->neighborhood->head;
	while (node != NULL) {
		Particle *pj = node->v;
		//printf("%lf %lf\n", pj->normal->x, pj->normal->y);
		xy *grad_W = grad_kernel(pi->pos, pj->pos, setup->kh, setup->kernel);
		//xy *corrected_gradient_W = correct_grad(grad_W, pi, setup);//Here it would biais the curvature I think
		denom += sqrt(squared(pi->pos->x - pj->pos->x) + squared(pi->pos->y - pj->pos->y)) *
			(pj->m / pj->rho) * norm(grad_W);
		free(grad_W);
		node = node->next;
	}
	return num / denom;
}*/

double compute_curvature(Particle *particle, Setup *setup, double epsilon) {
    double num = 0.0;
    // printf("%lf\n", compute_div(particle, Particle_get_normal, setup->kernel, setup->kh));
    Particle *pi = particle;
    xy *denom = xy_new(0,0);//grad_kernel(pi->pos, pj->pos, setup->kh, setup->kernel);
    //denom->x=0.0;
    //denom->y=0.0;
    ListNode *node = pi->neighborhood->head;
    double Csi=0.0;

    //Constriction of Csi just like the book
    while (node != NULL) {
        Particle *pj = node->v;
        double mrho= pj->m/pj->rho;
        Csi += mrho*eval_kernel(pi->pos,pj->pos,setup->kh,setup->kernel);
        node = node->next;
    }
    //Constriction of kappa just like the book
    node = pi->neighborhood->head;
    while (node != NULL) {
        Particle *pj = node->v;
        if(pj != pi){
	        //mass density ratio
	        double mrho = pj->m/pj->rho;

	        //nomalized position-> (X_i-X_j) / ||X_i-X_j||^2
	        xy *ij = xy_new(pi->pos->x - pj->pos->x, pi->pos->y - pj->pos->y);
	        double norm_ij = squared(norm(ij));
	        xy *unit = xy_new(ij->x/norm_ij, ij->y/norm_ij);
	        //double unit = (pi->pos-pj->pos)/squared(norm(pi->pos - pj->pos));

	        // Gradient of W
	        xy *grad_W = grad_kernel(pi->pos, pj->pos, setup->kh, setup->kernel);

	        //Vectorial decompsoition of denom
	        denom->x+=grad_W->x*mrho;
	        denom->y+=grad_W->y*mrho;

	        // Numerator assembly
	        num+=mrho*(Csi-1)*((unit->x*grad_W->x)+(unit->y*grad_W->y));

	        free(grad_W);
    	}
        node = node->next;
    }
    //Kappa assembly
    return (num / norm(denom))/10;
}

void compute_XSPH_correction(Particle *pi, Kernel kernel, double kh, double epsilon) {
	xy_reset(pi->XSPH_correction);
	ListNode *node = pi->neighborhood->head;
	while (node != NULL) {
		Particle *pj = node->v;
		pi->XSPH_correction->x += (pj->m / pj->rho) * (pi->v->x - pj->v->x) * eval_kernel(pi->pos, pj->pos, kh, kernel);
		pi->XSPH_correction->y += (pj->m / pj->rho) * (pi->v->y - pj->v->y) * eval_kernel(pi->pos, pj->pos, kh, kernel);
		//printf("%lf\n", pj->m);
		node = node->next;
	}
	pi->XSPH_correction->x *= epsilon;
	pi->XSPH_correction->y *= epsilon;
	//printf("pos = (%lf, %lf), Cs = %lf\n", particle->pos->x, particle->pos->y, particle->Cs);
}

/////////////
void update_positions_ellipse(Grid* grid, Particle** particles, Particle_derivatives** particles_derivatives, Residual** residuals, int n_p, Setup* setup) {

  	// Compute Cs, the XSPH correction on the velocity, and the divergence of the positions
	for (int i = 0; i < n_p; i++) {
		compute_Cs(particles[i], setup->kernel, setup->kh);
		if (setup->XSPH_epsilon != 0.0) compute_XSPH_correction(particles[i], setup->kernel, setup->kh,setup->XSPH_epsilon);
	}
	// Compute derivatives and residuals
	for (int i = 0; i < n_p; i++) {
		particles_derivatives[i]->div_v = compute_div(particles[i], Particle_get_v, setup->kernel, setup->kh);
		particles_derivatives[i]->lapl_v->x = compute_lapl(particles[i], Particle_get_v_x, setup->kernel, setup->kh);
		particles_derivatives[i]->lapl_v->y = compute_lapl(particles[i], Particle_get_v_y, setup->kernel, setup->kh);
		compute_grad(particles[i], Particle_get_P, setup->kernel, setup->kh, particles_derivatives[i]->grad_P);
		compute_grad(particles[i], Particle_get_Cs, setup->kernel, setup->kh, particles_derivatives[i]->grad_Cs);
		particles_derivatives[i]->lapl_Cs = compute_lapl(particles[i], Particle_get_Cs, setup->kernel, setup->kh);
		assemble_residual_NS(particles[i], particles_derivatives[i], residuals[i],setup);
	}
	int index_x_max, index_x_min, index_y_max, index_y_min;
	double pos_x_max = -INFINITY, pos_x_min = INFINITY, pos_y_max = -INFINITY, pos_y_min = INFINITY;
	// Integrate (new values (i.e. density, velocities) at time t+1)
	for (int i = 0; i < n_p; i++) {
		time_integrate(particles[i], residuals[i], setup->timestep);
		if (particles[i]->pos->x > pos_x_max) pos_x_max = particles[i]->pos->x, index_x_max = i;
		if (particles[i]->pos->x < pos_x_min) pos_x_min = particles[i]->pos->x, index_x_min = i;
		if (particles[i]->pos->y > pos_y_max) pos_y_max = particles[i]->pos->y, index_y_max = i;
		if (particles[i]->pos->y < pos_y_min) pos_y_min = particles[i]->pos->y, index_y_min = i;
	}
	// Compute semi-major axis of an ellipse
	double a_ellipse = particles[index_x_max]->pos->x - particles[index_x_min]->pos->x;
	double b_ellipse = particles[index_y_max]->pos->y - particles[index_y_min]->pos->y;
	// printf("a = %lf, b = %lf\n", a_ellipse * 0.5, b_ellipse * 0.5);
}

void update_positions_test_static_bubble(Grid* grid, Particle** particles, Particle_derivatives** particles_derivatives, Residual** residuals, int n_p, Setup* setup) {

	// Compute Cs, the XSPH correction on the velocity, and the divergence of the positions
	int index_x_max, index_x_min, index_y_max, index_y_min;
	double pos_x_max = -INFINITY, pos_x_min = INFINITY, pos_y_max = -INFINITY, pos_y_min = INFINITY;
	for (int i = 0; i < n_p; i++) {
		compute_Cs(particles[i], setup->kernel, setup->kh);
		if (setup->XSPH_epsilon != 0.0) compute_XSPH_correction(particles[i], setup->kernel, setup->kh,setup->XSPH_epsilon);
		// Compute radius of circle
		if (particles[i]->pos->x > pos_x_max) pos_x_max = particles[i]->pos->x, index_x_max = i;
		if (particles[i]->pos->x < pos_x_min) pos_x_min = particles[i]->pos->x, index_x_min = i;
		if (particles[i]->pos->y > pos_y_max) pos_y_max = particles[i]->pos->y, index_y_max = i;
		if (particles[i]->pos->y < pos_y_min) pos_y_min = particles[i]->pos->y, index_y_min = i;
	}
	double radius_circle = 0.5*(particles[index_x_max]->pos->x - particles[index_x_min]->pos->x);

	// Compute derivatives and normal
	for (int i = 0; i < n_p; i++) {
		particles_derivatives[i]->div_v = compute_div(particles[i], Particle_get_v, setup->kernel, setup->kh);
		particles_derivatives[i]->lapl_v->x = compute_lapl(particles[i], Particle_get_v_x, setup->kernel, setup->kh);
		particles_derivatives[i]->lapl_v->y = compute_lapl(particles[i], Particle_get_v_y, setup->kernel, setup->kh);
		compute_grad(particles[i], Particle_get_P, setup->kernel, setup->kh, particles_derivatives[i]->grad_P);
		compute_grad(particles[i], Particle_get_Cs, setup->kernel, setup->kh, particles_derivatives[i]->grad_Cs);
		particles_derivatives[i]->lapl_Cs = compute_lapl(particles[i], Particle_get_Cs, setup->kernel, setup->kh);
// 		assemble_residual_NS(particles[i], particles_derivatives[i], residuals[i], setup);
// 		assemble_residual_NS_test(particles[i], particles_derivatives[i], residuals[i], radius_circle, setup);
		compute_normal(particles[i], particles_derivatives[i]);
	}

	// Assemble residual and compute curvature
	for (int i = 0; i < n_p; i++) {
	    particles[i]->kappa = 2.0*compute_div(particles[i], Particle_get_normal, setup->kernel, setup->kh);
	    assemble_residual_NS_test(particles[i], particles_derivatives[i], residuals[i], radius_circle, setup);
	}

	// Integrate (obtain new values, i.e. density, velocities, pressure and positions, at time t+1)
	for (int i = 0; i < n_p; i++)
		time_integrate(particles[i], residuals[i], setup->timestep);
}

// Assemble the residual of the (incompressible) Navier-Stokes equations based on the derivatives available
void assemble_residual_NS_test(Particle* particle, Particle_derivatives* particle_derivatives, Residual* residual, double radius_circle, Setup* setup) {
	double mu_i = particle->param->dynamic_viscosity;

	double rho_i = particle->rho;
	double div_vel_i = particle_derivatives->div_v;
	xy* grad_P = particle_derivatives->grad_P;
	xy* lapl_v = particle_derivatives->lapl_v;


	xy *n = particle_derivatives->grad_Cs; // surface normal
	double norm_n = norm(n); // norm of n
	double lapl_Cs = particle_derivatives->lapl_Cs;
	double kappa = - lapl_Cs / norm_n; // curvature
	double kappa_2 = particle->kappa;

	// Exact values of normal and curvature for a circle centered in (0,0)
	xy* n_exact = xy_new(particle->pos->x, particle->pos->y);
	double norm_n_exact = norm(n_exact);
	double kappa_exact = 1.0 / radius_circle;

	double fs_x = 0; double fs_y = 0;
	// Apply surface tension only on particles in the vicinity the interface
	if (particle->on_free_surface) {
	    fs_x = - particle->param->sigma * kappa_exact * n->x / norm_n;
	    fs_y = - particle->param->sigma * kappa_exact * n->y / norm_n;
// 	  printf("pos = (%lf, %lf), n_exact = (%lf, %lf), n = (%lf, %lf), ||n|| = %lf, fs = (%lf, %lf), kappa_exact = %2.3f, kappa = %2.6f \n",
//	 particle->pos->x, particle->pos->y,-n_exact->x / norm_n_exact, -n_exact->y / norm_n_exact, n->x / norm_n, n->y / norm_n, norm_n, fs_x, fs_y, kappa_exact, kappa);
	  // printf("kappa_exact = %2.3f, kappa = %2.6f, kappa_div_n = %2.6f \n", kappa_exact, kappa, kappa_2);


	}

	residual->mass_eq = -rho_i * div_vel_i;
	residual->momentum_x_eq = (-1.0/rho_i) * grad_P->x + (mu_i/rho_i) * lapl_v->x + fs_x;
	residual->momentum_y_eq = (-1.0/rho_i) * grad_P->y + (mu_i/rho_i) * lapl_v->y + fs_y;

}

double compute_admissible_dt(double safety_param, double h_p, double c_0, double rho_0, double mu, double sigma) {
  // Relations from "Simulation of surface tension in 2D and 3D with smoothed particle hydrodynamics method", Zhang (2010)
  double dt_1 = 0.25 * h_p / c_0; // propagation of sound waves
  double dt_2 = INFINITY;
  if (mu > 0.0) dt_2 = 0.25 * (h_p*h_p) / (mu / rho_0); // viscous diffusion
  double dt_3 = INFINITY;
  if (sigma > 0.0) dt_3 = 0.25 * sqrt((rho_0*pow(h_p,3))/(2*M_PI*sigma)); // surface tension (capillary waves)

  double dt_min_interm = fmin(dt_1, dt_2);
  return safety_param * fmin(dt_min_interm, dt_3);

}

Boundary* Boundary_new(double xleft, double xright, double ybottom, double ytop,double CR, double CF){
	Boundary* boundary = (Boundary*) malloc(sizeof(Boundary));
	boundary->xleft = xleft;
	boundary->xright = xright;
	boundary->ybottom = ybottom;
	boundary->ytop = ytop;
	boundary->CR = CR;
	boundary->CF = CF;
	return boundary;
}
void Boundary_free(Boundary* boundary){
	free(boundary);
}

void center_reflection_right(Particle* pi, double CR, double Rp, double d){
		pi->pos->x -= (1+CR)*(Rp - d);
}
void center_reflection_left(Particle* pi, double CR, double Rp, double d){
		pi->pos->x += (1+CR)*(Rp - d);
}
void center_reflection_top(Particle* pi, double CR, double Rp, double d){
		pi->pos->y -= (1+CR)*(Rp - d);
}
void center_reflection_bottom(Particle* pi, double CR, double Rp, double d){
		pi->pos->y += (1+CR)*(Rp - d);
}
void velocity_reflection_vertical(Particle* pi, double CR, double CF){
	double vpN = pi->v->x;
	double vpT = pi->v->y;
	pi->v->x = -vpN*CR;
	pi->v->y = (1-CF)* vpT;
}
void velocity_reflection_horizontal(Particle* pi, double CR, double CF){
	double vpN = pi->v->y;
	double vpT = pi->v->x;
	pi->v->y = -vpN*CR;
	pi->v->x = (1-CF)* vpT;
}
void reflective_boundary(Particle** p, int n_p, Boundary* boundary, double Rp){
	// We have just computed the time integration. We correct the positions of the particles
	double CF = boundary->CF;
	double CR = boundary->CR;
	for(int i = 0; i < n_p ; i++){
		// For each particle we check if its position is close to a wall
		Particle* pi = p[i];
		double xright = boundary->xright;	double xleft = boundary->xleft;
		double ytop = boundary->ytop;			double ybottom = boundary->ybottom;
		double dxright,dxleft,dytop,dybottom;

		// Collision test
		int collision, collision_right, collision_left, collision_top, collision_bottom;
		collision = 0; collision_right = 0; collision_left = 0; collision_top = 0; collision_bottom = 0;

		int i = 0;
		while(pi->pos->x > boundary->xright - Rp){
			collision_right = 1;
			collision = 1;
			dxright = fabs(pi->pos->x - xright);
			center_reflection_right(pi,CR,Rp,dxright);
			velocity_reflection_vertical(pi,CR,CF);
			i++;
			// printf("i = %d\n",i);
			// printf("position : %f,%f\n", pi->pos->x, pi->pos->y);
		}
		while(pi->pos->x < boundary->xleft + Rp){
			collision_left = 1;
			collision = 1;
			dxleft = fabs(pi->pos->x - xleft);
			center_reflection_left(pi,CR,Rp,dxleft);
			velocity_reflection_vertical(pi,CR,CF);

		}
		while(pi->pos->y > boundary->ytop - Rp){
			collision_top = 1;
			collision = 1;
			dytop = fabs(pi->pos->y - ytop);
			center_reflection_top(pi,CR,Rp,dytop);
			velocity_reflection_horizontal(pi,CR,CF);
		}
		while(pi->pos->y < boundary->ybottom + Rp){
			collision_bottom = 1;
			collision = 1;
			dybottom = fabs(pi->pos->y - ybottom);
			center_reflection_bottom(pi,CR,Rp,dybottom);
			velocity_reflection_horizontal(pi,CR,CF);
		}
	}
}
