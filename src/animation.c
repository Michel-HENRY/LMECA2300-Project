#include "animation.h"

double P_max(Particle** p, int n_p){
	double Pmax = p[0]->fields->P;
	for(int i = 1; i < n_p ; i++){
		if(p[i]->fields->P > Pmax){
			Pmax = p[i]->fields->P;
		}
	}
	return Pmax;
}

double P_min(Particle** p, int n_p){
	double Pmin = p[0]->fields->P;
	for(int i = 1; i < n_p ; i++){
		if(p[i]->fields->P < Pmin){
			Pmin = p[i]->fields->P;
		}
	}
	return Pmin;
}

double rho_max(Particle** p, int n_p){
	double rmax = p[0]->fields->rho;
	for(int i = 1; i < n_p ; i++){
		if(p[i]->fields->rho > rmax){
			rmax = p[i]->fields->rho;
		}
	}
	return rmax;
}

double rho_min(Particle** p, int n_p){
	double rmin = p[0]->fields->rho;
	for(int i = 1; i < n_p ; i++){
		if(p[i]->fields->rho < rmin){
			rmin = p[i]->fields->rho;
		}
	}
	return rmin;
}

static double u_max(Particle** p, int n_p){
  double Pmax = sqrt(pow(p[0]->fields->u->X[0],2) + pow(p[0]->fields->u->X[1],2));
  for(int i = 1; i < n_p ; i++){
    // printf("i = %d/%d\n",i,n_p);
    if(sqrt(pow(p[i]->fields->u->X[0],2) + pow(p[i]->fields->u->X[1],2)) > Pmax){
      Pmax = sqrt(pow(p[i]->fields->u->X[0],2) + pow(p[i]->fields->u->X[1],2));
    }
  }
  return Pmax;
}
static double u_min(Particle** p, int n_p){
  double Pmin = sqrt(pow(p[0]->fields->u->X[0],2) + pow(p[0]->fields->u->X[1],2));
  for(int i = 1; i < n_p ; i++){
    // printf("i = %d/%d\n",i,n_p);
    if(sqrt(pow(p[i]->fields->u->X[0],2) + pow(p[i]->fields->u->X[1],2)) < Pmin){
      Pmin = sqrt(pow(p[i]->fields->u->X[0],2) + pow(p[i]->fields->u->X[1],2));
    }
  }
  return Pmin;
}

Animation* Animation_new(int n_p,double timeout, Grid* grid, double R_p, Vector** edges, int n_e){
  Animation* animation = (Animation*) malloc(sizeof(Animation));

	animation->timeout = timeout;
	animation->n_p = n_p;

	animation->window = bov_window_new(1024, 780, "Project: SPH");
	bov_window_set_color(animation->window, (GLfloat[]){0.0, 0.0f, 0.0f, 0.0f});
	bov_window_enable_help(animation->window);
	bov_window_set_zoom(animation->window, 0.8);

	GLfloat(*data)[8] = malloc(sizeof(data[0])*n_p);
	animation->bov_particles = bov_particles_new(data, n_p, GL_DYNAMIC_DRAW);
	bov_points_set_width(animation->bov_particles, R_p);
	bov_points_set_outline_width(animation->bov_particles, 0);

	GLfloat(*data2)[2] = malloc(sizeof(data2[0])*n_p*9);
	animation->shadow = bov_points_new(data2, n_p*9, GL_DYNAMIC_DRAW);
	bov_points_set_color(animation->shadow,(GLfloat[]){0.0, 0.0, 0.0, 0.2}); // 0.72
	bov_points_set_width(animation->shadow, 1e-15);

	GLfloat(*data3)[8] = malloc(sizeof(data3[0]));
	data3[0][0] = 2.0f; data3[0][1] = 1.5f;
	data3[0][2] = 0; data3[0][3] = 0; 
	data3[0][4] = 0.98f; data3[0][5] = 0.87f;
	data3[0][6] = 0.137f; data3[0][2] = 1; 
	animation->light = bov_particles_new(data3, 1, GL_STATIC_DRAW);
	bov_points_set_color(animation->light,(GLfloat[]){0.98f, 0.87f, 0.137f, 1.0});
	bov_points_set_outline_color(animation->light, (GLfloat[]){0.98f, 0.87f, 0.137f, 1.0});
	bov_points_set_marker(animation->light, 3);
	bov_points_set_width(animation->light, 3);
	bov_points_set_pos(animation->light, (GLfloat[]){0.0f, 0.0f});

	free(data);
	free(data2);
	free(data3);
	animation->n_e = n_e;
  	animation->domain = load_Domain(edges, n_e);
	bov_points_set_width(animation->domain, 0.01);
	bov_points_set_color(animation->domain,(GLfloat[]){0.0f, 0.0f, 0.0f, 1.0});

	if (grid != NULL){
		animation->grid = load_Grid(grid);
	} else{
		animation->grid = NULL;
	}

	animation->plot_none = bov_text_new((GLubyte[]) {
		"None\n"
	}, GL_STATIC_DRAW);

	animation->plot_pressure = bov_text_new((GLubyte[]) {
		"Pressure\n"
	}, GL_STATIC_DRAW);

	animation->plot_density = bov_text_new((GLubyte[]) {
		"Density\n"
	}, GL_STATIC_DRAW);
	animation->plot_velocity = bov_text_new((GLubyte[]) {
		"Velocity\n"
	}, GL_STATIC_DRAW);
	animation->plot_light = bov_text_new((GLubyte[]) {
		"Light\n"
	}, GL_STATIC_DRAW);
	animation->plot_liquid = bov_text_new((GLubyte[]) {
		"Liquid\n"
	}, GL_STATIC_DRAW);
	bov_text_set_space_type(animation->plot_none, PIXEL_SPACE);bov_text_set_space_type(animation->plot_pressure, PIXEL_SPACE);
	bov_text_set_space_type(animation->plot_density, PIXEL_SPACE);bov_text_set_space_type(animation->plot_velocity, PIXEL_SPACE);
	bov_text_set_space_type(animation->plot_light, PIXEL_SPACE);bov_text_set_space_type(animation->plot_liquid, PIXEL_SPACE);
	bov_text_set_fontsize(animation->plot_none, 16.0f);bov_text_set_fontsize(animation->plot_pressure, 16.0f);bov_text_set_fontsize(animation->plot_density, 16.0f);
	bov_text_set_fontsize(animation->plot_velocity, 16.0f);bov_text_set_fontsize(animation->plot_light, 16.0f);bov_text_set_fontsize(animation->plot_liquid, 16.0f);
	bov_text_set_pos(animation->plot_none, (GLfloat[2]){16.0f, 630.0f});
	bov_text_set_pos(animation->plot_pressure, (GLfloat[2]){16.0f, 650.0f});
	bov_text_set_pos(animation->plot_density, (GLfloat[2]){16.0f, 670.0f});
	bov_text_set_pos(animation->plot_velocity, (GLfloat[2]){16.0f, 690.0f});
	bov_text_set_pos(animation->plot_light, (GLfloat[2]){16.0f, 710.0f});
	bov_text_set_pos(animation->plot_liquid, (GLfloat[2]){16.0f, 730.0f});
	bov_text_set_color(animation->plot_none, (GLfloat[4]){0.0f, 0.8f, 0.0f, 1.0f});
  //bov_text_set_outline_width(window->indication, 0.5f);

	return animation;
}

void Animation_free(Animation* animation){
	bov_points_delete(animation->bov_particles);
	bov_points_delete(animation->shadow);
	bov_points_delete(animation->light);
	bov_text_delete(animation->plot_density);
	bov_text_delete(animation->plot_velocity);
	bov_text_delete(animation->plot_none);
	bov_text_delete(animation->plot_pressure);
	bov_text_delete(animation->plot_light);
	bov_text_delete(animation->plot_liquid);
	if(animation->domain = NULL){
		bov_points_delete(animation->domain);
	}
	if(animation->grid != NULL)
		bov_points_delete(animation->grid);
	bov_window_delete(animation->window);
	free(animation);
}
bov_points_t* load_Domain(Vector** edges, int n_e){
  GLfloat(*data)[2] = malloc(sizeof(data[0])*n_e);
  for(int i = 0 ; i < n_e;i ++){
    data[i][0] = edges[2*i]->X[0];
    data[i][1] = edges[2*i]->X[1];
    // printf("noeud %d :\t (%f,%f)\n", i, data[i][0], data[i][1]);
  }
  bov_points_t *points = bov_points_new(data, n_e, GL_STATIC_DRAW);
  free(data);
  return points;
}

bov_points_t* load_Grid(Grid* grid){
	int nLines = (grid->nCellx + 1) + (grid->nCelly + 1);
	GLfloat(*data)[2] = malloc(sizeof(data[0])*2*nLines);
	for (int i = 0;i < (grid->nCellx + 1);i++)
	{
		data[2 * i][0] = grid->left + i*grid->h;
		data[2 * i][1] = grid->bottom;
		data[2 * i + 1][0] = grid->left + i*grid->h;
		data[2 * i + 1][1] = grid->top;
	}
	int shift = 2 * (grid->nCellx + 1);
	for (int j = 0;j < (grid->nCelly + 1);j++)
	{
		data[shift + 2 * j][0] = grid->left;
		data[shift + 2 * j][1] = grid->bottom + j*grid->h;
		data[shift + 2 * j + 1][0] = grid->right;
		data[shift + 2 * j + 1][1] = grid->bottom + j*grid->h;
	}
	bov_points_t *points = bov_points_new(data, 2*nLines, GL_STATIC_DRAW);
	bov_points_set_width(points, 0.001);
	double L = grid->h*grid->nCellx;

	free(data);
	return points;
}

static void colormap(float v, float color[3]){
	float v1 = 3.5*(v-0.7);
	float v2 = 1.25*v;
	float v3 = fminf(0.5,v)*2.0;

  // color[0] = -v1*v1+1.0f;
  // color[1] = 6.0f*v2*v2*(1.0f-v2);
  // color[2] = 5.5f*v3*(1.0f-v3)*(1.0f-v3);

  /* alternative: classical jet colormap */
	color[0] = 1.5 - 4.0 * fabs(v - 0.75);
	color[1] = 1.5 - 4.0 * fabs(v - 0.5 );
	color[2] = 1.5 - 4.0 * fabs(v - 0.25);
}

static void fillData(GLfloat (*data)[8], GLfloat(*data2)[2], Particle** particles, int n_p, int mask){
	double light_X = 1.0f;
	double light_Y = 1.0f;

	double P_M = P_max(particles, n_p);
	double P_m = P_min(particles, n_p);
	double r_M = rho_max(particles, n_p);
	double r_m = rho_min(particles, n_p);
	double u_M = u_max(particles, n_p);
	double u_m = u_min(particles, n_p);
	for(int i=0; i<n_p; i++) {
		Particle* p = particles[i];
	    data[i][0] = p->fields->x->X[0]; // x (rand between -100 and 100)
	    data[i][1] = p->fields->x->X[1]; // y (rand between -100 and 100)
	    data[i][2] = p->fields->u->X[0]; // speed x (not used by default visualization)
	    data[i][3] = p->fields->u->X[1]; // speed y (not used by default visualization)
	    data[i][4] = 0;
	    data[i][5] = 0;
	    data[i][6] = 0;
	    // data[i][7] = 0;
	    if(mask%6 == 5){//Fluid
	    	//Set the marker for the fluid
	    	colormap(0.2, &data[i][4]);
	    	data[i][2] = -1000;
	    }
	    if(mask%6 == 4){//Shadow
	    	//Shadow
	    	//colormap(0.1, &data[i][4]);
	    	data[i][4] = 0.035;
		    data[i][5] = 0.267;
		    data[i][6] = 0.56;
		    double x = p->fields->x->X[0];
		    double y = p->fields->x->X[1];
		    data2[9*i][0] = x;
		    data2[9*i][1] = y;
		    double dist = sqrt((x-light_X)*(x-light_X) + (y-light_Y)*(y-light_Y));
		    double theta = atan((y-light_Y)/(x-light_X));
		    double alpha = atan(p->param->Rp / dist);
		    int mx = 1; 
		    if(light_X > x)
		    	mx = -1;
		    int my = 1;
		    if(light_Y > y)
		    	my = -1;
		    data2[9*i + 1][0] = 3*mx + light_X;
		    data2[9*i + 1][1] = 3*my*tan(theta-alpha) + light_Y;
		    data2[9*i + 2][0] = 3*mx + light_X;
		    data2[9*i + 2][1] = 3*my*tan(theta+alpha) + light_Y;

		    data2[9*i + 3][0] = mx*p->param->Rp*cos(theta-alpha+PI/2) + x;
		    data2[9*i + 3][1] = my*p->param->Rp*sin(theta-alpha+PI/2) + y;
		    data2[9*i + 4][0] = mx*p->param->Rp*cos(theta-alpha-PI/2) + x;
		    data2[9*i + 4][1] = my*p->param->Rp*sin(theta-alpha-PI/2) + y;
		    data2[9*i + 5][0] = data2[9*i + 1][0];
		    data2[9*i + 5][1] = data2[9*i + 1][1];

		    data2[9*i + 6][0] = data2[9*i + 3][0];
		    data2[9*i + 6][1] = data2[9*i + 3][1];
		    data2[9*i + 7][0] = data2[9*i + 1][0];
		    data2[9*i + 7][1] = data2[9*i + 1][1];
		    data2[9*i + 8][0] = data2[9*i + 2][0];
		    data2[9*i + 8][1] = data2[9*i + 2][1];

		    data[i][2] = -2000;
		} else if(mask%6 == 3){
		  colormap((p->fields->rho-r_m)/(r_M-r_m), &data[i][4]); // fill colormap//Useless --- to remove
		  data[i][5] = 0;
		  data[i][4] = 0;
		  data[i][6] = 0;
		  double u = sqrt(pow(p->fields->u->X[0],2) + pow(p->fields->u->X[1],2));
		  data[i][2] = 4000 + ((u-u_m)/(u_M-u_m));	
	    } else if(mask%6 == 2){
	      colormap((p->fields->rho-r_m)/(r_M-r_m), &data[i][4]); // fill colormap//Useless --- to remove
	      data[i][5] = 0;
	      data[i][4] = 0;
	      data[i][6] = 0;
	      data[i][2] = 4000 + ((p->fields->rho-r_m)/(r_M-r_m));
	  	} else if(mask%6 == 1){
	      colormap((p->fields->P-P_m)/(P_M-P_m), &data[i][4]); // fill colormap
	      data[i][5] = 0;
	      data[i][4] = 0;
	      data[i][6] = 0;
	      data[i][2] = 4000 + ((p->fields->P-P_m)/(P_M-P_m));
	  	} else{
	  	  colormap(0.2, &data[i][4]); // fill colormap
	  	}
	    data[i][7] = 0.8f; // transparency

	}
}
void show(Particle** particles, Animation* animation, int iter, bool wait, bool grid){
	int n_p = animation->n_p;
  // Update bov_particles
	bov_window_t* window = animation->window;
	int nbr = bov_window_get_counter(window);
	/*
	if(iter%600 <= 100){
		nbr = 1;
	} else if(iter%600 <= 200){
		nbr = 3;
	} else if(iter%600 <= 300){
		nbr = 2;
	} else if(iter%600 <= 400){
		nbr = 4;
	} else{
		nbr = 5;
	}*/
	bov_points_set_width(animation->bov_particles, particles[0]->param->Rp);
	bov_text_set_color(animation->plot_none, (GLfloat[4]){0.4f, 0.3f, 0.4f, 1.0f});
	bov_text_set_color(animation->plot_density, (GLfloat[4]){0.4f, 0.3f, 0.4f, 1.0f});
	bov_text_set_color(animation->plot_velocity, (GLfloat[4]){0.4f, 0.3f, 0.4f, 1.0f});
	bov_text_set_color(animation->plot_pressure, (GLfloat[4]){0.4f, 0.3f, 0.4f, 1.0f});
	bov_text_set_color(animation->plot_light, (GLfloat[4]){0.4f, 0.3f, 0.4f, 1.0f});
	bov_text_set_color(animation->plot_liquid, (GLfloat[4]){0.4f, 0.3f, 0.4f, 1.0f});
	if(nbr%6 == 5){
		bov_text_set_color(animation->plot_liquid, (GLfloat[4]){0.0f, 0.8f, 0.0f, 1.0f});
		bov_points_set_width(animation->bov_particles, particles[0]->param->Rp*4);
	} else if(nbr%6 == 4)
		bov_text_set_color(animation->plot_light, (GLfloat[4]){0.0f, 0.8f, 0.0f, 1.0f});
	else if(nbr%6 == 3){
		bov_text_set_color(animation->plot_velocity, (GLfloat[4]){0.0f, 0.8f, 0.0f, 1.0f});
		bov_points_set_width(animation->bov_particles, particles[0]->param->Rp*1.8);
	} else if(nbr%6 == 2){
		bov_text_set_color(animation->plot_density, (GLfloat[4]){0.0f, 0.8f, 0.0f, 1.0f});
		bov_points_set_width(animation->bov_particles, particles[0]->param->Rp*1.8);
	} else if(nbr%6 == 1){
		bov_text_set_color(animation->plot_pressure, (GLfloat[4]){0.0f, 0.8f, 0.0f, 1.0f});
		bov_points_set_width(animation->bov_particles, particles[0]->param->Rp*1.8);
	} else
		bov_text_set_color(animation->plot_none, (GLfloat[4]){0.0f, 0.8f, 0.0f, 1.0f});

	GLfloat(*data)[8] = malloc(sizeof(data[0])*n_p);
	GLfloat(*data2)[2] = malloc(sizeof(data2[0])*n_p*9);
	fillData(data, data2, particles, n_p, nbr);
	bov_points_t* bov_particles = animation->bov_particles;
	bov_particles = bov_particles_update(bov_particles, data,n_p);
	if(nbr%6 == 4){
		bov_points_t* shadow = animation->shadow;
		shadow = bov_points_update(shadow, data2, n_p*9);
	}
	free(data);
	free(data2);

  	// To make the screenshot

	char screenshot_name[64] = "animation_";
	char int_string[32];
	sprintf(int_string, "%d", iter);
	strcat(screenshot_name, int_string);
	double tbegin = bov_window_get_time(window);
	double timeout = animation->timeout;
	if(!wait){
		while(bov_window_get_time(window) - tbegin < timeout){
			if(animation->grid != NULL && grid)
				bov_lines_draw(window,animation->grid,0, BOV_TILL_END);
			//bov_fast_triangle_fan_draw(window, animation->domain, 0, BOV_TILL_END);
			if(nbr%6 == 4){
				bov_particles_draw(window, animation->light, 0, BOV_TILL_END);
			}
			bov_particles_draw(window, animation->bov_particles, 0, BOV_TILL_END);
			if(nbr%6 == 4){
				bov_triangles_draw(window,animation->shadow, 0, BOV_TILL_END);
			}
			bov_line_loop_draw(window, animation->domain,0,BOV_TILL_END);
			bov_text_draw(window, animation->plot_none);
			bov_text_draw(window, animation->plot_pressure);
			bov_text_draw(window, animation->plot_density);
			bov_text_draw(window, animation->plot_velocity);
			bov_text_draw(window, animation->plot_light);
			bov_text_draw(window, animation->plot_liquid);
			if (iter%10 == 0) {
				bov_window_screenshot(window, screenshot_name);
			}
			bov_window_update(window);
		}
	} else {
    // we want to keep the window open with everything displayed
		while (!bov_window_should_close(window)) {
			if(animation->grid != NULL && grid)
				bov_lines_draw(window,animation->grid,0, BOV_TILL_END);
			//bov_fast_triangle_fan_draw(window, animation->domain, 0, BOV_TILL_END);
			if(nbr%6 == 4){
				bov_particles_draw(window, animation->light, 0, BOV_TILL_END);
			}
			bov_particles_draw(window, animation->bov_particles, 0, BOV_TILL_END);
			if(nbr%6 == 4){
				bov_triangles_draw(window,animation->shadow, 0, BOV_TILL_END);
			}
			bov_line_loop_draw(window, animation->domain,0,BOV_TILL_END);
			bov_text_draw(window, animation->plot_none);
			bov_text_draw(window, animation->plot_pressure);
			bov_text_draw(window, animation->plot_density);
			bov_text_draw(window, animation->plot_velocity);
			bov_text_draw(window, animation->plot_light);
			bov_text_draw(window, animation->plot_liquid);
			bov_window_screenshot(window, screenshot_name);
			bov_window_update_and_wait_events(window);
		}
	}
}
/*
Pour la lumiere on va mettre une particule géante :o
=======
static void fillData(GLfloat (*data)[8], Particle** particles, int n_p){
  double Pmax = P_max(particles, n_p);
  double Pmin = P_min(particles, n_p);
  double umax = u_max(particles, n_p);
  double umin = u_min(particles, n_p);
	for(int i=0; i<n_p; i++) {
    // printf("i = %d/%d\n",i,n_p);
    Particle* p = particles[i];
		data[i][0] = p->fields->x->X[0]; // x (rand between -100 and 100)
		data[i][1] = p->fields->x->X[1]; // y (rand between -100 and 100)
		data[i][2] = p->fields->u->X[0]; // speed x (not used by default visualization)
		data[i][3] = p->fields->u->X[1]; // speed y (not used by default visualization)
		data[i][4] = 0;
		data[i][5] = 0;
		data[i][6] = 0;
		// data[i][7] = 0;
    // double P = p->fields->P;
    double u = sqrt(pow(p->fields->u->X[0],2) + pow(p->fields->u->X[1],2));
    double field;
    // if(Pmax == Pmin){
    //   field = 0.0;
    // }
    // else{
    //   field = (P - Pmin)/(Pmax - Pmin);
    // }
    if(umax == umin){
      field = 0.0;
    }
    else{
      field = (u - umin)/(umax - umin);
    }

		colormap(field, &data[i][4]); // fill color
		data[i][7] = 0.8f; // transparency
  }
}
void show(Particle** particles, Animation* animation, int iter, bool wait, bool grid){
  int n_p = animation->n_p;
  //------------------------ UPDATE PARTICLES ----------------------------------
  GLfloat(*data)[8] = malloc(sizeof(data[0])*n_p);
	fillData(data, particles, n_p);
  bov_points_t* bov_particles = animation->bov_particles;
  bov_particles = bov_particles_update(bov_particles, data,n_p);
  free(data);

  // --------------------------- SCREENSHOT ------------------------------------
  char screenshot_name[64] = "moving_circle_";
	char int_string[32];
	sprintf(int_string, "%d", iter);
	strcat(screenshot_name, int_string);
  int SCREENSHOT_STEP = 100;

  // -------------------------- CREATE THE WINDOW ------------------------------
  bov_window_t* window = animation->window;
  double tbegin = bov_window_get_time(window);
  double timeout = animation->timeout;

  // ------------------------------ SHOW ---------------------------------------
  if(!wait){
    while(bov_window_get_time(window) - tbegin < timeout){
      if(animation->grid != NULL && grid)
				bov_lines_draw(window,animation->grid,0, BOV_TILL_END);
    bov_particles_draw(window, animation->bov_particles, 0, BOV_TILL_END);
    bov_line_loop_draw(window, animation->domain,0,BOV_TILL_END);
    if (iter%SCREENSHOT_STEP == 0) {
      bov_window_screenshot(window, screenshot_name);
    }
    bov_window_update(window);
    }
	}
  // ----------------------- Keep it open if it is the end ---------------------
  else {
    while (!bov_window_should_close(window)) {
      if (animation->grid != NULL && grid)
        bov_lines_draw(window, animation->grid, 0, BOV_TILL_END);
      bov_particles_draw(window, animation->bov_particles, 0, BOV_TILL_END);
      bov_line_loop_draw(window, animation->domain,0,BOV_TILL_END);
      bov_window_screenshot(window, screenshot_name);
      bov_window_update_and_wait_events(window);
    }
  }
}
>>>>>>> a0c07f0b7ca9bf980f4894733c63d8008c661bee*/
