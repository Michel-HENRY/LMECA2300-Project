#ifndef ANIMATION_H
#define ANIMATION_H

#include "BOV.h"
#include <math.h>
#include "vector.h"
#include "particle.h"

// see stringification process
#define xstr(s) str(s)
#define str(s) #s
#define PI 3.1415926535897932

typedef struct Animation Animation;

struct Animation {
	bov_window_t* window;
	bov_points_t* bov_particles;
	double timeout;
	int n_p;
	bov_points_t* grid;
    bov_points_t* domain;
    bov_text_t* plot_none;
    bov_text_t* plot_density;
    bov_text_t* plot_pressure;
    bov_text_t* plot_light;
    bov_text_t* plot_liquid;
    bov_points_t* shadow;
    bov_points_t* light;
};

Animation* Animation_new(int n_p,double timeout,Grid* grid, double R_p, double domain[4]);
void Animation_free(Animation* animation);
bov_points_t* load_Grid(Grid* grid);
bov_points_t* load_Domain(double domain[4]);
static void colormap(float v, float color[3]);
static void fillData(GLfloat (*data)[8], GLfloat(*data2)[2], Particle** particles, int n_p, int mask);
void show(Particle** particles, Animation* animation, int iter, bool wait, bool grid);
#endif
