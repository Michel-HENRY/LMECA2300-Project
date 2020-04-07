#ifndef VECTOR_H
#define VECTOR_H

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>

typedef struct Vector Vector;

struct Vector {
  int DIM;
  double* X;
};

// Creator && Destructor
Vector* Vector_new(int DIM);
void Vector_initialise(Vector* v, double* X);
void Vector_free(Vector* v);

// Utils
void check_dim(Vector* v1, Vector* v2);
void Vector_print(Vector* v);

// Operation
Vector* sum(Vector* v1, Vector* v2);
void sum_into(Vector* v1, Vector* v2);
Vector* diff(Vector* v1, Vector* v2);
double dot(Vector* v1, Vector* v2);
double norm(Vector* v1);
double dist(Vector* v1, Vector* v2);
bool equal(Vector* v1, Vector* v2);
Vector* times(Vector* v, double a);
void times_into(Vector* v, double a);
Vector* prod(Vector* v1, Vector* v2);

void Vector_validation();

#endif
