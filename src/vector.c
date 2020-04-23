#include "vector.h"
#include <math.h>

// -------------------------------------------------------------------
// -------------------------------------------------------------------
// -------------------------------------------------------------------
// -------------------- Creator && Destructor ------------------------
// -------------------------------------------------------------------
// -------------------------------------------------------------------
// -------------------------------------------------------------------
Vector* Vector_new(int DIM){
  Vector* v = (Vector*) malloc(sizeof(Vector));
  v->DIM = DIM;
  v->X = (double*) malloc(DIM*sizeof(double));
  for(int i = 0; i < DIM ; i++){
    v->X[i] = 0.0;
  }
  return v;
}
void Vector_initialise(Vector* v, double* X){
  for(int i = 0; i < v->DIM ; i++){
    if (X != NULL){
      v->X[i] = X[i];
    }
  }
}
Vector* copy(Vector* v){
  Vector* u = (Vector*) malloc(sizeof(Vector));
  u->DIM = v->DIM;
  u->X = (double*) malloc(u->DIM*sizeof(double));
  for(int i = 0; i < v->DIM;i++){
      u->X[i] = v->X[i];
  }
  return u;
}
void Vector_free(Vector* v){
  if( v != NULL){
    if(v->X != NULL){
      free(v->X);
    }
    free(v);
  }
}
// -------------------------------------------------------------------
// -------------------------------------------------------------------
// -------------------------------------------------------------------
// ----------- Check dimension and print the Vector components -------
// -------------------------------------------------------------------
// -------------------------------------------------------------------
// -------------------------------------------------------------------
void check_dim(Vector* v1, Vector* v2){
  if(v1->DIM != v2->DIM){
    printf("ERROR : INVALID DIMENSION\n");
    exit(0);
  }
}
void Vector_print(Vector* v){
  printf("\nv = [");
  for(int i = 0; i < v->DIM; i++){
    printf("%e\t", v->X[i]);
  }
  printf("]\n");
}
// -------------------------------------------------------------------
// -------------------------------------------------------------------
// -------------------------------------------------------------------
// Operators : sum, diff, dot, norm, dist
// -------------------------------------------------------------------
// -------------------------------------------------------------------
// -------------------------------------------------------------------
Vector* sum(Vector* v1, Vector* v2){
  check_dim(v1,v2);
  Vector* w = Vector_new(v1->DIM);
  for(int i = 0; i < v1->DIM; i++){
    w->X[i] = v1->X[i] + v2->X[i];
  }
  return w;
}
void sum_into(Vector* v1, Vector* v2){
  check_dim(v1,v2);
  for(int i = 0; i < v1->DIM; i++){
    v1->X[i] += v2->X[i];
  }
}
Vector* diff(Vector* v1, Vector* v2){
  check_dim(v1,v2);
  Vector* w = Vector_new(v1->DIM);
  for(int i = 0; i < v1->DIM; i++){
    w->X[i] = v1->X[i] - v2->X[i];
  }
  return w;
}
double dot(Vector* v1, Vector* v2){
  check_dim(v1,v2);
  double result = 0;
  for(int i = 0; i < v1->DIM; i++){
    result += v1->X[i]*v2->X[i];
  }
  return result;
}
double norm(Vector* v1){
  double result = 0;
  for(int i = 0; i < v1->DIM; i++){
    result += v1->X[i]*v1->X[i];
  }
  return sqrt(result);
}
double dist(Vector* v1, Vector* v2){
  Vector* v3 = diff(v1,v2);
  double dist = norm(v3);
  Vector_free(v3);

  return dist;
}
bool equal(Vector* v1, Vector* v2){
  double eps = 1e-8;
  if(dist(v1,v2) < eps) return true;
  return false;
}
Vector* times(Vector* v, double a){
  Vector* result = Vector_new(v->DIM);
  for(int i = 0; i < v->DIM; i++){
    result->X[i] = a*v->X[i];
  }
  return result;
}
void times_into(Vector* v, double a){
  for(int i = 0; i < v->DIM; i++){
    v->X[i] *= a;
  }
}
Vector* prod(Vector* v1, Vector* v2){
  check_dim(v1,v2);
  Vector* result = Vector_new(v1->DIM);
  for(int i = 0; i < v1->DIM; i++){
    result->X[i] = v1->X[i]*v2->X[i];
  }
  return result;
}
void put_zeros(Vector* v){
  for(int i = 0; i < v->DIM; i++){
    v->X[i] = 0.0;
  }
}
// -------------------------------------------------------------------
// -------------------------------------------------------------------
// -------------------------------------------------------------------
// ---------------------- Validation testcase ------------------------
// -------------------------------------------------------------------
// -------------------------------------------------------------------
// -------------------------------------------------------------------
void Vector_validation(){
  int DIM = 2;
  Vector* v1 = Vector_new(DIM);
  Vector* v2 = Vector_new(DIM);

  double x1[2] = {3,4};
  double x2[2] = {2,1};

  Vector_initialise(v1,x1);
  Vector_initialise(v2,x2);

  Vector_print(v1);
  Vector_print(v2);

  double prod = dot(v1,v2);
  printf("dot = %f\n",prod);

  Vector* v3 = sum(v1,v2);
  Vector_print(v3);

  Vector* v4 = diff(v1,v2);
  Vector_print(v4);

  double norm_v1 = norm(v2);
  printf("norm = %f\n", norm_v1);

  double dist_v1v2 = dist(v1,v2);
  double dist_v2v1 = dist(v2,v1);

  printf("dist v1->v2 = %f\n", dist_v1v2);
  printf("dist v2->v1 = %f\n", dist_v2v1);

  Vector_free(v1);
  Vector_free(v2);
  Vector_free(v3);
  Vector_free(v4);
}
