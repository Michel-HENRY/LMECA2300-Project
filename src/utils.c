#include "utils.h"

List* List_new() {
	List *list = (List*)malloc(sizeof(List));
	list->n = 0;
	list->head = list->tail = NULL;
	return list;
}

void List_append(List *l, void *v) {
	ListNode* node = (ListNode*)malloc(sizeof(ListNode));
	node->v = v;
	node->next = NULL;

	if (l->n == 0)
		l->head = l->tail = node;
	else {
		l->tail->next = node;
		l->tail = node;
	}
	l->n++;
}

void List_free(List* list, destructor des) {
	ListNode* node1 = list->head;
	ListNode* node2 = NULL;
	for (int i = 0; i < list->n; i++) {
		if (des != NULL) des(node1->v);
		node2 = node1;
		node1 = node1->next;
		free(node2);
	}
	free(list);
}

xy* xy_new(double x, double y) {
	xy* pt = (xy*)malloc(sizeof(xy));
	pt->x = x, pt->y = y;
	return pt;
}

void xy_reset(xy *p) {
	p->x = 0;
	p->y = 0;
}

double rand_interval(double a, double b) {
	return (rand() / (double)RAND_MAX)*(b - a) + a;
}

double squared(double x) {
	return x*x;
}

double norm(xy *v) {
	return hypot(v->x, v->y);
}

xy* map_to_circle(xy* pos_square) {
    double pos_x_save = pos_square->x;
    double pos_y_save = pos_square->y;
    xy* pos_circle = xy_new(0.0, 0.0);
    pos_circle->x = pos_x_save*sqrt(1.0 - 0.5*pos_y_save*pos_y_save);
    pos_circle->y = pos_y_save*sqrt(1.0 - 0.5*pos_x_save*pos_x_save);
    
    return pos_circle;
}

xy* generate_circle(int k, int n, int nb, double radius) {
    double phi = (sqrt(5)+1)/2;
    double r;
    if (k > n-nb) r = radius;
    else r = radius*sqrt(k-1/2)/sqrt(n-(nb+1)/2);
    double theta = 2*M_PI*k/(phi*phi);
    xy* coord = xy_new(r*cos(theta), r*sin(theta));
    
    return coord;  
}