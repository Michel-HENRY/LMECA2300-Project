#include "list.h"

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

void List_validation(){
	List* list = List_new();

	for(int i = 0; i < 5; i++){
		List_append(list, &i); //Ca crée une liste : [5,5,5,5,5]
	}
	ListNode* node = list->head;
	while(node != NULL){
		printf("v = %d\n", *(int*)node->v); //La puissance des pointeurs
		node = node->next;
	}

	List_free(list,NULL);
}
