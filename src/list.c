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
