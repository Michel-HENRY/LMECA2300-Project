#ifndef LIST_H
#define LIST_H
#include <stdlib.h>
#include <stdio.h>


typedef struct List List; // Generic list
typedef struct ListNode ListNode;
typedef void(*destructor)(void *object);

struct List {
	int n; // list size
	ListNode *head; // pointer to first node for traversal
	ListNode *tail; // pointer to last node for append
};

struct ListNode {
	void *v; // node content
	ListNode *next; // pointer to next node
};


List* List_new();
void List_append(List *l, void *v);
void List_free(List*, destructor);

#endif
