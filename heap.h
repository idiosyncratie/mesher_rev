#ifndef INCLUSION_HEAP_H
#define INCLUSION_HEAP_H

#include <vector>
#include <iostream>

#define NOT_IN_HEAP -1

template <class T> class heap_node {
	public:
    double val;
	int loc;
    T *obj;
};

template <class T> class Heap {
public:
    Heap<T>() { size = 0;}
    ~Heap<T>() {}

    heap_node<T>* operator[] (int i) { return node[i]; }
    int heap_size() { return size; }

    void insert(T *obj,double v);
	void insert(heap_node<T>*);
    heap_node<T> *extract();
    heap_node<T> *top() { return size<1 ? NULL : node[0]; };
    heap_node<T> *kill(int i);
    void update(int,double);

private:
	std::vector<heap_node<T> *> node;
    int size;

    void swap(int i,int j);

    int parent(int i) { return (i-1)/2; }
    int left(int i) { return 2*i+1; }
    int right(int i) { return 2*i+2; }

    void upheap(int i);
    void downheap(int i);
};

// Heap::swap --
//
// Swaps two nodes in the heap.
//
template <class T> void Heap<T>::swap(int i,int j)
{
	heap_node<T>* tmp = node[i];

	node[i] = node[j];
	node[j] = tmp;

	node[i]->loc = i;
	node[j]->loc = j;
}


// Heap::upheap --
//
// The given node will be moved up in the heap, if necessary.
//
// NOTE: This function (as well as downheap) performs more swapping
// than is strictly necessary.
//
template <class T> void Heap<T>::upheap(int i)
{
	if( i==0 ) return;

	if( node[i]->val > node[parent(i)]->val ) {
		swap(i,parent(i));
		upheap(parent(i));
	}
}

// Heap::downheap --
//
// The given node is moved down through the heap, if necessary.
//
template <class T> void Heap<T>::downheap(int i)
{
	if (i>=size) return;	// perhaps just extracted the last

	int largest = i,
	l = left(i),
	r = right(i);

	if( l<size && node[l]->val > node[largest]->val ) largest = l;
	if( r<size && node[r]->val > node[largest]->val ) largest = r;

	if( largest != i ) {
		swap(i,largest);
		//Continue moving node i (now in position largest), down the heap
		downheap(largest);
	}
}

template <class T> void Heap<T>::insert(heap_node<T> *n) {
	int i = size++;
	
	n->loc = i;
	node.push_back(n);

	upheap(i);
}

// Heap::insert --
//
// Insert the given object into the heap using the specified key value.
//
template <class T> void Heap<T>::insert(T *obj,double v)
{
	heap_node<T> n;
	n->obj = obj;
	n->val = v;

	this->insert(n);
}

// Heap::extract --
//
// Extract the top element from the heap and return it.
//
template <class T> heap_node<T> *Heap<T>::extract()
{
	if( size<1 ) return 0;

	swap(0,size-1);
	size--;

	downheap(0);

	node[size]->loc = NOT_IN_HEAP;

	return node[size];
}

// Heap::kill --
//
// Kill a given node in the heap.
//
template <class T> heap_node<T>* Heap<T>::kill(int i)
{
	if( i>=size )
		cerr << "ATTEMPT TO DELETE OUTSIDE OF RANGE\n";

	swap(i,size-1);
	size--;

	if( node[i]->val < node[size]->val )
		downheap(i);
	else
		upheap(i);

	node[size]->loc = NOT_IN_HEAP;

	return node[size];
}


// Heap::update --
//
// This function is called when the key value of the given node has
// changed.  It will record this change and reorder the heap if
// necessary.
//
template <class T> void Heap<T>::update(int i,double v)
{
	double old=node[i]->val;
	node[i]->val = v;

	if( v<old ) {
		downheap(i);
	} else if (v > old) {
		upheap(i);
	}
}

#endif //#ifndef INCLUSION_HEAP_H
