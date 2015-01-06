#ifndef QUADEDGE_H
#define QUADEDGE_H

#include "pointdata.h"

class QuadEdge;
class Edge;

bool OnEdge(const Point2d& x, Edge* e);

class Edge {
	friend class QuadEdge;
private:
	Edge *next;
	void* data;
	unsigned int num;
public:
	Edge()			{ data = NULL; }
	Edge* rot();
	Edge* invRot();
	Edge* sym();
	Edge* oNext();
	Edge* oPrev();
	Edge* dNext();
	Edge* dPrev();
	Edge* lNext();
	Edge* lPrev();
	Edge* rNext();
	Edge* rPrev();
	void* org();
	void* dest();
	void* lDual();
	void* rDual();
	const PointData& org2d() const;
	const PointData& dest2d() const;
	void  endPoints(void*, void*);
	void setNext(Edge * next);
	QuadEdge* qEdge()	{ return (QuadEdge *)(this - num); }
	friend ostream& operator<<(ostream&, const Point2d&);
};

class QuadEdge {
private:
	Edge e[4];

public:
	QuadEdge();
	Edge * first();
};

//Return the first edge in the quadEdge
inline Edge * QuadEdge::first() {
	return e;
}
/************************* Edge Algebra *************************************/

// Return the dual of the current edge, directed from its right to its left.
inline Edge* Edge::rot()
{
	return (num < 3) ? this + 1 : this - 3;
}

// Return the dual of the current edge, directed from its left to its right.
inline Edge* Edge::invRot()
{
	return (num > 0) ? this - 1 : this + 3;
}

// Return the edge from the destination to the origin of the current edge.
inline Edge* Edge::sym()
{
	return (num < 2) ? this + 2 : this - 2;
}

// Return the next ccw edge around (from) the origin of the current edge.
inline Edge* Edge::oNext()
{
	return next;
}

// Return the next cw edge around (from) the origin of the current edge.
inline Edge* Edge::oPrev()
{
	return rot()->oNext()->rot();
}

// Return the next ccw edge around (into) the destination of the current edge.
inline Edge* Edge::dNext()
{
	return sym()->oNext()->sym();
}

// Return the next cw edge around (into) the destination of the current edge.
inline Edge* Edge::dPrev()
{
	return invRot()->oNext()->invRot();
}

// Return the ccw edge around the left face following the current edge.
inline Edge* Edge::lNext()
{
	return invRot()->oNext()->rot();
}

// Return the ccw edge around the left face before the current edge.
inline Edge* Edge::lPrev()
{
	return oNext()->sym();
}

// Return the edge around the right face ccw following the current edge.
inline Edge* Edge::rNext()
{
	return rot()->oNext()->invRot();
}

// Return the edge around the right face ccw before the current edge.
inline Edge* Edge::rPrev()
{
	return sym()->oNext();
}

inline ostream& operator<<(ostream& os, const Edge& e)
{
	os << '(' << e.org2d() << ", " << e.dest2d() << ')';
	return os;
}

/************** Access to data pointers *************************************/

//The origin of this edge
inline void* Edge::org()
{
	return data;
}

//The destination of this edge
inline void* Edge::dest()
{
	return sym()->data;
}

//The left dual point (usually the left face)
inline void* Edge::lDual() {
	return invRot()->data;
}

//The right dual point (usually the right face)
inline void* Edge::rDual() {
	return rot()->data;
}

inline const PointData& Edge::org2d() const
{
	return *(PointData*)data;
}

inline const PointData& Edge::dest2d() const
{
	return *(PointData*)((num < 2) ? ((this + 2)->data) : ((this - 2)->data));
}

inline void Edge::endPoints(void* org, void* de)
{
	data = org;
	sym()->data = de;
}

inline void Edge::setNext(Edge * n) {
	next = n;
}

//Returns true iff x is on the line defined by e
inline bool onLine(const Point2d& x, Edge* e) {
	return fabs(triArea(x, e->dest2d(), e->org2d())) < EPS;
}

//Returns true iff x is right of the directed edge e
inline bool rightOf(const Point2d& x, Edge* e)
{
	return ccw(x,  e->dest2d(),  e->org2d());
}

//Returns true iff x is left of the directed edge e
inline bool leftOf(const Point2d& x, Edge* e)
{
	return ccw(x,  e->org2d(),  e->dest2d());
}
#endif /* QUADEDGE_H */
