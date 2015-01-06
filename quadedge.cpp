#include "quadedge.h"

QuadEdge::QuadEdge()
{
	e[0].num = 0, e[1].num = 1, e[2].num = 2, e[3].num = 3;
	e[0].next = &(e[0]); e[1].next = &(e[3]);
	e[2].next = &(e[2]); e[3].next = &(e[1]);
}

//Returns true iff x is within EPS of e
bool OnEdge(const Point2d& x, Edge* e) {
	double t1, t2, t3;
	t1 = (x -  e->org2d()).norm();
	t2 = (x -  e->dest2d()).norm();
	if (t1 < EPS || t2 < EPS)
		return true;
	t3 = ( e->org2d() -  e->dest2d()).norm();
	if (t1 > t3 || t2 > t3)
		return false;
	Line line( e->org2d(),  e->dest2d());
	return (fabs(line.eval(x)) < EPS);
}
