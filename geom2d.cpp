#include "geom2d.h"
#include <algorithm>

using namespace std;

// Returns true iff the point d is inside the circle defined by the
// points a, b, c. See Guibas and Stolfi (1985) p.107.
bool inCircle(const Point2d& a, const Point2d& b,
			  const Point2d& c, const Point2d& d)
{
	//This is the equivalent of:
	//|a|^2 * triArea(b, c, d) - |b|^2 * triArea(a, c, d) + 
	//|c|^2 * triArea(a, b, d) - |d|^2 * triArea(a, b, c) > 0
	//Only we have shifted the coordinates so that a is at (0, 0), which simplifies
	//things a great deal, but gives the same result
	double ax = a.getX();
	double ay = a.getY();
	double bx = b.getX() - ax;
	double by = b.getY() - ay;
	double cx = c.getX() - ax;
	double cy = c.getY() - ay;
	double dx = d.getX() - ax;
	double dy = d.getY() - ay;
	return -(bx*bx + by*by) * (cx*dy - cy*dx) +
		(cx*cx + cy*cy) * (bx*dy - by*dx) -
		(dx*dx + dy*dy) * (bx*cy - by*cx) > 0;
}

//Check if a circle intersects a rectangle
//Straight out of Graphics Gems  by Glassner, p52, except an intersect with 0 area counts
//here
bool rectCircInersect(Point2d * ll, Point2d * ur, Point2d * center, double radius) {
	//Center the circle at (0,0)
	double llx = ll->getX() - center->getX();
	double lly = ll->getY() - center->getY();
	double urx = ur->getX() - center->getX();
	double ury = ur->getY() - center->getY();

	if (urx < 0) {
		if (ury < 0) {
			//The box is entirely below and left of the circle
			//The closest point is the upper right corner
			return urx*urx + ury*ury <= radius*radius;
		} else if (lly > 0) {
			//The box is entirely above and left of the circle
			//The closest point is the lower right corner
			return urx*urx + lly*lly <= radius*radius;
		} else {
			//The box is left of the circle, but crosses y=0
			//Closest approach is where the circle reaches to (-radius, 0)
			return -urx <= radius;
		}
	} else if (llx > 0) {
		//The same as above, only we are now to the right of the circle
		if (ury < 0) {
			return llx*llx + ury*ury <= radius*radius;
		} else if (lly > 0) {
			return llx*llx + lly*lly <= radius*radius;
		} else {
			return llx <= radius;
		}
	} else {
		//The rectangle intersects the y axis
		if (ury < 0) {
			//We are below the circle
			return -ury <= radius;
		} else if (lly > 0) {
			//We are above the cirlce
			return lly <= radius;
		} else {
			//We contain (0, 0)
			return true;
		}
	}
}

//Assuming sortedyVect is sorted by y value, return true iff p is in sortedyVect
bool pointInSortedVector(Point2d * p, vector<Point2d *> * sortedyVect) {
	struct LessThanY {
		bool operator()(Point2d * p1, Point2d * p2) {
			return p1->getY() < p2->getY();
		}
	} lty;

	vector<Point2d *>::iterator curPoint;
	curPoint = lower_bound(sortedyVect->begin(), sortedyVect->end(), p, lty);
	//if (curPoint == sortedyVect->end())
	//	return false;
	//else
	//	curPoint++;

	for (;curPoint != sortedyVect->end() && (*curPoint)->getY() == p->getY(); curPoint++) {
		if ((*curPoint)->getX() == p->getX())
			return true;
	}

	return false;
}
