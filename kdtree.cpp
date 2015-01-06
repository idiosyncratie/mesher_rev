#include "kdtree.h"
#include <algorithm>
#include <limits>

using namespace std;
const double KDTree::INFINITY = numeric_limits<double>::infinity();

//A function to sort points by x value
inline bool ltx(Point2d * p1, Point2d * p2) {
	return p1->getX() < p2->getX();
}

//A function to sort points by y value
inline bool lty(Point2d * p1, Point2d * p2) {
	return p1->getY() < p2->getY();
}

KDTree::KDTree(std::vector<Point2d *> * points) {
	init(points->begin(), points->end(), 0);
}

KDTree::KDTree(std::vector<Point2d *>::iterator start, std::vector<Point2d *>::iterator end, int depth) {
	init(start, end, depth);
}

//Find the median, and make a kd-tree out of the left and right sides of the median. This insures the tree is
//balenced
void KDTree::init(std::vector<Point2d *>::iterator start, std::vector<Point2d *>::iterator end, int depth) {
	subTrees = new KDTree*[2];

	//Sort the points
	std::vector<Point2d *>::iterator middle = start + ((end - start) / 2);

	switch (depth % AXES) {
		case X_AXIS:
			nth_element(start, middle, end, ltx);
			break;
		case Y_AXIS:
			nth_element(start, middle, end, lty);
			break;
	}

	point = *middle;

	if (start == middle) {
		subTrees[LEFT] = NULL;
	} else {
		subTrees[LEFT] = new KDTree(start, middle, depth + 1);
	}

	if (end == middle + 1) {
		subTrees[RIGHT] = NULL;
	} else {
		subTrees[RIGHT] = new KDTree(middle + 1, end, depth + 1);
	}

}

//Return the nearest neighbor to p
Point2d * KDTree::nearestNeighbor(Point2d * p) {
	double dist = INFINITY;
	Point2d start(*point);
	Point2d * pStart = &start;
	Point2d ** best = &pStart;
	Point2d * ur = new Point2d(INFINITY, INFINITY);
	Point2d * ll = new Point2d(-INFINITY, -INFINITY);
	nearestNeighbor(p, best, &dist, 0, ll, ur);

	delete ur;
	delete ll;
	return *best;
}

//Return the nearest neighbor to p
void KDTree::nearestNeighbor(Point2d * p, Point2d ** best, double * bestDist, int depth, Point2d * ll, Point2d * ur) {
	double distance = point->distance(*p);
	int nextSubdev;
	int otherSubdev;
	Point2d * newll = new Point2d(*ll);
	Point2d * newur = new Point2d(*ur);

	if (distance < *bestDist) {
		*best = point;
		*bestDist = distance;
	}

	switch (depth % AXES) {
		case X_AXIS:
			if (p->getX() < point->getX()) {
				nextSubdev = LEFT;
				otherSubdev = RIGHT;
				//We know all subpoints are to the left of our point
				newur->setX(point->getX());
			} else {
				nextSubdev = RIGHT;
				otherSubdev = LEFT;
				//We know all subpoints are to the right of our point
				newll->setX(point->getX());
			}
			break;
		case Y_AXIS:
			if (p->getY() < point->getY()) {
				otherSubdev = RIGHT;
				nextSubdev = LEFT;
				newur->setY(point->getY());
			} else {
				nextSubdev = RIGHT;
				otherSubdev = LEFT;
				newll->setY(point->getY());
			}
			break;
	}

	if (subTrees[nextSubdev]) {
		subTrees[nextSubdev]->nearestNeighbor(p, best, bestDist, depth + 1, newll, newur);
	}

	//Check if it is possible for the other tree to contain a point better than the one we have
	if (subTrees[otherSubdev]) {
		//Switch which way we have changed the bounding box, as we are on the other side of the tree
		switch (depth % AXES) {
		case X_AXIS:
			if (nextSubdev == LEFT) {
				newll->setX(point->getX());
				newur->setX(ur->getX());
			} else {
				newll->setX(ll->getX());
				newur->setX(point->getX());
			}
			break;
		case Y_AXIS:
			if (nextSubdev == LEFT) {
				newll->setY(point->getY());
				newur->setY(ur->getY());
			} else {
				newll->setY(ll->getY());
				newur->setY(point->getY());
			}
			break;
		}

		//Only search this tree if it's bounding rectangle intersects the circle
		//of bestDist around our target point
		if (rectCircInersect(newll, newur, p, *bestDist)) {
			subTrees[otherSubdev]->nearestNeighbor(p, best, bestDist, depth + 1, newll, newur);
		}
	}

	delete newll;
	delete newur;
}

//Remove the given point, if it is found
bool KDTree::remove(Point2d * p) {
	return remove(p, 0);
}

//Find a point in the tree
KDTree * KDTree::findPoint(Point2d * p, int * depth) {
	if (p == point)
		return this;

	int thisDepth = *depth;
	KDTree * tree;
	switch (thisDepth % AXES) {
		case X_AXIS:
			if (p->getX() == point->getX()) {
				//The point could either be on the left or the right
				(*depth) = thisDepth + 1;
				if (subTrees[LEFT] && (tree = subTrees[LEFT]->findPoint(p, depth)))
					return tree;
				else if (subTrees[RIGHT]) {
					(*depth) = thisDepth + 1;
					return subTrees[RIGHT]->findPoint(p, depth);
				} else
					return NULL;
			} else if (p->getX() < point->getX()) {
				(*depth) = thisDepth + 1;
				if(subTrees[LEFT])
					return subTrees[LEFT]->findPoint(p, depth);
				else
					return NULL;
			} else {
				(*depth) = thisDepth + 1;
				if(subTrees[RIGHT])
					return subTrees[RIGHT]->findPoint(p, depth);
				else
					return NULL;
			}
			break;
		case Y_AXIS:
			if (p->getY() == point->getY()) {
				//The point could either be on the left or the right
				(*depth) = thisDepth + 1;
				if (subTrees[LEFT] && (tree = subTrees[LEFT]->findPoint(p, depth)))
					return tree;
				else if (subTrees[RIGHT]) {
					(*depth) = thisDepth + 1;
					return subTrees[RIGHT]->findPoint(p, depth);
				} else
					return NULL;
			} else if (p->getY() < point->getY()) {
				(*depth) = thisDepth + 1;
				if(subTrees[LEFT])
					return subTrees[LEFT]->findPoint(p, depth);
				else
					return NULL;
			} else {
				(*depth) = thisDepth + 1;
				if(subTrees[RIGHT])
					return subTrees[RIGHT]->findPoint(p, depth);
				else
					return NULL;
			}
			break;
		default:
			return NULL;
	}
}

//Find the parent of a point in the tree
KDTree * KDTree::findParent(Point2d * p, int * depth) {
	if ((subTrees[LEFT] && p == subTrees[LEFT]->point) || (subTrees[RIGHT] && p == subTrees[RIGHT]->point))
		return this;

	int thisDepth = *depth;
	KDTree * tree;
	switch (thisDepth % AXES) {
		case X_AXIS:
			if (p->getX() == point->getX()) {
				//The point could either be on the left or the right
				(*depth) = thisDepth + 1;
				if (subTrees[LEFT] && (tree = subTrees[LEFT]->findParent(p, depth)))
					return tree;
				else if (subTrees[RIGHT]) {
					(*depth) = thisDepth + 1;
					return subTrees[RIGHT]->findParent(p, depth);
				} else
					return NULL;
			} else if (p->getX() < point->getX()) {
				(*depth) = thisDepth + 1;
				if(subTrees[LEFT])
					return subTrees[LEFT]->findParent(p, depth);
				else
					return NULL;
			} else {
				(*depth) = thisDepth + 1;
				if(subTrees[RIGHT])
					return subTrees[RIGHT]->findParent(p, depth);
				else
					return NULL;
			}
			break;
		case Y_AXIS:
			if (p->getY() == point->getY()) {
				//The point could either be on the left or the right
				(*depth) = thisDepth + 1;
				if (subTrees[LEFT] && (tree = subTrees[LEFT]->findParent(p, depth)))
					return tree;
				else if (subTrees[RIGHT]) {
					(*depth) = thisDepth + 1;
					return subTrees[RIGHT]->findParent(p, depth);
				} else
					return NULL;
			} else if (p->getY() < point->getY()) {
				(*depth) = thisDepth + 1;
				if(subTrees[LEFT])
					return subTrees[LEFT]->findParent(p, depth);
				else
					return NULL;
			} else {
				(*depth) = thisDepth + 1;
				if(subTrees[RIGHT])
					return subTrees[RIGHT]->findParent(p, depth);
				else
					return NULL;
			}
			break;
		default:
			return NULL;
	}
}

//Get the split axis based on the depth
double KDTree::getAxis(Point2d * p, int axis) {
	switch (axis) {
		case X_AXIS:
			return p->getX();
		default:
			return p->getY();
	}
}

//Remove a point, if it is found
bool KDTree::remove(Point2d * p, int depth) {
	int axis = depth % AXES;
	if (p != point) {
		int depthCopy = depth;
		KDTree * tree = findPoint(p, &depthCopy);
		if (tree) {
			tree->remove(p, depthCopy);
			if (tree->point == p) {
				KDTree * parent = findParent(p, &depth);
				if (tree == parent->subTrees[LEFT])
					parent->subTrees[LEFT] = NULL;
				else 
					parent->subTrees[RIGHT] = NULL;


				delete tree;
			}
			return true;
		} else {
			return false;
		}
	}
	//We found the point, now for the real work

	//Find the rightmost or topmost point in the left subtree, put it here
	//and remove it from it's old location
	double lbest = -INFINITY;
	double rbest = INFINITY;
	Point2d * a = NULL;
	Point2d * b = NULL;
	Point2d ** lp = &a;
	Point2d ** rp = &b;

	if (subTrees[LEFT])
		subTrees[LEFT]->findMax(&lbest, lp, axis, depth + 1);

	if (subTrees[RIGHT])
		subTrees[RIGHT]->findMin(&rbest, rp, axis, depth + 1);

	if (*lp && getAxis(*lp, axis) == getAxis(p, axis)) {
		//There is a point on the left which has the same value as the point we
		//removing, use it
		point = *lp;
		subTrees[LEFT]->remove(point, depth + 1);
		if (subTrees[LEFT]->point == point) {
			//The tree below us is a (now empty) leaf, delete it
			delete subTrees[LEFT];
			subTrees[LEFT] = NULL;
		}
	} else if (*rp && getAxis(*rp, axis) == getAxis(p, axis)) {
		point = *rp;
		//There is a point on the right which has the same value as the point we
		//removing, use it
		subTrees[RIGHT]->remove(point, depth + 1);
		if (subTrees[RIGHT]->point == point) {
			//The tree below us is a (now empty) leaf, delete it
			delete subTrees[RIGHT];
			subTrees[RIGHT] = NULL;
		}
	} else if (*lp) {
		//There is not a point with the same value as the point we are removing,
		//so use the left point
		point = *lp;
		subTrees[LEFT]->remove(point, depth + 1);
		if (subTrees[LEFT]->point == point) {
			//The tree below us is a (now empty) leaf, delete it
			delete subTrees[LEFT];
			subTrees[LEFT] = NULL;
		}
	} else if (*rp) {
		point = *rp;
		//There was no point to take from the left, so take from the right
		subTrees[RIGHT]->remove(point, depth + 1);
		if (subTrees[RIGHT]->point == point) {
			//The tree below us is a (now empty) leaf, delete it
			delete subTrees[RIGHT];
			subTrees[RIGHT] = NULL;
		}
	} else {
		//If we got here, we are removing from a leaf, so we are done
		return true;
	}

	return true;
}

//Find the point with maximum value with specified axis in the tree
void KDTree::findMax(double * best, Point2d ** bestPoint, int axis, int depth) {

	switch (axis) {
		case X_AXIS:
			if (point->getX() > *best) {
				*best = point->getX();
				*bestPoint = point;
			}
			break;
		case Y_AXIS:
			if (point->getY() > *best) {
				*best = point->getY();
				*bestPoint = point;
			}
			break;
	}

	if (axis != depth % AXES) {
		//We have no idea where the maximum point could be, so search both sides
		if (subTrees[LEFT])
			subTrees[LEFT]->findMax(best, bestPoint, axis, depth + 1);

		if (subTrees[RIGHT])
			subTrees[RIGHT]->findMax(best, bestPoint, axis, depth + 1);
	} else {
		//The best point can only be in the right tree if we are splitting along the axis of
		//interest.
		if (subTrees[RIGHT])
			subTrees[RIGHT]->findMax(best, bestPoint, axis, depth + 1);
	}
}

//Find the point with minimum value with specified axis in the tree
void KDTree::findMin(double * best, Point2d ** bestPoint, int axis, int depth) {
	switch (axis) {
		case X_AXIS:
			if (point->getX() < *best) {
				*best = point->getX();
				*bestPoint = point;
			}
			break;
		case Y_AXIS:
			if (point->getY() < *best) {
				*best = point->getY();
				*bestPoint = point;
			}
			break;
	}
	if (axis != depth % AXES) {
		//We have no idea where the maximum point could be, so search both sides
		if (subTrees[RIGHT])
			subTrees[RIGHT]->findMin(best, bestPoint, axis, depth + 1);

		if (subTrees[LEFT])
			subTrees[LEFT]->findMin(best, bestPoint, axis, depth + 1);
	} else {
		//The best point can only be in the left tree if we are splitting along the axis of
		//interest.
		if (subTrees[LEFT])
			subTrees[LEFT]->findMin(best, bestPoint, axis, depth + 1);
	}
}


//Do a nearest neighbor sort of the given points
void nnSort(vector<Point2d*> * vect) {
	KDTree tree(vect);
	int size = vect->size();
	Point2d * p = (*vect)[0];
	vect->clear();
	vect->push_back(p);
	int points = 0;
	
	tree.remove(p);
	double dist = 0;
	for (int i = 1; i < size; i++) {
		p = tree.nearestNeighbor(p);
		vect->push_back(p);
		dist += (*p - *(*vect)[i - 1]).norm();
		if (!tree.remove(p)) {
			cout << "error\n";
			exit(1);
		}
	}

	cout << dist << "\n";
}

//Return all points in the box
void KDTree::getPointsInBox(double llx, double lly, double urx, double ury, vector<Point2d *> * vect) {
	getPointsInBox(llx, lly, urx, ury, 0, vect);
}

//Return all points in the box
void KDTree::getPointsInBox(double llx, double lly, double urx, double ury, int depth, vector<Point2d *> * vect) {
	if (point->getX() >= llx && point->getY() >= lly && point->getX() <= urx && point->getY() <= ury)
		vect->push_back(point);

	switch (depth % AXES) {
		case X_AXIS:
			if (subTrees[LEFT] && point->getX() >= llx) {
				//points to the left of this one could be in the box
				subTrees[LEFT]->getPointsInBox(llx, lly, urx, ury, depth + 1, vect);
			}

			if (subTrees[RIGHT] && point->getX() <= urx) {
				//points to the right of this one could be in the box
				subTrees[RIGHT]->getPointsInBox(llx, lly, urx, ury, depth + 1, vect);
			}
			break;
		case Y_AXIS:
			if (subTrees[LEFT] && point->getY() >= lly) {
				//points to the left of this one could be in the box
				subTrees[LEFT]->getPointsInBox(llx, lly, urx, ury, depth + 1, vect);
			}

			if (subTrees[RIGHT] && point->getY() <= ury) {
				//points to the right of this one could be in the box
				subTrees[RIGHT]->getPointsInBox(llx, lly, urx, ury, depth + 1, vect);
			}
			break;
	}
}

KDTree::~KDTree() {
	if (subTrees[LEFT])
		delete subTrees[LEFT];

	if (subTrees[RIGHT])
		delete subTrees[RIGHT];

	delete subTrees;
}