#include "quadtree.h"
#include <stdio.h>
#include <stdlib.h>

QuadTree::QuadTree(double botLeftx, double botLefty, double width, double height) : 
botLeftx(botLeftx), botLefty(botLefty), width(width), height(height) {

	subdevs = (QuadTree **) malloc(sizeof(QuadTree *) * 4);
	subdevs[0] = NULL;
	subdevs[1] = NULL;
	subdevs[2] = NULL;
	subdevs[3] = NULL;
	//We have no mass, so put the CoM anywhere
	point = NULL;
	type = EMPTY;
}

QuadTree::QuadTree(double botLeftx, double botLefty, double width, double height, 
				   Point2d * point) : 
botLeftx(botLeftx), botLefty(botLefty), width(width), height(height) {

	subdevs = (QuadTree **) malloc(sizeof(QuadTree *) * 4);
	subdevs[0] = NULL;
	subdevs[1] = NULL;
	subdevs[2] = NULL;
	subdevs[3] = NULL;

	this->point = point;
	type = POINT;
}

QuadTree::~QuadTree() {
	for (int i = 0; i < 4; i++) {
		if (subdevs[i] != NULL) {
			delete(subdevs[i]);
		}
	}

	free(subdevs);
}


QuadTree * QuadTree::get00() {
	return subdevs[BOTLEFT];
}

QuadTree * QuadTree::get01() {
	return subdevs[TOPLEFT];
}

QuadTree * QuadTree::get10() {
	return subdevs[BOTRIGHT];
}

QuadTree * QuadTree::get11() {
	return subdevs[TOPRIGHT];
}

int QuadTree::findSubdev(Point2d * point) {
	if ((point->getX() - botLeftx) < width / 2 && 
		(point->getY() - botLefty) < height / 2) {
			return BOTLEFT;
	} else if ((point->getX() - botLeftx) < width / 2) {
		return TOPLEFT;
	} else if ((point->getY() - botLefty) < height / 2) {
		return BOTRIGHT;
	} else {
		return TOPRIGHT;
	}
}

double QuadTree::getsubx(int subdev) {
	if (subdev == BOTLEFT || subdev == TOPLEFT) {
		return botLeftx;
	} else {
		return botLeftx + width / 2;
	}
}

double QuadTree::getsuby(int subdev) {
	if (subdev == BOTLEFT || subdev == BOTRIGHT) {
		return botLefty;
	} else {
		return botLefty + height / 2;
	}
}

void QuadTree::addPoint(Point2d * newPoint) {
	if (type == EMPTY) {
		//The node is empty, put the point in it
		point = newPoint;
		type = POINT;
	} else if (type == POINT) {
		//There is a point here, make a subdevision for the existing point
		//and update the center of mass
		type = NODE;
		//Reinsert all points, assuming this is a node
		this->addPoint(point);
		this->addPoint(newPoint);
		point = NULL;
	} else {
		int subdev = findSubdev(newPoint);
		if (subdevs[subdev] == NULL) {
			subdevs[subdev] = new QuadTree(getsubx(subdev), getsuby(subdev),
				width / 2, height / 2, newPoint);
		} else {
			subdevs[subdev]->addPoint(newPoint);
		}
	}
}

Point2d * QuadTree::getNear(Point2d * point) {
	int subdev = findSubdev(point);

	if (type == EMPTY) {
		return NULL;
	} else if (type == POINT) {
		return this->point;
	} else {
		//Traverse the tree, trying to stay as close to the point as possible
		for (int i = 0; i < 4; i++) {
			if (subdevs[(subdev + i)%4]) {
				return subdevs[(subdev + i)%4]->getNear(point);
			}
		}
	}

	//We never get here
	return NULL;
}

void QuadTree::zOrder(vector<Point2d*> * vect) {
	switch (type) {
		case POINT:
			//Put this point on the back
			vect->push_back(point);
			break;
		case NODE:
			//Insert all the points from the subtrees, in order botleft, topleft, botright, topright
			for (int i = 0; i < 4; i++) {
				if (subdevs[i] != NULL) {
					subdevs[i]->zOrder(vect);
				}
			}
			break;
		default:
			break;
	}
}
