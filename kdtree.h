#ifndef INCLUSION_KDTREE_H
#define INCLUSION_KDTREE_H

#include <vector>
#include "geom2d.h"
using std::vector;

class KDTree {
private:
	static const int LEFT = 0;
	static const int RIGHT = 1;

	//The number of axes in the kd tree
	static const int AXES = 2;

	static const int X_AXIS = 0;
	static const int Y_AXIS = 1;

	static const double INFINITY;
	KDTree ** subTrees;
	Point2d * point;

	KDTree(std::vector<Point2d *>::iterator start, std::vector<Point2d *>::iterator end, int depth);
	
	void init(std::vector<Point2d *>::iterator start, std::vector<Point2d *>::iterator end, int depth);
	void nearestNeighbor(Point2d * p, Point2d ** best, double * bestDist, int depth, Point2d * ll, Point2d * ur);
	void findMin(double * best, Point2d ** bestPoint, int axis, int depth);
	void findMax(double * best, Point2d ** bestPoint, int axis, int depth);
	bool remove(Point2d * p, int depth);
	KDTree * findPoint(Point2d * p, int * depth);
	KDTree * findParent(Point2d * p, int * depth);
	double getAxis(Point2d * p, int axis);
	void getPointsInBox(double llx, double lly, double urx, double ury, int depth, vector<Point2d *> * vect);
	void print(int * points, int depth);

public:
	KDTree(vector<Point2d *> * points);
	Point2d * nearestNeighbor(Point2d * p);
	bool remove(Point2d * p);
	void getPointsInBox(double llx, double lly, double urx, double ury, vector<Point2d *> * vect);
	~KDTree();
};

void nnSort(vector<Point2d*> * vect);
#endif