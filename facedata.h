#ifndef INCLUSION_FACEDATA_H
#define INCLUSION_FACEDATA_H

class FaceData;
#define PI 3.1415926535

#include "quadedge.h"
#include "heap.h"
#include <cmath>

//A class to represent triangles in a subdivision
//Since we are using triangular decomposition, we can use this class
//to represent all faces.
class FaceData : public heap_node<FaceData*> {
private:
	//The edge which the face is anchored off of
	Edge * anchor;
	//The best point in this face to insert
	Point3d * bestPoint;
public:
	FaceData();
	FaceData(Edge *e);

	//Manipulation functions
	void setBest(double x, double y, double z);
	void setBest(Point3d * p);

	//Querying functions
	void longEdgeBis(PointData * p);
	void longEdgeBis(double * x, double * y);
	Edge * longEdge();
	bool isInterior(Point2d & p);
	void getPlane(double * a, double * b, double * c);
	void getPointsOrdY(double * ordx, double * ordy);
	void getPointsOrdY(int * ordx, int * ordy);
	void reAnchor(Edge * e);
	Edge * getAnchor();
	double area();
	double vol();
	double dihedralAngle();
	double getRho();
	void getCircumcenter(double * x, double * y);
	void getBest(Point3d * p);
	void printPoints();
	bool inSortedPoints(std::vector<Point3d *> * points);
	~FaceData();
};

inline FaceData::FaceData() {
	bestPoint = new Point3d();
}

inline FaceData::FaceData(Edge *e) {
	anchor = e;
	bestPoint = new Point3d();
}

inline FaceData::~FaceData() {
		delete bestPoint;
}

//Return the area of this face
inline double FaceData::area() {
	return triArea(this->anchor->lPrev()->org2d(), this->anchor->org2d(), this->anchor->dest2d()) / 2;
}

//Return the area under this face, referenced to z = 0
//Luck for us, this value is just the average of the zs times the projected x-y area.
inline double FaceData::vol() {
	return (this->anchor->lPrev()->org2d().getZ() + this->anchor->org2d().getZ() + this->anchor->dest2d().getZ()) * this->area() / 3;
}

//Determine whether a point is within this face (in x-y space)
inline bool FaceData::isInterior(Point2d & p) {
	return (!rightOf(p, anchor)) && (!rightOf(p, anchor->lNext())) && (!rightOf(p, anchor->lPrev()));
}

//Set the best point for this face
inline void FaceData::setBest(double x, double y, double z) {
	bestPoint->setX(x);
	bestPoint->setY(y);
	bestPoint->setZ(z);
}

//Set the best point for this face
inline void FaceData::setBest(Point3d * p) {
	bestPoint->set(p);
}

//Return the best point for this face
//Note that we don't return the pointer, so as not to give
//access to the internal data of the face
inline void FaceData::getBest(Point3d * p) {
	p->set(bestPoint);
}

//Return the anchor of this face
inline Edge* FaceData::getAnchor() {
	return anchor;
}

//reanchor the face
inline void FaceData::reAnchor(Edge * e) {
	anchor = e;
}

//Returns the magnitude of the dihedral angle between this face and the xy plane
//The dihedral angle is just the angle between the two normal vectors, which in
//case is the angle between (a, b, 1) and (0, 0, 1)
inline double FaceData::dihedralAngle() {
	double a = 0, b = 0, c = 0;
	getPlane(&a, &b, &c);
	return acos(pow(1 / (a * a + b * b + 1), 0.5)) * 180 / PI;
}

//Assuming points is sorted by y-value, return true iff all three points in this face are in the vector
inline bool FaceData::inSortedPoints(std::vector<Point3d*> *points) {
	return points && pointInSortedVector((PointData *) this->getAnchor()->org(), (std::vector<Point2d*>*) points)
		&& pointInSortedVector((PointData *) this->getAnchor()->dest(), (std::vector<Point2d*>*) points)
		&& pointInSortedVector((PointData *) this->getAnchor()->lPrev()->org(), (std::vector<Point2d*>*) points);
}

#endif
