#ifndef INCLUSION_POINTDATA_H
#define INCLUSION_POINTDATA_H

#include "geom2d.h"

//A class to represent points in a QuadEdge
class PointData : public Point3d {
private:
	//The index of this point, used in writing out .obj files
	unsigned int index;

	//static int numPointData;
	//static PointData * block;
public:
	//Constructors
	PointData()
	{this->x = 0; this->y = 0; this->z = 0; this->index = 0;}
	PointData(double x, double y)
		{ this->x = x; this->y = y; this->z = 0; this->index = 0;}
	PointData(double x, double y, double z, int index)
		{ this->x = x; this->y = y; this->z = z; this->index = index;}
	PointData(const PointData& p)	{ *this = p; }

	//Querying functions
	int getIndex() const;

	//set the index
	void setIndex(const int index);

	//Mathematical manipulation functions
	Vector2d operator-(const Point2d&) const;
	int operator==(const Point2d&) const;

	//Output function
	friend ostream& operator<<(ostream&, const PointData&);

	//The function from EdgeDataA
	void funct() {};

	//void * operator new(size_t bytes);
};

//Return the index of the point
inline int PointData::getIndex() const {
	return index;
}

//Return the index of the point
inline void PointData::setIndex(int index) {
	this->index = index;
}

//Return the vector from p to this point
inline Vector2d PointData::operator-(const Point2d& p) const
{
	return Vector2d(x - p.getX(), y - p.getY());
}

//Returns true iff p is within EPS of this point in the x-y plane
inline int PointData::operator==(const Point2d& p) const
{
	return ((*this - p).norm() < EPS);
}

//Print out the two dimentional coordinate of p
inline ostream& operator<<(ostream& os, const PointData& p)
{
	os << "(" << p.x << " " << p.y << ")";
	return os;
}
#endif //#ifndef INCLUSION_QUADDATA_H
