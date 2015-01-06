#ifndef GEOM2D_H
#define GEOM2D_H

#include <cmath>
#include <iostream>
#include <vector>

using std::istream;
using std::ostream;
using std::cerr;

#ifndef ABS
#define ABS(a)	((a) >= 0 ? (a) : -(a))
#endif

#ifndef MAX
#define MAX(a, b)	((a) >= (b) ? (a) : (b))
#define MIN(a, b)	((a) <= (b) ? (a) : (b))
#endif

#ifndef TRUE
#define FALSE 0
#define TRUE  1
#endif

#define EPS 1e-6

class Vector2d {
public:
	double x, y;
	Vector2d()					{ x = 0; y = 0; }
	Vector2d(double a, double b)	{ x = a; y = b; }
	double norm() const;
	void normalize();
	Vector2d operator+(const Vector2d&) const;
	Vector2d operator-(const Vector2d&) const;
	friend Vector2d operator*(double, const Vector2d&);
	friend double dot(const Vector2d&, const Vector2d&);
	friend istream& operator>>(istream&, Vector2d&);
	friend ostream& operator<<(ostream&, const Vector2d&);
};

//The standard implementation of a point
class Point2d {
protected:
	double x, y;
public:
	Point2d()					{ x = 0; y = 0; }
	Point2d(double a, double b)		{ x = a; y = b; }
	Point2d(const Point2d& p)	{ *this = p; }

	double getX() const;
	double getY() const;
	void setX(double xnew);
	void setY(double ynew);
	Point2d operator+(const Vector2d&) const;
	Vector2d operator-(const Point2d&) const;
	int operator==(const Point2d&) const;
	friend istream& operator>>(istream&, Point2d&);
	friend ostream& operator<<(ostream&, const Point2d&);
	double distance(const Point2d & p) const;
};

//Add a z value to a two dimentional point
class Point3d : public Point2d {
private:
protected:
	double z;
public:
	Point3d()			{x = 0; y = 0; z = 0;}
	Point3d(double a, double b, double c) {x = a; y = b; z = c;}
	Point3d(const Point3d& p) {*this = p;}
	double getZ() const {return z;}
	void setZ(double znew) {z = znew;}
	void set(double nx, double ny, double nz) {x = nx; y = ny; z = nz;}
	void set(Point3d * p) {x = p->x; y = p->y; z = p->z;}
	friend bool operator<(Point3d a, Point3d b);
	void move(double dx, double dy, double dz) {x += dx; y+=dy; z+=dz;}
};

//A 2d line
class Line {
public:
	Line()	{}
	Line(const Point2d&, const Point2d&);
	double eval(const Point2d&) const;
	int classify(const Point2d&) const;
private:
	double a, b, c;
};

// Vector2d:
inline double Vector2d::norm() const
{
	return sqrt(x * x + y * y);
}

inline void Vector2d::normalize()
{
	double len;

	if ((len = sqrt(x * x + y * y)) == 0.0)
		cerr << "Vector2d::normalize: Division by 0\n";
	else {
		x /= len;
		y /= len;
	}
}

inline Vector2d Vector2d::operator+(const Vector2d& v) const
{
	return Vector2d(x + v.x, y + v.y);
}

inline Vector2d Vector2d::operator-(const Vector2d& v) const
{
	return Vector2d(x - v.x, y - v.y);
}

inline Vector2d operator*(double c, const Vector2d& v)
{
	return Vector2d(c * v.x, c * v.y);
}

inline double dot(const Vector2d& u, const Vector2d& v)
{
	return u.x * v.x + u.y * v.y;
}

inline ostream& operator<<(ostream& os, const Vector2d& v)
{
	os << '(' << v.x << ", " << v.y << ')';
	return os;
}

inline istream& operator>>(istream& is, Vector2d& v)
{
	is >> v.x >> v.y;
	return is;
}

// Point2d:

inline double Point2d::getX() const {
	return x;
}

inline double Point2d::getY() const {
	return y;
}

inline void Point2d::setX(double xnew) {
	x = xnew;
}

inline void Point2d::setY(double ynew) {
	y = ynew;
}

inline double Point2d::distance(const Point2d & p) const {
	return pow(pow((x - p.x), 2) + pow((y - p.y), 2), 0.5);
}

inline Point2d Point2d::operator+(const Vector2d& v) const
{
	return Point2d(x + v.x, y + v.y);
}

inline Vector2d Point2d::operator-(const Point2d& p) const
{
	return Vector2d(x - p.getX(), y - p.getY());
}

inline int Point2d::operator==(const Point2d& p) const
{
	return ((*this - p).norm() < EPS);
}

inline istream& operator>>(istream& is, Point2d& p)
{
	is >> p.x >> p.y;
	return is;
}

inline ostream& operator<<(ostream& os, const Point2d& p)
{
	os << '(' << p.x << ", " << p.y << ')';
	return os;
}

// Line:

inline Line::Line(const Point2d& p, const Point2d& q)
// Computes the normalized line equation through the
// points p and q.
{
	Vector2d t = q - p;
	double len = t.norm();
	a =   t.y / len;
	b = - t.x / len;
	c = -(a*p.getX() + b*p.getY());
}

inline double Line::eval(const Point2d& p) const
// Plugs point p into the line equation.
{
	return (a * p.getX() + b* p.getY() + c);
}

inline int Line::classify(const Point2d& p) const
// Returns -1, 0, or 1, if p is to the left of, on,
// or right of the line, respectively.
{
	double d = eval(p);
	return (d < -EPS) ? -1 : (d > EPS ? 1 : 0);
}

inline bool operator<(Point3d a, Point3d b) {
	return a.y < b.y;
}

// Returns twice the area of the oriented triangle (a, b, c), i.e., the
// area is positive if the triangle is oriented counterclockwise.
inline double triArea(const Point2d& a, const Point2d& b, const Point2d& c)
{
	double ax = a.getX();
	double ay = a.getY();
	return (b.getX() - ax)*(c.getY() - ay) - (b.getY() - ay)*(c.getX() - ax);
}

// Returns twice the area of the oriented triangle (a, b, c), i.e., the
// area is positive if the triangle is oriented counterclockwise.
inline double triArea(const double ax, const double ay, 
					  const double bx, const double by, 
					  const double cx, const double cy) {
						  return (bx - ax)*(cy - ay) - (by - ay)*(cx - ax);
}

// Returns true iff the points a, b, c are in a counterclockwise order
inline bool ccw(const Point2d& a, const Point2d& b, const Point2d& c)
{
	return (triArea(a, b, c) > 0);
}

// Returns true iff the points a, b, c are in a clockwise order
inline bool cw(const Point2d& a, const Point2d& b, const Point2d& c) {
	return (triArea(a, b, c) < 0);
}


bool inCircle(const Point2d& a, const Point2d& b,
			  const Point2d& c, const Point2d& d);

bool rectCircInersect(Point2d * ll, Point2d * ur, Point2d * center, double radius);
bool pointInSortedVector(Point2d * p, std::vector<Point2d *> * sortedyVect);
#endif
