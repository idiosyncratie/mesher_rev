#ifndef INCLUSION_SUBDIVISION_H
#define INCLUSION_SUBDIVISION_H

#include "facedata.h"
#include "heap.h"
#include "rawField.h"
#include "quadedge.h"
#include "blockalloc.h"
#include <cstdlib>
#include <fstream>
#include <list>
#include <set>

using std::set;
using std::list;
using std::ofstream;
using std::cout;

#undef INFINITY

class Subdivision {
	friend class Edge;
	friend class QuadEdge;
public:

	enum inputtype_t {
		RASTER,
		XYZ_FIELD,
		PROGRAM
	};

	//Positive numbers indicate global reduction, negative represent some sort of local reduction
	enum inserttype_t {
		CICRCUMCENTER = -1,
		CICRCUMCENTER_VOL = -2,
		MAX_ERR = 1,
		MAX_ERR_VOL = 2,
		ALL = 0
	};

private:

	//A set of edited faces
	set<FaceData*> * editedFaces;
	//A list of all points in the mesh
	vector<PointData*> * points;
	//A heap of all the faces, sorted by error
	Heap<FaceData*> * faceHeap;
	//A list of all faces in the mesh
	vector<FaceData*> * faces;
	//A block allocation of all the edges
	BlockAlloc<QuadEdge> * edges;
	BlockAlloc<FaceData> * allocFaces;

	//The four corners of our mesh
	PointData * boundingPoints[4];

	//Utility pointers, so we don't have to create and destroy every time
	PointData * utilPd;
	Point3d * utilP3d;

	//Private constructors
	Subdivision(const PointData&, const PointData&, 
		const PointData&, const PointData&);

	Subdivision(const PointData&, const PointData&, 
		const PointData&, const PointData&, const char * filename);
	//Basic numerical properties
	int numPoints;
	int numFaces;

	//The scale to multiply x and y by
	double scale;

	//Point location functoins
	Edge* locate(Point2d&, Edge* e = NULL);

	//A face to represent deleted faces
	static FaceData * DELETED;

	//Functions to insert points, assuming the point data actually makes
	//sense in this mesh (ie. it's index is correct)
	bool insertSite(PointData*);
	bool insertSite(PointData*, Edge *, bool knownEdge = false);

	//the default edge to search from when locating points
	Edge *startingEdge;

	//The output stream of this mesh
	ofstream out;
	char outputFile[500];

	//Low level mesh manipulation functions
	FaceData * makeFace(Edge* edge);
	void deleteEdge(Edge* e);
	Edge* Connect(Edge* a, Edge* b);
	void swap(Edge * e);
	void swapEdges(Edge * e, PointData * x, Edge * end);
	bool deleteExternal(PointData *);
	void deleteFace(FaceData * face);
	void updateFaceWithoutScale(FaceData * f);
	void updateFaceWithScale(FaceData * f);
	void updateFace(FaceData * f);
	void updateHeap();
	void updateError(FaceData * f, Point3d * p, double a, double b, double c, double * maxErr);
	void updateError(FaceData * f, double x, double y, double z, double a, double b, double c, double * maxErr);
	void updateErrorWithScale(FaceData * f, Point3d * p, double a, double b, double c, double * maxErr);
	void updateErrorWithScale(FaceData * f, double x, double y, double z, double a, double b, double c, double * maxErr);
	void splice(Edge* a, Edge* b);
	Edge* makeEdge();
	void markEdited(FaceData * f);
	void reAnchorFace(FaceData * f, Edge * e);
	Edge* guessEdge(Point2d & p);
	bool isValid(FaceData * f);
	bool isValid(Edge * e);
	//The type of rawfield we are processing (raster vs. xyz field)
	inputtype_t intype;
	//The type of insertions we are using
	inserttype_t inserttype;
	bool fixSliver(FaceData * f, bool force = false);

	//The x values of the left and right, and the y values of the top and bottom
	//Note that there is not a demand for the subdivision to be a rectange, so this is not
	//necessarily useful. Currently, it is only used to apply the importance map correctly
	double leftEdge, rightEdge, topEdge, bottomEdge;
	double width, height;
public:
	//The raw field of points we are reducing
	RawField * raw;
	//The importance map
	//Allows for certain areas to be specified as more important than others
	RawField * impMap;

	static double INFINITY;

	//Constructor
	Subdivision(RawField * field, const double scale, const char * filename, inserttype_t itype, RawField * impMap);
	Subdivision(const Point2d &, const Point2d &, const Point2d &, const Point2d &, const char * filename);

	//The initializer function
	void init(const PointData&, const PointData&, 
		const PointData&, const PointData&);

	//Point insertion functions
	bool insertSite(Point3d & p);
	bool insertSite(const double x, const double y, const double z);
	bool insertSite(FaceData * f);
	bool insertSite(const double x, const double y, Edge * guess = NULL, bool force = false, bool knownEdge = false);
	//Mesh querying functions
	bool onBoundary(Edge * e) const;
	bool onBoundary(FaceData * f) const;
	bool onMeshBoundary(Edge * e) const;
	bool onMeshBoundary(FaceData * f) const;
	double vertDistTo(const double x, const double y, const double z);
	double vertDistTo(const PointData &);
	double vertDistTo(const Point3d *);
	void checkFaces() const;
	void checkDelauney() const;
	void printData() const { cout << numPoints << " points, " << numFaces << " faces\n";}
	double difference(Subdivision & mesh2);
	bool inBoundary(double x, double y) const;
	double errorScale(double x, double y) const;
	bool exists(PointData * p) const;

	//mesh output functions
	void outputFaces(vector<Point3d *> * pointGroup = NULL, bool fakeSurface = true);
	void outputPoints();
	void outputUV();

	//High level mesh manipulation functions
	bool insertBest();
	void insertAll();
	void clean(bool force = false);
	void smooth(int whichPoints = 0, int iterations = 1, vector<Point3d *> * skipPoints = NULL);
	void updatePointIndices();
	bool inMesh(double x, double y);
	bool inMesh(Point2d * p);
	void move(double dx, double dy, double dz);
	void stretch(double sx, double sy, double sz);
	void getFitPlane(double * a, double * b, double * c);
	void subtractPlane(double a, double b, double c);
	double getVariance();
	double getWidth() {return width;}
	double getHeight() {return height;}
	double getLeft() {return leftEdge;}
	double getRight() {return rightEdge;}
	double getTop() {return topEdge;}
	double getBottom() {return bottomEdge;}
	void findCenter(double * avex, double * avey, double * avez);
	void zOrderPoints();
	void nnOrderPoints();
	int removeSpikes(double threshold, double maxAngle = 90);
	int fixSlivers(double maxRho, bool force = false, vector<Point3d *> * skipPoints = NULL);
	//Deconstuctor
	~Subdivision();
};

//Returns true iff the edge is on the boundary of the bounding box
inline bool Subdivision::onBoundary(Edge * e) const {
	return (e->lDual() == NULL || e->rDual() == NULL);
}

//Returns ture iff the face has an edge on the boundary of the bounding box
inline bool Subdivision::onBoundary(FaceData * f) const {
	return (onBoundary(f->getAnchor()) || onBoundary(f->getAnchor()->lPrev()) || onBoundary(f->getAnchor()->lNext()));
}

//Returns ture iff the edge is on the boundary of the mesh
inline bool Subdivision::onMeshBoundary(Edge * e) const {
	return onBoundary(e) || e->lDual() == DELETED || e->rDual() == DELETED;
}

//Returns ture iff the face has an edge on the boundary of the mesh
inline bool Subdivision::onMeshBoundary(FaceData * f) const {
	return (onMeshBoundary(f->getAnchor()) || onMeshBoundary(f->getAnchor()->lPrev()) || onMeshBoundary(f->getAnchor()->lNext()));
}

//Returns true if the point is in the boundary of the box
inline bool Subdivision::inBoundary(double x, double y) const {
	Point2d p(x, y);
	return !cw(p, *boundingPoints[1], *boundingPoints[0]) &&
		!cw(p, *boundingPoints[2], *boundingPoints[1]) &&
		!cw(p, *boundingPoints[3], *boundingPoints[2]) &&
		!cw(p, *boundingPoints[0], *boundingPoints[3]);
}

// Inserts the best point of the supplied face into the subdivision
inline bool Subdivision::insertSite(FaceData * f) {
	f->getBest(utilP3d);
	if (utilP3d->getZ() != INFINITY)
		return insertSite(new PointData(utilP3d->getX(), utilP3d->getY(), utilP3d->getZ(), ++numPoints), f->getAnchor());
	else
		return insertSite(utilP3d->getX(), utilP3d->getY(), f->getAnchor(), false);
}

// Inserts a new point into the subdivision, checking to insure it is in the mesh first
inline bool Subdivision::insertSite(Point3d & p) {
	return (inBoundary(p.getX(), p.getY()) && this->insertSite(new PointData(p.getX(), p.getY(), p.getZ(), ++numPoints), NULL));
}

// Inserts a new point into the subdivision, checking to insure it is in the mesh first
inline bool Subdivision::insertSite(const double x, const double y, const double z) {
	return (inBoundary(x, y) && this->insertSite(new PointData(x, y, z, ++numPoints), NULL));
}

// Inserts a new point into the subdivision, checking to insure it is in the mesh first
inline bool Subdivision::insertSite(PointData * x) {
	return (inBoundary(x->getX(), x->getY()) && insertSite(x, NULL));
}

//Mark the given face as edited
inline void Subdivision::markEdited(FaceData * f) {
	if (inserttype != ALL && intype != PROGRAM) {
		editedFaces->insert(f);
	}
}

//Calculate the distance from p to the plane z = a*x + b*y + c, assuming f is described
//by this plane. If the error is greater than maxErr, update maxErr, and set the best
//point in f
inline void Subdivision::updateError(FaceData * f, double x, double y, double z, double a, double b, double c, double * maxErr) {
	if (z != INFINITY && fabs(z - a * x - b * y - c) > *maxErr) {
		f->setBest(x, y, z);
		*maxErr = fabs(z - a * x - b * y - c);
	}
}

//Calculate the distance from p to the plane z = a*x + b*y + c, assuming f is described
//by this plane. If the error is greater than maxErr, update maxErr, and set the best
//point in f
inline void Subdivision::updateError(FaceData * f, Point3d * p, double a, double b, double c, double * maxErr) {
	double z = p->getZ();
	if (fabs(z - a * p->getX() - b * p->getY() - c) > *maxErr && z != INFINITY) {
		f->setBest(p);
		*maxErr = fabs(z - a * p->getX() - b * p->getY() - c);
	}
}

//Same as updateError, except now we apply the scale specified by impMap to the error
inline void Subdivision::updateErrorWithScale(FaceData * f, double x, double y, double z, double a, double b, double c, double * maxErr) {
	if (z != INFINITY && fabs(z - a * x - b * y - c) * errorScale(x, y) > *maxErr) {
		f->setBest(x, y, z);
		*maxErr = fabs(z - a * x - b * y - c) * errorScale(x, y);
	}
}

//Same as updateError, except now we apply the scale specified by impMap to the error
inline void Subdivision::updateErrorWithScale(FaceData * f, Point3d * p, double a, double b, double c, double * maxErr) {
	double z = p->getZ();
	if (fabs(z - a * p->getX() - b * p->getY() - c) * errorScale(p->getX(), p->getY())> *maxErr && z != INFINITY) {
		f->setBest(p);
		*maxErr  = fabs(z - a * p->getX() - b * p->getY() - c) * errorScale(p->getX(), p->getY());
	}
}

//Return the vertical (z) difference to the specified point
//Note that an absolute value is not taken, so negative values will be returned
//if the point is below the mesh
inline double Subdivision::vertDistTo(const Point3d * p) {
	return vertDistTo(p->getX(), p->getY(), p->getZ());
}

//Return the vertical (z) difference to the specified point
//Note that an absolute value is not taken, so negative values will be returned
//if the point is below the mesh
inline double Subdivision::vertDistTo(const PointData & p) {
	return vertDistTo(p.getX(), p.getY(), p.getZ());
}

//Returns true if the given point is in boundary of the mesh
inline bool Subdivision::inMesh(Point2d * p) {
	return inMesh(p->getX(), p->getY());
}

//Reanchors a face
inline void Subdivision::reAnchorFace(FaceData * f, Edge * e) {
	if (isValid(f)) {
		f->reAnchor(e);
		markEdited(f);
	}
}

//Returns our best guess as to an edge near p
inline Edge* Subdivision::guessEdge(Point2d & p) {
	return startingEdge;
}

//Returns true iff f is a valid face in the mesh
inline bool Subdivision::isValid(FaceData * f) {
	return f && f != DELETED;
}

//Returns true iff e is a valid edge in the mesh
inline bool Subdivision::isValid(Edge * e) {
	return e && e->lDual() != DELETED;
}

//Treat a point as existing if it exists in the raw field, and we are not told to ignore it by the importance map
inline bool Subdivision::exists(PointData * p) const {
	return raw->exists(p) && this->errorScale(p->getX(), p->getY()) != 0;
}

#endif
