#include "subdivision.h"
#include "kdtree.h"
#include <limits>
#include <algorithm>
#include <vector>
#include <stdio.h>
#include <string.h>


#pragma warning(disable : 4244 4996)
#define PI 3.1415926535



using namespace std;

//A face to represent deleted faces
FaceData * Subdivision::DELETED = new FaceData();

//Store infinity, so we don't have to recalculate each time
double Subdivision::INFINITY = numeric_limits<double>::infinity();

/////////////////////////////////////////////////////////////////
//Constructors/initializer
/////////////////////////////////////////////////////////////////
//Create a new subdivision to represent the given RawField
//Scale x and y by scale, output results to filename, and use insertion scheme
//itype
Subdivision::Subdivision(RawField * field, const double scale, const char * filename, inserttype_t itype, RawField * impMap) {

	strcpy(outputFile, filename);
	this->impMap = impMap;
	if (field->type > 0) {
		intype = RASTER;
	} else if (field->type < 0) {
		intype = XYZ_FIELD;
	}

	inserttype = itype;

	this->scale = scale;
	if (intype == RASTER) {
		PointData p1(0, 0, field->getZ(0, 0), 1);
		PointData p2(0, field->height - 1, field->getZ(0, field->height - 1), 2);
		PointData p3(field->width - 1, field->height - 1, field->getZ(field->width - 1, field->height - 1), 3);
		PointData p4(field->width - 1, 0, field->getZ(field->width - 1, 0), 4);

		//The raw field may be telling us to leave these points out
		//We don't have that option, so insert them at 0
		if (p1.getZ() == INFINITY)
			p1.setZ(0);

		if (p2.getZ() == INFINITY)
			p2.setZ(0);

		if (p3.getZ() == INFINITY)
			p3.setZ(0);

		if (p4.getZ() == INFINITY)
			p4.setZ(0);

		this->raw = field;
		init(p1, p2, p3, p4);
	} else if (intype == XYZ_FIELD) {
		Point3d *p;
		PointData p1(field->xmin, field->ymin, field->avez, 1);
		p = field->getPointNear(field->xmin, field->ymin);
		//if (p->getX() == field->xmin && p->getY() == field->ymin)
			p1.setZ(p->getZ());

		PointData p2(field->xmin, field->ymax, field->avez, 2);
		p = field->getPointNear(field->xmin, field->ymax);
		//if (p->getX() == field->xmin && p->getY() == field->ymax)
			p2.setZ(p->getZ());

		PointData p3(field->xmax, field->ymax, field->avez, 3);
		p = field->getPointNear(field->xmax, field->ymax);
		//if (p->getX() == field->xmax && p->getY() == field->ymax)
			p3.setZ(p->getZ());

		PointData p4(field->xmax, field->ymin, field->avez, 4);
		p = field->getPointNear(field->xmax, field->ymin);
		//if (p->getX() == field->xmax && p->getY() == field->ymin)
			p4.setZ(p->getZ());

		this->raw = field;
		init(p1, p2, p3, p4);
	}
}

//Create a new subdivision of a quadrilateral defined by cw points
//a, b, c, d, and output file filename
Subdivision::Subdivision(const PointData& a, const PointData& b, 
						 const PointData& c, const PointData& d, 
						 const char * filename) {
							 impMap = NULL;
							 intype = PROGRAM;
							 init(a, b, c, d);
							 scale = 1;
							 strcpy(outputFile, filename);
}

//Just make a mesh with the four points a, b, c, d
//This mesh is not intended to be written to a file
Subdivision::Subdivision(const PointData& a, const PointData& b, 
						 const PointData& c, const PointData& d) {
							 impMap = NULL;
							 intype = PROGRAM;
							 scale = 1;
							 init(a, b, c, d);
}

//Create a new subdivision of a quadrilateral defined by cw points
//a, b, c, d, and output file filename
Subdivision::Subdivision(const Point2d& p1, const Point2d& p2, 
						 const Point2d& p3, const Point2d& p4, 
						 const char * filename) {
							 impMap = NULL;
							 intype = PROGRAM;
							 PointData a(p1.getX(), p1.getY());
							 PointData b(p2.getX(), p2.getY());
							 PointData c(p3.getX(), p3.getY());
							 PointData d(p4.getX(), p4.getY());
							 init(a, b, c, d);
							 scale = 1;
							 strcpy(outputFile, filename);

}

// Initialize a subdivision to the quadrilateral defined by the points a, b, c,d.
//a-b-c-d must be in cw order
void Subdivision::init(const PointData& a, const PointData& b, 
					   const PointData& c, const PointData& d)
{
	faceHeap = new Heap<FaceData*>();
	faces = new vector<FaceData*>();
	points = new vector<PointData*>();
	editedFaces = new set<FaceData*>();
	edges = new BlockAlloc<QuadEdge>(65536);
	allocFaces = new BlockAlloc<FaceData>(1024);

	utilPd = new PointData();
	utilP3d = new Point3d();

	numFaces = 0;
	numPoints = 4;

	//Copy the bounding points, insuring we have total control
	PointData *da, *db, *dc, *dd;
	da = new PointData(a), db = new PointData(b);
	dc = new PointData(c), dd = new PointData(d);

	points->push_back(da);
	points->push_back(db);
	points->push_back(dc);
	points->push_back(dd);

	//Save the bounding points
	boundingPoints[0] = da;
	boundingPoints[1] = db;
	boundingPoints[2] = dc;
	boundingPoints[3] = dd;

	//Construct the rectangle
	Edge* ea = makeEdge();
	ea->endPoints(da, db);
	Edge* eb = makeEdge();
	splice(ea->sym(), eb);
	eb->endPoints(db, dc);
	Edge* ec = makeEdge();
	splice(eb->sym(), ec);
	ec->endPoints(dc, dd);
	Edge *ed = makeEdge();
	splice(ec->sym(), ed);
	ed->endPoints(dd, da);
	splice(ed->sym(), ea);
	startingEdge = ea;

	//Construct the diagonal
	Edge *diag = makeEdge();
	splice(ed->sym(),diag);
	splice(eb->sym(),diag->sym());
	diag->endPoints(da,dc);

	FaceData *f1 = makeFace(ea->sym());
	FaceData *f2 = makeFace(ec->sym());

	ea->rot()->endPoints(f1, NULL);
	eb->rot()->endPoints(f1, NULL);
	ec->rot()->endPoints(f2, NULL);
	ed->rot()->endPoints(f2, NULL);
	diag->rot()->endPoints(f2, f1);

	leftEdge = da->getX();
	rightEdge = dc->getX();
	topEdge = dc->getY();
	bottomEdge = da->getY();
	width = rightEdge - leftEdge;
	height = topEdge - bottomEdge;

	if (inserttype == ALL) {
		insertAll();
	}
}


/////////////////////////////////////////////////////////////////
//Edge Level operations
//These functions edit the subdivision by changing the edges
//These are the only funtions that edit the underlying structure
//of the subdivision
/////////////////////////////////////////////////////////////////

// This operator affects the two edge rings around the origins of a and b,
// and, independently, the two edge rings around the left faces of a and b.
// In each case, (i) if the two rings are distinct, Splice will combine
// them into one; (ii) if the two are the same ring, Splice will break it
// into two separate pieces.
// Thus, Splice can be used both to attach the two edges together, and
// to break them apart. See Guibas and Stolfi (1985) p.96 for more details
// and illustrations.
void Subdivision::splice(Edge* a, Edge* b)
{
	Edge* alpha = a->oNext()->rot();
	Edge* beta  = b->oNext()->rot();

	Edge* t1 = b->oNext();
	Edge* t2 = a->oNext();
	Edge* t3 = beta->oNext();
	Edge* t4 = alpha->oNext();

	a->setNext(t1);
	b->setNext(t2);
	alpha->setNext(t3);
	beta->setNext(t4);
}

//Delete the specified edge from the mesh
void Subdivision::deleteEdge(Edge* e)
{
	splice(e, e->oPrev());
	splice(e->sym(), e->sym()->oPrev());
	edges->recycle(e->qEdge());
}

// Add a new edge e connecting the destination of a to the
// origin of b, in such a way that all three have the same
// left face after the connection is complete.
// Additionally, the data pointers of the new edge are set.
Edge* Subdivision::Connect(Edge* a, Edge* b)
{
	Edge* e = makeEdge();
	e->endPoints(a->dest(), b->org());
	e->rot()->endPoints(a->lNext()->lDual(), a->lDual());
	splice(e, a->lNext());
	splice(e->sym(), b);
	return e;
}

// Essentially turns edge e counterclockwise inside its enclosing
// quadrilateral. The data pointers are modified accordingly.
void Subdivision::swap(Edge* e)
{
	Edge* a = e->oPrev();
	Edge* b = e->sym()->oPrev();
	//The two faces which will change
	FaceData *fa = (FaceData*)a->lDual();
	FaceData *fb = (FaceData*)b->lDual();

	//update the anchors to edges whose faces will not change

	reAnchorFace(fa, e->dNext());
	reAnchorFace(fb, e->oNext()->sym());

	//Change the endpoints and faces
	e->endPoints(a->dest(), b->dest());
	a->rot()->endPoints(a->rDual(), fb);
	b->rot()->endPoints(b->rDual(), fa);

	//Resplice the data structure
	splice(e, a);
	splice(e->sym(), b);
	splice(e, a->lNext());
	splice(e->sym(), b->lNext());
}

//Make a new edge
Edge* Subdivision::makeEdge()
{
	QuadEdge * ql = edges->getNext();
	//edges->push_back(ql);
	return ql->first();
}

//Swaps all edges necessary to insure the Delauney condition is satisfied after the point x is
//inserted, starting with the edge e
void Subdivision::swapEdges(Edge * e, PointData * x, Edge * end) {
	// fixes the affected edges so that the result
	// is still a Delaunay triangulation. This is based on the
	// pseudocode from Guibas and Stolfi (1985) p.120, with slight
	// modifications and a bug fix.
	do {
		Edge* t = e->oPrev();
		if (!onMeshBoundary(e) &&
			inCircle(e->org2d(), t->dest2d(), e->dest2d(), *x)) {
				swap(e);
				e = e->oPrev();
		}
		else if (e->oNext() == end) {  // no more suspect edges
			break;
		}
		else  // pop a suspect edge
			e = e->oNext()->lPrev();
	} while (TRUE);
}

/////////////////////////////////////////////////////////////////
//Edge data manipulation functions
//These functions relate to the faces and points referenced by
//the underlying quadedge structure
/////////////////////////////////////////////////////////////////


//Update the error for this particular face, and set the best point for that
//face
inline void Subdivision::updateFace(FaceData * f) {
	impMap ? updateFaceWithScale(f): updateFaceWithoutScale(f);
}

//Update the face using the faster functions which do not check to see if we should be
//scaling the error. This is ~3x faster than the functions with scaling, when there is
//no scaling to be done
void Subdivision::updateFaceWithoutScale(FaceData * f) {
	double trueVal;
	double maxErr = 0;

	if (intype == RASTER ) {
		if (inserttype == CICRCUMCENTER || inserttype == CICRCUMCENTER_VOL) {
			double x = 0, y = 0;
			f->getCircumcenter(&x, &y);

			if (!inBoundary(x, y)) {
				x = (x < 0 ? 0 : (x > width ? width : x));
				y = (y < 0 ? 0 : (y > height ? height : y));
			}
			trueVal = raw->getZ((int) x, (int) y);
			maxErr = (trueVal == INFINITY ? 0 : fabs(vertDistTo((int) x, (int) y, trueVal)));
			if (inserttype == CICRCUMCENTER_VOL) {
				maxErr *= f->area();
			}

			f->setBest((int) x, (int) y, trueVal);
		} else if (inserttype == MAX_ERR || inserttype == MAX_ERR_VOL) {
			//scan convert the triangle

			//The points in acending order of y
			int ordx[3], ordy[3];
			f->getPointsOrdY(ordx, ordy);
			maxErr = 0;

			//The coefficients representing the plane of the triangle
			//z = a * x + b * y + c
			double a, b, c;
			f->getPlane(&a, &b, &c);

			//Scan convert from the bottom of the triangle to the middle vertical point
			if (ordy[0] != ordy[1]) {
				double dx1 = (double) (ordx[1] - ordx[0]) / (double) (ordy[1] - ordy[0]);
				double dx2 = (double) (ordx[2] - ordx[0]) / (double) (ordy[2] - ordy[0]);
				//Iterate over all grid points in the triangle
				for (int y = ordy[0]; y < ordy[1]; y++) {
					int xmin = ordx[0] + ceil(min(dx1 * (y - ordy[0]), dx2 * (y - ordy[0])));
					int xmax = ordx[0] + floor(max(dx1 * (y - ordy[0]), dx2 * (y - ordy[0])));
					for (int x = xmin; x <= xmax; x++) {
						updateError(f, x, y, raw->getZ(x, y), a, b, c, &maxErr);
					}
				}
			} else {
				//The y values are equal, just walk on the horizontal edge of the triangle
				for (int x = min(ordx[0], ordx[1]); x <= max(ordx[0], ordx[1]); x++) {
					updateError(f, x, ordy[0], raw->getZ(x, ordy[0]), a, b, c, &maxErr);
				}
			}

			//Scan convert from the middle of the triangle to the top point
			if (ordy[1] != ordy[2]) {
				double dx0 = (double) (ordx[0] - ordx[2]) / (double) (ordy[0] - ordy[2]);
				double dx1 = (double) (ordx[1] - ordx[2]) / (double) (ordy[1] - ordy[2]);
				//Iterate over all grid points in the triangle
				for (int y = ordy[1]; y <= ordy[2]; y++) {
					int xmin = ordx[2] + ceil(min(dx0 * (y - ordy[2]), dx1 * (y - ordy[2])));
					int xmax = ordx[2] + floor(max(dx0 * (y - ordy[2]), dx1 * (y - ordy[2])));
					for (int x = xmin; x <= xmax; x++) {
						updateError(f, x, y, raw->getZ(x, y), a, b, c, &maxErr);
					}
				}
			} else {
				//The y values are equal, just walk on the horizontal edge of the triangle
				utilP3d->setY(ordy[2]);
				for (int x = min(ordx[2], ordx[1]); x <= max(ordx[2], ordx[1]); x++) {
					updateError(f, x, ordy[2], raw->getZ(x, ordy[2]), a, b, c, &maxErr);
				}
			}

			if (inserttype == MAX_ERR_VOL) {
				maxErr *= f->area();
			}
		}
	} else if (intype == XYZ_FIELD) {
		if (inserttype == CICRCUMCENTER || inserttype == CICRCUMCENTER_VOL) {
			double x = 0, y = 0;
			f->getCircumcenter(&x, &y);
			Point3d * point = raw->getPointNear(x, y);
			maxErr = (point->getZ() == INFINITY ? 
				0 : fabs(vertDistTo(point->getX(), point->getY(), point->getZ())));
			if (inserttype == this->CICRCUMCENTER_VOL) {
				maxErr *= f->area();
			}
			f->setBest(point);
		} else if (inserttype == MAX_ERR || inserttype == MAX_ERR_VOL) {
			//Find all points in the triangle, and evaluate errors

			//The verticies of the triangle, in accending order of y
			double ordx[3], ordy[3];
			f->getPointsOrdY(ordx, ordy);
			maxErr = 0;

			//The coefficients representing the plane of the triangle
			//z = a * x + b * y + c
			double a, b, c;
			f->getPlane(&a, &b, &c);

			double lowx = ordx[0];
			double highx = ordx[0];
			lowx = min(ordx[1], lowx);
			lowx = min(ordx[2], lowx);
			highx = max(ordx[1], highx);
			highx = max(ordx[2], highx);
			vector<Point3d*> points;;
			raw->getPointsInBox(lowx, ordy[0], highx, ordy[2], &points);
			//Iterate over all points
			for(vector<Point3d*>::iterator it = points.begin(); it != points.end(); ++it) {
				double x, y;
				Point3d * p = *it;
				x = p->getX();
				y = p->getY();
				//Determine if this point is in the triangle
				//If so, update the error
				if (y < ordy[1]) {
					//ordy[0] <= y < ordy[1]
					double dx1 = (ordx[1] - ordx[0]) / (ordy[1] - ordy[0]);
					double dx2 = (ordx[2] - ordx[0]) / (ordy[2] - ordy[0]);
					double x1 = ordx[0] + dx1 * (y - ordy[0]);
					double x2 = ordx[0] + dx2 * (y - ordy[0]);
					if (x >= min(x1, x2) && x <= max(x1, x2)) {
						updateError(f, p, a, b, c, &maxErr);
					}
				} else if (y > ordy[1]) {
					//ordy[1] < y <= ordy[2]
					double dx0 = (ordx[0] - ordx[2]) / (ordy[0] - ordy[2]);
					double dx1 = (ordx[1] - ordx[2]) / (ordy[1] - ordy[2]);

					double x0 = ordx[2] + dx0 * (y - ordy[2]);
					double x1 = ordx[2] + dx1 * (y - ordy[2]);

					if (x >= min(x0, x1) && x <= max(x0, x1)) {
						updateError(f, p, a, b, c, &maxErr);
					}
				} else {
					//y = ordy[1]
					double dx0 = (ordx[2] - ordx[0]) / (ordy[2] - ordy[0]);
					double x0 = ordx[0] + dx0 * (y - ordy[0]);
					if (x >= min(ordx[1], x0) && x <= max(ordx[1], x0)) {
						updateError(f, p, a, b, c, &maxErr);
					}
				}
			}

			if (inserttype == MAX_ERR_VOL) {
				maxErr *= f->area();
			}
		}
	}

	//We know the error of this triangle
	//Update it's position in the heap
	faceHeap->update(f->loc, maxErr);
}

//Update the face using the slower functions which allow the error to be scaled based on
//position. This allows for the user to specify a certain region as more important than
//the rest, at a relativly high runtime cost, as to do the scaling, we now make two calls
//to RawField::getZ(), as well as being forced to do the computation to determine which
//image scale pixel we should use
void Subdivision::updateFaceWithScale(FaceData * f) {
	double trueVal;
	double maxErr = 0;
	if (intype == RASTER ) {
		if (inserttype == CICRCUMCENTER || inserttype == CICRCUMCENTER_VOL) {
			double x = 0, y = 0;
			f->getCircumcenter(&x, &y);

			if (!inBoundary(x, y)) {
				x = (x < 0 ? 0 : (x > width ? width : x));
				y = (y < 0 ? 0 : (y > height ? height : y));
			}
			trueVal = raw->getZ((int) x, (int) y);
			maxErr = (trueVal == INFINITY ? 0 : fabs(vertDistTo((int) x, (int) y, trueVal)));
			if (inserttype == CICRCUMCENTER_VOL) {
				maxErr *= f->area();
			}

			maxErr *= errorScale(x, y);
			f->setBest((int) x, (int) y, trueVal);
		} else if (inserttype == MAX_ERR || inserttype == MAX_ERR_VOL) {
			//scan convert the triangle

			//The points in acending order of y
			int ordx[3], ordy[3];
			f->getPointsOrdY(ordx, ordy);
			maxErr = 0;

			//The coefficients representing the plane of the triangle
			//z = a * x + b * y + c
			double a, b, c;
			f->getPlane(&a, &b, &c);

			//Scan convert from the bottom of the triangle to the middle vertical point
			if (ordy[0] != ordy[1]) {
				double dx1 = (double) (ordx[1] - ordx[0]) / (double) (ordy[1] - ordy[0]);
				double dx2 = (double) (ordx[2] - ordx[0]) / (double) (ordy[2] - ordy[0]);
				//Iterate over all grid points in the triangle
				for (int y = ordy[0]; y < ordy[1]; y++) {
					int xmin = ordx[0] + ceil(min(dx1 * (y - ordy[0]), dx2 * (y - ordy[0])));
					int xmax = ordx[0] + floor(max(dx1 * (y - ordy[0]), dx2 * (y - ordy[0])));
					for (int x = xmin; x <= xmax; x++) {
						updateErrorWithScale(f, x, y, raw->getZ(x, y), a, b, c, &maxErr);
					}
				}
			} else {
				//The y values are equal, just walk on the horizontal edge of the triangle
				for (int x = min(ordx[0], ordx[1]); x <= max(ordx[0], ordx[1]); x++) {
					updateErrorWithScale(f, x, ordy[0], raw->getZ(x, ordy[0]), a, b, c, &maxErr);
				}
			}

			//Scan convert from the middle of the triangle to the top point
			if (ordy[1] != ordy[2]) {
				double dx0 = (double) (ordx[0] - ordx[2]) / (double) (ordy[0] - ordy[2]);
				double dx1 = (double) (ordx[1] - ordx[2]) / (double) (ordy[1] - ordy[2]);
				//Iterate over all grid points in the triangle
				for (int y = ordy[1]; y <= ordy[2]; y++) {
					int xmin = ordx[2] + ceil(min(dx0 * (y - ordy[2]), dx1 * (y - ordy[2])));
					int xmax = ordx[2] + floor(max(dx0 * (y - ordy[2]), dx1 * (y - ordy[2])));
					for (int x = xmin; x <= xmax; x++) {
						updateErrorWithScale(f, x, y, raw->getZ(x, y), a, b, c, &maxErr);
					}
				}
			} else {
				//The y values are equal, just walk on the horizontal edge of the triangle
				utilP3d->setY(ordy[2]);
				for (int x = min(ordx[2], ordx[1]); x <= max(ordx[2], ordx[1]); x++) {
					updateErrorWithScale(f, x, ordy[2], raw->getZ(x, ordy[2]), a, b, c, &maxErr);
				}
			}

			if (inserttype == MAX_ERR_VOL) {
				maxErr *= f->area();
			}
		}
	} else if (intype == XYZ_FIELD) {
		if (inserttype == CICRCUMCENTER || inserttype == CICRCUMCENTER_VOL) {
			double x = 0, y = 0;
			f->getCircumcenter(&x, &y);

			Point3d * point = raw->getPointNear(x, y);
			maxErr = (point->getZ() == INFINITY ? 
				0 : fabs(vertDistTo(point->getX(), point->getY(), point->getZ())));
			if (inserttype == this->CICRCUMCENTER_VOL) {
				maxErr *= f->area();
			}
			maxErr *= errorScale(point->getX(), point->getY());
			f->setBest(point);
		} else if (inserttype == MAX_ERR || inserttype == MAX_ERR_VOL) {
			//Find all points in the triangle, and evaluate errors

			//The verticies of the triangle, in accending order of y
			double ordx[3], ordy[3];
			f->getPointsOrdY(ordx, ordy);
			maxErr = 0;

			//The coefficients representing the plane of the triangle
			//z = a * x + b * y + c
			double a, b, c;
			f->getPlane(&a, &b, &c);

			double lowx = ordx[0];
			double highx = ordx[0];
			lowx = min(ordx[1], lowx);
			lowx = min(ordx[2], lowx);
			highx = max(ordx[1], highx);
			highx = max(ordx[2], highx);
			vector<Point3d*> * points = new vector<Point3d*>();
			raw->getPointsInBox(lowx, ordy[0], highx, ordy[2], points);
			//Iterate over all points in the bounding box of the triangle
			for(vector<Point3d*>::iterator it = points->begin(); it!= points->end(); ++it) {
				double x, y;
				Point3d * p = *it;
				x = p->getX();
				y = p->getY();
				//Determine if this point is in the triangle
				//If so, update the error
				if (y < ordy[1]) {
					//ordy[0] <= y < ordy[1]
					double dx1 = (ordx[1] - ordx[0]) / (ordy[1] - ordy[0]);
					double dx2 = (ordx[2] - ordx[0]) / (ordy[2] - ordy[0]);
					double x1 = ordx[0] + dx1 * (y - ordy[0]);
					double x2 = ordx[0] + dx2 * (y - ordy[0]);
					if (x >= min(x1, x2) && x <= max(x1, x2)) {
						updateErrorWithScale(f, p, a, b, c, &maxErr);
					}
				} else if (y > ordy[1]) {
					//ordy[1] < y <= ordy[2]
					double dx0 = (ordx[0] - ordx[2]) / (ordy[0] - ordy[2]);
					double dx1 = (ordx[1] - ordx[2]) / (ordy[1] - ordy[2]);

					double x0 = ordx[2] + dx0 * (y - ordy[2]);
					double x1 = ordx[2] + dx1 * (y - ordy[2]);

					if (x >= min(x0, x1) && x <= max(x0, x1)) {
						updateErrorWithScale(f, p, a, b, c, &maxErr);
					}
				} else {
					//y = ordy[1]
					double dx0 = (ordx[2] - ordx[0]) / (ordy[2] - ordy[0]);
					double x0 = ordx[0] + dx0 * (y - ordy[0]);
					if (x >= min(ordx[1], x0) && x <= max(ordx[1], x0)) {
						updateErrorWithScale(f, p, a, b, c, &maxErr);
					}
				}
			}
			delete points;

			if (inserttype == MAX_ERR_VOL) {
				maxErr *= f->area();
			}
		}
	}

	//We know the error of this triangle
	//Update it's position in the heap
	faceHeap->update(f->loc, maxErr);
}

//Update the heap of faces to reflect the new error
//If we are doing a MAX_ERR reduction, only faces in editedFaces are updated
//Otherwise, faces in editedFaces and neigboring editedFaces are updated
void Subdivision::updateHeap() {
	if (inserttype == ALL || intype == PROGRAM) {
		editedFaces->clear();
		return;
	}

	for(set<FaceData*>::iterator it = editedFaces->begin(); 
		it != editedFaces->end(); 
		it++) {
			FaceData * f = *it;
			updateFace(f);
	}

	editedFaces->clear();
}

//Make a new face in the subdivision, and update all data pointers
//properly
FaceData *Subdivision::makeFace(Edge *e)
{
	numFaces++;
	FaceData *f = allocFaces->getNext();
	f->reAnchor(e);
	//Set the left face of the anchoring edge to this face
	e->rot()->endPoints(e->rDual(), f);

	//Put the face on the front of the list if it is on the boundary
	//This makes deleting points (usually on the edges), faster

	faces->push_back(f);
	f->val = -1;
	faceHeap->insert(f);
	markEdited(f);
	return f;
}

//"DELETE" a face. The edges still remain, so we can still traverse the mesh
//effectively, but the program now knows that this face is not actually in the
//completed mesh
void Subdivision::deleteFace(FaceData * face) {
	//Delete all references to the face
	updateHeap();
	//Make sure there is a face to delete
	if (!isValid(face)) {
		return;
	}

	Edge * e;
	for(e = face->getAnchor(); e != face->getAnchor()->lPrev(); e = e->lNext()) {
		e->rot()->endPoints(e->rDual(), DELETED);
	}

	e->rot()->endPoints(e->rDual(), DELETED);
	//vector<FaceData *>::iterator it = find(faces->begin(), faces->end(), face);

	//Erase the face from the vector of faces. This is an O(n) operation, which is sad, but
	//nothing can be done about it without using a hash map. This can be made a little faster
	//by using a list, but that carries the associated price of memory, and the find operation
	//is still O(n)
	faces->erase(find(faces->begin(), faces->end(), face));
	faceHeap->kill(face->loc);
	allocFaces->recycle(face);
	numFaces--;
}

//Delete a point on the boudary
//This is used to create a large boundary to put a mesh in, then delete it
//The deletion removes the point, and all faces from the point and face lists
//The deletion does not, however, delete the edges. They are preserved, so that
//the edges can still be walked along
//The pointers from the edges to their adjacent faces are switched such that they
//point to the DELETED face, indicating the face is not really there
bool Subdivision::deleteExternal(PointData * point) {
	Edge * e = locate(*point);
	Edge * e2;
	if (!(*point == e->org2d()) && !(*point == e->dest2d())) { // point is not in mesh
		return false;
	}

	if (*point == e->dest2d()) {
		e = e->sym();
	}

	//Make sure the point is on the boudary. If the point is on the boundary, one
	//of its edges is on the boundary
	for (e2 = e->oNext();;e2 = e2->oNext()) {
		if (onBoundary(e2)) {
			e = e2;
			break;
		} else if (e2 == e) {
			return false;
		}
	}

	//Don't try to delete the last face in the mesh, that is bad
	if (numFaces < 2) {
		return false;
	}

	if (e->lDual() == NULL) {
		e = e->oNext();
		e2 = e;
	}

	deleteFace((FaceData *) e->lDual());

	for (e = e2->oNext(); e != e2; e = e->oNext()) {
		deleteFace((FaceData*) e->lDual());


		//Don't try to delete the last face in the mesh, that is bad
		if (numFaces < 2) {
			startingEdge = (*faces->begin())->getAnchor();
			return false;
		}
	}

	startingEdge = (*faces->begin())->getAnchor();
	points->erase(find(points->begin(), points->end(), point));
	point->setIndex(0);
	numPoints--;
	return true;
}


/////////////////////////////////////////////////////////////////
//Insertion functions
/////////////////////////////////////////////////////////////////

// Inserts a new point near the specified x,y value, with a z value based on the raw field
bool Subdivision::insertSite(const double x, const double y, Edge * guess, bool force, bool knownEdge) {
	if (!inBoundary(x, y))
		return false;

	if (intype == RASTER) {
		bool inserted = false;
		if (raw->getZ((int) x, (int) y) != INFINITY) {
			//We have a good value to insert, do so
			inserted = this->insertSite(new PointData((int) x, (int) y, raw->getZ((int) x, (int) y), ++numPoints), guess);
		} else {
			//We have no data, just put the point at a z value that fits
			Point2d p((int) x, (int) y);
			Edge * e = locate(p, guess);
			if (isValid(e)) {
				double a = 0, b = 0, c = 0;
				((FaceData*) e->lDual())->getPlane(&a, &b, &c);
				inserted = this->insertSite(new PointData((int) x, (int) y, a * (int) x + b * (int) y + c, ++numPoints), e, true);
			} else {
				inserted = false;
			}
		}

		if (inserted)
			return true;

		if (force) {
			Point2d p(x, y);
			Edge * e = (knownEdge ? guess : locate(p, guess));
			if (isValid(e)) {
				double a = 0, b = 0, c = 0;
				((FaceData*) e->lDual())->getPlane(&a, &b, &c);
				inserted = this->insertSite(new PointData(x, y, a * x + b * y + c, ++numPoints), e, true);
			} else {
				inserted = false;
			}
		}

		return inserted;
	} else if (intype == XYZ_FIELD) {
		Point3d * p = raw->getPointNear(x, y);
		if (!force) {
			return this->insertSite(new PointData(p->getX(), p->getY(), p->getZ(), ++numPoints), guess);
		} else {
			//We have been told to force the point into the xyz-filed, if possible
			if (this->insertSite(new PointData(p->getX(), p->getY(), p->getZ(), ++numPoints), guess)) {
				return true;
			} else {
				Point2d p(x, y);
				Edge * e = (knownEdge ? guess : locate(p, guess));

				if (isValid(e)) {
					double a = 0, b = 0, c = 0;
					((FaceData*) e->lDual())->getPlane(&a, &b, &c);
					return this->insertSite(new PointData(x, y, a * x + b * y + c, ++numPoints), e, true);
				} else {
					return false;
				}
			}
		}
	} else if (intype == PROGRAM) {
		Point2d p(x, y);
		Edge * e = (knownEdge ? guess : locate(p, guess));
		
		if (isValid(e)) {
			double a = 0, b = 0, c = 0;
			((FaceData*) e->lDual())->getPlane(&a, &b, &c);
			return this->insertSite(new PointData(x, y, a * x + b * y + c, ++numPoints), e, true);
		} else {
			return false;
		}
	}

	return false;
}

// Inserts a new point into the subdivision, taking guess as the probable
//edge near the point
//Does not check to insure the point is in the mesh, but Locate will find out if it is not
//If it is unknown whether a point is in the mesh, call one of the overloads which does a quick check first
bool Subdivision::insertSite(PointData * x, Edge * guess, bool knownEdge) {

	//Don't locate the point if this is a known edge
	Edge* e = (knownEdge ? guess : locate(*x, guess));
	//Edge * e = locate(*x, guess);
	if (!isValid(e)) {
		//cout << "Error: point outside of mesh boundary\n";
		//cout << *x;
		numPoints--;
		return false;
	}

	Edge* base;
	//The first edge we inserted
	Edge* start;
	int connections = 2;

	if ((*x == e->org2d()) || (*x == e->dest2d())) { // point is already in
		numPoints--;
		//cout << "Point already in mesh\n";
		//cout << *x;
		return false;
	}

	points->push_back(x);

	//The block could probably be combined to handle the edge case and interior case together
	//If the point is on a boundary edge, the point is on the boundary
	if (onBoundary(e) && OnEdge(*x, e)) {
		if (!isValid((FaceData *) e->lDual())) {
			//Make the face on the left the inside
			e = e->sym();
		}

		Edge * oldEdge = e;

		//Make one new face, anchor the one exising face on an edge
		reAnchorFace((FaceData *) e->lDual(), e->lPrev());

		//The starting edge is the edge we are about to remove from the mesh
		//fix that
		startingEdge = e->lPrev();

		makeFace(e->lNext());

		base = makeEdge();
		start = base;
		startingEdge = base;
		base->endPoints(e->org(), x);
		base->rot()->endPoints(NULL, e->lPrev()->lDual());
		splice(base, e);

		//The edge that is going to be deleted
		e->rot()->endPoints(NULL, NULL);

		for (;connections > 0; connections--) {
			base = Connect(e, base->sym());
			e = base->oPrev();
		}

		//Delete the old edge
		deleteEdge(oldEdge);
		//Make sure the Delaunay swapping doesn't touch the boundary and get lost
		start = start->dNext();

	} else {

		if (OnEdge(*x, e)) {
			if (e->rDual() == DELETED) {
				e = e->sym();
			}
			e = e->oPrev();
			deleteEdge(e->oNext());
			connections = 3;
			reAnchorFace((FaceData*) e->lPrev()->lDual(), e->lPrev());
		}

		//Make two new faces, anchor the one exising face on an edge
		reAnchorFace((FaceData*) e->lDual(), e);
		makeFace(e->lNext());
		if (e->lNext()->lNext()->lDual() != DELETED) {
			makeFace(e->lNext()->lNext());
		}

		base = makeEdge();
		start = base;
		startingEdge = base;
		base->endPoints(e->org(), x);
		base->rot()->endPoints(e->lDual(), e->lPrev()->lDual());
		splice(base, e);

		// Connect the new point to the vertices of the containing
		// triangle (or quadrilateral, if the new point fell on an
		// existing edge.)
		for (;connections > 0; connections--) {
			base = Connect(e, base->sym());
			e = base->oPrev();
		}

	}

	// Examine suspect edges to ensure that the Delaunay condition
	// is satisfied.
	swapEdges(e, x, start);
	return true;
}


//Insert the best point (ie, the point with the maximum error based on our metric)
bool Subdivision::insertBest() {
	updateHeap();
	//Get the best face off the top of the heap
	FaceData * f = (FaceData *) faceHeap->top();
	bool inserted;
	//Insert the best point in that face
	inserted = this->insertSite(f);
	if (!inserted) {
		//Set the error of this point to -1
		faceHeap->update(0, -1);
	}
	return inserted;
}

//Insert every point in the raw field
void Subdivision::insertAll() {
	//Make sure the inserttype is ALL
	//This insures the error heap is not continually updated
	//Obviously, we don't care about the error heap when we insert all points
	this->inserttype = this->ALL;
	int i = 0;
	if (intype == RASTER) {
		for (int y = 0; y < raw->height; y++) {
			for (int x = 0; x < raw->width; x++) {
				if (raw->getZ(x, y) != INFINITY && insertSite(x, y, raw->getZ(x, y)) && !(++i%100000))
					cout << "Inserted: " << i << " points" << "\n";
				updateHeap();
			}
		}
	} else if (intype == XYZ_FIELD) {
		for (vector<Point3d*>::iterator id = raw->getPoints()->begin();
			id != raw->getPoints()->end(); ++id) {
				if ((*id)->getZ() != INFINITY && this->insertSite(*(*id)) && !(++i%100000)) {
					cout << "Inserted: " << i << " points" << "\n";	
				}
				updateHeap();
		}
	}
}

// Returns an edge e, s.t. the triangle to the left of e is interior to the
// subdivision and either x is on e (inclusive of endpoints) or x lies in the
// interior of the triangle to the left of e.
// The search starts from either hintedge, if it is not NULL, else
// startingEdge, and proceeds in the general direction of x.
Edge *Subdivision::locate(Point2d& x, Edge * guess) {
	// Algorithm is a variant of Green and Sibson's walking method for
	// point location, as described by Guibas and Stolfi (ACM Trans. on Graphics,
	// Apr. 1985, p.121), but modified in three ways:
	//	* Supports queries on perimeter of subdivision,
	//	  provided perimeter is convex.
	//	* Uses two area computations per step, not three.

	Edge * e = (guess ? guess : guessEdge(x)), *eo, *ed;
	double t, to, td;


	t = triArea(x, e->dest2d(), e->org2d());
	if (t>0) {			// x is to the right of edge e
		t = -t;
		e = e->sym();
	}

	// x is on e or to the left of e

	// edges e, eo, ed point upward in the diagram below:
	//
	//         /|
	//     ed / |
	//       /  |
	//      /   |
	//     /    |
	//     \    | e
	//      \   |
	//       \  |
	//     eo \ |
	//         \|

	double xx = x.getX();
	double xy = x.getY();
	double ox = ((PointData *) e->org())->getX();
	double oy = ((PointData *) e->org())->getY();
	double dx = ((PointData *) e->dest())->getX();
	double dy = ((PointData *) e->dest())->getY();

	while (TRUE) {
		eo = e->oNext();

		double eox = ((PointData *)eo->dest())->getX();
		double eoy = ((PointData *)eo->dest())->getY();
		double edx = eox;
		double edy = eoy;

		if (e->lDual() == NULL) {
			//The point is to the left of e, and left of e is outside the boudary
			//The point must be outside the boundary
			if (t < 0)
				return NULL;

			//The origin of ed and the destination of eo may not be the same
			ed = e->dPrev();
			edx = ((PointData *)ed->org())->getX();
			edy = ((PointData *)ed->org())->getY();
		}

		//There is a lot of unwrapping and chaching here, but the point is to get
		//the following
		//to = triArea(x, eo->dest2d(), eo->org2d());
		//td = triArea(x, ed->dest2d(), ed->org2d());
		to = triArea(xx, xy, eox, eoy, ox, oy);
		td = triArea(xx, xy, dx, dy, edx, edy);

		//Debug code, incase locate gets into an infinite loop
		//it happens, usually due to a slient error elsewhere
		/*if (numPoints >= 292) {
			cout << "\n";
			cout << x << " " << numFaces << "\n";
			cout << to << " " << td << "\n";
			cout << "e " << *((PointData*) e->org()) << " " << *((PointData*) e->dest()) << "\n";
			cout << "eo " << *((PointData*) eo->org()) << " " << *((PointData*) eo->dest()) << "\n";
			cout << "\n";
		}*/
		if (td>0)			// x is below ed
			if (to>0 || to==0 && t==0) {// x is interior, or origin endpoint
				startingEdge = e;
				return e;
			}
			else {			// x is below ed, below eo
				t = to;
				e = eo;
				dx = eox;
				dy = eoy;
			}
		else				// x is on or above ed
			if (to>0)			// x is above eo
				if (td==0 && t==0) {	// x is destination endpoint
					startingEdge = e;
					return e;
				}
				else {			// x is on or above ed and above eo
					t = td;
					e = e->dPrev();
					ox = edx;
					oy = edy;
				}
			else			// x is on or below eo
				if (t==0 && this->onBoundary(e)) {
					// x on e but subdiv. is to right
					e = e->sym();
					edx = ox;
					edy = oy;
					ox = dx;
					oy = dy;
					dx = edx;
					dy = edy;
				}
				else {
					t = td;
					e = e->dPrev();
					ox = edx;
					oy = edy;
				}
	}
}

/////////////////////////////////////////////////////////////////
//Global manipulation/high level querry functions
//These are functions which affect the entire subdivision,
//and are the funtions likely to be called by code outside of the
//subdivision class.
/////////////////////////////////////////////////////////////////

//There is redudent data stored in the structure in regards to the face data
//Make sure it is consistent
void Subdivision::checkFaces() const {
	if (faces->size() != numFaces) {
		cout << "List of faces disagrees in size with numfaces\n";
		exit(1);
	}

	if (faces->size() != faceHeap->heap_size()) {
		cout << "List of faces is a different size than the heap\n";
		exit(1);
	}

	int problems = 0;
	FaceData * f;
	for(vector<FaceData*>::iterator it = faces->begin(); it != faces->end(); ++it) {
		f = *it;

		if (f->loc == NOT_IN_HEAP || (*faceHeap)[f->loc] != f) {
			cout << "ERROR, face not in heap\n";
			problems++;
		}
		
		if (((FaceData *) f->getAnchor()->lDual()) != f) {
			cout << "Check Faces ERROR: inconsistent face data\n";
			problems++;
		} else if (((FaceData *) f->getAnchor()->lNext()->lDual()) != f ||
			((FaceData *) f->getAnchor()->lPrev()->lDual()) != f) {
				cout << "Check Faces ERROR: inconsistent face data\n";
				problems++;
		}
	}

	if (!problems)
		cout << "All face data is consitent\n";
	else
		cout << "Found " << problems << " incosistencies in the face data\n";
}

//Check to make sure that the triangluation is infact Delauney
//Note that this check does not insure the triangluation is correct, but will catch most errors
void Subdivision::checkDelauney() const {
	int problems = 0;
	FaceData * f;
	for (vector<FaceData*>::iterator it = faces->begin(); it != faces->end(); ++it) {
		f = *it;
		if (inCircle(f->getAnchor()->org2d(), f->getAnchor()->dest2d(), f->getAnchor()->lNext()->dest2d(), f->getAnchor()->oPrev()->dest2d()) ||
			inCircle(f->getAnchor()->org2d(), f->getAnchor()->dest2d(), f->getAnchor()->lNext()->dest2d(), f->getAnchor()->lNext()->oPrev()->dest2d()) ||
			inCircle(f->getAnchor()->org2d(), f->getAnchor()->dest2d(), f->getAnchor()->lNext()->dest2d(), f->getAnchor()->lPrev()->oPrev()->dest2d())) {
				cout << "Delaunay condition violated\n";
				problems++;
		}
	}

	if (!problems)
		cout << "Found no violations of the Deluanay condition\n";
	else
		cout << "Delaunay condition violated " << problems << " times\n";
}

//Delete External points, if necessary
void Subdivision::clean(bool force) {
	for (int i = 0; i < 4; i++) {
		if (!exists(boundingPoints[i]) || force)
			deleteExternal(boundingPoints[i]);
	}
}

//Fix all triangles with a Rho greater than maxRho
int Subdivision::fixSlivers(double maxRho, bool force, vector<Point3d *> * skipPoints) {
	//Iterate over all faces
	//For each face, if the mesh is actually changed, iterate over all affected faces until they are
	//all fixed. Then, proceed to the next face, until all faces are fixed.
	//This strategy has the effect of fixing all the slivers in a single area before moving on to the
	//next area, which makes a more efficent use of points, as well as improving runtime
	//Faces are still checked more than once, so some improvement could probably be made.

	struct LessThanY {
		bool operator()(Point2d * p1, Point2d * p2) {
			return p1->getY() < p2->getY();
		}
	} lty;

	if (skipPoints)
		sort(skipPoints->begin(), skipPoints->end(), lty);

	updateHeap();
	FaceData * f;
	set<FaceData *> suspects;
	int insertedPoints = 0;
	for (unsigned int i = 0; i < faces->size(); i++) {
		suspects.insert((*faces)[i]);
		while (suspects.size() != 0) {
			f = *suspects.begin();
			while (f->getRho() > maxRho && !f->inSortedPoints(skipPoints) && fixSliver(f, force)) {
				insertedPoints++;
				if (!(insertedPoints % 100000))
					cout << "Fixing slivers: Inserted " << insertedPoints << " points\n";
				suspects.insert(editedFaces->begin(), editedFaces->end());
				updateHeap();

			}
			suspects.erase(f);
		}
	}
	return insertedPoints;
}

//Fix a face, assuming it is a sliver
//Insert a point at the circumcenter of the face
bool Subdivision::fixSliver(FaceData * f, bool force) {
	double x = 0, y = 0;
	f->getCircumcenter(&x, &y);
	//Put our point in the boudary
	x = (x < leftEdge ? leftEdge : (x > rightEdge ? rightEdge : x));
	y = (y < bottomEdge ? bottomEdge : (y > topEdge ? topEdge : y));
	bool inserted = insertSite(x, y);
	if (!inserted) {
		//We didn't insert a point, try a bit harder, by allowing freedom on the edge
		if (onMeshBoundary(f->longEdge())) {
			f->longEdgeBis(&x, &y);
			inserted = insertSite(x, y, f->longEdge(), force, true);
		} else if (force) {
			//Try to force the point in
			inserted = insertSite(x, y, NULL, force);
		}
	}

	return inserted;
}

//Do a smoothing iteration
//The z value of every point is set to the average of itself, and any point it
//shares a face with

//the paramater whichPoints specifies which points to smooth
//0 indicates smooth everything
//1 indicates only the manufactured data
//Do smoothing iterations times, skip the points given in skip points
void Subdivision::smooth(int whichPoints, int iterations, vector<Point3d *> * skipPoints) {
	updatePointIndices();

	struct LessThanY {
		bool operator()(Point2d * p1, Point2d * p2) {
			return p1->getY() < p2->getY();
		}
	} lty;

	set<int> skipPointIndices;

	if (skipPoints) {
		sort(skipPoints->begin(), skipPoints->end(), lty);
		for (vector<PointData*>::iterator it = points->begin(); it != points->end(); it++) {
			if (pointInSortedVector(*it, (vector<Point2d *> *) skipPoints)) {
				skipPointIndices.insert((*it)->getIndex());
			}
		}
	}

	vector<double> zvals(points->size() + 1);

	zvals[0] = 0;
	vector<bool> exists(points->size() + 1, true);
	//Indicies in the mesh start at 1, so a point with index 0 is not in the mesh
	exists[0] = false;
	for (vector<PointData*>::iterator it = points->begin(); it != points->end(); it++) {
		if (!this->exists(*it)) {
			exists[(*it)->getIndex()] = false;
		}
	}

	for (int iter = 0; iter < iterations; iter++) {
		double varience = getVariance();

		for (vector<double>::iterator it = zvals.begin(); it != zvals.end(); ++it)
			*it = INFINITY;

		//Find all the new z-values, iterating over faces
		for (vector<FaceData*>::iterator it = faces->begin(); it != faces->end(); ++it) {
			Edge * faceE = (*it)->getAnchor()->lPrev();

			for (int i = 0; i < 3; faceE = faceE->lNext(), i++) {
				PointData * p = (PointData *) faceE->org();

				if (zvals[p->getIndex()] == INFINITY && ((!whichPoints && exists[p->getIndex()]) || (whichPoints == 1 && !exists[p->getIndex()]))) {
					//We have not yet updated this point, and we are supposed to
					int neighbors = 1;
					double avez;

					avez = p->getZ();

					//if(raw->exists((PointData *) faceE->dest())) {
					if (exists[((PointData *) faceE->dest())->getIndex()] || 
						(!exists[p->getIndex()] && ((PointData *) faceE->dest())->getIndex())) {
						//The neighboring point exists in the raw field, or this point does not exist,
						//in which case we just want it to look nice
						neighbors++;
						avez += ((PointData *) faceE->dest())->getZ();
					}

					for (Edge * eit = faceE->oNext(); eit != faceE; eit = eit->oNext()) {
						if (exists[((PointData *) eit->dest())->getIndex()] 
						|| (!exists[p->getIndex()] && ((PointData *) eit->dest())->getIndex())) {
							//if(raw->exists((PointData *) eit->dest())) {
							neighbors++;
							avez += ((PointData *) eit->dest())->getZ();
						}
					}

					avez /= neighbors;
					zvals[p->getIndex()] = avez;
				}
			}
		}

		int i = 1;
		//Replace all the z-values that we are supposed to replace
		for (vector<PointData*>::iterator it = points->begin(); it != points->end(); i++, ++it) {
			PointData * p = *it;
			//Don't change the values of points we are told not to
			if (((!whichPoints && exists[i]) || (whichPoints == 1 && !exists[i])) && skipPointIndices.find(i) == skipPointIndices.end()) {
				p->setZ(zvals[i]);
			}
		}

		//Try to keep true to the original scale by keepint the varience relative to the fit plane
		//the same. Only do this if we are smoothing real data
		if (!whichPoints)
			this->stretch(1, 1, pow(varience / this->getVariance(), 0.5));
	}
}

//Removes spikes from the mesh. A spike is defined as a point which has a y-value
//greater than or less than all of its neighboring points, and a dihedreal angle to the x-y plane of
//greater than threshold on all faces connected to the point, or 1 face with diheral angle greater than
//maxAngle
int Subdivision::removeSpikes(double threshold, double maxAngle) {

	//Change the threshold to reflect the fact that the internal representation of the mesh
	//does not have the scale in it
	threshold = atan( tan(threshold * PI / 180) * scale) * 180 / PI;

	int removedSpikes = 1;
	updatePointIndices();

	vector<double> zvals(points->size() + 1, INFINITY);
	vector<bool> exists(points->size() + 1, true);
	//Indicies in the mesh start at 1, so a point with index 0 is not in the mesh
	exists[0] = false;
	for (vector<PointData*>::iterator it = points->begin(); it != points->end(); it++) {
		if (!this->exists(*it)) {
			exists[(*it)->getIndex()] = false;
		}
	}

	vector<bool> possibleSpike(exists.begin(), exists.end());
	vector<PointData*> neighbors;

	while (removedSpikes != 0) {
		removedSpikes = 0;
		//Find all the new z-values, iterating over faces
		for (vector<FaceData*>::iterator it = faces->begin(); it != faces->end(); ++it) {
			Edge * faceE = (*it)->getAnchor()->lPrev();

			for (int i = 0; i < 3; faceE = faceE->lNext(), i++) {
				PointData * p = (PointData *) faceE->org();

				if (possibleSpike[p->getIndex()]) {
					bool greater;
					bool maxFace = false;
					double avez;

					double z = p->getZ();

					//Find out if if we are checking if this point is greater than or less than all points
					//around it

					neighbors.clear();
					avez = 0;
					if(exists[((PointData *) faceE->dest())->getIndex()]) {
						greater = (z > ((PointData *) faceE->dest())->getZ());
						if (isValid((FaceData *) faceE->lDual()) && ((FaceData *) faceE->lDual())->dihedralAngle() < threshold) {
							zvals[p->getIndex()] = -INFINITY;
							continue;
						} else if (isValid((FaceData *) faceE->lDual()) && ((FaceData *) faceE->lDual())->dihedralAngle() > maxAngle) {
							maxFace = true;
						}
						avez = ((PointData *) faceE->dest())->getZ();
						neighbors.push_back((PointData *) faceE->dest());
					}

					for (Edge * eit = faceE->oNext(); eit != faceE; eit = eit->oNext()) {
						if(exists[((PointData *) eit->dest())->getIndex()]) {
							if (!neighbors.size())
								//We just found the first valid point
								greater = (z > ((PointData *) faceE->dest())->getZ());

							neighbors.push_back((PointData *) eit->dest());
							avez += ((PointData *) eit->dest())->getZ();

							if (greater && z <= ((PointData *) eit->dest())->getZ()) {
								//The point is not greater than all its neighbors
								possibleSpike[p->getIndex()] = false;
								break;
							} else if (!greater && z >= ((PointData *) eit->dest())->getZ()) {
								//The point is not less than all its neighbors
								possibleSpike[p->getIndex()] = false;
								break;
							}

							double dihedralAngle = (isValid((FaceData *) eit->lDual()) ? ((FaceData *) eit->lDual())->dihedralAngle() : -1);
							if (isValid((FaceData *) eit->lDual()) && !maxFace && dihedralAngle < threshold) {
								possibleSpike[p->getIndex()] = false;
								//break;
							} else if (isValid((FaceData *) eit->lDual()) && dihedralAngle > maxAngle) {
								possibleSpike[p->getIndex()] = true;
								maxFace = true;
							}
						}
					}

					if (possibleSpike[p->getIndex()]) {
						//we found a spike, fix this one, and mark all neighbors as possible spikes
						avez /= neighbors.size();
						p->setZ(avez);
						removedSpikes++;
						for (vector<PointData*>::iterator it = neighbors.begin(); it != neighbors.end(); ++it) {
							possibleSpike[(*it)->getIndex()] = true;
						}
					}
				}
			}
		}
		if (removedSpikes)
			cout << "Removed " << removedSpikes << " spikes\n";
	}
	return 0;
}

//Update all the point indicies in the list, to insure they are correct
//They can become offset if a point is deleted
void Subdivision::updatePointIndices() {
	int i = 1;
	for (vector<PointData*>::iterator it = points->begin(); it != points->end(); ++it) {
		PointData * p = *it;
		p->setIndex(i);
		i++;
	}
}

//Return a metric representing the difference between this mesh and the input mesh
//This could use some work, and is currently not in active use
double Subdivision::difference(Subdivision & mesh) {
	PointData * bpcopy[4];

	//Copy the bounding points of this mesh, and miake them the
	//bounding points of the new, difference mesh
	bpcopy[0] = new PointData(*boundingPoints[0]);
	bpcopy[1] = new PointData(*boundingPoints[1]);
	bpcopy[2] = new PointData(*boundingPoints[2]);
	bpcopy[3] = new PointData(*boundingPoints[3]);

	//Set the Z corrdinate to the value in this mesh minus the value in the supplied mesh
	bpcopy[0]->setZ(mesh.vertDistTo(*bpcopy[0]));
	bpcopy[1]->setZ(mesh.vertDistTo(*bpcopy[1]));
	bpcopy[2]->setZ(mesh.vertDistTo(*bpcopy[2]));
	bpcopy[3]->setZ(mesh.vertDistTo(*bpcopy[3]));

	Subdivision * tempMesh = new Subdivision(*bpcopy[0], *bpcopy[1], *bpcopy[2], *bpcopy[3], "temp.obj");
	double thisVar = 0, otherVar = 0;

	//Insert all points in this mesh which overlap with the supplied mesh
	for(vector<PointData*>::iterator it = this->points->begin()
		;it != this->points->end(); ++it) {
			PointData * p = *it;
			if (mesh.inMesh(p))
				tempMesh->insertSite(p->getX(), p->getY(), p->getZ());
	}

	//Delete the external points if they are not really in this mesh, or if they have
	//infinite distance to the mesh we are comparing to
	for (int i = 0; i < 4; i++) {
		if (find(points->begin(), points->end(), boundingPoints[i]) == points->end() ||
			fabs(bpcopy[i]->getZ()) == INFINITY) {
				//If the delete fails, it means we ran out of polygons in the overlap mesh
				//These two meshes have no significant overlap

				if (!tempMesh->deleteExternal(tempMesh->boundingPoints[i])) {
					cout << "No significant overlap\n";
					return -1;
				}
		}
	}

	thisVar = tempMesh->getVariance();

	delete tempMesh;
	tempMesh = new Subdivision(*bpcopy[0], *bpcopy[1], *bpcopy[2], *bpcopy[3]);

	//Insert all points in the supplied mesh that overlap with this mesh
	for (vector<PointData*>::iterator it = mesh.points->begin();
		it != mesh.points->end();
		++it) {
			PointData * p = *it;
			if (this->inMesh(p))
				tempMesh->insertSite(p->getX(), p->getY(), p->getZ());

	}

	//Delete the external points if they are not really in this mesh, or if they have
	//infinite distance to the mesh we are comparing to
	for (int i = 0; i < 4; i++) {
		if (find(points->begin(), points->end(), boundingPoints[i]) == points->end() ||
			fabs(bpcopy[i]->getZ()) == INFINITY) {
				//If the delete fails, it means we ran out of polygons in the overlap mesh
				//These two meshes have no significant overlap
				if (!tempMesh->deleteExternal(tempMesh->boundingPoints[i])) {
					cout << "No significant overlap\n";
					return -1;
				}
		}
	}

	otherVar = tempMesh->getVariance();

	mesh.stretch(1, 1, pow(thisVar / otherVar, 0.5));
	delete tempMesh;

	//The difference between the two meshes := this - mesh
	Subdivision difMesh(*bpcopy[0], *bpcopy[1], *bpcopy[2], *bpcopy[3], "difference.obj");

	//Insert all points in this mesh which overlap with the supplied mesh
	for(vector<PointData*>::iterator it = this->points->begin()
		;it != this->points->end(); ++it) {
			PointData * p = *it;
			double z = mesh.vertDistTo(*p);
			if (fabs(z) != INFINITY)
				difMesh.insertSite(p->getX(), p->getY(), z);
	}

	//Insert all points in the supplied mesh that overlap with this mesh
	for (vector<PointData*>::iterator it = mesh.points->begin();
		it != mesh.points->end();
		++it) {
			PointData * p = *it;
			double z = -(this->vertDistTo(*p));
			if (fabs(z) != INFINITY)
				difMesh.insertSite(p->getX(), p->getY(), z);

	}

	//Delete the external points if they are not really in this mesh, or if they have
	//infinite distance to the mesh we are comparing to
	for (int i = 0; i < 4; i++) {
		if (find(points->begin(), points->end(), boundingPoints[i]) == points->end() ||
			fabs(bpcopy[i]->getZ()) == INFINITY) {
				//If the delete fails, it means we ran out of polygons in the overlap mesh
				//These two meshes have no significant overlap
				if (!difMesh.deleteExternal(difMesh.boundingPoints[i])) {
					cout << "No significant overlap\n";
					return -1;
				}
		}
	}


	//Return the mesh back to its original state
	mesh.stretch(1, 1, pow(otherVar / thisVar, 0.5));

	//difMesh now represents this - mesh. Do some calculations
	double volVar = difMesh.getVariance();
	cout << "R: " << volVar / thisVar << "\n";
	//double a = 0, b = 0, c = 0;
	//difMesh.getFitPlane(&a, &b, &c);
	//cout << "Fit Plane: " << a << " " << b << " " << c << "\n";
	//difMesh.subtractPlane(a, b, c);
	//difMesh.outputPoints();
	//difMesh.outputFaces();

	delete bpcopy[0];
	delete bpcopy[1];
	delete bpcopy[2];
	delete bpcopy[3];

	return volVar / thisVar;
}

//Move the entire subdivision the specified amount
void Subdivision::move(double x, double y, double z) {

	for (vector<PointData*>::iterator it = points->begin(); it != points->end(); ++it) {
		PointData * p = *it;
		p->move(x, y, z);
	}

	for (int i = 0; i < 4; i++) {
		if (find(points->begin(), points->end(), boundingPoints[i]) == points->end())
			boundingPoints[i]->move(x, y, z);
	}
}

//Returns the geometric center of all the points in the mesh
void Subdivision::findCenter(double * avex, double * avey, double * avez) {
	*avex = 0;
	*avey = 0;
	*avez = 0;

	int numPoints = 0;

	for (vector<PointData*>::iterator it = points->begin(); it != points->end(); ++it) {
		PointData * p = *it;
		*avex += p->getX();
		*avey += p->getY();
		*avez += p->getZ();
		numPoints++;
	}

	*avex /= numPoints;
	*avey /= numPoints;
	*avez /= numPoints;
}


//Scale the entire mesh from its center of mass
void Subdivision::stretch(double sx, double sy, double sz) {

	double avex = 0;
	double avey = 0;
	double avez = 0;

	findCenter(&avex, &avey, &avez);


	for (vector<PointData*>::iterator it = points->begin(); it != points->end(); ++it) {
		PointData * p = *it;

		p->setX(sx * p->getX() - avex * (sx - 1));
		p->setY(sy * p->getY() - avey * (sy - 1));
		p->setZ(sz * p->getZ() - avez * (sz - 1));
	}

	for (int i = 0; i < 4; i++) {
		if (find(points->begin(), points->end(), boundingPoints[i]) == points->end()) {
			boundingPoints[i]->setX(sx * boundingPoints[i]->getX() - avex * (sx - 1));
			boundingPoints[i]->setY(sy * boundingPoints[i]->getY() - avey * (sy - 1));
			boundingPoints[i]->setZ(sz * boundingPoints[i]->getZ() - avez * (sz - 1));
		}
	}
}

//Return the z-varience in the mesh, relative to the fit plane, weighted by area of triangle
double Subdivision::getVariance() {
	double a = 0, b = 0, c = 0;
	getFitPlane(&a, &b, &c);
	subtractPlane(a, b, c);

	double totArea = 0;
	double volVar = 0;

	for (vector<FaceData*>::iterator it = faces->begin(); it != faces->end(); ++it) {
		totArea += (*it)->area();
		volVar += pow((*it)->vol(), 2) / (*it)->area();
	}

	subtractPlane(-a, -b, -c);
	return volVar / totArea;
}

//Calculate the fit plane for all points in this subdivision
//So this using standard linear regression: Z = X.B + e
//Z is a column vector of all z values
//X is a nx3 vector with columns (1 xi yi)
//B is a 3x1 column vector of the best fit z = c + a*x + b*y
//e is a column vector of errors
//By standard regression, if x^t is the transpose of X:
//B = (X^tX)^-1 . X^t . Z
void Subdivision::getFitPlane(double * a, double * b, double * c) {
	double xtx[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
	double xtxinv[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};

	//The answer we are looking for {c, a, b}
	double beta[3] = {0, 0, 0};

	double determinant = 0;
	double detinv;
	//Calculate x^t . x
	//This is a symmetric 3x3 matrix
	for(vector<PointData*>::iterator it = points->begin(); it != points->end(); ++it) {
		PointData * p = *it;

		double x = p->getX();
		double y = p->getY();

		//xtx[0][0] is the number of elements
		xtx[0][0]++;
		//xtx[1][0] is the sum of the x values
		xtx[1][0] += x;
		//xtx[2][0] is the sum of the y values
		xtx[2][0] += y;
		//xtx[1][1] is the sum of the squares of the x values
		xtx[1][1] += x * x;
		//xtx[2][2] is the sum of the squares of the y values
		xtx[2][2] += y * y;
		//xtx[2][1] is the sum of the products of x and y
		xtx[2][1] += x * y;
	}

	//Make xtx symmetric
	xtx[0][1] = xtx[1][0];
	xtx[0][2] = xtx[2][0];
	xtx[1][2] = xtx[2][1];

	//invert xtx

	//calculate the determinant, using the permutations method
	determinant += xtx[0][0] * xtx[1][1] * xtx[2][2];
	determinant -= xtx[0][0] * xtx[2][1] * xtx[1][2];
	determinant -= xtx[1][0] * xtx[0][1] * xtx[2][2];
	determinant += xtx[1][0] * xtx[2][1] * xtx[0][2];
	determinant += xtx[2][0] * xtx[0][1] * xtx[1][2];
	determinant -= xtx[2][0] * xtx[1][1] * xtx[0][2];
	if (determinant == 0) {
		//This is not possible, as the plane can't be vertical in a 2.5D surface
		cout << "Fit plane is vertical, and can't be calculated\n";
		exit(1);
	}

	detinv = 1/determinant;

	xtxinv[0][0] = detinv * (xtx[1][1] * xtx[2][2] - xtx[2][1] * xtx[1][2]);
	xtxinv[1][0] = detinv * (xtx[1][2] * xtx[2][0] - xtx[2][2] * xtx[1][0]);
	xtxinv[2][0] = detinv * (xtx[1][0] * xtx[2][1] - xtx[2][0] * xtx[1][1]);
	xtxinv[0][1] = detinv * (xtx[0][2] * xtx[2][1] - xtx[2][2] * xtx[0][1]);
	xtxinv[1][1] = detinv * (xtx[0][0] * xtx[2][2] - xtx[2][0] * xtx[0][2]);
	xtxinv[2][1] = detinv * (xtx[0][1] * xtx[2][0] - xtx[2][1] * xtx[0][0]);
	xtxinv[0][2] = detinv * (xtx[0][1] * xtx[1][2] - xtx[1][1] * xtx[0][2]);
	xtxinv[1][2] = detinv * (xtx[0][2] * xtx[1][0] - xtx[1][2] * xtx[0][0]);
	xtxinv[2][2] = detinv * (xtx[0][0] * xtx[1][1] - xtx[1][0] * xtx[0][1]);


	for(vector<PointData*>::iterator it = points->begin(); it != points->end(); ++it) {
		PointData * p = *it;

		double x = p->getX();
		double y = p->getY();
		double z = p->getZ();

		beta[0] += (xtxinv[0][0] + xtxinv[1][0] * x + xtxinv[2][0] * y) * z;
		beta[1] += (xtxinv[0][1] + xtxinv[1][1] * x + xtxinv[2][1] * y) * z;
		beta[2] += (xtxinv[0][2] + xtxinv[1][2] * x + xtxinv[2][2] * y) * z;
	}

	*a = beta[1];
	*b = beta[2];
	*c = beta[0];
}

//Subtract the given plane from the mesh
void Subdivision::subtractPlane(double a, double b, double c) {
	for (vector<PointData*>::iterator it = points->begin(); it != points->end(); ++it) {
		PointData * p = *it;
		p->move(0, 0, - (a * p->getX() + b * p->getY() + c));
	}

	for (int i = 0; i < 4; i++) {
		if (find(points->begin(), points->end(), boundingPoints[i]) == points->end())
			boundingPoints[i]->move(0, 0, - (a * boundingPoints[i]->getX() + b * boundingPoints[i]->getY() + c));
	}
}


//Order the points in z-order
void Subdivision::zOrderPoints() {
	QuadTree tree(leftEdge, bottomEdge, width, height);
	for (vector<PointData*>::iterator it = points->begin(); it != points->end(); ++it) {
		tree.addPoint(*it);
	}

	points->clear();
	tree.zOrder((vector<Point2d *> *)points);
	updatePointIndices();
}

//Sort the points by nearest neighbor
void Subdivision::nnOrderPoints() {
	nnSort((vector<Point2d *> *) points);
}

//Print out all the faces for the mesh
//This should only be called once the meshing is done, and once outputPoints has
//been called
//Faces containing only points in the pointGroup are not output, allowing holes in the mesh
//to be defined according to thier boundaries. 
void Subdivision::outputFaces(vector<Point3d *> * pointGroup, bool fakeSurface) {

	struct LessThanY {
		bool operator()(Point2d * p1, Point2d * p2) {
			return p1->getY() < p2->getY();
		}
	} lty;

	vector<bool> pointExists(points->size() + 1, true);

	if (fakeSurface) {
		for (vector<PointData*>::iterator it = points->begin(); it != points->end(); it++) {
			pointExists[(*it)->getIndex()] = exists(*it);
		}
	}

	FaceData * f;

	if (!out.is_open()) {
		cout << "Output file not open";
		exit(1);
	}

	//A vector to hold the faces we have made up
	vector<FaceData*> fakeFaces;
	set<int> pointGroupIndices;


	vector<Point3d *>::iterator curPoint;
	if (pointGroup) {
		sort(pointGroup->begin(), pointGroup->end(), lty);
		for (vector<PointData*>::iterator it = points->begin(); it != points->end(); it++) {
			if (pointInSortedVector(*it, (vector<Point2d *> *) pointGroup)) {
				pointGroupIndices.insert((*it)->getIndex());
			}
		}
	}

	out << "usemtl Default\n";
	for(vector<FaceData*>::iterator it = faces->begin();
		it != faces->end();
		++it) {

			//Does this face actaully exist in the mesh
			bool exists = true;
			bool pointGroupFace = pointGroup ? true : false;
			f = *it;

			if (pointGroupFace && pointGroupIndices.find(((PointData *) f->getAnchor()->org())->getIndex()) == pointGroupIndices.end()) {
				pointGroupFace = false;
			}

			if (fakeSurface && !pointExists[((PointData *) f->getAnchor()->org())->getIndex()]) {
				exists = false;
			}

			for (Edge * e = f->getAnchor()->lPrev(); e != f->getAnchor(); e = e->lPrev()) {
				if (pointGroupFace && pointGroupIndices.find(((PointData *) e->org())->getIndex()) == pointGroupIndices.end()) {
					pointGroupFace = false;
				}
				if (fakeSurface && !pointExists[((PointData *) e->org())->getIndex()]) {
					exists = false;
				}
			}

			if (exists && !pointGroupFace) {
				out << "f " << f->getAnchor()->org2d().getIndex() << "/" << f->getAnchor()->org2d().getIndex();
				for (Edge * e = f->getAnchor()->lPrev(); e != f->getAnchor(); e = e->lPrev()) {
					out << " " << e->org2d().getIndex() << "/" << e->org2d().getIndex();
				}
				out << "\n";
			} else if (!exists && !pointGroupFace) {
				fakeFaces.push_back(f);
			}
	}

	if (fakeFaces.size()) {
		out << "usemtl Fake\n";
		for (vector<FaceData*>::iterator it = fakeFaces.begin(); it != fakeFaces.end(); ++it) {
			f = *it;
			out << "f " << f->getAnchor()->org2d().getIndex() << "/" << f->getAnchor()->org2d().getIndex();
			for (Edge * e = f->getAnchor()->lPrev(); e != f->getAnchor(); e = e->lPrev()) {
				out << " " << e->org2d().getIndex() << "/" << e->org2d().getIndex();
			}
			out << "\n";
		}
	}
}


//Print out all points to the output file
//This should only be called after the meshing is done
void Subdivision::outputPoints() {
	if (!out.is_open()) {
		out.open(outputFile, ios::out);
		if (!out.is_open()) {
			cout << "ERROR: unable to open output file\n";
			cout << "Check file permissions\n";
			exit(1);
		}
	}

	PointData * p;

	int outPoints = 0;

	out.setf(ios::scientific);
	out.precision(8);

	//Write points out to a file, store the index to which they have been written
	for (vector<PointData*>::iterator it = points->begin();
		it != points->end(); ++it) {
			outPoints++;
			p = *it;
			out << "v " << p->getX() * scale << " " << p->getZ() << " " << p->getY() * scale << "\n";
			p->setIndex(outPoints);
	}
}

//Output a UV map, based on a top down projection
void Subdivision::outputUV() {
	PointData * p;

	if (!out.is_open()) {
		cout << "Output file not open";
		exit(1);
	}

	if (intype == RASTER) {
		for (vector<PointData*>::iterator it = points->begin();
			it != points->end(); ++it) {
				p = *it;
				out << "vt " << p->getX() / (raw->width  - 1) << " " << 1 - p->getY() / (raw->height - 1) << "\n";
		}
	} else {
		cout << "Unable to write UV map for non-raster input\n";
	}
}

//Returns true iff the specified point is in the valid region of the mesh
bool Subdivision::inMesh(double x, double y) {
	//Do a quick check as to whether the point is in the mesh
	if (!inBoundary(x, y)) {
		return false;
	}

	utilP3d->setX(x);
	utilP3d->setY(y);
	//Find the point

	Edge * e = locate(*utilP3d);

	if (e == NULL || e->lDual() == DELETED) { 
		//The point is outside the mesh
		return false;
	}

	return true;
}

//Return the vertical (z) difference to the specified point
//Note that an absolute value is not taken, so negative values will be returned
//if the point is below the mesh
double Subdivision::vertDistTo(const double x, const double y, const double z) {

	//Do a quick check as to whether the point is in the mesh
	if (!inBoundary(x, y)) {
		return INFINITY;
	}

	utilP3d->setX(x);
	utilP3d->setY(y);
	//Find the point
	Edge * e = locate(*utilP3d);

	if (e == NULL || e->lDual() == DELETED) { 
		//The point is outside the mesh
		return INFINITY;
	}

	double a, b, c;
	((FaceData *) e->lDual())->getPlane(&a, &b, &c);

	return z - a * x - b * y - c;
}

//return the error scale based on the importance map, if it exists
inline double Subdivision::errorScale(double x, double y) const {
	return impMap ? impMap->getZ((x - leftEdge) / width * (impMap->width - 1), 
		(y - bottomEdge) / height * (impMap->height - 1)) : 1;
}

//The deconstructor
//We have been careful to keep track of all points, faces, and edges that we have
//created, as well as insure they are not used by any other class
//delete them all
Subdivision::~Subdivision() {

	out.flush();
	out.close();

	//Remove the bounding points from the vector, if they exist
	vector<PointData *>::iterator it;
	it = find(points->begin(), points->end(), boundingPoints[0]);
	if (it != points->end())
		points->erase(it);
	it = find(points->begin(), points->end(), boundingPoints[1]);
	if (it != points->end())
		points->erase(it);
	it = find(points->begin(), points->end(), boundingPoints[2]);
	if (it != points->end())
		points->erase(it);
	it = find(points->begin(), points->end(), boundingPoints[3]);
	if (it != points->end())
		points->erase(it);

	//Delete the boundingPoints
	delete boundingPoints[0];
	delete boundingPoints[1];
	delete boundingPoints[2];
	delete boundingPoints[3];

	//Delete all points in the mesh
	for (vector<PointData*>::iterator it = points->begin(); it != points->end(); ++it) {
		PointData * p = *it;
		delete p;
	}

	//Delete all the edges
	delete edges;

	delete faces;
	delete points;
	delete faceHeap;
	delete editedFaces;
	delete allocFaces;
}
