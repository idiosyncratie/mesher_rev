#ifndef INCLUSION_RAWFIELD_H
#define INCLUSION_RAWFIELD_H


//Comment this out if the libtiff library is not abvailable
#define TIFF_SUPPORT

#ifdef TIFF_SUPPORT
#include "tiffio.h"
#include "tiff.h"
#endif

#include "quadtree.h"
#include "geom2d.h"
#include "heap.h"
#include <fstream>
#include "kdtree.h"
#include <limits>

#ifdef _WIN32
#include <windows.h>
#endif

//The maximum value allowed to be read from an image
#define MAX_VALUE 10000000000000000000000000000000.0

#undef INFINITY

class ImageBuffer {
private:
	int * indices;
	void ** lines;
	void ** allocLines;
	int curLine;

public:
	const int width;
	const int height;
	const int rowsPerStrip;
	const int maxLines;
	const int lineSize;
	ImageBuffer(int maxLines, int width, int height, int rowsPerStrip, int lineSize);
	~ImageBuffer();

	void * getIndex(int index);
	void * writeIndex(int index);
};

//A class to represent all raw fields
//Since it has to represent many kinds of data fields, it is necessarily messy
//The idea is to put the mess here, and give a few access functions which are
//clean
//This code also has some OS dependent functions, as file handling of files
//greater than 2 GB is not standardized accross OSs.
class RawField {
public:
	enum intype_t
	{
		GEOTIFF = 1,
		VICAR = 2,
		CUBE = 3,
		OBJ = -1,
	};

	enum format_t
	{
		UNKNOWN = 0, 
		FORMAT_UINT = 1,
		FORMAT_INT = 2,
		FORMAT_IEEEFP = 3,
	};


private:
	//Data used for xyz fileds
	std::vector<Point3d*> points;

	//general data

	//The input file (of type void because it could either be a Windows HANDLE, or
	//an ifstream
	void * input;
	//Do we expect sequencial access to this field, changes caching priorities
	const bool seqAccess;
	//What value should we replace with the nodata flag
	const double replaceVal;
	//The amount of memory we are allowed to use
	const unsigned int imageMem;

	//The actuall input
	std::ifstream in;

	//Data used for raster files

	//The image buffer, containting the portion of image in memory
	ImageBuffer * buf;
	//The entire image, stored as cleaned doubles, if it fits in memory
	double * image;


	//Is the image in memory?
	bool inMemory;
	//The size of the header, used for vicar and cube files
	int headersize;
	//Is the file greater than 2GB?
	bool bigFile;
	//The number of bytes per pixel in the image
	int wordSize;
	//The format of the data in the raster (float, integer, unsigned integer etc.)
	format_t wordFormat;
	//The size of the file
	unsigned long long fileSize;
	//Should we switch the endian of the data?
	bool switchEndian;



	//Functions to read headers
	void readVicHeader();
	void readCubeHeader();

	//Read a raster file
	double readRaster(int xpix, int ypix);

	//Functions to convert from raw raster data into the true data

	//Remove pixels with values > MAX_VALUE, and replace values which are specified to be
	//the nodata flag
	double cleanPoint(double point);
	double cleanPoint(int point);
	double cleanPoint(unsigned int point);

	//Get a cleaned double from the specified pixel in the given line, assuming the line is
	//of the format wordFormat, and size wordSize.
	double getDouble(void * line, int pix);

public:
	double xmin, ymin, xmax, ymax, avez;
	KDTree * tree;
	RawField(const char * infile, const intype_t type, const bool seqAccess = true, 
		bool zeqnodata = false, double replaceVal = 0.0, double replaceWith = std::numeric_limits<double>::infinity(), double imageMem = (1<<20));
	~RawField();

	//The width and height of the image
	int width, height;
	//The type of input
	intype_t type;
	//What should we replace nodata with (infinity is a flag to ignore)
	const double replaceWith;
	//The number of points in the field
	unsigned long long numPoints;
	//Should we be replacing values with the nodata flag?
	const bool zeqnodata;


	//For a raster, do what is necessary to get the value at the desired pixel
	//This is the function which should be called from outside this class
	double getZ(int xpix, int ypix);

	void getPointsInBox(double llx, double lly, double urx, double ury, vector<Point3d *> * vect);

	Point3d * getPointNear(double x, double y);
	Point3d * getPointNear(Point2d * p);

	bool exists(Point2d * p);
	vector<Point3d *> * getPoints();
#ifdef TIFF_SUPPORT
	double ReadTIFF(TIFF* input, uint32 x_pixel, uint32 y_pixel, ImageBuffer * buffer);
#endif
};

//Returns the line at the specified index, if it exists, otherwise NULL
inline void * ImageBuffer::getIndex(int index) {
	return lines[index];
}

//Return the z value from a raster at the given location
inline double RawField::getZ(int xpix, int ypix) {
#ifdef _DEBUG
	if (xpix >= width || ypix >= height) {
		std::cout << "ERROR: asked for pixel outside of range\n";
		std::cout << "Asked for " << xpix << " " << ypix << " with box " << width << " " << height << "\n";
		exit(1);
	}
#endif
	return (inMemory ? image[ypix * width + xpix] : readRaster(xpix, ypix));
}

//Clean up the z value, detect missing datapoints, and replace them appropriately
inline double RawField::cleanPoint(double point) {
	return (fabs(point) > MAX_VALUE || (zeqnodata && point == replaceVal)) ? replaceWith : point;
}

//Clean up the z value, detect missing datapoints, and replace them appropriately
inline double RawField::cleanPoint(int point) {
	return (zeqnodata && point == replaceVal) ? replaceWith : (double) point;
}

//Clean up the z value, detect missing datapoints, and replace them appropriately
inline double RawField::cleanPoint(unsigned int point) {
	return (zeqnodata && point == replaceVal) ? replaceWith : (double) point;
}

//Return all the points in the raw field
//This should really be replaced with a call to get points in box, so that the internal vector cannot
//be edited
inline vector<Point3d *> * RawField::getPoints() {
	return &points;
}
#endif //#ifndef INCLUSION_RAWFIELD_H
