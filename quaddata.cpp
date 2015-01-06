#include "facedata.h"

#define POINTDATA_ALLOCSIZE 65536


/*int PointData::numPointData = 0;
PointData * PointData::block = new PointData[POINTDATA_ALLOCSIZE];

void * PointData::operator new (size_t bytes) {
	//Check that we are allocating an amount of memory which fits into our blocks
	if (bytes == sizeof(PointData)) {
		//Check to see if we have run out of allocated blocks
		if (numPointData == POINTDATA_ALLOCSIZE) {
			numPointData = 0;
			block = new PointData[POINTDATA_ALLOCSIZE];
		}

		return &(block[numPointData++]);
	} else {
		cout << "Using default new\n";
		return new char[bytes];
	}
}*/
