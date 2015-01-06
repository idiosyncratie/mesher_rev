#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <ctime>
#include <limits>

#include "rawField.h"
#include "subdivision.h"


using namespace std;

#ifdef TIFF_SUPPORT
double RawField::ReadTIFF(TIFF* input, uint32 x_pixel, uint32 y_pixel, ImageBuffer * buffer) {

	int pixelPos;
	int index;

	void * line;

	//This is the core of the program (Every getZ() comes to here for TIFF input)
	//Don't do an unnessesary division and % if we don't have to
	if (buffer->rowsPerStrip != 1) {
		index = y_pixel / buffer->rowsPerStrip;
		pixelPos = (y_pixel % buffer->rowsPerStrip) * buffer->width + x_pixel;
	} else {
		index = y_pixel;
		pixelPos = x_pixel;
	}

	line = buffer->getIndex(index);

	if (!line) {
		//The line is not in our buffer
		line = buffer->writeIndex(index);
		if (TIFFReadEncodedStrip(input, index, line, -1) == -1) {
			cout << "ERROR: unable to read image strip\n";	
			exit(1);
		}
	}

	return getDouble(line, pixelPos);
}
#endif

void WriteError(char * outputName, Subdivision * mesh) {
	if (mesh->raw->type < 0) {
		cout << "Cannot output error for non raster types\n";
		return;
	}

#ifdef TIFF_SUPPORT
	TIFF * output = TIFFOpen(outputName, "w");
	TIFFSetField(output, TIFFTAG_IMAGEWIDTH, (uint32) mesh->raw->width);
	TIFFSetField(output, TIFFTAG_IMAGELENGTH, (uint32) mesh->raw->height);
	TIFFSetField(output, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
	TIFFSetField(output, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
	TIFFSetField(output, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
	TIFFSetField(output, TIFFTAG_SAMPLESPERPIXEL, 1);
	TIFFSetField(output, TIFFTAG_BITSPERSAMPLE, 32);
	TIFFSetField(output, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP);
	TIFFSetField(output, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
	TIFFSetField(output, TIFFTAG_ROWSPERSTRIP, 1);

	tdata_t * outputBuf = (tdata_t *) _TIFFmalloc(TIFFStripSize(output));

	double maxError = 0;
	double error = 0;
	double rms = 0;
	for (int curRow = 0; curRow < mesh->raw->height; curRow++) {
		for (int curCol = 0; curCol < mesh->raw->width; curCol++) {
			error = fabs(mesh->vertDistTo(curCol, curRow, mesh->raw->getZ(curCol, curRow)));
			if (error == Subdivision::INFINITY)
				error = 0;

			if (error > maxError) {
				maxError = error;
			}
			rms += error*error;
			((float *) outputBuf)[curCol] = (float) error;
		}
		TIFFWriteEncodedStrip(output, curRow, outputBuf, sizeof(float) * mesh->raw->width);
	}

	cout << "Max Error: " << maxError << " RMS: " << pow(rms / mesh->raw->height / mesh->raw->width, 0.5) << "\n";

	TIFFFlush(output);
	TIFFClose(output);
#else
	cout << "ERROR: program not built with TIFF support, cannout output error map\n";
#endif
}
