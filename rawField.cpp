#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include "rawField.h"
#include <cmath>
#include <algorithm>
#include <string.h>


using namespace std;

//Don't tell me I'm using unsafe c functions, don't tell me I'm converting int to double
#pragma warning(disable:4996 4244)

RawField::RawField(const char * infile, const intype_t type, const bool seqAccess, bool zeqnodata, double replaceVal, double replaceWith, double imageMem) : seqAccess(seqAccess), zeqnodata (zeqnodata), replaceVal(replaceVal), replaceWith(zeqnodata ? replaceWith : numeric_limits<double>::infinity()), imageMem((unsigned int) imageMem) {

	this->type = type;
	//Make sure the view starts at a spot greater than anywhere in the file
	inMemory = false;
	switchEndian = false;

	int i = 0;
	int rowsPerStrip = 0;

	int samplesPerPixel = 0;
	this->height = 0;
	this->width = 0;
	this->wordSize = 0;
	int format = 0;


	switch ( type ) {
		case GEOTIFF:
#ifdef TIFF_SUPPORT
			input = TIFFOpen(infile, "r");
			if (!input) {
				cout << "Unable to open tiff file: " << infile << "\n";
				exit(1);
			}

			//Initialize variables
			wordFormat = UNKNOWN;
			wordSize = 0;
			width = 0;
			height = 0;
			rowsPerStrip = 0;

			TIFFGetField((TIFF*) input, TIFFTAG_IMAGEWIDTH, &width);
			TIFFGetField((TIFF*) input, TIFFTAG_IMAGELENGTH, &height);
			TIFFGetField((TIFF*) input, TIFFTAG_SAMPLEFORMAT, &format);
			TIFFGetField((TIFF*) input, TIFFTAG_BITSPERSAMPLE, &wordSize);
			TIFFGetField((TIFF*) input, TIFFTAG_ROWSPERSTRIP, &rowsPerStrip);
			TIFFGetField((TIFF*) input, TIFFTAG_SAMPLESPERPIXEL, &samplesPerPixel);

			wordSize /= 8;

			if (samplesPerPixel != 1) {
				cout << "ERROR: image must be grayscale, given " <<  samplesPerPixel << " samples per pixel\n";
				exit(1);
			}

			//Safely convert from libtiffs format enum to ours
			if (format == 0 || format == SAMPLEFORMAT_UINT) {
				wordFormat = FORMAT_UINT;
			} else if (format == SAMPLEFORMAT_INT) {
				wordFormat = FORMAT_INT;
			} else if (format == SAMPLEFORMAT_IEEEFP) {
				wordFormat = FORMAT_IEEEFP;
			} else {
				cout << "Unsupported TIFF format " << format << "\n";
				exit(1);
			}

			if (!rowsPerStrip) {
				cout << "ERROR: read 0 rows per strip\n";
				cout << "Perhaps the tiff image is in a tile oriented format?\n";
				cout << "Tile formats are not supported\n";
				exit(1);
			}

			numPoints = (unsigned long long) width * height;

			xmin = 0;
			ymin = 0;
			xmax = width;
			ymax = height;

			if (width * height * sizeof(double) < imageMem) {
				buf = new ImageBuffer(1, width, height, rowsPerStrip, TIFFStripSize((TIFF *) input));

				//Just store the image in memory
				cout << "Entire image fits in memory, reading it in\n";
				image = new double[width*height];
				for (int i = 0; i < width*height; i++) {
					image[i] = ReadTIFF((TIFF *) input, i % (int) width, i / (int) width, buf);

				}
				cout << "Done\n";
				delete buf;
				inMemory = true;
			} else {
				//The image does not fit in our alloted memory, make a buffer
				int maxLines = min((int) this->imageMem / TIFFStripSize((TIFF *) input), (int) height);
				buf = new ImageBuffer(maxLines, width, height, rowsPerStrip, TIFFStripSize((TIFF *) input));
				inMemory = false;
			}

			break;
#else
			cout << "ERROR, program not built with TIFF support.\n";
			cout << "Rebuild with the libtiff library and TIFF support\n";
			exit(1);
#endif
		case VICAR:
		case CUBE:
			char* view;
			in.open(infile, ios::in | ios::binary);
			if (!in.is_open()) {
				cout << "Unable to open file: " << infile << "\n";
				exit(1);
			}

			if (type == VICAR) 
				readVicHeader();
			else if (type == CUBE)
				readCubeHeader();

			if (width * height * sizeof(double) < imageMem) {
				//Just store the image in memory
				cout << "Entire image fits in memory, reading it in\n";

				inMemory = true;
				bigFile = false;
				view = new char[wordSize];
				char * bytes = new char[wordSize];
				image = new double[width*height];
				in.seekg(headersize, ios::beg);

				for (int i = 0; i < width * height; i++) {
					in.read(view, wordSize);

					if (switchEndian) {
						for (int j = 0; j < wordSize; j++) {
							bytes[j] = view[j];
						}

						for (int j = 0; j < wordSize; j++) {
							view[wordSize - 1 - j] = bytes[j];
						}
					}

					image[i] = getDouble(view, 0);
				}

				delete view;
				delete bytes;
			} else {
				in.close();
				inMemory = false;



#ifdef _WIN32
				////Use the win API file handling, so we can use files larger than 2GB
				input = CreateFile(infile, GENERIC_READ, FILE_SHARE_READ, NULL, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL);
				if (input == INVALID_HANDLE_VALUE) {
					cout << "ERROR: unable to open file " << infile << "\n";
					exit(1);
				}
				unsigned long fileSizeHigh;
				unsigned long fileSizeLow;
				if((fileSizeLow = GetFileSize(input, (LPDWORD) &fileSizeHigh)) < 1L<<31 && fileSizeHigh == 0) {
					bigFile = false;
				} else {
					bigFile = true;
				}
				fileSize = ((unsigned long long) fileSizeHigh<<32) + (unsigned long long) fileSizeLow;
#else

				////General Code Does not support files > 2GB in size
				input = new ifstream();
				((ifstream *)input)->open(infile, ios::in | ios::binary);
				((ifstream *)input)->seekg(0, ios::end);
				fileSize = ((ifstream *)input)->tellg();
				bigFile = false;
#endif

				cout << "Read File with size " << (fileSize>>20) << "MB\n";

				if (!seqAccess) {
					image = new double[1];
				} else {
					//The image does not fit in our alloted memory, make a buffer
					int maxLines = min((int) this->imageMem / width / wordSize, (int) height);
					buf = new ImageBuffer(maxLines, width, height, rowsPerStrip, width * wordSize);
				}
			}

			numPoints = (unsigned long long) width * height;
			break;
		case OBJ:
			in.open(infile, ios::in);
			if (!in.is_open()) {
				cout << "Unable to open file: " << infile << "\n";
				exit(1);
			}

			char line[500];

			double x, y, z;
			xmin = 0; xmax = 0; ymin = 0; ymax = 0, avez = 0;
			i = 0;
			for(in.getline(line, 500); !in.eof(); in.getline(line, 500)) {
				if (line[0] == 'v' && line[1] == ' ') {
					i++;
					sscanf(line, "v %lg %lg %lg\n", &x, &z, &y);
					Point3d * p = new Point3d(x, y, z);
					points.push_back(p);
					if (x < xmin)
						xmin = x;
					if (x > xmax)
						xmax = x;
					if (y < ymin)
						ymin = y;
					if (y > ymax)
						ymax = y;
					if (!(i%100000)) {
						cout << i << " points read\n";
					}
				}
			}

			cout << "Read " << i << " points\n";
			//Sort the points in acending order by y coordinate
			//The program would probably be (perhaps much) faster if this was replaced by a 
			//spiffy quadtree type sorting
			numPoints = i;
			width = xmax - xmin;
			height = ymax - ymin;
			cout << "Building kdtree...";
			tree = new KDTree((vector<Point2d *> *) &points);
			cout << "Done\n";
			break;
		default:
			//This will only occur if someone doesn't use the enum type supplied
			cout << "ERROR: unsupported file type\n";
			exit(1);
	}
}

//Read the actual raster to get the z value at the specified location
double RawField::readRaster(int xpix, int ypix) {
	switch(type) {
case GEOTIFF:
#ifdef TIFF_SUPPORT
	return ReadTIFF((TIFF*) input, xpix, ypix, buf);
#endif
	break;
case VICAR:
case CUBE:
	if(seqAccess) {
		void * line = buf->getIndex(ypix);
		// the line is not in the buffer
		if (!line) {
			line = buf->writeIndex(ypix);
			unsigned long long desPos = (unsigned long long) headersize + (unsigned long long) wordSize * ypix * width;
#ifdef _WIN32
			int charsRead;
			LARGE_INTEGER li;
			li.QuadPart = desPos;
			SetFilePointerEx(input, li, NULL, FILE_BEGIN);
			ReadFile(input, line, width * wordSize, (LPDWORD) &charsRead, NULL);
#else
			////General Code: no support for files > 2GB
			((ifstream*)input)->seekg(desPos, ios::beg);
			((ifstream*)input)->read((char *) line, width * wordSize);
#endif

			double bytes[sizeof(double)];

			if (switchEndian) {
				for (int i = 0; i < width * wordSize; i+=wordSize) {
					for (int j = 0; j < wordSize; j++) {
						bytes[j] = ((char *) line)[i + j];
					}

					for (int j = 0; j < wordSize; j++) {
						((char *) line)[i - j + wordSize - 1] = bytes[j];
					}
				}
			}
		}

		return getDouble(line, xpix);
	} else {
		////WINDOWS Code: Supports files > 2GB
#ifdef _WIN32
		if (!bigFile) {
			long seekto = headersize + wordSize * (ypix * width + xpix);
			SetFilePointer(input, seekto, NULL, FILE_BEGIN);
		} else {
			LARGE_INTEGER li;
			li.QuadPart = (long long) headersize + (long long) wordSize * (ypix * width + xpix);
			SetFilePointerEx(input, li, NULL, FILE_BEGIN);
		}
		long charsRead;
		ReadFile(input, image, wordSize, (LPDWORD) &charsRead, NULL);
#else
		////General Code, does not support files > 2GB
		((ifstream*) input)->seekg(headersize + wordSize * (ypix * width + xpix), ios::beg);
		((ifstream*) input)->read((char *) image, wordSize);
#endif
		return getDouble(image, 0);
	}
default:
	cout << "ERROR: this function is not compatitble with the input type\n";
	exit(1);
	}
}

//Get the information we need from a vicar header
void RawField::readVicHeader() {
	char * header;
	char * pc;
	char firstRead[100];
	in.get(firstRead, 100);
	sscanf(firstRead, "LBLSIZE=%d", &headersize);
	header = new char[headersize + 1];
	if (!headersize) {
		cout << "ERROR reading vicar header, headersize = 1\n";
		exit(1);
	}

	in.seekg(0, ios::beg);
	in.read(header, headersize + 1);

	height = 0;
	pc = header;
	while (height == 0) {
		//Find the first occurence of the flag
		pc = strstr(pc, "NL=");
		if (!pc) {
			//We didn't find the flag in the header
			cout << "ERROR reading vicar header, unable to find height\n";
			exit(1);
		} else {
			//Make sure the preceding character is whitespace
			pc--;
			sscanf(pc, "%*[\r\n\t ]NL=%d", &height);
			//Move past this flag
			pc+=2;
		}
	}

	width = 0;
	pc = header;
	while (width == 0) {
		//Find the first occurence of the flag
		pc = strstr(pc, "NS=");
		if (!pc) {
			//We didn't find the flag in the header
			cout << "ERROR reading vicar header, unable to find width\n";
			exit(1);
		} else {
			//Make sure the preceding character is whitespace
			pc--;
			sscanf(pc, "%*[\r\n\t ]NS=%d", &width);
			//Move past this flag
			pc+=2;
		}
	}

	pc = header;
	char format[100];
	wordSize = 0;
	while (wordSize == 0) {
		//Find the first occurence of the flag
		pc = strstr(pc, "FORMAT=");
		if (!pc) {
			//We didn't find the flag in the header
			cout << "ERROR reading vicar header, format is not DOUB or REAL, INT, or HALF\n";
			exit(1);
		} else {
			//Make sure the preceding character is whitespace
			pc--;
			sscanf(pc, "%*[\r\n\t ]FORMAT=%s", format);
			if (!strcmp(format, "'DOUB'")) {
				wordSize = sizeof(double);
				wordFormat = FORMAT_IEEEFP;
				cout << "Reading double precision floating point vicar file\n";
			} else if (!strcmp(format, "'REAL'")) {
				wordSize = sizeof(float);
				wordFormat = FORMAT_IEEEFP;
				cout << "Reading single precision floating point vicar file\n";
			//} else if (!stricmp(format, "'INT'")) {
			} else if (!strcmp(format, "'INT'")) {
				wordSize = 4;
				wordFormat = FORMAT_INT;
			//} else if (!stricmp(format, "'HALF'")) {
			} else if (!strcmp(format, "'HALF'")) {
				wordSize = 2;
				wordFormat = FORMAT_INT;
			}

			//Move past this flag
			pc+=2;
			strcpy(format, "");
		}
	}

	pc = header;
	while (pc) {
		//Find the first occurence of the flag
		pc = strstr(pc, "HOST=");
		if (!pc) {
			//do nothing
		} else {
			//Make sure the preceding character is whitespace
			pc--;
			sscanf(pc, "%*[\r\n\t ]HOST=%s", format);
			if (!strcmp(format, "'SUN-SOLR'")) {
				switchEndian = true;
				cout << "Switching endian\n";
			}
			//Move past this flag
			pc+=2;
			strcpy(format, "");
		}
	}
	delete header;
}

//Get the information we need from a cube header
void RawField::readCubeHeader() {
	char line[500];
	char format[100];
	int headRecords = 0;
	int recordSize = 0;
	char * pc;
	char * header;
	int bands;

	wordFormat = FORMAT_IEEEFP;

	while(!in.eof() && (headRecords == 0 || recordSize == 0)) {
		in.getline(line, 500);
		pc = strstr(line, "RECORD_BYTES = ");
		if (pc) {
			sscanf(pc, "RECORD_BYTES = %d", &recordSize);
		}

		pc = strstr(line, "^QUBE = ");
		if (pc) {
			sscanf(pc, "^QUBE = %d", &headRecords);
		}
	}

	if (!headRecords) {
		cout << "ERROR, unable to find QUBE pointer in file\n";
		exit(1);
	}

	if (!recordSize) {
		cout << "ERROR, unable to find RECORD_SIZE label\n";
		exit(1);
	}

	this->headersize = (headRecords - 1) * recordSize;
	header = new char[headersize + 1];
	in.seekg(0, ios::beg);
	in.read(header, headersize + 1);

	pc = header;
	while (height == 0) {
		//Find the first occurence of the flag
		pc = strstr(pc, "CORE_ITEMS = ");
		if (!pc) {
			//We didn't find the flag in the header
			cout << "ERROR reading cube header, unable to find cube dimentions\n";
			exit(1);
		} else {
			//Make sure the preceding character is whitespace
			pc--;
			sscanf(pc, "%*[\r\n\t ]CORE_ITEMS = (%d,%d,%d)", &width, &height, &bands);
			//Move past this flag
			pc+=2;
		}
	}

	if (bands != 1) {
		cout << "ERROR, a DEM must be a single banded image\n";
		cout << "Found " << bands << " bands\n";
		exit(1);
	}

	pc = header;
	while (wordSize == 0) {
		//Find the first occurence of the flag
		pc = strstr(pc, "CORE_ITEM_BYTES = ");
		if (!pc) {
			//We didn't find the flag in the header
			cout << "ERROR reading cube header, unable to find wordSize\n";
			exit(1);
		} else {
			//Make sure the preceding character is whitespace
			pc--;
			sscanf(pc, "%*[\r\n\t ]CORE_ITEM_BYTES = %d", &wordSize);
			//Move past this flag
			pc+=2;
		}
	}

	pc = header;
	while (pc) {
		//Find the first occurence of the flag
		pc = strstr(pc, "CORE_ITEM_TYPE = ");
		if (!pc) {
			//Do nothing
		} else {
			//Make sure the preceding character is whitespace
			pc--;
			sscanf(pc, "%*[\r\n\t ]CORE_ITEM_TYPE = %s", format);
			if (!strcmp(format, "SUN_REAL")) {
				switchEndian = true;
				cout << "Switching endian\n";
			}
			//Move past this flag
			pc+=2;
			strcpy(format, "");
		}
	}

	pc = header;
	while (pc) {
		//Find the first occurence of the flag
		pc = strstr(pc, "SUFFIX_ITEMS = ");
		if (!pc) {
			//Do nothing
		} else {
			//Make sure the preceding character is whitespace
			pc--;
			int a = 0, b = 0, c = 0;
			sscanf(pc, "%*[\r\n\t ]SUFFIX_ITEMS = (%d, %d, %d)", &a, &b, &c);
			if (a || b) {
				cout << "ERROR: unable to handle suffixes in cube files\n";
			}
			//Move past this flag
			pc+=2;
		}
	}
}

//Return the point in the raw field nearest to the specified point
Point3d * RawField::getPointNear(Point2d * p) {
	if (type == OBJ) {
		return (Point3d *) tree->nearestNeighbor(p);
	} else {
		cout << "Incompatible type for get point near\n";
		exit(1);
	}
}

//Return the point in the raw field nearest to the specified point
Point3d * RawField::getPointNear(double x, double y) {
	if (type == OBJ) {
		Point2d * p = new Point2d(x, y);
		Point3d * output = (Point3d *) tree->nearestNeighbor(p);
		delete p;
		return output;
	} else {
		cout << "Incompatible type for get point near\n";
		exit(1);
	}
}

//Appends the points in the specified box to the given vector
void RawField::getPointsInBox(double llx, double lly, double urx, double ury, vector<Point3d *> * vect) {
	if (this->type == OBJ)
		tree->getPointsInBox(llx, lly, urx, ury, (vector<Point2d *> *) vect);
	else {
		cout << "ERROR: getPointsInBox not implemented for non OBJ types\n";
		exit(1);
	}
}

//Return the correct value based on the sample format and size
double RawField::getDouble(void * line, int pix) {
	switch (wordFormat) {
				case FORMAT_IEEEFP:
					if (wordSize == 4) 
						return cleanPoint((double) ((float *)line)[pix]);
					else if (wordSize == 8)
						return cleanPoint(((double *)line)[pix]);

					break;
				case FORMAT_UINT:
					if (wordSize == 4)
						return cleanPoint(((uint32 *)line)[pix]);
					else if (wordSize == 2)
						return cleanPoint(((uint16 *)line)[pix]);
					else if (wordSize == 1)
						return cleanPoint(((uint8 *)line)[pix]);

					break;
				case FORMAT_INT:
					if (wordSize == 4)
						return cleanPoint(((int32 *)line)[pix]);
					else if (wordSize == 2)
						return cleanPoint(((int16 *)line)[pix]);
					else if (wordSize == 1)
						return cleanPoint(((int8 *)line)[pix]);

					break;
				default:
					cout << "ERROR: unrecognized sample type\n";
					cout << wordFormat << "\n";
					exit(1);
	}

	//Never executed
	return 0;
}

//Check to see if this point actually exits in the raw data
bool RawField::exists(Point2d * p) {
	if (type > 0) {
		//raster type, return true if the point is in the box, and has a real data value
		return (p->getX() >= 0 && p->getX() < width && p->getY() >= 0 && p->getY() < height
			&& (!zeqnodata || getZ((int) p->getX(), (int) p->getY()) != replaceWith));
	} else {
		//xyz field, return true if the nearest point to p has the same xy value
		Point3d * p3d = getPointNear(p);
		return (p->getX() == p3d->getX()) && (p->getY() == p3d->getY());
	}
}

RawField::~RawField() {
	if (type == GEOTIFF) {
#ifdef TIFF_SUPPORT
		TIFFClose((TIFF*) input);
#endif
	} else if (type == VICAR || type == CUBE) {
		////WINDOWS Code: Supports files > 2GB
#ifdef _WIN32
		CloseHandle(input);
#else
		////General Code: Does not support files >2GB
		if (!inMemory) {
			((ifstream*)input)->close();
			delete ((ifstream*)input);
		}
#endif
	}

	if (type > 0) {
		//raster
		if (inMemory) {
			delete image;
		} else {
			delete buf;
		}
	}
}


ImageBuffer::ImageBuffer(int maxLines, int width, int height, int rowsPerStrip, int lineSize) :
width(width), height(height), rowsPerStrip(rowsPerStrip), maxLines(maxLines), lineSize(lineSize) {
	lines = new void*[height];
	allocLines = new void*[maxLines];
	indices = new int[maxLines];

	curLine = 0;

	for (int i = 0; i < maxLines; i++) {
		indices[i] = -1;
		allocLines[i] = (void *) malloc(lineSize);
	}

	for (int i = 0; i < height; i++) {
		lines[i] = NULL;
	}
}

ImageBuffer::~ImageBuffer() {
	for (int i = 0; i < maxLines; i++) {
		free(allocLines[i]);
	}

	delete indices;
	delete lines;
	delete allocLines;
}

//Puts the specified index into the hashmap, returns a block of memory to
//store the line
void * ImageBuffer::writeIndex(int index) {
	if (indices[curLine] != -1) {
		lines[indices[curLine]] = NULL;
	}

	lines[index] = allocLines[curLine];
	indices[curLine] = index;

	curLine++;
	curLine%=maxLines;
	return lines[index];
}
