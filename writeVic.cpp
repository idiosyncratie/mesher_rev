#include "vicDefs.h"

//Code to write out vicar files
//not in active use

#pragma warning(disable : 4996)

using namespace std;
void openVic(char * filename, int width, int height, char * type, ofstream & out) {
	if (out.is_open()) {
		cout << "ERROR: cannot open an open stream\n";
		exit(1);
	}

	out.open(filename, ios::out | ios::binary);

	char label[500] = "";
	char endLabel[500] = "";
	//The length of everything except the string which specifies the length
	int prelimLength = 0;
	int length;

	sprintf(endLabel, " FORMAT='%s' NL=%d NS=%d\n", type, height, width);
	prelimLength = strlen(endLabel) + 8;
	length = prelimLength + (int) log10((double) prelimLength) + 1;
	

	if ((int) log10((double) prelimLength) != (int) log10((double) length)) {
		//length has one more digit than prelimLength does
		length++;
	}

	sprintf(label, "LBLSIZE=%d", length);
	strcat(label, endLabel);
	out << label;
}

void writeVic(RawField * raw, char * filename) {
	ofstream out;
	openVic(filename, raw->width, raw->height, "DOUB", out);

	for (int y = 0; y < raw->height; y++) {
		for (int x = 0; x < raw->width; x++) {
			double val = raw->getZ(x, y);
			out.write((char *) &val, sizeof(double));
		}
	}
	out.close();
}
