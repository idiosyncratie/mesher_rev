#ifndef INCLUSION_VICDEFS_H
#define INCLUSION_VICDEFS_H

#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include "rawField.h"

void openVic(char * filename, int width, int height, char * type, std::ofstream & out);
void writeVic(RawField * raw, char * filename);

#endif
