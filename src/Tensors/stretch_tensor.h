/*
 * stretch_tensor.h
 *
 *  Created on: Apr 6, 2015
 *      Author: olga
 */

#ifndef STRETCH_TENSOR_H_
#define STRETCH_TENSOR_H_

#include "../structures.h"
#include "../functions.h"

void createStretchTensor();
void initStretch();
void computeStretch();
void readStretch();
void writeStretch();
void printPDBStretch();
void printDATStretch();
void destroyStretchTensor();

extern Tensors stretchTensor;

extern FILE* tnsrStretch;

extern char outStretchFilename[256];

extern char pdbStretchFilename[256];


#endif /* STRETCH_TENSOR_H_ */
