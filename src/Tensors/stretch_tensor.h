/*
 * stretch_tensor.h
 *
 *  Created on: Apr 6, 2015
 *      Author: olga
 */

#ifndef STRETCH_TENSOR_H_
#define STRETCH_TENSOR_H_

#include "../structures.h"

void createStretchTensor();
void initStretch();
inline void computeStretch();
inline void readStretch();
inline void writeStretch();
inline void printPDBStretch();
inline void printDATStretch();
inline void destroyStretchTensor();

Tensors stretchTensor;

double** atomStretch;

char outStretchFilename[256];
FILE* outStretch;

char pdbStretchFilename[256];
FILE* pdbStretch;

char tnsrStretchFilename[256];
FILE* tnsrStretch;

double stretchScale;

#endif /* STRETCH_TENSOR_H_ */
