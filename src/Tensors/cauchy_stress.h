/*
 * cauchy_tensor.h
 *
 *  Created on: Mar 30, 2015
 *      Author: olga
 */

#ifndef CAUCHY_TENSOR_H_
#define CAUCHY_TENSOR_H_

#include "../structures.h"

void createCauchyStressTensor();
void initCauchyStress();
inline void computeCauchyStress();
inline void readCauchyStress();
inline void writeCauchyStress();
inline void printPDBCauchyStress();
inline void printDATCauchyStress();
inline void destroyCauchyStressTensor();

Tensors cauchyStressTensor;

double** atomCauchyStress;

char outCauchyStressFilename[256];
FILE* outCauchyStress;

char pdbCauchyStressFilename[256];
FILE* pdbCauchyStress;

char tnsrCauchyStressFilename[256];
FILE* tnsrCauchyStress;

double cauchyScale;

#endif /* CAUCHY_TENSOR_H_ */
