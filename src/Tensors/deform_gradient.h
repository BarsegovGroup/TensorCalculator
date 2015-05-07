/*
 * deform_gradient.h
 *
 *  Created on: Apr 1, 2015
 *      Author: olga
 */

#ifndef DEFORM_GRADIENT_H_
#define DEFORM_GRADIENT_H_

#include "../structures.h"

void createDeformGradientTensor();
void initDeformGradient();
inline void computeDeformGradient();
inline void readDeformGradient();
inline void writeDeformGradient();
inline void printPDBDeformGradient();
inline void printDATDeformGradient();
inline void destroyDeformGradientTensor();

Tensors deformGradientTensor;

double** atomDeformGradient;
Vector* atomRefCoord;

char outDeformGradientFilename[256];
FILE* outDeformGradient;

char pdbDeformGradientFilename[256];
FILE* pdbDeformGradient;

char tnsrDeformGradientFilename[256];
FILE* tnsrDeformGradient;

double weightCov;
double weightNat;
double weightRep;

#endif /* DEFORM_GRADIENT_H_ */
