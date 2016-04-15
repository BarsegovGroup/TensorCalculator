/*
 * deform_gradient.h
 *
 *  Created on: Apr 1, 2015
 *      Author: olga
 */

#ifndef DEFORM_GRADIENT_H_
#define DEFORM_GRADIENT_H_

#include "../structures.h"
#include "../functions.h"

void createDeformGradientTensor();
void initDeformGradient();
void computeDeformGradient();
void readDeformGradient();
void writeDeformGradient();
void printPDBDeformGradient();
void printDATDeformGradient();
void destroyDeformGradientTensor();

extern Tensors deformGradientTensor;

extern double** atomDeformGradient;

extern FILE* tnsrDeformGradient;

extern char outDeformGradientFilename[256];

extern char pdbDeformGradientFilename[256];

#endif /* DEFORM_GRADIENT_H_ */
