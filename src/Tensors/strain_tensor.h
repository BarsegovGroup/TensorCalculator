/*
 * strain_tensot.h
 *
 *  Created on: Apr 2, 2015
 *      Author: olga
 */

#ifndef STRAIN_TENSOT_H_
#define STRAIN_TENSOT_H_

#include "../structures.h"
#include "deform_gradient.h"
#include "stretch_tensor.h"

void createStrainTensor();
void initStrain();
inline void computeStrain();
inline void readStrain();
inline void writeStrain();
inline void printPDBStrain();
inline void printDATStrain();
inline void destroyStrainTensor();

Tensors strainTensor;

double** atomGreenStrain;
double** atomStretch;

char outGreenStrainFilename[256];
FILE* outGreenStrain;

char outStretchFilename[256];
FILE* outStretch;

char pdbGreenStrainFilename[256];
FILE* pdbGreenStrain;

char pdbStretchFilename[256];
FILE* pdbStretch;

char tnsrGreenStrainFilename[256];
FILE* tnsrGreenStrain;

char tnsrStretchFilename[256];
FILE* tnsrStretch;

double strainScale;

#endif /* STRAIN_TENSOT_H_ */
