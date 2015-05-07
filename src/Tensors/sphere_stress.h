/*
 * sphere_tensor.h
 *
 *  Created on: Mar 31, 2015
 *      Author: olga
 */

#ifndef SPHERE_TENSOR_H_
#define SPHERE_TENSOR_H_

#include "../structures.h"

void createSphereStressTensor();
void initSphereStress();
inline void computeSphereStress();
inline void readSphereStress();
inline void writeSphereStress();
inline void printPDBSphereStress();
inline void printDATSphereStress();
inline void destroySphereStressTensor();

Tensors sphereStressTensor;

Vector* atomSphereStress; //x - lateral, y - shear, z - radial

Vector ePhi, eTheta, eRad;

char outSphereStressFilename[256];
FILE* outSphereStress;

char pdbSphereStressFilename[256];
FILE* pdbSphereStress;

char tnsrSphereStressFilename[256];
FILE* tnsrSphereStress;

double sphereScale;

#endif /* SPHERE_TENSOR_H_ */
