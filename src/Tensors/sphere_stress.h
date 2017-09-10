/*
 * sphere_tensor.h
 *
 *  Created on: Mar 31, 2015
 *      Author: olga
 */

#ifndef SPHERE_TENSOR_H_
#define SPHERE_TENSOR_H_

#include "../structures.h"
#include "../functions.h"

void createSphereStressTensor();
void initSphereStress();
inline void computeSphereStress();
inline void readSphereStress();
inline void writeSphereStress();
inline void printPDBSphereStress();
inline void printDATSphereStress();
inline void destroySphereStressTensor();

#endif /* SPHERE_TENSOR_H_ */
