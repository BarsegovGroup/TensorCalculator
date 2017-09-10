/*
 * cauchy_tensor.h
 *
 *  Created on: Mar 30, 2015
 *      Author: olga
 */

#ifndef CAUCHY_TENSOR_H_
#define CAUCHY_TENSOR_H_

#include "../structures.h"
#include "../functions.h"

void createCauchyStressTensor();
void initCauchyStress();
inline void computeCauchyStress();
inline void readCauchyStress();
inline void writeCauchyStress();
inline void printPDBCauchyStress();
inline void printDATCauchyStress();
inline void destroyCauchyStressTensor();
inline void averageCauchyStressTensor();

#endif /* CAUCHY_TENSOR_H_ */
