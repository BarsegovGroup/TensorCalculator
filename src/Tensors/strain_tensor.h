/*
 * strain_tensot.h
 *
 *  Created on: Apr 2, 2015
 *      Author: olga
 */

#ifndef STRAIN_TENSOT_H_
#define STRAIN_TENSOT_H_

#include "../structures.h"
#include "../functions.h"

void createStrainTensor();
void initStrain();
void computeStrain();
void readStrain();
void writeStrain();
void printPDBStrain();
void printDATStrain();
void destroyStrainTensor();

#endif /* STRAIN_TENSOT_H_ */
