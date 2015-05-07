/*
 * energy.h
 *
 *  Created on: Mar 31, 2015
 *      Author: olga
 */

#ifndef ENERGY_H_
#define ENERGY_H_

#include "structures.h"
#include "functions.h"

void createEnergy();
void initEnergy();
inline void computeEnergy();
inline void readEnergy();
inline void writeEnergy();
inline void printPDBEnergy();
inline void printDATEnergy();
void destroyEnergy();

Tensors energy;

double* atomEnergy;

char outEnergyFilename[256];
FILE* outEnergy;

char pdbEnergyFilename[256];
FILE* pdbEnergy;

char tnsrEnergyFilename[256];
FILE* tnsrEnergy;

double energyScale;

#endif /* ENERGY_H_ */
