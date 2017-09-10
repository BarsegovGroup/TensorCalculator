/*
 * tensorio.h
 *
 *  Created on: Feb 12, 2013
 *      Author: olga
 */

#ifndef TENSORIO_H_
#define TENSORIO_H_

void printAtomValueToFile(FILE* file, Atom atomData, double atomEnergy);
int readAtomValueFrame(FILE* file, double* atomEnergy);

void printAtomVectorToFile(FILE* file, Atom atomData, Vector atomVector);
int readAtomVectorFrame(FILE* file, Vector* atomVector);

void printAtomTensorToFile(FILE* file, Atom atomData, double* atomTensor);
int readAtomTensorFrame(FILE* file, double** atomTensorFrame);

#endif /* TENSORIO_H_ */
