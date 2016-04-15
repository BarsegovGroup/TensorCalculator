/*
 * structures.h
 *
 *  Created on: Mar 30, 2015
 *      Author: olga
 */

#ifndef STRUCTURES_H_
#define STRUCTURES_H_

//#include "IO/pdbio.h"
#include "parameters.h"
#include "IO/pdbio.h"

typedef struct{
	char name[100];
	void (*compute)();
	void (*write)();
	void (*read)();
	void (*printPDB)();
	void (*printDAT)();
	void (*destroy)();
} Tensors;

typedef struct{
	double x;
	double y;
	double z;
	double value;
} Vector;

typedef struct{
	float refDist;
	float Eh;

} FuncParam;

extern int frame;
extern float* X;
extern float* Y;
extern float* Z;

extern Tensors **tensors;
extern int tensorsCount;

extern int computeOn;
extern int printPDBOn;
extern int printDatOn;

extern int stretchOn;

extern FuncParam funcParam;

extern int* bonds;
extern int* bondsCount;
extern double* bondsR0;

extern int* native;
extern int* nativeCount;
extern double* nativeR0;
extern double* nativeEh;

extern int* pairs;
extern int* pairsCount;

extern int CONECT;

extern Vector normalVec;

//long int totalPairs;

extern int cutoffAveOn;
extern int cutoffSumOn;
extern int segmentAveOn;
extern int segmentSumOn;

#endif /* STRUCTURES_H_ */
