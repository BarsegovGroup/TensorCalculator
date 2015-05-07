/*
 * structures.h
 *
 *  Created on: Mar 30, 2015
 *      Author: olga
 */

#ifndef STRUCTURES_H_
#define STRUCTURES_H_

#include "IO/pdbio.h"
#include "parameters.h"

typedef struct{
	char name[100];
	void (*compute)();
	void (*write)();
	void (*read)();
	void (*printPDB)();
	void (*printDAT)();
	void (*destroy)();
} Tensors;

Tensors **tensors;
int tensorsCount;

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

FuncParam funcParam;

int* bonds;
int* bondsCount;
double* bondsR0;

int* native;
int* nativeCount;
double* nativeR0;
double* nativeEh;

int* pairs;
int* pairsCount;

PDB pdbData;

FILE* dcdFile;
char dcdFilenameOrig[BUF_SIZE];
int N;
int NFILE;
int NSAVC;
int NPRIV;
float DELTA;
int frameCount;
int initFrame;
int finFrame;
float* X;
float* Y;
float* Z;

int CONECT;

Vector normalVec;

long int totalPairs;

int computeOn;
int printPDBOn;
int printDatOn;

int cauchyOn;
int sphereOn;
int energyOn;
int strainOn;
int stretchOn;
int cutoffAveOn;
int cutoffSumOn;
int segmentAveOn;
int segmentSumOn;

#endif /* STRUCTURES_H_ */
