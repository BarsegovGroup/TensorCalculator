/*
 * tensorio.c
 *
 *  Created on: Feb 12, 2013
 *      Author: olga
 */

#include "../functions.h"
#include "string.h"
#include "stdlib.h"

void printAtomValueToFile(FILE* file, Atom atomData, double atomEnergy);
int readAtomValueFrame(FILE* file, double* atomEnergy);

void printAtomVectorToFile(FILE* file, Atom atomData, Vector atomVector);
int readAtomVectorFrame(FILE* file, Vector* atomVector);

void printAtomTensorToFile(FILE* file, Atom atomData, double* atomTensor);
int readAtomTensorFrame(FILE* file, double** atomTensorFrame);

void printAtomTensorToFile(FILE* file, Atom atomData, double* atomTensor){

	fprintf(file, "ATOM  %5d %3s %4d %c %4s ",
			atomData.id, atomData.resName, atomData.resid, atomData.chain, atomData.segment);
	int i;
	for (i = 0; i < 9; i++)
		fprintf(file, "%10.5f", atomTensor[i]);
	fprintf(file, "\n");

}

int readAtomTensorFrame(FILE* file, double** atomTensorFrame){
	int bufSize = 1024;
	char buffer[bufSize];

	int atomCount = 0;
	char* pch;
	int k;

	while ( (fgets(buffer, bufSize, file) != NULL) && (strncmp(buffer, "END", 3) != 0) ){
		pch = strtok(buffer, " ");
		for (k = 0; k < 5; k++){
			pch = strtok(NULL, " ");
		}
		atomTensorFrame[atomCount][0] = atof(pch);
		for (k = 1; k < 9; k++){
			pch = strtok(NULL, " ");
			if ( (strcmp(pch,"nan") == 0) || (strcmp(pch,"nan") == 0))
				atomTensorFrame[atomCount][k] = 0.0;
			else
				atomTensorFrame[atomCount][k] = atof(pch);
		}
		atomCount ++;
	}

	if (feof(file) == 0) return 0;
	else return -1;
}

void printAtomVectorToFile(FILE* file, Atom atomData, Vector atomVector){

	fprintf(file, "ATOM  %5d %3s %4d %c %4s %10.5f %10.5f %10.5f %10.5f\n",
			atomData.id, atomData.resName, atomData.resid, atomData.chain, atomData.segment, atomVector.x, atomVector.y, atomVector.z, atomVector.value);;
}

int readAtomVectorFrame(FILE* file, Vector* atomVector){
	int bufSize = 1024;
	char buffer[bufSize];

	int atomCount = 0;
	char* pch;
	int k;

	while ( (fgets(buffer, bufSize, file) != NULL) && (strncmp(buffer, "END", 3) != 0) ){
		pch = strtok(buffer, " ");
		for (k = 0; k < 6; k++)
			pch = strtok(NULL, " ");
		atomVector[atomCount].x = atof(pch);
		pch = strtok(NULL, " ");
		atomVector[atomCount].y = atof(pch);
		pch = strtok(NULL, " ");
		atomVector[atomCount].z = atof(pch);
		atomCount ++;
	}

	if (feof(file) == 0) return 0;
	else return -1;
}

void printAtomValueToFile(FILE* file, Atom atomData, double atomEnergy){
	fprintf(file, "ATOM  %5d %3s %4d %c %4s %10.5f\n",
			atomData.id, atomData.resName, atomData.resid, atomData.chain, atomData.segment, atomEnergy);
}
int readAtomValueFrame(FILE* file, double* atomEnergy){
	int bufSize = 1024;
	char buffer[bufSize];

	int atomCount = 0;

	char* pch;
	int k;

	while ( (fgets(buffer, bufSize, file) != NULL) && (strncmp(buffer, "END", 3) != 0) ){
		pch = strtok(buffer, " ");
		for (k = 0; k < 6; k++){
			pch = strtok(NULL, " ");
		}
		atomEnergy[atomCount] = atof(pch);
		//printf("%s %f\n", pch, atomEnergy[atomCount]);
		atomCount ++;
	}
	//printf("\n");

	//int i;
	//for (i = 0; i < 5; i++) printf("%d: %f\n", i, atomEnergy[i]);

	if (feof(file) == 0){
		return 0;
	}else{
		return -1;
	}
}
