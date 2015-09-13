/*
 * stretch_tensor.c
 *
 *  Created on: Apr 6, 2015
 *      Author: olga
 */

#include "stretch_tensor.h"
#include "../functions.h"
#include "deform_gradient.h"

void createStretchTensor(){
	printf("Creating stretch tensor...\n");
	sprintf(stretchTensor.name, "stretch");
	stretchTensor.compute = &computeStretch;
	stretchTensor.read = &readStretch;
	stretchTensor.write = &writeStretch;
	stretchTensor.printPDB = &printPDBStretch;
	stretchTensor.printDAT = &printDATStretch;
	stretchTensor.destroy = &destroyStretchTensor;
	tensors[tensorsCount] = &stretchTensor;
	tensorsCount ++;
	initStretch();
}

void initStretch(){
	printf("Initializing stretch tensor...\n");

	int i;
	atomStretch = (double**)calloc(pdbData.atomCount, sizeof(double*));
	for (i = 0; i < pdbData.atomCount; i++){
		atomStretch[i] = (double*)calloc(9, sizeof(double));
	}

	//if (stretchOn){
		if (computeOn){
			getMaskedParameter(outStretchFilename, "outputStretch", "", 1);
			outStretch = fopen(outStretchFilename, "w");
			fclose(outStretch);
		}

		if (printPDBOn){
			if (!computeOn){
				getMaskedParameter(tnsrStretchFilename, "inputStretch", "", 1);
				tnsrStretch = fopen(tnsrStretchFilename, "r");
			}
			getMaskedParameter(pdbStretchFilename, "outputPDBStretch", "", 1);
			pdbStretch = fopen(pdbStretchFilename, "w");
			fclose(pdbStretch);
			stretchScale = getDoubleParameter("stretchScaleFactor", 1.0, 1);
		}
	//}
}

inline void computeStretch(){
	printf("Computing stretch tensor...\n");

	int atom, j;
	double eigenVec[9];
	double eigenVal[9];
	double eigenVecT[9];
	double strain[9];

	for (atom = 0; atom < pdbData.atomCount; atom++){
		for (j = 0; j < 9; j++)
			atomStretch[atom][j] = 0.0;
		transposeMatrix(atomDeformGradient[atom], strain);
		multMatrix(strain, atomDeformGradient[atom], strain);
		eigenDecompos3d(strain, eigenVal, eigenVec);
		eigenVal[0] = sqrt(eigenVal[0]);
		eigenVal[4] = sqrt(eigenVal[4]);
		eigenVal[8] = sqrt(eigenVal[8]);
		multMatrix(eigenVec, eigenVal, atomStretch[atom]);
		transposeMatrix(eigenVec, eigenVecT);
		multMatrix(atomStretch[atom], eigenVecT, atomStretch[atom]);
	}
}

inline void readStretch(){
	if (tnsrStretch != NULL){
		printf("Reading input stretch tensor values...\n");
		readAtomTensorFrame(tnsrStretch, atomStretch);
	}
}

inline void writeStretch(){
	//if (strcpy(outStretchFilename, "") != 0){
		printf("Printing stretch tensor...\n");
		int i,j;

		outStretch = fopen(outStretchFilename, "a");
		//fprintf(outStretch, "lalala\n");
		for (i = 0; i < pdbData.atomCount; i++){
			//pdbData.atoms[i].segment[0] = pdbData.atoms[i].chain;
			//pdbData.atoms[i].segment[1] = '\0';
			if ( strcmp(pdbData.atoms[i].segment, "") == 0)
				pdbData.atoms[i].segment[0] = pdbData.atoms[i].chain;
			printAtomTensorToFile(outStretch, pdbData.atoms[i], atomStretch[i]);

			//printf("ATOM  %5d %3s %4d %c %4s ",
			//		pdbData.atoms[i].id, pdbData.atoms[i].resName, pdbData.atoms[i].resid, pdbData.atoms[i].chain, pdbData.atoms[i].segment);
			//for (j = 0; j < 9; j++)
			//	printf("%10.5f", atomStretch[j]);
			//printf("\n");

		}
		//if (outStretch == NULL) printf("File does not exist\n");
		fprintf(outStretch, "END\n");
		fclose(outStretch);
	//}
}

inline void printPDBStretch(){
	printf("Printing PDB for stretch stress tensor...\n");
	//if (strcpy(pdbStretchFilename, "") != 0){
		double atomNormal;
		double atomShear;

		int i;

		for (i = 0; i < pdbData.atomCount; i++){
			calculateNormalShearComponent(atomStretch[i], normalVec, &atomNormal, &atomShear);
			pdbData.atoms[i].x = X[i];
			pdbData.atoms[i].y = Y[i];
			pdbData.atoms[i].z = Z[i];
			if (isnan(atomShear)) atomShear = 0.0;
			if (isnan(atomNormal)) atomNormal = 0.0;
			pdbData.atoms[i].beta = atomShear * stretchScale;
			pdbData.atoms[i].occupancy = atomNormal * stretchScale;
			pdbData.atoms[i].charge = 0.0;
			if ( strcmp(pdbData.atoms[i].segment, "") == 0)
				pdbData.atoms[i].segment[0] = pdbData.atoms[i].chain;
		}
	//}
	if (cutoffAveOn || cutoffSumOn) makeCutOffAveraged();
	if (segmentAveOn || segmentSumOn) makeSegmentAveraged();
	appendPDB(pdbStretchFilename, &pdbData, CONECT);
}

inline void printDATStretch(){

}

inline void destroyStretchTensor(){
	if (tnsrStretch != NULL) fclose(tnsrStretch);
	free(atomStretch);
}


