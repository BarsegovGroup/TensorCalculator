/*
 * strain_tensor.c
 *
 *  Created on: Apr 2, 2015
 *      Author: olga
 */

#include "strain_tensor.h"
#include "../functions.h"

void computeStretch();
void writeStretch();
void readStretch();
void printPDBStretch();

void createStrainTensor(){
	printf("Creating strain stress tensor...\n");
	sprintf(strainTensor.name, "strain");
	strainTensor.compute = &computeStrain;
	strainTensor.read = &readStrain;
	strainTensor.write = &writeStrain;
	strainTensor.printPDB = &printPDBStrain;
	strainTensor.printDAT = &printDATStrain;
	strainTensor.destroy = &destroyStrainTensor;
	tensors[tensorsCount] = &strainTensor;
	tensorsCount ++;
	initStrain();
}

void initStrain(){
	printf("Initializing strain tensor...\n");

	//createDeformGradientTensor();

	initDeformGradient();

	if (stretchOn) initStretch();

	atomGreenStrain = (double**)calloc(pdbData.atomCount, sizeof(double*));
	int i;
	for (i = 0; i < pdbData.atomCount; i++){
		atomGreenStrain[i] = (double*)calloc(9, sizeof(double));
	}

	if (computeOn){
		getMaskedParameter(outGreenStrainFilename, "outputGreenStrain", "", 0);
		outGreenStrain = fopen(outGreenStrainFilename, "w");
		fclose(outGreenStrain);
	}

	if (printPDBOn){
		if (!computeOn){
			getMaskedParameter(tnsrGreenStrainFilename, "inputGreenStrain", "", 0);
			tnsrGreenStrain = fopen(tnsrGreenStrainFilename, "r");
		}
		getMaskedParameter(pdbGreenStrainFilename, "outputPDBGreenStrain", "", 0);
		pdbGreenStrain = fopen(pdbGreenStrainFilename, "w");
		fclose(pdbGreenStrain);
		strainScale = getDoubleParameter("strainScaleFactor", 1.0, 1);
	}

}

inline void computeStrain(){
	printf("Computing Green strain tensor...\n");

	computeDeformGradient();

	int j, atom;

	double matrI[9];

	for (atom = 0; atom < pdbData.atomCount; atom++){

		for (j = 0; j < 9; j++){
			atomGreenStrain[atom][j] = 0.0;
		}

		transposeMatrix(atomDeformGradient[atom], matrI);
		multMatrix(matrI, atomDeformGradient[atom], atomGreenStrain[atom]);

		for(j = 0; j < 9; j++)
			matrI[j] = 0.0;
		matrI[0] = -1.0;
		matrI[4] = -1.0;
		matrI[8] = -1.0;

		addMatrix(atomGreenStrain[atom], matrI);
		multScalarMatrix(atomGreenStrain[atom], 0.5);

	}
	if(stretchOn)
		computeStretch();
}

inline void readStrain(){

	if (tnsrDeformGradient != NULL) readDeformGradient();

	printf("Reading input Green's strain tensor values...\n");
	readAtomTensorFrame(tnsrGreenStrain, atomGreenStrain);

	if (tnsrStretch != NULL) readStretch();
}

inline void writeStrain(){

	if (strcmp(outDeformGradientFilename, "") != 0) writeDeformGradient();

	printf("Printing Green's strain tensor...\n");
	int i;

	outGreenStrain = fopen(outGreenStrainFilename, "a");
	for (i = 0; i < pdbData.atomCount; i++){
		//pdbData.atoms[i].segment[0] = pdbData.atoms[i].chain;
		//pdbData.atoms[i].segment[1] = '\0';
		if ( strcmp(pdbData.atoms[i].segment, "") == 0)
			pdbData.atoms[i].segment[0] = pdbData.atoms[i].chain;
		printAtomTensorToFile(outGreenStrain, pdbData.atoms[i], atomGreenStrain[i]);
	}
	fprintf(outGreenStrain, "END\n");
	fclose(outGreenStrain);

	if (strcmp(outStretchFilename, "") != 0) writeStretch();
}

inline void printPDBStrain(){

	printf("Printing PDB for strain stress tensor...\n");

	if (strcmp(pdbDeformGradientFilename, "") != 0) printPDBDeformGradient();

	double atomNormal;
	double atomShear;

	int i;

	for (i = 0; i < pdbData.atomCount; i++){
		//if (Z[i] >= 0) normalVec.z = 1.0;
		//else normalVec.z = -1.0;
		calculateNormalShearComponent(atomGreenStrain[i], normalVec, &atomNormal, &atomShear);
		pdbData.atoms[i].x = X[i];
		pdbData.atoms[i].y = Y[i];
		pdbData.atoms[i].z = Z[i];
		pdbData.atoms[i].beta = atomShear * strainScale;
		pdbData.atoms[i].occupancy = atomNormal * strainScale;
		if ( strcmp(pdbData.atoms[i].segment, "") == 0)
			pdbData.atoms[i].segment[0] = pdbData.atoms[i].chain;
	}
	if (cutoffAveOn || cutoffSumOn) makeCutOffAveraged();
	if (segmentAveOn || segmentSumOn) makeSegmentAveraged();
	appendPDB(pdbGreenStrainFilename, &pdbData, CONECT);

	if (strcmp(pdbStretchFilename, "") != 0) printPDBStretch();
}

inline void printDATStrain(){

}

inline void destroyStrainTensor(){
	free(atomGreenStrain);
	if (tnsrGreenStrain != NULL) fclose(tnsrGreenStrain);
	destroyDeformGradientTensor();
	destroyStretchTensor();
}
/*
void computeStretch(){
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

void writeStretch(){
	printf("Printing stretch tensor...\n");
	int i;

	outStretch = fopen(outStretchFilename, "a");
	for (i = 0; i < pdbData.atomCount; i++){
		pdbData.atoms[i].segment[0] = pdbData.atoms[i].chain;
		pdbData.atoms[i].segment[1] = '\0';
		printAtomTensorToFile(outStretch, pdbData.atoms[i], atomStretch[i]);
	}
	fprintf(outStretch, "END\n");
	fclose(outStretch);
}

void readStretch(){
	printf("Reading input stretch tensor values...\n");
	readAtomTensorFrame(tnsrStretch, atomStretch);
}

void printPDBStretch(){
	double atomNormal;
	double atomShear;

	int i;

	for (i = 0; i < pdbData.atomCount; i++){
		calculateNormalShearComponent(atomStretch[i], normalVec, &atomNormal, &atomShear);
		pdbData.atoms[i].x = X[i];
		pdbData.atoms[i].y = Y[i];
		pdbData.atoms[i].z = Z[i];
		pdbData.atoms[i].beta = atomShear;
		pdbData.atoms[i].occupancy = atomNormal;
		pdbData.atoms[i].charge = 0.0;
	}

	appendPDB(pdbStretchFilename, &pdbData, CONECT);
}
*/

