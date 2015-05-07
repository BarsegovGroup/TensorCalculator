/*
 * deform_gradient.c
 *
 *  Created on: Apr 1, 2015
 *      Author: olga
 */

#include "deform_gradient.h"
#include "../functions.h"

void createDeformGradientTensor(){
	printf("Creating deformation gradient tensor...\n");
	sprintf(deformGradientTensor.name, "deform_gradient");
	deformGradientTensor.compute = &computeDeformGradient;
	deformGradientTensor.read = &readDeformGradient;
	deformGradientTensor.write = &writeDeformGradient;
	deformGradientTensor.printPDB = &printPDBDeformGradient;
	deformGradientTensor.printDAT = &printDATDeformGradient;
	deformGradientTensor.destroy = &destroyDeformGradientTensor;
	tensors[tensorsCount] = &deformGradientTensor;
	tensorsCount ++;
	initDeformGradient();
}

void initDeformGradient(){
	printf("Initializing deformation gradient tensor...\n");

	int i;
	atomDeformGradient = (double**)calloc(pdbData.atomCount, sizeof(double*));
	for (i = 0; i < pdbData.atomCount; i++){
		atomDeformGradient[i] = (double*)calloc(9, sizeof(double));
	}

	atomRefCoord = (Vector*)calloc(pdbData.atomCount, sizeof(Vector));
	for (i = 0; i < pdbData.atomCount; i++){
		atomRefCoord[i].x = pdbData.atoms[i].x;
		atomRefCoord[i].y = pdbData.atoms[i].y;
		atomRefCoord[i].z = pdbData.atoms[i].z;
	}

	weightCov = getFloatParameter("deformGradientWeightCovalent", 1.0, 1);
	weightNat = getFloatParameter("deformGradientWeightNative", 1.0, 1);
	weightRep = getFloatParameter("deformGradientWeightRepulsive", 1.0, 1);

	int printDeformGradient = getYesNoParameter("printDeformGrad", 0, 1);
	if (printDeformGradient){
		if (computeOn){
			getMaskedParameter(outDeformGradientFilename, "outputDeformGradient", "", 0);
			if (strcmp(outDeformGradientFilename, "") != 0)
				outDeformGradient = fopen(outDeformGradientFilename, "w");
			fclose(outDeformGradient);
		}
		if (printPDBOn){
			getMaskedParameter(pdbDeformGradientFilename, "outputPDBDeformGradient", "", 0);
			if (strcmp(pdbDeformGradientFilename, "") != 0)
				pdbDeformGradient = fopen(pdbDeformGradientFilename, "w");
			fclose(pdbDeformGradient);

			if (!computeOn){
				getMaskedParameter(tnsrDeformGradientFilename, "inputDeformGradient", "", 0);
				if (strcmp(tnsrDeformGradientFilename, "") != 0)
					tnsrDeformGradient = fopen(tnsrDeformGradientFilename, "r");
			}
		}
	}
}

inline void computeDeformGradient(){
	printf("Computing gradient deformation tensor...\n");

	int atom, id1, id2, j;
	Vector atomDistance, atomRefDistance;

	double matrD[9];
	double matrA[9];
	double matrI[9];

	for (atom = 0; atom < pdbData.atomCount; atom++){

		for (j = 0; j < 9; j++){
			matrA[j] = 0.0;
			matrD[j] = 0.0;
			matrI[j] = 0.0;
			atomDeformGradient[atom][j] = 0.0;
		}

		for (j = 0; j < bondsCount[atom]; j++){
			id1 = atom, id2 = bonds[atom * MAX_BONDS + j];
			copyVector(&atomDistance, DCDtoVector(X[id2] - X[id1], Y[id2] - Y[id1], Z[id2] - Z[id1]));
			copyVector(&atomRefDistance, DCDtoVector(atomRefCoord[id2].x - atomRefCoord[id1].x,atomRefCoord[id2].y - atomRefCoord[id1].y, atomRefCoord[id2].z - atomRefCoord[id1].z));

			multVectors(atomDistance, atomDistance, matrI);
			multScalarMatrix(matrI, weightCov);
			//if (atom < 5){
			//	printf("%d-%d (%f): %f %f %f\n", id1, id2, weightCov, matrI[0], matrI[3], matrI[8]);
			//}
			addMatrix(matrD, matrI);

			multVectors(atomRefDistance, atomDistance, matrI);
			multScalarMatrix(matrI, weightCov);
			addMatrix(matrA, matrI);
		}

		for (j = 0; j < nativeCount[atom]; j++){
			id1 = atom; id2 = native[atom * MAX_NATIVE + j];
			copyVector(&atomDistance, DCDtoVector(X[id2] - X[id1], Y[id2] - Y[id1], Z[id2] - Z[id1]));
			copyVector(&atomRefDistance, DCDtoVector(atomRefCoord[id2].x - atomRefCoord[id1].x,atomRefCoord[id2].y - atomRefCoord[id1].y, atomRefCoord[id2].z - atomRefCoord[id1].z));

			multVectors(atomDistance, atomDistance, matrI);
			multScalarMatrix(matrI, weightNat);
			addMatrix(matrD, matrI);

			multVectors(atomRefDistance, atomDistance, matrI);
			multScalarMatrix(matrI, weightNat);
			addMatrix(matrA, matrI);
		}

		for (j = 0; j < pairsCount[atom]; j++){
			id1 = atom; id2 = pairs[atom * MAX_PAIRS + j];
			copyVector(&atomDistance, DCDtoVector(X[id1] - X[id2], Y[id1] - Y[id2], Z[id1] - Z[id2]));
			if (atomDistance.value  < pairs_cutoff && id1 != id2 ){
				copyVector(&atomDistance, DCDtoVector(X[id2] - X[id1], Y[id2] - Y[id1], Z[id2] - Z[id1]));
				copyVector(&atomRefDistance, DCDtoVector(atomRefCoord[id2].x - atomRefCoord[id1].x,atomRefCoord[id2].y - atomRefCoord[id1].y, atomRefCoord[id2].z - atomRefCoord[id1].z));

				multVectors(atomDistance, atomDistance, matrI);
				multScalarMatrix(matrI, weightRep);
				addMatrix(matrD, matrI);

				multVectors(atomRefDistance, atomDistance, matrI);
				multScalarMatrix(matrI, weightRep);
				addMatrix(matrA, matrI);
			}

		}

		inverseMatrix3d(matrD, matrI);
		multMatrix(matrA, matrI, atomDeformGradient[atom]);
	}
}

inline void readDeformGradient(){
	if (tnsrDeformGradient != NULL){
		printf("Reading input deformation gradient values...\n");
		readAtomTensorFrame(tnsrDeformGradient, atomDeformGradient);
	}
}

inline void writeDeformGradient(){
	//if (strcmp(outDeformGradientFilename, "") != 0){
		printf("Printing deformation gradient tensor...\n");
		int i, j;

		outDeformGradient = fopen(outDeformGradientFilename, "a");
		for (i = 0; i < pdbData.atomCount; i++){
			if ( strcmp(pdbData.atoms[i].segment, "") == 0)
				pdbData.atoms[i].segment[0] = pdbData.atoms[i].chain;
			//pdbData.atoms[i].segment[0] = pdbData.atoms[i].chain;
			//pdbData.atoms[i].segment[1] = '\0';
			printAtomTensorToFile(outDeformGradient, pdbData.atoms[i], atomDeformGradient[i]);
		}
		fprintf(outDeformGradient, "END\n");
		fclose(outDeformGradient);
	//}
}

inline void printPDBDeformGradient(){

	//if (strcmp(pdbDeformGradientFilename, "") != 0){
		double atomNormal;
		double atomShear;

		int i;

		for (i = 0; i < pdbData.atomCount; i++){
			calculateNormalShearComponent(atomDeformGradient[i], normalVec, &atomNormal, &atomShear);
			pdbData.atoms[i].x = X[i];
			pdbData.atoms[i].y = Y[i];
			pdbData.atoms[i].z = Z[i];
			pdbData.atoms[i].beta = atomShear;
			pdbData.atoms[i].occupancy = atomNormal;
			if ( strcmp(pdbData.atoms[i].segment, "") == 0)
				pdbData.atoms[i].segment[0] = pdbData.atoms[i].chain;
		}
		appendPDB(pdbDeformGradientFilename, &pdbData, CONECT);
	//}
}

inline void printDATDeformGradient(){

}

inline void destroyDeformGradientTensor(){
	free(atomDeformGradient);
	free(atomRefCoord);
	if (tnsrDeformGradient != NULL) fclose(tnsrDeformGradient);
}

