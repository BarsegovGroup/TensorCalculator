/*
 * cauchy_tensor.c
 *
 *  Created on: Mar 30, 2015
 *      Author: olga
 */

#include "cauchy_stress.h"
#include "../functions.h"

void createCauchyStressTensor(){
	printf("Creating Cauchy stress tensor...\n");
	sprintf(cauchyStressTensor.name, "cauchy");
	cauchyStressTensor.compute = &computeCauchyStress;
	cauchyStressTensor.read = &readCauchyStress;
	cauchyStressTensor.write = &writeCauchyStress;
	cauchyStressTensor.printPDB = &printPDBCauchyStress;
	cauchyStressTensor.printDAT = &printDATCauchyStress;
	cauchyStressTensor.destroy = &destroyCauchyStressTensor;
	tensors[tensorsCount] = &cauchyStressTensor;
	tensorsCount ++;
	initCauchyStress();
}

void initCauchyStress(){
	printf("Initializing Cauchy stress tensor...\n");

	int i;
	atomCauchyStress = (double**)calloc(pdbData.atomCount, sizeof(double*));
	for (i = 0; i < pdbData.atomCount; i++){
		atomCauchyStress[i] = (double*)calloc(9, sizeof(double));
	}

	if (computeOn){
		getMaskedParameter(outCauchyStressFilename, "outputCauchyStress", "", 0);
		outCauchyStress = fopen(outCauchyStressFilename, "w");
		fclose(outCauchyStress);
	}
	if (printPDBOn){
		if (!computeOn){
			getMaskedParameter(tnsrCauchyStressFilename, "inputCauchyStress", "", 0);
			tnsrCauchyStress = fopen(tnsrCauchyStressFilename, "r");
		}
		getMaskedParameter(pdbCauchyStressFilename, "outputPDBCauchyStress", "", 0);
		pdbCauchyStress = fopen(pdbCauchyStressFilename, "w");
		fclose(pdbCauchyStress);
		cauchyScale = getDoubleParameter("cauchyScaleFactor", 1.0, 1);
	}

}

inline void computeCauchyStress(){
	printf("Computing Cauchy stress tensor...\n");

	int atom, id1, id2, j;
	Vector atomDistance, atomForce;
	double r, r2;
	double atomVolume;

	for (atom = 0; atom < pdbData.atomCount; atom++){
		for (j = 0; j < 9; j++)
			atomCauchyStress[atom][j] = 0.0;
	}

	for (atom = 0; atom < pdbData.atomCount; atom++){
		r = 0.0;
		r2 = 0.0;
		copyVector(&atomDistance, nullVector());
		copyVector(&atomForce, nullVector());
		for (j = 0; j < bondsCount[atom]; j++){
			id1 = atom, id2 = bonds[atom * MAX_BONDS + j];
			funcParam.refDist = bondsR0[atom * MAX_BONDS + j];
			copyVector(&atomDistance, DCDtoVector(X[id1] - X[id2], Y[id1] - Y[id2], Z[id1] - Z[id2]));
			copyVector(&atomForce, getFENEDerivative_1(atomDistance));
			atomCauchyStress[atom][0] += atomForce.x * atomDistance.x;
			atomCauchyStress[atom][1] += atomForce.x * atomDistance.y;
			atomCauchyStress[atom][2] += atomForce.x * atomDistance.z;
			atomCauchyStress[atom][3] += atomForce.y * atomDistance.x;
			atomCauchyStress[atom][4] += atomForce.y * atomDistance.y;
			atomCauchyStress[atom][5] += atomForce.y * atomDistance.z;
			atomCauchyStress[atom][6] += atomForce.z * atomDistance.x;
			atomCauchyStress[atom][7] += atomForce.z * atomDistance.y;
			atomCauchyStress[atom][8] += atomForce.z * atomDistance.z;
			r += 1.0/atomDistance.value;
			r2 += 1.0/(atomDistance.value * atomDistance.value);
		}

		for (j = 0; j < nativeCount[atom]; j++){
			id1 = atom; id2 = native[atom * MAX_NATIVE + j];
			funcParam.refDist = nativeR0[atom * MAX_NATIVE + j];
			funcParam.Eh = nativeEh[atom * MAX_NATIVE + j];
			copyVector(&atomDistance, DCDtoVector(X[id1] - X[id2], Y[id1] - Y[id2], Z[id1] - Z[id2]));
			copyVector(&atomForce, getNativeDerivative_1(atomDistance));
			atomCauchyStress[atom][0] += atomForce.x * atomDistance.x;
			atomCauchyStress[atom][1] += atomForce.x * atomDistance.y;
			atomCauchyStress[atom][2] += atomForce.x * atomDistance.z;
			atomCauchyStress[atom][3] += atomForce.y * atomDistance.x;
			atomCauchyStress[atom][4] += atomForce.y * atomDistance.y;
			atomCauchyStress[atom][5] += atomForce.y * atomDistance.z;
			atomCauchyStress[atom][6] += atomForce.z * atomDistance.x;
			atomCauchyStress[atom][7] += atomForce.z * atomDistance.y;
			atomCauchyStress[atom][8] += atomForce.z * atomDistance.z;
			r += 1.0/atomDistance.value;
			r2 += 1.0/(atomDistance.value * atomDistance.value);
		}

		for (j = 0; j < pairsCount[atom]; j++){
			id1 = atom; id2 = pairs[atom * MAX_PAIRS + j];
			copyVector(&atomDistance, DCDtoVector(X[id1] - X[id2], Y[id1] - Y[id2], Z[id1] - Z[id2]));
			if (atomDistance.value  < pairs_cutoff && id1 != id2 ){
				copyVector(&atomDistance, DCDtoVector(X[id1] - X[id2], Y[id1] - Y[id2], Z[id1] - Z[id2]));
				copyVector(&atomForce, getRepulsiveDerivative_1(atomDistance));
				atomCauchyStress[atom][0] += atomForce.x * atomDistance.x;
				atomCauchyStress[atom][1] += atomForce.x * atomDistance.y;
				atomCauchyStress[atom][2] += atomForce.x * atomDistance.z;
				atomCauchyStress[atom][3] += atomForce.y * atomDistance.x;
				atomCauchyStress[atom][4] += atomForce.y * atomDistance.y;
				atomCauchyStress[atom][5] += atomForce.y * atomDistance.z;
				atomCauchyStress[atom][6] += atomForce.z * atomDistance.x;
				atomCauchyStress[atom][7] += atomForce.z * atomDistance.y;
				atomCauchyStress[atom][8] += atomForce.z * atomDistance.z;
				r += 1.0/atomDistance.value;
				r2 += 1.0/(atomDistance.value * atomDistance.value);
			}

		}

		atomVolume = 4.0/3.0 * Pi * (0.5*r/r2) * (0.5*r/r2) * (0.5*r/r2);
		for (j = 0; j < 9; j++)
			atomCauchyStress[atom][j] *= 0.5 / atomVolume;

	}
}

inline void readCauchyStress(){
	printf("Reading input Cauchy stress values...\n");
	readAtomTensorFrame(tnsrCauchyStress, atomCauchyStress);
}

inline void writeCauchyStress(){
	printf("Printing Cauchy stress tensor...\n");
	outCauchyStress = fopen(outCauchyStressFilename, "a");
	int i;
	for (i = 0; i < pdbData.atomCount; i++){
		if ( strcmp(pdbData.atoms[i].segment, "") == 0)
			pdbData.atoms[i].segment[0] = pdbData.atoms[i].chain;
		//pdbData.atoms[i].segment[0] = '\0';
		printAtomTensorToFile(outCauchyStress, pdbData.atoms[i], atomCauchyStress[i]);
	}
	fprintf(outCauchyStress, "END\n");
	fclose(outCauchyStress);
}

inline void printPDBCauchyStress(){
	double atomNormal;
	double atomShear;

	int i;

	for (i = 0; i < pdbData.atomCount; i++){
		calculateNormalShearComponent(atomCauchyStress[i], normalVec, &atomNormal, &atomShear);
		pdbData.atoms[i].x = X[i];
		pdbData.atoms[i].y = Y[i];
		pdbData.atoms[i].z = Z[i];
		pdbData.atoms[i].beta = atomShear * cauchyScale;
		pdbData.atoms[i].occupancy = atomNormal *  cauchyScale;
		if ( strcmp(pdbData.atoms[i].segment, "") == 0)
			pdbData.atoms[i].segment[0] = pdbData.atoms[i].chain;
	}
	if (cutoffAveOn || cutoffSumOn) makeCutOffAveraged();
	if (segmentAveOn || segmentSumOn) makeSegmentAveraged();
	appendPDB(pdbCauchyStressFilename, &pdbData, CONECT);
}

inline void printDATCauchyStress(){

}

inline void destroyCauchyStressTensor(){
	free(atomCauchyStress);
	if (tnsrCauchyStress != NULL) fclose(tnsrCauchyStress);
	//if (outCauchyStress != NULL) fclose(outCauchyStress);
	//if (pdbCauchyStress != NULL) fclose(pdbCauchyStress);
}

