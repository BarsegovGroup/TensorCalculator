/*
 * cauchy_tensor.c
 *
 *  Created on: Mar 30, 2015
 *      Author: olga
 */

#include "cauchy_stress.h"

#include "../IO/config_reader.h"
#include "../IO/tensorio.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "math.h"

Tensors cauchyStressTensor;

double** atomCauchyStress;
double* atomRadius;

double** atomCauchyStressCov;
double* atomRadiusCov;

double** atomCauchyStressNat;
double* atomRadiusNat;

double** atomCauchyStressRep;
double* atomRadiusRep;

char outCauchyStressFilename[256];
FILE* outCauchyStress;

char pdbCauchyStressFilename[256];
FILE* pdbCauchyStress;

char tnsrCauchyStressFilename[256];
FILE* tnsrCauchyStress;

char datCauchyStressFilename[256];
FILE* datCauchyStress;

double cauchyNormalScale;
double cauchyShearScale;

int averageCauchyStress;

Vector normalVec;

void createCauchyStressTensor(){
	printf("Creating Cauchy stress tensor...\n");
	sprintf(cauchyStressTensor.name, "cauchy");
	cauchyStressTensor.compute = &computeCauchyStress;
	cauchyStressTensor.read = &readCauchyStress;
	cauchyStressTensor.write = &writeCauchyStress;
	cauchyStressTensor.printPDB = &printPDBCauchyStress;
	cauchyStressTensor.printDAT = &printDATCauchyStress;
	cauchyStressTensor.destroy = &destroyCauchyStressTensor;
	cauchyStressTensor.runningAverage = &averageCauchyStressTensor;
	tensors[tensorsCount] = &cauchyStressTensor;
	tensorsCount ++;
	initCauchyStress();
}

void initCauchyStress(){
	printf("Initializing Cauchy stress tensor...\n");

	int i;
	atomCauchyStress = (double**)calloc(pdbData.atomCount, sizeof(double*));
	atomCauchyStressCov = (double**)calloc(pdbData.atomCount, sizeof(double*));
	atomCauchyStressNat = (double**)calloc(pdbData.atomCount, sizeof(double*));
	atomCauchyStressRep = (double**)calloc(pdbData.atomCount, sizeof(double*));
	for (i = 0; i < pdbData.atomCount; i++){
		atomCauchyStress[i] = (double*)calloc(9, sizeof(double));
		atomCauchyStressCov[i] = (double*)calloc(9, sizeof(double));
		atomCauchyStressNat[i] = (double*)calloc(9, sizeof(double));
		atomCauchyStressRep[i] = (double*)calloc(9, sizeof(double));
	}

	atomRadius = (double*)calloc(pdbData.atomCount, sizeof(double));
	atomRadiusCov = (double*)calloc(pdbData.atomCount, sizeof(double));
	atomRadiusNat = (double*)calloc(pdbData.atomCount, sizeof(double));
	atomRadiusRep = (double*)calloc(pdbData.atomCount, sizeof(double));

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
		cauchyNormalScale = getDoubleParameter("cauchyNormalScaleFactor", 1.0, 1);
		cauchyShearScale = getDoubleParameter("cauchyShearScaleFactor", 1.0, 1);
		averageCauchyStress = getIntegerParameter("averageCauchyStress", 1, 1);
	}

	if (printDatOn){
		if (!computeOn){
			getMaskedParameter(tnsrCauchyStressFilename, "inputCauchyStress", "", 0);
			tnsrCauchyStress = fopen(tnsrCauchyStressFilename, "r");
		}
		getMaskedParameter(datCauchyStressFilename, "datCauchyStress", "", 0);
		datCauchyStress = fopen(datCauchyStressFilename, "w");
		fclose(datCauchyStress);
		cauchyNormalScale = getDoubleParameter("cauchyNormalScaleFactor", 1.0, 1);
		cauchyShearScale = getDoubleParameter("cauchyShearScaleFactor", 1.0, 1);
	}

}

inline void computeCauchyStress(){
	printf("Computing Cauchy stress tensor...\n");

	int atom, id1, id2, j;
	Vector atomDistance, atomForce;
	double r, r2;
	double atomVolume;

	for (atom = 0; atom < pdbData.atomCount; atom++){
		for (j = 0; j < 9; j++){
			atomCauchyStress[atom][j] = 0.0;
			atomCauchyStressCov[atom][j] = 0.0;
			atomCauchyStressNat[atom][j] = 0.0;
			atomCauchyStressRep[atom][j] = 0.0;			
		}
	}

	for (atom = 0; atom < pdbData.atomCount; atom++){
		atomRadius[atom] = 0.0;
		atomRadiusCov[atom] = 0.0;
		atomRadiusNat[atom] = 0.0;
		atomRadiusRep[atom] = 0.0;			
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
		if (r2 == 0) atomRadiusCov[atom] = 0.0;
		else atomRadiusCov[atom] = r/r2;

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
		if (r2 == 0) atomRadiusNat[atom] = 0.0;
		else atomRadiusNat[atom] = r/r2;

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
		if (r2 == 0) atomRadiusRep[atom] = 0.0;
		else atomRadiusRep[atom] = r/r2;
		
		atomVolume = 4.0/3.0 * Pi * (r/r2) * (r/r2) * (r/r2);

		for (j = 0; j < 9; j++){
			atomCauchyStress[atom][j] *= 1000.0 * 0.5 / atomVolume;
		}
		atomRadius[atom] = r/r2;
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

		printAtomTensorToFile(outCauchyStress, pdbData.atoms[i], atomCauchyStress[i]);
	}
	fprintf(outCauchyStress, "END\n");
	fclose(outCauchyStress);
}

inline void printPDBCauchyStress(){
	printf("Printing PDB for Cauchy stress tensor...\n");
	double atomNormal;
	double atomShear;
	Vector invar;
	Vector prncp;

	int i;

	for (i = 0; i < pdbData.atomCount; i++){
		if (CONECT == 1){
			pdbData.atoms[i].x = X[i];
			pdbData.atoms[i].y = Y[i];
			pdbData.atoms[i].z = Z[i];
		}
		if (normalVec.value > -1){
			calculateNormalShearComponent(atomCauchyStress[i], normalVec, &atomNormal, &atomShear);
			pdbData.atoms[i].beta = atomShear * cauchyShearScale;
			pdbData.atoms[i].occupancy = atomNormal *  cauchyNormalScale;
			pdbData.atoms[i].charge = 0.0;
		}else{
			getInvariants(atomCauchyStress[i], &invar);
			getPrincipalStresses(atomCauchyStress[i], &prncp);
			pdbData.atoms[i].occupancy = invar.x*cauchyNormalScale;
			pdbData.atoms[i].beta = prncp.x*cauchyNormalScale;
			pdbData.atoms[i].charge = getVonMisses(atomCauchyStress[i])*cauchyShearScale;
		}
		if ( strcmp(pdbData.atoms[i].segment, "") == 0)
			pdbData.atoms[i].segment[0] = pdbData.atoms[i].chain;
	}
	if (cutoffAveOn || cutoffSumOn) makeCutOffAveraged();
	if (segmentAveOn || segmentSumOn) makeSegmentAveraged();
	appendPDB(pdbCauchyStressFilename, &pdbData, CONECT);
}

inline void printDATCauchyStress(){
	datCauchyStress = fopen(datCauchyStressFilename, "a");
	int i;

	Vector invar;
	Vector prncpl;

	double atomNormal;
	double atomShear;
	float totalNormal = 0.0;
	float totalShear = 0.0;
	float totalInvar1 = 0.0;
	float totalInvar2 = 0.0;
	float totalInvar3 = 0.0;
	float totalPrincipal1 = 0.0;
	float totalPrincipal2 = 0.0;
	float totalPrincipal3 = 0.0;
	float totalVM = 0.0;
	float vm = 0.0;
	float totalRadius = 0.0;
	float totalTresca = 0.0;


	for (i = 0; i < pdbData.atomCount; i++){

		if (normalVec.value > -1){
			calculateNormalShearComponent(atomCauchyStress[i], normalVec, &atomNormal, &atomShear);
			totalNormal += atomNormal;
			totalShear += atomShear;
		}
		getInvariants(atomCauchyStress[i], &invar);
		totalInvar1 += invar.x;
		totalInvar2 += invar.y;
		totalInvar3 += invar.z;
		
		vm = getVonMisses(atomCauchyStress[i]);
		totalVM += vm;
		
		getPrincipalStresses(atomCauchyStress[i], &prncpl);
		
		totalPrincipal1 += prncpl.x;
		totalPrincipal2 += prncpl.y;
		totalPrincipal3 += prncpl.z;

		totalRadius += atomRadius[i];

		totalTresca += 0.5 * (prncpl.x - prncpl.z);

	}

	if (normalVec.value > -1){
		fprintf(datCauchyStress, "%8.5f\t%8.5f\t%8.5f\t%8.5f\t%8.5f\t%8.5f\t%8.5f\t%8.5f\t%8.5f\t%8.5f\t%8.5f\n", totalInvar1/pdbData.atomCount, totalInvar2/pdbData.atomCount, totalInvar3/pdbData.atomCount, totalPrincipal1/pdbData.atomCount, totalPrincipal2/pdbData.atomCount, totalPrincipal3/pdbData.atomCount, totalVM/pdbData.atomCount, totalTresca/pdbData.atomCount, totalRadius/pdbData.atomCount, totalNormal/pdbData.atomCount, totalShear/pdbData.atomCount);
	}else{
		fprintf(datCauchyStress, "%8.5f\t%8.5f\t%8.5f\t%8.5f\t%8.5f\t%8.5f\t%8.5f\t%8.5f\t%8.5f\n", totalInvar1/pdbData.atomCount, totalInvar2/pdbData.atomCount, totalInvar3/pdbData.atomCount, totalPrincipal1/pdbData.atomCount, totalPrincipal2/pdbData.atomCount, totalPrincipal3/pdbData.atomCount, totalVM/pdbData.atomCount, totalTresca/pdbData.atomCount, totalRadius/pdbData.atomCount);
	}

	fclose(datCauchyStress);
}

inline void destroyCauchyStressTensor(){
	free(atomCauchyStress);
	if (tnsrCauchyStress != NULL) fclose(tnsrCauchyStress);
}

inline void averageCauchyStressTensor(){
	if (averageCauchyStress > 1)
		makeRunningAveraged(pdbCauchyStressFilename, averageCauchyStress);
}

