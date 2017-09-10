/*
 * sphere_tensor.c
 *
 *  Created on: Mar 31, 2015
 *      Author: olga
 */

#include "sphere_stress.h"
#include "../IO/config_reader.h"
#include "../IO/tensorio.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

Tensors sphereStressTensor;

Vector* atomSphereStress; //x - lateral, y - shear, z - radial

Vector ePhi, eTheta, eRad;

char outSphereStressFilename[256];
FILE* outSphereStress;

char pdbSphereStressFilename[256];
FILE* pdbSphereStress;

char tnsrSphereStressFilename[256];
FILE* tnsrSphereStress;

double sphereLateralScale;
double sphere45ShearScale;
double sphereRadialScale;

void createSphereStressTensor(){
	printf("Creating sphere stress tensor...\n");
	sprintf(sphereStressTensor.name, "sphere");
	sphereStressTensor.compute = &computeSphereStress;
	sphereStressTensor.read = &readSphereStress;
	sphereStressTensor.write = &writeSphereStress;
	sphereStressTensor.printPDB = &printPDBSphereStress;
	sphereStressTensor.printDAT = &printDATSphereStress;
	sphereStressTensor.destroy = &destroySphereStressTensor;
	tensors[tensorsCount] = &sphereStressTensor;
	tensorsCount ++;
	initSphereStress();
}

void initSphereStress(){
	printf("Initializing sphere stress tensor...\n");

	atomSphereStress = (Vector*)calloc(pdbData.atomCount, sizeof(Vector));

	if (computeOn){
		getMaskedParameter(outSphereStressFilename, "outputSphereStress", "", 0);
		outSphereStress = fopen(outSphereStressFilename, "w");
		fclose(outSphereStress);
	}

	if (printPDBOn) {
		if (!computeOn){
			getMaskedParameter(tnsrSphereStressFilename, "inputSphereStress", "", 0);
			tnsrSphereStress = fopen(tnsrSphereStressFilename, "r");
		}
		getMaskedParameter(pdbSphereStressFilename, "outputPDBSphereStress", "", 0);
		pdbSphereStress = fopen(pdbSphereStressFilename, "w");
		fclose(pdbSphereStress);
		sphereLateralScale = getDoubleParameter("sphereScaleFactor", 1.0, 1);
		sphere45ShearScale = getDoubleParameter("sphereScaleFactor", 1.0, 1);
		sphereRadialScale = getDoubleParameter("sphereScaleFactor", 1.0, 1);
	}

}

inline void computeSphereStress(){
	printf("Computing sphere stress tensor...\n");

	int atom, id1, id2, j;
	Vector atomDistance, atomForce;
	double r, r2;
	double atomArea, atomVolume;
	double sp_1, sp_2, sp_3;

	for (atom = 0; atom < pdbData.atomCount; atom++){
		atomSphereStress[atom].x = 0.0;
		atomSphereStress[atom].y = 0.0;
		atomSphereStress[atom].z = 0.0;
		atomSphereStress[atom].value = 0.0;
	}

	for (atom = 0; atom < pdbData.atomCount; atom++){
		r = 0.0;
		r2 = 0.0;
		copyVector(&atomDistance, nullVector());
		copyVector(&atomForce, nullVector());

		copyVector(&ePhi, DCDtoEphi(X[atom], Y[atom], Z[atom]));
		copyVector(&eTheta, DCDtoEtheta(X[atom], Y[atom], Z[atom]));
		copyVector(&eRad, DCDtoErad(X[atom], Y[atom], Z[atom]));

		for (j = 0; j < bondsCount[atom]; j++){
			id1 = atom, id2 = bonds[atom * MAX_BONDS + j];
			funcParam.refDist = bondsR0[atom * MAX_BONDS + j];
			copyVector(&atomDistance, DCDtoVector(X[id1] - X[id2], Y[id1] - Y[id2], Z[id1] - Z[id2]));
			copyVector(&atomForce, getFENEDerivative_1(atomDistance));
			sp_1 = eTheta.x * atomDistance.x + eTheta.y * atomDistance.y + eTheta.z * atomDistance.z;
			sp_2 = ePhi.x * atomDistance.x + ePhi.y * atomDistance.y + ePhi.z * atomDistance.z;
			sp_3 = eRad.x * atomDistance.x + eRad.y * atomDistance.y + eRad.z * atomDistance.z;
			atomSphereStress[atom].x += 0.5 * atomForce.value/atomDistance.value * (sp_1 * sp_1 + sp_2 * sp_2); //lateral
			atomSphereStress[atom].y += atomForce.value/atomDistance.value * sp_1 * sp_2; //45-shear
			atomSphereStress[atom].z += atomForce.value/atomDistance.value * sp_3 * sp_3; //radial
			r += 1.0/atomDistance.value;
			r2 += 1.0/(atomDistance.value * atomDistance.value);
		}

		for (j = 0; j < nativeCount[atom]; j++){
			id1 = atom; id2 = native[atom * MAX_NATIVE + j];
			funcParam.refDist = nativeR0[atom * MAX_NATIVE + j];
			funcParam.Eh = nativeEh[atom * MAX_NATIVE + j];
			copyVector(&atomDistance, DCDtoVector(X[id1] - X[id2], Y[id1] - Y[id2], Z[id1] - Z[id2]));
			copyVector(&atomForce, getNativeDerivative_1(atomDistance));
			sp_1 = eTheta.x * atomDistance.x + eTheta.y * atomDistance.y + eTheta.z * atomDistance.z;
			sp_2 = ePhi.x * atomDistance.x + ePhi.y * atomDistance.y + ePhi.z * atomDistance.z;
			sp_3 = eRad.x * atomDistance.x + eRad.y * atomDistance.y + eRad.z * atomDistance.z;
			atomSphereStress[atom].x += 0.5 * atomForce.value/atomDistance.value * (sp_1 * sp_1 + sp_2 * sp_2); //lateral
			atomSphereStress[atom].y += atomForce.value/atomDistance.value * sp_1 * sp_2; //45-shear
			atomSphereStress[atom].z += atomForce.value/atomDistance.value * sp_3 * sp_3; //radial
			r += 1.0/atomDistance.value;
			r2 += 1.0/(atomDistance.value * atomDistance.value);
		}

		for (j = 0; j < pairsCount[atom]; j++){
			id1 = atom; id2 = pairs[atom * MAX_PAIRS + j];
			copyVector(&atomDistance, DCDtoVector(X[id1] - X[id2], Y[id1] - Y[id2], Z[id1] - Z[id2]));
			if (atomDistance.value  < pairs_cutoff && id1 != id2 ){
				copyVector(&atomDistance, DCDtoVector(X[id1] - X[id2], Y[id1] - Y[id2], Z[id1] - Z[id2]));
				copyVector(&atomForce, getRepulsiveDerivative_1(atomDistance));
				sp_1 = eTheta.x * atomDistance.x + eTheta.y * atomDistance.y + eTheta.z * atomDistance.z;
				sp_2 = ePhi.x * atomDistance.x + ePhi.y * atomDistance.y + ePhi.z * atomDistance.z;
				sp_3 = eRad.x * atomDistance.x + eRad.y * atomDistance.y + eRad.z * atomDistance.z;
				atomSphereStress[atom].x += 0.5 * atomForce.value/atomDistance.value * (sp_1 * sp_1 + sp_2 * sp_2); //lateral
				atomSphereStress[atom].y += atomForce.value/atomDistance.value * sp_1 * sp_2; //45-shear
				atomSphereStress[atom].z += atomForce.value/atomDistance.value * sp_3 * sp_3; //radial
				r += 1.0/atomDistance.value;
				r2 += 1.0/(atomDistance.value * atomDistance.value);
			}

		}

		atomArea = Pi * (0.5*r/r2) * (0.5*r/r2);
		atomVolume = 4.0/3.0 * Pi * (0.5*r/r2) * (0.5*r/r2) * (0.5*r/r2);
		atomSphereStress[atom].x /= atomArea;
		atomSphereStress[atom].y /= atomArea;
		atomSphereStress[atom].z /= atomVolume;
	}
}

inline void readSphereStress(){
	printf("Reading input sphere stress values...\n");
	readAtomVectorFrame(tnsrSphereStress, atomSphereStress);
}

inline void writeSphereStress(){
	printf("Print sphere stress tensor...\n");
	outSphereStress = fopen(outSphereStressFilename, "a");
	int i;
	for (i = 0; i < pdbData.atomCount; i++){
		//strcpy(pdbData.atoms[i].segment, pdbData.atoms[i].chain);
		if ( strcmp(pdbData.atoms[i].segment, "") == 0)
			pdbData.atoms[i].segment[0] = pdbData.atoms[i].chain;
		printAtomVectorToFile(outSphereStress, pdbData.atoms[i], atomSphereStress[i]);
	}
	fprintf(outSphereStress, "END\n");
	fclose(outSphereStress);
}

inline void printPDBSphereStress(){
	printf("Printing PDB file with sphere stress distribution...\n");
	int i;
	for (i = 0; i < pdbData.atomCount; i++){
		pdbData.atoms[i].x = X[i];
		pdbData.atoms[i].y = Y[i];
		pdbData.atoms[i].z = Z[i];
		//strcpy(pdbData.atoms[i].segment, pdbData.atoms[i].chain);
		if ( strcmp(pdbData.atoms[i].segment, "") == 0)
			pdbData.atoms[i].segment[0] = pdbData.atoms[i].chain;
		pdbData.atoms[i].occupancy = atomSphereStress[i].x * sphereLateralScale;
		pdbData.atoms[i].beta = atomSphereStress[i].y * sphere45ShearScale;
		pdbData.atoms[i].charge = atomSphereStress[i].z * sphereRadialScale;
	}
	if (cutoffAveOn || cutoffSumOn) makeCutOffAveraged();
	if (segmentAveOn || segmentSumOn) makeSegmentAveraged();
	appendPDB(pdbSphereStressFilename, &pdbData, CONECT);

}

inline void printDATSphereStress(){

}

inline void destroySphereStressTensor(){
	free(atomSphereStress);
	if (tnsrSphereStress != NULL) fclose(tnsrSphereStress);
	//if (outSphereStress != NULL) fclose(outSphereStress);
	//if (pdbSphereStress != NULL) fclose(pdbSphereStress);
}
