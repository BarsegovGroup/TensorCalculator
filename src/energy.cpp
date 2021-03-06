/*
 * energy.c
 *
 *  Created on: Mar 31, 2015
 *      Author: olga
 */

#include "energy.h"
#include "IO/config_reader.h"
#include "IO/tensorio.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

Tensors energy;

double* atomEnergy;

char outEnergyFilename[256];
FILE* outEnergy;

char pdbEnergyFilename[256];
FILE* pdbEnergy;

char datEnergyFilename[256];
FILE* datEnergy;

char tnsrEnergyFilename[256];
FILE* tnsrEnergy;

double energyScale;

int CONECT;

int avEnergy;


void createEnergy(){
	printf("Creating energy structure...\n");
	sprintf(energy.name, "energy");
	energy.compute = &computeEnergy;
	energy.read = &readEnergy;
	energy.write = &writeEnergy;
	energy.printPDB = &printPDBEnergy;
	energy.printDAT = &printDATEnergy;
	energy.destroy = &destroyEnergy;
	energy.runningAverage = &averageEnergy;
	tensors[tensorsCount] = &energy;
	tensorsCount ++;
	initEnergy();
}

void initEnergy(){
	printf("Initializing energy structure...\n");

	atomEnergy = (double*)calloc(pdbData.atomCount, sizeof(double));

	if (computeOn){
		getMaskedParameter(outEnergyFilename, "outputEnergy", "", 0);
		outEnergy = fopen(outEnergyFilename, "w");
		fclose(outEnergy);
	}

	if (printPDBOn){
		if (!computeOn) {
			getMaskedParameter(tnsrEnergyFilename, "inputEnergy", "", 0);
			tnsrEnergy = fopen(tnsrEnergyFilename, "r");
		}
		getMaskedParameter(pdbEnergyFilename, "outputPDBEnergy", "", 0);
		pdbEnergy = fopen(pdbEnergyFilename, "w");
		fclose(pdbEnergy);
		energyScale = getDoubleParameter("energyScaleFactor", 1.0, 1);
		avEnergy = getIntegerParameter("averageEnergy", 1, 1);
	}

	if (printDatOn){
		if (!computeOn) {
			getMaskedParameter(tnsrEnergyFilename, "inputEnergy", "", 0);
			tnsrEnergy = fopen(tnsrEnergyFilename, "r");
		}
		getMaskedParameter(datEnergyFilename, "datEnergy", "", 0);
		datEnergy = fopen(datEnergyFilename, "w");
		fclose(datEnergy);
		energyScale = getDoubleParameter("energyScaleFactor", 1.0, 1);
	}

}

inline void computeEnergy(){
	printf("Computing energy...\n");

	int atom, j, id1, id2;
	Vector atomDistance;
	double totalFENEEnergy, totalNativeEnergy, totalRepulEnergy;
	double tot = 0.0;

	for (atom = 0; atom < pdbData.atomCount; atom++)
		atomEnergy[atom] = 0.0;

	for (atom = 0; atom < pdbData.atomCount; atom++){

		totalFENEEnergy = 0.0;
		totalNativeEnergy = 0.0;
		totalRepulEnergy = 0.0;

		for (j = 0; j < bondsCount[atom]; j++){
			id1 = atom, id2 = bonds[atom * MAX_BONDS + j];
			funcParam.refDist = bondsR0[atom * MAX_BONDS + j];
			copyVector(&atomDistance, DCDtoVector(X[id1] - X[id2], Y[id1] - Y[id2], Z[id1] - Z[id2]));
			totalFENEEnergy += getFENEEnergy(atomDistance.value)/2.0;
		}

		for (j = 0; j < nativeCount[atom]; j++){
			id1 = atom; id2 = native[atom * MAX_NATIVE + j];
			funcParam.refDist = nativeR0[atom * MAX_NATIVE + j];
			funcParam.Eh = nativeEh[atom * MAX_NATIVE + j];
			copyVector(&atomDistance, DCDtoVector(X[id1] - X[id2], Y[id1] - Y[id2], Z[id1] - Z[id2]));
			totalNativeEnergy += getNativeEnergy(atomDistance.value)/2.0;
		}

		for (j = 0; j < pairsCount[atom]; j++){
			id1 = atom; id2 = pairs[atom * MAX_PAIRS + j];
			copyVector(&atomDistance, DCDtoVector(X[id1] - X[id2], Y[id1] - Y[id2], Z[id1] - Z[id2]));
			if (atomDistance.value  < pairs_cutoff && id1 != id2 ){
				totalRepulEnergy += getRepulsiveEnergy(atomDistance.value)/2.0;
			}

		}
		atomEnergy[atom] = totalFENEEnergy + totalNativeEnergy + totalRepulEnergy;

		tot += atomEnergy[atom];
	}

}

inline void readEnergy(){
	printf("Reading input energies...\n");
	readAtomValueFrame(tnsrEnergy, atomEnergy);
}

inline void writeEnergy(){
	printf("Printing energies:\n");
	outEnergy = fopen(outEnergyFilename, "a");
	int i;
	for (i = 0; i < pdbData.atomCount; i++){
		if ( strcmp(pdbData.atoms[i].segment, "") == 0)
			pdbData.atoms[i].segment[0] = pdbData.atoms[i].chain;
		printAtomValueToFile(outEnergy, pdbData.atoms[i], atomEnergy[i]);
	}
	fprintf(outEnergy, "END\n");
	fclose(outEnergy);
}

inline void printPDBEnergy(){
	printf("Printing PDB file with energy distribution...\n");
	int i;
	for (i = 0; i < pdbData.atomCount; i++){
		pdbData.atoms[i].x = X[i];
		pdbData.atoms[i].y = Y[i];
		pdbData.atoms[i].z = Z[i];
		if ( strcmp(pdbData.atoms[i].segment, "") == 0)
			pdbData.atoms[i].segment[0] = pdbData.atoms[i].chain;
		pdbData.atoms[i].beta = 0.0;
		pdbData.atoms[i].charge = 0.0;
		pdbData.atoms[i].occupancy = atomEnergy[i] * energyScale;

	}
	if (cutoffAveOn || cutoffSumOn) makeCutOffAveraged();
	if (segmentAveOn || segmentSumOn) makeSegmentAveraged();
	appendPDB(pdbEnergyFilename, &pdbData, CONECT);
}

inline void printDATEnergy(){
	printf("Printing DAT file with energy distribution...\n");
	int i;
	datEnergy = fopen(datEnergyFilename, "a");

	double totalEnergy = 0.0;
	for (i = 0; i < pdbData.atomCount; i++){
		totalEnergy += atomEnergy[i];
	}

	fprintf(datEnergy, "%f\n", totalEnergy);
	fclose(datEnergy);
}

void destroyEnergy(){
	free(atomEnergy);
	if (tnsrEnergy != NULL) fclose(tnsrEnergy);
}

inline void averageEnergy(){
	if (avEnergy > 1)
		makeRunningAveraged(pdbEnergyFilename, avEnergy);
}

