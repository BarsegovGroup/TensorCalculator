/*
 * tensors.c
 *
 *  Created on: Mar 30, 2015
 *      Author: olga
 */

#include "Tensors/cauchy_stress.h"
#include "Tensors/sphere_stress.h"
#include "Tensors/strain_tensor.h"
#include "energy.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "math.h"

#include "IO/config_reader.h"
#include "IO/topio.h"
#include "IO/dcdio.h"
#include "IO/pdbio.h"
#include "structures.h"

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
int frame;
float* X;
float* Y;
float* Z;

int computeOn;
int printPDBOn;
int printDatOn;

int cauchyOn;
int sphereOn;
int energyOn;
int strainOn;
//int stretchOn;

int cutoffAveOn;
int cutoffSumOn;
int segmentAveOn;
int segmentSumOn;

Tensors **tensors;
int tensorsCount;

void makeCoarseGraining();

int frameInit, frameStride;

int main(int argc, char *argv[]){

	parseFile(argv[1]);

	char pdbFilename[256];
	char topFilename[256];
	char dcdFilename[256];

	int allAtomOn = getYesNoParameter("allAtom", 0, 1);
	if (allAtomOn){
		makeCoarseGraining(); //if input structure is all-atom, extracting Calpha's and creating PDB and DCD
		getMaskedParameter(dcdFilename,"dcdCa", "", 0); //reading created DCD
	} else{
		getMaskedParameter(pdbFilename, "structure", "", 0); //otherwise read coarse-grained DCD and PDB
		readPDB(pdbFilename, &pdbData);  //!!!
		getMaskedParameter(dcdFilename,"dcd", "", 0);
	}

	getMaskedParameter(topFilename, "topology", "", 0); //topology file should be created using sop-top utility
	loadTOP(topFilename); //!!!

	dcdFile = fopen(dcdFilename, "r");//!!!
	dcd_read_header(dcdFile, dcdFilenameOrig, &N, &NFILE, &NPRIV, &NSAVC, &DELTA);
	X = (float*)calloc(pdbData.atomCount, sizeof(float));
	Y = (float*)calloc(pdbData.atomCount, sizeof(float));
	Z = (float*)calloc(pdbData.atomCount, sizeof(float));

	int i;
	tensors = (Tensors**)calloc(MAX_TENSORS, sizeof(Tensors*));
	for (i = 0; i < MAX_TENSORS; i++)
		tensors[i] = (Tensors*)malloc(sizeof(Tensors));

	frameInit = getIntegerParameter("frameInitial", 0, 1);
	frameCount = getIntegerParameter("frameCount", 1, 1);
	frameStride = getIntegerParameter("frameStride", 1, 1);

	computeOn = getYesNoParameter("compute", 1, 1);
	printPDBOn = getYesNoParameter("printPDB", 1, 1);
	printDatOn = getYesNoParameter("printDAT", 0, 1);

	if (computeOn)
		printf("Compute tensors.\n");

	fillContactLists();

	char average[256]; //choice of average values printed in PDB
	cutoffAveOn = 0;
	cutoffSumOn = 0;
	segmentAveOn = 0;
	segmentSumOn = 0;
	if (printPDBOn) {
		printf("Print PDB files.\n");
		getMaskedParameter(average, "average", "no", 1);
		if ( strcmp(average, "no") == 0){
			printf("Printed values will not be averaged.\n");
		}
		if ( strcmp(average, "cutoffAv") == 0){
			cutoffAveOn = 1;
			printf("Printed values will be averaged in cut-off %f A.\n", pairs_cutoff);
		}
		if ( strcmp(average, "cutoffSum") == 0){
			cutoffAveOn = 1;
			cutoffSumOn = 1;
			printf("Printed values will be summed in cut-off %f A.\n", pairs_cutoff);
		}
		if ( strcmp(average, "segmentAv") == 0){
			segmentAveOn = 1;
			printf("Printed values will be averaged for each segment. Make sure, that segment IDs in PDB are specified correctly!\n");
		}
		if ( strcmp(average, "segmentSum") == 0){
			segmentAveOn = 1;
			segmentSumOn = 1;
			printf("Printed values will be summed in chain. Make sure, that segment IDs in PDB are specified correctly!\n");
		}
	}

	if (printDatOn) printf("Print DAT files.\n");

	energyOn = getYesNoParameter("energy", 0, 1);
	cauchyOn = getYesNoParameter("cauchyStress", 0, 1);
	sphereOn = getYesNoParameter("sphereStress", 0, 1);
	strainOn = getYesNoParameter("strain", 0, 1);
	stretchOn = getYesNoParameter("stretch", 0, 1);

	if (energyOn) createEnergy();
	if (cauchyOn) createCauchyStressTensor();
	if (sphereOn) createSphereStressTensor();
	if (strainOn) createStrainTensor(); //makes initialization of DeformationGradient tensor and Stretch tensor if indicated

	char vector[256];
	if (printPDBOn || printDatOn){
		getMaskedParameter(vector, "normalVector", "&", 1);
		if (strcmp(vector, "endtoend") == 0){
			int end1 = getIntegerParameter("end1", 0, 0);
			int end2 = getIntegerParameter("end2", 0, 0);
			copyVector(&normalVec, DCDtoVector(pdbData.atoms[end2].x - pdbData.atoms[end2].x, pdbData.atoms[end2].y - pdbData.atoms[end1].y, pdbData.atoms[end2].z - pdbData.atoms[end1].z));
			printf("Using end-to-end direction for normal vector: (%f, %f, %f)\n", normalVec.x, normalVec.y, normalVec.z);
		}
		else if (strcmp(vector, "vector") == 0){
			float vecX, vecY, vecZ;
			getVectorParameter("normalDirection", &vecX, &vecY, &vecZ, 0, 0, 0, 0);
			normalVec.x = (double)vecX;
			normalVec.y = (double)vecY;
			normalVec.z = (double)vecZ;
			normalVec.value = sqrt(normalVec.x*normalVec.x + normalVec.y*normalVec.y + normalVec.z*normalVec.z);
			printf("Normal vector direction: (%f, %f, %f)\n", normalVec.x, normalVec.y, normalVec.z);
		}
		else{
			printf("First invariant, first principal and von Misses stresses will be calculated.\n");
			normalVec.value = -1.0;
		}
	}

	frame = 0;
	int eof = 0;
	int t = 0;
	CONECT = 1;
	while (!eof && frame < frameInit + frameCount) {
		eof = dcd_read_frame(dcdFile, pdbData.atomCount, X, Y, Z);
		if ((frame >= frameInit) && frame % frameStride == 0){
			printf("Frame %d:\n", frame);
			for (t = 0; t < tensorsCount; t++){
				if (computeOn){
					tensors[t]->compute();
					tensors[t]->write();
				}
				
				if (!computeOn) tensors[t]->read();
				if (printPDBOn){
					tensors[t]->printPDB();
				}
				if (printDatOn){
					tensors[t]->printDAT();
				}
			}
			CONECT = 0.0;
		}

		frame++;
	}
	for (t = 0; t < tensorsCount; t++)
		tensors[t]->destroy();
	
	for (t = 0; t < tensorsCount; t++)
		tensors[t]->runningAverage();

	printf("Done!\n");

	return 0;
}


void makeCoarseGraining(){
	printf("Creating coarse-grained structure...\n");
	char dcdFilename[256];
	char pdbFilename[256];

	PDB pdbAllAtom;
	getMaskedParameter(pdbFilename,"structure", "", 0);
	readPDB(pdbFilename, &pdbAllAtom);

	getMaskedParameter(dcdFilename,"dcd", "", 0);
	FILE* dcdAllAtomFile = fopen(dcdFilename, "r");//!!!
	dcd_read_header(dcdAllAtomFile, dcdFilenameOrig, &N, &NFILE, &NPRIV, &NSAVC, &DELTA);

	X = (float*)calloc(pdbAllAtom.atomCount, sizeof(float));
	Y = (float*)calloc(pdbAllAtom.atomCount, sizeof(float));
	Z = (float*)calloc(pdbAllAtom.atomCount, sizeof(float));

	int i;
	int caCount = 0;
	for (i = 0; i < pdbAllAtom.atomCount; i++)
		if (strcmp(pdbAllAtom.atoms[i].name, "CA") == 0) caCount++;

	pdbData.atoms = (Atom*)calloc(caCount, sizeof(Atom));
	pdbData.atomCount = caCount;

	printf("Found %d amino acid residues.\n", caCount);

	float* caX = (float*)calloc(caCount, sizeof(float));
	float* caY = (float*)calloc(caCount, sizeof(float));
	float* caZ = (float*)calloc(caCount, sizeof(float));

	char dcdOutFilename[256];
	getMaskedParameter(dcdOutFilename,"dcdCa", "", 0);
	FILE* dcdOutFile = fopen(dcdOutFilename, "w");
	dcd_write_header(dcdOutFile, dcdFilenameOrig, pdbData.atomCount, NFILE, NPRIV, NSAVC, DELTA);

	int eof = 0;
	while(!eof){
		eof = dcd_read_frame(dcdAllAtomFile, pdbAllAtom.atomCount, X, Y, Z);

		caCount = 0;
		for (i = 0; i < pdbAllAtom.atomCount; i++){
			if (strcmp(pdbAllAtom.atoms[i].name, "CA") == 0){
				caX[caCount] = X[i];
				caY[caCount] = Y[i];
				caZ[caCount] = Z[i];
				memcpy(&pdbData.atoms[caCount], &pdbAllAtom.atoms[i], sizeof(Atom));
				caCount++;
			}
		}
		dcd_write_frame(dcdOutFile, pdbData.atomCount, caX, caY, caZ);

	}
	fclose(dcdAllAtomFile);
	fclose(dcdOutFile);
	printf("Done creating coarse-grained structure.\n");
}
