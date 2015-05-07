/*
 * pdbio.h
 *
 *  Created on: Jul 5, 2009
 *      Author: zhmurov
 */

#ifndef PDBIO_H_
#define PDBIO_H_

#include <stdio.h>

/*
 * Structures
 */
typedef struct{

  int id;
  char   name[5], chain, resName[4], altLoc;
  int    resid;
  float x, y, z;

  float occupancy;
  float beta;
  float charge;

  char segment[5];

} Atom;

typedef struct{

	int resid1;
	char chain1;

	int resid2;
	char chain2;

} SSBond;

typedef struct{
	int i;
	int j;
	float r0;
} CovalentBond;

typedef struct{
	int i;
	int j;
	float r0;
	float eh;
} NativeContact;

typedef struct{
	int i;
	int j;
} PossiblePair;

typedef struct{
	int aminoCount;
	int bondCount;
	int nativeCount;
	int pairCount;
	Atom* aminos;
	CovalentBond* bonds;
	NativeContact* natives;
	PossiblePair* pairs;

	int additionalAminosCount;
	Atom* additionalAminos;
} SOP;

typedef struct {
	int atomCount;
	Atom* atoms;
	int ssCount;
	SSBond* ssbonds;
} PDB;

SOP sop;
/*
 * Public methods
 */
void readPDB(const char* filename, PDB* pdbData);
void writePDB(const char* filename, PDB* pdbData);
void appendPDB(const char* filename, PDB* pdbData, int connect);
void printAtom(Atom atomData);
void printAtomToFile(FILE* file, Atom atomData);


#endif /* PDBIO_H_ */
