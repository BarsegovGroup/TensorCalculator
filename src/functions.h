/*
 * functions.h
 *
 *  Created on: Mar 30, 2015
 *      Author: olga
 */

#ifndef FUNCTIONS_H_
#define FUNCTIONS_H_

#include "structures.h"
//#include "IO/config_reader.h"
//#include "IO/tensorio.h"

extern Vector getFENEDerivative_1(Vector dist);
extern Vector getNativeDerivative_1(Vector dist);
extern Vector getRepulsiveDerivative_1(Vector dist);

extern double getFENEEnergy(double dist);
extern double getNativeEnergy(double dist);
extern double getRepulsiveEnergy(double dist);

extern void copyVector(Vector* v1, const Vector v2);
extern Vector nullVector();
extern void printVector(Vector v);
extern double scalarProduct(Vector v1, Vector v2);
extern double getVectorComponent(const Vector v, int num);

extern Vector DCDtoVector(float x, float y, float z);
extern Vector DCDtoErad(float x, float y, float z);
extern Vector DCDtoEtheta(float x, float y, float z);
extern Vector DCDtoEphi(float x, float y, float z);

//void fillContactLists();
//void calculatePairs();

void calculateNormalShearComponent(const double* tensor, Vector n, double* normal, double* shear);

extern void inverseMatrix3d(const double* m_init, double* m_fin);
extern void multVectors(const Vector A, const Vector B, double* C);
extern void multMatrix(const double* A, const double* B, double* C);
extern void multScalarMatrix(double* A, const double B);
extern void addMatrix(double* m_init, const double* m_fin);
extern void transposeMatrix(const double* m_init, double* m_fin);

extern void eigenDecompos3d(const double* m_init, double* e_val, double* e_vec);

extern void makeCutOffAveraged();
extern void makeSegmentAveraged();

#endif /* FUNCTIONS_H_ */
