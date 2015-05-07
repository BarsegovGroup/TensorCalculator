/*
 * functions.h
 *
 *  Created on: Mar 30, 2015
 *      Author: olga
 */

#ifndef FUNCTIONS_H_
#define FUNCTIONS_H_

#include "structures.h"
#include "IO/config_reader.h"
#include "math.h"
#include "float.h"
#include <string.h>

extern Vector getFENEDerivative_1(Vector dist);
extern Vector getNativeDerivative_1(Vector dist);
extern Vector getRepulsiveDerivative_1(Vector dist);

extern double getFENEEnergy(double dist);
extern double getNativeEnergy(double dist);
extern double getRepulsiveEnergy(double dist);

extern inline void copyVector(Vector* v1, const Vector v2);
extern inline Vector nullVector();
extern inline void printVector(Vector v);
extern inline double scalarProduct(Vector v1, Vector v2);
extern inline double getVectorComponent(const Vector v, int num);

inline Vector DCDtoVector(float x, float y, float z);
inline Vector DCDtoErad(float x, float y, float z);
inline Vector DCDtoEtheta(float x, float y, float z);
inline Vector DCDtoEphi(float x, float y, float z);

void fillContactLists();
void calculatePairs();

void calculateNormalShearComponent(const double* tensor, Vector n, double* normal, double* shear);

inline void inverseMatrix3d(const double* m_init, double* m_fin);
inline void multVectors(const Vector A, const Vector B, double* C);
inline void multMatrix(const double* A, const double* B, double* C);
inline void multScalarMatrix(double* A, const double B);
inline void addMatrix(double* m_init, const double* m_fin);
inline void transposeMatrix(const double* m_init, double* m_fin);

extern void eigenDecompos3d(const double* m_init, double* e_val, double* e_vec);

extern void makeCutOffAveraged();
extern void makeSegmentAveraged();

#endif /* FUNCTIONS_H_ */
