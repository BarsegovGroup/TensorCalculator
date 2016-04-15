/*
 * functions.c
 *
 *  Created on: Mar 30, 2015
 *      Author: olga
 */

#include "functions.h"

#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "math.h"
#include "float.h"

FuncParam funcParam;

Vector getFENEDerivative_1(Vector dist){
	Vector force;
	copyVector(&force, nullVector());
	double ref = (dist.value - funcParam.refDist) / R_limit;
	force.value = kspring_cov * ref * R_limit/ (1.0 - ref * ref);
	force.x = force.value * dist.x / dist.value;
	force.y = force.value * dist.y / dist.value;
	force.z = force.value * dist.z / dist.value;

	return force;
}

double getFENEEnergy(double dist){
	double energy = 0.0;
	double ref = (dist - funcParam.refDist) / R_limit;
	energy = -0.5 * kspring_cov * R_limit * R_limit * log(1.0 - ref * ref);

	return energy;
}

Vector getNativeDerivative_1(Vector dist){
	Vector force;
	copyVector(&force, nullVector());
	double ref = funcParam.refDist / dist.value; // (r0/r)^2
	ref = ref * ref * ref * ref * ref *ref;
	force.value = -12.0 * funcParam.Eh / dist.value * (ref * ref - ref);
	force.x = force.value * dist.x / dist.value;
	force.y = force.value * dist.y / dist.value;
	force.z = force.value * dist.z / dist.value;
	return force;
}

double getNativeEnergy(double dist){
	double energy = 0.0;
	double ref = funcParam.refDist / dist; // (r0/r)
	ref = ref * ref * ref * ref * ref *ref;
	energy = funcParam.Eh  * (ref * ref - 2.0 * ref);
	return energy;
}

Vector getRepulsiveDerivative_1(Vector dist){
	Vector force;
	copyVector(&force, nullVector());
	double ref = Sigma * Sigma / (dist.value * dist.value);
	force.value = -6.0 * El * ref * ref * ref/ dist.value;
	force.x = force.value * dist.x / dist.value;
	force.y = force.value * dist.y / dist.value;
	force.z = force.value * dist.z / dist.value;
	return force;
}

double getRepulsiveEnergy(double dist){
	double energy = 0.0;
	double ref = Sigma * Sigma / (dist * dist);
	energy = El * ref * ref * ref;
	return energy;
}

void copyVector(Vector* v1, const Vector v2){
	v1->x = v2.x;
	v1->y = v2.y;
	v1->z = v2.z;
	v1->value = v2.value;
}

Vector nullVector(){
	Vector null;
	null.x = 0.0;
	null.y = 0.0;
	null.z = 0.0;
	null.value = 0.0;
	return null;
}

void printVector(Vector v){
	printf("%5.3f  %5.3f  %5.3f\n", v.x, v.y, v.z);
}

double scalarProduct(Vector v1, Vector v2){
	return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

Vector DCDtoVector(float x, float y, float z){
	Vector r;
	r.x = (double)x;
	r.y = (double)y;
	r.z = (double)z;
	r.value = sqrt(x * x + y * y + z * z);
	return r;
}

double getVectorComponent(const Vector v, int num){
	if (num == 0)
		return v.x;
	if (num == 1)
		return v.y;
	if (num == 2)
		return v.z;
	if (num > 2 && num < 0){
		printf("Wrong component number\n");
		exit(0);
	}
}

Vector DCDtoEtheta(float x, float y, float z){
	Vector e_theta;
	double rad = 0;
	double theta = 0, psi = 0;
	rad = sqrt(x * x + y * y + z * z);
	theta = acos(z/rad);
	psi = atan(y/x);
	e_theta.x = cos(theta) * cos(psi);
	e_theta.y = cos(theta) * sin(psi);
	e_theta.z = -sin(theta);
	e_theta.value = sqrt(e_theta.x * e_theta.x + e_theta.y * e_theta.y + e_theta.z * e_theta.z);
	return e_theta;
}

Vector DCDtoEphi(float x, float y, float z){
	Vector e_phi;
	//double rad = 0;
	//double theta = 0
	double phi = 0;
	//rad = sqrtf(x * x + y * y + z * z);
	//theta = acos(z/rad);
	phi = atan(y/x);
	e_phi.x = -sin(phi);
	e_phi.y = cos(phi);
	e_phi.z = 0;
	e_phi.value = sqrt(e_phi.x * e_phi.x + e_phi.y * e_phi.y + e_phi.z * e_phi.z);
	return e_phi;
}

Vector DCDtoErad(float x, float y, float z){
	Vector e_rad;
	double rad = 0;
	double theta = 0, psi = 0;
	rad = sqrt(x * x + y * y + z * z);
	theta = acos(z/rad);
	psi = atan(y/x);
	e_rad.x = sin(theta) * cos(psi);
	e_rad.y = sin(theta) * sin(psi);
	e_rad.z = cos(theta);
	e_rad.value = sqrt(e_rad.x * e_rad.x + e_rad.y * e_rad.y + e_rad.z * e_rad.z);
	return e_rad;

}

void calculateNormalShearComponent(const double* tensor, Vector n, double* normal, double* shear){
	double Px, Py, Pz;

	Px = n.x * tensor[0] + n.y * tensor[3] + n.z * tensor[6];
	Py = n.x * tensor[1] + n.y * tensor[4] + n.z * tensor[7];
	Pz = n.x * tensor[2] + n.y * tensor[5] + n.z * tensor[8];

	*normal = n.x * Px + n.y * Py + n.z * Pz;

	double sqShear;
	sqShear = Px * Px + Py * Py  + Pz * Pz - (*normal) * (*normal);
	if (sqShear <= 0) *shear = 0.0f;
	else *shear = sqrtf(sqShear);
}

void inverseMatrix3d(const double* m_init, double* m_fin){
	double detM = 0.0f;
	int i;
	for (i = 0; i < 9; i++)
		m_fin[i] = 0.0f;

	double A = m_init[4]*m_init[8] - m_init[5] * m_init[7];
	double B = m_init[5]*m_init[6] - m_init[3] * m_init[8];
	double C = m_init[3]*m_init[7] - m_init[4] * m_init[6];
	double D = m_init[2]*m_init[7] - m_init[1] * m_init[8];
	double E = m_init[0]*m_init[8] - m_init[2] * m_init[6];
	double F = m_init[6]*m_init[1] - m_init[0] * m_init[7];
	double G = m_init[1]*m_init[5] - m_init[2] * m_init[4];
	double H = m_init[2]*m_init[3] - m_init[0] * m_init[5];
	double K = m_init[0]*m_init[4] - m_init[1] * m_init[3];

	detM = m_init[0] * A + m_init[1] * B + m_init[2] * C;

	m_fin[0] = 1.0f/detM * A;
	m_fin[1] = 1.0f/detM * D;
	m_fin[2] = 1.0f/detM * G;
	m_fin[3] = 1.0f/detM * B;
	m_fin[4] = 1.0f/detM * E;
	m_fin[5] = 1.0f/detM * H;
	m_fin[6] = 1.0f/detM * C;
	m_fin[7] = 1.0f/detM * F;
	m_fin[8] = 1.0f/detM * K;
}

void multVectors(const Vector A, const Vector B, double* C){
	C[0] = A.x * B.x;
	C[1] = A.x * B.y;
	C[2] = A.x * B.z;
	C[3] = A.y * B.x;
	C[4] = A.y * B.y;
	C[5] = A.y * B.z;
	C[6] = A.z * B.x;
	C[7] = A.z * B.y;
	C[8] = A.z * B.z;
}

void addMatrix(double* m_fin, const double* m_init){
	int i;
	for (i = 0; i < 9; i++)
		m_fin[i] += m_init[i];
}

void multMatrix(const double* A, const double* B, double* C){
	int i, j, k;
	double Cprime[9];
	for (i = 0; i < 9; i++)
		Cprime[i] = 0.0f;
	for (i = 0; i < 3; i++){
		for (j = 0; j < 3; j++){
			for (k = 0; k < 3; k++){
				Cprime[i * 3 + j] += A[i * 3 + k] * B[k * 3 + j];
			}
		}
	}
	for (i = 0; i < 9; i++)
		C[i] = Cprime[i];
}

void transposeMatrix(const double* m_init, double* m_fin){
	int i = 0;
	for (i = 0; i < 9; i++)
		m_fin[i] = 0.0f;
	m_fin[0] = m_init[0];
	m_fin[1] = m_init[3];
	m_fin[2] = m_init[6];
	m_fin[3] = m_init[1];
	m_fin[4] = m_init[4];
	m_fin[5] = m_init[7];
	m_fin[6] = m_init[2];
	m_fin[7] = m_init[5];
	m_fin[8] = m_init[8];
}

void multScalarMatrix(double* A, const double B){
	int i;
	for (i = 0; i < 9; i++)
		A[i] *= B;
}

void eigenDecompos3d(const double* m_init, double* e_val, double* e_vec){
	double a, b, c, d, e, f;
	a = m_init[0];
	b = m_init[4];
	c = m_init[8];
	d = m_init[1];
	e = m_init[2];
	f = m_init[5];

	int i = 0;
	for (i = 0; i < 9; i++){
		e_vec[i] = 0.0f;
		e_val[i] = 0.0f;
	}

	double AA, BB, CC;
	AA = -(a + b + c);
	BB = -(e * e + f * f + d * d - a * b - b * c - a * c);
	CC = -(a * b * c + 2.0f * d * f * e - e * e * b - f * f * a - d * d * c);

	//printf("Characteristic polynomial: l^3 + %f*l^2 + %f*l + %f = 0.\n", AA, BB, CC);

	double QQ, RR;
	QQ = AA * AA/9.0 - BB/3.0;
	RR = AA * AA * AA/27.0 - AA * BB/6.0 + CC/2.0;

	if (QQ < 0)	printf("Trigonometric constants: Q = %f, R = %f.\n", QQ, RR);

	double phi, l1, l2, l3;

	if (QQ < 0) phi = 0.0;
	else phi = acos(RR/sqrt(QQ * QQ * QQ))/3.0;
	l1 = -2.0 * sqrt(QQ) * cos(phi) - AA/3.0;
	l2 = -2.0 * sqrt(QQ) * cos(phi + 2.0 * Pi/3.0) - AA/3.0;
	l3 = -2.0 * sqrt(QQ) * cos(phi - 2.0 * Pi/3.0) - AA/3.0;
	//printf("Roots of characteristic polynomial: (%f, %f, %f).\n", l1, l2, l3);

	double MM = d/(a - l1);
	double LL = e/(a - l1);
	double NN = (f/d - LL)/((b - l1)/d - MM);

	double x1, y1, z1, norm;
	z1 = 1.0f;
	y1 = -1.0f * NN * z1;
	x1 = -1.0f * MM * y1 - LL * z1;
	norm = sqrtf(x1 * x1 + y1 * y1 + z1 * z1);
	x1 /= norm;
	y1 /= norm;
	z1 /= norm;
	//printf("Eigenvalue = %f, eigenvector = (%f, %f, %f)\n", l1, x1, y1, z1);

	MM = d/(a - l2);
	LL = e/(a - l2);
	NN = (f/d - LL)/((b - l2)/d - MM);

	double x2, y2, z2;
	z2 = 1.0f;
	y2 = -1.0f * NN * z2;
	x2 = -1.0f * MM * y2 - LL * z2;
	norm = sqrtf(x2 * x2 + y2 * y2 + z2 * z2);
	x2 /= norm;
	y2 /= norm;
	z2 /= norm;
	//printf("Eigenvalue = %f, eigenvector = (%f, %f, %f)\n", l2, x2, y2, z2);

	MM = d/(a - l3);
	LL = e/(a - l3);
	NN = (f/d - LL)/((b - l3)/d - MM);

	double x3, y3, z3;
	z3 = 1.0f;
	y3 = -1.0f * NN * z3;
	x3 = -1.0f * MM * y3 - LL * z3;
	norm = sqrtf(x3 * x3 + y3 * y3 + z3 * z3);
	x3 /= norm;
	y3 /= norm;
	z3 /= norm;
	//printf("Eigenvalue = %f, eigenvector = (%f, %f, %f)\n", l3, x3, y3, z3);

	e_vec[0] = x1;
	e_vec[1] = x2;
	e_vec[2] = x3;
	e_vec[3] = y1;
	e_vec[4] = y2;
	e_vec[5] = y3;
	e_vec[6] = z1;
	e_vec[7] = z2;
	e_vec[8] = z3;

	e_val[0] = l1;
	e_val[4] = l2;
	e_val[8] = l3;

}

void makeCutOffAveraged(){
	double* occupancy = (double*)calloc(pdbData.atomCount, sizeof(double));
	double* beta = (double*)calloc(pdbData.atomCount, sizeof(double));
	double* charge = (double*)calloc(pdbData.atomCount, sizeof(double));
	int i, j, id;

	for (i = 0; i < pdbData.atomCount; i++){
		occupancy[i] = pdbData.atoms[i].occupancy;
		pdbData.atoms[i].occupancy = 0.0;
		beta[i] = pdbData.atoms[i].beta;
		pdbData.atoms[i].beta = 0.0;
		charge[i] = pdbData.atoms[i].charge;
		pdbData.atoms[i].charge = 0.0;
	}
	int count;
	for (i = 0; i < pdbData.atomCount; i++){
		count = 0;
		for (j = 0; j < bondsCount[i]; j++){
			id = bonds[i * MAX_BONDS + j];
			pdbData.atoms[i].occupancy += occupancy[id];
			pdbData.atoms[i].beta += beta[id];
			pdbData.atoms[i].charge += charge[id];
			count++;
		}
		for(j = 0; j < nativeCount[i]; j++){
			id = native[i * MAX_NATIVE + j];
			pdbData.atoms[i].occupancy += occupancy[id];
			pdbData.atoms[i].beta += beta[id];
			pdbData.atoms[i].charge += charge[id];
			count++;
		}
		for (j = 0; j < pairsCount[i]; j++){
			id = pairs[i * MAX_PAIRS + j];
			pdbData.atoms[i].occupancy += occupancy[id];
			pdbData.atoms[i].beta += beta[id];
			pdbData.atoms[i].charge += charge[id];
			count++;
		}
		if (!cutoffSumOn){
			pdbData.atoms[i].occupancy /= count;
			pdbData.atoms[i].beta /= count;
			pdbData.atoms[i].charge /= count;
		}
	}
}

void makeSegmentAveraged(){
	double* occupancy = (double*)calloc(pdbData.atomCount, sizeof(double));
	double* beta = (double*)calloc(pdbData.atomCount, sizeof(double));
	double* charge = (double*)calloc(pdbData.atomCount, sizeof(double));
	int i, j;
	
	for (i = 0; i < pdbData.atomCount; i++){
		occupancy[i] = pdbData.atoms[i].occupancy;
		pdbData.atoms[i].occupancy = 0.0;
		beta[i] = pdbData.atoms[i].beta;
		pdbData.atoms[i].beta = 0.0;
		charge[i] = pdbData.atoms[i].charge;
		pdbData.atoms[i].charge = 0.0;
	}
	int count;
	for (i = 0; i < pdbData.atomCount; i++){
		count = 0;
		for (j = 0; j < pdbData.atomCount; j++){
			if (strcmp(pdbData.atoms[i].segment, pdbData.atoms[j].segment) == 0){
				pdbData.atoms[i].occupancy += occupancy[j];
				pdbData.atoms[i].beta += beta[j];
				pdbData.atoms[i].charge += charge[j];
				count++;
			}
		}
		if (!segmentSumOn){
			pdbData.atoms[i].occupancy /= count;
			pdbData.atoms[i].beta /= count;
			pdbData.atoms[i].charge /= count;
		}
	}
}
