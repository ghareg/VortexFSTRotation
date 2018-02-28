#ifndef MATRIX_H_
#define MATRIX_H_
#include "basis.h"

struct SparseMat
{
	SparseMat(): mat(nullptr), ia(nullptr), ja(nullptr), start(0), end(0) {}
	MatType* mat;
	Count* ia;
	Count* ja;
	Count start;
	Count end;
};

struct InitVal
{
	InitVal(): BesselZeros(nullptr), normBes(nullptr), soiInt(nullptr), supInt(nullptr) {}
	double* BesselZeros;
	double* normBes;
	double* soiInt;
	double* supInt;
};

template <typename arrayT>
void resize(arrayT*& mat, Count currSize, Count newSize)
{
	arrayT* oldmat = mat;
	mat = new arrayT[newSize];
	for (Count i = 0; i < currSize; ++i) {
		mat[i] = oldmat[i];
	}
	delete[] oldmat;
}

void calcMat(SparseMat& baseHam, const State* basis, Count bSize, const InitVal& inVal, double mu);
void initStartEnd(SparseMat* baseHamMat, Count bSize);
void joinMat(SparseMat& MatFull, const SparseMat* MatArray, Count bSize, bool zeroBased);
void calcFullMat(SparseMat& baseHam, const State* basis, Count bSize, const InitVal& inVal, double mu, bool zeroBased);
void changeFullUMat(SparseMat& baseHam, double dmu, bool zeroBased);
#endif
