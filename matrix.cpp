#include "matrix.h"
#include "constants.h"
#include "operator.h"
#include <thread>
#include <algorithm>
#include <mkl_types.h>
#undef MKL_Complex16
#define MKL_Complex16 MatType
#include <mkl.h>
#include <cstdlib>
#include <fstream>

void calcMat(SparseMat& baseHam, const State* basis, Count bSize, const InitVal& inVal, double mu, double rot)
{
	Count start =  baseHam.start;
	Count end = baseHam.end;
	Count valueSize  = 3 * bSize;
	baseHam.mat = new MatType[valueSize];
	baseHam.ia = new Count[end - start + 1];
	baseHam.ja = new Count[valueSize];

	MatType matel = 0.0f;
	Count currHNum = 0;
	for (Count i = start; i < end; ++i) {
		for (Count j = i; j < bSize; ++j) {
			matel = MatElement(i, j, basis, inVal, mu, rot);
			if (j == i) {
				baseHam.ia[i - start] = currHNum;
			}
			if (j == i || std::abs(matel) > thresh) {
				if (currHNum == valueSize) {
					valueSize *= 2;
					resize(baseHam.mat, currHNum, valueSize);
					resize(baseHam.ja, currHNum, valueSize);
				}
				baseHam.mat[currHNum] = matel;
				baseHam.ja[currHNum] = j;
				++currHNum;
			}
		}
	}
	baseHam.ia[end - start] = currHNum;
}

void initStartEnd(SparseMat* baseHamMat, Count bSize)
{
	double total = 0.5 * (1.0 * bSize * bSize + bSize);	
	Count end0 = (Count) (bSize + 0.5 - 0.5 * sqrt((2.0 * bSize + 1.0) * (2.0 * bSize + 1.0) - 
			   8.0 * total / nCore));
	baseHamMat[0].start = 0;
	baseHamMat[0].end = end0;
	for (int i = 1; i < nCore - 1; ++i) {
		Count oldstart = baseHamMat[i - 1].start;
		Count oldend = baseHamMat[i - 1].end;
		baseHamMat[i].start = oldend;
		Count endi = (Count) (bSize + 0.5 - 0.5 * sqrt((2.0 * bSize + 1.0) * (2.0 * bSize + 1.0) -
				16.0 * bSize * oldend - 8.0 * oldend + 8.0 * oldend * oldend + 8.0 * bSize * 
				oldstart - 4.0 * oldstart * oldstart + 4.0 * oldstart));
		baseHamMat[i].end = endi;
	}
	if (nCore > 1) {
		baseHamMat[nCore - 1].start = baseHamMat[nCore - 2].end;
	}
	baseHamMat[nCore - 1].end = bSize;
}

void joinMat(SparseMat& MatFull, const SparseMat* MatArray, Count bSize, bool zeroBased)
{
	Count valuesSize = 0;
	for (int i = 0; i < nCore; ++i) {
		valuesSize += MatArray[i].ia[MatArray[i].end - MatArray[i].start];
	}
	MatFull.mat = new MatType[valuesSize];
	MatFull.ja = new Count[valuesSize];
	MatFull.ia = new Count[bSize + 1];
	Count jadd = zeroBased ? 0 : 1;
	Count indCount = 0;
	Count iaIndCount = 0;
	Count currTot = zeroBased ? 0 : 1;
	for (int i = 0; i < nCore; ++i) {
		for (Count j = 0; j < MatArray[i].ia[MatArray[i].end - MatArray[i].start]; ++j) {
			MatFull.mat[indCount] = MatArray[i].mat[j];
			MatFull.ja[indCount] = MatArray[i].ja[j] + jadd;
			indCount++;
		}
		for (Count j = 0; j < MatArray[i].end - MatArray[i].start; ++j) {
			MatFull.ia[iaIndCount] = currTot + MatArray[i].ia[j];
			++iaIndCount;
		}
		currTot += MatArray[i].ia[MatArray[i].end - MatArray[i].start];
		
		delete[] MatArray[i].mat;
		delete[] MatArray[i].ja;
		delete[] MatArray[i].ia;
	}
	
	MatFull.ia[bSize] = currTot; 
	MatFull.start = 0;
	MatFull.end = bSize;
}

void calcFullMat(SparseMat& baseHam, const State* basis, Count bSize, const InitVal& inVal, double mu, double rot, bool zeroBased)
{
	SparseMat baseHamMat[nCore];
	initStartEnd(baseHamMat, bSize);
	std::thread t[nCore - 1];
	for (int i = 0; i < nCore - 1; ++i) {
		t[i] = std::thread(calcMat, std::ref(baseHamMat[i]), basis, bSize, std::ref(inVal), mu, rot);
	}
	
	calcMat(baseHamMat[nCore - 1], basis, bSize, std::ref(inVal), mu, rot);

	for (int i = 0; i < nCore - 1; ++i) {
		t[i].join();
	}
	
	joinMat(baseHam, baseHamMat, bSize, zeroBased);
}

void changeFullUMat(SparseMat& baseHam, double dmu, bool zeroBased)
{
	Count Size = baseHam.end - baseHam.start;
	Count add = zeroBased ? 0 : 1;
	#pragma omp parallel for
	for (Count i = 0; i < Size / 2; ++i) {
		baseHam.mat[baseHam.ia[i] - add] += MatType(-dmu, 0.0);
	}
	#pragma omp parallel for
	for (Count i = Size / 2; i < Size; ++i) {
		baseHam.mat[baseHam.ia[i] - add] += MatType(dmu, 0.0);
	}
}

