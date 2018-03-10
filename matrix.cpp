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

void calcMat(SparseMat& baseHam, const State* basis, Count bSize, const InitVal& inVal, double mu)
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
			matel = MatElement(i, j, basis, inVal, mu);
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

void calcFullMat(SparseMat& baseHam, const State* basis, Count bSize, const InitVal& inVal, double mu, bool zeroBased)
{
	SparseMat baseHamMat[nCore];
	initStartEnd(baseHamMat, bSize);
	std::thread t[nCore - 1];
	for (int i = 0; i < nCore - 1; ++i) {
		t[i] = std::thread(calcMat, std::ref(baseHamMat[i]), basis, bSize, std::ref(inVal), mu);
	}
	
	calcMat(baseHamMat[nCore - 1], basis, bSize, std::ref(inVal), mu);

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

void diagtoFullMat(SparseMat& MatFull, SparseMat& MatUFull, bool zeroBased)
{
	const int add = zeroBased ? 0 : 1;
	const Count nimbsize = MatUFull.end - MatUFull.start;
	const Count USize = MatUFull.ia[nimbsize] - add;
	SparseMat MatDFull;
	MatDFull.mat = new MatType[USize];
	MatDFull.ja = new Count[USize];
	MatDFull.ia = new Count[nimbsize + 1];
	MKL_INT job[6];
	job[0] = 1;
	job[1] = add;
	job[2] = add;
	job[5] = 1;
	const MKL_INT n = nimbsize;
	MKL_INT info;
	mkl_zcsrcsc (job, &n , MatDFull.mat, MatDFull.ja , MatDFull.ia , MatUFull.mat, MatUFull.ja , MatUFull.ia , &info);
	for (Count i = 0; i < USize; ++i) {
		MatDFull.mat[i] = std::conj(MatDFull.mat[i]);
	}
	for (Count i = 0; i < nimbsize; ++i) {
		MatUFull.mat[MatUFull.ia[i] - add] = MatType(0.0, 0.0);
	}
	Count Size = 2 * (USize - nimbsize) + nimbsize;
	MatFull.mat = new MatType[Size];
	MatFull.ja = new Count[Size];
	MatFull.ia = new Count[nimbsize + 1];
	MatFull.start = MatUFull.start;
	MatFull.end = MatUFull.end;
	const char trans = 'N';
	MKL_INT request = 0;
	MatType beta = MatType(1.0, 0.0);
	const MKL_INT nzmax = Size;
	MKL_INT sort = 0;
	mkl_zcsradd (&trans, &request, &sort, &n , &n , MatUFull.mat, MatUFull.ja, MatUFull.ia, &beta, MatDFull.mat, MatDFull.ja, MatDFull.ia, MatFull.mat, MatFull.ja, MatFull.ia, &nzmax, &info);
	/*delete[] MatUFull.mat;
	delete[] MatUFull.ia;
	delete[] MatUFull.ja;*/
	delete[] MatDFull.mat;
	delete[] MatDFull.ia;
	delete[] MatDFull.ja;
}

void writeMatrixToFile(const SparseMat& baseHam, bool zeroBased)
{
	std::ofstream matfile;
	std::ofstream iafile;
	std::ofstream jafile;
	matfile.open("mat.dat");
	iafile.open("ia.dat");
	jafile.open("ja.dat");
	Count add = zeroBased ? 0 : 1;
	for (Count i = 0; i < baseHam.ia[baseHam.end - baseHam.start] - add; ++i) {
		if (baseHam.mat[i].imag() < 0) {
			matfile << baseHam.mat[i].real() << baseHam.mat[i].imag() << "j" << std::endl;
		}
		else {
			matfile << baseHam.mat[i].real() << "+" << baseHam.mat[i].imag() << "j" << std::endl;
		}
		jafile << baseHam.ja[i] - add << std::endl;
	}

	for (Count i = 0; i <= baseHam.end - baseHam.start; ++i) {
		iafile << baseHam.ia[i] - add << std::endl;
	}
	matfile.close();
	iafile.close();
	jafile.close();
}

Count DiagValueIndex(Count rowNum, const SparseMat& baseHam, bool zeroBased)
{
	Count add = zeroBased ? 0 : 1;
	Count rowStart = baseHam.ia[rowNum] - add;
	Count rowEnd = baseHam.ia[rowNum + 1] - add;
	Count* posStart = &baseHam.ja[rowStart];
	Count* posEnd = &baseHam.ja[rowEnd];
	Count* pos = std::lower_bound(posStart, posEnd, rowNum + add);
	if (*pos == rowNum + add && pos != posEnd) {
		return rowStart + (pos - posStart);
	}

	return -1;
}

void changeFullMat(SparseMat& baseHam, double dmu, bool zeroBased)
{
	Count Size = baseHam.end - baseHam.start;
	#pragma omp parallel for
	for (Count i = 0; i < Size / 2; ++i) {
		baseHam.mat[DiagValueIndex(i, baseHam, zeroBased)] += MatType(-dmu, 0.0);
	}
	#pragma omp parallel for
	for (Count i = Size / 2; i < Size; ++i) {
		baseHam.mat[DiagValueIndex(i, baseHam, zeroBased)] += MatType(dmu, 0.0);
	}
}

