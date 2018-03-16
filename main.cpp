#include <cmath>
#include <gsl/gsl_sf_bessel.h>
#include "matrix.h"
#include "constants.h"
#include "operator.h"
#include "diag.h"
#include <thread>
#include <stdio.h>
#include <algorithm>
#include <gsl/gsl_errno.h>

struct eigen {
	double value;
	const MatType *vector;
	
	bool operator<(eigen const &other) const {
		return value < other.value;
	}
};

void GenerateBesselZeros(double*& BesselZeros, Count& BesselSize);
void GenerateNormBes(double*& normBes, Count& NormBesSize, const double* BesselZeros);
void GenerateSupIntParallel(double*& supInt, Count& supIntSize, const double* BesselZeros, int stage);
void GenerateSupInt(double* supInt, int jStart, int jEnd, Count& lSupIntSize, const double* BesselZeros);
void sortEValues(eigen* eval, const double* evalues, const MatType* evecs, Count bSize);
void CalcLz(const int* eign, int nsize, double* Lzn, const eigen* eval, const State* basis, Count bSize, const InitVal& inVal);
void free(SparseMat& baseHam, InitVal& inVal, double* evalues, MatType* evecs);
void freeN(SparseMat& baseNHam);

int main(void)
{
	InitVal inVal;
	Count BesselSize;
	GenerateBesselZeros(inVal.BesselZeros, BesselSize);

	State* basis = nullptr;
	Count bSize;
	GenerateBasis(basis, bSize);
	
	Count NormBesSize;
	GenerateNormBes(inVal.normBes, NormBesSize, inVal.BesselZeros);

	//gsl_set_error_handler_off();
	Count supIntSize;
	GenerateSupIntParallel(inVal.supInt, supIntSize, inVal.BesselZeros, 2);

	SparseMat baseHam;
	SparseMat baseFHam;
	bool zeroBased = false;
	bool matrixCons = false;
	double* evalues = new double[neigs];
	MatType* evecs = new MatType[neigs * bSize];
	eigen* eval = new eigen[neigs];
	GenMatProd op;
	FILE* evalFile;
	evalFile = fopen("EndepMul10j300R5000xi3C4miLz.dat", "w");
	double mu = 48.15;
	double muEnd = 48.15;
	double dmu = 0.025;
	while (mu <= muEnd) {
		fprintf(evalFile, "%.6f\t", mu);
		if (!matrixCons) {
			calcFullMat(baseHam, basis, bSize, inVal, mu, zeroBased);
			//diagtoFullMat(baseFHam, baseHam, zeroBased);
			//writeMatrixToFile(baseFHam, zeroBased);
			op.init(baseHam);
			matrixCons = true;
			//fprintf(evalFile, "Construction done");
			fflush(evalFile);
		}
		else {
			changeFullUMat(baseHam, dmu, zeroBased);
			op.restart(baseHam);
		}
		
		calcEValues(baseHam, op, evalues, evecs);
		sortEValues(eval, evalues, evecs, bSize);

		const int nsize = 6;
		const int eign[nsize] = {neigs / 2 - 3, neigs / 2 - 2, neigs / 2 - 1, neigs / 2, neigs / 2 + 1, neigs / 2 + 2};
		double Lzn[nsize];
		CalcLz(eign, nsize, Lzn, eval, basis, bSize, inVal);


		for (int i = 0; i < neigs; ++i) {
			fprintf(evalFile, "\t%.6f", eval[i].value);
		}
		fprintf(evalFile, "\n");

		for (int i = 0; i < nsize; ++i) {
			fprintf(evalFile, "%.6f\t", eval[eign[i]].value);
		}
		fprintf(evalFile, "\n");

		for (int i = 0; i < nsize; ++i) {
			fprintf(evalFile, "%.6f\t", Lzn[i]);
		}
		fprintf(evalFile, "\n");

		mu += dmu;
		fflush(evalFile);
	}
	
	fclose(evalFile);
	free(baseHam, inVal, evalues, evecs);
	freeN(baseFHam);
	delete[] eval;
}

void sortEValues(eigen* eval, const double* evalues, const MatType* evecs, Count bSize)
{
	for (int i = 0; i < neigs; ++i) {
		eval[i].value = evalues[i];
		eval[i].vector = &evecs[i * bSize];
	}
	std::sort(&eval[0], &eval[0] + neigs);
}

void CalcLz(const int* eign, int nsize, double* Lzn, const eigen* eval, const State* basis, Count bSize, const InitVal& inVal)
{
	for (int i = 0; i < nsize; ++i) {
		Lzn[i] = 0.0;
	}
	int Lz[16] = {1, 2, 0, 1, 0, 1, -1, 0, 0, 1, -1, 0, -1, 0, -2, -1};

	for (int k = 0; k < nsize; ++k) {
		for (Count i = 0; i < bSize; i += 1) {
			Occup ph = basis[i].ph();
			Occup s = basis[i].s();
			Occup orb = basis[i].orb();
			Occup l = basis[i].l();
			Lzn[k] += (l + Lz[(ph + 1) * 4 + (s + 1) * 2 + orb]) * std::norm(eval[eign[k]].vector[i]);
		}
	}
}

void GenerateBesselZeros(double*& BesselZeros, Count& BesselSize)
{
	BesselSize = 0;
	int j0[3] = {558, 3150, 4916};
	BesselZeros = new double[(lACount + 1) * jmax];
	for (int l = 0; l <= lmax + 1; ++l) {
		for (int ind = 0; ind < 3; ++ind) {
			for (int j = j0[ind] - jsmax / 2; j < j0[ind] + jsmax / 2; ++j) {
				BesselZeros[BesselSize] = gsl_sf_bessel_zero_Jnu(l, j);
				BesselSize++;
			}
		}
	}
}

void GenerateNormBes(double*& normBes, Count& NormBesSize, const double* BesselZeros)
{
	NormBesSize = 0;
	normBes = new double[lACount * jmax];
	for (int l = 0; l <= lmax; ++l) {
		for (int j = 0; j < jmax; ++j) {
			normBes[l * jmax + j] = std::abs(std::sqrt(2.0f) / gsl_sf_bessel_Jn(l + 1, BesselZeros[l * jmax + j]));
			NormBesSize++;
		}
	}
}

void GenerateSupIntParallel(double*& supInt, Count& supIntSize, const double* BesselZeros, int stage)
{
	if (stage == 0 || stage == 1) {
		supIntSize = 0;
		supInt = new double[3 * lACount * jmax * jmax];
		Count lSupIntSize[nCore];
		std::thread t[nCore - 1];
		for (int i = 0; i < nCore - 1; ++i) {
			t[i] = std::thread(GenerateSupInt, supInt, jmax * i / nCore, jmax * (i + 1) / nCore,  std::ref(lSupIntSize[i]), BesselZeros);
		}

		GenerateSupInt(supInt, jmax * (nCore - 1) / nCore, jmax, lSupIntSize[nCore - 1], BesselZeros);

		for (int i = 0; i < nCore - 1; ++i) {
			t[i].join();
		}

		for (int i = 0; i < nCore; ++i) {
			supIntSize += lSupIntSize[i];
		}

		if (stage == 1) {
			FILE* supfile;
			supfile = fopen("SupIntl10j500R5000xi3RC.dat", "w");
			int istart = 0;
			for (int i = istart; i < istart + supIntSize; ++i) {
				fprintf(supfile, "%.8f\t", supInt[i]);
			}
			fclose(supfile);
		}
	}
	else {
		const int mlACount = lSupMax + 1;
		supIntSize = 3 * mlACount * jmax * jmax;
		supInt = new double[3 * mlACount * jmax * jmax];
		FILE* supfile;
		supfile = fopen("/scratch/areg.ghazaryan/VortexMat/SupIntl10j300R5000xi3T.dat", "r");
		int rCount = 0;
		while(rCount < supIntSize && fscanf(supfile, "%lf", &supInt[rCount]) != EOF) {
			rCount++;
		}
		supIntSize = rCount;
		fclose(supfile);
	}
}

void GenerateSupInt(double* supInt, int jStart, int jEnd, Count& lSupIntSize, const double* BesselZeros)
{
	lSupIntSize = 0;
	for (int l = 0; l <= lmax; ++l) {
		for (int j = jStart; j < jEnd; ++j) {
			for (int jp = 0; jp < jmax; ++jp) {
				int lp = l + 1;
				supInt[3 * l * jmax * jmax + 3 * j * jmax + 3 * jp] = IntTan(l - 1, lp - 1, BesselZeros[l * jmax + j], BesselZeros[lp * jmax + jp]);
				supInt[3 * l * jmax * jmax + 3 * j * jmax + 3 * jp + 1] = IntTan(l, lp, BesselZeros[l * jmax + j], BesselZeros[lp * jmax + jp]);
				supInt[3 * l * jmax * jmax + 3 * j * jmax + 3 * jp + 2] = IntTan(l + 1, lp + 1, BesselZeros[l * jmax + j], BesselZeros[lp * jmax + jp]);
				lSupIntSize += 3;
			}
		}
	}
}

void free(SparseMat& baseHam, InitVal& inVal, double* evalues, MatType* evecs)
{
	delete[] baseHam.mat;
	delete[] baseHam.ia;
	delete[] baseHam.ja;

	delete[] inVal.BesselZeros;
	delete[] inVal.normBes;
	delete[] inVal.supInt;

	delete[] evalues;
	delete[] evecs;
}

void freeN(SparseMat& baseNHam)
{
	delete[] baseNHam.mat;
	delete[] baseNHam.ia;
	delete[] baseNHam.ja;
}


