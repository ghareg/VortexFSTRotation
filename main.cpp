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

void GenerateBesselZeros(double*& BesselZeros, Count& BesselSize);
void GenerateNormBes(double*& normBes, Count& NormBesSize, const double* BesselZeros);
void GenerateSupIntParallel(double*& supInt, Count& supIntSize, const double* BesselZeros, int stage);
void GenerateSupInt(double* supInt, int jStart, int jEnd, Count& lSupIntSize, const double* BesselZeros);
void free(SparseMat& baseHam, InitVal& inVal, double* evalues, MatType* evecs);
void freeN(SparseMat& baseNHam);

struct eigen {
	double value;
	const MatType *vector;
	
	bool operator<(eigen const &other) const {
		return value < other.value;
	}
};

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
	GenMatProd op;
	FILE* evalFile;
	evalFile = fopen("EndepMul10j300R5000xi3C4mi.dat", "w");
	double mu = 40.0;
	double muEnd = 51.0;
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

		std::sort(&evalues[0], &evalues[0] + neigs);

		for (int i = 0; i < neigs; ++i) {
			fprintf(evalFile, "\t%.6f", evalues[i]);
		}
		fprintf(evalFile, "\n");
		mu += dmu;
		fflush(evalFile);
	}
	
	fclose(evalFile);
	free(baseHam, inVal, evalues, evecs);
	freeN(baseFHam);
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


