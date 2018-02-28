#include "diag.h"

extern "C" {

void znaupd_(int* ido, char* bmat, int* n, char* which, int* nev, double* tol, void* resid, 
		int* ncv, void* v, int* ldv, int* iparam, int* ipntr, void* workd,
	   	void* workl, int* lworkl, double* rwork, int* info);

void zneupd_(int* rvec, char* A, int* select, void* d, void* z, int* ldz, void* sigma,
		void* workev, char* bmat, int* n, char* which, int* nev, double* tol, void* resid, int* ncv, 
		void* v, int* ldv, int* iparam, int* ipntr, void* workd, void* workl, int* lworkl,
		double* rwork, int* ierr);
}

GenMatProd::GenMatProd()  
{
	pt_ = (_MKL_DSS_HANDLE_t) new MKL_INT[128];
	mtype_ = -4;
	iparm_ = new MKL_INT[128];
	perm_ = new MKL_INT[baseHam_.end - baseHam_.start];
	pardisoinit(pt_, &mtype_, iparm_);
	iparm_[34] = 0;
	iparm_[26] = 1;
	iparm_[27] = 0;
	iparm_[23] = 1;
	set_shift(shift);
}

void GenMatProd::init(SparseMat baseHam)
{
	baseHam_ = baseHam;
	qsize_ = baseHam.end - baseHam.start;
	phase_ = 13;
}

void GenMatProd::restart(SparseMat baseHam)
{
	baseHam_ = baseHam;
	qsize_ = baseHam.end - baseHam.start;
	phase_ = 23;
}
	

void GenMatProd::perform_op(const MatType* xIn, MatType* yOut)
{
	MKL_INT maxfct = 1;
	MKL_INT mnum = 1;
	MKL_INT n = qsize_;
	MKL_INT nrhs = 1;
	MKL_INT msglvl = 0;
	MKL_INT error = 0;

	pardiso(pt_, &maxfct, &mnum, &mtype_, &phase_, &n, (void*) baseHam_.mat, baseHam_.ia, baseHam_.ja, 
			perm_, &nrhs, iparm_, &msglvl, (void*) xIn, (void*) yOut, &error);
	phase_ = 33;
}

GenMatProd::~GenMatProd()
{
	phase_ = -1;
	perform_op(NULL, NULL);
	delete[] perm_;
	delete[] (MKL_INT* )pt_;
	delete[] iparm_;

}

void calcEValues(const SparseMat& baseHam, GenMatProd& op, double* evalues, MatType* evecs)
{
	Count maxn = baseHam.end - baseHam.start; 	
	int maxncv = 4 * neigs + 1;
	int ldv = maxn;
	int iparam[11] = {0};
	int ipntr[14];
	int select[maxncv];
	MatType* d = new MatType[maxncv];
	MatType* Z = new MatType[ldv * maxncv];
	MatType* v = new MatType[ldv * maxncv];
	MatType* workd = new MatType[3 * maxn];
	MatType* workev = new MatType[2 * maxncv];
	MatType* resid = new MatType[maxn];
	MatType* workl = new MatType[3 * maxncv * maxncv + 5 * maxncv];
	double rwork[maxn];
	char bmat[2] = "I";
	char which[3] = "LM";
	int ido, n, nx, nev, ncv, lworkl, info, ierr, ishfts, maxitr, mode1;
	MatType sigma = MatType(shift, 0.0);
	double tol;
	int rvec;

	nx = maxn;
	n = nx;
	nev = neigs;
	ncv = 3 * neigs + 1;
	lworkl  = 3 * ncv * ncv + 5 * ncv;
	tol    = 1E-8;
	ido    = 0;
	info   = 0;
	
	ishfts = 1;
	maxitr = 500;
	mode1 = 3;
	iparam[0] = ishfts;
	iparam[2] = maxitr;
	iparam[3] = 1;
	iparam[6] = mode1;

	while (true) {
		znaupd_(&ido, bmat, &n, which, &nev, &tol, (void*) resid, &ncv,
				(void*) v, &ldv, iparam, ipntr, (void*) workd, (void*) workl, &lworkl, rwork, &info);


		if (ido == -1 || ido == 1) {
			op.perform_op(&workd[ipntr[0] - 1], &workd[ipntr[1] - 1]);
		}
		else {
			break;
		}
	}

	if (info >= 0) {
		rvec = 1;
		char A = 'A';

		zneupd_(&rvec, &A, select, (void*) d, (void*) Z, &ldv, &sigma, (void*) workev, bmat, &n, which, &nev, 
				&tol, (void*) resid, &ncv, (void*) v, &ldv, iparam, ipntr, (void*) workd, (void*) workl, &lworkl,
				rwork, &ierr);

		for (int i = neigs - 1; i >= 0; --i) {
			evalues[neigs - i - 1] = d[i].real();
			for (int j = 0; j < n; ++j) {
				evecs[(neigs - i - 1) * n + j] = Z[i * n + j];
			}
		}
	}

	delete[] d;
	delete[] Z;
	delete[] v;
	delete[] workd;
	delete[] workev;
	delete[] resid;
	delete[] workl;
}

