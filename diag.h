#ifndef DIAG_H_
#define DIAG_H_
#include "matrix.h"
#include <mkl_types.h>
#undef MKL_Complex16
#define MKL_Complex16 MatType
#include <mkl.h>
#include "slu_zdefs.h"

class GenMatProd
{
public:
	GenMatProd();
	~GenMatProd();

	void init(SparseMat baseHam);
	void restart(SparseMat baseHam);

	Count rows() {return qsize_;}
	Count cols() {return qsize_;}
	void set_shift(double sigma) {sigma_ = sigma;}
	void perform_op(const MatType* xIn, MatType* yOut);
private:
	SparseMat baseHam_;
	Count qsize_;
	double sigma_;
	_MKL_DSS_HANDLE_t pt_;
	MKL_INT mtype_;
	MKL_INT* iparm_;
	MKL_INT* perm_;
	MKL_INT phase_;
};

void calcEValues(const SparseMat& baseHam, GenMatProd& op, double* evalues, MatType* evecs);
#endif
