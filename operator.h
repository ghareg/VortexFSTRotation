#ifndef OPERATOR_H_
#define OPERATOR_H_
#include "matrix.h"

MatType MatElement(Count i, Count j, const State* basis, const InitVal& inVal, double mu, double rot);
double EnDiag(int orb, int orbp, int l, int j, int s, double mu, const InitVal& inVal);
double EnNonDiag(int orb, int orbp, int l, int lp, int j, int jp, const InitVal& inVal);
double SocNonDiag(int orb, int orbp, int s);
double Super(int ph, int orb, int orbp, int l, int lp, int jc, int jcp, const InitVal& inVal);

MatType RotationDiag(int orb, int orbp, int l, int jc, double rot, const InitVal& inVal);
MatType RotationNonDiag(int orb, int orbp, int l, int lp, int jc, int jcp, double rot, const InitVal& inVal);
MatType RotationSocNonDiag(int s, int orb, int orbp, int l, int lp, int jc, int jcp, double rot, const InitVal& inVal);

struct intBparams {int l; int lp; double alpha; double alphap;};
double IntB (double r, void *paramsIntB);
double IntTanSlow(int l, int lp, double alpha, double alphap, int& err);
double IntTanFast(int l, int lp, double alpha, double alphap, int& err);
double IntTan(int l, int lp, double alpha, double alphap);
double FR(int l, int lp, int j, int jp, double alpha, double alphap);

#endif
