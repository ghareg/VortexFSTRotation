#include "operator.h"
#include "constants.h"
#include <cmath>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_integration.h>

MatType MatElement(Count i, Count j, const State* basis, const InitVal& inVal, double mu, double rot)
{
	Occup ph = basis[i].ph();
	Occup s = basis[i].s();
	Occup orb = basis[i].orb();
	Occup l = basis[i].l();
	Occup jc = basis[i].j();
	Occup php = basis[j].ph();
	Occup sp = basis[j].s();
	Occup orbp = basis[j].orb();
	Occup lp = basis[j].l();
	Occup jcp = basis[j].j();
	MatType result = MatType(0.0, 0.0);
	if (ph == php) {
		if (s == sp) {
			if (l == lp && jc == jcp) {
					result += MatType(-ph * EnDiag(orb, orbp, l, jc, s, mu, inVal), 0.0);
					result += -ph * 1.0 * RotationDiag(orb, orbp, l, jc, rot, inVal);
			}
			result += MatType(-ph * EnNonDiag(orb, orbp, l, lp, jc, jcp, inVal), 0.0);
			result += -ph * 1.0 * RotationNonDiag(orb, orbp, l, lp, jc, jcp, rot, inVal);
		}
		else {
			if (l == lp && jc == jcp) {
				result += MatType(0.0, -ph * SocNonDiag(orb, orbp, s));
			}
			result += -ph * 1.0 * RotationSocNonDiag(s, orb, orbp, l, lp, jc, jcp, rot, inVal);
		}
	}
	else {
		if (s == sp) {
			result += MatType(Super(ph, orb, orbp, l, lp, jc, jcp, inVal), 0.0);
		}
	}

	return result;	
}


double EnDiag(int orb, int orbp, int l, int j, int s, double mu, const InitVal& inVal)
{
	double EnC = 0.0;
	if (orb == orbp) {
		int la = std::abs(l);
		double alphap = inVal.BesselZeros[la * jmax + j]; 
		if (orb == 0) {
			EnC +=  E1 + alphap * alphap * m1xm1 - mu;
		}
		else if (orb == 1) {
			EnC +=  E2 + alphap * alphap * m2xm1 - mu;
			EnC += s * lbd1; 
		}
		else if (orb == 2) {
			EnC +=  E2 + alphap * alphap * m2xm1 - mu;
			EnC += -s * lbd1; 
		}
		else if (orb == 3) {
			EnC +=  E4 + alphap * alphap * m4xm1 - mu;
		}
	}
	
	return EnC;
}

double EnNonDiag(int orb, int orbp, int l, int lp, int j, int jp, const InitVal& inVal)
{
	double EnC = 0.0;
	if (lp == l + 1) {
		if (orb == 1 && orbp == 3) {
			int la = std::abs(l);
			int lpa = std::abs(lp);
			double alphap = inVal.BesselZeros[lpa * jmax + jp];
			double alpha = inVal.BesselZeros[la * jmax + j];
			EnC += deltaR * alphap * inVal.normBes[la * jmax + j] * inVal.normBes[lpa * jmax + jp] * FR(l, lp, j, jp, alpha, alphap);
		}
		else if (orb == 3 && orbp == 2) {
			int la = std::abs(l);
			int lpa = std::abs(lp);
			double alphap = inVal.BesselZeros[lpa * jmax + jp];
			double alpha = inVal.BesselZeros[la * jmax + j];
			EnC += -deltaR * alphap * inVal.normBes[la * jmax + j] * inVal.normBes[lpa * jmax + jp] * FR(l, lp, j, jp, alpha, alphap);
		}
	}
	else if (l == lp + 1) {
		if (orb == 3 && orbp == 1) {
			int la = std::abs(l);
			int lpa = std::abs(lp);
			double alphap = inVal.BesselZeros[lpa * jmax + jp];
			double alpha = inVal.BesselZeros[la * jmax + j];
			EnC += deltaR * alphap * inVal.normBes[la * jmax + j] * inVal.normBes[lpa * jmax + jp] * FR(l, lp, j, jp, alpha, alphap);
		}
		else if (orb == 2 && orbp == 3) {
			int la = std::abs(l);
			int lpa = std::abs(lp);
			double alphap = inVal.BesselZeros[lpa * jmax + jp];
			double alpha = inVal.BesselZeros[la * jmax + j];
			EnC += -deltaR * alphap * inVal.normBes[la * jmax + j] * inVal.normBes[lpa * jmax + jp] * FR(l, lp, j, jp, alpha, alphap);
		}
	}
	else if (lp == l + 2) {
		if (orb == 1 && orbp == 2) {
			int la = std::abs(l);
			int lpa = std::abs(lp);
			double alphap = inVal.BesselZeros[lpa * jmax + jp];
			double alpha = inVal.BesselZeros[la * jmax + j];
			EnC += betaR * alphap * alphap * inVal.normBes[la * jmax + j] * inVal.normBes[lpa * jmax + jp] * FR(l, lp, j, jp, alpha, alphap);
		}
	}
	else if (l == lp + 2) {
		if (orb == 2 && orbp == 1) {
			int la = std::abs(l);
			int lpa = std::abs(lp);
			double alphap = inVal.BesselZeros[lpa * jmax + jp];
			double alpha = inVal.BesselZeros[la * jmax + j];
			EnC += betaR * alphap * alphap * inVal.normBes[la * jmax + j] * inVal.normBes[lpa * jmax + jp] * FR(l, lp, j, jp, alpha, alphap);
		}
	}

	return EnC;
}

double SocNonDiag(int orb, int orbp, int s)
{
	if (s == -1) {
		if (orb == 0 && orbp == 1) {
			return 	sqrt(2.0) * lbd2;
		}
		else if (orb == 2 && orbp == 0) {
			return -sqrt(2.0) * lbd2;
		}
	}
	else {
		if (orb == 1 && orbp == 0) {
			return -sqrt(2.0) * lbd2;
		}
		else if (orb == 0 && orbp == 2) {
			return sqrt(2.0) * lbd2;
		}
	}

	return 0.0;
}

double Super(int ph, int orb, int orbp, int l, int lp, int jc, int jcp, const InitVal& inVal)
{
	if (orb == orbp) {
		if ((ph == -1 && lp == l + 1) || (ph == 1 && l == lp + 1)) {
			double DeltaL = (orb == 0 || orb == 3) ? Delta1 : Delta2;
			int la = std::abs(l);
			int lpa = std::abs(lp);
			double sInt = 0.0;
			if (lp > l) {
				if (l < 0) {
					sInt = -inVal.supInt[3 * lpa * jmax * jmax + 3 * jcp * jmax + 3 * jc + 1];
				}
				else {
					sInt = inVal.supInt[3 * la * jmax * jmax + 3 * jc * jmax + 3 * jcp + 1];
				}
			}
			else {
				if (lp < 0) {
					sInt = -inVal.supInt[3 * la * jmax * jmax + 3 * jc * jmax + 3 * jcp + 1];
				}
				else {
					sInt = inVal.supInt[3 * lpa * jmax * jmax + 3 * jcp * jmax + 3 * jc + 1];
				}
			}

			return DeltaL * inVal.normBes[la * jmax + jc] * inVal.normBes[lpa * jmax + jcp] * sInt;
		}
	}

	return 0.0;
}

MatType RotationDiag(int orb, int orbp, int l, int jc, double rot, const InitVal& inVal)
{
	if (orb == orbp) {
		int la = std::abs(l);
		double alphap = inVal.BesselZeros[la * jmax + j];
		if (orb == 0) {
			return MatType(0.25 * t1z * c * c * rot * rot * alphap * alphap * Rm2, 0.0);
		}
		else if (orb == 1 || orb == 2) {
			return MatType(0.25 * t2z * c * c * rot * rot * alphap * alphap * Rm2, 0.0);
		}
		else {
			return MatType(0.25 * t4z * c * c * rot * rot * alphap * alphap * Rm2, 0.0);
		}
	}
	else if ((orb == 0 && orbp == 1) || (orb == 2 && orbp == 0)) {
		int la = std::abs(l);
		double alphap = inVal.BesselZeros[la * jmax + j];
		return MatType(0.0, 0.5 * Gamma * c * rot * alphap * alphap * Rm2);
	}
	else if ((orb == 0 && orbp == 2) || (orb == 1 && orbp == 0)) {
		int la = std::abs(l);
		double alphap = inVal.BesselZeros[la * jmax + j];
		return MatType(0.0, -0.5 * Gamma * c * rot * alphap * alphap * Rm2);
	}

	return MatType(0.0, 0.0);
}

MatType RotationNonDiag(int orb, int orbp, int l, int lp, int jc, int jcp, double rot, const InitVal& inVal)
{
	if (orb == orbp) {
		if (l == lp + 2 || lp == l + 2) {
			int la = std::abs(l);
			int lpa = std::abs(lp);
			double alpha = inVal.BesselZeros[la * jmax + j];
			double alphap = inVal.BesselZeros[lpa * jmax + j];
			double normB = inVal.normBes[la * jmax + jc] * inVal.normBes[lpa * jmax + jcp]; 
			if (orb == 0) {
				return MatType(0.125 * t1z * c * c * rot * rot * alphap * alphap * Rm2 * normB * FR(l, lp, jc, jcp, alpha, alphap), 0.0);
			}
			else if (orb == 1 || orb == 2) {
				return MatType(0.125 * t2z * c * c * rot * rot * alphap * alphap * Rm2 * normB * FR(l, lp, jc, jcp, alpha, alphap), 0.0);
			}
			else {
				return MatType(0.125 * t4z * c * c * rot * rot * alphap * alphap * Rm2 * normB * FR(l, lp, jc, jcp, alpha, alphap), 0.0);
			}
		}
	}
	else if ((orb == 0 && orbp == 1) || (orb == 2 && orbp == 0)) {
		if (lp == l + 2) {
			int la = std::abs(l);
			int lpa = std::abs(lp);
			double alpha = inVal.BesselZeros[la * jmax + j];
			double alphap = inVal.BesselZeros[lpa * jmax + j];
			double normB = inVal.normBes[la * jmax + jc] * inVal.normBes[lpa * jmax + jcp];
			return MatType(0.0, 0.5 * Gamma * c * rot * alphap * alphap * Rm2 * normB * FR(l, lp, jc, jcp, alpha, alphap));
		}
	}
	else if ((orb == 0 && orbp == 2) || (orb == 1 && orbp == 0)) {
		if (l == lp + 2) {
			int la = std::abs(l);
			int lpa = std::abs(lp);
			double alpha = inVal.BesselZeros[la * jmax + j];
			double alphap = inVal.BesselZeros[lpa * jmax + j];
			double normB = inVal.normBes[la * jmax + jc] * inVal.normBes[lpa * jmax + jcp];
			return MatType(0.0, -0.5 * Gamma * c * rot * alphap * alphap * Rm2 * normB * FR(l, lp, jc, jcp, alpha, alphap));
		}
	}

	return MatType(0.0, 0.0);
}

MatType RotationSocNonDiag(int s, int orb, int orbp, int l, int lp, int jc, int jcp, double rot, const InitVal& inVal)
{
	if (l == lp + 1 || lp == l + 1) {
		if (s == -1) {
			if ((orb == 1 && orbp == 3) || (orb == 3 && orbp == 2)) {
				int la = std::abs(l);
				int lpa = std::abs(lp);
				double alpha = inVal.BesselZeros[la * jmax + j];
				double alphap = inVal.BesselZeros[lpa * jmax + j];
				double normB = inVal.normBes[la * jmax + jc] * inVal.normBes[lpa * jmax + jcp];
				return MatType(-c * lbd3 * rot * sqrt2 * alphap * Rm1 * normB * FR(l, lp, jc, jcp, alpha, alphap), 0.0);
			}
		}
		else {
			if ((orb == 2 && orbp == 3) || (orb == 3 && orbp == 1)) {
				int la = std::abs(l);
				int lpa = std::abs(lp);
				double alpha = inVal.BesselZeros[la * jmax + j];
				double alphap = inVal.BesselZeros[lpa * jmax + j];
				double normB = inVal.normBes[la * jmax + jc] * inVal.normBes[lpa * jmax + jcp];
				return MatType(-c * lbd3 * rot * sqrt2 * alphap * Rm1 * normB * FR(l, lp, jc, jcp, alpha, alphap), 0.0);
			}
		}
	}

	return MatType(0.0, 0.0);
}

double IntB (double r, void *paramsIntB)
{
	intBparams *intBprm = (intBparams *)paramsIntB;
	int l = (intBprm -> l);
	int lp = (intBprm -> lp);
	double alpha = (intBprm -> alpha);
	double alphap = (intBprm -> alphap);
	
	return r * std::tanh(r * R / xi) * gsl_sf_bessel_Jn(l, alpha * r) * gsl_sf_bessel_Jn(lp, alphap * r);
}	

// Superconductivity integral with Bessel
double IntTanSlow(int l, int lp, double alpha, double alphap, int& err)
{
	gsl_integration_cquad_workspace* w = gsl_integration_cquad_workspace_alloc (10000);
	double min = 0;
	double max = 1;
	double eabs = 1.0e-7;
	double erel = 1.0e-7;
	double result = 0.0;
	double error;

	intBparams intBprm = {l, lp, alpha, alphap};

	gsl_function F;
	
	F.function = &IntB;
	F.params = &intBprm;

	err = gsl_integration_cquad(&F, min, max, eabs, erel, w, &result, &error, NULL);
	gsl_integration_cquad_workspace_free(w);
	return result;
}

double IntTanFast(int l, int lp, double alpha, double alphap, int& err)
{
	gsl_integration_workspace* w = gsl_integration_workspace_alloc (10000);
	double min = 0;
	double max = 1;
	double eabs = 1.0e-7;
	double erel = 1.0e-7;
	double result = 0.0;
	double error;

	intBparams intBprm = {l, lp, alpha, alphap};

	gsl_function F;
	
	F.function = &IntB;
	F.params = &intBprm;

	err = gsl_integration_qag(&F, min, max, eabs, erel, 5000, 6, w, &result, &error);
	gsl_integration_workspace_free(w);
	return result;
}

double IntTan(int l, int lp, double alpha, double alphap)
{
	double result = 0.0;
	int err = 0;
	result = IntTanFast(l, lp, alpha, alphap, err);
	if (err) {
		result = IntTanSlow(l, lp, alpha, alphap, err);
	}

	return result;
}

double FR(int l, int lp, int j, int jp, double alpha, double alphap)
{
	if (abs(l) == abs(lp)) {
		if (j == jp) {
			double BesselJ = gsl_sf_bessel_Jn(l + 1, alpha);
			return 0.5 * BesselJ * BesselJ;
		}

		return 0.0;
	}

	return alpha * gsl_sf_bessel_Jn(l, alphap) * gsl_sf_bessel_Jn(l + 1, alpha) / (alpha * alpha - alphap * alphap);
}
