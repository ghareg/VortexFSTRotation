#ifndef BASIS_H_
#define BASIS_H_
#include <stdint.h>
#include <complex>
#include "constants.h"

typedef int16_t Occup;
typedef int Count;
typedef std::complex<double> MatType;

class State
{
public:
	State(int ph, int s, int orb, int l, int j):
		mph(ph), ms(s), morb(orb), ml(l), mj(j) {}

	State(): mph(-1), ms(-1), morb(0), ml(-lmax), mj(0) {} 

	Occup ph() const {return mph;}
	Occup s() const {return ms;}
	Occup orb() const {return morb;}
	Occup l() const {return ml;}
	Occup j() const {return mj;}
	
private:
	Occup mph;
	Occup ms;
	Occup morb;
	Occup ml;
	Occup mj;
};

void GenerateBasis(State*& basis, Count& bSize);

#endif
