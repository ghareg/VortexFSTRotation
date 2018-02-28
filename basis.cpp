#include "basis.h"
#include "matrix.h"

void GenerateBasis(State*& basis, Count& bSize)
{
	bSize = 0;
	Count maxbSize = lCount * jmax * 4 * 2 * 2;
	basis = new State[maxbSize];
	int ls = 0;
	int le = 0;
	for (int ph = -1; ph <= 1; ph += 2) {
		ls = (ph == -1) ? -lmax : -lmax + 1;
		le = (ph == -1) ? lmax - 1 : lmax;
		for (int s = -1; s <= 1; s += 2) {
			for (int orb = 0; orb < 4; ++orb) {
				for (int l = ls; l <= le; l += 1) {
					for (int j = 0; j < jmax; ++j) {
						if (bSize == maxbSize) {
							maxbSize = 3 * maxbSize / 2;
							resize(basis, bSize, maxbSize);
						}
						basis[bSize++] = State(ph, s, orb, l, j);
					}
				}
			}
		}
	}
}
