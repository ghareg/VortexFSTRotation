#include "basis.h"
#include "matrix.h"

void GenerateBasis(State*& basis, Count& bSize)
{
	bSize = 0;
	Count maxbSize = lCount * jmax * 4 * 2 * 2;
	basis = new State[maxbSize];
	int ls = 0;
	int le = 0;
	//int les[16] = {-10, -9, -7, -8, -7, -10, -8, -9, -9, -8, -6, -7, -6, -9, -7, -8}; // Eig C4 1
	//int lee[16] = {6, 7, 9, 8, 9, 6, 8, 7, 7, 8, 10, 9, 10, 7, 9, 8}; // Eig C4 1
	//int les[16] = {-9, -8, -10, -7, -10, -9, -7, -8, -8, -7, -9, -6, -9, -8, -6, -7}; // Eig C4 i
	//int lee[16] = {7, 8, 6, 9, 6, 7, 9, 8, 8, 9, 7, 10, 7, 8, 10, 9}; // Eig C4 i
	//int les[16] = {-8, -7, -9, -10, -9, -8, -10, -7, -7, -6, -8, -9, -8, -7, -9, -6}; // Eig C4 -1
	//int lee[16] = {8, 9, 7, 6, 7, 8, 6, 9, 9, 10, 8, 7, 8, 9, 7, 10}; // Eig C4 -1
	int les[16] = {-7, -10, -8, -9, -8, -7, -9, -10, -6, -9, -7, -8, -7, -6, -8, -9}; // Eig C4 -i
	int lee[16] = {9, 6, 8, 7, 8, 9, 7, 6, 10, 7, 9, 8, 9, 10, 8, 7}; // Eig C4 -i

	for (int ph = -1; ph <= 1; ph += 2) {
		ls = (ph == -1) ? -lmax : -lmax + 1;
		le = (ph == -1) ? lmax - 1 : lmax;
		for (int s = -1; s <= 1; s += 2) {
			for (int orb = 0; orb < 4; ++orb) {
				//for (int l = ls; l <= le; l += 1) {
				for (int l = les[(ph + 1) * 4 + (s + 1) * 2 + orb]; l <= lee[(ph + 1) * 4 + (s + 1) * 2 + orb]; l += 4) {
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
