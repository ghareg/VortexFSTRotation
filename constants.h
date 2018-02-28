#ifndef CONSTANTS_H_
#define CONSTANTS_H_


const double Pi = 3.14159;
const double E1 = 226.0; //in meV
const double E2 = 120.0; //in meV
const double E4 = 275.0; //in meV
const double m1x = -0.0271; //in meV^-1 x nm^-2
const double m2x = -0.0151; // in meV^-1 x nm^-2
const double m4x = 0.0131; // in meV^-1 x nm^-2
const double t1z = -4.0; // in meV
const double t2z = 76.0; // in meV
const double t4z = -426.0; // in meV
const double c = 0.59552; // in nm
const double beta = 15.24; // in meV x nm^2
const double Gamma = 0.3; // in meV x nm
const double delta = 244.8; // in meV x nm
const double lbd1 = 50.0; //in meV
const double lbd2 = 25.0; //in meV
const double lbd3 = 8.0; //in meV
const double Delta1 = 2.5; // in meV
const double Delta2 = 1.7; // in meV
const double R = 5000; // in nm
const double xi = 3; // in nm

//Auxiliary definitions
const double m1xm1 = 1.0f / (2.0f * m1x * R * R);
const double m2xm1 = 1.0f / (2.0f * m2x * R * R);
const double m4xm1 = 1.0f / (2.0f * m4x * R * R);
const double betaR = beta / (R * R);
const double deltaSR = delta * delta / (R * R);
const double deltaR = delta / R;
const double Rm2 = 1.0 / (R * R);
const double Rm1 = 1.0 / R;
const double sqrt2 = 1.0 / sqrt(2.0);


const int lmax = 10;
const int lCount = 2 * lmax + 1;
const int lACount = lmax + 1;
const int lSupMax = 10;
const int jsmax = 300;
const int jmax = 3 * jsmax;
const double thresh = 1E-6;
const int neigs = 40;
const int nCore = 10;
const double shift = 0.0;
#endif
