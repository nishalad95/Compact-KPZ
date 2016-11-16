//
// Created by nisha on 03/11/2016.
//

#include<cstdlib>
#include <cmath>
#include <random>
#include <cstdio>
#include <array>
#include <math.h>
#include "energy.h"

using namespace std;
typedef std::array<std::array<double,S>,S> mat;

// handle periodic boundary conditions
inline int pmod(int i, int n = S) {
	return (i % n + n) % n;
}

double energy::calcEnergy(mat &theta, double D[]) {
	double energy = 0.0;
	for (int j = 0; j < S; j++) {
		for (int k = 0; k < S; k++) {
			energy += (-D[0]*(cos(theta[j][k] - theta[j][pmod(k+1)]) + cos(theta[j][k] - theta[j][pmod(k-1)]))
					   -D[1]*(cos(theta[j][k] - theta[pmod(j+1)][k]) + cos(theta[j][k] - theta[pmod(j-1)][k]))) / (S*S);
		}
	}
	return energy;
}

