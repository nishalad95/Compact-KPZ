/*
 * Nisha Lad 20/10/16
 * Supervised by Dr Marzena H. Syzmanska & Dr Alejandro Zamora
 * University College London - Master's Project
 * Monte Carlo simulation of the dynamics of the phase of Polaritons
 * described by the compact KPZ equation & periodic boundary conditions
 */

// original compact eq - 4 nearest neighbours

#include <cstdlib>
#include <cmath>
#include <random>
#include <cstdio>
#include <array>

double D[2] = {1.0,1.0};			// {Dx, Dy}
double LAMBDA[2] = {1.0, 1.0};		// {Lx, Ly}
const int N = 1;					// monte-carlo iterations
const int S = 4;					// size of lattice - currently code assumes square S x S matrix
double DT = 0.5;					// time increment

typedef std::array<std::array<double,S>,S> mat;		// matrix definition

void initializeTheta(mat &theta) {
	//std::random_device rd;
	//std::mt19937 gen(rd());
	//std::uniform_real_distribution<> dis(0, 1);
	for(int i = 0; i < S; i++) {
		for (int j = 0; j < S; j++) {
			//theta[i][j] = dis(gen);
			theta[i][j] = 0.5;
		}
	}
}

void initializeNoise(mat &noise) {
	//std::random_device rd;
	//std::mt19937 gen(rd());
	//std::uniform_real_distribution<> dis(0, 1);
	for(int i = 0; i < S; i++) {
		for(int j = 0; j < S; j++) {
			//noise[i][j] = dis(gen);
			noise[i][j] = 0.0;
		}
	}
}

// handle periodic boundary conditions
inline int pmod(int i, int n = S) {
	return (i % n + n) % n;
}

void calcPhase(mat &theta) {
	mat noise, current;
	for (int i=0; i<N; i++) {
		initializeNoise(noise);
		for (int j = 0; j < S; j++) {
			for (int k = 0; k < S; k++) {
				double diffX = -D[0]*(sin(theta[j][k] - theta[j][pmod(k+1)])
									  + sin(theta[j][k] - theta[j][pmod(k-1)]));
				double diffY = -D[1]*(sin(theta[j][k] - theta[pmod(j+1)][k])
									  + sin(theta[j][k] - theta[pmod(j-1)][k]));
				double nonLinX = -LAMBDA[0]*0.5*(cos(theta[j][k] - theta[j][pmod(k+1)])
												 + cos(theta[j][k] - theta[j][pmod(k-1)]));
				double nonLinY = -LAMBDA[1]*0.5*(cos(theta[j][k] - theta[pmod(j+1)][k])
												 + cos(theta[j][k] - theta[pmod(j-1)][k]));

				current[j][k] = DT*(theta[j][k] + diffX + diffY + nonLinX + nonLinY + noise[j][k]);
			}
		}
		std::swap(current, theta);
	}
}

int main() {
	printf("Running compact KPZ: \n");
	mat theta;
	initializeTheta(theta);

	printf("Initial theta values: \n");
	for (int i = 0; i < S; i++) {
		for (int j = 0; j < S; j++) {
			printf("%.5f ", theta[i][j]);
		}
		printf("\n");
	}

	calcPhase(theta);

	printf("Final theta values: \n");
	for (int i = 0; i < S; i++) {
		for (int j = 0; j < S; j++) {
			printf("%.5f ", theta[i][j]);
		}
		printf("\n");
	}

	// TODO: plot graph of time against phase (thetaF)

	return EXIT_SUCCESS;
}
