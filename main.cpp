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
#include <math.h>
#include "energy.h"

double D[2] = {1.0,1.0};			// {Dx, Dy}
double LAMBDA[2] = {1.0, 1.0};		// {Lx, Ly}
const int N = 1;					// monte-carlo iterations
double DT = 0.05;					// time increment
const double CL = 0.0;					// "temperature" of the system
const double PI  =3.141592653589793238463;

typedef std::array<std::array<double,S>,S> mat;		// matrix definition

void initializeMatrix(mat &theta, bool n) {
	//std::random_device rd;
	//std::mt19937 gen(rd());
	//std::uniform_real_distribution<> dis(0, 1);
	for(int i = 0; i < S; i++) {
		for (int j = 0; j < S; j++) {
			theta[i][j] = 0.5;
			//theta[i][j] = dis(gen);
			if (n) { theta[i][j] *= 2.0*PI*CL; }
		}
	}
}

// handle periodic boundary conditions
inline int pmod(int i, int n = S) {
	return (i % n + n) % n;
}

void calcPhase(mat &theta, double ener[]) {
	mat noise, current;
	for (int i=0; i<N; i++) {
		initializeMatrix(noise, true);
		for (int j = 0; j < S; j++) {
			for (int k = 0; k < S; k++) {
				double XPlus = theta[j][k] - theta[j][pmod(k+1)];
				double XMinus = theta[j][k] - theta[j][pmod(k-1)];
				double YPlus = theta[j][k] - theta[pmod(j+1)][k];
				double YMinus = theta[j][k] - theta[pmod(j-1)][k];

				double diffX = D[0]*(sin(XPlus) + sin(XMinus));
				double diffY = D[1]*(sin(YPlus) + sin(YMinus));
				double nonLinX = LAMBDA[0]*(cos(XPlus) + cos(XMinus));
				double nonLinY = LAMBDA[1]*(cos(YPlus) + cos(YMinus));

				current[j][k] = theta[j][k] - DT*(diffX + diffY + nonLinX + nonLinY + noise[j][k] - LAMBDA[0] - LAMBDA[1]);
				current[j][k] = fabs(fmod(current[j][k],2.0*PI));
			}
		}
		energy e; // create an object for my energy class
		ener[i] = e.calcEnergy(current, D);
		std::swap(current, theta);
	}
}

int main() {
	printf("Running compact KPZ: \n");
	double ener[N];
	mat theta;
	initializeMatrix(theta, false);

	printf("Number of iterations = %d \n", N);
	printf("Noise = 0.0 \n");
	printf("Initial theta values: \n");
	for (int i = 0; i < S; i++) {
		for (int j = 0; j < S; j++) {
			printf("%.5f ", theta[i][j]);
		}
		printf("\n");
	}

	calcPhase(theta, ener);

	printf("Final theta values: \n");
	for (int i = 0; i < S; i++) {
		for (int j = 0; j < S; j++) {
			printf("%.5f ", theta[i][j]);
		}
		printf("\n");
	}

	printf("Energy values: \n");
	for (int j = 0; j < N; j++) {
		printf("%.5f ", ener[j]);
	}
	printf("\n");

	// TODO: plot graph of time against phase (thetaF)

	return EXIT_SUCCESS;
}
