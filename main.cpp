/*
 * Nisha Lad 20/10/16
 * Supervised by Dr Marzena H. Syzmanska & Dr Alejandro Zamora
 * University College London - Master's Project
 * Monte Carlo simulation of the dynamics of the phase of Polaritons
 * described by the compact KPZ equation & periodic boundary conditions
 */

// original compact eq - 4 nearest neighbours

#define _USE_MATH_DEFINES

#include <math.h>
#include <string>
#include <sstream>
#include<cstdlib>
#include <cmath>
#include <random>
#include <cstdio>
#include <array>
#include <fstream>

using namespace std;

const int S = 64;
double Dx = 1.0;
double Dy = 1.0;
double Lx = 0.0;
double Ly = 0.0;
const int N = 100;					// monte-carlo iterations
const int R = 1;					// number of stochastic realisations
double dt = 0.05;					// time increment
const double CL = 0.0;				// "temperature" of the system

random_device rd;
mt19937 gen(rd());

typedef std::array<std::array<double,S>,S> mat;		// matrix definition

void initializeMatrix(mat &theta, bool n, uniform_real_distribution<> &dis,
					  mt19937 &gen, uniform_real_distribution<> &r_i) {
	for(int i = 0; i < S; i++) {
		for (int j = 0; j < S; j++) {
			theta[i][j] = dis(gen);
			if (n) { theta[i][j] *= CL*r_i(gen); }
		}
	}
}

// handle periodic boundary conditions
inline int pmod(int i, int n = S) {
	return (i % n + n) % n;
}

double calcEnergy(mat &theta) {
	double energy = 0.0;
	for (int j = 0; j < S; j++) {
		for (int k = 0; k < S; k++) {
			energy += -Dx*(cos(theta[j][k] - theta[j][pmod(k+1)]) + cos(theta[j][k] - theta[j][pmod(k-1)]))
					  -Dy*(cos(theta[j][k] - theta[pmod(j+1)][k]) + cos(theta[j][k] - theta[pmod(j-1)][k]));
		}
	}
	return energy;
}

double calcNumVortices(mat &theta) {
	double num = 0.0;
	for (int j = 0; j < S; j++) {
		for (int k = 0; k < S; k++) {
			num += fabs(theta[j][k] - theta[pmod(j-1)][k]
						+ theta[pmod(j-1)][pmod(k)] - theta[pmod(j-1)][pmod(k-1)]
						+ theta[pmod(j-1)][pmod(k-1)] - theta[pmod(j)][pmod(k-1)]
						+ theta[pmod(j)][pmod(k-1)] - theta[j][k]);
		}
	}
	return num;
}

void calcPhase(mat &theta, double *ener, double *vortexNum, uniform_real_distribution<> &dis,
			   mt19937 &gen, uniform_real_distribution<> &r_i) {
	mat noise, current;
	for (int i=0; i<N; i++) {
		initializeMatrix(noise, true, dis, gen, r_i);
		for (int j = 0; j < S; j++) {
			for (int k = 0; k < S; k++) {
				double XPlus = theta[j][k] - theta[j][pmod(k+1)];
				double XMinus = theta[j][k] - theta[j][pmod(k-1)];
				double YPlus = theta[j][k] - theta[pmod(j+1)][k];
				double YMinus = theta[j][k] - theta[pmod(j-1)][k];

				double diffX = Dx*(sin(XPlus) + sin(XMinus));
				double diffY = Dy*(sin(YPlus) + sin(YMinus));
				double nonLinX = Lx*(cos(XPlus) + cos(XMinus));
				double nonLinY = Ly*(cos(YPlus) + cos(YMinus));

				current[j][k] = theta[j][k] - dt*(diffX + diffY + nonLinX + nonLinY + noise[j][k] - Lx - Ly);
				current[j][k] = fmod(current[j][k],2.0*M_PI);
			}
		}
		ener[i] = calcEnergy(current) / (S*S*1.0);
		vortexNum[i] = calcNumVortices(current) / (2.0*M_PI);
		std::swap(current, theta);
	}
}

// for matrices
void outputMatToFile(mat &theta, string filename) {
	ofstream energyVals(filename); // opening output stream for file
	if (energyVals.is_open()) {
		for (int i=0; i<S; i++) {
			for (int j=0; j<S; j++) {
				energyVals << theta[i][j] << " "; //contiguous memory output
			}
		}
	} else { printf("File could not be opened.\n"); }
}

// for arrays
void outputToFile(double arr[], string filename) {
	ofstream energyVals(filename); // opening output stream for file
	if (energyVals.is_open()) {
		for (int i=0; i<N; i++) {
			energyVals << arr[i] << " ";
		}
	} else { printf("File could not be opened.\n"); }
}


void runKPZEquation() {
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dis(0, 2.0*M_PI);
	std::uniform_real_distribution<> r_i(-0.5, 0.5);

	for (int i = 1; i <= R; i++) {
		double ener[N];
		double vortexNum[N];
		mat theta;
		initializeMatrix(theta, false, dis, gen, r_i);        // initialize theta init
		//printf("Initial theta values: \n");
		//for (int i = 0; i < S; i++) {
		//	for (int j = 0; j < S; j++) {
		//		printf("%.5f ", theta[i][j]);
		//	}
		//	printf("\n");
		//}

		// Writing to file
		stringstream thetaInitFilename;
		thetaInitFilename << "CL_" << + CL << "/thetaInit_CL" << CL << "_N" << N << "_R" << i << ".txt";
		outputMatToFile(theta, thetaInitFilename.str());
		calcPhase(theta, ener, vortexNum, dis, gen, r_i);

//		printf("Final theta values: \n");
//		for (int i = 0; i < S; i++) {
//			for (int j = 0; j < S; j++) {
//				printf("%.5f ", theta[i][j]);
//			}
//			printf("\n");
//		}

//		printf("Energy values: \n");
//		for (int i = 0; i < N; i++) {
//			printf("%.5f ", ener[i] * S*S);
//		}
//		printf("\n");

//		printf("Energy density values: \n");
//		for (int i = 0; i < N; i++) {
//			printf("%.5f ", ener[i]);
//		}

		stringstream thetaFinFilename;
		thetaFinFilename << "CL_" << + CL << "/thetaFin_CL" << CL << "_N" << N << "_R" << i << ".txt";
		outputMatToFile(theta, thetaFinFilename.str());
		stringstream energyValFilename;
		energyValFilename << "CL_" << + CL << "/Energy_CL" << CL << "_N" << N << "_R" << i << ".txt";
		outputToFile(ener, energyValFilename.str());
		stringstream vortexNumFilename;
		vortexNumFilename << "CL_" << + CL <<  "/VortexNum_CL" << CL << "_N" << N << "_R" << i << ".txt";
		outputToFile(vortexNum, vortexNumFilename.str());
	}
}

int main() {
	printf("Running compact KPZ: \n");

	runKPZEquation();

	return EXIT_SUCCESS;

}