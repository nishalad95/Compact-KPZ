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

const int S = 128;					// lattice size
double Dx = 1.0;					// diffusion
double Dy = 1.0;
double Lx = -0.5;					// non-linearity
double Ly = 0.75;
const int N = 800;					// monte-carlo iterations
const int R = 50;					// stochastic realisations
double dt = 0.05;					// time increment
//double CL = 3.8;					// "temperature" of the system

random_device rd;
mt19937 gen(rd());					// twister to generate random values

typedef std::array<std::array<double,S>,S> mat;		// matrix definition

void initializeMatrix(mat &theta, mt19937 &gen, double CL, bool noise = false) {
	uniform_real_distribution<> dis(0,2*M_PI);
	uniform_real_distribution<> disNoise(-0.5,0.5);
	for(int i = 0; i < S; i++) {
		for (int j = 0; j < S; j++) {
			if (noise) { theta[i][j] = 2.0*M_PI*CL*disNoise(gen); }
			else { theta[i][j] = dis(gen); }
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
	double diff[4];

	for (int j = 0; j < S; j++) {
		for (int k = 0; k < S; k++) {

			diff[0] = theta[j][k] - theta[pmod(j-1)][k];
			diff[1] = theta[pmod(j-1)][pmod(k)] - theta[pmod(j-1)][pmod(k-1)];
			diff[2] = theta[pmod(j-1)][pmod(k-1)] - theta[pmod(j)][pmod(k-1)];
			diff[3] = theta[pmod(j)][pmod(k-1)] - theta[j][k];

			for (int l = 0; l < 4; l++) {
				if (diff[l] > M_PI) {
					diff[l] -= 2*M_PI;
				} else if (diff[l] < -M_PI) {
					diff[l] += 2*M_PI;
				}
			}
			num += fabs(diff[0] + diff[1] + diff[2] + diff[3]);
		}
	}
	return num;
}

void calcPhase(mat &theta, double *ener, double *vortexNum, mt19937 &gen, double CL) {
	mat noise, current;
	for (int i=0; i<N; i++) {
		initializeMatrix(noise, gen, CL, true);
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
				if (current[j][k] < 0) {
					current[j][k] += 2.0*M_PI;
				}
			}
		}
		ener[i] = calcEnergy(current) / (S*S*1.0);
		vortexNum[i] = calcNumVortices(current) / (2.0*M_PI);
		std::swap(current, theta);
	}
}

// for arrays
void outputToFile(double arr[], string filename) {
	ofstream energyVals(filename); // opening output stream for file
	if (energyVals.is_open()) {
		for (int i=0; i<N; i++) {
			energyVals << arr[i] << " ";
		}
	} else { printf("File could not be opened or not found.\n"); }
}


void runKPZEquation(double CL) {
	for (int i = 1; i <= R; i++) {
		double ener[N];
		double vortexNum[N];
		mat theta;
		initializeMatrix(theta, gen, CL);

		calcPhase(theta, ener, vortexNum, gen, CL);

		// Writing to file
		stringstream energyValFilename;
		energyValFilename << "data/S=128/Lx=" <<Lx<< ",Ly=" <<Ly<< "/CL_"
						  << CL << "/Energy_CL" << CL << "_N" << N << "_R" << i << ".txt";
		outputToFile(ener, energyValFilename.str());
		stringstream vortexNumFilename;
		vortexNumFilename << "data/S=128/Lx=" <<Lx<< ",Ly=" <<Ly<< "/CL_"
						  << CL << "/VortexNum_CL" << CL << "_N" << N << "_R" << i << ".txt";
		outputToFile(vortexNum, vortexNumFilename.str());
	}
}

int main() {
	printf("Running compact KPZ: \n");
	double CL[11] = {3.1, 3.15, 3.2, 3.25, 3.3, 3.35, 3.4, 3.45, 3.55, 3.6, 3.65};
	for (int a=0; a<11; a++) {
		runKPZEquation(CL[a]);
		printf("Completed simulation for %f", CL[a]);
		printf("\n");
	}

	return EXIT_SUCCESS;
}