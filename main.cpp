/*
 * Nisha Lad 20/10/16
 * Supervised by Dr Marzena H. Syzmanska & Dr Alejandro Zamora
 * University College London - Master's Project
 * Monte Carlo simulation of the dynamics of the phase of Polaritons
 * described by the compact KPZ equation & periodic boundary conditions
 */

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

const int S = 64;					// lattice size
double Dx = 1.0;					// diffusion
double Dy = 1.0;
double Lx = 0.0;					// non-linearity
double Ly = 0.0;
double dt = 0.05;					// time increment
const double tolerance = 1e-5;		//convergence threshold

random_device rd;
mt19937 gen(rd());					// twister to generate random values

typedef std::array<std::array<double,S>,S> mat;		// matrix definition
typedef std::vector<double> vec;

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

int calcPhase(mat &theta, vec &ener, vec &vortexNum, mt19937 &gen, double CL) {
	mat noise, current;
	double diff = 100.0;
	int i = 0;
	while (diff > tolerance) {
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
				if (current[j][k] < 0) { current[j][k] += 2.0*M_PI; }
			}
		}
		double energyDensity = calcEnergy(current) / (S*S*1.0);
		double numVortices = calcNumVortices(current) / (2.0*M_PI);
		ener.push_back(energyDensity);
		vortexNum.push_back(numVortices);
		if (i != 0) { diff = fabs(ener[i] - ener[i-1]); }
		i++;
		std::swap(current, theta);
	}
	return i;
}

void outputToFile(int count, vec &array, string filename) {
	ofstream energyVals(filename); // opening output stream for file
	if (energyVals.is_open()) {
		for (int i=0; i<count; i++) {
			energyVals << array[i] << " ";
		}
	} else { printf("File could not be opened or not found.\n"); }
}


void runKPZEquation(double CL) {
	int R;
	if (CL == 0.0) { R = 1; }
	else if (CL <= 2.5) { R = 10; }
	else if (CL <= 3.0) { R = 30; }
	else if (CL <= 4.0) { R = 50; }
	else if (CL <= 5.5) { R = 80; }
	else { R = 100; }
	for (int i = 1; i <= R; i++) {
		vec ener;
		vec vortexNum;
		mat theta;
		initializeMatrix(theta, gen, CL);
		int count = calcPhase(theta, ener, vortexNum, gen, CL);
		// Writing to file
		stringstream energyValFilename;
		energyValFilename << "data/S=" << S << "/Lx=" <<Lx<< ",Ly=" <<Ly<< "/CL_"
						  << CL << "/Energy_CL" << CL << "_R" << i << ".txt";
		outputToFile(count, ener, energyValFilename.str());
		stringstream vortexNumFilename;
		vortexNumFilename << "data/S=" << S << "/Lx=" <<Lx<< ",Ly=" <<Ly<< "/CL_"
						  << CL << "/VortexNum_CL" << CL << "_R" << i << ".txt";
		outputToFile(count, vortexNum, vortexNumFilename.str());
	}
}

int main() {
	printf("Running compact KPZ: \n");
	double CL[30] = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 2.95, 3.0, 3.05, 3.1, 3.15, 3.2, 3.25,
					 3.3, 3.35, 3.4, 3.45, 3.5, 3.55, 3.6, 3.65, 3.7, 3.75, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0};
	for (int a=0; a<30; a++) {
		runKPZEquation(CL[a]);
		printf("Completed simulation for %f", CL[a]);
		printf("\n");
	}
	return EXIT_SUCCESS;
}