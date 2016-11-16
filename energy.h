//
// Created by nisha on 03/11/2016.
//

#ifndef KPZ_TEST_ENERGY_H
#define KPZ_TEST_ENERGY_H
#endif //KPZ_TEST_ENERGY_H

#include<cstdlib>
#include <cmath>
#include <random>
#include <cstdio>
#include <array>

using namespace std;
const int S = 4; 		// size of lattice - currently code assumes square S x S matrix

typedef std::array<std::array<double,S>,S> mat;
class energy
{
public:
	double calcEnergy(mat &theta, double D[]);
};