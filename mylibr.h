#pragma once

#ifndef HEADER_H
#define HEADER_H

#include <stdlib.h>
#include <string.h>
#include <malloc.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <exception>
#include <string>
#include <iomanip>

#include "randomc.h"

using namespace std;

struct Particle {
	double x;
	double y;
	double z;
};

double DotP(Particle l, Particle r);

const string Today = "27.09.21_old";

const int Nstep = 1600; // should be divisible by 5 and 4

const double PI = 3.141592653589793;

const double kB = 1.3806503e-23;

const double Rcore = 6e-9; // radius

const double Tshell = 5e-9; // tolshina obolochki

const int Nside = 13; // number of particles per cube side

const int iSeed = 13; // for randomc

const double dMoment = 3e-19; //dipol moment

const double M0 = 1.2566370614e-6; 

// const Particle E_ext = { 0, 0, 8000 };
const double kA = 30000;
const double Vol = 4.0 / 3.0 * PI * static_cast<double>(pow(Rcore, 3));

const double SideLength = Nside * (2 * (Tshell + Rcore));
const double HalfSideLength = SideLength / 2;

const int Tmax = 400;
const int Tmin = 150;
const int DeltaT = 25;
//const double kBT = T * kB;

const double Particle_diam = 2 * (Tshell + Rcore); // find particle diameter
const double Particle_radius = Particle_diam / 2;

const int Nmb_particles = Nside * Nside * Nside;




#endif

