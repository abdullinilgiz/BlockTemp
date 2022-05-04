#pragma once

#include <cstdlib>
#include <string>
#include <malloc.h>
#include <cstdio>
#include <cmath>
#include <ctime>
#include <vector>
#include <iostream>
#include <fstream>
#include <exception>
#include <string>
#include <iomanip>

#include "randomc.h"
#include "profile.h"

using namespace std;

struct Particle {
    double x;
    double y;
    double z;
};

const string Today = "27.10.21";

const int Nstep = 5000; // should be divisible by 5 and 4

const double PI = 3.141592653589793;

const double kB = 1.3806503e-23;

const double Rcore = 6e-9; // radius

const double Tshell = 5e-9; // tolshina obolochki `

const int iSeed = 13; // for randomc

const double dMoment = 3e-19; //dipol moment

const double M0 = 1.2566370614e-6;

// const Particle E_ext = { 0, 0, 8000 };
const double kA = 30000;

const double Vol = 4.0 / 3.0 * PI * pow(Rcore, 3);

const double Particle_radius = (Tshell + Rcore);

const double Particle_diam = 2 * Particle_radius; // find particle diameter

extern double XSideLength;
extern double XHalfSideLength;
extern double YSideLength;
extern double YHalfSideLength;
extern double ZSideLength;
extern double ZHalfSideLength;

extern double ShortestHalfSideLength;

const int Tmax = 350;
const int Tmin = 150;
const int DeltaT = 50;
//const double kBT = T * kB;

const int Nmb_particles = 11 * 11 * 11;

