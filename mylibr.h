#pragma once

#ifndef HEADER_H
#define HEADER_H

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

double DotP(Particle l, Particle r);

const string Today = "27.09.21";

const int Nstep = 1600; // should be divisible by 5 and 4

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

const int Nside_X = 13; // number of particles per cube side
const double Xstep = Particle_diam;
const double XSideLength = Nside_X * Xstep;
const double XHalfSideLength = XSideLength / 2.0;

const double Ystep = Particle_radius * sqrt(3.0);
const int Nside_Y = round(XSideLength / (Ystep * (Nside_X - 1) + Particle_diam) * static_cast<int>(Nside_X));
const double YSideLength = (Nside_Y - 1) * Ystep + Particle_diam;
const double YHalfSideLength = YSideLength / 2.0;

const double Zstep = Particle_radius * 2 * sqrt(6) / 3;
const int Nside_Z = round(XSideLength / (Zstep * (Nside_X-1) + Particle_diam) * static_cast<int>(Nside_X));
const double ZSideLength = (Nside_Z - 1) * Zstep + Particle_diam;
const double ZHalfSideLength = ZSideLength / 2.0;

const double ShortestHalfSideLength = min(XSideLength, min(YSideLength, ZSideLength)) / 2.0;

const int Tmax = 400;
const int Tmin = 150;
const int DeltaT = 25;
//const double kBT = T * kB;

const int Nmb_particles = Nside_X * Nside_Y * Nside_Z;




#endif

