#pragma once

#include "mylibr.h"

void FillVectorWithCoords(vector <Particle>& particle_coords);

Particle FullSysMoment(
	const vector<Particle>& moments,
	const int& v_size);

Particle GetRandomDir(TRandomMersenne* rg);

void RotateVector(
	Particle& vector,
	const Particle& axis,
	const double& phi);