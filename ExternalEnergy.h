#pragma once

#include "mylibr.h"

double GetOneExtEnergy(
	const vector<Particle>& moments,
	const int& particle_index,
	const Particle& E_ext);

double ExtSysEnergy(
	const vector<Particle>& moments,
	const int& v_size,
	const Particle& E_ext);