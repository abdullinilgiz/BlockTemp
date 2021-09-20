#pragma once

#include "ExternalEnergy.h"

double GetOneExtEnergy(
	const vector<Particle>& moments,
	const int& particle_index,
	const Particle& E_ext) 
{
	double energy = -DotP(E_ext, moments[particle_index]);

	return M0 * energy;
}

double ExtSysEnergy(
	const vector<Particle>& moments,
	const int& v_size,
	const Particle& E_ext)
{
	double energy = 0;
	for (int i = 0; i < v_size; ++i) {
		energy += GetOneExtEnergy(moments, i, E_ext);
	}

	return energy;
}

