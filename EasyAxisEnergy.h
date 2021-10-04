#pragma once

#include "mylibr.h"

double OneEasyAxis(
	const Particle& moments,
	const Particle& particle_easy_axis);

double EasySysEnergy(
	const vector<Particle>& moments,
	vector<double>& easy_nrg,
	const vector<Particle>& easy_axis_dir,
	const int& nmb_of_particles);