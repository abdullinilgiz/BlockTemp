#pragma once

#include "EasyAxisEnergy.h"

double OneEasyAxis(
	const Particle& moment,
	const Particle& particle_easy_axis)
{
	double nrg = kA * Vol * (1 - pow(DotP(moment, particle_easy_axis),2)
		/ DotP(moment, moment));
	return nrg;
}

double EasySysEnergy(
	const vector<Particle>& moments, 
	vector<double>& easy_nrg,
	const vector<Particle>& easy_axis_dir,
	const int& nmb_of_particles) 
{
	double nrg = 0;
	for (int i = 0; i < nmb_of_particles; ++i) {
		easy_nrg[i] = OneEasyAxis(moments[i], easy_axis_dir[i]);
		nrg += easy_nrg[i];
	}
	return nrg; 
}
