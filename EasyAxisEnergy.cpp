#pragma once

#include "EasyAxisEnergy.h"

double OneEasyAxis(
	const vector<Particle>& moments,
	const vector<Particle>& particle_easy_axis, 
	const int& index)
{
	double nrg = kA * Vol * (1 - pow(DotP(moments[index], particle_easy_axis[index]),2)
		/ DotP(moments[index], moments[index]));
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
		easy_nrg[i] = OneEasyAxis(moments, easy_axis_dir, i);
		nrg += easy_nrg[i];
	}
	return nrg; 
}
