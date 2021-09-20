#pragma once

#include "mylibr.h"
#include "EasyAxisEnergy.h"
#include "ExternalEnergy.h"
#include "InteractionEnergy.h"
#include "Moments.h"
#include "operations.h"
#include "randomc.h"

bool MKIteration(
	const int& Temper,
	TRandomMersenne* rg,
	vector<double>& particles_easy_nrg,
	vector<vector<double>>& int_nrg_matrix,
	const vector<vector<Particle>>& v_distances,
	const vector<vector<double>>& distances,
	vector <Particle>& particle_moments,
	double& int_sys_nrg,
	double& ext_sys_nrg,
	double& easy_sys_nrg,
	double& full_sys_nrg,
	const double& phi_max,
	const vector <Particle>& easy_axis_dir,
	const Particle& E_ext);

void AngleChanger(const int& apply_counter, double& phi_max);