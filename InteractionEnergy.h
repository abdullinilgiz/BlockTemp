#pragma once

#include "mylibr.h"

double CalculateDipol_Dipol(
	const vector<Particle>& moments,
	const vector<vector<double>>& distances_matrix,
	const vector<vector<Particle>>& v_distances_matrix,
	const int& index_1,
	const int& index_2);

double IntParticleEnergy(
	const vector<Particle>& moments,
	vector<vector<double>>& energy_matrix,
	const vector<vector<double>>& distances_matrix,
	const vector<vector<Particle>>& v_distances_matrix,
	const int& particle_index);

double IntSysEnergy(
	const vector<Particle>& moments,
	vector<vector<double>>& matrix,
	const vector<vector<double>>& distances_matrix,
	const vector<vector<Particle>>& v_distances_matrix,
	const	int& v_size);

double GetOneIntEnergy(
	const vector<vector<double>>& matrix,
	const int& particle_index,
	const int& v_size);