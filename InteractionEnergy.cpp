#pragma once

#include "InteractionEnergy.h"

double CalculateDipol_Dipol(
	const vector<Particle>& moments, 
	const vector<vector<double>>& distances_matrix,
	const vector<vector<Particle>>& v_distances_matrix,
	const int& index_1, 
	const int& index_2) 
{
	double energy = 0;
	double mm; 
	double mrmr3;
	Particle r_v = {0,0,0};
	double dis_r_v;

	// abs r_v
	dis_r_v = distances_matrix[index_1][index_2];
	if (dis_r_v > HalfSideLength) {
		return energy;
	}
	// first summand
	mm = DotP(moments[index_1], moments[index_2]);
	//get distance from saved matrix
	r_v = v_distances_matrix[index_1][index_2];
	//second summand
	mrmr3 = 3.0 * DotP(moments[index_1], r_v) * DotP(moments[index_2], r_v);
	

	energy = M0 / (4.0 * PI) * (dis_r_v * dis_r_v * mm - mrmr3) / pow(dis_r_v, 5);

	return energy;
}

double IntParticleEnergy(
	const vector<Particle>& moments,
	vector<vector<double>>& energy_matrix,
	const vector<vector<double>>& distances_matrix,
	const vector<vector<Particle>>& v_distances_matrix,
	const int& particle_index,
	const int& v_size) 
{
	double energy = 0;
	//double dis_r_v = 0;

	for (int i = 0; i < v_size; ++i) {
		if (i != particle_index) {
			energy_matrix[i][particle_index] = CalculateDipol_Dipol(moments, distances_matrix,
																		v_distances_matrix, i, particle_index);
			energy_matrix[particle_index][i] = energy_matrix[i][particle_index];
			energy += energy_matrix[i][particle_index];
		}
	}
	return energy;
}


double IntSysEnergy(
	const vector<Particle>& moments,
	vector<vector<double>>& matrix,
	const vector<vector<double>>& distances_matrix,
	const vector<vector<Particle>>& v_distances_matrix,
	const int& v_size)
{
	double energy = 0;
	// calculate one interaction only once
	for (int i = 0; i < v_size; ++i) {
		for (int j = i + 1; j < v_size; ++j) {
			matrix[i][j] = CalculateDipol_Dipol(moments, distances_matrix, v_distances_matrix, i, j);
			matrix[j][i] = matrix[i][j];
			energy += matrix[i][j];
		}
	}

	return energy;
}

double GetOneIntEnergy(
	const vector<vector<double>>& matrix,
	const int& particle_index,
	const int& v_size)
{
	double energy = 0;

	for (int i = 0; i < v_size; ++i) {
		energy += matrix[i][particle_index];
	}

	return energy;
}