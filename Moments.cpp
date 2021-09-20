#pragma once
#include "Moments.h"

void FillVectorWithCoords(vector <Particle>& particle_coords){
    for (int i = 0; i < Nside; ++i) {
        for (int j = 0; j < Nside; ++j) {
            for (int k = 0; k < Nside; ++k) {
                particle_coords[i * Nside * Nside + j * Nside + k].x = Particle_radius + i * Particle_diam;
                particle_coords[i * Nside * Nside + j * Nside + k].y = Particle_radius + j * Particle_diam;
                particle_coords[i * Nside * Nside + j * Nside + k].z = Particle_radius + k * Particle_diam;
            }
        }
    }
}

Particle FullSysMoment(
	const vector<Particle>& moments, 
	const int& v_size) 
{
	Particle full_moment = { 0, 0, 0 };
	for (int i = 0; i < v_size; ++i) {
		full_moment.x += moments[i].x;
		full_moment.y += moments[i].y;
		full_moment.z += moments[i].z;
	}

	return full_moment;
}

//������� ���������� ��������� ������ ��������� �����, �.�. ���������� ��������� �����������
Particle GetRandomDir(TRandomMersenne* rg)
{
	Particle axis;
	double dPhi = 2.0 * PI * rg->Random();
	double dCosTheta = 2.0 * rg->Random() - 1.0;

	axis.x = sqrt(1 - dCosTheta * dCosTheta) * cos(dPhi);
	axis.y = sqrt(1 - dCosTheta * dCosTheta) * sin(dPhi);
	axis.z = dCosTheta;

	return axis;
}


//������� ������������ ������� ������� vector �� ���� phi ������������� ��� axis
void RotateVector(
	Particle& vector,
	const Particle& axis, 
	const double& phi)
{
	if (DotP(axis, axis) - 1 > 1e-15) {
		exit(1);
	}
	/*
	//��������, ��� ������, ������������ �������� ������������ ������ ��������� �����
	double buffer = axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2];

	if ((buffer - 1) * (buffer - 1) > 0.000000000001)
	{
		buffer = sqrt(buffer);
		axis[0] = axis[0] / buffer;
		axis[1] = axis[1] / buffer;
		axis[2] = axis[2] / buffer;
	}
	//����� ��������
	*/

	double CosPhi = cos(phi);
	double SinPhi = sin(phi);
	double OneMinCosPhi = 1 - CosPhi;

	Particle tempvector;//double tempvector[3];

	tempvector.x = (CosPhi + OneMinCosPhi * axis.x * axis.x) * vector.x +
		(axis.x * axis.y * OneMinCosPhi - axis.z * SinPhi) * vector.y +
		(axis.x * axis.z * OneMinCosPhi + axis.y * SinPhi) * vector.z;

	tempvector.y = (axis.x * axis.y * OneMinCosPhi + axis.z * SinPhi) * vector.x +
		(CosPhi + OneMinCosPhi * axis.y * axis.y) * vector.y +
		(axis.y * axis.z * OneMinCosPhi - axis.x * SinPhi) * vector.z;

	tempvector.z = (axis.x * axis.z * OneMinCosPhi - axis.y * SinPhi) * vector.x +
		(axis.y * axis.z * OneMinCosPhi + axis.x * SinPhi) * vector.y +
		(CosPhi + OneMinCosPhi * axis.z * axis.z) * vector.z;

	vector.x = tempvector.x;
	vector.y = tempvector.y;
	vector.z = tempvector.z;
}
