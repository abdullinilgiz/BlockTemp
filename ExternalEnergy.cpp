#include "ExternalEnergy.h"
#include "operations.h"


double GetOneExtEnergy(
	const Particle& moment,
	const Particle& E_ext) 
{
	double energy = -DotP(E_ext, moment);

	return M0 * energy;
}

double ExtSysEnergy(
	const vector<Particle>& moments,
	const int& v_size,
	const Particle& E_ext)
{
	double energy = 0;
	for (int i = 0; i < v_size; ++i) {
		energy += GetOneExtEnergy(moments[i], E_ext);
	}

	return energy;
}

