#pragma once

#include "mylibr.h"

double GetOneExtEnergy(
        const Particle& moment,
        const Particle& E_ext);

double ExtSysEnergy(
	const vector<Particle>& moments,
	const int& v_size,
	const Particle& E_ext);