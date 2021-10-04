#include "MonteCarlo.h"

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
	const Particle& E_ext) 
{

	//we should choose random particle
	int rng_particle = static_cast<int>(Nmb_particles * rg->Random());

	//calculate old moment energy
    const double old_int_particle_nrg = GetOneIntEnergy(int_nrg_matrix, rng_particle, Nmb_particles);
    const double old_ext_particle_nrg = GetOneExtEnergy(particle_moments[rng_particle], E_ext);
    const double old_easy_particle_nrg = particles_easy_nrg[rng_particle];

	//save old interaction energy and external field energy
    //const double old_full_sys_nrg = int_sys_nrg + ext_sys_nrg + easy_sys_nrg;

	//save matrix energies for random_particle_number
	vector<double> old_particle_nrgs(Nmb_particles);
	for (int index3 = 0; index3 < Nmb_particles; ++index3) {
		old_particle_nrgs[index3] = int_nrg_matrix[index3][rng_particle];
	}

	Particle old_moment;
	//save old moment value
	old_moment.x = particle_moments[rng_particle].x;
	old_moment.y = particle_moments[rng_particle].y;
	old_moment.z = particle_moments[rng_particle].z;

	// new moment for randomed number particle
	{
		double d_phi = (2 * rg->Random() - 1) * phi_max;
		Particle axis = GetRandomDir(rg);
		RotateVector(particle_moments[rng_particle], axis, d_phi);
	}

	//new moment energy
	const double new_easy_particle_nrg = OneEasyAxis(particle_moments[rng_particle], easy_axis_dir[rng_particle]);

	//new energies
	const double delta_int_sys_nrg =
            IntParticleEnergy(particle_moments, int_nrg_matrix, distances, v_distances, rng_particle)
            - old_int_particle_nrg;
	const double delta_ext_sys_nrg = GetOneExtEnergy(particle_moments[rng_particle], E_ext)
            - old_ext_particle_nrg;
	const double delta_easy_sys_nrg = new_easy_particle_nrg - old_easy_particle_nrg;

    //delta energy
    const double delta_full_sys_nrg = delta_int_sys_nrg + delta_ext_sys_nrg + delta_easy_sys_nrg;

	//apply or not
	double ksi = rg->Random();
	double porog = exp(-delta_full_sys_nrg / kB / Temper);
	//cout << "Porog" << porog << endl;
	if (ksi < porog) {
		full_sys_nrg += delta_full_sys_nrg;
		int_sys_nrg += delta_int_sys_nrg;
		ext_sys_nrg += delta_ext_sys_nrg;
		easy_sys_nrg += delta_easy_sys_nrg;
		particles_easy_nrg[rng_particle] = new_easy_particle_nrg;
		return true;
	}
	else {
		particle_moments[rng_particle].x = old_moment.x;
		particle_moments[rng_particle].y = old_moment.y;
		particle_moments[rng_particle].z = old_moment.z;
		//old energy matrix
		for (int j = 0; j < Nmb_particles; ++j) {
			int_nrg_matrix[j][rng_particle] = old_particle_nrgs[j];
			int_nrg_matrix[rng_particle][j] = old_particle_nrgs[j];
		}
		return false;
	}
}

void AngleChanger(const int& apply_counter, double& phi_max) {

	double koef = static_cast<double>(apply_counter) / (static_cast<double>(Nstep));
	//cout << "Koef: " << koef << endl;
	if (koef > 0.6) {
		if (phi_max < 1.1) {
			phi_max *= 1.1;
			cout << "Action: increase" << endl;
		}
	}
	if (koef < 0.4) {
		if (phi_max > 0.001) {
			phi_max *= 0.8;
			cout << "Action: decrease" << endl;
		}
	}
}