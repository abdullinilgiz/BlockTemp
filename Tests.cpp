#include "test_runner.h"
#include "EasyAxisEnergy.h"
#include "ExternalEnergy.h"
#include "InteractionEnergy.h"
#include "Moments.h"
#include "operations.h"

void TestInteraction(){
    int particles_number = 4;

    vector<Particle> moments(particles_number);

    moments[0].x = 2;
    moments[0].y = 0;
    moments[0].z = 3;

    moments[1].x = 1;
    moments[1].y = 2;
    moments[1].z = 1;

    moments[2].x = 2;
    moments[2].y = 1;
    moments[2].z = 2;

    moments[3].x = 3;
    moments[3].y = 2;
    moments[3].z = 1;


    vector<vector<double>> matrix(particles_number, vector<double>(particles_number, 0));
    vector<vector<double>> distances(particles_number, vector<double>(particles_number, 0));

    distances[0][1] = 2e-8;
    distances[1][0] = 2e-8;
    distances[0][2] = 2e-8;
    distances[2][0] = 2e-8;
    distances[0][3] = 2e-8;
    distances[3][0] = 2e-8;

    distances[1][2] = 2e-8;
    distances[2][1] = 2e-8;
    distances[1][3] = 2e-8;
    distances[3][1] = 2e-8;

    distances[3][2] = 2e-8;
    distances[2][3] = 2e-8;

    vector<vector<Particle>> v_distances(particles_number, vector<Particle>(particles_number));

    v_distances[0][0] = {0, 0, 0};
    v_distances[1][1] = {0, 0, 0};
    v_distances[2][2] =  {0, 0, 0};
    v_distances[3][3] =  {0, 0, 0};

    v_distances[0][1] = { 4, 1, 4 };
    v_distances[1][0] = v_distances[1][0];
    v_distances[2][0] = { 2, 3, 2 };
    v_distances[0][2] = v_distances[2][0];
    v_distances[3][0] = { 2, 3, 2 };
    v_distances[0][3] = v_distances[3][0];

    v_distances[1][2] = { 3, 2, 3 };
    v_distances[2][1] = v_distances[1][2];
    v_distances[1][3] = { 3, 2, 3 };
    v_distances[3][1] = v_distances[1][2];

    v_distances[2][3] = { 3, 2, 3 };
    v_distances[3][2] = v_distances[1][2];

    double int_sys_energy = IntSysEnergy(moments, matrix, distances, v_distances, 3);
    double sum_energy = 0;
    double one_inter;
    for (int i = 0; i < 3; ++i){
        for (int j = i + 1; j < 3; ++j){
            one_inter = M0 / 4 / PI * (DotP(moments[i],moments[j]) / pow(distances[i][j],3)
                    - 3 * DotP(moments[i], v_distances[i][j]) * DotP(moments[j], v_distances[i][j]) / pow(distances[i][j],5));

            ASSERT_EQUAL(abs((one_inter - matrix[i][j]) / one_inter) < 1e-15, true);
            sum_energy += one_inter;
        }
    }
    ASSERT_EQUAL(sum_energy, int_sys_energy);
}

void TestExt() {
	vector<Particle> moments(3);
	moments[1].x = 5;
	moments[1].y = 5;
	moments[1].z = 0;

	moments[0].x = 0;
	moments[0].y = 0;
	moments[0].z = 5;

	moments[2].x = 1;
	moments[2].y = 2;
	moments[2].z = 3;

	Particle E_ext = {100, 100, 100};

	double expected = 0;
	for (int i = 0; i < 3; ++i){
	    expected += - M0 * DotP(E_ext,moments[i]);
	}

	double ext_sys_nrg = ExtSysEnergy(moments, 3, E_ext);
    ASSERT_EQUAL(expected, ext_sys_nrg);
}

void TestEasy() {

	vector<Particle> moments(3);
	moments[1].x = 2;
	moments[1].y = 1;
	moments[1].z = 2;
	moments[0].x = 1;
	moments[0].y = 2;
	moments[0].z = 1;
	moments[2].x = 3;
	moments[2].y = 0;
	moments[2].z = 1;

	vector<double> easy_nrg(3, 0);
	vector<Particle> easy_axis_dir(3);

    TRandomMersenne* Test = new TRandomMersenne(iSeed);
	for (int i = 0; i < 3; ++i){
	    easy_axis_dir[i] = GetRandomDir(Test);
	}

	double expected = EasySysEnergy(moments, easy_nrg, easy_axis_dir, 3);
	double cal_sum = 0;
	double one_easy;
	for (int i = 0; i < 3; ++i){
	    one_easy = kA * Vol * (1 - pow(DotP(moments[i], easy_axis_dir[i]),2)/DotP(moments[i],moments[i]));
	    ASSERT_EQUAL(one_easy, easy_nrg[i]);
	    cal_sum += one_easy;
	}
	ASSERT_EQUAL(cal_sum, expected);
}
void TestRotateVector() {
	TRandomMersenne* Test = new TRandomMersenne(iSeed);
	Particle vector = GetRandomDir(Test) * 1;
	Particle axis = GetRandomDir(Test);
	const double d_phi = PI / 2;
	Particle first = vector;
	for (int i = 0; i < 8; ++i) {
		RotateVector(vector, axis, d_phi);
	}
	double zero = 0;
	ASSERT_EQUAL(first.x - vector.x < 1e-15, true);
	ASSERT_EQUAL(first.y - vector.y < 1e-15, true);
	ASSERT_EQUAL(first.z - vector.z < 1e-15, true);




}
void TestRandomDir() {
    TRandomMersenne* rg = new TRandomMersenne(iSeed);

	Particle fi = GetRandomDir(rg);
	Particle s = GetRandomDir(rg);
	Particle th = GetRandomDir(rg);
	Particle fo = GetRandomDir(rg);
    double first = DotP(fi, fi) - 1;
    double second = DotP(s, s) - 1;
    double third = DotP(th, th) - 1;
    double fourth = DotP(fo, fo) - 1;
    double zero = 0;
	ASSERT_EQUAL(first < 1e-16, true);
	ASSERT_EQUAL(second < 1e-16, true);
	ASSERT_EQUAL(third < 1e-16, true);
	ASSERT_EQUAL(fourth < 1e-16, true);
}
double Trans_rg(TRandomMersenne* rg){
    return rg->Random();
}

void TestRandomTransmission (){

    TRandomMersenne* rg = new TRandomMersenne(iSeed);
    double first = Trans_rg(rg);
    double second = Trans_rg(rg);
    double third = Trans_rg(rg);

    ASSERT_EQUAL((first == second || first == third || second == third), false);
}

void TestMK_interaction_nrg(const double& int_sys_nrg,
                            const vector<vector<double>>& int_nrg_matrix,
                            const vector<Particle>& particle_moments,
                            const vector<vector<double>>& distances,
                            const vector<vector<Particle>>& v_distances,
                            const int& Nmb_particles){

    vector<vector<double>> test_int_nrg_matrix(Nmb_particles, vector<double>(Nmb_particles, 0));
    double test_int_sys_nrg = IntSysEnergy(particle_moments, test_int_nrg_matrix, distances, v_distances, Nmb_particles);
    bool net = false;
    ASSERT_EQUAL(test_int_sys_nrg - int_sys_nrg > 1e-30, net);
    for (int i=0 ; i < Nmb_particles; ++i){
        for (int j = 0; j < Nmb_particles; ++j){
            ASSERT_EQUAL(test_int_nrg_matrix[i][j] - int_nrg_matrix[i][j] > 1e-35, net );
        }

    }
}

void TestMK_Ext_nrg(const double& ext_sys_nrg,
                    const vector<Particle>& particle_moments,
                    const int& Nmb_particles,
                    const Particle& E_ext){

    double test_ext_system_nrg = ExtSysEnergy(particle_moments, Nmb_particles, E_ext);
    ASSERT_EQUAL(test_ext_system_nrg, ext_sys_nrg);

}

void TestMK_Easy_nrg(const double& easy_sys_nrg,
                      const vector<double>& particles_easy_nrg,
                      const vector<Particle>& particle_moments,
                      const vector<Particle>&easy_axis_dir,
                      const int& Nmb_particles){
    vector<double> test_particles_easy_nrg(Nmb_particles, 0);
    double test_easy_sys_nrg = EasySysEnergy(particle_moments,test_particles_easy_nrg, easy_axis_dir, Nmb_particles);
    double zero = 0;
    ASSERT_EQUAL(test_easy_sys_nrg - easy_sys_nrg > 1e-31, false);
    for (int i=0; i < Nmb_particles; ++i){
        ASSERT_EQUAL(test_particles_easy_nrg[i], particles_easy_nrg[i]);
    }
}


void OutPutParam(ostream& stream) {
    OUTPUT(stream, Nstep);
    OUTPUT(stream, Rcore);
    OUTPUT(stream, Tshell);
    OUTPUT(stream, Nside_X);
    OUTPUT(stream, Nside_Y);
    OUTPUT(stream, Nside_Z);
    OUTPUT(stream, iSeed);
    OUTPUT(stream, dMoment);
    OUTPUT(stream, kA);
    OUTPUT(stream, Tmax);
    OUTPUT(stream, Tmin);
}

	
