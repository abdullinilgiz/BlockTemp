#include "EasyAxisEnergy.h"
#include "ExternalEnergy.h"
#include "InteractionEnergy.h"
#include "Moments.h"
#include "operations.h"
#include "randomc.h"
#include "MonteCarlo.h"
#include "PhiThetaSurface.h"
#include "test_runner.h"

int main() {

	string date = Today + "_kA=" + to_string(static_cast<int>(kA));
	cout << "kA = " << kA << endl;

    ofstream stream;
    stream.open(date + "_Constants.txt" , ios::app);
    OutPutParam(stream);
    stream.close();

	cout.precision(12);
//	{
//		TestRunner tr;
//		RUN_TEST(tr, TestExt);
//		RUN_TEST(tr, TestInteraction);
//		RUN_TEST(tr, TestRandomTransmission);
//		RUN_TEST(tr, TestRotateVector);
//		RUN_TEST(tr, TestRandomDir);
//		RUN_TEST(tr, TestEasy);
//	}

	// create a vector of coords
	vector <Particle> particle_coords(Nmb_particles);

	//put coords inside vector
    FillVectorWithCoords(particle_coords);

	//create matrix of distances
	Particle r_v = { 0, 0, 0 };
	vector<vector<Particle>> v_distances(Nmb_particles, vector<Particle>(Nmb_particles));
	vector<vector<double>> distances(Nmb_particles, vector<double>(Nmb_particles));
	for (int i = 0; i < Nmb_particles; ++i) {
		for (int j = i; j < Nmb_particles; ++j) {
			if (i == j) {
				r_v = { 0, 0, 0 };
				distances[i][j] = 0;
			}
			else {
				//x summand
				r_v.x = particle_coords[i].x - particle_coords[j].x;
				if (r_v.x > XHalfSideLength) {
					r_v.x = particle_coords[i].x - particle_coords[j].x - XSideLength;
				}
				if (r_v.x < -XHalfSideLength) {
					r_v.x = particle_coords[i].x - particle_coords[j].x + XSideLength;
				}
				//y summand
				r_v.y = particle_coords[i].y - particle_coords[j].y;
				if (r_v.y > YHalfSideLength) {
					r_v.y = particle_coords[i].y - particle_coords[j].y - YSideLength;
				}
				if (r_v.y < -YHalfSideLength) {
					r_v.y = particle_coords[i].y - particle_coords[j].y + YSideLength;
				}
				//z summand 
				r_v.z = particle_coords[i].z - particle_coords[j].z;
				if (r_v.z > ZHalfSideLength) {
					r_v.z = particle_coords[i].z - particle_coords[j].z - ZSideLength;
				}
				if (r_v.z < -ZHalfSideLength) {
					r_v.z = particle_coords[i].z - particle_coords[j].z + ZSideLength;
				}
			}
			v_distances[i][j] = r_v;
			v_distances[j][i] = v_distances[i][j];
			distances[i][j] = sqrt(DotP(r_v, r_v));
			distances[j][i] = distances[i][j];
		}
	}
	

	//create vector of moments and easy axis
	vector <Particle> particle_moments(Nmb_particles);
	vector <Particle> easy_axis_dir(Nmb_particles);

	//put moments inside vector
	//put directions inside easy axis vector
	
	TRandomMersenne* rg = new TRandomMersenne(iSeed);
	Particle E_ext = { 0, 0, 0 };
	string marker = "ZFC";
	for (int i1 = 0; i1 < 2; i1++) {

		for (int i2 = 0; i2 < Nmb_particles; ++i2) {
			particle_moments[i2] = GetRandomDir(rg) * dMoment;
			easy_axis_dir[i2] = GetRandomDir(rg);
		}
	
		//full system moment
		Particle sys_moment = FullSysMoment(particle_moments, Nmb_particles);
		cout << " *** " << "Sys moment" << '\n';
		cout << sys_moment.x << '\n';
		cout << sys_moment.y << '\n';
		cout << sys_moment.z << endl;


		//create energy matrix
		vector<vector<double>> int_nrg_matrix(Nmb_particles, vector<double>(Nmb_particles, 0));
		vector<double> particles_easy_nrg(Nmb_particles, 0);

		//Metropolis algorithm

		double phi_max = PI / 3;
		int total_apply_counter = 0;
		int apply_counter = 0;
		double int_sys_nrg, ext_sys_nrg, easy_sys_nrg;
			 
		ofstream zfc_fc;
        zfc_fc.open(date + "_" + marker + ".txt" , ios::app);

		int decreaser = 4;
		int down_up = 2;
		int Temper = Tmax;
		int dT = -DeltaT; //should be negative
		Particle sys_moment_avr = { 0, 0, 0 };
		double moment_error = 0;

		if (marker == "FC"){

			E_ext = { 0, 0, 8000 };
			decreaser = 1;
			down_up = 1; //2 - zfc; 1 - fc
		}

		for (int i3 = 0; i3 < down_up; i3++) {



			ofstream convergence;
			convergence.open(date + "_" + marker + "_convergence" + ".txt", ios::app);

            ofstream koef_stream;
            koef_stream.open(date + "_" + marker + "_koef" + ".txt", ios::app);

            ofstream moments;
            moments.open(date + "_" + marker + "_moments" + ".txt", ios::app);

            ofstream delta_nrg;
            delta_nrg.open(date + "_" + marker + "_delta_nrg" + ".txt", ios::app);
            delta_nrg << "Temper" << "Intr - " << "Extr - " << "Easy - " << '\n';

            PhiThetaSurface(particle_coords, v_distances, distances, particle_moments,
                            easy_axis_dir, E_ext, marker + "_init");

			for (Temper; Temper > (Tmin - 1) && Temper < (Tmax + 1); Temper += dT) {

				cout << "                           New Temperature: " << Temper << '\n'
				<< "                           Field " << E_ext << '\n';

				vector<Particle> sys_moment_vector;

				//need to recalculate because of "double" type accuracy 
				//calculate energy of the system and save interactions matrix
				int_sys_nrg = IntSysEnergy(particle_moments, int_nrg_matrix, distances, v_distances, Nmb_particles);

				//calculate energy of interaction with E_external
				ext_sys_nrg = ExtSysEnergy(particle_moments, Nmb_particles, E_ext);

				//easy axis sum energy
				easy_sys_nrg = EasySysEnergy(particle_moments, particles_easy_nrg, easy_axis_dir, Nmb_particles);

				//full system energy
				double full_sys_nrg = int_sys_nrg + ext_sys_nrg + easy_sys_nrg;
				cout << " ***  START NRGS" << '\n';
				cout << "Intr: " << int_sys_nrg << '\n';
				cout << "Extr: " << ext_sys_nrg << '\n';
				cout << "Easy: " << easy_sys_nrg << '\n';
				cout << "Full: " << full_sys_nrg << '\n';
				cout << " *** " << endl;

				//phi_max = PI / 4;
				total_apply_counter = 0;
				apply_counter = 0;

				for (int index1 = 0; index1 < Nstep / decreaser ; ++index1) {
					for (int index2 = 0; index2 < Nstep; ++index2) {
						if (MKIteration(
                                Temper,
                                rg,
                                particles_easy_nrg,
                                int_nrg_matrix,
                                v_distances,
                                distances,
                                particle_moments,
                                int_sys_nrg,
                                ext_sys_nrg,
                                easy_sys_nrg,
                                full_sys_nrg,
                                phi_max,
                                easy_axis_dir,
                                E_ext)) {
							++apply_counter;
						}
					}
					AngleChanger(apply_counter, phi_max);
					total_apply_counter += apply_counter;

                    sys_moment = FullSysMoment(particle_moments, Nmb_particles);

                    if (decreaser == 1) {
                        koef_stream << static_cast<double>(apply_counter) / Nstep << " " << phi_max << '\n';
                        convergence << Temper << " " << index1 << " "
							 << sqrt(DotP(sys_moment, sys_moment)) << '\n';
                        for (size_t id = 0; id < 1001; id += 200) {
                            moments << Temper << " " << id << " " << particle_moments[id] << "\n";
                        }
                        if (Nstep - Nstep / 5 <= index1) {
                            sys_moment_vector.push_back(sys_moment);
                        }

                    }

                    delta_nrg << Temper << " " << int_sys_nrg
                    << " " << ext_sys_nrg << " " << easy_sys_nrg << " " << '\n';

                    apply_counter = 0;
                }

                cout << " ***  END NRGS" << '\n';
                cout << "Intr: " << int_sys_nrg << '\n';
                cout << "Extr: " << ext_sys_nrg << '\n';
                cout << "Easy: " << easy_sys_nrg << '\n';
                cout << "Full: " << full_sys_nrg << '\n';
                cout << " *** " << '\n';


				//sys_moment = FullSysMoment(particle_moments, Nmb_particles);
				if (decreaser == 1) {
					sys_moment_avr = ParticleVectorAvr(sys_moment_vector);
					moment_error = ParticleVectorStDiv(sys_moment_avr, sys_moment_vector);
					zfc_fc << Temper << " " << sys_moment_avr << " " << sqrt(DotP(sys_moment_avr, sys_moment_avr));
					zfc_fc << " " << moment_error << '\n';
				}

				cout << "total_apply_counter: " << total_apply_counter << '\n';
				cout << "phi_max: " << phi_max << '\n';
				cout << "Koef: " << static_cast<double>(total_apply_counter) / Nstep / Nstep / decreaser << '\n';
				cout << " .*.*.*.*.*.*.*.*.*. " << endl;
                sys_moment_vector.clear();

                if (decreaser == 1) {
                    PhiThetaSurface(particle_coords, v_distances, distances, particle_moments,
                                    easy_axis_dir, E_ext, marker + to_string(Temper));
                }
			}
			dT *= -1;
 			Temper = Tmin;
			if (marker == "ZFC") {
				E_ext = { 0, 0, 8000 };
			}
			else {
				E_ext = { 0, 0, 0 };
			}
			decreaser = 1;

			convergence.close();
			koef_stream.close();
			delta_nrg.close();
			moments.close();
		}
		zfc_fc.close();
		marker = "FC";
	}

	delete rg;
	
	return 0;
}