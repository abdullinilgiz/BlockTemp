//
// Created by Ya on 02.10.2021.
//
#include "mylibr.h"
#include "Moments.h"
#include "operations.h"
#include "ExternalEnergy.h"
#include "EasyAxisEnergy.h"
#include "InteractionEnergy.h"


void PhiThetaSurface(const vector<Particle>& coords,
                     const vector<vector<Particle>>& v_distances,
                     const vector<vector<double>>& distances,
                     vector <Particle>& moments,
                     const vector <Particle>& easy_axis_dir,
                     const Particle& E_ext,
                     const string& marker);

