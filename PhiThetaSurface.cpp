//
// Created by Ya on 02.10.2021.
//

#include "PhiThetaSurface.h"

void PhiThetaSurface(const vector<Particle>& coords,
                     const vector<vector<Particle>>& v_distances,
                     const vector<vector<double>>& distances,
                     vector <Particle>& moments,
                     const vector <Particle>& easy_axis_dir,
                     const Particle& E_ext,
                     const string& marker){
    vector<int> inds = {0, 200, 400, 600, 800};
    size_t points = 101;
    vector<double> phi(points);
    for (size_t i = 0; i < points; ++i){
        phi[i] = 2 * PI / points * i;
    }
    vector<double> theta(points);
    for (size_t i = 0; i < points; ++i){
        theta[i] = PI / points * i;
    }
    for (auto item : inds) {

        ofstream stream;
        stream.open(Today + "_phitheta_" + marker + "_index" + to_string(item) + ".txt" , ios::app);

        Particle old_moment = moments[item];

        for (size_t i2 = 0; i2 < points; ++i2) {
            for (size_t i3 = 0; i3 < points; ++i3) {
                moments[item] = GetDir(phi[i2], theta[i3]) * dMoment;
                double ext_nrg = GetOneExtEnergy(moments[item], E_ext);
                double easy_nrg = OneEasyAxis(moments[item], easy_axis_dir[item]);

                double int_nrg = 0;
                for (int i = 0; i < moments.size(); ++i) {
                    if (i != item) {
                        int_nrg += CalculateDipol_Dipol(moments, distances, v_distances, i, item);
                    }
                }

                stream << phi[i2] << " " << theta[i3] << " "
                       << ext_nrg << " " << easy_nrg << " " << int_nrg << '\n';
            }
        }

        stream.close();

        moments[item] = old_moment;
    }

}
