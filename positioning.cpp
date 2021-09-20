#include <iostream>
#include <fstream>
#include <cmath>
#include "profile.h"

using namespace std;

int main() {

    ofstream data;
    data.open("positioning.txt");

    int n = 5;

    double R = 1;
    double x_step = 2 * R;
    double y_step = R * sqrt(3);
    double z_step = 2 *  R * sqrt(6) / 3;
    double dy_z;
    double dx_y;
    double dx_z;

    {
        LOG_DURATION("Zapolnenie")
        for (int z = 0; z < n; z++) {
            if (z % 2 == 1) {
                dy_z = R / sqrt(3);
                dx_z = R;
            } else {
                dx_z = 0;
                dy_z = 0;
            }
            for (int y = 0; y < n; y++) {
                if (y % 2 == 1) {
                    dx_y = -R;
                } else {
                    dx_y = 0;
                }
                for (int x = 0; x < n; x++) {
                    data << x * x_step + dx_y + dx_z<< " " << y * y_step + dy_z << " " << z * z_step << '\n';
                }
            }
        }
    }

    return 0;
}