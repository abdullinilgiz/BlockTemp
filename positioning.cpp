#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

int main() {

    ofstream data;
    data.open("positioning.txt", ios::app);

    int n = 10;
    int m = 10;
    int l = 10;

    double R = 1;
    double d = 2*R;
    double h = R * sqrt(3);
    double y_d;
    double z_h;
    double z_d;

    for (int z = 0; z < l; z++) {
        if (z % 2 == 1){
            z_h = R / sqrt(3);
            z_d = R;
        }
        else {
            z_h = 0;
            z_d = 0;
        }
            for (int y = 0; y < m; y++) {
                if (y % 2 == 0) {
                    y_d = 0;
                } else {
                    y_d = R;
                }
                for (int x = 0; x < n; x++) {
                    data << z_d + y_d + x * d << " " << z_h + y * h  << " " << z * h << '\n';
                }
            }
        }

    return 0;
}