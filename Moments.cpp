#include "Moments.h"
#include "operations.h"

int ClosestParticle(vector <Particle>& particle_coords, const vector<bool>& pressed_particles){
    double min_dist = 10000;
    double dist;
    int min_ind = 0;
    for (int i = 1; i < particle_coords.size(); ++i){
        dist = DotP(particle_coords[i], particle_coords[i]);
        if (min_dist > dist && !pressed_particles[i]){
            min_dist = dist;
            min_ind = i;
        }
    }
    return min_ind;
}

void ToCenter(vector<Particle>& particle_coords, int ind){
    Particle norm_dir = particle_coords[ind] / -sqrt(DotP(particle_coords[ind], particle_coords[ind]));
    double koef = 0.05;
    double min_dist = 10000, dist;
    double range_min = 1.95 * Particle_radius;
    Particle new_coords = particle_coords[ind];
    while (min_dist > range_min * range_min){
        particle_coords[ind] = new_coords;
        new_coords = particle_coords[ind] + norm_dir * Particle_radius * koef;
        min_dist = 10000;
        for (int i = 0; i < particle_coords.size(); ++i){
            dist = DotP(particle_coords[i] - new_coords, particle_coords[i] - new_coords);
            if (dist < min_dist && i != ind){
                min_dist = dist;
            }
        }
    }
}

void ParticlesPressing(vector <Particle>& particle_coords){
    vector<bool> pressed_particles(Nmb_particles, false);
    pressed_particles[0] = true;
    for (int i = 1; i < pressed_particles.size(); ++i) {
        int ind = ClosestParticle(particle_coords, pressed_particles);
        ToCenter(particle_coords, ind);
        pressed_particles[ind] = true;
    }
}

void ParticlesShift(vector <Particle>& particle_coords){
    double max_x = -100, max_y = -100, max_z = -100;
    double min_x = 100, min_y = 100, min_z = 100;
    for(const auto item : particle_coords){
        if (max_x < item.x) max_x = item.x;
        if (max_y < item.y) max_y = item.y;
        if (max_z < item.z) max_z = item.z;
        if (min_x > item.x) min_x = item.x;
        if (min_y > item.y) min_y = item.y;
        if (min_z > item.z) min_z = item.z;
    }
    for(auto& item : particle_coords){
        item.x -= min_x;
        item.y -= min_y;
        item.z -= min_z;
    }
    XSideLength = max_x - min_x - 2 * Particle_radius;
    XHalfSideLength = XSideLength / 2;
    YSideLength = max_y - min_y - 2 * Particle_radius;
    YHalfSideLength = YSideLength / 2;
    ZSideLength = max_z - min_z - 2 * Particle_radius;
    ZHalfSideLength = ZSideLength / 2;
    ShortestHalfSideLength = min(XSideLength, min(YSideLength, ZSideLength)) / 2.0;
    cout << XSideLength << " " << YSideLength << " " << ZSideLength << endl;
}

void FillVectorWithCoords(vector <Particle>& particle_coords){
    double L = Particle_diam * 22;
    double x,y,z;
    bool marker;
    TRandomMersenne* rg = new TRandomMersenne(iSeed);
    particle_coords[0] = Particle{0,0, 0};
    for(int i = 1; i < Nmb_particles; ++i){
        cout << i << endl;
        x = (1 - 2 * rg->Random()) * L / 2;
        y = (1 - 2 * rg->Random()) * L / 2;
        z = (1 - 2 * rg->Random()) * L / 2;
        marker = true;
        for (int j = 0; j <= i; ++j){
            if (sqrt(pow(particle_coords[j].x - x, 2) + pow(particle_coords[j].y - y, 2) +
                     pow(particle_coords[j].z - z, 2)) < 2 * Particle_radius){
                marker = false;
                --i;
                cout << "FAIL" << endl;
                break;
            }
        }
        if (marker){
            particle_coords[i] = Particle{x,y,z};
        }
    }

    ofstream stream;
    stream.open("initial positions.txt");
    for (const auto& item : particle_coords){
        stream << item << endl;
    }
    stream.close();

    ParticlesPressing(particle_coords);
    ParticlesShift(particle_coords);
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

Particle GetDir(const double phi, const double theta)
{
    Particle axis;
    axis.x = sin(theta) * cos(phi);
    axis.y = sin(theta) * sin(phi);
    axis.z = cos(theta);

    return axis;
}

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
