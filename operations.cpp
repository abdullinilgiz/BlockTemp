#include "operations.h"

double DotP(Particle l, Particle r) {
	return l.x * r.x + l.y * r.y + l.z * r.z;
}

Particle operator+ (const Particle& l, const Particle& r) {
	return { l.x + r.x, l.y + r.y, l.z + r.z };
}
Particle operator- (const Particle& l, const Particle& r) {
	return { l.x - r.x, l.y - r.y, l.z - r.z };
}
Particle operator* (const Particle& l, const double& r) {
	return { l.x * r , l.y * r, l.z * r };
}
Particle operator/ (const Particle& l, const double& r) {
	return { l.x / r , l.y / r, l.z / r };
}
ostream& operator<< (ostream& os, const Particle& l) {
	os << l.x << " " << l.y << " " << l.z;
	return os;
}

Particle ParticleVectorAvr(const vector<Particle>& input_vector) {
	
	Particle vector_avr = { 0,0,0 };
	for (const auto& item : input_vector) {
		vector_avr = vector_avr + item;
	}
	vector_avr = vector_avr / input_vector.size();

	return vector_avr;
}

double ParticleVectorStDiv(const Particle& avr, const vector<Particle>& vector) {
	
	Particle delta;
	double error = 0;
	for (const auto& item : vector) {
		delta = item - avr;
		error += DotP(delta, delta);
	}
	error = sqrt(error / vector.size());

	return error;
}