#pragma once

#include <sstream>
#include <stdexcept>
#include <map>
#include <set>
#include "mylibr.h"

using namespace std;

template <class T>
ostream& operator << (ostream& os, const vector<T>& s) {
	os << "{";
	bool first = true;
	for (const auto& x : s) {
		if (!first) {
			os << ", ";
		}
		first = false;
		os << x;
	}
	return os << "}";
}

template <class T>
ostream& operator << (ostream& os, const set<T>& s) {
	os << "{";
	bool first = true;
	for (const auto& x : s) {
		if (!first) {
			os << ", ";
		}
		first = false;
		os << x;
	}
	return os << "}";
}

template <class K, class V>
ostream& operator << (ostream& os, const map<K, V>& m) {
	os << "{";
	bool first = true;
	for (const auto& kv : m) {
		if (!first) {
			os << ", ";
		}
		first = false;
		os << kv.first << ": " << kv.second;
	}
	return os << "}";
}

template<class T, class U>
void AssertEqual(const T& t, const U& u, const string& hint = {}) {
	if (!(t == u)) {
		ostringstream os;
		os << "Assertion failed: " << t << " != " << u;
		if (!hint.empty()) {
			os << " hint: " << hint;
		}
		throw runtime_error(os.str());
	}
}

inline void Assert(bool b, const string& hint) {
	AssertEqual(b, true, hint);
}

class TestRunner {
public:
	template <class TestFunc>
	void RunTest(TestFunc func, const string& test_name) {
		try {
			func();
			cerr << test_name << " OK" << endl;
		}
		catch (exception & e) {
			++fail_count;
			cerr << test_name << " fail: " << e.what() << endl;
		}
		catch (...) {
			++fail_count;
			cerr << "Unknown exception caught" << endl;
		}
	}

	~TestRunner() {
		if (fail_count > 0) {
			cerr << fail_count << " unit tests failed. Terminate" << endl;
			//exit(1);
		}
	}

private:
	int fail_count = 0;
};

#define ASSERT_EQUAL(x, y) {            \
  ostringstream os;                     \
  os << #x << " != " << #y << ", "      \
    << __FILE__ << ":" << __LINE__;     \
  AssertEqual(x, y, os.str());          \
}

#define ASSERT(x) {                     \
  ostringstream os;                     \
  os << #x << " is false, "             \
    << __FILE__ << ":" << __LINE__;     \
  Assert(x, os.str());                  \
}

#define RUN_TEST(tr, func) \
  tr.RunTest(func, #func)


void TestInteraction();
void TestExt();
void TestEasy();
void TestRotateVector();
void TestRandomDir();
void TestRandomTransmission ();
double Trans_rg(TRandomMersenne* rg);
void TestMK_interaction_nrg(const double& int_sys_nrg,
                            const vector<vector<double>>& int_nrg_matrix,
                            const vector<Particle>& particle_moments,
                            const vector<vector<double>>& distances,
                            const vector<vector<Particle>>& v_distances,
                            const int& Nmb_particles);

void TestMK_Ext_nrg(const double& ext_sys_nrg,
                    const vector<Particle>& particle_moments,
                    const int& Nmb_particles,
                    const Particle& E_ext);

void TestMK_Easy_nrg(const double& easy_sys_nrg,
                     const vector<double>& particles_easy_nrg,
                     const vector<Particle>& particle_moments,
                     const vector<Particle>&easy_axis_dir,
                     const int& Nmb_particles);

#define OUTPUT(stream, param) 	{stream << #param << ": " << param << '\n';}  

void OutPutParam(ostream& stream);
