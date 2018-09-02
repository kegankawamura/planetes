#include <iostream>
#include <string>
#include <math.h>
#include <complex>
#include <new>
#include <Eigen/Dense>
#include "../tools/NumIntegrators.h"

using namespace std;
using namespace Eigen;

// x = t^2 --> dx/dt = 2t
double parabola (double x, double t) {
    return 2.0*t;
}

// ma + cv + kx = 0
Vector2d oscillation(Vector2d x, double t) {
    Matrix2d sys;
    double k = 1;
    double c = 1;
    double m = 1;
    sys << 0, 1,
	   -k/m,-c/m;
    return sys * x;
}

void oscillation_analytical(const Vector2d & x0, double t, Vector2d * x) {
    double k = 1;
    double c = 1;
    double m = 1;
    double nat_freq = sqrt(k/m);
    double damp = c/(2*sqrt(m*k));
    double damp_freq = nat_freq * sqrt(1 - pow(damp,2));
    if (pow(c,2) < 4*k*m) {
        double c1, c2, decay;
        decay = c/(2*m);
        c1 = x0(0);
        c2 = 1/damp_freq * (x0(1) + c1*decay);
        (*x)(0) = c1*exp(-decay*t)*cos(damp_freq*t) 
                + c2*exp(-decay*t)*sin(damp_freq*t);
        (*x)(1) = -decay * c1*exp(-decay*t)*cos(damp_freq*t)
                - damp_freq * c1*exp(-decay*t)*sin(damp_freq*t) 
                - decay*c2*exp(-decay*t)*sin(damp_freq*t)
                + damp_freq * c2*exp(-decay*t)*cos(damp_freq*t);
    } else if (pow(c,2) > 4*k*m) {
        double r1, r2, c1, c2;
        r1 = (-c + sqrt(pow(c,2) - 4*k*m))/(2*m);
        r2 = (-c - sqrt(pow(c,2) - 4*k*m))/(2*m);
        c1 = (x0(1) - x0(0)*r2) / (2*damp_freq);
        c2 = x0(0) - c1;
        (*x)(0) = c1 * exp(r1*t) + c2 * exp(r2*t); 
        (*x)(1) = c1 * r1 * exp(r1*t) + c2 * r2 * exp(r2*t);
    } else {
        double r, c1, c2;
        r = -c/(2*m);
        c1 = x0(0);
        c2 = x0(1) - c1 * r;
        (*x)(0) = exp(r * t) * (c1 + c2*t);
        (*x)(1) = exp(r * t) * (r*c1 + r*c2*t + c2); 
    }
}

bool runge_kutta_t1 (bool print) {
    double dt = 0.01;
    const double t0 = 0;
    const double x0 = 1.0;
    const int num_steps = 1000;
    int info_steps = num_steps/10;
    double * x = new (nothrow) double[num_steps];
    if ( x == NULL ) {
	   cout << "Memory could not be assigned dynamically for Integrator_test test1." << endl;
	   return;
    }
    double t = t0;
    x[0] = x0;
    unsigned int err_cnt = 0;
    const unsigned int err_thrsh = 2;
    for (int i = 0; i < num_steps-1; i++) {
        x[i+1] = runge_kutta_4th<double, double>(x[i], t, dt, parabola);
	    t += dt;
        if (fabs(x[i+1] - x0 - pow(t, 2.0)) > 1e-6) {
            err_cnt++;
        }
    	if (print && i % info_steps == 0 || i == num_steps - 2) {
            
    	    cout << "Timestep:\t" << i+1 << endl;
    	    cout << "Time:\t" << t << endl;
    	    cout << "Estimated:\t" << x[i+1] << endl;
    	    cout << "Actual:\t" << x0 + pow(t, 2.0) << endl;
    	    cout << "Error:\t" << fabs(x[i+1] - x0 - pow(t, 2.0)) << endl;
	    }

        if (err_cnt >= err_thrsh) {
            delete[] x;
            return false;
        }
    }
    delete[] x;
    return true;
}

bool runge_kutta_t2(bool print) {
    double dt = 0.02;
    double len_t = 10;
    const double t0 = 0;
    Vector2d x0(1.0, 0);
    const int num_steps = floor(len_t/dt);
    int info_steps = floor(num_steps/4);
//    Matrix<double, num_steps, 2> x;
    double t = t0;
//    x << x0(0), x0(1);
    Vector2d xi = x0;
    Vector2d x_analytical;
    unsigned int err_cnt = 0;
    const unsigned int err_thrsh = 4;
    for (int i = 0; i < num_steps - 1; i++) {
	    xi = runge_kutta_4th<Vector2d, double> (xi, t, dt, oscillation);
	    t += dt;
        oscillation_analytical(x0, t, &x_analytical);
        if (fabs(xi(0) - x_analytical(0)) > 1e-6) {
            err_cnt++;
        }
        if (fabs(xi(1) - x_analytical(1)) > 1e-6) {
            err_cnt++;
        }
        if (print && i % info_steps == 0 || i == num_steps - 2) {
            cout << "Timestep:\t" << i+1 << endl;
            cout << "Time:\t" << t << endl;
            cout << "Estimated:\t" << '\n' << xi << endl;
            cout << "Actual:\t" << '\n' << x_analytical << endl;
            cout << "Error:\t" << '\n' << xi - x_analytical << endl;
        }
        if (err_cnt >= err_thrsh) {
            return false;
        }
    }

    return true;
}


int main(int argv, char* argc[]) {
    char p_flag[] = "-p";
    bool print = false;
    if (argv == 2 && strcmp(argc[1], p_flag) == 0) {
        print = true;
    }
    unsigned int fail = 0;
    short tests[] = {0, 1};
    if(tests[0] && !runge_kutta_t1(print))
        fail++;
    if (tests[1] && !runge_kutta_t2(print))
        fail++;
    cout << fail << " test(s) have failed" << endl;
    return 0;
}
