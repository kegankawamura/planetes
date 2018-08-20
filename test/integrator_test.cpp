#include <iostream>
#include <string>
#include <math.h>
#include <new>
#include <Eigen/Dense>
#include "../tools/NumIntegrators.h"

using namespace std;
using namespace Eigen;

// x = t^2 --> dx/dt = 2t
double parabola (double x, double t) {
    return 2.0*t;
}

void runge_kutta_t1 () {
    double dt = 0.01;
    const double t0 = 0;
    const double x0 = 1.0;
    int num_steps = 1000;
    int info_steps = num_steps/10;
    double * x = new (nothrow) double[num_steps];
    if ( x == NULL ) {
	   cout << "Memory could not be assigned dynamically." << endl;
	   return;
    }
    double t = t0;
    x[0] = x0;
    for (int i = 0; i < num_steps-1; i++) {
        x[i+1] = runge_kutta_4th<double, double>(x[i], t, dt, parabola);
	   t += dt;
    	if (i % info_steps == 0 || i == num_steps - 2) {
    	    cout << "Timestep:\t" << i+1 << endl;
    	    cout << "Time:\t" << t << endl;
    	    cout << "Estimated:\t" << x[i+1] << endl;
    	    cout << "Actual:\t" << x0 + pow(t, 2.0) << endl;
    	    cout << "Error:\t" << fabs(x[i+1] - x0 - pow(t, 2.0)) << endl;
    	    cout << endl;
    	}
    } 
    delete[] x;
}


int main() {
    runge_kutta_t1();
    return 0;
}
