#ifndef NUM_INTEGRATORS_H__
#define NUM_INTEGRATORS_H__

/*
*	Includes numerical integration techniques.
*	Runge-Kutta
*/

#include <Eigen/Dense>


// Eigen::VectorXd runge_kutta_4th(const Eigen::VectorXd & x, double t, double dt, Eigen::VectorXd (*f)(Eigen::VectorXd, double));

// double runge_kutta_4th(const double & x, double t, double dt, double (*f)(double, double));

template <class X, class T>
X runge_kutta_4th(X x, T t, T dt, X (*f)(X, T)) {
	X k1, k2, k3, k4;
	k1 = f(x, t);
	k2 = f(x + k1*dt/2.0, t + dt/2.0);
	k3 = f(x + k2*dt/2.0, t + dt/2.0);
	k4 = f(x + k3*dt, t + dt);
	return x + (k1/6.0 + k2/3.0 + k3/3.0 + k4/6.0)*dt; 
}
#endif