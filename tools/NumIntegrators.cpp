#include "NumIntegrators.h"

using namespace Eigen;

/*VectorXd runge_kutta_4th(VectorXd x, double t, double dt, VectorXd (*f)(VectorXd, double)) {
	VectorXd k1, k2, k3, k4;
	k1 = f(x, t);
	k2 = f(x + k1*dt/2.0, t + dt/2.0);
	k3 = f(x + k2*dt/2.0, t + dt/2.0);
	k4 = f(x + k3*dt, t + dt);
	return x + (k1/6.0 + k2/3.0 + k3/3.0 + k4/6.0)*dt; 
}*/


