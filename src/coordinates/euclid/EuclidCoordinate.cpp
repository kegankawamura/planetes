#include <cstdlib>
#include <cmath>
#include "EuclidCoordinate.h"



EuclidCoordinate::EuclidCoordinate() {
	setQ(0,0,0);
	setQDot(0,0,0);
}





void EuclidCoordinate::setQ (double rad, double theta, double z) {
	q[0] = abs(rad);
	if (rad < 0) {
		if (theta >= 0) {
			theta += M_PI;
		} else {
			theta -= M_PI;
		}
	}
	q[1] = theta;
	q[2] = z;
	computeR();
}

void EuclidCoordinate::setQDot (double raddot, double thetadot, double zdot) {
	qdot[0] = raddot;
	qdot[1] = thetadot;
	qdot[2] = zdot;
	computeRDot();
}


void EuclidCoordinate::computeR() {
	r <<	q[0]*cos(q[1]),
			q[0]*sin(q[1]), 
			q[2];
}

void EuclidCoordinate::computeR() {
	rdot <<	qdot[0],
			qdot[1], 
			qdot[2];
}


void EuclidCooordinate::computeCovariant() {
	covariant[0] <<	1,
					0,
					0;
	covariant[1] <<	0,
					1,
					0;					
	covariant[2] <<	0,
					0,
					1;
}

void EuclidCoordiante::computeChristoffel_first() {
	Christoffel_first[0][0] = 0;
	Christoffel_first[1][0] = 0;
	Christoffel_first[2][0] = 0;
	Christoffel_first[3][0] = 0;
	Christoffel_first[4][0] = 0;
	Christoffel_first[5][0] = 0;
	Christoffel_first[0][1] = 0;
	Christoffel_first[1][1] = 0;
	Christoffel_first[2][1] = 0;
	Christoffel_first[3][1] = 0;
	Christoffel_first[4][1] = 0;
	Christoffel_first[5][1] = 0;
	Christoffel_first[0][2] = 0;
	Christoffel_first[1][2] = 0;
	Christoffel_first[2][2] = 0;
	Christoffel_first[3][2] = 0;
	Christoffel_first[4][2] = 0;
	Christoffel_first[5][2] = 0;
}