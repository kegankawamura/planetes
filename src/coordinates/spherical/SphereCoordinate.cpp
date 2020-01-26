#include <cstdlib>
#include <cmath>
#include "SphereCoordinate.h"



SphereCoordinate::SphereCoordinate() {
	setQ(0,0,0);
	setQDot(0,0,0);
}





void SphereCoordinate::setQ (double Rad, double theta, double phi) {
	if (phi < 0 || phi > M_PI) {
		// TODO: error
		phi = M_PI/2;
	}
	q[0] = abs(Rad);
	theta %= 2*M_PI;
	if (Rad < 0) {
		if (theta >= 0) {
			theta -= M_PI;
		} else {
			theta += M_PI;
		}
		phi = M_PI - phi;
	}
	q[1] = theta;
	sinTheta = sin(theta);
	cosTheta = cos(theta);
	q[2] = phi;
	sinPhi = sin(phi);
	cosPhi = cos(phi);
	computeR();
}

void CylinderCoordinate::setQDot (double Raddot, double thetadot, double phidot) {
	qdot[0] = Raddot;
	qdot[1] = thetadot;
	qdot[2] = phidot;
	computeRDot();
}


void CylinderCoordinate::computeR() {
	sinTheta = sin(q[1]);
	cosTheta = cos(q[1]);
	sinPhi = sin(q[2]);
	cosPhi = cos(q[2]);
	r <<	q[0]*cosTheta*sinPhi,
			q[0]*sinTheta*sinPhi, 
			q[0]*cosPhi;
}

void CylinderCoordinate::computeR() {
	sinTheta = sin(q[1]);
	cosTheta = cos(q[1]);
	sinPhi = sin(q[2]);
	cosPhi = cos(q[2]);
	rdot <<	qdot[0]*cosTheta*sinPhi + q[0]*(-qdot[1]*sinTheta*sinPhi + qdot[2]*cosTheta*cosPhi),
			qdot[0]*sinTheta*sinPhi + q[0]*(qdot[1]*cosTheta*sinPhi + qdot[2]*sinTheta*cosPhi), 
			qdot[0]*cosPhi - q[0]*qdot[2]*sinPhi;
}

// a_1 = e_r, a_2 = rad*e_theta, a_3 = E_3
void EuclidCooordinate::computeCovariant() {
	sinPhi = sin(q[2]);
	covariant[0] <<	1,
					0,
					0;
	covariant[1] <<	0,
					q[0]*q[0]*sinPhi*sinPhi,
					0;					
	covariant[2] <<	0,
					0,
					q[0]*q[0];
}

void EuclidCoordiante::computeChristoffel_first() {
	sinPhi = sin(q[2]);
	cosPhi = cos(q[2]);
	Christoffel_first[0][0] = 0;
	Christoffel_first[1][0] = 0;
	Christoffel_first[2][0] = 0;
	Christoffel_first[3][0] = -q[0]*sinPhi*sinPhi;
	Christoffel_first[4][0] = 0;
	Christoffel_first[5][0] = -q[0];
	Christoffel_first[0][1] = 0;
	Christoffel_first[1][1] = q[0]*sinPhi*sinPhi;
	Christoffel_first[2][1] = 0;
	Christoffel_first[3][1] = 0;
	Christoffel_first[4][1] = q[0]*q[0]*sinPhi*cosPhi;
	Christoffel_first[5][1] = 0;
	Christoffel_first[0][2] = 0;
	Christoffel_first[1][2] = 0;
	Christoffel_first[2][2] = q[0];
	Christoffel_first[3][2] = -q[0]*sinPhi*cosPhi;
	Christoffel_first[4][2] = 0;
	Christoffel_first[5][2] = 0;
}