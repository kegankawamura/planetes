#include <cstdlib>
#include <cmath>
#include "CylinderCoordinate.h"



CylinderCoordinate::CylidnerCoordinate() {
	setQ(0,0,0);
	setQDot(0,0,0);
}





void CylinderCoordinate::setQ (double rad, double theta, double z) {
	q[0] = abs(rad);
	theta %= 2*M_PI;
	if (rad < 0) {
		if (theta >= 0) {
			theta -= M_PI;
		} else {
			theta += M_PI;
		}
	}
	q[1] = theta;
	sinTheta = sin(theta);
	cosTheta = cos(theta);
	q[2] = z;
	computeR();
}

void CylinderCoordinate::setQDot (double raddot, double thetadot, double zdot) {
	qdot[0] = raddot;
	qdot[1] = thetadot;
	qdot[2] = zdot;
	computeRDot();
}


void CylinderCoordinate::computeR() {
	sinTheta = sin(q[1]);
	cosTheta = cos(q[1]);
	r <<	q[0]*cosTheta,
			q[0]*sinTheta, 
			q[2];
}

void CylinderCoordinate::computeR() {
	sinTheta = sin(q[1]);
	cosTheta = cos(q[1]);
	rdot <<	qdot[0]*cosTheta - q[0]*qdot[1]*sinTheta,
			qdot[0]*sinTheta + q[0]*qdot[1]*cosTheta, 
			qdot[2];
}

// a_1 = e_r, a_2 = rad*e_theta, a_3 = E_3
void CylinderCooordinate::computeCovariant() {
	/**
	covariant[0] <<	1,
					0,
					0;
	covariant[1] <<	0,
					q[0]*q[0],
					0;					
	covariant[2] <<	0,
					0,
					1;
	*/
	computeCovariant(q, covariant);
}


void CylinderCooordinate::computeCovariant(Vector3d Q, Vector3d * cov) {
	cov[0] <<	1,
				0,
				0;
	cov[1] <<	0,
				Q[0]*Q[0],
				0;					
	cov[2] <<	0,
				0,
				1;
}


void CylinderCoordiante::computeChristoffel_first() {
	/**Christoffel_first[0][0] = 0;
	Christoffel_first[1][0] = 0;
	Christoffel_first[2][0] = 0;
	Christoffel_first[3][0] = -q[0];
	Christoffel_first[4][0] = 0;
	Christoffel_first[5][0] = 0;
	Christoffel_first[0][1] = 0;
	Christoffel_first[1][1] = q[0];
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
	*/
	computeChristoffel_first(q, Christoffel_first);

}

void CylinderCoordiante::computeChristoffel_first(Vector3d Q, double ** Gamma_1st) {
	Gamma_1st[0][0] = 0;
	Gamma_1st[1][0] = 0;
	Gamma_1st[2][0] = 0;
	Gamma_1st[3][0] = -Q[0];
	Gamma_1st[4][0] = 0;
	Gamma_1st[5][0] = 0;
	Gamma_1st[0][1] = 0;
	Gamma_1st[1][1] = Q[0];
	Gamma_1st[2][1] = 0;
	Gamma_1st[3][1] = 0;
	Gamma_1st[4][1] = 0;
	Gamma_1st[5][1] = 0;
	Gamma_1st[0][2] = 0;
	Gamma_1st[1][2] = 0;
	Gamma_1st[2][2] = 0;
	Gamma_1st[3][2] = 0;
	Gamma_1st[4][2] = 0;
	Gamma_1st[5][2] = 0;
}