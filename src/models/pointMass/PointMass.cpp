#include "PointMass.h"

PointMass::PointMass(unsigned int CoordSys) : CelestialObject(CoordSys) {
	
}


void PointMass::computeQddot(Vector3d Q, Vector3d Qdot, Vector3d & Qddot) {

	Eigen3d::Vector3d * covariant;

	coordSys.computeCovariant(Q, covariant);

	Eigen::Matrix3d covMatrix;
	computeCovariantMatrix(covMatrix, covariant);

	Eigen::Vector3d forceComp;
	computeForceComponents(forceComp, covariant);

	Eigen::Vector3d christoffelSum;
	computeChristoffelSum(christoffelSum);

	// Solve Ax = b for x
	// A = covMatrix
	// b = 1/m*forceComp - christoffelSum
	// LDLT method used for speed and accuracy

	Qddot = covMatrix.ldlt().solve(1/m*forceComp - christoffelSum);
}


void PointMass::computeQddot() {
	Eigen3d::Vector3d * covariant = coordSys.getCovariant();

	Eigen::Matrix3d covMatrix;
	computeCovariantMatrix(covMatrix, covariant);

	Eigen::Vector3d forceComp;
	computeForceComponents(forceComp, covariant);

	Eigen::Vector3d christoffelSum;
	computeChristoffelSum(christoffelSum);

	// Solve Ax = b for x
	// A = covMatrix
	// b = 1/m*forceComp - christoffelSum
	// LDLT method used for speed and accuracy

	qddot = covMatrix.ldlt().solve(1/m*forceComp - christoffelSum);
}

void PointMass::computeCovariantMatrix(Eigen::Matrix3d & covMatrix, const Eigen3d::Vector3d * covariant) {
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j <= i; j++) {
			covMatrix(i, j) = covariant[i].dot(covariant[j]);
			if (i != j) {
				covMatrix(j, i) = covMatrix(i, j);
			}
		}
	}
}

void PointMass::computeForceComponents(Eigen::Vector3d & forceComp, const Eigen3d::Vector3d * covariant) {
	for (int i = 0; i < 3; i++) {
		forceComp[i] = totalForce.dot(covariant[i]);
	}
}

void PointMass::computeChristoffelSum(Eigen::Vector3d & christoffelSum) {
	double christoffel_firsts[18];

	for (int k = 1; k <= 3; k++) {
		for (int i = 1; i <= 3; i++) {
			for (int j = i; j <= 3; j++) {
				if (i == 1) {
					christoffel_firsts[6*(k-1) + j-i] = coordSys.getChristoffel_first(i, j, k);
				} else if (i == 2) {
					christoffel_firsts[6*(k-1) + 3 + j - i] = coordSys.getChristoffel_first(i, j, k);
				} else if (i == 3) {
					christoffel_firsts[6*(k-1) + 3 + 2 + j - i] = coordSys.getChristoffel_first(i, j, k);
				}
			}
		}
	}

	Eigen::Matrix3d mat_1 << 	christoffel_firsts[0], christoffel_firsts[1], christoffel_firsts[2],
								christoffel_firsts[6], christoffel_firsts[7], christoffel_firsts[8],
								christoffel_firsts[12],christoffel_firsts[13],christoffel_firsts[14];

	Eigen::Matrix3d mat_2 << 	christoffel_firsts[1], christoffel_firsts[3], christoffel_firsts[4],
								christoffel_firsts[7], christoffel_firsts[9], christoffel_firsts[10],
								christoffel_firsts[13],christoffel_firsts[15],christoffel_firsts[16];

	Eigen::Matrix3d mat_3 << 	christoffel_firsts[2], christoffel_firsts[4], christoffel_firsts[5],
								christoffel_firsts[8], christoffel_firsts[10], christoffel_firsts[11],
								christoffel_firsts[14],christoffel_firsts[16],christoffel_firsts[17];

	double * qdot = coordSys.getQdot();

	Eigen::Matrix3d mat = mat_1 * qdot[0] + mat_2 * qdot[1] + mat_3 * qdot[2];

	Vector3d QDot(qdot[0], qdot[1], qdot[2]);

	christoffelSum = mat * QDot;

}

void PointMass::addForce(const Eigen::Vector3d & force) {
	totalForce += force;
}

void PointMass::propagate(double del_t) {
	
	Vector3d Qd_1, Qdd_1, Qd_2, Qdd_2, Qd_3, Qdd_3, Qd_4, Qdd_4;


	double * q = coordSys.getQ();
	double * qdot = coordSys.getQdot();
	Vector3d Q(q[0], q[1], q[2]);
	Vector3d QDot(qdot[0], qdot[1], qdot[2]);



	Qd_1 = QDot;
	computeQddot(Q, QDot, Qdd_1);

	Qd_2 = QDot + Qdd_1 * del_t/2.0;
	computeQddot(Q + Qd_1 * del_t/2.0, QDot + Qdd_1 * del_t/2.0, Qdd_2);

	Qd_3 = QDot + Qdd_2 * del_t/2.0;
	computeQddot(Q + Qd_2 * del_t/2.0, QDot + Qdd_2 * del_t/2.0, Qdd_3);

	Qd_4 = QDot + Qdd_3 * del_t/2.0;
	computeQddot(Q + Qd_3 * del_t, QDot + Qdd_3 * del_t, Qdd_4);

	// NEED TO APPLY CONSTRAINTS
	// Constraints defined by Psi(r,t) = 0, or f*v + e = 0
	// With this setup, our inputs become r and t, or f, v and e
	// These are functions of q, qdot, and coordinate vectors
	// Ideally each constraint is restricted to one component of q each.


	coordSys.setQ(Q + (Qd_1/6.0 + Qd_2/3.0 + Qd_3/3.0 + Qd_4/6.0) * del_t);
	coordSys.setQdot(QDot + (Qd_1/6.0 + Qd_2/3.0 + Qd_3/3.0 + Qd_4/6.0) * del_t);

/**
	X k1, k2, k3, k4;
	k1 = f(x, t);
	k2 = f(x + k1*dt/2.0, t + dt/2.0);
	k3 = f(x + k2*dt/2.0, t + dt/2.0);
	k4 = f(x + k3*dt, t + dt);
	return x + (k1/6.0 + k2/3.0 + k3/3.0 + k4/6.0)*dt; 
*/

}