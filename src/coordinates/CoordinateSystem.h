#ifndef COORDINATESYSTEM_H__
#define COORDINATESYSTEM_H__

#include <math.h>
#include <Eigen/Dense>


using namespace Eigen;


class CoordinateSystem {
	protected:
		double q[3];
		double qdot[3];

		Vector3d r, rdot;
		Vector3d covariant[3];

		// Since [ij,k] = [ji,k], six combinations
		double Christoffel_first[6][3];

		// Defines r(q1, q2, q3)
		virtual void computeR() = 0;

		// Defines rdot(q1, q2, q3, q1dot, q2dot, q3dot)
		virtual void computeRDot() = 0;

		// Defines a1(q1,q2,q3), a2(q1,q2,q3), a3(q1,q2,q3)
		virtual void computeCovariant() = 0;


		// Defines Christofell Symbols of the 1st kind
		virtual void computeChristoffel_first() = 0;
		virtual void computeChristoffel_first(Vector3d Q, Vector3d * cov) = 0;

	
	public:

		virtual ~CoordinateSystem() {};

		virtual void setQ (double, double, double) = 0;
		virtual void setQDot (double, double, double) = 0;

		virtual const double * getQ {
			return q;
		}
		virtual const double * getQDot() {
			return qdot;
		}

		virtual const Vector3d & getR() {
			return r;
		}
		virtual const Vector3d & getRDot() {
			return rdot;
		}
		virtual const Vector3d * getCovariant() {
			return covariant;
		}

		virtual double getChristoffel_first(int i, int j, int k) {
			if (i > 3 || j > 3 || k > 3) {
				return 0;
			}

			if (i > j) {
				int tmp = i;
				i = j;
				j = tmp;
			}

			if (i == 1) {
				return Christoffel_first[j-i][k-1];
			} else if (i == 2) {
				return Christoffel_first[3 + j - i][k-1];
			} else if (i == 3) {
				return Christoffel_first[3 + 2 + j - i][k-1];
			}
			return 0;
		}

		virtual void computeCovariant(Vector3d Q, Vector3d * cov) = 0;

		virtual void computeChristoffel_first(Vector3d Q, Vector3d * cov) = 0;
		

}

#endif