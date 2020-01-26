#ifndef EUCLIDCOORDINATE_H__
#define EUCLIDCOORDINATE_H__

#include "CoordinateSystem.h"

class EuclidCoordinate : public Coordinate {
	public:
		EuclidCoordinate();

	protected:
		double q[3];
		double qdot[3];

		Vector3d r, rdot;
		Vector3d covariant[3];
		double Christoffel_first[6][3];

		// Defines r(q1, q2, q3)
		virtual void computeR();

		// Defines rdot(q1, q2, q3, q1dot, q2dot, q3dot)
		virtual void computeRDot();

		// Defines a1(q1,q2,q3), a2(q1,q2,q3), a3(q1,q2,q3)
		virtual void computeCovariant();

		// Defines Christofell Symbols of the 1st kind
		virtual void computeChristoffel_first();

	
	public:

		virtual ~CoordinateSystem() {}

		virtual void setQ (double x, double y, double z);
		virtual void setQDot (double xdot, double ydot, double zdot);

		virtual const double * getQ();
		virtual const double * getQDot();

		virtual const Vector3d & getR();
		virtual const Vector3d & getRDot();
		virtual const Vector3d * getCovariant();
		virtual double getChristoffel_first(int i, int j, int k);


}


#endif