#ifndef CELESTIALOBJECT_H__
#define CELESTIALOBJECT_H__

#include <Eigen/Dense>
#include "CoordinateSystem.h"
#include "EuclidCoordinate.h"
#include "CylinderCoordinate.h"
#include "SphereCoordinate.h"


namespace CoordType {
	unsigned int EUCLID = 0;
	unsigned int CYLIND = 1;
	unsigned int SPHERE = 2;
}


class CelestialObject {
	protected:

		// The formula for this is general given the
		// qdot, Christoffel symbol of 2nd kind, total force, contravariant vectors, and mass 
		// TODO: How to incoporate constraints... handle in propagate?
		virtual void computeQddot(int num);

		Eigen::Vector3d qddot;

		CoordinateSystem coordSys;

		double m;

		Eigen::Vector3d totalForce;

	public:

		CelestialObject(unsigned int CoordSys);

		virtual ~CelestialObject() {};

		virtual void addForce(const Eigen::Vector3d & force);

		virtual void propagate() = 0;


}



#end