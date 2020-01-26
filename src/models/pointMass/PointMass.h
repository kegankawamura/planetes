#ifndef POINTMASS_H__
#define POINTMASS_H__

#include <Eigen/Dense>
#include "CoordinateSystem.h"
#include "tools/NumIntegrators.h"


class PointMass : public CelestialObject {
	protected:

		// The formula for this is general given the
		// qdot, Christoffel symbol of 2nd kind, total force, contravariant vectors, and mass 
		// TODO: How to incoporate constraints... handle in propagate?
		virtual void computeQddot();

		virtual void computeCovariantMatrix(Eigen::Matrix3d & covMatrix);



		Eigen::Vector3d qddot;

		CoordinateSystem coordSys;

		double m;

		Eigen::Vector3d totalForce;



	public:

		PointMass(unsigned int CoordSys);

		virtual ~PointMass() {};

		virtual void addForce(const Eigen::Vector3d & force);

		virtual void propagate(double del_t);


	private:

		void computeForceComponents(Eigen::Vector3d & forceComp, const Eigen3d::Vector3d * covariant);
		void computeChristoffelSum(Eigen::Vector3d & christoffelSum);



}



#end