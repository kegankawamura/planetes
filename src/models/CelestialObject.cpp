#include "CelestialObject.h"

CelestialObject::CelestialObject(unsigned int CoordSys) {
	switch(CoordSys) {
		case CoordType::EUCLID:

			coordSys = EuclidCoordinate();
			
			break;
		case CoordType::CYLIND:

			coordSys = CylinderCoordinate();

			break;

		case CoordType::SPHERE:

			coordSys = SphereCoordinate();

			break;

	}
}