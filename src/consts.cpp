#include "consts.h"

std::string _pushname(const pushmethod_t PUSH_TYPE)
{
	std::string name;

	switch (PUSH_TYPE)
	{
		case PUSH_RUNGE_KUTTA:
			name = "Runge-Kutta";
			break;
		case PUSH_SIMPLE:
			name = "Simple";
			break;
		case PUSH_FIELD:
			name = "Field";
			break;
	}

	return name;
}

PushMethod::PushMethod(const pushmethod_t PUSH_TYPE) :
	_PUSH_TYPE(PUSH_TYPE),
	name(_pushname(PUSH_TYPE))
{

}

bool PushMethod::operator==(const PushMethod &other) const
{
	return (this->_PUSH_TYPE == other._PUSH_TYPE);
}

std::string _zdistname(const zdist_t zdist)
{
	std::string name;

	switch (zdist)
	{
		case Z_DIST_FLAT:
			name = "Flat";
			break;
		case Z_DIST_GAUSS:
			name = "Gauss";
			break;
	}

	return name;
}

zDist::zDist(const zdist_t Z_DIST) :
	_Z_DIST(Z_DIST),
	name(_zdistname(Z_DIST))
{
}
