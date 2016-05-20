#include "consts.h"

PushMethod::PushMethod(const pushmethod_t PUSH_TYPE) :
	_PUSH_TYPE(PUSH_TYPE)
{
	switch (PUSH_TYPE)
	{
		case 1:
			name = "Runge-Kutta";
			break;
		case 2:
			name = "Simple";
			break;
		case 3:
			name = "Field";
			break;
	}

}

bool PushMethod::operator==(const PushMethod &other) const
{
	return (this->_PUSH_TYPE == other._PUSH_TYPE);
}
