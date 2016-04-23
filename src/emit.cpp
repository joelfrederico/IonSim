#include "emit.h"
#include "support_func.h"

// ==============================
// Emit
// ==============================
Emit::Emit() 
{
}

void Emit::set_emit(double emit, double E_GeV)
{
	_emit = emit;
}

void Emit::set_emit_n(double emit_n, double E_GeV)
{
	_emit = emit_n / ionsim::GeV2gamma(E_GeV);
}

double Emit::emit() const
{
	return _emit;
}
