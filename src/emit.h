#ifndef __EMIT_H_INCLUDED__
#define __EMIT_H_INCLUDED__

class Emit
{
	private:
		double _emit;
	public:
		Emit();

		void set_emit_n(double emit_n, double E_GeV);
		void set_emit(double emit, double E_GeV);

		double emit();
};

#endif
