#ifndef __BASECLASS_H_INCLUDED__
#define __BASECLASS_H_INCLUDED__

class Field;

class Field
{
	private:
		double *_x_data;
		double *_y_data;
		long _x_pts;
		long _y_pts;
		long _n_pts;
		int _init(long x_pts, long y_pts);
		int _copy(const Field &rhs);
		bool _samedim(const Field &rhs);

	public:
		Field(long x_pts, long y_pts);
		Field(const Field &rhs);
		~Field();

		Field &operator=(const Field &rhs);
		void operator()(long i, long j);
		Field &operator+=(const Field &rhs);
		Field &operator-=(const Field &rhs);
		template <class T>
		Field &operator*=(const T rhs)
		{
			if ( (*this)._samedim(rhs) )
			{
				*(*this)._x_data *= rhs;
				*(*this)._y_data *= rhs;
			} else {
				throw "Cannot subtract fields of different sizes";
			}
			return *this;
		}
		template <class T>
		Field &operator/=(const T rhs)
		{
			if ( (*this)._samedim(rhs) )
			{
				*(*this)._x_data /= rhs;
				*(*this)._y_data /= rhs;
			} else {
				throw "Cannot subtract fields of different sizes";
			}
			return *this;
		}

		const Field operator+(const Field &rhs);
		const Field operator-(const Field &rhs);

		template <class T>
		const Field operator*(T rhs)
		{
			return Field(*this) *= rhs;
		}
		template <class T>
		const Field operator/(T rhs)
		{
			return Field(*this) /= rhs;
		}
};

#endif
