///////////////////////////////////////////////////////////////////////////////
///
///	\file    VerticalStretch.h
///	\author  Paul Ullrich
///	\version July 13, 2014
///
///	<remarks>
///		Copyright 2000-2010 Paul Ullrich
///
///		This file is distributed as part of the Tempest source code package.
///		Permission is granted to use, copy, modify and distribute this
///		source code and its documentation under the terms of the GNU General
///		Public License.  This software is provided "as is" without express
///		or implied warranty.
///	</remarks>

#ifndef _VERTICALSTRETCH_H_
#define _VERTICALSTRETCH_H_

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A base class defining functions for stretching the vertical
///		coordinate.
///	</summary>
class VerticalStretchFunction {

public:
	///	<summary>
	///		Virtual destructor.
	///	</summary>
	virtual ~VerticalStretchFunction()
	{ }

public:
	///	<summary>
	///		Stretching function, returns the new xi value after stretch and
	///		the derivative of the stretching function at that point.
	///	</summary>
	virtual void operator()(
		double dREta,
		double & dREtaStretch,
		double & dDxREtaStretch
	) = 0;
};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Uniform vertical stretching.
///	</summary>
class VerticalStretchUniform : public VerticalStretchFunction {

public:
	///	<summary>
	///		Stretching function, returns the new xi value after stretch and
	///		the derivative of the stretching function at that point.
	///	</summary>
	virtual void operator()(
		double dREta,
		double & dREtaStretch,
		double & dDxREtaStretch
	) {
		dREtaStretch = dREta;
		dDxREtaStretch = 1.0;
	}
};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Cubic vertical stretching.
///	</summary>
class VerticalStretchCubic : public VerticalStretchFunction {

public:
	///	<summary>
	///		Stretching function, returns the new xi value after stretch and
	///		the derivative of the stretching function at that point.
	///	</summary>
	virtual void operator()(
		double dREta,
		double & dREtaStretch,
		double & dDxREtaStretch
	) {
		const double dConstS1 = 0.5;
		const double dConstS2 = 2.0;

		dREtaStretch =
			dConstS1 * dREta
			+ ( 3.0 - 2.0 * dConstS1 - dConstS2) * dREta * dREta
			+ (-2.0 +       dConstS1 + dConstS2) * dREta * dREta * dREta;

		dDxREtaStretch =
			dConstS1
			+ 2.0 * ( 3.0 - 2.0 * dConstS1 - dConstS2) * dREta
			+ 3.0 * (-2.0 +       dConstS1 + dConstS2) * dREta * dREta;
	}

};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Piecewise linear vertical stretching.
///	</summary>
class VerticalStretchPiecewiseLinear : public VerticalStretchFunction {

public:
	///	<summary>
	///		Stretching function, returns the new xi value after stretch and
	///		the derivative of the stretching function at that point.
	///	</summary>
	virtual void operator()(
		double dREta,
		double & dREtaStretch,
		double & dDxREtaStretch
	) {
		if (dREta < 2.0/3.0) {
			dREtaStretch = 0.5 * dREta;
			dDxREtaStretch = 0.5;
		} else {
			dREtaStretch = 2.0 * (dREta - 2.0/3.0) + 1.0/3.0;
			dDxREtaStretch = 2.0;
		}
	}

};

///////////////////////////////////////////////////////////////////////////////

#endif

