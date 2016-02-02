///////////////////////////////////////////////////////////////////////////////
///
///	\file    TimestepSchemeARK232.h
///	\author  Paul Ullrich
///	\version April 22, 2014
///
///	<remarks>
///		Copyright 2000-2010 Paul Ullrich, Jorge Guerra
///
///		This file is distributed as part of the Tempest source code package.
///		Permission is granted to use, copy, modify and distribute this
///		source code and its documentation under the terms of the GNU General
///		Public License.  This software is provided "as is" without express
///		or implied warranty.
///	</remarks>

#ifndef _TIMESTEPSCHEMEARK232_H_
#define _TIMESTEPSCHEMEARK232_H_

#include "TimestepScheme.h"
#include "Exception.h"
#include "DataArray1D.h"

///////////////////////////////////////////////////////////////////////////////

class Model;
class Time;

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Adaptive Fourth-Order Runge-Kutta time stepping.
///	</summary>
class TimestepSchemeARK232 : public TimestepScheme {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	TimestepSchemeARK232(
		Model & model
	);

public:
	///	<summary>
	///		Get the number of component data instances.
	///	</summary>
	virtual int GetComponentDataInstances() const {
		return 8;
	}

	///	<summary>
	///		Get the number of tracer data instances.
	///	</summary>
	virtual int GetTracerDataInstances() const {
		return 8;
	}

protected:
	///	<summary>
	///		Perform one time step.
	///	</summary>
	virtual void Step(
		bool fFirstStep,
		bool fLastStep,
		const Time & time,
		double dDeltaT
	);

private:
        ///	<summary>
	///		ARK232 parameter gamma
	///	</summary>
	static const double m_dgamma;

        ///	<summary>
	///		ARK232 parameter alpha
	///	</summary>
	static const double m_dalpha;

        ///	<summary>
	///		ARK232 parameter delta
	///	</summary>
	static const double m_ddelta;

        ///	<summary>
	///		Coefficients for the explicit ARK232.
	///	</summary>
	static const double m_dExpCf[3][3];

	///	<summary>
	///		Coefficients for the explicit ARK232.
	///	</summary>
	static const double m_dImpCf[3][3];

	///		Explicit evaluation at the 2nd substage
	///	</summary>
	DataArray1D<double> m_du2fCombo;

	///	<summary>
	///		Explicit evaluation at the 3rd substage
	///	</summary>
	DataArray1D<double> m_du3fCombo;
};

///////////////////////////////////////////////////////////////////////////////

#endif


