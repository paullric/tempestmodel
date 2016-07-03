///////////////////////////////////////////////////////////////////////////////
///
///	\file    TimestepSchemeSSP3332.h
///	\author  Jorge Guerra
///	\version January 27, 2016
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

#ifndef _TIMESTEPSCHEMESSP3332_H_
#define _TIMESTEPSCHEMESSP3332_H_

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
class TimestepSchemeSSP3332 : public TimestepScheme {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	TimestepSchemeSSP3332(
		Model & model
	);

public:
	///	<summary>
	///		Get the number of component data instances.
	///	</summary>
	virtual int GetComponentDataInstances() const {
		return 9;
	}

	///	<summary>
	///		Get the number of tracer data instances.
	///	</summary>
	virtual int GetTracerDataInstances() const {
		return 9;
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
	///		SSP3(332) parameter gamma
	///	</summary>
	static const double m_dgamma;

	///	<summary>H
	///		Coefficients for the time increment SSP3(332).
	///	</summary>
	static const double m_dTimeCf[4];

	///	<summary>
	///		Coefficients for the explicit SSP3(332).
	///	</summary>
	static const double m_dExpCf[4][4];

	///	<summary>
	///		Coefficients for the explicit SSP3(332).
	///	</summary>
	static const double m_dImpCf[4][4];

	///		Explicit evaluation at the 2nd substage
	///	</summary>
	DataArray1D<double> m_du2fCombo;

	///	<summary>
	///		Explicit evaluation at the 3rd substage
	///	</summary>
	DataArray1D<double> m_du3fCombo;

	///	<summary>
	///		Explicit evaluation at the 4th substage
	///	</summary>
	DataArray1D<double> m_du4fCombo;
};

///////////////////////////////////////////////////////////////////////////////

#endif


