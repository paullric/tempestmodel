///////////////////////////////////////////////////////////////////////////////
///
///	\file    TimestepSchemeARS232.h
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

#ifndef _TIMESTEPSCHEMEARS232_H_
#define _TIMESTEPSCHEMEARS232_H_

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
class TimestepSchemeARS232 : public TimestepScheme {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	TimestepSchemeARS232(
		Model & model
	);

public:
	///	<summary>
	///		Get the number of component data instances.
	///	</summary>
	virtual int GetComponentDataInstances() const {
		return 5;
	}

	///	<summary>
	///		Get the number of tracer data instances.
	///	</summary>
	virtual int GetTracerDataInstances() const {
		return 5;
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
	///		Diagonal explicit coefficients.
	///	</summary>
	DataArray1D<double> m_dDiagExpCf;

	///	<summary>
	///		Diagonal implicit coefficients.
	///	</summary>
	DataArray1D<double> m_dDiagImpCf;

	///	<summary>
	///		Linear combination at the 2nd substage.
	///	</summary>
	DataArray1D<double> m_dU2fCombo;

	///	<summary>
	///		Linear combination at the 3rd substage.
	///	</summary>
	DataArray1D<double> m_dU3fCombo;
};

///////////////////////////////////////////////////////////////////////////////

#endif
