///////////////////////////////////////////////////////////////////////////////
///
///	\file    TimestepSchemeARS222.h
///	\author  Jorge Guerra
///	\version January 29, 2016
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

#ifndef _TIMESTEPSCHEMEARS222_H_
#define _TIMESTEPSCHEMEARS222_H_

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
class TimestepSchemeARS222 : public TimestepScheme {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	TimestepSchemeARS222(
		Model & model
	);

public:
	///	<summary>
	///		Get the number of component data instances.
	///	</summary>
	virtual int GetComponentDataInstances() const {
		return 4;
	}

	///	<summary>
	///		Get the number of tracer data instances.
	///	</summary>
	virtual int GetTracerDataInstances() const {
		return 4;
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
	///		ARS222 parameter gamma
	///	</summary>
	static const double m_dgamma;
	///	<summary>
	///		ARS222 parameter delta
	///	</summary>
	static const double m_ddelta;

	///	<summary>
	///		Coefficients for the time increment ARS222.
	///	</summary>
	static const double m_dTimeCf[2];

	///	<summary>
	///		Coefficients for the explicit ARS222.
	///	</summary>
	static const double m_dExpCf[2][2];

	///	<summary>
	///		Coefficients for the explicit ARS222.
	///	</summary>
	static const double m_dImpCf[2][2];

	///	<summary>
	///		Explicit evaluation at the 2nd substage
	///	</summary>
	DataArray1D<double> m_du2fCombo;
};

///////////////////////////////////////////////////////////////////////////////

#endif


