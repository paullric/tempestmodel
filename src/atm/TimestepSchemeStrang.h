///////////////////////////////////////////////////////////////////////////////
///
///	\file    TimestepSchemeStrang.h
///	\author  Paul Ullrich
///	\version June 18, 2013
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

#ifndef _TIMESTEPSCHEMESTRANG_H_
#define _TIMESTEPSCHEMESTRANG_H_

#include "TimestepScheme.h"
#include "Exception.h"

///////////////////////////////////////////////////////////////////////////////

class Model;
class Time;

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Adaptive Fourth-Order Runge-Kutta time stepping.
///	</summary>
class TimestepSchemeStrang : public TimestepScheme {

public:
	///	<summary>
	///		Explicit time discretization.
	///	</sumamry>
	enum ExplicitDiscretization {
		RungeKutta4,
		RungeKuttaSSP3
	};

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	TimestepSchemeStrang(
		Model & model,
		double dOffCentering = 0.0,
		ExplicitDiscretization eExplicitDiscretization = RungeKuttaSSP3
	) :
		TimestepScheme(model),
		m_dOffCentering(dOffCentering),
		m_eExplicitDiscretization(eExplicitDiscretization)
	{
		// Check bounds of OffCentering parameter
		if ((m_dOffCentering < 0.0) || (m_dOffCentering > 1.0)) {
			_EXCEPTIONT("OffCentering parameter out of range [0,1]");
		}
	}

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
	///		Off-centering parameter.
	///	</summary>
	double m_dOffCentering;

	///	<summary>
	///		Explicit time discretization to use.
	///	</summary>
	ExplicitDiscretization m_eExplicitDiscretization;
};

///////////////////////////////////////////////////////////////////////////////

#endif

