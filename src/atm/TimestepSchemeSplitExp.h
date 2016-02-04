///////////////////////////////////////////////////////////////////////////////
///
///	\file    TimestepSchemeSplitExp.h
///	\author  Paul Ullrich, Jorge Guerra
///	\version January 11, 2016
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

#ifndef _TIMESTEPSCHEMESplitExp_H_
#define _TIMESTEPSCHEMESplitExp_H_

#include "TimestepScheme.h"
#include "Exception.h"
#include "DataArray1D.h"

///////////////////////////////////////////////////////////////////////////////

class Model;
class Time;

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Multistep Runge-Kutta 3 (Split Explicit) time stepping.
///	</summary>
class TimestepSchemeSplitExp : public TimestepScheme {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	TimestepSchemeSplitExp(
		Model & model
	);

public:
	///	<summary>
	///		Get the number of component data instances.
	///	</summary>
	virtual int GetComponentDataInstances() const {
		return 6;
	}

	///	<summary>
	///		Get the number of tracer data instances.
	///	</summary>
	virtual int GetTracerDataInstances() const {
		return 6;
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
	///		Linear combination coefficients used by KGU53.
	///	</summary>
	DataArray1D<double> m_dKinnmarkGrayUllrichCombination;

	///	<summary>
	///		Linear combination coefficients used by SSPRK3 (combination A).
	///	</summary>
	DataArray1D<double> m_dSSPRK3CombinationA;

	///	<summary>
	///		Linear combination coefficients used by SSPRK3 (combination B).
	///	</summary>
	DataArray1D<double> m_dSSPRK3CombinationB;
};

///////////////////////////////////////////////////////////////////////////////

#endif


