///////////////////////////////////////////////////////////////////////////////
///
///	\file    TimestepSchemeARS343b.h
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

#ifndef _TIMESTEPSCHEMEARS343b_H_
#define _TIMESTEPSCHEMEARS343b_H_

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
class TimestepSchemeARS343b : public TimestepScheme {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	TimestepSchemeARS343b(
		Model & model
	);

public:
	///	<summary>
	///		Get the number of component data instances.
	///	</summary>
	virtual int GetComponentDataInstances() const {
		return 7;
	}

	///	<summary>
	///		Get the number of tracer data instances.
	///	</summary>
	virtual int GetTracerDataInstances() const {
		return 7;
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
	///		Data instances for 2nd stage combination.
	///	</summary>
	DataArray1D<int> m_iU2fCombo;

	///	<summary>
	///		Coefficients for 2nd stage combination.
	///	</summary>
	DataArray1D<double> m_dU2fCombo;

	///	<summary>
	///		Data instances for 3rd stage combination.
	///	</summary>
	DataArray1D<int> m_iU3fCombo;

	///	<summary>
	///		Coefficients for 3rd stage combination.
	///	</summary>
	DataArray1D<double> m_dU3fCombo;

	///	<summary>
	///		Data instances for 4th stage combination.
	///	</summary>
	DataArray1D<int> m_iU4fCombo;

	///	<summary>
	///		Coefficients for 4th stage combination.
	///	</summary>
	DataArray1D<double> m_dU4fCombo;

	///	<summary>
	///		Data instances for copying data from individual data instances.
	///	</summary>
	DataArray1D<int> m_iCopyInst[7];

	///	<summary>
	///		Coefficients for copying data from individual data instances.
	///	</summary>
	DataArray1D<double> m_dCopyInst[7];
};

///////////////////////////////////////////////////////////////////////////////

#endif


