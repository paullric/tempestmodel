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
#include "DataArray1D.h"

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
		ForwardEuler,
		RungeKutta4,
		RungeKuttaSSP3,
		KinnmarkGrayUllrich35,
		RungeKuttaSSPRK53
	};

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	TimestepSchemeStrang(
		Model & model,
		double dOffCentering = 0.0,
		ExplicitDiscretization eExplicitDiscretization = KinnmarkGrayUllrich35
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

	///	<summary>
	///		Get the number of substeps.
	///	</summary>
	virtual int GetSubStepCount() const;

	///	<summary>
	///		Get the maximum stable Courant number for the explicit part of the
	///		Timesteps scheme.
	///	</summary>
	virtual double GetMaximumStableCourantNumber(
		TimestepScheme::MixedMethodPart eMixedMethodPart,
		int nOrder
	) const;

protected:
	///	<summary>
	///		Perform one time substep.
	///	</summary>
	virtual int SubStep(
		bool fFirstStep,
		bool fLastStep,
		const Time & time,
		double dDeltaT,
		int iSubStep
	);

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

	///	<summary>
	///		Number of sub-steps in explicit discretization.
	///	</summary>
	int m_nExplicitSubSteps;

	///	<summary>
	///		Carryover combination.
	///	</summary>
	DataArray1D<double> m_dCarryoverCombination;

	///	<summary>
	///		Off-centering combination.
	///	</summary>
	DataArray1D<double> m_dOffCenteringCombination;

	///	<summary>
	///		Final carryover combination.
	///	</summary>
	DataArray1D<double> m_dCarryoverFinal;

	///	<summary>
	///		Linear combination coefficients used by RK4.
	///	</summary>
	DataArray1D<double> m_dRK4Combination;

	///	<summary>
	///		Linear combination coefficients used by SSPRK3 (combination A).
	///	</summary>
	DataArray1D<double> m_dSSPRK3CombinationA;

	///	<summary>
	///		Linear combination coefficients used by SSPRK3 (combination B).
	///	</summary>
	DataArray1D<double> m_dSSPRK3CombinationB;

	///	<summary>
	///		Linear combination coefficients used by KGU35.
	///	</summary>
	DataArray1D<double> m_dKinnmarkGrayUllrichCombination;

	///	<summary>
	///		Linear combination coefficients used by SSPRK53 (combination A).
	///	</summary>
	DataArray1D<double> m_dSSPRK53CombinationA;

	///	<summary>
	///		Linear combination coefficients used by SSPRK53 (combination B).
	///	</summary>
	DataArray1D<double> m_dSSPRK53CombinationB;

	///	<summary>
	///		Linear combination coefficients used by SSPRK53 (combination C).
	///	</summary>
	DataArray1D<double> m_dSSPRK53CombinationC;

};

///////////////////////////////////////////////////////////////////////////////

#endif

