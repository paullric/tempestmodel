///////////////////////////////////////////////////////////////////////////////
///
///	\file    TimestepSchemeStrang.cpp
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

#include "TimestepSchemeStrang.h"
#include "Model.h"
#include "Grid.h"
#include "HorizontalDynamics.h"
#include "VerticalDynamics.h"

///////////////////////////////////////////////////////////////////////////////

TimestepSchemeStrang::TimestepSchemeStrang(
	Model & model,
	double dOffCentering,
	ExplicitDiscretization eExplicitDiscretization
) :
	TimestepScheme(model),
	m_dOffCentering(dOffCentering),
	m_eExplicitDiscretization(eExplicitDiscretization)
{
	// Check bounds of OffCentering parameter
	if ((m_dOffCentering < 0.0) || (m_dOffCentering > 1.0)) {
		_EXCEPTIONT("OffCentering parameter out of range [0,1]");
	}

	// Carryover combination
	m_dCarryoverCombination.Initialize(2);
	m_dCarryoverCombination[0] = 1.0;
	m_dCarryoverCombination[1] = 1.0;

	// Off-centering combination
	m_dOffCenteringCombination.Initialize(2);
	m_dOffCenteringCombination[0] = (2.0 - m_dOffCentering) / 2.0;
	m_dOffCenteringCombination[1] = m_dOffCentering / 2.0;

	// Final carryover combination
	m_dCarryoverFinal.Initialize(2);
	m_dCarryoverFinal[0] = +1.0;
	m_dCarryoverFinal[1] = -1.0;

	// RK4 combination
	m_dRK4Combination.Initialize(5);
	m_dRK4Combination[0] = - 1.0 / 3.0;
	m_dRK4Combination[1] = + 1.0 / 3.0;
	m_dRK4Combination[2] = + 2.0 / 3.0;
	m_dRK4Combination[3] = + 1.0 / 3.0;

	// SSPRK3 combination A
	m_dSSPRK3CombinationA.Initialize(3);
	m_dSSPRK3CombinationA[0] = 3.0 / 4.0;
	m_dSSPRK3CombinationA[1] = 1.0 / 4.0;
	m_dSSPRK3CombinationA[2] = 0.0;

	// SSPRK3 combination B
	m_dSSPRK3CombinationB.Initialize(5);
	m_dSSPRK3CombinationB[0] = 1.0 / 3.0;
	m_dSSPRK3CombinationB[1] = 0.0;
	m_dSSPRK3CombinationB[2] = 2.0 / 3.0;
	m_dSSPRK3CombinationB[3] = 0.0;
	m_dSSPRK3CombinationB[4] = 0.0;

	// KGU 3-5 combination
	m_dKinnmarkGrayUllrichCombination.Initialize(5);
	m_dKinnmarkGrayUllrichCombination[0] = - 1.0 / 4.0;
	m_dKinnmarkGrayUllrichCombination[1] =   5.0 / 4.0;
	m_dKinnmarkGrayUllrichCombination[2] =   0.0;
	m_dKinnmarkGrayUllrichCombination[3] =   0.0;
	m_dKinnmarkGrayUllrichCombination[4] =   0.0;
}

///////////////////////////////////////////////////////////////////////////////

void TimestepSchemeStrang::Step(
	bool fFirstStep,
	bool fLastStep,
	const Time & time,
	double dDeltaT
) {
	// Get a copy of the grid
	Grid * pGrid = m_model.GetGrid();

	// Get a copy of the HorizontalDynamics
	HorizontalDynamics * pHorizontalDynamics =
		m_model.GetHorizontalDynamics();

	// Get a copy of the VerticalDynamics
	VerticalDynamics * pVerticalDynamics =
		m_model.GetVerticalDynamics();

	// Half a time step
	double dHalfDeltaT = 0.5 * dDeltaT;

	// Vertical timestep
	if (fFirstStep) {
		pGrid->CopyData(0, 1, DataType_State);
		pVerticalDynamics->StepImplicit(0, 1, time, dHalfDeltaT);

	} else {
		pGrid->LinearCombineData(m_dCarryoverCombination, 0, DataType_State);
	}

	// Explicit fourth-order Runge-Kutta
	if (m_eExplicitDiscretization == RungeKutta4) {
		pGrid->CopyData(0, 1, DataType_State);
		pHorizontalDynamics->StepExplicit(0, 1, time, dHalfDeltaT);
		pVerticalDynamics->StepExplicit(0, 1, time, dHalfDeltaT);
		pGrid->PostProcessSubstage(1, DataType_State);
		pGrid->PostProcessSubstage(1, DataType_Tracers);

		pGrid->CopyData(0, 2, DataType_State);
		pHorizontalDynamics->StepExplicit(1, 2, time, dHalfDeltaT);
		pVerticalDynamics->StepExplicit(1, 2, time, dHalfDeltaT);
		pGrid->PostProcessSubstage(2, DataType_State);
		pGrid->PostProcessSubstage(2, DataType_Tracers);

		pGrid->CopyData(0, 3, DataType_State);
		pHorizontalDynamics->StepExplicit(2, 3, time, dDeltaT);
		pVerticalDynamics->StepExplicit(2, 3, time, dDeltaT);
		pGrid->PostProcessSubstage(3, DataType_State);
		pGrid->PostProcessSubstage(3, DataType_Tracers);

		pGrid->LinearCombineData(m_dRK4Combination, 4, DataType_State);

		pHorizontalDynamics->StepExplicit(3, 4, time, dDeltaT / 6.0);
		pVerticalDynamics->StepExplicit(3, 4, time, dDeltaT / 6.0);
		pGrid->PostProcessSubstage(4, DataType_State);
		pGrid->PostProcessSubstage(4, DataType_Tracers);

	// Explicit strong stability preserving third-order Runge-Kutta
	} else if (m_eExplicitDiscretization == RungeKuttaSSP3) {

		pGrid->CopyData(0, 1, DataType_State);
		pHorizontalDynamics->StepExplicit(0, 1, time, dDeltaT);
		pVerticalDynamics->StepExplicit(0, 1, time, dDeltaT);
		pGrid->PostProcessSubstage(1, DataType_State);
		pGrid->PostProcessSubstage(1, DataType_Tracers);

		pGrid->LinearCombineData(m_dSSPRK3CombinationA, 2, DataType_State);
		pHorizontalDynamics->StepExplicit(1, 2, time, 0.25 * dDeltaT);
		pVerticalDynamics->StepExplicit(1, 2, time, 0.25 * dDeltaT);
		pGrid->PostProcessSubstage(2, DataType_State);
		pGrid->PostProcessSubstage(2, DataType_Tracers);

		pGrid->LinearCombineData(m_dSSPRK3CombinationB, 4, DataType_State);
		pHorizontalDynamics->StepExplicit(2, 4, time, (2.0/3.0) * dDeltaT);
		pVerticalDynamics->StepExplicit(2, 4, time, (2.0/3.0) * dDeltaT);
		pGrid->PostProcessSubstage(4, DataType_State);
		pGrid->PostProcessSubstage(4, DataType_Tracers);

	// Explicit Kinnmark, Gray and Ullrich third-order five-stage Runge-Kutta
	} else if (m_eExplicitDiscretization == KinnmarkGrayUllrich35) {

		pGrid->CopyData(0, 1, DataType_State);
		pHorizontalDynamics->StepExplicit(0, 1, time, dDeltaT / 5.0);
		pVerticalDynamics->StepExplicit(0, 1, time, dDeltaT / 5.0);
		pGrid->PostProcessSubstage(1, DataType_State);
		pGrid->PostProcessSubstage(1, DataType_Tracers);

		pGrid->CopyData(0, 2, DataType_State);
		pHorizontalDynamics->StepExplicit(1, 2, time, dDeltaT / 5.0);
		pVerticalDynamics->StepExplicit(1, 2, time, dDeltaT / 5.0);
		pGrid->PostProcessSubstage(2, DataType_State);
		pGrid->PostProcessSubstage(2, DataType_Tracers);

		pGrid->CopyData(0, 3, DataType_State);
		pHorizontalDynamics->StepExplicit(2, 3, time, dDeltaT / 3.0);
		pVerticalDynamics->StepExplicit(2, 3, time, dDeltaT / 3.0);
		pGrid->PostProcessSubstage(3, DataType_State);
		pGrid->PostProcessSubstage(3, DataType_Tracers);

		pGrid->CopyData(0, 2, DataType_State);
		pHorizontalDynamics->StepExplicit(3, 2, time, 2.0 * dDeltaT / 3.0);
		pVerticalDynamics->StepExplicit(3, 2, time, 2.0 * dDeltaT / 3.0);
		pGrid->PostProcessSubstage(2, DataType_State);
		pGrid->PostProcessSubstage(2, DataType_Tracers);

		pGrid->LinearCombineData(
			m_dKinnmarkGrayUllrichCombination, 4, DataType_State);
		pHorizontalDynamics->StepExplicit(2, 4, time, 3.0 * dDeltaT / 4.0);
		pVerticalDynamics->StepExplicit(2, 4, time, 3.0 * dDeltaT / 4.0);
		pGrid->PostProcessSubstage(4, DataType_State);
		pGrid->PostProcessSubstage(4, DataType_Tracers);

	// Invalid explicit discretization
	} else {
		_EXCEPTIONT("Invalid explicit discretization");
	}

	// Apply hyperdiffusion
	pGrid->CopyData(4, 1, DataType_State);
	pHorizontalDynamics->StepAfterSubCycle(4, 1, 2, time, dDeltaT);

	// Vertical timestep
	double dOffCenterDeltaT = 0.5 * (1.0 + m_dOffCentering) * dDeltaT;

	pGrid->CopyData(1, 0, DataType_State);
	pVerticalDynamics->StepImplicit(0, 0, time, dOffCenterDeltaT);

	pGrid->LinearCombineData(m_dOffCenteringCombination, 0, DataType_State);

	if (!fLastStep) {
		pGrid->LinearCombineData(m_dCarryoverFinal, 1, DataType_State);
	}

	//pGrid->CopyData(0, 1, DataType_State);
	//pVerticalDynamics->StepExplicit(0, 1, time, dDeltaT);
	//pVerticalDynamics->StepImplicit(0, 0, time, dDeltaT);
/*
	if (0) {
	//if (fLastStep) {
		DataVector<double> dDifference;
		dDifference.Initialize(2);
		dDifference[0] = -1.0;
		dDifference[1] = 1.0;
		pGrid->LinearCombineData(dDifference, 0, DataType_State);
	} else {
		pGrid->CopyData(4, 0, DataType_State);
	}
*/
}

///////////////////////////////////////////////////////////////////////////////

