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

void TimestepSchemeStrang::Step(
	bool fFirstStep,
	bool fLastStep,
	const Time & time,
	double dDeltaT
) {
	// Get a copy of the grid
	Grid * pGrid = m_model.GetGrid();

	// Get a copy of the HorizontalDynamics
	HorizontalDynamics * pHorizontalDynamics = m_model.GetHorizontalDynamics();

	// Get a copy of the VerticalDynamics
	VerticalDynamics * pVerticalDynamics = m_model.GetVerticalDynamics();

	// Half a time step
	double dHalfDeltaT = 0.5 * dDeltaT;

	// Vertical timestep
	if (fFirstStep) {
		pGrid->CopyData(0, 1, DataType_State);
		pVerticalDynamics->StepImplicit(0, 1, time, dHalfDeltaT);

	} else {
		DataVector<double> dCarryoverCombination;
		dCarryoverCombination.Initialize(2);
		dCarryoverCombination[0] = 1.0;
		dCarryoverCombination[1] = 1.0;
		pGrid->LinearCombineData(dCarryoverCombination, 0, DataType_State);
	}

	// Explicit fourth-order Runge-Kutta
	if (m_eExplicitDiscretization == RungeKutta4) {
		pGrid->CopyData(0, 1, DataType_State);
		pHorizontalDynamics->StepExplicit(0, 1, time, dHalfDeltaT);
		pVerticalDynamics->StepExplicit(0, 1, time, dHalfDeltaT);

		pGrid->CopyData(0, 2, DataType_State);
		pHorizontalDynamics->StepExplicit(1, 2, time + dHalfDeltaT, dHalfDeltaT);
		pVerticalDynamics->StepExplicit(1, 2, time + dHalfDeltaT, dHalfDeltaT);

		pGrid->CopyData(0, 3, DataType_State);
		pHorizontalDynamics->StepExplicit(2, 3, time + dDeltaT, dDeltaT);
		pVerticalDynamics->StepExplicit(2, 3, time + dDeltaT, dDeltaT);

		DataVector<double> dRK4Combination;
		dRK4Combination.Initialize(5);
		dRK4Combination[0] = - 1.0 / 3.0;
		dRK4Combination[1] = + 1.0 / 3.0;
		dRK4Combination[2] = + 2.0 / 3.0;
		dRK4Combination[3] = + 1.0 / 3.0;
		pGrid->LinearCombineData(dRK4Combination, 4, DataType_State);

		pHorizontalDynamics->StepExplicit(
			3, 4, time + dDeltaT / 6.0, dDeltaT / 6.0);
		pVerticalDynamics->StepExplicit(
			3, 4, time + dDeltaT / 6.0, dDeltaT / 6.0);

	// Explicit strong stability preserving third-order Runge-Kutta
	} else if (m_eExplicitDiscretization == RungeKuttaSSP3) {
		pGrid->CopyData(0, 1, DataType_State);
		pHorizontalDynamics->StepExplicit(0, 1, time, dDeltaT);
		pVerticalDynamics->StepExplicit(0, 1, time, dDeltaT);

		DataVector<double> dSSPRK3CombinationA;
		dSSPRK3CombinationA.Initialize(3);
		dSSPRK3CombinationA[0] = 3.0 / 4.0;
		dSSPRK3CombinationA[1] = 1.0 / 4.0;
		dSSPRK3CombinationA[2] = 0.0;
		pGrid->LinearCombineData(dSSPRK3CombinationA, 2, DataType_State);
		pHorizontalDynamics->StepExplicit(1, 2, time + dDeltaT, 0.25 * dDeltaT);
		pVerticalDynamics->StepExplicit(1, 2, time + dDeltaT, 0.25 * dDeltaT);

		DataVector<double> dSSPRK3CombinationB;
		dSSPRK3CombinationB.Initialize(5);
		dSSPRK3CombinationB[0] = 1.0 / 3.0;
		dSSPRK3CombinationB[1] = 0.0;
		dSSPRK3CombinationB[2] = 2.0 / 3.0;
		dSSPRK3CombinationB[3] = 0.0;
		dSSPRK3CombinationB[4] = 0.0;
		pGrid->LinearCombineData(dSSPRK3CombinationB, 4, DataType_State);
		pHorizontalDynamics->StepExplicit(
			2, 4, time + 0.5 * dDeltaT, (2.0/3.0) * dDeltaT);
		pVerticalDynamics->StepExplicit(
			2, 4, time + 0.5 * dDeltaT, (2.0/3.0) * dDeltaT);

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

	DataVector<double> dOffCenteringCombination;
	dOffCenteringCombination.Initialize(2);
	dOffCenteringCombination[0] = (2.0 - m_dOffCentering) / 2.0;
	dOffCenteringCombination[1] = m_dOffCentering / 2.0;
	pGrid->LinearCombineData(dOffCenteringCombination, 0, DataType_State);

	if (!fLastStep) {
		DataVector<double> dCarryoverFinal;
		dCarryoverFinal.Initialize(2);
		dCarryoverFinal[0] = +1.0;
		dCarryoverFinal[1] = -1.0;
		pGrid->LinearCombineData(dCarryoverFinal, 1, DataType_State);
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

