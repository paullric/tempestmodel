///////////////////////////////////////////////////////////////////////////////
///
///	\file	TimestepSchemeSplitExp.cpp
///	\author  Paul Ullrich, Jorge Guerra
///	\version January 11, 2016
///
///	<remarks>
///		Copyright 2000-2010 Paul Ullrich and Jorge Guerra
///
///		This file is distributed as part of the Tempest source code package.
///		Permission is granted to use, copy, modify and distribute this
///		source code and its documentation under the terms of the GNU General
///		Public License.  This software is provided "as is" without express
///		or implied warranty.
///	</remarks>

#include "TimestepSchemeSplitExp.h"
#include "Model.h"
#include "Grid.h"
#include "HorizontalDynamics.h"
#include "VerticalDynamicsFEM.h"

///////////////////////////////////////////////////////////////////////////////
// Implements that time split RK3 scheme of Skamarock 2002
TimestepSchemeSplitExp::TimestepSchemeSplitExp(
	Model & model
) :
	TimestepScheme(model)
{
	// Allocate memory for the function evaluations
	// KGU 3-5 combination
	m_dKinnmarkGrayUllrichCombination.Allocate(5);
	m_dKinnmarkGrayUllrichCombination[0] = - 1.0 / 4.0;
	m_dKinnmarkGrayUllrichCombination[1] =   5.0 / 4.0;
	m_dKinnmarkGrayUllrichCombination[2] =   0.0;
	m_dKinnmarkGrayUllrichCombination[3] =   0.0;
	m_dKinnmarkGrayUllrichCombination[4] =   0.0;

	// SSPRK3 combination A
	m_dSSPRK3CombinationA.Allocate(3);
	m_dSSPRK3CombinationA[0] = 3.0 / 4.0;
	m_dSSPRK3CombinationA[1] = 1.0 / 4.0;
	m_dSSPRK3CombinationA[2] = 0.0;

	// SSPRK3 combination B
	m_dSSPRK3CombinationB.Allocate(5);
	m_dSSPRK3CombinationB[0] = 1.0 / 3.0;
	m_dSSPRK3CombinationB[1] = 0.0;
	m_dSSPRK3CombinationB[2] = 2.0 / 3.0;
	m_dSSPRK3CombinationB[3] = 0.0;
	m_dSSPRK3CombinationB[4] = 0.0;
}

///////////////////////////////////////////////////////////////////////////////

void TimestepSchemeSplitExp::Step(
	bool fFirstStep,
	bool fLastStep,
	const Time & time,
	double dDeltaT
) {
	// Get a copy of the grid
	Grid * pGrid = m_model.GetGrid();

	// Pointer to grid
	//GridGLL * pGridGLL = dynamic_cast<GridGLL *>(m_model.GetGrid());
	//if (pGridGLL == NULL) {
	//	_EXCEPTIONT("Invalid grid in Split Explicit -- expected GridGLL");
	//}

	// Get a copy of the HorizontalDynamics
	HorizontalDynamics * pHorizontalDynamics = m_model.GetHorizontalDynamics();

	// Get a copy of the VerticalDynamics
	VerticalDynamicsFEM * pVerticalDynamics = (VerticalDynamicsFEM *) m_model.GetVerticalDynamics();

	// Check that vertical dynamics is set to full explicit evaluation
	if (pVerticalDynamics->IsFullyExplicit() == false) {
		_EXCEPTIONT("Set --explicitvertical in the command line for"  
				" Split Explicit --spex time integration!");
	}

	// Set up the small step partition based on vertical sound speed
	// EVERYTHING IS HARD CODED FOR NOW -- CLEAN UP TO USE PHYSICAL CONSTANTS
	int nRElements = pGrid->GetRElements();
	double dZtop = pGrid->GetZtop();
	double dZGridSpace = dZtop / nRElements;
	// Compute the stiff time step based on CFL = 1 and c = 350.0 m/s
	double dStiffDT = dZGridSpace / 350.0;
	int ns = int (2.0 * dDeltaT / dStiffDT);
	//ns = 20;
	//std::cout << "Number of small steps: " << ns << std::endl;

	// Take the full horizontal step with KGU53
	pGrid->CopyData(0, 1, DataType_State);
	pGrid->CopyData(0, 1, DataType_Tracers);
	pHorizontalDynamics->StepExplicit(0, 1, time, dDeltaT / 5.0);
	pGrid->PostProcessSubstage(1, DataType_State);
	pGrid->PostProcessSubstage(1, DataType_Tracers);

	pGrid->CopyData(0, 2, DataType_State);
	pGrid->CopyData(0, 2, DataType_Tracers);
	pHorizontalDynamics->StepExplicit(1, 2, time, dDeltaT / 5.0);
	pGrid->PostProcessSubstage(2, DataType_State);
	pGrid->PostProcessSubstage(2, DataType_Tracers);

	pGrid->CopyData(0, 3, DataType_State);
	pGrid->CopyData(0, 3, DataType_Tracers);
	pHorizontalDynamics->StepExplicit(2, 3, time, dDeltaT / 3.0);
	pGrid->PostProcessSubstage(3, DataType_State);
	pGrid->PostProcessSubstage(3, DataType_Tracers);

	pGrid->CopyData(0, 2, DataType_State);
	pGrid->CopyData(0, 2, DataType_Tracers);
	pHorizontalDynamics->StepExplicit(3, 2, time, 2.0 * dDeltaT / 3.0);
	pGrid->PostProcessSubstage(2, DataType_State);
	pGrid->PostProcessSubstage(2, DataType_Tracers);

	pGrid->LinearCombineData(
		m_dKinnmarkGrayUllrichCombination, 4, DataType_State);
	pGrid->LinearCombineData(
		m_dKinnmarkGrayUllrichCombination, 4, DataType_Tracers);
	pHorizontalDynamics->StepExplicit(2, 4, time, 3.0 * dDeltaT / 4.0);
	pGrid->PostProcessSubstage(4, DataType_State);
	pGrid->PostProcessSubstage(4, DataType_Tracers);

/*	// Take the full horizontal step with SSPRK3
	pHorizontalDynamics->StepExplicit(0, 1, time, dDeltaT);
	pGrid->PostProcessSubstage(1, DataType_State);
	pGrid->PostProcessSubstage(1, DataType_Tracers);

	pGrid->LinearCombineData(m_dSSPRK3CombinationA, 2, DataType_State);
	pGrid->LinearCombineData(m_dSSPRK3CombinationA, 2, DataType_Tracers);
	pHorizontalDynamics->StepExplicit(1, 2, time, 0.25 * dDeltaT);
	pGrid->PostProcessSubstage(2, DataType_State);
	pGrid->PostProcessSubstage(2, DataType_Tracers);

	pGrid->LinearCombineData(m_dSSPRK3CombinationB, 4, DataType_State);
	pGrid->LinearCombineData(m_dSSPRK3CombinationB, 4, DataType_Tracers);
	pHorizontalDynamics->StepExplicit(2, 4, time, (2.0/3.0) * dDeltaT);
	pGrid->PostProcessSubstage(4, DataType_State);
	pGrid->PostProcessSubstage(4, DataType_Tracers);

	// Apply hyperdiffusion
	pGrid->CopyData(4, 1, DataType_State);
	pGrid->CopyData(4, 1, DataType_Tracers);
	pHorizontalDynamics->StepAfterSubCycle(4, 1, 2, time, dDeltaT);
	pGrid->CopyData(1, 0, DataType_State);
	pGrid->CopyData(1, 0, DataType_Tracers);
*/
	pGrid->CopyData(4, 0, DataType_State);
	pGrid->CopyData(4, 0, DataType_Tracers);

	//std::cout << "Entering substages at the timestep... \n";
	// Compute the small step loop for stiff vertical terms

	for (int n = 0; n < ns; n++) {
		pGrid->CopyData(0, 1, DataType_State);
		pGrid->CopyData(0, 1, DataType_Tracers);
/*
		pVerticalDynamics->StepExplicit(0, 1, time, dDeltaT / 5.0 / ns);
		pGrid->PostProcessSubstage(1, DataType_State);
		pGrid->PostProcessSubstage(1, DataType_Tracers);

		pGrid->CopyData(0, 2, DataType_State);
		pGrid->CopyData(0, 2, DataType_Tracers);
		pVerticalDynamics->StepExplicit(1, 2, time, dDeltaT / 5.0 / ns);
		pGrid->PostProcessSubstage(2, DataType_State);
		pGrid->PostProcessSubstage(2, DataType_Tracers);

		pGrid->CopyData(0, 3, DataType_State);
		pGrid->CopyData(0, 3, DataType_Tracers);
		pVerticalDynamics->StepExplicit(2, 3, time, dDeltaT / 3.0 / ns);
		pGrid->PostProcessSubstage(3, DataType_State);
		pGrid->PostProcessSubstage(3, DataType_Tracers);

		pGrid->CopyData(0, 2, DataType_State);
		pGrid->CopyData(0, 2, DataType_Tracers);
		pVerticalDynamics->StepExplicit(3, 2, time, 2.0 * dDeltaT / 3.0 / ns);
		pGrid->PostProcessSubstage(2, DataType_State);
		pGrid->PostProcessSubstage(2, DataType_Tracers);

		pGrid->LinearCombineData(
		m_dKinnmarkGrayUllrichCombination, 4, DataType_State);
		pGrid->LinearCombineData(
		m_dKinnmarkGrayUllrichCombination, 4, DataType_Tracers);
		pVerticalDynamics->StepExplicit(2, 4, time, 3.0 * dDeltaT / 4.0 / ns);
		pGrid->PostProcessSubstage(4, DataType_State);
		pGrid->PostProcessSubstage(4, DataType_Tracers);
*/
		pVerticalDynamics->StepExplicit(0, 1, time, dDeltaT / ns);
		pGrid->PostProcessSubstage(1, DataType_State);
		pGrid->PostProcessSubstage(1, DataType_Tracers);

		pGrid->LinearCombineData(m_dSSPRK3CombinationA, 2, DataType_State);
		pGrid->LinearCombineData(m_dSSPRK3CombinationA, 2, DataType_Tracers);
		pVerticalDynamics->StepExplicit(1, 2, time, 0.25 * dDeltaT / ns);
		pGrid->PostProcessSubstage(2, DataType_State);
		pGrid->PostProcessSubstage(2, DataType_Tracers);

		pGrid->LinearCombineData(m_dSSPRK3CombinationB, 4, DataType_State);
		pGrid->LinearCombineData(m_dSSPRK3CombinationB, 4, DataType_Tracers);
		pVerticalDynamics->StepExplicit(2, 4, time, (2.0/3.0) * dDeltaT / ns);
		pGrid->PostProcessSubstage(4, DataType_State);
		pGrid->PostProcessSubstage(4, DataType_Tracers);

		if (n == ns - 1) {
		// Apply hyperdiffusion
		pGrid->CopyData(4, 1, DataType_State);
		pGrid->CopyData(4, 1, DataType_Tracers);
		pHorizontalDynamics->StepAfterSubCycle(4, 1, 2, time, dDeltaT);
		pGrid->CopyData(1, 0, DataType_State);
		pGrid->CopyData(1, 0, DataType_Tracers);
		}
		else {
		pGrid->CopyData(4, 0, DataType_State);
		pGrid->CopyData(4, 0, DataType_Tracers);
		}
	}
}

///////////////////////////////////////////////////////////////////////////////


