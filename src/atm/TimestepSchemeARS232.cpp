///
///	\file    TimestepSchemeARS232.cpp
///	\author  Paul Ullrich
///	\version May 29, 2014
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

#include "TimestepSchemeARS232.h"
#include "Model.h"
#include "Grid.h"
#include "HorizontalDynamics.h"
#include "VerticalDynamics.h"

///////////////////////////////////////////////////////////////////////////////

TimestepSchemeARS232::TimestepSchemeARS232(
	Model & model
) :
	TimestepScheme(model)
{
	m_dU2fCombo.Allocate(5);
	m_dU3fCombo.Allocate(5);

	m_dDiagExpCf.Allocate(3);
	m_dDiagImpCf.Allocate(3);

	///////////////////////////////////////////////////////////////////////////////
	// THESE COEFFICIENTS ARE COMPUTED FROM THE ORIGINAL TABLEAUX
	// IMPLEMENTS ARS(2,3,2) FROM ASCHER ET AL. 1997 Section 2.5
	// Time step coefficients
	const double dGamma = 1.0 - 1.0 / std::sqrt(2.0);
	const double dDelta = -(2.0 * std::sqrt(2.0)) / 3.0;

	// Implicit stage coefficients
	const double dImpCf[3][3] = {
	  {dGamma, 0., 0.},
	  {1.0 - dGamma, dGamma, 0.},
	  {1.0 - dGamma, dGamma, 0.}};

	// Explicit stage coefficients
	const double dExpCf[3][3] = {
	  {dGamma, 0., 0.},
	  {dDelta, 1.0 - dDelta, 0.},
	  {0., 1.0 - dGamma, dGamma}};

	// Diagnoal explicit coefficients
  	m_dDiagExpCf[0] = dExpCf[0][0];
  	m_dDiagExpCf[1] = dExpCf[1][1];
  	m_dDiagExpCf[2] = dExpCf[2][2];

  	// Diagnoal implicit coefficients
  	m_dDiagImpCf[0] = dImpCf[0][0];
  	m_dDiagImpCf[1] = dImpCf[1][1];
  	m_dDiagImpCf[2] = dImpCf[2][2];

	// u2 explicit evaluation combination
	m_dU2fCombo[0] = 1.0 - dExpCf[1][0] / dExpCf[0][0];
	m_dU2fCombo[1] = dExpCf[1][0] / dExpCf[0][0] -
					 dImpCf[1][0] / dImpCf[0][0];
	m_dU2fCombo[2] = dImpCf[1][0] / dImpCf[0][0];
	m_dU2fCombo[3] = 0.0;
	m_dU2fCombo[4] = 0.0;

	// u3 explicit evaluation combination
	m_dU3fCombo[0] = 1.0 - dExpCf[2][0] / dExpCf[0][0];
	m_dU3fCombo[1] = dExpCf[2][0] / dExpCf[0][0] -
					 dImpCf[2][0] / dImpCf[0][0];
	m_dU3fCombo[2] = dImpCf[2][0] / dImpCf[0][0];
	m_dU3fCombo[3] = dExpCf[2][1] / dExpCf[1][1] -
					 dImpCf[2][1] / dImpCf[1][1];
	m_dU3fCombo[4] = dImpCf[2][1] / dImpCf[1][1];

	const double dU3fCombo5 = -dExpCf[2][1] / dExpCf[1][1];

	// Recombination terms
	m_dU3fCombo[0] += dU3fCombo5 * m_dU2fCombo[0];
	m_dU3fCombo[1] += dU3fCombo5 * m_dU2fCombo[1];
	m_dU3fCombo[2] += dU3fCombo5 * m_dU2fCombo[2];
}

///////////////////////////////////////////////////////////////////////////////

void TimestepSchemeARS232::Step(
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

	// STAGE 1
	// Compute uf1 into index 1
	pGrid->CopyData(0, 1, DataType_State);
	pGrid->CopyData(0, 1, DataType_Tracers);
	pHorizontalDynamics->StepExplicit(
		0, 1, time, m_dDiagExpCf[0] * dDeltaT);
	pVerticalDynamics->StepExplicit(
		0, 1, time, m_dDiagExpCf[0] * dDeltaT);
	pGrid->PostProcessSubstage(1, DataType_State);
	pGrid->PostProcessSubstage(1, DataType_Tracers);

	// Compute u1 into index 2
	pGrid->CopyData(1, 2, DataType_State);
	pGrid->CopyData(1, 2, DataType_Tracers);
	pVerticalDynamics->StepImplicit(
		2, 2, time, m_dDiagImpCf[0] * dDeltaT);
	pHorizontalDynamics->StepImplicit(
		2, 2, time, m_dDiagImpCf[0] * dDeltaT);

	// STAGE 2
	// Compute uf2 from u1 (index 2) into index 3
	pGrid->LinearCombineData(m_dU2fCombo, 3, DataType_State);
	pGrid->LinearCombineData(m_dU2fCombo, 3, DataType_Tracers);
	pHorizontalDynamics->StepExplicit(
		2, 3, time, m_dDiagExpCf[1] * dDeltaT);
	pVerticalDynamics->StepExplicit(
		2, 3, time, m_dDiagExpCf[1] * dDeltaT);
	pGrid->PostProcessSubstage(3, DataType_State);
	pGrid->PostProcessSubstage(3, DataType_Tracers);

	// Compute u2 from uf2 (index 3) into index 4
	pGrid->CopyData(3, 4, DataType_State);
	pGrid->CopyData(3, 4, DataType_Tracers);
	pVerticalDynamics->StepImplicit(
		4, 4, time, m_dDiagImpCf[1] * dDeltaT);
	pHorizontalDynamics->StepImplicit(
		4, 4, time, m_dDiagImpCf[1] * dDeltaT);

	// STAGE 3
	// Compute uf3 from u2 (index 4) into index 1
	pGrid->LinearCombineData(m_dU3fCombo, 1, DataType_State);
	pGrid->LinearCombineData(m_dU3fCombo, 1, DataType_Tracers);
	pHorizontalDynamics->StepExplicit(
		4, 1, time, m_dDiagExpCf[2] * dDeltaT);
	pVerticalDynamics->StepExplicit(
		4, 1, time, m_dDiagExpCf[2] * dDeltaT);
	pGrid->PostProcessSubstage(1, DataType_State);
	pGrid->PostProcessSubstage(1, DataType_Tracers);

	// NO IMPLICIT STEP ON THE LAST STAGE

	// Apply hyperdiffusion at the end of the explicit substep
	pGrid->CopyData(1, 0, DataType_State);
	pGrid->CopyData(1, 0, DataType_Tracers);
	pHorizontalDynamics->StepAfterSubCycle(1, 0, 2, time, dDeltaT);
}

///////////////////////////////////////////////////////////////////////////////
