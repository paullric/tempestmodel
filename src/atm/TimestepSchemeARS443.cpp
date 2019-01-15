///
///	\file    TimestepSchemeARS443.cpp
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

#include "TimestepSchemeARS443.h"
#include "Model.h"
#include "Grid.h"
#include "HorizontalDynamics.h"
#include "VerticalDynamics.h"

TimestepSchemeARS443::TimestepSchemeARS443(
	Model & model
) :
	TimestepScheme(model)
{
	m_dU2fCombo.Allocate(7);
	m_dU3fCombo.Allocate(7);
	m_dU4fCombo.Allocate(7);

	m_dDiagExpCf.Allocate(4);
	m_dDiagImpCf.Allocate(4);

	///////////////////////////////////////////////////////////////////////////////
	// THESE COEFFICIENTS ARE COMPUTED FROM THE ORIGINAL TABLEAUX
	// IMPLEMENTS ARK(4,4,3) FROM ASCHER ET AL. 1997 PG. 9
	// Time step coefficients
	const double dTimeCf[4] = {1./2., 2./3., 1./2., 1.};
	// Implicit stage coefficients - Converted to U from Ascher et al. 1997 Pg 9
	const double dImpCf[4][4] = {
		{1./2., 0., 0., 0.},
		{1./6., 1./2., 0., 0.},
		{-1./2., 1./2., 1./2., 0.},
		{3./2., -3./2., 1./2., 1./2.}};

	// Explicit stage coefficients - Converted to U from Ascher et al. 1997 Pg 9
	const double dExpCf[4][4] = {
		{1./2., 0., 0., 0.},
		{11./18., 1./18., 0., 0.},
		{5./6., -5./6., 1./2., 0.},
		{1./4., 7./4., 3./4., -7./4.}};

	// Diagnoal explicit coefficients
	m_dDiagExpCf[0] = dExpCf[0][0];
	m_dDiagExpCf[1] = dExpCf[1][1];
	m_dDiagExpCf[2] = dExpCf[2][2];
	m_dDiagExpCf[3] = dExpCf[3][3];

	// Diagnoal implicit coefficients
	m_dDiagImpCf[0] = dImpCf[0][0];
	m_dDiagImpCf[1] = dImpCf[1][1];
	m_dDiagImpCf[2] = dImpCf[2][2];
	m_dDiagImpCf[3] = dImpCf[3][3];

	// u2 explicit evaluation combination
	m_dU2fCombo[0] = 1.0 - dExpCf[1][0] / dExpCf[0][0];
	m_dU2fCombo[1] = dExpCf[1][0] / dExpCf[0][0] -
					 dImpCf[1][0] / dImpCf[0][0];
	m_dU2fCombo[2] = dImpCf[1][0] / dImpCf[0][0];
	m_dU2fCombo[3] = 0.0;
	m_dU2fCombo[4] = 0.0;
	m_dU2fCombo[5] = 0.0;
	m_dU2fCombo[6] = 0.0;

	// u3 explicit evaluation combination
	m_dU3fCombo[0] = 1.0 - dExpCf[2][0] / dExpCf[0][0];
	m_dU3fCombo[1] = dExpCf[2][0] / dExpCf[0][0] -
					 dImpCf[2][0] / dImpCf[0][0];
	m_dU3fCombo[2] = dImpCf[2][0] / dImpCf[0][0];
	m_dU3fCombo[3] = dExpCf[2][1] / dExpCf[1][1] -
					 dImpCf[2][1] / dImpCf[1][1];
	m_dU3fCombo[4] = dImpCf[2][1] / dImpCf[1][1];
	m_dU3fCombo[5] = 0.0;
	m_dU3fCombo[6] = 0.0;

	const double dU3fCombo7 = -dExpCf[2][1] / dExpCf[1][1];

	// u4 explicit evaluation combination
	m_dU4fCombo[0] = 1.0 - dExpCf[3][0] / dExpCf[0][0];
	m_dU4fCombo[1] = dExpCf[3][0] / dExpCf[0][0] -
					 dImpCf[3][0] / dImpCf[0][0];
	m_dU4fCombo[2] = dImpCf[3][0] / dImpCf[0][0];
	m_dU4fCombo[3] = dExpCf[3][1] / dExpCf[1][1] -
					 dImpCf[3][1] / dImpCf[1][1];
	m_dU4fCombo[4] = dImpCf[3][1] / dImpCf[1][1];
	m_dU4fCombo[5] = dExpCf[3][2] / dExpCf[2][2] -
					 dImpCf[3][2] / dImpCf[2][2];
	m_dU4fCombo[6] = dImpCf[3][2] / dImpCf[2][2];

	const double dU4fCombo7 = -dExpCf[3][1] / dExpCf[1][1];
	const double dU4fCombo8 = -dExpCf[3][2] / dExpCf[2][2];

	// Recombination terms
	m_dU3fCombo[0] += dU3fCombo7 * m_dU2fCombo[0];
	m_dU3fCombo[1] += dU3fCombo7 * m_dU2fCombo[1];
	m_dU3fCombo[2] += dU3fCombo7 * m_dU2fCombo[2];

	m_dU4fCombo[0] +=
		  dU4fCombo7 * m_dU2fCombo[0]
		+ dU4fCombo8 * m_dU3fCombo[0];
	m_dU4fCombo[1] +=
		  dU4fCombo7 * m_dU2fCombo[1]
		+ dU4fCombo8 * m_dU3fCombo[1];
	m_dU4fCombo[2] +=
		  dU4fCombo7 * m_dU2fCombo[2]
		+ dU4fCombo8 * m_dU3fCombo[2];

	m_dU4fCombo[3] += dU4fCombo8 * m_dU3fCombo[3];
	m_dU4fCombo[4] += dU4fCombo8 * m_dU3fCombo[4];
}

///////////////////////////////////////////////////////////////////////////////

void TimestepSchemeARS443::Step(
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
	pGrid->CopyData(1, 2, DataType_State);
	pVerticalDynamics->StepImplicit(
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
	// Compute uf3 from u2 (index 4) into index 8
	pGrid->LinearCombineData(m_dU3fCombo, 5, DataType_State);
	pGrid->LinearCombineData(m_dU3fCombo, 5, DataType_Tracers);
	pHorizontalDynamics->StepExplicit(
		4, 5, time, m_dDiagExpCf[2] * dDeltaT);
	pVerticalDynamics->StepExplicit(
		4, 5, time, m_dDiagExpCf[2] * dDeltaT);
	pGrid->PostProcessSubstage(5, DataType_State);
	pGrid->PostProcessSubstage(5, DataType_Tracers);

	// Compute u3 from uf3 (index 5) into index 6
	pGrid->CopyData(5, 6, DataType_State);
	pGrid->CopyData(5, 6, DataType_Tracers);
	pVerticalDynamics->StepImplicit(
		6, 6, time, m_dDiagImpCf[2] * dDeltaT);
	pHorizontalDynamics->StepImplicit(
		6, 6, time, m_dDiagImpCf[2] * dDeltaT);

	// STAGE 4
	// Compute uf4 from u3 (index 6) into index 1
	pGrid->LinearCombineData(m_dU4fCombo, 1, DataType_State);
	pGrid->LinearCombineData(m_dU4fCombo, 1, DataType_Tracers);
	pHorizontalDynamics->StepExplicit(
		6, 1, time, m_dDiagExpCf[3] * dDeltaT);
	pVerticalDynamics->StepExplicit(
		6, 1, time, m_dDiagExpCf[3] * dDeltaT);
	pGrid->PostProcessSubstage(1, DataType_State);
	pGrid->PostProcessSubstage(1, DataType_Tracers);

	// Compute u4 from uf4 (index 6) into index 1
	pVerticalDynamics->StepImplicit(
		1, 1, time, m_dDiagImpCf[3] * dDeltaT);
	pHorizontalDynamics->StepImplicit(
		1, 1, time, m_dDiagImpCf[3] * dDeltaT);

	// Apply hyperdiffusion at the end of the explicit substep
	pGrid->CopyData(1, 0, DataType_State);
	pGrid->CopyData(1, 0, DataType_Tracers);
	pHorizontalDynamics->StepAfterSubCycle(1, 0, 2, time, dDeltaT);
}

///////////////////////////////////////////////////////////////////////////////
