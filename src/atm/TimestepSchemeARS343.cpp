///////////////////////////////////////////////////////////////////////////////
///
///	\file    TimestepSchemeARS343.h
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

#include "TimestepSchemeARS343.h"
#include "Model.h"
#include "Grid.h"
#include "HorizontalDynamics.h"
#include "VerticalDynamics.h"

///////////////////////////////////////////////////////////////////////////////
// COEFFICIENTS COMPUTED FROM THE ORIGINAL TABLEAUX
// IMPLEMENTS ARS(3,4,3) FROM ASCHER ET AL. 1997 PG. 9
//
const double TimestepSchemeARS343::m_dgamma = 0.4358665215084590;
const double TimestepSchemeARS343::m_db1 = -1.5 * m_dgamma * m_dgamma + 
											4.0 * m_dgamma - 0.25;
const double TimestepSchemeARS343::m_db2 = 1.5 * m_dgamma * m_dgamma -
											5.0 * m_dgamma + 1.2;
// Implicit stage coefficients - Converted to U from Ascher et al. 1997 Pg 9
const double TimestepSchemeARS343::m_dImpCf[4][4] = {
	{m_dgamma, 0., 0., 0.},
	{0.5 * (1.0 - m_dgamma), m_dgamma, 0., 0.},
	{m_db1, m_db2, m_dgamma, 0.},
	{1.208496649176010, -0.6443631706844690, m_dgamma, 0.}};

// Explicit stage coefficients - Converted to U from Ascher et al. 1997 Pg 9
const double TimestepSchemeARS343::m_dExpCf[4][4] = {
	{m_dgamma, 0., 0., 0.},
	{0.3212788860286278, 0.3966543747256017, 0., 0.},
	{-0.1058582960718797, 0.5529291480359398, 0.5529291480359398, 0.},
	{0., 1.208496649176010, -0.6443631706844690, m_dgamma}};
/*
// MODIFIED ARS343 COEFFICIENTS COMPUTED BY (Boscarino, 2009)
const double TimestepSchemeARS343::m_dgamma = 0.4358665215084590;
const double TimestepSchemeARS343::m_db1 = -1.5 * m_dgamma * m_dgamma + 
											4.0 * m_dgamma - 0.25;
const double TimestepSchemeARS343::m_db2 = 1.5 * m_dgamma * m_dgamma -
											5.0 * m_dgamma + 1.2;

// Implicit stage coefficients - Boscarino, 2009
const double TimestepSchemeARS343::m_dImpCf[4][4] = {
	{m_dgamma, 0., 0., 0.},
	{0.5 * (1.0 - m_dgamma), m_dgamma, 0., 0.},
	{m_db1, m_db2, m_dgamma, 0.},
	{1.208496649176010, -0.6443631706844690, m_dgamma, 0.}};

// Explicit stage coefficients - Boscarino, 2009
const double TimestepSchemeARS343::m_dExpCf[4][4] = {
	{m_dgamma, 0., 0., 0.},
	{0.5353965403073543, 0.1825367204468752, 0., 0.},
	{0.6304125581528727, -0.8319339010630845, 1.201521342910212, 0.},
	{0., 1.208496649176010, -0.6443631706844690, m_dgamma}};
*/
TimestepSchemeARS343::TimestepSchemeARS343(
	Model & model
) :
	TimestepScheme(model)
{
	m_du2fCombo.Allocate(8);
	m_du3fCombo.Allocate(9);
	m_du4fCombo.Allocate(10);
}

///////////////////////////////////////////////////////////////////////////////

void TimestepSchemeARS343::Step(
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

	// u2 explicit evaluation combination
	m_du2fCombo[0] = 1.0 - m_dExpCf[1][0] / m_dExpCf[0][0];
	m_du2fCombo[1] = m_dExpCf[1][0] / m_dExpCf[0][0] - 
					 m_dImpCf[1][0] / m_dImpCf[0][0];
	m_du2fCombo[2] = m_dImpCf[1][0] / m_dImpCf[0][0];
	m_du2fCombo[3] = 0.0;
	m_du2fCombo[4] = 0.0;
	m_du2fCombo[5] = 0.0;
	m_du2fCombo[6] = 0.0;
	m_du2fCombo[7] = 0.0;

	// u3 explicit evaluation combination
	m_du3fCombo[0] = 1.0 - m_dExpCf[2][0] / m_dExpCf[0][0];
	m_du3fCombo[1] = m_dExpCf[2][0] / m_dExpCf[0][0] - 
					 m_dImpCf[2][0] / m_dImpCf[0][0];
	m_du3fCombo[2] = m_dImpCf[2][0] / m_dImpCf[0][0];
	m_du3fCombo[3] = m_dExpCf[2][1] / m_dExpCf[1][1] - 
					 m_dImpCf[2][1] / m_dImpCf[1][1];
	m_du3fCombo[4] = m_dImpCf[2][1] / m_dImpCf[1][1];
	m_du3fCombo[5] = 0.0;
	m_du3fCombo[6] = 0.0;
	m_du3fCombo[7] = -m_dExpCf[2][1] / m_dExpCf[1][1];
	m_du3fCombo[8] = 0.0;
	
	// u4 explicit evaluation combination
	m_du4fCombo[0] = 1.0 - m_dExpCf[3][0] / m_dExpCf[0][0];
	m_du4fCombo[1] = m_dExpCf[3][0] / m_dExpCf[0][0] - 
					 m_dImpCf[3][0] / m_dImpCf[0][0];
	m_du4fCombo[2] = m_dImpCf[3][0] / m_dImpCf[0][0];
	m_du4fCombo[3] = m_dExpCf[3][1] / m_dExpCf[1][1] - 
					 m_dImpCf[3][1] / m_dImpCf[1][1];
	m_du4fCombo[4] = m_dImpCf[3][1] / m_dImpCf[1][1];
	m_du4fCombo[5] = m_dExpCf[3][2] / m_dExpCf[2][2] - 
					 m_dImpCf[3][2] / m_dImpCf[2][2];
	m_du4fCombo[6] = m_dImpCf[3][2] / m_dImpCf[2][2];
	m_du4fCombo[7] = -m_dExpCf[3][1] / m_dExpCf[1][1];
	m_du4fCombo[8] = -m_dExpCf[3][2] / m_dExpCf[2][2];
	m_du4fCombo[9] = 0.0;

	// STAGE 1
	// Compute uf1 into index 1
	pGrid->CopyData(0, 1, DataType_State);
	pGrid->CopyData(0, 1, DataType_Tracers);
	pHorizontalDynamics->StepExplicit(
		0, 1, time, m_dExpCf[0][0] * dDeltaT);
	pVerticalDynamics->StepExplicit(
		0, 1, time, m_dExpCf[0][0] * dDeltaT);
	pGrid->PostProcessSubstage(1, DataType_State);
	pGrid->PostProcessSubstage(1, DataType_Tracers);

	// Compute u1 into index 2
	pGrid->CopyData(1, 2, DataType_State);
	pGrid->CopyData(1, 2, DataType_Tracers);
	pVerticalDynamics->StepImplicit(
		2, 2, time, m_dImpCf[0][0] * dDeltaT);
	pGrid->PostProcessSubstage(2, DataType_State);
	pGrid->PostProcessSubstage(2, DataType_Tracers);

	// STAGE 2
	// Compute uf2 from u1 (index 2) into index 7
	pGrid->LinearCombineData(m_du2fCombo, 7, DataType_State);
	pGrid->LinearCombineData(m_du2fCombo, 7, DataType_Tracers);
	pGrid->CopyData(7, 3, DataType_State);
	pGrid->CopyData(7, 3, DataType_Tracers);
	pHorizontalDynamics->StepExplicit(
		2, 3, time, m_dExpCf[1][1] * dDeltaT);
	pVerticalDynamics->StepExplicit(
		2, 3, time, m_dExpCf[1][1] * dDeltaT);
	pGrid->PostProcessSubstage(3, DataType_State);
	pGrid->PostProcessSubstage(3, DataType_Tracers);

	// Compute u2 from uf2 (index 3) into index 4
	pGrid->CopyData(3, 4, DataType_State);
	pGrid->CopyData(3, 4, DataType_Tracers);
	pVerticalDynamics->StepImplicit(
		4, 4, time, m_dImpCf[1][1] * dDeltaT);
	pGrid->PostProcessSubstage(4, DataType_State);
	pGrid->PostProcessSubstage(4, DataType_Tracers);

	// STAGE 3
	// Compute uf3 from u2 (index 4) into index 8
	pGrid->LinearCombineData(m_du3fCombo, 8, DataType_State);
	pGrid->LinearCombineData(m_du3fCombo, 8, DataType_Tracers);
	pGrid->CopyData(8, 5, DataType_State);
	pGrid->CopyData(8, 5, DataType_Tracers);
	pHorizontalDynamics->StepExplicit(
		4, 5, time, m_dExpCf[2][2] * dDeltaT);
	pVerticalDynamics->StepExplicit(
		4, 5, time, m_dExpCf[2][2] * dDeltaT);
	pGrid->PostProcessSubstage(5, DataType_State);
	pGrid->PostProcessSubstage(5, DataType_Tracers);

	// Compute u3 from uf3 (index 3) into index 6
	pGrid->CopyData(5, 6, DataType_State);
	pGrid->CopyData(5, 6, DataType_Tracers);
	pVerticalDynamics->StepImplicit(
		6, 6, time, m_dImpCf[2][2] * dDeltaT);
	pGrid->PostProcessSubstage(6, DataType_State);
	pGrid->PostProcessSubstage(6, DataType_Tracers);

	// STAGE 4
	// Compute uf4 from u3 (index 6) into index 9
	pGrid->LinearCombineData(m_du4fCombo, 9, DataType_State);
	pGrid->LinearCombineData(m_du4fCombo, 9, DataType_Tracers);
	pHorizontalDynamics->StepExplicit(
		6, 9, time, m_dExpCf[3][3] * dDeltaT);
	pVerticalDynamics->StepExplicit(
		6, 9, time, m_dExpCf[3][3] * dDeltaT);
	pGrid->PostProcessSubstage(9, DataType_State);
	pGrid->PostProcessSubstage(9, DataType_Tracers);

	// NO IMPLICIT STEP ON THE LAST STAGE

	// Apply hyperdiffusion at the end of the explicit substep (ask Paul)
	pGrid->CopyData(9, 2, DataType_State);
	pGrid->CopyData(9, 2, DataType_Tracers);
	pHorizontalDynamics->StepAfterSubCycle(2, 1, 9, time, dDeltaT);
	pGrid->CopyData(1, 0, DataType_State);
	pGrid->CopyData(1, 0, DataType_Tracers);
}

///////////////////////////////////////////////////////////////////////////////

