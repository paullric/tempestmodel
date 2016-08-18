///////////////////////////////////////////////////////////////////////////////
///
///	\file	TimestepSchemeGARK2.cpp
///	\author  Jorge Guerra
///	\version January 29, 2016
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

#include "TimestepSchemeGARK2.h"
#include "Model.h"
#include "Grid.h"
#include "HorizontalDynamics.h"
#include "VerticalDynamics.h"

///////////////////////////////////////////////////////////////////////////////
// THESE COEFFICIENTS ARE COMPUTED FROM THE ORIGINAL TABLEAUX
// IMPLEMENTS 2nd ORDER IMEX GARK (Example 7, Sandu & Gunther, 2013)
// Tableaux parameters
const double TimestepSchemeGARK2::m_dgamma = 1.0 - 0.5 * std::sqrt(2.0);
const double TimestepSchemeGARK2::m_dalpha = 0.5;

// Implicit stage coefficients
const double TimestepSchemeGARK2::m_dImpCf[2][2] = {
	{m_dgamma, 0.0},
	{1.0 - m_dgamma, m_dgamma}};

// Explicit stage coefficients
const double TimestepSchemeGARK2::m_dExpCf[2][2] = {
	{0.0, 0.0},
	{1.0, 0.0}};

// Explicit-Implicit stage coupling coefficients
const double TimestepSchemeGARK2::m_dEICf[2][2] = {
	{0.0, 0.0},
	{1.0, 0.0}};

// Implicit-Explicit stage coupling coefficients
const double TimestepSchemeGARK2::m_dIECf[2][2] = {
	{m_dgamma, 0.0},
	{m_dalpha, 1.0 - m_dalpha}};

TimestepSchemeGARK2::TimestepSchemeGARK2(
	Model & model
) :
	TimestepScheme(model)
{
	m_du2fCombo.Allocate(4);
	m_du3fCombo.Allocate(5);
}

///////////////////////////////////////////////////////////////////////////////

void TimestepSchemeGARK2::Step(
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

	// Y2 evaluation combination
	m_du2fCombo[0] = 1.0 - m_dExpCf[1][0] / m_dIECf[0][0];
	m_du2fCombo[1] = m_dExpCf[1][0] / m_dIECf[0][0] - 
					 m_dEICf[1][0] / m_dImpCf[0][0];
	m_du2fCombo[2] = m_dEICf[1][0] / m_dImpCf[0][0];
	m_du2fCombo[3] = 0.0;

	// Z2 evaluation combination
	m_du3fCombo[0] = 1.0 - m_dIECf[1][0] / m_dIECf[0][0];
	m_du3fCombo[1] = m_dIECf[1][0] / m_dIECf[0][0] - 
					 m_dImpCf[1][0] / m_dImpCf[0][0];
	m_du3fCombo[2] = m_dImpCf[1][0] / m_dImpCf[0][0];
	m_du3fCombo[3] = 0.0;
	m_du3fCombo[4] = 0.0;

	//std::cout << "Entering substages at the timestep... \n";

	// STAGE 1
	// Compute uf1 into index 1
	pGrid->CopyData(0, 1, DataType_State);
	pGrid->CopyData(0, 1, DataType_Tracers);
	pHorizontalDynamics->StepExplicit(
		0, 1, time, m_dIECf[0][0] * dDeltaT);
	pVerticalDynamics->StepExplicit(
		0, 1, time, m_dIECf[0][0] * dDeltaT);
	pGrid->PostProcessSubstage(1, DataType_State);
	pGrid->PostProcessSubstage(1, DataType_Tracers);

	// Compute u1 into index 2
	pGrid->CopyData(1, 2, DataType_State);
	pGrid->CopyData(1, 2, DataType_State);
	pVerticalDynamics->StepImplicit(
		2, 2, time, m_dImpCf[0][0] * dDeltaT);
	pGrid->PostProcessSubstage(2, DataType_State);
	pGrid->PostProcessSubstage(2, DataType_Tracers);

	// Compute uf2 from u1 (index 2) into index 3
	pGrid->LinearCombineData(m_du2fCombo, 3, DataType_State);
	pGrid->LinearCombineData(m_du2fCombo, 3, DataType_Tracers);

	// STAGE 2
	// Compute uf2 from u1 (index 2) into index 3
	pGrid->LinearCombineData(m_du3fCombo, 4, DataType_State);
	pGrid->LinearCombineData(m_du3fCombo, 4, DataType_Tracers);
	pHorizontalDynamics->StepExplicit(
		3, 4, time, m_dIECf[1][1] * dDeltaT);
	pVerticalDynamics->StepExplicit(
		3, 4, time, m_dIECf[1][1] * dDeltaT);
	pGrid->PostProcessSubstage(4, DataType_State);
	pGrid->PostProcessSubstage(4, DataType_Tracers);

	// Compute u2 from uf2 (index 3) into index 2
	pVerticalDynamics->StepImplicit(
		4, 4, time, m_dImpCf[1][1] * dDeltaT);
	pGrid->PostProcessSubstage(4, DataType_State);
	pGrid->PostProcessSubstage(4, DataType_Tracers);

	// Apply hyperdiffusion at the end of the explicit substep (ask Paul)
	pGrid->CopyData(4, 2, DataType_State);
	pGrid->CopyData(4, 2, DataType_Tracers);
	pHorizontalDynamics->StepAfterSubCycle(2, 1, 4, time, dDeltaT);
	pGrid->CopyData(1, 0, DataType_State);
	pGrid->CopyData(1, 0, DataType_Tracers);
}

///////////////////////////////////////////////////////////////////////////////


