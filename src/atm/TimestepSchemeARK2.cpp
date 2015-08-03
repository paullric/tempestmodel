///////////////////////////////////////////////////////////////////////////////
///
///	\file	TimestepSchemeARK2.cpp
///	\author  Paul Ullrich
///	\version April 22, 2014
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

#include "TimestepSchemeARK2.h"
#include "Model.h"
#include "Grid.h"
#include "HorizontalDynamics.h"
#include "VerticalDynamics.h"

#pragma message "File tagged for cleanup"

///////////////////////////////////////////////////////////////////////////////
// THESE COEFFICIENTS ARE COMPUTED FROM THE ORIGINAL TABLEAUX
// IMPLEMENTS ARK(2,2,2) FROM ASCHER ET AL. 1997 PG. 9
// Tableaux parameters
const double TimestepSchemeARK2::m_dgamma = 1.0 - 0.5 * std::sqrt(2.0);
const double TimestepSchemeARK2::m_ddelta = 1.0 - 1.0 / (2.0 * m_dgamma);
// Time step coefficients
const double TimestepSchemeARK2::m_dTimeCf[2] = {m_dgamma, 1.0};
// Implicit stage coefficients - Converted to U from Ascher et al. 1997 Pg 9
const double TimestepSchemeARK2::m_dImpCf[2][2] = {
  {m_dgamma, 0.0},
  {1.0 - m_dgamma, m_dgamma}};

// Explicit stage coefficients - Converted to U from Ascher et al. 1997 Pg 9
const double TimestepSchemeARK2::m_dExpCf[3][3] = {
  {0.0, 0.0, 0.0},
  {m_dgamma, 0.0, 0.0},
  {m_ddelta, 1.0 - m_ddelta, 0.0}};

// Final stage coefficients
const double TimestepSchemeARK2::m_dBCoeff1[2] = {1.0 - m_dgamma, m_dgamma};
const double TimestepSchemeARK2::m_dBCoeff2[3] =
			 {m_ddelta, 1.0 - m_ddelta, 0.0};

TimestepSchemeARK2::TimestepSchemeARK2(
	Model & model
) :
	TimestepScheme(model)
{
	m_dKh1Combo.Allocate(5);
	m_dK1Combo.Allocate(5);
	m_du2fCombo.Allocate(5);
}

///////////////////////////////////////////////////////////////////////////////

void TimestepSchemeARK2::Step(
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

	// Kh1 combination
	m_dKh1Combo[0] = -1.0 / (dDeltaT * m_dExpCf[1][0]);
	m_dKh1Combo[1] = 1.0 / (dDeltaT * m_dExpCf[1][0]);
	m_dKh1Combo[2] = 0.0;
	m_dKh1Combo[3] = 0.0;
	m_dKh1Combo[4] = 0.0;

	// K1 combination
	m_dK1Combo[0] = 0.0;
	m_dK1Combo[1] = -1.0 / (dDeltaT * m_dImpCf[0][0]);
	m_dK1Combo[2] = 1.0 / (dDeltaT * m_dImpCf[0][0]);
	m_dK1Combo[3] = 0.0;
	m_dK1Combo[4] = 0.0;

	// u2 explicit evaluation combination
	m_du2fCombo[0] = 1.0;
	m_du2fCombo[1] = 0.0;
	m_du2fCombo[2] = 0.0;
	m_du2fCombo[3] = dDeltaT * m_dImpCf[1][0];
	m_du2fCombo[4] = dDeltaT * m_dExpCf[2][0];

	//std::cout << "Entering substages at the timestep... \n";

	// SUBSTAGE 1
	// Compute uf1
	pGrid->CopyData(0, 1, DataType_State);
	pHorizontalDynamics->StepExplicit(
		0, 1, time /*+ m_dTimeCf[0] * dDeltaT*/, m_dExpCf[1][0] * dDeltaT);
	pVerticalDynamics->StepExplicit(
		0, 1, time /*+ m_dTimeCf[0] * dDeltaT*/, m_dExpCf[1][0] * dDeltaT);
	pGrid->PostProcessSubstage(1, DataType_State);
	pGrid->PostProcessSubstage(1, DataType_Tracers);

	// Store the evaluation Kh1 to index 4
	pGrid->LinearCombineData(m_dKh1Combo, 4, DataType_State);

	// Compute u1
	pGrid->CopyData(1, 2, DataType_State);
	pVerticalDynamics->StepImplicit(
		1, 2, time /*+ m_dTimeCf[0] * dDeltaT*/, m_dImpCf[0][0] * dDeltaT);
	pGrid->PostProcessSubstage(2, DataType_State);
	pGrid->PostProcessSubstage(2, DataType_Tracers);

	// Store the evaluation K1 to index 3
	pGrid->LinearCombineData(m_dK1Combo, 3, DataType_State);

	// SUBSTAGE 2
	// Compute uf2 from u1 and store it to index 1 (over uf1)
	pGrid->LinearCombineData(m_du2fCombo, 1, DataType_State);
	pHorizontalDynamics->StepExplicit(
		2, 1, time /*+ m_dTimeCf[1] * dDeltaT*/, m_dExpCf[2][1] * dDeltaT);
	pVerticalDynamics->StepExplicit(
		2, 1, time /*+ m_dTimeCf[1] * dDeltaT*/, m_dExpCf[2][1] * dDeltaT);
	pGrid->PostProcessSubstage(1, DataType_State);
	pGrid->PostProcessSubstage(1, DataType_Tracers);

	// Compute u2 from uf2 and u2 and store it to index 0 (over u0)
	pGrid->CopyData(1, 2, DataType_State);
	pVerticalDynamics->StepImplicit(
		1, 2, time /*+ m_dTimeCf[1] * dDeltaT*/, m_dImpCf[1][1] * dDeltaT);
	pGrid->PostProcessSubstage(2, DataType_State);
	pGrid->PostProcessSubstage(2, DataType_Tracers);

	// Apply hyperdiffusion at the end of the explicit substep (ask Paul)
	pGrid->CopyData(2, 1, DataType_State);
	pHorizontalDynamics->StepAfterSubCycle(2, 1, 3, time, dDeltaT);
	pGrid->CopyData(1, 0, DataType_State);
}

///////////////////////////////////////////////////////////////////////////////


