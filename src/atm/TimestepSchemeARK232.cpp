///
///	\file    TimestepSchemeARK232.cpp
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

#include "TimestepSchemeARK232.h"
#include "Model.h"
#include "Grid.h"
#include "HorizontalDynamics.h"
#include "VerticalDynamics.h"

///////////////////////////////////////////////////////////////////////////////
// THESE COEFFICIENTS ARE COMPUTED FROM THE ORIGINAL TABLEAUX
// IMPLEMENTS ARK(2,3,2) FROM WELLER (Giraldo) ET AL. 2013 Fig. 2
// Time step coefficients
const double TimestepSchemeARK232::m_dgamma = 1.0 - 1.0 / std::sqrt(2.0);
const double TimestepSchemeARK232::m_ddelta = 1.0 / (2.0 * std::sqrt(2.0));
const double TimestepSchemeARK232::m_dalpha = 1.0 / 6.0 * 
                                              (3.0 + 2.0 * std::sqrt(2.0));
// Implicit stage coefficients
const double TimestepSchemeARK232::m_dImpCf[3][3] = {
  {m_dgamma, m_dgamma, 0.},
  {m_ddelta, m_ddelta, m_dgamma},
  {m_ddelta, m_ddelta, m_dgamma}};

// Explicit stage coefficients
const double TimestepSchemeARK232::m_dExpCf[3][3] = {
  {2.0 * m_dgamma, 0., 0.},
  {1.0 - m_dalpha, m_dalpha, 0.},
  {m_ddelta, m_ddelta, m_dgamma}};

TimestepSchemeARK232::TimestepSchemeARK232(
	Model & model
) :
	TimestepScheme(model)
{
    m_du2fCombo.Allocate(7);
    m_du3fCombo.Allocate(8);
}

///////////////////////////////////////////////////////////////////////////////

void TimestepSchemeARK232::Step(
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
	m_du2fCombo[2] = m_dImpCf[1][0] / m_dImpCf[0][0] - 
                         m_dImpCf[1][1] / m_dImpCf[0][1];
        m_du2fCombo[3] = m_dImpCf[1][1] / m_dImpCf[0][1];
        m_du2fCombo[4] = 0.0;
        m_du2fCombo[5] = 0.0;
        m_du2fCombo[6] = 0.0;

	// u3 explicit evaluation combination TO DO
	m_du3fCombo[0] = 1.0 - m_dExpCf[2][0] / m_dExpCf[0][0];
	m_du3fCombo[1] = m_dExpCf[2][0] / m_dExpCf[0][0] - 
                         m_dImpCf[2][0] / m_dImpCf[0][0];
	m_du3fCombo[2] = m_dImpCf[2][0] / m_dImpCf[0][0] - 
                         m_dImpCf[2][1] / m_dImpCf[0][1];
	m_du3fCombo[3] = m_dImpCf[2][1] / m_dImpCf[0][1];
	m_du3fCombo[4] = m_dExpCf[2][1] / m_dExpCf[1][1] - 
                         m_dImpCf[2][2] / m_dImpCf[1][2];
	m_du3fCombo[5] = m_dImpCf[2][2] / m_dImpCf[1][2];
        m_du3fCombo[6] = -m_dExpCf[2][1] / m_dExpCf[1][1];
        m_du3fCombo[7] = 0.0;

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

        // Compute ug1 into index 2 (NO IMPLICIT SOLVE NEEDED HERE)
        pGrid->CopyData(1, 2, DataType_State);
        pGrid->CopyData(1, 2, DataType_State);
        pVerticalDynamics->StepImplicitTermsExplicitly(
		1, 2, time, m_dImpCf[0][0] * dDeltaT);
//        pVerticalDynamics->StepImplicit(
//		0, 2, time, m_dImpCf[0][0] * dDeltaT);
        pGrid->PostProcessSubstage(2, DataType_State);
        pGrid->PostProcessSubstage(2, DataType_Tracers);

	// Compute u1 into index 3
	pGrid->CopyData(2, 3, DataType_State);
        pGrid->CopyData(2, 3, DataType_State);        
        pVerticalDynamics->StepImplicit(
		3, 3, time, m_dImpCf[0][1] * dDeltaT);
        pGrid->PostProcessSubstage(3, DataType_State);
	pGrid->PostProcessSubstage(3, DataType_Tracers);

	// STAGE 2
	// Compute uf2 from u1 (index 2) into index 6
	pGrid->LinearCombineData(m_du2fCombo, 6, DataType_State);
        pGrid->LinearCombineData(m_du2fCombo, 6, DataType_Tracers);
        pGrid->CopyData(6, 4, DataType_State);
        pGrid->CopyData(6, 4, DataType_Tracers);
	pHorizontalDynamics->StepExplicit(
		3, 4, time, m_dExpCf[1][1] * dDeltaT);
	pVerticalDynamics->StepExplicit(
		3, 4, time, m_dExpCf[1][1] * dDeltaT);
	pGrid->PostProcessSubstage(4, DataType_State);
	pGrid->PostProcessSubstage(4, DataType_Tracers);

	// Compute u2 from uf2 (index 3) into index 4
        pGrid->CopyData(4, 5, DataType_State);
        pGrid->CopyData(4, 5, DataType_State);
	pVerticalDynamics->StepImplicit(
		5, 5, time, m_dImpCf[1][2] * dDeltaT);
	pGrid->PostProcessSubstage(5, DataType_State);
	pGrid->PostProcessSubstage(5, DataType_Tracers);

	// STAGE 3
	// Compute uf3 from u2 (index 4) into index 8
	pGrid->LinearCombineData(m_du3fCombo, 7, DataType_State);
        pGrid->LinearCombineData(m_du3fCombo, 7, DataType_Tracers);
	pHorizontalDynamics->StepExplicit(
		5, 7, time, m_dExpCf[2][2] * dDeltaT);
	pVerticalDynamics->StepExplicit(
		5, 7, time, m_dExpCf[2][2] * dDeltaT);
	pGrid->PostProcessSubstage(7, DataType_State);
	pGrid->PostProcessSubstage(7, DataType_Tracers);

	// Apply hyperdiffusion at the end of the explicit substep (ask Paul)
	pGrid->CopyData(7, 2, DataType_State);
        pGrid->CopyData(7, 2, DataType_Tracers);
	pHorizontalDynamics->StepAfterSubCycle(7, 1, 3, time, dDeltaT);
	pGrid->CopyData(1, 0, DataType_State);
        pGrid->CopyData(1, 0, DataType_Tracers);
}

///////////////////////////////////////////////////////////////////////////////

