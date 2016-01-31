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

///////////////////////////////////////////////////////////////////////////////
// THESE COEFFICIENTS ARE COMPUTED FROM THE ORIGINAL TABLEAUX
// IMPLEMENTS ARK(4,4,3) FROM ASCHER ET AL. 1997 PG. 9
// Time step coefficients
const double TimestepSchemeARS443::m_dTimeCf[4] = {1./2., 2./3., 1./2., 1.};
// Implicit stage coefficients - Converted to U from Ascher et al. 1997 Pg 9
const double TimestepSchemeARS443::m_dImpCf[4][4] = {
  {1./2., 0., 0., 0.},
  {1./6., 1./2., 0., 0.},
  {-1./2., 1./2., 1./2., 0.},
  {3./2., -3./2., 1./2., 1./2.}};

// Explicit stage coefficients - Converted to U from Ascher et al. 1997 Pg 9
const double TimestepSchemeARS443::m_dExpCf[4][4] = {
  {1./2., 0., 0., 0.},
  {11./18., 1./18., 0., 0.},
  {5./6., -5./6., 1./2., 0.},
  {1./4., 7./4., 3./4., -7./4.}};

TimestepSchemeARS443::TimestepSchemeARS443(
	Model & model
) :
	TimestepScheme(model)
{
    m_du2fCombo.Allocate(8);
    m_du3fCombo.Allocate(9);
    m_du4fCombo.Allocate(10);
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
        pGrid->CopyData(1, 2, DataType_State);        
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
        pGrid->CopyData(3, 4, DataType_State);
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
        pGrid->CopyData(5, 6, DataType_State);
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

	// Compute u2 from uf2 (index 7) into index 2
	pVerticalDynamics->StepImplicit(
		9, 9, time, m_dImpCf[3][3] * dDeltaT);
	pGrid->PostProcessSubstage(9, DataType_State);
	pGrid->PostProcessSubstage(9, DataType_Tracers);

	// Apply hyperdiffusion at the end of the explicit substep (ask Paul)
	pGrid->CopyData(9, 2, DataType_State);
        pGrid->CopyData(9, 2, DataType_Tracers);
	pHorizontalDynamics->StepAfterSubCycle(2, 1, 9, time, dDeltaT);
	pGrid->CopyData(1, 0, DataType_State);
        pGrid->CopyData(1, 0, DataType_Tracers);
}

///////////////////////////////////////////////////////////////////////////////

