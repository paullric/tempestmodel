///////////////////////////////////////////////////////////////////////////////
///
///	\file    TimestepSchemeHARK332.h
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

#include "TimestepSchemeHARK332.h"
#include "Model.h"
#include "Grid.h"
#include "HorizontalDynamics.h"
#include "VerticalDynamics.h"

///////////////////////////////////////////////////////////////////////////////
// IMPLEMENTATION OF SSP2(3,3,2) FROM Kupka et. al. 2012 (Scheme by Higueras)
// Weight coefficients (THESE ARE IDENTICAL FOR BOTH METHODS)
const double TimestepSchemeHARK332::m_dbCf[3] = {1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0};
// Implicit stage coefficients
const double TimestepSchemeHARK332::m_dImpCf[3][3] = {
  {0.25, 0.0, 0.0},
  {0.0, 0.25, 0.0},
  {1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0}};

// Explicit stage coefficients
const double TimestepSchemeHARK332::m_dExpCf[3][3] = {
  {0.0, 0.0, 0.0},
  {0.5, 0.0, 0.0},
  {0.5, 0.5, 0.0}};

TimestepSchemeHARK332::TimestepSchemeHARK332(
	Model & model
) :
	TimestepScheme(model)
{
	m_dKv1Combo.Allocate(8);
        m_dKv2Combo.Allocate(8);
        m_dKv3Combo.Allocate(8);
        m_dKh1Combo.Allocate(8);
        m_dKh2Combo.Allocate(8);
	m_du2Combo.Allocate(8);
	m_du3Combo.Allocate(8);
        m_dufCombo.Allocate(8);
}

///////////////////////////////////////////////////////////////////////////////

void TimestepSchemeHARK332::Step(
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

        // Kv1 combination
	m_dKv1Combo[0] = -1.0 / (dDeltaT * m_dImpCf[0][0]);
	m_dKv1Combo[1] = 1.0 / (dDeltaT * m_dImpCf[0][0]);
	m_dKv1Combo[2] = 0.0;
	m_dKv1Combo[3] = 0.0;
	m_dKv1Combo[4] = 0.0;
        m_dKv2Combo[5] = 0.0;
	m_dKv1Combo[6] = 0.0;
        m_dKv2Combo[7] = 0.0;

        // Kv2 combination
	m_dKv2Combo[0] = -1.0 / (dDeltaT * m_dImpCf[1][1]);
	m_dKv2Combo[1] = 0.0;
	m_dKv2Combo[2] = 1.0 / (dDeltaT * m_dImpCf[1][1]);
	m_dKv2Combo[3] = -m_dImpCf[1][0] / m_dImpCf[1][1];
	m_dKv2Combo[4] = -m_dExpCf[1][0] / m_dImpCf[1][1];
        m_dKv2Combo[5] = 0.0;
        m_dKv2Combo[6] = 0.0;
        m_dKv2Combo[7] = 0.0;

        // Kv3 combination
	m_dKv3Combo[0] = -1.0 / (dDeltaT * m_dImpCf[2][2]);
        m_dKv3Combo[1] = 1.0 / (dDeltaT * m_dImpCf[2][2]);
	m_dKv3Combo[2] = 0.0;
	m_dKv3Combo[3] = -m_dImpCf[2][0] / m_dImpCf[2][2];
	m_dKv3Combo[4] = -m_dExpCf[2][0] / m_dImpCf[2][2];
        m_dKv3Combo[5] = -m_dImpCf[2][1] / m_dImpCf[2][2];
        m_dKv3Combo[6] = -m_dExpCf[2][1] / m_dImpCf[2][2];
        m_dKv3Combo[7] = 0.0;

        // Kh1 combination
	m_dKh1Combo[0] = -1.0 / (dDeltaT * 0.1);
	m_dKh1Combo[1] = 0.0;
	m_dKh1Combo[2] = 1.0 / (dDeltaT * 0.1);
	m_dKh1Combo[3] = 0.0;
	m_dKh1Combo[4] = 0.0;
        m_dKh1Combo[5] = 0.0;
	m_dKh1Combo[6] = 0.0;
        m_dKh1Combo[7] = 0.0;

        // Kh2 combination
	m_dKh2Combo[0] = -1.0 / (dDeltaT * m_dExpCf[2][1]);
        m_dKh2Combo[1] = 1.0 / (dDeltaT * m_dExpCf[2][1]);
	m_dKh2Combo[2] = 0.0;
	m_dKh2Combo[3] = -m_dImpCf[2][0] / m_dExpCf[2][1];
	m_dKh2Combo[4] = -m_dExpCf[2][0] / m_dExpCf[2][1];
        m_dKh2Combo[5] = -m_dImpCf[2][1] / m_dExpCf[2][1];
        m_dKh2Combo[6] = 0.0;
        m_dKh2Combo[7] = 0.0;

	// u2 evaluation combination
	m_du2Combo[0] = 1.0;
	m_du2Combo[1] = 0.0;
	m_du2Combo[2] = 0.0;
	m_du2Combo[3] = dDeltaT * m_dImpCf[1][0];
	m_du2Combo[4] = dDeltaT * m_dExpCf[1][0];
        m_du2Combo[5] = 0.0;
	m_du2Combo[6] = 0.0;
        m_du2Combo[7] = 0.0;

        // u3 evaluation combination
	m_du3Combo[0] = 1.0;
	m_du3Combo[1] = 0.0;
	m_du3Combo[2] = 0.0;
	m_du3Combo[3] = dDeltaT * m_dImpCf[2][0];
	m_du3Combo[4] = dDeltaT * m_dExpCf[2][0];
        m_du3Combo[5] = dDeltaT * m_dImpCf[2][1];
	m_du3Combo[6] = 0.0;
        m_du3Combo[7] = 0.0;

	// uf evaluation combination
	m_dufCombo[0] = 1.0;
	m_dufCombo[1] = 0.0;
	m_dufCombo[2] = 0.0;
	m_dufCombo[3] = dDeltaT * m_dbCf[0];
	m_dufCombo[4] = dDeltaT * m_dbCf[0];
        m_dufCombo[5] = dDeltaT * m_dbCf[1];
	m_dufCombo[6] = dDeltaT * m_dbCf[1];
        m_dufCombo[7] = dDeltaT * m_dbCf[2];

	//std::cout << "Entering substages at the timestep... \n";

	// STAGE 1
	// Compute u1
        pGrid->CopyData(0, 1, DataType_State);
        pGrid->CopyData(0, 1, DataType_Tracers);
        pVerticalDynamics->StepImplicit(
		1, 1, time, m_dImpCf[0][0] * dDeltaT);
        pGrid->PostProcessSubstage(1, DataType_State);
	pGrid->PostProcessSubstage(1, DataType_Tracers);

	// Store the evaluation Kv1 to index 3
	pGrid->LinearCombineData(m_dKv1Combo, 3, DataType_State);
        pGrid->LinearCombineData(m_dKv1Combo, 3, DataType_Tracers);

        // This is needed to make the first stage explicit evaluation
        pGrid->CopyData(0, 2, DataType_State);
        pGrid->CopyData(0, 2, DataType_Tracers);
        pHorizontalDynamics->StepExplicit(
		0, 2, time, 0.1 * dDeltaT);
	pVerticalDynamics->StepExplicit(
		0, 2, time, 0.1 * dDeltaT);
	pGrid->PostProcessSubstage(2, DataType_State);
	pGrid->PostProcessSubstage(2, DataType_Tracers);

        // Store the evaluation Kh1 to index 4
	pGrid->LinearCombineData(m_dKh1Combo, 4, DataType_State);
        pGrid->LinearCombineData(m_dKh1Combo, 4, DataType_Tracers);

        // STAGE 2
	// Compute u2
        pGrid->LinearCombineData(m_du2Combo, 2, DataType_State);
        pGrid->LinearCombineData(m_du2Combo, 2, DataType_Tracers);

        // Last implicit solve in stage 2
        pVerticalDynamics->StepImplicit(
		2, 2, time, m_dImpCf[1][1] * dDeltaT);
        pGrid->PostProcessSubstage(2, DataType_State);
	pGrid->PostProcessSubstage(2, DataType_Tracers);

        // Store the evaluation Kv2 to index 5
	pGrid->LinearCombineData(m_dKv2Combo, 5, DataType_State);
        pGrid->LinearCombineData(m_dKv2Combo, 5, DataType_Tracers);

        // STAGE 3
        // Compute u3
        pGrid->LinearCombineData(m_du3Combo, 2, DataType_State);
        pGrid->LinearCombineData(m_du3Combo, 2, DataType_Tracers);
        pGrid->CopyData(2, 1, DataType_State);
        pGrid->CopyData(2, 1, DataType_Tracers);
	pHorizontalDynamics->StepExplicit(
		2, 1, time, m_dExpCf[2][1] * dDeltaT);
	pVerticalDynamics->StepExplicit(
		2, 1, time, m_dExpCf[2][1] * dDeltaT);
	pGrid->PostProcessSubstage(1, DataType_State);
	pGrid->PostProcessSubstage(1, DataType_Tracers);

        // Store the evaluation Kh2 to index 6
	pGrid->LinearCombineData(m_dKh2Combo, 6, DataType_State);
        pGrid->LinearCombineData(m_dKh2Combo, 6, DataType_Tracers);

        // Implicit solve in stage 2
        pVerticalDynamics->StepImplicit(
		1, 1, time, m_dImpCf[2][2] * dDeltaT);

        // Store the evaluation Kv2 to index 7
	pGrid->LinearCombineData(m_dKv3Combo, 7, DataType_State);
        pGrid->LinearCombineData(m_dKv3Combo, 7, DataType_Tracers);

	// FINAL STAGE
	pGrid->LinearCombineData(m_dufCombo, 1, DataType_State);
        pGrid->LinearCombineData(m_dufCombo, 1, DataType_Tracers);
        pGrid->CopyData(1, 2, DataType_State);
        pGrid->CopyData(1, 2, DataType_Tracers);
        pHorizontalDynamics->StepExplicit(
		1, 2, time, m_dbCf[1] * dDeltaT);
	pVerticalDynamics->StepExplicit(
		1, 2, time, m_dbCf[1] * dDeltaT);
	pGrid->PostProcessSubstage(2, DataType_State);
	pGrid->PostProcessSubstage(2, DataType_Tracers);

	// Apply hyperdiffusion at the end of the explicit substep (ask Paul)
        pGrid->CopyData(1, 2, DataType_State);
        pGrid->CopyData(1, 2, DataType_Tracers);
	pHorizontalDynamics->StepAfterSubCycle(2, 1, 3, time, dDeltaT);
	pGrid->CopyData(1, 0, DataType_State);
        pGrid->CopyData(1, 0, DataType_Tracers);
}

///////////////////////////////////////////////////////////////////////////////


