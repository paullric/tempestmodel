///////////////////////////////////////////////////////////////////////////////
///
///	\file    TimestepSchemeARK343.h
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

#include "TimestepSchemeARK343.h"
#include "Model.h"
#include "Grid.h"
#include "HorizontalDynamics.h"
#include "VerticalDynamics.h"

///////////////////////////////////////////////////////////////////////////////
// THESE COEFFICIENTS ARE COMPUTED FROM THE ORIGINAL TABLEAUX
// IMPLEMENTS ARK(3,4,3) FROM ASCHER ET AL. 1997 PG. 9
const double TimestepSchemeARK343::m_dgamma = 0.4358665215;
const double TimestepSchemeARK343::m_db1 = -1.5 * m_dgamma * m_dgamma + 
                                           4.0 * m_dgamma - 0.25;
const double TimestepSchemeARK343::m_db2 = 1.5 * m_dgamma * m_dgamma -
                                           5.0 * m_dgamma + 1.2;
// Time step coefficients
const double TimestepSchemeARK343::m_dbCf[3] = {m_db1, m_db2, m_dgamma};
// Implicit stage coefficients - Converted to U from Ascher et al. 1997 Pg 9
const double TimestepSchemeARK343::m_dImpCf[3][3] = {
  {m_dgamma, 0., 0.},
  {0.5 * (1.0 - m_dgamma), m_dgamma, 0.},
  {m_db1, m_db2, m_dgamma}};

// Explicit stage coefficients - Converted to U from Ascher et al. 1997 Pg 9
const double TimestepSchemeARK343::m_dExpCf[4][4] = {
  {0., 0., 0., 0.},
  {m_dgamma, 0., 0., 0.},
  {0.3212788860, 0.3966543747, 0., 0.},
  {-0.105858296, 0.5529291479, 0.5529291479, 0.}};

TimestepSchemeARK343::TimestepSchemeARK343(
	Model & model
) :
	TimestepScheme(model)
{
    m_dKh1Combo.Allocate(5);
    m_dKh2Combo.Allocate(7);
    m_dKh3Combo.Allocate(9);
    m_dK1Combo.Allocate(4);
    m_dK2Combo.Allocate(7);
    m_dK3Combo.Allocate(9);
    m_du2fCombo.Allocate(5);
    m_du3fCombo.Allocate(7);
    m_du4fCombo.Allocate(9);
}

///////////////////////////////////////////////////////////////////////////////

void TimestepSchemeARK343::Step(
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

	// u2 implicit evaluation combination
	m_du2fCombo[0] = 1.0;
	m_du2fCombo[1] = 0.0;
	m_du2fCombo[2] = 0.0;
	m_du2fCombo[3] = dDeltaT * m_dImpCf[1][0];
	m_du2fCombo[4] = dDeltaT * m_dExpCf[2][0];

	// K2 combination
	m_dK2Combo[0] = -1.0 / (dDeltaT * m_dImpCf[1][1]);
	m_dK2Combo[1] = 0.0;
	m_dK2Combo[2] = 1.0 / (dDeltaT * m_dImpCf[1][1]);
	m_dK2Combo[3] = -m_dImpCf[1][0] / m_dImpCf[1][1];
	m_dK2Combo[4] = -m_dExpCf[2][0] / m_dImpCf[1][1];
	m_dK2Combo[5] = 0.0;
	m_dK2Combo[6] = -m_dExpCf[2][1] / m_dImpCf[1][1];

	// Kh2 combination
	m_dKh2Combo[0] = -1.0 / (dDeltaT * m_dExpCf[2][1]);
	m_dKh2Combo[1] = 1.0 / (dDeltaT * m_dExpCf[2][1]);
	m_dKh2Combo[2] = 0.0;
	m_dKh2Combo[3] = -m_dImpCf[1][0] / m_dExpCf[2][1];
	m_dKh2Combo[4] = -m_dExpCf[2][0] / m_dExpCf[2][1];
	m_dKh2Combo[5] = 0.0;
	m_dKh2Combo[6] = 0.0;

	// u3 implicit evaluation combination
	m_du3fCombo[0] = 1.0;
	m_du3fCombo[1] = 0.0;
	m_du3fCombo[2] = 0.0;
	m_du3fCombo[3] = dDeltaT * m_dImpCf[2][0];
	m_du3fCombo[4] = dDeltaT * m_dExpCf[3][0];
	m_du3fCombo[5] = dDeltaT * m_dImpCf[2][1];
	m_du3fCombo[6] = dDeltaT * m_dExpCf[3][1];

	// K3 combination
	m_dK3Combo[0] = -1.0 / (dDeltaT * m_dImpCf[2][2]);
	m_dK3Combo[1] = 0.0;
	m_dK3Combo[2] = 1.0 / (dDeltaT * m_dImpCf[2][2]);
	m_dK3Combo[3] = -m_dImpCf[2][0] / m_dImpCf[2][2];
	m_dK3Combo[4] = -m_dExpCf[3][0] / m_dImpCf[2][2];
	m_dK3Combo[5] = -m_dImpCf[2][1] / m_dImpCf[2][2];
	m_dK3Combo[6] = -m_dExpCf[3][1] / m_dImpCf[2][2];
	m_dK3Combo[7] = 0.0;
	m_dK3Combo[8] = -m_dExpCf[3][2] / m_dImpCf[2][2];

	// Kh3 combination
	m_dKh3Combo[0] = -1.0 / (dDeltaT * m_dExpCf[3][2]);
	m_dKh3Combo[1] = 1.0 / (dDeltaT * m_dExpCf[3][2]);
	m_dKh3Combo[2] = 0.0;
	m_dKh3Combo[3] = -m_dImpCf[2][0] / m_dExpCf[3][2];
	m_dKh3Combo[4] = -m_dExpCf[3][0] / m_dExpCf[3][2];
	m_dKh3Combo[5] = -m_dImpCf[2][1] / m_dExpCf[3][2];
	m_dKh3Combo[6] = -m_dExpCf[3][1] / m_dExpCf[3][2];
	m_dKh3Combo[7] = 0.0;
	m_dKh3Combo[8] = 0.0;

	// uf4 explicit evaluation combination
	m_du4fCombo[0] = 1.0;
	m_du4fCombo[1] = 0.0;
	m_du4fCombo[2] = 0.0;
	m_du4fCombo[3] = dDeltaT * m_dbCf[0];
	m_du4fCombo[4] = 0.0;
	m_du4fCombo[5] = dDeltaT * m_dbCf[1];
	m_du4fCombo[6] = dDeltaT * m_dbCf[0];
	m_du4fCombo[7] = dDeltaT * m_dbCf[2];
	m_du4fCombo[8] = dDeltaT * m_dbCf[1];

	// SUBSTAGE 1
	// Compute the explicit step to index 1
	Time timeSub1 = time;
	double dtSub1 = m_dExpCf[1][0] * dDeltaT;
	pGrid->CopyData(0, 1, DataType_State);
        pGrid->CopyData(0, 1, DataType_Tracers);
	pHorizontalDynamics->StepExplicit(0, 1, timeSub1, dtSub1);
	pVerticalDynamics->StepExplicit(0, 1, timeSub1, dtSub1);
	pGrid->PostProcessSubstage(1, DataType_State);
	pGrid->PostProcessSubstage(1, DataType_Tracers);

	// Store the evaluation Kh1 to index 4
	pGrid->LinearCombineData(m_dKh1Combo, 4, DataType_State);
        pGrid->LinearCombineData(m_dKh1Combo, 4, DataType_Tracers);

	// Compute implicit step based on known data to index 2
	pGrid->CopyData(1, 2, DataType_State);
        pGrid->CopyData(1, 2, DataType_Tracers);
	dtSub1 = m_dImpCf[0][0] * dDeltaT;
	pVerticalDynamics->StepImplicit(1, 2, timeSub1, dtSub1);
	pGrid->PostProcessSubstage(2, DataType_State);
	pGrid->PostProcessSubstage(2, DataType_Tracers);

	// Store the evaluation K1 to index 3
	pGrid->LinearCombineData(m_dK1Combo, 3, DataType_State);

	//std::cout << "Substage 1 done ... \n";

	// SUBSTAGE 2
	// Compute uf2
	Time timeSub2 = time;
	double dtSub2 = m_dExpCf[2][1] * dDeltaT;
	pGrid->LinearCombineData(m_du2fCombo, 1, DataType_State);
        pGrid->LinearCombineData(m_du2fCombo, 1, DataType_Tracers);
	pHorizontalDynamics->StepExplicit(2, 1, timeSub2, dtSub2);
	pVerticalDynamics->StepExplicit(2, 1, timeSub2, dtSub2);
	pGrid->PostProcessSubstage(1, DataType_State);
	pGrid->PostProcessSubstage(1, DataType_Tracers);

	// Store the evaluation Kh2 to index 6
	pGrid->LinearCombineData(m_dKh2Combo, 6, DataType_State);

	// Compute u2 from uf2 and store it to index 2 (over u1)
	dtSub2 = m_dImpCf[1][1] * dDeltaT;
	pGrid->CopyData(1, 2, DataType_State);
        pGrid->CopyData(1, 2, DataType_Tracers);
	pVerticalDynamics->StepImplicit(1, 2, timeSub2, dtSub2);
	pGrid->PostProcessSubstage(2, DataType_State);
	pGrid->PostProcessSubstage(2, DataType_Tracers);

	// Store the evaluation K2 to index 5
	pGrid->LinearCombineData(m_dK2Combo, 5, DataType_State);
        pGrid->LinearCombineData(m_dK2Combo, 5, DataType_Tracers);

	//std::cout << "Substage 2 done ... \n";

	// SUBSTAGE 3
	// Compute uf3
	Time timeSub3 = time;
	double dtSub3 = m_dExpCf[3][2] * dDeltaT;
	pGrid->LinearCombineData(m_du3fCombo, 1, DataType_State);
        pGrid->LinearCombineData(m_du3fCombo, 1, DataType_Tracers);
	pHorizontalDynamics->StepExplicit(2, 1, timeSub2, dtSub2);
	pVerticalDynamics->StepExplicit(2, 1, timeSub2, dtSub2);
	pGrid->PostProcessSubstage(1, DataType_State);
	pGrid->PostProcessSubstage(1, DataType_Tracers);

	// Store the evaluation Kh3 to index 8
	pGrid->LinearCombineData(m_dKh3Combo, 8, DataType_State);
        pGrid->LinearCombineData(m_dKh3Combo, 8, DataType_Tracers);

	// Compute u3 from uf3 and store it to index 2 (over u2)
	dtSub3 = m_dImpCf[2][2] * dDeltaT;
	pGrid->CopyData(1, 2, DataType_State);
        pGrid->CopyData(1, 2, DataType_Tracers);
	pVerticalDynamics->StepImplicit(1, 2, timeSub3, dtSub3);
	pGrid->PostProcessSubstage(2, DataType_State);
	pGrid->PostProcessSubstage(2, DataType_Tracers);

	// Store the evaluation K3 to index 7
	pGrid->LinearCombineData(m_dK3Combo, 7, DataType_State);
        pGrid->LinearCombineData(m_dK3Combo, 7, DataType_Tracers);

	//std::cout << "Substage 3 done ... \n";

	// FINAL STAGE
	// Compute uf4 from u3 and store it to index 2
	Time timeSub4 = time;
	double dtSub4 = m_dbCf[2] * dDeltaT;
	pGrid->LinearCombineData(m_du4fCombo, 1, DataType_State);
        pGrid->LinearCombineData(m_du4fCombo, 1, DataType_Tracers);
	pHorizontalDynamics->StepExplicit(2, 1, timeSub4, dtSub4);
	pVerticalDynamics->StepExplicit(2, 1, timeSub4, dtSub4);
	pGrid->PostProcessSubstage(1, DataType_State);
	pGrid->PostProcessSubstage(1, DataType_Tracers);

	// Apply hyperdiffusion at the end of the explicit substep (ask Paul)
	pHorizontalDynamics->StepAfterSubCycle(2, 1, 3, time, dDeltaT);
	pGrid->CopyData(1, 0, DataType_State);
        pGrid->CopyData(1, 0, DataType_Tracers);
}

///////////////////////////////////////////////////////////////////////////////

