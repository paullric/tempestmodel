///
///	\file    TimestepSchemeARK3.cpp
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

#include "TimestepSchemeARK3.h"
#include "Model.h"
#include "Grid.h"
#include "HorizontalDynamics.h"
#include "VerticalDynamics.h"

#pragma message "File tagged for cleanup"

///////////////////////////////////////////////////////////////////////////////
// THESE COEFFICIENTS ARE COMPUTED FROM THE ORIGINAL TABLEAUX
// IMPLEMENTS ARK(4,4,3) FROM ASCHER ET AL. 1997 PG. 9
// Time step coefficients
const double TimestepSchemeARK3::m_dTimeCf[4] = {1./2., 2./3., 1./2., 1.};
// Implicit stage coefficients - Converted to U from Ascher et al. 1997 Pg 9
const double TimestepSchemeARK3::m_dImpCf[4][4] = {
  {1./2., 0., 0., 0.},
  {1./6., 1./2., 0., 0.},
  {-1./2., 1./2., 1./2., 0.},
  {3./2., -3./2., 1./2., 1./2.}};

// Explicit stage coefficients - Converted to U from Ascher et al. 1997 Pg 9
const double TimestepSchemeARK3::m_dExpCf[5][5] = {
  {0., 0., 0., 0., 0.},
  {1./2., 0., 0., 0., 0.},
  {11./18., 1./18., 0., 0., 0.},
  {5./6., -5./6., 1./2., 0., 0.},
  {1./4., 7./4., 3./4., -7./4., 0.}};

TimestepSchemeARK3::TimestepSchemeARK3(
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

void TimestepSchemeARK3::Step(
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
    m_du4fCombo[3] = dDeltaT * m_dImpCf[3][0];
    m_du4fCombo[4] = dDeltaT * m_dExpCf[4][0];
    m_du4fCombo[5] = dDeltaT * m_dImpCf[3][1];
    m_du4fCombo[6] = dDeltaT * m_dExpCf[4][1];
    m_du4fCombo[7] = dDeltaT * m_dImpCf[3][2];
    m_du4fCombo[8] = dDeltaT * m_dExpCf[4][2];

    // SUBSTAGE 1
    // Compute the explicit step to index 1
    Time timeSub1 = time;// + m_dTimeCf[0] * dDeltaT;
    double dtSub1 = m_dExpCf[1][0] * dDeltaT;
    pGrid->CopyData(0, 1, DataType_State);
    pHorizontalDynamics->StepExplicit(0, 1, timeSub1, dtSub1);
	pVerticalDynamics->StepExplicit(0, 1, timeSub1, dtSub1);
    pGrid->PostProcessSubstage(1, DataType_State);
    pGrid->PostProcessSubstage(1, DataType_Tracers);

    // Store the evaluation Kh1 to index 4
    pGrid->LinearCombineData(m_dKh1Combo, 4, DataType_State);

    // Compute implicit step based on known data to index 2
    pGrid->CopyData(1, 2, DataType_State);
    dtSub1 = m_dImpCf[0][0] * dDeltaT;
    pVerticalDynamics->StepImplicit(1, 2, timeSub1, dtSub1);
    pGrid->PostProcessSubstage(2, DataType_State);
    pGrid->PostProcessSubstage(2, DataType_Tracers);

    // Store the evaluation K1 to index 3
    pGrid->LinearCombineData(m_dK1Combo, 3, DataType_State);

    //std::cout << "Substage 1 done ... \n";

    // SUBSTAGE 2
    // Compute uf2
    Time timeSub2 = time; // + m_dTimeCf[1] * dDeltaT;
    double dtSub2 = m_dExpCf[2][1] * dDeltaT;
    pGrid->LinearCombineData(m_du2fCombo, 1, DataType_State);
    pHorizontalDynamics->StepExplicit(2, 1, timeSub2, dtSub2);
    pVerticalDynamics->StepExplicit(2, 1, timeSub2, dtSub2);
    pGrid->PostProcessSubstage(1, DataType_State);
    pGrid->PostProcessSubstage(1, DataType_Tracers);

    // Store the evaluation Kh2 to index 6
    pGrid->LinearCombineData(m_dKh2Combo, 6, DataType_State);

    // Compute u2 from uf2 and store it to index 2 (over u1)
    dtSub2 = m_dImpCf[1][1] * dDeltaT;
    pGrid->CopyData(1, 2, DataType_State);
    pVerticalDynamics->StepImplicit(1, 2, timeSub2, dtSub2);
    pGrid->PostProcessSubstage(2, DataType_State);
    pGrid->PostProcessSubstage(2, DataType_Tracers);

    // Store the evaluation K2 to index 5
    pGrid->LinearCombineData(m_dK2Combo, 5, DataType_State);

    //std::cout << "Substage 2 done ... \n";

    // SUBSTAGE 3
    // Compute uf3
    Time timeSub3 = time; // + m_dTimeCf[2] * dDeltaT;
    double dtSub3 = m_dExpCf[3][2] * dDeltaT;
    pGrid->LinearCombineData(m_du3fCombo, 1, DataType_State);
    pHorizontalDynamics->StepExplicit(2, 1, timeSub2, dtSub2);
    pVerticalDynamics->StepExplicit(2, 1, timeSub2, dtSub2);
    pGrid->PostProcessSubstage(1, DataType_State);
    pGrid->PostProcessSubstage(1, DataType_Tracers);

    // Store the evaluation Kh3 to index 8
    pGrid->LinearCombineData(m_dKh3Combo, 8, DataType_State);

    // Compute u3 from uf3 and store it to index 2 (over u2)
    dtSub3 = m_dImpCf[2][2] * dDeltaT;
    pGrid->CopyData(1, 2, DataType_State);
    pVerticalDynamics->StepImplicit(1, 2, timeSub3, dtSub3);
    pGrid->PostProcessSubstage(2, DataType_State);
    pGrid->PostProcessSubstage(2, DataType_Tracers);

    // Store the evaluation K3 to index 7
    pGrid->LinearCombineData(m_dK3Combo, 7, DataType_State);

    //std::cout << "Substage 3 done ... \n";

    // SUBSTAGE 4
    // Compute uf4 from u3 and store it to index 2
    Time timeSub4 = time; // + m_dTimeCf[3] * dDeltaT;
    double dtSub4 = m_dExpCf[4][3] * dDeltaT;
    //pGrid->CopyData(2, 1, DataType_State);
    pGrid->LinearCombineData(m_du4fCombo, 1, DataType_State);
    pHorizontalDynamics->StepExplicit(2, 1, timeSub4, dtSub4);
    pVerticalDynamics->StepExplicit(2, 1, timeSub4, dtSub4);
    pGrid->PostProcessSubstage(1, DataType_State);
    pGrid->PostProcessSubstage(1, DataType_Tracers);

    // Compute u4 from uf4 and store it to index 2 (over u3)
    dtSub4 = m_dImpCf[3][3] * dDeltaT;
    pGrid->CopyData(1, 2, DataType_State);
    pVerticalDynamics->StepImplicit(1, 2, timeSub4, dtSub4);
    pGrid->PostProcessSubstage(2, DataType_State);
    pGrid->PostProcessSubstage(2, DataType_Tracers);

    // Apply hyperdiffusion at the end of the explicit substep (ask Paul)
	pGrid->CopyData(2, 1, DataType_State);
	pHorizontalDynamics->StepAfterSubCycle(2, 1, 3, time, dDeltaT);
	pGrid->CopyData(1, 0, DataType_State);
	//*/
}

///////////////////////////////////////////////////////////////////////////////

