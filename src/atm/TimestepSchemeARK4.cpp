///
///	\file    TimestepSchemeARK4.cpp
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

#include "TimestepSchemeARK4.h"
#include "Model.h"
#include "Grid.h"
#include "HorizontalDynamics.h"
#include "VerticalDynamics.h"

#pragma message "File tagged for cleanup"

///////////////////////////////////////////////////////////////////////////////
// THESE COEFFICIENTS ARE COMPUTED FROM THE ORIGINAL TABLEAUX
// IMPLEMENTS RK.4.A.1 from Liu et al. 2006 Pg 87
// HAS A ZERO EXPLICIT EVALUATION AT THE 4TH STAGE - HAD TO MODIFY - JEG
// Time step coefficients
const double TimestepSchemeARK4::m_dTimeCf[6] = {1./3., 1./3., 
												 1./2., 1./2., 1., 1.};
// Implicit stage coefficients
const double TimestepSchemeARK4::m_dImpCf[7][7] = {
  {1./10., 0., 0., 0., 0., 0., 0.},
  {-1./6., 1./2., 0., 0., 0., 0., 0.},
  {1./6., -1./3., 1./2., 0., 0., 0., 0.},
  {3./8., -3./8., 0., 1./2., 0., 0., 0.},
  {1./8., 0., 3./8., -1./2., 1./2., 0., 0.},
  {-1./2., 0., 3., -2., 0., 1./2., 0.},
  {1./6., 0., 0., 0., 2./3., -1./2., 2./3.}};

// Explicit stage coefficients
const double TimestepSchemeARK4::m_dExpCf[7][7] = {
  {0., 0., 0., 0., 0., 0., 0.},
  {1./3., 0., 0., 0., 0., 0., 0.},
  {1./6., 1./6., 0., 0., 0., 0., 0.},
  {1./8., 0., 3./8., 0., 0., 0., 0.},
  {1./8., 0., 3./8., 1., 0., 0., 0.},
  {1./2., 0., -3./2., 1., 1., 0., 0.},
  {1./6., 0., 0., 0., 2./3., 1./6., 0.}};

TimestepSchemeARK4::TimestepSchemeARK4(
	Model & model
) :
	TimestepScheme(model)
{
    m_dKh1Combo.Allocate(5);
    m_dKh2Combo.Allocate(7);
    m_dKh3Combo.Allocate(9);
	m_dKh4Combo.Allocate(11);
	m_dKh5Combo.Allocate(13);
	m_dK0Combo.Allocate(4);
    m_dK1Combo.Allocate(6);
    m_dK2Combo.Allocate(8);
    m_dK3Combo.Allocate(10);
	m_dK4Combo.Allocate(12);
	m_dK5Combo.Allocate(14);
	m_du1fCombo.Allocate(4);
    m_du2fCombo.Allocate(6);
    m_du3fCombo.Allocate(8);
    m_du4fCombo.Allocate(10);
	m_du5fCombo.Allocate(12);
	m_du6fCombo.Allocate(14);
}

///////////////////////////////////////////////////////////////////////////////

void TimestepSchemeARK4::Step(
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
	
	// K0 combination
	m_dK0Combo[0] = -1.0 / (dDeltaT * m_dImpCf[0][0]);;
    m_dK0Combo[1] = 1.0 / (dDeltaT * m_dImpCf[0][0]);;
    m_dK0Combo[2] = 0.0;
	m_dK0Combo[3] = 0.0;

	// u1 implicit evaluation combination
    m_du1fCombo[0] = 1.0;
    m_du1fCombo[1] = 0.0;
    m_du1fCombo[2] = 0.0;
    m_du1fCombo[3] = dDeltaT * m_dImpCf[1][0];

    // Kh1 combination
    m_dKh1Combo[0] = -1.0 / (dDeltaT * m_dExpCf[1][0]);
    m_dKh1Combo[1] = 1.0 / (dDeltaT * m_dExpCf[1][0]);
	m_dKh1Combo[2] = 0.0;
	m_dKh1Combo[3] = -m_dImpCf[0][0] / m_dExpCf[1][0];
	m_dKh1Combo[4] = 0.0;

    // K1 combination
    m_dK1Combo[0] = -1.0 / (dDeltaT * m_dImpCf[1][1]);
    m_dK1Combo[1] = 0.0;
    m_dK1Combo[2] = 1.0 / (dDeltaT * m_dImpCf[1][1]);
	m_dK1Combo[3] = -m_dImpCf[0][0] / m_dImpCf[1][1];
	m_dK1Combo[4] = -m_dExpCf[1][0] / m_dImpCf[1][1];
	m_dK1Combo[5] = 0.0;

    // u2 implicit evaluation combination
    m_du2fCombo[0] = 1.0;
    m_du2fCombo[1] = 0.0;
    m_du2fCombo[2] = 0.0;
    m_du2fCombo[3] = dDeltaT * m_dImpCf[2][0];
    m_du2fCombo[4] = dDeltaT * m_dExpCf[2][0];
	m_du2fCombo[5] = dDeltaT * m_dImpCf[2][1];

    // K2 combination
    m_dK2Combo[0] = -1.0 / (dDeltaT * m_dImpCf[2][2]);
    m_dK2Combo[1] = 0.0;
    m_dK2Combo[2] = 1.0 / (dDeltaT * m_dImpCf[2][2]);
    m_dK2Combo[3] = -m_dImpCf[1][0] / m_dImpCf[2][2];
    m_dK2Combo[4] = -m_dExpCf[2][0] / m_dImpCf[2][2];
    m_dK2Combo[5] = -m_dImpCf[2][1] / m_dImpCf[2][2];
    m_dK2Combo[6] = -m_dExpCf[2][1] / m_dImpCf[2][2];
	m_dK2Combo[7] = 0.0;

    // Kh2 combination
    m_dKh2Combo[0] = -1.0 / (dDeltaT * m_dExpCf[2][1]);
    m_dKh2Combo[1] = 1.0 / (dDeltaT * m_dExpCf[2][1]);
    m_dKh2Combo[2] = 0.0;
    m_dKh2Combo[3] = -m_dImpCf[2][0] / m_dExpCf[2][1];
    m_dKh2Combo[4] = -m_dExpCf[2][0] / m_dExpCf[2][1];
	m_dKh2Combo[5] = -m_dImpCf[2][1] / m_dExpCf[2][1];
	m_dKh2Combo[6] = 0.0;

    // u3 implicit evaluation combination
    m_du3fCombo[0] = 1.0;
    m_du3fCombo[1] = 0.0;
    m_du3fCombo[2] = 0.0;
    m_du3fCombo[3] = dDeltaT * m_dImpCf[3][0];
    m_du3fCombo[4] = dDeltaT * m_dExpCf[3][0];
    m_du3fCombo[5] = dDeltaT * m_dImpCf[3][1];
    m_du3fCombo[6] = dDeltaT * m_dExpCf[3][1];
	m_du3fCombo[7] = dDeltaT * m_dImpCf[3][2];

    // K3 combination
    m_dK3Combo[0] = -1.0 / (dDeltaT * m_dImpCf[3][3]);
    m_dK3Combo[1] = 0.0;
    m_dK3Combo[2] = 1.0 / (dDeltaT * m_dImpCf[3][3]);
    m_dK3Combo[3] = -m_dImpCf[3][0] / m_dImpCf[3][3];
    m_dK3Combo[4] = -m_dExpCf[3][0] / m_dImpCf[3][3];
    m_dK3Combo[5] = -m_dImpCf[3][1] / m_dImpCf[3][3];
    m_dK3Combo[6] = -m_dExpCf[3][1] / m_dImpCf[3][3];
    m_dK3Combo[7] = -m_dImpCf[3][2] / m_dImpCf[3][3];
    m_dK3Combo[8] = -m_dExpCf[3][2] / m_dImpCf[3][3];
	m_dK3Combo[9] = 0.0;

    // Kh3 combination
    m_dKh3Combo[0] = -1.0 / (dDeltaT * m_dExpCf[3][2]);
    m_dKh3Combo[1] = 1.0 / (dDeltaT * m_dExpCf[3][2]);
    m_dKh3Combo[2] = 0.0;
    m_dKh3Combo[3] = -m_dImpCf[3][0] / m_dExpCf[3][2];
    m_dKh3Combo[4] = -m_dExpCf[3][0] / m_dExpCf[3][2];
    m_dKh3Combo[5] = -m_dImpCf[3][1] / m_dExpCf[3][2];
    m_dKh3Combo[6] = -m_dExpCf[3][1] / m_dExpCf[3][2];
	m_dKh3Combo[7] = -m_dImpCf[3][2] / m_dExpCf[3][2];
	m_dKh3Combo[8] = 0.0;

    // uf4 explicit evaluation combination
    m_du4fCombo[0] = 1.0;
    m_du4fCombo[1] = 0.0;
    m_du4fCombo[2] = 0.0;
    m_du4fCombo[3] = dDeltaT * m_dImpCf[4][0];
    m_du4fCombo[4] = dDeltaT * m_dExpCf[4][0];
    m_du4fCombo[5] = dDeltaT * m_dImpCf[4][1];
    m_du4fCombo[6] = dDeltaT * m_dExpCf[4][1];
    m_du4fCombo[7] = dDeltaT * m_dImpCf[4][2];
    m_du4fCombo[8] = dDeltaT * m_dExpCf[4][2];
	m_du4fCombo[9] = dDeltaT * m_dImpCf[4][3];

	// K4 combination
    m_dK4Combo[0] = -1.0 / (dDeltaT * m_dImpCf[4][4]);
    m_dK4Combo[1] = 0.0;
    m_dK4Combo[2] = 1.0 / (dDeltaT * m_dImpCf[4][4]);
    m_dK4Combo[3] = -m_dImpCf[4][0] / m_dImpCf[4][4];
    m_dK4Combo[4] = -m_dExpCf[4][0] / m_dImpCf[4][4];
    m_dK4Combo[5] = -m_dImpCf[4][1] / m_dImpCf[4][4];
    m_dK4Combo[6] = -m_dExpCf[4][1] / m_dImpCf[4][4];
    m_dK4Combo[7] = -m_dImpCf[4][2] / m_dImpCf[4][4];
    m_dK4Combo[8] = -m_dExpCf[4][2] / m_dImpCf[4][4];
	m_dK4Combo[9] = -m_dImpCf[4][3] / m_dImpCf[4][4];
    m_dK4Combo[10] = -m_dExpCf[4][3] / m_dImpCf[4][4];
	m_dK4Combo[11] = 0.0;

    // Kh4 combination
    m_dKh4Combo[0] = -1.0 / (dDeltaT * m_dExpCf[4][3]);
    m_dKh4Combo[1] = 1.0 / (dDeltaT * m_dExpCf[4][3]);
    m_dKh4Combo[2] = 0.0;
    m_dKh4Combo[3] = -m_dImpCf[4][0] / m_dExpCf[4][3];
    m_dKh4Combo[4] = -m_dExpCf[4][0] / m_dExpCf[4][3];
    m_dKh4Combo[5] = -m_dImpCf[4][1] / m_dExpCf[4][3];
    m_dKh4Combo[6] = -m_dExpCf[4][1] / m_dExpCf[4][3];
	m_dKh4Combo[7] = -m_dImpCf[4][2] / m_dExpCf[4][3];
	m_dKh4Combo[8] = -m_dExpCf[4][2] / m_dExpCf[4][3];
	m_dKh4Combo[9] = -m_dImpCf[4][3] / m_dExpCf[4][3];
	m_dKh4Combo[10] = 0.0;

	// uf5 explicit evaluation combination
    m_du5fCombo[0] = 1.0;
    m_du5fCombo[1] = 0.0;
    m_du5fCombo[2] = 0.0;
    m_du5fCombo[3] = dDeltaT * m_dImpCf[5][0];
    m_du5fCombo[4] = dDeltaT * m_dExpCf[5][0];
    m_du5fCombo[5] = dDeltaT * m_dImpCf[5][1];
    m_du5fCombo[6] = dDeltaT * m_dExpCf[5][1];
    m_du5fCombo[7] = dDeltaT * m_dImpCf[5][2];
    m_du5fCombo[8] = dDeltaT * m_dExpCf[5][2];
	m_du5fCombo[9] = dDeltaT * m_dImpCf[5][3];
	m_du5fCombo[10] = dDeltaT * m_dExpCf[5][3];
	m_du5fCombo[11] = dDeltaT * m_dImpCf[5][4];

	// K5 combination
    m_dK5Combo[0] = -1.0 / (dDeltaT * m_dImpCf[5][5]);
    m_dK5Combo[1] = 0.0;
    m_dK5Combo[2] = 1.0 / (dDeltaT * m_dImpCf[5][5]);
    m_dK5Combo[3] = -m_dImpCf[5][0] / m_dImpCf[5][5];
    m_dK5Combo[4] = -m_dExpCf[5][0] / m_dImpCf[5][5];
    m_dK5Combo[5] = -m_dImpCf[5][1] / m_dImpCf[5][5];
    m_dK5Combo[6] = -m_dExpCf[5][1] / m_dImpCf[5][5];
    m_dK5Combo[7] = -m_dImpCf[5][2] / m_dImpCf[5][5];
    m_dK5Combo[8] = -m_dExpCf[5][2] / m_dImpCf[5][5];
	m_dK5Combo[9] = -m_dImpCf[5][3] / m_dImpCf[5][5];
    m_dK5Combo[10] = -m_dExpCf[5][3] / m_dImpCf[5][5];
	m_dK5Combo[11] = -m_dImpCf[5][4] / m_dImpCf[5][5];
    m_dK5Combo[12] = -m_dExpCf[5][4] / m_dImpCf[5][5];
	m_dK5Combo[13] = 0.0;

    // Kh5 combination
    m_dKh5Combo[0] = -1.0 / (dDeltaT * m_dExpCf[5][4]);
    m_dKh5Combo[1] = 1.0 / (dDeltaT * m_dExpCf[5][4]);
    m_dKh5Combo[2] = 0.0;
    m_dKh5Combo[3] = -m_dImpCf[5][0] / m_dExpCf[5][4];
    m_dKh5Combo[4] = -m_dExpCf[5][0] / m_dExpCf[5][4];
    m_dKh5Combo[5] = -m_dImpCf[5][1] / m_dExpCf[5][4];
    m_dKh5Combo[6] = -m_dExpCf[5][1] / m_dExpCf[5][4];
	m_dKh5Combo[7] = -m_dImpCf[5][2] / m_dExpCf[5][4];
	m_dKh5Combo[8] = -m_dExpCf[5][2] / m_dExpCf[5][4];
	m_dKh5Combo[9] = -m_dImpCf[5][3] / m_dExpCf[5][4];
	m_dKh5Combo[10] = -m_dExpCf[5][3] / m_dExpCf[5][4];
	m_dKh5Combo[11] = -m_dImpCf[5][4] / m_dExpCf[5][4];
	m_dKh5Combo[12] = 0.0;

	// uf6 explicit evaluation combination
    m_du6fCombo[0] = 1.0;
    m_du6fCombo[1] = 0.0;
    m_du6fCombo[2] = 0.0;
    m_du6fCombo[3] = dDeltaT * m_dImpCf[6][0];
    m_du6fCombo[4] = dDeltaT * m_dExpCf[6][0];
    m_du6fCombo[5] = dDeltaT * m_dImpCf[6][1];
    m_du6fCombo[6] = dDeltaT * m_dExpCf[6][1];
    m_du6fCombo[7] = dDeltaT * m_dImpCf[6][2];
    m_du6fCombo[8] = dDeltaT * m_dExpCf[6][2];
	m_du6fCombo[9] = dDeltaT * m_dImpCf[6][3];
	m_du6fCombo[10] = dDeltaT * m_dExpCf[6][3];
	m_du6fCombo[11] = dDeltaT * m_dImpCf[6][4];
	m_du6fCombo[12] = dDeltaT * m_dExpCf[6][4];
	m_du6fCombo[13] = dDeltaT * m_dImpCf[6][5];

	// PRE-STAGE 1 - Implicit solve to start
	Time timeSub0 = time;
	double dtSub0 = m_dImpCf[0][0] * dDeltaT;
	pGrid->CopyData(0, 1, DataType_State);
	pHorizontalDynamics->StepImplicit(0, 1, timeSub0, dtSub0);
	pVerticalDynamics->StepImplicit(0, 1, timeSub0, dtSub0);
	pGrid->PostProcessSubstage(1, DataType_State);
	pGrid->PostProcessSubstage(1, DataType_Tracers);

	// Store the evaluation K1 to index 3
	pGrid->LinearCombineData(m_dK0Combo, 3, DataType_State);

	// SUBSTAGE 1
	// Compute the explicit step to index 1
	Time timeSub1 = time;
	double dtSub1 = m_dExpCf[1][0] * dDeltaT;
	pGrid->LinearCombineData(m_du1fCombo, 1, DataType_State);
	pHorizontalDynamics->StepExplicit(0, 1, timeSub1, dtSub1);
	pVerticalDynamics->StepExplicit(0, 1, 0, timeSub1, dtSub1, false);
	pGrid->PostProcessSubstage(1, DataType_State);
	pGrid->PostProcessSubstage(1, DataType_Tracers);

	// Store the evaluation Kh1 to index 4
	pGrid->LinearCombineData(m_dKh1Combo, 4, DataType_State);

	// Compute implicit step based on known data to index 2
	pGrid->CopyData(1, 2, DataType_State);
	dtSub1 = m_dImpCf[1][1] * dDeltaT;
	pVerticalDynamics->StepImplicit(1, 2, timeSub1, dtSub1);
	pGrid->PostProcessSubstage(2, DataType_State);
	pGrid->PostProcessSubstage(2, DataType_Tracers);

	// Store the evaluation K1 to index 3
	pGrid->LinearCombineData(m_dK1Combo, 5, DataType_State);

	//std::cout << "Substage 1 done ... \n";

	// SUBSTAGE 2
	// Compute uf2
	Time timeSub2 = time;
	double dtSub2 = m_dExpCf[2][1] * dDeltaT;
	pGrid->LinearCombineData(m_du2fCombo, 1, DataType_State);
	pHorizontalDynamics->StepExplicit(2, 1, timeSub2, dtSub2);
	pVerticalDynamics->StepExplicit(2, 1, 0, timeSub2, dtSub2, false);
	pGrid->PostProcessSubstage(1, DataType_State);
	pGrid->PostProcessSubstage(1, DataType_Tracers);

	// Store the evaluation Kh2 to index 6
	pGrid->LinearCombineData(m_dKh2Combo, 6, DataType_State);

	// Compute u2 from uf2 and store it to index 2 (over u1)
	dtSub2 = m_dImpCf[2][2] * dDeltaT;
	pGrid->CopyData(1, 2, DataType_State);
	pVerticalDynamics->StepImplicit(1, 2, timeSub2, dtSub2);
	pGrid->PostProcessSubstage(2, DataType_State);
	pGrid->PostProcessSubstage(2, DataType_Tracers);

	// Store the evaluation K2 to index 7
	pGrid->LinearCombineData(m_dK2Combo, 7, DataType_State);

	//std::cout << "Substage 2 done ... \n";

	// SUBSTAGE 3
	// Compute uf3
	Time timeSub3 = time;
	double dtSub3 = m_dExpCf[3][2] * dDeltaT;
	pGrid->LinearCombineData(m_du3fCombo, 1, DataType_State);
	pHorizontalDynamics->StepExplicit(2, 1, timeSub2, dtSub2);
	pVerticalDynamics->StepExplicit(2, 1, 0, timeSub2, dtSub2, false);
	pGrid->PostProcessSubstage(1, DataType_State);
	pGrid->PostProcessSubstage(1, DataType_Tracers);

    // Store the evaluation Kh3 to index 8
    pGrid->LinearCombineData(m_dKh3Combo, 8, DataType_State);

    // Compute u3 from uf3 and store it to index 2 (over u2)
    dtSub3 = m_dImpCf[3][3] * dDeltaT;
    pGrid->CopyData(1, 2, DataType_State);
    pVerticalDynamics->StepImplicit(1, 2, timeSub3, dtSub3);
    pGrid->PostProcessSubstage(2, DataType_State);
    pGrid->PostProcessSubstage(2, DataType_Tracers);

    // Store the evaluation K3 to index 9
    pGrid->LinearCombineData(m_dK3Combo, 9, DataType_State);

	//std::cout << "Substage 3 done ... \n";

    // SUBSTAGE 4 - MODIFIED DUE TO EXPLICIT ZERO EVALUATION
    // Compute uf4 from u3 and store it to index 2
    Time timeSub4 = time;
    double dtSub4 = m_dExpCf[4][3] * dDeltaT;
    pGrid->LinearCombineData(m_du4fCombo, 1, DataType_State);
    pHorizontalDynamics->StepExplicit(2, 1, timeSub4, dtSub4);
    pVerticalDynamics->StepExplicit(2, 1, 0, timeSub4, dtSub4, false);
    pGrid->PostProcessSubstage(1, DataType_State);
    pGrid->PostProcessSubstage(1, DataType_Tracers);

	// Store the evaluation Kh4 to index 10
    pGrid->LinearCombineData(m_dKh4Combo, 10, DataType_State);

    // Compute u4 from uf4 and store it to index 2 (over u3)
    dtSub4 = m_dImpCf[4][4] * dDeltaT;
    //pGrid->CopyData(1, 2, DataType_State);
	// Start with the previous data sum here because of the explicit zero
	pGrid->LinearCombineData(m_du4fCombo, 2, DataType_State);
    pVerticalDynamics->StepImplicit(2, 2, timeSub4, dtSub4);
    pGrid->PostProcessSubstage(2, DataType_State);
    pGrid->PostProcessSubstage(2, DataType_Tracers);

	// Store the evaluation K4 to index 11
    pGrid->LinearCombineData(m_dK4Combo, 11, DataType_State);

	//std::cout << "Substage 4 done ... \n";

    // SUBSTAGE 5
    // Compute uf5 from u4 and store it to index 2
    Time timeSub5 = time;
    double dtSub5 = m_dExpCf[5][4] * dDeltaT;
    pGrid->LinearCombineData(m_du5fCombo, 1, DataType_State);
    pHorizontalDynamics->StepExplicit(2, 1, timeSub5, dtSub5);
    pVerticalDynamics->StepExplicit(2, 1, 0, timeSub5, dtSub5, false);
    pGrid->PostProcessSubstage(1, DataType_State);
    pGrid->PostProcessSubstage(1, DataType_Tracers);

	// Store the evaluation Kh5 to index 12
    pGrid->LinearCombineData(m_dKh5Combo, 12, DataType_State);

    // Compute u5 from uf5 and store it to index 2 (over u4)
    dtSub5 = m_dImpCf[5][5] * dDeltaT;
    pGrid->CopyData(1, 2, DataType_State);
    pVerticalDynamics->StepImplicit(1, 2, timeSub5, dtSub5);
    pGrid->PostProcessSubstage(2, DataType_State);
    pGrid->PostProcessSubstage(2, DataType_Tracers);

	// Store the evaluation K5 to index 13
    pGrid->LinearCombineData(m_dK5Combo, 13, DataType_State);

	//std::cout << "Substage 5 done ... \n";

    // SUBSTAGE 6
    // Compute uf6 from u5 and store it to index 2
    Time timeSub6 = time;
    double dtSub6 = m_dExpCf[6][5] * dDeltaT;
    pGrid->LinearCombineData(m_du6fCombo, 1, DataType_State);
    pHorizontalDynamics->StepExplicit(2, 1, timeSub6, dtSub6);
    pVerticalDynamics->StepExplicit(2, 1, 0, timeSub6, dtSub6, false);
    pGrid->PostProcessSubstage(1, DataType_State);
    pGrid->PostProcessSubstage(1, DataType_Tracers);

    // Compute u5 from uf5 and store it to index 2 (over u4)
    dtSub5 = m_dImpCf[6][6] * dDeltaT;
    pGrid->CopyData(1, 2, DataType_State);
    pVerticalDynamics->StepImplicit(1, 2, timeSub6, dtSub6);
    pGrid->PostProcessSubstage(2, DataType_State);
    pGrid->PostProcessSubstage(2, DataType_Tracers);

    // Apply hyperdiffusion at the end of the explicit substep (ask Paul)
	pGrid->CopyData(2, 1, DataType_State);
	pHorizontalDynamics->StepAfterSubCycle(2, 1, 3, time, dDeltaT);
	pGrid->CopyData(1, 0, DataType_State);
	//*/
}

///////////////////////////////////////////////////////////////////////////////

