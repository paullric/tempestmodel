///////////////////////////////////////////////////////////////////////////////
///
///	\file    TimestepSchemeStrang.cpp
///	\author  Paul Ullrich
///	\version June 18, 2013
///
///	<remarks>
///		Copyright 2000-2010 Paul Ullrich
///
///		This file is distributed as part of the Tempest source code package.
///		Permission is granted to use, copy, modify and distribute this
///		source code and its documentation under the terms of the GNU General
///		Public License.  This software is provided "as is" without express
///		or implied warranty.
///	</remarks>

#include "TimestepSchemeStrang.h"
#include "Model.h"
#include "Grid.h"
#include "HorizontalDynamics.h"
#include "VerticalDynamics.h"
#include "Announce.h"
#include "Defines.h"

///////////////////////////////////////////////////////////////////////////////

TimestepSchemeStrang::TimestepSchemeStrang(
	Model & model,
	double dOffCentering,
	ExplicitDiscretization eExplicitDiscretization
) :
	TimestepScheme(model),
	m_dOffCentering(dOffCentering),
	m_eExplicitDiscretization(eExplicitDiscretization)
{
	// Check bounds of OffCentering parameter
	if ((m_dOffCentering < 0.0) || (m_dOffCentering > 1.0)) {
		_EXCEPTIONT("OffCentering parameter out of range [0,1]");
	}

	if (m_eExplicitDiscretization == ForwardEuler) {
		m_nExplicitSubSteps = 1;
	} else if (m_eExplicitDiscretization == RungeKutta4) {
		m_nExplicitSubSteps = 4;
	} else if (m_eExplicitDiscretization == RungeKuttaSSP3) {
		m_nExplicitSubSteps = 3;
	} else if (m_eExplicitDiscretization == KinnmarkGrayUllrich35) {
		m_nExplicitSubSteps = 5;
	} else if (m_eExplicitDiscretization == RungeKuttaSSPRK53) {
		m_nExplicitSubSteps = 5;
	} else {
		_EXCEPTIONT("Invalid ExplicitDiscretization");
	}

	// Carryover combination
	m_dCarryoverCombination.Allocate(2);
	m_dCarryoverCombination[0] = 1.0;
	m_dCarryoverCombination[1] = 1.0;

	// Off-centering combination
	m_dOffCenteringCombination.Allocate(2);
	m_dOffCenteringCombination[0] = (2.0 - m_dOffCentering) / 2.0;
	m_dOffCenteringCombination[1] = m_dOffCentering / 2.0;

	// Final carryover combination
	m_dCarryoverFinal.Allocate(2);
	m_dCarryoverFinal[0] = +1.0;
	m_dCarryoverFinal[1] = -1.0;

	// RK4 combination
	m_dRK4Combination.Allocate(5);
	m_dRK4Combination[0] = - 1.0 / 3.0;
	m_dRK4Combination[1] = + 1.0 / 3.0;
	m_dRK4Combination[2] = + 2.0 / 3.0;
	m_dRK4Combination[3] = + 1.0 / 3.0;

	// SSPRK3 combination A
	m_dSSPRK3CombinationA.Allocate(3);
	m_dSSPRK3CombinationA[0] = 3.0 / 4.0;
	m_dSSPRK3CombinationA[1] = 1.0 / 4.0;
	m_dSSPRK3CombinationA[2] = 0.0;

	// SSPRK3 combination B
	m_dSSPRK3CombinationB.Allocate(5);
	m_dSSPRK3CombinationB[0] = 1.0 / 3.0;
	m_dSSPRK3CombinationB[1] = 0.0;
	m_dSSPRK3CombinationB[2] = 2.0 / 3.0;
	m_dSSPRK3CombinationB[3] = 0.0;
	m_dSSPRK3CombinationB[4] = 0.0;

	// KGU 3-5 combination
	m_dKinnmarkGrayUllrichCombination.Allocate(5);
	m_dKinnmarkGrayUllrichCombination[0] = - 1.0 / 4.0;
	m_dKinnmarkGrayUllrichCombination[1] =   5.0 / 4.0;
	m_dKinnmarkGrayUllrichCombination[2] =   0.0;
	m_dKinnmarkGrayUllrichCombination[3] =   0.0;
	m_dKinnmarkGrayUllrichCombination[4] =   0.0;

	// SSPRK53 combination A
	m_dSSPRK53CombinationA.Allocate(4);
	m_dSSPRK53CombinationA[0] = 0.355909775063327;
	m_dSSPRK53CombinationA[1] = 0.0;
	m_dSSPRK53CombinationA[2] = 0.644090224936674;
	m_dSSPRK53CombinationA[3] = 0.0;

	// SSPRK53 combination B
	m_dSSPRK53CombinationB.Allocate(4);
	m_dSSPRK53CombinationB[0] = 0.367933791638137;
	m_dSSPRK53CombinationB[1] = 0.0;
	m_dSSPRK53CombinationB[2] = 0.0;
	m_dSSPRK53CombinationB[3] = 0.632066208361863;

	// SSPRK53 combination C
	m_dSSPRK53CombinationC.Allocate(5);
	m_dSSPRK53CombinationC[0] = 0.762406163401431;
	m_dSSPRK53CombinationC[1] = 0.0;
	m_dSSPRK53CombinationC[2] = 0.237593836598569;
	m_dSSPRK53CombinationC[3] = 0.0;
	m_dSSPRK53CombinationC[4] = 0.0;

	m_dResCombine.Allocate(5);
}

///////////////////////////////////////////////////////////////////////////////
/*
int TimestepSchemeStrang::GetSubStepCount() const {

	HorizontalDynamics * pHorizontalDynamics =
		m_model.GetHorizontalDynamics();

	if (pHorizontalDynamics == NULL) {
		_EXCEPTIONT("HorizontalDynamics has not been initialized");
	}

	int nHorizontalDynamicsSubSteps =
		pHorizontalDynamics->GetSubStepAfterSubCycleCount();

	return (m_nExplicitSubSteps + nHorizontalDynamicsSubSteps + 1);
}
*/
///////////////////////////////////////////////////////////////////////////////

double TimestepSchemeStrang::GetMaximumStableCourantNumber(
	TimestepScheme::MixedMethodPart eMixedMethodPart,
	int nOrder
) const {
	if (m_eExplicitDiscretization == KinnmarkGrayUllrich35) {
		if (eMixedMethodPart == TimestepScheme::ContinuousPart) {
			if (nOrder == 2) {
				return 3.873077;
			} else if (nOrder == 3) {
				return 2.582184;
			} else if (nOrder == 4) {
				return 2.121307;
			} else if (nOrder == 5) {
				return 1.851593;
			} else if (nOrder == 6) {
				return 1.653839;
			} else if (nOrder == 7) {
				return 1.491241;
			} else if (nOrder == 8) {
				return 1.353363;
			} else if (nOrder == 9) {
				return 1.234161;
			} else if (nOrder == 10) {
				return 1.130890;
			} else {
				return 0.0;
			}

		} else {
			if (nOrder == 1) {
				return 1.524200;
			} else if (nOrder == 2) {
				return 2.824432;
			} else if (nOrder == 3) {
				return 1.686798;
			} else if (nOrder == 4) {
				return 1.263824;
			} else if (nOrder == 5) {
				return 1.034760;
			} else if (nOrder == 6) {
				return 0.888092;
			} else if (nOrder == 7) {
				return 0.784271;
			} else if (nOrder == 8) {
				return 0.705719;
			} else if (nOrder == 9) {
				return 0.644745;
			} else if (nOrder == 10) {
				return 0.594757;
			} else {
				return 0.0;
			}
		}

	} else {
		return 0.0;
	}
}

///////////////////////////////////////////////////////////////////////////////

void TimestepSchemeStrang::Step(
	bool fFirstStep,
	bool fLastStep,
	const Time & time,
	double dDeltaT
) {
	//Announce("Top of the time Step... %s");
	// Get a copy of the grid
	Grid * pGrid = m_model.GetGrid();

	// Get a copy of the HorizontalDynamics
	HorizontalDynamics * pHorizontalDynamics =
		m_model.GetHorizontalDynamics();

	// Get a copy of the VerticalDynamics
	VerticalDynamics * pVerticalDynamics =
		m_model.GetVerticalDynamics();

	// Half a time step
	double dHalfDeltaT = 0.5 * dDeltaT;

	// Vertical timestep
	if (fFirstStep) {
#if defined(RESIDUAL_DIFFUSION)
		pGrid->ZeroData(0, DataType_Residual);
		pGrid->ZeroData(1, DataType_Residual);
		pGrid->ZeroData(2, DataType_Residual);

		pGrid->CopyData(0, 4, DataType_State);
		pGrid->CopyData(0, 4, DataType_Tracers);
		pVerticalDynamics->StepImplicit(0, 4, time, dHalfDeltaT);
		pHorizontalDynamics->StepImplicit(0, 4, time, dHalfDeltaT);

		//std::cout << "Before computing G function in TimestepSchemeStrang";
		// Send the function G (implicit stuff) to index 0 in residual storage
		//DataArray1D<double> m_dResCombine;
		//m_dResCombine.Allocate(5);
		m_dResCombine[0] = -1.0 / dDeltaT;
		m_dResCombine[1] = 0.0;
		m_dResCombine[2] = 0.0;
		m_dResCombine[3] = 0.0;
		m_dResCombine[4] = 1.0 / dDeltaT;
		pGrid->LinearCombineData2Residual(m_dResCombine, 0, DataType_Residual);

		pGrid->CopyData(4, 0, DataType_State);
		pGrid->CopyData(4, 0, DataType_Tracers);
#else
		pVerticalDynamics->StepImplicit(0, 0, time, dHalfDeltaT);
		pHorizontalDynamics->StepImplicit(0, 0, time, dHalfDeltaT);
#endif

	} else {
		pGrid->LinearCombineData(m_dCarryoverCombination, 0, DataType_State);
		pGrid->LinearCombineData(m_dCarryoverCombination, 0, DataType_Tracers);

		pVerticalDynamics->FilterNegativeTracers(0);
	}

	// Forward Euler
	if (m_eExplicitDiscretization == ForwardEuler) {
		pGrid->CopyData(0, 4, DataType_State);
		pGrid->CopyData(0, 4, DataType_Tracers);
		pHorizontalDynamics->StepExplicit(0, 4, time, dDeltaT);
		pVerticalDynamics->StepExplicit(0, 4, time, dDeltaT);
		pGrid->PostProcessSubstage(4, DataType_State);
		pGrid->PostProcessSubstage(4, DataType_Tracers);

#if defined(RESIDUAL_DIFFUSION)
		// Send the function F (explicit stuff) to index 1 in residual storage
		//DataArray1D<double> m_dResCombine;
		//m_dResCombine.Initialize(5);
		m_dResCombine[0] = -1.0 / dDeltaT;
		m_dResCombine[1] = 0.0;
		m_dResCombine[2] = 0.0;
		m_dResCombine[3] = 0.0;
		m_dResCombine[4] = 1.0 / dDeltaT;
		pGrid->LinearCombineData2Residual(m_dResCombine, 1, DataType_Residual);
#endif

	// Explicit fourth-order Runge-Kutta
	} else if (m_eExplicitDiscretization == RungeKutta4) {
		pGrid->CopyData(0, 1, DataType_State);
		pGrid->CopyData(0, 1, DataType_Tracers);
		pHorizontalDynamics->StepExplicit(0, 1, time, dHalfDeltaT);
		pVerticalDynamics->StepExplicit(0, 1, time, dHalfDeltaT);
		pGrid->PostProcessSubstage(1, DataType_State);
		pGrid->PostProcessSubstage(1, DataType_Tracers);

		pGrid->CopyData(0, 2, DataType_State);
		pGrid->CopyData(0, 2, DataType_Tracers);
		pHorizontalDynamics->StepExplicit(1, 2, time, dHalfDeltaT);
		pVerticalDynamics->StepExplicit(1, 2, time, dHalfDeltaT);
		pGrid->PostProcessSubstage(2, DataType_State);
		pGrid->PostProcessSubstage(2, DataType_Tracers);

		pGrid->CopyData(0, 3, DataType_State);
		pGrid->CopyData(0, 3, DataType_Tracers);
		pHorizontalDynamics->StepExplicit(2, 3, time, dDeltaT);
		pVerticalDynamics->StepExplicit(2, 3, time, dDeltaT);
		pGrid->PostProcessSubstage(3, DataType_State);
		pGrid->PostProcessSubstage(3, DataType_Tracers);

		pGrid->LinearCombineData(m_dRK4Combination, 4, DataType_State);
		pGrid->LinearCombineData(m_dRK4Combination, 4, DataType_Tracers);

		pHorizontalDynamics->StepExplicit(3, 4, time, dDeltaT / 6.0);
		pVerticalDynamics->StepExplicit(3, 4, time, dDeltaT / 6.0);
		pGrid->PostProcessSubstage(4, DataType_State);
		pGrid->PostProcessSubstage(4, DataType_Tracers);

#if defined(RESIDUAL_DIFFUSION)
		// Send the function F (explicit stuff) to index 1 in residual storage
		//DataArray1D<double> m_dResCombine;
		//m_dResCombine.Initialize(5);
		m_dResCombine[0] = +2.0 / dDeltaT;
		m_dResCombine[1] = -2.0 / dDeltaT;
		m_dResCombine[2] = -4.0 / dDeltaT;
		m_dResCombine[3] = -2.0 / dDeltaT;
		m_dResCombine[4] = +6.0 / dDeltaT;
		pGrid->LinearCombineData2Residual(m_dResCombine, 1, DataType_Residual);
#endif

	// Explicit strong stability preserving third-order Runge-Kutta
	} else if (m_eExplicitDiscretization == RungeKuttaSSP3) {

		pGrid->CopyData(0, 1, DataType_State);
		pGrid->CopyData(0, 1, DataType_Tracers);
		pHorizontalDynamics->StepExplicit(0, 1, time, dDeltaT);
		pVerticalDynamics->StepExplicit(0, 1, time, dDeltaT);
		pGrid->PostProcessSubstage(1, DataType_State);
		pGrid->PostProcessSubstage(1, DataType_Tracers);

		pGrid->LinearCombineData(m_dSSPRK3CombinationA, 2, DataType_State);
		pGrid->LinearCombineData(m_dSSPRK3CombinationA, 2, DataType_Tracers);
		pHorizontalDynamics->StepExplicit(1, 2, time, 0.25 * dDeltaT);
		pVerticalDynamics->StepExplicit(1, 2, time, 0.25 * dDeltaT);
		pGrid->PostProcessSubstage(2, DataType_State);
		pGrid->PostProcessSubstage(2, DataType_Tracers);

		pGrid->LinearCombineData(m_dSSPRK3CombinationB, 4, DataType_State);
		pGrid->LinearCombineData(m_dSSPRK3CombinationB, 4, DataType_Tracers);
		pHorizontalDynamics->StepExplicit(2, 4, time, (2.0/3.0) * dDeltaT);
		pVerticalDynamics->StepExplicit(2, 4, time, (2.0/3.0) * dDeltaT);
		pGrid->PostProcessSubstage(4, DataType_State);
		pGrid->PostProcessSubstage(4, DataType_Tracers);

#if defined(RESIDUAL_DIFFUSION)
		// Send the function F (explicit stuff) to index 1 in residual storage
		//DataArray1D<double> m_dResCombine;
		//m_dResCombine.Initialize(5);
		m_dResCombine[0] = -1.0 / (2.0 * dDeltaT);
		m_dResCombine[1] = 0.0;
		m_dResCombine[2] = -1.0 / dDeltaT;
		m_dResCombine[3] = 0.0;
		m_dResCombine[4] = +3.0 / (2.0 * dDeltaT);
		pGrid->LinearCombineData2Residual(m_dResCombine, 1, DataType_Residual);
#endif

	// Explicit Kinnmark, Gray and Ullrich third-order five-stage Runge-Kutta
	} else if (m_eExplicitDiscretization == KinnmarkGrayUllrich35) {

		pGrid->CopyData(0, 1, DataType_State);
		pGrid->CopyData(0, 1, DataType_Tracers);
		pHorizontalDynamics->StepExplicit(0, 1, time, dDeltaT / 5.0);
		pVerticalDynamics->StepExplicit(0, 1, time, dDeltaT / 5.0);
		pGrid->PostProcessSubstage(1, DataType_State);
		pGrid->PostProcessSubstage(1, DataType_Tracers);

		pGrid->CopyData(0, 2, DataType_State);
		pGrid->CopyData(0, 2, DataType_Tracers);
		pHorizontalDynamics->StepExplicit(1, 2, time, dDeltaT / 5.0);
		pVerticalDynamics->StepExplicit(1, 2, time, dDeltaT / 5.0);
		pGrid->PostProcessSubstage(2, DataType_State);
		pGrid->PostProcessSubstage(2, DataType_Tracers);

		pGrid->CopyData(0, 3, DataType_State);
		pGrid->CopyData(0, 3, DataType_Tracers);
		pHorizontalDynamics->StepExplicit(2, 3, time, dDeltaT / 3.0);
		pVerticalDynamics->StepExplicit(2, 3, time, dDeltaT / 3.0);
		pGrid->PostProcessSubstage(3, DataType_State);
		pGrid->PostProcessSubstage(3, DataType_Tracers);

		pGrid->CopyData(0, 2, DataType_State);
		pGrid->CopyData(0, 2, DataType_Tracers);
		pHorizontalDynamics->StepExplicit(3, 2, time, 2.0 * dDeltaT / 3.0);
		pVerticalDynamics->StepExplicit(3, 2, time, 2.0 * dDeltaT / 3.0);
		pGrid->PostProcessSubstage(2, DataType_State);
		pGrid->PostProcessSubstage(2, DataType_Tracers);

		pGrid->LinearCombineData(
			m_dKinnmarkGrayUllrichCombination, 4, DataType_State);
		pGrid->LinearCombineData(
			m_dKinnmarkGrayUllrichCombination, 4, DataType_Tracers);
		pHorizontalDynamics->StepExplicit(2, 4, time, 3.0 * dDeltaT / 4.0);
		pVerticalDynamics->StepExplicit(2, 4, time, 3.0 * dDeltaT / 4.0);
		pGrid->PostProcessSubstage(4, DataType_State);
		pGrid->PostProcessSubstage(4, DataType_Tracers);

#if defined(RESIDUAL_DIFFUSION)
		// Send the function F (explicit stuff) to index 1 in residual storage
		//DataArray1D<double> m_dResCombine;
		//m_dResCombine.Initialize(5);
		m_dResCombine[0] = +1.0 / (3.0 * dDeltaT);
		m_dResCombine[1] = -5.0 / (3.0 * dDeltaT);
		m_dResCombine[2] = 0.0;
		m_dResCombine[3] = 0.0;
		m_dResCombine[4] = +4.0 / (3.0 * dDeltaT);
		pGrid->LinearCombineData2Residual(m_dResCombine, 1, DataType_Residual);
#endif

	// Explicit strong stability preserving five-stage third-order Runge-Kutta
	} else if (m_eExplicitDiscretization == RungeKuttaSSPRK53) {

		const double dStepOne = 0.377268915331368;

		pGrid->CopyData(0, 1, DataType_State);
		pGrid->CopyData(0, 1, DataType_Tracers);
		pHorizontalDynamics->StepExplicit(0, 1, time, dStepOne * dDeltaT);
		pVerticalDynamics->StepExplicit(0, 1, time, dStepOne * dDeltaT);
		pGrid->PostProcessSubstage(1, DataType_State);
		pGrid->PostProcessSubstage(1, DataType_Tracers);

		pGrid->CopyData(1, 2, DataType_State);
		pGrid->CopyData(1, 2, DataType_Tracers);
		pHorizontalDynamics->StepExplicit(1, 2, time, dStepOne * dDeltaT);
		pVerticalDynamics->StepExplicit(1, 2, time, dStepOne * dDeltaT);
		pGrid->PostProcessSubstage(2, DataType_State);
		pGrid->PostProcessSubstage(2, DataType_Tracers);

		const double dStepThree = 0.242995220537396;

		pGrid->LinearCombineData(m_dSSPRK53CombinationA, 3, DataType_State);
		pGrid->LinearCombineData(m_dSSPRK53CombinationA, 3, DataType_Tracers);
		pHorizontalDynamics->StepExplicit(2, 3, time, dStepThree * dDeltaT);
		pVerticalDynamics->StepExplicit(2, 3, time, dStepThree * dDeltaT);
		pGrid->PostProcessSubstage(3, DataType_State);
		pGrid->PostProcessSubstage(3, DataType_Tracers);

		const double dStepFour = 0.238458932846290;

		pGrid->LinearCombineData(m_dSSPRK53CombinationB, 0, DataType_State);
		pGrid->LinearCombineData(m_dSSPRK53CombinationB, 0, DataType_Tracers);
		pHorizontalDynamics->StepExplicit(3, 0, time, dStepFour * dDeltaT);
		pVerticalDynamics->StepExplicit(3, 0, time, dStepFour * dDeltaT);
		pGrid->PostProcessSubstage(0, DataType_State);
		pGrid->PostProcessSubstage(0, DataType_Tracers);

		const double dStepFive = 0.287632146308408;

		pGrid->LinearCombineData(m_dSSPRK53CombinationC, 4, DataType_State);
		pGrid->LinearCombineData(m_dSSPRK53CombinationC, 4, DataType_Tracers);
		pHorizontalDynamics->StepExplicit(0, 4, time, dStepFive * dDeltaT);
		pVerticalDynamics->StepExplicit(0, 4, time, dStepFive * dDeltaT);
		pGrid->PostProcessSubstage(4, DataType_State);
		pGrid->PostProcessSubstage(4, DataType_Tracers);

#if defined(RESIDUAL_DIFFUSION)
		// Send the function F (explicit stuff) to index 1 in residual storage
		//DataArray1D<double> m_dResCombine;
		//m_dResCombine.Initialize(5);
		m_dResCombine[0] = -m_dSSPRK53CombinationC[0] / (dStepFive * dDeltaT);
		m_dResCombine[1] = 0.0;
		m_dResCombine[2] = -m_dSSPRK53CombinationC[2] / (dStepFive * dDeltaT);
		m_dResCombine[3] = 0.0;
		m_dResCombine[4] = +1.0 / (dStepFive * dDeltaT);
		pGrid->LinearCombineData2Residual(m_dResCombine, 1, DataType_Residual);
#endif

	// Invalid explicit discretization
	} else {
		_EXCEPTIONT("Invalid explicit discretization");
	}

#if defined(RESIDUAL_DIFFUSION)
		// Estimate the total time derivative and send to index 2 residual storage
		//DataArray1D<double> m_dResCombine;
		//m_dResCombine.Initialize(5);
		m_dResCombine[0] = -1.0 / dDeltaT;
		m_dResCombine[1] = 0.0;
		m_dResCombine[2] = 0.0;
		m_dResCombine[3] = 0.0;
		m_dResCombine[4] = +1.0 / dDeltaT;
		pGrid->LinearCombineData2Residual(m_dResCombine, 2, DataType_Residual);

		// Estimate the residual and send to index 2 of residual storage
		//DataArray1D<double> m_dResCombine;
		//m_dResCombine.Initialize(3);
		m_dResCombine[0] = +1.0;
		m_dResCombine[1] = +1.0;
		m_dResCombine[2] = +1.0;
		m_dResCombine[3] = 0.0;
		m_dResCombine[4] = 0.0;
		pGrid->LinearCombineData(m_dResCombine, 2, DataType_Residual);
#endif

	// Apply hyperdiffusion and Rayleigh
	pGrid->CopyData(4, 1, DataType_State);
	pGrid->CopyData(4, 1, DataType_Tracers);
	pHorizontalDynamics->StepAfterSubCycle(4, 1, 2, time, dDeltaT);

	// Vertical timestep
	double dOffCenterDeltaT = 0.5 * (1.0 + m_dOffCentering) * dDeltaT;

	pGrid->CopyData(1, 0, DataType_State);
	pGrid->CopyData(1, 0, DataType_Tracers);
	pVerticalDynamics->StepImplicit(0, 0, time, dOffCenterDeltaT);
	pHorizontalDynamics->StepImplicit(0, 0, time, dOffCenterDeltaT);

	pGrid->LinearCombineData(m_dOffCenteringCombination, 0, DataType_State);
	pGrid->LinearCombineData(m_dOffCenteringCombination, 0, DataType_Tracers);

	if (!fLastStep) {
		pGrid->LinearCombineData(m_dCarryoverFinal, 1, DataType_State);
		pGrid->LinearCombineData(m_dCarryoverFinal, 1, DataType_Tracers);
	}

	//pGrid->CopyData(0, 1, DataType_State);
	//pVerticalDynamics->StepExplicit(0, 1, time, dDeltaT);
	//pVerticalDynamics->StepImplicit(0, 0, time, dDeltaT);
/*
	if (0) {
	//if (fLastStep) {
		DataArray1D<double> dDifference;
		dDifference.Initialize(2);
		dDifference[0] = -1.0;
		dDifference[1] = 1.0;
		pGrid->LinearCombineData(dDifference, 0, DataType_State);
	} else {
		pGrid->CopyData(4, 0, DataType_State);
	}
*/
}

///////////////////////////////////////////////////////////////////////////////
