///////////////////////////////////////////////////////////////////////////////
///
///	\file    TimestepSchemeARK4.cpp
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

#include "TimestepSchemeARK4.h"
#include "Model.h"
#include "Grid.h"
#include "HorizontalDynamics.h"
#include "VerticalDynamics.h"

///////////////////////////////////////////////////////////////////////////////
/*
// Time coefficients for each stage
const double TimestepSchemeARK4::m_timeCoeff[]
	= { 0.0, 0.5, 0.332, 0.62, 0.85, 1.0 };

// Stage coefficients - each row is for that stage
const double TimestepSchemeARK4::m_dExplicitCoeff[][] = {
  {0., 0., 0., 0., 0., 0.},
  {0.5, 0., 0., 0., 0., 0.},
  {0.221776, 0.110224, 0., 0., 0., 0.},
  {-0.04884659515311857, -0.17772065232640102, 0.8465672474795197, 0., 0., 0.},
  {-0.15541685842491548, -0.3567050098221991, 1.0587258798684427, 0.30339598837867193, 0., 0.},
  { 0.2014243506726763, 0.008742057842904185, 0.15993995707168115, 0.4038290605220775, 0.22606457389066084, 0.}
};

// Implicit stage coefficients
const double TimestepSchemeARK4::m_dImplicitCoeff[][] = {
  {0., 0., 0., 0., 0., 0.},
  {0.25, 0.25, 0., 0., 0., 0.},
  {0.137776, -0.055776, 0.25, 0., 0., 0.},
  {0.14463686602698217, -0.22393190761334475, 0.4492950415863626, 0.25, 0., 0.},
  {0.09825878328356477, -0.5915442428196704, 0.8101210538282996, 0.283164405707806, 0.25, 0.},
  {0.15791629516167136, 0., 0.18675894052400077, 0.6805652953093346, -0.27524053099500667, 0.25}
};

// Final stage coefficients
const double TimestepSchemeARK4::m_dBCoeff[][] = {
  {0.15791629516167136, 0., 0.18675894052400077, 0.6805652953093346, -0.27524053099500667, 0.25};

// Coefficients for dense ouput, 4th-order interpolation
const double TimestepSchemeARK4::m_dInterpCoeff[][] = {
  {0.961753400252887, 0., 0.787405595186356, -2.74544192086633, 3.70351728061223, -1.70723435518514},
  {-1.76418754019038, 0., -0.774504669155511, 9.64023584441292, -12.544886411271, 5.44334277620397},
  {0.960350435099165, 0., 0.173858014493155, -6.21422862823726, 8.56612859966376, -3.48610842101883}
};
*/
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
	HorizontalDynamics * pHorizontalDynamics =
		m_model.GetHorizontalDynamics();

	// Get a copy of the VerticalDynamics
	VerticalDynamics * pVerticalDynamics =
		m_model.GetVerticalDynamics();

	// Half a time step
	double dHalfDeltaT = 0.5 * dDeltaT;

	// Vertical timestep
	if (fFirstStep) {
		pVerticalDynamics->StepImplicit(0, 0, time, dHalfDeltaT);
	} else {
		DataVector<double> dCarryoverCombination;
		dCarryoverCombination.Initialize(2);
		dCarryoverCombination[0] = 1.0;
		dCarryoverCombination[1] = 1.0;
		pGrid->LinearCombineData(dCarryoverCombination, 0, DataType_State);
		dCarryoverCombination.Deinitialize();
	}
/*
	// Explicit RK4
	pGrid->CopyData(0, 1, DataType_State);
	pHorizontalDynamics->StepExplicit(0, 1, time, dHalfDeltaT);
	pVerticalDynamics->StepExplicit(0, 1, time, dHalfDeltaT);

	pGrid->CopyData(0, 2, DataType_State);
	pHorizontalDynamics->StepExplicit(1, 2, time, dHalfDeltaT);
	pVerticalDynamics->StepExplicit(1, 2, time, dHalfDeltaT);

	pGrid->CopyData(0, 3, DataType_State);
	pHorizontalDynamics->StepExplicit(2, 3, time, dDeltaT);
	pVerticalDynamics->StepExplicit(2, 3, time, dDeltaT);

	DataVector<double> dRK4Combination;
	dRK4Combination.Initialize(5);
	dRK4Combination[0] = - 1.0 / 3.0;
	dRK4Combination[1] = + 1.0 / 3.0;
	dRK4Combination[2] = + 2.0 / 3.0;
	dRK4Combination[3] = + 1.0 / 3.0;
	pGrid->LinearCombineData(dRK4Combination, 4, DataType_State);

	pHorizontalDynamics->StepExplicit(
		3, 4, time + dDeltaT / 6.0, dDeltaT / 6.0);
	pVerticalDynamics->StepExplicit(
		3, 4, time + dDeltaT / 6.0, dDeltaT / 6.0);
*/

	// SSP RK3
	pGrid->CopyData(0, 1, DataType_State);
	pHorizontalDynamics->StepExplicit(0, 1, time, dDeltaT);
	pVerticalDynamics->StepExplicit(0, 1, time, dDeltaT);

	DataVector<double> dSSPRK3CombinationA;
	dSSPRK3CombinationA.Initialize(3);
	dSSPRK3CombinationA[0] = 3.0 / 4.0;
	dSSPRK3CombinationA[1] = 1.0 / 4.0;
	dSSPRK3CombinationA[2] = 0.0;
	pGrid->LinearCombineData(dSSPRK3CombinationA, 2, DataType_State);
	pHorizontalDynamics->StepExplicit(1, 2, time, 0.25 * dDeltaT);
	pVerticalDynamics->StepExplicit(1, 2, time, 0.25 * dDeltaT);
	dSSPRK3CombinationA.Deinitialize();

	DataVector<double> dSSPRK3CombinationB;
	dSSPRK3CombinationB.Initialize(5);
	dSSPRK3CombinationB[0] = 1.0 / 3.0;
	dSSPRK3CombinationB[1] = 0.0;
	dSSPRK3CombinationB[2] = 2.0 / 3.0;
	dSSPRK3CombinationB[3] = 0.0;
	dSSPRK3CombinationB[4] = 0.0;
	pGrid->LinearCombineData(dSSPRK3CombinationB, 4, DataType_State);
	pHorizontalDynamics->StepExplicit(2, 4, time, (2.0/3.0) * dDeltaT);
	pVerticalDynamics->StepExplicit(2, 4, time, (2.0/3.0) * dDeltaT);
    dSSPRK3CombinationB.Deinitialize();

	// Apply hyperdiffusion
	pGrid->CopyData(4, 1, DataType_State);
	pHorizontalDynamics->StepAfterSubCycle(4, 1, 2, time, dDeltaT);

	// Vertical timestep
	pGrid->CopyData(1, 0, DataType_State);
	pVerticalDynamics->StepImplicit(0, 0, time, dHalfDeltaT);

	if (!fLastStep) {
		DataVector<double> dCarryoverFinal;
		dCarryoverFinal.Initialize(2);
		dCarryoverFinal[0] = +1.0;
		dCarryoverFinal[1] = -1.0;
		pGrid->LinearCombineData(dCarryoverFinal, 1, DataType_State);
		dCarryoverFinal.Deinitialize();
	}

	//pGrid->CopyData(0, 1, DataType_State);
	//pVerticalDynamics->StepExplicit(0, 1, time, dDeltaT);
	//pVerticalDynamics->StepImplicit(0, 0, time, dDeltaT);
/*
	if (0) {
	//if (fLastStep) {
		DataVector<double> dDifference;
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

