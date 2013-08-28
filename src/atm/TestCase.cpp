///////////////////////////////////////////////////////////////////////////////
///
///	\file    TestCase.cpp
///	\author  Paul Ullrich
///	\version February 25, 2013
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

#include "TestCase.h"
/*
///////////////////////////////////////////////////////////////////////////////

void TestCase::PrepareOutput(
	const TestCaseInitialConditions & state,
	GridType eGridType,
	int nResolution,
	int nOrder,
	std::string strOutputFile
) {
}

///////////////////////////////////////////////////////////////////////////////

void TestCase::InitializeData(
	const TestCaseInitialConditions & state,
	PhysicalConstants & phys,
	GridData4D & data
) {
	// Evaluate physical constants
	state.EvaluatePhysicalConstants(phys);
}

///////////////////////////////////////////////////////////////////////////////


#include "TestCase.h"
#include "SystemState.h"
#include "PhysicalConstants.h"
#include "Geometry.h"

#include "MathHelper.h"
#include "GaussianQuadrature.h"

///////////////////////////////////////////////////////////////////////////////

void TestCase::EvaluatePointwiseVelocity(
	const SystemStateBlock & block,
	const PhysicalConstants & phys,
	double dTime,
	int iPanel,
	double dR,
	double dAlpha,
	double dBeta,
	double * dU
) {
	_EXCEPTIONT("No pointwise velocity function specified:\n"
		"Prescribed velocities unavailable.");
}

///////////////////////////////////////////////////////////////////////////////

void TestCase::EvaluateElementwiseState(
	const SystemStateBlock & block,
	const PhysicalConstants & phys,
	int iPanel,
	int iR,
	int iA,
	int iB,
	double * dState,
	double * dTracer
) {

	int p;
	int q;

	int c;

	double dA;
	double dB;
	double dR;

	double dREta;

	double dSqrtG;
	double dVolume;

	// Geometry
	const Geometry & geometry = block.GetGeometry();

	// Initialize parameters for Gaussian quadrature
	const int PointCount =
		GaussianQuadrature::GetGaussPointCount(
			GAUSSIAN_QUADRATURE_ORDER);

	DataVector<double> dG;
	dG.Initialize(PointCount);

	DataVector<double> dW;
	dW.Initialize(PointCount);

	GaussianQuadrature::GetGaussPointsWeights(
		GAUSSIAN_QUADRATURE_ORDER, dG, dW);

	DataVector<double> dStateSample;
	dStateSample.Initialize(block.GetComponents());

	DataVector<double> dTracerSample;
	dTracerSample.Initialize(Max(1, block.GetTracerCount()));

	double dWeight;

	// Set state to zero
	for (c = 0; c < block.GetComponents(); c++) {
		dState[c] = 0.0;
	}
	for (c = 0; c < block.GetTracerCount(); c++) {
		dTracer[c] = 0.0;
	}

	// REta coordinate at this level
	dREta = block.GetCentroidREta(iR);

	// Sample at all Gaussian quadrature points
	for (p = 0; p < PointCount; p++) {
	for (q = 0; q < PointCount; q++) {

		dA = block.GetCentroidA(iA) + 0.5 * block.GetDeltaA() * dG[p];
		dB = block.GetCentroidB(iB) + 0.5 * block.GetDeltaA() * dG[q];

		// Convert radial coordinate to true radius
		dR = geometry.ConvertREtaToRadius(
			block.GetPanel(), dREta, dA, dB);

		dWeight = 0.25 * dW[p] * dW[q];

		dStateSample.Zero();
		dTracerSample.Zero();

		EvaluatePointwiseState(
			block, phys, block.GetPanel(),
			dR, dA, dB,
			&(dStateSample[0]), &(dTracerSample[0]));

		dSqrtG = geometry.EvaluateSqrtG(iPanel, iR, dA, dB);

		for (c = 0; c < block.GetComponents(); c++) {
			dState[c] += dWeight * dStateSample[c] * dSqrtG;
		}
		for (c = 0; c < block.GetTracerCount(); c++) {
			dTracer[c] += dWeight * dTracerSample[c] * dSqrtG;
		}
	}
	}

	// Average
	dVolume = geometry.GetElementVolume(iR, iA, iB);
	for (c = 0; c < block.GetComponents(); c++) {
		dState[c] *= block.GetDeltaA() * block.GetDeltaA() / dVolume;
	}
	for (c = 0; c < block.GetTracerCount(); c++) {
		dTracer[c] *= block.GetDeltaA() * block.GetDeltaA() / dVolume;
	}

	// Subtract off hydrostatic background
	const SimulationParameters & params =
		block.GetSystemState().GetSimulationParameters();

	const BackgroundState & bgstate =
		block.GetBackgroundState();

	if (params.m_eEquationSet == SimulationParameters::EquationSet_Euler) {
		dState[0] -= bgstate.GetHydroRho(iR, iA, iB);
		dState[4] -= bgstate.GetHydroRhoTheta(iR, iA, iB);
	}
}

///////////////////////////////////////////////////////////////////////////////

void TestCase::InitialConditions(
	const Preferences & aPrefs,
	SystemState & state
) {
	int n;
	int k;
	int i;
	int j;

	int p;
	int q;

	int c;

	double dA;
	double dB;
	double dR;

	double dREta;

	double dSqrtG;
	double dVolume;

	// Initialize parameters for Gaussian quadrature
	const int PointCount =
		GaussianQuadrature::GetGaussPointCount(
			GAUSSIAN_QUADRATURE_ORDER);

	DataVector<double> dG;
	dG.Initialize(PointCount);

	DataVector<double> dW;
	dW.Initialize(PointCount);

	GaussianQuadrature::GetGaussPointsWeights(
		GAUSSIAN_QUADRATURE_ORDER, dG, dW);

	DataVector<double> dStateSample;
	dStateSample.Initialize(state.GetComponents());

	DataVector<double> dTracerSample;
	dTracerSample.Initialize(Max(1, state.GetTracerCount()));

	double dWeight;

	// Loop over all blocks
	for (n = 0; n < state.GetBlockCount(); n++) {
		SystemStateBlock & block = state.GetSystemStateBlock(n);

		const Geometry & geometry = block.GetGeometry();

		const PhysicalConstants & phys = state.GetPhysicalConstants();

		// Initialize data
		SystemStateData & data = block.GetData(0);
		data.Zero();

		SystemStateData * tracerdata;
		if (block.GetTracerCount() != 0) {
			tracerdata = &(block.GetTracerData(0));
			tracerdata->Zero();
		}

		// Loop over all elements within a block
		for (k = 0; k < block.GetRElements(); k++) {
		for (i = 0; i < block.GetATotalElements(); i++) {
		for (j = 0; j < block.GetBTotalElements(); j++) {

			// REta coordinate at this level
			dREta = block.GetCentroidREta(k);

			// Sample at all Gaussian quadrature points
			for (p = 0; p < PointCount; p++) {
			for (q = 0; q < PointCount; q++) {

				dA = block.GetCentroidA(i) + 0.5 * block.GetDeltaA() * dG[p];
				dB = block.GetCentroidB(j) + 0.5 * block.GetDeltaA() * dG[q];

				// Convert REta coordinate to radius
				dR = geometry.ConvertREtaToRadius(
					block.GetPanel(), dREta, dA, dB);

				dWeight = 0.25 * dW[p] * dW[q];

				EvaluatePointwiseState(
					block, phys, block.GetPanel(),
					dR, dA, dB,
					&(dStateSample[0]), &(dTracerSample[0]));

				dSqrtG = geometry.EvaluateSqrtG(block.GetPanel(), k, dA, dB);

				for (c = 0; c < state.GetComponents(); c++) {
					data[k][i][j][c] +=
						dWeight * dStateSample[c] * dSqrtG;
				}

				if (state.GetTracerCount() != 0) {
					for (c = 0; c < state.GetTracerCount(); c++) {
						(*tracerdata)[k][i][j][c] +=
							dWeight * dTracerSample[c] * dSqrtG;
					}
				}
			}
			}

			// Average
			dVolume = geometry.GetElementVolume(k, i, j);

			for (c = 0; c < state.GetComponents(); c++) {
				data[k][i][j][c] *=
					block.GetDeltaA() * block.GetDeltaA() / dVolume;
			}
			if (state.GetTracerCount() != 0) {
				for (c = 0; c < state.GetTracerCount(); c++) {
					(*tracerdata)[k][i][j][c] *=
						block.GetDeltaA() * block.GetDeltaA() / dVolume;
				}
			}

			// Subtract off hydrostatic background
			const SimulationParameters & params =
				block.GetSystemState().GetSimulationParameters();

			const BackgroundState & bgstate =
				block.GetBackgroundState();

			if (params.m_eEquationSet ==
				SimulationParameters::EquationSet_Euler
			) {
				data[k][i][j][0] -= bgstate.GetHydroRho(k, i, j);
				data[k][i][j][4] -= bgstate.GetHydroRhoTheta(k, i, j);
			}
		}
		}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void TestCase::PrescribeHorizontalVelocities(
	const SystemState & state,
	DataMatrix4D<double> & dAlphaU,
	DataMatrix4D<double> & dBetaU,
	double dTime
) {
	int n;
	int k;
	int i;
	int j;

	// Array of velocity components
	double dU[3];

	// Physical constants
	const PhysicalConstants & phys = state.GetPhysicalConstants();

	// Loop over all blocks
	for (n = 0; n < state.GetBlockCount(); n++) {
		SystemStateBlock & block = state.GetSystemStateBlock(n);

		const Geometry & geometry = block.GetGeometry();

		// Find velocities at alpha edge centerpoints
		for (k = 0; k < block.GetRElements(); k++) {
		for (i = 0; i <= block.GetATotalElements(); i++) {
		for (j = 0; j < block.GetBTotalElements(); j++) {

			double dA = block.GetEdgeA(i);
			double dB = block.GetCentroidB(j);
			double dREta = block.GetCentroidREta(k);

			double dR = geometry.ConvertREtaToRadius(
				block.GetPanel(), dREta, dA, dB);

			EvaluatePointwiseVelocity(
				block, phys, dTime,
				block.GetPanel(), dR, dA, dB,
				dU);

			// Determine component of U perpendicular to edge
			geometry.OrthonormalizeA(k, i, j, dU[0], dU[1]);

			dAlphaU[n][k][i][j] = dU[0];
		}
		}
		}

		// Find velocities at beta edge centerpoints
		for (k = 0; k < block.GetRElements(); k++) {
		for (i = 0; i < block.GetATotalElements(); i++) {
		for (j = 0; j <= block.GetBTotalElements(); j++) {

			double dA = block.GetCentroidA(i);
			double dB = block.GetEdgeB(j);
			double dREta = block.GetCentroidREta(k);

			double dR = geometry.ConvertREtaToRadius(
				block.GetPanel(), dREta, dA, dB);

			EvaluatePointwiseVelocity(
				block, phys, dTime,
				block.GetPanel(), dR, dA, dB,
				dU);

			// Determine component of U perpendicular to edge
			geometry.OrthonormalizeB(k, i, j, dU[0], dU[1]);

			dBetaU[n][k][i][j] = dU[0];
		}
		}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void TestCase::PrescribeVerticalVelocities(
	const SystemStateBlock & block,
	int iA,
	int iB,
	DataVector<double> & dVelocity,
	double dTime
) {
	int k;

	// Array of velocity components
	double dU[3];

	// Physical constants
	const PhysicalConstants & phys =
		block.GetSystemState().GetPhysicalConstants();

	const Geometry & geometry = block.GetGeometry();

	// Set velocity to zero at the top and bottom of model
	dVelocity[0] = 0.0;
	dVelocity[block.GetRElements()] = 0.0;

	// Find velocities at Reta edge centerpoints
	for (k = 1; k < block.GetRElements(); k++) {

		double dA = block.GetCentroidA(iA);
		double dB = block.GetCentroidB(iB);
		double dREta = block.GetEdgeREta(k);

		double dR = geometry.ConvertREtaToRadius(
			block.GetPanel(), dREta, dA, dB);

		EvaluatePointwiseVelocity(
			block, phys, dTime,
			block.GetPanel(), dR, dA, dB,
			dU);

		// Determine component of U perpendicular to edge
		dVelocity[k] =
			  geometry.GetROrthoCoeffs(iA, iB, k, 0) * dU[0]
			+ geometry.GetROrthoCoeffs(iA, iB, k, 1) * dU[1]
			+ geometry.GetROrthoCoeffs(iA, iB, k, 2) * dU[2];
	}
}
*/
///////////////////////////////////////////////////////////////////////////////

