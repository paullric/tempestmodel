///////////////////////////////////////////////////////////////////////////////
///
///	\file    HeldSuarezPhysics.cpp
///	\author  Paul Ullrich
///	\version December 15, 2011
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

#include "HeldSuarezPhysics.h"

#include "Model.h"

#include "Announce.h"

///////////////////////////////////////////////////////////////////////////////

// Boundary layer thickness (Sigma)
static const double ParamBoundarySigma = 0.7;

// Coefficient of friction in boundary layer (1/s)
static const double ParamKFriction = 1.0 / 86400.0;

// Polar coefficient of temperature diffusion (1/s)
static const double ParamKA = (1.0 / 40.0) / 86400.0;

// Equatorial coefficient of temperature diffusion (1/s)
static const double ParamKS = (1.0 / 4.0)  / 86400.0;

// Meridional temperature scale
static const double ParamDeltaTy = 60.0;

// Vertical potential temperature scale
static const double ParamDeltaThetaZ = 10.0;

// Lower bound on temperature
static const double ParamMinimumT = 200.0;

// Upper bound on temperature
static const double ParamMaximumT = 315.0;

///////////////////////////////////////////////////////////////////////////////

HeldSuarezPhysics::HeldSuarezPhysics(
	Model & model,
	const Time & timeFrequency
) :
	WorkflowProcess(
		model,
		timeFrequency)
{ }

///////////////////////////////////////////////////////////////////////////////

void HeldSuarezPhysics::Perform(
	const Time & time
) {
	// Indices of EquationSet variables
	const int UIx = 0;
	const int VIx = 1;
	const int TIx = 2;
	const int WIx = 3;
	const int RIx = 4;

	// Get DeltaT
	double dDeltaT = m_timeFrequency.GetSeconds();

	// Get a copy of the GLL grid
	Grid * pGrid = m_model.GetGrid();

	// Physical constants
	const PhysicalConstants & phys = m_model.GetPhysicalConstants();

	// Perform local update
	for (int n = 0; n < pGrid->GetActivePatchCount(); n++) {
		GridPatch * pPatch = pGrid->GetActivePatch(n);

		const PatchBox & box = pPatch->GetPatchBox();

		// Get latitude
		const DataArray2D<double> & dataLatitude = pPatch->GetLatitude();

		// Grid data
		DataArray4D<double> & dataNode =
			pPatch->GetDataState(0, DataLocation_Node);

		DataArray4D<double> & dataREdge =
			pPatch->GetDataState(0, DataLocation_REdge);

		// Perform interpolations as required due to vertical staggering
		if (pGrid->GetVarsAtLocation(DataLocation_REdge) != 0) {
			if (pGrid->GetVarLocation(TIx) == DataLocation_REdge) {
				pPatch->InterpolateREdgeToNode(TIx, 0);
				pPatch->InterpolateNodeToREdge(RIx, 0);
			}
		}

		// Loop over all horizontal nodes in GridPatch
		for (int i = box.GetAInteriorBegin(); i < box.GetAInteriorEnd(); i++) {
		for (int j = box.GetBInteriorBegin(); j < box.GetBInteriorEnd(); j++) {

			// Calculate surface pressure
			double dSurfacePressure =
				phys.PressureFromRhoTheta(
					dataREdge[RIx][0][i][j] * dataREdge[TIx][0][i][j]);

			// Loop over all levels in column
			for (int k = 0; k < pGrid->GetRElements(); k++) {

				// Calculate pressure
				double dPressure =
					phys.PressureFromRhoTheta(
						dataNode[RIx][k][i][j] * dataNode[TIx][k][i][j]);

				double dSigma = dPressure / dSurfacePressure;

				double dBoundaryScale =
					(dSigma - ParamBoundarySigma) / (1.0 - ParamBoundarySigma);

				if (dBoundaryScale < 0.0) {
					dBoundaryScale = 0.0;
				}

				// Apply velocity diffusion using backward Euler
				dataNode[UIx][k][i][j] /=
					(1.0 + ParamKFriction * dBoundaryScale * dDeltaT);
				dataNode[VIx][k][i][j] /=
					(1.0 + ParamKFriction * dBoundaryScale * dDeltaT);

			}

			// Theta on model levels
			if (pGrid->GetVarLocation(TIx) == DataLocation_Node) {
				for (int k = 0; k < pGrid->GetRElements(); k++) {

					// Calculate pressure
					double dPressure =
						phys.PressureFromRhoTheta(
							dataNode[RIx][k][i][j]
							* dataNode[TIx][k][i][j]);

					double dSigma = dPressure / dSurfacePressure;

					double dBoundaryScale =
						(dSigma - ParamBoundarySigma)
						/ (1.0 - ParamBoundarySigma);

					if (dBoundaryScale < 0.0) {
						dBoundaryScale = 0.0;
					}

					// Pointwise temperature
					double dT = dPressure
						/ (dataNode[RIx][k][i][j] * phys.GetR());

					// Get latitude
					double dLat = dataLatitude[i][j];

					// Temperature diffusion rate
					double dSinLat = sin(dLat);
					double dCosLat = cos(dLat);
					double dCos4Lat = dCosLat * dCosLat * dCosLat * dCosLat;
					double dKT = ParamKA
						+ (ParamKS - ParamKA) * dBoundaryScale * dCos4Lat;

					// Equilibrium temperature
					double dTeq =
						ParamMaximumT
						- ParamDeltaTy * dSinLat * dSinLat
						- ParamDeltaThetaZ * log(dPressure / phys.GetP0())
							* dCosLat * dCosLat;

					dTeq *= pow(dPressure / phys.GetP0(), phys.GetKappa());

					if (dTeq < ParamMinimumT) {
						dTeq = ParamMinimumT;
					}
/*
					// Apply temperature diffusion via backward Euler
					double dTnew =
						(dT + dDeltaT * dKT * dTeq) / (1.0 + dDeltaT * dKT);

					dataNode[TIx][k][i][j] =
						dTnew * pow(phys.GetP0() / dPressure, phys.GetKappa());
*/
					double dDH = - dKT / phys.GetGamma()
						* (1.0 + (phys.GetGamma() - 1.0) * dTeq / dT);

					double dH = - dKT / phys.GetGamma() * (1.0 - dTeq / dT);

					dataNode[TIx][k][i][j] *=
						1.0 + dDeltaT / (1.0 - dDeltaT * dDH) * dH;
				}
			}

			// Theta on model interfaces
			if (pGrid->GetVarLocation(TIx) == DataLocation_REdge) {
				for (int k = 0; k <= pGrid->GetRElements(); k++) {

					// Calculate pressure
					double dPressure =
						phys.PressureFromRhoTheta(
							dataREdge[RIx][k][i][j]
							* dataREdge[TIx][k][i][j]);

					double dSigma = dPressure / dSurfacePressure;

					double dBoundaryScale =
						(dSigma - ParamBoundarySigma)
						/ (1.0 - ParamBoundarySigma);

					if (dBoundaryScale < 0.0) {
						dBoundaryScale = 0.0;
					}

					// Pointwise temperature
					double dT = dPressure
						/ (dataREdge[RIx][k][i][j] * phys.GetR());

					// Get latitude
					double dLat = dataLatitude[i][j];

					// Temperature diffusion rate
					double dSinLat = sin(dLat);
					double dCosLat = cos(dLat);
					double dCos4Lat = dCosLat * dCosLat * dCosLat * dCosLat;
					double dKT = ParamKA
						+ (ParamKS - ParamKA) * dBoundaryScale * dCos4Lat;

					// Equilibrium temperature
					double dTeq =
						ParamMaximumT
						- ParamDeltaTy * dSinLat * dSinLat
						- ParamDeltaThetaZ * log(dPressure / phys.GetP0())
							* dCosLat * dCosLat;

					dTeq *= pow(dPressure / phys.GetP0(), phys.GetKappa());

					if (dTeq < ParamMinimumT) {
						dTeq = ParamMinimumT;
					}
/*
					// Apply temperature diffusion via backward Euler
					double dTnew =
						(dT + dDeltaT * dKT * dTeq) / (1.0 + dDeltaT * dKT);

					dataREdge[TIx][k][i][j] =
						dTnew * pow(phys.GetP0() / dPressure, phys.GetKappa());
*/
					double dDH = - dKT / phys.GetGamma()
						* (1.0 + (phys.GetGamma() - 1.0) * dTeq / dT);

					double dH = - dKT / phys.GetGamma() * (1.0 - dTeq / dT);

					dataREdge[TIx][k][i][j] *=
						1.0 + dDeltaT / (1.0 - dDeltaT * dDH) * dH;
				}
			}
/*
			// Vertical Velocity diffusion
			dataREdge[WIx][k][i][j] /=
				(1.0 + ParamKFriction * dBoundaryScale * dDeltaT);
*/
		}
		}
	}

	// Call up the stack to update performance time
	WorkflowProcess::Perform(time);
}

///////////////////////////////////////////////////////////////////////////////

