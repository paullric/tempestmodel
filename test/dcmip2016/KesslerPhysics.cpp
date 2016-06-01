///////////////////////////////////////////////////////////////////////////////
///
///	\file    KesslerPhysics.cpp
///	\author  Antonin Verlet-Banide
///	\version May 3, 2016
///
///	<remarks>
///		Copyright 2000-2016 Paul Ullrich
///
///		This file is distributed as part of the Tempest source code package.
///		Permission is granted to use, copy, modify and distribute this
///		source code and its documentation under the terms of the GNU General
///		Public License.  This software is provided "as is" without express
///		or implied warranty.
///	</remarks>

#include "KesslerPhysics.h"

#include "Model.h"
#include "GridGLL.h"

#include "Announce.h" 

///////////////////////////////////////////////////////////////////////////////

extern "C" {
	void kessler_(
		double * t,
		double * qv,
		double * qc,
		double * qr,
		double * rho,
		double * pk,
		double * dt,
		double * z,
		int * nz,
		double * rainnc);
}

///////////////////////////////////////////////////////////////////////////////

KesslerPhysics::KesslerPhysics(
	Model & model,
	const Time & timeFrequency
) :
WorkflowProcess(
	model,
	timeFrequency)
{ }

//////////////////////////////////////////////////////////////////////////////

	void KesslerPhysics::Initialize(
		const Time & timeStart
) {
		// Indices of EquationSet variables
		const int UIx = 0;
		const int VIx = 1;
		const int TIx = 2;
		const int WIx = 3;
		const int RIx = 4;

		// Get a copy of the GLL grid
		Grid * pGrid = m_model.GetGrid();

		// Check position of variables
		int nRElements = pGrid->GetRElements();

		m_dQv.Allocate(nRElements);
		m_dQc.Allocate(nRElements);
		m_dQr.Allocate(nRElements);
		m_dRho.Allocate(nRElements);
		m_dZc.Allocate(nRElements);
		m_dPk.Allocate(nRElements);
		m_dTheta.Allocate(nRElements);

		m_dThetaVNode.Allocate(nRElements);
		m_dThetaVREdge.Allocate(nRElements+1);
}

	
///////////////////////////////////////////////////////////////////////////////
	
void KesslerPhysics::Perform(
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
	GridGLL * pGridGLL = dynamic_cast<GridGLL*>(m_model.GetGrid());
	if (pGridGLL == NULL) {
		_EXCEPTIONT("Not implemented for a general Grid");
	}

	// Physical constants
	const PhysicalConstants & phys = m_model.GetPhysicalConstants();

	// Number of radial elements
	int nRElements = pGridGLL->GetRElements();

	// Perform local update
	for (int n = 0; n < pGridGLL->GetActivePatchCount(); n++) {
		GridPatch * pPatch = pGridGLL->GetActivePatch(n);

		const PatchBox & box = pPatch->GetPatchBox();

		// Get latitude
		const DataArray2D<double> & dataLatitude = pPatch->GetLatitude();

		// Grid data
		DataArray4D<double> & dataNode =
			pPatch->GetDataState(0, DataLocation_Node);

		DataArray4D<double> & dataREdge =
			pPatch->GetDataState(0, DataLocation_REdge);

		DataArray4D<double> & dataTracer =
			pPatch->GetDataTracers(0);

		DataArray3D<double> & dataUserData2D =
			pPatch->GetUserData2D();

		if (dataUserData2D.GetRows() == 0) {
			_EXCEPTIONT("Insufficient entries in UserData2D");
		}

		const DataArray3D<double> & dataZLevels = pPatch->GetZLevels();

		// Loop over all horizontal nodes in GridPatch
		for (int i = box.GetAInteriorBegin(); i < box.GetAInteriorEnd(); i++) {
		for (int j = box.GetBInteriorBegin(); j < box.GetBInteriorEnd(); j++) {

			// Remap virtual potential temperature tendency to levels
			if (pGridGLL->GetVarLocation(TIx) == DataLocation_REdge) {
				for (int k = 0; k < nRElements+1; k++) {
					m_dThetaVREdge[k] = dataREdge[TIx][k][i][j];
				}

				pGridGLL->InterpolateREdgeToNode(
					&(m_dThetaVREdge[0]),
					&(m_dThetaVNode[0]));

			// Store virtual potential temperature on levels
			} else {
				for (int k = 0; k < nRElements; k++) {
					m_dThetaVNode[k] = dataNode[TIx][k][i][j];
				}
			}

			// Loop over all levels in column
			for (int k = 0; k < nRElements; k++) {

				// Total density
				double dRho = dataNode[RIx][k][i][j];

				// Dry air density
				double dRhoD =
					dRho
					- dataTracer[0][k][i][j]
					- dataTracer[1][k][i][j]
					- dataTracer[2][k][i][j];

				// Calculate moist pressure
				double dPressure =
					phys.PressureFromRhoTheta(dRho * m_dThetaVNode[k]);

				// Pointwise virtual temperature
				double dTv = dPressure / (dRho * phys.GetR());

				// Water vapor (RhoQv / Rho)
				m_dQv[k] = dataTracer[0][k][i][j] / dataNode[RIx][k][i][j];
				if (m_dQv[k] < 0.0) {
					m_dQv[k] = 0.0;
				}

				// Cloud water (RhoQc / Rho)
				m_dQc[k] = dataTracer[1][k][i][j] / dataNode[RIx][k][i][j];
				if (m_dQc[k] < 0.0) {
					m_dQc[k] = 0.0;
				}

				// Rain water (RhoQr / Rho)
				m_dQr[k] = dataTracer[2][k][i][j] / dataNode[RIx][k][i][j];
				if (m_dQr[k] < 0.0) {
					m_dQr[k] = 0.0;
				}

				// Potential temperature
				m_dTheta[k] = m_dThetaVNode[k] / (1.0 + 0.61 * m_dQv[k]);

				// Dry air density
				m_dRho[k] = dRhoD;
 
				// Heights of each level
				m_dZc[k] = dataZLevels[k][i][j];
	
				// Exner function (Tv/thetav)
				m_dPk[k] = dTv / m_dThetaVNode[k];
			}

			double dRainNc = 0.0;
			kessler_(
				&(m_dTheta[0]),
				&(m_dQv[0]),
				&(m_dQc[0]),
				&(m_dQr[0]),
				&(m_dRho[0]),
				&(m_dPk[0]),
				&(dDeltaT),
				&(m_dZc[0]),
				&(nRElements),
				&(dRainNc));

			// Store accumualted precipitation
			dataUserData2D[0][i][j] += dRainNc;

			// Remap virtual potential temperature tendency to interfaces
			if (pGridGLL->GetVarLocation(TIx) == DataLocation_REdge) {

				// Store negative virtual potential tendency
				for (int k = 0; k < nRElements; k++) {
					m_dThetaVNode[k] -=
						m_dTheta[k] * (1.0 + 0.61 * m_dQv[k]);
				}

				// Interpolate negative tendency
				pGridGLL->InterpolateNodeToREdge(
					&(m_dThetaVNode[0]),
					&(m_dThetaVREdge[0]));

				// Update thetaV on interfaces
				for (int k = 0; k <= nRElements; k++) {
					dataREdge[TIx][k][i][j] -= m_dThetaVREdge[k];
				}

			// Update thetaV on nodes
			} else {
				for (int k = 0; k < nRElements; k++) {
					dataNode[TIx][k][i][j] =
						m_dTheta[k] * (1.0 + 0.61 * m_dQv[k]);
				}
			}

			// Update mass variables
			for (int k = 0; k < nRElements; k++) {

				// Moist density
				dataNode[RIx][k][i][j] =
					m_dRho[k] / (1.0 - m_dQv[k] - m_dQc[k] - m_dQr[k]);

				// Virtual potential temperature
				dataTracer[0][k][i][j] = m_dQv[k] * dataNode[RIx][k][i][j];
				dataTracer[1][k][i][j] = m_dQc[k] * dataNode[RIx][k][i][j];
				dataTracer[2][k][i][j] = m_dQr[k] * dataNode[RIx][k][i][j];
			}
		}
		}
	}

	// Call up the stack to update performance time
	WorkflowProcess::Perform(time);
}

