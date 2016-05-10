///////////////////////////////////////////////////////////////////////////////
///
///	\file    DCMIPPhysics.cpp
///	\author  Paul Ullrich
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

#include "DCMIPPhysics.h"

#include "CubedSphereTrans.h"
#include "Model.h"
#include "GridGLL.h"

#include "Announce.h" 

///////////////////////////////////////////////////////////////////////////////

extern "C" {
	void kessler(
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

	void simple_physics(
		int * pcols,
		int * pver,
		double * dtime,
		double * lat,
		double * t,
		double * q,
		double * u,
		double * v,
		double * pmid,
		double * pint,
		double * pdel,
		double * rpdel,
		double * ps,
		double * precl,
		int * test,
		int * RJ2012_precip,
		int * TC_PBL_mod);
}

///////////////////////////////////////////////////////////////////////////////

DCMIPPhysics::DCMIPPhysics(
	Model & model,
	const Time & timeFrequency
) :
WorkflowProcess(
	model,
	timeFrequency)
{ }

//////////////////////////////////////////////////////////////////////////////

	void DCMIPPhysics::Initialize(
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
		m_dPmid.Allocate(nRElements);
		m_dPint.Allocate(nRElements+1);
		m_dPdel.Allocate(nRElements);
		m_dRPdel.Allocate(nRElements);
		m_dT.Allocate(nRElements);
		m_dTheta.Allocate(nRElements);

		m_dU.Allocate(nRElements);
		m_dV.Allocate(nRElements);

		m_dTvNode.Allocate(nRElements);
		m_dTvREdge.Allocate(nRElements+1);
}

	
///////////////////////////////////////////////////////////////////////////////
	
void DCMIPPhysics::Perform(
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

	// Interpolate rho and theta to both model interfaces and levels
	((Grid*)pGridGLL)->InterpolateNodeToREdge(RIx, 0);

	if (pGridGLL->GetVarLocation(TIx) == DataLocation_Node) {
		((Grid*)pGridGLL)->InterpolateNodeToREdge(TIx, 0);
	} else {
		((Grid*)pGridGLL)->InterpolateREdgeToNode(TIx, 0);
	}

	// Perform local update
	for (int n = 0; n < pGridGLL->GetActivePatchCount(); n++) {
		GridPatch * pPatch = pGridGLL->GetActivePatch(n);

		const PatchBox & box = pPatch->GetPatchBox();

		// Grid data
		DataArray4D<double> & dataNode =
			pPatch->GetDataState(0, DataLocation_Node);

		DataArray4D<double> & dataREdge =
			pPatch->GetDataState(0, DataLocation_REdge);

		DataArray4D<double> & dataTracer =
			pPatch->GetDataTracers(0);

		const DataArray2D<double> & dataLatitude =
			pPatch->GetLatitude();

		const DataArray3D<double> & dataZLevels = pPatch->GetZLevels();

		// Loop over all horizontal nodes in GridPatch
		for (int i = box.GetAInteriorBegin(); i < box.GetAInteriorEnd(); i++) {
		for (int j = box.GetBInteriorBegin(); j < box.GetBInteriorEnd(); j++) {

			// Latitude
			double dLat = dataLatitude[i][j];

			// Calculate pressure on interfaces
			for (int k = 0; k <= nRElements; k++) {
				m_dPint[k] =
					phys.PressureFromRhoTheta(
						dataREdge[RIx][k][i][j]
						* dataREdge[TIx][k][i][j]);
			}

			// Loop over all levels in column
			for (int k = 0; k < nRElements; k++) {

				// Total density
				double dRho = dataNode[RIx][k][i][j];

				// Dry air density
				double dRhoD =
					dRho
					- dataTracer[0][k][i][j];
					//- dataTracer[1][k][i][j]
					//- dataTracer[2][k][i][j];

				// Calculate moist pressure
				m_dPmid[k] =
					phys.PressureFromRhoTheta(
						dataNode[RIx][k][i][j]
						* dataNode[TIx][k][i][j]);

				// Layer thickness
				m_dPdel[k] = m_dPint[k+1] - m_dPint[k];

				// Reciprocal of layer thickness
				m_dRPdel[k] = 1.0 / m_dPdel[k];

				// Water vapor (RhoQv / Rho)
				m_dQv[k] = dataTracer[0][k][i][j] / dataNode[RIx][k][i][j];

				// Cloud water (RhoQc / Rho)
				//m_dQc[k] = dataTracer[1][k][i][j] / dataNode[RIx][k][i][j];

				// Rain water (RhoQr / Rho)
				//m_dQr[k] = dataTracer[2][k][i][j] / dataNode[RIx][k][i][j];

				// Pointwise virtual temperature
				m_dTvNode[k] = m_dPmid[k] / (dRho * phys.GetR());

				// Temperature
				m_dT[k] = m_dTvNode[k] / (1.0 + 0.61 * m_dQv[k]);

				// Dry air density
				m_dRho[k] = dRhoD;
 
				// Heights of each level
				m_dZc[k] = dataZLevels[k][i][j];

				// Zonal velocity
				double dUalpha =
					dataNode[UIx][k][i][j] / phys.GetEarthRadius();
				double dUbeta =
					dataNode[VIx][k][i][j] / phys.GetEarthRadius();

				CubedSphereTrans::CoVecTransRLLFromABP(
					tan(pPatch->GetANode(i)),
					tan(pPatch->GetBNode(j)),
					pPatch->GetPatchBox().GetPanel(),
					dUalpha,
					dUbeta,
					m_dU[k],
					m_dV[k]);
			}

			// Number of columns
			int nPcols = 1;

			// Surface pressure
			double dPS = m_dPint[0];

			// Precipitation
			double dPRECL = 0.0;

			// Test (0 = tropical cyclone, 1 = moist baroclinic wave)
			int iTest = 0;

			// Activate Reed Jablonowski precipitation
			int iRJ2012_precip = 1;

			// Boundary layer modification
			int iTC_PBL_mod = 0;


			// Call simple_physics
			simple_physics(
				&(nPcols),
				&(nRElements),
				&(dDeltaT),
				&(dLat),
				&(m_dT[0]),
				&(m_dQv[0]),
				&(m_dU[0]),
				&(m_dV[0]),
				&(m_dPmid[0]),
				&(m_dPint[0]),
				&(m_dPdel[0]),
				&(m_dRPdel[0]),
				&(dPS),
				&(dPRECL),
				&(iTest),
				&(iRJ2012_precip),
				&(iTC_PBL_mod));

/*
			double dRainNc = 0.0;
			kessler(
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
*/

			// Update thetaV on interfaces
			if (pGridGLL->GetVarLocation(TIx) == DataLocation_REdge) {

				// Store negative virtual temperature tendency
				for (int k = 0; k < nRElements; k++) {
					m_dTvNode[k] -= m_dT[k] * (1.0 + 0.61 * m_dQv[k]);
				}

				// Interpolate negative tendency
				pGridGLL->InterpolateNodeToREdge(
					&(m_dTvNode[0]),
					&(m_dTvREdge[0]));

				// Update thetaV on interfaces
				for (int k = 0; k <= nRElements; k++) {
					dataREdge[TIx][k][i][j] -=
						m_dTvREdge[k]
						* pow(phys.GetP0() / m_dPint[k], phys.GetKappa());
				}

			// Update thetaV on nodes
			} else {
				for (int k = 0; k < nRElements; k++) {
					dataNode[TIx][k][i][j] =
						m_dT[k] * (1.0 + 0.61 * m_dQv[k])
						* pow(phys.GetP0() / m_dPmid[k], phys.GetKappa());
				}
			}

			// Update mass variables
			for (int k = 0; k < nRElements; k++) {

				// Moist density
				dataNode[RIx][k][i][j] =
					m_dRho[k] / (1.0 - m_dQv[k]); // - m_dQc[k] - m_dQr[k]);

				// Water vapor density
				dataTracer[0][k][i][j] = m_dQv[k] * dataNode[RIx][k][i][j];
				//dataTracer[1][k][i][j] = m_dQc[k] * dataNode[RIx][k][i][j];
				//dataTracer[2][k][i][j] = m_dQr[k] * dataNode[RIx][k][i][j];

				CubedSphereTrans::CoVecTransABPFromRLL(
					tan(pPatch->GetANode(i)),
					tan(pPatch->GetBNode(i)),
					pPatch->GetPatchBox().GetPanel(),
					m_dU[k],
					m_dV[k],
					dataNode[UIx][k][i][j],
					dataNode[VIx][k][i][j]);

				dataNode[UIx][k][i][j] *= phys.GetEarthRadius();
				dataNode[VIx][k][i][j] *= phys.GetEarthRadius();
			}
		}
		}
	}

	// Call up the stack to update performance time
	WorkflowProcess::Perform(time);
}

