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

		m_dEddyStateNode.Allocate(5, nRElements);
		m_dEddyStateREdge.Allocate(5, nRElements+1);

		m_dEddyTracerNode.Allocate(3, nRElements);
		m_dEddyTracerREdge.Allocate(3, nRElements+1);
}

		
///////////////////////////////////////////////////////////////////////////////
/*
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

		// Metric components
		const DataArray3D<double> & dContraMetric2DA =
			pPatch->GetContraMetric2DA();
		const DataArray3D<double> & dContraMetric2DB =
			pPatch->GetContraMetric2DB();

		// Grid data
		DataArray4D<double> & dataNode =
			pPatch->GetDataState(0, DataLocation_Node);
		DataArray4D<double> & dataREdge =
			pPatch->GetDataState(0, DataLocation_REdge);
		DataArray4D<double> & dataTracer =
			pPatch->GetDataTracers(0);

		const DataArray2D<double> & dataLatitude =
			pPatch->GetLatitude();
		const DataArray3D<double> & dataZLevels =
			pPatch->GetZLevels();

		// Loop over all horizontal nodes in GridPatch
		for (int i = box.GetAInteriorBegin(); i < box.GetAInteriorEnd(); i++) {
		for (int j = box.GetBInteriorBegin(); j < box.GetBInteriorEnd(); j++) {

			// Contravariant velocity in lowest model level
			double dCovUa = dataNode[UIx][0][i][j];
			double dCovUb = dataNode[VIx][0][i][j];

			// Contravariant velocity in lowest model level
			double dConUa =
				  dContraMetric2DA[i][j][0] * dCovUa
				+ dContraMetric2DA[i][j][1] * dCovUb;

			double dConUb =
				  dContraMetric2DB[i][j][0] * dCovUa
				+ dContraMetric2DB[i][j][1] * dCovUb;

			// Absolute velocity in lowest model level
			double dAbsU = sqrt(dCovUa * dConUa + dCovUb * dConUb);

			if (dAbsU < 0.0) {
				_EXCEPTIONT("Logic error");
			}
			if (dAbsU > 200.0) {
				_EXCEPTIONT("Logic error");
			}

			// Tropical cyclone surface temperatures
			double dTsurf = 302.15;

			// Mixing coefficients for surface fluxes
			double dCd;
			if (dAbsU > 20.0) {
				dCd = 0.002;
			} else {
				dCd = 7.0e-4 + 6.5e-5 * dAbsU;
			}

			double dCH = 0.0011;
			double dCE = 0.0011;

			// Pressure and temperature in lowest model level
			double dPa =
				phys.PressureFromRhoTheta(
					dataNode[RIx][0][i][j] * dataNode[TIx][0][i][j]);

			double dTva = dPa / (phys.GetR() * dataNode[RIx][0][i][j]);

			double dQa = dataTracer[0][0][i][j] / dataNode[RIx][0][i][j];

			double dTa = dTva / (1.0 + 0.61 * dQa);

			// Surface pressure
			double dPsurf =
				phys.PressureFromRhoTheta(
					dataREdge[RIx][0][i][j] * dataREdge[TIx][0][i][j]);

			// Saturation specific humidity at surface
			double dQsatsurf =
				0.622 / dPsurf * 610.78
				* exp(-(2.5e6 / 461.5) * ((1.0 / dTsurf) - (1.0 / 273.16)));

			// Altitude of lowest model level
			double dZa = dataZLevels[0][i][j];

			// Eddy fluxes at each node
			for (int k = 0; k < nRElements; k++) {
				m_dU[k] = dataNode[UIx][k][i][j];
				m_dV[k] = dataNode[VIx][k][i][j];
				m_dTheta[k] = dataNode[TIx][k][i][j];

				m_dQv[k] = dataTracer[0][k][i][j];
			}

			pGridGLL->DifferentiateNodeToREdge(
				&(m_dU[0]), &(m_dEddyStateREdge[UIx][0]));
			pGridGLL->DifferentiateNodeToREdge(
				&(m_dV[0]), &(m_dEddyStateREdge[VIx][0]));
			pGridGLL->DifferentiateNodeToREdge(
				&(m_dTheta[0]), &(m_dEddyStateREdge[TIx][0]));

			pGridGLL->DifferentiateNodeToREdge(
				&(m_dQv[0]), &(m_dEddyTracerREdge[0][0]));

			for (int k = 0; k <= nRElements; k++) {
				double dRho = dataREdge[RIx][k][i][j];
				double dThetaV = dataREdge[TIx][k][i][j];

				m_dPint[k] =
					phys.PressureFromRhoTheta(dRho * dThetaV);

				double dKm;
				double dKE;

				const double dPtop = 85000.0;
				const double dPstrato = 10000.0;

				if (m_dPint[k] > 85000.0) {
					dKm = dCd * dAbsU * dZa;
					dKE = dCE * dAbsU * dZa;
				} else {
					double dDecay = (dPtop - m_dPint[k]) / dPstrato;

					dDecay = exp(- dDecay * dDecay);

					dKm = dCd * dAbsU * dZa * dDecay;
					dKE = dCE * dAbsU * dZa * dDecay;
				}

				m_dEddyStateREdge[UIx][k] *= - dRho * dKm;
				m_dEddyStateREdge[VIx][k] *= - dRho * dKm;
				m_dEddyStateREdge[TIx][k] *= - dRho * dKE;
				m_dEddyTracerREdge[0][k] *= - dRho * dKE;
			}

			pGridGLL->DifferentiateREdgeToNode(
				&(m_dEddyStateREdge[UIx][0]), &(m_dEddyStateNode[UIx][0]));
			pGridGLL->DifferentiateREdgeToNode(
				&(m_dEddyStateREdge[VIx][0]), &(m_dEddyStateNode[VIx][0]));
			pGridGLL->DifferentiateREdgeToNode(
				&(m_dEddyStateREdge[TIx][0]), &(m_dEddyStateNode[TIx][0]));

			pGridGLL->DifferentiateREdgeToNode(
				&(m_dEddyTracerREdge[0][0]), &(m_dEddyTracerNode[0][0]));

			// Apply surface fluxes
			dataNode[UIx][0][i][j] -=
				dDeltaT * dCd * dAbsU * dCovUa / dZa;
			dataNode[VIx][0][i][j] -=
				dDeltaT * dCd * dAbsU * dCovUb / dZa;

			dataNode[TIx][0][i][j] +=
				dDeltaT * dCH * dAbsU * (dTsurf - dTa) / dZa
				* pow(phys.GetP0() / dPa, phys.GetKappa());

			dataTracer[0][0][i][j] +=
				dDeltaT * dCE * dAbsU * (dQsatsurf - dQa) / dZa;

			// Apply boundary layer parameterization
			for (int k = 0; k < nRElements; k++) {

				if ((dAbsU > 10.0) && (k < 10)) {
					printf("%i %1.15e\n", k,
						dataNode[UIx][k][i][j]
						/ (dDeltaT * m_dEddyStateNode[UIx][k] / dataNode[RIx][k][i][j]));
				}

				dataNode[UIx][k][i][j] -=
					dDeltaT * m_dEddyStateNode[UIx][k] / dataNode[RIx][k][i][j];
				dataNode[VIx][k][i][j] -=
					dDeltaT * m_dEddyStateNode[VIx][k] / dataNode[RIx][k][i][j];
				dataNode[TIx][k][i][j] -=
					dDeltaT * m_dEddyStateNode[TIx][k] / dataNode[RIx][k][i][j];

				dataTracer[0][k][i][j] -=
					dDeltaT * m_dEddyTracerNode[0][k] / dataNode[RIx][k][i][j];
			}
		}
		}
	}

	_EXCEPTION();
	// Call up the stack to update performance time
	WorkflowProcess::Perform(time);
}
*/
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

		const DataArray3D<double> & dataZLevels =
			pPatch->GetZLevels();

		// Loop over all horizontal nodes in GridPatch
		for (int i = box.GetAInteriorBegin(); i < box.GetAInteriorEnd(); i++) {
		for (int j = box.GetBInteriorBegin(); j < box.GetBInteriorEnd(); j++) {

			// Latitude
			double dLat = dataLatitude[i][j];

			// Calculate pressure on interfaces
			for (int k = 0; k <= nRElements; k++) {
				int kswap = nRElements - k;

				m_dPint[kswap] =
					phys.PressureFromRhoTheta(
						dataREdge[RIx][k][i][j]
						* dataREdge[TIx][k][i][j]);
			}

			// Loop over all levels in column
			for (int k = 0; k < nRElements; k++) {

				int kswap = nRElements - k - 1;

				// Total density
				double dRho = dataNode[RIx][k][i][j];

				// Dry air density
				double dRhoD =
					dRho
					- dataTracer[0][k][i][j];
					//- dataTracer[1][k][i][j]
					//- dataTracer[2][k][i][j];

				// Calculate moist pressure
				m_dPmid[kswap] =
					phys.PressureFromRhoTheta(
						dataNode[RIx][k][i][j]
						* dataNode[TIx][k][i][j]);

				// Layer thickness
				m_dPdel[kswap] = m_dPint[kswap+1] - m_dPint[kswap];

				if (m_dPdel[kswap] <= 0.0) {
					_EXCEPTIONT("Logic error");
				}

				// Reciprocal of layer thickness
				m_dRPdel[kswap] = 1.0 / m_dPdel[kswap];

				// Water vapor (RhoQv / Rho)
				m_dQv[kswap] = dataTracer[0][k][i][j] / dataNode[RIx][k][i][j];

				// Cloud water (RhoQc / Rho)
				//m_dQc[k] = dataTracer[1][k][i][j] / dataNode[RIx][k][i][j];

				// Rain water (RhoQr / Rho)
				//m_dQr[k] = dataTracer[2][k][i][j] / dataNode[RIx][k][i][j];

				// Pointwise virtual temperature (bottom to top)
				m_dTvNode[k] = m_dPmid[kswap] / (dRho * phys.GetR());

				// Temperature
				m_dT[kswap] = m_dTvNode[k] / (1.0 + 0.61 * m_dQv[kswap]);

				// Dry air density
				m_dRho[kswap] = dRhoD;
 
				// Heights of each level
				m_dZc[kswap] = dataZLevels[k][i][j];

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
					m_dU[kswap],
					m_dV[kswap]);
			}

			// Number of columns
			int nPcols = 1;

			// Surface pressure
			double dPS = m_dPint[nRElements];

			// Precipitation
			double dPRECL = 0.0;

			// Test (0 = tropical cyclone, 1 = moist baroclinic wave)
			int iTest = 0;

			// Activate Reed Jablonowski precipitation
			int iRJ2012_precip = 0;

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

			// Update thetaV on interfaces
			if (pGridGLL->GetVarLocation(TIx) == DataLocation_REdge) {

				// Store negative virtual temperature tendency
				for (int k = 0; k < nRElements; k++) {
					int kswap = nRElements - k - 1;

					m_dTvNode[k] -= m_dT[kswap] * (1.0 + 0.61 * m_dQv[kswap]);
				}

				// Interpolate negative tendency
				pGridGLL->InterpolateNodeToREdge(
					&(m_dTvNode[0]),
					&(m_dTvREdge[0]));

				// Update thetaV on interfaces
				for (int k = 0; k <= nRElements; k++) {
					int kswap = nRElements - k;

					dataREdge[TIx][k][i][j] -=
						m_dTvREdge[k]
						* pow(phys.GetP0() / m_dPint[kswap], phys.GetKappa());
				}

			// Update thetaV on nodes
			} else {
				for (int k = 0; k < nRElements; k++) {
					int kswap = nRElements - k - 1;

					dataNode[TIx][k][i][j] =
						m_dT[kswap] * (1.0 + 0.61 * m_dQv[kswap])
						* pow(phys.GetP0() / m_dPmid[kswap], phys.GetKappa());
				}
			}

			// Update mass variables
			for (int k = 0; k < nRElements; k++) {
				int kswap = nRElements - k - 1;

				// Moist density
				dataNode[RIx][k][i][j] =
					m_dRho[kswap] / (1.0 - m_dQv[kswap]); // - m_dQc[k] - m_dQr[k]);

				// Water vapor density
				dataTracer[0][k][i][j] = m_dQv[kswap] * dataNode[RIx][k][i][j];
				//dataTracer[1][k][i][j] = m_dQc[k] * dataNode[RIx][k][i][j];
				//dataTracer[2][k][i][j] = m_dQr[k] * dataNode[RIx][k][i][j];

				CubedSphereTrans::CoVecTransABPFromRLL(
					tan(pPatch->GetANode(i)),
					tan(pPatch->GetBNode(j)),
					pPatch->GetPatchBox().GetPanel(),
					m_dU[kswap],
					m_dV[kswap],
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


