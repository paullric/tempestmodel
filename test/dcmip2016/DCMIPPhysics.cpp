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

//#define DCMIP_KESSLER

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

	void dcmip2016_physics_(
		int * test,
		double * u,
		double * v,
		double * p,
		double * qv,
		double * qc,
		double * qr,
		double * rho,
		double * dt,
		double * z,
		double * zi,
		double * lat,
		int * nz,
		double * precl,
		int * pbl_type,
		int * prec_type);
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
		m_dZi.Allocate(nRElements+1);
		m_dPmid.Allocate(nRElements);
		m_dPint.Allocate(nRElements+1);
		m_dPdel.Allocate(nRElements);
		m_dRPdel.Allocate(nRElements);
		m_dT.Allocate(nRElements);
		m_dThetaVNode.Allocate(nRElements);
		m_dTheta.Allocate(nRElements);

		m_dU.Allocate(nRElements);
		m_dV.Allocate(nRElements);

		m_dTvNode.Allocate(nRElements);
		m_dTvREdge.Allocate(nRElements+1);

		m_dZNode.Allocate(nRElements);
		m_dZREdge.Allocate(nRElements+1);

		m_dEddyStateNode.Allocate(5, nRElements);
		m_dEddyStateREdge.Allocate(5, nRElements+1);

		m_dEddyTracerNode.Allocate(3, nRElements);
		m_dEddyTracerREdge.Allocate(3, nRElements+1);

		m_dKm.Allocate(nRElements+1);
		m_dKE.Allocate(nRElements+1);

		m_dEm.Allocate(nRElements);
		m_dEE.Allocate(nRElements);

		m_dF.Allocate(6, nRElements);

		m_dRhoREdge.Allocate(nRElements);
		m_dThetaREdge.Allocate(nRElements+1);
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

	// Not implemented for CPH
	if (pGridGLL->GetVarLocation(TIx) == DataLocation_REdge) {
		_EXCEPTIONT("Not implemented");
	}

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
		const DataArray3D<double> & dataZInterfaces = pPatch->GetZInterfaces();

		// Loop over all horizontal nodes in GridPatch
		for (int i = box.GetAInteriorBegin(); i < box.GetAInteriorEnd(); i++) {
		for (int j = box.GetBInteriorBegin(); j < box.GetBInteriorEnd(); j++) {
/*
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
*/
			// Latitude
			double dLat = dataLatitude[i][j];

			// Calculate quantities on model levels
			for (int k = 0; k < nRElements; k++) {

				// Virtual potential temperature
				m_dThetaVNode[k] = dataNode[TIx][k][i][j];

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

				// Total density
				double dRho = dataNode[RIx][k][i][j];

				// Dry air density
				double dRhoD =
					dRho
					- dataTracer[0][k][i][j]
					- dataTracer[1][k][i][j]
					- dataTracer[2][k][i][j];

				// Calculate moist pressure
				m_dPmid[k] =
					phys.PressureFromRhoTheta(dRho * m_dThetaVNode[k]);

				// Pointwise virtual temperature
				double dTv = m_dPmid[k] / (dRho * phys.GetR());

				// Water vapor mixing ratio (RhoQv / Rho)
				m_dQv[k] = dataTracer[0][k][i][j] / dataNode[RIx][k][i][j];
				if (m_dQv[k] < 0.0) {
					m_dQv[k] = 0.0;
				}

				// Cloud water mixing ratio (RhoQc / Rho)
				m_dQc[k] = dataTracer[1][k][i][j] / dataNode[RIx][k][i][j];
				if (m_dQc[k] < 0.0) {
					m_dQc[k] = 0.0;
				}

				// Rain water mixing ratio (RhoQr / Rho)
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
			}

			// Calculate quantities on model interfaces
			for (int k = 0; k <= nRElements; k++) {
				// Heights of each interface
				m_dZi[k] = dataZInterfaces[k][i][j];
			}

			int test = 1;
			int pbl_type = 0;
			int prec_type = 1;

			double dPrecL = 0.0;

			dcmip2016_physics_(
				&(test),
				&(m_dU[0]),
				&(m_dV[0]),
				&(m_dPmid[0]),
				&(m_dQv[0]),
				&(m_dQc[0]),
				&(m_dQr[0]),
				&(m_dRho[0]),
				&(dDeltaT),
				&(m_dZc[0]),
				&(m_dZi[0]),
				&(dLat),
				&(nRElements),
				&(dPrecL),
				&(pbl_type),
				&(prec_type));

			// Store precipitation
			dataUserData2D[0][i][j] += dPrecL;
/*
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
*/
			// Update mass variables
			for (int k = 0; k < nRElements; k++) {

				// Moist density
				dataNode[RIx][k][i][j] =
					m_dRho[k] / (1.0 - m_dQv[k] - m_dQc[k] - m_dQr[k]);

				// Virtual temperature
				double dTv =
					m_dPmid[k] / (dataNode[RIx][k][i][j] * phys.GetR());

				// Virtual potential temperature
				dataNode[TIx][k][i][j] =
					dTv * pow(phys.GetP0() / m_dPmid[k], phys.GetKappa());

				// Tracers
				dataTracer[0][k][i][j] = m_dQv[k] * dataNode[RIx][k][i][j];
				dataTracer[1][k][i][j] = m_dQc[k] * dataNode[RIx][k][i][j];
				dataTracer[2][k][i][j] = m_dQr[k] * dataNode[RIx][k][i][j];

				// Velocities
				CubedSphereTrans::CoVecTransABPFromRLL(
					tan(pPatch->GetANode(i)),
					tan(pPatch->GetBNode(j)),
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
	const int QIx = 5;

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

	if (pGridGLL->GetVarLocation(TIx) == DataLocation_REdge) {
		_EXCEPTIONT("CPH staggering not implemented (yet)");
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
		const DataArray3D<double> & dataZNode =
			pPatch->GetZLevels();
		const DataArray3D<double> & dataZREdge =
			pPatch->GetZInterfaces();

		// Loop over all horizontal nodes in GridPatch
		for (int i = box.GetAInteriorBegin(); i < box.GetAInteriorEnd(); i++) {
		for (int j = box.GetBInteriorBegin(); j < box.GetBInteriorEnd(); j++) {

			///////////////////////////////////////////////////////////////////
			// Large-scale precipitation
			double dPrecL = 0.0;

			for (int k = 0; k < nRElements; k++) {
				double dRho = dataNode[RIx][k][i][j];
				double dTheta = dataNode[TIx][k][i][j];

				double dP = phys.PressureFromRhoTheta(dRho * dTheta);

				double dTv = dP / (phys.GetR() * dRho);

				double dQv = dataTracer[0][k][i][j] / dRho;

				double dT = dTv / (1.0 + 0.61 * dQv);

				double dQsat =
					0.622 * 610.78 / dP
					* exp(- (phys.GetLvap() / phys.GetRvap())
						* ((1.0 / dT) - (1.0 / 273.16)));

				if (dQv > dQsat) {
					double dDeltaQv = (dQv - dQsat)
						/ (1.0 + (phys.GetLvap() / phys.GetCp())
							* (0.622 * phys.GetLvap() * dQsat
								/ (phys.GetR() * dT * dT)));

					dQv -= dDeltaQv;

					double dDeltaZ =
						(dataZREdge[k+1][i][j] - dataZREdge[k][i][j]);

					dPrecL +=
						dDeltaQv * dataNode[RIx][k][i][j] * dDeltaZ / 1000.0;

					double dDeltaTemp =
						phys.GetLvap() / phys.GetCp() * dDeltaQv;

					dTv = (dT + dDeltaTemp) * (1.0 + 0.61 * dQv);

					dP = dataNode[RIx][k][i][j] * phys.GetR() * dTv;

					dataTracer[0][k][i][j] = dQv * dRho;

					dataNode[TIx][k][i][j] =
						phys.RhoThetaFromPressure(dP) / dRho;
				}
			}

			///////////////////////////////////////////////////////////////////
			// Surface fluxes and Boundary layer

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
					dataNode[RIx][0][i][j]
					* dataNode[TIx][0][i][j]);

			double dTva = dPa / (phys.GetR() * dataNode[RIx][0][i][j]);

			double dQva = dataTracer[0][0][i][j] / dataNode[RIx][0][i][j];

			double dTa = dTva / (1.0 + 0.61 * dQva);

			// Surface pressure
			double dPsurf =
				phys.PressureFromRhoTheta(
					dataREdge[RIx][0][i][j] * dataREdge[TIx][0][i][j]);

			// Saturation specific humidity at surface
			double dQsatsurf =
				0.622 / dPsurf * 610.78
				* exp(-(phys.GetLvap() / phys.GetRvap())
					* ((1.0 / dTsurf) - (1.0 / 273.16)));

			// Store primitive quantities
			for (int k = 0; k < nRElements; k++) {
				m_dU[k] = dataNode[UIx][k][i][j];
				m_dV[k] = dataNode[VIx][k][i][j];
				m_dRho[k] = dataNode[RIx][k][i][j];

				m_dQv[k] = dataTracer[0][k][i][j] / m_dRho[k];

				m_dTheta[k] = dataNode[TIx][k][i][j] / (1.0 + 0.61 * m_dQv[k]);

				m_dZNode[k] = dataZNode[k][i][j];
			}

			// Altitude of lowest model level
			double dZa = m_dZNode[0];

			// Quantities on interfaces
			for (int k = 0; k <= nRElements; k++) {
				m_dRhoREdge[k] = dataREdge[RIx][k][i][j];

				m_dZREdge[k] = dataZREdge[k][i][j];

				double dPint =
					phys.PressureFromRhoTheta(
						dataREdge[RIx][k][i][j]
						* dataREdge[TIx][k][i][j]);

				// Mixing coefficients
				const double dPtop = 85000.0;
				const double dPstrato = 10000.0;

				if (dPint > 85000.0) {
					m_dKm[k] = dCd * dAbsU * dZa;
					m_dKE[k] = dCE * dAbsU * dZa;
				} else {
					double dDecay = (dPtop - dPint) / dPstrato;

					dDecay = exp(- dDecay * dDecay);

					m_dKm[k] = dCd * dAbsU * dZa * dDecay;
					m_dKE[k] = dCE * dAbsU * dZa * dDecay;
				}

				// DEBUG
				if (m_dKm[k] < 0.0) {
					_EXCEPTIONT("Negative Km detected");
				}
				if (m_dKE[k] < 0.0) {
					_EXCEPTIONT("Negative KE detected");
				}
			}

			// Apply surface fluxes (implicitly)
			m_dU[0] /=
				(1.0 + dDeltaT * dCd * dAbsU / dZa);
			m_dV[0] /=
				(1.0 + dDeltaT * dCd * dAbsU / dZa);
			m_dQv[0] =
				(m_dQv[0] + dCH * dAbsU * dQsatsurf * dDeltaT / dZa)
				/ (1.0 + dCH * dAbsU * dDeltaT / dZa);
			m_dTheta[0] =
				(dTa + dCH * dAbsU * dTsurf * dDeltaT / dZa)
				/ (1.0 + dCH * dAbsU * dDeltaT / dZa)
				* pow(phys.GetP0() / dPa, phys.GetKappa());

			// Calculate boundary layer parameterization matrix coefficients
			for (int k = 0; k < nRElements; k++) {

				double dA = 0.0;
				double dC = 0.0;

				if (k != 0) {
					dA = m_dRhoREdge[k] / m_dRho[k] * dDeltaT
						/ (m_dZNode[k] - m_dZNode[k-1])
						/ (m_dZREdge[k+1] - m_dZREdge[k]);
				}

				if (k != nRElements-1) {
					dC = m_dRhoREdge[k+1] / m_dRho[k] * dDeltaT
						/ (m_dZNode[k+1] - m_dZNode[k])
						/ (m_dZREdge[k+1] - m_dZREdge[k]);
				}

				double dAm = dA * m_dKm[k];
				double dCm = dC * m_dKm[k+1];

				double dAE = dA * m_dKE[k];
				double dCE = dC * m_dKE[k+1];

				double dBm = 1.0 + dAm + dCm;
				double dBE = 1.0 + dAE + dCE;

				if (k == 0) {
					m_dEm[0] = dCm / dBm;
					m_dEE[0] = dCE / dBE;

					m_dF[UIx][0] = m_dU[0] / dBm;
					m_dF[VIx][0] = m_dV[0] / dBm;
					m_dF[TIx][0] = m_dTheta[0] / dBE;
					m_dF[QIx][0] = m_dQv[0] / dBE;

				} else {
					m_dEm[k] = dCm / (dBm - dAm * m_dEm[k-1]);
					m_dEE[k] = dCE / (dBE - dAE * m_dEE[k-1]);

					m_dF[UIx][k] = (m_dU[k] + dAm * m_dF[UIx][k-1])
						/ (dBm - dAm * m_dEm[k-1]);
					m_dF[VIx][k] = (m_dV[k] + dAm * m_dF[VIx][k-1])
						/ (dBm - dAm * m_dEm[k-1]);
					m_dF[TIx][k] = (m_dTheta[k] + dAE * m_dF[TIx][k-1])
						/ (dBE - dAE * m_dEE[k-1]);
					m_dF[QIx][k] = (m_dQv[k] + dAE * m_dF[QIx][k-1])
						/ (dBE - dAE * m_dEE[k-1]);
				}
			}

			dataNode[UIx][nRElements-1][i][j] = m_dF[UIx][nRElements-1];
			dataNode[VIx][nRElements-1][i][j] = m_dF[VIx][nRElements-1];
			dataNode[TIx][nRElements-1][i][j] = m_dF[TIx][nRElements-1];
			dataTracer[0][nRElements-1][i][j] = m_dF[QIx][nRElements-1];

			for (int k = nRElements-2; k >= 0; k--) {
				dataNode[UIx][k][i][j] =
					m_dEm[k] * dataNode[UIx][k+1][i][j] + m_dF[UIx][k];
				dataNode[VIx][k][i][j] =
					m_dEm[k] * dataNode[VIx][k+1][i][j] + m_dF[VIx][k];

				dataTracer[0][k][i][j] =
					m_dEE[k] * dataTracer[0][k+1][i][j] + m_dF[QIx][k];

				dataNode[TIx][k][i][j] =
					m_dEE[k] * dataNode[TIx][k+1][i][j] + m_dF[TIx][k];
			}

			for (int k = 0; k < nRElements; k++) {
				m_dTheta[k] = dataNode[TIx][k][i][j];

				dataNode[TIx][k][i][j] *=
					(1.0 + 0.61 * dataTracer[0][k][i][j]);

				dataTracer[0][k][i][j] *= dataNode[RIx][k][i][j];
			}
		}
		}
	}

	//_EXCEPTION();
	// Call up the stack to update performance time
	WorkflowProcess::Perform(time);
}
*/
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

#if defined(DCMIP_KESSLER)
				// Dry air density
				double dRhoD =
					dRho
					- dataTracer[0][k][i][j]
					- dataTracer[1][k][i][j]
					- dataTracer[2][k][i][j];
#else
				// Dry air density
				double dRhoD =
					dRho - dataTracer[0][k][i][j];
#endif

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
			int iRJ2012_precip = 1;

			// Boundary layer modification
			int iTC_PBL_mod = 1;

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
*/
