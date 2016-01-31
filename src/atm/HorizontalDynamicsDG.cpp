///////////////////////////////////////////////////////////////////////////////
///
///	\file    HorizontalDynamicsDG.cpp
///	\author  Paul Ullrich
///	\version September 19, 2013
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

#include "Defines.h"
#include "HorizontalDynamicsDG.h"
#include "PhysicalConstants.h"
#include "Model.h"
#include "Grid.h"

#include "GridGLL.h"
#include "GridPatchGLL.h"

#define DIFFERENTIAL_FORM

//#define ADVECTION_ONLY

#define PENALIZE_DISCONTINUITY

#ifdef DIFFERENTIAL_FORM
#pragma message "WARNING: DIFFERENTIAL_FORM will lose mass over topography"
#endif

///////////////////////////////////////////////////////////////////////////////

static const double ParamHypervisScaling = 3.2;

///////////////////////////////////////////////////////////////////////////////

HorizontalDynamicsDG::HorizontalDynamicsDG(
	Model & model,
	int nHorizontalOrder,
	int nHyperviscosityOrder,
	double dNuScalar,
	double dNuDiv,
	double dNuVort,
	double dInstepNuDiv
) :
	HorizontalDynamicsFEM(
		model,
		nHorizontalOrder,
		nHyperviscosityOrder,
		dNuScalar,
		dNuDiv,
		dNuVort,
		dInstepNuDiv)
{
}

///////////////////////////////////////////////////////////////////////////////

void HorizontalDynamicsDG::StepShallowWater(
	int iDataInitial,
	int iDataUpdate,
	const Time & time,
	double dDeltaT
) {
/*
	// Get a copy of the GLL grid
	GridGLL * pGrid = dynamic_cast<GridGLL*>(m_model.GetGrid());

	// Physical constants
	const PhysicalConstants & phys = m_model.GetPhysicalConstants();

	// Get indices of variables to update
	const int UIx = pGrid->GetVarIndex(0);
	const int VIx = pGrid->GetVarIndex(1);
	const int HIx = pGrid->GetVarIndex(2);

	// GLL Weights, used by DG dynamics
	const DataArray1D<double> & dGLLWeights1D = pGrid->GetGLLWeights1D();

	// Perform local update
	for (int n = 0; n < pGrid->GetActivePatchCount(); n++) {
		GridPatchGLL * pPatch =
			dynamic_cast<GridPatchGLL*>(pGrid->GetActivePatch(n));

		const PatchBox & box = pPatch->GetPatchBox();

		const DataArray3D<double> & dJacobian =
			pPatch->GetJacobian();
		const DataArray4D<double> & dContraMetricA =
			pPatch->GetContraMetricA();
		const DataArray4D<double> & dContraMetricB =
			pPatch->GetContraMetricB();
		const DataArray3D<double> & dChristoffelA =
			pPatch->GetChristoffelA();
		const DataArray3D<double> & dChristoffelB =
			pPatch->GetChristoffelB();
		const DataArray2D<double> & dCoriolisF =
			pPatch->GetCoriolisF();
		const DataArray2D<double> & dTopography =
			pPatch->GetTopography();

		// Connectivity for patch
		Connectivity & connect = pPatch->GetConnectivity();

		// Data
		const DataArray4D<double> & dataInitialNode =
			pPatch->GetDataState(iDataInitial, DataLocation_Node);

		DataArray4D<double> & dataUpdateNode =
			pPatch->GetDataState(iDataUpdate, DataLocation_Node);

		// Element grid spacing and derivative coefficients
		const double dElementDeltaA = pPatch->GetElementDeltaA();
		const double dElementDeltaB = pPatch->GetElementDeltaB();

		const DataArray2D<double> & dDxBasis1D = pGrid->GetDxBasis1D();
		const DataArray2D<double> & dStiffness1D = pGrid->GetStiffness1D();

		// Get number of finite elements in each coordinate direction
		int nElementCountA = pPatch->GetElementCountA();
		int nElementCountB = pPatch->GetElementCountB();

		for (int k = 0; k < pGrid->GetRElements(); k++) {
		for (int a = 0; a < nElementCountA; a++) {
		for (int b = 0; b < nElementCountB; b++) {

			// Pointwise fluxes within spectral element
			for (int i = 0; i < m_nHorizontalOrder; i++) {
			for (int j = 0; j < m_nHorizontalOrder; j++) {

				int iA = a * m_nHorizontalOrder + i + box.GetHaloElements();
				int iB = b * m_nHorizontalOrder + j + box.GetHaloElements();

				// Height flux
				m_dAlphaFlux[i][j] =
					dJacobian[k][iA][iB]
					* (dataInitialNode[HIx][k][iA][iB] - dTopography[iA][iB])
					* dataInitialNode[UIx][k][iA][iB];

				m_dBetaFlux[i][j] =
					dJacobian[k][iA][iB]
					* (dataInitialNode[HIx][k][iA][iB] - dTopography[iA][iB])
					* dataInitialNode[VIx][k][iA][iB];

				// Pointwise pressure
				m_dPressure[i][j] =
					phys.GetG() * dataInitialNode[HIx][k][iA][iB];
			}
			}

			// Pointwise update of horizontal velocities
			for (int i = 0; i < m_nHorizontalOrder; i++) {
			for (int j = 0; j < m_nHorizontalOrder; j++) {

				int iA = a * m_nHorizontalOrder + i + box.GetHaloElements();
				int iB = b * m_nHorizontalOrder + j + box.GetHaloElements();

				int iElementA = a * m_nHorizontalOrder + box.GetHaloElements();
				int iElementB = b * m_nHorizontalOrder + box.GetHaloElements();

				// Derivatives of the velocity field
				double dDaUa = 0.0;
				double dDaUb = 0.0;
				double dDbUa = 0.0;
				double dDbUb = 0.0;

				// Derivatives of the pressure field
				double dDaP = 0.0;
				double dDbP = 0.0;

				// Aliases for alpha and beta velocities
				double dUa = dataInitialNode[UIx][k][iA][iB];
				double dUb = dataInitialNode[VIx][k][iA][iB];

				// Update density
				double dLocalUpdateHa = 0.0;
				for (int s = 0; s < m_nHorizontalOrder; s++) {
#ifdef DIFFERENTIAL_FORM
					// Update density: Differential formulation
					dLocalUpdateHa -=
						m_dAlphaFlux[s][j]
						* dDxBasis1D[s][i];
#else
					// Update density: Variational formulation
					dLocalUpdateHa +=
						m_dAlphaFlux[s][j]
						* dStiffness1D[i][s];
#endif
					// Derivative of alpha velocity with respect to alpha
					dDaUa +=
						dataInitialNode[UIx][k][iElementA+s][iB]
						* dDxBasis1D[s][i];

					// Derivative of beta velocity with respect to alpha
					dDaUb +=
						dataInitialNode[VIx][k][iElementA+s][iB]
						* dDxBasis1D[s][i];

					// Derivative of pressure with respect to alpha
					dDaP +=
						m_dPressure[s][j]
						* dDxBasis1D[s][i];
				}

				double dLocalUpdateHb = 0.0;
				for (int s = 0; s < m_nHorizontalOrder; s++) {
#ifdef DIFFERENTIAL_FORM
					// Update density: Differential formulation
					dLocalUpdateHb -=
						m_dBetaFlux[i][s]
						* dDxBasis1D[s][j];
#else
					// Update density: Variational formulation
					dLocalUpdateHb +=
						m_dBetaFlux[i][s]
						* dStiffness1D[j][s];
#endif
					// Derivative of alpha velocity with respect to beta
					dDbUa +=
						dataInitialNode[UIx][k][iA][iElementB+s]
						* dDxBasis1D[s][j];

					// Derivative of beta velocity with respect to beta
					dDbUb +=
						dataInitialNode[VIx][k][iA][iElementB+s]
						* dDxBasis1D[s][j];

					// Derivative of pressure with respect to beta
					dDbP +=
						m_dPressure[i][s]
						* dDxBasis1D[s][j];
				}

				// Scale derivatives
				dDaUa /= dElementDeltaA;
				dDaUb /= dElementDeltaA;
				dDaP  /= dElementDeltaA;

				dDbUa /= dElementDeltaB;
				dDbUb /= dElementDeltaB;
				dDbP  /= dElementDeltaB;

				dLocalUpdateHa /= (dJacobian[k][iA][iB] * dElementDeltaA);
				dLocalUpdateHb /= (dJacobian[k][iA][iB] * dElementDeltaB);

				// Local update for momentum
				double dLocalUpdateUa = 0.0;
				double dLocalUpdateUb = 0.0;

				// Curvature terms
				dLocalUpdateUa -= 
						+ dChristoffelA[iA][iB][0] * dUa * dUa
						+ dChristoffelA[iA][iB][1] * dUa * dUb
						+ dChristoffelA[iA][iB][2] * dUb * dUb;

				dLocalUpdateUb -=
						+ dChristoffelB[iA][iB][0] * dUa * dUa
						+ dChristoffelB[iA][iB][1] * dUa * dUb
						+ dChristoffelB[iA][iB][2] * dUb * dUb;

				// Coriolis forces
				dLocalUpdateUa -=
					dCoriolisF[iA][iB] * dJacobian[k][iA][iB] * (
						+ dContraMetricA[k][iA][iB][1] * dUa
						- dContraMetricA[k][iA][iB][0] * dUb);

				dLocalUpdateUb -=
					dCoriolisF[iA][iB] * dJacobian[k][iA][iB] * (
						+ dContraMetricB[k][iA][iB][1] * dUa
						- dContraMetricB[k][iA][iB][0] * dUb);

				// Interior dynamics
				double dLocalUpdateUaDiff = 0.0;
				double dLocalUpdateUbDiff = 0.0;

				dLocalUpdateUaDiff +=
					- dUa * dDaUa - dContraMetricA[k][iA][iB][0] * dDaP
					- dUb * dDbUa - dContraMetricA[k][iA][iB][1] * dDbP;

				dLocalUpdateUbDiff +=
					- dUa * dDaUb - dContraMetricB[k][iA][iB][0] * dDaP
					- dUb * dDbUb - dContraMetricB[k][iA][iB][1] * dDbP;

#ifndef ADVECTION_ONLY
				// Apply update
				dataUpdateNode[UIx][k][iA][iB] +=
					dDeltaT * (dLocalUpdateUa + dLocalUpdateUaDiff);
				dataUpdateNode[VIx][k][iA][iB] +=
					dDeltaT * (dLocalUpdateUb + dLocalUpdateUbDiff);
#endif

				// Update free surface height
				dataUpdateNode[HIx][k][iA][iB] +=
					dDeltaT * (dLocalUpdateHa + dLocalUpdateHb);

			}
			}
		}
		}
		}
	}
*/
}

///////////////////////////////////////////////////////////////////////////////

void HorizontalDynamicsDG::ElementFluxesShallowWater(
	int iDataInitial,
	int iDataUpdate,
	const Time & time,
	double dDeltaT
) {
	// Get a copy of the GLL grid
	GridGLL * pGrid = dynamic_cast<GridGLL*>(m_model.GetGrid());

	// Physical constants
	const PhysicalConstants & phys = m_model.GetPhysicalConstants();

	// Indices of EquationSet variables
	const int UIx = 0;
	const int VIx = 1;
	const int HIx = 2;

	// Perform a global exchange
	pGrid->Exchange(DataType_State, iDataInitial);

	// Perform local update
	for (int n = 0; n < pGrid->GetActivePatchCount(); n++) {
		GridPatchGLL * pPatch =
			dynamic_cast<GridPatchGLL*>(pGrid->GetActivePatch(n));

		const PatchBox & box = pPatch->GetPatchBox();

		const DataArray3D<double> & dJacobian =
			pPatch->GetJacobian();
		const DataArray4D<double> & dContraMetricA =
			pPatch->GetContraMetricA();
		const DataArray4D<double> & dContraMetricB =
			pPatch->GetContraMetricB();
		const DataArray2D<double> & dTopography =
			pPatch->GetTopography();

		const DataArray1D<double> & dGLLWeights1D = pGrid->GetGLLWeights1D();
		const DataArray1D<double> & dFluxDeriv1D = pGrid->GetFluxDeriv1D();

		// Data
		DataArray4D<double> & dataInitialNode =
			pPatch->GetDataState(iDataInitial, DataLocation_Node);
		DataArray4D<double> & dataUpdateNode =
			pPatch->GetDataState(iDataUpdate, DataLocation_Node);

		// Element grid spacing and derivative coefficients
		const double dElementDeltaA = pPatch->GetElementDeltaA();
		const double dElementDeltaB = pPatch->GetElementDeltaB();

		// Time over grid spacing ratio
		double dCourantA = dDeltaT / dElementDeltaA;
		double dCourantB = dDeltaT / dElementDeltaB;

		// Get number of finite elements in each coordinate direction
		int nElementCountA = pPatch->GetElementCountA();
		int nElementCountB = pPatch->GetElementCountB();

		// Post-process velocities received during exchange
		pPatch->TransformHaloVelocities(iDataInitial);
/*
		// Flux reconstruction update coefficient
		double dUpdateDeriv =
			  dDeltaT
			* dFluxDeriv1D[m_nHorizontalOrder-1]
			/ dElementDeltaA;
*/
		// Loop over edges of constant alpha
		for (int k = 0; k < pGrid->GetRElements(); k++) {
		for (int a = 0; a <= nElementCountA; a++) {

			int i = box.GetAInteriorBegin() + a * m_nHorizontalOrder;
			int j = box.GetBInteriorBegin();
			for (; j < box.GetBInteriorEnd(); j++) {

				double dUaL = dataInitialNode[UIx][k][i-1][j];
				double dUbL = dataInitialNode[VIx][k][i-1][j];
				double dHL  = dataInitialNode[HIx][k][i-1][j];

				double dUaR = dataInitialNode[UIx][k][i][j];
				double dUbR = dataInitialNode[VIx][k][i][j];
				double dHR  = dataInitialNode[HIx][k][i][j];

				// Calculate pointwise height flux
				double dPtJacobian;
				double dZs;
				if (a == nElementCountA) {
					dPtJacobian = dJacobian[k][i-1][j];
					dZs         = dTopography[i-1][j];
				} else {
					dPtJacobian = dJacobian[k][i][j];
					dZs         = dTopography[i][j];
				}

				double dHFL = (dHL - dZs) * dUaL * dPtJacobian;
				double dHFR = (dHR - dZs) * dUaR * dPtJacobian;

				double dHF = 0.5 * (dHFL + dHFR);

				// Nodal velocities and pressure
				double dUa = 0.5 * (dUaL + dUaR);
				double dUb = 0.5 * (dUbL + dUbR);

				double dPL = phys.GetG() * dHL;
				double dPR = phys.GetG() * dHR;
				double dP  = 0.5 * (dPL + dPR);

				// Upwinding
				double dA = fabs(dUa)
					+ sqrt(phys.GetG() * 0.5 * (dHL + dHR))
						/ phys.GetEarthRadius();

				double dUaF = 0.0;
				double dUbF = 0.0;

#ifdef PENALIZE_DISCONTINUITY
				dHF -= 0.5 * dA * dPtJacobian * (dHR - dHL);
				dUaF = - 0.5 * dA * (dUaR - dUaL);
				dUbF = - 0.5 * dA * (dUbR - dUbL);
#endif

#ifndef DIFFERENTIAL_FORM
				// Update the height field
				dataUpdateNode[HIx][k][i-1][j] +=
					- dDeltaT
					/ dGLLWeights1D[m_nHorizontalOrder-1]
					/ dElementDeltaA
					* dHF / dPtJacobian;

				dataUpdateNode[HIx][k][i][j] +=
					+ dDeltaT
					/ dGLLWeights1D[m_nHorizontalOrder-1]
					/ dElementDeltaA
					* dHF / dPtJacobian;
#endif

				for (int s = 0; s < m_nHorizontalOrder; s++) {
					double dUpdateDeriv =
						  dDeltaT
						* dFluxDeriv1D[m_nHorizontalOrder - 1 - s]
						/ dElementDeltaA;

					// Calculate modified derivatives in alpha
					if (a != 0) {
#ifdef DIFFERENTIAL_FORM 
						dataUpdateNode[HIx][k][i-1-s][j] +=
							- dUpdateDeriv * (dHF - dHFL)
							/ dJacobian[k][i-1-s][j];
#endif
#ifndef ADVECTION_ONLY
						dataUpdateNode[UIx][k][i-1-s][j] +=
							- dUpdateDeriv * (
								dUaL * (dUa - dUaL)
								+ dContraMetricA[k][i-1][j][0] * (dP - dPL))
							- dUpdateDeriv * dUaF;

						dataUpdateNode[VIx][k][i-1-s][j] +=
							- dUpdateDeriv * (
								dUaL * (dUb - dUbL)
								+ dContraMetricB[k][i-1][j][0] * (dP - dPL))
							- dUpdateDeriv * dUbF;
#endif
					}

					if (a != nElementCountA) {
#ifdef DIFFERENTIAL_FORM
						dataUpdateNode[HIx][k][i+s][j] +=
							+ dUpdateDeriv * (dHF - dHFR)
							/ dJacobian[k][i+s][j];
#endif
#ifndef ADVECTION_ONLY
						dataUpdateNode[UIx][k][i+s][j] +=
							+ dUpdateDeriv * (
								dUaR * (dUa - dUaR)
								+ dContraMetricA[k][i][j][0] * (dP - dPR))
							+ dUpdateDeriv * dUaF;

						dataUpdateNode[VIx][k][i+s][j] +=
							+ dUpdateDeriv * (
								dUaR * (dUb - dUbR)
								+ dContraMetricB[k][i][j][0] * (dP - dPR))
							+ dUpdateDeriv * dUbF;
#endif
					}
				}
			}
		}
		}

		// Loop over edges of constant beta
		for (int k = 0; k < pGrid->GetRElements(); k++) {
		for (int b = 0; b <= nElementCountB; b++) {

			int i = box.GetAInteriorBegin();
			int j = box.GetBInteriorBegin() + b * m_nHorizontalOrder;
			for (; i < box.GetAInteriorEnd(); i++) {

				double dUaL = dataInitialNode[UIx][k][i][j-1];
				double dUbL = dataInitialNode[VIx][k][i][j-1];
				double dHL  = dataInitialNode[HIx][k][i][j-1];

				double dUaR = dataInitialNode[UIx][k][i][j];
				double dUbR = dataInitialNode[VIx][k][i][j];
				double dHR  = dataInitialNode[HIx][k][i][j];

				// Calculate pointwise height flux
				double dPtJacobian;
				double dZs;
				if (b == nElementCountB) {
					dPtJacobian = dJacobian[k][i][j-1];
					dZs = dTopography[i][j-1];
				} else {
					dPtJacobian = dJacobian[k][i][j];
					dZs = dTopography[i][j];
				}

				double dHFL = (dHL - dZs) * dUbL * dPtJacobian;
				double dHFR = (dHR - dZs) * dUbR * dPtJacobian;

				double dHF = 0.5 * (dHFL + dHFR);

				// Nodal velocities and pressure
				double dUa = 0.5 * (dUaL + dUaR);
				double dUb = 0.5 * (dUbL + dUbR);

				double dPL = phys.GetG() * dHL;
				double dPR = phys.GetG() * dHR;
				double dP  = 0.5 * (dPL + dPR);

				// Upwinding
				double dA = fabs(dUb)
					+ sqrt(phys.GetG() * 0.5 * (dHL + dHR))
						/ phys.GetEarthRadius();

				double dUaF = 0.0;
				double dUbF = 0.0;

#ifdef PENALIZE_DISCONTINUITY
				dHF -= 0.5 * dA * dPtJacobian * (dHR - dHL);
				dUaF = - 0.5 * dA * (dUaR - dUaL);
				dUbF = - 0.5 * dA * (dUbR - dUbL);
#endif

#ifndef DIFFERENTIAL_FORM
				// Update the height field
				dataUpdateNode[HIx][k][i][j-1] +=
					- dDeltaT
					/ dGLLWeights1D[m_nHorizontalOrder-1]
					/ dElementDeltaA
					* dHF / dPtJacobian;

				dataUpdateNode[HIx][k][i][j] +=
					+ dDeltaT
					/ dGLLWeights1D[m_nHorizontalOrder-1]
					/ dElementDeltaA
					* dHF / dPtJacobian;
#endif

				for (int s = 0; s < m_nHorizontalOrder; s++) {
					double dUpdateDeriv =
						  dDeltaT
						* dFluxDeriv1D[m_nHorizontalOrder - 1 - s]
						/ dElementDeltaB;

					// Calculate modified derivatives in beta
					if (b != 0) {
#ifdef DIFFERENTIAL_FORM 
						dataUpdateNode[HIx][k][i][j-1-s] +=
							- dUpdateDeriv * (dHF - dHFL)
							/ dJacobian[k][i][j-1-s];
#endif
#ifndef ADVECTION_ONLY
						dataUpdateNode[UIx][k][i][j-1-s] +=
							- dUpdateDeriv * (
								dUbL * (dUa - dUaL)
							 	+ dContraMetricA[k][i][j-1][1] * (dP - dPL))
							- dUpdateDeriv * dUaF;

						dataUpdateNode[VIx][k][i][j-1-s] +=
							- dUpdateDeriv * (
								dUbL * (dUb - dUbL)
								+ dContraMetricB[k][i][j-1][1] * (dP - dPL))
							- dUpdateDeriv * dUbF;
#endif
					}

					if (b != nElementCountB) {
#ifdef DIFFERENTIAL_FORM
						dataUpdateNode[HIx][k][i][j+s] +=
							+ dUpdateDeriv * (dHF - dHFR)
							/ dJacobian[k][i][j+s];
#endif
#ifndef ADVECTION_ONLY
						dataUpdateNode[UIx][k][i][j+s] +=
							+ dUpdateDeriv * (
								dUbR * (dUa - dUaR)
								+ dContraMetricA[k][i][j][1] * (dP - dPR))
							+ dUpdateDeriv * dUaF;

						dataUpdateNode[VIx][k][i][j+s] +=
							+ dUpdateDeriv * (
								dUbR * (dUb - dUbR)
							 	+ dContraMetricB[k][i][j][1] * (dP - dPR))
							+ dUpdateDeriv * dUbF;
#endif
					}
				}
			}
		}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void HorizontalDynamicsDG::StepNonhydrostaticPrimitive(
	int iDataInitial,
	int iDataUpdate,
	const Time & time,
	double dDeltaT
) {
/*
	// Get a copy of the GLL grid
	GridGLL * pGrid = dynamic_cast<GridGLL*>(m_model.GetGrid());

	// Physical constants
	const PhysicalConstants & phys = m_model.GetPhysicalConstants();

	// Indices of EquationSet variables
	const int UIx = 0;
	const int VIx = 1;
	const int TIx = 2;
	const int WIx = 3;
	const int RIx = 4;

	//std::cout << "Inside the horizontal step! \n";

	// Perform local update
	for (int n = 0; n < pGrid->GetActivePatchCount(); n++) {
		GridPatchGLL * pPatch =
			dynamic_cast<GridPatchGLL*>(pGrid->GetActivePatch(n));

		const PatchBox & box = pPatch->GetPatchBox();

		const DataArray2D<double> & dJacobian2D =
			pPatch->GetJacobian2D();
		const DataArray3D<double> & dJacobian =
			pPatch->GetJacobian();
		const DataArray4D<double> & dContraMetricA =
			pPatch->GetContraMetricA();
		const DataArray4D<double> & dContraMetricB =
			pPatch->GetContraMetricB();
		const DataArray4D<double> & dContraMetricXi =
			pPatch->GetContraMetricXi();
		const DataArray3D<double> & dChristoffelA =
			pPatch->GetChristoffelA();
		const DataArray3D<double> & dChristoffelB =
			pPatch->GetChristoffelB();
		const DataArray4D<double> & dChristoffelXi =
			pPatch->GetChristoffelXi();
		const DataArray2D<double> & dCoriolisF =
			pPatch->GetCoriolisF();
		const DataArray2D<double> & dTopography =
			pPatch->GetTopography();

		const double dZtop = pGrid->GetZtop();

		// Data
		DataArray4D<double> & dataInitialNode =
			pPatch->GetDataState(iDataInitial, DataLocation_Node);

		DataArray4D<double> & dataUpdateNode =
			pPatch->GetDataState(iDataUpdate, DataLocation_Node);

		DataArray4D<double> & dataInitialREdge =
			pPatch->GetDataState(iDataInitial, DataLocation_REdge);

		DataArray4D<double> & dataUpdateREdge =
			pPatch->GetDataState(iDataUpdate, DataLocation_REdge);

		// Perform interpolations as required due to vertical staggering
		if (pGrid->GetVarsAtLocation(DataLocation_REdge) != 0) {

			// Interpolate Theta to model levels
			if (pGrid->GetVarLocation(TIx) == DataLocation_REdge) {
				pPatch->InterpolateREdgeToNode(TIx, iDataInitial);
			}
			if ((pGrid->GetVarLocation(TIx) == DataLocation_Node) &&
				(pGrid->GetVarLocation(WIx) == DataLocation_REdge)
			) {
				pPatch->InterpolateNodeToREdge(TIx, iDataInitial);
			}

			// Interpolate U, V and Rho to model interfaces
			if ((pGrid->GetVarLocation(RIx) == DataLocation_Node) &&
				(pGrid->GetVarLocation(WIx) == DataLocation_REdge)
			) {
				pPatch->InterpolateNodeToREdge(UIx, iDataInitial);
				pPatch->InterpolateNodeToREdge(VIx, iDataInitial);
				pPatch->InterpolateNodeToREdge(RIx, iDataInitial);
			}
		}

		// Element grid spacing and derivative coefficients
		const double dElementDeltaA = pPatch->GetElementDeltaA();
		const double dElementDeltaB = pPatch->GetElementDeltaB();

		const DataArray2D<double> & dDxBasis1D = pGrid->GetDxBasis1D();
		const DataArray2D<double> & dStiffness1D = pGrid->GetStiffness1D();

		// Time over grid spacing ratio
		double dCourantA = dDeltaT / dElementDeltaA;
		double dCourantB = dDeltaT / dElementDeltaB;

		// Get number of finite elements in each coordinate direction
		int nElementCountA = pPatch->GetElementCountA();
		int nElementCountB = pPatch->GetElementCountB();

		// Pressure data
		DataArray3D<double> & dataPressure = pPatch->GetDataPressure();
		DataArray3D<double> & dataDxPressure = pPatch->GetDataDxPressure();

#pragma message "This can be optimized at finite element edges"
		// Loop over all nodes and compute pressure
		for (int i = box.GetAInteriorBegin(); i < box.GetAInteriorEnd(); i++) {
		for (int j = box.GetBInteriorBegin(); j < box.GetBInteriorEnd(); j++) {
			for (int k = 0; k < pGrid->GetRElements(); k++) {
				m_dColumnPressure[k] =
					phys.ExnerPressureFromRhoTheta(
						  dataInitialNode[RIx][k][i][j]
						* dataInitialNode[TIx][k][i][j]);

				dataPressure[k][i][j] = m_dColumnPressure[k];
			}

			// Differentiate pressures
			pGrid->DifferentiateNodeToNode(
				m_dColumnPressure,
				m_dColumnDxPressure);

			for (int k = 0; k < pGrid->GetRElements(); k++) {
				dataDxPressure[k][i][j] =
					m_dColumnDxPressure[k];
			}
		}
		}

		// Loop over all nodes
		for (int k = 0; k < pGrid->GetRElements(); k++) {
		for (int a = 0; a < nElementCountA; a++) {
		for (int b = 0; b < nElementCountB; b++) {

			// Pointwise fluxes and pressure within spectral element
			for (int i = 0; i < m_nHorizontalOrder; i++) {
			for (int j = 0; j < m_nHorizontalOrder; j++) {

				int iA = a * m_nHorizontalOrder + i + box.GetHaloElements();
				int iB = b * m_nHorizontalOrder + j + box.GetHaloElements();

				// Density flux
				m_dAlphaFlux[i][j] =
					dJacobian[k][iA][iB]
					* dataInitialNode[RIx][k][iA][iB]
					* dataInitialNode[UIx][k][iA][iB];

				m_dBetaFlux[i][j] =
					dJacobian[k][iA][iB]
					* dataInitialNode[RIx][k][iA][iB]
					* dataInitialNode[VIx][k][iA][iB];

				// Pointwise pressure
				m_dPressure[i][j] = dataPressure[k][iA][iB];
			}
			}

			// Pointwise update of quantities on model levels
			for (int i = 0; i < m_nHorizontalOrder; i++) {
			for (int j = 0; j < m_nHorizontalOrder; j++) {

				int iA = a * m_nHorizontalOrder + i + box.GetHaloElements();
				int iB = b * m_nHorizontalOrder + j + box.GetHaloElements();

				int iElementA = a * m_nHorizontalOrder + box.GetHaloElements();
				int iElementB = b * m_nHorizontalOrder + box.GetHaloElements();

				// Derivatives of the velocity field
				double dDaUa = 0.0;
				double dDaUb = 0.0;
				double dDbUa = 0.0;
				double dDbUb = 0.0;

				// Derivatives of the pressure field
				double dDaP = 0.0;
				double dDbP = 0.0;
				double dDxP = 0.0;

				// Aliases for alpha and beta velocities
				const double dUa = dataInitialNode[UIx][k][iA][iB];
				const double dUb = dataInitialNode[VIx][k][iA][iB];

				// Calculate derivatives in the alpha direction
				double dLocalUpdateRhoA = 0.0;
				for (int s = 0; s < m_nHorizontalOrder; s++) {
#ifdef DIFFERENTIAL_FORM
					// Update density: Differential formulation
					dLocalUpdateRhoA -=
						m_dAlphaFlux[s][j]
						* dDxBasis1D[s][i];
#else
					// Update density: Variational formulation
					dLocalUpdateRhoA +=
						m_dAlphaFlux[s][j]
						* dStiffness1D[i][s];
#endif
					// Derivative of alpha velocity with respect to alpha
					dDaUa +=
						dataInitialNode[UIx][k][iElementA+s][iB]
						* dDxBasis1D[s][i];

					// Derivative of beta velocity with respect to alpha
					dDaUb +=
						dataInitialNode[VIx][k][iElementA+s][iB]
						* dDxBasis1D[s][i];

					// Derivative of pressure with respect to alpha
					dDaP +=
						m_dPressure[s][j]
						* dDxBasis1D[s][i];
				}

				// Calculate derivatives in the beta direction
				double dLocalUpdateRhoB = 0.0;
				for (int s = 0; s < m_nHorizontalOrder; s++) {
#ifdef DIFFERENTIAL_FORM
					// Update density: Differential formulation
					dLocalUpdateRhoB -=
						m_dBetaFlux[i][s]
						* dDxBasis1D[s][j];
#else
					// Update density: Variational formulation
					dLocalUpdateRhoB +=
						m_dBetaFlux[i][s]
						* dStiffness1D[j][s];
#endif
					// Derivative of alpha velocity with respect to beta
					dDbUa +=
						dataInitialNode[UIx][k][iA][iElementB+s]
						* dDxBasis1D[s][j];

					// Derivative of beta velocity with respect to beta
					dDbUb +=
						dataInitialNode[VIx][k][iA][iElementB+s]
						* dDxBasis1D[s][j];

					// Derivative of pressure with respect to beta
					dDbP +=
						m_dPressure[i][s]
						* dDxBasis1D[s][j];
				}

				// Scale derivatives
				dDaUa /= dElementDeltaA;
				dDaUb /= dElementDeltaA;
				dDaP  /= dElementDeltaA;

				dDbUa /= dElementDeltaB;
				dDbUb /= dElementDeltaB;
				dDbP  /= dElementDeltaB;

#pragma message "Reference state?"
				dDxP  = dataDxPressure[k][iA][iB];

				// Momentum advection terms
				double dLocalUpdateUa = 0.0;
				double dLocalUpdateUb = 0.0;

				dLocalUpdateUa -= dUa * dDaUa + dUb * dDbUa;

				dLocalUpdateUb -= dUa * dDaUb + dUb * dDbUb;

				// Curvature terms
				dLocalUpdateUa -= 
						+ dChristoffelA[iA][iB][0] * dUa * dUa
						+ dChristoffelA[iA][iB][1] * dUa * dUb
						+ dChristoffelA[iA][iB][2] * dUb * dUb;

				dLocalUpdateUb -=
						+ dChristoffelB[iA][iB][0] * dUa * dUa
						+ dChristoffelB[iA][iB][1] * dUa * dUb
						+ dChristoffelB[iA][iB][2] * dUb * dUb;

				// Pressure derivatives
				dLocalUpdateUa -=
						( dContraMetricA[k][iA][iB][0] * dDaP
						+ dContraMetricA[k][iA][iB][1] * dDbP
						+ dContraMetricA[k][iA][iB][2] * dDxP)
							* dataInitialNode[TIx][k][iA][iB];

				dLocalUpdateUb -=
						( dContraMetricB[k][iA][iB][0] * dDaP
						+ dContraMetricB[k][iA][iB][1] * dDbP
						+ dContraMetricB[k][iA][iB][2] * dDxP)
							* dataInitialNode[TIx][k][iA][iB];

				// Coriolis forces
				dLocalUpdateUa -=
					dCoriolisF[iA][iB]
					* dJacobian2D[iA][iB]
					* ( + dContraMetricA[k][iA][iB][1] * dUa
						- dContraMetricA[k][iA][iB][0] * dUb);

				dLocalUpdateUb -=
					dCoriolisF[iA][iB]
					* dJacobian2D[iA][iB]
					* ( + dContraMetricB[k][iA][iB][1] * dUa
						- dContraMetricB[k][iA][iB][0] * dUb);

 				// Apply update to horizontal velocity on model levels
				dataUpdateNode[UIx][k][iA][iB] += dDeltaT * dLocalUpdateUa;
				dataUpdateNode[VIx][k][iA][iB] += dDeltaT * dLocalUpdateUb;

				// Update density on model levels
				dataUpdateNode[RIx][k][iA][iB] +=
					(dCourantA * dLocalUpdateRhoA
						+ dCourantB * dLocalUpdateRhoB)
					/ dJacobian[k][iA][iB];

				// Update the vertical velocity (on model levels)
				if (pGrid->GetVarLocation(WIx) == DataLocation_Node) {
					const double dUx = dataInitialNode[WIx][k][iA][iB];

					double dDaUx = 0.0;
					double dDbUx = 0.0;

					for (int s = 0; s < m_nHorizontalOrder; s++) {
						// Derivative of xi velocity with respect to alpha
						dDaUx +=
							dataInitialNode[WIx][k][iElementA+s][iB]
							* dDxBasis1D[s][i];

						// Derivative of xi velocity with respect to beta
						dDbUx +=
							dataInitialNode[WIx][k][iA][iElementB+s]
							* dDxBasis1D[s][j];
					}

					// Scale derivatives
					dDaUx /= dElementDeltaA;
					dDbUx /= dElementDeltaB;

					// Update vertical velocity
					dataUpdateNode[WIx][k][iA][iB] -=
						dDeltaT * (dUa * dDaUx + dUb * dDbUx);
				}

				// Update the potential temperature (on model levels)
				if (pGrid->GetVarLocation(TIx) == DataLocation_Node) {
					double dDaTheta = 0.0;
					double dDbTheta = 0.0;

					for (int s = 0; s < m_nHorizontalOrder; s++) {
						// Derivative of theta with respect to alpha
						dDaTheta +=
							dataInitialNode[TIx][k][iElementA+s][iB]
							* dDxBasis1D[s][i];

						// Derivative of theta with respect to beta
						dDbTheta +=
							dataInitialNode[TIx][k][iA][iElementB+s]
							* dDxBasis1D[s][j];
					}

					// Scale derivatives
					dDaTheta /= dElementDeltaA;
					dDbTheta /= dElementDeltaB;

					// Update potential temperature
					dataUpdateNode[TIx][k][iA][iB] -=
						dDeltaT * (dUa * dDaTheta + dUb * dDbTheta);
				}

				// Update tracers
			}
			}
		}
		}
		}

		// Update quantities on model interfaces
		for (int a = 0; a < nElementCountA; a++) {
		for (int b = 0; b < nElementCountB; b++) {

			// Pointwise update of horizontal velocities
			for (int i = 0; i < m_nHorizontalOrder; i++) {
			for (int j = 0; j < m_nHorizontalOrder; j++) {

				int iA = a * m_nHorizontalOrder + i + box.GetHaloElements();
				int iB = b * m_nHorizontalOrder + j + box.GetHaloElements();

				int iElementA = a * m_nHorizontalOrder + box.GetHaloElements();
				int iElementB = b * m_nHorizontalOrder + box.GetHaloElements();

				// Update the vertical velocity (on model interfaces)
				for (int k = 1; k < pGrid->GetRElements(); k++) {
					if (pGrid->GetVarLocation(WIx) != DataLocation_REdge) {
						break;
					}

					double dDaUx = 0.0;
					double dDbUx = 0.0;

					double dUaREdge = dataInitialREdge[UIx][k][iA][iB];
					double dUbREdge = dataInitialREdge[VIx][k][iA][iB];
					double dUxREdge = dataInitialREdge[WIx][k][iA][iB];

					for (int s = 0; s < m_nHorizontalOrder; s++) {
						// Derivative of xi velocity with respect to alpha
						dDaUx +=
							dataInitialREdge[WIx][k][iElementA+s][iB]
							* dDxBasis1D[s][i];

						// Derivative of xi velocity with respect to beta
						dDbUx +=
							dataInitialREdge[WIx][k][iA][iElementB+s]
							* dDxBasis1D[s][j];
					}

					// Scale derivatives
					dDaUx /= dElementDeltaA;
					dDbUx /= dElementDeltaB;

					// Update vertical velocity
					dataUpdateREdge[WIx][k][iA][iB] -=
						dDeltaT * (dUaREdge * dDaUx + dUbREdge * dDbUx);
				}

				// Update the potential temperature (on model interfaces)
				for (int k = 0; k <= pGrid->GetRElements(); k++) {
					if (pGrid->GetVarLocation(TIx) != DataLocation_REdge) {
						break;
					}

					double dDaTheta = 0.0;
					double dDbTheta = 0.0;

					const double dUaREdge = dataInitialREdge[UIx][k][iA][iB];
					const double dUbREdge = dataInitialREdge[VIx][k][iA][iB];

					for (int s = 0; s < m_nHorizontalOrder; s++) {
						// Derivative of theta with respect to alpha
						dDaTheta +=
							dataInitialREdge[TIx][k][iElementA+s][iB]
							* dDxBasis1D[s][i];

						// Derivative of theta with respect to beta
						dDbTheta +=
							dataInitialREdge[TIx][k][iA][iElementB+s]
							* dDxBasis1D[s][j];
					}

					// Scale derivatives
					dDaTheta /= dElementDeltaA;
					dDbTheta /= dElementDeltaB;

					// Update vertical velocity
					dataUpdateREdge[TIx][k][iA][iB] -=
						dDeltaT * (dUaREdge * dDaTheta + dUbREdge * dDbTheta);
				}
			}
			}
		}
		}
	}
*/
}

///////////////////////////////////////////////////////////////////////////////

void HorizontalDynamicsDG::StepExplicit(
	int iDataInitial,
	int iDataUpdate,
	const Time & time,
	double dDeltaT
) {
	if (iDataInitial == iDataUpdate) {
		_EXCEPTIONT(
			"HorizontalDynamics Step must have iDataInitial != iDataUpdate");
	}

	// Equation set
	const EquationSet & eqn = m_model.GetEquationSet();

	// Step the primitive nonhydrostatic equations
	if (eqn.GetType() == EquationSet::PrimitiveNonhydrostaticEquations) {
		StepNonhydrostaticPrimitive(iDataInitial, iDataUpdate, time, dDeltaT);

	// Step the shallow water equations
	} else if (eqn.GetType() == EquationSet::ShallowWaterEquations) {
		StepShallowWater(iDataInitial, iDataUpdate, time, dDeltaT);

		ElementFluxesShallowWater(iDataInitial, iDataUpdate, time, dDeltaT);

	// Invalid EquationSet
	} else {
		_EXCEPTIONT("Invalid EquationSet");
	}
}

///////////////////////////////////////////////////////////////////////////////

void HorizontalDynamicsDG::InitializeApplyHyperdiffusionToBoundary(
	int iDataInitial
) {
	// Get a copy of the GLL grid
	GridGLL * pGrid = dynamic_cast<GridGLL*>(m_model.GetGrid());

	// Perform exchange
	pGrid->Exchange(DataType_State, iDataInitial);

	// Transform halo velocities
	for (int n = 0; n < pGrid->GetActivePatchCount(); n++) {
		GridPatchGLL * pPatch =
			dynamic_cast<GridPatchGLL*>(pGrid->GetActivePatch(n));

		pPatch->TransformHaloVelocities(iDataInitial);
	}
}

///////////////////////////////////////////////////////////////////////////////

void HorizontalDynamicsDG::ApplyScalarHyperdiffusionToBoundary(
	int iDataState,
	int iDataUpdate,
	double dDeltaT,
	double dNu,
	bool fScaleNuLocally
) {
	// Get a copy of the GLL grid
	GridGLL * pGrid = dynamic_cast<GridGLL*>(m_model.GetGrid());

	// Number of components
	int nComponents = m_model.GetEquationSet().GetComponents();

	// Derivatives of the flux reconstruction function
	const DataArray1D<double> & dFluxDeriv1D = pGrid->GetFluxDeriv1D();

	// Loop over all patches
	for (int n = 0; n < pGrid->GetActivePatchCount(); n++) {
		GridPatchGLL * pPatch =
			dynamic_cast<GridPatchGLL*>(pGrid->GetActivePatch(n));

		const PatchBox & box = pPatch->GetPatchBox();

		const DataArray4D<double> & dContraMetricA =
			pPatch->GetContraMetricA();
		const DataArray4D<double> & dContraMetricB =
			pPatch->GetContraMetricB();

		// Connectivity for patch
		Connectivity & connect = pPatch->GetConnectivity();

		// Data
		DataArray4D<double> & dataStateNode =
			pPatch->GetDataState(iDataState, DataLocation_Node);
		DataArray4D<double> & dataUpdateNode =
			pPatch->GetDataState(iDataUpdate, DataLocation_Node);

		DataArray4D<double> & dataStateREdge =
			pPatch->GetDataState(iDataState, DataLocation_REdge);
		DataArray4D<double> & dataUpdateREdge =
			pPatch->GetDataState(iDataUpdate, DataLocation_REdge);

		// Element grid spacing and derivative coefficients
		const double dElementDeltaA = pPatch->GetElementDeltaA();
		const double dElementDeltaB = pPatch->GetElementDeltaB();

		const DataArray2D<double> & dDxBasis1D = pGrid->GetDxBasis1D();
		const DataArray2D<double> & dStiffness1D = pGrid->GetStiffness1D();
		const DataArray1D<double> & dGLLWeights1D = pGrid->GetGLLWeights1D();

		// Flux reconstruction update coefficient
		double dUpdateDerivA =
			  dFluxDeriv1D[m_nHorizontalOrder-1] / dElementDeltaA;
		double dUpdateDerivB =
			  dFluxDeriv1D[m_nHorizontalOrder-1] / dElementDeltaB;

		// Get number of finite elements in each coordinate direction
		int nElementCountA = pPatch->GetElementCountA();
		int nElementCountB = pPatch->GetElementCountB();

		// Compute new hyperviscosity coefficient
		double dLocalNu  = dNu;

		if (fScaleNuLocally) {
			double dReferenceLength = pGrid->GetReferenceLength();
			dLocalNu *=
				pow(dElementDeltaA / dReferenceLength, ParamHypervisScaling);
		}

		// Loop over all scalar components
		for (int c = 2; c < nComponents; c++) {

			int nElementCountR;

			DataArray4D<double> * pDataState;
			DataArray4D<double> * pDataUpdate;
			if (pGrid->GetVarLocation(c) == DataLocation_Node) {
				pDataState = &dataStateNode;
				pDataUpdate = &dataUpdateNode;
				nElementCountR = pGrid->GetRElements();

			} else if (pGrid->GetVarLocation(c) == DataLocation_REdge) {
				pDataState = &dataStateREdge;
				pDataUpdate = &dataUpdateREdge;
				nElementCountR = pGrid->GetRElements()+1;

			} else {
				_EXCEPTIONT("UNIMPLEMENTED");
			}

			// Loop over perimeter of all elements
			for (int k = 0; k < nElementCountR; k++) {
			for (int a = 0; a < nElementCountA; a++) {
			for (int b = 0; b < nElementCountB; b++) {

				// Pointwise update of scalar quantities
				for (int i = 0; i < m_nHorizontalOrder; i++) {
				for (int j = 0; j < m_nHorizontalOrder; j++) {

					// Only perform computation on perimeter elements
					if ((i != 0) && (i != m_nHorizontalOrder-1) &&
						(j != 0) && (j != m_nHorizontalOrder-1)
					) {
						continue;
					}

					// Local indices
					int iA = a * m_nHorizontalOrder + i + box.GetHaloElements();
					int iB = b * m_nHorizontalOrder + j + box.GetHaloElements();

					int iElementA =
						a * m_nHorizontalOrder + box.GetHaloElements();
					int iElementB =
						b * m_nHorizontalOrder + box.GetHaloElements();

					// Calculate local derivatives
					double dDaPsi = 0.0;
					double dDbPsi = 0.0;

					for (int s = 0; s < m_nHorizontalOrder; s++) {
						// Derivative with respect to alpha
						dDaPsi +=
							(*pDataState)[c][k][iElementA+s][iB]
							* dDxBasis1D[s][i];

						// Derivative with respect to beta
						dDbPsi +=
							(*pDataState)[c][k][iA][iElementB+s]
							* dDxBasis1D[s][j];
					}

					dDaPsi /= dElementDeltaA;
					dDbPsi /= dElementDeltaB;

					// Add contribution due to boundaries
					if (i == 0) {
						double dPsiL = (*pDataState)[c][k][iA-1][iB];
						double dPsiR = (*pDataState)[c][k][iA  ][iB];

						dDaPsi += 0.5 * dUpdateDerivA * (dPsiR - dPsiL);
					}
					if (i == m_nHorizontalOrder-1) {
						double dPsiL = (*pDataState)[c][k][iA  ][iB];
						double dPsiR = (*pDataState)[c][k][iA+1][iB];

						dDaPsi += 0.5 * dUpdateDerivA * (dPsiR - dPsiL);
					}
					if (j == 0) {
						double dPsiL = (*pDataState)[c][k][iA][iB-1];
						double dPsiR = (*pDataState)[c][k][iA][iB  ];

						dDbPsi += 0.5 * dUpdateDerivB * (dPsiR - dPsiL);
					}
					if (j == m_nHorizontalOrder-1) {
						double dPsiL = (*pDataState)[c][k][iA][iB  ];
						double dPsiR = (*pDataState)[c][k][iA][iB+1];

						dDbPsi += 0.5 * dUpdateDerivB * (dPsiR - dPsiL);
					}

					// Calculate contravariant derivative
					double dGradDaPsi;
					double dGradDbPsi;

					dGradDaPsi =
						  dContraMetricA[k][iA][iB][0] * dDaPsi
						+ dContraMetricA[k][iA][iB][1] * dDbPsi;

					dGradDbPsi =
						  dContraMetricB[k][iA][iB][0] * dDaPsi
						+ dContraMetricB[k][iA][iB][1] * dDbPsi;

					// Calculate update
					double dUpdateA =
						0.5 * dDeltaT * dLocalNu * dGradDaPsi
							/ (dGLLWeights1D[i] * dElementDeltaA);

					double dUpdateB =
						0.5 * dDeltaT * dLocalNu * dGradDbPsi
							/ (dGLLWeights1D[j] * dElementDeltaB);

					// Either set up communication with neighbor or apply
					// scalar hyperdiffusion along right edge.
					if (i == m_nHorizontalOrder-1) {
						if (a == nElementCountA-1) {
							connect.SetSendBuffer(
								Direction_Right, c, k, iB, dUpdateA);

						} else {
							(*pDataUpdate)[c][k][iA+1][iB] -= dUpdateA;
						}
						(*pDataUpdate)[c][k][iA][iB] += dUpdateA;
					}

					// Either set up communication with neighbor or apply
					// scalar hyperdiffusion along top edge.
					if (j == m_nHorizontalOrder-1) {
						if (b == nElementCountB-1) {
							connect.SetSendBuffer(
								Direction_Top, c, k, iA, dUpdateB);

						} else {
							(*pDataUpdate)[c][k][iA][iB+1] -= dUpdateB;
						}
						(*pDataUpdate)[c][k][iA][iB] += dUpdateB;
					}

					// Either set up communication with neighbor or apply
					// scalar hyperdiffusion along left edge.
					if (i == 0) {
						if (a == 0) {
							connect.SetSendBuffer(
								Direction_Left, c, k, iB, dUpdateA);
						} else {
							(*pDataUpdate)[c][k][iA-1][iB] += dUpdateA;
						}
						(*pDataUpdate)[c][k][iA][iB] -= dUpdateA;
					}

					// Either set up communication with neighbor or apply
					// scalar hyperdiffusion along bottom edge.
					if (j == 0) {
						if (b == 0) {
							connect.SetSendBuffer(
								Direction_Bottom, c, k, iA, dUpdateB);

						} else {
							(*pDataUpdate)[c][k][iA][iB-1] += dUpdateB;
						}
						(*pDataUpdate)[c][k][iA][iB] -= dUpdateB;
					}
				}
				}
			}
			}
			}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void HorizontalDynamicsDG::ApplyVectorHyperdiffusionToBoundary(
	int iDataState,
	int iDataUpdate,
	double dDeltaT,
	double dNuDiv,
	double dNuVort,
	bool fScaleNuLocally
) {
	// Variable indices
	const int UIx = 0;
	const int VIx = 1;

	// Get a copy of the GLL grid
	GridGLL * pGrid = dynamic_cast<GridGLL*>(m_model.GetGrid());

	// Number of components
	int nComponents = m_model.GetEquationSet().GetComponents();

	// Loop over all patches
	for (int n = 0; n < pGrid->GetActivePatchCount(); n++) {
		GridPatchGLL * pPatch =
			dynamic_cast<GridPatchGLL*>(pGrid->GetActivePatch(n));

		const PatchBox & box = pPatch->GetPatchBox();

		const DataArray2D<double> & dJacobian =
			pPatch->GetJacobian2D();
		const DataArray4D<double> & dContraMetricA =
			pPatch->GetContraMetricA();
		const DataArray4D<double> & dContraMetricB =
			pPatch->GetContraMetricB();

		// Connectivity for patch
		Connectivity & connect = pPatch->GetConnectivity();

		// Data
		DataArray4D<double> & dataStateNode =
			pPatch->GetDataState(iDataState, DataLocation_Node);
		DataArray4D<double> & dataUpdateNode =
			pPatch->GetDataState(iDataUpdate, DataLocation_Node);

		DataArray4D<double> & dataStateREdge =
			pPatch->GetDataState(iDataState, DataLocation_REdge);
		DataArray4D<double> & dataUpdateREdge =
			pPatch->GetDataState(iDataUpdate, DataLocation_REdge);

		// Element grid spacing and derivative coefficients
		const double dElementDeltaA = pPatch->GetElementDeltaA();
		const double dElementDeltaB = pPatch->GetElementDeltaB();

		const DataArray2D<double> & dDxBasis1D = pGrid->GetDxBasis1D();
		const DataArray2D<double> & dStiffness1D = pGrid->GetStiffness1D();
		const DataArray1D<double> & dGLLWeights1D = pGrid->GetGLLWeights1D();

		// Get number of finite elements in each coordinate direction
		int nElementCountA = pPatch->GetElementCountA();
		int nElementCountB = pPatch->GetElementCountB();

		// Get curl and divergence
		const DataArray3D<double> & dataCurl = pPatch->GetDataVorticity();
		const DataArray3D<double> & dataDiv  = pPatch->GetDataDivergence();

		// Compute new hyperviscosity coefficient
		double dLocalNuDiv  = dNuDiv;
		double dLocalNuVort = dNuVort;

		if (fScaleNuLocally) {
			double dReferenceLength = pGrid->GetReferenceLength();
			dLocalNuDiv =
				dLocalNuDiv
				* pow(dElementDeltaA / dReferenceLength, ParamHypervisScaling);
			dLocalNuVort =
				dLocalNuVort
				* pow(dElementDeltaA / dReferenceLength, ParamHypervisScaling);
		}

		// Pointers to data
		int nElementCountR;

		DataArray3D<double> dataStateU;
		DataArray3D<double> dataStateV;

		DataArray3D<double> dataUpdateU;
		DataArray3D<double> dataUpdateV;

		if (pGrid->GetVarLocation(UIx) == DataLocation_Node) {
			dataStateU.SetSize(
				dataStateNode.GetSize(1),
				dataStateNode.GetSize(2),
				dataStateNode.GetSize(3));

			dataStateV.SetSize(
				dataStateNode.GetSize(1),
				dataStateNode.GetSize(2),
				dataStateNode.GetSize(3));

			dataUpdateU.SetSize(
				dataStateNode.GetSize(1),
				dataStateNode.GetSize(2),
				dataStateNode.GetSize(3));

			dataUpdateV.SetSize(
				dataStateNode.GetSize(1),
				dataStateNode.GetSize(2),
				dataStateNode.GetSize(3));

			dataStateU.AttachTo(dataStateNode[UIx]);
			dataStateV.AttachTo(dataStateNode[VIx]);
			dataUpdateU.AttachTo(dataStateNode[UIx]);
			dataUpdateV.AttachTo(dataStateNode[VIx]);
			nElementCountR = dataStateNode.GetSize(1);

		} else if (pGrid->GetVarLocation(UIx) == DataLocation_REdge) {
			dataStateU.SetSize(
				dataStateREdge.GetSize(1),
				dataStateREdge.GetSize(2),
				dataStateREdge.GetSize(3));

			dataStateV.SetSize(
				dataStateREdge.GetSize(1),
				dataStateREdge.GetSize(2),
				dataStateREdge.GetSize(3));

			dataUpdateU.SetSize(
				dataStateREdge.GetSize(1),
				dataStateREdge.GetSize(2),
				dataStateREdge.GetSize(3));

			dataUpdateV.SetSize(
				dataStateREdge.GetSize(1),
				dataStateREdge.GetSize(2),
				dataStateREdge.GetSize(3));

			dataStateU.AttachTo(dataStateREdge[UIx]);
			dataStateV.AttachTo(dataStateREdge[VIx]);
			dataUpdateU.AttachTo(dataStateREdge[UIx]);
			dataUpdateV.AttachTo(dataStateREdge[VIx]);
			nElementCountR = dataStateREdge.GetSize(1);

		} else {
			_EXCEPTIONT("UNIMPLEMENTED");
		}

		// Loop over perimeter of all elements
		for (int k = 0; k < nElementCountR; k++) {
		for (int a = 0; a < nElementCountA; a++) {
		for (int b = 0; b < nElementCountB; b++) {

			// Pointwise update of scalar quantities
			for (int i = 0; i < m_nHorizontalOrder; i++) {
			for (int j = 0; j < m_nHorizontalOrder; j++) {

				// Only perform computation on perimeter elements
				if ((i != 0) && (i != m_nHorizontalOrder-1) &&
					(j != 0) && (j != m_nHorizontalOrder-1)
				) {
					continue;
				}

				// Local indices
				int iA = a * m_nHorizontalOrder + i + box.GetHaloElements();
				int iB = b * m_nHorizontalOrder + j + box.GetHaloElements();

				int iElementA =
					a * m_nHorizontalOrder + box.GetHaloElements();
				int iElementB =
					b * m_nHorizontalOrder + box.GetHaloElements();

				// Calculate update due to divergence of velocity field
				double dDivUpdateAA =
						0.5 * dDeltaT * dLocalNuDiv
						* dContraMetricA[k][iA][iB][0]
						* dataDiv[k][iA][iB]
						/ (dGLLWeights1D[i] * dElementDeltaA);

				double dDivUpdateAB = 
						0.5 * dDeltaT * dLocalNuDiv
						* dContraMetricA[k][iA][iB][1]
						* dataDiv[k][iA][iB]
						/ (dGLLWeights1D[j] * dElementDeltaB);

				double dDivUpdateBA =
						0.5 * dDeltaT * dLocalNuDiv
						* dContraMetricB[k][iA][iB][0]
						* dataDiv[k][iA][iB]
						/ (dGLLWeights1D[i] * dElementDeltaA);

				double dDivUpdateBB =
						0.5 * dDeltaT * dLocalNuDiv
						* dContraMetricB[k][iA][iB][1]
						* dataDiv[k][iA][iB]
						/ (dGLLWeights1D[j] * dElementDeltaB);

				// Calculate update due to curl of velocity field
				double dVortUpdateA =
						0.5 * dDeltaT * dLocalNuVort
						* dataCurl[k][iA][iB]
						/ dJacobian[iA][iB]
						/ (dGLLWeights1D[j] * dElementDeltaB);

				double dVortUpdateB =
						0.5 * dDeltaT * dLocalNuVort
						* dataCurl[k][iA][iB]
						/ dJacobian[iA][iB]
						/ (dGLLWeights1D[i] * dElementDeltaA);

				// Either set up communication with neighbor or apply
				// vector hyperdiffusion along right edge.
				if (i == m_nHorizontalOrder-1) {
					double dUpdateU = dDivUpdateAA;
					double dUpdateV = dDivUpdateBA + dVortUpdateB;

					if (a == nElementCountA-1) {
						connect.SetSendBuffer(
							Direction_Right, UIx, k, iB, dUpdateU);
						connect.SetSendBuffer(
							Direction_Right, VIx, k, iB, dUpdateV);

					} else {
						dataUpdateU[k][iA+1][iB] -= dUpdateU;
						dataUpdateV[k][iA+1][iB] -= dUpdateV;
					}
					dataUpdateU[k][iA][iB] += dUpdateU;
					dataUpdateV[k][iA][iB] += dUpdateV;
				}

				// Either set up communication with neighbor or apply
				// vector hyperdiffusion along top edge.
				if (j == m_nHorizontalOrder-1) {
					double dUpdateU = dDivUpdateAB - dVortUpdateA;
					double dUpdateV = dDivUpdateBB;

					if (b == nElementCountB-1) {
						connect.SetSendBuffer(
							Direction_Top, UIx, k, iA, dUpdateU);
						connect.SetSendBuffer(
							Direction_Top, VIx, k, iA, dUpdateV);

					} else {
						dataUpdateU[k][iA][iB+1] -= dUpdateU;
						dataUpdateV[k][iA][iB+1] -= dUpdateV;
					}
					dataUpdateU[k][iA][iB] += dUpdateU;
					dataUpdateV[k][iA][iB] += dUpdateV;
				}

				// Either set up communication with neighbor or apply
				// vector hyperdiffusion along left edge.
				if (i == 0) {
					double dUpdateU = dDivUpdateAA;
					double dUpdateV = dDivUpdateBA + dVortUpdateB;

					if (a == 0) {
						connect.SetSendBuffer(
							Direction_Left, UIx, k, iB, dUpdateU);
						connect.SetSendBuffer(
							Direction_Left, VIx, k, iB, dUpdateV);

					} else {
						dataUpdateU[k][iA-1][iB] += dUpdateU;
						dataUpdateV[k][iA-1][iB] += dUpdateV;
					}
					dataUpdateU[k][iA][iB] -= dUpdateU;
					dataUpdateV[k][iA][iB] -= dUpdateV;
				}

				// Either set up communication with neighbor or apply
				// vector hyperdiffusion along bottom edge.
				if (j == 0) {
					double dUpdateU = dDivUpdateAB - dVortUpdateA;
					double dUpdateV = dDivUpdateBB;

					if (b == 0) {
						connect.SetSendBuffer(
							Direction_Bottom, UIx, k, iA, dUpdateU);
						connect.SetSendBuffer(
							Direction_Bottom, VIx, k, iA, dUpdateV);

					} else {
						dataUpdateU[k][iA][iB-1] += dUpdateU;
						dataUpdateV[k][iA][iB-1] += dUpdateV;
					}
					dataUpdateU[k][iA][iB] -= dUpdateU;
					dataUpdateV[k][iA][iB] -= dUpdateV;
				}
			}
			}
		}
		}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void HorizontalDynamicsDG::FinalizeApplyHyperdiffusionToBoundary(
	int iDataState,
	int iDataUpdate
) {
	// Get a copy of the GLL grid
	GridGLL * pGrid = dynamic_cast<GridGLL*>(m_model.GetGrid());

	// Number of components
	int nComponents = m_model.GetEquationSet().GetComponents();

	// Perform a global exchange
	pGrid->ExchangeBuffersAndUnpack(DataType_State, iDataUpdate);

	// Loop over all patches
	for (int n = 0; n < pGrid->GetActivePatchCount(); n++) {
		GridPatchGLL * pPatch =
			dynamic_cast<GridPatchGLL*>(pGrid->GetActivePatch(n));

		const PatchBox & box = pPatch->GetPatchBox();

		// Post-process velocities received during exchange
		pPatch->TransformHaloVelocities(iDataUpdate);

		// Connectivity for patch
		Connectivity & connect = pPatch->GetConnectivity();

		// Data
		DataArray4D<double> & dataUpdateNode =
			pPatch->GetDataState(iDataUpdate, DataLocation_Node);
		DataArray4D<double> & dataUpdateREdge =
			pPatch->GetDataState(iDataUpdate, DataLocation_REdge);

		// Loop over all components
		for (int c = 0; c < nComponents; c++) {

			int nElementCountR;

			DataArray4D<double> * pDataState;
			DataArray4D<double> * pDataUpdate;
			if (pGrid->GetVarLocation(c) == DataLocation_Node) {
				pDataUpdate = &dataUpdateNode;
				nElementCountR = pGrid->GetRElements();

			} else if (pGrid->GetVarLocation(c) == DataLocation_REdge) {
				pDataUpdate = &dataUpdateREdge;
				nElementCountR = pGrid->GetRElements()+1;

			} else {
				_EXCEPTIONT("UNIMPLEMENTED");
			}

			// Apply hyperviscosity to right edge
			for (int k = 0; k < nElementCountR; k++) {
			for (int j = box.GetBInteriorBegin();
				j < box.GetBInteriorEnd(); j++
			) {
				int i = box.GetAInteriorEnd();
				double dUpdateA = (*pDataUpdate)[c][k][i][j];

				if (connect.IsCoordinateFlipped(Direction_Right, j)) {
					dUpdateA *= -1.0;
				}

				(*pDataUpdate)[c][k][i-1][j] += dUpdateA;
			}
			}

			// Apply hyperviscosity to top edge
			for (int k = 0; k < nElementCountR; k++) {
			for (int i = box.GetAInteriorBegin();
				i < box.GetAInteriorEnd(); i++
			) {
				int j = box.GetBInteriorEnd();
				double dUpdateB = (*pDataUpdate)[c][k][i][j];

				if (connect.IsCoordinateFlipped(Direction_Top, i)) {
					dUpdateB *= -1.0;
				}

				(*pDataUpdate)[c][k][i][j-1] += dUpdateB;
			}
			}

			// Apply hyperviscosity to left edge
			for (int k = 0; k < nElementCountR; k++) {
			for (int j = box.GetBInteriorBegin();
				j < box.GetBInteriorEnd(); j++
			) {
				int i = box.GetAInteriorBegin()-1;
				double dUpdateA = (*pDataUpdate)[c][k][i][j];

				if (connect.IsCoordinateFlipped(Direction_Left, j)) {
					dUpdateA *= -1.0;
				}

				(*pDataUpdate)[c][k][i+1][j] -= dUpdateA;
			}
			}

			// Apply hyperviscosity to bottom edge
			for (int k = 0; k < nElementCountR; k++) {
			for (int i = box.GetAInteriorBegin();
				i < box.GetAInteriorEnd(); i++
			) {
				int j = box.GetBInteriorBegin()-1;
				double dUpdateB = (*pDataUpdate)[c][k][i][j];

				if (connect.IsCoordinateFlipped(Direction_Bottom, i)) {
					dUpdateB *= -1.0;
				}

				(*pDataUpdate)[c][k][i][j+1] -= dUpdateB;
			}
			}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void HorizontalDynamicsDG::StepAfterSubCycle(
	int iDataInitial,
	int iDataUpdate,
	int iDataWorking,
	const Time & time,
	double dDeltaT
) {

	// Check indices
	if (iDataInitial == iDataWorking) {
		_EXCEPTIONT("Invalid indices "
			"-- initial and working data must be distinct");
	}
	if (iDataUpdate == iDataWorking) {
		_EXCEPTIONT("Invalid indices "
			"-- working and update data must be distinct");
	}

	// Get the GLL grid
	GridGLL * pGrid = dynamic_cast<GridGLL*>(m_model.GetGrid());

	// No hyperdiffusion
	if ((m_dNuScalar == 0.0) && (m_dNuDiv == 0.0) && (m_dNuVort == 0.0)) {

	// Apply hyperdiffusion
	} else {
		// Apply scalar and vector hyperdiffusion (first application)
		pGrid->ZeroData(iDataWorking, DataType_State);
		//pGrid->CopyData(iDataInitial, iDataWorking, DataType_State);

		InitializeApplyHyperdiffusionToBoundary(iDataInitial);

		ApplyScalarHyperdiffusion(
			iDataInitial, iDataWorking, 1.0, 1.0, false);
		ApplyVectorHyperdiffusion(
			iDataInitial, iDataWorking, 1.0, 1.0, 1.0, false);

		ApplyScalarHyperdiffusionToBoundary(
			iDataInitial, iDataWorking, 1.0, 1.0, false);
		ApplyVectorHyperdiffusionToBoundary(
			iDataInitial, iDataWorking, 1.0, 1.0, 1.0, false);
		FinalizeApplyHyperdiffusionToBoundary(
			iDataInitial, iDataWorking);

		// Apply scalar and vector hyperdiffusion (second application)
		pGrid->CopyData(iDataInitial, iDataUpdate, DataType_State);

		InitializeApplyHyperdiffusionToBoundary(iDataWorking);

		ApplyScalarHyperdiffusion(
			iDataWorking, iDataUpdate, -dDeltaT, m_dNuScalar, true);
		ApplyVectorHyperdiffusion(
			iDataWorking, iDataUpdate, -dDeltaT, m_dNuDiv, m_dNuVort, true);

		ApplyScalarHyperdiffusionToBoundary(
			iDataWorking, iDataUpdate, -dDeltaT, m_dNuScalar, true);
		ApplyVectorHyperdiffusionToBoundary(
			iDataWorking, iDataUpdate, -dDeltaT, m_dNuDiv, m_dNuVort, true);
		FinalizeApplyHyperdiffusionToBoundary(
			iDataWorking, iDataUpdate);
	}


#ifdef APPLY_RAYLEIGH_WITH_HYPERVIS
	// Apply Rayleigh damping
	if (pGrid->HasRayleighFriction()) {
		ApplyRayleighFriction(iDataUpdate, dDeltaT);
	}
#endif
}

///////////////////////////////////////////////////////////////////////////////

