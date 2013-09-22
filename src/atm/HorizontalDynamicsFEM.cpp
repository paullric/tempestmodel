///////////////////////////////////////////////////////////////////////////////
///
///	\file    HorizontalDynamicsFEM.cpp
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

#include "HorizontalDynamicsFEM.h"
#include "PhysicalConstants.h"
#include "Model.h"
#include "Grid.h"
#include "GaussLobattoQuadrature.h"
#include "FluxReconstructionFunction.h"

#include "GridGLL.h"
#include "GridPatchGLL.h"

//#define DIFFERENTIAL_FORM

///////////////////////////////////////////////////////////////////////////////

HorizontalDynamicsFEM::HorizontalDynamicsFEM(
	Model & model,
	int nHorizontalOrder,
	HorizontalDynamicsFEM::Type eHorizontalDynamicsType,
	bool fNoHyperdiffusion
) :
	HorizontalDynamics(model),
	m_nHorizontalOrder(nHorizontalOrder),
	m_eHorizontalDynamicsType(eHorizontalDynamicsType),
	m_fNoHyperdiffusion(fNoHyperdiffusion),
	m_dNuScalar(1.0e15),
	m_dNuDiv(1.0e15),
	m_dNuVort(1.0e15)
{
	// Initialize the alpha and beta fluxes
	m_dAlphaFlux.Initialize(
		m_nHorizontalOrder,
		m_nHorizontalOrder);

	m_dBetaFlux.Initialize(
		m_nHorizontalOrder,
		m_nHorizontalOrder);

	m_dPressure.Initialize(
		m_nHorizontalOrder,
		m_nHorizontalOrder);

	// Get quadrature points for Gauss-Lobatto quadrature
	DataVector<double> dG;
	DataVector<double> dW;

	GaussLobattoQuadrature::GetPoints(m_nHorizontalOrder, 0.0, 1.0, dG, dW);

	// Get the derivatives of the flux reconstruction function
	m_dFluxDeriv1D.Initialize(m_nHorizontalOrder);
	FluxReconstructionFunction::GetDerivatives(
		2, m_nHorizontalOrder, dG, m_dFluxDeriv1D);
}

///////////////////////////////////////////////////////////////////////////////

void HorizontalDynamicsFEM::StepShallowWater(
	int iDataInitial,
	int iDataUpdate,
	double dTime,
	double dDeltaT
) {

	// Get a copy of the GLL grid
	GridGLL * pGrid = dynamic_cast<GridGLL*>(m_model.GetGrid());

	// Physical constants
	const PhysicalConstants & phys = m_model.GetPhysicalConstants();

	// Get indices of variables to update
	const int UIx = pGrid->GetVarIndex(0);
	const int VIx = pGrid->GetVarIndex(1);
	const int HIx = pGrid->GetVarIndex(2);

	// Perform local update
	for (int n = 0; n < pGrid->GetActivePatchCount(); n++) {
		GridPatchGLL * pPatch =
			dynamic_cast<GridPatchGLL*>(pGrid->GetActivePatch(n));

		const PatchBox & box = pPatch->GetPatchBox();

		const DataMatrix3D<double> & dJacobian =
			pPatch->GetJacobian();
		const DataMatrix4D<double> & dContraMetricA =
			pPatch->GetContraMetricA();
		const DataMatrix4D<double> & dContraMetricB =
			pPatch->GetContraMetricB();
		const DataMatrix3D<double> & dChristoffelA =
			pPatch->GetChristoffelA();
		const DataMatrix3D<double> & dChristoffelB =
			pPatch->GetChristoffelB();
		const DataMatrix<double> & dLatitude =
			pPatch->GetLatitude();
		const DataMatrix<double> & dCoriolisF =
			pPatch->GetCoriolisF();

		// Data
		const GridData4D & dataInitialNode =
			pPatch->GetDataState(iDataInitial, DataLocation_Node);

		GridData4D & dataUpdateNode =
			pPatch->GetDataState(iDataUpdate, DataLocation_Node);

		// Element grid spacing and derivative coefficients
		const double dElementDeltaA = pPatch->GetElementDeltaA();
		const double dElementDeltaB = pPatch->GetElementDeltaB();

		const DataMatrix<double> & dDxBasis1D = pGrid->GetDxBasis1D();
		const DataMatrix<double> & dStiffness1D = pGrid->GetStiffness1D();

		// Time over grid spacing ratio
		double dCourantA = dDeltaT / dElementDeltaA;
		double dCourantB = dDeltaT / dElementDeltaB;

		// Check interior domain size
		int nAPatchInteriorWidth = pPatch->GetPatchBox().GetAInteriorWidth();
		int nBPatchInteriorWidth = pPatch->GetPatchBox().GetBInteriorWidth();

		int nAElements = nAPatchInteriorWidth / m_nHorizontalOrder;
		int nBElements = nBPatchInteriorWidth / m_nHorizontalOrder;

		if ((nAPatchInteriorWidth % m_nHorizontalOrder) != 0) {
			_EXCEPTIONT("Logic Error: Invalid PatchBox alpha spacing");
		}
		if ((nBPatchInteriorWidth % m_nHorizontalOrder) != 0) {
			_EXCEPTIONT("Logic Error: Invalid PatchBox beta spacing");
		}

		for (int k = 0; k < pGrid->GetRElements(); k++) {
		for (int a = 0; a < nAElements; a++) {
		for (int b = 0; b < nBElements; b++) {

			// Pointwise fluxes within spectral element
			for (int i = 0; i < m_nHorizontalOrder; i++) {
			for (int j = 0; j < m_nHorizontalOrder; j++) {

				int iA = a * m_nHorizontalOrder + i + box.GetHaloElements();
				int iB = b * m_nHorizontalOrder + j + box.GetHaloElements();

				// Density flux
				m_dAlphaFlux[i][j] =
					dJacobian[k][iA][iB]
					* dataInitialNode[HIx][k][iA][iB]
					* dataInitialNode[UIx][k][iA][iB];

				m_dBetaFlux[i][j] =
					dJacobian[k][iA][iB]
					* dataInitialNode[HIx][k][iA][iB]
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

				int iAElement = a * m_nHorizontalOrder + box.GetHaloElements();
				int iBElement = b * m_nHorizontalOrder + box.GetHaloElements();

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
						dataInitialNode[UIx][k][iAElement+s][iB]
						* dDxBasis1D[s][i];

					// Derivative of beta velocity with respect to alpha
					dDaUb +=
						dataInitialNode[VIx][k][iAElement+s][iB]
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
						dataInitialNode[UIx][k][iA][iBElement+s]
						* dDxBasis1D[s][j];

					// Derivative of beta velocity with respect to beta
					dDbUb +=
						dataInitialNode[VIx][k][iA][iBElement+s]
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
						+ dContraMetricA[k][iA][iB][0] * dDaP
						+ dContraMetricA[k][iA][iB][1] * dDbP;

				dLocalUpdateUb -=
						+ dContraMetricB[k][iA][iB][0] * dDaP
						+ dContraMetricB[k][iA][iB][1] * dDbP;

				// Coriolis forces
				dLocalUpdateUa -=
					dCoriolisF[iA][iB] * dJacobian[k][iA][iB] * (
						+ dContraMetricA[k][iA][iB][1] * dUa
						- dContraMetricA[k][iA][iB][0] * dUb);

				dLocalUpdateUb -=
					dCoriolisF[iA][iB] * dJacobian[k][iA][iB] * (
						+ dContraMetricB[k][iA][iB][1] * dUa
						- dContraMetricB[k][iA][iB][0] * dUb);

				// Apply update
				dataUpdateNode[UIx][k][iA][iB] += dDeltaT * dLocalUpdateUa;
				dataUpdateNode[VIx][k][iA][iB] += dDeltaT * dLocalUpdateUb;

				// Update free surface height
				dataUpdateNode[HIx][k][iA][iB] +=
					(dCourantA * dLocalUpdateHa + dCourantB * dLocalUpdateHb)
					/ dJacobian[k][iA][iB];
			}
			}
		}
		}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void HorizontalDynamicsFEM::ElementFluxesShallowWater(
	int iDataInitial,
	int iDataUpdate,
	double dTime,
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

		const DataMatrix3D<double> & dJacobian =
			pPatch->GetJacobian();
		const DataMatrix4D<double> & dContraMetricA =
			pPatch->GetContraMetricA();
		const DataMatrix4D<double> & dContraMetricB =
			pPatch->GetContraMetricB();

		// Data
		GridData4D & dataInitialNode =
			pPatch->GetDataState(iDataInitial, DataLocation_Node);

		GridData4D & dataUpdateNode =
			pPatch->GetDataState(iDataUpdate, DataLocation_Node);

		GridData4D & dataInitialREdge =
			pPatch->GetDataState(iDataInitial, DataLocation_REdge);

		GridData4D & dataUpdateREdge =
			pPatch->GetDataState(iDataUpdate, DataLocation_REdge);

		// Element grid spacing and derivative coefficients
		const double dElementDeltaA = pPatch->GetElementDeltaA();
		const double dElementDeltaB = pPatch->GetElementDeltaB();

		// Time over grid spacing ratio
		double dCourantA = dDeltaT / dElementDeltaA;
		double dCourantB = dDeltaT / dElementDeltaB;

		// Check interior domain size
		int nAPatchInteriorWidth = pPatch->GetPatchBox().GetAInteriorWidth();
		int nBPatchInteriorWidth = pPatch->GetPatchBox().GetBInteriorWidth();

		int nAElements = nAPatchInteriorWidth / m_nHorizontalOrder;
		int nBElements = nBPatchInteriorWidth / m_nHorizontalOrder;

		if ((nAPatchInteriorWidth % m_nHorizontalOrder) != 0) {
			_EXCEPTIONT("Logic Error: Invalid PatchBox alpha spacing");
		}
		if ((nBPatchInteriorWidth % m_nHorizontalOrder) != 0) {
			_EXCEPTIONT("Logic Error: Invalid PatchBox beta spacing");
		}

		// Post-process velocities received during exchange
		pPatch->TransformHaloVelocities(iDataInitial);

		// Flux reconstruction update coefficient
		double dUpdateDeriv =
			  dDeltaT
			* m_dFluxDeriv1D[m_nHorizontalOrder-1]
			/ dElementDeltaA;

		// Loop over edges of constant alpha
		for (int k = 0; k < pGrid->GetRElements(); k++) {
		for (int a = 0; a <= nAElements; a++) {

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
				double dHFL = dHL * dUaL;
				double dHFR = dHR * dUaR;

				double dHF = 0.5 * (dHFL + dHFR);

#ifdef DIFFERENTIAL_FORM 
				_EXCEPTIONT("Not implemented.");
#else
				dataUpdateNode[HIx][k][i-1][j] -=
					  dUpdateDeriv * dHF;

				dataUpdateNode[HIx][k][i][j] +=
					  dUpdateDeriv * dHF;
#endif

				// Nodal pressure
				double dUa = 0.5 * (dUaL + dUaR);
				double dUb = 0.5 * (dUbL + dUbR);

				double dPL = phys.GetG() * dHL;
				double dPR = phys.GetG() * dHR;
				double dP  = 0.5 * (dPL + dPR);

				// Calculate modified derivatives in alpha
				dataUpdateNode[UIx][k][i-1][j] -=
					dUpdateDeriv * dUaL * (dUa - dUaL);

				dataUpdateNode[UIx][k][i][j] +=
					dUpdateDeriv * dUaR * (dUa - dUaR);

				dataUpdateNode[UIx][k][i-1][j] -=
					dUpdateDeriv * dContraMetricA[k][i-1][j][0] * (dP - dPL);

				dataUpdateNode[UIx][k][i][j] +=
					dUpdateDeriv * dContraMetricA[k][i][j][0] * (dP - dPR);

				dataUpdateNode[VIx][k][i-1][j] -=
					dUpdateDeriv * dUaL * (dUb - dUbL);

				dataUpdateNode[VIx][k][i][j] +=
					dUpdateDeriv * dUaR * (dUb - dUbR);

				dataUpdateNode[VIx][k][i-1][j] -=
					dUpdateDeriv * dContraMetricB[k][i-1][j][0] * (dP - dPL);

				dataUpdateNode[VIx][k][i][j] +=
					dUpdateDeriv * dContraMetricB[k][i][j][0] * (dP - dPR);
			}
		}
		}

		// Loop over edges of constant beta
		for (int k = 0; k < pGrid->GetRElements(); k++) {
		for (int b = 0; b <= nBElements; b++) {

			int i = box.GetAInteriorBegin();
			int j = box.GetBInteriorBegin() + b * m_nHorizontalOrder;
			for (; i < box.GetBInteriorEnd(); i++) {

				double dUaL = dataInitialNode[UIx][k][i][j-1];
				double dUbL = dataInitialNode[VIx][k][i][j-1];
				double dHL  = dataInitialNode[HIx][k][i][j-1];

				double dUaR = dataInitialNode[UIx][k][i][j];
				double dUbR = dataInitialNode[VIx][k][i][j];
				double dHR  = dataInitialNode[HIx][k][i][j];

				// Calculate pointwise height flux
				double dHFL = dHL * dUbL;
				double dHFR = dHR * dUbR;

				double dHF = 0.5 * (dHFL + dHFR);

#ifdef DIFFERENTIAL_FORM 
				_EXCEPTIONT("Not implemented.");
#else
				dataUpdateNode[HIx][k][i][j-1] -=
					  dUpdateDeriv * dHF;

				dataUpdateNode[HIx][k][i][j] +=
					  dUpdateDeriv * dHF;
#endif

				// Nodal pressure
				double dUa = 0.5 * (dUaL + dUaR);
				double dUb = 0.5 * (dUbL + dUbR);

				double dPL = phys.GetG() * dHL;
				double dPR = phys.GetG() * dHR;
				double dP  = 0.5 * (dPL + dPR);

				// Calculate modified derivatives in beta
				dataUpdateNode[UIx][k][i][j-1] -=
					dUpdateDeriv * dUbL * (dUa - dUaL);

				dataUpdateNode[UIx][k][i][j] +=
					dUpdateDeriv * dUbR * (dUa - dUaR);

				dataUpdateNode[UIx][k][i][j-1] -=
					dUpdateDeriv * dContraMetricA[k][i][j-1][1] * (dP - dPL);

				dataUpdateNode[UIx][k][i][j] +=
					dUpdateDeriv * dContraMetricA[k][i][j][1] * (dP - dPR);

				dataUpdateNode[VIx][k][i][j-1] -=
					dUpdateDeriv * dUbL * (dUb - dUbL);

				dataUpdateNode[VIx][k][i][j] +=
					dUpdateDeriv * dUbR * (dUb - dUbR);

				dataUpdateNode[VIx][k][i][j-1] -=
					dUpdateDeriv * dContraMetricB[k][i][j-1][1] * (dP - dPL);

				dataUpdateNode[VIx][k][i][j] +=
					dUpdateDeriv * dContraMetricB[k][i][j][1] * (dP - dPR);
			}
		}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void HorizontalDynamicsFEM::StepNonhydrostaticPrimitive(
	int iDataInitial,
	int iDataUpdate,
	double dTime,
	double dDeltaT
) {
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

	// Perform local update
	for (int n = 0; n < pGrid->GetActivePatchCount(); n++) {
		GridPatchGLL * pPatch =
			dynamic_cast<GridPatchGLL*>(pGrid->GetActivePatch(n));

		const PatchBox & box = pPatch->GetPatchBox();

		const DataMatrix3D<double> & dJacobian =
			pPatch->GetJacobian();
		const DataMatrix4D<double> & dContraMetricA =
			pPatch->GetContraMetricA();
		const DataMatrix4D<double> & dContraMetricB =
			pPatch->GetContraMetricB();
		const DataMatrix4D<double> & dContraMetricXi =
			pPatch->GetContraMetricXi();
		const DataMatrix3D<double> & dChristoffelA =
			pPatch->GetChristoffelA();
		const DataMatrix3D<double> & dChristoffelB =
			pPatch->GetChristoffelB();
		const DataMatrix4D<double> & dChristoffelXi =
			pPatch->GetChristoffelXi();
		const DataMatrix<double> & dLatitude =
			pPatch->GetLatitude();
		const DataMatrix<double> & dCoriolisF =
			pPatch->GetCoriolisF();
		const DataMatrix<double> & dTopography =
			pPatch->GetTopography();

		const double dZtop = pGrid->GetZtop();

		// Data
		GridData4D & dataInitialNode =
			pPatch->GetDataState(iDataInitial, DataLocation_Node);

		GridData4D & dataUpdateNode =
			pPatch->GetDataState(iDataUpdate, DataLocation_Node);

		GridData4D & dataInitialREdge =
			pPatch->GetDataState(iDataInitial, DataLocation_REdge);

		GridData4D & dataUpdateREdge =
			pPatch->GetDataState(iDataUpdate, DataLocation_REdge);

		// Perform interpolations as required due to vertical staggering
		if (pGrid->GetVarsAtLocation(DataLocation_REdge) != 0) {

			// Interpolate U and V to model interfaces
			pPatch->InterpolateNodeToREdge(UIx, iDataInitial);
			pPatch->InterpolateNodeToREdge(VIx, iDataInitial);

			// Interpolate Theta to model levels
			if (pGrid->GetVarLocation(TIx) == DataLocation_REdge) {
				pPatch->InterpolateREdgeToNode(TIx, iDataInitial);
			}
		}

		// Element grid spacing and derivative coefficients
		const double dElementDeltaA = pPatch->GetElementDeltaA();
		const double dElementDeltaB = pPatch->GetElementDeltaB();

		const DataMatrix<double> & dDxBasis1D = pGrid->GetDxBasis1D();
		const DataMatrix<double> & dStiffness1D = pGrid->GetStiffness1D();

		// Time over grid spacing ratio
		double dCourantA = dDeltaT / dElementDeltaA;
		double dCourantB = dDeltaT / dElementDeltaB;

		// Calculate pointwise fluxes
		int nAPatchInteriorWidth = pPatch->GetPatchBox().GetAInteriorWidth();
		int nBPatchInteriorWidth = pPatch->GetPatchBox().GetBInteriorWidth();

		int nAElements = nAPatchInteriorWidth / m_nHorizontalOrder;
		int nBElements = nBPatchInteriorWidth / m_nHorizontalOrder;

		if ((nAPatchInteriorWidth % m_nHorizontalOrder) != 0) {
			_EXCEPTIONT("Logic Error: Invalid PatchBox alpha spacing");
		}
		if ((nBPatchInteriorWidth % m_nHorizontalOrder) != 0) {
			_EXCEPTIONT("Logic Error: Invalid PatchBox beta spacing");
		}

		for (int k = 0; k < pGrid->GetRElements(); k++) {
		for (int a = 0; a < nAElements; a++) {
		for (int b = 0; b < nBElements; b++) {

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
				m_dPressure[i][j] =
					phys.PressureFromRhoTheta(
						dataInitialNode[RIx][k][iA][iB]
						* dataInitialNode[TIx][k][iA][iB]);
			}
			}

			// Pointwise update of quantities on model levels
			for (int i = 0; i < m_nHorizontalOrder; i++) {
			for (int j = 0; j < m_nHorizontalOrder; j++) {

				int iA = a * m_nHorizontalOrder + i + box.GetHaloElements();
				int iB = b * m_nHorizontalOrder + j + box.GetHaloElements();

				int iAElement = a * m_nHorizontalOrder + box.GetHaloElements();
				int iBElement = b * m_nHorizontalOrder + box.GetHaloElements();

				// Derivatives of the velocity field
				double dDaUa = 0.0;
				double dDaUb = 0.0;
				double dDbUa = 0.0;
				double dDbUb = 0.0;

				// Derivatives of the pressure field
				double dDaP = 0.0;
				double dDbP = 0.0;

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
						dataInitialNode[UIx][k][iAElement+s][iB]
						* dDxBasis1D[s][i];

					// Derivative of beta velocity with respect to alpha
					dDaUb +=
						dataInitialNode[VIx][k][iAElement+s][iB]
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
						dataInitialNode[UIx][k][iA][iBElement+s]
						* dDxBasis1D[s][j];

					// Derivative of beta velocity with respect to beta
					dDbUb +=
						dataInitialNode[VIx][k][iA][iBElement+s]
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
#pragma message "What about metric terms affected by vertical pressure derivative?"
				dLocalUpdateUa -=
						( dContraMetricA[k][iA][iB][0] * dDaP
						+ dContraMetricA[k][iA][iB][1] * dDbP)
							/ dataInitialNode[RIx][k][iA][iB];

				dLocalUpdateUb -=
						( dContraMetricB[k][iA][iB][0] * dDaP
						+ dContraMetricB[k][iA][iB][1] * dDbP)
							/ dataInitialNode[RIx][k][iA][iB];

				// Coriolis forces
				double dDomainHeight = dZtop - dTopography[iA][iB];

				dLocalUpdateUa -=
					dCoriolisF[iA][iB]
					* dJacobian[k][iA][iB]
					/ dDomainHeight * (
						+ dContraMetricA[k][iA][iB][1] * dUa
						- dContraMetricA[k][iA][iB][0] * dUb);

				dLocalUpdateUb -=
					dCoriolisF[iA][iB]
					* dJacobian[k][iA][iB]
					/ dDomainHeight * (
						+ dContraMetricB[k][iA][iB][1] * dUa
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
							dataInitialNode[WIx][k][iAElement+s][iB]
							* dDxBasis1D[s][i];

						// Derivative of xi velocity with respect to beta
						dDbUx +=
							dataInitialNode[WIx][k][iA][iBElement+s]
							* dDxBasis1D[s][j];
					}

					// Scale derivatives
					dDaUx /= dElementDeltaA;
					dDbUx /= dElementDeltaA;

					// Update vertical velocity
					dataUpdateNode[WIx][k][iA][iB] -=
						dDeltaT * (dUa * dDaUx + dUb * dDbUx);

					// Curvature terms
					double dCurvatureXi =
						+ dChristoffelXi[k][iA][iB][0] * dUa * dUa
						+ dChristoffelXi[k][iA][iB][1] * dUa * dUb
						+ dChristoffelXi[k][iA][iB][2] * dUa * dUx
						+ dChristoffelXi[k][iA][iB][3] * dUb * dUb
						+ dChristoffelXi[k][iA][iB][4] * dUb * dUx
						+ dChristoffelXi[k][iA][iB][5] * dUx * dUx;

					dataUpdateNode[WIx][k][iA][iB] -=
						dDeltaT * dCurvatureXi;

				}

				// Update the potential temperature (on model levels)
				if (pGrid->GetVarLocation(TIx) == DataLocation_Node) {
					double dDaTheta = 0.0;
					double dDbTheta = 0.0;

					for (int s = 0; s < m_nHorizontalOrder; s++) {
						// Derivative of theta with respect to alpha
						dDaTheta +=
							dataInitialNode[TIx][k][iAElement+s][iB]
							* dDxBasis1D[s][i];

						// Derivative of theta with respect to beta
						dDbTheta +=
							dataInitialNode[TIx][k][iA][iBElement+s]
							* dDxBasis1D[s][j];
					}

					// Scale derivatives
					dDaTheta /= dElementDeltaA;
					dDbTheta /= dElementDeltaA;

					// Update vertical velocity
					dataUpdateNode[TIx][k][iA][iB] -=
						dDeltaT * (dUa * dDaTheta + dUb * dDbTheta);
				}
			}
			}
		}
		}
		}

		// Update quantities on model interfaces
		for (int k = 0; k <= pGrid->GetRElements(); k++) {
		for (int a = 0; a < nAElements; a++) {
		for (int b = 0; b < nBElements; b++) {

			// Pointwise update of horizontal velocities
			for (int i = 0; i < m_nHorizontalOrder; i++) {
			for (int j = 0; j < m_nHorizontalOrder; j++) {

				int iA = a * m_nHorizontalOrder + i + box.GetHaloElements();
				int iB = b * m_nHorizontalOrder + j + box.GetHaloElements();

				int iAElement = a * m_nHorizontalOrder + box.GetHaloElements();
				int iBElement = b * m_nHorizontalOrder + box.GetHaloElements();

				// Update the vertical velocity (on model interfaces)
				if (pGrid->GetVarLocation(WIx) == DataLocation_REdge) {
					double dDaUx = 0.0;
					double dDbUx = 0.0;

					const double dUaREdge = dataInitialREdge[UIx][k][iA][iB];
					const double dUbREdge = dataInitialREdge[VIx][k][iA][iB];
					const double dUxREdge = dataInitialREdge[WIx][k][iA][iB];

					for (int s = 0; s < m_nHorizontalOrder; s++) {
						// Derivative of xi velocity with respect to alpha
						dDaUx +=
							dataInitialREdge[WIx][k][iAElement+s][iB]
							* dDxBasis1D[s][i];

						// Derivative of xi velocity with respect to beta
						dDbUx +=
							dataInitialREdge[WIx][k][iA][iBElement+s]
							* dDxBasis1D[s][j];
					}

					// Scale derivatives
					dDaUx /= dElementDeltaA;
					dDbUx /= dElementDeltaA;

					// Update vertical velocity
					dataUpdateNode[WIx][k][iA][iB] -=
						dDeltaT * (dUaREdge * dDaUx + dUbREdge * dDbUx);

					// Curvature terms
					double dCurvatureXi =
						+ dChristoffelXi[k][iA][iB][0] * dUaREdge * dUaREdge
						+ dChristoffelXi[k][iA][iB][1] * dUaREdge * dUbREdge
						+ dChristoffelXi[k][iA][iB][2] * dUaREdge * dUxREdge
						+ dChristoffelXi[k][iA][iB][3] * dUbREdge * dUbREdge
						+ dChristoffelXi[k][iA][iB][4] * dUbREdge * dUxREdge
						+ dChristoffelXi[k][iA][iB][5] * dUxREdge * dUxREdge;

					dataUpdateNode[WIx][k][iA][iB] -=
						dDeltaT * dCurvatureXi;
				}

				// Update the potential temperature (on model interfaces)
				if (pGrid->GetVarLocation(TIx) == DataLocation_REdge) {
					double dDaTheta = 0.0;
					double dDbTheta = 0.0;

					const double dUaREdge = dataInitialREdge[UIx][k][iA][iB];
					const double dUbREdge = dataInitialREdge[VIx][k][iA][iB];

					for (int s = 0; s < m_nHorizontalOrder; s++) {
						// Derivative of theta with respect to alpha
						dDaTheta +=
							dataInitialREdge[TIx][k][iAElement+s][iB]
							* dDxBasis1D[s][i];

						// Derivative of theta with respect to beta
						dDbTheta +=
							dataInitialREdge[TIx][k][iA][iBElement+s]
							* dDxBasis1D[s][j];
					}

					// Scale derivatives
					dDaTheta /= dElementDeltaA;
					dDbTheta /= dElementDeltaA;

					// Update vertical velocity
					dataUpdateREdge[TIx][k][iA][iB] -=
						dDeltaT * (dUaREdge * dDaTheta + dUbREdge * dDbTheta);
				}
			}
			}
		}
		}
		}

	}
}

///////////////////////////////////////////////////////////////////////////////

void HorizontalDynamicsFEM::StepExplicit(
	int iDataInitial,
	int iDataUpdate,
	double dTime,
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
		StepNonhydrostaticPrimitive(iDataInitial, iDataUpdate, dTime, dDeltaT);

	// Step the shallow water equations
	} else if (eqn.GetType() == EquationSet::ShallowWaterEquations) {
		StepShallowWater(iDataInitial, iDataUpdate, dTime, dDeltaT);

		if (m_eHorizontalDynamicsType == DiscontinuousGalerkin) {
			ElementFluxesShallowWater(iDataInitial, iDataUpdate, dTime, dDeltaT);
		}

	// Invalid EquationSet
	} else {
		_EXCEPTIONT("Invalid EquationSet");
	}

	// Apply Direct Stiffness Summation (DSS) procedure
	if (m_eHorizontalDynamicsType == SpectralElement) {
		GridGLL * pGrid = dynamic_cast<GridGLL*>(m_model.GetGrid());
		pGrid->ApplyDSS(iDataUpdate);
	}
}

///////////////////////////////////////////////////////////////////////////////
// OUTPUT INSTRUCTIONS FOR DEBUGGING

/*
				if (((n == 0) && (iA == box.GetAInteriorEnd()-1) && (iB == box.GetBInteriorEnd()-1)) ||
					((n == 1) && (iA == box.GetAInteriorBegin()) && (iB == box.GetBInteriorEnd()-1)) ||
					((n == 4) && (iA == box.GetAInteriorEnd()-1) && (iB == box.GetBInteriorBegin()))
				) {
					printf(" DP %i: %1.10e %1.10e\n", n, dDaP, dDbP);
				}
*/
/*
				if ((n == 0) && (iA == 3) && (iB == 3)) {
					printf(" AB: %1.10e %1.10e\n",
						box.GetANode(iA), box.GetBNode(iB));
					printf("  H: %1.10e\n", dataInitialNode[0][k][iA][iB] / dJacobian[k][iA][iB]);
					printf("  U: %1.10e %1.10e\n", dUa, dUb);
					printf("DaU: %1.10e %1.10e\n", dDaUa, dDaUb);
					printf("DbU: %1.10e %1.10e\n", dDbUa, dDbUb);
					printf(" DP: %1.10e %1.10e\n", dDaP, dDbP);
					printf("GammaA: %1.10e %1.10e %1.10e\n",
						dChristoffel[iA][iB][0],
						dChristoffel[iA][iB][1],
						dChristoffel[iA][iB][2]);
					printf("GammaB: %1.10e %1.10e %1.10e\n",
						dChristoffel[iA][iB][3],
						dChristoffel[iA][iB][4],
						dChristoffel[iA][iB][5]);
					printf("Mome: %1.10e %1.10e\n",
						(dUa * dDaUa + dUb * dDbUa), 
						(dUa * dDaUb + dUb * dDbUb));
					printf("Curv: %1.10e %1.10e\n",
						( dChristoffel[iA][iB][0] * dUa * dUa
						+ dChristoffel[iA][iB][1] * dUa * dUb
						+ dChristoffel[iA][iB][2] * dUb * dUb),
						( dChristoffel[iA][iB][3] * dUa * dUa
						+ dChristoffel[iA][iB][4] * dUa * dUb
						+ dChristoffel[iA][iB][5] * dUb * dUb));
					printf("Pres: %1.10e %1.10e\n",
						( dContraMetric[iA][iB][0] * dDaP
						+ dContraMetric[iA][iB][1] * dDbP),
						( dContraMetric[iA][iB][1] * dDaP
						+ dContraMetric[iA][iB][2] * dDbP));
					printf("Cori: %1.10e %1.10e\n",
						dF * dJacobian[k][iA][iB] * (
						+ dContraMetric[iA][iB][1] * dUa
						- dContraMetric[iA][iB][0] * dUb),
						dF * dJacobian[k][iA][iB] * (
						+ dContraMetric[iA][iB][2] * dUa
						- dContraMetric[iA][iB][1] * dUb));

				}
*/

///////////////////////////////////////////////////////////////////////////////

void HorizontalDynamicsFEM::ApplyScalarHyperdiffusion(
	int iDataInitial,
	int iDataUpdate,
	double dDeltaT,
	double dNu,
	bool fScaleNuLocally
) {
	// Get a copy of the GLL grid
	GridGLL * pGrid = dynamic_cast<GridGLL*>(m_model.GetGrid());

	// Physical constants
	const PhysicalConstants & phys = m_model.GetPhysicalConstants();

	// Components of the gradient at each point
	DataMatrix<double> dJGradientA;
	dJGradientA.Initialize(m_nHorizontalOrder, m_nHorizontalOrder);
	DataMatrix<double> dJGradientB;
	dJGradientB.Initialize(m_nHorizontalOrder, m_nHorizontalOrder);

	// Perform local update
	for (int n = 0; n < pGrid->GetActivePatchCount(); n++) {
		GridPatchGLL * pPatch =
			dynamic_cast<GridPatchGLL*>(pGrid->GetActivePatch(n));

		const PatchBox & box = pPatch->GetPatchBox();

		const DataMatrix<double> & dJacobian =
			pPatch->GetJacobian2D();
		const DataMatrix3D<double> & dContraMetricA =
			pPatch->GetContraMetric2DA();
		const DataMatrix3D<double> & dContraMetricB =
			pPatch->GetContraMetric2DB();
		const DataMatrix<double> & dTopography =
			pPatch->GetTopography();

		const double dZtop = pGrid->GetZtop();

		// Grid data
		GridData4D & dataInitialNode =
			pPatch->GetDataState(iDataInitial, DataLocation_Node);

		GridData4D & dataInitialREdge =
			pPatch->GetDataState(iDataInitial, DataLocation_REdge);

		GridData4D & dataUpdateNode =
			pPatch->GetDataState(iDataUpdate, DataLocation_Node);

		GridData4D & dataUpdateREdge =
			pPatch->GetDataState(iDataUpdate, DataLocation_REdge);

		// Element grid spacing and derivative coefficients
		const double dElementDeltaA = pPatch->GetElementDeltaA();
		const double dElementDeltaB = pPatch->GetElementDeltaB();

		const DataMatrix<double> & dDxBasis1D = pGrid->GetDxBasis1D();
		const DataMatrix<double> & dStiffness1D = pGrid->GetStiffness1D();

		// Number of finite elements
		int nAElements =
			box.GetAInteriorWidth() / m_nHorizontalOrder;
		int nBElements =
			box.GetBInteriorWidth() / m_nHorizontalOrder;

		// Compute new hyperviscosity coefficient
		double dLocalNu  = dNu;

		if (fScaleNuLocally) {
			double dReferenceLength = pGrid->GetReferenceLength();
			dLocalNu *= pow(dElementDeltaA / dReferenceLength, 3.2);
		}

		// Loop over all components
		int nComponents = m_model.GetEquationSet().GetComponents();
		for (int c = 2; c < nComponents; c++) {

			int nRElements;

			double *** pDataInitial;
			double *** pDataUpdate;
			if (pGrid->GetVarLocation(c) == DataLocation_Node) {
				pDataInitial = dataInitialNode[c];
				pDataUpdate  = dataUpdateNode[c];
				nRElements = dataInitialNode.GetRElements();

			} else if (pGrid->GetVarLocation(c) == DataLocation_REdge) {
				pDataInitial = dataInitialREdge[c];
				pDataUpdate  = dataUpdateREdge[c];
				nRElements = dataInitialREdge.GetRElements();

			} else {
				_EXCEPTIONT("UNIMPLEMENTED");
			}

			// Loop over all finite elements
			for (int k = 0; k < nRElements; k++) {
			for (int a = 0; a < nAElements; a++) {
			for (int b = 0; b < nBElements; b++) {

				int iAElement = a * m_nHorizontalOrder + box.GetHaloElements();
				int iBElement = b * m_nHorizontalOrder + box.GetHaloElements();

				// Calculate the pointwise gradient of the scalar field
				for (int i = 0; i < m_nHorizontalOrder; i++) {
				for (int j = 0; j < m_nHorizontalOrder; j++) {
					int iA = iAElement + i;
					int iB = iBElement + j;

					double dDaPsi = 0.0;
					double dDbPsi = 0.0;
					for (int s = 0; s < m_nHorizontalOrder; s++) {
						dDaPsi +=
							pDataInitial[k][iAElement+s][iB]
							* dDxBasis1D[s][i];

						dDbPsi +=
							pDataInitial[k][iA][iBElement+s]
							* dDxBasis1D[s][j];
					}

					dDaPsi /= dElementDeltaA;
					dDbPsi /= dElementDeltaB;

#pragma message "There should probably be a better method than using DomainHeight"
					double dDomainHeight = dZtop - dTopography[iA][iB];

					dJGradientA[i][j] = dJacobian[iA][iB] * dDomainHeight * (
						+ dContraMetricA[iA][iB][0] * dDaPsi
						+ dContraMetricA[iA][iB][1] * dDbPsi);

					dJGradientB[i][j] = dJacobian[iA][iB] * dDomainHeight * (
						+ dContraMetricB[iA][iB][0] * dDaPsi
						+ dContraMetricB[iA][iB][1] * dDbPsi);
				}
				}

				// Pointwise updates
				for (int i = 0; i < m_nHorizontalOrder; i++) {
				for (int j = 0; j < m_nHorizontalOrder; j++) {
					int iA = iAElement + i;
					int iB = iBElement + j;

					// Compute integral term
					double dUpdateA = 0.0;
					double dUpdateB = 0.0;

					for (int s = 0; s < m_nHorizontalOrder; s++) {
						dUpdateA +=
							dJGradientA[s][j]
							* dStiffness1D[i][s];

						dUpdateB +=
							dJGradientB[i][s]
							* dStiffness1D[j][s];
					}

					dUpdateA /= dElementDeltaA;
					dUpdateB /= dElementDeltaB;

					// Apply update
					double dDomainHeight = dZtop - dTopography[iA][iB];

					double dInvJacobian =
						1.0 / (dJacobian[iA][iB] * dDomainHeight);

					pDataUpdate[k][iA][iB] +=
						dDeltaT * dInvJacobian * dLocalNu
							* (dUpdateA + dUpdateB);
				}
				}
			}
			}
			}
		}
/*
			// Pointwise fluxes within spectral element
			for (int i = 0; i < m_nHorizontalOrder; i++) {
			for (int j = 0; j < m_nHorizontalOrder; j++) {

				double dDaGradient = 0.0;
				double dDbGradient = 0.0;

				for (int s = 0; s < m_nHorizontalOrder; s++) {
					dDaGradient +=
						m_dGradient[0][s][j]
						* m_dDxBasis1D[s][i];

					dDbGradient +=
						m_dGradient[1][i][s]
						* m_dDxBasis1D[s][j];
				}

				dDaGradient /= dElementDeltaA;
				dDbGradient /= dElementDeltaB;

				//printf("%1.10e %1.10e\n", dDaGradient, dDbGradient);

				// Update this variable
				int iA = a * m_nHorizontalOrder + i + box.GetHaloElements();
				int iB = b * m_nHorizontalOrder + j + box.GetHaloElements();

				dataUpdate[iC][k][iA][iB] +=
					dCoeff * (dDaGradient + dDbGradient) / dJacobian[k][iA][iB];
			}
			}

		}
		}
		}
*/
	}
}

///////////////////////////////////////////////////////////////////////////////

void HorizontalDynamicsFEM::ApplyVectorHyperdiffusion(
	int iDataInitial,
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

	// Physical constants
	const PhysicalConstants & phys = m_model.GetPhysicalConstants();

	// Loop over all patches
	for (int n = 0; n < pGrid->GetActivePatchCount(); n++) {
		GridPatchGLL * pPatch =
			dynamic_cast<GridPatchGLL*>(pGrid->GetActivePatch(n));

		const PatchBox & box = pPatch->GetPatchBox();

		const DataMatrix<double> & dJacobian =
			pPatch->GetJacobian2D();
		const DataMatrix3D<double> & dContraMetricA =
			pPatch->GetContraMetric2DA();
		const DataMatrix3D<double> & dContraMetricB =
			pPatch->GetContraMetric2DB();

		GridData4D & dataInitial =
			pPatch->GetDataState(iDataInitial, DataLocation_Node);

		GridData4D & dataUpdate =
			pPatch->GetDataState(iDataUpdate, DataLocation_Node);

		// Element grid spacing and derivative coefficients
		const double dElementDeltaA = pPatch->GetElementDeltaA();
		const double dElementDeltaB = pPatch->GetElementDeltaB();

		const DataMatrix<double> & dDxBasis1D = pGrid->GetDxBasis1D();
		const DataMatrix<double> & dStiffness1D = pGrid->GetStiffness1D();

		// Compute curl and divergence of U on the grid
		GridData3D dataUa;
		GridData3D dataUb;

		dataInitial.GetAsGridData3D(UIx, dataUa);
		dataInitial.GetAsGridData3D(VIx, dataUb);

		pPatch->ComputeCurlAndDiv(dataUa, dataUb);

		// Get curl and divergence
		const GridData3D & dataCurl = pPatch->GetDataVorticity();
		const GridData3D & dataDiv  = pPatch->GetDataDivergence();

		// Compute new hyperviscosity coefficient
		double dLocalNuDiv  = dNuDiv;
		double dLocalNuVort = dNuVort;

		if (fScaleNuLocally) {
			double dReferenceLength = pGrid->GetReferenceLength();
			dLocalNuDiv =
				dLocalNuDiv  * pow(dElementDeltaA / dReferenceLength, 3.2);
			dLocalNuVort =
				dLocalNuVort * pow(dElementDeltaA / dReferenceLength, 3.2);
		}

		// Number of finite elements
		int nAElements = box.GetAInteriorWidth() / m_nHorizontalOrder;
		int nBElements = box.GetBInteriorWidth() / m_nHorizontalOrder;

		// Loop over all finite elements
		for (int k = 0; k < pGrid->GetRElements(); k++) {
		for (int a = 0; a < nAElements; a++) {
		for (int b = 0; b < nBElements; b++) {

			int iAElement = a * m_nHorizontalOrder + box.GetHaloElements();
			int iBElement = b * m_nHorizontalOrder + box.GetHaloElements();

			// Pointwise update of horizontal velocities
			for (int i = 0; i < m_nHorizontalOrder; i++) {
			for (int j = 0; j < m_nHorizontalOrder; j++) {

				int iA = iAElement + i;
				int iB = iBElement + j;

				// Compute hyperviscosity sums
				double dAlphaDivTermA = 0.0;
				double dAlphaDivTermB = 0.0;
				double dAlphaCurlTerm = 0.0;

				double dBetaDivTermA = 0.0;
				double dBetaDivTermB = 0.0;
				double dBetaCurlTerm = 0.0;

				for (int s = 0; s < m_nHorizontalOrder; s++) {
					double dAlphaDiv =
						dJacobian[iAElement+s][iB]
						* dStiffness1D[i][s]
						* dataDiv[k][iAElement+s][iB];

					double dBetaDiv =
						dJacobian[iA][iBElement+s]
						* dStiffness1D[j][s]
						* dataDiv[k][iA][iBElement+s];

					dAlphaDivTermA +=
						dAlphaDiv * dContraMetricA[iAElement+s][iB][0];
					dAlphaDivTermB +=
						dBetaDiv  * dContraMetricB[iA][iBElement+s][0];

					dAlphaCurlTerm +=
						dStiffness1D[j][s]
						* dataCurl[k][iA][iBElement+s];

					dBetaDivTermA +=
						dAlphaDiv * dContraMetricA[iAElement+s][iB][1];
					dBetaDivTermB +=
						dBetaDiv  * dContraMetricB[iA][iBElement+s][1];

					dBetaCurlTerm +=
						dStiffness1D[i][s]
						* dataCurl[k][iAElement+s][iB];
				}

				dAlphaDivTermA /= dElementDeltaA;
				dAlphaDivTermB /= dElementDeltaB;
				dAlphaCurlTerm /= dElementDeltaB;

				dBetaDivTermA /= dElementDeltaA;
				dBetaDivTermB /= dElementDeltaB;
				dBetaCurlTerm /= dElementDeltaA;

				// Apply update
				double dInvJacobian = 1.0 / dJacobian[iA][iB];

				dataUpdate[UIx][k][iA][iB] += dDeltaT * dInvJacobian * (
					+ dLocalNuDiv * (
						+ dAlphaDivTermA
						+ dAlphaDivTermB)
					- dLocalNuVort * dAlphaCurlTerm);

				dataUpdate[VIx][k][iA][iB] += dDeltaT * dInvJacobian * (
					+ dLocalNuDiv * (
						+ dBetaDivTermA
						+ dBetaDivTermB)
					+ dLocalNuVort * dBetaCurlTerm);

			}
			}
		}
		}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void HorizontalDynamicsFEM::StepAfterSubCycle(
	int iDataInitial,
	int iDataUpdate,
	int iDataWorking,
	double dTime,
	double dDeltaT
) {
	// Only proceed if hyperdiffusion is applied
	if (m_fNoHyperdiffusion) {
		return;
	}

	// Check indices
	if (iDataInitial == iDataWorking) {
		_EXCEPTIONT("Invalid indices "
			"-- initial and working data must be distinct");
	}
	if (iDataUpdate == iDataWorking) {
		_EXCEPTIONT("Invalid indices "
			"-- working and update data must be distinct");
	}

#pragma message "Altering the horizontal velocities should also modify the vertical velocities"

	// Apply Direct Stiffness Summation (DSS) procedure
	GridGLL * pGrid = dynamic_cast<GridGLL*>(m_model.GetGrid());

	// Apply vector Laplacian (first application)
	pGrid->ZeroData(iDataWorking, DataType_State);

	ApplyScalarHyperdiffusion(
		iDataInitial, iDataWorking, 1.0, 1.0, false);
	ApplyVectorHyperdiffusion(
		iDataInitial, iDataWorking, 1.0, 1.0, 1.0, false);

	pGrid->ApplyDSS(iDataWorking);

	// Apply vector Laplacian (second application)
	ApplyScalarHyperdiffusion(
		iDataWorking, iDataUpdate, -dDeltaT, m_dNuScalar, true);
	ApplyVectorHyperdiffusion(
		iDataWorking, iDataUpdate, -dDeltaT, m_dNuDiv, m_dNuVort, true);

	pGrid->ApplyDSS(iDataUpdate);
}

///////////////////////////////////////////////////////////////////////////////

