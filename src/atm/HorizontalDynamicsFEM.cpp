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

#ifdef DIFFERENTIAL_FORM
#pragma message "WARNING: DIFFERENTIAL_FORM will lose mass over topography"
#endif

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

}

///////////////////////////////////////////////////////////////////////////////

void HorizontalDynamicsFEM::Initialize() {

	int nRElements = m_model.GetGrid()->GetRElements();

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

	// Column storage
	m_dZeroColumn.Initialize(nRElements + 1);

	m_dColumnPressure.Initialize(nRElements);

	m_dColumnDaPressure.Initialize(nRElements);
	m_dColumnDbPressure.Initialize(nRElements);
	m_dColumnDxPressure.Initialize(nRElements);

	m_dColumnDaPressureREdge.Initialize(nRElements + 1);
	m_dColumnDbPressureREdge.Initialize(nRElements + 1);

	// Initialize buffers for derivatives of Jacobian
	m_dJGradientA.Initialize(
		m_nHorizontalOrder,
		m_nHorizontalOrder);

	m_dJGradientB.Initialize(
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
	const Time & time,
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

		// Get number of finite elements in each coordinate direction
		int nElementCountA = pPatch->GetElementCountA();
		int nElementCountB = pPatch->GetElementCountB();

		// Post-process velocities received during exchange
		pPatch->TransformHaloVelocities(iDataInitial);

		// Flux reconstruction update coefficient
		double dUpdateDeriv =
			  dDeltaT
			* m_dFluxDeriv1D[m_nHorizontalOrder-1]
			/ dElementDeltaA;

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
		for (int b = 0; b <= nElementCountB; b++) {

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
	const int TIx = 2;
	const int WIx = 3;
	const int RIx = 4;

	//std::cout << "Inside the horizontal step! \n";

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

		const DataMatrix<double> & dDxBasis1D = pGrid->GetDxBasis1D();
		const DataMatrix<double> & dStiffness1D = pGrid->GetStiffness1D();

		// Time over grid spacing ratio
		double dCourantA = dDeltaT / dElementDeltaA;
		double dCourantB = dDeltaT / dElementDeltaB;

		// Get number of finite elements in each coordinate direction
		int nElementCountA = pPatch->GetElementCountA();
		int nElementCountB = pPatch->GetElementCountB();

		// Pressure data
		GridData3D & dataPressure = pPatch->GetDataPressure();
		GridData3D & dataDaPressure = pPatch->GetDataDaPressure();
		GridData3D & dataDbPressure = pPatch->GetDataDbPressure();
		GridData3D & dataDxPressure = pPatch->GetDataDxPressure();

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

				// Store local horizontal pressure derivatives
				dataDaPressure[k][iA][iB] = dDaP;
				dataDbPressure[k][iA][iB] = dDbP;

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

				// Store pressure derivatives in column
				for (int k = 0; k < pGrid->GetRElements(); k++) {
					m_dColumnDaPressure[k] = dataDaPressure[k][iA][iB];
					m_dColumnDbPressure[k] = dataDbPressure[k][iA][iB];
				}

				pGrid->InterpolateNodeToREdge(
					m_dColumnDaPressure,
					m_dZeroColumn,
					m_dColumnDaPressureREdge,
					m_dZeroColumn);

				pGrid->InterpolateNodeToREdge(
					m_dColumnDbPressure,
					m_dZeroColumn,
					m_dColumnDbPressureREdge,
					m_dZeroColumn);

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
}

///////////////////////////////////////////////////////////////////////////////

void HorizontalDynamicsFEM::StepExplicit(
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

		if (m_eHorizontalDynamicsType == DiscontinuousGalerkin) {
			ElementFluxesShallowWater(iDataInitial, iDataUpdate, time, dDeltaT);
		}

	// Invalid EquationSet
	} else {
		_EXCEPTIONT("Invalid EquationSet");
	}

/*
	// Apply Direct Stiffness Summation (DSS) procedure
	if (m_eHorizontalDynamicsType == SpectralElement) {
		GridGLL * pGrid = dynamic_cast<GridGLL*>(m_model.GetGrid());
		pGrid->ApplyDSS(iDataUpdate);
	}
*/
}

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
		int nElementCountA = pPatch->GetElementCountA();
		int nElementCountB = pPatch->GetElementCountB();

		// Compute new hyperviscosity coefficient
		double dLocalNu  = dNu;

		if (fScaleNuLocally) {
			double dReferenceLength = pGrid->GetReferenceLength();
			dLocalNu *= pow(dElementDeltaA / dReferenceLength, 3.2);
		}

		// Loop over all components
		int nComponents = m_model.GetEquationSet().GetComponents();
		for (int c = 2; c < nComponents; c++) {

			int nElementCountR;

			double *** pDataInitial;
			double *** pDataUpdate;
			if (pGrid->GetVarLocation(c) == DataLocation_Node) {
				pDataInitial = dataInitialNode[c];
				pDataUpdate  = dataUpdateNode[c];
				nElementCountR = dataInitialNode.GetRElements();

			} else if (pGrid->GetVarLocation(c) == DataLocation_REdge) {
				pDataInitial = dataInitialREdge[c];
				pDataUpdate  = dataUpdateREdge[c];
				nElementCountR = dataInitialREdge.GetRElements();

			} else {
				_EXCEPTIONT("UNIMPLEMENTED");
			}

			// Loop over all finite elements
			for (int k = 0; k < nElementCountR; k++) {
			for (int a = 0; a < nElementCountA; a++) {
			for (int b = 0; b < nElementCountB; b++) {

				int iElementA = a * m_nHorizontalOrder + box.GetHaloElements();
				int iElementB = b * m_nHorizontalOrder + box.GetHaloElements();

				// Calculate the pointwise gradient of the scalar field
				for (int i = 0; i < m_nHorizontalOrder; i++) {
				for (int j = 0; j < m_nHorizontalOrder; j++) {
					int iA = iElementA + i;
					int iB = iElementB + j;

					double dDaPsi = 0.0;
					double dDbPsi = 0.0;
					for (int s = 0; s < m_nHorizontalOrder; s++) {
						dDaPsi +=
							pDataInitial[k][iElementA+s][iB]
							* dDxBasis1D[s][i];

						dDbPsi +=
							pDataInitial[k][iA][iElementB+s]
							* dDxBasis1D[s][j];
					}

					dDaPsi /= dElementDeltaA;
					dDbPsi /= dElementDeltaB;

#pragma message "There should probably be a better method than using DomainHeight"
					double dDomainHeight = dZtop - dTopography[iA][iB];

					m_dJGradientA[i][j] = dJacobian[iA][iB] * dDomainHeight * (
						+ dContraMetricA[iA][iB][0] * dDaPsi
						+ dContraMetricA[iA][iB][1] * dDbPsi);

					m_dJGradientB[i][j] = dJacobian[iA][iB] * dDomainHeight * (
						+ dContraMetricB[iA][iB][0] * dDaPsi
						+ dContraMetricB[iA][iB][1] * dDbPsi);
				}
				}

				// Pointwise updates
				for (int i = 0; i < m_nHorizontalOrder; i++) {
				for (int j = 0; j < m_nHorizontalOrder; j++) {
					int iA = iElementA + i;
					int iB = iElementB + j;

					// Compute integral term
					double dUpdateA = 0.0;
					double dUpdateB = 0.0;

					for (int s = 0; s < m_nHorizontalOrder; s++) {
						dUpdateA +=
							m_dJGradientA[s][j]
							* dStiffness1D[i][s];

						dUpdateB +=
							m_dJGradientB[i][s]
							* dStiffness1D[j][s];
					}

					dUpdateA /= dElementDeltaA;
					dUpdateB /= dElementDeltaB;

					// Apply update
					double dDomainHeight = dZtop - dTopography[iA][iB];

					double dInvJacobian =
						1.0 / (dJacobian[iA][iB] * dDomainHeight);

					pDataUpdate[k][iA][iB] -=
						dDeltaT * dInvJacobian * dLocalNu
							* (dUpdateA + dUpdateB);
				}
				}
			}
			}
			}
		}
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
		int nElementCountA = pPatch->GetElementCountA();
		int nElementCountB = pPatch->GetElementCountB();

		// Loop over all finite elements
		for (int k = 0; k < pGrid->GetRElements(); k++) {
		for (int a = 0; a < nElementCountA; a++) {
		for (int b = 0; b < nElementCountB; b++) {

			int iElementA = a * m_nHorizontalOrder + box.GetHaloElements();
			int iElementB = b * m_nHorizontalOrder + box.GetHaloElements();

			// Pointwise update of horizontal velocities
			for (int i = 0; i < m_nHorizontalOrder; i++) {
			for (int j = 0; j < m_nHorizontalOrder; j++) {

				int iA = iElementA + i;
				int iB = iElementB + j;

				// Compute hyperviscosity sums
				double dAlphaDivTermA = 0.0;
				double dAlphaDivTermB = 0.0;
				double dAlphaCurlTerm = 0.0;

				double dBetaDivTermA = 0.0;
				double dBetaDivTermB = 0.0;
				double dBetaCurlTerm = 0.0;

				for (int s = 0; s < m_nHorizontalOrder; s++) {
					double dAlphaDiv =
						dJacobian[iElementA+s][iB]
						* dStiffness1D[i][s]
						* dataDiv[k][iElementA+s][iB];

					double dBetaDiv =
						dJacobian[iA][iElementB+s]
						* dStiffness1D[j][s]
						* dataDiv[k][iA][iElementB+s];

					dAlphaDivTermA +=
						dAlphaDiv * dContraMetricA[iElementA+s][iB][0];
					dAlphaDivTermB +=
						dBetaDiv  * dContraMetricB[iA][iElementB+s][0];

					dAlphaCurlTerm +=
						dStiffness1D[j][s]
						* dataCurl[k][iA][iElementB+s];

					dBetaDivTermA +=
						dAlphaDiv * dContraMetricA[iElementA+s][iB][1];
					dBetaDivTermB +=
						dBetaDiv  * dContraMetricB[iA][iElementB+s][1];

					dBetaCurlTerm +=
						dStiffness1D[i][s]
						* dataCurl[k][iElementA+s][iB];
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

void HorizontalDynamicsFEM::ApplyScalarHyperdiffusionToBoundary(
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

	// Physical constants
	const PhysicalConstants & phys = m_model.GetPhysicalConstants();

	// Loop over all patches
	for (int n = 0; n < pGrid->GetActivePatchCount(); n++) {
		GridPatchGLL * pPatch =
			dynamic_cast<GridPatchGLL*>(pGrid->GetActivePatch(n));

		const PatchBox & box = pPatch->GetPatchBox();

		const DataMatrix4D<double> & dContraMetricA =
			pPatch->GetContraMetricA();
		const DataMatrix4D<double> & dContraMetricB =
			pPatch->GetContraMetricB();
		const DataMatrix<double> & dTopography =
			pPatch->GetTopography();

		const double dZtop = pGrid->GetZtop();

		// Connectivity for patch
		Connectivity & connect = pPatch->GetConnectivity();

		// Data
		GridData4D & dataStateNode =
			pPatch->GetDataState(iDataState, DataLocation_Node);
		GridData4D & dataUpdateNode =
			pPatch->GetDataState(iDataUpdate, DataLocation_Node);

		GridData4D & dataStateREdge =
			pPatch->GetDataState(iDataState, DataLocation_REdge);
		GridData4D & dataUpdateREdge =
			pPatch->GetDataState(iDataUpdate, DataLocation_REdge);

		// Element grid spacing and derivative coefficients
		const double dElementDeltaA = pPatch->GetElementDeltaA();
		const double dElementDeltaB = pPatch->GetElementDeltaB();

		const DataMatrix<double> & dDxBasis1D = pGrid->GetDxBasis1D();
		const DataMatrix<double> & dStiffness1D = pGrid->GetStiffness1D();
		const DataVector<double> & dGLLWeights1D = pGrid->GetGLLWeights1D();

		// Get number of finite elements in each coordinate direction
		int nElementCountA = pPatch->GetElementCountA();
		int nElementCountB = pPatch->GetElementCountB();

		// Compute new hyperviscosity coefficient
		double dLocalNu  = dNu;

		if (fScaleNuLocally) {
			double dReferenceLength = pGrid->GetReferenceLength();
			dLocalNu *= pow(dElementDeltaA / dReferenceLength, 3.2);
		}

		// Loop over all scalar components
		int nComponents = m_model.GetEquationSet().GetComponents();
		for (int c = 2; c < nComponents; c++) {

			int nElementCountR;

			double *** pDataState;
			double *** pDataUpdate;
			if (pGrid->GetVarLocation(c) == DataLocation_Node) {
				pDataState = dataStateNode[c];
				pDataUpdate = dataUpdateNode[c];
				nElementCountR = dataStateNode.GetRElements();

			} else if (pGrid->GetVarLocation(c) == DataLocation_REdge) {
				pDataState = dataStateREdge[c];
				pDataUpdate = dataUpdateREdge[c];
				nElementCountR = dataStateREdge.GetRElements();

			} else {
				_EXCEPTIONT("UNIMPLEMENTED");
			}

			// Loop over perimeter of all elements
			for (int k = 0; k < nElementCountR; k++) {
			for (int a = 0; a < nElementCountA; a++) {
			for (int b = 0; b < nElementCountB; b++) {

				// Pointwise update of horizontal velocities
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
							pDataState[k][iElementA+s][iB]
							* dDxBasis1D[s][i];

						// Derivative with respect to beta
						dDbPsi +=
							pDataState[k][iA][iElementB+s]
							* dDxBasis1D[s][j];
					}

					dDaPsi /= dElementDeltaA;
					dDbPsi /= dElementDeltaB;

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
							pDataUpdate[k][iA+1][iB] -= dUpdateA;
						}
						pDataUpdate[k][iA][iB] += dUpdateA;
					}

					// Either set up communication with neighbor or apply
					// scalar hyperdiffusion along top edge.
					if (j == m_nHorizontalOrder-1) {
						if (b == nElementCountB-1) {
							connect.SetSendBuffer(
								Direction_Top, c, k, iA, dUpdateB);

						} else {
							pDataUpdate[k][iA][iB+1] -= dUpdateB;
						}
						pDataUpdate[k][iA][iB] += dUpdateB;
					}

					// Either set up communication with neighbor or apply
					// scalar hyperdiffusion along left edge.
					if (i == 0) {
						if (a == 0) {
							connect.SetSendBuffer(
								Direction_Left, c, k, iB, dUpdateA);
						} else {
							pDataUpdate[k][iA-1][iB] += dUpdateA;
						}
						pDataUpdate[k][iA][iB] -= dUpdateA;
					}

					// Either set up communication with neighbor or apply
					// scalar hyperdiffusion along bottom edge.
					if (j == 0) {
						if (b == 0) {
							connect.SetSendBuffer(
								Direction_Bottom, c, k, iA, dUpdateB);

						} else {
							pDataUpdate[k][iA][iB-1] += dUpdateB;
						}
						pDataUpdate[k][iA][iB] -= dUpdateB;
					}
				}
				}
			}
			}
			}
		}
	}

	// Perform a global exchange
	pGrid->ExchangeBuffers();
	pGrid->ExchangeBuffers();

	// Loop over all patches
	for (int n = 0; n < pGrid->GetActivePatchCount(); n++) {
		GridPatchGLL * pPatch =
			dynamic_cast<GridPatchGLL*>(pGrid->GetActivePatch(n));

		const PatchBox & box = pPatch->GetPatchBox();

		// Connectivity for patch
		Connectivity & connect = pPatch->GetConnectivity();

		// Data
		GridData4D & dataUpdateNode =
			pPatch->GetDataState(iDataUpdate, DataLocation_Node);
		GridData4D & dataUpdateREdge =
			pPatch->GetDataState(iDataUpdate, DataLocation_REdge);

		// Loop over all scalar components
		int nComponents = m_model.GetEquationSet().GetComponents();
		for (int c = 2; c < nComponents; c++) {


			int nElementCountR;

			double *** pDataState;
			double *** pDataUpdate;
			if (pGrid->GetVarLocation(c) == DataLocation_Node) {
				pDataUpdate = dataUpdateNode[c];
				nElementCountR = dataUpdateNode.GetRElements();

			} else if (pGrid->GetVarLocation(c) == DataLocation_REdge) {
				pDataUpdate = dataUpdateREdge[c];
				nElementCountR = dataUpdateREdge.GetRElements();

			} else {
				_EXCEPTIONT("UNIMPLEMENTED");
			}

			// Apply scalar hyperviscosity to right edge
			for (int k = 0; k < nElementCountR; k++) {
			for (int j = box.GetBInteriorBegin(); j < box.GetBInteriorEnd(); j++) {
				double dUpdateA =
					connect.GetRecvBuffer(Direction_Right, c, k, j);

				if (connect.IsCoordinateFlipped(Direction_Right, j)) {
					dUpdateA *= -1.0;
				}

				pDataUpdate[k][box.GetAInteriorEnd() - 1][j] += dUpdateA;
			}
			}

			// Apply scalar hyperviscosity to top edge
			for (int k = 0; k < nElementCountR; k++) {
			for (int i = box.GetAInteriorBegin(); i < box.GetAInteriorEnd(); i++) {
				double dUpdateB =
					connect.GetRecvBuffer(Direction_Top, c, k, i);

				if (connect.IsCoordinateFlipped(Direction_Top, i)) {
					dUpdateB *= -1.0;
				}

				pDataUpdate[k][i][box.GetBInteriorEnd()-1] += dUpdateB;
			}
			}

			// Apply scalar hyperviscosity to left edge
			for (int k = 0; k < nElementCountR; k++) {
			for (int j = box.GetBInteriorBegin(); j < box.GetBInteriorEnd(); j++) {
				double dUpdateA =
					connect.GetRecvBuffer(Direction_Left, c, k, j);

				if (connect.IsCoordinateFlipped(Direction_Left, j)) {
					dUpdateA *= -1.0;
				}

				pDataUpdate[k][box.GetAInteriorBegin()][j] -= dUpdateA;
			}
			}

			// Apply scalar hyperviscosity to bottom edge
			for (int k = 0; k < nElementCountR; k++) {
			for (int i = box.GetAInteriorBegin(); i < box.GetAInteriorEnd(); i++) {
				double dUpdateB =
					connect.GetRecvBuffer(Direction_Bottom, c, k, i);

				if (connect.IsCoordinateFlipped(Direction_Bottom, i)) {
					dUpdateB *= -1.0;
				}

				pDataUpdate[k][i][box.GetBInteriorBegin()] -= dUpdateB;
			}
			}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void HorizontalDynamicsFEM::ApplyVectorHyperdiffusionToBoundary(
	int iDataState,
	int iDataUpdate,
	double dDeltaT,
	double dNuDiff,
	double dNuVort,
	bool fScaleNuLocally
) {
}

///////////////////////////////////////////////////////////////////////////////

void HorizontalDynamicsFEM::ApplyRayleighFriction(
	int iDataUpdate,
	double dDeltaT
) {
	// Get a copy of the GLL grid
	GridGLL * pGrid = dynamic_cast<GridGLL*>(m_model.GetGrid());

	// Physical constants
	const PhysicalConstants & phys = m_model.GetPhysicalConstants();

	// Perform local update
	for (int n = 0; n < pGrid->GetActivePatchCount(); n++) {
		GridPatchGLL * pPatch =
			dynamic_cast<GridPatchGLL*>(pGrid->GetActivePatch(n));

		const PatchBox & box = pPatch->GetPatchBox();

		// Grid data
		GridData4D & dataUpdateNode =
			pPatch->GetDataState(iDataUpdate, DataLocation_Node);

		GridData4D & dataUpdateREdge =
			pPatch->GetDataState(iDataUpdate, DataLocation_REdge);

		// Reference state
		const GridData4D & dataReferenceNode =
			pPatch->GetReferenceState(DataLocation_Node);

		const GridData4D & dataReferenceREdge =
			pPatch->GetReferenceState(DataLocation_REdge);

		// Rayleigh friction strength
		const GridData3D & dataRayleighStrengthNode =
			pPatch->GetRayleighStrength(DataLocation_Node);

		const GridData3D & dataRayleighStrengthREdge =
			pPatch->GetRayleighStrength(DataLocation_REdge);

		// Loop over all nodes in patch
		for (int i = box.GetAInteriorBegin(); i < box.GetAInteriorEnd(); i++) {
		for (int j = box.GetBInteriorBegin(); j < box.GetBInteriorEnd(); j++) {

			// Rayleigh damping on nodes
			for (int k = 0; k < pGrid->GetRElements(); k++) {

				double dNu = dataRayleighStrengthNode[k][i][j];

				// Backwards Euler
				if (dNu == 0.0) {
					continue;
				}

				double dNuNode = 1.0 / (1.0 + dDeltaT * dNu);

				// Loop over all components
				int nComponents = m_model.GetEquationSet().GetComponents();
				for (int c = 0; c < nComponents; c++) {
					if (pGrid->GetVarLocation(c) == DataLocation_Node) {
						dataUpdateNode[c][k][i][j] =
							dNuNode * dataUpdateNode[c][k][i][j]
							+ (1.0 - dNuNode)
								* dataReferenceNode[c][k][i][j];
					}
				}
			}

			// Rayleigh damping on interfaces
			for (int k = 0; k <= pGrid->GetRElements(); k++) {

				double dNu = dataRayleighStrengthREdge[k][i][j];

				// Backwards Euler
				if (dNu == 0.0) {
					continue;
				}

				double dNuREdge = 1.0 / (1.0 + dDeltaT * dNu);

				// Loop over all components
				int nComponents = m_model.GetEquationSet().GetComponents();
				for (int c = 0; c < nComponents; c++) {
					if (pGrid->GetVarLocation(c) == DataLocation_REdge) {
						dataUpdateREdge[c][k][i][j] =
							dNuREdge * dataUpdateREdge[c][k][i][j]
							+ (1.0 - dNuREdge)
								* dataReferenceREdge[c][k][i][j];
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

	// Apply hyperdiffusion
	if (!m_fNoHyperdiffusion) {

		// Apply scalar and vector hyperdiffusion (first application)
		pGrid->ZeroData(iDataWorking, DataType_State);

		ApplyScalarHyperdiffusion(
			iDataInitial, iDataWorking, 1.0, 1.0, false);
		ApplyVectorHyperdiffusion(
			iDataInitial, iDataWorking, 1.0, 1.0, 1.0, false);

		// For the spectral element method apply DSS
		if (m_eHorizontalDynamicsType == SpectralElement) {
			pGrid->ApplyDSS(iDataWorking);

		// For Discontinous Galerkin apply hyperdiffusive fluxes
		// along boundaries
		} else if (m_eHorizontalDynamicsType == DiscontinuousGalerkin) {
			ApplyScalarHyperdiffusionToBoundary(
				iDataInitial, iDataWorking, 1.0, 1.0, false);
			ApplyVectorHyperdiffusionToBoundary(
				iDataInitial, iDataWorking, 1.0, 1.0, 1.0, false);
		}
		//pGrid->CopyData(iDataWorking, iDataUpdate, DataType_State);

#pragma message "Check sign of VectorHyperdiffusion"

		// Apply scalar and vector hyperdiffusion (second application)
		ApplyScalarHyperdiffusion(
			iDataWorking, iDataUpdate, dDeltaT, m_dNuScalar, true);
		ApplyVectorHyperdiffusion(
			iDataWorking, iDataUpdate, -dDeltaT, m_dNuDiv, m_dNuVort, true);

		// For the spectral element method apply DSS
		if (m_eHorizontalDynamicsType == SpectralElement) {
			pGrid->ApplyDSS(iDataUpdate);

		// For Discontinous Galerkin apply hyperdiffusive fluxes
		// along boundaries
		} else if (m_eHorizontalDynamicsType == DiscontinuousGalerkin) {
			ApplyScalarHyperdiffusionToBoundary(
				iDataWorking, iDataUpdate, dDeltaT, m_dNuScalar, true);
			ApplyVectorHyperdiffusionToBoundary(
				iDataWorking, iDataUpdate, -dDeltaT, m_dNuDiv, m_dNuVort, true);
		}
	}

	// Apply Rayleigh damping
	if (pGrid->HasRayleighFriction()) {
		ApplyRayleighFriction(iDataUpdate, dDeltaT);
	}
}

///////////////////////////////////////////////////////////////////////////////

