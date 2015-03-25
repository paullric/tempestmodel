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

#include "Defines.h"
#include "HorizontalDynamicsFEM.h"
#include "PhysicalConstants.h"
#include "Model.h"
#include "Grid.h"
#include "GaussLobattoQuadrature.h"

#include "Announce.h"
#include "GridGLL.h"
#include "GridPatchGLL.h"

#include "iomanip"

//#define DIFFERENTIAL_FORM

#ifdef DIFFERENTIAL_FORM
#pragma message "WARNING: DIFFERENTIAL_FORM will lose mass over topography"
#endif

///////////////////////////////////////////////////////////////////////////////

HorizontalDynamicsFEM::HorizontalDynamicsFEM(
	Model & model,
	int nHorizontalOrder,
	double dNuScalar,
	double dNuDiv,
	double dNuVort
) :
	HorizontalDynamics(model),
	m_nHorizontalOrder(nHorizontalOrder),
	m_dNuScalar(dNuScalar),
	m_dNuDiv(dNuDiv),
	m_dNuVort(dNuVort)
{

}

///////////////////////////////////////////////////////////////////////////////

void HorizontalDynamicsFEM::Initialize() {

	int nRElements = m_model.GetGrid()->GetRElements();

	// Initialize the alpha and beta mass fluxes
	m_dAlphaMassFlux.Initialize(
		m_nHorizontalOrder,
		m_nHorizontalOrder);

	m_dBetaMassFlux.Initialize(
		m_nHorizontalOrder,
		m_nHorizontalOrder);

	// Initialize the alpha and beta pressure fluxes
	m_dAlphaPressureFlux.Initialize(
		m_nHorizontalOrder,
		m_nHorizontalOrder);

	m_dBetaPressureFlux.Initialize(
		m_nHorizontalOrder,
		m_nHorizontalOrder);

#pragma message "Obtain number of auxiliary data elements from EquationSet type"
	// Auxiliary data
	m_dAuxDataNode.Initialize(
		9,
		nRElements,
		m_nHorizontalOrder,
		m_nHorizontalOrder);
/*

	m_dPressure.Initialize(
		m_nHorizontalOrder,
		m_nHorizontalOrder);

	m_dUx.Initialize(
		m_nHorizontalOrder,
		m_nHorizontalOrder);

	m_dCovUa.Initialize(
		m_nHorizontalOrder,
		m_nHorizontalOrder);

	m_dCovUb.Initialize(
		m_nHorizontalOrder,
		m_nHorizontalOrder);

	m_dCovUx.Initialize(
		m_nHorizontalOrder,
		m_nHorizontalOrder);

	m_dEnergy.Initialize(
		m_nHorizontalOrder,
		m_nHorizontalOrder);

	// Column storage
	m_dZeroColumn.Initialize(nRElements + 1);

	m_dColumnPressure.Initialize(nRElements);

	m_dColumnDxPressure.Initialize(nRElements);

	m_dColumnKineticEnergy.Initialize(nRElements);

	m_dColumnDxKineticEnergy.Initialize(nRElements);
*/
	// Initialize buffers for derivatives of Jacobian
	m_dJGradientA.Initialize(
		m_nHorizontalOrder,
		m_nHorizontalOrder);

	m_dJGradientB.Initialize(
		m_nHorizontalOrder,
		m_nHorizontalOrder);

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
/*
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
		const DataMatrix<double> & dCoriolisF =
			pPatch->GetCoriolisF();
		const DataMatrix<double> & dTopography =
			pPatch->GetTopography();

		const DataMatrix4D<double> & dCovMetricA =
			pPatch->GetCovMetricA();
		const DataMatrix4D<double> & dCovMetricB =
			pPatch->GetCovMetricB();

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
				m_dAlphaMassFlux[i][j] =
					dJacobian[k][iA][iB]
					* (dataInitialNode[HIx][k][iA][iB] - dTopography[iA][iB])
					* dataInitialNode[UIx][k][iA][iB];

				m_dBetaMassFlux[i][j] =
					dJacobian[k][iA][iB]
					* (dataInitialNode[HIx][k][iA][iB] - dTopography[iA][iB])
					* dataInitialNode[VIx][k][iA][iB];

				double dUa = dataInitialNode[UIx][k][iA][iB];
				double dUb = dataInitialNode[VIx][k][iA][iB];

				// Covariant velocities
				m_dCovUa[i][j] =
					  dCovMetricA[k][iA][iB][0] * dUa
					+ dCovMetricA[k][iA][iB][1] * dUb;

				m_dCovUb[i][j] =
					  dCovMetricB[k][iA][iB][0] * dUa
					+ dCovMetricB[k][iA][iB][1] * dUb;

				// Kinetic energy (KE)
				double dKE =
					0.5 * (m_dCovUa[i][j] * dUa + m_dCovUb[i][j] * dUb);

				// Pointwise pressure plus KE
				m_dPressure[i][j] =
					phys.GetG() * dataInitialNode[HIx][k][iA][iB] + dKE;
			}
			}

			// Pointwise update of horizontal velocities
			for (int i = 0; i < m_nHorizontalOrder; i++) {
			for (int j = 0; j < m_nHorizontalOrder; j++) {

				int iA = a * m_nHorizontalOrder + i + box.GetHaloElements();
				int iB = b * m_nHorizontalOrder + j + box.GetHaloElements();

				int iElementA = a * m_nHorizontalOrder + box.GetHaloElements();
				int iElementB = b * m_nHorizontalOrder + box.GetHaloElements();

				// Derivatives of the contravariant velocity field
				double dDaUa = 0.0;
				double dDaUb = 0.0;
				double dDbUa = 0.0;
				double dDbUb = 0.0;

				// Derivatives of the covariant velocity field
				double dCovDaUb = 0.0;
				double dCovDbUa = 0.0;

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
						m_dAlphaMassFlux[s][j]
						* dDxBasis1D[s][i];
#else
					// Update density: Variational formulation
					dLocalUpdateHa +=
						m_dAlphaMassFlux[s][j]
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

					// Derivative of covariant beta velocity wrt alpha
					dCovDaUb +=
						m_dCovUb[s][j]
						* dDxBasis1D[s][i];
				}

				double dLocalUpdateHb = 0.0;
				for (int s = 0; s < m_nHorizontalOrder; s++) {
#ifdef DIFFERENTIAL_FORM
					// Update density: Differential formulation
					dLocalUpdateHb -=
						m_dBetaMassFlux[i][s]
						* dDxBasis1D[s][j];
#else
					// Update density: Variational formulation
					dLocalUpdateHb +=
						m_dBetaMassFlux[i][s]
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

					// Derivative of alpha velocity wrt beta
					dCovDbUa +=
						m_dCovUa[i][s]
						* dDxBasis1D[s][j];
				}

				// Scale derivatives
				dDaUa /= dElementDeltaA;
				dDaUb /= dElementDeltaA;
				dDaP  /= dElementDeltaA;

				dDbUa /= dElementDeltaB;
				dDbUb /= dElementDeltaB;
				dDbP  /= dElementDeltaB;

				dCovDaUb /= dElementDeltaA;
				dCovDbUa /= dElementDeltaB;

				dLocalUpdateHa /= (dJacobian[k][iA][iB] * dElementDeltaA);
				dLocalUpdateHb /= (dJacobian[k][iA][iB] * dElementDeltaB);

				// Local update for momentum
				double dLocalUpdateUa = 0.0;
				double dLocalUpdateUb = 0.0;

				// Relative vorticity
				double dZeta = (dCovDaUb - dCovDbUa);

				// Rotational terms
				dLocalUpdateUa -=
					(dZeta + dCoriolisF[iA][iB] * dJacobian[k][iA][iB]) * (
						+ dContraMetricA[k][iA][iB][1] * dUa
						- dContraMetricA[k][iA][iB][0] * dUb);

				dLocalUpdateUb -=
					(dZeta + dCoriolisF[iA][iB] * dJacobian[k][iA][iB]) * (
						+ dContraMetricB[k][iA][iB][1] * dUa
						- dContraMetricB[k][iA][iB][0] * dUb);

				// Gradient of pressure
				dLocalUpdateUa -=
					+ dContraMetricA[k][iA][iB][0] * dDaP
					+ dContraMetricA[k][iA][iB][1] * dDbP;

				dLocalUpdateUb -=
					+ dContraMetricB[k][iA][iB][0] * dDaP
					+ dContraMetricB[k][iA][iB][1] * dDbP;

				// Apply update
				dataUpdateNode[UIx][k][iA][iB] +=
					dDeltaT * dLocalUpdateUa;
				dataUpdateNode[VIx][k][iA][iB] +=
					dDeltaT * dLocalUpdateUb;

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

void HorizontalDynamicsFEM::StepNonhydrostaticPrimitive(
	int iDataInitial,
	int iDataUpdate,
	const Time & time,
	double dDeltaT
) {
	// Get a copy of the GLL grid
	GridGLL * pGrid = dynamic_cast<GridGLL*>(m_model.GetGrid());

	// Number of vertical elements
	const int nRElements = pGrid->GetRElements();

	// Physical constants
	const PhysicalConstants & phys = m_model.GetPhysicalConstants();

	// Indices of EquationSet variables
	const int UIx = 0;
	const int VIx = 1;
	const int PIx = 2;
	const int WIx = 3;
	const int RIx = 4;

	// Indices of auxiliary data variables
	const int ConUaIx = 0;
	const int ConUbIx = 1;
	const int ConUxIx = 2;
	const int CovUxIx = 3;
	const int KIx = 4;
	const int PFluxIx = 5;
	const int RFluxIx = 6;

	// Vertical level stride in local data arrays
	const int nVerticalElementStride =
		m_nHorizontalOrder * m_nHorizontalOrder;

	//std::cout << "Inside the horizontal step! \n";

	// Perform local update
	for (int n = 0; n < pGrid->GetActivePatchCount(); n++) {
		GridPatchGLL * pPatch =
			dynamic_cast<GridPatchGLL*>(pGrid->GetActivePatch(n));

		const PatchBox & box = pPatch->GetPatchBox();

		const DataMatrix3D<double> & dElementArea =
			pPatch->GetElementArea();
		const DataMatrix<double> & dJacobian2D =
			pPatch->GetJacobian2D();
		const DataMatrix3D<double> & dJacobian =
			pPatch->GetJacobian();
		const DataMatrix4D<double> & dContraMetricA =
			pPatch->GetContraMetricA();
		const DataMatrix4D<double> & dContraMetricB =
			pPatch->GetContraMetricB();
		const DataMatrix4D<double> & dContraMetricXi =
			pPatch->GetContraMetricXi();
		const DataMatrix4D<double> & dDerivRNode =
			pPatch->GetDerivRNode();

		const DataMatrix4D<double> & dOrthonormNode =
			pPatch->GetOrthonormNode();
		const DataMatrix4D<double> & dOrthonormREdge =
			pPatch->GetOrthonormREdge();

		const DataMatrix3D<double> & dCovMetric2DA =
			pPatch->GetCovMetric2DA();
		const DataMatrix3D<double> & dCovMetric2DB =
			pPatch->GetCovMetric2DB();

		const DataMatrix4D<double> & dCovMetricA =
			pPatch->GetCovMetricA();
		const DataMatrix4D<double> & dCovMetricB =
			pPatch->GetCovMetricB();
		const DataMatrix4D<double> & dCovMetricXi =
			pPatch->GetCovMetricXi();

		const DataMatrix<double> & dCoriolisF =
			pPatch->GetCoriolisF();
		const DataMatrix<double> & dTopography =
			pPatch->GetTopography();

		const double dZtop = pGrid->GetZtop();
		
		// Get the coordinates to check force balance externally
		const DataMatrix<double>   & dLongitude  = pPatch->GetLongitude();
		const DataMatrix<double>   & dLatitude   = pPatch->GetLatitude();
		const DataMatrix3D<double> & dZLevels    = pPatch->GetZLevels();

		// Data
		GridData4D & dataInitialNode =
			pPatch->GetDataState(iDataInitial, DataLocation_Node);

		GridData4D & dataUpdateNode =
			pPatch->GetDataState(iDataUpdate, DataLocation_Node);

		GridData4D & dataInitialREdge =
			pPatch->GetDataState(iDataInitial, DataLocation_REdge);

		GridData4D & dataUpdateREdge =
			pPatch->GetDataState(iDataUpdate, DataLocation_REdge);

		const int nVerticalStateStride =
			dataInitialNode.GetAElements() * dataInitialNode.GetBElements();

		// Perform interpolations as required due to vertical staggering
		if (pGrid->GetVarsAtLocation(DataLocation_REdge) != 0) {

			_EXCEPTIONT("Not implemented");

			// Interpolate W to model levels
			if (pGrid->GetVarLocation(WIx) == DataLocation_REdge) {
				pPatch->InterpolateREdgeToNode(WIx, iDataInitial);
			}

			// Interpolate Theta to model levels
			if (pGrid->GetVarLocation(PIx) == DataLocation_REdge) {
				pPatch->InterpolateREdgeToNode(PIx, iDataInitial);
			}
			if ((pGrid->GetVarLocation(PIx) == DataLocation_Node) &&
				(pGrid->GetVarLocation(WIx) == DataLocation_REdge)
			) {
				pPatch->InterpolateNodeToREdge(PIx, iDataInitial);
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

		// Get number of finite elements in each coordinate direction
		int nElementCountA = pPatch->GetElementCountA();
		int nElementCountB = pPatch->GetElementCountB();

		// Kinetic energy data
		GridData3D & dataKineticEnergy = pPatch->GetDataKineticEnergy();
		GridData3D & dataDxKineticEnergy = pPatch->GetDataDxKineticEnergy();

		// Loop over all elements
		for (int a = 0; a < nElementCountA; a++) {
		for (int b = 0; b < nElementCountB; b++) {
/*
			// Apply boundary conditions to W
			if (pGrid->GetVerticalStaggering() ==
				Grid::VerticalStaggering_Interfaces
			) {
				for (int i = 0; i < m_nHorizontalOrder; i++) {
				for (int j = 0; j < m_nHorizontalOrder; j++) {

					int iA = a * m_nHorizontalOrder + i + box.GetHaloElements();
					int iB = b * m_nHorizontalOrder + j + box.GetHaloElements();

					dataInitialNode[WIx][0][iA][iB] =
						pPatch->CalculateNoFlowUrNode(
							0, iA, iB,
							dataInitialNode[UIx][0][iA][iB],
							dataInitialNode[VIx][0][iA][iB]);

					dataInitialNode[WIx][nRElements-1][iA][iB] =
						pPatch->CalculateNoFlowUrNode(
							nRElements-1, iA, iB,
							dataInitialNode[UIx][nRElements-1][iA][iB],
							dataInitialNode[VIx][nRElements-1][iA][iB]);

					dataUpdateNode[WIx][0][iA][iB] =
						pPatch->CalculateNoFlowUrNode(
							0, iA, iB,
							dataUpdateNode[UIx][0][iA][iB],
							dataUpdateNode[VIx][0][iA][iB]);

					dataUpdateNode[WIx][nRElements-1][iA][iB] =
						pPatch->CalculateNoFlowUrNode(
							nRElements-1, iA, iB,
							dataUpdateNode[UIx][nRElements-1][iA][iB],
							dataUpdateNode[VIx][nRElements-1][iA][iB]);

					dataInitialNode[WIx][0][iA][iB] =
						- ( dContraMetricXi[0][iA][iB][0]
							* dataInitialNode[UIx][0][iA][iB]
						+ dContraMetricXi[0][iA][iB][1]
							* dataInitialNode[VIx][0][iA][iB])
						/ dContraMetricXi[0][iA][iB][2]
						/ dDerivRNode[0][iA][iB][2];

					dataInitialNode[WIx][nRElements-1][iA][iB] =
						- ( dContraMetricXi[nRElements-1][iA][iB][0]
							* dataInitialNode[UIx][nRElements-1][iA][iB]
						+ dContraMetricXi[nRElements-1][iA][iB][1]
							* dataInitialNode[VIx][nRElements-1][iA][iB])
						/ dContraMetricXi[nRElements-1][iA][iB][2]
						/ dDerivRNode[nRElements-1][iA][iB][2];

					dataUpdateNode[WIx][0][iA][iB] =
						- ( dContraMetricXi[0][iA][iB][0]
							* dataUpdateNode[UIx][0][iA][iB]
						+ dContraMetricXi[0][iA][iB][1]
							* dataUpdateNode[VIx][0][iA][iB])
						/ dContraMetricXi[0][iA][iB][2]
						/ dDerivRNode[0][iA][iB][2];

					dataUpdateNode[WIx][nRElements-1][iA][iB] =
						- ( dContraMetricXi[nRElements-1][iA][iB][0]
							* dataUpdateNode[UIx][nRElements-1][iA][iB]
						+ dContraMetricXi[nRElements-1][iA][iB][1]
							* dataUpdateNode[VIx][nRElements-1][iA][iB])
						/ dContraMetricXi[nRElements-1][iA][iB][2]
						/ dDerivRNode[nRElements-1][iA][iB][2];

				}
				}
			}
*/
			// Compute auxiliary data in element
			for (int k = 0; k < nRElements; k++) {
			for (int i = 0; i < m_nHorizontalOrder; i++) {
			for (int j = 0; j < m_nHorizontalOrder; j++) {

				int iA = a * m_nHorizontalOrder + i + box.GetHaloElements();
				int iB = b * m_nHorizontalOrder + j + box.GetHaloElements();

				int iElementA = a * m_nHorizontalOrder + box.GetHaloElements();
				int iElementB = b * m_nHorizontalOrder + box.GetHaloElements();

				// Contravariant velocities
				double dCovUa = dataInitialNode[UIx][k][iA][iB];
				double dCovUb = dataInitialNode[VIx][k][iA][iB];
				double dCovUr = dataInitialNode[WIx][k][iA][iB];

				// Calculate Ux
				double dCovUx =
					  dataInitialNode[WIx][k][iA][iB]
					* dDerivRNode[k][iA][iB][2];

				// Covariant xi velocity
				m_dAuxDataNode[CovUxIx][k][i][j] = dCovUx;

				// Contravariant velocities
				m_dAuxDataNode[ConUaIx][k][i][j] =
					  dContraMetricA[k][iA][iB][0] * dCovUa
					+ dContraMetricA[k][iA][iB][1] * dCovUb
					+ dContraMetricA[k][iA][iB][2] * dCovUx;

				m_dAuxDataNode[ConUbIx][k][i][j] =
					  dContraMetricB[k][iA][iB][0] * dCovUa
					+ dContraMetricB[k][iA][iB][1] * dCovUb
					+ dContraMetricB[k][iA][iB][2] * dCovUx;

				m_dAuxDataNode[ConUxIx][k][i][j] =
					  dContraMetricXi[k][iA][iB][0] * dCovUa
					+ dContraMetricXi[k][iA][iB][1] * dCovUb
					+ dContraMetricXi[k][iA][iB][2] * dCovUx;

				// Specific kinetic energy
				m_dAuxDataNode[KIx][k][i][j] = 0.5 * (
					  m_dAuxDataNode[ConUaIx][k][i][j] * dCovUa
					+ m_dAuxDataNode[ConUbIx][k][i][j] * dCovUb
					+ m_dAuxDataNode[ConUxIx][k][i][j] * dCovUx);
			}
			}
			}

			// Update all elements
			for (int k = 0; k < nRElements; k++) {

				// Pointwise fluxes and pressure within spectral element
				for (int i = 0; i < m_nHorizontalOrder; i++) {
				for (int j = 0; j < m_nHorizontalOrder; j++) {

					int iA = a * m_nHorizontalOrder + i + box.GetHaloElements();
					int iB = b * m_nHorizontalOrder + j + box.GetHaloElements();

					// Base fluxes (area times velocity)
					double dAlphaBaseFlux =
						dJacobian[k][iA][iB]
						* m_dAuxDataNode[ConUaIx][k][i][j];

					double dBetaBaseFlux =
						dJacobian[k][iA][iB]
						* m_dAuxDataNode[ConUbIx][k][i][j];

					// Density flux
					m_dAlphaMassFlux[i][j] =
						  dAlphaBaseFlux
						* dataInitialNode[RIx][k][iA][iB];

					m_dBetaMassFlux[i][j] =
						  dBetaBaseFlux
						* dataInitialNode[RIx][k][iA][iB];

					// Pressure flux
					m_dAlphaPressureFlux[i][j] =
						  dAlphaBaseFlux
						* phys.GetGamma()
						* dataInitialNode[PIx][k][iA][iB];

					m_dBetaPressureFlux[i][j] =
						  dBetaBaseFlux
						* phys.GetGamma()
						* dataInitialNode[PIx][k][iA][iB];
				}
				}

				double dElementDeltaAE = 0.0;
				double dElementDeltaTE = 0.0;
				double dElementDeltaKE = 0.0;
				double dElementDeltaIE = 0.0;

				// Pointwise update of quantities on model levels
				for (int i = 0; i < m_nHorizontalOrder; i++) {
				for (int j = 0; j < m_nHorizontalOrder; j++) {

					int iA = a * m_nHorizontalOrder + i + box.GetHaloElements();
					int iB = b * m_nHorizontalOrder + j + box.GetHaloElements();

					int iElementA = a * m_nHorizontalOrder + box.GetHaloElements();
					int iElementB = b * m_nHorizontalOrder + box.GetHaloElements();

					// Derivatives of the covariant velocity field
					double dCovDaUb = 0.0;
					double dCovDaUx = 0.0;

					double dCovDbUa = 0.0;
					double dCovDbUx = 0.0;

					// Derivative of the kinetic energy
					double dDaKE = 0.0;
					double dDbKE = 0.0;

					// Derivatives of the pressure field
					double dDaP = 0.0;
					double dDbP = 0.0;

					// Aliases for alpha and beta velocities
					const double dConUa = m_dAuxDataNode[ConUaIx][k][i][j];
					const double dConUb = m_dAuxDataNode[ConUbIx][k][i][j];

					const double dCovUx = m_dAuxDataNode[CovUxIx][k][i][j];

					const double dConUx =
						  dContraMetricXi[k][iA][iB][0]
							* dataInitialNode[UIx][0][iA][iB]
						+ dContraMetricXi[k][iA][iB][1]
							* dataInitialNode[VIx][0][iA][iB]
						+ dContraMetricXi[k][iA][iB][2]
							* m_dAuxDataNode[CovUxIx][k][i][j];

					// Calculate derivatives in the alpha direction
					double dDaRhoFluxA = 0.0;
					double dDaPressureFluxA = 0.0;

					for (int s = 0; s < m_nHorizontalOrder; s++) {
#ifdef DIFFERENTIAL_FORM
						// Update density: Differential formulation
						dDaRhoFluxA +=
							m_dAlphaMassFlux[s][j]
							* dDxBasis1D[s][i];

						// Update pressure: Differential formulation
						dDaPressureFluxA +=
							m_dAlphaPressureFlux[s][j]
							* dDxBasis1D[s][i];

#else
						// Update density: Variational formulation
						dDaRhoFluxA -=
							m_dAlphaMassFlux[s][j]
							* dStiffness1D[i][s];

						// Update pressure: Variational formulation
						dDaPressureFluxA -=
							m_dAlphaPressureFlux[s][j]
							* dStiffness1D[i][s];
#endif

						// Derivative of pressure with respect to alpha
						dDaP +=
							dataInitialNode[PIx][k][iElementA+s][iB]
							* dDxBasis1D[s][i];

						// Derivative of covariant beta velocity wrt alpha
						dCovDaUb +=
							dataInitialNode[VIx][k][iElementA+s][iB]
							* dDxBasis1D[s][i];

						// Derivative of covariant xi velocity wrt alpha
						dCovDaUx +=
							m_dAuxDataNode[CovUxIx][k][s][j]
							* dDxBasis1D[s][i];

						// Derivative of specific kinetic energy wrt alpha
						dDaKE +=
							m_dAuxDataNode[KIx][k][s][j]
							* dDxBasis1D[s][i];
					}

					// Calculate derivatives in the beta direction
					double dDbRhoFluxB = 0.0;
					double dDbPressureFluxB = 0.0;

					for (int s = 0; s < m_nHorizontalOrder; s++) {
#ifdef DIFFERENTIAL_FORM
						// Update density: Differential formulation
						dDbRhoFluxB +=
							m_dBetaMassFlux[i][s]
							* dDxBasis1D[s][j];

						// Update pressure: Differential formulation
						dDbPressureFluxB +=
							m_dBetaPressureFlux[i][s]
							* dDxBasis1D[s][j];

#else
						// Update density: Variational formulation
						dDbRhoFluxB -=
							m_dBetaMassFlux[i][s]
							* dStiffness1D[j][s];

						// Update pressure: Variational formulation
						dDbPressureFluxB -=
							m_dBetaPressureFlux[i][s]
							* dStiffness1D[j][s];
#endif

						// Derivative of pressure with respect to beta
						dDbP +=
							dataInitialNode[PIx][k][iA][iElementB+s]
							* dDxBasis1D[s][j];

						// Derivative of covariant alpha velocity wrt beta
						dCovDbUa +=
							dataInitialNode[VIx][k][iA][iElementB+s]
							* dDxBasis1D[s][j];

						// Derivative of covariant xi velocity wrt beta
						dCovDbUx +=
							m_dAuxDataNode[CovUxIx][k][i][s]
							* dDxBasis1D[s][j];

						// Derivative of specific kinetic energy wrt beta
						dDbKE +=
							m_dAuxDataNode[KIx][k][i][s]
							* dDxBasis1D[s][j];
					}

					// Scale derivatives
					dDaP  /= dElementDeltaA;
					dDbP  /= dElementDeltaB;

					dCovDaUb /= dElementDeltaA;
					dCovDaUx /= dElementDeltaA;
					dDaKE    /= dElementDeltaA;

					dCovDbUa /= dElementDeltaB;
					dCovDbUx /= dElementDeltaB;
					dDbKE    /= dElementDeltaB;

					dDaRhoFluxA /= dElementDeltaA;
					dDbRhoFluxB /= dElementDeltaB;

					dDaPressureFluxA /= dElementDeltaA;
					dDbPressureFluxB /= dElementDeltaB;

					// Vertical derivatives
					double dCovDxUa =
						pGrid->DifferentiateNodeToNode(
							&(dataInitialNode[UIx][0][iA][iB]),
							k, nVerticalStateStride);

					double dCovDxUb =
						pGrid->DifferentiateNodeToNode(
							&(dataInitialNode[VIx][0][iA][iB]),
							k, nVerticalStateStride);
/*
					double dDxKE =
						pGrid->DifferentiateNodeToNode(
							&(m_dAuxDataNode[KIx][0][i][j]),
							k, nVerticalElementStride);

					double dDxRhoFluxXi =
						pGrid->DifferentiateNodeToNode(
							&(m_dAuxDataNode[RFluxIx][0][i][j]),
							k, nVerticalElementStride);

					double dDxPressureFluxXi =
						pGrid->DifferentiateNodeToNode(
							&(m_dAuxDataNode[PFluxIx][0][i][j]),
							k, nVerticalElementStride);

					double dDxP =
						pGrid->DifferentiateNodeToNode(
							&(dataInitialNode[PIx][0][iA][iB]),
							k, nVerticalStateStride);
*/
					// Pointwise horizontal momentum update
					double dLocalUpdateUa = 0.0;
					double dLocalUpdateUb = 0.0;
					double dLocalUpdateUx = 0.0;

					// Relative vorticity
					double dZetaA = (dCovDbUx - dCovDxUb);
					double dZetaB = (dCovDxUa - dCovDaUx);
					double dZetaXi = (dCovDaUb - dCovDbUa);

					// Rotational terms (covariant)
					double dCovUCrossZetaA  = dConUb * dZetaXi - dConUx * dZetaB;
					double dCovUCrossZetaB  = dConUx * dZetaA  - dConUa * dZetaXi;
					double dCovUCrossZetaXi = dConUa * dZetaB  - dConUb * dZetaA;

					// u^xi must be zero at boundaries
					if ((k == 0) || (k == nRElements-1)) {
						dCovUCrossZetaXi = 0.0;
					}

					// Updates due to rotational terms
					dLocalUpdateUa += dCovUCrossZetaA;

					dLocalUpdateUb += dCovUCrossZetaB;

					dLocalUpdateUx += dCovUCrossZetaXi;
/*
					// Updates due to rotational terms
					dLocalUpdateUa -=
						(dZetaXi + dCoriolisF[iA][iB] * dJacobian2D[iA][iB]) * (
							+ dContraMetricA[k][iA][iB][1] * dUa
							- dContraMetricA[k][iA][iB][0] * dUb);

					dLocalUpdateUb -=
						(dZetaXi + dCoriolisF[iA][iB] * dJacobian2D[iA][iB]) * (
							+ dContraMetricB[k][iA][iB][1] * dUa
							- dContraMetricB[k][iA][iB][0] * dUb);
*/
					// Gravity
					double dDaPhi = phys.GetG() * dDerivRNode[k][iA][iB][0];
					double dDbPhi = phys.GetG() * dDerivRNode[k][iA][iB][1];
					double dDxPhi = phys.GetG() * dDerivRNode[k][iA][iB][2];

					// Horizontal updates
					double dDaUpdate =
						  dDaP / dataInitialNode[RIx][k][iA][iB]
						+ dDaKE + dDaPhi;

					double dDbUpdate =
						  dDbP / dataInitialNode[RIx][k][iA][iB]
						+ dDbKE + dDbPhi;

					dLocalUpdateUa -= dDaUpdate;
					dLocalUpdateUb -= dDbUpdate;

					// Vertical velocity update
					double dLocalUpdateUr = 0.0;

					// Apply horizontal updates
					if (k != 0) {
						dLocalUpdateUr =
							dLocalUpdateUx / dDerivRNode[k][iA][iB][2];

					// Apply horizontal updates at lower boundary
					} else {
						//dLocalUpdateUr =
						//	dLocalUpdateUx / dDerivRNode[k][iA][iB][2];
/*
						if (((iA > 20) && (iA < 30)) && (k == 0)) {
							printf("%1.15e %1.15e\n",
								dLocalUpdateUx / dDerivRNode[k][iA][iB][2],
								- ( dContraMetricXi[0][iA][iB][0] * dDaUpdate
								  + dContraMetricXi[0][iA][iB][1] * dDbUpdate)
								/ dContraMetricXi[0][iA][iB][2]
								/ dDerivRNode[0][iA][iB][2]);

						}
*/

						dLocalUpdateUr =
							- ( dContraMetricXi[0][iA][iB][0] * dLocalUpdateUa
							  + dContraMetricXi[0][iA][iB][1] * dLocalUpdateUb)
							/ dContraMetricXi[0][iA][iB][2]
							/ dDerivRNode[0][iA][iB][2];
					}

					// Ensure no vertical velocity update in topmost element
					if (k == nRElements - 1) {
						dLocalUpdateUr = 0.0;
					}

 					// Apply update to horizontal velocity on model levels
					dataUpdateNode[UIx][k][iA][iB] +=
						dDeltaT * dLocalUpdateUa;
					//dataUpdateNode[VIx][k][iA][iB] +=
					//	dDeltaT * dLocalUpdateUb;
					dataUpdateNode[WIx][k][iA][iB] +=
						dDeltaT * dLocalUpdateUr;

					// Check boundary condition
					if (k == 0) {
						double dConUx =
							dContraMetricXi[0][iA][iB][0]
								* dataUpdateNode[UIx][k][iA][iB]
							+ dContraMetricXi[0][iA][iB][1]
								* dataUpdateNode[VIx][k][iA][iB]
							+ dContraMetricXi[0][iA][iB][2]
								* dDerivRNode[0][iA][iB][2]
								* dataUpdateNode[WIx][k][iA][iB];

						if (fabs(dConUx) > 1.0e-14) {
							printf("%1.15e\n", dConUx);
							_EXCEPTIONT("Boundary condition failure");
						}
					}

					// Update density on model levels
					dataUpdateNode[RIx][k][iA][iB] -=
						dDeltaT / dJacobian[k][iA][iB] * (
							  dDaRhoFluxA
							+ dDbRhoFluxB);
							//+ dDxRhoFluxXi);

					// Update pressure on model levels
					dataUpdateNode[PIx][k][iA][iB] +=
						dDeltaT * (phys.GetGamma() - 1.0)
						* (dConUa * dDaP + dConUb * dDbP);
							//+ dUx * dDxP);

					dataUpdateNode[PIx][k][iA][iB] -=
						dDeltaT / dJacobian[k][iA][iB] * (
							  dDaPressureFluxA
							+ dDbPressureFluxB);
							//+ dDxPressureFluxXi);

					if ((iA == 35) && (iB == 1)) {
/*
						printf("%i %1.5e %1.5e %1.5e %1.5e %1.5e\n",
							k,
							dataUpdateNode[UIx][k][iA][iB],
							dataUpdateNode[VIx][k][iA][iB],
							dataUpdateNode[WIx][k][iA][iB],
							dataUpdateNode[RIx][k][iA][iB],
							dataUpdateNode[PIx][k][iA][iB]);
*/
/*
						printf("%i %1.5e %1.5e %1.5e\n",
							k,
							dDxKE / dDerivRNode[k][iA][iB][2],
							dDxUpdate / dDerivRNode[k][iA][iB][2],
							dLocalUpdateUr);
*/
						if (k == nRElements-1) {
							//_EXCEPTION();
						}
					}

					// Update tracers
				}
				}
			}
/*
			// Apply boundary conditions to W
			if (pGrid->GetVerticalStaggering() ==
				Grid::VerticalStaggering_Interfaces
			) {
				for (int i = 0; i < m_nHorizontalOrder; i++) {
				for (int j = 0; j < m_nHorizontalOrder; j++) {

					int iA = a * m_nHorizontalOrder + i + box.GetHaloElements();
					int iB = b * m_nHorizontalOrder + j + box.GetHaloElements();

					// 2D Covariant velocities
					double dUa = dataUpdateNode[UIx][0][iA][iB];
					double dUb = dataUpdateNode[VIx][0][iA][iB];
					double dUr = dataUpdateNode[WIx][0][iA][iB];

					double d2DCovUa =
						  dCovMetric2DA[iA][iB][0] * dUa
						+ dCovMetric2DA[iA][iB][1] * dUb;

					double d2DCovUb =
						  dCovMetric2DB[iA][iB][0] * dUa
						+ dCovMetric2DB[iA][iB][1] * dUb;

					// Initial Kinetic energy
					double dKE =
						0.5 * (d2DCovUa * dUa + d2DCovUb * dUb + dUr * dUr);

					// Solutions preserving kinetic energy
					double dUaF = sqrt(2.0 * dKE / dCovMetricA[0][iA][iB][0]);

					if (dUaF * dUa > 0.0) {
						dataUpdateNode[UIx][0][iA][iB] = dUaF;
					} else {
						dataUpdateNode[UIx][0][iA][iB] = - dUaF;
					}

					dataUpdateNode[WIx][0][iA][iB] =
						dDerivRNode[0][iA][iB][0]
						* dataUpdateNode[UIx][0][iA][iB];
				}
				}
			}
*/
		}
		}

/*
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
						if (pGrid->GetVarLocation(PIx) != DataLocation_REdge) {
							break;
						}

						_EXCEPTION();

						double dDaTheta = 0.0;
						double dDbTheta = 0.0;

						const double dUaREdge = dataInitialREdge[UIx][k][iA][iB];
						const double dUbREdge = dataInitialREdge[VIx][k][iA][iB];

						for (int s = 0; s < m_nHorizontalOrder; s++) {
							// Derivative of theta with respect to alpha
							dDaTheta +=
								dataInitialREdge[PIx][k][iElementA+s][iB]
								* dDxBasis1D[s][i];

							// Derivative of theta with respect to beta
							dDbTheta +=
								dataInitialREdge[PIx][k][iA][iElementB+s]
								* dDxBasis1D[s][j];
						}

						// Scale derivatives
						dDaTheta /= dElementDeltaA;
						dDbTheta /= dElementDeltaB;

						// Update potential temperature
						dataUpdateREdge[PIx][k][iA][iB] -=
							dDeltaT * (dUaREdge * dDaTheta + dUbREdge * dDbTheta);
					}
				}
				}
			}
			}
*/
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

	// Perform local update
	for (int n = 0; n < pGrid->GetActivePatchCount(); n++) {
		GridPatchGLL * pPatch =
			dynamic_cast<GridPatchGLL*>(pGrid->GetActivePatch(n));

		const PatchBox & box = pPatch->GetPatchBox();

		const DataMatrix3D<double> & dJacobianNode =
			pPatch->GetJacobian();
		const DataMatrix3D<double> & dJacobianREdge =
			pPatch->GetJacobianREdge();
		const DataMatrix3D<double> & dContraMetricA =
			pPatch->GetContraMetric2DA();
		const DataMatrix3D<double> & dContraMetricB =
			pPatch->GetContraMetric2DB();

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
			double *** pJacobian;
			if (pGrid->GetVarLocation(c) == DataLocation_Node) {
				pDataInitial = dataInitialNode[c];
				pDataUpdate  = dataUpdateNode[c];
				nElementCountR = dataInitialNode.GetRElements();
				pJacobian = dJacobianNode;

			} else if (pGrid->GetVarLocation(c) == DataLocation_REdge) {
				pDataInitial = dataInitialREdge[c];
				pDataUpdate  = dataUpdateREdge[c];
				nElementCountR = dataInitialREdge.GetRElements();
				pJacobian = dJacobianREdge;

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

					m_dJGradientA[i][j] = pJacobian[k][iA][iB] * (
						+ dContraMetricA[iA][iB][0] * dDaPsi
						+ dContraMetricA[iA][iB][1] * dDbPsi);

					m_dJGradientB[i][j] = pJacobian[k][iA][iB] * (
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
					double dInvJacobian = 1.0 / pJacobian[k][iA][iB];

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
	const int WIx = 3;

	// Get a copy of the GLL grid
	GridGLL * pGrid = dynamic_cast<GridGLL*>(m_model.GetGrid());

	// Loop over all patches
	for (int n = 0; n < pGrid->GetActivePatchCount(); n++) {
		GridPatchGLL * pPatch =
			dynamic_cast<GridPatchGLL*>(pGrid->GetActivePatch(n));

		const PatchBox & box = pPatch->GetPatchBox();

		const DataMatrix4D<double> & dContraMetricA =
			pPatch->GetContraMetricA();
		const DataMatrix4D<double> & dContraMetricB =
			pPatch->GetContraMetricB();
		const DataMatrix4D<double> & dContraMetricXi =
			pPatch->GetContraMetricXi();
		const DataMatrix4D<double> & dDerivRNode =
			pPatch->GetDerivRNode();

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

		// Compute curl and divergence of U on the grid
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
				double dDaDiv = 0.0;
				double dDbDiv = 0.0;

				double dDaCurl = 0.0;
				double dDbCurl = 0.0;

				for (int s = 0; s < m_nHorizontalOrder; s++) {
					dDaDiv -=
						  dStiffness1D[i][s]
						* dataDiv[k][iElementA+s][iB];

					dDbDiv -=
						  dStiffness1D[j][s]
						* dataDiv[k][iA][iElementB+s];

					dDaCurl -=
						  dStiffness1D[i][s]
						* dataCurl[k][iElementA+s][iB];

					dDbCurl -=
						  dStiffness1D[j][s]
						* dataCurl[k][iA][iElementB+s];
				}

				dDaDiv /= dElementDeltaA;
				dDbDiv /= dElementDeltaB;

				dDaCurl /= dElementDeltaA;
				dDbCurl /= dElementDeltaB;

				// Apply update
				double dUpdateUa =
					+ dLocalNuDiv * dDaDiv
					+ dLocalNuVort * (
						  dContraMetricA[k][iA][iB][0] * dDaCurl
						+ dContraMetricA[k][iA][iB][1] * dDbCurl);

				double dUpdateUb =
					+ dLocalNuDiv * dDbDiv
					+ dLocalNuVort * (
						  dContraMetricB[k][iA][iB][0] * dDaCurl
						+ dContraMetricB[k][iA][iB][1] * dDbCurl);

				dataUpdate[UIx][k][iA][iB] -= dDeltaT * dUpdateUa;

				//dataUpdate[VIx][k][iA][iB] -= dDeltaT * dUpdateUb;

				if (k == 0) {
					double dUpdateUr =
						- ( dContraMetricXi[0][iA][iB][0] * dUpdateUa
						  + dContraMetricXi[0][iA][iB][1] * dUpdateUb)
						/ dContraMetricXi[0][iA][iB][2]
						/ dDerivRNode[0][iA][iB][2];

					dataUpdate[WIx][k][iA][iB] -= dDeltaT * dUpdateUr;
				}
			}
			}
		}
		}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void HorizontalDynamicsFEM::ApplyRayleighFriction(
	int iDataUpdate,
	double dDeltaT
) {
	// Get a copy of the GLL grid
	GridGLL * pGrid = dynamic_cast<GridGLL*>(m_model.GetGrid());

	// Number of components to hit with friction
	int nComponents = m_model.GetEquationSet().GetComponents();

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

		// Apply boundary conditions
		pPatch->ApplyBoundaryConditions(iDataUpdate);
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

	// No hyperdiffusion
	if ((m_dNuScalar == 0.0) && (m_dNuDiv == 0.0) && (m_dNuVort == 0.0)) {

	// Apply hyperdiffusion
	} else {

		// Apply scalar and vector hyperdiffusion (first application)
		pGrid->ZeroData(iDataWorking, DataType_State);
		//pGrid->CopyData(iDataInitial, iDataWorking, DataType_State);

		ApplyScalarHyperdiffusion(
			iDataInitial, iDataWorking, 1.0, 1.0, false);
		ApplyVectorHyperdiffusion(
			iDataInitial, iDataWorking, 1.0, 1.0, 1.0, false);

		// Apply Direct Stiffness Summation
		pGrid->ApplyDSS(iDataWorking);

		// Apply scalar and vector hyperdiffusion (second application)
		pGrid->CopyData(iDataInitial, iDataUpdate, DataType_State);

		ApplyScalarHyperdiffusion(
			iDataWorking, iDataUpdate, -dDeltaT, m_dNuScalar, true);
		ApplyVectorHyperdiffusion(
			iDataWorking, iDataUpdate, -dDeltaT, m_dNuDiv, m_dNuVort, true);

		// Apply Direct Stiffness Summation
		pGrid->ApplyDSS(iDataUpdate);
	}

#ifdef APPLY_RAYLEIGH_WITH_HYPERVIS
	// Apply Rayleigh damping
	if (pGrid->HasRayleighFriction()) {
		ApplyRayleighFriction(iDataUpdate, dDeltaT);
	}
#endif
}

///////////////////////////////////////////////////////////////////////////////

