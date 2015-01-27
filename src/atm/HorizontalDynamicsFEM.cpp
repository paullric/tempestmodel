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

//#define VECTOR_INVARIANT_FORM

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

	m_dCovUa.Initialize(
		m_nHorizontalOrder,
		m_nHorizontalOrder);

	m_dCovUb.Initialize(
		m_nHorizontalOrder,
		m_nHorizontalOrder);

	m_dEnergy.Initialize(
		m_nHorizontalOrder,
		m_nHorizontalOrder);

	// Column storage
	m_dZeroColumn.Initialize(nRElements + 1);

	m_dColumnPressure.Initialize(nRElements);

	m_dColumnDxPressure.Initialize(nRElements);

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

#ifdef VECTOR_INVARIANT_FORM
		const DataMatrix4D<double> & dCovMetricA =
			pPatch->GetCovMetricA();
		const DataMatrix4D<double> & dCovMetricB =
			pPatch->GetCovMetricB();
#else
		const DataMatrix3D<double> & dChristoffelA =
			pPatch->GetChristoffelA();
		const DataMatrix3D<double> & dChristoffelB =
			pPatch->GetChristoffelB();
#endif

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
				m_dAlphaFlux[i][j] =
					dJacobian[k][iA][iB]
					* (dataInitialNode[HIx][k][iA][iB] - dTopography[iA][iB])
					* dataInitialNode[UIx][k][iA][iB];

				m_dBetaFlux[i][j] =
					dJacobian[k][iA][iB]
					* (dataInitialNode[HIx][k][iA][iB] - dTopography[iA][iB])
					* dataInitialNode[VIx][k][iA][iB];

#ifdef VECTOR_INVARIANT_FORM
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
#else
				// Pointwise pressure
				m_dPressure[i][j] =
					phys.GetG() * dataInitialNode[HIx][k][iA][iB];
#endif
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

#ifdef VECTOR_INVARIANT_FORM
				// Derivatives of the covariant velocity field
				double dCovDaUb = 0.0;
				double dCovDbUa = 0.0;
#endif

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

#ifdef VECTOR_INVARIANT_FORM
					// Derivative of covariant beta velocity wrt alpha
					dCovDaUb +=
						m_dCovUb[s][j]
						* dDxBasis1D[s][i];
#endif				
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

#ifdef VECTOR_INVARIANT_FORM
					// Derivative of alpha velocity wrt beta
					dCovDbUa +=
						m_dCovUa[i][s]
						* dDxBasis1D[s][j];
#endif
				}

				// Scale derivatives
				dDaUa /= dElementDeltaA;
				dDaUb /= dElementDeltaA;
				dDaP  /= dElementDeltaA;

				dDbUa /= dElementDeltaB;
				dDbUb /= dElementDeltaB;
				dDbP  /= dElementDeltaB;

#ifdef VECTOR_INVARIANT_FORM
				dCovDaUb /= dElementDeltaA;
				dCovDbUa /= dElementDeltaB;
#endif

				dLocalUpdateHa /= (dJacobian[k][iA][iB] * dElementDeltaA);
				dLocalUpdateHb /= (dJacobian[k][iA][iB] * dElementDeltaB);

				// Local update for momentum
				double dLocalUpdateUa = 0.0;
				double dLocalUpdateUb = 0.0;

#ifdef VECTOR_INVARIANT_FORM
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

#else
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

				// Advection
				dLocalUpdateUa -= dUa * dDaUa + dUb * dDbUa;
				dLocalUpdateUb -= dUa * dDaUb + dUb * dDbUb;
#endif

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

		const DataMatrix<double> & dJacobian2D =
			pPatch->GetJacobian2D();
		const DataMatrix3D<double> & dJacobian =
			pPatch->GetJacobian();
		const DataMatrix4D<double> & dContraMetricA =
			pPatch->GetContraMetricA();
		const DataMatrix4D<double> & dContraMetricB =
			pPatch->GetContraMetricB();

#ifdef VECTOR_INVARIANT_FORM
		const DataMatrix4D<double> & dCovMetricA =
			pPatch->GetCovMetricA();
		const DataMatrix4D<double> & dCovMetricB =
			pPatch->GetCovMetricB();
#else
		const DataMatrix3D<double> & dChristoffelA =
			pPatch->GetChristoffelA();
		const DataMatrix3D<double> & dChristoffelB =
			pPatch->GetChristoffelB();
#endif

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

			// Differentiate column pressure
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

#ifdef VECTOR_INVARIANT_FORM
				// Add kinetic energy to pressure
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
				m_dEnergy[i][j] =
					0.5 * (m_dCovUa[i][j] * dUa + m_dCovUb[i][j] * dUb);
#endif
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

#ifdef VECTOR_INVARIANT_FORM
				// Derivatives of the covariant velocity field
				double dCovDaUb = 0.0;
				double dCovDbUa = 0.0;

				// Derivative of the kinetic energy
				double dDaKE = 0.0;
				double dDbKE = 0.0;
#endif

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

#ifdef VECTOR_INVARIANT_FORM
					// Derivative of covariant beta velocity wrt alpha
					dCovDaUb +=
						m_dCovUb[s][j]
						* dDxBasis1D[s][i];

					// Derivative of kinetic energy with respect to alpha
					dDaKE +=
						m_dEnergy[s][j]
						* dDxBasis1D[s][i];
#endif				
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

#ifdef VECTOR_INVARIANT_FORM
					// Derivative of alpha velocity wrt beta
					dCovDbUa +=
						m_dCovUa[i][s]
						* dDxBasis1D[s][j];

					// Derivative of kinetic energy with respect to beta
					dDbKE +=
						m_dEnergy[i][s]
						* dDxBasis1D[s][j];

#endif
				}

				// Scale derivatives
				dDaUa /= dElementDeltaA;
				dDaUb /= dElementDeltaA;
				dDaP  /= dElementDeltaA;

				dDbUa /= dElementDeltaB;
				dDbUb /= dElementDeltaB;
				dDbP  /= dElementDeltaB;

#ifdef VECTOR_INVARIANT_FORM
				dCovDaUb /= dElementDeltaA;
				dDaKE    /= dElementDeltaA;

				dCovDbUa /= dElementDeltaB;
				dDbKE    /= dElementDeltaB;
#endif

#pragma message "Reference state?"
				dDxP  = dataDxPressure[k][iA][iB];

				// Store local horizontal pressure derivatives
				dataDaPressure[k][iA][iB] = dDaP;
				dataDbPressure[k][iA][iB] = dDbP;

				// Pointwise horizontal momentum update
				double dLocalUpdateUa = 0.0;
				double dLocalUpdateUb = 0.0;

#ifdef VECTOR_INVARIANT_FORM
				// Relative vorticity
				double dZeta = (dCovDaUb - dCovDbUa);

				// Rotational terms
				dLocalUpdateUa -=
					(dZeta + dCoriolisF[iA][iB] * dJacobian2D[iA][iB]) * (
						+ dContraMetricA[k][iA][iB][1] * dUa
						- dContraMetricA[k][iA][iB][0] * dUb);

				dLocalUpdateUb -=
					(dZeta + dCoriolisF[iA][iB] * dJacobian2D[iA][iB]) * (
						+ dContraMetricB[k][iA][iB][1] * dUa
						- dContraMetricB[k][iA][iB][0] * dUb);

				// Energy derivatives
				dLocalUpdateUa -=
						( dContraMetricA[k][iA][iB][0] * dDaKE
						+ dContraMetricA[k][iA][iB][1] * dDbKE);

				dLocalUpdateUb -=
						( dContraMetricB[k][iA][iB][0] * dDaKE
						+ dContraMetricB[k][iA][iB][1] * dDbKE);

#else
				// Advection of momentum
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

				double dCorForce = 
					dCoriolisF[iA][iB]
					* dJacobian2D[iA][iB]
					* ( + dContraMetricB[k][iA][iB][1] * dUa
						- dContraMetricB[k][iA][iB][0] * dUb);
#endif
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
				
				double dPresForce = 
						( dContraMetricB[k][iA][iB][0] * dDaP
						+ dContraMetricB[k][iA][iB][1] * dDbP
						+ dContraMetricB[k][iA][iB][2] * dDxP)
							* dataInitialNode[TIx][k][iA][iB];
				/*
				// OUTPUT THE CORIOLIS AND PRESSURE UPDATES AT THEIR LOCATIONS
				std::cout << std::fixed;
				std::cout << std::setprecision(12) << dLatitude[iA][iB] 
						  << " " << dLongitude[iA][iB] 
  						  << " " << dZLevels[k][iA][iB]
        				  << " " << dCorForce << " " << dPresForce << "\n";
				*/
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

	// Get a copy of the GLL grid
	GridGLL * pGrid = dynamic_cast<GridGLL*>(m_model.GetGrid());

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

				dataUpdate[UIx][k][iA][iB] -= dDeltaT * dInvJacobian * (
					+ dLocalNuDiv * (
						+ dAlphaDivTermA
						+ dAlphaDivTermB)
					- dLocalNuVort * dAlphaCurlTerm);

				dataUpdate[VIx][k][iA][iB] -= dDeltaT * dInvJacobian * (
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

