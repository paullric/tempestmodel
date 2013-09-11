///////////////////////////////////////////////////////////////////////////////
///
///	\file    HorizontalDynamicsFEM.cpp
///	\author  Paul Ullrich
///	\version June 18, 2013
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
#include "PolynomialInterp.h"
#include "CubedSphereTrans.h"

#include "GridCSGLL.h"

//#define DIFFERENTIAL_FORM

///////////////////////////////////////////////////////////////////////////////

HorizontalDynamicsFEM::HorizontalDynamicsFEM(
	Model & model,
	int nHorizontalOrder,
	bool fUseHyperdiffusion
) :
	HorizontalDynamics(model),
	m_nHorizontalOrder(nHorizontalOrder),
	m_fUseHyperdiffusion(fUseHyperdiffusion),
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

	GaussLobattoQuadrature::GetPoints(m_nHorizontalOrder, dG, dW);

	// Store the nodal weights in the reference element
	m_dGLLWeight = dW;

	// Derivatives of the 1D basis functions at each point on the reference
	// element [0, 1]
#pragma message "Pull this information from the Grid"
	m_dDxBasis1D.Initialize(m_nHorizontalOrder, m_nHorizontalOrder);
	m_dStiffness1D.Initialize(m_nHorizontalOrder, m_nHorizontalOrder);

	DataVector<double> dCoeffs;
	dCoeffs.Initialize(m_nHorizontalOrder);

	for (int i = 0; i < m_nHorizontalOrder; i++) {
		PolynomialInterp::DiffLagrangianPolynomialCoeffs(
			m_nHorizontalOrder, dG, dCoeffs, dG[i]);

		for (int m = 0; m < m_nHorizontalOrder; m++) {
			m_dDxBasis1D[m][i] = 2.0 * dCoeffs[m];

			m_dStiffness1D[m][i] = m_dDxBasis1D[m][i] * dW[i] / dW[m];
		}
	}

	// Generate hyperdiffusion matrix
	GenerateHyperdiffusionMatrix();

	m_dGradient.Initialize(
		2,
		m_nHorizontalOrder,
		m_nHorizontalOrder);
}

///////////////////////////////////////////////////////////////////////////////

void HorizontalDynamicsFEM::GenerateHyperdiffusionMatrix() {
/*
	// Initialize the matrix
	m_dHyperdiffusion.Initialize(
		m_nHorizontalOrder * m_nHorizontalOrder,
		m_nHorizontalOrder * m_nHorizontalOrder);

	// Get quadrature points for Gauss-Lobatto quadrature
	DataVector<double> dGX;
	DataVector<double> dWX;

	DataVector<double> dGY;
	DataVector<double> dWY;

	GaussLobattoQuadrature::GetPoints(m_nHorizontalOrder, dGX, dWX);
	GaussLobattoQuadrature::GetPoints(m_nHorizontalOrder, dGY, dWY);

	// Derivative matrix
	// Element [i,j] stores the derivative of basis function j at point i
	DataMatrix<double> dDxBasis;
	DataMatrix<double> dDyBasis;

	dDxBasis.Initialize(m_nHorizontalOrder, m_nHorizontalOrder);
	dDyBasis.Initialize(m_nHorizontalOrder, m_nHorizontalOrder);

	for (int i = 0; i < m_nHorizontalOrder; i++) {
		PolynomialInterp::DiffLagrangianPolynomialCoeffs(
			m_nHorizontalOrder, dGX, &(dDxBasis[i][0]), dGX[i]);
		PolynomialInterp::DiffLagrangianPolynomialCoeffs(
			m_nHorizontalOrder, dGY, &(dDyBasis[i][0]), dGY[i]);
	}

	// Loop through each element of the diffusion matrix and integrate
	for (int n = 0; n < m_nHorizontalOrder; n++) {
	for (int m = 0; m < m_nHorizontalOrder; m++) {
		int nbasis = m * m_nHorizontalOrder + n;

		for (int s = 0; s < m_nHorizontalOrder; s++) {
		for (int t = 0; t < m_nHorizontalOrder; t++) {
			int mbasis = t * m_nHorizontalOrder + s;

			if (m == t) {
				for (int p = 0; p < m_nHorizontalOrder; p++) {
					m_dHyperdiffusion[nbasis][mbasis] +=
						dWX[p] * dWY[m] * dDxBasis[p][n] * dDxBasis[p][s];
				}
			}
			if (n == s) {
				for (int q = 0; q < m_nHorizontalOrder; q++) {
					m_dHyperdiffusion[nbasis][mbasis] +=
						dWX[n] * dWY[q] * dDyBasis[q][m] * dDyBasis[q][t];
				}
			}
		}
		}
	}
	}

	// Multiply by inverse lumped mass matrix
	for (int n = 0; n < m_nHorizontalOrder; n++) {
	for (int m = 0; m < m_nHorizontalOrder; m++) {
		int nbasis = m * m_nHorizontalOrder + n;

		for (int s = 0; s < m_nHorizontalOrder; s++) {
		for (int t = 0; t < m_nHorizontalOrder; t++) {
			int mbasis = t * m_nHorizontalOrder + s;

			m_dHyperdiffusion[nbasis][mbasis] *= 1.0 / (dWX[n] * dWY[m]);
		}
		}
	}
	}
*/
}

///////////////////////////////////////////////////////////////////////////////

void HorizontalDynamicsFEM::StepShallowWater(
	int iDataInitial,
	int iDataUpdate,
	double dTime,
	double dDeltaT
) {

	// Get a copy of the grid
	Grid * pGrid = m_model.GetGrid();

	// Physical constants
	const PhysicalConstants & phys = m_model.GetPhysicalConstants();

	// Get indices of variables to update
	const int UIx = pGrid->GetVarIndex(0);
	const int VIx = pGrid->GetVarIndex(1);
	const int HIx = pGrid->GetVarIndex(2);

	// Perform local update
	for (int n = 0; n < pGrid->GetActivePatchCount(); n++) {
		GridPatch * pPatch = pGrid->GetActivePatch(n);

		const PatchBox & box = pPatch->GetPatchBox();

		const DataMatrix3D<double> & dJacobian =
			pPatch->GetJacobian();
		const DataMatrix4D<double> & dContraMetricA =
			pPatch->GetContraMetricA();
		const DataMatrix4D<double> & dContraMetricB =
			pPatch->GetContraMetricB();
		const DataMatrix4D<double> & dChristoffelA =
			pPatch->GetChristoffelA();
		const DataMatrix4D<double> & dChristoffelB =
			pPatch->GetChristoffelB();
		const DataMatrix<double> & dLatitude =
			pPatch->GetLatitude();

		// Data
		const GridData4D & dataInitialNode =
			pPatch->GetDataState(iDataInitial, DataLocation_Node);

		GridData4D & dataUpdateNode =
			pPatch->GetDataState(iDataUpdate, DataLocation_Node);

		// Element grid spacing
		double dElementDeltaA =
			  box.GetAEdge(box.GetHaloElements() + m_nHorizontalOrder)
			- box.GetAEdge(box.GetHaloElements());

		double dElementDeltaB =
			  box.GetBEdge(box.GetHaloElements() + m_nHorizontalOrder)
			- box.GetBEdge(box.GetHaloElements());

		// Time over grid spacing ratio
		double dCourantA = dDeltaT / dElementDeltaA;
		double dCourantB = dDeltaT / dElementDeltaB;

		// Calculate pointwise fluxes
		int nAElements =
			pPatch->GetPatchBox().GetAInteriorWidth() / m_nHorizontalOrder;
		int nBElements =
			pPatch->GetPatchBox().GetBInteriorWidth() / m_nHorizontalOrder;

		if ((pPatch->GetPatchBox().GetAInteriorWidth() % m_nHorizontalOrder) != 0) {
			_EXCEPTIONT("Logic Error: Invalid PatchBox alpha spacing");
		}
		if ((pPatch->GetPatchBox().GetBInteriorWidth() % m_nHorizontalOrder) != 0) {
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
					dLocalUpdateHa -= m_dAlphaFlux[s][j] * m_dDxBasis1D[s][i];
#else
					// Update density: Variational formulation
					dLocalUpdateHa += m_dAlphaFlux[s][j] * m_dStiffness1D[i][s];
#endif
					// Derivative of alpha velocity with respect to alpha
					dDaUa +=
						dataInitialNode[UIx][k][iAElement+s][iB]
						* m_dDxBasis1D[s][i];

					// Derivative of beta velocity with respect to alpha
					dDaUb +=
						dataInitialNode[VIx][k][iAElement+s][iB]
						* m_dDxBasis1D[s][i];

					// Derivative of pressure with respect to alpha
					dDaP +=
						m_dPressure[s][j]
						* m_dDxBasis1D[s][i];
				}

				double dLocalUpdateHb = 0.0;
				for (int s = 0; s < m_nHorizontalOrder; s++) {
#ifdef DIFFERENTIAL_FORM
					// Update density: Differential formulation
					dLocalUpdateHb -= m_dBetaFlux[i][s] * m_dDxBasis1D[s][j];
#else
					// Update density: Variational formulation
					dLocalUpdateHb += m_dBetaFlux[i][s] * m_dStiffness1D[j][s];
#endif
					// Derivative of alpha velocity with respect to beta
					dDbUa +=
						dataInitialNode[UIx][k][iA][iBElement+s]
						* m_dDxBasis1D[s][j];

					// Derivative of beta velocity with respect to beta
					dDbUb +=
						dataInitialNode[VIx][k][iA][iBElement+s]
						* m_dDxBasis1D[s][j];

					// Derivative of pressure with respect to beta
					dDbP +=
						m_dPressure[i][s]
						* m_dDxBasis1D[s][j];
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
						+ dChristoffelA[k][iA][iB][0] * dUa * dUa
						+ dChristoffelA[k][iA][iB][1] * dUa * dUb
						+ dChristoffelA[k][iA][iB][2] * dUb * dUb;

				dLocalUpdateUb -=
						+ dChristoffelB[k][iA][iB][0] * dUa * dUa
						+ dChristoffelB[k][iA][iB][1] * dUa * dUb
						+ dChristoffelB[k][iA][iB][2] * dUb * dUb;

				// Pressure derivatives
				dLocalUpdateUa -=
						+ dContraMetricA[k][iA][iB][0] * dDaP
						+ dContraMetricA[k][iA][iB][1] * dDbP;

				dLocalUpdateUb -=
						+ dContraMetricB[k][iA][iB][0] * dDaP
						+ dContraMetricB[k][iA][iB][1] * dDbP;

				// Coriolis forces
				double dF = 2.0 * phys.GetOmega() * sin(dLatitude[iA][iB]);

				dLocalUpdateUa -=
					dF * dJacobian[k][iA][iB] * (
						+ dContraMetricA[k][iA][iB][1] * dUa
						- dContraMetricA[k][iA][iB][0] * dUb);

				dLocalUpdateUb -=
					dF * dJacobian[k][iA][iB] * (
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

void HorizontalDynamicsFEM::StepNonhydrostaticPrimitive(
	int iDataInitial,
	int iDataUpdate,
	double dTime,
	double dDeltaT
) {
	// Get a copy of the grid
	Grid * pGrid = m_model.GetGrid();

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
		GridPatch * pPatch = pGrid->GetActivePatch(n);

		const PatchBox & box = pPatch->GetPatchBox();

		const DataMatrix3D<double> & dJacobian =
			pPatch->GetJacobian();
		const DataMatrix4D<double> & dContraMetricA =
			pPatch->GetContraMetricA();
		const DataMatrix4D<double> & dContraMetricB =
			pPatch->GetContraMetricB();
		const DataMatrix4D<double> & dContraMetricXi =
			pPatch->GetContraMetricXi();
		const DataMatrix4D<double> & dChristoffelA =
			pPatch->GetChristoffelA();
		const DataMatrix4D<double> & dChristoffelB =
			pPatch->GetChristoffelB();
		const DataMatrix4D<double> & dChristoffelXi =
			pPatch->GetChristoffelXi();
		const DataMatrix<double> & dLatitude =
			pPatch->GetLatitude();

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

		// Element grid spacing
		double dElementDeltaA =
			  box.GetAEdge(box.GetHaloElements() + m_nHorizontalOrder)
			- box.GetAEdge(box.GetHaloElements());

		double dElementDeltaB =
			  box.GetBEdge(box.GetHaloElements() + m_nHorizontalOrder)
			- box.GetBEdge(box.GetHaloElements());

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
				double dUa = dataInitialNode[UIx][k][iA][iB];
				double dUb = dataInitialNode[VIx][k][iA][iB];

				// Calculate derivatives in the alpha direction
				double dLocalUpdateRhoA = 0.0;
				for (int s = 0; s < m_nHorizontalOrder; s++) {
#ifdef DIFFERENTIAL_FORM
					// Update density: Differential formulation
					dLocalUpdateRhoA -=
						m_dAlphaFlux[s][j] * m_dDxBasis1D[s][i];
#else
					// Update density: Variational formulation
					dLocalUpdateRhoA +=
						m_dAlphaFlux[s][j] * m_dStiffness1D[i][s];
#endif
					// Derivative of alpha velocity with respect to alpha
					dDaUa +=
						dataInitialNode[UIx][k][iAElement+s][iB]
						* m_dDxBasis1D[s][i];

					// Derivative of beta velocity with respect to alpha
					dDaUb +=
						dataInitialNode[VIx][k][iAElement+s][iB]
						* m_dDxBasis1D[s][i];

					// Derivative of pressure with respect to alpha
					dDaP +=
						m_dPressure[s][j]
						* m_dDxBasis1D[s][i];
				}

				// Calculate derivatives in the beta direction
				double dLocalUpdateRhoB = 0.0;
				for (int s = 0; s < m_nHorizontalOrder; s++) {
#ifdef DIFFERENTIAL_FORM
					// Update density: Differential formulation
					dLocalUpdateRhoB -=
						m_dBetaFlux[i][s] * m_dDxBasis1D[s][j];
#else
					// Update density: Variational formulation
					dLocalUpdateRhoB +=
						m_dBetaFlux[i][s] * m_dStiffness1D[j][s];
#endif
					// Derivative of alpha velocity with respect to beta
					dDbUa +=
						dataInitialNode[UIx][k][iA][iBElement+s]
						* m_dDxBasis1D[s][j];

					// Derivative of beta velocity with respect to beta
					dDbUb +=
						dataInitialNode[VIx][k][iA][iBElement+s]
						* m_dDxBasis1D[s][j];

					// Derivative of pressure with respect to beta
					dDbP +=
						m_dPressure[i][s]
						* m_dDxBasis1D[s][j];
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
						+ dChristoffelA[k][iA][iB][0] * dUa * dUa
						+ dChristoffelA[k][iA][iB][1] * dUa * dUb
						+ dChristoffelA[k][iA][iB][2] * dUb * dUb;

				dLocalUpdateUb -=
						+ dChristoffelB[k][iA][iB][0] * dUa * dUa
						+ dChristoffelB[k][iA][iB][1] * dUa * dUb
						+ dChristoffelB[k][iA][iB][2] * dUb * dUb;

				// Pressure derivatives
				dLocalUpdateUa -=
						( dContraMetricA[k][iA][iB][0] * dDaP
						+ dContraMetricA[k][iA][iB][1] * dDbP)
							/ dataInitialNode[RIx][k][iA][iB];

				dLocalUpdateUb -=
						( dContraMetricB[k][iA][iB][0] * dDaP
						+ dContraMetricB[k][iA][iB][1] * dDbP)
							/ dataInitialNode[RIx][k][iA][iB];

				// Coriolis forces
				double dF = 2.0 * phys.GetOmega() * sin(dLatitude[iA][iB]);

				dLocalUpdateUa -=
					dF * dJacobian[k][iA][iB] * (
						+ dContraMetricA[k][iA][iB][1] * dUa
						- dContraMetricA[k][iA][iB][0] * dUb);

				dLocalUpdateUb -=
					dF * dJacobian[k][iA][iB] * (
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
					double dDaUx = 0.0;
					double dDbUx = 0.0;

					for (int s = 0; s < m_nHorizontalOrder; s++) {
						// Derivative of xi velocity with respect to alpha
						dDaUx +=
							dataInitialNode[WIx][k][iAElement+s][iB]
							* m_dDxBasis1D[s][i];

						// Derivative of xi velocity with respect to beta
						dDbUx +=
							dataInitialNode[WIx][k][iA][iBElement+s]
							* m_dDxBasis1D[s][j];
					}

					// Scale derivatives
					dDaUx /= dElementDeltaA;
					dDbUx /= dElementDeltaA;

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
							dataInitialNode[TIx][k][iAElement+s][iB]
							* m_dDxBasis1D[s][i];

						// Derivative of theta with respect to beta
						dDbTheta +=
							dataInitialNode[TIx][k][iA][iBElement+s]
							* m_dDxBasis1D[s][j];
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

					double dUaREdge = dataInitialREdge[UIx][k][iA][iB];
					double dUbREdge = dataInitialREdge[VIx][k][iA][iB];

					for (int s = 0; s < m_nHorizontalOrder; s++) {
						// Derivative of xi velocity with respect to alpha
						dDaUx +=
							dataInitialREdge[WIx][k][iAElement+s][iB]
							* m_dDxBasis1D[s][i];

						// Derivative of xi velocity with respect to beta
						dDbUx +=
							dataInitialREdge[WIx][k][iA][iBElement+s]
							* m_dDxBasis1D[s][j];
					}

					// Scale derivatives
					dDaUx /= dElementDeltaA;
					dDbUx /= dElementDeltaA;

					// Update vertical velocity
					dataUpdateNode[WIx][k][iA][iB] -=
						dDeltaT * (dUaREdge * dDaUx + dUbREdge * dDbUx);
				}

				// Update the potential temperature (on model interfaces)
				if (pGrid->GetVarLocation(TIx) == DataLocation_REdge) {
					double dDaTheta = 0.0;
					double dDbTheta = 0.0;

					double dUaREdge = dataInitialREdge[UIx][k][iA][iB];
					double dUbREdge = dataInitialREdge[VIx][k][iA][iB];

					for (int s = 0; s < m_nHorizontalOrder; s++) {
						// Derivative of theta with respect to alpha
						dDaTheta +=
							dataInitialREdge[TIx][k][iAElement+s][iB]
							* m_dDxBasis1D[s][i];

						// Derivative of theta with respect to beta
						dDbTheta +=
							dataInitialREdge[TIx][k][iA][iBElement+s]
							* m_dDxBasis1D[s][j];
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

	// Invalid EquationSet
	} else {
		_EXCEPTIONT("Invalid EquationSet");
	}

	// Apply Direct Stiffness Summation (DSS) procedure
	GridCSGLL * pGridCSGLL = dynamic_cast<GridCSGLL*>(m_model.GetGrid());
	pGridCSGLL->ApplyDSS(iDataUpdate);

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
	int iC,
	double dDeltaT,
	bool fUseHyperdiffusionCoeff
) {
	// Get a copy of the grid
	Grid * pGrid = m_model.GetGrid();

	// Physical constants
	const PhysicalConstants & phys = m_model.GetPhysicalConstants();

	// Perform local update
	for (int n = 0; n < pGrid->GetActivePatchCount(); n++) {
		GridPatch * pPatch = pGrid->GetActivePatch(n);

		const PatchBox & box = pPatch->GetPatchBox();

		const DataMatrix3D<double> & dJacobian =
			pPatch->GetJacobian();
		const DataMatrix4D<double> & dContraMetricA =
			pPatch->GetContraMetricA();
		const DataMatrix4D<double> & dContraMetricB =
			pPatch->GetContraMetricB();

		GridData4D & dataInitial =
			pPatch->GetDataState(iDataInitial, DataLocation_Node);

		GridData4D & dataUpdate =
			pPatch->GetDataState(iDataUpdate, DataLocation_Node);

		// Element grid spacing
		double dElementDeltaA =
			  box.GetAEdge(box.GetHaloElements() + m_nHorizontalOrder)
			- box.GetAEdge(box.GetHaloElements());

		double dElementDeltaB =
			  box.GetBEdge(box.GetHaloElements() + m_nHorizontalOrder)
			- box.GetBEdge(box.GetHaloElements());

		// Calculate local hyperdiffusion coefficient
		double dCoeff;

		if (fUseHyperdiffusionCoeff) {
			const double dElementRefDeltaA = 0.5 * M_PI / 30.0;

			dCoeff =
				- dDeltaT * m_dNuScalar
					* pow(dElementDeltaA / dElementRefDeltaA, 3.2);
		} else {
			dCoeff = 1.0;
		}

		// Number of finite elements
		int nAElements =
			box.GetAInteriorWidth() / m_nHorizontalOrder;
		int nBElements =
			box.GetBInteriorWidth() / m_nHorizontalOrder;

		// Loop over all finite elements
		for (int k = 0; k < dataInitial.GetRElements(); k++) {
		for (int a = 0; a < nAElements; a++) {
		for (int b = 0; b < nBElements; b++) {

			int iAElement = a * m_nHorizontalOrder + box.GetHaloElements();
			int iBElement = b * m_nHorizontalOrder + box.GetHaloElements();

			// Calculate the gradient within each element
			m_dGradient.Zero();

			// Pointwise update of horizontal velocities
			for (int i = 0; i < m_nHorizontalOrder; i++) {
			for (int j = 0; j < m_nHorizontalOrder; j++) {

				int iA = iAElement + i;
				int iB = iBElement + j;

				// Calculate pointwise derivatives
				double dDaPhi = 0.0;
				double dDbPhi = 0.0;

				for (int s = 0; s < m_nHorizontalOrder; s++) {
					dDaPhi +=
						dataInitial[iC][k][iAElement+s][iB]
						* m_dDxBasis1D[s][i];

					dDbPhi +=
						dataInitial[iC][k][iA][iBElement+s]
						* m_dDxBasis1D[s][j];
				}

				dDaPhi /= dElementDeltaA;
				dDbPhi /= dElementDeltaB;

				//printf("%1.10e %1.10e\n", dDaPhi, dDbPhi);

				// Compute gradient terms
				m_dGradient[0][i][j] = dJacobian[k][iA][iB] * (
					+ dContraMetricA[k][iA][iB][0] * dDaPhi
					+ dContraMetricA[k][iA][iB][1] * dDbPhi);

				m_dGradient[1][i][j] = dJacobian[k][iA][iB] * (
					+ dContraMetricB[k][iA][iB][0] * dDaPhi
					+ dContraMetricB[k][iA][iB][1] * dDbPhi);
			}
			}

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
/*
				if ((n == 0) && (a == 0) && (b == 0)) {
					printf("%1.10e %1.10e\n", dCoeff, (dDaGradient + dDbGradient) / dJacobian[k][iA][iB]);
				}
*/
			}
			}
/*
			if (dCoeff != 1.0) {
			printf("%1.10e %1.10e\n", dElementDeltaA, dElementDeltaB);
			for (int i = 0; i < m_nHorizontalOrder; i++) {
			for (int j = 0; j < m_nHorizontalOrder; j++) {
				int iA = a * m_nHorizontalOrder + i + box.GetHaloElements();
				int iB = b * m_nHorizontalOrder + j + box.GetHaloElements();

				//printf("%1.10e %1.10e\n", box.GetANode(iA), box.GetBNode(iB));
				printf("[%1.10e, %1.10e],\n", dataInitial[iC][k][iA][iB], dataUpdate[iC][k][iA][iB]);
			}
			}
			_EXCEPTION();
			}
*/
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

	// Get a copy of the grid
	Grid * pGrid = m_model.GetGrid();

	// Physical constants
	const PhysicalConstants & phys = m_model.GetPhysicalConstants();

	// Loop over all patches
	for (int n = 0; n < pGrid->GetActivePatchCount(); n++) {
		GridPatch * pPatch = pGrid->GetActivePatch(n);

		const PatchBox & box = pPatch->GetPatchBox();

		const DataMatrix3D<double> & dJacobian =
			pPatch->GetJacobian();
		const DataMatrix4D<double> & dContraMetricA =
			pPatch->GetContraMetricA();
		const DataMatrix4D<double> & dContraMetricB =
			pPatch->GetContraMetricB();

		GridData4D & dataInitial =
			pPatch->GetDataState(iDataInitial, DataLocation_Node);

		GridData4D & dataUpdate =
			pPatch->GetDataState(iDataUpdate, DataLocation_Node);

		// Element grid spacing
		double dElementDeltaA =
			  box.GetAEdge(box.GetHaloElements() + m_nHorizontalOrder)
			- box.GetAEdge(box.GetHaloElements());

		double dElementDeltaB =
			  box.GetBEdge(box.GetHaloElements() + m_nHorizontalOrder)
			- box.GetBEdge(box.GetHaloElements());

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
			dLocalNuDiv =
				dLocalNuDiv  * pow(dElementDeltaA / (0.5 * M_PI / 30.0), 3.2);
			dLocalNuVort =
				dLocalNuVort * pow(dElementDeltaA / (0.5 * M_PI / 30.0), 3.2);
		}

		// Number of finite elements
		int nAElements =
			box.GetAInteriorWidth() / m_nHorizontalOrder;
		int nBElements =
			box.GetBInteriorWidth() / m_nHorizontalOrder;

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
						dJacobian[k][iAElement+s][iB]
						* m_dStiffness1D[i][s]
						* dataDiv[k][iAElement+s][iB];

					double dBetaDiv =
						dJacobian[k][iA][iBElement+s]
						* m_dStiffness1D[j][s]
						* dataDiv[k][iA][iBElement+s];

					dAlphaDivTermA +=
						dAlphaDiv * dContraMetricA[k][iAElement+s][iB][0];
					dAlphaDivTermB +=
						dBetaDiv  * dContraMetricB[k][iA][iBElement+s][0];

					dAlphaCurlTerm +=
						m_dStiffness1D[j][s]
						* dataCurl[k][iA][iBElement+s];

					dBetaDivTermA +=
						dAlphaDiv * dContraMetricA[k][iAElement+s][iB][1];
					dBetaDivTermB +=
						dBetaDiv  * dContraMetricB[k][iA][iBElement+s][1];

					dBetaCurlTerm +=
						m_dStiffness1D[i][s]
						* dataCurl[k][iAElement+s][iB];
				}

				dAlphaDivTermA /= dElementDeltaA;
				dAlphaDivTermB /= dElementDeltaB;
				dAlphaCurlTerm /= dElementDeltaB;

				dBetaDivTermA /= dElementDeltaA;
				dBetaDivTermB /= dElementDeltaB;
				dBetaCurlTerm /= dElementDeltaA;

				// Apply update
				double dInvJacobian = 1.0 / dJacobian[k][iA][iB];

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
	double dTime,
	double dDeltaT
) {
	// Only proceed if hyperdiffusion is applied
	if (!m_fUseHyperdiffusion) {
		return;
	}

#pragma message "Altering the horizontal velocities should also modify the vertical velocities"

	// Apply Direct Stiffness Summation (DSS) procedure
	GridCSGLL * pGridCSGLL = dynamic_cast<GridCSGLL*>(m_model.GetGrid());

	// Zero the target state data
	pGridCSGLL->ZeroData(1, DataType_State);

	// Apply vector Laplacian (first application)
	ApplyVectorHyperdiffusion(0, 1, 1.0, 1.0, 1.0, false);

	pGridCSGLL->ApplyDSS(1);

	//pGridCSGLL->CopyData(1, 0, DataType_State);

	// Apply vector Laplacian (second application)
	ApplyVectorHyperdiffusion(1, 0, -dDeltaT, m_dNuDiv, m_dNuVort, true);

	pGridCSGLL->ApplyDSS(0);

/*
	// Variable indices
	const int UIx = 0;
	const int VIx = 1;
	const int HIx = 2;

	// Perform local update
	for (int n = 0; n < pGridCSGLL->GetActivePatchCount(); n++) {
		GridPatch * pPatch = pGridCSGLL->GetActivePatch(n);

		GridData4D & dataInitial =
			pPatch->GetDataState(1, DataLocation_Node);

		dataInitial.Zero();
	}

	// Apply the fourth-order hyperdiffusion operator
	ApplyScalarHyperdiffusion(iDataInitial, 1, UIx, dDeltaT, false);
	ApplyScalarHyperdiffusion(iDataInitial, 1, VIx, dDeltaT, false);

	// Apply Direct Stiffness Summation (DSS) procedure
	pGridCSGLL->ApplyDSS(iDataUpdate);

	pGridCSGLL->ApplyDSS(1);

	// Apply the fourth-order hyperdiffusion operator
	ApplyScalarHyperdiffusion(1, iDataUpdate, UIx, dDeltaT, true);
	ApplyScalarHyperdiffusion(1, iDataUpdate, VIx, dDeltaT, true);

	// Apply Direct Stiffness Summation (DSS) procedure
	pGridCSGLL->ApplyDSS(iDataUpdate);
*/
}

///////////////////////////////////////////////////////////////////////////////

