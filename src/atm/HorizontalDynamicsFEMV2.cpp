///////////////////////////////////////////////////////////////////////////////
///
///	\file    HorizontalDynamicsFEMV2.cpp
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
#include "HorizontalDynamicsFEMV2.h"
#include "PhysicalConstants.h"
#include "Model.h"
#include "Grid.h"
#include "FunctionTimer.h"

#include "Announce.h"
#include "GridGLL.h"
#include "GridPatchGLL.h"

//#define DIFFERENTIAL_FORM

#ifdef DIFFERENTIAL_FORM
#pragma message "WARNING: DIFFERENTIAL_FORM will lose mass over topography"
#endif

//#define INSTEP_DIVERGENCE_DAMPING

#define FIX_ELEMENT_MASS_NONHYDRO

///////////////////////////////////////////////////////////////////////////////

HorizontalDynamicsFEMV2::HorizontalDynamicsFEMV2(
	Model & model,
	int nHorizontalOrder,
	int nHyperviscosityOrder,
	double dNuScalar,
	double dNuDiv,
	double dNuVort,
	double dInstepNuDiv
) :
	HorizontalDynamics(model),
	m_nHorizontalOrder(nHorizontalOrder),
	m_nHyperviscosityOrder(nHyperviscosityOrder),
	m_dNuScalar(dNuScalar),
	m_dNuDiv(dNuDiv),
	m_dNuVort(dNuVort),
	m_dInstepNuDiv(dInstepNuDiv)
{
}

///////////////////////////////////////////////////////////////////////////////

void HorizontalDynamicsFEMV2::Initialize() {

	// Get a copy of the GLL grid
	GridGLL * pGrid = dynamic_cast<GridGLL*>(m_model.GetGrid());
	if (pGrid == NULL) {
		_EXCEPTIONT("Grid must be of type GridGLL");
	}

	// Number of vertical levels
	int nRElements = pGrid->GetRElements();

	// Number of tracers
	int nTracerCount = m_model.GetEquationSet().GetTracers();

	// Initialize the alpha and beta mass fluxes
	m_dAlphaMassFlux.Allocate(
		m_nHorizontalOrder,
		m_nHorizontalOrder);

	m_dBetaMassFlux.Allocate(
		m_nHorizontalOrder,
		m_nHorizontalOrder);

	// Initialize the alpha and beta mass fluxes for mass fix
	m_dAlphaElMassFlux.Allocate(
		m_nHorizontalOrder,
		m_nHorizontalOrder);

	m_dBetaElMassFlux.Allocate(
		m_nHorizontalOrder,
		m_nHorizontalOrder);

	// Initialize the alpha and beta pressure fluxes
	m_dAlphaPressureFlux.Allocate(
		m_nHorizontalOrder,
		m_nHorizontalOrder);

	m_dBetaPressureFlux.Allocate(
		m_nHorizontalOrder,
		m_nHorizontalOrder);

	// Initialize tracer fluxes
	if (nTracerCount > 0) {
		m_dAlphaTracerFlux.Allocate(
			nTracerCount,
			m_nHorizontalOrder,
			m_nHorizontalOrder);

		m_dBetaTracerFlux.Allocate(
			nTracerCount,
			m_nHorizontalOrder,
			m_nHorizontalOrder);
	}

	// Auxiliary data
	m_dAuxDataNode.Allocate(
		9,
		nRElements,
		m_nHorizontalOrder,
		m_nHorizontalOrder);

	m_dAuxDataREdge.Allocate(
		9,
		nRElements+1,
		m_nHorizontalOrder,
		m_nHorizontalOrder);

	m_dDivergence.Allocate(
		nRElements,
		m_nHorizontalOrder,
		m_nHorizontalOrder);

	// Contravariant metric terms
	m_dLocalCoriolisF.Allocate(
		m_nHorizontalOrder,
		m_nHorizontalOrder);

	m_dLocalJacobian2D.Allocate(
		m_nHorizontalOrder,
		m_nHorizontalOrder);

	m_dLocalJacobian.Allocate(
		nRElements,
		m_nHorizontalOrder,
		m_nHorizontalOrder);

	m_dLocalDerivR.Allocate(
		nRElements,
		m_nHorizontalOrder,
		m_nHorizontalOrder,
		3);

	m_dLocalContraMetric.Allocate(
		nRElements,
		m_nHorizontalOrder,
		m_nHorizontalOrder,
		6);

	// Initialize buffers for derivatives of Jacobian
	m_dJGradientA.Allocate(
		m_nHorizontalOrder,
		m_nHorizontalOrder);

	m_dJGradientB.Allocate(
		m_nHorizontalOrder,
		m_nHorizontalOrder);

	// Buffer state
	m_dBufferState.Allocate(
		m_nHorizontalOrder,
		m_nHorizontalOrder);
}

///////////////////////////////////////////////////////////////////////////////

void HorizontalDynamicsFEMV2::FilterNegativeTracers(
	int iDataUpdate
) {
#ifdef POSITIVE_DEFINITE_FILTER_TRACERS
	// Get a copy of the GLL grid
	GridGLL * pGrid = dynamic_cast<GridGLL*>(m_model.GetGrid());

	// Number of vertical elements
	const int nRElements = pGrid->GetRElements();

	// Perform local update
	for (int n = 0; n < pGrid->GetActivePatchCount(); n++) {
		GridPatchGLL * pPatch =
			dynamic_cast<GridPatchGLL*>(pGrid->GetActivePatch(n));

		const PatchBox & box = pPatch->GetPatchBox();

		const DataArray3D<double> & dElementAreaNode =
			pPatch->GetElementAreaNode();

		DataArray4D<double> & dataUpdateTracer =
			pPatch->GetDataTracers(iDataUpdate);

		// Number of tracers
		const int nTracerCount = dataUpdateTracer.GetSize(0);

		// Get number of finite elements in each coordinate direction
		int nElementCountA = pPatch->GetElementCountA();
		int nElementCountB = pPatch->GetElementCountB();

		// Loop over all elements
		for (int a = 0; a < nElementCountA; a++) {
		for (int b = 0; b < nElementCountB; b++) {

			// Loop overall tracers and vertical levels
			for (int c = 0; c < nTracerCount; c++) {
			for (int k = 0; k < nRElements; k++) {

				// Compute total mass and non-negative mass
				double dTotalInitialMass = 0.0;
				double dTotalMass = 0.0;
				double dNonNegativeMass = 0.0;

				for (int i = 0; i < m_nHorizontalOrder; i++) {
				for (int j = 0; j < m_nHorizontalOrder; j++) {

					int iA = a * m_nHorizontalOrder + i + box.GetHaloElements();
					int iB = b * m_nHorizontalOrder + j + box.GetHaloElements();

					double dPointwiseMass =
						  dataUpdateTracer[c][k][iA][iB]
						* dElementAreaNode[k][iA][iB];

					dTotalMass += dPointwiseMass;

					if (dataUpdateTracer[c][k][iA][iB] >= 0.0) {
						dNonNegativeMass += dPointwiseMass;
					}
				}
				}
/*
				// Check for negative total mass
				if (dTotalMass < 1.0e-14) {
					printf("%i %i %i %1.15e %1.15e\n", n, a, b, dTotalInitialMass, dTotalMass);
					_EXCEPTIONT("Negative element mass detected in filter");
				}
*/
				// Apply scaling ratio to points with non-negative mass
				double dR = dTotalMass / dNonNegativeMass;

				for (int i = 0; i < m_nHorizontalOrder; i++) {
				for (int j = 0; j < m_nHorizontalOrder; j++) {

					int iA = a * m_nHorizontalOrder + i + box.GetHaloElements();
					int iB = b * m_nHorizontalOrder + j + box.GetHaloElements();

					if (dataUpdateTracer[c][k][iA][iB] > 0.0) {
						dataUpdateTracer[c][k][iA][iB] *= dR;
					} else {
						dataUpdateTracer[c][k][iA][iB] = 0.0;
					}
				}
				}
			}
			}
		}
		}
	}
#endif
}

///////////////////////////////////////////////////////////////////////////////

template<int n>
inline double ApplyStiffness(
	double d0,
	double d1,
	double d2,
	double d3
) { }

template<>
inline double ApplyStiffness<0>(
	double d0,
	double d1,
	double d2,
	double d3
) {
	return (
		+ 6.0 * d0
		+ 8.09016994374947451262869435595E0 * d1
		- 3.09016994374947451262869435595E0 * d2
		+ d3);
}

template<>
inline double ApplyStiffness<1>(
	double d0,
	double d1,
	double d2,
	double d3
) {
	return (
		- 1.61803398874989490252573887119E0 * d0
		+ 2.23606797749978980505147774238E0 * d2
		- 6.18033988749894902525738871191E-1 * d3);
}

template<>
inline double ApplyStiffness<2>(
	double d0,
	double d1,
	double d2,
	double d3
) {
	return (
		+ 6.18033988749894902525738871191E-1 * d0
		- 2.23606797749978980505147774238E0 * d1
		+ 1.61803398874989490252573887119E0 * d3);
}

template<>
inline double ApplyStiffness<3>(
	double d0,
	double d1,
	double d2,
	double d3
) {
	return (
		- d0
		+ 3.09016994374947451262869435595E0 * d1
		- 8.09016994374947451262869435595E0 * d2
		- 6.0 * d3);
}

/*

		m_dStiffness1D[0][0] = -6.0;
		m_dStiffness1D[1][0] =  1.61803398874989490252573887119E0;
		m_dStiffness1D[2][0] = -6.18033988749894902525738871191E-1;
		m_dStiffness1D[3][0] =  1.0;

		m_dStiffness1D[0][1] = -8.09016994374947451262869435595E0;
		m_dStiffness1D[1][1] =  0.0;
		m_dStiffness1D[2][1] =  2.23606797749978980505147774238E0;
		m_dStiffness1D[3][1] = -3.09016994374947451262869435595E0;

		m_dStiffness1D[0][2] =  3.09016994374947451262869435595E0;
		m_dStiffness1D[1][2] = -2.23606797749978980505147774238E0;
		m_dStiffness1D[2][2] =  0.0;
		m_dStiffness1D[3][2] =  8.09016994374947451262869435595E0;

		m_dStiffness1D[0][3] = -1.0;
		m_dStiffness1D[1][3] =  6.18033988749894902525738871191E-1;
		m_dStiffness1D[2][3] = -1.61803398874989490252573887119E0;
		m_dStiffness1D[3][3] =  6.0;
*/

///////////////////////////////////////////////////////////////////////////////

void HorizontalDynamicsFEMV2::StepNonhydrostaticPrimitive(
	int iDataInitial,
	int iDataUpdate,
	const Time & time,
	double dDeltaT
) {

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
	const int KIx = 4;
	const int UCrossZetaAIx = 5;
	const int UCrossZetaBIx = 6;
	const int UCrossZetaXIx = 7;
	const int ExnerIx = 8;

#if defined(FORMULATION_THETA_FLUX)
	_EXCEPTIONT("Not implemented for FORMULATION_THETA_FLUX");
#endif
#if defined(FORMULATION_PRESSURE)
	_EXCEPTIONT("Not implemented for FORMULATION_PRESSURE");
#endif

	// Start the function timer
	FunctionTimer timer("HorizontalStepNonhydrostaticPrimitive");

	// Get a copy of the GLL grid
	GridGLL * pGrid = dynamic_cast<GridGLL*>(m_model.GetGrid());

	// Check variable positions
	if ((pGrid->GetVarLocation(WIx) != DataLocation_REdge) ||
	    (pGrid->GetVarLocation(UIx) != DataLocation_Node) ||
	    (pGrid->GetVarLocation(VIx) != DataLocation_Node) ||
	    (pGrid->GetVarLocation(UIx) != DataLocation_Node) ||
	    (pGrid->GetVarLocation(VIx) != DataLocation_Node)
	) {
		_EXCEPTIONT("Only Lorenz staggering is supported");
	}

	if (m_nHorizontalOrder != 4) {
		_EXCEPTIONT("Only 4th order currently supported");
	}

	if (pGrid->HasUniformDiffusion()) {
		_EXCEPTIONT("Uniform diffusion not supported");
	}

	// Number of vertical elements
	const int nRElements = pGrid->GetRElements();

	// Physical constants
	const PhysicalConstants & phys = m_model.GetPhysicalConstants();

	// Vertical level stride in local data arrays
	const int nVerticalElementStride =
		m_nHorizontalOrder * m_nHorizontalOrder;

	// Perform local update
	for (int n = 0; n < pGrid->GetActivePatchCount(); n++) {
		GridPatchGLL * pPatch =
			dynamic_cast<GridPatchGLL*>(pGrid->GetActivePatch(n));

		const PatchBox & box = pPatch->GetPatchBox();

		const DataArray3D<double> & dElementAreaNode =
			pPatch->GetElementAreaNode();
		const DataArray2D<double> & dJacobian2D =
			pPatch->GetJacobian2D();
		const DataArray3D<double> & dJacobian =
			pPatch->GetJacobian();
		const DataArray3D<double> & dJacobianREdge =
			pPatch->GetJacobianREdge();
		const DataArray4D<double> & dContraMetricA =
			pPatch->GetContraMetricA();
		const DataArray4D<double> & dContraMetricB =
			pPatch->GetContraMetricB();
		const DataArray4D<double> & dContraMetricXi =
			pPatch->GetContraMetricXi();
		const DataArray4D<double> & dContraMetricAREdge =
			pPatch->GetContraMetricAREdge();
		const DataArray4D<double> & dContraMetricBREdge =
			pPatch->GetContraMetricBREdge();
		const DataArray4D<double> & dContraMetricXiREdge =
			pPatch->GetContraMetricXiREdge();

		const DataArray4D<double> & dDerivRNode =
			pPatch->GetDerivRNode();
		const DataArray4D<double> & dDerivRREdge =
			pPatch->GetDerivRREdge();

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

		DataArray4D<double> & dataInitialTracer =
			pPatch->GetDataTracers(iDataInitial);

		DataArray4D<double> & dataUpdateTracer =
			pPatch->GetDataTracers(iDataUpdate);

		// Number of tracers
		const int nTracerCount = dataInitialTracer.GetSize(0);

		// Spacing between vertical levels in dataInitialNode
		const int nVerticalStateStride =
			  dataInitialNode.GetSize(2)
			* dataInitialNode.GetSize(3);

		// Perform interpolations as required due to vertical staggering
		if (pGrid->GetVarsAtLocation(DataLocation_REdge) != 0) {

			// Interpolate W to model levels
			if (pGrid->GetVarLocation(WIx) == DataLocation_REdge) {
				pPatch->InterpolateREdgeToNode(WIx, iDataInitial);

				pPatch->InterpolateNodeToREdge(UIx, iDataInitial);
				pPatch->InterpolateNodeToREdge(VIx, iDataInitial);
			}

			// Interpolate Theta to model levels
			if (pGrid->GetVarLocation(PIx) == DataLocation_REdge) {
				pPatch->InterpolateREdgeToNode(PIx, iDataInitial);
			}
		}

		// Element grid spacing and derivative coefficients
		const double dElementDeltaA = pPatch->GetElementDeltaA();
		const double dElementDeltaB = pPatch->GetElementDeltaB();

		const double dInvElementDeltaA = 1.0 / dElementDeltaA;
		const double dInvElementDeltaB = 1.0 / dElementDeltaB;

		const DataArray2D<double> & dDxBasis1D = pGrid->GetDxBasis1D();
		const DataArray2D<double> & dStiffness1D = pGrid->GetStiffness1D();

		// Get number of finite elements in each coordinate direction
		const int nElementCountA = pPatch->GetElementCountA();
		const int nElementCountB = pPatch->GetElementCountB();

		// Loop over all elements
		for (int a = 0; a < nElementCountA; a++) {
		for (int b = 0; b < nElementCountB; b++) {

			const int iElementA = a * m_nHorizontalOrder + box.GetHaloElements();
			const int iElementB = b * m_nHorizontalOrder + box.GetHaloElements();

			// Store 2D Jacobian
			for (int i = 0; i < m_nHorizontalOrder; i++) {
			for (int j = 0; j < m_nHorizontalOrder; j++) {
				const int iA = iElementA + i;
				const int iB = iElementB + j;

				m_dLocalCoriolisF[i][j] = dCoriolisF[iA][iB];
				m_dLocalJacobian2D[i][j] = dJacobian2D[iA][iB];
			}
			}

			// Compute auxiliary data in element
			for (int k = 0; k < nRElements; k++) {
			for (int i = 0; i < m_nHorizontalOrder; i++) {
			for (int j = 0; j < m_nHorizontalOrder; j++) {
				const int iA = iElementA + i;
				const int iB = iElementB + j;

				// Contravariant velocities
				const double dCovUa = dataInitialNode[UIx][k][iA][iB];
				const double dCovUb = dataInitialNode[VIx][k][iA][iB];
				const double dCovUx = dataInitialNode[WIx][k][iA][iB];

				// Store metric quantities
				m_dLocalJacobian[k][i][j] = dJacobian[k][iA][iB];

				m_dLocalDerivR[k][i][j][0] = dDerivRNode[k][iA][iB][0];
				m_dLocalDerivR[k][i][j][1] = dDerivRNode[k][iA][iB][1];
				m_dLocalDerivR[k][i][j][2] = dDerivRNode[k][iA][iB][2];

				m_dLocalContraMetric[k][i][j][0] =
					dContraMetricA[k][iA][iB][0];
				m_dLocalContraMetric[k][i][j][1] =
					dContraMetricA[k][iA][iB][1];
				m_dLocalContraMetric[k][i][j][2] =
					dContraMetricA[k][iA][iB][2];
				m_dLocalContraMetric[k][i][j][3] =
					dContraMetricB[k][iA][iB][1];
				m_dLocalContraMetric[k][i][j][4] =
					dContraMetricB[k][iA][iB][2];
				m_dLocalContraMetric[k][i][j][5] =
					dContraMetricXi[k][iA][iB][2];

				// Contravariant velocities
				const double dConUa =
					  m_dLocalContraMetric[k][i][j][0] * dCovUa
					+ m_dLocalContraMetric[k][i][j][1] * dCovUb
					+ m_dLocalContraMetric[k][i][j][2] * dCovUx;

				const double dConUb =
					  m_dLocalContraMetric[k][i][j][1] * dCovUa
					+ m_dLocalContraMetric[k][i][j][3] * dCovUb
					+ m_dLocalContraMetric[k][i][j][4] * dCovUx;

				const double dConUx =
					  m_dLocalContraMetric[k][i][j][2] * dCovUa
					+ m_dLocalContraMetric[k][i][j][4] * dCovUb
					+ m_dLocalContraMetric[k][i][j][5] * dCovUx;

				m_dAuxDataNode[ConUaIx][k][i][j] = dConUa;
				m_dAuxDataNode[ConUbIx][k][i][j] = dConUb;
				m_dAuxDataNode[ConUxIx][k][i][j] = dConUx;

				// Specific kinetic energy
				m_dAuxDataNode[KIx][k][i][j] = 0.5 * (
					    dConUa * dCovUa
					  + dConUb * dCovUb
					  + dConUx * dCovUx);

#ifdef FORMULATION_RHOTHETA_P
				// Pressure
				m_dAuxDataNode[ExnerIx][k][i][j] =
					phys.PressureFromRhoTheta(
						dataInitialNode[PIx][k][iA][iB]);
#endif
#ifdef FORMULATION_RHOTHETA_PI
				// Exner pressure
				m_dAuxDataNode[ExnerIx][k][i][j] =
					phys.ExnerPressureFromRhoTheta(
						dataInitialNode[PIx][k][iA][iB]);
#endif
#if defined(FORMULATION_THETA) || defined(FORMULATION_THETA_FLUX)
				// Exner pressure
				m_dAuxDataNode[ExnerIx][k][i][j] =
					phys.ExnerPressureFromRhoTheta(
						  dataInitialNode[RIx][k][iA][iB]
						* dataInitialNode[PIx][k][iA][iB]);
#endif

				// Derivatives of the covariant velocity field
				double dCovDaUb = 0.0;
				double dCovDaUx = 0.0;
				double dCovDbUa = 0.0;
				double dCovDbUx = 0.0;

				// Vertical derivatives of the covariant velocity field			
				const double dCovDxUa =
					pGrid->DifferentiateNodeToNode(
						&(dataUpdateNode[UIx][0][iA][iB]),
						k, nVerticalStateStride);

				const double dCovDxUb =
					pGrid->DifferentiateNodeToNode(
						&(dataUpdateNode[VIx][0][iA][iB]),
						k, nVerticalStateStride);

				// Derivative needed for calculating relative vorticity
				for (int s = 0; s < m_nHorizontalOrder; s++) {

					// Derivative of covariant beta velocity wrt alpha
					dCovDaUb +=
						dataInitialNode[VIx][k][iElementA+s][iB]
						* dDxBasis1D[s][i];

					// Derivative of covariant xi velocity wrt alpha
					dCovDaUx +=
						dataInitialNode[WIx][k][iElementA+s][iB]
						* dDxBasis1D[s][i];

					// Derivative of covariant alpha velocity wrt beta
					dCovDbUa +=
						dataInitialNode[UIx][k][iA][iElementB+s]
						* dDxBasis1D[s][j];

					// Derivative of covariant xi velocity wrt beta
					dCovDbUx +=
						dataInitialNode[WIx][k][iA][iElementB+s]
						* dDxBasis1D[s][j];
				}

				dCovDaUb *= dInvElementDeltaA;
				dCovDaUx *= dInvElementDeltaA;
				dCovDbUa *= dInvElementDeltaB;
				dCovDbUx *= dInvElementDeltaB;

				// Relative vorticity (contravariant)
				const double dJZetaA = (dCovDbUx - dCovDxUb);
				const double dJZetaB = (dCovDxUa - dCovDaUx);
				const double dJZetaX = (dCovDaUb - dCovDbUa);

				// U cross Relative Vorticity (contravariant)
				m_dAuxDataNode[UCrossZetaAIx][k][i][j] =
					dConUb * dJZetaX - dConUx * dJZetaB;

				m_dAuxDataNode[UCrossZetaBIx][k][i][j] =
					dConUx * dJZetaA - dConUa * dJZetaX;

				m_dAuxDataNode[UCrossZetaXIx][k][i][j] =
					- dConUa * dCovDaUx - dConUb * dCovDbUx;
			}
			}
			}

			// Update quantities on nodes
			for (int k = 0; k < nRElements; k++) {

				for (int i = 0; i < m_nHorizontalOrder; i++) {
				for (int j = 0; j < m_nHorizontalOrder; j++) {

					const int iA = iElementA + i;
					const int iB = iElementB + j;

					// Base fluxes (area times velocity)
					const double dAlphaBaseFlux =
						m_dAuxDataNode[ConUaIx][k][i][j]
						* m_dLocalJacobian[k][i][j];

					const double dBetaBaseFlux =
						m_dAuxDataNode[ConUbIx][k][i][j]
						* m_dLocalJacobian[k][i][j];

					// Density flux
					m_dAlphaMassFlux[i][j] =
						  dAlphaBaseFlux
						* dataInitialNode[RIx][k][iA][iB];

					m_dBetaMassFlux[i][j] =
						  dBetaBaseFlux
						* dataInitialNode[RIx][k][iA][iB];

#if defined(FORMULATION_RHOTHETA_PI) \
 || defined(FORMULATION_RHOTHETA_P)
					// RhoTheta flux
					m_dAlphaPressureFlux[i][j] =
						  dAlphaBaseFlux
						* dataInitialNode[PIx][k][iA][iB];

					m_dBetaPressureFlux[i][j] =
						  dBetaBaseFlux
						* dataInitialNode[PIx][k][iA][iB];
#endif
					for (int c = 0; c < nTracerCount; c++) {
						m_dAlphaTracerFlux[c][i][j] =
							dAlphaBaseFlux
							* dataInitialTracer[c][k][iA][iB];

						m_dBetaTracerFlux[c][i][j] =
							dBetaBaseFlux
							* dataInitialTracer[c][k][iA][iB];
					}
				}
				}

				for (int i = 0; i < m_nHorizontalOrder; i++) {
				for (int j = 0; j < m_nHorizontalOrder; j++) {

					const int iA = iElementA + i;
					const int iB = iElementB + j;

					// Inverse Jacobian
					const double dInvJacobian =
						1.0 / m_dLocalJacobian[k][i][j];

					// Aliases for alpha and beta velocities
					const double dConUa = m_dAuxDataNode[ConUaIx][k][i][j];
					const double dConUb = m_dAuxDataNode[ConUbIx][k][i][j];
					const double dConUx = m_dAuxDataNode[ConUxIx][k][i][j];

					// Derivative of the kinetic energy
					double dDaKE = 0.0;
					double dDbKE = 0.0;

					// Derivatives of the pressure field
					double dDaP = 0.0;
					double dDbP = 0.0;

					// Gradient of the divergence
					double dDaDiv = 0.0;
					double dDbDiv = 0.0;

					// Calculate derivatives in the alpha direction
					double dDaRhoFluxA = 0.0;
					double dDaPressureFluxA = 0.0;

					// Calculate derivatives in the beta direction
					double dDbRhoFluxB = 0.0;
					double dDbPressureFluxB = 0.0;

					for (int s = 0; s < m_nHorizontalOrder; s++) {

						// Update density: Variational formulation
						dDaRhoFluxA -=
							m_dAlphaMassFlux[s][j]
							* dStiffness1D[i][s];

						// Update pressure: Variational formulation
						dDaPressureFluxA -=
							m_dAlphaPressureFlux[s][j]
							* dStiffness1D[i][s];

						// Derivative of (Exner) pressure with respect to alpha
						dDaP +=
							m_dAuxDataNode[ExnerIx][k][s][j]
							* dDxBasis1D[s][i];

						// Derivative of specific kinetic energy wrt alpha
						dDaKE +=
							m_dAuxDataNode[KIx][k][s][j]
							* dDxBasis1D[s][i];

						// Update density: Variational formulation
						dDbRhoFluxB -=
							m_dBetaMassFlux[i][s]
							* dStiffness1D[j][s];

						// Update pressure: Variational formulation
						dDbPressureFluxB -=
							m_dBetaPressureFlux[i][s]
							* dStiffness1D[j][s];
 
						// Derivative of (Exner) pressure with respect to beta
						dDbP +=
							m_dAuxDataNode[ExnerIx][k][i][s]
							* dDxBasis1D[s][j];

						// Derivative of specific kinetic energy wrt beta
						dDbKE +=
							m_dAuxDataNode[KIx][k][i][s]
							* dDxBasis1D[s][j];
					}

					// Scale derivatives
					dDaRhoFluxA *= dInvElementDeltaA;
					dDbRhoFluxB *= dInvElementDeltaB;

					dDaPressureFluxA *= dInvElementDeltaA;
					dDbPressureFluxB *= dInvElementDeltaB;

					dDaP *= dInvElementDeltaA;
					dDbP *= dInvElementDeltaB;

					dDaKE *= dInvElementDeltaA;
					dDbKE *= dInvElementDeltaB;

					// Pointwise momentum updates
					double dLocalUpdateUa = 0.0;
					double dLocalUpdateUb = 0.0;

					// Updates due to rotational terms
					dLocalUpdateUa += m_dAuxDataNode[UCrossZetaAIx][k][i][j];
					dLocalUpdateUb += m_dAuxDataNode[UCrossZetaBIx][k][i][j];

					// Coriolis terms
					dLocalUpdateUa +=
						m_dLocalCoriolisF[i][j]
						* m_dLocalJacobian2D[i][j]
						* dConUb;

					dLocalUpdateUb -=
						m_dLocalCoriolisF[i][j]
						* m_dLocalJacobian2D[i][j]
						* dConUa;

					// Pressure gradient force
#if defined(FORMULATION_RHOTHETA_P)
					const double dPressureGradientForceUa =
						dDaP / dataInitialNode[RIx][k][iA][iB];
					const double dPressureGradientForceUb =
						dDbP / dataInitialNode[RIx][k][iA][iB];
#endif
#if defined(FORMULATION_RHOTHETA_PI)
					const double dPressureGradientForceUa =
						dDaP * dataInitialNode[PIx][k][iA][iB]
						/ dataInitialNode[RIx][k][iA][iB];
					const double dPressureGradientForceUb =
						dDbP * dataInitialNode[PIx][k][iA][iB]
						/ dataInitialNode[RIx][k][iA][iB];
#endif
#if defined(FORMULATION_THETA)
					const double dPressureGradientForceUa =
						dDaP * dataInitialNode[PIx][k][iA][iB];
					const double dPressureGradientForceUb =
						dDbP * dataInitialNode[PIx][k][iA][iB];
#endif

					// Gravity
					const double dDaPhi = phys.GetG() * m_dLocalDerivR[k][i][j][0];
					const double dDbPhi = phys.GetG() * m_dLocalDerivR[k][i][j][1];

					// Horizontal updates due to gradient terms
					const double dDaUpdate =
						dPressureGradientForceUa + dDaKE + dDaPhi;

					const double dDbUpdate =
						dPressureGradientForceUb + dDbKE + dDbPhi;

					// Apply gradient term update to total update
					dLocalUpdateUa -= dDaUpdate;
					dLocalUpdateUb -= dDbUpdate;

 					// Apply update to horizontal velocity on model levels
					dataUpdateNode[UIx][k][iA][iB] +=
						dDeltaT * dLocalUpdateUa;

					dataUpdateNode[VIx][k][iA][iB] +=
					dDeltaT * dLocalUpdateUb;

					// Update density on model levels
					dataUpdateNode[RIx][k][iA][iB] -=
						dDeltaT * dInvJacobian * (
							  dDaRhoFluxA
							+ dDbRhoFluxB);

#if defined(FORMULATION_RHOTHETA_PI) || defined(FORMULATION_RHOTHETA_P)
					// Update RhoTheta on model levels
					dataUpdateNode[PIx][k][iA][iB] -=
						dDeltaT * dInvJacobian * (
							  dDaPressureFluxA
							+ dDbPressureFluxB);
#endif

#ifdef FORMULATION_THETA
					// Derivatives of the theta field
					double dDaTheta = 0.0;
					double dDbTheta = 0.0;

					for (int s = 0; s < m_nHorizontalOrder; s++) {
						dDaTheta +=
							dataInitialNode[PIx][k][iElementA+s][iB]
							* dDxBasis1D[s][i];

						dDbTheta +=
							dataInitialNode[PIx][k][iA][iElementB+s]
							* dDxBasis1D[s][j];
					}

					dDaTheta *= dInvElementDeltaA;
					dDbTheta *= dInvElementDeltaB;

					// Update Theta on model levels
					dataUpdateNode[PIx][k][iA][iB] -=
						dDeltaT * (dConUa * dDaTheta + dConUb * dDbTheta);
#endif

					// Update tracers
					for (int c = 0; c < nTracerCount; c++) {
						double dDaTracerFluxA = 0.0;
						double dDbTracerFluxB = 0.0;

						for (int s = 0; s < m_nHorizontalOrder; s++) {
							dDaTracerFluxA -=
								m_dAlphaTracerFlux[c][s][j]
								* dStiffness1D[i][s];

							dDbTracerFluxB -=
								m_dBetaTracerFlux[c][i][s]
								* dStiffness1D[j][s];
						}

						dDaTracerFluxA *= dInvElementDeltaA;
						dDbTracerFluxB *= dInvElementDeltaB;

						dataUpdateTracer[c][k][iA][iB] -=
							dDeltaT * dInvJacobian * (
								  dDaTracerFluxA
								+ dDbTracerFluxB);
					}
				}
				}
			}

			// Update vertical velocity on interfaces
			{
				// Update vertical velocity at bottom boundary
				for (int i = 0; i < m_nHorizontalOrder; i++) {
				for (int j = 0; j < m_nHorizontalOrder; j++) {

					const int iA = iElementA + i;
					const int iB = iElementB + j;

					// Interpolate horizontal velocity to bottom boundary
					const double dU0 =
						pGrid->InterpolateNodeToREdge(
							&(dataUpdateNode[UIx][0][iA][iB]),
							NULL, 0, 0.0,
							nVerticalStateStride);

					const double dV0 =
						pGrid->InterpolateNodeToREdge(
							&(dataUpdateNode[VIx][0][iA][iB]),
							NULL, 0, 0.0,
							nVerticalStateStride);

					// Update vertical velocity on boundary
					dataUpdateREdge[WIx][0][iA][iB] =
						- ( dContraMetricXiREdge[0][iA][iB][0] * dU0
						  + dContraMetricXiREdge[0][iA][iB][1] * dV0)
							/ dContraMetricXiREdge[0][iA][iB][2];
				}
				}

				// Update interior interfaces
				for (int k = 1; k < nRElements; k++) {
				for (int i = 0; i < m_nHorizontalOrder; i++) {
				for (int j = 0; j < m_nHorizontalOrder; j++) {

					const int iA = iElementA + i;
					const int iB = iElementB + j;

					double dUCrossZetaXIx =
						pGrid->InterpolateNodeToREdge(
							&(m_dAuxDataNode[UCrossZetaXIx][0][i][j]),
							NULL,
							k,
							0.0,
							nVerticalElementStride);

					dataUpdateREdge[WIx][k][iA][iB] +=
						dDeltaT * dUCrossZetaXIx;
				}
				}
				}
			}
		}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void HorizontalDynamicsFEMV2::StepExplicit(
	int iDataInitial,
	int iDataUpdate,
	const Time & time,
	double dDeltaT
) {
	if (iDataInitial == iDataUpdate) {
		_EXCEPTIONT(
			"HorizontalDynamics Step must have iDataInitial != iDataUpdate");
	}

	// Get a copy of the GLL grid
	GridGLL * pGrid = dynamic_cast<GridGLL*>(m_model.GetGrid());

	// Equation set
	const EquationSet & eqn = m_model.GetEquationSet();

	// Step the primitive nonhydrostatic equations
	if (eqn.GetType() == EquationSet::PrimitiveNonhydrostaticEquations) {
		StepNonhydrostaticPrimitive(iDataInitial, iDataUpdate, time, dDeltaT);

	// Step the shallow water equations
	} else if (eqn.GetType() == EquationSet::ShallowWaterEquations) {
		_EXCEPTIONT("Not implemented");

	// Invalid EquationSet
	} else {
		_EXCEPTIONT("Invalid EquationSet");
	}

	// Apply positive definite filter to tracers
	FilterNegativeTracers(iDataUpdate);
}

///////////////////////////////////////////////////////////////////////////////

void HorizontalDynamicsFEMV2::ApplyScalarHyperdiffusion(
	int iDataInitial,
	int iDataUpdate,
	double dDeltaT,
	double dNu,
	bool fScaleNuLocally,
	int iComponent,
	bool fRemoveRefState
) {
	// Indices of EquationSet variables
	const int UIx = 0;
	const int VIx = 1;
	const int PIx = 2;
	const int WIx = 3;
	const int RIx = 4;

	// Get a copy of the GLL grid
	GridGLL * pGrid = dynamic_cast<GridGLL*>(m_model.GetGrid());

	// Number of radial elements in grid
	int nRElements = pGrid->GetRElements();

	// Check argument
	if (iComponent < (-1)) {
		_EXCEPTIONT("Invalid component index");
	}

	// Perform local update
	for (int n = 0; n < pGrid->GetActivePatchCount(); n++) {
		GridPatchGLL * pPatch =
			dynamic_cast<GridPatchGLL*>(pGrid->GetActivePatch(n));

		const PatchBox & box = pPatch->GetPatchBox();

		const DataArray3D<double> & dElementAreaNode =
			pPatch->GetElementAreaNode();
		const DataArray3D<double> & dJacobianNode =
			pPatch->GetJacobian();
		const DataArray3D<double> & dJacobianREdge =
			pPatch->GetJacobianREdge();
		const DataArray3D<double> & dContraMetricA =
			pPatch->GetContraMetric2DA();
		const DataArray3D<double> & dContraMetricB =
			pPatch->GetContraMetric2DB();

		// Grid data
		DataArray4D<double> & dataInitialNode =
			pPatch->GetDataState(iDataInitial, DataLocation_Node);

		DataArray4D<double> & dataInitialREdge =
			pPatch->GetDataState(iDataInitial, DataLocation_REdge);

		DataArray4D<double> & dataUpdateNode =
			pPatch->GetDataState(iDataUpdate, DataLocation_Node);

		DataArray4D<double> & dataUpdateREdge =
			pPatch->GetDataState(iDataUpdate, DataLocation_REdge);

		const DataArray4D<double> & dataRefNode =
			pPatch->GetReferenceState(DataLocation_Node);

		const DataArray4D<double> & dataRefREdge =
			pPatch->GetReferenceState(DataLocation_REdge);

		// Tracer data
		DataArray4D<double> & dataInitialTracer =
			pPatch->GetDataTracers(iDataInitial);

		DataArray4D<double> & dataUpdateTracer =
			pPatch->GetDataTracers(iDataUpdate);

		// Element grid spacing and derivative coefficients
		const double dElementDeltaA = pPatch->GetElementDeltaA();
		const double dElementDeltaB = pPatch->GetElementDeltaB();

		const double dInvElementDeltaA = 1.0 / dElementDeltaA;
		const double dInvElementDeltaB = 1.0 / dElementDeltaB;

		const DataArray2D<double> & dDxBasis1D = pGrid->GetDxBasis1D();
		const DataArray2D<double> & dStiffness1D = pGrid->GetStiffness1D();

		// Number of finite elements
		int nElementCountA = pPatch->GetElementCountA();
		int nElementCountB = pPatch->GetElementCountB();

		// Compute new hyperviscosity coefficient
		double dLocalNu = dNu;

		if (fScaleNuLocally) {
			double dReferenceLength = pGrid->GetReferenceLength();
			if (dReferenceLength != 0.0) {
				dLocalNu *= pow(dElementDeltaA / dReferenceLength, 3.2);
			}
		}

		// State variables and tracers
		for (int iType = 0; iType < 2; iType++) {

			int nComponentStart;
			int nComponentEnd;

			if (iType == 0) {
				nComponentStart = 2;
				nComponentEnd = m_model.GetEquationSet().GetComponents();

				if (iComponent != (-1)) {
					if (iComponent >= nComponentEnd) {
						_EXCEPTIONT("Invalid component index");
					}

					nComponentStart = iComponent;
					nComponentEnd = iComponent+1;
				}

			} else {
				nComponentStart = 0;
				nComponentEnd = m_model.GetEquationSet().GetTracers();

				if (iComponent != (-1)) {
					continue;
				}
			}

			// Loop over all components
			for (int c = nComponentStart; c < nComponentEnd; c++) {

				int nElementCountR;

				const DataArray4D<double> * pDataInitial;
				DataArray4D<double> * pDataUpdate;
				const DataArray4D<double> * pDataRef;
				const DataArray3D<double> * pJacobian;

				if (iType == 0) {
					if (pGrid->GetVarLocation(c) == DataLocation_Node) {
						pDataInitial = &dataInitialNode;
						pDataUpdate  = &dataUpdateNode;
						pDataRef = &dataRefNode;
						nElementCountR = nRElements;
						pJacobian = &dJacobianNode;

					} else if (pGrid->GetVarLocation(c) == DataLocation_REdge) {
						pDataInitial = &dataInitialREdge;
						pDataUpdate  = &dataUpdateREdge;
						pDataRef = &dataRefREdge;
						nElementCountR = nRElements + 1;
						pJacobian = &dJacobianREdge;

					} else {
						_EXCEPTIONT("UNIMPLEMENTED");
					}

				} else {
					pDataInitial = &dataInitialTracer;
					pDataUpdate = &dataUpdateTracer;
					nElementCountR = nRElements;
					pJacobian = &dJacobianNode;
				}

				// Loop over all finite elements
				for (int k = 0; k < nElementCountR; k++) {
				for (int a = 0; a < nElementCountA; a++) {
				for (int b = 0; b < nElementCountB; b++) {

					int iElementA =
						a * m_nHorizontalOrder + box.GetHaloElements();
					int iElementB =
						b * m_nHorizontalOrder + box.GetHaloElements();

					// Store the buffer state
					for (int i = 0; i < m_nHorizontalOrder; i++) {
					for (int j = 0; j < m_nHorizontalOrder; j++) {
						int iA = iElementA + i;
						int iB = iElementB + j;

						m_dBufferState[i][j] = (*pDataInitial)[c][k][iA][iB];
					}
					}

					// Remove the reference state from the buffer state
					if (fRemoveRefState) {
						for (int i = 0; i < m_nHorizontalOrder; i++) {
						for (int j = 0; j < m_nHorizontalOrder; j++) {
							int iA = iElementA + i;
							int iB = iElementB + j;

							m_dBufferState[i][j] -=
								(*pDataRef)[c][k][iA][iB];
						}
						}
					}

					// Calculate the pointwise gradient of the scalar field
					for (int i = 0; i < m_nHorizontalOrder; i++) {
					for (int j = 0; j < m_nHorizontalOrder; j++) {
						int iA = iElementA + i;
						int iB = iElementB + j;

						double dDaPsi = 0.0;
						double dDbPsi = 0.0;
						for (int s = 0; s < m_nHorizontalOrder; s++) {
							dDaPsi +=
								m_dBufferState[s][j]
								* dDxBasis1D[s][i];

							dDbPsi +=
								m_dBufferState[i][s]
								* dDxBasis1D[s][j];
						}

						dDaPsi *= dInvElementDeltaA;
						dDbPsi *= dInvElementDeltaB;

						m_dJGradientA[i][j] = (*pJacobian)[k][iA][iB] * (
							+ dContraMetricA[iA][iB][0] * dDaPsi
							+ dContraMetricA[iA][iB][1] * dDbPsi);

						m_dJGradientB[i][j] = (*pJacobian)[k][iA][iB] * (
							+ dContraMetricB[iA][iB][0] * dDaPsi
							+ dContraMetricB[iA][iB][1] * dDbPsi);
					}
					}

					// Pointwise updates
					for (int i = 0; i < m_nHorizontalOrder; i++) {
					for (int j = 0; j < m_nHorizontalOrder; j++) {
						int iA = iElementA + i;
						int iB = iElementB + j;

						// Inverse Jacobian and Jacobian
						const double dInvJacobian = 1.0 / (*pJacobian)[k][iA][iB];

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

						dUpdateA *= dInvElementDeltaA;
						dUpdateB *= dInvElementDeltaB;

						// Apply update
						(*pDataUpdate)[c][k][iA][iB] -=
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
}

///////////////////////////////////////////////////////////////////////////////

void HorizontalDynamicsFEMV2::ApplyVectorHyperdiffusion(
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

	// Apply viscosity to reference state
	bool fApplyToRefState = false;
	if (iDataInitial == DATA_INDEX_REFERENCE) {
		iDataInitial = 0;
		fApplyToRefState = true;
	}

	// Loop over all patches
	for (int n = 0; n < pGrid->GetActivePatchCount(); n++) {
		GridPatchGLL * pPatch =
			dynamic_cast<GridPatchGLL*>(pGrid->GetActivePatch(n));

		const PatchBox & box = pPatch->GetPatchBox();

		const DataArray2D<double> & dJacobian2D =
			pPatch->GetJacobian2D();
		const DataArray3D<double> & dContraMetric2DA =
			pPatch->GetContraMetric2DA();
		const DataArray3D<double> & dContraMetric2DB =
			pPatch->GetContraMetric2DB();

		DataArray4D<double> & dataInitial =
			pPatch->GetDataState(iDataInitial, DataLocation_Node);

		DataArray4D<double> & dataUpdate =
			pPatch->GetDataState(iDataUpdate, DataLocation_Node);

		DataArray4D<double> & dataRef =
			pPatch->GetReferenceState(DataLocation_Node);

		// Element grid spacing and derivative coefficients
		const double dElementDeltaA = pPatch->GetElementDeltaA();
		const double dElementDeltaB = pPatch->GetElementDeltaB();

		const double dInvElementDeltaA = 1.0 / dElementDeltaA;
		const double dInvElementDeltaB = 1.0 / dElementDeltaB;

		const DataArray2D<double> & dDxBasis1D = pGrid->GetDxBasis1D();
		const DataArray2D<double> & dStiffness1D = pGrid->GetStiffness1D();

		// Compute curl and divergence of U on the grid
		DataArray3D<double> dataUa;
		dataUa.SetSize(
			dataInitial.GetSize(1),
			dataInitial.GetSize(2),
			dataInitial.GetSize(3));

		DataArray3D<double> dataUb;
		dataUb.SetSize(
			dataInitial.GetSize(1),
			dataInitial.GetSize(2),
			dataInitial.GetSize(3));

		if (fApplyToRefState) {
			dataUa.AttachToData(&(dataRef[UIx][0][0][0]));
			dataUb.AttachToData(&(dataRef[VIx][0][0][0]));
		} else {
			dataUa.AttachToData(&(dataInitial[UIx][0][0][0]));
			dataUb.AttachToData(&(dataInitial[VIx][0][0][0]));
		}

		// Compute curl and divergence of U on the grid
		pPatch->ComputeCurlAndDiv(dataUa, dataUb);

		// Get curl and divergence
		const DataArray3D<double> & dataCurl = pPatch->GetDataVorticity();
		const DataArray3D<double> & dataDiv  = pPatch->GetDataDivergence();

		// Compute new hyperviscosity coefficient
		double dLocalNuDiv  = dNuDiv;
		double dLocalNuVort = dNuVort;

		if (fScaleNuLocally) {
			double dReferenceLength = pGrid->GetReferenceLength();
			if (dReferenceLength != 0.0) {
				dLocalNuDiv =
					dLocalNuDiv  * pow(dElementDeltaA / dReferenceLength, 3.2);
				dLocalNuVort =
					dLocalNuVort * pow(dElementDeltaA / dReferenceLength, 3.2);
			}
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

				dDaDiv *= dInvElementDeltaA;
				dDbDiv *= dInvElementDeltaB;

				dDaCurl *= dInvElementDeltaA;
				dDbCurl *= dInvElementDeltaB;

				// Apply update
				double dUpdateUa =
					+ dLocalNuDiv * dDaDiv
					- dLocalNuVort * dJacobian2D[iA][iB] * (
						  dContraMetric2DB[iA][iB][0] * dDaCurl
						+ dContraMetric2DB[iA][iB][1] * dDbCurl);

				double dUpdateUb =
					+ dLocalNuDiv * dDbDiv
					+ dLocalNuVort * dJacobian2D[iA][iB] * (
						  dContraMetric2DA[iA][iB][0] * dDaCurl
						+ dContraMetric2DA[iA][iB][1] * dDbCurl);

				dataUpdate[UIx][k][iA][iB] -= dDeltaT * dUpdateUa;

				if (pGrid->GetIsCartesianXZ() == false) {
					dataUpdate[VIx][k][iA][iB] -= dDeltaT * dUpdateUb;
				}
			}
			}
		}
		}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

//#pragma message "jeguerra: Clean up this function"

void HorizontalDynamicsFEMV2::ApplyRayleighFriction(
	int iDataUpdate,
	double dDeltaT
) {
	// Get a copy of the GLL grid
	GridGLL * pGrid = dynamic_cast<GridGLL*>(m_model.GetGrid());

	// Number of components to hit with friction
	int nComponents = m_model.GetEquationSet().GetComponents();

	// Equation set being solved
	int nEqSet = m_model.GetEquationSet().GetType();

	bool fCartXZ = pGrid->GetIsCartesianXZ();

	int nEffectiveC[nComponents];
	// 3D primitive nonhydro models with no density treatment
	if ((nEqSet == EquationSet::PrimitiveNonhydrostaticEquations) && !fCartXZ) {
		nEffectiveC[0] = 0; nEffectiveC[1] = 1;
		nEffectiveC[2] = 2; nEffectiveC[3] = 3;
		nEffectiveC[nComponents - 1] = 0;
		nComponents = nComponents - 1;
	}
	// 2D Cartesian XZ primitive nonhydro models with no density treatment
	else if ((nEqSet == EquationSet::PrimitiveNonhydrostaticEquations) && fCartXZ) {
		nEffectiveC[0] = 0;
		nEffectiveC[1] = 2;
		nEffectiveC[2] = 3;
		nEffectiveC[nComponents - 2] = 0;
		nEffectiveC[nComponents - 1] = 0;
		nComponents = nComponents - 2;
	}
	// Other model types (advection, shallow water, mass coord)
	else {
		for (int nc = 0; nc < nComponents; nc++) {
			nEffectiveC[nc] = nc;
		}
	}

	// Subcycle the rayleigh update
	int nRayleighCycles = 10;
	double dRayleighFactor = 1.0 / nRayleighCycles;

	// Perform local update
	for (int n = 0; n < pGrid->GetActivePatchCount(); n++) {
		GridPatchGLL * pPatch =
			dynamic_cast<GridPatchGLL*>(pGrid->GetActivePatch(n));

		const PatchBox & box = pPatch->GetPatchBox();

		// Grid data
		DataArray4D<double> & dataUpdateNode =
			pPatch->GetDataState(iDataUpdate, DataLocation_Node);

		DataArray4D<double> & dataUpdateREdge =
			pPatch->GetDataState(iDataUpdate, DataLocation_REdge);

		// Reference state
		const DataArray4D<double> & dataReferenceNode =
			pPatch->GetReferenceState(DataLocation_Node);

		const DataArray4D<double> & dataReferenceREdge =
			pPatch->GetReferenceState(DataLocation_REdge);

		// Rayleigh friction strength
		const DataArray3D<double> & dataRayleighStrengthNode =
			pPatch->GetRayleighStrength(DataLocation_Node);

		const DataArray3D<double> & dataRayleighStrengthREdge =
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

				// Loop over all effective components
				for (int c = 0; c < nComponents; c++) {
					if (pGrid->GetVarLocation(nEffectiveC[c]) ==
						DataLocation_Node) {
						for (int si = 0; si < nRayleighCycles; si++) {
							dNuNode = 1.0 / (1.0 + dRayleighFactor * dDeltaT * dNu);
							dataUpdateNode[nEffectiveC[c]][k][i][j] =
								dNuNode * dataUpdateNode[nEffectiveC[c]][k][i][j]
								+ (1.0 - dNuNode)
								* dataReferenceNode[nEffectiveC[c]][k][i][j];
						}
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

				// Loop over all effective components
				for (int c = 0; c < nComponents; c++) {
					if (pGrid->GetVarLocation(nEffectiveC[c]) ==
						DataLocation_REdge) {
						for (int si = 0; si < nRayleighCycles; si++) {
							dNuREdge = 1.0 / (1.0 + dRayleighFactor * dDeltaT * dNu);
							dataUpdateREdge[nEffectiveC[c]][k][i][j] =
							dNuREdge * dataUpdateREdge[nEffectiveC[c]][k][i][j]
							+ (1.0 - dNuREdge)
							* dataReferenceREdge[nEffectiveC[c]][k][i][j];
						}
					}
				}
			}
		}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

int HorizontalDynamicsFEMV2::GetSubStepAfterSubCycleCount() {
	int iSubStepCount = m_nHyperviscosityOrder / 2;

	return iSubStepCount;
}

///////////////////////////////////////////////////////////////////////////////

int HorizontalDynamicsFEMV2::SubStepAfterSubCycle(
	int iDataInitial,
	int iDataUpdate,
	int iDataWorking,
	const Time & time,
	double dDeltaT,
	int iSubStep
) {
	// Get the GLL grid
	GridGLL * pGrid = dynamic_cast<GridGLL*>(m_model.GetGrid());

	// First calculation of Laplacian
	if (iSubStep == 0) {

		// Apply scalar and vector hyperviscosity (first application)
		pGrid->ZeroData(iDataWorking, DataType_State);
		pGrid->ZeroData(iDataWorking, DataType_Tracers);

		ApplyScalarHyperdiffusion(
			iDataInitial, iDataWorking, 1.0, 1.0, false);
		ApplyVectorHyperdiffusion(
			iDataInitial, iDataWorking, 1.0, 1.0, 1.0, false);

		return iDataWorking;

	// Second calculation of Laplacian
	} else if (iSubStep == 1) {

		// Copy initial data to updated data
		pGrid->CopyData(iDataInitial, iDataUpdate, DataType_State);
		pGrid->CopyData(iDataInitial, iDataUpdate, DataType_Tracers);

		// Apply scalar and vector hyperviscosity (second application)
		ApplyScalarHyperdiffusion(
			iDataWorking, iDataUpdate, -dDeltaT, m_dNuScalar, true);
		ApplyVectorHyperdiffusion(
			iDataWorking, iDataUpdate, -dDeltaT, m_dNuDiv, m_dNuVort, true);

		// Apply positive definite filter to tracers
		FilterNegativeTracers(iDataUpdate);

#ifdef APPLY_RAYLEIGH_WITH_HYPERVIS
		// Apply Rayleigh damping
		if (pGrid->HasRayleighFriction()) {
			ApplyRayleighFriction(iDataUpdate, dDeltaT);
		}
#endif
		return iDataUpdate;
	}

	_EXCEPTIONT("Invalid iSubStep");
}

///////////////////////////////////////////////////////////////////////////////

void HorizontalDynamicsFEMV2::StepAfterSubCycle(
	int iDataInitial,
	int iDataUpdate,
	int iDataWorking,
	const Time & time,
	double dDeltaT
) {
	// Start the function timer
	FunctionTimer timer("StepAfterSubCycle");

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

	// Copy initial data to updated data
	pGrid->CopyData(iDataInitial, iDataUpdate, DataType_State);
	pGrid->CopyData(iDataInitial, iDataUpdate, DataType_Tracers);

	// No hyperdiffusion
	if ((m_dNuScalar == 0.0) && (m_dNuDiv == 0.0) && (m_dNuVort == 0.0)) {

	// Apply hyperdiffusion
	} else if (m_nHyperviscosityOrder == 0) {

	// Apply viscosity
	} else if (m_nHyperviscosityOrder == 2) {

		// Apply scalar and vector viscosity
		ApplyScalarHyperdiffusion(
			iDataInitial, iDataUpdate, dDeltaT, m_dNuScalar, false);
		ApplyVectorHyperdiffusion(
			iDataInitial, iDataUpdate, -dDeltaT, m_dNuDiv, m_dNuVort, false);

		// Apply positive definite filter to tracers
		FilterNegativeTracers(iDataUpdate);

		// Apply Direct Stiffness Summation
		pGrid->ApplyDSS(iDataUpdate, DataType_State);
		pGrid->ApplyDSS(iDataUpdate, DataType_Tracers);

	// Apply hyperviscosity
	} else if (m_nHyperviscosityOrder == 4) {

		// Apply scalar and vector hyperviscosity (first application)
		pGrid->ZeroData(iDataWorking, DataType_State);
		pGrid->ZeroData(iDataWorking, DataType_Tracers);

		ApplyScalarHyperdiffusion(
			iDataInitial, iDataWorking, 1.0, 1.0, false);
		ApplyVectorHyperdiffusion(
			iDataInitial, iDataWorking, 1.0, 1.0, 1.0, false);

		// Apply Direct Stiffness Summation
		pGrid->ApplyDSS(iDataWorking, DataType_State);
		pGrid->ApplyDSS(iDataWorking, DataType_Tracers);

		// Apply scalar and vector hyperviscosity (second application)
		ApplyScalarHyperdiffusion(
			iDataWorking, iDataUpdate, -dDeltaT, m_dNuScalar, true);
		ApplyVectorHyperdiffusion(
			iDataWorking, iDataUpdate, -dDeltaT, m_dNuDiv, m_dNuVort, true);

		// Apply positive definite filter to tracers
		FilterNegativeTracers(iDataUpdate);

		// Apply Direct Stiffness Summation
		pGrid->ApplyDSS(iDataUpdate, DataType_State);
		pGrid->ApplyDSS(iDataUpdate, DataType_Tracers);

	// Invalid viscosity order
	} else {
		_EXCEPTIONT("Invalid viscosity order");
	}

#ifdef APPLY_RAYLEIGH_WITH_HYPERVIS
		// Apply Rayleigh damping
		if (pGrid->HasRayleighFriction()) {
			ApplyRayleighFriction(iDataUpdate, dDeltaT);
		}
#endif
}

///////////////////////////////////////////////////////////////////////////////
