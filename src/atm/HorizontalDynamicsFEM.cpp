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

#define RESIDUAL_DIFFUSION_THERMO
//#define RESIDUAL_DIFFUSION_RHO

///////////////////////////////////////////////////////////////////////////////

HorizontalDynamicsFEM::HorizontalDynamicsFEM(
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

void HorizontalDynamicsFEM::Initialize() {

#if defined(PROGNOSTIC_CONTRAVARIANT_MOMENTA)
	_EXCEPTIONT("Prognostic contrvariant velocities not supported");
#endif

	// Get a copy of the GLL grid
	GridGLL * pGrid = dynamic_cast<GridGLL*>(m_model.GetGrid());
	if (pGrid == NULL) {
		_EXCEPTIONT("Grid must be of type GridGLL");
	}

#if defined(FIXED_HORIZONTAL_ORDER)
	const int nHorizontalOrder = FIXED_HORIZONTAL_ORDER;
	if (nHorizontalOrder != m_nHorizontalOrder) {
		_EXCEPTIONT("Command line order must match FIXED_HORIZONTAL_ORDER");
	}
#else
	const int nHorizontalOrder = m_nHorizontalOrder;
#endif

#if defined(FIXED_RELEMENTS)
	const int nRElements = FIXED_RELEMENTS;
	if (nRElements != pGrid->GetRElements()) {
		_EXCEPTIONT("Command line levels must match FIXED_RELEMENTS");
	}
#else
	const int nRElements = pGrid->GetRElements();
#endif

	// Number of tracers
	const int nTracerCount = m_model.GetEquationSet().GetTracers();

	// Initialize the alpha and beta mass fluxes
	m_dAlphaMassFlux.Allocate(
		nHorizontalOrder,
		nHorizontalOrder,
		nRElements);

	m_dBetaMassFlux.Allocate(
		nHorizontalOrder,
		nHorizontalOrder,
		nRElements);

	// Initialize the alpha and beta mass fluxes for mass fix
	m_dAlphaElMassFlux.Allocate(
		nHorizontalOrder,
		nHorizontalOrder,
		nRElements);

	m_dBetaElMassFlux.Allocate(
		nHorizontalOrder,
		nHorizontalOrder,
		nRElements);

	// Initialize element-wise arrays used for mass fix
	m_dElementMassFluxA.Allocate(
		nRElements);

	m_dElementMassFluxB.Allocate(
		nRElements);

	m_dElTotalArea.Allocate(
		nRElements);

	// Initialize the alpha and beta pressure fluxes
	m_dAlphaPressureFlux.Allocate(
		nHorizontalOrder,
		nHorizontalOrder,
		nRElements);

	m_dBetaPressureFlux.Allocate(
		nHorizontalOrder,
		nHorizontalOrder,
		nRElements);

	// Initialize tracer fluxes
	if (nTracerCount > 0) {
		m_dAlphaTracerFlux.Allocate(
			nTracerCount,
			nHorizontalOrder,
			nHorizontalOrder,
			nRElements);

		m_dBetaTracerFlux.Allocate(
			nTracerCount,
			nHorizontalOrder,
			nHorizontalOrder,
			nRElements);
	}

	// Auxiliary data
	m_dAuxDataNode.Allocate(
		9,
		nHorizontalOrder,
		nHorizontalOrder,
		nRElements);

	m_dAuxDataREdge.Allocate(
		9,
		nHorizontalOrder,
		nHorizontalOrder,
		nRElements+1);

	m_dDivergence.Allocate(
		nHorizontalOrder,
		nHorizontalOrder,
		nRElements);

	// Contravariant metric terms
	m_dLocalCoriolisF.Allocate(
		nHorizontalOrder,
		nHorizontalOrder);

	m_dLocalJacobian2D.Allocate(
		nHorizontalOrder,
		nHorizontalOrder);

	m_dLocalJacobian.Allocate(
		nHorizontalOrder,
		nHorizontalOrder,
		nRElements);

	m_dLocalDerivR.Allocate(
		nHorizontalOrder,
		nHorizontalOrder,
		nRElements,
		3);

	m_dLocalContraMetric.Allocate(
		nHorizontalOrder,
		nHorizontalOrder,
		nRElements,
		6);

	// Initialize buffers for derivatives of Jacobian
	m_dJGradientA.Allocate(
		nHorizontalOrder,
		nHorizontalOrder);

	m_dJGradientB.Allocate(
		nHorizontalOrder,
		nHorizontalOrder);

	// Buffer state
	m_dBufferState.Allocate(
		nHorizontalOrder,
		nHorizontalOrder);
}

///////////////////////////////////////////////////////////////////////////////

void HorizontalDynamicsFEM::FilterNegativeTracers(
	int iDataUpdate
) {
#ifdef POSITIVE_DEFINITE_FILTER_TRACERS
	// Get a copy of the GLL grid
	GridGLL * pGrid = dynamic_cast<GridGLL*>(m_model.GetGrid());

#if defined(FIXED_HORIZONTAL_ORDER)
	const int nHorizontalOrder = FIXED_HORIZONTAL_ORDER;
	if (nHorizontalOrder != m_nHorizontalOrder) {
		_EXCEPTIONT("Command line order must match FIXED_HORIZONTAL_ORDER");
	}
#else
	const int nHorizontalOrder = m_nHorizontalOrder;
#endif

#if defined(FIXED_RELEMENTS)
	const int nRElements = FIXED_RELEMENTS;
	if (nRElements != pGrid->GetRElements()) {
		_EXCEPTIONT("Command line levels must match FIXED_RELEMENTS");
	}
#else
	const int nRElements = pGrid->GetRElements();
#endif

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

				for (int i = 0; i < nHorizontalOrder; i++) {
				for (int j = 0; j < nHorizontalOrder; j++) {

					int iA = a * nHorizontalOrder + i + box.GetHaloElements();
					int iB = b * nHorizontalOrder + j + box.GetHaloElements();

					double dPointwiseMass =
						  dataUpdateTracer(c,iA,iB,k)
						* dElementAreaNode(iA,iB,k);

					dTotalMass += dPointwiseMass;

					if (dataUpdateTracer(c,iA,iB,k) >= 0.0) {
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

				for (int i = 0; i < nHorizontalOrder; i++) {
				for (int j = 0; j < nHorizontalOrder; j++) {

					int iA = a * nHorizontalOrder + i + box.GetHaloElements();
					int iB = b * nHorizontalOrder + j + box.GetHaloElements();

					if (dataUpdateTracer(c,iA,iB,k) > 0.0) {
						dataUpdateTracer(c,iA,iB,k) *= dR;
					} else {
						dataUpdateTracer(c,iA,iB,k) = 0.0;
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

void HorizontalDynamicsFEM::StepShallowWater(
	int iDataInitial,
	int iDataUpdate,
	const Time & time,
	double dDeltaT
) {

	// Get a copy of the GLL grid
	GridGLL * pGrid = dynamic_cast<GridGLL*>(m_model.GetGrid());

#if defined(FIXED_HORIZONTAL_ORDER)
	const int nHorizontalOrder = FIXED_HORIZONTAL_ORDER;
	if (nHorizontalOrder != m_nHorizontalOrder) {
		_EXCEPTIONT("Command line order must match FIXED_HORIZONTAL_ORDER");
	}
#else
	const int nHorizontalOrder = m_nHorizontalOrder;
#endif

#if defined(FIXED_RELEMENTS)
	const int nRElements = FIXED_RELEMENTS;
	if (nRElements != pGrid->GetRElements()) {
		_EXCEPTIONT("Command line levels must match FIXED_RELEMENTS");
	}
#else
	const int nRElements = pGrid->GetRElements();
#endif

	// Physical constants
	const PhysicalConstants & phys = m_model.GetPhysicalConstants();

	// Get indices of variables to update
	const int UIx = 0;
	const int VIx = 1;
	const int HIx = 2;

	// Indices of auxiliary data variables
	const int ConUaIx = 0;
	const int ConUbIx = 1;
	const int KIx = 4;

	// Perform local update
	for (int n = 0; n < pGrid->GetActivePatchCount(); n++) {
		GridPatchGLL * pPatch =
			dynamic_cast<GridPatchGLL*>(pGrid->GetActivePatch(n));

		const PatchBox & box = pPatch->GetPatchBox();

		const DataArray3D<double> & dElementAreaNode =
			pPatch->GetElementAreaNode();
		const DataArray2D<double> & dJacobian2D =
			pPatch->GetJacobian2D();
		const DataArray3D<double> & dContraMetric2DA =
			pPatch->GetContraMetric2DA();
		const DataArray3D<double> & dContraMetric2DB =
			pPatch->GetContraMetric2DB();
		const DataArray2D<double> & dCoriolisF =
			pPatch->GetCoriolisF();
		const DataArray2D<double> & dTopography =
			pPatch->GetTopography();

		// Data
		const DataArray4D<double> & dataInitialNode =
			pPatch->GetDataState(iDataInitial, DataLocation_Node);

		DataArray4D<double> & dataUpdateNode =
			pPatch->GetDataState(iDataUpdate, DataLocation_Node);

		// Element grid spacing and derivative coefficients
		const double dElementDeltaA = pPatch->GetElementDeltaA();
		const double dElementDeltaB = pPatch->GetElementDeltaB();

		const double dInvElementDeltaA = 1.0 / dElementDeltaA;
		const double dInvElementDeltaB = 1.0 / dElementDeltaB;

		const DataArray2D<double> & dDxBasis1D = pGrid->GetDxBasis1D();
		const DataArray2D<double> & dStiffness1D = pGrid->GetStiffness1D();

		// Get number of finite elements in each coordinate direction
		int nElementCountA = pPatch->GetElementCountA();
		int nElementCountB = pPatch->GetElementCountB();

		// Loop over all elements
		for (int a = 0; a < nElementCountA; a++) {
		for (int b = 0; b < nElementCountB; b++) {

			int iElementA = a * nHorizontalOrder + box.GetHaloElements();
			int iElementB = b * nHorizontalOrder + box.GetHaloElements();

			// Compute auxiliary data in element
			for (int k = 0; k < nRElements; k++) {
			for (int i = 0; i < nHorizontalOrder; i++) {
			for (int j = 0; j < nHorizontalOrder; j++) {

				const int iA = iElementA + i;
				const int iB = iElementB + j;

				// Contravariant velocities
				const double dCovUa = dataInitialNode(UIx,iA,iB,k);
				const double dCovUb = dataInitialNode(VIx,iA,iB,k);

				// Contravariant velocities
				m_dAuxDataNode(ConUaIx,i,j,k) =
					  dContraMetric2DA(iA,iB,0) * dCovUa
					+ dContraMetric2DA(iA,iB,1) * dCovUb;

				m_dAuxDataNode(ConUbIx,i,j,k) =
					  dContraMetric2DB(iA,iB,0) * dCovUa
					+ dContraMetric2DB(iA,iB,1) * dCovUb;

				// Specific kinetic energy plus pointwise pressure
				m_dAuxDataNode(KIx,i,j,k) = 0.5 * (
					  m_dAuxDataNode(ConUaIx,i,j,k) * dCovUa
					+ m_dAuxDataNode(ConUbIx,i,j,k) * dCovUb);

				m_dAuxDataNode(KIx,i,j,k) +=
					phys.GetG() * dataInitialNode(HIx,iA,iB,k);
			}
			}
			}

			// Update all elements
			for (int k = 0; k < nRElements; k++) {

				// Pointwise fluxes and pressure within spectral element
				for (int i = 0; i < nHorizontalOrder; i++) {
				for (int j = 0; j < nHorizontalOrder; j++) {

					const int iA = iElementA + i;
					const int iB = iElementB + j;

					// Height flux
					m_dAlphaMassFlux(i,j,k) =
						dJacobian2D(iA,iB)
						* (dataInitialNode(HIx,iA,iB,k)
							- dTopography(iA,iB))
						* m_dAuxDataNode(ConUaIx,i,j,k);

					m_dBetaMassFlux(i,j,k) =
						dJacobian2D(iA,iB)
						* (dataInitialNode(HIx,iA,iB,k)
							- dTopography(iA,iB))
						* m_dAuxDataNode(ConUbIx,i,j,k);

				}
				}

				// Pointwise update of quantities on model levels
				for (int i = 0; i < nHorizontalOrder; i++) {
				for (int j = 0; j < nHorizontalOrder; j++) {

					const int iA = iElementA + i;
					const int iB = iElementB + j;

					// Inverse Jacobian
					double dInvJacobian2D = 1.0 / dJacobian2D(iA,iB);

					// Derivatives of the covariant velocity field
					double dCovDaUb = 0.0;
					double dCovDbUa = 0.0;

					// Derivative of the kinetic energy
					double dDaKE = 0.0;
					double dDbKE = 0.0;

					// Aliases for alpha and beta velocities
					const double dConUa = m_dAuxDataNode(ConUaIx,i,j,k);
					const double dConUb = m_dAuxDataNode(ConUbIx,i,j,k);

					// Calculate derivatives in the alpha direction
					double dDaMassFluxA = 0.0;

#pragma unroll
					for (int s = 0; s < nHorizontalOrder; s++) {
#ifdef DIFFERENTIAL_FORM
						// Update density: Differential formulation
						dDaMassFluxA +=
							m_dAlphaMassFlux(s,j,k)
							* dDxBasis1D(s,i);
#else
						// Update density: Variational formulation
						dDaMassFluxA -=
							m_dAlphaMassFlux(s,j,k)
							* dStiffness1D(i,s);
#endif
						// Derivative of covariant beta velocity wrt alpha
						dCovDaUb +=
							dataInitialNode(VIx,iElementA+s,iB,k)
							* dDxBasis1D(s,i);

						// Derivative of specific kinetic energy wrt alpha
						dDaKE +=
							m_dAuxDataNode(KIx,s,j,k)
							* dDxBasis1D(s,i);
					}

					// Calculate derivatives in the beta direction
					double dDbMassFluxB = 0.0;

#pragma unroll
					for (int s = 0; s < nHorizontalOrder; s++) {
#ifdef DIFFERENTIAL_FORM
						// Update density: Differential formulation
						dDbMassFluxB +=
							m_dBetaMassFlux(i,s,k)
							* dDxBasis1D(s,j);

#else
						// Update density: Variational formulation
						dDbMassFluxB -=
							m_dBetaMassFlux(i,s,k)
							* dStiffness1D(j,s);
#endif
						// Derivative of covariant alpha velocity wrt beta
						dCovDbUa +=
							dataInitialNode(UIx,iA,iElementB+s,k)
							* dDxBasis1D(s,j);

						// Derivative of specific kinetic energy wrt beta
						dDbKE +=
							m_dAuxDataNode(KIx,i,s,k)
							* dDxBasis1D(s,j);
					}

					// Scale derivatives
					dDaMassFluxA *= dInvElementDeltaA;
					dCovDaUb     *= dInvElementDeltaA;
					dDaKE        *= dInvElementDeltaA;

					dDbMassFluxB *= dInvElementDeltaB;
					dCovDbUa     *= dInvElementDeltaB;
					dDbKE        *= dInvElementDeltaB;

					// Pointwise horizontal momentum update
					double dLocalUpdateUa = 0.0;
					double dLocalUpdateUb = 0.0;

					// Relative vorticity
					double dZetaXi = (dCovDaUb - dCovDbUa);

					// Rotational terms (covariant)
					double dCovUCrossZetaA  =   dConUb * dZetaXi;
					double dCovUCrossZetaB  = - dConUa * dZetaXi;

					// Coriolis terms
					dLocalUpdateUa +=
						dCoriolisF(iA,iB)
						* dJacobian2D(iA,iB)
						* dConUb;

					dLocalUpdateUb -=
						dCoriolisF(iA,iB)
						* dJacobian2D(iA,iB)
						* dConUa;

					// Horizontal updates
					dLocalUpdateUa += - dDaKE + dCovUCrossZetaA;
					dLocalUpdateUb += - dDbKE + dCovUCrossZetaB;

 					// Apply update to horizontal velocity
					dataUpdateNode(UIx,iA,iB,k) +=
						dDeltaT * dLocalUpdateUa;
					dataUpdateNode(VIx,iA,iB,k) +=
						dDeltaT * dLocalUpdateUb;

					// Update height
					dataUpdateNode(HIx,iA,iB,k) -=
						dDeltaT * dInvJacobian2D * (
							  dDaMassFluxA
							+ dDbMassFluxB);
				}
				}
			}
		}
		}
	}
}
/*
///////////////////////////////////////////////////////////////////////////////

inline void Diff4(
	double * const dataIn,
	const size_t jStrideIn,
	const size_t iStrideIn,
	double * dataOut,
	const size_t jStrideOut,
	const size_t iStrideOut
) {
	// Alpha derivatives
	for (int j = 0; j < 4; j++) {
		double * const dataInOff = dataIn + jStride*j;

		dataOut[0] =
			- 6.0 * dataInOff[jStride*j]
			- 8.09016994374947451262869435595E0 * dataInOff[iStride+jStride*j]
			+ 3.09016994374947451262869435595E0 * dataInOff[2*iStride]
			- 1.0 * dataInOff[3*iStride];

		dataOut[iStrideOut] =
			+ 1.61803398874989490252573887119E0 * dataInOff[0]
			- 2.23606797749978980505147774238E0 * dataInOff[2*iStride]
			+ 6.18033988749894902525738871191E-1 * dataInOff[3*iStride];

		dataOut[2*iStrideOut] =
			- 6.18033988749894902525738871191E-1 * dataInOff[0]

		m_dStiffness1D(0,0) = -6.0;
		m_dStiffness1D(1,0) =  1.61803398874989490252573887119E0;
		m_dStiffness1D(2,0) = -6.18033988749894902525738871191E-1;
		m_dStiffness1D(3,0) =  1.0;

		m_dStiffness1D(0,1) = -8.09016994374947451262869435595E0;
		m_dStiffness1D(1,1) =  0.0;
		m_dStiffness1D(2,1) =  2.23606797749978980505147774238E0;
		m_dStiffness1D(3,1) = -3.09016994374947451262869435595E0;

		m_dStiffness1D(0,2) =  3.09016994374947451262869435595E0;
		m_dStiffness1D(1,2) = -2.23606797749978980505147774238E0;
		m_dStiffness1D(2,2) =  0.0;
		m_dStiffness1D(3,2) =  8.09016994374947451262869435595E0;

		m_dStiffness1D(0,3) = -1.0;
		m_dStiffness1D(1,3) =  6.18033988749894902525738871191E-1;
		m_dStiffness1D(2,3) = -1.61803398874989490252573887119E0;
		m_dStiffness1D(3,3) =  6.0;

}
*/
///////////////////////////////////////////////////////////////////////////////

void HorizontalDynamicsFEM::StepNonhydrostaticPrimitive(
	int iDataInitial,
	int iDataUpdate,
	const Time & time,
	double dDeltaT
) {
	// Start the function timer
	FunctionTimer timer("HorizontalStepNonhydrostaticPrimitive");

	// Get a copy of the GLL grid
	GridGLL * pGrid = dynamic_cast<GridGLL*>(m_model.GetGrid());

	// Number of vertical elements

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
	const int UCrossZetaAIx = 5;
	const int UCrossZetaBIx = 6;
	const int UCrossZetaXIx = 7;
	const int ExnerIx = 8;

#if defined(FIXED_HORIZONTAL_ORDER)
	const int nHorizontalOrder = FIXED_HORIZONTAL_ORDER;
	if (nHorizontalOrder != m_nHorizontalOrder) {
		_EXCEPTIONT("Command line order must match FIXED_HORIZONTAL_ORDER");
	}
#else
	const int nHorizontalOrder = m_nHorizontalOrder;
#endif

#if defined(FIXED_RELEMENTS)
	const int nRElements = FIXED_RELEMENTS;
	if (nRElements != pGrid->GetRElements()) {
		_EXCEPTIONT("Command line levels must match FIXED_RELEMENTS");
	}
#else
	const int nRElements = pGrid->GetRElements();
#endif

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
		const DataArray4D<double> & dataInitialNode =
			pPatch->GetDataState(iDataInitial, DataLocation_Node);

		DataArray4D<double> & dataUpdateNode =
			pPatch->GetDataState(iDataUpdate, DataLocation_Node);

		const DataArray4D<double> & dataInitialREdge =
			pPatch->GetDataState(iDataInitial, DataLocation_REdge);

		DataArray4D<double> & dataUpdateREdge =
			pPatch->GetDataState(iDataUpdate, DataLocation_REdge);

		const DataArray4D<double> & dataInitialTracer =
			pPatch->GetDataTracers(iDataInitial);

		DataArray4D<double> & dataUpdateTracer =
			pPatch->GetDataTracers(iDataUpdate);

		DataArray4D<double> & m_dataXiDiffNode =
			pPatch->GetXiDiffNode(DataLocation_Node);

		// Number of tracers
		const int nTracerCount = dataInitialTracer.GetSize(0);

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

			const int iElementA = a * nHorizontalOrder + box.GetHaloElements();
			const int iElementB = b * nHorizontalOrder + box.GetHaloElements();

#if defined(FIX_ELEMENT_MASS_NONHYDRO)
			// Zero mass fixer arrays
			m_dElementMassFluxA.Zero();
			m_dElementMassFluxB.Zero();
			m_dElTotalArea.Zero();
#endif

			// Store 2D Jacobian
			for (int i = 0; i < nHorizontalOrder; i++) {
			for (int j = 0; j < nHorizontalOrder; j++) {

				const int iA = iElementA + i;
				const int iB = iElementB + j;

				m_dLocalCoriolisF(i,j) =
					dCoriolisF(iA,iB);
				m_dLocalJacobian2D(i,j) =
					dJacobian2D(iA,iB);
			}
			}

			// Compute auxiliary data in element
			for (int i = 0; i < nHorizontalOrder; i++) {
			for (int j = 0; j < nHorizontalOrder; j++) {
			for (int k = 0; k < nRElements; k++) {

				const int iA = iElementA + i;
				const int iB = iElementB + j;

				// Contravariant velocities
				const double dCovUa = dataInitialNode(UIx,iA,iB,k);
				const double dCovUb = dataInitialNode(VIx,iA,iB,k);
				const double dCovUx = dataInitialNode(WIx,iA,iB,k);

				// Store metric quantities
				m_dLocalJacobian(i,j,k) =
					dJacobian(iA,iB,k);

				m_dLocalDerivR(i,j,k,0) =
					dDerivRNode(iA,iB,k,0);
				m_dLocalDerivR(i,j,k,1) =
					dDerivRNode(iA,iB,k,1);
				m_dLocalDerivR(i,j,k,2) =
					dDerivRNode(iA,iB,k,2);

				m_dLocalContraMetric(i,j,k,0) =
					dContraMetricA(iA,iB,k,0);
				m_dLocalContraMetric(i,j,k,1) =
					dContraMetricA(iA,iB,k,1);
				m_dLocalContraMetric(i,j,k,2) =
					dContraMetricA(iA,iB,k,2);
				m_dLocalContraMetric(i,j,k,3) =
					dContraMetricB(iA,iB,k,1);
				m_dLocalContraMetric(i,j,k,4) =
					dContraMetricB(iA,iB,k,2);
				m_dLocalContraMetric(i,j,k,5) =
					dContraMetricXi(iA,iB,k,2);

				// Calculate covariant xi velocity and store
				m_dAuxDataNode(CovUxIx,i,j,k) = dCovUx;

				// Contravariant velocities
				m_dAuxDataNode(ConUaIx,i,j,k) =
					  m_dLocalContraMetric(i,j,k,0) * dCovUa
					+ m_dLocalContraMetric(i,j,k,1) * dCovUb
					+ m_dLocalContraMetric(i,j,k,2) * dCovUx;

				m_dAuxDataNode(ConUbIx,i,j,k) =
					  m_dLocalContraMetric(i,j,k,1) * dCovUa
					+ m_dLocalContraMetric(i,j,k,3) * dCovUb
					+ m_dLocalContraMetric(i,j,k,4) * dCovUx;

				m_dAuxDataNode(ConUxIx,i,j,k) =
					  m_dLocalContraMetric(i,j,k,2) * dCovUa
					+ m_dLocalContraMetric(i,j,k,4) * dCovUb
					+ m_dLocalContraMetric(i,j,k,5) * dCovUx;

				// Specific kinetic energy
				m_dAuxDataNode(KIx,i,j,k) = 0.5 * (
					  m_dAuxDataNode(ConUaIx,i,j,k) * dCovUa
					+ m_dAuxDataNode(ConUbIx,i,j,k) * dCovUb
					+ m_dAuxDataNode(ConUxIx,i,j,k) * dCovUx);

#ifdef FORMULATION_RHOTHETA_P
				// NOTE: For some reason using parenthetical notation on the
				//       LHS of these expressions crash the Intel compiler
				//       (observed with icpc 16.0.3)

				// Pressure
				m_dAuxDataNode[ExnerIx][i][j][k] =
					phys.PressureFromRhoTheta(
						dataInitialNode(PIx,iA,iB,k));
#endif
#ifdef FORMULATION_RHOTHETA_PI
				// Exner pressure
				m_dAuxDataNode[ExnerIx][i][j][k] =
					phys.ExnerPressureFromRhoTheta(
						dataInitialNode(PIx,iA,iB,k));
#endif
#if defined(FORMULATION_THETA) || defined(FORMULATION_THETA_FLUX)
				// Exner pressure
				m_dAuxDataNode[ExnerIx][i][j][k] =
					phys.ExnerPressureFromRhoTheta(
						  dataInitialNode(RIx,iA,iB,k)
						* dataInitialNode(PIx,iA,iB,k));
#endif

			}
			}
			}

			// Compute U cross Relative vorticity
			for (int i = 0; i < nHorizontalOrder; i++) {
			for (int j = 0; j < nHorizontalOrder; j++) {
			for (int k = 0; k < nRElements; k++) {

				const int iA = iElementA + i;
				const int iB = iElementB + j;

				// Vertical derivatives of the covariant velocity field
				const double dCovDxUa =
					pGrid->DifferentiateNodeToNode(
						&(dataInitialNode(UIx,iA,iB,0)),
						k, 1);

				const double dCovDxUb =
					pGrid->DifferentiateNodeToNode(
						&(dataInitialNode(VIx,iA,iB,0)),
						k, 1);

				// Derivative needed for calculating relative vorticity
				double dCovDaUb = 0.0;
				double dCovDaUx = 0.0;
				double dCovDbUa = 0.0;
				double dCovDbUx = 0.0;

				for (int s = 0; s < nHorizontalOrder; s++) {

					// Derivative of covariant beta velocity wrt alpha
					dCovDaUb +=
						dataInitialNode(VIx,iElementA+s,iB,k)
						* dDxBasis1D(s,i);

					// Derivative of covariant xi velocity wrt alpha
					dCovDaUx +=
						m_dAuxDataNode(CovUxIx,s,j,k)
						* dDxBasis1D(s,i);

					// Derivative of covariant alpha velocity wrt beta
					dCovDbUa +=
						dataInitialNode(UIx,iA,iElementB+s,k)
						* dDxBasis1D(s,j);

					// Derivative of covariant xi velocity wrt beta
					dCovDbUx +=
						m_dAuxDataNode(CovUxIx,i,s,k)
						* dDxBasis1D(s,j);
				}

				dCovDaUb *= dInvElementDeltaA;
				dCovDaUx *= dInvElementDeltaA;
				dCovDbUa *= dInvElementDeltaB;
				dCovDbUx *= dInvElementDeltaB;

				// Contravariant velocities
				const double dConUa = m_dAuxDataNode(ConUaIx,i,j,k);
				const double dConUb = m_dAuxDataNode(ConUbIx,i,j,k);
				const double dConUx = m_dAuxDataNode(ConUxIx,i,j,k);

				// Relative vorticity (contravariant)
				const double dJZetaA = (dCovDbUx - dCovDxUb);
				const double dJZetaB = (dCovDxUa - dCovDaUx);
				const double dJZetaX = (dCovDaUb - dCovDbUa);

				// U cross Relative Vorticity (contravariant)
				m_dAuxDataNode(UCrossZetaAIx,i,j,k) =
					dConUb * dJZetaX - dConUx * dJZetaB;

				m_dAuxDataNode(UCrossZetaBIx,i,j,k) =
					dConUx * dJZetaA - dConUa * dJZetaX;

				m_dAuxDataNode(UCrossZetaXIx,i,j,k) =
					- dConUa * dCovDaUx - dConUb * dCovDbUx;

				// Store the vertical derivative bits...
				m_dataXiDiffNode(0, iA, iB, k) =
					- dConUx * dCovDxUa;
				m_dataXiDiffNode(1, iA, iB, k) =
					  dConUx * dCovDxUb;
			}
			}
			}

			// Pointwise fluxes and pressure within spectral element
			for (int i = 0; i < nHorizontalOrder; i++) {
			for (int j = 0; j < nHorizontalOrder; j++) {
			for (int k = 0; k < nRElements; k++) {

				const int iA = iElementA + i;
				const int iB = iElementB + j;

				// Base fluxes (area times velocity)
				const double dAlphaBaseFlux =
					m_dLocalJacobian(i,j,k)
					* m_dAuxDataNode(ConUaIx,i,j,k);

				const double dBetaBaseFlux =
					m_dLocalJacobian(i,j,k)
					* m_dAuxDataNode(ConUbIx,i,j,k);

				// Density flux
				m_dAlphaMassFlux(i,j,k) =
					  dAlphaBaseFlux
					* dataInitialNode(RIx,iA,iB,k);

				m_dBetaMassFlux(i,j,k) =
					  dBetaBaseFlux
					* dataInitialNode(RIx,iA,iB,k);

#ifdef FORMULATION_PRESSURE
				// Pressure flux
				m_dAlphaPressureFlux(i,j,k) =
					  dAlphaBaseFlux
					* phys.GetGamma()
					* dataInitialNode(PIx,iA,iB,k);

				m_dBetaPressureFlux(i,j,k) =
					  dBetaBaseFlux
					* phys.GetGamma()
					* dataInitialNode(PIx,iA,iB,k);
#endif
#if defined(FORMULATION_RHOTHETA_PI) \
 || defined(FORMULATION_RHOTHETA_P)
				// RhoTheta flux
				m_dAlphaPressureFlux(i,j,k) =
					  dAlphaBaseFlux
					* dataInitialNode(PIx,iA,iB,k);

				m_dBetaPressureFlux(i,j,k) =
					  dBetaBaseFlux
					* dataInitialNode(PIx,iA,iB,k);
#endif
/*
#pragma unroll
				for (int c = 0; c < nTracerCount; c++) {
					m_dAlphaTracerFlux(c,i,j,k) =
						dAlphaBaseFlux
						* dataInitialTracer(c,iA,iB,k);

					m_dBetaTracerFlux(c,i,j,k) =
						dBetaBaseFlux
						* dataInitialTracer(c,iA,iB,k);
				}
*/
#if !defined(DISABLE_UNIFORM_DIFFUSION_CHECKS)
				////////////////////////////////////////////////////////
				// Apply uniform diffusion to tracers
				if (pGrid->HasUniformDiffusion()) {

					for (int c = 0; c < nTracerCount; c++) {

						// Derivatives of tracer mixing ratio
						double dCovDaQ = 0.0;
						double dCovDbQ = 0.0;

#pragma unroll
						for (int s = 0; s < nHorizontalOrder; s++) {
							dCovDaQ +=
								dataInitialTracer(c,iElementA+s,iB,k)
								/ dataInitialNode(RIx,iElementA+s,iB,k)
								* dDxBasis1D(s,i);

							dCovDbQ +=
								dataInitialTracer(c,iA,iElementB+s,k)
								/ dataInitialNode(RIx,iA,iElementB+s,k)
								* dDxBasis1D(s,j);
						}

						dCovDaQ *= dInvElementDeltaA;
						dCovDbQ *= dInvElementDeltaB;

						// Gradient of tracer mixing ratio
						double dConDaQ =
							  m_dLocalContraMetric(i,j,k,0) * dCovDaQ
							+ m_dLocalContraMetric(i,j,k,1) * dCovDbQ;

						double dConDbQ =
							  m_dLocalContraMetric(i,j,k,1) * dCovDaQ
							+ m_dLocalContraMetric(i,j,k,3) * dCovDbQ;

						m_dAlphaTracerFlux(c,i,j,k) -=
							pGrid->GetScalarUniformDiffusionCoeff()
							* m_dLocalJacobian(i,j,k)
							* dataInitialNode(RIx,iA,iB,k)
							* dConDaQ;

						m_dBetaTracerFlux(c,i,j,k) -=
							pGrid->GetScalarUniformDiffusionCoeff()
							* m_dLocalJacobian(i,j,k)
							* dataInitialNode(RIx,iA,iB,k)
							* dConDbQ;
					}
				}
#endif

#if defined(INSTEP_DIVERGENCE_DAMPING)
				// Derivatives of J U^i
				double dDaJUa = 0.0;
				double dDbJUb = 0.0;

				for (int s = 0; s < nHorizontalOrder; s++) {
					// Alpha derivative of J U^a
					dDaJUa +=
						m_dLocalJacobian(s,j,k)
						* m_dAuxDataNode(ConUaIx,s,j,k)
						* dDxBasis1D(s,i);

					// Beta derivative of J U^b
					dDbJUb +=
						m_dLocalJacobian(i,s,k)
						* m_dAuxDataNode(ConUbIx,i,s,k)
						* dDxBasis1D(s,j);
				}

				dDaJUa *= dInvElementDeltaA;
				dDbJUb *= dInvElementDeltaB;

				m_dDivergence(i,j,k) =
					(dDaJUa + dDbJUb) / m_dLocalJacobian(i,j,k);
#endif
			}
			}
			}

			// Pointwise update of quantities on model levels
			for (int i = 0; i < nHorizontalOrder; i++) {
			for (int j = 0; j < nHorizontalOrder; j++) {
			for (int k = 0; k < nRElements; k++) {

				const int iA = iElementA + i;
				const int iB = iElementB + j;

				// Inverse Jacobian
				const double dInvJacobian =
					1.0 / m_dLocalJacobian(i,j,k);

				// Aliases for alpha and beta velocities
				const double dConUa = m_dAuxDataNode(ConUaIx,i,j,k);
				const double dConUb = m_dAuxDataNode(ConUbIx,i,j,k);
				const double dConUx = m_dAuxDataNode(ConUxIx,i,j,k);

				const double dCovUx = m_dAuxDataNode(CovUxIx,i,j,k);

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

				// Alpha derivatives
				for (int s = 0; s < nHorizontalOrder; s++) {

					// Update density: Variational formulation
					dDaRhoFluxA -=
						m_dAlphaMassFlux(s,j,k)
						* dStiffness1D(i,s);

					// Update pressure: Variational formulation
					dDaPressureFluxA -=
						m_dAlphaPressureFlux(s,j,k)
						* dStiffness1D(i,s);

#if defined(FORMULATION_PRESSURE)
					// Derivative of pressure with respect to alpha
					dDaP +=
						dataInitialNode(PIx,iElementA+s,iB,k)
						* dDxBasis1D(s,i);
#endif
#if defined(FORMULATION_RHOTHETA_PI) \
 || defined(FORMULATION_RHOTHETA_P) \
 || defined(FORMULATION_THETA) \
 || defined(FORMULATION_THETA_FLUX)
					// Derivative of (Exner) pressure with respect to alpha
					dDaP +=
						m_dAuxDataNode(ExnerIx,s,j,k)
						* dDxBasis1D(s,i);
#endif

					// Derivative of specific kinetic energy wrt alpha
					dDaKE +=
						m_dAuxDataNode(KIx,s,j,k)
						* dDxBasis1D(s,i);

#if defined(INSTEP_DIVERGENCE_DAMPING)
					dDaDiv -=
						m_dDivergence(s,j,k)
						* dStiffness1D(i,s);
#endif
				}

				// Beta derivatives
				for (int s = 0; s < nHorizontalOrder; s++) {

					// Update density: Variational formulation
					dDbRhoFluxB -=
						m_dBetaMassFlux(i,s,k)
						* dStiffness1D(j,s);

					// Update pressure: Variational formulation
					dDbPressureFluxB -=
						m_dBetaPressureFlux(i,s,k)
						* dStiffness1D(j,s);

#if defined(FORMULATION_PRESSURE)
					// Derivative of pressure with respect to beta
					dDbP +=
						dataInitialNode(PIx,iA,iElementB+s,k)
						* dDxBasis1D(s,j);
#endif
#if defined(FORMULATION_RHOTHETA_PI) \
 || defined(FORMULATION_RHOTHETA_P) \
 || defined(FORMULATION_THETA) \
 || defined(FORMULATION_THETA_FLUX)
					// Derivative of (Exner) pressure with respect to beta
					dDbP +=
						m_dAuxDataNode(ExnerIx,i,s,k)
						* dDxBasis1D(s,j);
#endif

					// Derivative of specific kinetic energy wrt beta
					dDbKE +=
						m_dAuxDataNode(KIx,i,s,k)
						* dDxBasis1D(s,j);

#if defined(INSTEP_DIVERGENCE_DAMPING)
					dDbDiv -=
						m_dDivergence(i,s,k)
						* dStiffness1D(j,s);
#endif
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

#if defined(INSTEP_DIVERGENCE_DAMPING)
				dDaDiv *= dInvElementDeltaA;
				dDbDiv *= dInvElementDeltaB;
#endif

				// Pointwise momentum updates
				double dLocalUpdateUa = 0.0;
				double dLocalUpdateUb = 0.0;

				// Updates due to rotational terms
				dLocalUpdateUa += m_dAuxDataNode(UCrossZetaAIx,i,j,k);
				dLocalUpdateUb += m_dAuxDataNode(UCrossZetaBIx,i,j,k);

				// Coriolis terms
				dLocalUpdateUa +=
					m_dLocalCoriolisF(i,j)
					* m_dLocalJacobian2D(i,j)
					* dConUb;

				dLocalUpdateUb -=
					m_dLocalCoriolisF(i,j)
					* m_dLocalJacobian2D(i,j)
					* dConUa;

				// Pressure gradient force
#if defined(FORMULATION_PRESSURE) || defined(FORMULATION_RHOTHETA_P)
				const double dPressureGradientForceUa =
					dDaP / dataInitialNode(RIx,iA,iB,k);
				const double dPressureGradientForceUb =
					dDbP / dataInitialNode(RIx,iA,iB,k);
#endif
#if defined(FORMULATION_RHOTHETA_PI)
				const double dPressureGradientForceUa =
					dDaP * dataInitialNode(PIx,iA,iB,k)
					/ dataInitialNode(RIx,iA,iB,k);
				const double dPressureGradientForceUb =
					dDbP * dataInitialNode(PIx,iA,iB,k)
					/ dataInitialNode(RIx,iA,iB,k);
#endif
#if defined(FORMULATION_THETA) || defined(FORMULATION_THETA_FLUX)
				const double dPressureGradientForceUa =
					dDaP * dataInitialNode(PIx,iA,iB,k);
				const double dPressureGradientForceUb =
					dDbP * dataInitialNode(PIx,iA,iB,k);
#endif

				// Gravity
				const double dDaPhi = phys.GetG() * m_dLocalDerivR(i,j,k,0);
				const double dDbPhi = phys.GetG() * m_dLocalDerivR(i,j,k,1);

				// Horizontal updates due to gradient terms
				const double dDaUpdate =
					dPressureGradientForceUa + dDaKE + dDaPhi;

				const double dDbUpdate =
					dPressureGradientForceUb + dDbKE + dDbPhi;

				// Apply gradient term update to total update
				dLocalUpdateUa -= dDaUpdate;
				dLocalUpdateUb -= dDbUpdate;

 				// Apply update to horizontal velocity on model levels
				dataUpdateNode(UIx,iA,iB,k) +=
					dDeltaT * dLocalUpdateUa;

				// Omit beta update for XZ 2D models
				if (pGrid->GetIsCartesianXZ() == false) {
					dataUpdateNode(VIx,iA,iB,k) +=
						dDeltaT * dLocalUpdateUb;
				}

#if defined(INSTEP_DIVERGENCE_DAMPING)
				// Apply instep divergence
				dataUpdateNode(UIx,iA,iB,k) +=
					dDeltaT * m_dInstepNuDiv * dDaDiv;
				if (pGrid->GetIsCartesianXZ() == false) {
					dataUpdateNode(VIx,iA,iB,k) +=
						dDeltaT * m_dInstepNuDiv * dDbDiv;
				}
#endif

#if !defined(FIX_ELEMENT_MASS_NONHYDRO)
				// Update density on model levels
				dataUpdateNode(RIx,iA,iB,k) -=
					dDeltaT * dInvJacobian * (
						  dDaRhoFluxA
						+ dDbRhoFluxB);
#endif

#if defined(FORMULATION_PRESSURE)
				// Update pressure on model levels
				dataUpdateNode(PIx,iA,iB,k) +=
					dDeltaT * (phys.GetGamma() - 1.0)
					* (dConUa * dDaP + dConUb * dDbP);

				dataUpdateNode(PIx,iA,iB,k) -=
					dDeltaT * dInvJacobian * (
						  dDaPressureFluxA
						+ dDbPressureFluxB);
#endif
#if defined(FORMULATION_RHOTHETA_PI) || defined(FORMULATION_RHOTHETA_P)
				// Update RhoTheta on model levels
				dataUpdateNode(PIx,iA,iB,k) -=
					dDeltaT * dInvJacobian * (
						  dDaPressureFluxA
						+ dDbPressureFluxB);
#endif

				// Update vertical velocity on nodes
				if (pGrid->GetVarLocation(WIx) == DataLocation_Node) {

					// Calculate vertical velocity update
					double dLocalUpdateUr =
						m_dAuxDataNode(UCrossZetaXIx,i,j,k);

					if (k == 0) {
						dLocalUpdateUr =
							- ( m_dLocalContraMetric(iA,iB,0,2)
									* dLocalUpdateUa
							  + m_dLocalContraMetric(iA,iB,0,4)
							  		* dLocalUpdateUb)
							/ m_dLocalContraMetric(iA,iB,0,5);

					} else if (k == nRElements-1) {
						dLocalUpdateUr = 0.0;
					}

					// Update vertical velocity
					dataUpdateNode(WIx,iA,iB,k) +=
						dDeltaT * dLocalUpdateUr;
				}

#if defined(FORMULATION_THETA)
				// Update thermodynamic variable on nodes
				if (pGrid->GetVarLocation(PIx) == DataLocation_Node) {

					// Derivatives of the theta field
					double dDaTheta = 0.0;
					double dDbTheta = 0.0;

#pragma unroll
					for (int s = 0; s < nHorizontalOrder; s++) {
						dDaTheta +=
							dataInitialNode(PIx,iElementA+s,iB,k)
							* dDxBasis1D(s,i);

						dDbTheta +=
							dataInitialNode(PIx,iA,iElementB+s,k)
							* dDxBasis1D(s,j);
					}

					dDaTheta *= dInvElementDeltaA;
					dDbTheta *= dInvElementDeltaB;

					// Update Theta on model levels
					dataUpdateNode(PIx,iA,iB,k) -=
						dDeltaT * (dConUa * dDaTheta + dConUb * dDbTheta);
				}
#endif
#if defined(FORMULATION_THETA_FLUX)
				// Update thermodynamic variable on nodes
				if (pGrid->GetVarLocation(PIx) == DataLocation_Node) {

					// Derivatives of the theta field
					double dDaJUa = 0.0;
					double dDbJUb = 0.0;

					double dDaJThetaUa = 0.0;
					double dDbJThetaUb = 0.0;

#pragma unroll
					for (int s = 0; s < nHorizontalOrder; s++) {
						dDaJUa +=
							m_dLocalJacobian(s,j,k)
							* m_dAuxDataNode(ConUaIx,s,j,k)
							* dDxBasis1D(s,i);

						dDbJUb +=
							m_dLocalJacobian(i,s,k)
							* m_dAuxDataNode(ConUbIx,i,s,k)
							* dDxBasis1D(s,j);

						dDaJThetaUa +=
							m_dLocalJacobian(s,j,k)
							* dataInitialNode(PIx,iElementA+s,iB,k)
							* m_dAuxDataNode(ConUaIx,s,j,k)
							* dDxBasis1D(s,i);

						dDbJThetaUb +=
							m_dLocalJacobian(i,s,k)
							* dataInitialNode(PIx,iA,iElementB+s,k)
							* m_dAuxDataNode(ConUbIx,i,s,k)
							* dDxBasis1D(s,j);
					}

					dDaJUa *= dInvElementDeltaA;
					dDbJUb *= dInvElementDeltaB;

					dDaJThetaUa *= dInvElementDeltaA;
					dDbJThetaUb *= dInvElementDeltaB;

					// Update Theta on model levels
					const double dUpdateTheta =
						(dDaJThetaUa + dDbJThetaUb)
						- dataInitialNode(PIx,iA,iB,k)
							* (dDaJUa + dDbJUb);

					dataUpdateNode(PIx,iA,iB,k) -=
						dDeltaT
						* dUpdateTheta
						* dInvJacobian;
				}
#endif

				// Update tracers
				for (int c = 0; c < nTracerCount; c++) {
					double dDaTracerFluxA = 0.0;
					double dDbTracerFluxB = 0.0;

#pragma unroll
					for (int s = 0; s < nHorizontalOrder; s++) {
						dDaTracerFluxA -=
							m_dAlphaTracerFlux(c,s,j,k)
							* dStiffness1D(i,s);

						dDbTracerFluxB -=
							m_dBetaTracerFlux(c,i,s,k)
							* dStiffness1D(j,s);
					}

					dDaTracerFluxA *= dInvElementDeltaA;
					dDbTracerFluxB *= dInvElementDeltaB;

					dataUpdateTracer(c,iA,iB,k) -=
						dDeltaT * dInvJacobian * (
							  dDaTracerFluxA
							+ dDbTracerFluxB);
				}

#if defined(FIX_ELEMENT_MASS_NONHYDRO)
				// Integrate element mass fluxes
				m_dElementMassFluxA[k] +=
					dInvJacobian
					* dDaRhoFluxA
					* dElementAreaNode(iA,iB,k);

				m_dElementMassFluxB[k] +=
					dInvJacobian
					* dDbRhoFluxB
					* dElementAreaNode(iA,iB,k);

				m_dElTotalArea[k] +=
					dElementAreaNode(iA,iB,k);

				// Store the local element fluxes
				m_dAlphaElMassFlux(i,j,k) = dDaRhoFluxA;
				m_dBetaElMassFlux(i,j,k) = dDbRhoFluxB;
#endif
			}
			}
			}
#if defined(FIX_ELEMENT_MASS_NONHYDRO)
			// Apply mass fixer
			for (int k = 0; k < nRElements; k++) {
				m_dElementMassFluxA[k] /= m_dElTotalArea[k];
				m_dElementMassFluxB[k] /= m_dElTotalArea[k];
			}

			for (int i = 0; i < nHorizontalOrder; i++) {
			for (int j = 0; j < nHorizontalOrder; j++) {
			for (int k = 0; k < nRElements; k++) {

				const int iA = iElementA + i;
				const int iB = iElementB + j;

				const double dInvJacobian =
					1.0 / m_dLocalJacobian(i,j,k);
				const double dJacobian =
					m_dLocalJacobian(i,j,k);

				m_dAlphaElMassFlux(i,j,k) -=
					dJacobian * m_dElementMassFluxA[k];
				m_dBetaElMassFlux(i,j,k) -=
					dJacobian * m_dElementMassFluxB[k];

				// Update density on model levels
				dataUpdateNode(RIx,iA,iB,k) -=
					dDeltaT * dInvJacobian * (
						  m_dAlphaElMassFlux(i,j,k)
						+ m_dBetaElMassFlux(i,j,k));
			}
			}
			}
#endif

			// Update vertical velocity on interfaces
			if (pGrid->GetVarLocation(WIx) == DataLocation_REdge) {

				// Update vertical velocity at bottom boundary
				for (int i = 0; i < nHorizontalOrder; i++) {
				for (int j = 0; j < nHorizontalOrder; j++) {

					const int iA = iElementA + i;
					const int iB = iElementB + j;

					// Interpolate horizontal velocity to bottom boundary
					const double dU0 =
						pGrid->InterpolateNodeToREdge(
							&(dataUpdateNode(UIx,iA,iB,0)),
							NULL, 0, 0.0, 1);

					const double dV0 =
						pGrid->InterpolateNodeToREdge(
							&(dataUpdateNode(VIx,iA,iB,0)),
							NULL, 0, 0.0, 1);

					// Update vertical velocity on boundary
					dataUpdateREdge(WIx,iA,iB,0) =
						- ( dContraMetricXiREdge(iA,iB,0,0) * dU0
						  + dContraMetricXiREdge(iA,iB,0,1) * dV0)
							/ dContraMetricXiREdge(iA,iB,0,2);

				}
				}

				// Update interior interfaces
				for (int i = 0; i < nHorizontalOrder; i++) {
				for (int j = 0; j < nHorizontalOrder; j++) {
				for (int k = 1; k < nRElements; k++) {

					const int iA = iElementA + i;
					const int iB = iElementB + j;

					// Interpolate U cross Zeta to interfaces
					const double dUCrossZetaX =
						pGrid->InterpolateNodeToREdge(
							&(m_dAuxDataNode(UCrossZetaXIx,i,j,0)),
							NULL, k, 0.0, 1);

					// Calculate vertical velocity update
					dataUpdateREdge(WIx,iA,iB,k) +=
						dDeltaT * dUCrossZetaX;
				}
				}
				}
			}

#if defined(FORMULATION_THETA) || defined(FORMULATION_THETA_FLUX)
			// Update thermodynamic variable on interfaces
			if (pGrid->GetVarLocation(PIx) == DataLocation_REdge) {

				for (int i = 0; i < nHorizontalOrder; i++) {
				for (int j = 0; j < nHorizontalOrder; j++) {
				for (int k = 0; k <= nRElements; k++) {

					const int iA = iElementA + i;
					const int iB = iElementB + j;

					// Contravariant velocities
					const double dCovUa = dataInitialREdge(UIx,iA,iB,k);
					const double dCovUb = dataInitialREdge(VIx,iA,iB,k);
					const double dCovUx = dataInitialREdge(WIx,iA,iB,k);

					// Contravariant velocities on interfaces
					m_dAuxDataREdge(ConUaIx,i,j,k) =
						  dContraMetricAREdge(iA,iB,k,0) * dCovUa
						+ dContraMetricAREdge(iA,iB,k,1) * dCovUb
						+ dContraMetricAREdge(iA,iB,k,2) * dCovUx;

					m_dAuxDataREdge(ConUbIx,i,j,k) =
						  dContraMetricBREdge(iA,iB,k,0) * dCovUa
						+ dContraMetricBREdge(iA,iB,k,1) * dCovUb
						+ dContraMetricBREdge(iA,iB,k,2) * dCovUx;
				}
				}
				}

				for (int i = 0; i < nHorizontalOrder; i++) {
				for (int j = 0; j < nHorizontalOrder; j++) {
				for (int k = 0; k <= nRElements; k++) {

					const int iA = iElementA + i;
					const int iB = iElementB + j;

#if defined(FORMULATION_THETA)
					// Derivatives of the theta field on interfaces
					double dDaTheta = 0.0;
					double dDbTheta = 0.0;

#pragma unroll
					for (int s = 0; s < nHorizontalOrder; s++) {
						dDaTheta +=
							dataInitialREdge(PIx,iElementA+s,iB,k)
							* dDxBasis1D(s,i);

						dDbTheta +=
							dataInitialREdge(PIx,iA,iElementB+s,k)
							* dDxBasis1D(s,j);
					}

					dDaTheta *= dInvElementDeltaA;
					dDbTheta *= dInvElementDeltaB;

					// Update Theta on interfaces
					const double dConUa = m_dAuxDataREdge(ConUaIx,i,j,k);
					const double dConUb = m_dAuxDataREdge(ConUbIx,i,j,k);

					dataUpdateREdge(PIx,iA,iB,k) -=
						dDeltaT * (dConUa * dDaTheta + dConUb * dDbTheta);
#endif
#if defined(FORMULATION_THETA_FLUX)
					// Derivatives of the theta field
					double dDaJUa = 0.0;
					double dDbJUb = 0.0;

					double dDaJThetaUa = 0.0;
					double dDbJThetaUb = 0.0;

#pragma unroll
					for (int s = 0; s < nHorizontalOrder; s++) {
						dDaJUa +=
							dJacobianREdge(iElementA+s,iB,k)
							* m_dAuxDataREdge(ConUaIx,s,j,k)
							* dDxBasis1D(s,i);

						dDbJUb +=
							dJacobianREdge(iA,iElementB+s,k)
							* m_dAuxDataREdge(ConUbIx,i,s,k)
							* dDxBasis1D(s,j);

						dDaJThetaUa +=
							dJacobianREdge(iElementA+s,iB,k)
							* dataInitialREdge(PIx,iElementA+s,iB,k)
							* m_dAuxDataREdge(ConUaIx,s,j,k)
							* dDxBasis1D(s,i);

						dDbJThetaUb +=
							dJacobianREdge(iA,iElementB+s,k)
							* dataInitialREdge(PIx,iA,iElementB+s,k)
							* m_dAuxDataREdge(ConUbIx,i,s,k)
							* dDxBasis1D(s,j);
					}

					dDaJUa *= dInvElementDeltaA;
					dDbJUb *= dInvElementDeltaB;

					dDaJThetaUa *= dInvElementDeltaA;
					dDbJThetaUb *= dInvElementDeltaB;

					// Update Theta on model levels
					double dUpdateTheta =
						(dDaJThetaUa + dDbJThetaUb)
						- dataInitialREdge(PIx,iA,iB,k)
							* (dDaJUa + dDbJUb);

					dataUpdateREdge(PIx,iA,iB,k) -=
						dDeltaT
						* dUpdateTheta
						/ dJacobianREdge(iA,iB,k);
#endif
				}
				}
				}
			}
#endif
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

	// Get a copy of the GLL grid
	GridGLL * pGrid = dynamic_cast<GridGLL*>(m_model.GetGrid());

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

#if !defined(DISABLE_UNIFORM_DIFFUSION_CHECKS)
	// Uniform diffusion of U and V with vector diffusion coeff
	if (pGrid->HasUniformDiffusion()) {
		ApplyVectorHyperdiffusion(
			iDataInitial,
			iDataUpdate,
			dDeltaT,
			- pGrid->GetVectorUniformDiffusionCoeff(),
			- pGrid->GetVectorUniformDiffusionCoeff(),
			false);

		ApplyVectorHyperdiffusion(
			DATA_INDEX_REFERENCE,
			iDataUpdate,
			dDeltaT,
			pGrid->GetVectorUniformDiffusionCoeff(),
			pGrid->GetVectorUniformDiffusionCoeff(),
			false);

		if (eqn.GetType() == EquationSet::PrimitiveNonhydrostaticEquations) {

			// Uniform diffusion of Theta with scalar diffusion coeff
			ApplyScalarHyperdiffusion(
				iDataInitial,
				iDataUpdate,
				dDeltaT,
				pGrid->GetScalarUniformDiffusionCoeff(),
				false,
				2,
				true);

			// Uniform diffusion of W with vector diffusion coeff
			ApplyScalarHyperdiffusion(
				iDataInitial,
				iDataUpdate,
				dDeltaT,
				pGrid->GetVectorUniformDiffusionCoeff(),
				false,
				3,
				true);
		}
	}
#endif

#if defined(RESIDUAL_DIFFUSION_THERMO)
	// Residual diffusion of Theta with residual based diffusion coeff
	ApplyScalarHyperdiffusionResidual(
		iDataInitial,
		iDataUpdate,
		2,
		pGrid->GetScalarUniformDiffusionCoeff(),
		dDeltaT,
		2,
		true);
#endif
#if defined(RESIDUAL_DIFFUSION_RHO)
	// Residual diffusion of Rho with residual based diffusion coeff
	ApplyScalarHyperdiffusionResidual(
		iDataInitial,
		iDataUpdate,
		2,
		pGrid->GetScalarUniformDiffusionCoeff(),
		dDeltaT,
		4,
		true);
#endif

	// Apply positive definite filter to tracers
	FilterNegativeTracers(iDataUpdate);
}

///////////////////////////////////////////////////////////////////////////////

void HorizontalDynamicsFEM::ApplyScalarHyperdiffusion(
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

#if defined(FIXED_HORIZONTAL_ORDER)
	const int nHorizontalOrder = FIXED_HORIZONTAL_ORDER;
	if (nHorizontalOrder != m_nHorizontalOrder) {
		_EXCEPTIONT("Command line order must match FIXED_HORIZONTAL_ORDER");
	}
#else
	const int nHorizontalOrder = m_nHorizontalOrder;
#endif

#if defined(FIXED_RELEMENTS)
	const int nRElements = FIXED_RELEMENTS;
	if (nRElements != pGrid->GetRElements()) {
		_EXCEPTIONT("Command line levels must match FIXED_RELEMENTS");
	}
#else
	const int nRElements = pGrid->GetRElements();
#endif

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
				for (int a = 0; a < nElementCountA; a++) {
				for (int b = 0; b < nElementCountB; b++) {
				for (int k = 0; k < nElementCountR; k++) {

					int iElementA =
						a * nHorizontalOrder + box.GetHaloElements();
					int iElementB =
						b * nHorizontalOrder + box.GetHaloElements();

					// Store the buffer state
					for (int i = 0; i < nHorizontalOrder; i++) {
					for (int j = 0; j < nHorizontalOrder; j++) {
						int iA = iElementA + i;
						int iB = iElementB + j;

						m_dBufferState(i,j) =
							(*pDataInitial)(c,iA,iB,k);
					}
					}

					// Remove the reference state from the buffer state
					if (fRemoveRefState) {
						for (int i = 0; i < nHorizontalOrder; i++) {
						for (int j = 0; j < nHorizontalOrder; j++) {
							int iA = iElementA + i;
							int iB = iElementB + j;

							m_dBufferState(i,j) -=
								(*pDataRef)(c,iA,iB,k);
						}
						}
					}

					// Calculate the pointwise gradient of the scalar field
					for (int i = 0; i < nHorizontalOrder; i++) {
					for (int j = 0; j < nHorizontalOrder; j++) {
						int iA = iElementA + i;
						int iB = iElementB + j;

						double dDaPsi = 0.0;
						double dDbPsi = 0.0;

#pragma unroll
						for (int s = 0; s < nHorizontalOrder; s++) {
							dDaPsi +=
								m_dBufferState(s,j)
								* dDxBasis1D(s,i);

							dDbPsi +=
								m_dBufferState(i,s)
								* dDxBasis1D(s,j);
						}

						dDaPsi *= dInvElementDeltaA;
						dDbPsi *= dInvElementDeltaB;

						m_dJGradientA(i,j) =
							(*pJacobian)(iA,iB,k) * (
								+ dContraMetricA(iA,iB,0) * dDaPsi
								+ dContraMetricA(iA,iB,1) * dDbPsi);

						m_dJGradientB(i,j) =
							(*pJacobian)(iA,iB,k) * (
								+ dContraMetricB(iA,iB,0) * dDaPsi
								+ dContraMetricB(iA,iB,1) * dDbPsi);
					}
					}

					// Pointwise updates
					double dElementMassFluxA = 0.0;
					double dElementMassFluxB = 0.0;
					double dMassFluxPerNodeA = 0.0;
					double dMassFluxPerNodeB = 0.0;
					double dElTotalArea = 0.0;
					for (int i = 0; i < nHorizontalOrder; i++) {
					for (int j = 0; j < nHorizontalOrder; j++) {
						int iA = iElementA + i;
						int iB = iElementB + j;

						// Inverse Jacobian and Jacobian
						const double dInvJacobian = 1.0 / (*pJacobian)(iA,iB,k);

						// Compute integral term
						double dUpdateA = 0.0;
						double dUpdateB = 0.0;

#pragma unroll
						for (int s = 0; s < nHorizontalOrder; s++) {
							dUpdateA +=
								m_dJGradientA(s,j)
								* dStiffness1D(i,s);

							dUpdateB +=
								m_dJGradientB(i,s)
								* dStiffness1D(j,s);
						}

						dUpdateA *= dInvElementDeltaA;
						dUpdateB *= dInvElementDeltaB;

#ifdef FIX_ELEMENT_MASS_NONHYDRO
							if (c == RIx) {
		 						double dInvJacobian = 1.0 / (*pJacobian)(iA,iB,k);
								// Integrate element mass fluxes
								dElementMassFluxA += dInvJacobian *
									dUpdateA * dElementAreaNode(iA,iB,k);
								dElementMassFluxB += dInvJacobian *
									dUpdateB * dElementAreaNode(iA,iB,k);
								dElTotalArea += dElementAreaNode(iA,iB,k);

								// Store the local element fluxes
								m_dAlphaElMassFlux(i,j,k) = dUpdateA;
								m_dBetaElMassFlux(i,j,k) = dUpdateB;
							} else {
								// Apply update
								(*pDataUpdate)(c,iA,iB,k) -=
									dDeltaT * dInvJacobian * dLocalNu
										* (dUpdateA + dUpdateB);
							}
#else
						// Apply update
						(*pDataUpdate)(c,iA,iB,k) -=
							dDeltaT * dInvJacobian * dLocalNu
								* (dUpdateA + dUpdateB);
#endif
					}
					}
#ifdef FIX_ELEMENT_MASS_NONHYDRO
					if (c == RIx) {
						double dMassFluxPerNodeA = dElementMassFluxA / dElTotalArea;
						double dMassFluxPerNodeB = dElementMassFluxB / dElTotalArea;

						// Compute the fixed mass update
						for (int i = 0; i < nHorizontalOrder; i++) {
						for (int j = 0; j < nHorizontalOrder; j++) {

							int iA = iElementA + i;
							int iB = iElementB + j;

							// Inverse Jacobian
							const double dInvJacobian = 1.0 / (*pJacobian)(iA,iB,k);
							const double dJacobian = (*pJacobian)(iA,iB,k);

							m_dAlphaElMassFlux(i,j,k) -= dJacobian * dMassFluxPerNodeA;
							m_dBetaElMassFlux(i,j,k) -= dJacobian * dMassFluxPerNodeB;

							// Update fixed density on model levels
							(*pDataUpdate)(c,iA,iB,k) -=
								dDeltaT * dInvJacobian * dLocalNu *
									(m_dAlphaElMassFlux(i,j,k)
									+ m_dBetaElMassFlux(i,j,k));
						}
						}
					}
#endif
				}
				}
				}
			}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void HorizontalDynamicsFEM::ApplyScalarHyperdiffusionResidual(
	int iDataInitial,
	int iDataUpdate,
	int iDataResidual,
	double dNuRayleigh,
	double dDeltaT,
	int iComponent,
	bool fRemoveRefState
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
	const int CovUxIx = 3;
	const int KIx = 4;
	const int UCrossZetaAIx = 5;
	const int UCrossZetaBIx = 6;
	const int UCrossZetaXIx = 7;
	const int ExnerIx = 8;

	// Get a copy of the GLL grid
	GridGLL * pGrid = dynamic_cast<GridGLL*>(m_model.GetGrid());

#if defined(FIXED_HORIZONTAL_ORDER)
	const int nHorizontalOrder = FIXED_HORIZONTAL_ORDER;
	if (nHorizontalOrder != m_nHorizontalOrder) {
		_EXCEPTIONT("Command line order must match FIXED_HORIZONTAL_ORDER");
	}
#else
	const int nHorizontalOrder = m_nHorizontalOrder;
#endif

#if defined(FIXED_RELEMENTS)
	const int nRElements = FIXED_RELEMENTS;
	if (nRElements != pGrid->GetRElements()) {
		_EXCEPTIONT("Command line levels must match FIXED_RELEMENTS");
	}
#else
	const int nRElements = pGrid->GetRElements();
#endif

	// Physical constants
	const PhysicalConstants & phys = m_model.GetPhysicalConstants();

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
		const DataArray3D<double> & dContraMetric2DA =
			pPatch->GetContraMetric2DA();
		const DataArray3D<double> & dContraMetric2DB =
			pPatch->GetContraMetric2DB();

		const DataArray4D<double> & dDerivRNode =
			pPatch->GetDerivRNode();
		const DataArray4D<double> & dDerivRREdge =
			pPatch->GetDerivRREdge();

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

		// Grid data
		DataArray4D<double> & dataInitialNode =
			pPatch->GetDataState(iDataInitial, DataLocation_Node);

		DataArray4D<double> & dataInitialREdge =
			pPatch->GetDataState(iDataInitial, DataLocation_REdge);

		DataArray4D<double> & dataResidualNode =
			pPatch->GetDataResidual(iDataResidual, DataLocation_Node);

		DataArray4D<double> & dataResidualREdge =
			pPatch->GetDataResidual(iDataResidual, DataLocation_REdge);

		// Get the residual storage to update for output_SGS
		DataArray4D<double> & dataDynSGSNode =
			pPatch->GetDataResidualSGS(DataLocation_Node);

		DataArray4D<double> & dataDynSGSREdge =
			pPatch->GetDataResidualSGS(DataLocation_REdge);

		DataArray4D<double> & dataUpdateNode =
			pPatch->GetDataState(iDataUpdate, DataLocation_Node);

		DataArray4D<double> & dataUpdateREdge =
			pPatch->GetDataState(iDataUpdate, DataLocation_REdge);

		const DataArray4D<double> & dataRefNode =
			pPatch->GetReferenceState(DataLocation_Node);

		const DataArray4D<double> & dataRefREdge =
			pPatch->GetReferenceState(DataLocation_REdge);

		// Tracer data
#pragma Residual diffusion on tracers NOT IMPLEMENTED yet.
		DataArray4D<double> & dataInitialTracer =
			pPatch->GetDataTracers(iDataInitial);

		DataArray4D<double> & dataUpdateTracer =
			pPatch->GetDataTracers(iDataUpdate);

		DataArray4D<double> & dataResidualTracer =
			pPatch->GetDataTracers(iDataResidual);

		// Element grid spacing and derivative coefficients
		const double dElementDeltaA = pPatch->GetElementDeltaA();
		const double dElementDeltaB = pPatch->GetElementDeltaB();

		const double dElementLength = 0.5 *
			(dElementDeltaA + dElementDeltaB);

		const double dInvElementDeltaA = 1.0 / dElementDeltaA;
		const double dInvElementDeltaB = 1.0 / dElementDeltaB;

		const DataArray2D<double> & dDxBasis1D = pGrid->GetDxBasis1D();
		const DataArray2D<double> & dStiffness1D = pGrid->GetStiffness1D();

		// Number of finite elements
		int nElementCountA = pPatch->GetElementCountA();
		int nElementCountB = pPatch->GetElementCountB();

		// Initialize new hyperviscosity coefficient
		double dResNu[nHorizontalOrder][nHorizontalOrder];

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
			 	const DataArray4D<double> * pDataResidual;
				const DataArray4D<double> * pDataRef;
				const DataArray3D<double> * pJacobian;
				const DataArray4D<double> * pDerivR;
				const DataArray4D<double> * pContraMetricA;
				const DataArray4D<double> * pContraMetricB;
				const DataArray4D<double> * pContraMetricXi;

				if (iType == 0) {
					if (pGrid->GetVarLocation(c) == DataLocation_Node) {
						pDataInitial = &dataInitialNode;
						pDataUpdate  = &dataUpdateNode;
						pDataResidual = &dataResidualNode;
						pDataRef = &dataRefNode;
						nElementCountR = nRElements;
						pJacobian = &dJacobianNode;
						pDerivR = &dDerivRNode;
						pContraMetricA = &dContraMetricA;
						pContraMetricB = &dContraMetricB;
						pContraMetricXi = &dContraMetricXi;

					} else if (pGrid->GetVarLocation(c) == DataLocation_REdge) {
						pDataInitial = &dataInitialREdge;
						pDataUpdate  = &dataUpdateREdge;
						pDataResidual = &dataResidualREdge;
						pDataRef = &dataRefREdge;
						nElementCountR = nRElements + 1;
						pJacobian = &dJacobianREdge;
						pDerivR = &dDerivRREdge;
						pContraMetricA = &dContraMetricAREdge;
						pContraMetricB = &dContraMetricBREdge;
						pContraMetricXi = &dContraMetricXiREdge;

					} else {
						_EXCEPTIONT("UNIMPLEMENTED");
					}

				} else {
					pDataInitial = &dataInitialTracer;
					pDataUpdate = &dataUpdateTracer;
					pDataResidual = &dataResidualTracer;
					nElementCountR = nRElements;
					pJacobian = &dJacobianNode;
				}

	// Loop over all finite elements
	for (int a = 0; a < nElementCountA; a++) {
	for (int b = 0; b < nElementCountB; b++) {
	for (int k = 0; k < nElementCountR; k++) {

		int iElementA =
			a * nHorizontalOrder + box.GetHaloElements();
		int iElementB =
			b * nHorizontalOrder + box.GetHaloElements();

		// Store the buffer state
		for (int i = 0; i < nHorizontalOrder; i++) {
		for (int j = 0; j < nHorizontalOrder; j++) {
			int iA = iElementA + i;
			int iB = iElementB + j;

			m_dBufferState(i,j) =
				(*pDataInitial)(c,iA,iB,k);
		}
		}

		// Remove the reference state from the buffer state
		if (fRemoveRefState) {
			for (int i = 0; i < nHorizontalOrder; i++) {
			for (int j = 0; j < nHorizontalOrder; j++) {
				int iA = iElementA + i;
				int iB = iElementB + j;

				m_dBufferState(i,j) -=
					(*pDataRef)(c,iA,iB,k);
			}
			}
		}

		// Compute maximum local wave speed in an element
		for (int i = 0; i < nHorizontalOrder; i++) {
		for (int j = 0; j < nHorizontalOrder; j++) {
			int iA = iElementA + i;
			int iB = iElementB + j;

			// Contravariant velocities
			double dCovUa = (*pDataInitial)(UIx,iA,iB,k);
			double dCovUb = (*pDataInitial)(VIx,iA,iB,k);

			// Contravariant velocities
			m_dAuxDataNode(ConUaIx,i,j,k) =
				dContraMetric2DA(iA,iB,0) * dCovUa
			      + dContraMetric2DA(iA,iB,1) * dCovUb;

			m_dAuxDataNode(ConUbIx,i,j,k) =
				dContraMetric2DB(iA,iB,0) * dCovUa
			      + dContraMetric2DB(iA,iB,1) * dCovUb;

			double dRhoTheta = 0.0;
#if defined(FORMULATION_THETA)
			dRhoTheta = (*pDataInitial)(PIx,iA,iB,k) *
				  (*pDataInitial)(RIx,iA,iB,k);
#endif

#if defined(FORMULATION_RHOTHETA_P) || defined(FORMULATION_RHOTHETA_PI)
			dRhoTheta = (*pDataInitial)(PIx,iA,iB,k);
#endif

			// Maximum wave speed in this element
			m_dAuxDataNode(KIx,i,j,k) = sqrt(
			      m_dAuxDataNode(ConUaIx,i,j,k) * dCovUa
			    + m_dAuxDataNode(ConUbIx,i,j,k) * dCovUb)
			    + sqrt(phys.GetGamma() *
				       phys.PressureFromRhoTheta(dRhoTheta) /
				       (*pDataInitial)(RIx,iA,iB,k));
		}
		}

		double dAvgA = 0.0;
		double dAvgB = 0.0;
		double dEAvgU = 0.0;
		double dEAvgV = 0.0;
		double dEAvgW = 0.0;
		double dEAvgP = 0.0;
		double dEAvgR = 0.0;
		// Calculate the element average of the scalar field
		for (int i = 0; i < nHorizontalOrder; i++) {
		for (int j = 0; j < nHorizontalOrder; j++) {
			int iA = iElementA + i;
			int iB = iElementB + j;

			for (int o = 0; o < m_model.GetEquationSet().GetComponents(); o++) {
				for (int s = 0; s < nHorizontalOrder; s++) {
					dAvgA +=
					  (*pDataInitial)(o,s,iB,k)
					  * dStiffness1D(i,s);

					dAvgB +=
					  (*pDataInitial)(o,iA,s,k)
					  * dStiffness1D(j,s);
				}

				dAvgA *= dInvElementDeltaA;
				dAvgB *= dInvElementDeltaB;

				switch (o) {
				  case 0: dEAvgU = 0.5 * (dAvgA + dAvgB); break;
				  case 1: dEAvgV = 0.5 * (dAvgA + dAvgB); break;
				  case 2: dEAvgW = 0.5 * (dAvgA + dAvgB); break;
				  case 3: dEAvgP = 0.5 * (dAvgA + dAvgB); break;
				  case 4: dEAvgR = 0.5 * (dAvgA + dAvgB); break;
				}
			}
		}
		}

		double dResU = 0.0;
		double dResV = 0.0;
		double dResW = 0.0;
		double dResP = 0.0;
		double dResR = 0.0;
		double dResMax = 0.0;
		double dResCoeff = 0.0;
		double dNuMax = 0.0;
		double dGridLength = dElementLength / (nHorizontalOrder - 1);

		// Compute the local diffusion coefficient
		for (int i = 0; i < nHorizontalOrder; i++) {
		for (int j = 0; j < nHorizontalOrder; j++) {
			int iA = iElementA + i;
			int iB = iElementB + j;

			// Compute the local diffusion coefficient
			//
			dResU = 0.0;
			//dResU = fabs((*pDataResidual)(UIx,iA,iB,k));//  / fabs(
	                             //(*pDataInitial)(UIx,iA,iB,k) - dEAvgU);
			dResV = 0.0;
			//dResV = fabs((*pDataResidual)(VIx,iA,iB,k));//  / fabs(
	                             //(*pDataInitial)(VIx,iA,iB,k) - dEAvgV);
			dResW = 0.0;
			//dResW = fabs((*pDataResidual)(WIx,iA,iB,k));//  / fabs(
	                             //(*pDataInitial)(WIx,iA,iB,k) - dEAvgW);
			dResP = 1.0;
			//dResP = fabs((*pDataResidual)(PIx,iA,iB,k));// / fabs(
				     //(*pDataInitial)(PIx,iA,iB,k) - dEAvgP);
			dResR = 0.0;
			//dResR = fabs((*pDataResidual)(RIx,iA,iB,k));// / fabs(
				     //(*pDataInitial)(RIx,iA,iB,k) - dEAvgR);
			/*
			dResU = fabs((*pDataResidual)(UIx,iA,iB,k))  / fabs(50.0);
			dResV = fabs((*pDataResidual)(VIx,iA,iB,k))  / fabs(50.0);
			dResW = fabs((*pDataResidual)(WIx,iA,iB,k))  / fabs(4.0);
			dResP = fabs((*pDataResidual)(PIx,iA,iB,k)) / fabs(0.1);
			dResR = fabs((*pDataResidual)(RIx,iA,iB,k)) / fabs(0.1);
			*/
			// Select the maximum residual
			dResMax = std::max(dResU, dResV);
			dResMax = std::max(dResMax, dResW);
			dResMax = std::max(dResMax, dResP);
			dResMax = std::max(dResMax, dResR);

			if (dResMax > 0.0) {
			    dResCoeff = dResMax;
			} else {
			    dResCoeff = 0.0;
			}

			// Scale to the average element length
			double dABLength = dGridLength
				/ pGrid->GetReferenceLength();
			dResNu[i][j] = dABLength * dABLength * fabs(dResCoeff);

			// Apply Pr number of 0.7 to thermodynamic stress
			if (c == PIx) {
				dResNu[i][j] *= 1.75;
			}

			// Get the maximum possible coefficient (upwind)
			//dNuMax = 0.5 * m_dAuxDataNode(KIx,i,j,k); // /
				//pGrid->GetReferenceLength();

			// Limit the coefficients to the upwind value
			//if (dResNu[i][j] >= dNuMax) {
			//	dResNu[i][j] = dNuMax;
			//}
			// Check for Inf or NaN and adjust
			if (!std::isfinite(dResNu[i][j])) {
				dResNu[i][j] = dNuMax;
			}
		}
		}

		// Calculate the pointwise gradient of the scalar field
		for (int i = 0; i < nHorizontalOrder; i++) {
		for (int j = 0; j < nHorizontalOrder; j++) {
			int iA = iElementA + i;
			int iB = iElementB + j;

			double dDaPsi = 0.0;
			double dDbPsi = 0.0;
			for (int s = 0; s < nHorizontalOrder; s++) {
				dDaPsi +=
					//(*pDataInitial)(c,s,iB,k)
					//* dDxBasis1D(s,i);
					m_dBufferState(s,j)
					* dDxBasis1D(s,i);

				dDbPsi +=
					//(*pDataInitial)(c,iA,s,k)
					//* dDxBasis1D(s,j);
					m_dBufferState(i,s)
					* dDxBasis1D(s,j);
			}

				dDaPsi *= dInvElementDeltaA;
				dDbPsi *= dInvElementDeltaB;

				dDaPsi *= dResNu[i][j];
				dDbPsi *= dResNu[i][j];

				m_dJGradientA(i,j) =
					(*pJacobian)(iA,iB,k) * (
						+ dContraMetric2DA(iA,iB,0) * dDaPsi
						+ dContraMetric2DA(iA,iB,1) * dDbPsi);

				m_dJGradientB(i,j) =
					(*pJacobian)(iA,iB,k) * (
						+ dContraMetric2DB(iA,iB,0) * dDaPsi
						+ dContraMetric2DB(iA,iB,1) * dDbPsi);

				// Update the DynSGS stress field for output
				dataDynSGSNode(0,iA,iB,k) = m_dJGradientA(i,j);
				dataDynSGSNode(1,iA,iB,k) = m_dJGradientB(i,j);
		}
		}

		// Pointwise updates
		for (int i = 0; i < nHorizontalOrder; i++) {
		for (int j = 0; j < nHorizontalOrder; j++) {
			int iA = iElementA + i;
			int iB = iElementB + j;

			// Inverse Jacobian and Jacobian
			const double dInvJacobian = 1.0 / (*pJacobian)(iA,iB,k);
			const double dInvRho = 1.0 / (*pDataInitial)(RIx,iA,iB,k);

			// Compute integral term
			double dUpdateA = 0.0;
			double dUpdateB = 0.0;

			for (int s = 0; s < nHorizontalOrder; s++) {
				dUpdateA +=
					m_dJGradientA(s,j)
					* dStiffness1D(i,s);

				dUpdateB +=
					m_dJGradientB(i,s)
					* dStiffness1D(j,s);
			}

			dUpdateA *= dInvElementDeltaA;
			dUpdateB *= dInvElementDeltaB;

			// Apply update
			(*pDataUpdate)(c,iA,iB,k) +=
				dDeltaT * dInvRho
					* dInvJacobian
					* (dUpdateA + dUpdateB);

			//Announce("%1.10e %1.10e %1.10e %1.10e", dResNu[i][j], dUpdateA, dUpdateB, (*pDataUpdate)(c,iA,iB,k));
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

void HorizontalDynamicsFEM::ApplyVectorHyperdiffusion(
	int iDataInitial,
	int iDataUpdate,
	double dDeltaT,
	double dNuDiv,
	double dNuVort,
	bool fScaleNuLocally
) {
	// Indices of auxiliary data variables
	const int ConUaIx = 0;
	const int ConUbIx = 1;

	// Variable indices
	const int UIx = 0;
	const int VIx = 1;
	const int WIx = 3;

	// Density or height index
	int RIx = 4;

	const EquationSet & eqn = m_model.GetEquationSet();
	if (eqn.GetType() == EquationSet::ShallowWaterEquations) {
		RIx = 2;
	}

	// Get a copy of the GLL grid
	GridGLL * pGrid = dynamic_cast<GridGLL*>(m_model.GetGrid());

#if defined(FIXED_HORIZONTAL_ORDER)
	const int nHorizontalOrder = FIXED_HORIZONTAL_ORDER;
	if (nHorizontalOrder != m_nHorizontalOrder) {
		_EXCEPTIONT("Command line order must match FIXED_HORIZONTAL_ORDER");
	}
#else
	const int nHorizontalOrder = m_nHorizontalOrder;
#endif

#if defined(FIXED_RELEMENTS)
	const int nRElements = FIXED_RELEMENTS;
	if (nRElements != pGrid->GetRElements()) {
		_EXCEPTIONT("Command line levels must match FIXED_RELEMENTS");
	}
#else
	const int nRElements = pGrid->GetRElements();
#endif

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

		DataArray3D<double> dataRho;
		dataRho.SetSize(
			dataInitial.GetSize(1),
			dataInitial.GetSize(2),
			dataInitial.GetSize(3));

		if (fApplyToRefState) {
			dataUa.AttachToData(&(dataRef(UIx,0,0,0)));
			dataUb.AttachToData(&(dataRef(VIx,0,0,0)));
			dataRho.AttachToData(&(dataRef(RIx,0,0,0)));
		} else {
			dataUa.AttachToData(&(dataInitial(UIx,0,0,0)));
			dataUb.AttachToData(&(dataInitial(VIx,0,0,0)));
			dataRho.AttachToData(&(dataInitial(RIx,0,0,0)));
		}

		// Compute curl and divergence of U on the grid
		pPatch->ComputeCurlAndDiv(dataUa, dataUb, dataRho);

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
		for (int a = 0; a < nElementCountA; a++) {
		for (int b = 0; b < nElementCountB; b++) {
		for (int k = 0; k < nRElements; k++) {

			const int iElementA = a * nHorizontalOrder + box.GetHaloElements();
			const int iElementB = b * nHorizontalOrder + box.GetHaloElements();

			// Pointwise update of horizontal velocities
			for (int i = 0; i < nHorizontalOrder; i++) {
			for (int j = 0; j < nHorizontalOrder; j++) {

				const int iA = iElementA + i;
				const int iB = iElementB + j;

				// Compute hyperviscosity sums
				double dDaDiv = 0.0;
				double dDbDiv = 0.0;

				double dDaCurl = 0.0;
				double dDbCurl = 0.0;

#pragma unroll
				for (int s = 0; s < nHorizontalOrder; s++) {
					dDaDiv -=
						  dStiffness1D(i,s)
						* dataDiv(iElementA+s,iB,k);

					dDbDiv -=
						  dStiffness1D(j,s)
						* dataDiv(iA,iElementB+s,k);

					dDaCurl -=
						  dStiffness1D(i,s)
						* dataCurl(iElementA+s,iB,k);

					dDbCurl -=
						  dStiffness1D(j,s)
						* dataCurl(iA,iElementB+s,k);
				}

				dDaDiv *= dInvElementDeltaA;
				dDbDiv *= dInvElementDeltaB;

				dDaCurl *= dInvElementDeltaA;
				dDbCurl *= dInvElementDeltaB;

				// Apply update
				double dUpdateUa =
					+ dLocalNuDiv * dDaDiv
					- dLocalNuVort * dJacobian2D(iA,iB) * (
						  dContraMetric2DB(iA,iB,0) * dDaCurl
						+ dContraMetric2DB(iA,iB,1) * dDbCurl);

				double dUpdateUb =
					+ dLocalNuDiv * dDbDiv
					+ dLocalNuVort * dJacobian2D(iA,iB) * (
						  dContraMetric2DA(iA,iB,0) * dDaCurl
						+ dContraMetric2DA(iA,iB,1) * dDbCurl);

				dataUpdate(UIx,iA,iB,k) -= dDeltaT * dUpdateUa;

				if (pGrid->GetIsCartesianXZ() == false) {
					dataUpdate(VIx,iA,iB,k) -= dDeltaT * dUpdateUb;
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
	int iDataInitial,
	int iDataUpdate,
	double dDeltaT
) {
	// Get a copy of the GLL grid
	GridGLL * pGrid = dynamic_cast<GridGLL*>(m_model.GetGrid());

#if defined(FIXED_HORIZONTAL_ORDER)
	const int nHorizontalOrder = FIXED_HORIZONTAL_ORDER;
	if (nHorizontalOrder != m_nHorizontalOrder) {
		_EXCEPTIONT("Command line order must match FIXED_HORIZONTAL_ORDER");
	}
#else
	const int nHorizontalOrder = m_nHorizontalOrder;
#endif

#if defined(FIXED_RELEMENTS)
	const int nRElements = FIXED_RELEMENTS;
	if (nRElements != pGrid->GetRElements()) {
		_EXCEPTIONT("Command line levels must match FIXED_RELEMENTS");
	}
#else
	const int nRElements = pGrid->GetRElements();
#endif

	// Number of components to hit with friction
	int nComponents = m_model.GetEquationSet().GetComponents();

	// Equation set being solved
	const int nEqSet = m_model.GetEquationSet().GetType();

	const bool fCartXZ = pGrid->GetIsCartesianXZ();

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
	int nRayCycles = 2;
	double dRayFactor = 1.0 / nRayCycles;

	// Perform local update
	for (int n = 0; n < pGrid->GetActivePatchCount(); n++) {
		GridPatchGLL * pPatch =
			dynamic_cast<GridPatchGLL*>(pGrid->GetActivePatch(n));

		const PatchBox & box = pPatch->GetPatchBox();

		// Grid data
		const DataArray4D<double> & dataInitialNode =
			pPatch->GetDataState(iDataInitial, DataLocation_Node);

		const DataArray4D<double> & dataInitialREdge =
			pPatch->GetDataState(iDataInitial, DataLocation_REdge);

		DataArray4D<double> & dataUpdateNode =
			pPatch->GetDataState(iDataUpdate, DataLocation_Node);

		DataArray4D<double> & dataUpdateREdge =
			pPatch->GetDataState(iDataUpdate, DataLocation_REdge);

		// Reference state
		const DataArray4D<double> & dataReferenceNode =
			pPatch->GetReferenceState(DataLocation_Node);

		const DataArray4D<double> & dataReferenceREdge =
			pPatch->GetReferenceState(DataLocation_REdge);

		// Top PML layer strength
		const DataArray3D<double> & dataRayStrengthNode =
			pPatch->GetRayleighStrength(DataLocation_Node);

		const DataArray3D<double> & dataRayStrengthREdge =
			pPatch->GetRayleighStrength(DataLocation_REdge);

		const DataArray4D<double> & m_dataXiDiffNode =
			pPatch->GetXiDiffNode(DataLocation_Node);

		// Loop over all nodes in patch
		for (int i = box.GetAInteriorBegin(); i < box.GetAInteriorEnd(); i++) {
		for (int j = box.GetBInteriorBegin(); j < box.GetBInteriorEnd(); j++) {

			// Rayleigh damping on nodes
			for (int k = 0; k < nRElements; k++) {
				double dNu = dataRayStrengthNode(i,j,k);

				// Backwards Euler
				if (dNu == 0.0) {
					continue;
				}

				double dNuInv = 0.0;

				// Loop over all effective components
				for (int c = 0; c < nComponents; c++) {
					if (pGrid->GetVarLocation(nEffectiveC[c]) ==
						DataLocation_Node) {
						// Get the PML part of the update
						double dPhi = (1.0 / dDeltaT) *
							(dataUpdateNode(nEffectiveC[c],i,j,k)
							 - dataInitialNode(nEffectiveC[c],i,j,k));

						if (nEffectiveC[0] == 0) {
							dPhi -= m_dataXiDiffNode(0, i, j, k);
						}

						if (nEffectiveC[1] == 1) {
							dPhi -= m_dataXiDiffNode(1, i, j, k);
						}

						double dPML = 0.0;
						for (int si = 0; si < nRayCycles; si++) {
							dNuInv = 1.0 / (1.0 + dRayFactor * dDeltaT * dNu);
							dPML = dataUpdateNode(nEffectiveC[c],i,j,k)
							- dataReferenceNode(nEffectiveC[c],i,j,k)
							+ dRayFactor * dDeltaT * dNu * dPhi;

							dataUpdateNode(nEffectiveC[c],i,j,k) =
							dNuInv * dataUpdateNode(nEffectiveC[c],i,j,k)
							+ (1.0 - dNuInv)
							* dataReferenceNode(nEffectiveC[c],i,j,k);
							//- dNuInv * dRayFactor * dDeltaT * dPML;
						}
					}
				}
			}

			// Rayleigh damping on interfaces
			for (int k = 0; k <= nRElements; k++) {
				double dNu = dataRayStrengthREdge(i,j,k);

				// Backwards Euler
				if (dNu == 0.0) {
					continue;
				}

				double dNuInv = 0.0;

				// Loop over all effective components
				for (int c = 0; c < nComponents; c++) {
					if (pGrid->GetVarLocation(nEffectiveC[c]) ==
						DataLocation_REdge) {
						// Get the PML part of the update
						double dPhi = (1.0 / dDeltaT) *
							(dataUpdateREdge(nEffectiveC[c],i,j,k)
							 - dataInitialREdge(nEffectiveC[c],i,j,k));

						if (nEffectiveC[0] == 0) {
 							dPhi -= m_dataXiDiffNode(0, i, j, k);
 						}

 						if (nEffectiveC[1] == 1) {
 							dPhi -= m_dataXiDiffNode(1, i, j, k);
 						}

						double dPML = 0.0;
						for (int si = 0; si < nRayCycles; si++) {
							dNuInv = 1.0 / (1.0 + dRayFactor * dDeltaT * dNu);
							dPML = dataUpdateREdge(nEffectiveC[c],i,j,k)
							- dataReferenceREdge(nEffectiveC[c],i,j,k)
							+ dRayFactor * dDeltaT * dNu * dPhi;

							dataUpdateREdge(nEffectiveC[c],i,j,k) =
							dNuInv * dataUpdateREdge(nEffectiveC[c],i,j,k)
							+ (1.0 - dNuInv)
							* dataReferenceREdge(nEffectiveC[c],i,j,k);
							//- dNuInv * dRayFactor * dDeltaT * dPML;
						}
					}
				}
			}
		}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

int HorizontalDynamicsFEM::GetSubStepAfterSubCycleCount() {
	int iSubStepCount = m_nHyperviscosityOrder / 2;

	return iSubStepCount;
}

///////////////////////////////////////////////////////////////////////////////

int HorizontalDynamicsFEM::SubStepAfterSubCycle(
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
			ApplyRayleighFriction(iDataInitial, iDataUpdate, dDeltaT);
		}
#endif
		return iDataUpdate;
	}

	_EXCEPTIONT("Invalid iSubStep");
}

///////////////////////////////////////////////////////////////////////////////

void HorizontalDynamicsFEM::StepAfterSubCycle(
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
			ApplyRayleighFriction(iDataInitial, iDataUpdate, dDeltaT);

			// Apply Direct Stiffness Summation
			//pGrid->ApplyDSS(iDataUpdate, DataType_State);
			//pGrid->ApplyDSS(iDataUpdate, DataType_Tracers);
		}
#endif
}

///////////////////////////////////////////////////////////////////////////////
