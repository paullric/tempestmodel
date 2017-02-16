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
	// Only horizontal order 4 supported here
	if (nHorizontalOrder != 4) {
		_EXCEPTIONT("Only horizontal order 4 supported");
	}
}

///////////////////////////////////////////////////////////////////////////////

void HorizontalDynamicsFEMV2::Initialize() {

#if defined(PROGNOSTIC_CONTRAVARIANT_MOMENTA)
	_EXCEPTIONT("Prognostic contrvariant velocities not supported");
#endif

	// Get a copy of the GLL grid
	GridGLL * pGrid = dynamic_cast<GridGLL*>(m_model.GetGrid());
	if (pGrid == NULL) {
		_EXCEPTIONT("Grid must be of type GridGLL");
	}

	// Only horizontal order 4 supported
	const int nHorizontalOrder = 4;

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

void HorizontalDynamicsFEMV2::FilterNegativeTracers(
	int iDataUpdate
) {
#ifdef POSITIVE_DEFINITE_FILTER_TRACERS
	// Get a copy of the GLL grid
	GridGLL * pGrid = dynamic_cast<GridGLL*>(m_model.GetGrid());

	// Horizontal order 4
	const int nHorizontalOrder = 4;

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

void HorizontalDynamicsFEMV2::StepShallowWater(
	int iDataInitial,
	int iDataUpdate,
	const Time & time,
	double dDeltaT
) {

	// Get a copy of the GLL grid
	GridGLL * pGrid = dynamic_cast<GridGLL*>(m_model.GetGrid());

	// Horizontal order 4
	const int nHorizontalOrder = 4;

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

///////////////////////////////////////////////////////////////////////////////

void HorizontalDynamicsFEMV2::StepNonhydrostaticPrimitive(
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

	// Horizontal order 4
	const int nHorizontalOrder = 4;

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

				// NOTE: For some reason using parenthetical notation on the
				//       LHS of these expressions crash the Intel compiler
				//       (observed with icpc 16.0.3)

				// Exner pressure
				m_dAuxDataNode[ExnerIx][i][j][k] =
					phys.ExnerPressureFromRhoTheta(
						dataInitialNode(PIx,iA,iB,k));

			}
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
/*
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
*/
				// Contravariant velocities
				const double dConUa =
					  dContraMetricA(iA,iB,k,0) * dCovUa
					+ dContraMetricA(iA,iB,k,1) * dCovUb
					+ dContraMetricA(iA,iB,k,2) * dCovUx;

				const double dConUb =
					  dContraMetricB(iA,iB,k,0) * dCovUa
					+ dContraMetricB(iA,iB,k,1) * dCovUb
					+ dContraMetricB(iA,iB,k,2) * dCovUx;

				const double dConUx =
					  dContraMetricXi(iA,iB,k,0) * dCovUa
					+ dContraMetricXi(iA,iB,k,1) * dCovUb
					+ dContraMetricXi(iA,iB,k,2) * dCovUx;

				m_dAuxDataNode(ConUaIx,i,j,k) = dConUa;
				m_dAuxDataNode(ConUbIx,i,j,k) = dConUb;
				m_dAuxDataNode(ConUxIx,i,j,k) = dConUx;

				// Specific kinetic energy
				m_dAuxDataNode(KIx,i,j,k) = 0.5 * (
					  dConUa * dCovUa
					+ dConUb * dCovUb
					+ dConUx * dCovUx);
			}
			}
			}

			// Compute curl term
			for (int i = 0; i < nHorizontalOrder; i++) {
			for (int j = 0; j < nHorizontalOrder; j++) {
			for (int k = 0; k < nRElements; k++) {

				const int iA = iElementA + i;
				const int iB = iElementB + j;

				// Vertical derivatives of the covariant velocity field
				const double dConUa = m_dAuxDataNode(ConUaIx,i,j,k);
				const double dConUb = m_dAuxDataNode(ConUbIx,i,j,k);
				const double dConUx = m_dAuxDataNode(ConUxIx,i,j,k);

				const double dCovDxUa =
					pGrid->DifferentiateNodeToNode(
						&(dataInitialNode(UIx,iA,iB,0)),
						k, 1);

				const double dCovDxUb =
					pGrid->DifferentiateNodeToNode(
						&(dataInitialNode(VIx,iA,iB,0)),
						k, 1);

				const double dCovDaUb = dInvElementDeltaA * (
					  dataInitialNode(VIx,iElementA+0,iB,k) * dDxBasis1D(0,i)
					+ dataInitialNode(VIx,iElementA+1,iB,k) * dDxBasis1D(1,i)
					+ dataInitialNode(VIx,iElementA+2,iB,k) * dDxBasis1D(2,i)
					+ dataInitialNode(VIx,iElementA+3,iB,k) * dDxBasis1D(3,i));

				const double dCovDaUx = dInvElementDeltaA * (
					  dataInitialNode(WIx,iElementA+0,iB,k) * dDxBasis1D(0,i)
					+ dataInitialNode(WIx,iElementA+1,iB,k) * dDxBasis1D(1,i)
					+ dataInitialNode(WIx,iElementA+2,iB,k) * dDxBasis1D(2,i)
					+ dataInitialNode(WIx,iElementA+3,iB,k) * dDxBasis1D(3,i));

				const double dCovDbUa = dInvElementDeltaB * (
					  dataInitialNode(UIx,iA,iElementB+0,k) * dDxBasis1D(0,j)
					+ dataInitialNode(UIx,iA,iElementB+1,k) * dDxBasis1D(1,j)
					+ dataInitialNode(UIx,iA,iElementB+2,k) * dDxBasis1D(2,j)
					+ dataInitialNode(UIx,iA,iElementB+3,k) * dDxBasis1D(3,j));

				const double dCovDbUx = dInvElementDeltaB * (
					  dataInitialNode(WIx,iA,iElementB+0,k) * dDxBasis1D(0,j)
					+ dataInitialNode(WIx,iA,iElementB+1,k) * dDxBasis1D(1,j)
					+ dataInitialNode(WIx,iA,iElementB+2,k) * dDxBasis1D(2,j)
					+ dataInitialNode(WIx,iA,iElementB+3,k) * dDxBasis1D(3,j));

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

			}
			}
			}

			// Compute fluxes within spectral element
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

				// RhoTheta flux
				m_dAlphaPressureFlux(i,j,k) =
					  dAlphaBaseFlux
					* dataInitialNode(PIx,iA,iB,k);

				m_dBetaPressureFlux(i,j,k) =
					  dBetaBaseFlux
					* dataInitialNode(PIx,iA,iB,k);
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

				// Derivatives
				const double dDaRhoFluxA = dInvElementDeltaA * (
					- m_dAlphaMassFlux(0,j,k) * dStiffness1D(i,0)
					- m_dAlphaMassFlux(1,j,k) * dStiffness1D(i,1)
					- m_dAlphaMassFlux(2,j,k) * dStiffness1D(i,2)
					- m_dAlphaMassFlux(3,j,k) * dStiffness1D(i,3));

				const double dDaPressureFluxA = dInvElementDeltaA * (
					- m_dAlphaPressureFlux(0,j,k) * dStiffness1D(i,0)
					- m_dAlphaPressureFlux(1,j,k) * dStiffness1D(i,1)
					- m_dAlphaPressureFlux(2,j,k) * dStiffness1D(i,2)
					- m_dAlphaPressureFlux(3,j,k) * dStiffness1D(i,3));

				const double dDbRhoFluxB = dInvElementDeltaB * (
					- m_dBetaMassFlux(i,0,k) * dStiffness1D(j,0)
					- m_dBetaMassFlux(i,1,k) * dStiffness1D(j,1)
					- m_dBetaMassFlux(i,2,k) * dStiffness1D(j,2)
					- m_dBetaMassFlux(i,3,k) * dStiffness1D(j,3));

				const double dDbPressureFluxB = dInvElementDeltaB * (
					- m_dBetaPressureFlux(i,0,k) * dStiffness1D(j,0)
					- m_dBetaPressureFlux(i,1,k) * dStiffness1D(j,1)
					- m_dBetaPressureFlux(i,2,k) * dStiffness1D(j,2)
					- m_dBetaPressureFlux(i,3,k) * dStiffness1D(j,3));

				const double dDaP = dInvElementDeltaA * (
					  m_dAuxDataNode(ExnerIx,0,j,k) * dDxBasis1D(0,i)
					+ m_dAuxDataNode(ExnerIx,1,j,k) * dDxBasis1D(1,i)
					+ m_dAuxDataNode(ExnerIx,2,j,k) * dDxBasis1D(2,i)
					+ m_dAuxDataNode(ExnerIx,3,j,k) * dDxBasis1D(3,i));

				const double dDbP = dInvElementDeltaB * (
					  m_dAuxDataNode(ExnerIx,i,0,k) * dDxBasis1D(0,j)
					+ m_dAuxDataNode(ExnerIx,i,1,k) * dDxBasis1D(1,j)
					+ m_dAuxDataNode(ExnerIx,i,2,k) * dDxBasis1D(2,j)
					+ m_dAuxDataNode(ExnerIx,i,3,k) * dDxBasis1D(3,j));

				const double dDaKE = dInvElementDeltaA * (
					  m_dAuxDataNode(KIx,0,j,k) * dDxBasis1D(0,i)
					+ m_dAuxDataNode(KIx,1,j,k) * dDxBasis1D(1,i)
					+ m_dAuxDataNode(KIx,2,j,k) * dDxBasis1D(2,i)
					+ m_dAuxDataNode(KIx,3,j,k) * dDxBasis1D(3,i));

				const double dDbKE = dInvElementDeltaB * (
					  m_dAuxDataNode(KIx,i,0,k) * dDxBasis1D(0,j)
					+ m_dAuxDataNode(KIx,i,1,k) * dDxBasis1D(1,j)
					+ m_dAuxDataNode(KIx,i,2,k) * dDxBasis1D(2,j)
					+ m_dAuxDataNode(KIx,i,3,k) * dDxBasis1D(3,j));

				// Pressure gradient force
				const double dInvRho = 1.0 / dataInitialNode(RIx,iA,iB,k);

				const double dPressureGradientForceUa =
					dInvRho * dataInitialNode(PIx,iA,iB,k) * dDaP;

				const double dPressureGradientForceUb =
					dInvRho * dataInitialNode(PIx,iA,iB,k) * dDbP;

				// Gravity
				const double dDaPhi = phys.GetG() * m_dLocalDerivR(i,j,k,0);
				const double dDbPhi = phys.GetG() * m_dLocalDerivR(i,j,k,1);

 				// Apply update to horizontal velocity
				dataUpdateNode(UIx,iA,iB,k) +=
					dDeltaT * (
						- dPressureGradientForceUa
						- dDaKE
						- dDaPhi
						+ m_dLocalCoriolisF(i,j) * m_dLocalJacobian2D(i,j) * dConUb
						+ m_dAuxDataNode(UCrossZetaAIx,i,j,k));

				dataUpdateNode(VIx,iA,iB,k) +=
					dDeltaT * (
						- dPressureGradientForceUb
						- dDbKE
						- dDbPhi
						- m_dLocalCoriolisF(i,j) * m_dLocalJacobian2D(i,j) * dConUa
						+ m_dAuxDataNode(UCrossZetaBIx,i,j,k));

				// Apply update to density
				dataUpdateNode(RIx,iA,iB,k) -=
					dDeltaT * dInvJacobian * (
						  dDaRhoFluxA
						+ dDbRhoFluxB);

				// Apply update to rhotheta
				dataUpdateNode(PIx,iA,iB,k) -=
					dDeltaT * dInvJacobian * (
						  dDaPressureFluxA
						+ dDbPressureFluxB);
/*
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
*/
			}
			}
			}

			// Update vertical velocity on interfaces
			{

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

	// Horizontal order 4
	const int nHorizontalOrder = 4;

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

	// Density or height index
	int RIx = 4;

	const EquationSet & eqn = m_model.GetEquationSet();
	if (eqn.GetType() == EquationSet::ShallowWaterEquations) {
		RIx = 2;
	}

	// Get a copy of the GLL grid
	GridGLL * pGrid = dynamic_cast<GridGLL*>(m_model.GetGrid());

	// Horizontal order 4
	const int nHorizontalOrder = 4;

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

//#pragma message "jeguerra: Clean up this function"

void HorizontalDynamicsFEMV2::ApplyRayleighFriction(
	int iDataUpdate,
	double dDeltaT
) {
	// Get a copy of the GLL grid
	GridGLL * pGrid = dynamic_cast<GridGLL*>(m_model.GetGrid());

	// Horizontal order 4
	const int nHorizontalOrder = 4;

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
			for (int k = 0; k < nRElements; k++) {

				double dNu = dataRayleighStrengthNode(i,j,k);

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
							dataUpdateNode(nEffectiveC[c],i,j,k) =
								dNuNode * dataUpdateNode(nEffectiveC[c],i,j,k)
								+ (1.0 - dNuNode)
								* dataReferenceNode(nEffectiveC[c],i,j,k);
						}
					}
				}
			}

			// Rayleigh damping on interfaces
			for (int k = 0; k <= nRElements; k++) {

				double dNu = dataRayleighStrengthREdge(i,j,k);

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
							dataUpdateREdge(nEffectiveC[c],i,j,k) =
							dNuREdge * dataUpdateREdge(nEffectiveC[c],i,j,k)
							+ (1.0 - dNuREdge)
							* dataReferenceREdge(nEffectiveC[c],i,j,k);
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
