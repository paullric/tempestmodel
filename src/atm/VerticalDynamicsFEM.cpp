///////////////////////////////////////////////////////////////////////////////
///
///	\file    VerticalDynamicsFEM.cpp
///	\author  Paul Ullrich
///	\version May 20, 2013
///
///	<remarks>
///		Copyright 2000-2013 Paul Ullrich
///
///		This file is distributed as part of the Tempest source code package.
///		Permission is granted to use, copy, modify and distribute this
///		source code and its documentation under the terms of the GNU General
///		Public License.  This software is provided "as is" without express
///		or implied warranty.
///	</remarks>

#include "Defines.h"
#include "VerticalDynamicsFEM.h"
#include "TimestepScheme.h"

#include "Announce.h"
#include "Model.h"
#include "Grid.h"
#include "GridCSGLL.h"
#include "EquationSet.h"
#include "TimeObj.h"
#include "PolynomialInterp.h"
#include "LinearAlgebra.h"

///////////////////////////////////////////////////////////////////////////////

//#define VERTICAL_HYPERVISCOSITY
//#define VERTICAL_UPWINDING

//#define DIFFUSE_HORIZONTAL_VELOCITIES
//#define DIFFUSE_THERMO
//#define DIFFUSE_VERTICAL_VELOCITY
//#define DIFFUSE_RHO

//#define DETECT_CFL_VIOLATION
//#define CAP_VERTICAL_VELOCITY

//#define EXPLICIT_THERMO
//#define EXPLICIT_VERTICAL_VELOCITY_ADVECTION

//#define VERTICAL_VELOCITY_ADVECTION_CLARK

///////////////////////////////////////////////////////////////////////////////

#define DEBUG

///////////////////////////////////////////////////////////////////////////////

VerticalDynamicsFEM::VerticalDynamicsFEM(
	Model & model,
	int nHorizontalOrder,
	int nVerticalOrder,
	int nHypervisOrder,
	bool fFullyExplicit,
	bool fUseReferenceState,
	bool fForceMassFluxOnLevels
) :
	VerticalDynamics(model),
	m_nHorizontalOrder(nHorizontalOrder),
	m_nVerticalOrder(nVerticalOrder),
	m_fFullyExplicit(fFullyExplicit),
	m_fUseReferenceState(fUseReferenceState),
	m_fForceMassFluxOnLevels(fForceMassFluxOnLevels),
	m_nHypervisOrder(nHypervisOrder),
	m_dHypervisCoeff(0.0)
{
	if (nHypervisOrder % 2 == 1) {
		_EXCEPTIONT("Vertical hyperdiffusion order must be even.");
	}

	if (nHypervisOrder < 0) {
		_EXCEPTIONT("Vertical hyperdiffusion order must be nonnegative.");
	}
}

///////////////////////////////////////////////////////////////////////////////

VerticalDynamicsFEM::~VerticalDynamicsFEM() {
#ifdef USE_JFNK_PETSC
	SNESDestroy(&m_snes);
	VecDestroy(&m_vecX);
	VecDestroy(&m_vecR);
	MatDestroy(&m_matJ);
#endif
}

///////////////////////////////////////////////////////////////////////////////

void VerticalDynamicsFEM::Initialize() {

	// Indices of EquationSet variables
	const int UIx = 0;
	const int VIx = 1;
	const int PIx = 2;
	const int WIx = 3;
	const int RIx = 4;

	// Pointer to grid
	GridGLL * pGrid = dynamic_cast<GridGLL *>(m_model.GetGrid());
	if (pGrid == NULL) {
		_EXCEPTIONT("Invalid grid -- expected GridGLL");
	}

	// Physical constants
	const PhysicalConstants & phys = m_model.GetPhysicalConstants();

	// Number of radial elements
	int nRElements = pGrid->GetRElements();

	// Number of radial elements in the column
	m_nRElements = nRElements;

	// Number of degrees of freedom per column in u/v/rho/w/theta
	m_nColumnStateSize = FTot * (nRElements + 1);

#ifdef USE_JFNK_PETSC
	// Initialize the PetSc solver context
	SNESCreate(PETSC_COMM_SELF, &m_snes);

	// Create vectors
	VecCreate(PETSC_COMM_SELF, &m_vecX);
	VecSetSizes(m_vecX, PETSC_DECIDE, m_nColumnStateSize);
	VecSetFromOptions(m_vecX);
	VecDuplicate(m_vecX, &m_vecR);

	// Set tolerances
	SNESSetTolerances(
		m_snes,
		1.0e-8,
		1.0e-8,
		1.0e-8,
		1,
		50);

	// Set the function
	SNESSetFunction(
		m_snes,
		m_vecR,
		VerticalDynamicsFEM_FormFunction,
		(void*)(this));

	MatCreateSNESMF(m_snes, &m_matJ);

	SNESSetJacobian(m_snes, m_matJ, m_matJ, MatMFFDComputeJacobian, NULL);

	// Set the SNES context from options
	SNESSetFromOptions(m_snes);
#endif
#ifdef USE_JFNK_GMRES
	// Initialize JFNK
	InitializeJFNK(m_nColumnStateSize, m_nColumnStateSize, 1.0e-5);
#endif
#ifdef USE_DIRECTSOLVE_APPROXJ
	// Initialize Jacobian matrix
	m_matJacobianF.Allocate(m_nColumnStateSize, m_nColumnStateSize);

	// Initialize pivot vector
	m_vecIPiv.Allocate(m_nColumnStateSize);
#endif
#ifdef USE_DIRECTSOLVE
	// Initialize Jacobian matrix
	m_matJacobianF.Allocate(m_nColumnStateSize, m_nColumnStateSize);

	// Initialize pivot vector
	m_vecIPiv.Allocate(m_nColumnStateSize);
#endif
#if defined(USE_DIRECTSOLVE_APPROXJ) || defined(USE_DIRECTSOLVE)
#ifdef USE_JACOBIAN_DIAGONAL
	if (m_nHypervisOrder > 2) {
		_EXCEPTIONT("Diagonal Jacobian only implemented for "
			"Hypervis order <= 2");
	}
	if (m_nVerticalOrder == 1) {
		m_nJacobianFKL = 4;
		m_nJacobianFKU = 4;
	} else if (m_nVerticalOrder == 2) {
		m_nJacobianFKL = 9;
		m_nJacobianFKU = 9;
	} else if (m_nVerticalOrder == 3) {
		m_nJacobianFKL = 15;
		m_nJacobianFKU = 15;
	} else if (m_nVerticalOrder == 4) {
		m_nJacobianFKL = 22;
		m_nJacobianFKU = 22;
	} else if (m_nVerticalOrder == 5) {
		m_nJacobianFKL = 30;
		m_nJacobianFKU = 30;
	} else {
		_EXCEPTIONT("UNIMPLEMENTED: At this vertical order");
	}
#endif
#endif

	// Announce vertical dynamics configuration
	AnnounceStartBlock("Configuring VerticalDynamicsFEM");

#if defined(VERTICAL_HYPERVISCOSITY)
	Announce("Hyperviscosity enabled");
#elif defined(VERTICAL_UPWINDING)
	Announce("Upwinding enabled");
#else
	Announce("Vertical diffusion disabled");
#endif

#if defined(VERTICAL_HYPERVISCOSITY) || defined(VERTICAL_UPWINDING)
#if defined(DIFFUSE_HORIZONTAL_VELOCITIES)
	Announce("Diffuse horizontal velocities");
#endif
#if defined(DIFFUSE_THERMO)
	Announce("Diffuse thermodynamic variable");
#endif
#if defined(DIFFUSE_VERTICAL_VELOCITY)
	Announce("Diffuse vertical velocity");
#endif
#if defined(DIFFUSE_RHO)
	Announce("Diffuse density");
#endif
#endif

#if defined(EXPLICIT_THERMO)
	Announce("Explicit thermodynamic advection");
#else
	Announce("Implicit thermodynamic advection");
#endif
#if defined(EXPLICIT_VERTICAL_VELOCITY_ADVECTION)
	Announce("Explicit vertical velocity advection");
#else
	Announce("Implicit vertical velocity advection");
#endif
#if defined(VERTICAL_VELOCITY_ADVECTION_CLARK)
	Announce("Vertical velocity advection in Clark form");
#else
	Announce("Vertical velocity advection in advective form");
#endif

	// End block
	AnnounceEndBlock(NULL);

	// Upwind weights
	m_dUpwindWeights.Allocate(nRElements / m_nVerticalOrder - 1);

	// Allocate column for JFNK
	m_dColumnState.Allocate(m_nColumnStateSize);

	// Allocation reference column
	m_dStateRefNode.Allocate(5, nRElements);
	m_dStateRefREdge.Allocate(5, nRElements+1);

	// Solution vector from JFNK
	m_dSoln.Allocate(m_nColumnStateSize);

	// State vector at levels
	m_dStateNode.Allocate(
		m_model.GetEquationSet().GetComponents(),
		nRElements);

	// State vector at interfaces
	m_dStateREdge.Allocate(
		m_model.GetEquationSet().GetComponents(),
		nRElements+1);
/*
	// Auxiliary variables at interfaces
	int nFiniteElements = nRElements / m_nVerticalOrder;
	if (nRElements % m_nVerticalOrder != 0) {
		_EXCEPTIONT("Logic error: Vertical order must divide RElements");
	}
*/
	// Auxiliary variables
	m_dStateAux.Allocate(nRElements+1);
	m_dStateAuxDiff.Allocate(nRElements+1);

	m_dXiDotNode.Allocate(nRElements);
	m_dXiDotREdge.Allocate(nRElements+1);

	m_dDiffUa.Allocate(nRElements+1);
	m_dDiffUb.Allocate(nRElements+1);

	m_dDiffPNode.Allocate(nRElements);
	m_dDiffPREdge.Allocate(nRElements+1);

#ifdef UNIFORM_DIFFUSION
	m_dDiffDiffState.Allocate(
		m_model.GetEquationSet().GetComponents(),
		nRElements+1);
#endif

	m_dHyperDiffState.Allocate(
		m_model.GetEquationSet().GetComponents(),
		nRElements+1);

	m_dDiffThetaNode.Allocate(nRElements);
	m_dDiffThetaREdge.Allocate(nRElements+1);

	m_dDiffWNode.Allocate(nRElements);
	m_dDiffWREdge.Allocate(nRElements+1);

	m_dHorizKineticEnergyNode.Allocate(nRElements);
	m_dKineticEnergyNode.Allocate(nRElements);
	m_dDiffKineticEnergyNode.Allocate(nRElements);
	m_dDiffKineticEnergyREdge.Allocate(nRElements+1);

	m_dMassFluxNode.Allocate(nRElements);
	m_dDiffMassFluxNode.Allocate(nRElements);
	m_dMassFluxREdge.Allocate(nRElements+1);
	m_dDiffMassFluxREdge.Allocate(nRElements+1);

	m_dPressureFluxNode.Allocate(nRElements);
	m_dDiffPressureFluxNode.Allocate(nRElements);
	m_dPressureFluxREdge.Allocate(nRElements+1);
	m_dDiffPressureFluxREdge.Allocate(nRElements+1);

	m_dExnerNode.Allocate(nRElements);
	m_dExnerREdge.Allocate(nRElements+1);

	m_dTracerDensityNode.Allocate(nRElements);
	m_dTracerDensityREdge.Allocate(nRElements+1);

	m_dInitialDensityNode.Allocate(nRElements);
	m_dInitialDensityREdge.Allocate(nRElements+1);

	m_dUpdateDensityNode.Allocate(nRElements);
	m_dUpdateDensityREdge.Allocate(nRElements+1);

	m_vecTracersF.Allocate(nRElements);
	m_matTracersLUDF.Allocate(nRElements, nRElements);
	m_vecTracersIPiv.Allocate(nRElements);
/*
	m_dExnerNode.Allocate(nRElements);
	m_dExnerRefNode.Allocate(nRElements);

	m_dDiffExnerNode.Allocate(nRElements);
	m_dDiffExnerRefNode.Allocate(nRElements);

	m_dExnerREdge.Allocate(nRElements+1);
	m_dExnerRefREdge.Allocate(nRElements+1);

	m_dDiffExnerREdge.Allocate(nRElements+1);
	m_dDiffExnerRefREdge.Allocate(nRElements+1);
*/
	// Variables to upwind
	m_fUpwind.Allocate(m_model.GetEquationSet().GetComponents());

#if defined(VERTICAL_HYPERVISCOSITY) || defined(VERTICAL_UPWINDING)
#ifdef DIFFUSE_HORIZONTAL_VELOCITIES
	m_fUpwind[UIx] = true;
	m_fUpwind[VIx] = true;
#endif
#ifdef DIFFUSE_THERMO
	m_fUpwind[PIx] = true;
#endif
#ifdef DIFFUSE_VERTICAL_VELOCITY
	m_fUpwind[WIx] = true;
#endif
#ifdef DIFFUSE_RHO
	m_fUpwind[RIx] = true;
#endif
#endif

	// Compute hyperviscosity coefficient
	if (m_nHypervisOrder == 0) {
		m_dHypervisCoeff = 0.0;

	} else if (m_nHypervisOrder == 2) {
		m_dHypervisCoeff = (1.0 / 2.0)
			* pow(1.0 / static_cast<double>(nRElements), 1.0);

	} else if (m_nHypervisOrder == 4) {
		m_dHypervisCoeff = - (1.0 / 6.0) //(1.0 / 12.0)
			* pow(1.0 / static_cast<double>(nRElements), 3.0);

	} else if (m_nHypervisOrder == 6) {
		m_dHypervisCoeff = (1.0 / 60.0)
			* pow(1.0 / static_cast<double>(nRElements), 5.0);

	} else if (m_nHypervisOrder == 8) {
		m_dHypervisCoeff = - (3.0 / 840.0)
			* pow(1.0 / static_cast<double>(nRElements), 7.0);

	} else {
		_EXCEPTIONT("UNIMPLEMENTED: Vertical hyperdiffusion order > 8");
	}

	// Metric quantities
	m_dColumnJacobianNode.Allocate(nRElements);
	m_dColumnJacobianREdge.Allocate(nRElements+1);
	m_dColumnElementArea.Allocate(nRElements);
	m_dColumnInvJacobianNode.Allocate(nRElements);
	m_dColumnInvJacobianREdge.Allocate(nRElements+1);
	m_dColumnDerivRNode.Allocate(nRElements, 3);
	m_dColumnDerivRREdge.Allocate(nRElements+1, 3);
	m_dColumnContraMetricA.Allocate(nRElements, 3);
	m_dColumnContraMetricB.Allocate(nRElements, 3);
	m_dColumnContraMetricXi.Allocate(nRElements, 3);
	m_dColumnContraMetricAREdge.Allocate(nRElements+1, 3);
	m_dColumnContraMetricBREdge.Allocate(nRElements+1, 3);
	m_dColumnContraMetricXiREdge.Allocate(nRElements+1, 3);
}

///////////////////////////////////////////////////////////////////////////////
void VerticalDynamicsFEM::ForceStepExplicit(
	int iDataInitial,
	int iDataUpdate,
	const Time & time,
	double dDeltaT
) {
        bool fOldFullyExplicit = m_fFullyExplicit;
        m_fFullyExplicit = true;
        StepExplicit(iDataInitial, iDataUpdate, time, dDeltaT);
        m_fFullyExplicit = fOldFullyExplicit;
}

///////////////////////////////////////////////////////////////////////////////
void VerticalDynamicsFEM::StepImplicitTermsExplicitly(
	int iDataInitial,
	int iDataUpdate,
	const Time & time,
	double dDeltaT
) {
	// Get a copy of the grid
	GridGLL * pGrid = dynamic_cast<GridGLL *>(m_model.GetGrid());

	// Physical constants
	const PhysicalConstants & phys = m_model.GetPhysicalConstants();

	// Indices of EquationSet variables
	const int UIx = 0;
	const int VIx = 1;
	const int PIx = 2;
	const int WIx = 3;
	const int RIx = 4;

	// Number of elements
	const int nRElements = pGrid->GetRElements();

	// Number of finite elements in the vertical
	const int nFiniteElements = nRElements / m_nVerticalOrder;

	// Store timestep size
	m_dDeltaT = dDeltaT;

	// Reset the reference state
	memset(m_dStateRefNode[WIx],  0,  nRElements   *sizeof(double));
	memset(m_dStateRefREdge[WIx], 0, (nRElements+1)*sizeof(double));

	// Perform local update
	for (int n = 0; n < pGrid->GetActivePatchCount(); n++) {
		GridPatch * pPatch = pGrid->GetActivePatch(n);

		const PatchBox & box = pPatch->GetPatchBox();

		// State Data
		const DataArray4D<double> & dataRefNode =
			pPatch->GetReferenceState(DataLocation_Node);

		const DataArray4D<double> & dataInitialNode =
			pPatch->GetDataState(iDataInitial, DataLocation_Node);

		DataArray4D<double> & dataUpdateNode =
			pPatch->GetDataState(iDataUpdate, DataLocation_Node);

		const DataArray4D<double> & dataRefREdge =
			pPatch->GetReferenceState(DataLocation_REdge);

		const DataArray4D<double> & dataInitialREdge =
			pPatch->GetDataState(iDataInitial, DataLocation_REdge);

		DataArray4D<double> & dataUpdateREdge =
			pPatch->GetDataState(iDataUpdate, DataLocation_REdge);

		DataArray4D<double> & dataInitialTracer =
			pPatch->GetDataTracers(iDataInitial);

		DataArray4D<double> & dataUpdateTracer =
			pPatch->GetDataTracers(iDataUpdate);

		// Metric quantities
		const DataArray4D<double> & dDerivRNode =
			pPatch->GetDerivRNode();

		const DataArray4D<double> & dContraMetricXi =
			pPatch->GetContraMetricXi();

		const DataArray4D<double> & dDerivRREdge =
			pPatch->GetDerivRREdge();

		const DataArray4D<double> & dContraMetricXiREdge =
			pPatch->GetContraMetricXiREdge();

#if defined(EXPLICIT_VERTICAL_VELOCITY_ADVECTION) && \
    defined(VERTICAL_VELOCITY_ADVECTION_CLARK)
		const DataArray4D<double> & dContraMetricA =
			pPatch->GetContraMetricA();

		const DataArray4D<double> & dContraMetricB =
			pPatch->GetContraMetricB();

		const DataArray4D<double> & dContraMetricAREdge =
			pPatch->GetContraMetricAREdge();

		const DataArray4D<double> & dContraMetricBREdge =
			pPatch->GetContraMetricBREdge();
#endif
#if defined(VERTICAL_UPWINDING) \
 || defined(EXPLICIT_VERTICAL_VELOCITY_ADVECTION)
		// Interpolate U and V to interfaces
		pPatch->InterpolateNodeToREdge(UIx, iDataInitial);
		pPatch->InterpolateNodeToREdge(VIx, iDataInitial);
#endif
#if defined(EXPLICIT_THERMO) \
 || (defined(EXPLICIT_VERTICAL_VELOCITY_ADVECTION) \
     && defined(VERTICAL_VELOCITY_ADVECTION_CLARK))
		// Interpolate W to levels
		pPatch->InterpolateREdgeToNode(WIx, iDataInitial);
#endif

		// Loop over all nodes
		for (int i = box.GetAInteriorBegin(); i < box.GetAInteriorEnd(); i++) {
		for (int j = box.GetBInteriorBegin(); j < box.GetBInteriorEnd(); j++) {

                        //UPDATE ONLY THE TERMS TREATED IMPLICITLY BUT EXPLICITLY
			int iA = i;
			int iB = j;

			SetupReferenceColumn(
				pPatch, iA, iB,
				dataRefNode,
				dataInitialNode,
				dataRefREdge,
				dataInitialREdge);

			Evaluate(m_dColumnState, m_dSoln);

			// Apply update to P
			if (pGrid->GetVarLocation(PIx) == DataLocation_REdge) {
				for (int k = 0; k <= nRElements; k++) {
					dataUpdateREdge[PIx][k][iA][iB] -=
						dDeltaT * m_dSoln[VecFIx(FPIx, k)];
				}
			} else {
				for (int k = 0; k < nRElements; k++) {
					dataUpdateNode[PIx][k][iA][iB] -=
						dDeltaT * m_dSoln[VecFIx(FPIx, k)];
				}
			}

			// Apply update to W
			if (pGrid->GetVarLocation(WIx) == DataLocation_REdge) {
				for (int k = 0; k <= nRElements; k++) {
					dataUpdateREdge[WIx][k][iA][iB] -=
						dDeltaT * m_dSoln[VecFIx(FWIx, k)];
				}

			} else {
				for (int k = 0; k < nRElements; k++) {
					dataUpdateNode[WIx][k][iA][iB] -=
						dDeltaT * m_dSoln[VecFIx(FWIx, k)];
				}
			}

			// Apply update to Rho
			if (pGrid->GetVarLocation(RIx) == DataLocation_REdge) {
				for (int k = 0; k <= nRElements; k++) {
					dataUpdateREdge[RIx][k][iA][iB] -=
						dDeltaT * m_dSoln[VecFIx(FRIx, k)];
				}
			} else {
				for (int k = 0; k < nRElements; k++) {
					dataUpdateNode[RIx][k][iA][iB] -=
						dDeltaT * m_dSoln[VecFIx(FRIx, k)];
				}
			}

			// Update tracers in column
			UpdateColumnTracers(
				dDeltaT,
				dataInitialNode,
				dataUpdateNode,
				dataInitialREdge,
				dataUpdateREdge,
				dataInitialTracer,
				dataUpdateTracer);
                }
                }
        }             
}

///////////////////////////////////////////////////////////////////////////////
void VerticalDynamicsFEM::StepExplicit(
	int iDataInitial,
	int iDataUpdate,
	const Time & time,
	double dDeltaT
) {
	// Get a copy of the grid
	GridGLL * pGrid = dynamic_cast<GridGLL *>(m_model.GetGrid());

	// Physical constants
	const PhysicalConstants & phys = m_model.GetPhysicalConstants();

	// Indices of EquationSet variables
	const int UIx = 0;
	const int VIx = 1;
	const int PIx = 2;
	const int WIx = 3;
	const int RIx = 4;

	// Number of elements
	const int nRElements = pGrid->GetRElements();

	// Number of finite elements in the vertical
	const int nFiniteElements = nRElements / m_nVerticalOrder;

	// Store timestep size
	m_dDeltaT = dDeltaT;

	// Reset the reference state
	memset(m_dStateRefNode[WIx],  0,  nRElements   *sizeof(double));
	memset(m_dStateRefREdge[WIx], 0, (nRElements+1)*sizeof(double));

	// Perform local update
	for (int n = 0; n < pGrid->GetActivePatchCount(); n++) {
		GridPatch * pPatch = pGrid->GetActivePatch(n);

		const PatchBox & box = pPatch->GetPatchBox();

		// State Data
		const DataArray4D<double> & dataRefNode =
			pPatch->GetReferenceState(DataLocation_Node);

		const DataArray4D<double> & dataInitialNode =
			pPatch->GetDataState(iDataInitial, DataLocation_Node);

		DataArray4D<double> & dataUpdateNode =
			pPatch->GetDataState(iDataUpdate, DataLocation_Node);

		const DataArray4D<double> & dataRefREdge =
			pPatch->GetReferenceState(DataLocation_REdge);

		const DataArray4D<double> & dataInitialREdge =
			pPatch->GetDataState(iDataInitial, DataLocation_REdge);

		DataArray4D<double> & dataUpdateREdge =
			pPatch->GetDataState(iDataUpdate, DataLocation_REdge);

		DataArray4D<double> & dataInitialTracer =
			pPatch->GetDataTracers(iDataInitial);

		DataArray4D<double> & dataUpdateTracer =
			pPatch->GetDataTracers(iDataUpdate);

		// Metric quantities
		const DataArray4D<double> & dDerivRNode =
			pPatch->GetDerivRNode();

		const DataArray4D<double> & dContraMetricXi =
			pPatch->GetContraMetricXi();

		const DataArray4D<double> & dDerivRREdge =
			pPatch->GetDerivRREdge();

		const DataArray4D<double> & dContraMetricXiREdge =
			pPatch->GetContraMetricXiREdge();

#if defined(EXPLICIT_VERTICAL_VELOCITY_ADVECTION) && \
    defined(VERTICAL_VELOCITY_ADVECTION_CLARK)
		const DataArray4D<double> & dContraMetricA =
			pPatch->GetContraMetricA();

		const DataArray4D<double> & dContraMetricB =
			pPatch->GetContraMetricB();

		const DataArray4D<double> & dContraMetricAREdge =
			pPatch->GetContraMetricAREdge();

		const DataArray4D<double> & dContraMetricBREdge =
			pPatch->GetContraMetricBREdge();
#endif
#if defined(VERTICAL_UPWINDING) \
 || defined(EXPLICIT_VERTICAL_VELOCITY_ADVECTION)
		// Interpolate U and V to interfaces
		pPatch->InterpolateNodeToREdge(UIx, iDataInitial);
		pPatch->InterpolateNodeToREdge(VIx, iDataInitial);
#endif
#if defined(EXPLICIT_THERMO) \
 || (defined(EXPLICIT_VERTICAL_VELOCITY_ADVECTION) \
     && defined(VERTICAL_VELOCITY_ADVECTION_CLARK))
		// Interpolate W to levels
		pPatch->InterpolateREdgeToNode(WIx, iDataInitial);
#endif

		// Loop over all nodes
		for (int i = box.GetAInteriorBegin(); i < box.GetAInteriorEnd(); i++) {
		for (int j = box.GetBInteriorBegin(); j < box.GetBInteriorEnd(); j++) {

			// Update thermodynamic variables
			if (m_fFullyExplicit) {

				int iA = i;
				int iB = j;

				SetupReferenceColumn(
					pPatch, iA, iB,
					dataRefNode,
					dataInitialNode,
					dataRefREdge,
					dataInitialREdge);

				Evaluate(m_dColumnState, m_dSoln);

				// Apply update to P
				if (pGrid->GetVarLocation(PIx) == DataLocation_REdge) {
					for (int k = 0; k <= nRElements; k++) {
						dataUpdateREdge[PIx][k][iA][iB] -=
							dDeltaT * m_dSoln[VecFIx(FPIx, k)];
					}
				} else {
					for (int k = 0; k < nRElements; k++) {
						dataUpdateNode[PIx][k][iA][iB] -=
							dDeltaT * m_dSoln[VecFIx(FPIx, k)];
					}
				}

				// Apply update to W
				if (pGrid->GetVarLocation(WIx) == DataLocation_REdge) {
					for (int k = 0; k <= nRElements; k++) {
						dataUpdateREdge[WIx][k][iA][iB] -=
							dDeltaT * m_dSoln[VecFIx(FWIx, k)];
					}

				} else {
					for (int k = 0; k < nRElements; k++) {
						dataUpdateNode[WIx][k][iA][iB] -=
							dDeltaT * m_dSoln[VecFIx(FWIx, k)];
					}
				}

				// Apply update to Rho
				if (pGrid->GetVarLocation(RIx) == DataLocation_REdge) {
					for (int k = 0; k <= nRElements; k++) {
						dataUpdateREdge[RIx][k][iA][iB] -=
							dDeltaT * m_dSoln[VecFIx(FRIx, k)];
					}
				} else {
					for (int k = 0; k < nRElements; k++) {
						dataUpdateNode[RIx][k][iA][iB] -=
							dDeltaT * m_dSoln[VecFIx(FRIx, k)];
					}
				}

				// Update tracers in column
				UpdateColumnTracers(
					dDeltaT,
					dataInitialNode,
					dataUpdateNode,
					dataInitialREdge,
					dataUpdateREdge,
					dataInitialTracer,
					dataUpdateTracer);
/*
				// Check boundary condition
				{
					const DataArray4D<double> & dContraMetricXi =
						pPatch->GetContraMetricXi();
					const DataArray4D<double> & dDerivRNode =
						pPatch->GetDerivRNode();

					double dConUx =
						dContraMetricXi[0][iA][iB][0]
							* dataUpdateNode[UIx][0][iA][iB]
						+ dContraMetricXi[0][iA][iB][1]
							* dataUpdateNode[VIx][0][iA][iB]
						+ dContraMetricXi[0][iA][iB][2]
							* dDerivRNode[0][iA][iB][2]
							* dataUpdateNode[WIx][0][iA][iB];

					if (fabs(dConUx) > 1.0e-12) {
						printf("%1.15e\n", dConUx);
						_EXCEPTION();
					}
				}
*/
			}

#if defined(EXPLICIT_VERTICAL_VELOCITY_ADVECTION)
#if defined(VERTICAL_VELOCITY_ADVECTION_CLARK)
			// Explicit vertical velocity advection (Clark form)
			{
				// Calculate specific kinetic energy on model levels
				for (int k = 0; k < nRElements; k++) {
					m_dStateNode[UIx][k] = dataInitialNode[UIx][k][i][j];
					m_dStateNode[VIx][k] = dataInitialNode[VIx][k][i][j];

					double dCovUa = dataInitialNode[UIx][k][i][j];
					double dCovUb = dataInitialNode[VIx][k][i][j];
					double dCovUx = dataInitialNode[WIx][k][i][j]
						* dDerivRNode[k][i][j][2];

					double dConUa =
						  dContraMetricA[k][i][j][0] * dCovUa
						+ dContraMetricA[k][i][j][1] * dCovUb
						+ dContraMetricA[k][i][j][2] * dCovUx;

					double dConUb =
						  dContraMetricB[k][i][j][0] * dCovUa
						+ dContraMetricB[k][i][j][1] * dCovUb
						+ dContraMetricB[k][i][j][2] * dCovUx;

					double dConUx =
						  dContraMetricXi[k][i][j][0] * dCovUa
						+ dContraMetricXi[k][i][j][1] * dCovUb
						+ dContraMetricXi[k][i][j][2] * dCovUx;

					m_dKineticEnergyNode[k] = 0.5 * (
						  dConUa * dCovUa
						+ dConUb * dCovUb
						+ dConUx * dCovUx);
				}

				// Stride through state matrix
				int nVerticalStateStride =
					dataInitialNode.GetSize(2)
					* dataInitialNode.GetSize(3);

				// Update vertical velocity on levels
				if (pGrid->GetVarLocation(WIx) == DataLocation_Node) {
					pGrid->DifferentiateNodeToNode(
						m_dKineticEnergyNode,
						m_dDiffKineticEnergyNode);

					pGrid->DifferentiateNodeToNode(
						m_dStateNode[UIx],
						m_dDiffUa);

					pGrid->DifferentiateNodeToNode(
						m_dStateNode[VIx],
						m_dDiffUb);

					for (int k = 0; k < nRElements; k++) {
						double dCovUa = dataInitialNode[UIx][k][i][j];
						double dCovUb = dataInitialNode[VIx][k][i][j];
						double dCovUx = dataInitialNode[WIx][k][i][j]
							* dDerivRNode[k][i][j][2];

						double dConUa =
							  dContraMetricA[k][i][j][0] * dCovUa
							+ dContraMetricA[k][i][j][1] * dCovUb
							+ dContraMetricA[k][i][j][2] * dCovUx;
	
						double dConUb =
							  dContraMetricB[k][i][j][0] * dCovUa
							+ dContraMetricB[k][i][j][1] * dCovUb
							+ dContraMetricB[k][i][j][2] * dCovUx;

						double dCurlTerm =
							- dConUa * m_dDiffUa[k]
							- dConUb * m_dDiffUb[k];

						dataUpdateNode[WIx][k][i][j] -=
							dDeltaT
							* (m_dDiffKineticEnergyNode[k] + dCurlTerm)
							/ dDerivRNode[k][i][j][2];
					}

				} else {
					pGrid->DifferentiateNodeToREdge(
						m_dKineticEnergyNode,
						m_dDiffKineticEnergyREdge);
	
					pGrid->DifferentiateNodeToREdge(
						m_dStateNode[UIx],
						m_dDiffUa);

					pGrid->DifferentiateNodeToREdge(
						m_dStateNode[VIx],
						m_dDiffUb);

					for (int k = 1; k < nRElements; k++) {
						double dCovUa = dataInitialREdge[UIx][k][i][j];
						double dCovUb = dataInitialREdge[VIx][k][i][j];
						double dCovUx = dataInitialREdge[WIx][k][i][j]
							* dDerivRREdge[k][i][j][2];

						double dConUa =
							  dContraMetricAREdge[k][i][j][0] * dCovUa
							+ dContraMetricAREdge[k][i][j][1] * dCovUb
							+ dContraMetricAREdge[k][i][j][2] * dCovUx;

						double dConUb =
							  dContraMetricBREdge[k][i][j][0] * dCovUa
							+ dContraMetricBREdge[k][i][j][1] * dCovUb
							+ dContraMetricBREdge[k][i][j][2] * dCovUx;

						double dCurlTerm =
							- dConUa * m_dDiffUa[k]
							- dConUb * m_dDiffUb[k];
/*
						// DEBUG
						double dValue =
							(m_dDiffKineticEnergyREdge[k] + dCurlTerm)
							/ dDerivRREdge[k][i][j][2];
						if (fabs(dValue) > 1.0e-5) {
							printf("V: %1.5e\n", dValue);

							double dDxW =
								pGrid->DifferentiateREdgeToREdge(
									&(dataInitialREdge[WIx][0][i][j]),
									k,
									nVerticalStateStride);

							double dXiDotREdge =
								  dContraMetricXiREdge[k][i][j][0] * dCovUa
								+ dContraMetricXiREdge[k][i][j][1] * dCovUb
								+ dContraMetricXiREdge[k][i][j][2] * dCovUx;

							printf("UxDxW: %1.5e\n",
								dXiDotREdge * dDxW);
						}
*/
						dataUpdateREdge[WIx][k][i][j] -=
							dDeltaT
							* (m_dDiffKineticEnergyREdge[k] + dCurlTerm)
							/ dDerivRREdge[k][i][j][2];
					}
				}
			}
#else
			// Explicit vertical velocity advection (advective form)
			{
				// Stride through state matrix
				int nVerticalStateStride =
					dataInitialNode.GetSize(2)
					* dataInitialNode.GetSize(3);

				// Update vertical velocity on levels
				if (pGrid->GetVarLocation(WIx) == DataLocation_Node) {
					_EXCEPTIONT("Not implemented");

				// Update vertical velocity on interfaces
				} else {
					for (int k = 1; k < nRElements; k++) {
						double dCovUa = dataInitialREdge[UIx][k][i][j];
						double dCovUb = dataInitialREdge[VIx][k][i][j];

						double dCovUx =
							  dataInitialREdge[WIx][k][i][j]
							* dDerivRREdge[k][i][j][2];

						double dDxW =
							pGrid->DifferentiateREdgeToREdge(
								&(dataInitialREdge[WIx][0][i][j]),
								k,
								nVerticalStateStride);

						double dXiDotREdge =
							  dContraMetricXiREdge[k][i][j][0] * dCovUa
							+ dContraMetricXiREdge[k][i][j][1] * dCovUb
							+ dContraMetricXiREdge[k][i][j][2] * dCovUx;

						dataUpdateREdge[WIx][k][i][j] -=
							dDeltaT * dXiDotREdge * dDxW;
					}
				}
			}
#endif
#endif

#if defined(EXPLICIT_THERMO) && defined(FORMULATION_THETA)
			// Explicit update of thermodynamic equation
			{
				// Calculate u^xi on model levels
				if (!m_fFullyExplicit) {
					for (int k = 0; k < nRElements; k++) {
						double dCovUa = dataInitialNode[UIx][k][i][j];
						double dCovUb = dataInitialNode[VIx][k][i][j];

						double dCovUx =
							  dataInitialNode[WIx][k][i][j]
							* dDerivRNode[k][i][j][2];

						m_dXiDotNode[k] =
							  dContraMetricXi[k][i][j][0] * dCovUa
							+ dContraMetricXi[k][i][j][1] * dCovUb
							+ dContraMetricXi[k][i][j][2] * dCovUx;

						m_dStateNode[UIx][k] = dCovUa;
						m_dStateNode[VIx][k] = dCovUb;
					}
				}

				int nUpwindStride =
					dataInitialNode.GetSize(2)
					* dataInitialNode.GetSize(3);

				// Differentiate theta
				const LinearColumnDiffFEM & opDiffNodeToNode =
					pGrid->GetOpDiffNodeToNode();

				opDiffNodeToNode.Apply(
					&(dataInitialNode[PIx][0][i][j]),
					&(m_dDiffThetaNode[0]),
					nUpwindStride,
					1);

				// Calculate update to theta
				for (int k = 0; k < nRElements; k++) {
					dataUpdateNode[PIx][k][i][j] -=
						dDeltaT * m_dXiDotNode[k] * m_dDiffThetaNode[k];
				}
			}
#endif
#ifdef VERTICAL_UPWINDING
			// Apply upwinding (discontinuous penalization)
			{

				// Calculate u^xi on model interfaces (if not done above)
				if (!m_fFullyExplicit) {
					for (int k = 0; k <= nRElements; k++) {
						double dCovUa = dataInitialREdge[UIx][k][i][j];
						double dCovUb = dataInitialREdge[VIx][k][i][j];

						double dCovUx =
							  dataInitialREdge[WIx][k][i][j]
							* dDerivRREdge[k][i][j][2];

						m_dXiDotREdge[k] =
							  dContraMetricXiREdge[k][i][j][0] * dCovUa
							+ dContraMetricXiREdge[k][i][j][1] * dCovUb
							+ dContraMetricXiREdge[k][i][j][2] * dCovUx;
					}
				}

				// Calculate weights
				for (int a = 0; a < nFiniteElements - 1; a++) {
					int k = (a+1) * m_nVerticalOrder;
					m_dUpwindWeights[a] =
						0.5 * dDeltaT * fabs(m_dXiDotREdge[k]);
				}

				// Apply upwinding
				const LinearColumnDiscPenaltyFEM & opPenalty =
					pGrid->GetOpPenaltyNodeToNode();

				int nUpwindStride =
					dataInitialNode.GetSize(2)
					* dataInitialNode.GetSize(3);

				// Apply upwinding to U and V
				if (m_fUpwind[UIx]) {
					opPenalty.Apply(
						&(m_dUpwindWeights[0]),
						&(dataInitialNode[UIx][0][i][j]),
						&(dataUpdateNode[UIx][0][i][j]),
						nUpwindStride,
						nUpwindStride);

					opPenalty.Apply(
						&(m_dUpwindWeights[0]),
						&(dataInitialNode[VIx][0][i][j]),
						&(dataUpdateNode[VIx][0][i][j]),
						nUpwindStride,
						nUpwindStride);
				}

#ifdef EXPLICIT_THERMO
				// Apply upwinding to thermodynamic variable
				if ((m_fUpwind[PIx]) && (!m_fFullyExplicit)) {
					if (pGrid->GetVarLocation(PIx) == DataLocation_REdge) {
						_EXCEPTIONT("Not implemented: Upwinding of thermo"
							   " on interfaces");
					} else {
						opPenalty.Apply(
							&(m_dUpwindWeights[0]),
							&(dataInitialNode[PIx][0][i][j]),
							&(dataUpdateNode[PIx][0][i][j]),
							nUpwindStride,
							nUpwindStride);
					}
				}
#endif
			}
#endif

#if defined(VERTICAL_HYPERVISCOSITY) || defined(UNIFORM_DIFFUSION)
			// Apply hyperviscosity to U and V
			if (m_fUpwind[UIx] && (m_nHypervisOrder > 0)) {

				// Calculate u^xi on model levels (if not done above)
				if (!m_fFullyExplicit) {
					for (int k = 0; k < nRElements; k++) {
						double dCovUa = dataInitialNode[UIx][k][i][j];
						double dCovUb = dataInitialNode[VIx][k][i][j];

						double dCovUx =
							  dataInitialNode[WIx][k][i][j]
							* dDerivRNode[k][i][j][2];

						m_dXiDotNode[k] =
							  dContraMetricXi[k][i][j][0] * dCovUa
							+ dContraMetricXi[k][i][j][1] * dCovUb
							+ dContraMetricXi[k][i][j][2] * dCovUx;

						m_dStateNode[UIx][k] = dCovUa;
						m_dStateNode[VIx][k] = dCovUb;
					}
				}

				// Second derivatives of horizontal velocity on model levels
				pGrid->DiffDiffNodeToNode(
					m_dStateNode[UIx],
					m_dHyperDiffState[UIx]);

				pGrid->DiffDiffNodeToNode(
					m_dStateNode[VIx],
					m_dHyperDiffState[VIx]);

#if defined(UNIFORM_DIFFUSION)
				// Apply uniform diffusion in the vertical
				double dZtop = pGrid->GetZtop();

				double dUniformDiffusionCoeff = 75.0 / (dZtop * dZtop);

				for (int k = 0; k < nRElements; k++) {
					dataUpdateNode[UIx][k][i][j] +=
						dDeltaT
						* dUniformDiffusionCoeff
						* m_dHyperDiffState[UIx][k];

					dataUpdateNode[VIx][k][i][j] +=
						dDeltaT
						* dUniformDiffusionCoeff
						* m_dHyperDiffState[VIx][k];
				}
#endif
				// Compute higher derivatives of u and v used for
				// hyperviscosity
				for (int h = 2; h < m_nHypervisOrder; h += 2) {
					memcpy(
						m_dStateAux,
						m_dHyperDiffState[UIx],
						nRElements * sizeof(double));

					pGrid->DiffDiffNodeToNode(
						m_dStateAux,
						m_dHyperDiffState[UIx]
					);

					memcpy(
						m_dStateAux,
						m_dHyperDiffState[VIx],
						nRElements * sizeof(double));

					pGrid->DiffDiffNodeToNode(
						m_dStateAux,
						m_dHyperDiffState[VIx]
					);
				}

				// Apply hyperviscosity
				for (int k = 0; k < nRElements; k++) {
					dataUpdateNode[UIx][k][i][j] +=
						dDeltaT
						* m_dHypervisCoeff
						* fabs(m_dXiDotNode[k])
						* m_dHyperDiffState[UIx][k];

					dataUpdateNode[VIx][k][i][j] +=
						dDeltaT
						* m_dHypervisCoeff
						* fabs(m_dXiDotNode[k])
						* m_dHyperDiffState[VIx][k];
				}
			}
#endif
		}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void VerticalDynamicsFEM::BootstrapJacobian() {

	static const double Epsilon = 1.0e-5;

	int nDim = m_dColumnState.GetRows();

	DataArray2D<double> dJacobian(nDim, nDim);
	DataArray1D<double> dJC(nDim);
	DataArray1D<double> dG(nDim);
	DataArray1D<double> dJCref(nDim);

	Evaluate(m_dColumnState, dJCref);

	for (int i = 0; i < m_dColumnState.GetRows(); i++) {
		dG = m_dColumnState;
		dG[i] = dG[i] + Epsilon;

		Evaluate(dG, dJC);

		for (int j = 0; j < m_dColumnState.GetRows(); j++) {
			dJacobian[i][j] = (dJC[j] - dJCref[j]) / Epsilon;
		}
	}

	std::cout << "DeltaT: " << m_dDeltaT << std::endl;

	FILE * fp;
	fp = fopen("DGRef.txt", "w");
	for (int i = 0; i < nDim; i++) {
		for (int j = 0; j < nDim; j++) {
			fprintf(fp, "%1.15e", dJacobian[i][j]);
			if (j != nDim-1) {
				fprintf(fp, " ");
			}
		}
		fprintf(fp, "\n");
	}
	fclose(fp);

	fp = fopen("G.txt", "w");
	for (int i = 0; i < nDim; i++) {
		fprintf(fp, "%1.15e", dJCref[i]);
		if (i != nDim-1) {
			fprintf(fp, "\n");
		}
	}
	fclose(fp);

	BuildJacobianF(m_dSoln, &(dJacobian[0][0]));

	fp = fopen("DG.txt", "w");
	for (int i = 0; i < nDim; i++) {
		for (int j = 0; j < nDim; j++) {
			fprintf(fp, "%1.15e", dJacobian[i][j]);
			if (j != nDim-1) {
				fprintf(fp, " ");
			}
		}
		fprintf(fp, "\n");
	}
	fclose(fp);

	_EXCEPTION();
}

///////////////////////////////////////////////////////////////////////////////

void VerticalDynamicsFEM::StepImplicit(
	int iDataInitial,
	int iDataUpdate,
	const Time & time,
	double dDeltaT
) {
	// If fully explicit do nothing
	if (m_fFullyExplicit) {
		return;
	}

	// Get a copy of the grid
	Grid * pGrid = m_model.GetGrid();

	// Physical constants
	const PhysicalConstants & phys = m_model.GetPhysicalConstants();

	// Indices of EquationSet variables
	const int UIx = 0;
	const int VIx = 1;
	const int PIx = 2;
	const int WIx = 3;
	const int RIx = 4;

	// Store timestep size
	m_dDeltaT = dDeltaT;

	// Perform local update
	for (int n = 0; n < pGrid->GetActivePatchCount(); n++) {
		GridPatch * pPatch = pGrid->GetActivePatch(n);

		const PatchBox & box = pPatch->GetPatchBox();

		// Contravariant metric components
		const DataArray4D<double> & dContraMetricA =
			pPatch->GetContraMetricA();
		const DataArray4D<double> & dContraMetricB =
			pPatch->GetContraMetricB();

		// State Data
		const DataArray4D<double> & dataRefNode =
			pPatch->GetReferenceState(DataLocation_Node);

		const DataArray4D<double> & dataInitialNode =
			pPatch->GetDataState(iDataInitial, DataLocation_Node);

		DataArray4D<double> & dataUpdateNode =
			pPatch->GetDataState(iDataUpdate, DataLocation_Node);

		const DataArray4D<double> & dataRefREdge =
			pPatch->GetReferenceState(DataLocation_REdge);

		const DataArray4D<double> & dataInitialREdge =
			pPatch->GetDataState(iDataInitial, DataLocation_REdge);

		DataArray4D<double> & dataUpdateREdge =
			pPatch->GetDataState(iDataUpdate, DataLocation_REdge);

		// Tracer Data
		DataArray4D<double> & dataInitialTracer =
			pPatch->GetDataTracers(iDataInitial);

		DataArray4D<double> & dataUpdateTracer =
			pPatch->GetDataTracers(iDataUpdate);

		// Number of tracers
		const int nTracerCount = dataInitialTracer.GetSize(0);

		// Number of finite elements
		int nAElements =
			box.GetAInteriorWidth() / m_nHorizontalOrder;
		int nBElements =
			box.GetBInteriorWidth() / m_nHorizontalOrder;

		// Loop over all nodes, but only perform calculation on shared
		// nodes once
		for (int a = 0; a < nAElements; a++) {
		for (int b = 0; b < nBElements; b++) {

			int iEnd;
			int jEnd;

			if (a == nAElements-1) {
				iEnd = m_nHorizontalOrder;
			} else {
				iEnd = m_nHorizontalOrder-1;
			}

			if (b == nBElements-1) {
				jEnd = m_nHorizontalOrder;
			} else {
				jEnd = m_nHorizontalOrder-1;
			}

		for (int i = 0; i < iEnd; i++) {
		for (int j = 0; j < jEnd; j++) {

			int iA = box.GetAInteriorBegin() + a * m_nHorizontalOrder + i;
			int iB = box.GetBInteriorBegin() + b * m_nHorizontalOrder + j;

			SetupReferenceColumn(
				pPatch, iA, iB,
				dataRefNode,
				dataInitialNode,
				dataRefREdge,
				dataInitialREdge);

#ifdef USE_JACOBIAN_DEBUG
			BootstrapJacobian();
#endif
#ifdef USE_JFNK_PETSC
			// Use PetSc to solve
			double * dX;
			VecGetArray(m_vecX, &dX);
			memcpy(dX, m_dColumnState, m_nColumnStateSize * sizeof(double));
			VecRestoreArray(m_vecX, &dX);

			// Solve
			PetscErrorCode ierr;
			SNESSolve(m_snes, NULL, m_vecX);

			SNESConvergedReason reason;
			SNESGetConvergedReason(m_snes, &reason);
			if ((reason < 0) && (reason != (-5))) {
				_EXCEPTION1("PetSc solver failed to converge (%i)", reason);
			}
/*
			{
				// Evaluate
				DataArray1D<double> dEval;
				dEval.Initialize(m_dColumnState.GetRows());
				VecGetArray(m_vecX, &dX);
				Evaluate(dX, dEval);

				nAvgValues++;
				double dEvalNorm = 0.0;
				for (int n = 0; n < dEval.GetRows(); n++) {
					dEvalNorm += dEval[n] * dEval[n];
					std::cout << dEval[n] << std::endl;
				}
				dAvgNorm += sqrt(dEvalNorm / dEval.GetRows());

				VecRestoreArray(m_vecX, &dX);
				_EXCEPTION();
			}
*/
			VecGetArray(m_vecX, &dX);
			memcpy(m_dSoln, dX, m_nColumnStateSize * sizeof(double));
			VecRestoreArray(m_vecX, &dX);
#endif
#ifdef USE_JFNK_GMRES
			// Use Jacobian-Free Newton-Krylov to solve
			m_dSoln = m_dColumnState;

			double dError =
				PerformJFNK_NewtonStep_Safe(
				//PerformBICGSTAB_NewtonStep_Safe(
					m_dSoln,
					m_dSoln.GetRows(),
					1.0e-8);

/*
			bool fConverged =
				PerformJFNK(
					m_dSoln,
					m_dSoln.GetRows(),
					1.0e-12,
					m_dSoln.GetRows(),
					1.0e-12);

			if (!fConverged) {
				_EXCEPTIONT("Convergence failure");
			}
*/
			// DEBUG (check for NANs in output)
			if (!(m_dSoln[0] == m_dSoln[0])) {
                DataArray1D<double> dEval;
                dEval.Allocate(m_dColumnState.GetRows());
                Evaluate(m_dSoln, dEval);

                for (int p = 0; p < dEval.GetRows(); p++) {
                    printf("%1.15e %1.15e %1.15e\n",
						dEval[p], m_dSoln[p] - m_dColumnState[p], m_dColumnState[p]);
                }
				for (int p = 0; p < m_dExnerRefREdge.GetRows(); p++) {
					printf("%1.15e %1.15e\n",
						m_dExnerRefREdge[p], dataRefREdge[RIx][p][iA][iB]);
				}
                _EXCEPTIONT("Inversion failure");
            }

#endif
#ifdef USE_DIRECTSOLVE_APPROXJ
			static const double Epsilon = 1.0e-5;

			// Prepare the column
			PrepareColumn(m_dColumnState);

			// Build the F vector
			BuildF(m_dColumnState, m_dSoln);

			DataArray1D<double> dJC;
			dJC.Allocate(m_dColumnState.GetRows());

			DataArray1D<double> dG;
			dG.Allocate(m_dColumnState.GetRows());

			DataArray1D<double> dJCref;
			dJCref.Allocate(m_dColumnState.GetRows());

			Evaluate(m_dColumnState, dJCref);

			for (int i = 0; i < m_dColumnState.GetRows(); i++) {
				dG = m_dColumnState;
				dG[i] = dG[i] + Epsilon;

				Evaluate(dG, dJC);

				for (int j = 0; j < m_dColumnState.GetRows(); j++) {
					m_matJacobianF[i][j] = (dJC[j] - dJCref[j]) / Epsilon;
				}
			}

			// Use direct solver
			LAPACK::DGESV(m_matJacobianF, m_dSoln, m_vecIPiv);

			for (int k = 0; k < m_dSoln.GetRows(); k++) {
				m_dSoln[k] = m_dColumnState[k] - m_dSoln[k];
			}
#endif
#ifdef USE_DIRECTSOLVE
			// Prepare the column
			PrepareColumn(m_dColumnState);

			// Build the F vector
			BuildF(m_dColumnState, m_dSoln);

			// Build the Jacobian
			BuildJacobianF(m_dColumnState, &(m_matJacobianF[0][0]));

#ifdef USE_JACOBIAN_GENERAL
/*
			FILE * fpDF = fopen("StateDF.txt", "w");
			for (int k = 0; k < 3*(m_nRElements+1); k++) {
				for (int i = 0; i < 3*(m_nRElements+1); i++) {
					fprintf(fpDF, "%1.15e\t", m_matJacobianF[k][i]);
				}
				fprintf(fpDF, "\n");
			}
			fclose(fpDF);

			FILE * fpF = fopen("StateF.txt", "w");
			for (int k = 0; k < 3*(m_nRElements+1); k++) {
				fprintf(fpF, "%1.15e\n", m_dSoln[k]);
			}
			fclose(fpF);
*/
			// Use direct solver
			LAPACK::DGESV(m_matJacobianF, m_dSoln, m_vecIPiv);
#endif
#ifdef USE_JACOBIAN_DIAGONAL
			// Use diagonal solver
			int iInfo = LAPACK::DGBSV(
				m_matJacobianF, m_dSoln, m_vecIPiv,
				m_nJacobianFKL, m_nJacobianFKU);

			if (iInfo != 0) {
				_EXCEPTION1("Solution failed: %i", iInfo);
			}

			// DEBUG (check for NANs in output)
			if (!(m_dSoln[0] == m_dSoln[0])) {
				DataArray1D<double> dEval;
				dEval.Allocate(m_dColumnState.GetRows());
				Evaluate(m_dSoln, dEval);

				for (int p = 0; p < dEval.GetRows(); p++) {
					printf("%1.15e %1.15e %1.15e\n",
						dEval[p], m_dSoln[p] - m_dColumnState[p], m_dColumnState[p]);
				}
				for (int p = 0; p < m_dExnerRefREdge.GetRows(); p++) {
					printf("%1.15e %1.15e\n",
						m_dExnerRefREdge[p], dataRefREdge[RIx][p][iA][iB]);
				}
				_EXCEPTIONT("Inversion failure");
			}

#endif

			for (int k = 0; k < m_dSoln.GetRows(); k++) {
				m_dSoln[k] = m_dColumnState[k] - m_dSoln[k];
			}
#endif

#ifdef EXPLICIT_THERMO
			// Verify thermodynamic closure is untouched by update
			if (pGrid->GetVarLocation(PIx) == DataLocation_REdge) {
				for (int k = 0; k <= pGrid->GetRElements(); k++) {
					if (fabs(m_dSoln[VecFIx(FPIx, k)] - dataInitialREdge[PIx][k][iA][iB]) > 1.0e-12) {
						_EXCEPTIONT("Logic error");
					}
				}

			} else {
				for (int k = 0; k < pGrid->GetRElements(); k++) {
					if (fabs(m_dSoln[VecFIx(FPIx, k)] - dataInitialNode[PIx][k][iA][iB]) > 1.0e-12) {
						_EXCEPTIONT("Logic error");
					}
				}
			}

#else
			// Apply updated state to thermodynamic closure
			if (pGrid->GetVarLocation(PIx) == DataLocation_REdge) {
				for (int k = 0; k <= pGrid->GetRElements(); k++) {
					dataUpdateREdge[PIx][k][iA][iB] =
						m_dSoln[VecFIx(FPIx, k)];
				}
			} else {
				for (int k = 0; k < pGrid->GetRElements(); k++) {
					dataUpdateNode[PIx][k][iA][iB] =
						m_dSoln[VecFIx(FPIx, k)];
				}
			}
#endif
			// Copy over W
			if (pGrid->GetVarLocation(WIx) == DataLocation_REdge) {
				for (int k = 0; k <= pGrid->GetRElements(); k++) {
					dataUpdateREdge[WIx][k][iA][iB] =
						m_dSoln[VecFIx(FWIx, k)];
				}
			} else {
				for (int k = 0; k < pGrid->GetRElements(); k++) {
					dataUpdateNode[WIx][k][iA][iB] =
						m_dSoln[VecFIx(FWIx, k)];
				}
			}

			// Copy over Rho
			if (pGrid->GetVarLocation(RIx) == DataLocation_REdge) {
				for (int k = 0; k <= pGrid->GetRElements(); k++) {
					dataUpdateREdge[RIx][k][iA][iB] =
						m_dSoln[VecFIx(FRIx, k)];
				}
			} else {
				for (int k = 0; k < pGrid->GetRElements(); k++) {
					dataUpdateNode[RIx][k][iA][iB] =
						m_dSoln[VecFIx(FRIx, k)];
				}
			}

			// Update tracers in column
			UpdateColumnTracers(
				dDeltaT,
				dataInitialNode,
				dataUpdateNode,
				dataInitialREdge,
				dataUpdateREdge,
				dataInitialTracer,
				dataUpdateTracer);
		}
		}

		}
		}

		// Copy over new state on shared nodes (edges of constant alpha)
		for (int a = 1; a < nAElements; a++) {
			int iA = box.GetAInteriorBegin() + a * m_nHorizontalOrder - 1;

			for (int b = 0; b < nBElements; b++) {

				// Top element contains more information
				int jEnd;
				if (b == nBElements-1) {
					jEnd = m_nHorizontalOrder;
				} else {
					jEnd = m_nHorizontalOrder-1;
				}

				// Loop along edges of constant alpha
				for (int j = 0; j < jEnd; j++) {

					int iB = box.GetBInteriorBegin() + b * m_nHorizontalOrder + j;

					for (int k = 0; k < pGrid->GetRElements(); k++) {
						//dataUpdateNode[UIx][k][iA][iB]
						//	= dataUpdateNode[UIx][k][iA+1][iB];
						//dataUpdateNode[VIx][k][iA][iB]
						//	= dataUpdateNode[VIx][k][iA+1][iB];
						dataUpdateNode[PIx][k][iA][iB]
							= dataUpdateNode[PIx][k][iA+1][iB];
						dataUpdateNode[WIx][k][iA][iB]
							= dataUpdateNode[WIx][k][iA+1][iB];
						dataUpdateNode[RIx][k][iA][iB]
							= dataUpdateNode[RIx][k][iA+1][iB];

						for (int c = 0; c < nTracerCount; c++) {
							dataUpdateTracer[c][k][iA][iB]
								= dataUpdateTracer[c][k][iA+1][iB];
						}
					}

					for (int k = 0; k <= pGrid->GetRElements(); k++) {
						//dataUpdateREdge[UIx][k][iA][iB]
						//	= dataUpdateREdge[UIx][k][iA+1][iB];
						//dataUpdateREdge[VIx][k][iA][iB]
						//	= dataUpdateREdge[VIx][k][iA+1][iB];
						dataUpdateREdge[PIx][k][iA][iB]
							= dataUpdateREdge[PIx][k][iA+1][iB];
						dataUpdateREdge[WIx][k][iA][iB]
							= dataUpdateREdge[WIx][k][iA+1][iB];
						dataUpdateREdge[RIx][k][iA][iB]
							= dataUpdateREdge[RIx][k][iA+1][iB];
					}
				}
			}
		}

		// Copy over new state on shared nodes (edges of constant beta)
		for (int b = 1; b < nBElements; b++) {
		for (int i = box.GetAInteriorBegin(); i < box.GetAInteriorEnd(); i++) {
			int iB = box.GetBInteriorBegin() + b * m_nHorizontalOrder - 1;

			for (int k = 0; k < pGrid->GetRElements(); k++) {
				//dataUpdateNode[UIx][k][i][iB]
				//	= dataUpdateNode[UIx][k][i][iB+1];
				//dataUpdateNode[VIx][k][i][iB]
				//	= dataUpdateNode[VIx][k][i][iB+1];
				dataUpdateNode[PIx][k][i][iB]
					= dataUpdateNode[PIx][k][i][iB+1];
				dataUpdateNode[WIx][k][i][iB]
					= dataUpdateNode[WIx][k][i][iB+1];
				dataUpdateNode[RIx][k][i][iB]
					= dataUpdateNode[RIx][k][i][iB+1];

				for (int c = 0; c < nTracerCount; c++) {
					dataUpdateTracer[c][k][i][iB]
						= dataUpdateTracer[c][k][i][iB+1];
				}

			}

			for (int k = 0; k <= pGrid->GetRElements(); k++) {
				//dataUpdateREdge[UIx][k][i][iB]
				//	= dataUpdateREdge[UIx][k][i][iB+1];
				//dataUpdateREdge[VIx][k][i][iB]
				//	= dataUpdateREdge[VIx][k][i][iB+1];
				dataUpdateREdge[PIx][k][i][iB]
					= dataUpdateREdge[PIx][k][i][iB+1];
				dataUpdateREdge[WIx][k][i][iB]
					= dataUpdateREdge[WIx][k][i][iB+1];
				dataUpdateREdge[RIx][k][i][iB]
					= dataUpdateREdge[RIx][k][i][iB+1];
			}
		}
		}
	}

	// Filter negative tracers
	FilterNegativeTracers(iDataUpdate);
}


///////////////////////////////////////////////////////////////////////////////

void VerticalDynamicsFEM::SetupReferenceColumn(
	GridPatch * pPatch,
	int iA,
	int iB,
	const DataArray4D<double> & dataRefNode,
	const DataArray4D<double> & dataInitialNode,
	const DataArray4D<double> & dataRefREdge,
	const DataArray4D<double> & dataInitialREdge
) {

	// Indices of EquationSet variables
	const int UIx = 0;
	const int VIx = 1;
	const int PIx = 2;
	const int WIx = 3;
	const int RIx = 4;

	// Get a copy of the grid
	GridGLL * pGrid = dynamic_cast<GridGLL *>(m_model.GetGrid());

	// Number of radial elements
	const int nRElements = pGrid->GetRElements();

	// Store active patch in index
	m_pPatch = pPatch;
	m_iA = iA;
	m_iB = iB;

	// Store U in State structure
	if (pGrid->GetVarLocation(UIx) == DataLocation_Node) {
		for (int k = 0; k < nRElements; k++) {
			m_dStateNode[UIx][k] = dataInitialNode[UIx][k][iA][iB];
		}

		if (pGrid->GetVarsAtLocation(DataLocation_REdge) != 0) {
			pGrid->InterpolateNodeToREdge(
				m_dStateNode[UIx],
				m_dStateREdge[UIx]);
		}

	} else {
		for (int k = 0; k <= nRElements; k++) {
			m_dStateREdge[UIx][k] = dataInitialREdge[UIx][k][iA][iB];
		}

		pGrid->InterpolateREdgeToNode(
			m_dStateREdge[UIx],
			m_dStateNode[UIx]);
	}

	// Store V in State structure
	if (pGrid->GetVarLocation(VIx) == DataLocation_Node) {
		for (int k = 0; k < nRElements; k++) {
			m_dStateNode[VIx][k] = dataInitialNode[VIx][k][iA][iB];
		}

		if (pGrid->GetVarsAtLocation(DataLocation_REdge) != 0) {
			pGrid->InterpolateNodeToREdge(
				m_dStateNode[VIx],
				m_dStateREdge[VIx]);
		}

	} else {
		for (int k = 0; k <= nRElements; k++) {
			m_dStateREdge[VIx][k] = dataInitialREdge[VIx][k][iA][iB];
		}

		pGrid->InterpolateREdgeToNode(
			m_dStateREdge[VIx],
			m_dStateNode[VIx]);
	}

#if !defined(EXPLICIT_VERTICAL_VELOCITY_ADVECTION) && \
	defined(VERTICAL_VELOCITY_ADVECTION_CLARK)
	// Calculate vertical derivatives of horizontal velocity
	if (pGrid->GetVarLocation(WIx) == DataLocation_Node) {
		pGrid->DifferentiateNodeToNode(
			m_dStateNode[UIx],
			m_dDiffUa);

		pGrid->DifferentiateNodeToNode(
			m_dStateNode[VIx],
			m_dDiffUb);

	} else {
		pGrid->DifferentiateNodeToREdge(
			m_dStateNode[UIx],
			m_dDiffUa);

		pGrid->DifferentiateNodeToREdge(
			m_dStateNode[VIx],
			m_dDiffUb);
	}
#endif

	// Copy over Theta
	if (pGrid->GetVarLocation(PIx) == DataLocation_REdge) {
		for (int k = 0; k <= pGrid->GetRElements(); k++) {
			m_dColumnState[VecFIx(FPIx, k)] =
				dataInitialREdge[PIx][k][iA][iB];
		}
	} else {
		for (int k = 0; k < pGrid->GetRElements(); k++) {
			m_dColumnState[VecFIx(FPIx, k)] =
				dataInitialNode[PIx][k][iA][iB];
		}
	}

	// Copy over W
	if (pGrid->GetVarLocation(WIx) == DataLocation_REdge) {
		for (int k = 0; k <= pGrid->GetRElements(); k++) {
			m_dColumnState[VecFIx(FWIx, k)] =
				dataInitialREdge[WIx][k][iA][iB];
		}
	} else {
		for (int k = 0; k < pGrid->GetRElements(); k++) {
			m_dColumnState[VecFIx(FWIx, k)] =
				dataInitialNode[WIx][k][iA][iB];
		}
	}

	// Copy over rho
	if (pGrid->GetVarLocation(RIx) == DataLocation_REdge) {
		for (int k = 0; k <= pGrid->GetRElements(); k++) {
			m_dColumnState[VecFIx(FRIx, k)] =
				dataInitialREdge[RIx][k][iA][iB];
		}
	} else {
		for (int k = 0; k < pGrid->GetRElements(); k++) {
			m_dColumnState[VecFIx(FRIx, k)] =
				dataInitialNode[RIx][k][iA][iB];
		}
	}
/*
	// Construct reference column
	if (m_fUseReferenceState) {
		for (int k = 0; k < pGrid->GetRElements(); k++) {
			m_dStateRefNode[RIx][k] = dataRefNode[RIx][k][iA][iB];
			m_dStateRefNode[PIx][k] = dataRefNode[PIx][k][iA][iB];
		}
		for (int k = 0; k <= pGrid->GetRElements(); k++) {
			m_dStateRefREdge[RIx][k] = dataRefREdge[RIx][k][iA][iB];
			m_dStateRefREdge[PIx][k] = dataRefREdge[PIx][k][iA][iB];
		}

		// Build the Exner pressure reference
		for (int k = 0; k < pGrid->GetRElements(); k++) {
			m_dExnerRefNode[k] = dataExnerNode[k][iA][iB];
			m_dDiffExnerRefNode[k] = dataDiffExnerNode[k][iA][iB];
		}

		for (int k = 0; k <= pGrid->GetRElements(); k++) {
			m_dExnerRefREdge[k] = dataExnerREdge[k][iA][iB];
			m_dDiffExnerRefREdge[k] = dataDiffExnerREdge[k][iA][iB];
		}
	}
*/

	// Metric terms
	const DataArray3D<double> & dJacobian =
		m_pPatch->GetJacobian();
	const DataArray3D<double> & dElementArea =
		m_pPatch->GetElementArea();
	const DataArray3D<double> & dJacobianREdge =
		m_pPatch->GetJacobianREdge();
	const DataArray4D<double> & dDerivRNode =
		m_pPatch->GetDerivRNode();
	const DataArray4D<double> & dDerivRREdge =
		m_pPatch->GetDerivRREdge();
	const DataArray4D<double> & dContraMetricA =
		m_pPatch->GetContraMetricA();
	const DataArray4D<double> & dContraMetricB =
		m_pPatch->GetContraMetricB();
	const DataArray4D<double> & dContraMetricXi =
		m_pPatch->GetContraMetricXi();
	const DataArray4D<double> & dContraMetricAREdge =
		m_pPatch->GetContraMetricAREdge();
	const DataArray4D<double> & dContraMetricBREdge =
		m_pPatch->GetContraMetricBREdge();
	const DataArray4D<double> & dContraMetricXiREdge =
		m_pPatch->GetContraMetricXiREdge();

	for (int k = 0; k < pGrid->GetRElements(); k++) {
		m_dColumnJacobianNode[k] = dJacobian[k][iA][iB];
		m_dColumnElementArea[k] = dElementArea[k][iA][iB];
		m_dColumnInvJacobianNode[k] = 1.0 / m_dColumnJacobianNode[k];

		m_dColumnDerivRNode[k][0] = dDerivRNode[k][iA][iB][0];
		m_dColumnDerivRNode[k][1] = dDerivRNode[k][iA][iB][1];
		m_dColumnDerivRNode[k][2] = dDerivRNode[k][iA][iB][2];

		m_dColumnContraMetricA[k][0] = dContraMetricA[k][iA][iB][0];
		m_dColumnContraMetricA[k][1] = dContraMetricA[k][iA][iB][1];
		m_dColumnContraMetricA[k][2] = dContraMetricA[k][iA][iB][2];

		m_dColumnContraMetricB[k][0] = dContraMetricB[k][iA][iB][0];
		m_dColumnContraMetricB[k][1] = dContraMetricB[k][iA][iB][1];
		m_dColumnContraMetricB[k][2] = dContraMetricB[k][iA][iB][2];

		m_dColumnContraMetricXi[k][0] = dContraMetricXi[k][iA][iB][0];
		m_dColumnContraMetricXi[k][1] = dContraMetricXi[k][iA][iB][1];
		m_dColumnContraMetricXi[k][2] = dContraMetricXi[k][iA][iB][2];
	}

	for (int k = 0; k <= pGrid->GetRElements(); k++) {
		m_dColumnJacobianREdge[k] = dJacobianREdge[k][iA][iB];
		m_dColumnInvJacobianREdge[k] = 1.0 / m_dColumnJacobianREdge[k];

		m_dColumnDerivRREdge[k][0] = dDerivRREdge[k][iA][iB][0];
		m_dColumnDerivRREdge[k][1] = dDerivRREdge[k][iA][iB][1];
		m_dColumnDerivRREdge[k][2] = dDerivRREdge[k][iA][iB][2];

		m_dColumnContraMetricAREdge[k][0] = dContraMetricAREdge[k][iA][iB][0];
		m_dColumnContraMetricAREdge[k][1] = dContraMetricAREdge[k][iA][iB][1];
		m_dColumnContraMetricAREdge[k][2] = dContraMetricAREdge[k][iA][iB][2];

		m_dColumnContraMetricBREdge[k][0] = dContraMetricBREdge[k][iA][iB][0];
		m_dColumnContraMetricBREdge[k][1] = dContraMetricBREdge[k][iA][iB][1];
		m_dColumnContraMetricBREdge[k][2] = dContraMetricBREdge[k][iA][iB][2];

		m_dColumnContraMetricXiREdge[k][0] = dContraMetricXiREdge[k][iA][iB][0];
		m_dColumnContraMetricXiREdge[k][1] = dContraMetricXiREdge[k][iA][iB][1];
		m_dColumnContraMetricXiREdge[k][2] = dContraMetricXiREdge[k][iA][iB][2];
	}
}

///////////////////////////////////////////////////////////////////////////////

void VerticalDynamicsFEM::PrepareColumn(
	const double * dX
) {
	// Indices of EquationSet variables
	const int UIx = 0;
	const int VIx = 1;
	const int PIx = 2;
	const int WIx = 3;
	const int RIx = 4;

	// Finite element grid
	const GridGLL * pGrid = dynamic_cast<const GridGLL *>(m_model.GetGrid());

	// Physical constants
	const PhysicalConstants & phys = m_model.GetPhysicalConstants();

	// Number of radial elements
	const int nRElements = pGrid->GetRElements();

	// Store state data from this column in state vector
	for (int k = 0; k < nRElements; k++) {

		// Store pressure
		if (pGrid->GetVarLocation(PIx) == DataLocation_Node) {
			m_dStateNode[PIx][k] = dX[VecFIx(FPIx, k)];
		} else {
			m_dStateREdge[PIx][k] = dX[VecFIx(FPIx, k)];
		}

		// Store vertical velocity
		if (pGrid->GetVarLocation(WIx) == DataLocation_Node) {
			m_dStateNode[WIx][k] = dX[VecFIx(FWIx, k)];
		} else {
			m_dStateREdge[WIx][k] = dX[VecFIx(FWIx, k)];
		}

		// Store density
		m_dStateNode[RIx][k] = dX[VecFIx(FRIx, k)];
	}

	if (pGrid->GetVarLocation(WIx) == DataLocation_REdge) {
		m_dStateREdge[WIx][nRElements] = dX[VecFIx(FWIx, nRElements)];
	}
	if (pGrid->GetVarLocation(PIx) == DataLocation_REdge) {
		m_dStateREdge[PIx][nRElements] = dX[VecFIx(FPIx, nRElements)];
	}

	// Vertical velocity on model levels
	if (pGrid->GetVarLocation(WIx) == DataLocation_Node) {

#ifdef FORMULATION_PRESSURE
		// Calculate derivative of P at nodes
		pGrid->DifferentiateNodeToNode(
			m_dStateNode[PIx],
			m_dDiffPNode);
#endif
#ifdef FORMULATION_RHOTHETA_PI
		// Calculate Exner pressure at nodes
		for (int k = 0; k < nRElements; k++) {
			m_dExnerNode[k] =
				phys.ExnerPressureFromRhoTheta(m_dStateNode[PIx][k]);
		}

		// Calculate derivative of Exner pressure at nodes
		pGrid->DifferentiateNodeToNode(
			m_dExnerNode,
			m_dDiffPNode);
#endif
#ifdef FORMULATION_RHOTHETA_P
		// Calculate pressure at nodes
		for (int k = 0; k < nRElements; k++) {
			m_dExnerNode[k] =
				phys.PressureFromRhoTheta(m_dStateNode[PIx][k]);
		}

		// Calculate derivative of Exner pressure at nodes
		pGrid->DifferentiateNodeToNode(
			m_dExnerNode,
			m_dDiffPNode);
#endif
#if defined(FORMULATION_THETA) || defined(FORMULATION_THETA_FLUX)
		// Calculate Exner pressure at nodes
		for (int k = 0; k < nRElements; k++) {
			m_dExnerNode[k] =
				phys.ExnerPressureFromRhoTheta(
					m_dStateNode[RIx][k] * m_dStateNode[PIx][k]);
		}

		// Calculate derivative of Exner pressure at nodes
		pGrid->DifferentiateNodeToNode(
			m_dExnerNode,
			m_dDiffPNode);

		// Theta derivatives on model levels
		if (pGrid->GetVarLocation(PIx) == DataLocation_Node) {
			pGrid->DifferentiateNodeToNode(
				m_dStateNode[PIx],
				m_dDiffThetaNode);
		} else {
			pGrid->DifferentiateREdgeToREdge(
				m_dStateREdge[PIx],
				m_dDiffThetaREdge);
		}

#endif

	// Vertical velocity on model interfaces
	} else {

		// U, V already interpolated in SetupReferenceColumn

		// W is needed on model levels
		pGrid->InterpolateREdgeToNode(
			m_dStateREdge[WIx],
			m_dStateNode[WIx]);

		// Rho are needed on model interfaces
		pGrid->InterpolateNodeToREdge(
			m_dStateNode[RIx],
			m_dStateREdge[RIx]);

#ifdef FORMULATION_PRESSURE
		// Interpolate P to edges
		pGrid->InterpolateNodeToREdge(
			m_dStateNode[PIx],
			m_dStateREdge[PIx]);

		// Calculate derivative of P at edges
		pGrid->DifferentiateNodeToREdge(
			m_dStateNode[PIx],
			m_dDiffPREdge);

		// Calculate derivative of P at nodes
		if (pGrid->GetVarLocation(PIx) == DataLocation_Node) {
			pGrid->DifferentiateNodeToNode(
				m_dStateNode[PIx],
				m_dDiffPNode);
		}
#endif
#ifdef FORMULATION_RHOTHETA_PI
		// Interpolate RhoTheta to edges
		pGrid->InterpolateNodeToREdge(
			m_dStateNode[PIx],
			m_dStateREdge[PIx]);

		// Calculate Exner pressure at nodes
		for (int k = 0; k < nRElements; k++) {
			m_dExnerNode[k] =
				phys.ExnerPressureFromRhoTheta(
					m_dStateNode[PIx][k]);
		}

		// Calculate derivative of Exner pressure at interfaces
		pGrid->DifferentiateNodeToREdge(
			m_dExnerNode,
			m_dDiffPREdge);
#endif
#ifdef FORMULATION_RHOTHETA_P
		// Interpolate RhoTheta to edges
		pGrid->InterpolateNodeToREdge(
			m_dStateNode[PIx],
			m_dStateREdge[PIx]);

		// Calculate pressure at nodes
		for (int k = 0; k < nRElements; k++) {
			m_dExnerNode[k] =
				phys.PressureFromRhoTheta(
					m_dStateNode[PIx][k]);
		}

		// Calculate derivative of pressure at interfaces
		pGrid->DifferentiateNodeToREdge(
			m_dExnerNode,
			m_dDiffPREdge);
#endif
#if defined(FORMULATION_THETA) || defined(FORMULATION_THETA_FLUX)

		// Theta on model levels
		if (pGrid->GetVarLocation(PIx) == DataLocation_Node) {	
			// Theta is needed on model interfaces
			pGrid->InterpolateNodeToREdge(
				m_dStateNode[PIx],
				m_dStateREdge[PIx]);

			// Theta derivatives on model levels
			pGrid->DifferentiateNodeToNode(
				m_dStateNode[PIx],
				m_dDiffThetaNode);

		// Theta on model interfaces
		} else {
			// Theta is needed on model levels
			pGrid->InterpolateREdgeToNode(
				m_dStateREdge[PIx],
				m_dStateNode[PIx]);

			// Theta derivatives on model interfaces
			pGrid->DifferentiateREdgeToREdge(
				m_dStateREdge[PIx],
				m_dDiffThetaREdge);
		}

		// Calculate Exner pressure at nodes
		for (int k = 0; k < nRElements; k++) {
			m_dExnerNode[k] =
				phys.ExnerPressureFromRhoTheta(
					m_dStateNode[RIx][k] * m_dStateNode[PIx][k]);
		}

		// Calculate derivative of Exner pressure at interfaces
		pGrid->DifferentiateNodeToREdge(
			m_dExnerNode,
			m_dDiffPREdge);
#endif
	}

	// Calculate u^xi on model levels
	for (int k = 0; k < nRElements; k++) {
		double dCovUx =
			m_dStateNode[WIx][k] * m_dColumnDerivRNode[k][2];

		m_dXiDotNode[k] =
			  m_dColumnContraMetricXi[k][0] * m_dStateNode[UIx][k]
			+ m_dColumnContraMetricXi[k][1] * m_dStateNode[VIx][k]
			+ m_dColumnContraMetricXi[k][2] * dCovUx;
	}

	// Calculate u^xi on model interfaces
	if ((pGrid->GetVarLocation(WIx) == DataLocation_REdge) ||
	    (pGrid->GetVarLocation(PIx) == DataLocation_REdge)
	) {
		for (int k = 1; k < nRElements; k++) {
			double dCovUx =
				m_dStateREdge[WIx][k] * m_dColumnDerivRREdge[k][2];

			m_dXiDotREdge[k] =
				  m_dColumnContraMetricXiREdge[k][0]
					* m_dStateREdge[UIx][k]
				+ m_dColumnContraMetricXiREdge[k][1]
					* m_dStateREdge[VIx][k]
				+ m_dColumnContraMetricXiREdge[k][2]
					* dCovUx;
		}

		m_dXiDotREdge[0] = 0.0;
		m_dXiDotREdge[nRElements] = 0.0;
	}

#if !defined(EXPLICIT_VERTICAL_VELOCITY_ADVECTION) \
 && !defined(VERTICAL_VELOCITY_ADVECTION_CLARK)
	if (pGrid->GetVarLocation(WIx) == DataLocation_Node) {
		pGrid->DifferentiateNodeToNode(
			m_dStateNode[WIx],
			m_dDiffWNode);
	} else {
		pGrid->DifferentiateREdgeToREdge(
			m_dStateREdge[WIx],
			m_dDiffWREdge);
	}
#endif

#if defined(VERTICAL_HYPERVISCOSITY) || defined(UNIFORM_DIFFUSION)
	// Calculate second derivatives of other state variables
	if (m_nHypervisOrder > 0) {

		// Do not upwind horizontal velocity here
		for (int c = 2; c < 5; c++) {

			// Only upwind select variables
			if (!m_fUpwind[c]) {
				continue;
			}

			// High-order derivatives on interfaces
			if (pGrid->GetVarLocation(c) == DataLocation_REdge) {
				pGrid->DiffDiffREdgeToREdge(
					m_dStateREdge[c],
					m_dHyperDiffState[c]);

#ifdef UNIFORM_DIFFUSION
				// Uniform diffusion needs second derivatives of the state
				for (int k = 0; k <= nRElements; k++) {
					m_dDiffDiffState[c][k] = m_dHyperDiffState[c][k];
				}
#endif

				for (int h = 2; h < m_nHypervisOrder; h += 2) {
					memcpy(
						m_dStateAux,
						m_dHyperDiffState[c],
						(nRElements+1) * sizeof(double));

					pGrid->DiffDiffREdgeToREdge(
						m_dStateAux,
						m_dHyperDiffState[c]
					);
				}

			// High-order derivatives on levels
			} else {
				pGrid->DiffDiffNodeToNode(
					m_dStateNode[c],
					m_dHyperDiffState[c]);

#ifdef UNIFORM_DIFFUSION
				// Uniform diffusion needs second derivatives of the state
				for (int k = 0; k < nRElements; k++) {
					m_dDiffDiffState[c][k] = m_dHyperDiffState[c][k];
				}
#endif

				for (int h = 2; h < m_nHypervisOrder; h += 2) {
					memcpy(
						m_dStateAux,
						m_dHyperDiffState[c],
						nRElements * sizeof(double));

					pGrid->DiffDiffNodeToNode(
						m_dStateAux,
						m_dHyperDiffState[c]
					);
				}
			}
		}
	}
#endif
}

///////////////////////////////////////////////////////////////////////////////

void VerticalDynamicsFEM::BuildF(
	const double * dX,
	double * dF
) {
	// Indices of EquationSet variables
	const int UIx = 0;
	const int VIx = 1;
	const int PIx = 2;
	const int WIx = 3;
	const int RIx = 4;

	// Finite element grid
	const GridGLL * pGrid = dynamic_cast<const GridGLL *>(m_model.GetGrid());

	// Physical constants
	const PhysicalConstants & phys = m_model.GetPhysicalConstants();

	// Number of radial elements
	const int nRElements = pGrid->GetRElements();

	// Under this configuration, set fluxes at boundaries to zero
	bool fZeroBoundaries =
		(pGrid->GetVarLocation(WIx) == DataLocation_REdge);

	// Mass flux on levels
	bool fMassFluxOnLevels = false;
	if (m_fForceMassFluxOnLevels) {
		fMassFluxOnLevels = true;
	} else if (pGrid->GetVarLocation(WIx) == DataLocation_Node) {
		fMassFluxOnLevels = true;
	}

	// Zero F
	memset(dF, 0, m_nColumnStateSize * sizeof(double));

	// Mass flux on model interfaces
	if (!fMassFluxOnLevels) {
		for (int k = 1; k < nRElements; k++) {
			m_dMassFluxREdge[k] =
				m_dColumnJacobianREdge[k]
				* m_dStateREdge[RIx][k]
				* m_dXiDotREdge[k];
		}

		pGrid->DifferentiateREdgeToNode(
			m_dMassFluxREdge,
			m_dDiffMassFluxNode);

	// Mass flux on model levels
	} else {
		for (int k = 0; k < nRElements; k++) {
			m_dMassFluxNode[k] =
				m_dColumnJacobianNode[k]
				* m_dStateNode[RIx][k]
				* m_dXiDotNode[k];
		}

		pGrid->DifferentiateNodeToNode(
			m_dMassFluxNode,
			m_dDiffMassFluxNode,
			fZeroBoundaries);
	}

	// Change in density on model levels
	double dSum = 0.0;
	for (int k = 0; k < nRElements; k++) {
		dSum += m_dColumnElementArea[k] * m_dDiffMassFluxNode[k];

		dF[VecFIx(FRIx, k)] =
			m_dDiffMassFluxNode[k]
			* m_dColumnInvJacobianNode[k];
	}

#ifndef EXPLICIT_THERMO
#ifdef FORMULATION_PRESSURE
	// Pressure flux calculated on model interfaces
	if (!fMassFluxOnLevels) {
		for (int k = 1; k < nRElements; k++) {
			m_dPressureFluxREdge[k] =
				m_dColumnJacobianREdge[k]
				* phys.GetGamma()
				* m_dStateREdge[PIx][k]
				* m_dXiDotREdge[k];
		}

		pGrid->DifferentiateREdgeToNode(
			m_dPressureFluxREdge,
			m_dDiffPressureFluxNode);

	// Pressure flux calculated on model levels
	} else {
		for (int k = 0; k < nRElements; k++) {
			m_dPressureFluxNode[k] =
				m_dColumnJacobianNode[k]
				* phys.GetGamma()
				* m_dStateNode[PIx][k]
				* m_dXiDotNode[k];
		}

		pGrid->DifferentiateNodeToNode(
			m_dPressureFluxNode,
			m_dDiffPressureFluxNode,
			fZeroBoundaries);
	}

	// Change in pressure on model levels
	for (int k = 0; k < nRElements; k++) {
		dF[VecFIx(FPIx, k)] =
			- (phys.GetGamma() - 1.0)
			* m_dXiDotNode[k]
			* m_dDiffPNode[k];

		dF[VecFIx(FPIx, k)] +=
			m_dDiffPressureFluxNode[k]
			* m_dColumnInvJacobianNode[k];
	}
#endif
#if defined(FORMULATION_RHOTHETA_PI) || defined(FORMULATION_RHOTHETA_P)
	// RhoTheta flux on model interfaces
	if (!fMassFluxOnLevels) {
		for (int k = 1; k < nRElements; k++) {
			m_dPressureFluxREdge[k] =
				m_dColumnJacobianREdge[k]
				* m_dStateREdge[PIx][k]
				* m_dXiDotREdge[k];
		}

		pGrid->DifferentiateREdgeToNode(
			m_dPressureFluxREdge,
			m_dDiffPressureFluxNode);

	// RhoTheta flux on model levels
	} else {
		for (int k = 0; k < nRElements; k++) {
			m_dPressureFluxNode[k] =
				m_dColumnJacobianNode[k]
				* m_dStateNode[PIx][k]
				* m_dXiDotNode[k];
		}

		pGrid->DifferentiateNodeToNode(
			m_dPressureFluxNode,
			m_dDiffPressureFluxNode,
			fZeroBoundaries);
	}

	// Change in RhoTheta on model levels
	for (int k = 0; k < nRElements; k++) {
		dF[VecFIx(FPIx, k)] +=
			m_dDiffPressureFluxNode[k]
			* m_dColumnInvJacobianNode[k];
	}

#endif
#ifdef FORMULATION_THETA
	// Update theta on model levels
	if (pGrid->GetVarLocation(PIx) == DataLocation_Node) {

		// Change in Theta on model levels
		for (int k = 0; k < nRElements; k++) {
			dF[VecFIx(FPIx, k)] +=
				m_dXiDotNode[k] * m_dDiffThetaNode[k];
		}

	// Update theta on model interfaces
	} else {
		// Change in Theta on model interfaces
		for (int k = 0; k <= nRElements; k++) {
			dF[VecFIx(FPIx, k)] +=
				m_dXiDotREdge[k] * m_dDiffThetaREdge[k];
		}
	}
#endif
#ifdef FORMULATION_THETA_FLUX
	// Update theta on model levels
	if (pGrid->GetVarLocation(PIx) == DataLocation_Node) {
		// Pressure flux on model levels
		for (int k = 0; k < nRElements; k++) {
			m_dPressureFluxNode[k] =
				m_dColumnJacobianNode[k]
				* m_dStateNode[PIx][k]
				* m_dXiDotNode[k];
		}

		// Xidot derivatives on model levels
		pGrid->DifferentiateNodeToNode(
			m_dXiDotNode,
			m_dStateAuxDiff);

		// Theta flux derivatives on model levels
		pGrid->DifferentiateNodeToNode(
			m_dPressureFluxNode,
			m_dDiffPressureFluxNode);

		// Change in Theta on model levels
		for (int k = 0; k < nRElements; k++) {
			dF[VecFIx(FPIx, k)] +=
				m_dDiffPressureFluxNode[k]
			  * m_dColumnInvJacobianNode[k];

			dF[VecFIx(FPIx, k)] -=
				m_dStateNode[PIx][k] * m_dStateAuxDiff[k];
		}

	// Update theta on model interfaces
	} else {
		// Pressure flux on model levels
		for (int k = 0; k <= nRElements; k++) {
			m_dPressureFluxREdge[k] =
				m_dColumnJacobianREdge[k]
				* m_dStateREdge[PIx][k]
				* m_dXiDotREdge[k];

			m_dStateAux[k] =
				m_dColumnJacobianREdge[k]
				* m_dXiDotREdge[k];
		}

		// Xidot divergence on model interfaces
		pGrid->DifferentiateREdgeToREdge(
			m_dStateAux,
			m_dStateAuxDiff);

		// Theta flux derivatives on model interfaces
		pGrid->DifferentiateREdgeToREdge(
			m_dPressureFluxREdge,
			m_dDiffPressureFluxREdge);

		// Change in Theta on model interfaces
		for (int k = 0; k <= nRElements; k++) {
			dF[VecFIx(FPIx, k)] +=
				(m_dDiffPressureFluxREdge[k]
				- m_dStateREdge[PIx][k] * m_dStateAuxDiff[k])
				* m_dColumnInvJacobianREdge[k];
		}
	}
#endif
#endif

#if !defined(EXPLICIT_VERTICAL_VELOCITY_ADVECTION)
#if defined(VERTICAL_VELOCITY_ADVECTION_CLARK)
	// Kinetic energy on model levels
	for (int k = 0; k < nRElements; k++) {
		double dCovUa = m_dStateNode[UIx][k];
		double dCovUb = m_dStateNode[VIx][k];
		double dCovUx = m_dStateNode[WIx][k] * m_dColumnDerivRNode[k][2];

		double dConUa =
			  m_dColumnContraMetricA[k][0] * dCovUa
			+ m_dColumnContraMetricA[k][1] * dCovUb
			+ m_dColumnContraMetricA[k][2] * dCovUx;

		double dConUb =
			  m_dColumnContraMetricB[k][0] * dCovUa
			+ m_dColumnContraMetricB[k][1] * dCovUb
			+ m_dColumnContraMetricB[k][2] * dCovUx;

		double dConUx =
			  m_dColumnContraMetricXi[k][0] * dCovUa
			+ m_dColumnContraMetricXi[k][1] * dCovUb
			+ m_dColumnContraMetricXi[k][2] * dCovUx;

		// Specific kinetic energy
		m_dKineticEnergyNode[k] =
			  0.5 * (dConUa * dCovUa + dConUb * dCovUb + dConUx * dCovUx);
	}

	if (pGrid->GetVarLocation(WIx) == DataLocation_Node) {
		pGrid->DifferentiateNodeToNode(
			m_dKineticEnergyNode,
			m_dDiffKineticEnergyNode);
	} else {
		pGrid->DifferentiateNodeToREdge(
			m_dKineticEnergyNode,
			m_dDiffKineticEnergyREdge);
	}
#endif
#endif

	// Update equation for vertical velocity on levels
	if (pGrid->GetVarLocation(WIx) == DataLocation_Node) {

		for (int k = 1; k < nRElements-1; k++) {

			// Vertical pressure gradient force
#if defined(FORMULATION_PRESSURE) \
 || defined(FORMULATION_RHOTHETA_P)
			double dPressureGradientForce =
				m_dDiffPNode[k] / m_dStateNode[RIx][k];
#endif
#ifdef FORMULATION_RHOTHETA_PI
			double dPressureGradientForce =
				  m_dDiffPNode[k]
				* m_dStateNode[PIx][k]
				/ m_dStateNode[RIx][k];
#endif
#if defined(FORMULATION_THETA) || defined(FORMULATION_THETA_FLUX)
			double dPressureGradientForce =
				  m_dDiffPNode[k]
				* m_dStateNode[PIx][k];
#endif

			dF[VecFIx(FWIx, k)] =
				dPressureGradientForce / m_dColumnDerivRNode[k][2];

			dF[VecFIx(FWIx, k)] +=
				phys.GetG();

#if !defined(EXPLICIT_VERTICAL_VELOCITY_ADVECTION)
#if defined(VERTICAL_VELOCITY_ADVECTION_CLARK)
			// Vertical advection of vertical velocity
			double dCovUa = m_dStateNode[UIx][k];
			double dCovUb = m_dStateNode[VIx][k];
			double dCovUx = m_dStateNode[WIx][k] * m_dColumnDerivRNode[k][2];

			double dConUa =
				  m_dColumnContraMetricA[k][0] * dCovUa
				+ m_dColumnContraMetricA[k][1] * dCovUb
				+ m_dColumnContraMetricA[k][2] * dCovUx;

			double dConUb =
				  m_dColumnContraMetricB[k][0] * dCovUa
				+ m_dColumnContraMetricB[k][1] * dCovUb
				+ m_dColumnContraMetricB[k][2] * dCovUx;

			double dCurlTerm =
				- dConUa * m_dDiffUa[k]
				- dConUb * m_dDiffUb[k];

			dF[VecFIx(FWIx, k)] +=
				(m_dDiffKineticEnergyNode[k] + dCurlTerm)
					/ m_dColumnDerivRNode[k][2];

#else // VERTICAL VELOCITY ADVECTION (ADVECTIVE FORM)
			dF[VecFIx(FWIx, k)] +=
				m_dXiDotNode[k] * m_dDiffWNode[k];
#endif
#endif
		}

	// Update equation for vertical velocity on interfaces
	} else {
		for (int k = 1; k < nRElements; k++) {

			// Pressure gradient force
#if defined(FORMULATION_PRESSURE) \
 || defined(FORMULATION_RHOTHETA_P)
			double dPressureGradientForce =
				m_dDiffPREdge[k] / m_dStateREdge[RIx][k];
#endif
#ifdef FORMULATION_RHOTHETA_PI
			double dPressureGradientForce =
				  m_dDiffPREdge[k]
				* m_dStateREdge[PIx][k]
				/ m_dStateREdge[RIx][k];
#endif
#if defined(FORMULATION_THETA) || defined(FORMULATION_THETA_FLUX)
			double dPressureGradientForce =
				  m_dDiffPREdge[k]
				* m_dStateREdge[PIx][k];
#endif

			dF[VecFIx(FWIx, k)] =
				dPressureGradientForce / m_dColumnDerivRREdge[k][2];

			dF[VecFIx(FWIx, k)] +=
				phys.GetG();

#if !defined(EXPLICIT_VERTICAL_VELOCITY_ADVECTION)
#if defined(VERTICAL_VELOCITY_ADVECTION_CLARK)
			// Vertical advection of vertical velocity
			double dCovUa = m_dStateREdge[UIx][k];
			double dCovUb = m_dStateREdge[VIx][k];
			double dCovUx = m_dStateREdge[WIx][k] * m_dColumnDerivRREdge[k][2];

			double dConUa =
				  m_dColumnContraMetricAREdge[k][0] * dCovUa
				+ m_dColumnContraMetricAREdge[k][1] * dCovUb
				+ m_dColumnContraMetricAREdge[k][2] * dCovUx;

			double dConUb =
				  m_dColumnContraMetricBREdge[k][0] * dCovUa
				+ m_dColumnContraMetricBREdge[k][1] * dCovUb
				+ m_dColumnContraMetricBREdge[k][2] * dCovUx;

			double dCurlTerm =
				- dConUa * m_dDiffUa[k]
				- dConUb * m_dDiffUb[k];

			dF[VecFIx(FWIx, k)] +=
				(m_dDiffKineticEnergyREdge[k] + dCurlTerm)
					/ m_dColumnDerivRREdge[k][2];

#else // VERTICAL VELOCITY ADVECTION (ADVECTIVE FORM)
			dF[VecFIx(FWIx, k)] +=
				m_dXiDotREdge[k] * m_dDiffWREdge[k];

#endif
#endif
		}
	}

#ifdef UNIFORM_DIFFUSION
	// Apply uniform diffusion to theta and vertical velocity
	double dZtop = pGrid->GetZtop();

	double dUniformDiffusionCoeff = 75.0 / (dZtop * dZtop);

	// NOTE: Do not apply diffusion to density
	for (int c = 2; c < 4; c++) {

		// Uniform diffusion on interfaces
		if (pGrid->GetVarLocation(c) == DataLocation_REdge) {
			for (int k = 0; k <= nRElements; k++) {
				dF[VecFIx(FIxFromCIx(c), k)] -=
					dUniformDiffusionCoeff
					* m_dDiffDiffState[c][k];
			}

		// Uniform diffusion on levels
		} else {
			for (int k = 0; k < nRElements; k++) {
				dF[VecFIx(FIxFromCIx(c), k)] -=
					dUniformDiffusionCoeff
					* m_dDiffDiffState[c][k];
			}
		}
	}
#endif

#ifdef VERTICAL_UPWINDING
	{
		if (m_fUpwind[3]) {
			_EXCEPTIONT("Not implemented: Vertical upwinding of vertical velocity");
		}

		// Get penalty operator
		const LinearColumnDiscPenaltyFEM & opPenalty =
			pGrid->GetOpPenaltyNodeToNode();

		int nFiniteElements = nRElements / m_nVerticalOrder;

		// Calculate weights
		for (int a = 0; a < nFiniteElements - 1; a++) {
			int k = (a+1) * m_nVerticalOrder;
			m_dUpwindWeights[a] =
				0.5 * fabs(m_dXiDotREdge[k]);
		}

		// Loop through all variables
		for (int c = 2; c < 5; c++) {

			// Only upwind select variables
			if (!m_fUpwind[c]) {
				continue;
			}

#ifdef EXPLICIT_THERMO
			// Don't upwind thermodynamic variable if explicit
			if (c == PIx) {
				continue;
			}
#endif

			// Upwinding on interfaces (not implemented)
			if (pGrid->GetVarLocation(c) == DataLocation_REdge) {
				_EXCEPTIONT("Upwinding on interfaces not implemented");

			// Upwinding on levels
			} else {

				// Apply upwinding
				m_dStateAux.Zero();
				opPenalty.Apply(
					&(m_dUpwindWeights[0]),
					&(m_dStateNode[c][0]),
					&(m_dStateAux[0]),
					1,
					1);

				for (int k = 0; k < nRElements; k++) {
					dF[VecFIx(FIxFromCIx(c), k)] -= m_dStateAux[k];
				}
			}
		}
	}
#endif

#if defined(VERTICAL_HYPERVISCOSITY) || defined(UNIFORM_DIFFUSION)
	// Apply flow-dependent hyperviscosity
	if (m_nHypervisOrder > 0) {
		for (int c = 2; c < 5; c++) {

			// Only upwind select variables
			if (!m_fUpwind[c]) {
				continue;
			}

			// Flow-dependent hyperviscosity on interfaces
			if (pGrid->GetVarLocation(c) == DataLocation_REdge) {

				for (int k = 0; k <= nRElements; k++) {
					dF[VecFIx(FIxFromCIx(c), k)] -=
						m_dHypervisCoeff
						* fabs(m_dXiDotREdge[k])
						* m_dHyperDiffState[c][k];
				}

			// Flow-dependent hyperviscosity on levels
			} else {

				for (int k = 0; k < nRElements; k++) {
					dF[VecFIx(FIxFromCIx(c), k)] -=
						m_dHypervisCoeff
						* fabs(m_dXiDotNode[k])
						* m_dHyperDiffState[c][k];
				}
			}
		}
	}
#endif

#ifdef DEBUG
	if (dF[VecFIx(FWIx, 0)] != 0.0) {
		_EXCEPTIONT("No updates to W at bottom boundary allowed");
	}
	if (pGrid->GetVarLocation(WIx) == DataLocation_REdge) {
		if (dF[VecFIx(FWIx, nRElements)] != 0.0) {
			_EXCEPTIONT("No updates to W at top boundary allowed");
		}
	} else {
		if (dF[VecFIx(FWIx, nRElements-1)] != 0.0) {
			_EXCEPTIONT("No updates to W at top boundary allowed");
		}
	}
#endif

	// Construct the time-dependent component of the RHS
	double dInvDeltaT = 1.0 / m_dDeltaT;
	for (int i = 0; i < m_nColumnStateSize; i++) {
		dF[i] += (dX[i] - m_dColumnState[i]) * dInvDeltaT;
	}
/*
	// Apply boundary conditions to W
	if (pGrid->GetVarLocation(WIx) == DataLocation_REdge) {
		dF[VecFIx(FWIx, 0)] =
			  dContraMetricXiREdge[0][m_iA][m_iB][0] * m_dStateREdge[UIx][0]
			+ dContraMetricXiREdge[0][m_iA][m_iB][1] * m_dStateREdge[VIx][0]
			+ dContraMetricXiREdge[0][m_iA][m_iB][2]
				* dDerivRREdge[0][m_iA][m_iB][2] * m_dStateREdge[WIx][0];

		int k = nRElements;

		dF[VecFIx(FWIx, k)] =
			  dContraMetricXiREdge[k][m_iA][m_iB][0] * m_dStateREdge[UIx][k]
			+ dContraMetricXiREdge[k][m_iA][m_iB][1] * m_dStateREdge[VIx][k]
			+ dContraMetricXiREdge[k][m_iA][m_iB][2]
				* dDerivRREdge[k][m_iA][m_iB][2] * m_dStateREdge[WIx][k];

	} else
*/
/*
	if (
		pGrid->GetVerticalStaggering() ==
			Grid::VerticalStaggering_Interfaces
	) {
		dF[VecFIx(FWIx, 0)] =
			  m_dColumnContraMetricXi[0][0] * m_dStateNode[UIx][0]
			+ m_dColumnContraMetricXi[0][1] * m_dStateNode[VIx][0]
			+ m_dColumnContraMetricXi[0][2]
				* m_dColumnDerivRNode[0][2] * m_dStateNode[WIx][0];

		int k = nRElements-1;

		dF[VecFIx(FWIx, k)] =
			  m_dColumnContraMetricXi[k][0] * m_dStateNode[UIx][k]
			+ m_dColumnContraMetricXi[k][1] * m_dStateNode[VIx][k]
			+ m_dColumnContraMetricXi[k][2]
				* m_dColumnDerivRNode[k][2] * m_dStateNode[WIx][k];

	} else {
		_EXCEPTIONT("UNIMPLEMENTED");
	}
*/
}

///////////////////////////////////////////////////////////////////////////////

void VerticalDynamicsFEM::BuildJacobianF(
	const double * dX,
	double * dDG
) {
	// Indices of EquationSet variables
	const int UIx = 0;
	const int VIx = 1;
	const int PIx = 2;
	const int WIx = 3;
	const int RIx = 4;

	// Finite element grid
	const GridGLL * pGrid = dynamic_cast<const GridGLL *>(m_model.GetGrid());

	// Get the column interpolation and differentiation coefficients
	const LinearColumnInterpFEM & opInterpNodeToREdge =
		pGrid->GetOpInterpNodeToREdge();
	const LinearColumnInterpFEM & opInterpREdgeToNode =
		pGrid->GetOpInterpREdgeToNode();
	const LinearColumnDiffFEM & opDiffNodeToNode =
		pGrid->GetOpDiffNodeToNode();
	const LinearColumnDiffFEM & opDiffNodeToREdge =
		pGrid->GetOpDiffNodeToREdge();
	const LinearColumnDiffFEM & opDiffREdgeToNode =
		pGrid->GetOpDiffREdgeToNode();
	const LinearColumnDiffFEM & opDiffREdgeToREdge =
		pGrid->GetOpDiffREdgeToREdge();

	const DataArray2D<double> & dInterpNodeToREdge =
		opInterpNodeToREdge.GetCoeffs();
	const DataArray2D<double> & dInterpREdgeToNode =
		opInterpREdgeToNode.GetCoeffs();
	const DataArray2D<double> & dDiffNodeToNode =
		opDiffNodeToNode.GetCoeffs();
	const DataArray2D<double> & dDiffNodeToREdge =
		opDiffNodeToREdge.GetCoeffs();
	const DataArray2D<double> & dDiffREdgeToNode =
		opDiffREdgeToNode.GetCoeffs();
	const DataArray2D<double> & dDiffREdgeToREdge =
		opDiffREdgeToREdge.GetCoeffs();

	const DataArray1D<int> & iInterpNodeToREdgeBegin =
		opInterpNodeToREdge.GetIxBegin();
	const DataArray1D<int> & iInterpREdgeToNodeBegin =
		opInterpREdgeToNode.GetIxBegin();
	const DataArray1D<int> & iDiffNodeToNodeBegin =
		opDiffNodeToNode.GetIxBegin();
	const DataArray1D<int> & iDiffNodeToREdgeBegin =
		opDiffNodeToREdge.GetIxBegin();
	const DataArray1D<int> & iDiffREdgeToNodeBegin =
		opDiffREdgeToNode.GetIxBegin();
	const DataArray1D<int> & iDiffREdgeToREdgeBegin =
		opDiffREdgeToREdge.GetIxBegin();

	const DataArray1D<int> & iInterpNodeToREdgeEnd =
		opInterpNodeToREdge.GetIxEnd();
	const DataArray1D<int> & iInterpREdgeToNodeEnd =
		opInterpREdgeToNode.GetIxEnd();
	const DataArray1D<int> & iDiffNodeToNodeEnd =
		opDiffNodeToNode.GetIxEnd();
	const DataArray1D<int> & iDiffNodeToREdgeEnd =
		opDiffNodeToREdge.GetIxEnd();
	const DataArray1D<int> & iDiffREdgeToNodeEnd =
		opDiffREdgeToNode.GetIxEnd();
	const DataArray1D<int> & iDiffREdgeToREdgeEnd =
		opDiffREdgeToREdge.GetIxEnd();

#ifdef VERTICAL_UPWINDING
	// Get the column penalization coefficients
	const LinearColumnDiscPenaltyFEM & opPenalty =
		pGrid->GetOpPenaltyNodeToNode();
	const LinearColumnOperator & opPenaltyLeft =
		opPenalty.GetLeftOp();
	const LinearColumnOperator & opPenaltyRight =
		opPenalty.GetRightOp();

	const DataArray2D<double> & dPenaltyLeft =
		opPenaltyLeft.GetCoeffs();
	const DataArray2D<double> & dPenaltyRight =
		opPenaltyRight.GetCoeffs();

	const DataArray1D<int> & iPenaltyLeftBegin =
		opPenaltyLeft.GetIxBegin();
	const DataArray1D<int> & iPenaltyRightBegin =
		opPenaltyRight.GetIxBegin();

	const DataArray1D<int> & iPenaltyLeftEnd =
		opPenaltyLeft.GetIxEnd();
	const DataArray1D<int> & iPenaltyRightEnd =
		opPenaltyRight.GetIxEnd();
#endif

	// Physical constants
	const PhysicalConstants & phys = m_model.GetPhysicalConstants();

	// Number of radial elements
	const int nRElements = pGrid->GetRElements();

	// Number of radial finite elements
	const int nFiniteElements = nRElements / m_nVerticalOrder;

	// Zero DG
	memset(dDG, 0,
		m_nColumnStateSize * m_nColumnStateSize * sizeof(double));

	// Mass flux on levels
	bool fMassFluxOnLevels = false;
	if (m_fForceMassFluxOnLevels) {
		fMassFluxOnLevels = true;
	} else if (pGrid->GetVarLocation(WIx) == DataLocation_Node) {
		fMassFluxOnLevels = true;
	}

	if (fMassFluxOnLevels) {
		_EXCEPTIONT("Mass flux on levels -- not implemented");
	}

#ifdef VERTICAL_HYPERVISCOSITY
	// Check upwinding
	for (int c = 2; c < 5; c++) {
		if (m_fUpwind[c]) {
			_EXCEPTIONT("Hyperviscosity not implemented in BuildJacobian");
		}
	}
#endif

//////////////////////////////////////////////
// Prognostic thermodynamic variable pressure
#ifdef FORMULATION_PRESSURE

	// Vertical velocity on interfaces (CPH or LOR staggering)
	if (pGrid->GetVarLocation(WIx) == DataLocation_REdge) {
		_EXCEPTIONT("Not implemented");

	// Vertical velocity on nodes (LEV or INT staggering)
	} else {

#ifndef EXPLICIT_THERMO
		// dP_k/dP_n
		for (int k = 0; k < nRElements; k++) {

			int n = iDiffNodeToNodeBegin[k];
			for (; n < iDiffNodeToNodeEnd[k]; n++) {

				// Pressure flux
				dDG[MatFIx(FPIx, n, FPIx, k)] +=
					phys.GetGamma()
					* m_dColumnJacobianNode[n]
					* m_dColumnInvJacobianNode[k]
					* dDiffNodeToNode[k][n]
					* m_dXiDotNode[n];

				// Correction terms
				dDG[MatFIx(FPIx, n, FPIx, k)] +=
					- (phys.GetGamma() - 1.0)
					* m_dXiDotNode[k]
					* dDiffNodeToNode[k][n];
			}
		}

		// dP_k/dW_n
		for (int k = 0; k < nRElements; k++) {

			int n = iDiffNodeToNodeBegin[k];
			for (; n < iDiffNodeToNodeEnd[k]; n++) {
				if (pGrid->GetVerticalStaggering() ==
				    Grid::VerticalStaggering_Interfaces
				) {
					if ((n == 0) || (n == nRElements-1)) {
						continue;
					}
				}

				// Pressure flux
				dDG[MatFIx(FWIx, n, FPIx, k)] +=
					dDiffNodeToNode[k][n]
					* m_dColumnInvJacobianNode[k]
					* m_dColumnJacobianNode[n]
					* phys.GetGamma()
					* m_dStateNode[PIx][n]
					* m_dColumnContraMetricXi[n][2]
					* m_dColumnDerivRNode[n][2];
			}

			// Correction terms
			if (pGrid->GetVerticalStaggering() ==
			    Grid::VerticalStaggering_Interfaces
			) {
				if ((k == 0) || (k == nRElements-1)) {
					continue;
				}
			}

			dDG[MatFIx(FWIx, k, FPIx, k)] +=
				- (phys.GetGamma() - 1.0)
				* m_dColumnContraMetricXi[k][2]
				* m_dColumnDerivRNode[k][2]
				* m_dDiffPNode[k];
		}
#endif

		// Account for interfaces
		int kBegin = 0;
		int kEnd = nRElements;

		if (pGrid->GetVerticalStaggering() ==
		    Grid::VerticalStaggering_Interfaces
		) {
			kBegin = 1;
			kEnd = nRElements-1;
		}

		// dW_k/dP_n and dW_k/dR_k
		for (int k = kBegin; k < kEnd; k++) {

			int n = iDiffNodeToNodeBegin[k];
			for (; n < iDiffNodeToNodeEnd[k]; n++) {
				dDG[MatFIx(FPIx, n, FWIx, k)] +=
					dDiffNodeToNode[k][n]
					/ m_dStateNode[RIx][k]
					/ m_dColumnDerivRNode[k][2];
			}

			dDG[MatFIx(FRIx, k, FWIx, k)] +=
				- m_dDiffPNode[k]
				/ m_dColumnDerivRNode[k][2]
				/ (m_dStateNode[RIx][k] * m_dStateNode[RIx][k]);
		}
	}

#endif
#ifdef FORMULATION_RHOTHETA_PI
	_EXCEPTIONT("Not implemented");
#endif
#ifdef FORMULATION_RHOTHETA_P
	_EXCEPTIONT("Not implemented");
#endif

//////////////////////////////////////////////
// Prognostic thermodynamic variable theta
#ifdef FORMULATION_THETA

	// Lorenz staggering
	if (pGrid->GetVarLocation(PIx) == DataLocation_Node) {
		if (pGrid->GetVarLocation(WIx) != DataLocation_REdge) {
			_EXCEPTIONT("Not implemented");
		}

#if !defined(EXPLICIT_THERMO)
		// dT_k/dW_l
		for (int k = 0; k < nRElements; k++) {
			int l = iInterpREdgeToNodeBegin[k];
			for (; l < iInterpREdgeToNodeEnd[k]; l++) {
				dDG[MatFIx(FWIx, l, FPIx, k)] +=
					m_dDiffThetaNode[k]
					* dInterpREdgeToNode[k][l]
					* m_dColumnContraMetricXi[k][2]
					* m_dColumnDerivRNode[k][2];
			}
		}

		// dT_k/dT_l
		for (int k = 0; k < nRElements; k++) {
			int l = iDiffNodeToNodeBegin[k];
			for (; l < iDiffNodeToNodeEnd[k]; l++) {
				dDG[MatFIx(FPIx, l, FPIx, k)] +=
					dDiffNodeToNode[k][l]
					* m_dXiDotNode[k];
			}
		}
#endif

		// dW_k/dT_m and dW_k/dR_m
		for (int k = 1; k < nRElements; k++) {

			double dRHSWCoeff = 
				1.0 / m_dColumnDerivRNode[k][2]
				* m_dStateREdge[PIx][k]
				* phys.GetR()
				/ phys.GetCv();

			int m = iDiffNodeToREdgeBegin[k];
			for (; m < iDiffNodeToREdgeEnd[k]; m++) {

				dDG[MatFIx(FPIx, m, FWIx, k)] +=
					dRHSWCoeff 
					* dDiffNodeToREdge[k][m]
					* m_dExnerNode[m]
					/ m_dStateNode[PIx][m];

				dDG[MatFIx(FRIx, m, FWIx, k)] +=
					dRHSWCoeff
					* dDiffNodeToREdge[k][m]
					* m_dExnerNode[m]
					/ m_dStateNode[RIx][m];
			}
		}

		// dW_k/dT_k (first theta in RHS)
		for (int k = 1; k < nRElements; k++) {
			int l = iInterpNodeToREdgeBegin[k];
			for (; l < iInterpNodeToREdgeEnd[k]; l++) {
				dDG[MatFIx(FPIx, l, FWIx, k)] +=
					 1.0 / m_dColumnDerivRREdge[k][2]
					 * dInterpNodeToREdge[k][l]
					 * m_dDiffPREdge[k];
			}
		}

	// Charney-Phillips staggering
	} else {
		if (pGrid->GetVarLocation(WIx) != DataLocation_REdge) {
			_EXCEPTIONT("Not implemented");
		}

#ifndef EXPLICIT_THERMO
		// dT_k/dW_k
		for (int k = 1; k < nRElements; k++) {
			dDG[MatFIx(FWIx, k, FPIx, k)] +=
				m_dDiffThetaREdge[k]
				* m_dColumnContraMetricXiREdge[k][2]
				* m_dColumnDerivRREdge[k][2];
		}

		// dT_k/dT_l
		for (int k = 0; k <= nRElements; k++) {
			int l = iDiffREdgeToREdgeBegin[k];
			for (; l < iDiffREdgeToREdgeEnd[k]; l++) {
				dDG[MatFIx(FPIx, l, FPIx, k)] +=
					dDiffREdgeToREdge[k][l]
					* m_dXiDotREdge[k];
			}
		}
#endif

		// dW_k/dT_l and dW_k/dR_m
		for (int k = 1; k < nRElements; k++) {

			double dRHSWCoeff = 
				1.0 / m_dColumnDerivRNode[k][2]
				* m_dStateREdge[PIx][k]
				* phys.GetR()
				/ phys.GetCv();

			int m = iDiffNodeToREdgeBegin[k];
			for (; m < iDiffNodeToREdgeEnd[k]; m++) {

				double dTEntry =
					dRHSWCoeff 
					* dDiffNodeToREdge[k][m]
					* m_dExnerNode[m]
					/ m_dStateNode[PIx][m];

				int l = iInterpREdgeToNodeBegin[m];
				for (; l < iInterpREdgeToNodeEnd[m]; l++) {
					dDG[MatFIx(FPIx, l, FWIx, k)] +=
						dTEntry * dInterpREdgeToNode[m][l];
				}

				dDG[MatFIx(FRIx, m, FWIx, k)] +=
					dRHSWCoeff
					* dDiffNodeToREdge[k][m]
					* m_dExnerNode[m]
					/ m_dStateNode[RIx][m];
			}
		}

		// dW_k/dT_k (first theta in RHS)
		for (int k = 1; k < nRElements; k++) {
			dDG[MatFIx(FPIx, k, FWIx, k)] +=
				 1.0 / m_dColumnDerivRREdge[k][2]
				 * m_dDiffPREdge[k];
		}
	}
#endif

	// Vertical velocity on interfaces (CPH or LOR staggering)
	if (pGrid->GetVarLocation(WIx) == DataLocation_REdge) {
		for (int k = 0; k < nRElements; k++) {

			// dRho_k/dW_l and dRho_k/dRho_n
			int m = iDiffREdgeToNodeBegin[k];
			for (; m < iDiffREdgeToNodeEnd[k]; m++) {

				if ((m != 0) && (m != nRElements)) {
					dDG[MatFIx(FWIx, m, FRIx, k)] +=
						dDiffREdgeToNode[k][m]
						* m_dColumnInvJacobianREdge[m]
						* m_dColumnJacobianNode[k]
						* m_dStateREdge[RIx][m]
						* m_dColumnContraMetricXiREdge[m][2]
						* m_dColumnDerivRREdge[m][2];
				}

				int n = iInterpNodeToREdgeBegin[m];
				for (; n < iInterpNodeToREdgeEnd[m]; n++) {

					dDG[MatFIx(FRIx, n, FRIx, k)] +=
						dDiffREdgeToNode[k][m]
						* m_dColumnJacobianREdge[m]
						* m_dColumnInvJacobianNode[k]
						* dInterpNodeToREdge[m][n]
						* m_dXiDotREdge[m];
				}
			}
		}

#if !defined(EXPLICIT_VERTICAL_VELOCITY_ADVECTION)
#if defined(VERTICAL_VELOCITY_ADVECTION_CLARK)
		// dW_k/dW_m
		for (int k = 1; k < nRElements; k++) {
			int l = iDiffNodeToREdgeBegin[k];
			for (; l < iDiffNodeToREdgeEnd[k]; l++) {

				int m = iInterpREdgeToNodeBegin[l];
				for (; m < iInterpREdgeToNodeEnd[l]; m++) {
					dDG[MatFIx(FWIx, m, FWIx, k)] +=
						dInterpREdgeToNode[l][m]
						* dDiffNodeToREdge[k][l]
						/ m_dColumnDerivRREdge[k][2]
						* m_dColumnDerivRNode[l][2]
						* m_dXiDotNode[l];
				}
			}
		}
#else

		// dW_k/dW_m
		for (int k = 1; k < nRElements; k++) {
			int m = iDiffREdgeToREdgeBegin[k];
			for (; m < iDiffREdgeToREdgeEnd[k]; m++) {
				dDG[MatFIx(FWIx, m, FWIx, k)] +=
					m_dXiDotREdge[k]
					* dDiffREdgeToREdge[k][m];
			}

			dDG[MatFIx(FWIx, k, FWIx, k)] +=
				m_dDiffWREdge[k]
				* m_dColumnContraMetricXiREdge[k][2]
				* m_dColumnDerivRREdge[k][2];
		}

#endif
#endif

	// Vertical velocity on nodes (LEV or INT staggering)
	} else {

		for (int k = 0; k < nRElements; k++) {

			int n = iDiffNodeToNodeBegin[k];
			for (; n < iDiffNodeToNodeEnd[k]; n++) {

				// dRho_k/dRho_n
				dDG[MatFIx(FRIx, n, FRIx, k)] +=
					dDiffNodeToNode[k][n]
					* m_dColumnJacobianNode[n]
					* m_dColumnInvJacobianNode[k]
					* m_dXiDotNode[n];

				// Boundary conditions
				if (pGrid->GetVerticalStaggering() ==
				    Grid::VerticalStaggering_Interfaces
				) {
					if ((n == 0) || (n == nRElements-1)) {
						continue;
					}
				}

				// dRho_k/dW_n
				dDG[MatFIx(FWIx, n, FRIx, k)] +=
					dDiffNodeToNode[k][n]
					* m_dColumnJacobianNode[n]
					* m_dColumnInvJacobianNode[k]
					* m_dStateNode[RIx][n]
					* m_dColumnDerivRNode[n][2]
					* m_dColumnContraMetricXi[n][2];
			}

			// Boundary conditions
			if (pGrid->GetVerticalStaggering() ==
			    Grid::VerticalStaggering_Interfaces
			) {
				if ((k == 0) || (k == nRElements-1)) {
					continue;
				}
			}

#if !defined(EXPLICIT_VERTICAL_VELOCITY_ADVECTION)
#if defined(VERTICAL_VELOCITY_ADVECTION_CLARK)
			// dW_k/dW_n
			n = iDiffNodeToNodeBegin[k];
			for (; n < iDiffNodeToNodeEnd[k]; n++) {
				dDG[MatFIx(FWIx, n, FWIx, k)] +=
					dDiffNodeToNode[k][n]
					/ m_dColumnDerivRNode[k][2]
					* m_dColumnDerivRNode[n][2]
					* m_dXiDotNode[n];

			}
#else
		// dW_k/dW_m
		for (int k = 0; k < nRElements; k++) {
			int m = iDiffNodeToNodeBegin[k];
			for (; m < iDiffNodeToNodeEnd[k]; m++) {
				dDG[MatFIx(FWIx, m, FWIx, k)] +=
					m_dXiDotNode[k]
					* dDiffNodeToNode[k][m];
			}

			dDG[MatFIx(FWIx, k, FWIx, k)] +=
				m_dDiffWNode[k]
				* m_dColumnContraMetricXi[k][2];
		}
#endif
#endif

		}
	}

#ifdef VERTICAL_UPWINDING
	// Vertical upwinding
	for (int c = 2; c < 5; c++) {

#ifdef EXPLICIT_THERMO
		// Don't upwind thermodynamic variable if explicit
		if (c == PIx) {
			continue;
		}
#endif

		// Check upwinding
		if (!m_fUpwind[c]) {
			continue;
		}
		if (pGrid->GetVarLocation(c) != DataLocation_Node) {
			_EXCEPTIONT("Upwinding on interfaces not implemented");
		}
		if (pGrid->GetVarLocation(WIx) != DataLocation_REdge) {
			_EXCEPTIONT("Upwinding DIRECTSOLVE requires W on interfaces");
		}

		for (int a = 1; a < nFiniteElements; a++) {
			double dWeight = 0.5 * fabs(m_dXiDotREdge[a * m_nVerticalOrder]);
			double dSignWeight;
			if (m_dXiDotREdge[a * m_nVerticalOrder] > 0.0) {
				dSignWeight = 0.5;
			} else if (m_dXiDotREdge[a * m_nVerticalOrder] < 0.0) {
				dSignWeight = -0.5;
			} else {
				dSignWeight = 0.0;
			}

			int kLeftBegin = (a-1) * m_nVerticalOrder;
			int kLeftEnd = a * m_nVerticalOrder;

			int kRightBegin = a * m_nVerticalOrder;
			int kRightEnd = (a+1) * m_nVerticalOrder;

#pragma message "Need to account for horizontal flow contribution to xi_dot"
			// dC_k/dW_a (left operator)
			for (int k = kLeftBegin; k < kLeftEnd; k++) {
			for (int n = iPenaltyLeftBegin[k]; n < iPenaltyLeftEnd[k]; n++) {
				//printf("%1.15e %1.15e %1.15e\n",
				//	dSignWeight, dPenaltyLeft[k][n], m_dStateNode[c][n]);
				dDG[MatFIx(FWIx, kLeftEnd, FIxFromCIx(c), k)] -=
					dSignWeight
					/ m_dColumnDerivRREdge[kLeftEnd][2]
					* dPenaltyLeft[k][n]
					* m_dStateNode[c][n];
			}
			}

			// dC_k/dW_a (right operator)
			for (int k = kRightBegin; k < kRightEnd; k++) {
			for (int n = iPenaltyRightBegin[k]; n < iPenaltyRightEnd[k]; n++) {
				dDG[MatFIx(FWIx, kRightBegin, FIxFromCIx(c), k)] -=
					dSignWeight
					/ m_dColumnDerivRREdge[kRightBegin][2]
					* dPenaltyRight[k][n]
					* m_dStateNode[c][n];
			}
			}

			// dC_k/dC_n (left operator)
			for (int k = kLeftBegin; k < kLeftEnd; k++) {
			for (int n = iPenaltyLeftBegin[k]; n < iPenaltyLeftEnd[k]; n++) {
				dDG[MatFIx(FIxFromCIx(c), n, FIxFromCIx(c), k)] -=
					dWeight * dPenaltyLeft[k][n];
			}
			}

			// dC_k/dC_n (right operator)
			for (int k = kRightBegin; k < kRightEnd; k++) {
			for (int n = iPenaltyRightBegin[k]; n < iPenaltyRightEnd[k]; n++) {
				dDG[MatFIx(FIxFromCIx(c), n, FIxFromCIx(c), k)] -=
					dWeight * dPenaltyRight[k][n];
			}
			}
		}
	}
#endif

	// Add the identity components
	for (int k = 0; k <= nRElements; k++) {
		dDG[MatFIx(FPIx, k, FPIx, k)] += 1.0 / m_dDeltaT;
		dDG[MatFIx(FWIx, k, FWIx, k)] += 1.0 / m_dDeltaT;
		dDG[MatFIx(FRIx, k, FRIx, k)] += 1.0 / m_dDeltaT;
	}

	// Boundary conditions for W
	if (pGrid->GetVerticalStaggering() ==
	    Grid::VerticalStaggering_Interfaces
	) {
		dDG[MatFIx(FWIx, 0, FWIx, 0)] =
			m_dColumnDerivRNode[0][2]
			* m_dColumnContraMetricXi[0][2];
		dDG[MatFIx(FWIx, nRElements-1, FWIx, nRElements-1)] =
			m_dColumnDerivRNode[nRElements-1][2]
			* m_dColumnContraMetricXi[nRElements-1][2];
	}
}

///////////////////////////////////////////////////////////////////////////////

void VerticalDynamicsFEM::Evaluate(
	const double * dX,
	double * dF
) {
	// Prepare the column
	PrepareColumn(dX);

	// Evaluate the zero equations
	BuildF(dX, dF);
}

///////////////////////////////////////////////////////////////////////////////

void VerticalDynamicsFEM::UpdateColumnTracers(
	double dDeltaT,
	const DataArray4D<double> & dataInitialNode,
	const DataArray4D<double> & dataUpdateNode,
	const DataArray4D<double> & dataInitialREdge,
	const DataArray4D<double> & dataUpdateREdge,
	const DataArray4D<double> & dataInitialTracer,
	const DataArray4D<double> & dataUpdateTracer
) {
	// Indices of EquationSet variables
	const int UIx = 0;
	const int VIx = 1;
	const int PIx = 2;
	const int WIx = 3;
	const int RIx = 4;

	// Finite element grid
	const GridGLL * pGrid = dynamic_cast<const GridGLL *>(m_model.GetGrid());

	// Number of model levels
	const int nRElements = pGrid->GetRElements();

	// Number of tracer components
	const int nComponents = dataInitialTracer.GetSize(0);

	// If no tracer components, no update necessary
	if (nComponents == 0) {
		return;
	}

	// Get the column interpolation and differentiation coefficients
	const LinearColumnInterpFEM & opInterpNodeToREdge =
		pGrid->GetOpInterpNodeToREdge();
	const LinearColumnDiffFEM & opDiffNodeToNode =
		pGrid->GetOpDiffNodeToNode();
	const LinearColumnDiffFEM & opDiffREdgeToNode =
		pGrid->GetOpDiffREdgeToNode();

	const DataArray2D<double> & dInterpNodeToREdge =
		opInterpNodeToREdge.GetCoeffs();
	const DataArray2D<double> & dDiffNodeToNode =
		opDiffNodeToNode.GetCoeffs();
	const DataArray2D<double> & dDiffREdgeToNode =
		opDiffREdgeToNode.GetCoeffs();

	const DataArray1D<int> & iInterpNodeToREdgeBegin =
		opInterpNodeToREdge.GetIxBegin();
	const DataArray1D<int> & iDiffNodeToNodeBegin =
		opDiffNodeToNode.GetIxBegin();
	const DataArray1D<int> & iDiffREdgeToNodeBegin =
		opDiffREdgeToNode.GetIxBegin();

	const DataArray1D<int> & iInterpNodeToREdgeEnd =
		opInterpNodeToREdge.GetIxEnd();
	const DataArray1D<int> & iDiffNodeToNodeEnd =
		opDiffNodeToNode.GetIxEnd();
	const DataArray1D<int> & iDiffREdgeToNodeEnd =
		opDiffREdgeToNode.GetIxEnd();

	// Metric quantities
	const DataArray4D<double> & dContraMetricXi =
		m_pPatch->GetContraMetricXi();
	const DataArray4D<double> & dContraMetricXiREdge =
		m_pPatch->GetContraMetricXiREdge();
	const DataArray3D<double> & dElementArea =
		m_pPatch->GetElementArea();
	const DataArray3D<double> & dJacobianNode =
		m_pPatch->GetJacobian();
	const DataArray3D<double> & dJacobianREdge =
		m_pPatch->GetJacobianREdge();
	const DataArray4D<double> & dDerivRNode =
		m_pPatch->GetDerivRNode();
	const DataArray4D<double> & dDerivRREdge =
		m_pPatch->GetDerivRREdge();

	// Under this configuration, set fluxes at boundaries to zero
	bool fZeroBoundaries =
		(pGrid->GetVarLocation(WIx) == DataLocation_REdge);

	// Mass flux on levels
	bool fMassFluxOnLevels = false;
	if (m_fForceMassFluxOnLevels) {
		fMassFluxOnLevels = true;
	} else if (pGrid->GetVarLocation(WIx) == DataLocation_Node) {
		fMassFluxOnLevels = true;
	}

	// Zero the Jacobian
	m_matTracersLUDF.Zero();
	double * dTracersLUDF = &(m_matTracersLUDF[0][0]);

	// Only compute off-diagonal terms of the Jacobian if implicit advection
	// is being performed.
	if (!m_fFullyExplicit) {

		// Mass flux on model levels
		if (fMassFluxOnLevels) {

			// Calculate u^xi on model levels
			if (pGrid->GetVarLocation(WIx) == DataLocation_REdge) {
				for (int k = 0; k <= nRElements; k++) {
					m_dStateREdge[WIx][k] = m_dColumnState[VecFIx(FWIx, k)];
				}

				pGrid->InterpolateREdgeToNode(
					m_dStateREdge[WIx],
					m_dStateNode[WIx]);

			} else {
				for (int k = 0; k < nRElements; k++) {
					m_dStateNode[WIx][k] = m_dColumnState[VecFIx(FWIx, k)];
				}
			}

			for (int k = 1; k < nRElements-1; k++) {
				double dCovUx =
					m_dStateNode[WIx][k] * dDerivRNode[k][m_iA][m_iB][2];

				m_dXiDotNode[k] =
					  dContraMetricXi[k][m_iA][m_iB][0] * m_dStateNode[UIx][k]
					+ dContraMetricXi[k][m_iA][m_iB][1] * m_dStateNode[VIx][k]
					+ dContraMetricXi[k][m_iA][m_iB][2] * dCovUx;
			}

			m_dXiDotNode[0] = 0.0;
			m_dXiDotNode[nRElements-1] = 0.0;

			_EXCEPTION();
			// dRhoQ_k/dRhoQ_n
			for (int k = 0; k < nRElements; k++) {

				int n = iDiffNodeToNodeBegin[k];
				for (; n < iDiffNodeToNodeEnd[k]; n++) {

					dTracersLUDF[TracerMatFIx(n, k)] +=
						dDiffNodeToNode[k][n]
						* dJacobianNode[n][m_iA][m_iB]
						/ dJacobianNode[k][m_iA][m_iB]
						* m_dXiDotNode[n];
				}
			}

		// Mass flux on model interfaces
		} else {

			// Calculate u^xi on model interfaces
			if (pGrid->GetVarLocation(WIx) == DataLocation_REdge) {
				for (int k = 0; k <= nRElements; k++) {
					m_dStateREdge[WIx][k] = m_dColumnState[VecFIx(FWIx, k)];
				}

			} else {
				for (int k = 0; k < nRElements; k++) {
					m_dStateNode[WIx][k] = m_dColumnState[VecFIx(FWIx, k)];
				}

				pGrid->InterpolateNodeToREdge(
					m_dStateNode[WIx],
					m_dStateREdge[WIx]);
			}

			for (int k = 1; k < nRElements; k++) {
				double dCovUx =
					m_dStateREdge[WIx][k] * dDerivRREdge[k][m_iA][m_iB][2];

				m_dXiDotREdge[k] =
					  dContraMetricXiREdge[k][m_iA][m_iB][0]
						* m_dStateREdge[UIx][k]
					+ dContraMetricXiREdge[k][m_iA][m_iB][1]
						* m_dStateREdge[VIx][k]
					+ dContraMetricXiREdge[k][m_iA][m_iB][2]
						* dCovUx;
			}

			m_dXiDotREdge[0] = 0.0;
			m_dXiDotREdge[nRElements] = 0.0;
/*
			// dRhoQ_k/dRhoQ_n
			for (int k = 0; k < nRElements; k++) {

				int m = iDiffREdgeToNodeBegin[k];
				for (; m < iDiffREdgeToNodeEnd[k]; m++) {

					int n = iInterpNodeToREdgeBegin[m];
					for (; n < iInterpNodeToREdgeEnd[m]; n++) {

						dTracersLUDF[TracerMatFIx(n, k)] +=
							dDiffREdgeToNode[k][m]
							* dJacobianREdge[m][m_iA][m_iB]
							/ dJacobianNode[k][m_iA][m_iB]
							* dInterpNodeToREdge[m][n]
							* m_dXiDotREdge[m];
					}
				}
			}
*/
		}
	}

	// Add the identity components
	for (int k = 0; k < nRElements; k++) {
		dTracersLUDF[TracerMatFIx(k, k)] += 1.0 / m_dDeltaT;
	}
/*
		FILE * fpLUDF = fopen("TracersLUDF.txt", "w");
		for (int k = 0; k < nRElements; k++) {
			for (int i = 0; i < nRElements; i++) {
				fprintf(fpLUDF, "%1.15e\t", m_matTracersLUDF[k][i]);
			}
			fprintf(fpLUDF, "\n");
		}
		fclose(fpLUDF);
*/
#if defined(USE_JACOBIAN_GENERAL) || defined(USE_JACOBIAN_DEBUG)
	// LU Decomposition
	int iInfo =
		LAPACK::DGETRF(
			m_matTracersLUDF,
			m_vecTracersIPiv);
#elif defined(USE_JACOBIAN_DIAGONAL)
	// Banded diagonal LU decomposition
	int iInfo =
		LAPACK::DGBTRF(
			m_matTracersLUDF,
			m_vecTracersIPiv,
			2 * m_nVerticalOrder - 1,
			2 * m_nVerticalOrder - 1);
#else
	_EXCEPTIONT("Invalid Jacobian type");
#endif

	if (iInfo != 0) {
		_EXCEPTIONT("Triangulation failure");
	}

	// Calculate xi dot using the updated vertical velocity
	// with mass flux on levels
	if (fMassFluxOnLevels) {

		// Calculate u^xi on model levels
		if (pGrid->GetVarLocation(WIx) == DataLocation_REdge) {
			for (int k = 0; k <= nRElements; k++) {
				m_dStateREdge[WIx][k] = dataUpdateREdge[WIx][k][m_iA][m_iB];
			}

			pGrid->InterpolateREdgeToNode(
				m_dStateREdge[WIx],
				m_dStateNode[WIx]);

		} else {
			for (int k = 0; k < nRElements; k++) {
				m_dStateNode[WIx][k] = dataUpdateNode[WIx][k][m_iA][m_iB];
			}
		}

		for (int k = 1; k < nRElements-1; k++) {
			double dCovUx =
				m_dStateNode[WIx][k] * dDerivRNode[k][m_iA][m_iB][2];

			m_dXiDotNode[k] =
				  dContraMetricXi[k][m_iA][m_iB][0] * m_dStateNode[UIx][k]
				+ dContraMetricXi[k][m_iA][m_iB][1] * m_dStateNode[VIx][k]
				+ dContraMetricXi[k][m_iA][m_iB][2] * dCovUx;
		}

		m_dXiDotNode[0] = 0.0;
		m_dXiDotNode[nRElements-1] = 0.0;

	// Mass flux on model interfaces
	} else {

		// Calculate u^xi on model interfaces
		if (pGrid->GetVarLocation(WIx) == DataLocation_REdge) {
			for (int k = 0; k <= nRElements; k++) {
				m_dStateREdge[WIx][k] = dataUpdateREdge[WIx][k][m_iA][m_iB];
			}

		} else {
			for (int k = 0; k < nRElements; k++) {
				m_dStateNode[WIx][k] = dataUpdateNode[WIx][k][m_iA][m_iB];
			}

			pGrid->InterpolateNodeToREdge(
				m_dStateNode[WIx],
				m_dStateREdge[WIx]);
		}

		for (int k = 1; k < nRElements; k++) {
			double dCovUx =
				m_dStateREdge[WIx][k] * dDerivRREdge[k][m_iA][m_iB][2];

			m_dXiDotREdge[k] =
				  dContraMetricXiREdge[k][m_iA][m_iB][0]
					* m_dStateREdge[UIx][k]
				+ dContraMetricXiREdge[k][m_iA][m_iB][1]
					* m_dStateREdge[VIx][k]
				+ dContraMetricXiREdge[k][m_iA][m_iB][2]
					* dCovUx;
		}

		m_dXiDotREdge[0] = 0.0;
		m_dXiDotREdge[nRElements] = 0.0;
	}

	// Loop through all tracer species and apply update
	for (int c = 0; c < nComponents; c++) {

		// Calculate mass flux on nodes
		if (fMassFluxOnLevels) {
			for (int k = 0; k < nRElements; k++) {
				m_dMassFluxNode[k] =
					dJacobianNode[k][m_iA][m_iB]
					* dataInitialTracer[c][k][m_iA][m_iB]
					* m_dXiDotNode[k];
			}

			pGrid->DifferentiateNodeToNode(
				m_dMassFluxNode,
				m_dDiffMassFluxNode,
				fZeroBoundaries);

		// Calculate mass flux on interfaces
		} else {

			// Interpolate tracer density to interfaces
			for (int k = 0; k < nRElements; k++) {
				m_dTracerDensityNode[k] = dataInitialTracer[c][k][m_iA][m_iB];
			}

			pGrid->InterpolateNodeToREdge(
				m_dTracerDensityNode,
				m_dTracerDensityREdge);

			// Calculate mass flux
			for (int k = 0; k <= nRElements; k++) {
				m_dMassFluxREdge[k] =
					dJacobianREdge[k][m_iA][m_iB]
					* m_dTracerDensityREdge[k]
					* m_dXiDotREdge[k];
			}

			pGrid->DifferentiateREdgeToNode(
				m_dMassFluxREdge,
				m_dDiffMassFluxNode);
		}

		// Update tracers
		for (int k = 0; k < nRElements; k++) {
			m_vecTracersF[k] =
				m_dDiffMassFluxNode[k]
				/ dJacobianNode[k][m_iA][m_iB];
		}

		double dDiff = 0.0;
		double dSum = 0.0;
		for (int k = 0; k < nRElements; k++) {
			dDiff += m_vecTracersF[k] * dJacobianNode[k][m_iA][m_iB];
			dSum += m_dTracerDensityNode[k] * dJacobianNode[k][m_iA][m_iB];
		}

		if (c == 0) {
			//printf("%i %i %i %1.15e\n", m_pPatch->GetPatchIndex(), m_iA, m_iB, dSum);
			if ((fabs(dDiff) > 1.0e-12 * dSum) && (dSum > 1.0e-12)) {
				printf("%1.15e %1.15e\n", dDiff, dSum);
				for (int k = 0; k <= nRElements; k++) {
					printf("%i %1.15e\n", k, m_dMassFluxREdge[k]);
				}
				for (int k = 0; k < nRElements; k++) {
					printf("%i %1.15e\n", k, m_dDiffMassFluxNode[k]);
				}
				_EXCEPTION();
			}
		}

/*
		FILE * fpF = fopen("TracersF.txt", "w");
		for (int k = 0; k < nRElements; k++) {
			fprintf(fpF, "%1.15e\n", m_vecTracersF[k]);
		}
		fclose(fpF);
*/
/*
#if defined(USE_JACOBIAN_GENERAL) || defined(USE_JACOBIAN_DEBUG)
		// Solve the matrix system using LU decomposed matrix
		int iInfo =
			LAPACK::DGETRS(
				'N',
				m_matTracersLUDF,
				m_vecTracersF,
				m_vecTracersIPiv);

#elif defined(USE_JACOBIAN_DIAGONAL)
		// Solve the matrix system using banded LU decomposed matrix
		int iInfo =
			LAPACK::DGBTRS(
				'N',
				m_matTracersLUDF,
				m_vecTracersF,
				m_vecTracersIPiv,
				2 * m_nVerticalOrder - 1,
				2 * m_nVerticalOrder - 1);
#else
	_EXCEPTIONT("Invalid Jacobian type");
#endif
*/
		if (iInfo != 0) {
			_EXCEPTIONT("Inversion failure");
		}

/*
		FILE * fpX = fopen("TracersX.txt", "w");
		for (int k = 0; k < nRElements; k++) {
			fprintf(fpX, "%1.15e\n", m_vecTracersF[k]);
		}
		fclose(fpX);
*/

		// Update the state
		for (int k = 0; k < nRElements; k++) {
			dataUpdateTracer[c][k][m_iA][m_iB] -=
				m_vecTracersF[k];
		}
	}

/*
	for (int k = 0; k < nRElements; k++) {
		printf("%1.15e %1.15e\n",
			dataUpdateTracer[0][k][m_iA][m_iB],
			dataUpdateNode[RIx][k][m_iA][m_iB]);
	}
	_EXCEPTION();
*/
}

///////////////////////////////////////////////////////////////////////////////

void VerticalDynamicsFEM::FilterNegativeTracers(
	int iDataUpdate
) {
#ifdef POSITIVE_DEFINITE_FILTER_TRACERS
#pragma message "Apply limiter only to finite element?"
	// Get a copy of the GLL grid
	GridGLL * pGrid = dynamic_cast<GridGLL*>(m_model.GetGrid());

	// Number of vertical elements
	const int nRElements = pGrid->GetRElements();

	// Perform local update
	for (int n = 0; n < pGrid->GetActivePatchCount(); n++) {
		GridPatch * pPatch = pGrid->GetActivePatch(n);

		const PatchBox & box = pPatch->GetPatchBox();

		const DataArray3D<double> & dElementArea =
			pPatch->GetElementArea();

		DataArray4D<double> & dataUpdateTracer =
			pPatch->GetDataTracers(iDataUpdate);

		// Number of tracers
		const int nTracerCount = dataUpdateTracer.GetSize(0);

		for (int i = box.GetAInteriorBegin(); i < box.GetAInteriorEnd(); i++) {
		for (int j = box.GetBInteriorBegin(); j < box.GetBInteriorEnd(); j++) {

			for (int c = 0; c < nTracerCount; c++) {

				// Calculate total mass and non-negative mass
				double dTotalMass = 0.0;
				double dNonNegativeMass = 0.0;
				for (int k = 0; k < nRElements; k++) {
					double dPointwiseMass =
						  dataUpdateTracer[c][k][i][j]
						* dElementArea[k][i][j];

					dTotalMass += dPointwiseMass;

					if (dataUpdateTracer[c][k][i][j] >= 0.0) {
						dNonNegativeMass += dPointwiseMass;
					}
				}

				// Apply scaling ratio to points with non-negative mass
				double dR = dTotalMass / dNonNegativeMass;

				for (int k = 0; k < nRElements; k++) {
					if (dataUpdateTracer[c][k][i][j] > 0.0) {
						dataUpdateTracer[c][k][i][j] *= dR;
					} else {
						dataUpdateTracer[c][k][i][j] = 0.0;
					}
				}
			}
		}
		}
	}
#endif
}

///////////////////////////////////////////////////////////////////////////////
// GLOBAL CONTEXT
///////////////////////////////////////////////////////////////////////////////

#ifdef USE_JFNK_PETSC
PetscErrorCode VerticalDynamicsFEM_FormFunction(
	SNES snes,
	Vec x,
	Vec f,
	void * pDyn
) {
	// Pointers to the PetSc vector data
	const double * dX;
	double * dF;

	// Get pointers to the vector data
	VecGetArrayRead(x, &dX);
	VecGetArray(f, &dF);

	// Cast the context to VerticalDynamics and call Evaluate
	((VerticalDynamicsFEM*)(pDyn))->Evaluate(dX, dF);

	// Restore the array
	VecRestoreArrayRead(x, &dX);
	VecRestoreArray(f, &dF);

	return 0;
}
#endif
///////////////////////////////////////////////////////////////////////////////

