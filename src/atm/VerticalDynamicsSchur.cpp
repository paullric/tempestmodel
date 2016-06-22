///////////////////////////////////////////////////////////////////////////////
///
///	\file    VerticalDynamicsSchur.cpp
///	\author  Paul Ullrich
///	\version June 20, 2016
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
#include "VerticalDynamicsSchur.h"
#include "TimestepScheme.h"
#include "FunctionTimer.h"

#include "Announce.h"
#include "Model.h"
#include "Grid.h"
#include "GridCSGLL.h"
#include "EquationSet.h"
#include "TimeObj.h"
#include "PolynomialInterp.h"
#include "LinearAlgebra.h"

///////////////////////////////////////////////////////////////////////////////

//#define HYPERVISC_HORIZONTAL_VELOCITIES
//#define HYPERVISC_THERMO
//#define HYPERVISC_VERTICAL_VELOCITY

#define UPWIND_HORIZONTAL_VELOCITIES
#define UPWIND_THERMO
#define UPWIND_VERTICAL_VELOCITY
#define UPWIND_RHO_AND_TRACERS

#define UNIFORM_DIFFUSION_HORIZONTAL_VELOCITIES
#define UNIFORM_DIFFUSION_THERMO
#define UNIFORM_DIFFUSION_VERTICAL_VELOCITY
#define UNIFORM_DIFFUSION_TRACERS

//#define EXPLICIT_THERMO
#define EXPLICIT_VERTICAL_VELOCITY_ADVECTION

#define VERTICAL_VELOCITY_ADVECTION_CLARK

///////////////////////////////////////////////////////////////////////////////

#define DEBUG

///////////////////////////////////////////////////////////////////////////////

VerticalDynamicsSchur::VerticalDynamicsSchur(
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

VerticalDynamicsSchur::~VerticalDynamicsSchur() {
#ifdef USE_JFNK_PETSC
	SNESDestroy(&m_snes);
	VecDestroy(&m_vecX);
	VecDestroy(&m_vecR);
	MatDestroy(&m_matJ);
#endif
}

///////////////////////////////////////////////////////////////////////////////

void VerticalDynamicsSchur::Initialize() {

	// Indices of EquationSet variables
	const int UIx = 0;
	const int VIx = 1;
	const int PIx = 2;
	const int WIx = 3;
	const int RIx = 4;
	const int TracerIx = 5;

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
		VerticalDynamicsSchur_FormFunction,
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

#if defined(USE_DIRECTSOLVE_APPROXJ) || defined(USE_DIRECTSOLVE)
#if defined(USE_JACOBIAN_DIAGONAL)
	if (m_nHypervisOrder > 2) {
		_EXCEPTIONT("Diagonal Jacobian only implemented for "
			"Hypervis order <= 2");
	}

	// Upwind weights
	if (pGrid->GetVerticalDiscretization() ==
	    Grid::VerticalDiscretization_FiniteVolume
	) {
		if (m_nVerticalOrder <= 2) {
			m_nJacobianFOffD = 4;
		} else if (m_nVerticalOrder == 4) {
			m_nJacobianFOffD = 7;
		} else if (m_nVerticalOrder == 6) {
			m_nJacobianFOffD = 10;
		} else {
			_EXCEPTIONT("UNIMPLEMENTED: At this vertical order");
		}

	} else {
		if (m_nVerticalOrder == 1) {
			m_nJacobianFOffD = 4;
		} else if (m_nVerticalOrder == 2) {
			m_nJacobianFOffD = 9;
		} else if (m_nVerticalOrder == 3) {
			m_nJacobianFOffD = 15;
		} else if (m_nVerticalOrder == 4) {
			m_nJacobianFOffD = 22;
		} else if (m_nVerticalOrder == 5) {
			m_nJacobianFOffD = 30;
		} else {
			_EXCEPTIONT("UNIMPLEMENTED: At this vertical order");
		}
	}
#endif
#endif

	// Number of off-diagonals
	if (pGrid->GetVerticalDiscretization() ==
	    Grid::VerticalDiscretization_FiniteVolume
	) {
		m_nOffDiagonals = m_nVerticalOrder / 2;
		m_nJacobianFSchurOffD = m_nVerticalOrder + 1;
		m_nJacobianFBandwidth = 2 * m_nJacobianFOffD + 1;
		m_nJacobianFSchurWidth = 3 * m_nJacobianFSchurOffD + 1;

	} else {
		_EXCEPTIONT("Not implemented");
	}

#if defined(USE_JACOBIAN_DIAGONAL)
	// Initialize Jacobian matrix
	m_matJacobianF.Allocate(m_nColumnStateSize, m_nJacobianFBandwidth);

	// Initialize Schur complement Jacobian matrix
	m_matJacobianFSchur.Allocate(
		STot * m_nRElements,
		m_nJacobianFSchurWidth);

#elif defined(USE_JACOBIAN_GENERAL)
	// Initialize Jacobian matrix
	m_matJacobianF.Allocate(
		m_nColumnStateSize,
		m_nColumnStateSize);

	// Initialize Schur complement Jacobian matrix
	m_matJacobianFSchur.Allocate(
		STot * m_nRElements,
		STot * m_nRElements);

#else
	_EXCEPTIONT("Not implemented");
#endif

	// Initialize pivot vector
	m_vecIPiv.Allocate(m_nColumnStateSize);

	// Initialize solution in Schur complement
	m_dSolnSchur.Allocate(STot * m_nRElements);

	// Announce vertical dynamics configuration
	AnnounceStartBlock("Configuring VerticalDynamicsSchur");

	m_fHypervisVar.Allocate(6);
	m_fUpwindVar.Allocate(6);
	m_fUniformDiffusionVar.Allocate(6);

#if defined(HYPERVISC_HORIZONTAL_VELOCITIES)
	Announce("Hyperviscosity on horizontal velocities");
	m_fHypervisVar[UIx] = true;
	m_fHypervisVar[VIx] = true;
#endif
#if defined(HYPERVISC_THERMO)
	Announce("Hyperviscosity on thermodynamic variable");
	m_fHypervisVar[PIx] = true;
#endif
#if defined(HYPERVISC_VERTICAL_VELOCITY)
	Announce("Hyperviscosity on vertical velocity");
	m_fHypervisVar[WIx] = true;
#endif

#if defined(UPWIND_HORIZONTAL_VELOCITIES)
	Announce("Upwinding on horizontal velocities");
	m_fUpwindVar[UIx] = true;
	m_fUpwindVar[VIx] = true;
#endif
#if defined(UPWIND_THERMO)
	Announce("Upwinding on thermodynamic variable");
	m_fUpwindVar[PIx] = true;
#endif
#if defined(UPWIND_VERTICAL_VELOCITY)
	Announce("Upwinding on vertical velocity");
	m_fUpwindVar[WIx] = true;
#endif
#if defined(UPWIND_RHO_AND_TRACERS)
	Announce("Upwinding on rho and tracers");
	m_fUpwindVar[RIx] = true;
	m_fUpwindVar[TracerIx] = true;
#endif

#if defined(UNIFORM_DIFFUSION_HORIZONTAL_VELOCITIES)
	if (pGrid->HasUniformDiffusion()) {
		Announce("Uniform diffusion on horizontal velocities");
		m_fUniformDiffusionVar[UIx] = true;
		m_fUniformDiffusionVar[VIx] = true;
	}
#endif
#if defined(UNIFORM_DIFFUSION_THERMO)
	if (pGrid->HasUniformDiffusion()) {
		Announce("Uniform diffusion on thermodynamic variable");
		m_fUniformDiffusionVar[PIx] = true;
	}
#endif
#if defined(UNIFORM_DIFFUSION_VERTICAL_VELOCITY)
	if (pGrid->HasUniformDiffusion()) {
		Announce("Uniform diffusion on vertical velocity");
		m_fUniformDiffusionVar[WIx] = true;
	}
#endif
#if defined(UNIFORM_DIFFUSION_TRACERS)
	if (pGrid->HasUniformDiffusion()) {
		Announce("Uniform diffusion on tracers");
		m_fUniformDiffusionVar[TracerIx] = true;
	}
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
	if (pGrid->GetVerticalDiscretization() ==
	    Grid::VerticalDiscretization_FiniteVolume
	) {
		m_dUpwindWeights.Allocate(nRElements-1);

	} else {
		m_dUpwindWeights.Allocate(nRElements / m_nVerticalOrder - 1);
	}

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

	// Auxiliary variables
	m_dStateAux.Allocate(nRElements+1);
	m_dStateAuxDiff.Allocate(nRElements+1);

	m_dXiDotNode.Allocate(nRElements);
	m_dXiDotREdge.Allocate(nRElements+1);
	m_dXiDotREdgeInitial.Allocate(nRElements+1);

	m_dDiffUa.Allocate(nRElements+1);
	m_dDiffUb.Allocate(nRElements+1);

	m_dDiffPNode.Allocate(nRElements);
	m_dDiffPREdge.Allocate(nRElements+1);

	m_dDiffDiffStateUpwind.Allocate(
		m_model.GetEquationSet().GetComponents(),
		nRElements+1);

	m_dDiffDiffStateUniform.Allocate(
		m_model.GetEquationSet().GetComponents(),
		nRElements+1);

	m_dDiffDiffStateHypervis.Allocate(
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

	// Compute upwinding coefficient
	m_dUpwindCoeff = (1.0 / 2.0)
		* pow(1.0 / static_cast<double>(nRElements), 1.0);

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

void VerticalDynamicsSchur::StepImplicitTermsExplicitly(
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
/*
	for (std::ptrdiff_t e = 0; e < nRElements; ++e)
		m_dStateRefNode[WIx][e];
	for (std::ptrdiff_t e = 0; e < (nRElements+1); ++e)
		m_dStateRefREdge[WIx][e];
*/
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

		DataArray4D<double> & dataReferenceTracer =
			pPatch->GetReferenceTracers();

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
#if defined(EXPLICIT_VERTICAL_VELOCITY_ADVECTION)
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
				dataReferenceTracer,
				dataInitialTracer,
				dataUpdateTracer);
		}
		}
	}	     
}

///////////////////////////////////////////////////////////////////////////////

void VerticalDynamicsSchur::StepExplicit(
	int iDataInitial,
	int iDataUpdate,
	const Time & time,
	double dDeltaT
) {
	// Start the function timer
	FunctionTimer timer("VerticalStepExplicit");

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
	int nFiniteElements = nRElements / m_nVerticalOrder;
	int nNodesPerFiniteElement = m_nVerticalOrder;

	if (pGrid->GetVerticalDiscretization() ==
	    Grid::VerticalDiscretization_FiniteVolume
	) {
		nFiniteElements = nRElements;
		nNodesPerFiniteElement = 1;
	}

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

		DataArray4D<double> & dataReferenceTracer =
			pPatch->GetReferenceTracers();

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
/*
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
*/
		// Loop over all nodes
		for (int i = box.GetAInteriorBegin(); i < box.GetAInteriorEnd(); i++) {
		for (int j = box.GetBInteriorBegin(); j < box.GetBInteriorEnd(); j++) {

			// Setup the reference column:  Store U and V in m_dState arrays
			// and interpolate U and V from levels to interfaces.
			SetupReferenceColumn(
				pPatch, i, j,
				dataRefNode,
				dataInitialNode,
				dataRefREdge,
				dataInitialREdge);

			// Store W in m_dState structure on levels and interfaces
			if (pGrid->GetVarLocation(WIx) == DataLocation_Node) {
				for (int k = 0; k < nRElements; k++) {
					m_dStateNode[WIx][k] = dataInitialNode[WIx][k][i][j];
				}

				pGrid->InterpolateNodeToREdge(
					m_dStateNode[WIx],
					m_dStateREdge[WIx]);

			} else {
				for (int k = 0; k <= nRElements; k++) {
					m_dStateREdge[WIx][k] = dataInitialREdge[WIx][k][i][j];
				}

				pGrid->InterpolateREdgeToNode(
					m_dStateREdge[WIx],
					m_dStateNode[WIx]);
			}

			// Update thermodynamic variables
			if (m_fFullyExplicit) {

				// Evaluate the time tendency equations
				Evaluate(m_dColumnState, m_dSoln);

				// Apply update to P
				if (pGrid->GetVarLocation(PIx) == DataLocation_REdge) {
					for (int k = 0; k <= nRElements; k++) {
						dataUpdateREdge[PIx][k][i][j] -=
							dDeltaT * m_dSoln[VecFIx(FPIx, k)];
					}
				} else {
					for (int k = 0; k < nRElements; k++) {
						dataUpdateNode[PIx][k][i][j] -=
							dDeltaT * m_dSoln[VecFIx(FPIx, k)];
					}
				}

				// Apply update to W
				if (pGrid->GetVarLocation(WIx) == DataLocation_REdge) {
					for (int k = 0; k <= nRElements; k++) {
						dataUpdateREdge[WIx][k][i][j] -=
							dDeltaT * m_dSoln[VecFIx(FWIx, k)];
					}

				} else {
					for (int k = 0; k < nRElements; k++) {
						dataUpdateNode[WIx][k][i][j] -=
							dDeltaT * m_dSoln[VecFIx(FWIx, k)];
					}
				}

				// Apply update to Rho
				if (pGrid->GetVarLocation(RIx) == DataLocation_REdge) {
					for (int k = 0; k <= nRElements; k++) {
						dataUpdateREdge[RIx][k][i][j] -=
							dDeltaT * m_dSoln[VecFIx(FRIx, k)];
					}
				} else {
					for (int k = 0; k < nRElements; k++) {
						dataUpdateNode[RIx][k][i][j] -=
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
					dataReferenceTracer,
					dataInitialTracer,
					dataUpdateTracer);

			// Calculate xi dot on interfaces (this is done in Evaluate()
			// when m_fFullyExplicit is enabled)
			} else {
				for (int k = 0; k <= nRElements; k++) {
					double dCovUa = m_dStateREdge[UIx][k];
					double dCovUb = m_dStateREdge[VIx][k];
					double dCovUx =
						  m_dStateREdge[WIx][k]
						* m_dColumnDerivRREdge[k][2];

					m_dXiDotREdge[k] =
						  m_dColumnContraMetricXiREdge[k][0] * dCovUa
						+ m_dColumnContraMetricXiREdge[k][1] * dCovUb
						+ m_dColumnContraMetricXiREdge[k][2] * dCovUx;
				}

				m_dXiDotREdge[0] = 0.0;
				m_dXiDotREdge[nRElements] = 0.0;

				// Calculate u^xi on model levels (interpolated
				// from interfaces)
				pGrid->InterpolateREdgeToNode(
					m_dXiDotREdge,
					m_dXiDotNode);

/*
				for (int k = 0; k < nRElements; k++) {
					double dCovUa = m_dStateNode[UIx][k];
					double dCovUb = m_dStateNode[VIx][k];
					double dCovUx =
						  m_dStateNode[WIx][k]
						* m_dColumnDerivRNode[k][2];

					m_dXiDotNode[k] =
						  m_dColumnContraMetricXi[k][0] * dCovUa
						+ m_dColumnContraMetricXi[k][1] * dCovUb
						+ m_dColumnContraMetricXi[k][2] * dCovUx;
				}
*/
			}

			//////////////////////////////////////////////////////////////
			// Explicit vertical horizontal velocity advection
			for (int k = 0; k < nRElements; k++) {

				double dCovDxUa =
					pGrid->DifferentiateNodeToNode(
						m_dStateNode[UIx], k);

				double dCovDxUb =
					pGrid->DifferentiateNodeToNode(
						m_dStateNode[VIx], k);

				dataUpdateNode[UIx][k][i][j] -=
					dDeltaT * m_dXiDotNode[k] * dCovDxUa;

				dataUpdateNode[VIx][k][i][j] -=
					dDeltaT * m_dXiDotNode[k] * dCovDxUb;

			}

			//////////////////////////////////////////////////////////////
			// Explicit vertical velocity advection (Clark form)

#if defined(EXPLICIT_VERTICAL_VELOCITY_ADVECTION)
#if defined(VERTICAL_VELOCITY_ADVECTION_CLARK)

			{
				// Calculate specific kinetic energy on model levels
				for (int k = 0; k < nRElements; k++) {
					double dCovUa = m_dStateNode[UIx][k];
					double dCovUb = m_dStateNode[VIx][k];
					double dCovUx = m_dStateNode[WIx][k]
						* dDerivRNode[k][i][j][2];

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
						double dCovUa = m_dStateNode[UIx][k];
						double dCovUb = m_dStateNode[VIx][k];
						double dCovUx = m_dStateNode[WIx][k]
							* m_dColumnDerivRNode[k][2];

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

						dataUpdateNode[WIx][k][i][j] -=
							dDeltaT
							* (m_dDiffKineticEnergyNode[k] + dCurlTerm)
							/ m_dColumnDerivRNode[k][2];
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
						double dCovUa = m_dStateREdge[UIx][k];
						double dCovUb = m_dStateREdge[VIx][k];
						double dCovUx = m_dStateREdge[WIx][k]
							* m_dColumnDerivRREdge[k][2];

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

						dataUpdateREdge[WIx][k][i][j] -=
							dDeltaT
							* (m_dDiffKineticEnergyREdge[k] + dCurlTerm)
							/ m_dColumnDerivRREdge[k][2];
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
						double dDxW =
							pGrid->DifferentiateREdgeToREdge(
								&(m_dStateREdge[WIx][0]),
								k,
								nVerticalStateStride);

						dataUpdateREdge[WIx][k][i][j] -=
							dDeltaT * m_dXiDotREdge[k] * dDxW;
					}
				}
			}
#endif
#endif

#if defined(EXPLICIT_THERMO)
#if !defined(FORMULATION_THETA)
			_EXCEPTIONT("Not implemented: Explicit update of thermodynamic "
				"variable for RHOTHETA or P formulations");
#endif

			//////////////////////////////////////////////////////////////
			// Explicit update of thermodynamic equation
			{
				// Differentiate theta
				const LinearColumnDiffFEM & opDiffNodeToNode =
					pGrid->GetOpDiffNodeToNode();

				int nUpwindStride =
					dataInitialNode.GetSize(2)
					* dataInitialNode.GetSize(3);

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

			//////////////////////////////////////////////////////////////
			// Apply upwinding (discontinuous penalization)
			{
				// Calculate weights
				for (int a = 0; a < nFiniteElements - 1; a++) {
					int k = (a+1) * nNodesPerFiniteElement;
					m_dUpwindWeights[a] =
						dDeltaT * fabs(m_dXiDotREdge[k]);
				}

				// Apply upwinding
				const LinearColumnDiscPenaltyFEM & opPenalty =
					pGrid->GetOpPenaltyNodeToNode();

				int nUpwindStride =
					dataInitialNode.GetSize(2)
					* dataInitialNode.GetSize(3);

				// Apply upwinding to U and V
				if (m_fUpwindVar[UIx]) {
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

#if defined(EXPLICIT_THERMO)
				// Apply upwinding to thermodynamic variable
				if (m_fUpwindVar[PIx]) {
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
#if defined(EXPLICIT_VERTICAL_VELOCITY_ADVECTION)
				// Apply upwinding to vertical velocity
				if (m_fUpwindVar[WIx]) {
					if (pGrid->GetVarLocation(WIx) == DataLocation_REdge) {

						pGrid->DiffDiffREdgeToREdge(
							m_dStateREdge[WIx],
							m_dDiffDiffStateUpwind[WIx]);

						for (int k = 1; k < nRElements; k++) {
							dataUpdateNode[WIx][k][i][j] +=
								m_dUpwindCoeff
								* fabs(m_dXiDotREdge[k])
								* m_dDiffDiffStateUpwind[WIx][k];
						}

					} else {
						opPenalty.Apply(
							&(m_dUpwindWeights[0]),
							&(dataInitialNode[WIx][0][i][j]),
							&(dataUpdateNode[WIx][0][i][j]),
							nUpwindStride,
							nUpwindStride);
					}
				}
#endif
			}

			//////////////////////////////////////////////////////////////
			// Apply hyperviscosity or uniform diffusion to U and V
			if ((m_fHypervisVar[UIx]) ||
			    (m_fUniformDiffusionVar[UIx])
			) {

				// Second derivatives of horizontal velocity on model levels
				pGrid->DiffDiffNodeToNode(
					m_dStateNode[UIx],
					m_dDiffDiffStateHypervis[UIx]);

				pGrid->DiffDiffNodeToNode(
					m_dStateNode[VIx],
					m_dDiffDiffStateHypervis[VIx]);

				// Apply uniform diffusion in the vertical
				if (m_fUniformDiffusionVar[UIx]) {
					double dZtop = pGrid->GetZtop();

					double dUniformDiffusionCoeff =
						pGrid->GetVectorUniformDiffusionCoeff()
						/ (dZtop * dZtop);

					for (int k = 0; k < nRElements; k++) {
						m_dStateRefNode[UIx][k] = dataRefNode[UIx][k][i][j];
						m_dStateRefNode[VIx][k] = dataRefNode[VIx][k][i][j];
					}

					pGrid->DiffDiffNodeToNode(
						m_dStateRefNode[UIx],
						m_dDiffDiffStateUniform[UIx]);

					pGrid->DiffDiffNodeToNode(
						m_dStateRefNode[VIx],
						m_dDiffDiffStateUniform[VIx]);

					for (int k = 0; k < nRElements; k++) {
						dataUpdateNode[UIx][k][i][j] +=
							dDeltaT
							* dUniformDiffusionCoeff
							* (m_dDiffDiffStateHypervis[UIx][k]
								- m_dDiffDiffStateUniform[UIx][k]);

						dataUpdateNode[VIx][k][i][j] +=
							dDeltaT
							* dUniformDiffusionCoeff
							* (m_dDiffDiffStateHypervis[VIx][k]
								- m_dDiffDiffStateUniform[VIx][k]);
					}
				}

				// Apply hyperviscosity in the vertical
				if (m_fHypervisVar[UIx]) {

					// No hyperviscosity from command line
					if (m_nHypervisOrder == 0) {
						continue;
					}

					// Compute higher derivatives of U and V used for
					// hyperviscosity
					for (int h = 2; h < m_nHypervisOrder; h += 2) {
						memcpy(
							m_dStateAux,
							m_dDiffDiffStateHypervis[UIx],
							nRElements * sizeof(double));

						pGrid->DiffDiffNodeToNode(
							m_dStateAux,
							m_dDiffDiffStateHypervis[UIx]
						);

						memcpy(
							m_dStateAux,
							m_dDiffDiffStateHypervis[VIx],
							nRElements * sizeof(double));

						pGrid->DiffDiffNodeToNode(
							m_dStateAux,
							m_dDiffDiffStateHypervis[VIx]
						);
					}

					// Apply hyperviscosity
					for (int k = 0; k < nRElements; k++) {
						dataUpdateNode[UIx][k][i][j] +=
							dDeltaT
							* m_dHypervisCoeff
							* fabs(m_dXiDotNode[k])
							* m_dDiffDiffStateHypervis[UIx][k];

						dataUpdateNode[VIx][k][i][j] +=
							dDeltaT
							* m_dHypervisCoeff
							* fabs(m_dXiDotNode[k])
							* m_dDiffDiffStateHypervis[VIx][k];
					}
				}
			}
		}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void VerticalDynamicsSchur::DebugJacobian() {

			double * dDG = &(m_matJacobianF[0][0]);

			double * dDGSchur = &(m_matJacobianFSchur[0][0]);

			{
				FILE * fp = fopen("F.txt", "w");
				for (int i = 0; i < m_nColumnStateSize; i++) {
					fprintf(fp, "%1.15e", m_dSoln[i]);
					if (i != m_nColumnStateSize-1) {
						fprintf(fp, "\n");
					}
				}
				fclose(fp);
			}

			{
				FILE * fp = fopen("DG.txt", "w");
				for (int i = 0; i < m_matJacobianF.GetRows(); i++) {
					for (int j = 0; j < m_matJacobianF.GetColumns(); j++) {
						fprintf(fp, "%1.15e", m_matJacobianF[i][j]);
						if (j != m_nColumnStateSize-1) {
							fprintf(fp, " ");
						}
					}
					fprintf(fp, "\n");
				}
				fclose(fp);
			}

			{
				FILE * fp = fopen("FSchur.txt", "w");
				for (int i = 0; i < m_dSolnSchur.GetRows(); i++) {
					fprintf(fp, "%1.15e", m_dSolnSchur[i]);
					if (i != m_dSolnSchur.GetRows()-1) {
						fprintf(fp, "\n");
					}
				}
				fclose(fp);
			}

			{
				FILE * fp = fopen("DGSchur.txt", "w");
				for (int i = 0; i < m_matJacobianFSchur.GetRows(); i++) {
					for (int j = 0; j < m_matJacobianFSchur.GetColumns(); j++) {
						fprintf(fp, "%1.15e", m_matJacobianFSchur[i][j]);
						if (j != m_matJacobianFSchur.GetColumns()-1) {
							fprintf(fp, " ");
						}
					}
					fprintf(fp, "\n");
				}
				fclose(fp);
			}

			// Solve for the solution with the general Jacobian
			DataArray2D<double> dBackupJacobian(
				m_matJacobianF.GetRows(),
				m_matJacobianF.GetColumns());

			DataArray1D<double> dBackupSoln(
				m_dSoln.GetRows());

			dBackupJacobian = m_matJacobianF;
			dBackupSoln = m_dSoln;

			LAPACK::DGESV(dBackupJacobian, dBackupSoln, m_vecIPiv);

			{
				FILE * fp = fopen("X.txt", "w");
				for (int i = 0; i < dBackupSoln.GetRows(); i++) {
					fprintf(fp, "%1.15e", dBackupSoln[i]);
					if (i != dBackupSoln.GetRows()-1) {
						fprintf(fp, "\n");
					}
				}
				fclose(fp);
			}

			// Solve the Schur complement
			LAPACK::DGESV(m_matJacobianFSchur, m_dSolnSchur, m_vecIPiv);

			// Back out the full solution
			memcpy(&(m_dSoln[0]), &(m_dSolnSchur[0]),
				2 * (m_nRElements+1) * sizeof(double));

			for (int i = 0; i <= m_nRElements; i++) {
				m_dSoln[VecFIx(FWIx,i)] /= dDG[MatFIx(FWIx, i, FWIx, i)];
			}

			for (int i = 0; i <= m_nRElements; i++) {
				for (int l = 0; l < m_nRElements; l++) {
					m_dSoln[VecFIx(FWIx,i)] -=
						dDG[MatFIx(FPIx, l, FWIx, i)]
						/ dDG[MatFIx(FWIx, i, FWIx, i)]
						* m_dSolnSchur[VecFIx(FPIx,l)];

					m_dSoln[VecFIx(FWIx,i)] -=
						dDG[MatFIx(FRIx, l, FWIx, i)]
						/ dDG[MatFIx(FWIx, i, FWIx, i)]
						* m_dSolnSchur[VecFIx(FRIx,l)];
				}
			}

			{
				FILE * fp = fopen("XSchur.txt", "w");
				for (int i = 0; i < m_dSoln.GetRows(); i++) {
					fprintf(fp, "%1.15e", m_dSoln[i]);
					if (i != m_dSoln.GetRows()-1) {
						fprintf(fp, "\n");
					}
				}
				fclose(fp);
			}

			_EXCEPTION();
}

///////////////////////////////////////////////////////////////////////////////

void VerticalDynamicsSchur::StepImplicit(
	int iDataInitial,
	int iDataUpdate,
	const Time & time,
	double dDeltaT
) {
	// Start the function timer
	FunctionTimer timer("VerticalStepImplicit");

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
		DataArray4D<double> & dataReferenceTracer =
			pPatch->GetReferenceTracers();

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
#ifdef USE_DIRECTSOLVE
			// Prepare the column
			PrepareColumn(m_dColumnState);

			// Build the F vector
			BuildF(m_dColumnState, m_dSoln);

			// Build the Jacobian
			m_matJacobianF.Zero();
			
			int nDim = m_dColumnState.GetRows();

			BuildJacobianF(m_dColumnState, &(m_matJacobianF[0][0]));

			// Build the Schur complement array
			double * dDG = &(m_matJacobianF[0][0]);

			double * dDGSchur = &(m_matJacobianFSchur[0][0]);

			m_matJacobianFSchur.Zero();

			for (int k = 0; k < m_nRElements; k++) {
				m_dSolnSchur[VecSIx(SPIx,k)] =
					m_dSoln[VecFIx(FPIx,k)];
				m_dSolnSchur[VecSIx(SRIx,k)] =
					m_dSoln[VecFIx(FRIx,k)];
			}

			for (int k = 0; k < m_nRElements; k++) {

				int ibegin = k;
				int iend = k + 2;
				if (ibegin < 1) {
					ibegin = 1;
				}
				if (iend > m_nRElements) {
					iend = m_nRElements;
				}

				for (int i = ibegin; i < iend; i++) {

					m_dSolnSchur[VecSIx(SPIx,k)] -=
						dDeltaT
						* dDG[MatFIx(FWIx, i, FPIx, k)]
						// dDG[MatFIx(FWIx, i, FWIx, i)]
						* m_dSoln[VecFIx(FWIx,i)];

					m_dSolnSchur[VecSIx(SRIx,k)] -=
						dDeltaT
						* dDG[MatFIx(FWIx, i, FRIx, k)]
						// dDG[MatFIx(FWIx, i, FWIx, i)]
						* m_dSoln[VecFIx(FWIx,i)];

					int lbegin = i - m_nOffDiagonals;
					int lend = i + m_nOffDiagonals + 1;
					if (lbegin < 0) {
						lbegin = 0;
					}
					if (lend > m_nRElements) {
						lend = m_nRElements;
					}

					if (i == ibegin) {
						for (int l = lbegin; l < lend; l++) {
							dDGSchur[MatSIx(SPIx, l, SPIx, k)] =
								dDG[MatFIx(FPIx, l, FPIx, k)];
							dDGSchur[MatSIx(SPIx, l, SRIx, k)] =
								dDG[MatFIx(FPIx, l, FRIx, k)];
							dDGSchur[MatSIx(SRIx, l, SPIx, k)] =
								dDG[MatFIx(FRIx, l, FPIx, k)];
							dDGSchur[MatSIx(SRIx, l, SRIx, k)] =
								dDG[MatFIx(FRIx, l, FRIx, k)];
						}
					}

					for (int l = lbegin; l < lend; l++) {
						dDGSchur[MatSIx(SPIx, l, SPIx, k)] -=
							dDeltaT
							* dDG[MatFIx(FWIx, i, FPIx, k)]
							// dDG[MatFIx(FWIx, i, FWIx, i)]
							* dDG[MatFIx(FPIx, l, FWIx, i)];

						dDGSchur[MatSIx(SPIx, l, SRIx, k)] -=
							dDeltaT
							* dDG[MatFIx(FWIx, i, FRIx, k)]
							// dDG[MatFIx(FWIx, i, FWIx, i)]
							* dDG[MatFIx(FPIx, l, FWIx, i)];

						dDGSchur[MatSIx(SRIx, l, SPIx, k)] -=
							dDeltaT
							* dDG[MatFIx(FWIx, i, FPIx, k)]
							// dDG[MatFIx(FWIx, i, FWIx, i)]
							* dDG[MatFIx(FRIx, l, FWIx, i)];

						dDGSchur[MatSIx(SRIx, l, SRIx, k)] -=
							dDeltaT
							* dDG[MatFIx(FWIx, i, FRIx, k)]
							// dDG[MatFIx(FWIx, i, FWIx, i)]
							* dDG[MatFIx(FRIx, l, FWIx, i)];
					}
				}
			}

			// Debug
			//DebugJacobian();
/*
			{
				FILE * fp = fopen("DG.txt", "w");
				for (int i = 0; i < m_matJacobianF.GetRows(); i++) {
					for (int j = 0; j < m_matJacobianF.GetColumns(); j++) {
						fprintf(fp, "%1.15e", m_matJacobianF[i][j]);
						if (j != m_matJacobianF.GetColumns()-1) {
							fprintf(fp, " ");
						}
					}
					fprintf(fp, "\n");
				}
				fclose(fp);
			}
			_EXCEPTION();
*/
/*
			{
				FILE * fp = fopen("DGSchur.txt", "w");
				for (int i = 0; i < m_matJacobianFSchur.GetRows(); i++) {
					for (int j = 0; j < m_matJacobianFSchur.GetColumns(); j++) {
						fprintf(fp, "%1.15e", m_matJacobianFSchur[i][j]);
						if (j != m_matJacobianFSchur.GetColumns()-1) {
							fprintf(fp, " ");
						}
					}
					fprintf(fp, "\n");
				}
				fclose(fp);
			}
			_EXCEPTION();
*/
#if defined(USE_JACOBIAN_GENERAL)
			// Solve using the general matrix solver
			int iInfo = LAPACK::DGESV(
				m_matJacobianFSchur, m_dSolnSchur, m_vecIPiv);

			if (iInfo != 0) {
				_EXCEPTION1("Solution failed: %i", iInfo);
			}
/*
			// Use direct solver
			LAPACK::DGESV(m_matJacobianF, m_dSoln, m_vecIPiv);

			if (iInfo != 0) {
				_EXCEPTION1("Solution failed: %i", iInfo);
			}
*/
#endif
#if defined(USE_JACOBIAN_DIAGONAL)
			// Solve using the diagonal matrix solver
			int iInfo = LAPACK::DGBSV(
				m_matJacobianFSchur, m_dSolnSchur, m_vecIPiv,
				m_nJacobianFSchurOffD, m_nJacobianFSchurOffD);

			if (iInfo != 0) {
				_EXCEPTION1("Solution failed: %i", iInfo);
			}
#endif

			// Back out the full solution
			for (int k = 0; k < m_nRElements; k++) {
				m_dSoln[VecFIx(FPIx,k)] = m_dSolnSchur[VecSIx(SPIx,k)];
				m_dSoln[VecFIx(FRIx,k)] = m_dSolnSchur[VecSIx(SRIx,k)];

				m_dSoln[VecFIx(FWIx,k)] *= dDeltaT;
					// 1.0 /dDG[MatFIx(FWIx, i, FWIx, i)];
			}
			m_dSoln[VecFIx(FWIx,m_nRElements)] = 0.0;

			for (int i = 1; i < m_nRElements; i++) {

				int lbegin = i - m_nOffDiagonals;
				int lend = i + m_nOffDiagonals + 1;
				if (lbegin < 0) {
					lbegin = 0;
				}
				if (lend > m_nRElements) {
					lend = m_nRElements;
				}

				for (int l = lbegin; l < lend; l++) {
					m_dSoln[VecFIx(FWIx,i)] -=
						dDeltaT
						* dDG[MatFIx(FPIx, l, FWIx, i)]
						// dDG[MatFIx(FWIx, k, FWIx, k)]
						* m_dSolnSchur[VecSIx(SPIx,l)];

					m_dSoln[VecFIx(FWIx,i)] -=
						dDeltaT
						* dDG[MatFIx(FRIx, l, FWIx, i)]
						// dDG[MatFIx(FWIx, k, FWIx, k)]
						* m_dSolnSchur[VecSIx(SRIx,l)];
				}
			}

/*
			for (int k = 0; k <= m_nRElements; k++) {
				printf("%1.10e %1.10e %1.10e\n",
					m_dSoln[VecFIx(FPIx,k)],
					m_dSoln[VecFIx(FRIx,k)],
					m_dSoln[VecFIx(FWIx,k)]);
			}
			_EXCEPTION();
*/
			// DEBUG (check for NANs in output)
			if (!(m_dSoln[0] == m_dSoln[0])) {
				DataArray1D<double> dEval;
				dEval.Allocate(m_dColumnState.GetRows());
				Evaluate(m_dSoln, dEval);

				for (int p = 0; p < dEval.GetRows(); p++) {
					printf("%1.15e %1.15e %1.15e\n",
						dEval[p], m_dSoln[p] - m_dColumnState[p],
						m_dColumnState[p]);
				}
				for (int p = 0; p < m_dExnerRefREdge.GetRows(); p++) {
					printf("%1.15e %1.15e\n",
						m_dExnerRefREdge[p], dataRefREdge[RIx][p][iA][iB]);
				}
				_EXCEPTIONT("Inversion failure");
			}

			// Update the solution
			for (int k = 0; k < m_dSoln.GetRows(); k++) {
				m_dSoln[k] = m_dColumnState[k] - m_dSoln[k];
			}
#endif

#if defined(EXPLICIT_THERMO)
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
				dataReferenceTracer,
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

void VerticalDynamicsSchur::SetupReferenceColumn(
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
	for (int k = 0; k < nRElements; k++) {
		m_dStateNode[UIx][k] = dataInitialNode[UIx][k][iA][iB];
	}

	if (pGrid->GetVarsAtLocation(DataLocation_REdge) != 0) {
		pGrid->InterpolateNodeToREdge(
			m_dStateNode[UIx],
			m_dStateREdge[UIx]);
	}

	// Store V in State structure
	for (int k = 0; k < nRElements; k++) {
		m_dStateNode[VIx][k] = dataInitialNode[VIx][k][iA][iB];
	}

	if (pGrid->GetVarsAtLocation(DataLocation_REdge) != 0) {
		pGrid->InterpolateNodeToREdge(
			m_dStateNode[VIx],
			m_dStateREdge[VIx]);
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

	// Construct reference column
	if ((pGrid->HasUniformDiffusion()) && (m_fUseReferenceState)) {
		for (int k = 0; k < pGrid->GetRElements(); k++) {
			m_dStateRefNode[RIx][k] = dataRefNode[RIx][k][iA][iB];
			m_dStateRefNode[WIx][k] = dataRefNode[WIx][k][iA][iB];
			m_dStateRefNode[PIx][k] = dataRefNode[PIx][k][iA][iB];
		}
		for (int k = 0; k <= pGrid->GetRElements(); k++) {
			m_dStateRefREdge[RIx][k] = dataRefREdge[RIx][k][iA][iB];
			m_dStateRefREdge[WIx][k] = dataRefREdge[WIx][k][iA][iB];
			m_dStateRefREdge[PIx][k] = dataRefREdge[PIx][k][iA][iB];
		}
	}

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

void VerticalDynamicsSchur::PrepareColumn(
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
	for (int k = 1; k < nRElements; k++) {
		double dCovUx =
			m_dStateREdge[WIx][k] * m_dColumnDerivRREdge[k][2];

		m_dXiDotREdge[k] =
			  m_dColumnContraMetricXiREdge[k][0] * m_dStateREdge[UIx][k]
			+ m_dColumnContraMetricXiREdge[k][1] * m_dStateREdge[VIx][k]
			+ m_dColumnContraMetricXiREdge[k][2] * dCovUx;
	}

	m_dXiDotREdge[0] = 0.0;
	m_dXiDotREdge[nRElements] = 0.0;

#if !defined(EXPLICIT_VERTICAL_VELOCITY_ADVECTION) \
 && !defined(VERTICAL_VELOCITY_ADVECTION_CLARK)
	// Calculate vertical derivatives of W needed for advective form of W
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

	// Calculate second derivatives needed for upwinding on interfaces.
	for (int c = 2; c < 5; c++) {
		if (!m_fUpwindVar[c]) {
			continue;
		}

		if (pGrid->GetVarLocation(c) == DataLocation_REdge) {
			pGrid->DiffDiffREdgeToREdge(
				m_dStateREdge[c],
				m_dDiffDiffStateUpwind[c]);
		}
	}

	// Calculate higher-order even derivatives needed for hyperviscosity
	// and uniform diffusion.
	for (int c = 2; c < 5; c++) {

		// Only apply hypervis or uniform diffusion to select variables
		if ((!m_fHypervisVar[c]) && (!m_fUniformDiffusionVar[c])) {
			continue;
		}

		// Variable on model interfaces
		if (pGrid->GetVarLocation(c) == DataLocation_REdge) {
			pGrid->DiffDiffREdgeToREdge(
				m_dStateREdge[c],
				m_dDiffDiffStateHypervis[c]);

			// Calculate second-derivative of variable minus reference
			if (m_fUniformDiffusionVar[c]) {
				pGrid->DiffDiffREdgeToREdge(
					m_dStateRefREdge[c],
					m_dDiffDiffStateUniform[c]);

				for (int k = 0; k <= nRElements; k++) {
					m_dDiffDiffStateUniform[c][k] =
						m_dDiffDiffStateHypervis[c][k]
						- m_dDiffDiffStateUniform[c][k];
				}
			}

			// Calculate higher-order derivatives for hyperviscosity
			if (m_fHypervisVar[c]) {
				for (int h = 2; h < m_nHypervisOrder; h += 2) {
					memcpy(
						m_dStateAux,
						m_dDiffDiffStateHypervis[c],
						(nRElements+1) * sizeof(double));
	
					pGrid->DiffDiffREdgeToREdge(
						m_dStateAux,
						m_dDiffDiffStateHypervis[c]
					);
				}
			}

		// Variable on model levels
		} else {
			pGrid->DiffDiffNodeToNode(
				m_dStateNode[c],
				m_dDiffDiffStateHypervis[c]);

			// Calculate second-derivative of variable minus reference
			if (m_fUniformDiffusionVar[c]) {
				pGrid->DiffDiffNodeToNode(
					m_dStateRefNode[c],
					m_dDiffDiffStateUniform[c]);

				for (int k = 0; k < nRElements; k++) {
					m_dDiffDiffStateUniform[c][k] =
						m_dDiffDiffStateHypervis[c][k]
						- m_dDiffDiffStateUniform[c][k];
				}
			}

			// Calculate higher-order derivatives for hyperviscosity
			if (m_fHypervisVar[c]) {
				for (int h = 2; h < m_nHypervisOrder; h += 2) {
					memcpy(
						m_dStateAux,
						m_dDiffDiffStateHypervis[c],
						nRElements * sizeof(double));

					pGrid->DiffDiffNodeToNode(
						m_dStateAux,
						m_dDiffDiffStateHypervis[c]
					);
				}
			}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void VerticalDynamicsSchur::BuildF(
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

#if !defined(EXPLICIT_THERMO)
#if defined(FORMULATION_PRESSURE)
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
#if defined(FORMULATION_THETA)
	// Update theta on model levels
	if (pGrid->GetVarLocation(PIx) == DataLocation_Node) {

		// Test using u^xi on interfaces interpolated to nodes for Lorenz
		pGrid->InterpolateREdgeToNode(
				m_dXiDotREdge,
				m_dXiDotNode);

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
#if defined(FORMULATION_THETA_FLUX)
	// Update theta on model levels
	if (pGrid->GetVarLocation(PIx) == DataLocation_Node) {

		// Test using u^xi on interfaces interpolated to nodes for Lorenz
		pGrid->InterpolateREdgeToNode(
				m_dXiDotREdge,
				m_dXiDotNode);

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

	///////////////////////////////////////////////////////////////////
	// Apply uniform diffusion to theta and vertical velocity
	for (int c = 2; c < 4; c++) {

		// Only upwind select variables
		if (!m_fUniformDiffusionVar[c]) {
			continue;
		}

		double dZtop = pGrid->GetZtop();

		// Do not diffusion vertical velocity on boundaries
		if (c == WIx) {
			m_dDiffDiffStateUniform[c][0] = 0.0;
			m_dDiffDiffStateUniform[c][nRElements] = 0.0;
		}

		// Uniform diffusion coefficient
		double dUniformDiffusionCoeff = 0.0;
		if (c == PIx) {
			dUniformDiffusionCoeff =
				pGrid->GetScalarUniformDiffusionCoeff() / (dZtop * dZtop);
		}
		if (c == WIx) {
			dUniformDiffusionCoeff =
				pGrid->GetVectorUniformDiffusionCoeff() / (dZtop * dZtop);
		}

		// Get fix
		const int fix = FIxFromCIx(c);

		// Uniform diffusion on interfaces
		if (pGrid->GetVarLocation(c) == DataLocation_REdge) {
			for (int k = 0; k <= nRElements; k++) {
				dF[VecFIx(fix, k)] -=
					dUniformDiffusionCoeff
					* m_dDiffDiffStateUniform[c][k];
			}

		// Uniform diffusion on levels
		} else {
			for (int k = 0; k < nRElements; k++) {
				dF[VecFIx(FIxFromCIx(c), k)] -=
					dUniformDiffusionCoeff
					* m_dDiffDiffStateUniform[c][k];
			}
		}
	}

	///////////////////////////////////////////////////////////////////
	// Apply vertical upwinding
	{
		// Get penalty operator
		const LinearColumnDiscPenaltyFEM & opPenalty =
			pGrid->GetOpPenaltyNodeToNode();

		// Number of finite elements in the vertical
		int nFiniteElements = nRElements / m_nVerticalOrder;
		int nNodesPerFiniteElement = m_nVerticalOrder;

		if (pGrid->GetVerticalDiscretization() ==
		    Grid::VerticalDiscretization_FiniteVolume
		) {
			nFiniteElements = nRElements;
			nNodesPerFiniteElement = 1;
		}

		// Calculate weights
		for (int a = 0; a < nFiniteElements - 1; a++) {
			int k = (a+1) * nNodesPerFiniteElement;
			m_dUpwindWeights[a] = fabs(m_dXiDotREdge[k]);
		}

		// Loop through all variables
		for (int c = 2; c < 5; c++) {

			// Only upwind select variables
			if (!m_fUpwindVar[c]) {
				continue;
			}

#if defined(EXPLICIT_THERMO)
			// Don't upwind thermodynamic variable if explicit
			if (c == PIx) {
				continue;
			}
#endif

#if defined(EXPLICIT_VERTICAL_VELOCITY_ADVECTION)
			// Don't upwind vertical velocity if explicit
			if (c == WIx) {
				continue;
			}
#endif

			// Upwinding on interfaces (use flow-dependent constant
			// coefficient hyperviscosity at second order)
			if (pGrid->GetVarLocation(c) == DataLocation_REdge) {
				//if (m_nVerticalOrder != 1) {
				//	_EXCEPTIONT("Upwinding on interfaces only at "
				//		"VO1 supported");
				//}

				// No upwinding on W on domain boundaries
				if (c == WIx) {
					m_dDiffDiffStateUpwind[c][0] = 0.0;
					m_dDiffDiffStateUpwind[c][nRElements] = 0.0;
				}

				for (int k = 0; k <= nRElements; k++) {
					dF[VecFIx(FIxFromCIx(c), k)] -=
						m_dUpwindCoeff
						* fabs(m_dXiDotREdge[k])
						* m_dDiffDiffStateUpwind[c][k];
				}

			// Upwinding on levels (discontinuous penalization)
			} else {
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

	// Apply flow-dependent hyperviscosity
	if (m_nHypervisOrder > 0) {
		for (int c = 2; c < 5; c++) {

			// Only upwind select variables
			if (!m_fHypervisVar[c]) {
				continue;
			}

			// Flow-dependent hyperviscosity on interfaces
			if (pGrid->GetVarLocation(c) == DataLocation_REdge) {

				for (int k = 0; k <= nRElements; k++) {
					dF[VecFIx(FIxFromCIx(c), k)] -=
						m_dHypervisCoeff
						* fabs(m_dXiDotREdge[k])
						* m_dDiffDiffStateHypervis[c][k];
				}

			// Flow-dependent hyperviscosity on levels
			} else {

				for (int k = 0; k < nRElements; k++) {
					dF[VecFIx(FIxFromCIx(c), k)] -=
						m_dHypervisCoeff
						* fabs(m_dXiDotNode[k])
						* m_dDiffDiffStateHypervis[c][k];
				}
			}
		}
	}

	if (dF[VecFIx(FWIx, 0)] != 0.0) {
		dF[VecFIx(FWIx, 0)] = 0.0;
	}
	if (pGrid->GetVarLocation(WIx) == DataLocation_REdge) {
		if (dF[VecFIx(FWIx, nRElements)] != 0.0) {
			dF[VecFIx(FWIx, nRElements)] = 0.0;
		}
	} else {
		if (dF[VecFIx(FWIx, nRElements-1)] != 0.0) {
			dF[VecFIx(FWIx, nRElements-1)] = 0.0;
		}
	}

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
}

///////////////////////////////////////////////////////////////////////////////

void VerticalDynamicsSchur::BuildJacobianF_Diffusion(
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

	// Physical constants
	const PhysicalConstants & phys = m_model.GetPhysicalConstants();

	// Number of radial elements
	const int nRElements = pGrid->GetRElements();

	// Number of finite elements in the vertical
	int nFiniteElements = nRElements / m_nVerticalOrder;
	int nNodesPerFiniteElement = m_nVerticalOrder;

	if (pGrid->GetVerticalDiscretization() ==
	    Grid::VerticalDiscretization_FiniteVolume
	) {
		nFiniteElements = nRElements;
		nNodesPerFiniteElement = 1;
	}

	// Get the column penalization coefficients (for vertical upwinding)
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

	// Second derivative coefficients (for vertical hyperviscosity)
	const LinearColumnDiffDiffFEM & opDiffDiffNodeToNode =
		pGrid->GetOpDiffDiffNodeToNode();
	const LinearColumnDiffDiffFEM & opDiffDiffREdgeToREdge =
		pGrid->GetOpDiffDiffREdgeToREdge();

	const DataArray2D<double> & dDiffDiffNodeToNode =
		opDiffDiffNodeToNode.GetCoeffs();
	const DataArray2D<double> & dDiffDiffREdgeToREdge =
		opDiffDiffREdgeToREdge.GetCoeffs();

	const DataArray1D<int> & iDiffDiffNodeToNodeBegin =
		opDiffDiffNodeToNode.GetIxBegin();
	const DataArray1D<int> & iDiffDiffREdgeToREdgeBegin =
		opDiffDiffREdgeToREdge.GetIxBegin();

	const DataArray1D<int> & iDiffDiffNodeToNodeEnd =
		opDiffDiffNodeToNode.GetIxEnd();
	const DataArray1D<int> & iDiffDiffREdgeToREdgeEnd =
		opDiffDiffREdgeToREdge.GetIxEnd();

	/////////////////////////////////////////////////////////////////////
	// Vertical upwinding
	for (int c = 2; c < 5; c++) {

#if defined(EXPLICIT_THERMO)
		// Don't upwind thermodynamic variable if explicit
		if (c == PIx) {
			continue;
		}
#endif

#if defined(EXPLICIT_VERTICAL_VELOCITY_ADVECTION)
		// Don't upwind vertical velocity if explicit
		if (c == WIx) {
			continue;
		}
#endif

		// Check upwinding
		if (!m_fUpwindVar[c]) {
			continue;
		}
		if (pGrid->GetVarLocation(WIx) != DataLocation_REdge) {
			_EXCEPTIONT("Upwinding DIRECTSOLVE requires W on interfaces");
		}

		// Convert to Jacobian index
		int fix = FIxFromCIx(c);

		// Upwinding on interfaces (velocity-weighted hyperviscosity)
		if (pGrid->GetVarLocation(c) == DataLocation_REdge) {

			// dC_k/dW_k
			for (int k = 0; k <= nRElements; k++) {
				double dSignWeight;
				if (m_dXiDotREdge[k] > 0.0) {
					dSignWeight = 1.0;
				} else if (m_dXiDotREdge[k] < 0.0) {
					dSignWeight = -1.0;
				} else {
					dSignWeight = 0.0;
				}

				dDG[MatFIx(FWIx, k, fix, k)] -=
					m_dUpwindCoeff
					* dSignWeight
					/ m_dColumnDerivRREdge[k][2]
					* m_dDiffDiffStateUpwind[c][k];
			}

			// dC_k/dC_n
			for (int k = 0; k <= nRElements; k++) {
				int n = iDiffDiffREdgeToREdgeBegin[k];
				for (; n < iDiffDiffREdgeToREdgeEnd[k]; n++) {
					dDG[MatFIx(fix, n, fix, k)] -=
						m_dUpwindCoeff
						* fabs(m_dXiDotREdge[k])
						* dDiffDiffREdgeToREdge[k][n];
				}
			}

		// Upwinding on levels
		} else {

			for (int a = 1; a < nFiniteElements; a++) {
				double dWeight =
					fabs(m_dXiDotREdge[a * nNodesPerFiniteElement]);

				double dSignWeight;
				if (m_dXiDotREdge[a * nNodesPerFiniteElement] > 0.0) {
					dSignWeight = 1.0;
				} else if (m_dXiDotREdge[a * nNodesPerFiniteElement] < 0.0) {
					dSignWeight = -1.0;
				} else {
					dSignWeight = 0.0;
				}

				int kLeftBegin = (a-1) * nNodesPerFiniteElement;
				int kLeftEnd = a * nNodesPerFiniteElement;

				int kRightBegin = a * nNodesPerFiniteElement;
				int kRightEnd = (a+1) * nNodesPerFiniteElement;

				// dC_k/dW_a (left operator)
				for (int k = kLeftBegin; k < kLeftEnd; k++) {
					int n = iPenaltyLeftBegin[k];
					for (; n < iPenaltyLeftEnd[k]; n++) {
						dDG[MatFIx(FWIx, kLeftEnd, fix, k)] -=
							dSignWeight
							/ m_dColumnDerivRREdge[kLeftEnd][2]
							* dPenaltyLeft[k][n]
							* m_dStateNode[c][n];
					}
				}

				// dC_k/dW_a (right operator)
				for (int k = kRightBegin; k < kRightEnd; k++) {
					int n = iPenaltyRightBegin[k];
					for (; n < iPenaltyRightEnd[k]; n++) {
						dDG[MatFIx(FWIx, kRightBegin, fix, k)] -=
							dSignWeight
							/ m_dColumnDerivRREdge[kRightBegin][2]
							* dPenaltyRight[k][n]
							* m_dStateNode[c][n];
					}
				}

				// dC_k/dC_n (left operator)
				for (int k = kLeftBegin; k < kLeftEnd; k++) {
					int n = iPenaltyLeftBegin[k];
					for (; n < iPenaltyLeftEnd[k]; n++) {
						dDG[MatFIx(fix, n, fix, k)] -=
							dWeight * dPenaltyLeft[k][n];
					}
				}

				// dC_k/dC_n (right operator)
				for (int k = kRightBegin; k < kRightEnd; k++) {
					int n = iPenaltyRightBegin[k];
					for (; n < iPenaltyRightEnd[k]; n++) {
						dDG[MatFIx(fix, n, fix, k)] -=
							dWeight * dPenaltyRight[k][n];
					}
				}
			}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void VerticalDynamicsSchur::BuildJacobianF(
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

	// Physical constants
	const PhysicalConstants & phys = m_model.GetPhysicalConstants();

	// Number of radial elements
	const int nRElements = pGrid->GetRElements();

	// Zero DG
	memset(dDG, 0,
		m_nColumnStateSize * m_nJacobianFBandwidth * sizeof(double));

#if defined(EXPLICIT_THERMO)
	// Only support implicit thermodynamics
	_EXCEPTIONT("Not implemented");
#endif
#if !defined(FORMULATION_RHOTHETA_PI)
	// Only support FORMULATION_RHOTHETA_PI
	_EXCEPTIONT("Invalid call to BuildJacobianF_LOR_RhoTheta_Pi");
#endif

	// Only support Lorenz staggering
	if (pGrid->GetVarLocation(PIx) == DataLocation_REdge) {
		_EXCEPTIONT("Invalid call to BuildJacobianF_LOR_RhoTheta_Pi");
	}

	// Jacobian terms for Rho and RhoTheta updates
	for (int k = 0; k < nRElements; k++) {

		// Conservative flux terms
		int m = iDiffREdgeToNodeBegin[k];
		for (; m < iDiffREdgeToNodeEnd[k]; m++) {

			if ((m != 0) && (m != nRElements)) {
				double dMassFluxCoeff =
					dDiffREdgeToNode[k][m]
					* m_dColumnJacobianREdge[m]
					* m_dColumnInvJacobianNode[k]
					* m_dColumnContraMetricXiREdge[m][2]
					* m_dColumnDerivRREdge[m][2];

				// dP_k/dW_l
				dDG[MatFIx(FWIx, m, FPIx, k)] +=
					dMassFluxCoeff * m_dStateREdge[PIx][m];

				// dR_k/dW_l
				dDG[MatFIx(FWIx, m, FRIx, k)] +=
					dMassFluxCoeff * m_dStateREdge[RIx][m];
			}

			// dR_k/dR_n
			int n = iInterpNodeToREdgeBegin[m];
			for (; n < iInterpNodeToREdgeEnd[m]; n++) {

				double dCoeffVerticalFlux = 
					dDiffREdgeToNode[k][m]
					* m_dColumnJacobianREdge[m]
					* m_dColumnInvJacobianNode[k]
					* dInterpNodeToREdge[m][n]
					* m_dXiDotREdge[m];

				dDG[MatFIx(FRIx, n, FRIx, k)] += dCoeffVerticalFlux;

				dDG[MatFIx(FPIx, n, FPIx, k)] += dCoeffVerticalFlux;
			}
		}
	}

	// Jacobian terms for W updates
	for (int k = 1; k < nRElements; k++) {

		// dW_k/dP_m (in pressure gradient)
		double dRHSWCoeffA = 
			m_dStateREdge[PIx][k]
			* phys.GetR()
			/ (m_dColumnDerivRNode[k][2]
				* m_dStateREdge[RIx][k]
				* phys.GetCv());

		int m = iDiffNodeToREdgeBegin[k];
		for (; m < iDiffNodeToREdgeEnd[k]; m++) {
			dDG[MatFIx(FPIx, m, FWIx, k)] +=
				dRHSWCoeffA 
				* dDiffNodeToREdge[k][m]
				* m_dExnerNode[m]
				/ m_dStateNode[PIx][m];
		}

		// Rhotheta/rho term in front of pressure gradient
		double dRHSWCoeffB = 
			1.0 / m_dColumnDerivRREdge[k][2]
			/ (m_dStateREdge[RIx][k] * m_dStateREdge[RIx][k])
			* m_dDiffPREdge[k];

		int n = iInterpNodeToREdgeBegin[k];
		for (; n < iInterpNodeToREdgeEnd[k]; n++) {

			double dRHSWCoeffC = dRHSWCoeffB * dInterpNodeToREdge[k][n];

			// dW_k/dP_n (first rhotheta in RHS)
			dDG[MatFIx(FPIx, n, FWIx, k)] +=
				dRHSWCoeffC * m_dStateREdge[RIx][k];

			// dW_k/dR_l (first rho in RHS)
			dDG[MatFIx(FRIx, n, FWIx, k)] +=
				- dRHSWCoeffC * m_dStateREdge[PIx][k];
		}
	}

#if !defined(EXPLICIT_VERTICAL_VELOCITY_ADVECTION)
	_EXCEPTIONT("Not implemented");
/*
	// dW_k/dW_m (Clark form)
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
*/
#endif

	// Add the diffusion terms
	BuildJacobianF_Diffusion(dX, dDG);

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

void VerticalDynamicsSchur::Evaluate(
	const double * dX,
	double * dF
) {
	// Prepare the column
	PrepareColumn(dX);

	// Evaluate the zero equations
	BuildF(dX, dF);
}

///////////////////////////////////////////////////////////////////////////////

void VerticalDynamicsSchur::UpdateColumnTracers(
	double dDeltaT,
	const DataArray4D<double> & dataInitialNode,
	const DataArray4D<double> & dataUpdateNode,
	const DataArray4D<double> & dataInitialREdge,
	const DataArray4D<double> & dataUpdateREdge,
	const DataArray4D<double> & dataRefTracer,
	const DataArray4D<double> & dataInitialTracer,
	DataArray4D<double> & dataUpdateTracer
) {
	// Indices of EquationSet variables
	const int UIx = 0;
	const int VIx = 1;
	const int PIx = 2;
	const int WIx = 3;
	const int RIx = 4;
	const int TracerIx = 5;

	// Finite element grid
	const GridGLL * pGrid = dynamic_cast<const GridGLL *>(m_model.GetGrid());

	// Number of model levels
	const int nRElements = pGrid->GetRElements();

	// Number of finite elements in the vertical
	int nFiniteElements = nRElements / m_nVerticalOrder;
	int nNodesPerFiniteElement = m_nVerticalOrder;

	if (pGrid->GetVerticalDiscretization() ==
	    Grid::VerticalDiscretization_FiniteVolume
	) {
		nFiniteElements = nRElements;
		nNodesPerFiniteElement = 1;
	}

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

	// Penalty operators
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

		// No implicit uniform diffusion of tracers (yet)
		if (m_fUniformDiffusionVar[TracerIx]) {
			_EXCEPTIONT("Not implemented");
		}

		// Calculate u^xi on model interfaces and nodes
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

			pGrid->InterpolateNodeToREdge(
				m_dStateNode[WIx],
				m_dStateREdge[WIx]);
		}

#pragma message "Replace with column arrays"
		// Calculate u^xi on interfaces
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

			m_dXiDotREdgeInitial[k] = m_dXiDotREdge[k];
		}

		m_dXiDotREdge[0] = 0.0;
		m_dXiDotREdge[nRElements] = 0.0;

		m_dXiDotREdgeInitial[0] = 0.0;
		m_dXiDotREdgeInitial[nRElements] = 0.0;

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

		// Apply upwinding to Jacobian
		if (m_fUpwindVar[TracerIx]) {
			for (int a = 1; a < nFiniteElements; a++) {

				int kLeftBegin = (a-1) * nNodesPerFiniteElement;
				int kLeftEnd = a * nNodesPerFiniteElement;

				int kRightBegin = a * nNodesPerFiniteElement;
				int kRightEnd = (a+1) * nNodesPerFiniteElement;

				double dWeight = fabs(m_dXiDotREdge[kLeftEnd]);

				double dSignWeight;
				if (m_dXiDotREdge[kLeftEnd] > 0.0) {
					dSignWeight = 1.0;
				} else if (m_dXiDotREdge[kLeftEnd] < 0.0) {
					dSignWeight = -1.0;
				} else {
					dSignWeight = 0.0;
				}

				// dRhoQ_k/dRhoQ_n (left operator)
				for (int k = kLeftBegin; k < kLeftEnd; k++) {
					int n = iPenaltyLeftBegin[k];
					for (; n < iPenaltyLeftEnd[k]; n++) {
						dTracersLUDF[TracerMatFIx(n, k)] -=
							dWeight * dPenaltyLeft[k][n];
					}
				}

				// dRhoQ_k/dRhoQ_n (right operator)
				for (int k = kRightBegin; k < kRightEnd; k++) {
					int n = iPenaltyRightBegin[k];
					for (; n < iPenaltyRightEnd[k]; n++) {
						dTracersLUDF[TracerMatFIx(n, k)] -=
							dWeight * dPenaltyRight[k][n];
					}
				}
			}
		}
	}

	// Add the identity components
	for (int k = 0; k < nRElements; k++) {
		dTracersLUDF[TracerMatFIx(k, k)] += 1.0 / m_dDeltaT;
	}

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

	// Calculate u^xi on model interfaces (using updated state if implicit)
	if (m_fFullyExplicit) {
		if (pGrid->GetVarLocation(WIx) == DataLocation_REdge) {
			for (int k = 0; k <= nRElements; k++) {
				m_dStateREdge[WIx][k] =
					dataInitialREdge[WIx][k][m_iA][m_iB];
			}

		} else {
			for (int k = 0; k < nRElements; k++) {
				m_dStateNode[WIx][k] =
					dataInitialNode[WIx][k][m_iA][m_iB];
			}
		}

	} else {
		if (pGrid->GetVarLocation(WIx) == DataLocation_REdge) {
			for (int k = 0; k <= nRElements; k++) {
				m_dStateREdge[WIx][k] =
					dataUpdateREdge[WIx][k][m_iA][m_iB];
			}

		} else {
			for (int k = 0; k < nRElements; k++) {
				m_dStateNode[WIx][k] =
					dataUpdateNode[WIx][k][m_iA][m_iB];
			}
		}
	}

	if (pGrid->GetVarLocation(WIx) == DataLocation_Node) {
		pGrid->InterpolateNodeToREdge(
			m_dStateNode[WIx],
			m_dStateREdge[WIx]);
	}

#pragma message "Change to column arrays"
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

	// Loop through all tracer species and apply update
	for (int c = 0; c < nComponents; c++) {

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

		////////////////////////////////////////////////////////////
		// Apply uniform diffusion
		if (m_fUniformDiffusionVar[TracerIx]) {
			for (int k = 0; k < nRElements; k++) {
				m_dStateAux[k] =
					m_dTracerDensityNode[k]
					/ m_dStateNode[RIx][k];

				m_dStateAux[k] -=
					dataRefTracer[c][k][m_iA][m_iB]
					/ m_dStateRefNode[RIx][k];
			}

			pGrid->DifferentiateNodeToREdge(
				m_dStateAux,
				m_dStateAuxDiff);

			for (int k = 1; k < nRElements; k++) {
				m_dMassFluxREdge[k] -=
					pGrid->GetScalarUniformDiffusionCoeff()
					* m_dStateREdge[RIx][k]
					* m_dStateAuxDiff[k];
			}
		}

		// Set boundary conditions and differentiate mass flux
		m_dMassFluxREdge[0] = 0.0;
		m_dMassFluxREdge[nRElements] = 0.0;

		pGrid->DifferentiateREdgeToNode(
			m_dMassFluxREdge,
			m_dDiffMassFluxNode);

		// Update tracers
		for (int k = 0; k < nRElements; k++) {
			m_vecTracersF[k] =
				m_dDiffMassFluxNode[k]
				/ dJacobianNode[k][m_iA][m_iB];
		}

		////////////////////////////////////////////////////////////
		// Apply upwinding to tracers
		if (m_fUpwindVar[TracerIx]) {

			// Get penalty operator
			const LinearColumnDiscPenaltyFEM & opPenalty =
				pGrid->GetOpPenaltyNodeToNode();

			// Calculate weights
			if (m_fFullyExplicit) {
				for (int a = 0; a < nFiniteElements - 1; a++) {
					int k = (a+1) * nNodesPerFiniteElement;
					m_dUpwindWeights[a] = fabs(m_dXiDotREdge[k]);
				}

			} else {
				for (int a = 0; a < nFiniteElements - 1; a++) {
					int k = (a+1) * nNodesPerFiniteElement;
					m_dUpwindWeights[a] = fabs(m_dXiDotREdgeInitial[k]);
				}
			}

			// Apply upwinding
			m_dStateAux.Zero();
			opPenalty.Apply(
				&(m_dUpwindWeights[0]),
				&(m_dTracerDensityNode[0]),
				&(m_dStateAux[0]),
				1,
				1);

			for (int k = 0; k < nRElements; k++) {
				m_vecTracersF[k] -= m_dStateAux[k];
			}

			// Apply implicit velocity correction
			if (!m_fFullyExplicit) {
				for (int a = 1; a < nFiniteElements; a++) {

					int kLeftBegin = (a-1) * nNodesPerFiniteElement;
					int kLeftEnd = a * nNodesPerFiniteElement;

					int kRightBegin = a * nNodesPerFiniteElement;
					int kRightEnd = (a+1) * nNodesPerFiniteElement;

					double dWeight =
						fabs(m_dXiDotREdgeInitial[kLeftEnd]);

					double dSignWeight;
					if (m_dXiDotREdgeInitial[kLeftEnd] > 0.0) {
						dSignWeight = 1.0;
					} else if (m_dXiDotREdgeInitial[kLeftEnd] < 0.0) {
						dSignWeight = -1.0;
					} else {
						dSignWeight = 0.0;
					}

					// dRhoQ_k/dW_a (left operator)
					double dLeftJumpConUx =
						dSignWeight
						* (dataUpdateREdge[WIx][kLeftEnd][m_iA][m_iB]
							- m_dColumnState[VecFIx(FWIx, kLeftEnd)])
						/ m_dColumnDerivRREdge[kLeftEnd][2];

					for (int k = kLeftBegin; k < kLeftEnd; k++) {
						int n = iPenaltyLeftBegin[k];
						for (; n < iPenaltyLeftEnd[k]; n++) {
							m_vecTracersF[k] -=
								dPenaltyLeft[k][n]
								* m_dTracerDensityNode[n]
								* dLeftJumpConUx;
						}
					}

					// dRhoQ_k/dW_a (right operator)
					double dRightJumpConUx =
						dSignWeight
						* (dataUpdateREdge[WIx][kRightBegin][m_iA][m_iB]
							- m_dColumnState[VecFIx(FWIx, kRightBegin)])
						/ m_dColumnDerivRREdge[kRightBegin][2];

					for (int k = kRightBegin; k < kRightEnd; k++) {
						int n = iPenaltyRightBegin[k];
						for (; n < iPenaltyRightEnd[k]; n++) {
							m_vecTracersF[k] -=
								dPenaltyRight[k][n]
								* m_dTracerDensityNode[n]
								* dRightJumpConUx;
						}
					}
				}
			}
		}

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

		if (iInfo != 0) {
			_EXCEPTIONT("Inversion failure");
		}

		// Update the state
		for (int k = 0; k < nRElements; k++) {
			dataUpdateTracer[c][k][m_iA][m_iB] -=
				m_vecTracersF[k];
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void VerticalDynamicsSchur::FilterNegativeTracers(
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
PetscErrorCode VerticalDynamicsSchur_FormFunction(
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
	((VerticalDynamicsSchur*)(pDyn))->Evaluate(dX, dF);

	// Restore the array
	VecRestoreArrayRead(x, &dX);
	VecRestoreArray(f, &dF);

	return 0;
}
#endif
///////////////////////////////////////////////////////////////////////////////

