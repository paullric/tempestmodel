///////////////////////////////////////////////////////////////////////////////
///
///	\file    VerticalDynamicsFullImp.cpp
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
#include "VerticalDynamicsFullImp.h"
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

//#define UPWIND_HORIZONTAL_VELOCITIES
//#define UPWIND_THERMO
//#define UPWIND_VERTICAL_VELOCITY
//#define UPWIND_RHO_AND_TRACERS

#define UNIFORM_DIFFUSION_HORIZONTAL_VELOCITIES
#define UNIFORM_DIFFUSION_THERMO
#define UNIFORM_DIFFUSION_VERTICAL_VELOCITY
#define UNIFORM_DIFFUSION_TRACERS

#define VERTICAL_VELOCITY_ADVECTION_CLARK

///////////////////////////////////////////////////////////////////////////////

#define DEBUG

///////////////////////////////////////////////////////////////////////////////

VerticalDynamicsFullImp::VerticalDynamicsFullImp(
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

VerticalDynamicsFullImp::~VerticalDynamicsFullImp() {
#ifdef USE_JFNK_PETSC
	SNESDestroy(&m_snes);
	VecDestroy(&m_vecX);
	VecDestroy(&m_vecR);
	MatDestroy(&m_matJ);
#endif
}

///////////////////////////////////////////////////////////////////////////////

void VerticalDynamicsFullImp::Initialize() {

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
		VerticalDynamicsFullImp_FormFunction,
		(void*)(this));

	MatCreateSNESMF(m_snes, &m_matJ);

	SNESSetJacobian(m_snes, m_matJ, m_matJ, MatMFFDComputeJacobian, NULL);

	// Set the SNES context from options
	SNESSetFromOptions(m_snes);
#endif
#ifdef USE_JFNK_GMRES
	// Initialize JFNK
	InitializeJFNK(m_nColumnStateSize, m_nColumnStateSize, 1.0e-5);

	m_nJacobianFWidth = 1;
#endif
#if defined(USE_DIRECTSOLVE_APPROXJ) || defined(USE_DIRECTSOLVE)
#ifdef USE_JACOBIAN_DIAGONAL
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

	// Total bandwidth
	m_nJacobianFWidth = FTot * m_nJacobianFOffD + 1;
#endif
#endif

#if defined(USE_JACOBIAN_DIAGONAL)
	// Initialize Jacobian matrix
	m_matJacobianF.Allocate(m_nColumnStateSize, m_nJacobianFWidth);
#else
	// Initialize Jacobian matrix
	m_matJacobianF.Allocate(m_nColumnStateSize, m_nColumnStateSize);
#endif
	
	// Initialize pivot vector
	m_vecIPiv.Allocate(m_nColumnStateSize);

	// Announce vertical dynamics configuration
	AnnounceStartBlock("Configuring VerticalDynamicsFullImp");

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

	m_dCurlTermAREdge.Allocate(nRElements+1);
	m_dCurlTermBREdge.Allocate(nRElements+1);

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

void VerticalDynamicsFullImp::StepImplicitTermsExplicitly(
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
		}
		}
	}	     
}

///////////////////////////////////////////////////////////////////////////////

void VerticalDynamicsFullImp::StepExplicit(
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

				// Apply update to U
				if (pGrid->GetVarLocation(UIx) == DataLocation_REdge) {
					for (int k = 0; k <= nRElements; k++) {
						dataUpdateREdge[UIx][k][i][j] -=
							dDeltaT * m_dSoln[VecFIx(FUIx, k)];
					}
				} else {
					for (int k = 0; k < nRElements; k++) {
						dataUpdateNode[UIx][k][i][j] -=
							dDeltaT * m_dSoln[VecFIx(FUIx, k)];
					}
				}

				// Apply update to V
				if (pGrid->GetVarLocation(VIx) == DataLocation_REdge) {
					for (int k = 0; k <= nRElements; k++) {
						dataUpdateREdge[VIx][k][i][j] -=
							dDeltaT * m_dSoln[VecFIx(FVIx, k)];
					}
				} else {
					for (int k = 0; k < nRElements; k++) {
						dataUpdateNode[VIx][k][i][j] -=
							dDeltaT * m_dSoln[VecFIx(FVIx, k)];
					}
				}

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
			}
		}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void VerticalDynamicsFullImp::BootstrapJacobian() {

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

	//BuildJacobianF(m_dSoln, &(dJacobian[0][0]));

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

void VerticalDynamicsFullImp::StepImplicit(
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
#ifdef USE_DIRECTSOLVE_APPROXJ
			static const double Epsilon = 1.0e-5;

			// Prepare the column
			PrepareColumn(m_dColumnState);

			// Build the F vector
			BuildF(m_dColumnState, m_dSoln);
/*
			{
				double dDeltaEnergy = 0.0;
				for (int k = 0; k < pGrid->GetRElements(); k++) {
					double dPressure =
						phys.PressureFromRhoTheta(
							dataInitialNode[PIx][k][i][j]);

					dDeltaEnergy += dPressure / (phys.GetGamma() - 1.0);
				}

				printf("%1.15e\n", dDeltaEnergy);
				_EXCEPTION();
			}
*/
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

			// Copy over U
			if (pGrid->GetVarLocation(UIx) == DataLocation_REdge) {
				for (int k = 0; k <= pGrid->GetRElements(); k++) {
					dataUpdateREdge[UIx][k][iA][iB] =
						m_dSoln[VecFIx(FUIx, k)];
				}
			} else {
				for (int k = 0; k < pGrid->GetRElements(); k++) {
					dataUpdateNode[UIx][k][iA][iB] =
						m_dSoln[VecFIx(FUIx, k)];
				}
			}

			// Copy over V
			if (pGrid->GetVarLocation(VIx) == DataLocation_REdge) {
				for (int k = 0; k <= pGrid->GetRElements(); k++) {
					dataUpdateREdge[VIx][k][iA][iB] =
						m_dSoln[VecFIx(FVIx, k)];
				}
			} else {
				for (int k = 0; k < pGrid->GetRElements(); k++) {
					dataUpdateNode[VIx][k][iA][iB] =
						m_dSoln[VecFIx(FVIx, k)];
				}
			}

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
						dataUpdateNode[UIx][k][iA][iB]
							= dataUpdateNode[UIx][k][iA+1][iB];
						dataUpdateNode[VIx][k][iA][iB]
							= dataUpdateNode[VIx][k][iA+1][iB];
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
						dataUpdateREdge[UIx][k][iA][iB]
							= dataUpdateREdge[UIx][k][iA+1][iB];
						dataUpdateREdge[VIx][k][iA][iB]
							= dataUpdateREdge[VIx][k][iA+1][iB];
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
				dataUpdateNode[UIx][k][i][iB]
					= dataUpdateNode[UIx][k][i][iB+1];
				dataUpdateNode[VIx][k][i][iB]
					= dataUpdateNode[VIx][k][i][iB+1];
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
				dataUpdateREdge[UIx][k][i][iB]
					= dataUpdateREdge[UIx][k][i][iB+1];
				dataUpdateREdge[VIx][k][i][iB]
					= dataUpdateREdge[VIx][k][i][iB+1];
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
}


///////////////////////////////////////////////////////////////////////////////

void VerticalDynamicsFullImp::SetupReferenceColumn(
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

	// Get a copy of the Grid
	GridGLL * pGrid = dynamic_cast<GridGLL *>(m_model.GetGrid());

	// Get a copy of the PhysicalConstants
	const PhysicalConstants & phys = m_model.GetPhysicalConstants();

	// Number of radial elements
	const int nRElements = pGrid->GetRElements();

	// Store active patch in index
	m_pPatch = pPatch;
	m_iA = iA;
	m_iB = iB;

	// Copy over U
	if (pGrid->GetVarLocation(UIx) == DataLocation_REdge) {
		for (int k = 0; k <= pGrid->GetRElements(); k++) {
			m_dColumnState[VecFIx(FUIx, k)] =
				dataInitialREdge[UIx][k][iA][iB];
		}
	} else {
		for (int k = 0; k < pGrid->GetRElements(); k++) {
			m_dColumnState[VecFIx(FUIx, k)] =
				dataInitialNode[UIx][k][iA][iB];
		}
	}

	// Copy over V
	if (pGrid->GetVarLocation(VIx) == DataLocation_REdge) {
		for (int k = 0; k <= pGrid->GetRElements(); k++) {
			m_dColumnState[VecFIx(FVIx, k)] =
				dataInitialREdge[VIx][k][iA][iB];
		}
	} else {
		for (int k = 0; k < pGrid->GetRElements(); k++) {
			m_dColumnState[VecFIx(FVIx, k)] =
				dataInitialNode[VIx][k][iA][iB];
		}
	}

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

void VerticalDynamicsFullImp::PrepareColumn(
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

	// Get the array of interfaces
	const DataArray1D<double> & dREtaInterfaces = pGrid->GetREtaInterfaces();

	// Store state data from this column in state vector
	for (int k = 0; k < nRElements; k++) {

		// Store horizontal velocity
		m_dStateNode[UIx][k] = dX[VecFIx(FUIx, k)];
		m_dStateNode[VIx][k] = dX[VecFIx(FVIx, k)];

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

	// Interpolate horizontal velocities to interfaces
	pGrid->InterpolateNodeToREdge(
		m_dStateNode[UIx],
		m_dStateREdge[UIx]);

	pGrid->InterpolateNodeToREdge(
		m_dStateNode[VIx],
		m_dStateREdge[VIx]);

#if defined(VERTICAL_VELOCITY_ADVECTION_CLARK)
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
		for (int k = 0; k < nRElements; k++) {
			m_dStateNode[PIx][k] /= m_dStateNode[RIx][k];
		}

		pGrid->InterpolateNodeToREdge(
			m_dStateNode[PIx],
			m_dStateREdge[PIx]);

		for (int k = 0; k < nRElements; k++) {
			m_dStateNode[PIx][k] *= m_dStateNode[RIx][k];
		}
		for (int k = 0; k <= nRElements; k++) {
			m_dStateREdge[PIx][k] *= m_dStateREdge[RIx][k];
		}

		// Calculate Exner pressure at nodes
		for (int k = 0; k < nRElements; k++) {
			m_dExnerNode[k] =
				phys.ExnerPressureFromRhoTheta(
					m_dStateNode[PIx][k]);
		}

		for (int k = 1; k < nRElements; k++) {
			m_dDiffPREdge[k] =
				m_dExnerNode[k] / (dREtaInterfaces[k+1] - dREtaInterfaces[k])
				- m_dExnerNode[k-1] / (dREtaInterfaces[k] - dREtaInterfaces[k-1]);
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

	// Calculate curl terms
	if (pGrid->GetVarLocation(WIx) == DataLocation_Node) {
		_EXCEPTIONT("Not implemented");
	} else {
		for (int k = 1; k < nRElements; k++) {
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

			m_dCurlTermAREdge[k] = - dConUa * m_dDiffUa[k];
			m_dCurlTermBREdge[k] = - dConUb * m_dDiffUb[k];
		}
	}

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

void VerticalDynamicsFullImp::BuildF(
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

	// Z on levels
	const DataArray3D<double> & dZLevels = m_pPatch->GetZLevels();

	// Under this configuration, set fluxes at boundaries to zero
	bool fZeroBoundaries =
		(pGrid->GetVarLocation(WIx) == DataLocation_REdge);

	// Zero F
	memset(dF, 0, m_nColumnStateSize * sizeof(double));

	// Vertical advection of horizontal velocity
	const LinearColumnInterpFEM & opInterpNodeToREdge =
		pGrid->GetOpInterpNodeToREdge();
	const DataArray2D<double> & dInterpNodeToREdge =
		opInterpNodeToREdge.GetCoeffs();

	for (int k = 0; k < nRElements; k++) {
/*
		// Original method
		double dCovDxUa =
			pGrid->DifferentiateNodeToNode(
				m_dStateNode[UIx], k);

		double dCovDxUb =
			pGrid->DifferentiateNodeToNode(
				m_dStateNode[VIx], k);

		dF[VecFIx(FUIx, k)] +=
			m_dXiDotNode[k] * dCovDxUa;

		dF[VecFIx(FVIx, k)] +=
			m_dXiDotNode[k] * dCovDxUb;
*/

   		// Updated method
		for (int l = 0; l < nRElements; l++) {

			double dCovDxUa =
				pGrid->DifferentiateNodeToREdge(
					m_dStateNode[UIx], l);

			double dCovDxUb =
				pGrid->DifferentiateNodeToREdge(
					m_dStateNode[VIx], l);

			dF[VecFIx(FUIx, k)] +=
				m_dStateREdge[RIx][l]
				* m_dColumnJacobianREdge[l]
				* m_dXiDotREdge[l]
				* dInterpNodeToREdge[l][k]
				* dCovDxUa;

			dF[VecFIx(FVIx, k)] +=
				m_dStateREdge[RIx][l]
				* m_dColumnJacobianREdge[l]
				* m_dXiDotREdge[l]
				* dInterpNodeToREdge[l][k]
				* dCovDxUb;
		}

		dF[VecFIx(FUIx, k)] /=
			(m_dStateNode[RIx][k] * m_dColumnJacobianNode[k]);
		dF[VecFIx(FVIx, k)] /=
			(m_dStateNode[RIx][k] * m_dColumnJacobianNode[k]);
	}

	// Calculate mass flux on model interfaces
	for (int k = 1; k < nRElements; k++) {
		m_dMassFluxREdge[k] =
			m_dColumnJacobianREdge[k]
			* m_dStateREdge[RIx][k]
			* m_dXiDotREdge[k];
	}

	pGrid->DifferentiateREdgeToNode(
		m_dMassFluxREdge,
		m_dDiffMassFluxNode);

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
	for (int k = 1; k < nRElements; k++) {
		m_dPressureFluxREdge[k] =
			m_dColumnJacobianREdge[k]
			* m_dStateREdge[PIx][k]
			* m_dXiDotREdge[k];
	}

	pGrid->DifferentiateREdgeToNode(
		m_dPressureFluxREdge,
		m_dDiffPressureFluxNode);

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

		// Specific kinetic energy plus potential energy
		m_dKineticEnergyNode[k] =
			  0.5 * (dConUa * dCovUa + dConUb * dCovUb + dConUx * dCovUx);

		m_dKineticEnergyNode[k] += phys.GetG() * dZLevels[k][m_iA][m_iB];
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

#if !defined(EXPLICIT_VERTICAL_VELOCITY_ADVECTION)
#if defined(VERTICAL_VELOCITY_ADVECTION_CLARK)

			dF[VecFIx(FWIx, k)] +=
				(m_dDiffKineticEnergyREdge[k]
				+ m_dCurlTermAREdge[k]
				+ m_dCurlTermBREdge[k]
					) / m_dColumnDerivRREdge[k][2];

#else // VERTICAL VELOCITY ADVECTION (ADVECTIVE FORM)
			dF[VecFIx(FWIx, k)] +=
				m_dXiDotREdge[k] * m_dDiffWREdge[k];

#endif
#endif
		}
	}
/*
	// Calculate Energy difference
	{
		double dDeltaKEh = 0.0;
		double dDeltaKEv = 0.0;
		double dDeltaKErho = 0.0;

		for (int k = 1; k < nRElements; k++) {
			dDeltaKEv -=
				m_dStateREdge[RIx][k]
				* m_dXiDotREdge[k]
				* dF[VecFIx(FWIx,k)]
				* m_dColumnDerivRREdge[k][2];
		}
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

			dDeltaKEh -=
				m_dStateNode[RIx][k]
				* dConUa
				* dF[VecFIx(FUIx,k)];

			dDeltaKEh -=
				m_dStateNode[RIx][k]
				* dConUb
				* dF[VecFIx(FVIx,k)];
		}
		for (int k = 0; k < nRElements; k++) {
			dDeltaKErho -= dF[VecFIx(FRIx,k)]
			* (m_dKineticEnergyNode[k] - phys.GetG() * dZLevels[k][m_iA][m_iB]);
		}

		double dDeltaIE = 0.0;
		for (int k = 0; k < nRElements; k++) {
			double dExner =
				phys.ExnerPressureFromRhoTheta(m_dStateNode[PIx][k]);

			dDeltaIE -=
				dExner * dF[VecFIx(FPIx,k)];
		}

		double dDeltaPE = 0.0;
		for (int k = 0; k < nRElements; k++) {
			double dExner =
				phys.ExnerPressureFromRhoTheta(m_dStateNode[PIx][k]);

			dDeltaPE -= 
				dF[VecFIx(FRIx,k)] * phys.GetG() * dZLevels[k][m_iA][m_iB];
		}

		//printf("dKE:    %1.15e %1.15e\n", dDeltaKEh, dDeltaKEv);
		//printf("dKErho: %1.15e\n", dDeltaKErho);
		//printf("dIE/PE: %1.15e (%1.15e %1.15e)\n", dDeltaIE + dDeltaPE, dDeltaIE, dDeltaPE);
		printf("Delta:  %1.15e\n", dDeltaIE + dDeltaPE + dDeltaKEh + dDeltaKEv + dDeltaKErho);
		//_EXCEPTION();
	}
*/
	///////////////////////////////////////////////////////////////////
	// Apply uniform diffusion to theta and vertical velocity
	for (int c = 0; c < 4; c++) {

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

		// Uniform diffusion on interfaces
		if (pGrid->GetVarLocation(c) == DataLocation_REdge) {
			for (int k = 0; k <= nRElements; k++) {
				dF[VecFIx(FIxFromCIx(c), k)] -=
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
					/*printf("%1.15e %1.15e %1.15e\n",
						m_dUpwindCoeff,
						fabs(m_dXiDotREdge[k]),
						m_dDiffDiffStateUpwind[c][k]);*/
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
		for (int c = 0; c < 5; c++) {

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

void VerticalDynamicsFullImp::Evaluate(
	const double * dX,
	double * dF
) {
	// Prepare the column
	PrepareColumn(dX);

	// Evaluate the zero equations
	BuildF(dX, dF);
}

///////////////////////////////////////////////////////////////////////////////

