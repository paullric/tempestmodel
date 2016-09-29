///////////////////////////////////////////////////////////////////////////////
///
///	\file    VerticalDynamicsFLL.cpp
///	\author  Paul Ullrich
///	\version June 29, 2016
///
///	<remarks>
///		Copyright 2000-2016 Paul Ullrich
///
///		This file is distributed as part of the Tempest source code package.
///		Permission is granted to use, copy, modify and distribute this
///		source code and its documentation under the terms of the GNU General
///		Public License.  This software is provided "as is" without express
///		or implied warranty.
///	</remarks>

#include "Defines.h"
#include "VerticalDynamicsFLL.h"
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
//#define EXPLICIT_VERTICAL_VELOCITY_ADVECTION

#define VERTICAL_VELOCITY_ADVECTION_CLARK

///////////////////////////////////////////////////////////////////////////////

#define DEBUG

///////////////////////////////////////////////////////////////////////////////

VerticalDynamicsFLL::VerticalDynamicsFLL(
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

VerticalDynamicsFLL::~VerticalDynamicsFLL() {
#ifdef USE_JFNK_PETSC
	SNESDestroy(&m_snes);
	VecDestroy(&m_vecX);
	VecDestroy(&m_vecR);
	MatDestroy(&m_matJ);
#endif
}

///////////////////////////////////////////////////////////////////////////////

void VerticalDynamicsFLL::Initialize() {

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
		VerticalDynamicsFLL_FormFunction,
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
	m_nJacobianFWidth = 3 * m_nJacobianFOffD + 1;
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
	AnnounceStartBlock("Configuring VerticalDynamicsFLL");

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

void VerticalDynamicsFLL::StepExplicit(
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

		// Loop over all nodes
		for (int i = box.GetAInteriorBegin(); i < box.GetAInteriorEnd(); i++) {
		for (int j = box.GetBInteriorBegin(); j < box.GetBInteriorEnd(); j++) {

			// Update thermodynamic variables
			if (m_fFullyExplicit) {

				// Setup the reference column:  Store U and V in m_dState
				// arrays and interpolate U and V from levels to interfaces.
				SetupReferenceColumn(
					pPatch, i, j,
					dataRefNode,
					dataInitialNode,
					dataRefREdge,
					dataInitialREdge);

				// Evaluate the time tendency equations
				Evaluate(m_dColumnState, m_dSoln);

				// Update density and rhotheta
				for (int k = 0; k < pGrid->GetRElements(); k++) {
					double dInitialDeltaZ =
						m_dColumnState[VecFIx(FZIx, k+1)]
						- m_dColumnState[VecFIx(FZIx, k)];

					double dUpdateDeltaZ =
						m_dSoln[VecFIx(FZIx, k+1)]
						- m_dSoln[VecFIx(FZIx, k)];

					dataUpdateNode[PIx][k][i][j] =
						dataInitialNode[PIx][k][i][j]
						* dInitialDeltaZ
						/ dUpdateDeltaZ;

					dataUpdateNode[RIx][k][i][j] =
						dataInitialNode[RIx][k][i][j]
						* dInitialDeltaZ
						/ dUpdateDeltaZ;
				}

/*
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
*/
			}
		}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void VerticalDynamicsFLL::BootstrapJacobian() {

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

void VerticalDynamicsFLL::StepImplicit(
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
			// Use direct solver
			int iInfo = LAPACK::DGESV(m_matJacobianF, m_dSoln, m_vecIPiv);

			if (iInfo != 0) {
				_EXCEPTION1("Solution failed: %i", iInfo);
			}
#endif
#ifdef USE_JACOBIAN_DIAGONAL
			// Use diagonal solver
			int iInfo = LAPACK::DGBSV(
				m_matJacobianF, m_dSoln, m_vecIPiv,
				m_nJacobianFOffD, m_nJacobianFOffD);

			if (iInfo != 0) {
				_EXCEPTION1("Solution failed: %i", iInfo);
			}
#endif

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

			for (int k = 0; k < m_dSoln.GetRows(); k++) {
				m_dSoln[k] = m_dColumnState[k] - m_dSoln[k];
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

			// Update density and rhotheta
			for (int k = 0; k < pGrid->GetRElements(); k++) {
				double dInitialDeltaZ =
					m_dColumnState[VecFIx(FZIx, k+1)]
					- m_dColumnState[VecFIx(FZIx, k)];

				double dUpdateDeltaZ =
					m_dSoln[VecFIx(FZIx, k+1)]
					- m_dSoln[VecFIx(FZIx, k)];

				dataUpdateNode[PIx][k][iA][iB] =
					dataInitialNode[PIx][k][iA][iB]
					* dInitialDeltaZ
					/ dUpdateDeltaZ;

				dataUpdateNode[RIx][k][iA][iB] =
					dataInitialNode[RIx][k][iA][iB]
					* dInitialDeltaZ
					/ dUpdateDeltaZ;
			}
/*
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
*/
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

void VerticalDynamicsFLL::SetupReferenceColumn(
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

	// Store V in State structure
	for (int k = 0; k < nRElements; k++) {
		m_dStateNode[VIx][k] = dataInitialNode[VIx][k][iA][iB];
	}
/*
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
*/
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
/*
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
*/
/*
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
*/
}

///////////////////////////////////////////////////////////////////////////////

void VerticalDynamicsFLL::PrepareColumn(
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
/*
		// Store pressure
		if (pGrid->GetVarLocation(PIx) == DataLocation_Node) {
			m_dStateNode[PIx][k] = dX[VecFIx(FPIx, k)];
		} else {
			m_dStateREdge[PIx][k] = dX[VecFIx(FPIx, k)];
		}
*/
		// Store vertical velocity
		if (pGrid->GetVarLocation(WIx) == DataLocation_Node) {
			m_dStateNode[WIx][k] = dX[VecFIx(FWIx, k)];
		} else {
			m_dStateREdge[WIx][k] = dX[VecFIx(FWIx, k)];
		}

		// Store density
		//m_dStateNode[RIx][k] = dX[VecFIx(FRIx, k)];
	}

	if (pGrid->GetVarLocation(WIx) == DataLocation_REdge) {
		m_dStateREdge[WIx][nRElements] = dX[VecFIx(FWIx, nRElements)];
	}
/*
	if (pGrid->GetVarLocation(PIx) == DataLocation_REdge) {
		m_dStateREdge[PIx][nRElements] = dX[VecFIx(FPIx, nRElements)];
	}
*/
}

///////////////////////////////////////////////////////////////////////////////

void VerticalDynamicsFLL::BuildF(
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

	// Zero F
	memset(dF, 0, m_nColumnStateSize * sizeof(double));

	// Vertical velocity update
	double dExnerCoeff =
		phys.GetCp()
			* exp(phys.GetR() / phys.GetCv()
				* log(phys.GetR() / phys.GetP0()));

	for (int k = 1; k < m_nRElements; k++) {

		// Future altitude of this interface
		double dZint = dX[VecFIx(FZIx,k)];

		// Future layer thicknesses bounding this interface
		double dDeltaZAboveN1 = dX[VecFIx(FZIx,k+1)] - dX[VecFIx(FZIx,k  )];
		double dDeltaZBelowN1 = dX[VecFIx(FZIx,k  )] - dX[VecFIx(FZIx,k-1)];

		// Current layer thicknesses bounding this interface
		double dDeltaZAboveN0 =
			m_dColumnState[VecFIx(FZIx,k+1)] - m_dColumnState[VecFIx(FZIx,k)];
		double dDeltaZBelowN0 =
			m_dColumnState[VecFIx(FZIx,k)] - m_dColumnState[VecFIx(FZIx,k-1)];

		// Current model levels
		double dZkAbove = 0.5 * (dX[VecFIx(FZIx,k+1)] + dX[VecFIx(FZIx,k  )]);
		double dZkBelow = 0.5 * (dX[VecFIx(FZIx,k  )] + dX[VecFIx(FZIx,k-1)]);

		double dInvDeltaZk = 1.0 / (dZkAbove - dZkBelow);

		// Theta at interface
		double dThetaBelow = m_dStateNode[PIx][k-1] / m_dStateNode[RIx][k-1];
		double dThetaAbove = m_dStateNode[PIx][k  ] / m_dStateNode[RIx][k  ];

		double dThetaInt =
			  (dZint - dZkBelow) * dInvDeltaZk * dThetaAbove
			+ (dZkAbove - dZint) * dInvDeltaZk * dThetaBelow;

		// Update to vertical velocity due to buoyancy
		dF[VecFIx(FWIx, k)] +=
			phys.GetG();

		double dExnerBelow =
			exp(phys.GetR() / phys.GetCv()
				* log(m_dStateNode[PIx][k-1]
					* dDeltaZBelowN0
					/ dDeltaZBelowN1));

		double dExnerAbove =
			exp(phys.GetR() / phys.GetCv()
				* log(m_dStateNode[PIx][k]
					* dDeltaZAboveN0
					/ dDeltaZAboveN1));

		dF[VecFIx(FWIx, k)] +=
			dThetaInt * dInvDeltaZk * dExnerCoeff * (dExnerAbove - dExnerBelow);

		// Update to interface altitudes
		dF[VecFIx(FZIx, k)] +=
			dX[VecFIx(FWIx, k)];
	}

	// Construct the time-dependent component of the RHS
	double dInvDeltaT = 1.0 / m_dDeltaT;
	for (int i = 0; i < m_nColumnStateSize; i++) {
		dF[i] += (dX[i] - m_dColumnState[i]) * dInvDeltaT;
	}
}

///////////////////////////////////////////////////////////////////////////////

void VerticalDynamicsFLL::BuildJacobianF(
	const double * dX,
	double * dDG
) {
	_EXCEPTIONT("Not implemented");
}

///////////////////////////////////////////////////////////////////////////////

void VerticalDynamicsFLL::Evaluate(
	const double * dX,
	double * dF
) {
	// Prepare the column
	PrepareColumn(dX);

	// Evaluate the zero equations
	BuildF(dX, dF);
}

///////////////////////////////////////////////////////////////////////////////

void VerticalDynamicsFLL::UpdateColumnTracers(
	double dDeltaT,
	const DataArray4D<double> & dataInitialNode,
	const DataArray4D<double> & dataUpdateNode,
	const DataArray4D<double> & dataInitialREdge,
	const DataArray4D<double> & dataUpdateREdge,
	const DataArray4D<double> & dataRefTracer,
	const DataArray4D<double> & dataInitialTracer,
	DataArray4D<double> & dataUpdateTracer
) {
	_EXCEPTIONT("Not implemented");
}

///////////////////////////////////////////////////////////////////////////////

void VerticalDynamicsFLL::FilterNegativeTracers(
	int iDataUpdate
) {
	_EXCEPTIONT("Not implemented");
}

///////////////////////////////////////////////////////////////////////////////

