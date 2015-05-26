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

#include "Model.h"
#include "Grid.h"
#include "GridCSGLL.h"
#include "EquationSet.h"
#include "TimeObj.h"
#include "PolynomialInterp.h"
#include "LinearAlgebra.h"

///////////////////////////////////////////////////////////////////////////////

//#define UPWIND_HORIZONTAL_VELOCITIES
//#define UPWIND_RHO
//#define UPWIND_THERMO

//#define DETECT_CFL_VIOLATION
//#define CAP_VERTICAL_VELOCITY

///////////////////////////////////////////////////////////////////////////////

#define DEBUG

//#define THREE_COMPONENT_SOLVE

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
	m_matJacobianF.Initialize(m_nColumnStateSize, m_nColumnStateSize);

	// Initialize pivot vector
	m_vecIPiv.Initialize(m_nColumnStateSize);
#endif
#ifdef USE_DIRECTSOLVE
	// Initialize Jacobian matrix
	m_matJacobianF.Initialize(m_nColumnStateSize, m_nColumnStateSize);

	// Initialize pivot vector
	m_vecIPiv.Initialize(m_nColumnStateSize);
#endif
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

	// Allocate column for JFNK
	m_dColumnState.Initialize(m_nColumnStateSize);

	// Allocation reference column
	m_dStateRefNode.Initialize(5, nRElements);
	m_dStateRefREdge.Initialize(5, nRElements+1);

	// Solution vector from JFNK
	m_dSoln.Initialize(m_nColumnStateSize);

	// State vector at levels
	m_dStateNode.Initialize(
		m_model.GetEquationSet().GetComponents(),
		nRElements);

	// State vector at interfaces
	m_dStateREdge.Initialize(
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
	m_dStateAux.Initialize(nRElements+1);
	m_dStateAuxDiff.Initialize(nRElements+1);

	m_dXiDotNode.Initialize(nRElements);
	m_dXiDotREdge.Initialize(nRElements+1);

	m_dDiffUa.Initialize(nRElements+1);
	m_dDiffUb.Initialize(nRElements+1);

	m_dDiffPNode.Initialize(nRElements);
	m_dDiffPREdge.Initialize(nRElements+1);
	m_dDiffDiffTheta.Initialize(nRElements+1);

	m_dDiffThetaNode.Initialize(nRElements);
	m_dDiffThetaREdge.Initialize(nRElements+1);

	m_dHorizKineticEnergyNode.Initialize(nRElements);
	m_dKineticEnergyNode.Initialize(nRElements);
	m_dDiffKineticEnergyNode.Initialize(nRElements);
	m_dDiffKineticEnergyREdge.Initialize(nRElements+1);

	m_dMassFluxNode.Initialize(nRElements);
	m_dDiffMassFluxNode.Initialize(nRElements);
	m_dMassFluxREdge.Initialize(nRElements+1);
	m_dDiffMassFluxREdge.Initialize(nRElements+1);

	m_dPressureFluxNode.Initialize(nRElements);
	m_dDiffPressureFluxNode.Initialize(nRElements);
	m_dPressureFluxREdge.Initialize(nRElements+1);
	m_dDiffPressureFluxREdge.Initialize(nRElements+1);

	m_dExnerNode.Initialize(nRElements);
	m_dExnerREdge.Initialize(nRElements+1);

	m_dTracerDensityNode.Initialize(nRElements);
	m_dTracerDensityREdge.Initialize(nRElements+1);

	m_dInitialDensityNode.Initialize(nRElements);
	m_dInitialDensityREdge.Initialize(nRElements+1);

	m_dUpdateDensityNode.Initialize(nRElements);
	m_dUpdateDensityREdge.Initialize(nRElements+1);
/*
	m_dExnerNode.Initialize(nRElements);
	m_dExnerRefNode.Initialize(nRElements);

	m_dDiffExnerNode.Initialize(nRElements);
	m_dDiffExnerRefNode.Initialize(nRElements);

	m_dExnerREdge.Initialize(nRElements+1);
	m_dExnerRefREdge.Initialize(nRElements+1);

	m_dDiffExnerREdge.Initialize(nRElements+1);
	m_dDiffExnerRefREdge.Initialize(nRElements+1);
*/
/*
	// Vertical elemental grid spacing
	double dElementDeltaXi =
		static_cast<double>(m_nVerticalOrder)
		/ static_cast<double>(nRElements);
*/
/*
	// Compute second differentiation coefficients from nodes to nodes
	m_dDiffDiffNodeToNode.Initialize(
		m_nVerticalOrder, m_nVerticalOrder);

	for (int n = 0; n < m_nVerticalOrder; n++) {
	for (int m = 0; m < m_nVerticalOrder; m++) {
		m_dDiffDiffNodeToNode[n][m] = 0.0;
		for (int s = 0; s < m_nVerticalOrder; s++) {
			m_dDiffDiffNodeToNode[n][m] -=
				  dDiffNodeToNode[s][n]
				* dDiffNodeToNode[s][m]
				* dW[s];
		}
		m_dDiffDiffNodeToNode[n][m] /= dW[n];
	}
	}
*/
/*
	// Scale by 1/dxi
	for (int n = 0; n < m_nVerticalOrder; n++) {
	for (int m = 0; m < m_nVerticalOrder; m++) {
		m_dDiffDiffNodeToNode[n][m] *= dElementDeltaXi;
	}
	}
*/
/*
	// Compute second differentiation coefficients from edges to edges
	m_dDiffDiffREdgeToREdge.Initialize(
		m_nVerticalOrder+1, m_nVerticalOrder+1);

	for (int n = 0; n <= m_nVerticalOrder; n++) {
	for (int m = 0; m <= m_nVerticalOrder; m++) {
		m_dDiffDiffREdgeToREdge[n][m] = 0.0;
		for (int s = 0; s <= m_nVerticalOrder; s++) {
			m_dDiffDiffREdgeToREdge[n][m] -=
				  dDiffREdgeToREdge[s][n]
				* dDiffREdgeToREdge[s][m]
				* dWL[s];
		}
		m_dDiffDiffREdgeToREdge[n][m] /= dWL[n];
	}
	}
*/
/*
	// Scale by 1/dxi
	for (int n = 0; n <= m_nVerticalOrder; n++) {
	for (int m = 0; m <= m_nVerticalOrder; m++) {
		m_dDiffDiffREdgeToREdge[n][m] *= dElementDeltaXi;
	}
	}
*/
/*
	// Compute hyperviscosity operator
	m_dHypervisREdgeToREdge.Initialize(
		nRElements+1, nRElements+1);

	// Compute second derivative operator over whole column
	for (int a = 0; a < nFiniteElements; a++) {
		for (int n = 0; n <= m_nVerticalOrder; n++) {
		for (int m = 0; m <= m_nVerticalOrder; m++) {

			int iCol = a * m_nVerticalOrder + n;
			int iRow = a * m_nVerticalOrder + m;

			m_dHypervisREdgeToREdge[iRow][iCol] +=
				m_dDiffDiffREdgeToREdge[n][m];
		}
		}
	}
*/
/*
	for (int n = 0; n < nRElements+1; n++) {
		m_dHypervisREdgeToREdge[n][0] = 0.0;
		m_dHypervisREdgeToREdge[n][nFiniteElements * m_nVerticalOrder] = 0.0;
	}
*/
/*
	for (int a = 1; a < nFiniteElements; a++) {
		for (int n = 0; n < nRElements+1; n++) {
			m_dHypervisREdgeToREdge[n][a * m_nVerticalOrder] *= 0.5;
		}
	}
*/
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

/*
	// Compute higher powers of the second derivative operator
	DataMatrix<double> dHypervisFirstBuffer;
	dHypervisFirstBuffer.Initialize(nRElements+1, nRElements+1);

	DataMatrix<double> dHypervisSecondBuffer;
	dHypervisSecondBuffer.Initialize(nRElements+1, nRElements+1);

	if (m_nHypervisOrder > 2) {
		dHypervisFirstBuffer = m_dHypervisREdgeToREdge;

		for (int h = 2; h < m_nHypervisOrder; h += 2) {
			dHypervisSecondBuffer = m_dHypervisREdgeToREdge;

			LAPACK::DGEMM(
				m_dHypervisREdgeToREdge,
				dHypervisFirstBuffer,
				dHypervisSecondBuffer,
				1.0, 0.0);
		}
	}
*/
/*
	// Calculate Exner pressure reference profile
	for (int n = 0; n < pGrid->GetActivePatchCount(); n++) {
		GridPatch * pPatch = pGrid->GetActivePatch(n);

		const PatchBox & box = pPatch->GetPatchBox();

		// Data
		const GridData4D & dataRefNode =
			pPatch->GetReferenceState(DataLocation_Node);

		const GridData4D & dataRefREdge =
			pPatch->GetReferenceState(DataLocation_REdge);

		GridData3D & dataExnerNode =
			pPatch->GetVerticalDynamicsAuxData(0, DataLocation_Node);

		GridData3D & dataDiffExnerNode =
			pPatch->GetVerticalDynamicsAuxData(1, DataLocation_Node);

		GridData3D & dataExnerREdge =
			pPatch->GetVerticalDynamicsAuxData(0, DataLocation_REdge);

		GridData3D & dataDiffExnerREdge =
			pPatch->GetVerticalDynamicsAuxData(1, DataLocation_REdge);

		// Loop over all nodes
		for (int i = box.GetAInteriorBegin(); i < box.GetAInteriorEnd(); i++) {
		for (int j = box.GetBInteriorBegin(); j < box.GetBInteriorEnd(); j++) {

			// Initialize reference Exner pressure at model levels
			for (int k = 0; k < nRElements; k++) {
				dataExnerNode[k][i][j] =
					phys.ExnerPressureFromRhoTheta(
						  dataRefNode[PIx][k][i][j]
						* dataRefNode[RIx][k][i][j]);

				m_dExnerRefNode[k] = dataExnerNode[k][i][j];
			}

			// Initialize reference Exner pressure at model interfaces
			for (int k = 0; k <= nRElements; k++) {
				dataExnerREdge[k][i][j] =
					phys.ExnerPressureFromRhoTheta(
						  dataRefREdge[PIx][k][i][j]
						* dataRefREdge[RIx][k][i][j]);

				m_dExnerRefREdge[k] = dataExnerREdge[k][i][j];
			}

			// Differentiate the reference Exner pressure
			if ((pGrid->GetVarLocation(RIx) == DataLocation_Node) &&
				(pGrid->GetVarLocation(WIx) == DataLocation_Node)
			) {
				pGrid->DifferentiateNodeToNode(
					m_dExnerRefNode,
					m_dDiffExnerRefNode);

				pGrid->DifferentiateREdgeToREdge(
					m_dExnerRefREdge,
					m_dDiffExnerRefREdge);

			} else if (
				(pGrid->GetVarLocation(RIx) == DataLocation_Node) &&
				(pGrid->GetVarLocation(WIx) == DataLocation_REdge)
			) {
				if (m_fExnerPressureOnLevels) {
					pGrid->DifferentiateNodeToNode(
						m_dExnerRefNode,
						m_dDiffExnerRefNode);

					pGrid->DifferentiateNodeToREdge(
						m_dExnerRefNode,
						m_dDiffExnerRefREdge);

				} else {
					pGrid->DifferentiateREdgeToNode(
						m_dExnerRefREdge,
						m_dDiffExnerRefNode);

					pGrid->DifferentiateREdgeToREdge(
						m_dExnerRefREdge,
						m_dDiffExnerRefREdge);
				}

			} else {
				_EXCEPTIONT("Invalid variable staggering"
					" (possibly UNIMPLEMENTED)");
			}

			// Store derivatives at this point
			for (int k = 0; k < nRElements; k++) {
				dataDiffExnerNode[k][i][j] = m_dDiffExnerRefNode[k];
			}

			for (int k = 0; k <= nRElements; k++) {
				dataDiffExnerREdge[k][i][j] = m_dDiffExnerRefREdge[k];
			}
		}
		}
	}
*/
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

		// Data
		const GridData4D & dataRefNode =
			pPatch->GetReferenceState(DataLocation_Node);

		const GridData4D & dataInitialNode =
			pPatch->GetDataState(iDataInitial, DataLocation_Node);

		GridData4D & dataUpdateNode =
			pPatch->GetDataState(iDataUpdate, DataLocation_Node);

		const GridData4D & dataRefREdge =
			pPatch->GetReferenceState(DataLocation_REdge);

		const GridData4D & dataInitialREdge =
			pPatch->GetDataState(iDataInitial, DataLocation_REdge);

		GridData4D & dataUpdateREdge =
			pPatch->GetDataState(iDataUpdate, DataLocation_REdge);

		GridData4D & dataInitialTracer =
			pPatch->GetDataTracers(iDataInitial);

		GridData4D & dataUpdateTracer =
			pPatch->GetDataTracers(iDataUpdate);

		// Auxiliary data storing Exner pressure
		const GridData3D & dataExnerNode =
			pPatch->GetVerticalDynamicsAuxData(0, DataLocation_Node);

		const GridData3D & dataExnerREdge =
			pPatch->GetVerticalDynamicsAuxData(0, DataLocation_REdge);

		const GridData3D & dataDiffExnerNode =
			pPatch->GetVerticalDynamicsAuxData(1, DataLocation_Node);

		const GridData3D & dataDiffExnerREdge =
			pPatch->GetVerticalDynamicsAuxData(1, DataLocation_REdge);

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
/*
#ifndef THREE_COMPONENT_SOLVE
				// Apply update to U
				if (pGrid->GetVarLocation(UIx) == DataLocation_REdge) {
					for (int k = 0; k <= pGrid->GetRElements(); k++) {
						dataUpdateREdge[UIx][k][iA][iB] -=
							dDeltaT * m_dSoln[VecFIx(FUIx, k)];
					}
				} else {
					for (int k = 0; k < pGrid->GetRElements(); k++) {
						dataUpdateNode[UIx][k][iA][iB] -=
							dDeltaT * m_dSoln[VecFIx(FUIx, k)];
					}
				}

				// Apply update to V
				if (pGrid->GetVarLocation(VIx) == DataLocation_REdge) {
					for (int k = 0; k <= pGrid->GetRElements(); k++) {
						dataUpdateREdge[VIx][k][iA][iB] -=
							dDeltaT * m_dSoln[VecFIx(FVIx, k)];
					}
				} else {
					for (int k = 0; k < pGrid->GetRElements(); k++) {
						dataUpdateNode[VIx][k][iA][iB] -=
							dDeltaT * m_dSoln[VecFIx(FVIx, k)];
					}
				}
#endif
*/
				// Apply update to P
				if (pGrid->GetVarLocation(PIx) == DataLocation_REdge) {
					for (int k = 0; k <= pGrid->GetRElements(); k++) {
						dataUpdateREdge[PIx][k][iA][iB] -=
							dDeltaT * m_dSoln[VecFIx(FPIx, k)];
					}
				} else {
					for (int k = 0; k < pGrid->GetRElements(); k++) {
						dataUpdateNode[PIx][k][iA][iB] -=
							dDeltaT * m_dSoln[VecFIx(FPIx, k)];
					}
				}

				// Apply update to W
				if (pGrid->GetVarLocation(WIx) == DataLocation_REdge) {
					for (int k = 0; k <= pGrid->GetRElements(); k++) {
						dataUpdateREdge[WIx][k][iA][iB] -=
							dDeltaT * m_dSoln[VecFIx(FWIx, k)];
					}

				} else {
					for (int k = 0; k < pGrid->GetRElements(); k++) {
						dataUpdateNode[WIx][k][iA][iB] -=
							dDeltaT * m_dSoln[VecFIx(FWIx, k)];
					}
				}

				// Apply update to Rho
				if (pGrid->GetVarLocation(RIx) == DataLocation_REdge) {
					for (int k = 0; k <= pGrid->GetRElements(); k++) {
						dataUpdateREdge[RIx][k][iA][iB] -=
							dDeltaT * m_dSoln[VecFIx(FRIx, k)];
					}
				} else {
					for (int k = 0; k < pGrid->GetRElements(); k++) {
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
					const DataMatrix4D<double> & dContraMetricXi =
						pPatch->GetContraMetricXi();
					const DataMatrix4D<double> & dDerivRNode =
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
		}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void VerticalDynamicsFEM::BootstrapJacobian() {

	static const double Epsilon = 1.0e-5;

	int nDim = m_dColumnState.GetRows();

	DataMatrix<double> dJacobian;
	dJacobian.Initialize(nDim, nDim);

	DataVector<double> dJC;
	dJC.Initialize(nDim);

	DataVector<double> dG;
	dG.Initialize(nDim);

	DataVector<double> dJCref;
	dJCref.Initialize(nDim);

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
		const DataMatrix4D<double> & dContraMetricA =
			pPatch->GetContraMetricA();
		const DataMatrix4D<double> & dContraMetricB =
			pPatch->GetContraMetricB();

		// State Data
		const GridData4D & dataRefNode =
			pPatch->GetReferenceState(DataLocation_Node);

		const GridData4D & dataInitialNode =
			pPatch->GetDataState(iDataInitial, DataLocation_Node);

		GridData4D & dataUpdateNode =
			pPatch->GetDataState(iDataUpdate, DataLocation_Node);

		const GridData4D & dataRefREdge =
			pPatch->GetReferenceState(DataLocation_REdge);

		const GridData4D & dataInitialREdge =
			pPatch->GetDataState(iDataInitial, DataLocation_REdge);

		GridData4D & dataUpdateREdge =
			pPatch->GetDataState(iDataUpdate, DataLocation_REdge);

		// Tracer Data
		GridData4D & dataInitialTracer =
			pPatch->GetDataTracers(iDataInitial);

		GridData4D & dataUpdateTracer =
			pPatch->GetDataTracers(iDataUpdate);

		// Number of tracers
		const int nTracerCount = dataInitialTracer.GetComponents();

		// Auxiliary data storing Exner pressure
		const GridData3D & dataExnerNode =
			pPatch->GetVerticalDynamicsAuxData(0, DataLocation_Node);

		const GridData3D & dataExnerREdge =
			pPatch->GetVerticalDynamicsAuxData(0, DataLocation_REdge);

		const GridData3D & dataDiffExnerNode =
			pPatch->GetVerticalDynamicsAuxData(1, DataLocation_Node);

		const GridData3D & dataDiffExnerREdge =
			pPatch->GetVerticalDynamicsAuxData(1, DataLocation_REdge);

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
				DataVector<double> dEval;
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
                DataVector<double> dEval;
                dEval.Initialize(m_dColumnState.GetRows());
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

			DataVector<double> dJC;
			dJC.Initialize(m_dColumnState.GetRows());

			DataVector<double> dG;
			dG.Initialize(m_dColumnState.GetRows());

			DataVector<double> dJCref;
			dJCref.Initialize(m_dColumnState.GetRows());

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
				DataVector<double> dEval;
				dEval.Initialize(m_dColumnState.GetRows());
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
/*
#ifndef THREE_COMPONENT_SOLVE
			// Apply updated state to U
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

			// Apply updated state to V
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
#endif
*/
			// Apply updated state to theta
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
}


///////////////////////////////////////////////////////////////////////////////

void VerticalDynamicsFEM::SetupReferenceColumn(
	GridPatch * pPatch,
	int iA,
	int iB,
	const GridData4D & dataRefNode,
	const GridData4D & dataInitialNode,
	const GridData4D & dataRefREdge,
	const GridData4D & dataInitialREdge
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

/*
#ifndef THREE_COMPONENT_SOLVE
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
#endif
*/
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

	// Metric terms
	const DataMatrix3D<double> & dJacobian =
		m_pPatch->GetJacobian();
	const DataMatrix4D<double> & dDerivRNode =
		m_pPatch->GetDerivRNode();
	const DataMatrix4D<double> & dDerivRREdge =
		m_pPatch->GetDerivRREdge();
	const DataMatrix4D<double> & dContraMetricXi =
		m_pPatch->GetContraMetricXi();
	const DataMatrix4D<double> & dContraMetricXiREdge =
		m_pPatch->GetContraMetricXiREdge();

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

#ifdef UPWIND_THERMO
			// Second derivatives of theta on model levels
			pGrid->DiffDiffNodeToNode(
				m_dStateNode[PIx],
				m_dDiffDiffTheta);
#endif

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

#ifdef UPWIND_THERMO
			// Second derivatives of theta on interfaces
			pGrid->DiffDiffREdgeToREdge(
				m_dStateREdge[PIx],
				m_dDiffDiffTheta);
#endif
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
			m_dStateNode[WIx][k] * dDerivRNode[k][m_iA][m_iB][2];

		m_dXiDotNode[k] =
			  dContraMetricXi[k][m_iA][m_iB][0] * m_dStateNode[UIx][k]
			+ dContraMetricXi[k][m_iA][m_iB][1] * m_dStateNode[VIx][k]
			+ dContraMetricXi[k][m_iA][m_iB][2] * dCovUx;
	}

	// Calculate u^xi on model interfaces
	if ((pGrid->GetVarLocation(WIx) == DataLocation_REdge) ||
	    (pGrid->GetVarLocation(PIx) == DataLocation_REdge)
	) {
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

	// Strictly enforce conservation
	if (pGrid->GetVarLocation(WIx) == DataLocation_Node) {
		m_dXiDotNode[0] = 0.0;
		m_dXiDotNode[nRElements-1] = 0.0;
	}
/*
#ifdef UPWIND_THETA
	// Compute second derivatives of theta
	DiffDiffREdgeToREdge(
		m_dStateREdge[PIx],
		m_dDiffDiffTheta
	);

	// Compute higher derivatives of theta used for hyperdiffusion
	if (m_nHypervisOrder > 0) {

		if (pGrid->GetVarLocation(PIx) != DataLocation_REdge) {
			_EXCEPTIONT("Not implemented");
		}

		for (int h = 2; h < m_nHypervisOrder; h += 2) {
			memcpy(
				m_dStateAux,
				m_dDiffDiffTheta,
				(nRElements + 1) * sizeof(double));

			DiffDiffREdgeToREdge(
				m_dStateAux,
				m_dDiffDiffTheta
			);
		}
	}
#endif
*/
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

	// Metric terms
	const DataMatrix3D<double> & dJacobian =
		m_pPatch->GetJacobian();
	const DataMatrix3D<double> & dJacobianREdge =
		m_pPatch->GetJacobianREdge();
	const DataMatrix4D<double> & dDerivRNode =
		m_pPatch->GetDerivRNode();
	const DataMatrix4D<double> & dDerivRREdge =
		m_pPatch->GetDerivRREdge();
	const DataMatrix4D<double> & dContraMetricA =
		m_pPatch->GetContraMetricA();
	const DataMatrix4D<double> & dContraMetricB =
		m_pPatch->GetContraMetricB();
	const DataMatrix4D<double> & dContraMetricXi =
		m_pPatch->GetContraMetricXi();
	const DataMatrix4D<double> & dContraMetricAREdge =
		m_pPatch->GetContraMetricAREdge();
	const DataMatrix4D<double> & dContraMetricBREdge =
		m_pPatch->GetContraMetricBREdge();
	const DataMatrix4D<double> & dContraMetricXiREdge =
		m_pPatch->GetContraMetricXiREdge();
	const DataMatrix3D<double> & dElementArea =
		m_pPatch->GetElementArea();

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
				dJacobianREdge[k][m_iA][m_iB]
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
				dJacobian[k][m_iA][m_iB]
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
		dSum += dElementArea[k][m_iA][m_iB] * m_dDiffMassFluxNode[k];

		dF[VecFIx(FRIx, k)] =
			m_dDiffMassFluxNode[k]
			/ dJacobian[k][m_iA][m_iB];
	}

#ifdef FORMULATION_PRESSURE
	// Pressure flux calculated on model interfaces
	if (!fMassFluxOnLevels) {
		for (int k = 1; k < nRElements; k++) {
			m_dPressureFluxREdge[k] =
				dJacobianREdge[k][m_iA][m_iB]
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
				dJacobian[k][m_iA][m_iB]
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
			/ dJacobian[k][m_iA][m_iB];
	}
#endif
#if defined(FORMULATION_RHOTHETA_PI) || defined(FORMULATION_RHOTHETA_P)
	// RhoTheta flux on model interfaces
	if (!fMassFluxOnLevels) {
		for (int k = 1; k < nRElements; k++) {
			m_dPressureFluxREdge[k] =
				dJacobianREdge[k][m_iA][m_iB]
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
				dJacobian[k][m_iA][m_iB]
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
			/ dJacobian[k][m_iA][m_iB];
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

#ifdef UPWIND_THERMO
		double dDeltaXi = 1.0 / static_cast<double>(nRElements);
		for (int k = 0; k < nRElements; k++) {
			dF[VecFIx(FPIx, k)] -=
				0.5 * fabs(m_dXiDotNode[k]) * dDeltaXi * m_dDiffDiffTheta[k];
		}
#endif

	// Update theta on model interfaces
	} else {
		// Change in Theta on model interfaces
		for (int k = 0; k <= nRElements; k++) {
			dF[VecFIx(FPIx, k)] +=
				m_dXiDotREdge[k] * m_dDiffThetaREdge[k];
		}

#ifdef UPWIND_THERMO
		double dDeltaXi = 1.0 / static_cast<double>(nRElements);
		for (int k = 0; k <= nRElements; k++) {
			dF[VecFIx(FPIx, k)] -=
				0.5 * fabs(m_dXiDotREdge[k]) * dDeltaXi * m_dDiffDiffTheta[k];
		}
#endif

	}
#endif
#ifdef FORMULATION_THETA_FLUX
	// Update theta on model levels
	if (pGrid->GetVarLocation(PIx) == DataLocation_Node) {
		// Pressure flux on model levels
		for (int k = 0; k < nRElements; k++) {
			m_dPressureFluxNode[k] =
				dJacobian[k][m_iA][m_iB]
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
				m_dDiffPressureFluxNode[k] / dJacobian[k][m_iA][m_iB];

			dF[VecFIx(FPIx, k)] -=
				m_dStateNode[PIx][k] * m_dStateAuxDiff[k];
		}

	// Update theta on model interfaces
	} else {
		// Pressure flux on model levels
		for (int k = 0; k <= nRElements; k++) {
			m_dPressureFluxREdge[k] =
				dJacobianREdge[k][m_iA][m_iB]
				* m_dStateREdge[PIx][k]
				* m_dXiDotREdge[k];

			m_dStateAux[k] =
				dJacobianREdge[k][m_iA][m_iB]
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
			   	/ dJacobianREdge[k][m_iA][m_iB];
		}
	}
#endif

	// Kinetic energy on model levels
	for (int k = 0; k < nRElements; k++) {
		double dCovUa = m_dStateNode[UIx][k];
		double dCovUb = m_dStateNode[VIx][k];
		double dCovUx = m_dStateNode[WIx][k] * dDerivRNode[k][m_iA][m_iB][2];

		double dConUa =
			  dContraMetricA[k][m_iA][m_iB][0] * dCovUa
			+ dContraMetricA[k][m_iA][m_iB][1] * dCovUb
			+ dContraMetricA[k][m_iA][m_iB][2] * dCovUx;

		double dConUb =
			  dContraMetricB[k][m_iA][m_iB][0] * dCovUa
			+ dContraMetricB[k][m_iA][m_iB][1] * dCovUb
			+ dContraMetricB[k][m_iA][m_iB][2] * dCovUx;

		double dConUx =
			  dContraMetricXi[k][m_iA][m_iB][0] * dCovUa
			+ dContraMetricXi[k][m_iA][m_iB][1] * dCovUb
			+ dContraMetricXi[k][m_iA][m_iB][2] * dCovUx;

		// Specific kinetic energy
		m_dKineticEnergyNode[k] =
			  0.5 * (dConUa * dCovUa + dConUb * dCovUb + dConUx * dCovUx);
	}

	// Update equation for vertical velocity on levels
	if (pGrid->GetVarLocation(WIx) == DataLocation_Node) {
		pGrid->DifferentiateNodeToNode(
			m_dKineticEnergyNode,
			m_dDiffKineticEnergyNode);

		for (int k = 1; k < nRElements-1; k++) {

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

			double dCovUa = m_dStateNode[UIx][k];
			double dCovUb = m_dStateNode[VIx][k];
			double dCovUx = m_dStateNode[WIx][k] * dDerivRNode[k][m_iA][m_iB][2];

			double dConUa =
				  dContraMetricA[k][m_iA][m_iB][0] * dCovUa
				+ dContraMetricA[k][m_iA][m_iB][1] * dCovUb
				+ dContraMetricA[k][m_iA][m_iB][2] * dCovUx;

			double dConUb =
				  dContraMetricB[k][m_iA][m_iB][0] * dCovUa
				+ dContraMetricB[k][m_iA][m_iB][1] * dCovUb
				+ dContraMetricB[k][m_iA][m_iB][2] * dCovUx;

			double dCurlTerm =
				- dConUa * m_dDiffUa[k]
				- dConUb * m_dDiffUb[k];

			dF[VecFIx(FWIx, k)] =
				(m_dDiffKineticEnergyNode[k]
				 	+ dCurlTerm
					+ dPressureGradientForce)
				/ dDerivRNode[k][m_iA][m_iB][2];

			dF[VecFIx(FWIx, k)] +=
				phys.GetG();
		}

	// Update equation for vertical velocity on interfaces
	} else {
		pGrid->DifferentiateNodeToREdge(
			m_dKineticEnergyNode,
			m_dDiffKineticEnergyREdge);

		pGrid->DifferentiateREdgeToREdge(
			m_dStateREdge[WIx],
			m_dStateAuxDiff);

		for (int k = 1; k < nRElements; k++) {
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

			double dCovUa = m_dStateREdge[UIx][k];
			double dCovUb = m_dStateREdge[VIx][k];
			double dCovUx = m_dStateREdge[WIx][k] * dDerivRREdge[k][m_iA][m_iB][2];

			double dConUa =
				  dContraMetricAREdge[k][m_iA][m_iB][0] * dCovUa
				+ dContraMetricAREdge[k][m_iA][m_iB][1] * dCovUb
				+ dContraMetricAREdge[k][m_iA][m_iB][2] * dCovUx;

			double dConUb =
				  dContraMetricBREdge[k][m_iA][m_iB][0] * dCovUa
				+ dContraMetricBREdge[k][m_iA][m_iB][1] * dCovUb
				+ dContraMetricBREdge[k][m_iA][m_iB][2] * dCovUx;

			double dCurlTerm =
				- dConUa * m_dDiffUa[k]
				- dConUb * m_dDiffUb[k];

			dF[VecFIx(FWIx, k)] =
				(m_dDiffKineticEnergyREdge[k]
					+ dCurlTerm
					+ dPressureGradientForce)
				/ dDerivRREdge[k][m_iA][m_iB][2];

			dF[VecFIx(FWIx, k)] +=
				phys.GetG();
		}

		//_EXCEPTION();
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
	for (int i = 0; i < m_nColumnStateSize; i++) {
		dF[i] += (dX[i] - m_dColumnState[i]) / m_dDeltaT;
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

	} else if (
		pGrid->GetVerticalStaggering() ==
			Grid::VerticalStaggering_Interfaces
	) {
		dF[VecFIx(FWIx, 0)] =
			  dContraMetricXi[0][m_iA][m_iB][0] * m_dStateNode[UIx][0]
			+ dContraMetricXi[0][m_iA][m_iB][1] * m_dStateNode[VIx][0]
			+ dContraMetricXi[0][m_iA][m_iB][2]
				* dDerivRNode[0][m_iA][m_iB][2] * m_dStateNode[WIx][0];

		int k = nRElements-1;

		dF[VecFIx(FWIx, k)] =
			  dContraMetricXi[k][m_iA][m_iB][0] * m_dStateNode[UIx][k]
			+ dContraMetricXi[k][m_iA][m_iB][1] * m_dStateNode[VIx][k]
			+ dContraMetricXi[k][m_iA][m_iB][2]
				* dDerivRNode[k][m_iA][m_iB][2] * m_dStateNode[WIx][k];

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

	const DataMatrix<double> & dInterpNodeToREdge =
		opInterpNodeToREdge.GetCoeffs();
	const DataMatrix<double> & dInterpREdgeToNode =
		opInterpREdgeToNode.GetCoeffs();
	const DataMatrix<double> & dDiffNodeToNode =
		opDiffNodeToNode.GetCoeffs();
	const DataMatrix<double> & dDiffNodeToREdge =
		opDiffNodeToREdge.GetCoeffs();
	const DataMatrix<double> & dDiffREdgeToNode =
		opDiffREdgeToNode.GetCoeffs();
	const DataMatrix<double> & dDiffREdgeToREdge =
		opDiffREdgeToREdge.GetCoeffs();

	const DataVector<int> & iInterpNodeToREdgeBegin =
		opInterpNodeToREdge.GetIxBegin();
	const DataVector<int> & iInterpREdgeToNodeBegin =
		opInterpREdgeToNode.GetIxBegin();
	const DataVector<int> & iDiffNodeToNodeBegin =
		opDiffNodeToNode.GetIxBegin();
	const DataVector<int> & iDiffNodeToREdgeBegin =
		opDiffNodeToREdge.GetIxBegin();
	const DataVector<int> & iDiffREdgeToNodeBegin =
		opDiffREdgeToNode.GetIxBegin();
	const DataVector<int> & iDiffREdgeToREdgeBegin =
		opDiffREdgeToREdge.GetIxBegin();

	const DataVector<int> & iInterpNodeToREdgeEnd =
		opInterpNodeToREdge.GetIxEnd();
	const DataVector<int> & iInterpREdgeToNodeEnd =
		opInterpREdgeToNode.GetIxEnd();
	const DataVector<int> & iDiffNodeToNodeEnd =
		opDiffNodeToNode.GetIxEnd();
	const DataVector<int> & iDiffNodeToREdgeEnd =
		opDiffNodeToREdge.GetIxEnd();
	const DataVector<int> & iDiffREdgeToNodeEnd =
		opDiffREdgeToNode.GetIxEnd();
	const DataVector<int> & iDiffREdgeToREdgeEnd =
		opDiffREdgeToREdge.GetIxEnd();

	// Metric components
	const DataMatrix3D<double> & dJacobianNode =
		m_pPatch->GetJacobian();
	const DataMatrix3D<double> & dJacobianREdge =
		m_pPatch->GetJacobianREdge();
	const DataMatrix4D<double> & dContraMetricXi =
		m_pPatch->GetContraMetricXi();
	const DataMatrix4D<double> & dDerivRNode =
		m_pPatch->GetDerivRNode();
	const DataMatrix4D<double> & dContraMetricXiREdge =
		m_pPatch->GetContraMetricXiREdge();
	const DataMatrix4D<double> & dDerivRREdge =
		m_pPatch->GetDerivRREdge();

	// Physical constants
	const PhysicalConstants & phys = m_model.GetPhysicalConstants();

	// Number of radial elements
	const int nRElements = pGrid->GetRElements();

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

//////////////////////////////////////////////
// Prognostic thermodynamic variable pressure
#ifdef FORMULATION_PRESSURE

	// Vertical velocity on interfaces (CPH or LOR staggering)
	if (pGrid->GetVarLocation(WIx) == DataLocation_REdge) {
		_EXCEPTIONT("Not implemented");

	// Vertical velocity on nodes (LEV or INT staggering)
	} else {

		// dP_k/dP_n
		for (int k = 0; k < nRElements; k++) {

			int n = iDiffNodeToNodeBegin[k];
			for (; n < iDiffNodeToNodeEnd[k]; n++) {

				// Pressure flux
				dDG[MatFIx(FPIx, n, FPIx, k)] +=
					phys.GetGamma()
					* dJacobianNode[n][m_iA][m_iB]
					/ dJacobianNode[k][m_iA][m_iB]
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
					/ dJacobianNode[k][m_iA][m_iB]
					* dJacobianNode[n][m_iA][m_iB]
					* phys.GetGamma()
					* m_dStateNode[PIx][n]
					* dContraMetricXi[n][m_iA][m_iB][2]
					* dDerivRNode[n][m_iA][m_iB][2];
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
				* dContraMetricXi[k][m_iA][m_iB][2]
				* dDerivRNode[k][m_iA][m_iB][2]
				* m_dDiffPNode[k];
		}

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
					/ dDerivRNode[k][m_iA][m_iB][2];
			}

			dDG[MatFIx(FRIx, k, FWIx, k)] +=
				- m_dDiffPNode[k]
				/ dDerivRNode[k][m_iA][m_iB][2]
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

		// dT_k/dW_l
		for (int k = 0; k < nRElements; k++) {
			int l = iInterpREdgeToNodeBegin[k];
			for (; l < iInterpREdgeToNodeEnd[k]; l++) {
				dDG[MatFIx(FWIx, l, FPIx, k)] +=
					m_dDiffThetaNode[k]
					* dInterpREdgeToNode[k][l]
					* dContraMetricXi[k][m_iA][m_iB][2]
					* dDerivRNode[k][m_iA][m_iB][2];
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

		// dW_k/dT_m and dW_k/dR_m
		for (int k = 1; k < nRElements; k++) {

			double dRHSWCoeff = 
				1.0 / dDerivRNode[k][m_iA][m_iB][2]
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
					 1.0 / dDerivRREdge[k][m_iA][m_iB][2]
					 * dInterpNodeToREdge[k][l]
					 * m_dDiffPREdge[k];
			}
		}

	// Charney-Phillips staggering
	} else {
		if (pGrid->GetVarLocation(WIx) != DataLocation_REdge) {
			_EXCEPTIONT("Not implemented");
		}

		// dT_k/dW_k
		for (int k = 1; k < nRElements; k++) {
			dDG[MatFIx(FWIx, k, FPIx, k)] +=
				m_dDiffThetaREdge[k]
				* dContraMetricXiREdge[k][m_iA][m_iB][2]
				* dDerivRREdge[k][m_iA][m_iB][2];
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

		// dW_k/dT_l and dW_k/dR_m
		for (int k = 1; k < nRElements; k++) {

			double dRHSWCoeff = 
				1.0 / dDerivRNode[k][m_iA][m_iB][2]
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
				 1.0 / dDerivRREdge[k][m_iA][m_iB][2]
				 * m_dDiffPREdge[k];
		}
	}
/*
#ifdef UPWIND_THETA
	if (pGrid->GetVarLocation(PIx) == DataLocation_REdge) {
		for (int k = 0; k <= nRElements; k++) {

			// dT_k/dW_k
			double dSignW = -1.0;

			if (m_dXiDotREdge[k] < 0.0) {
				dSignW = -1.0;
			} else if (m_dXiDotREdge[k] > 0.0) {
				dSignW = 1.0;
			}

			dDG[MatFIx(FWIx, k, FPIx, k)] -=
				m_dHypervisCoeff
				* dOrthonomREdge[k][m_iA][m_iB][2]
				* dSignW
				* m_dDiffDiffTheta[k];

			// dT_k/dT_m
			for (int m = 0; m <= nRElements; m++) {
				dDG[MatFIx(FPIx, m, FPIx, k)] -=
					m_dHypervisCoeff
					* fabs(m_dXiDotREdge[k])
					* m_dHypervisREdgeToREdge[m][k];
			}
		}
	}
#endif
*/
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
						/ dJacobianREdge[m][m_iA][m_iB]
						* dJacobianNode[k][m_iA][m_iB]
						* m_dStateREdge[RIx][m]
						* dContraMetricXiREdge[m][m_iA][m_iB][2]
						* dDerivRREdge[m][m_iA][m_iB][2];
				}

				int n = iInterpNodeToREdgeBegin[m];
				for (; n < iInterpNodeToREdgeEnd[m]; n++) {

					dDG[MatFIx(FRIx, n, FRIx, k)] +=
						dDiffREdgeToNode[k][m]
						* dJacobianREdge[m][m_iA][m_iB]
						/ dJacobianNode[k][m_iA][m_iB]
						* dInterpNodeToREdge[m][n]
						* m_dXiDotREdge[m];
				}
			}
		}

		// dW_k/dW_m
		for (int k = 1; k < nRElements; k++) {
			int l = iDiffNodeToREdgeBegin[k];
			for (; l < iDiffNodeToREdgeEnd[k]; l++) {

				int m = iInterpREdgeToNodeBegin[l];
				for (; m < iInterpREdgeToNodeEnd[l]; m++) {
					dDG[MatFIx(FWIx, m, FWIx, k)] +=
						dInterpREdgeToNode[l][m]
						* dDiffNodeToREdge[k][l]
						/ dDerivRREdge[k][m_iA][m_iB][2]
						* dDerivRNode[l][m_iA][m_iB][2]
						* m_dXiDotNode[l];
				}
			}

/*
#ifdef UPWIND_RHO
			// Diffusion of Rho (dRho_k/dW)
			if (pGrid->GetVarLocation(RIx) == DataLocation_Node) {
				int lBeginNext = lBegin + m_nVerticalOrder;

				double dSignWL = (m_dXiDotREdge[lBegin] > 0.0)?(1.0):(-1.0);
				double dSignWR = (m_dXiDotREdge[lBeginNext] > 0.0)?(1.0):(-1.0);

				dDG[MaxFIx(FWIx, lBegin, FRIx, k)] +=
					dDiffReconsPolyNode[l]

			} else {
				_EXCEPTIONT("Upwind Rho on interfaces: Unimplemented");
			}
#endif
*/
		}

	// Vertical velocity on nodes (LEV or INT staggering)
	} else {

		for (int k = 0; k < nRElements; k++) {

			int n = iDiffNodeToNodeBegin[k];
			for (; n < iDiffNodeToNodeEnd[k]; n++) {

				// dRho_k/dRho_n
				dDG[MatFIx(FRIx, n, FRIx, k)] +=
					dDiffNodeToNode[k][n]
					* dJacobianNode[n][m_iA][m_iB]
					/ dJacobianNode[k][m_iA][m_iB]
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
					* dJacobianNode[n][m_iA][m_iB]
					/ dJacobianNode[k][m_iA][m_iB]
					* m_dStateNode[RIx][n]
					* dDerivRNode[n][m_iA][m_iB][2]
					* dContraMetricXi[n][m_iA][m_iB][2];
			}

			// Boundary conditions
			if (pGrid->GetVerticalStaggering() ==
			    Grid::VerticalStaggering_Interfaces
			) {
				if ((k == 0) || (k == nRElements-1)) {
					continue;
				}
			}

			// dW_k/dW_n
			n = iDiffNodeToNodeBegin[k];
			for (; n < iDiffNodeToNodeEnd[k]; n++) {
				dDG[MatFIx(FWIx, n, FWIx, k)] +=
					dDiffNodeToNode[k][n]
					/ dDerivRNode[k][m_iA][m_iB][2]
					* dDerivRNode[n][m_iA][m_iB][2]
					* m_dXiDotNode[n];

			}
		}
	}

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
			dDerivRNode[0][m_iA][m_iB][2]
			* dContraMetricXi[0][m_iA][m_iB][2];
		dDG[MatFIx(FWIx, nRElements-1, FWIx, nRElements-1)] =
			dDerivRNode[nRElements-1][m_iA][m_iB][2]
			* dContraMetricXi[nRElements-1][m_iA][m_iB][2];
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
	const GridData4D & dataInitialNode,
	const GridData4D & dataUpdateNode,
	const GridData4D & dataInitialREdge,
	const GridData4D & dataUpdateREdge,
	const GridData4D & dataInitialTracer,
	const GridData4D & dataUpdateTracer
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
	const int nComponents = dataInitialTracer.GetComponents();

	// Metric quantities
	const DataMatrix4D<double> & dContraMetricXi =
		m_pPatch->GetContraMetricXi();
	const DataMatrix4D<double> & dContraMetricXiREdge =
		m_pPatch->GetContraMetricXiREdge();
	const DataMatrix3D<double> & dElementArea =
		m_pPatch->GetElementArea();
	const DataMatrix3D<double> & dJacobianNode =
		m_pPatch->GetJacobian();
	const DataMatrix3D<double> & dJacobianREdge =
		m_pPatch->GetJacobianREdge();
	const DataMatrix4D<double> & dDerivRNode =
		m_pPatch->GetDerivRNode();
	const DataMatrix4D<double> & dDerivRREdge =
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

	// Calculate u^xi on model levels
	if (fMassFluxOnLevels) {
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

		for (int k = 0; k < nRElements; k++) {
			double dCovUx =
				m_dStateNode[WIx][k] * dDerivRNode[k][m_iA][m_iB][2];

			m_dXiDotNode[k] =
				  dContraMetricXi[k][m_iA][m_iB][0] * m_dStateNode[UIx][k]
				+ dContraMetricXi[k][m_iA][m_iB][1] * m_dStateNode[VIx][k]
				+ dContraMetricXi[k][m_iA][m_iB][2] * dCovUx;
		}

	// Calculate u^xi on model interfaces
	} else {
		if (pGrid->GetVarLocation(WIx) == DataLocation_REdge) {
			for (int k = 1; k < nRElements; k++) {
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

	// Strictly enforce conservation
	if (pGrid->GetVarLocation(WIx) == DataLocation_Node) {
		m_dXiDotNode[0] = 0.0;
		m_dXiDotNode[nRElements-1] = 0.0;
	}

	// Loop through all tracer species
	for (int c = 0; c < nComponents; c++) {

		// Calculate mass flux on nodes
		if (fMassFluxOnLevels) {
			for (int k = 0; k < nRElements; k++) {
				m_dMassFluxNode[k] =
					dJacobianNode[k][m_iA][m_iB]
					* dataInitialTracer[c][k][m_iA][m_iB]
					* dataUpdateNode[RIx][k][m_iA][m_iB]
					/ dataInitialNode[RIx][k][m_iA][m_iB]
					* m_dXiDotNode[k];
			}

			pGrid->DifferentiateNodeToNode(
				m_dMassFluxNode,
				m_dDiffMassFluxNode,
				fZeroBoundaries);

		// Calculate mass flux on interfaces
		} else {

			// Interpolate initial density to interfaces
			for (int k = 0; k < nRElements; k++) {
				m_dInitialDensityNode[k] = dataInitialNode[RIx][k][m_iA][m_iB];
			}

			pGrid->InterpolateNodeToREdge(
				m_dInitialDensityNode,
				m_dInitialDensityREdge);

			// Interpolate updated density to interfaces
			for (int k = 0; k < nRElements; k++) {
				m_dUpdateDensityNode[k] = dataUpdateNode[RIx][k][m_iA][m_iB];
			}

			pGrid->InterpolateNodeToREdge(
				m_dUpdateDensityNode,
				m_dUpdateDensityREdge);

			// Interpolate tracer density to interfaces
			for (int k = 0; k < nRElements; k++) {
				m_dTracerDensityNode[k] = dataInitialTracer[c][k][m_iA][m_iB];
			}

			pGrid->InterpolateNodeToREdge(
				m_dTracerDensityNode,
				m_dTracerDensityREdge);

			// Calculate mass flux
			for (int k = 0; k < nRElements; k++) {
				m_dMassFluxREdge[k] =
					dJacobianREdge[k][m_iA][m_iB]
					* m_dTracerDensityREdge[k]
					* dataUpdateREdge[RIx][k][m_iA][m_iB]
					/ dataInitialREdge[RIx][k][m_iA][m_iB]
					* m_dXiDotREdge[k];
			}

			pGrid->DifferentiateREdgeToNode(
				m_dMassFluxREdge,
				m_dDiffMassFluxNode);
		}

		// Update tracers
		for (int k = 0; k < nRElements; k++) {
			dataUpdateTracer[c][k][m_iA][m_iB] -=
				dDeltaT
				* m_dDiffMassFluxNode[k]
				/ dJacobianNode[k][m_iA][m_iB];
		}

#pragma message "Apply limiter only to finite element?"
		// Apply positive definite limiter to the tracers
		double dTotalMass = 0.0;
		double dNonNegativeMass = 0.0;
		for (int k = 0; k < nRElements; k++) {
			double dPointwiseMass =
				  dataUpdateTracer[c][k][m_iA][m_iB]
				* dElementArea[k][m_iA][m_iB];

			dTotalMass += dPointwiseMass;

			if (dataUpdateTracer[c][k][m_iA][m_iB] >= 0.0) {
				dNonNegativeMass += dPointwiseMass;
			}
		}

		// Apply scaling ratio to points with non-negative mass
		double dR = dTotalMass / dNonNegativeMass;

		for (int k = 0; k < nRElements; k++) {
			if (dataUpdateTracer[c][k][m_iA][m_iB] > 0.0) {
				dataUpdateTracer[c][k][m_iA][m_iB] *= dR;
			} else {
				dataUpdateTracer[c][k][m_iA][m_iB] = 0.0;
			}
		}
	}
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

