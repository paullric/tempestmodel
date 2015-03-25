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
#include "GaussQuadrature.h"
#include "GaussLobattoQuadrature.h"

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
//#define UPWIND_THETA

//#define DETECT_CFL_VIOLATION
//#define CAP_VERTICAL_VELOCITY

///////////////////////////////////////////////////////////////////////////////

VerticalDynamicsFEM::VerticalDynamicsFEM(
	Model & model,
	int nHorizontalOrder,
	int nVerticalOrder,
	int nHypervisOrder,
	bool fFullyExplicit,
	bool fUseReferenceState,
	bool fExnerPressureOnLevels,
	bool fMassFluxOnLevels
) :
	VerticalDynamics(model),
	m_nHorizontalOrder(nHorizontalOrder),
	m_nVerticalOrder(nVerticalOrder),
	m_fFullyExplicit(fFullyExplicit),
	m_fUseReferenceState(fUseReferenceState),
	m_fExnerPressureOnLevels(fExnerPressureOnLevels),
	m_fMassFluxOnLevels(fMassFluxOnLevels),
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

	// Get the vertical interpolation and differentiation coefficients
	const DataMatrix<double> & dDiffREdgeToNode =
		pGrid->GetDiffREdgeToNode();
	const DataMatrix<double> & dDiffNodeToREdge =
		pGrid->GetDiffNodeToREdge();
	const DataMatrix<double> & dDiffNodeToNode =
		pGrid->GetDiffNodeToNode();
	const DataMatrix<double> & dDiffREdgeToREdge =
		pGrid->GetDiffREdgeToREdge();

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
	InitializeJFNK(m_nColumnStateSize, m_nColumnStateSize);
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
	if (pGrid->GetReconstructionPolyType() != 2) {
		_EXCEPTIONT("Diagonal Jacobian only implemented for "
			"ReconstructionPolyType == 2");
	}
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

	// Get points for Gaussian quadrature
	DataVector<double> dG;
	DataVector<double> dGL;
	DataVector<double> dW;
	DataVector<double> dWL;

	// Reference element [0,1] model levels
	GaussQuadrature::GetPoints(m_nVerticalOrder, 0.0, 1.0, dG, dW);

	// Reference element [0,1] model interfaces
	GaussLobattoQuadrature::GetPoints(m_nVerticalOrder+1, 0.0, 1.0, dGL, dWL);

	// State vector at levels
	m_dStateNode.Initialize(
		m_model.GetEquationSet().GetComponents(),
		nRElements);

	// State vector at interfaces
	m_dStateREdge.Initialize(
		m_model.GetEquationSet().GetComponents(),
		nRElements+1);

	// Auxiliary variables at interfaces
	int nFiniteElements = nRElements / m_nVerticalOrder;
	if (nRElements % m_nVerticalOrder != 0) {
		_EXCEPTIONT("Logic error: Vertical order must divide RElements");
	}

	m_dStateAux.Initialize(nRElements+1);
	m_dStateAuxDiff.Initialize(nRElements+1);

	m_dXiDotNode.Initialize(nRElements);
	m_dXiDotREdge.Initialize(nRElements+1);

	m_dDiffP.Initialize(nRElements+1);
	m_dDiffDiffP.Initialize(nRElements+1);

	m_dHorizKineticEnergyNode.Initialize(nRElements);
	m_dKineticEnergyNode.Initialize(nRElements);
	m_dDiffKineticEnergyNode.Initialize(nRElements);

	m_dMassFluxNode.Initialize(nRElements);
	m_dDiffMassFluxNode.Initialize(nRElements);
	m_dMassFluxREdge.Initialize(nRElements+1);
	m_dDiffMassFluxREdge.Initialize(nRElements+1);

	m_dPressureFluxNode.Initialize(nRElements);
	m_dDiffPressureFluxNode.Initialize(nRElements);
	m_dPressureFluxREdge.Initialize(nRElements+1);
	m_dDiffPressureFluxREdge.Initialize(nRElements+1);

	m_dExnerPertNode.Initialize(nRElements);
	m_dExnerRefNode.Initialize(nRElements);

	m_dDiffExnerPertNode.Initialize(nRElements);
	m_dDiffExnerRefNode.Initialize(nRElements);

	m_dExnerPertREdge.Initialize(nRElements+1);
	m_dExnerRefREdge.Initialize(nRElements+1);

	m_dDiffExnerPertREdge.Initialize(nRElements+1);
	m_dDiffExnerRefREdge.Initialize(nRElements+1);

	// Vertical elemental grid spacing
	double dElementDeltaXi =
		static_cast<double>(m_nVerticalOrder)
		/ static_cast<double>(nRElements);

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
/*
	// Scale by 1/dxi
	for (int n = 0; n < m_nVerticalOrder; n++) {
	for (int m = 0; m < m_nVerticalOrder; m++) {
		m_dDiffDiffNodeToNode[n][m] *= dElementDeltaXi;
	}
	}
*/
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
/*
	// Scale by 1/dxi
	for (int n = 0; n <= m_nVerticalOrder; n++) {
	for (int m = 0; m <= m_nVerticalOrder; m++) {
		m_dDiffDiffREdgeToREdge[n][m] *= dElementDeltaXi;
	}
	}
*/
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
/*
	for (int n = 0; n < nRElements+1; n++) {
		m_dHypervisREdgeToREdge[n][0] = 0.0;
		m_dHypervisREdgeToREdge[n][nFiniteElements * m_nVerticalOrder] = 0.0;
	}
*/
	for (int a = 1; a < nFiniteElements; a++) {
		for (int n = 0; n < nRElements+1; n++) {
			m_dHypervisREdgeToREdge[n][a * m_nVerticalOrder] *= 0.5;
		}
	}

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

	// Allocate data for column Jacobian
	m_dJacobian3DNode.Initialize(nRElements);
	m_dJacobian3DREdge.Initialize(nRElements+1);

	m_dDxRNode.Initialize(nRElements);
	m_dDxRREdge.Initialize(nRElements+1);
}

///////////////////////////////////////////////////////////////////////////////

void VerticalDynamicsFEM::DiffDiffREdgeToREdge(
	const double * dDataREdge,
	double * dDiffDiffREdge
) {
	// Number of radial elements
	int nRElements = m_model.GetGrid()->GetRElements();

	int nFiniteElements = nRElements / m_nVerticalOrder;

	// Zero the data
	memset(dDiffDiffREdge, 0, (nRElements+1) * sizeof(double));

	// Apply all interfaces values to all interfaces within element
	for (int a = 0; a < nFiniteElements; a++) {
	for (int l = 0; l <= m_nVerticalOrder; l++) {

		int lBegin = a * m_nVerticalOrder;

		for (int m = 0; m <= m_nVerticalOrder; m++) {
			dDiffDiffREdge[lBegin + l] +=
				  m_dDiffDiffREdgeToREdge[l][m]
				* dDataREdge[lBegin + m];
		}
	}
	}

	// Halve interior element interface values
	for (int a = 1; a < nFiniteElements; a++) {
		dDiffDiffREdge[a * m_nVerticalOrder] *= 0.5;
	}
}

///////////////////////////////////////////////////////////////////////////////

void VerticalDynamicsFEM::SetupReferenceColumn(
	GridPatch * pPatch,
	int iA,
	int iB,
	const DataMatrix<double> & dataJacobian2D,
	const DataMatrix3D<double> & dataJacobian3DNode,
	const DataMatrix3D<double> & dataJacobian3DREdge,
	const GridData4D & dataRefNode,
	const GridData4D & dataInitialNode,
	const GridData4D & dataRefREdge,
	const GridData4D & dataInitialREdge,
	const GridData3D & dataExnerNode,
	const GridData3D & dataDiffExnerNode,
	const GridData3D & dataExnerREdge,
	const GridData3D & dataDiffExnerREdge
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

	// Store the column Jacobian
	for (int k = 0; k < nRElements; k++) {
		m_dJacobian3DNode[k] = dataJacobian3DNode[k][iA][iB];
		m_dDxRNode[k] = m_dJacobian3DNode[k] / dataJacobian2D[iA][iB];
	}
	for (int k = 0; k <= nRElements; k++) {
		m_dJacobian3DREdge[k] = dataJacobian3DREdge[k][iA][iB];
		m_dDxRREdge[k] = m_dJacobian3DREdge[k] / dataJacobian2D[iA][iB];
	}

	// Store U in State structure
	if (pGrid->GetVarLocation(UIx) == DataLocation_Node) {
		for (int k = 0; k < nRElements; k++) {
			m_dStateNode[UIx][k] = dataInitialNode[UIx][k][iA][iB];
		}

		pGrid->InterpolateNodeToREdge(
			m_dStateNode[UIx],
			m_dStateRefNode[UIx],
			m_dStateREdge[UIx],
			m_dStateRefREdge[UIx]);

	} else {
		for (int k = 0; k <= nRElements; k++) {
			m_dStateREdge[UIx][k] = dataInitialREdge[UIx][k][iA][iB];
		}

		pGrid->InterpolateREdgeToNode(
			m_dStateREdge[UIx],
			m_dStateRefREdge[UIx],
			m_dStateNode[UIx],
			m_dStateRefNode[UIx]);
	}

	// Store V in State structure
	if (pGrid->GetVarLocation(VIx) == DataLocation_Node) {
		for (int k = 0; k < nRElements; k++) {
			m_dStateNode[VIx][k] = dataInitialNode[VIx][k][iA][iB];
		}

		pGrid->InterpolateNodeToREdge(
			m_dStateNode[VIx],
			m_dStateRefNode[VIx],
			m_dStateREdge[VIx],
			m_dStateRefREdge[VIx]);

	} else {
		for (int k = 0; k <= nRElements; k++) {
			m_dStateREdge[VIx][k] = dataInitialREdge[VIx][k][iA][iB];
		}

		pGrid->InterpolateREdgeToNode(
			m_dStateREdge[VIx],
			m_dStateRefREdge[VIx],
			m_dStateNode[VIx],
			m_dStateRefNode[VIx]);
	}
/*
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

		// Auxiliary data storing Exner pressure
		const GridData3D & dataExnerNode =
			pPatch->GetVerticalDynamicsAuxData(0, DataLocation_Node);

		const GridData3D & dataExnerREdge =
			pPatch->GetVerticalDynamicsAuxData(0, DataLocation_REdge);

		const GridData3D & dataDiffExnerNode =
			pPatch->GetVerticalDynamicsAuxData(1, DataLocation_Node);

		const GridData3D & dataDiffExnerREdge =
			pPatch->GetVerticalDynamicsAuxData(1, DataLocation_REdge);

		// Pointwise topography
		const DataMatrix<double> & dataJacobian2D = pPatch->GetJacobian2D();

		// Pointwise Jacobian on nodes and interfaces
		const DataMatrix3D<double> & dataJacobian3DNode =
			pPatch->GetJacobian();
		const DataMatrix3D<double> & dataJacobian3DREdge =
			pPatch->GetJacobianREdge();

		// Loop over all nodes
		for (int i = box.GetAInteriorBegin(); i < box.GetAInteriorEnd(); i++) {
		for (int j = box.GetBInteriorBegin(); j < box.GetBInteriorEnd(); j++) {

			// Update thermodynamic variables
			if (m_fFullyExplicit) {

				int iA = i;
				int iB = j;

				SetupReferenceColumn(
					pPatch, iA, iB,
					dataJacobian2D,
					dataJacobian3DNode,
					dataJacobian3DREdge,
					dataRefNode,
					dataInitialNode,
					dataRefREdge,
					dataInitialREdge,
					dataExnerNode,
					dataDiffExnerNode,
					dataExnerREdge,
					dataDiffExnerREdge);

				Evaluate(m_dColumnState, m_dSoln);
/*
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

	BuildJacobianF2(m_dSoln, &(dJacobian[0][0]));

	fp = fopen("DG2.txt", "w");
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

		// Auxiliary data storing Exner pressure
		const GridData3D & dataExnerNode =
			pPatch->GetVerticalDynamicsAuxData(0, DataLocation_Node);

		const GridData3D & dataExnerREdge =
			pPatch->GetVerticalDynamicsAuxData(0, DataLocation_REdge);

		const GridData3D & dataDiffExnerNode =
			pPatch->GetVerticalDynamicsAuxData(1, DataLocation_Node);

		const GridData3D & dataDiffExnerREdge =
			pPatch->GetVerticalDynamicsAuxData(1, DataLocation_REdge);

		// Pointwise topography
		const DataMatrix<double> & dataJacobian2D = pPatch->GetJacobian2D();

		// Contravariant metric components
		const DataMatrix3D<double> & dataJacobian3DNode =
			pPatch->GetJacobian();
		const DataMatrix3D<double> & dataJacobian3DREdge =
			pPatch->GetJacobianREdge();

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
				dataJacobian2D,
				dataJacobian3DNode,
				dataJacobian3DREdge,
				dataRefNode,
				dataInitialNode,
				dataRefREdge,
				dataInitialREdge,
				dataExnerNode,
				dataDiffExnerNode,
				dataExnerREdge,
				dataDiffExnerREdge);

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
			BuildJacobianF2(m_dColumnState, &(m_matJacobianF[0][0]));

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

void VerticalDynamicsFEM::PrepareColumn(
	const double * dX
) {
	// Indices of EquationSet variables
	const int UIx = 0;
	const int VIx = 1;
	const int PIx = 2;
	const int WIx = 3;
	const int RIx = 4;

	// Finite element grid spacing
	const GridGLL * pGrid = dynamic_cast<const GridGLL *>(m_model.GetGrid());

	// Physical constants
	const PhysicalConstants & phys = m_model.GetPhysicalConstants();

	// Number of radial elements
	const int nRElements = pGrid->GetRElements();
/*
	// U and V on model levels
	for (int k = 0; k < nRElements; k++) {
		m_dStateNode[UIx][k] = dX[VecFIx(FUIx, k)];
		m_dStateNode[VIx][k] = dX[VecFIx(FVIx, k)];
	}
*/
	// W on model interfaces
	if (pGrid->GetVarLocation(WIx) == DataLocation_REdge) {
		for (int k = 0; k <= nRElements; k++) {
			m_dStateREdge[WIx][k] = dX[VecFIx(FWIx, k)];
		}

		// W is needed on model levels
		pGrid->InterpolateREdgeToNode(
			m_dStateREdge[WIx],
			m_dStateRefREdge[WIx],
			m_dStateNode[WIx],
			m_dStateRefNode[WIx]);

	// W on model levels
	} else {
		for (int k = 0; k < nRElements; k++) {
			m_dStateNode[WIx][k] = dX[VecFIx(FWIx, k)];
		}

		// W is needed on model interfaces
		pGrid->InterpolateNodeToREdge(
			m_dStateNode[WIx],
			m_dStateRefNode[WIx],
			m_dStateREdge[WIx],
			m_dStateRefREdge[WIx]);
	}

	// Pressure on model interfaces
	if (pGrid->GetVarLocation(PIx) == DataLocation_REdge) {
		for (int k = 0; k <= nRElements; k++) {
			m_dStateREdge[PIx][k] = dX[VecFIx(FPIx, k)];
		}

		// Calculate update for T
		pGrid->DifferentiateREdgeToREdge(
			m_dStateREdge[PIx],
			m_dDiffP);

		// T is needed on model levels
		if ((m_fExnerPressureOnLevels) ||
			(pGrid->GetVarLocation(WIx) == DataLocation_Node)
		) {
			pGrid->InterpolateREdgeToNode(
				m_dStateREdge[PIx],
				m_dStateRefREdge[PIx],
				m_dStateNode[PIx],
				m_dStateRefNode[PIx]);
		}

	// Pressure on model levels
	} else {
		for (int k = 0; k < nRElements; k++) {
			m_dStateNode[PIx][k] = dX[VecFIx(FPIx, k)];
		}

		// Calculate update for T
		pGrid->DifferentiateNodeToNode(
			m_dStateNode[PIx],
			m_dDiffP);

		// Pressure is needed on model interfaces
		if ((!m_fExnerPressureOnLevels) ||
			(pGrid->GetVarLocation(WIx) == DataLocation_REdge)
		) {
			pGrid->InterpolateNodeToREdge(
				m_dStateNode[PIx],
				m_dStateRefNode[PIx],
				m_dStateREdge[PIx],
				m_dStateRefREdge[PIx]);
		}
	}

	// Rho on model interfaces
	if (pGrid->GetVarLocation(RIx) == DataLocation_REdge) {
		for (int k = 0; k <= nRElements; k++) {
			m_dStateREdge[RIx][k] = dX[VecFIx(FRIx, k)];
		}

		// Rho is needed on model levels
		if ((m_fMassFluxOnLevels) ||
			(m_fExnerPressureOnLevels)
		) {
			pGrid->InterpolateREdgeToNode(
				m_dStateREdge[RIx],
				m_dStateRefREdge[RIx],
				m_dStateNode[RIx],
				m_dStateRefNode[RIx]);
		}

	// Rho on model levels
	} else {
		for (int k = 0; k < nRElements; k++) {
			m_dStateNode[RIx][k] = dX[VecFIx(FRIx, k)];
		}

		// Rho is needed on model interfaces
		if ((!m_fMassFluxOnLevels) ||
			(!m_fExnerPressureOnLevels)
		) {
			pGrid->InterpolateNodeToREdge(
				m_dStateNode[RIx],
				m_dStateRefNode[RIx],
				m_dStateREdge[RIx],
				m_dStateRefREdge[RIx]);
		}
	}

#ifdef UPWIND_THETA
	// Compute second derivatives of theta
	DiffDiffREdgeToREdge(
		m_dStateREdge[PIx],
		m_dDiffDiffP
	);

	// Compute higher derivatives of theta used for hyperdiffusion
	if (m_nHypervisOrder > 0) {

		if (pGrid->GetVarLocation(PIx) != DataLocation_REdge) {
			_EXCEPTIONT("Not implemented");
		}

		for (int h = 2; h < m_nHypervisOrder; h += 2) {
			memcpy(
				m_dStateAux,
				m_dDiffDiffP,
				(nRElements + 1) * sizeof(double));

			DiffDiffREdgeToREdge(
				m_dStateAux,
				m_dDiffDiffP
			);
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

	// Finite element grid spacing
	const GridGLL * pGrid = dynamic_cast<const GridGLL *>(m_model.GetGrid());

	// Physical constants
	const PhysicalConstants & phys = m_model.GetPhysicalConstants();

	// Number of radial elements
	const int nRElements = pGrid->GetRElements();

	// Metric terms
	const DataMatrix4D<double> & dDerivRNode =
		m_pPatch->GetDerivRNode();

	// Zero F
	memset(dF, 0, m_nColumnStateSize * sizeof(double));

	// Compute Xi dot on model levels
	for (int k = 0; k < nRElements; k++) {
		const DataMatrix4D<double> & dContraMetricXi =
			m_pPatch->GetContraMetricXi();

		double dCovUx =
			m_dStateNode[WIx][k]* dDerivRNode[k][m_iA][m_iB][2];

		m_dXiDotNode[k] =
			  dContraMetricXi[k][m_iA][m_iB][0] * m_dStateNode[UIx][k]
			+ dContraMetricXi[k][m_iA][m_iB][1] * m_dStateNode[VIx][k]
			+ dContraMetricXi[k][m_iA][m_iB][2] * dCovUx;
	}

	// Check conservation
	if (fabs(m_dXiDotNode[0]) > 1.0e-14) {
		_EXCEPTIONT("Bottom boundary: Conservation error");
	}
	if (fabs(m_dXiDotNode[nRElements-1]) > 1.0e-14) {
		_EXCEPTIONT("Top boundary: Conservation error");
	}

	// Mass flux on model levels
	for (int k = 0; k < nRElements; k++) {
		m_dMassFluxNode[k] =
			m_dJacobian3DNode[k]
			* m_dStateNode[RIx][k]
			* m_dXiDotNode[k];
	}

	pGrid->DifferentiateNodeToNode(
		m_dMassFluxNode,
		m_dDiffMassFluxNode);

	// Change in density on model levels
	for (int k = 0; k < nRElements; k++) {
		dF[VecFIx(FRIx, k)] =
			m_dDiffMassFluxNode[k]
			/ m_dJacobian3DNode[k];
	}

	// Pressure flux on model levels
	for (int k = 0; k < nRElements; k++) {
		m_dPressureFluxNode[k] =
			m_dJacobian3DNode[k]
			* phys.GetGamma()
			* m_dStateNode[PIx][k]
			* m_dXiDotNode[k];
	}

	pGrid->DifferentiateNodeToNode(
		m_dPressureFluxNode,
		m_dDiffPressureFluxNode);

	// Change in pressure on model levels
	for (int k = 0; k < nRElements; k++) {
		dF[VecFIx(FPIx, k)] =
			- (phys.GetGamma() - 1.0)
			* m_dXiDotNode[k]
			* m_dDiffP[k];

		dF[VecFIx(FPIx, k)] +=
			m_dDiffPressureFluxNode[k]
			/ m_dJacobian3DNode[k];
	}

	// Kinetic energy on model levels
	for (int k = 0; k < nRElements; k++) {
		const DataMatrix4D<double> & dContraMetricA =
			m_pPatch->GetContraMetricA();
		const DataMatrix4D<double> & dContraMetricB =
			m_pPatch->GetContraMetricB();
		const DataMatrix4D<double> & dContraMetricXi =
			m_pPatch->GetContraMetricXi();

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

	pGrid->DifferentiateNodeToNode(
		m_dKineticEnergyNode,
		m_dDiffKineticEnergyNode);
/*
	// Change in alpha/beta velocity on model levels
	for (int k = 0; k < nRElements; k++) {
		const DataMatrix4D<double> & dContraMetricA =
			m_pPatch->GetContraMetricA();
		const DataMatrix4D<double> & dContraMetricB =
			m_pPatch->GetContraMetricB();

		double dDerivUpdate =
			  m_dDiffP[k] / m_dStateNode[RIx][k]
			+ m_dDiffKineticEnergyNode[k]
			+ phys.GetG() * dDerivRNode[k][m_iA][m_iB][2];

		if (k == 0) {
			dDerivUpdate = 0.0;
		}

		dF[VecFIx(FUIx, k)] =
			dContraMetricA[k][m_iA][m_iB][2] * dDerivUpdate;

		dF[VecFIx(FVIx, k)] =
			dContraMetricB[k][m_iA][m_iB][2] * dDerivUpdate;
	}
*/
	// Change in vertical velocity on model levels
	for (int k = 0; k < nRElements; k++) {

		if (k == 0) {
			dF[VecFIx(FWIx, k)] = 0.0;

		} else if (k == nRElements-1) {
			dF[VecFIx(FWIx, k)] = 0.0;

		} else {
			dF[VecFIx(FWIx, k)] =
				(m_dDiffKineticEnergyNode[k]
					+ m_dDiffP[k] / m_dStateNode[RIx][k])
				/ m_dDxRNode[k];

			dF[VecFIx(FWIx, k)] +=
				phys.GetG();
		}
/*
		if ((m_iA == 35) && (m_iB == 1)) {
			printf("%i %1.5e %1.5e\n",
				k,
				m_dDiffKineticEnergyNode[k] / m_dDxRNode[k],
				dF[VecFIx(FWIx, k)]);
		}
*/
	}

	// Construct the time-dependent component of the RHS
	for (int i = 0; i < m_nColumnStateSize; i++) {
		dF[i] += (dX[i] - m_dColumnState[i]) / m_dDeltaT;
	}
/*
	// Apply boundary conditions to W
	if (pGrid->GetVarLocation(WIx) == DataLocation_REdge) {
		dF[VecFIx(FWIx, 0)] =
			m_pPatch->CalculateXiDotREdge(
				0, m_iA, m_iB,
				m_dStateREdge[UIx][0],
				m_dStateREdge[VIx][0],
				m_dStateREdge[WIx][0]);

		dF[VecFIx(FWIx, nRElements)] =
			m_pPatch->CalculateXiDotREdge(
				nRElements, m_iA, m_iB,
				m_dStateREdge[UIx][nRElements],
				m_dStateREdge[VIx][nRElements],
				m_dStateREdge[WIx][nRElements]);

	} else if (
		pGrid->GetVerticalStaggering() ==
			Grid::VerticalStaggering_Interfaces
	) {
		dF[VecFIx(FWIx, 0)] =
			m_pPatch->CalculateXiDotREdge(
				0, m_iA, m_iB,
				m_dStateREdge[UIx][0],
				m_dStateREdge[VIx][0],
				m_dStateREdge[WIx][0]);

		dF[VecFIx(FWIx, nRElements-1)] =
			m_pPatch->CalculateXiDotREdge(
				nRElements, m_iA, m_iB,
				m_dStateREdge[UIx][nRElements-1],
				m_dStateREdge[VIx][nRElements-1],
				m_dStateREdge[WIx][nRElements-1]);
	} else {
		_EXCEPTIONT("UNIMPLEMENTED");
	}
*/
/*
#ifdef APPLY_RAYLEIGH_WITH_VERTICALDYN

	// Rayleigh friction strength
	const GridData3D & dataRayleighStrengthNode =
		m_pPatch->GetRayleighStrength(DataLocation_Node);

	const GridData3D & dataRayleighStrengthREdge =
		m_pPatch->GetRayleighStrength(DataLocation_REdge);

	// Rayleigh damping on nodes
	for (int k = 0; k < pGrid->GetRElements(); k++) {
		double dNu = dataRayleighStrengthNode[k][m_iA][m_iB];

		if (dNu == 0.0) {
			continue;
		}

		if (pGrid->GetVarLocation(PIx) == DataLocation_Node) {
			dF[VecFIx(FPIx, k)] =
				dNu * (m_dStateNode[VIx][k] - m_dStateRefNode[VIx][k]);
		}
		if (pGrid->GetVarLocation(WIx) == DataLocation_Node) {
			dF[VecFIx(FWIx, k)] =
				dNu * (m_dStateNode[WIx][k] - m_dStateRefNode[WIx][k]);
		}
		if (pGrid->GetVarLocation(RIx) == DataLocation_Node) {
			dF[VecFIx(FRIx, k)] =
				dNu * (m_dStateNode[RIx][k] - m_dStateRefNode[RIx][k]);
		}
	}

	// Rayleigh damping on interfaces
	for (int k = 0; k <= pGrid->GetRElements(); k++) {
		double dNu = dataRayleighStrengthREdge[k][m_iA][m_iB];

		if (dNu == 0.0) {
			continue;
		}

		if (pGrid->GetVarLocation(PIx) == DataLocation_REdge) {
			dF[VecFIx(FPIx, k)] =
				dNu * (m_dStateREdge[VIx][k] - m_dStateRefREdge[VIx][k]);
		}
		if (pGrid->GetVarLocation(WIx) == DataLocation_REdge) {
			dF[VecFIx(FWIx, k)] =
				dNu * (m_dStateREdge[WIx][k] - m_dStateRefREdge[WIx][k]);
		}
		if (pGrid->GetVarLocation(RIx) == DataLocation_REdge) {
			dF[VecFIx(FRIx, k)] =
				dNu * (m_dStateREdge[RIx][k] - m_dStateRefREdge[RIx][k]);
		}
	}
#endif
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

	if ((pGrid->GetVarLocation(UIx) != DataLocation_Node) ||
		(pGrid->GetVarLocation(VIx) != DataLocation_Node) ||
		(pGrid->GetVarLocation(PIx) != DataLocation_REdge) ||
		(pGrid->GetVarLocation(WIx) != DataLocation_REdge) ||
		(pGrid->GetVarLocation(RIx) != DataLocation_Node) ||
		(m_fMassFluxOnLevels) ||
		(!m_fExnerPressureOnLevels)
	) {
		_EXCEPTIONT("Not implemented");
	}

	// Get the vertical interpolation and differentiation coefficients
	const DataMatrix<double> & dInterpNodeToREdge =
		pGrid->GetOpInterpNodeToREdge().GetCoeffs();
	const DataMatrix<double> & dInterpREdgeToNode =
		pGrid->GetOpInterpREdgeToNode().GetCoeffs();

	const DataMatrix<double> & dDiffREdgeToNode =
		pGrid->GetDiffREdgeToNode();
	const DataMatrix<double> & dDiffNodeToREdge =
		pGrid->GetDiffNodeToREdge();
	const DataMatrix<double> & dDiffNodeToNode =
		pGrid->GetDiffNodeToNode();
	const DataMatrix<double> & dDiffREdgeToREdge =
		pGrid->GetDiffREdgeToREdge();
	const DataMatrix<double> & dDiffNodeToREdgeAmal =
		pGrid->GetDiffNodeToREdgeAmal();
	const DataMatrix<double> & dDiffNodeToREdgeLeft =
		pGrid->GetDiffNodeToREdgeLeft();
	const DataMatrix<double> & dDiffNodeToREdgeRight =
		pGrid->GetDiffNodeToREdgeRight();
	const DataVector<double> & dDiffReconsPolyREdge =
		pGrid->GetDiffReconsPolyREdge();
	const DataVector<double> & dDiffReconsPolyNode =
		pGrid->GetDiffReconsPolyNode();

	// Orthonormalization coefficients
	const DataMatrix4D<double> & dOrthonomREdge =
		m_pPatch->GetOrthonormREdge();

	// Physical constants
	const PhysicalConstants & phys = m_model.GetPhysicalConstants();

	// Number of radial elements
	const int nRElements = pGrid->GetRElements();

	// Zero DG
	memset(dDG, 0,
		m_nColumnStateSize * m_nColumnStateSize * sizeof(double));

	// dT_k/dT
	int nFiniteElements = nRElements / m_nVerticalOrder;
	for (int a = 0; a < nFiniteElements; a++) {
	for (int l = 0; l <= m_nVerticalOrder; l++) {
		int lBegin = a * m_nVerticalOrder;

		for (int k = 0; k <= m_nVerticalOrder; k++) {
			dDG[MatFIx(FPIx, lBegin+l, FPIx, lBegin+k)] +=
				dDiffREdgeToREdge[k][l]
				* m_dXiDotREdge[lBegin+k];
		}
	}
	}

	for (int a = 1; a < nFiniteElements; a++) {
		int k = a * m_nVerticalOrder;
		int lBegin = k - m_nVerticalOrder;
		int lEnd   = k + m_nVerticalOrder + 1;
		for (int l = lBegin; l < lEnd; l++) {
			dDG[MatFIx(FPIx, l, FPIx, k)] *= 0.5;
		}
	}

	// dT_k/dW
	for (int k = 0; k <= nRElements; k++) {
		dDG[MatFIx(FWIx, k, FPIx, k)] =
			m_dDiffP[k] * dOrthonomREdge[k][m_iA][m_iB][2];
	}

	if ((nFiniteElements == 1) || (nFiniteElements == 2)) {
		_EXCEPTIONT("UNIMPLEMENTED: At least three elements needed");
	}

	int kBegin;
	int kEnd;

	// dW_k/dT_l and dW_k/dR_m (bottom element)
	for (int k = 1; k < m_nVerticalOrder; k++) {
		double dRHSWCoeff = 
			1.0 / m_dDxRREdge[k]
			* m_dStateREdge[PIx][k]
			* phys.GetR() / phys.GetCv();

		for (int m = 0; m < 2 * m_nVerticalOrder; m++) {
			double dTEntry =
				dRHSWCoeff 
				* dDiffNodeToREdgeLeft[k][m]
				* (m_dExnerPertNode[m] + m_dExnerRefNode[m])
					/ m_dStateNode[PIx][m];

			int mx = m % m_nVerticalOrder;

			for (int l = 0; l <= m_nVerticalOrder; l++) {
				dDG[MatFIx(FPIx, l, FWIx, k)] +=
					dTEntry * dInterpREdgeToNode[mx][l];
			}

			dDG[MatFIx(FRIx, m, FWIx, k)] +=
				dRHSWCoeff
				* dDiffNodeToREdgeLeft[k][m]
				* (m_dExnerPertNode[m] + m_dExnerRefNode[m])
					/ m_dStateNode[RIx][m];
		}
	}

	// dW_k/dT_l and dW_k/dR_m (middle elements)
	int kLast = (nFiniteElements-1) * m_nVerticalOrder;
	for (int k = m_nVerticalOrder; k <= kLast; k++) {
		int kx = k % m_nVerticalOrder;
		int a = k / m_nVerticalOrder;
		if (k == kLast) {
			a--;
			kx = m_nVerticalOrder;
		}

		int lPrev = (a-1) * m_nVerticalOrder;
		int lBegin = lPrev + m_nVerticalOrder;

		double dRHSWCoeff = 
			1.0 / m_dDxRREdge[k]
			* m_dStateREdge[PIx][k]
			* phys.GetR() / phys.GetCv();

		for (int m = 0; m < 3 * m_nVerticalOrder; m++) {

			int mx = m % m_nVerticalOrder;
			int ma = (m / m_nVerticalOrder) * m_nVerticalOrder;

			double dTEntry =
				dRHSWCoeff 
				* dDiffNodeToREdgeAmal[kx][m]
				* (m_dExnerPertNode[lPrev + m] + m_dExnerRefNode[lPrev + m])
					/ m_dStateNode[PIx][lPrev + m];

			for (int l = 0; l <= m_nVerticalOrder; l++) {
				dDG[MatFIx(FPIx, lPrev+ma+l, FWIx, k)] +=
					dTEntry * dInterpREdgeToNode[lPrev+m][lPrev+ma+l];
			}

			dDG[MatFIx(FRIx, lPrev+m, FWIx, k)] +=
				dRHSWCoeff
				* dDiffNodeToREdgeAmal[kx][m]
				* (m_dExnerPertNode[lPrev + m] + m_dExnerRefNode[lPrev + m])
					/ m_dStateNode[RIx][lPrev + m];
		}
	}

	// dW_k/dT_l and dW_k/dR_m (top element)
	kBegin = (nFiniteElements-1) * m_nVerticalOrder;
	kEnd = kBegin + m_nVerticalOrder;
	for (int k = kBegin + 1; k < kEnd; k++) {
		int kx = k - kBegin;
		int lPrev = kBegin - m_nVerticalOrder;

		double dRHSWCoeff = 
			1.0 / m_dDxRREdge[k]
			* m_dStateREdge[PIx][k]
			* phys.GetR() / phys.GetCv();

		for (int m = 0; m < 2 * m_nVerticalOrder; m++) {
			double dTEntry =
				dRHSWCoeff 
				* dDiffNodeToREdgeRight[kx][m]
				* (m_dExnerPertNode[lPrev + m] + m_dExnerRefNode[lPrev + m])
					/ m_dStateNode[PIx][lPrev + m];

			int mx = m % m_nVerticalOrder;
			int ma = (m / m_nVerticalOrder) * m_nVerticalOrder;

			for (int l = 0; l <= m_nVerticalOrder; l++) {
				dDG[MatFIx(FPIx, lPrev+ma+l, FWIx, k)] +=
					dTEntry * dInterpREdgeToNode[lPrev+m][lPrev+ma+l];
			}

			dDG[MatFIx(FRIx, lPrev+m, FWIx, k)] +=
				dRHSWCoeff
				* dDiffNodeToREdgeRight[kx][m]
				* (m_dExnerPertNode[lPrev + m] + m_dExnerRefNode[lPrev + m])
					/ m_dStateNode[RIx][lPrev + m];
		}
	}

	// dW_k/dT_k (first theta in RHS)
	for (int k = 1; k < nRElements; k++) {
		dDG[MatFIx(FPIx, k, FWIx, k)] +=
			 1.0 / m_dDxRREdge[k]
			 * (m_dDiffExnerRefREdge[k] + m_dDiffExnerPertREdge[k]);
	}

	// dRho
	for (int k = 0; k < nRElements; k++) {
		int a = k / m_nVerticalOrder;
		int l = k % m_nVerticalOrder;
		int lBegin = a * m_nVerticalOrder;

		// dRho_k/dW
		for (int m = 0; m <= m_nVerticalOrder; m++) {
			dDG[MatFIx(FWIx, lBegin+m, FRIx, k)] +=
				dDiffREdgeToNode[l][m]
				* m_dJacobian3DREdge[lBegin+m]
				/ m_dJacobian3DNode[k]
				* m_dStateREdge[RIx][lBegin+m]
				* dOrthonomREdge[lBegin+m][m_iA][m_iB][2];
		}

#ifdef UPWIND_RHO
		// Diffusion of Rho (dRho_k/dW)
		if (pGrid->GetVarLocation(RIx) == DataLocation_Node) {
			int lBeginNext = lBegin + m_nVerticalOrder;

			double dSignWL = (m_dXiDotREdge[lBegin] > 0.0)?(1.0):(-1.0);
			double dSignWR = (m_dXiDotREdge[lBeginNext] > 0.0)?(1.0):(-1.0);

			dDG[MatFIx(FWIx, lBegin, FRIx, k)] +=
				dDiffReconsPolyNode[l];

		} else {
			_EXCEPTIONT("Upwind Rho on interfaces: Unimplemented");
		}
#endif

		// dRho_k/dRho
		if (a != 0) {
			int lPrev = lBegin - m_nVerticalOrder;
			for (int n = 0; n < m_nVerticalOrder; n++) {
				dDG[MatFIx(FRIx, lPrev+n, FRIx, k)] +=
					dDiffREdgeToNode[l][0]
					* m_dJacobian3DREdge[lBegin]
					/ m_dJacobian3DNode[k]
					* 0.5 * dInterpNodeToREdge[lBegin][lPrev+n]
					* m_dXiDotREdge[lBegin];
			}
		}

		for (int m = 0; m <= m_nVerticalOrder; m++) {
		for (int n = 0; n < m_nVerticalOrder; n++) {
			double dMult = 1.0;
			if ((m == 0) || (m == m_nVerticalOrder)) {
				dMult = 0.5;
			}
			dDG[MatFIx(FRIx, lBegin+n, FRIx, k)] +=
				dDiffREdgeToNode[l][m]
				* m_dJacobian3DREdge[lBegin+m]
				/ m_dJacobian3DNode[k]
				* dMult * dInterpNodeToREdge[lBegin+m][lBegin+n]
				* m_dXiDotREdge[lBegin+m];
		}
		}

		if (a != nFiniteElements-1) {
			int lNext = lBegin + m_nVerticalOrder;
			for (int n = 0; n < m_nVerticalOrder; n++) {
				dDG[MatFIx(FRIx, lNext+n, FRIx, k)] +=
					dDiffREdgeToNode[l][m_nVerticalOrder]
					* m_dJacobian3DREdge[lNext]
					/ m_dJacobian3DNode[k]
					* 0.5 * dInterpNodeToREdge[lNext][lNext+n]
					* m_dXiDotREdge[lNext];
			}
		}
	}

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
				* m_dDiffDiffP[k];

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

	// Add the identity components
	for (int k = 0; k <= nRElements; k++) {
		dDG[MatFIx(FPIx, k, FPIx, k)] += 1.0 / m_dDeltaT;
		dDG[MatFIx(FWIx, k, FWIx, k)] += 1.0 / m_dDeltaT;
		dDG[MatFIx(FRIx, k, FRIx, k)] += 1.0 / m_dDeltaT;
	}

	for (int k = 0; k <= nRElements ; k++) {
		dDG[MatFIx(FWIx, 0, FRIx, k)] = 0.0;
		dDG[MatFIx(FWIx, nRElements, FRIx, k)] = 0.0;
	}

	dDG[MatFIx(FWIx, 0, FWIx, 0)] =
		dOrthonomREdge[0][m_iA][m_iB][2];
	dDG[MatFIx(FWIx, nRElements, FWIx, nRElements)] =
		dOrthonomREdge[nRElements][m_iA][m_iB][2];
/*
#ifdef APPLY_RAYLEIGH_WITH_VERTICALDYN

	// Rayleigh friction strength
	const GridData3D & dataRayleighStrengthNode =
		m_pPatch->GetRayleighStrength(DataLocation_Node);

	const GridData3D & dataRayleighStrengthREdge =
		m_pPatch->GetRayleighStrength(DataLocation_REdge);

	// Rayleigh friction on nodes
	for (int k = 0; k < pGrid->GetRElements(); k++) {
		double dNu = dataRayleighStrengthNode[k][m_iA][m_iB];

		if (dNu == 0.0) {
			continue;
		}

		if (pGrid->GetVarLocation(PIx) == DataLocation_Node) {
			dDG[MatFIx(FPIx, k, FPIx, k)] += dNu;
		}
		if (pGrid->GetVarLocation(WIx) == DataLocation_Node) {
			dDG[MatFIx(FWIx, k, FWIx, k)] += dNu;
		}
		if (pGrid->GetVarLocation(RIx) == DataLocation_Node) {
			dDG[MatFIx(FRIx, k, FRIx, k)] += dNu;
		}
	}

	// Rayleigh friction on interfaces
	for (int k = 0; k <= pGrid->GetRElements(); k++) {
		double dNu = dataRayleighStrengthREdge[k][m_iA][m_iB];

		if (dNu == 0.0) {
			continue;
		}

		if (pGrid->GetVarLocation(PIx) == DataLocation_REdge) {
			dDG[MatFIx(FPIx, k, FPIx, k)] += dNu;
		}
		if (pGrid->GetVarLocation(WIx) == DataLocation_REdge) {
			dDG[MatFIx(FWIx, k, FWIx, k)] += dNu;
		}
		if (pGrid->GetVarLocation(RIx) == DataLocation_REdge) {
			dDG[MatFIx(FRIx, k, FRIx, k)] += dNu;
		}
	}

#endif
*/
}

///////////////////////////////////////////////////////////////////////////////

void VerticalDynamicsFEM::BuildJacobianF2(
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

	if ((pGrid->GetVarLocation(UIx) != DataLocation_Node) ||
		(pGrid->GetVarLocation(VIx) != DataLocation_Node) ||
		(pGrid->GetVarLocation(PIx) != DataLocation_REdge) ||
		(pGrid->GetVarLocation(WIx) != DataLocation_REdge) ||
		(pGrid->GetVarLocation(RIx) != DataLocation_Node) ||
		(m_fMassFluxOnLevels) ||
		(!m_fExnerPressureOnLevels)
	) {
		_EXCEPTIONT("Not implemented");
	}

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

	// Orthonormalization coefficients
	const DataMatrix4D<double> & dOrthonomREdge =
		m_pPatch->GetOrthonormREdge();

	// Physical constants
	const PhysicalConstants & phys = m_model.GetPhysicalConstants();

	// Number of radial elements
	const int nRElements = pGrid->GetRElements();

	// Zero DG
	memset(dDG, 0,
		m_nColumnStateSize * m_nColumnStateSize * sizeof(double));

	// dT_k/dT_l
	for (int k = 0; k <= nRElements; k++) {
		int l = iDiffREdgeToREdgeBegin[k];
		for (; l < iDiffREdgeToREdgeEnd[k]; l++) {
			dDG[MatFIx(FPIx, l, FPIx, k)] +=
				dDiffREdgeToREdge[k][l]
				* m_dXiDotREdge[k];
		}
	}

	// dT_k/dW_k
	for (int k = 0; k <= nRElements; k++) {
		dDG[MatFIx(FWIx, k, FPIx, k)] =
			m_dDiffP[k] * dOrthonomREdge[k][m_iA][m_iB][2];
	}

	// dW_k/dT_l and dW_k/dR_m
	for (int k = 1; k < nRElements; k++) {

		double dRHSWCoeff = 
			1.0 / m_dDxRREdge[k]
			* m_dStateREdge[PIx][k]
			* phys.GetR() / phys.GetCv();

		int m = iDiffNodeToREdgeBegin[k];
		for (; m < iDiffNodeToREdgeEnd[k]; m++) {

			double dTEntry =
				dRHSWCoeff 
				* dDiffNodeToREdge[k][m]
				* (m_dExnerPertNode[m] + m_dExnerRefNode[m])
					/ m_dStateNode[PIx][m];

			int l = iInterpREdgeToNodeBegin[m];
			for (; l < iInterpREdgeToNodeEnd[m]; l++) {
				dDG[MatFIx(FPIx, l, FWIx, k)] +=
					dTEntry * dInterpREdgeToNode[m][l];
			}

			dDG[MatFIx(FRIx, m, FWIx, k)] +=
				dRHSWCoeff
				* dDiffNodeToREdge[k][m]
				* (m_dExnerPertNode[m] + m_dExnerRefNode[m])
					/ m_dStateNode[RIx][m];
		}
	}

	// dW_k/dT_k (first theta in RHS)
	for (int k = 1; k < nRElements; k++) {
		dDG[MatFIx(FPIx, k, FWIx, k)] +=
			 1.0 / m_dDxRREdge[k]
			 * (m_dDiffExnerRefREdge[k] + m_dDiffExnerPertREdge[k]);
	}

	// dRho
	for (int k = 0; k < nRElements; k++) {

		// dRho_k/dW_l
		int l = iDiffREdgeToNodeBegin[k];
		for (; l < iDiffREdgeToNodeEnd[k]; l++) {
			dDG[MatFIx(FWIx, l, FRIx, k)] +=
				dDiffREdgeToNode[k][l]
				* m_dJacobian3DREdge[l]
				/ m_dJacobian3DNode[k]
				* m_dStateREdge[RIx][l]
				* dOrthonomREdge[l][m_iA][m_iB][2];
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

		// dRho_k/dRho_n
		int m = iDiffREdgeToNodeBegin[k];
		for (; m < iDiffREdgeToNodeEnd[k]; m++) {

			int n = iInterpNodeToREdgeBegin[m];
			for (; n < iInterpNodeToREdgeEnd[m]; n++) {

				dDG[MatFIx(FRIx, n, FRIx, k)] +=
					dDiffREdgeToNode[k][m]
					* m_dJacobian3DREdge[m]
					/ m_dJacobian3DNode[k]
					* dInterpNodeToREdge[m][n]
					* m_dXiDotREdge[m];
			}
		}
	}

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
				* m_dDiffDiffP[k];

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

	// Add the identity components
	for (int k = 0; k <= nRElements; k++) {
		dDG[MatFIx(FPIx, k, FPIx, k)] += 1.0 / m_dDeltaT;
		dDG[MatFIx(FWIx, k, FWIx, k)] += 1.0 / m_dDeltaT;
		dDG[MatFIx(FRIx, k, FRIx, k)] += 1.0 / m_dDeltaT;
	}

#pragma message "Boundary conditions for W"
	// Enforce that density can't be a source of vertical velocity at a boundary
	for (int k = 0; k <= nRElements ; k++) {
		dDG[MatFIx(FWIx, 0, FRIx, k)] = 0.0;
		dDG[MatFIx(FWIx, nRElements, FRIx, k)] = 0.0;
	}
	
	// Enforce no horizontal velocity at the top and bottom through a boundarys
	dDG[MatFIx(FWIx, 0, FWIx, 0)] =
		dOrthonomREdge[0][m_iA][m_iB][2];
	dDG[MatFIx(FWIx, nRElements, FWIx, nRElements)] =
		dOrthonomREdge[nRElements][m_iA][m_iB][2];
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

