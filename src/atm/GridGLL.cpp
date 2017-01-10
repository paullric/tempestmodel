///////////////////////////////////////////////////////////////////////////////
///
///	\file    GridGLL.cpp
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

#include "GridGLL.h"
#include "Model.h"
#include "GridSpacing.h"

#include "Direction.h"
#include "FluxCorrectionFunction.h"
#include "PolynomialInterp.h"
#include "GaussQuadrature.h"
#include "GaussLobattoQuadrature.h"

///////////////////////////////////////////////////////////////////////////////

GridGLL::GridGLL(
	Model & model
) :
	Grid(model)
{ }

///////////////////////////////////////////////////////////////////////////////

void GridGLL::DefineParameters() {

	if (m_dcGridParameters.IsAttached()) {
		_EXCEPTIONT("Attempting to recall DefineParameters");
	}

	// Initialize the GridParameters relevant to GridGLL
	m_dcGridParameters.PushDataChunk(&m_nHorizontalOrder);
	m_dcGridParameters.PushDataChunk(&m_nVerticalOrder);

	// Call up the stack
	Grid::DefineParameters();
}

///////////////////////////////////////////////////////////////////////////////

void GridGLL::SetParameters(
	int nRElements,
	int nMaxPatchCount,
	int nABaseResolution,
	int nBBaseResolution,
	int nRefinementRatio,
	int nHorizontalOrder,
	int nVerticalOrder,
	VerticalDiscretization eVerticalDiscretization,
	VerticalStaggering eVerticalStaggering
) {
	if (!m_dcGridParameters.IsAttached()) {
		_EXCEPTIONT("DefineParameters() must be called before SetParameters()");
	}

	// Call up the stack
	Grid::SetParameters(
		nRElements,
		nMaxPatchCount,
		nABaseResolution,
		nBBaseResolution,
		nRefinementRatio,
		eVerticalDiscretization,
		eVerticalStaggering
	);

	// Set the default parameters
	m_nHorizontalOrder = nHorizontalOrder;
	m_nVerticalOrder = nVerticalOrder;
}

///////////////////////////////////////////////////////////////////////////////

void GridGLL::Initialize() {

	// Call up the stack
	Grid::Initialize();

	// Quadrature points for Gauss and Gauss-Lobatto quadrature
	DataArray1D<double> dG;
	DataArray1D<double> dW;

	DataArray1D<double> dGL;
	DataArray1D<double> dWL;

	///////////////////////////////////////////////////////////////////////////
	// Get quadrature points for Gauss-Lobatto quadrature (horizontal)

	// Derivatives of the 1D basis functions at each point on the reference
	// element [0, 1]
	m_dDxBasis1D.Allocate(m_nHorizontalOrder, m_nHorizontalOrder);
	m_dStiffness1D.Allocate(m_nHorizontalOrder, m_nHorizontalOrder);

	// Exact fourth order coefficients to machine truncation
	if (m_nHorizontalOrder == 4) {
		dGL.Allocate(4);
		dWL.Allocate(4);

		dGL[0] = 0.0;
		dGL[1] = 2.76393202250021063903773210768E-1;
		dGL[2] = 7.23606797749978936096226789232E-1;
		dGL[3] = 1.0;

		dWL[0] = 8.33333333333333148296162562474E-2;
		dWL[1] = 4.16666666666666685170383743753E-1;
		dWL[2] = 4.16666666666666685170383743753E-1;
		dWL[3] = 8.33333333333333148296162562474E-2;

		m_dDxBasis1D[0][0] = -6.0;
		m_dDxBasis1D[1][0] =  8.09016994374947451262869435595E0;
		m_dDxBasis1D[2][0] = -3.09016994374947451262869435595E0;
		m_dDxBasis1D[3][0] =  1.0;

		m_dDxBasis1D[0][1] = -1.61803398874989490252573887119E0;
		m_dDxBasis1D[1][1] =  0.0;
		m_dDxBasis1D[2][1] =  2.23606797749978980505147774238E0;
		m_dDxBasis1D[3][1] = -6.18033988749894902525738871191E-1;

		m_dDxBasis1D[0][2] =  6.18033988749894902525738871191E-1;
		m_dDxBasis1D[1][2] = -2.23606797749978980505147774238E0;
		m_dDxBasis1D[2][2] =  0.0;
		m_dDxBasis1D[3][2] =  1.61803398874989490252573887119E0;

		m_dDxBasis1D[0][3] = -1.0;
		m_dDxBasis1D[1][3] =  3.09016994374947451262869435595E0;
		m_dDxBasis1D[2][3] = -8.09016994374947451262869435595E0;
		m_dDxBasis1D[3][3] =  6.0;

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

	// Arbitrary order coefficients
	} else {

		GaussLobattoQuadrature::GetPoints(
			m_nHorizontalOrder, 0.0, 1.0, dGL, dWL);

		DataArray1D<double> dCoeffs(m_nHorizontalOrder);

		for (int i = 0; i < m_nHorizontalOrder; i++) {

			PolynomialInterp::DiffLagrangianPolynomialCoeffs(
				m_nHorizontalOrder, dGL, dCoeffs, dGL[i]);

			for (int m = 0; m < m_nHorizontalOrder; m++) {
				m_dDxBasis1D[m][i] = dCoeffs[m];

				m_dStiffness1D[m][i] = m_dDxBasis1D[m][i] * dWL[i] / dWL[m];
			}
		}
	}

	// Store GLL weights
	m_dGLLWeights1D = dWL;

	// Get the derivatives of the flux correction function
	m_dFluxDeriv1D.Allocate(m_nHorizontalOrder);
	FluxCorrectionFunction::GetDerivatives(
		2, m_nHorizontalOrder, dGL, m_dFluxDeriv1D);

	// Initialize differentiation operators
	if (m_eVerticalDiscretization == VerticalDiscretization_FiniteVolume) {

		// Interpolation operators
		m_opInterpNodeToREdge.InitializeReconstructed(
			LinearColumnInterpFEM::InterpSource_Levels,
			LinearColumnInterpFEM::InterpSource_Interfaces,
			m_nVerticalOrder,
			m_dREtaLevels,
			m_dREtaInterfaces);

		m_opInterpREdgeToNode.InitializeReconstructed(
			LinearColumnInterpFEM::InterpSource_Interfaces,
			LinearColumnInterpFEM::InterpSource_Levels,
			m_nVerticalOrder,
			m_dREtaLevels,
			m_dREtaInterfaces);

		// Differentiation operators
		m_opDiffNodeToNode.InitializeReconstructed(
			LinearColumnDiffFEM::InterpSource_Levels,
			LinearColumnDiffFEM::InterpSource_Levels,
			m_nVerticalOrder,
			m_dREtaLevels,
			m_dREtaInterfaces);

		m_opDiffNodeToREdge.InitializeReconstructed(
			LinearColumnDiffFEM::InterpSource_Levels,
			LinearColumnDiffFEM::InterpSource_Interfaces,
			m_nVerticalOrder,
			m_dREtaLevels,
			m_dREtaInterfaces);

		m_opDiffREdgeToNode.InitializeReconstructed(
			LinearColumnDiffFEM::InterpSource_Interfaces,
			LinearColumnDiffFEM::InterpSource_Levels,
			m_nVerticalOrder,
			m_dREtaLevels,
			m_dREtaInterfaces);

		m_opDiffREdgeToREdge.InitializeReconstructed(
			LinearColumnDiffFEM::InterpSource_Interfaces,
			LinearColumnDiffFEM::InterpSource_Interfaces,
			m_nVerticalOrder,
			m_dREtaLevels,
			m_dREtaInterfaces);

		// Initialize second derivative operators
		m_opDiffDiffNodeToNode.Initialize(
			LinearColumnDiffDiffFEM::InterpSource_Levels,
			1,
			m_dREtaLevels,
			m_dREtaInterfaces);

		m_opDiffDiffREdgeToREdge.Initialize(
			LinearColumnDiffDiffFEM::InterpSource_Interfaces,
			1,
			m_dREtaLevels,
			m_dREtaInterfaces);

		// Initialize discontinuous penalty operator
		m_opPenaltyNodeToNode.InitializeReconstructed(
			m_nVerticalOrder,
			m_dREtaLevels,
			m_dREtaInterfaces);

	// Finite element discretization on interfaces
	} else if (m_eVerticalStaggering == VerticalStaggering_Interfaces) {

		// Special differentiation operator
		m_opDiffNodeToNode.InitializeGLLNodes(
			m_nVerticalOrder,
			m_dREtaLevels,
			m_dREtaLevels);

		// Special differentiation operator
		m_opDiffDiffNodeToNode.InitializeGLLNodes(
			m_nVerticalOrder,
			m_dREtaLevels);
/*
		// Initialize integral operator
		m_opIntNodeToNode.InitializeNodeToNodeInterfaceMethod(
			m_nVerticalOrder,
			m_dREtaLevels,
			m_dREtaInterfaces);
*/

	} else {
		// Interpolation operators
		m_opInterpNodeToREdge.Initialize(
			LinearColumnInterpFEM::InterpSource_Levels,
			m_nVerticalOrder,
			m_dREtaStretchLevels,
			m_dREtaStretchInterfaces,
			m_dREtaStretchInterfaces);

		m_opInterpREdgeToNode.Initialize(
			LinearColumnInterpFEM::InterpSource_Interfaces,
			m_nVerticalOrder,
			m_dREtaStretchLevels,
			m_dREtaStretchInterfaces,
			m_dREtaStretchLevels);

		// Differentiation operators
		m_opDiffNodeToNode.InitializeInterfaceMethod(
			LinearColumnDiffFEM::InterpSource_Levels,
			m_nVerticalOrder,
			m_dREtaLevels,
			m_dREtaInterfaces,
			m_dREtaLevels,
			false);

		m_opDiffNodeToNodeZeroBoundaries.InitializeInterfaceMethod(
			LinearColumnDiffFEM::InterpSource_Levels,
			m_nVerticalOrder,
			m_dREtaLevels,
			m_dREtaInterfaces,
			m_dREtaLevels,
			true);

		m_opDiffNodeToREdge.InitializeFluxCorrectionMethod(
			LinearColumnDiffFEM::InterpSource_Levels,
			m_nVerticalOrder,
			m_dREtaLevels,
			m_dREtaInterfaces,
			m_dREtaInterfaces,
			true);
/*
		m_opDiffNodeToREdge.InitializeVariationalNodeToREdge(
			m_nVerticalOrder,
			m_dREtaLevels,
			m_dREtaInterfaces);
*/
		m_opDiffREdgeToNode.InitializeInterfaceMethod(
			LinearColumnDiffFEM::InterpSource_Interfaces,
			m_nVerticalOrder,
			m_dREtaLevels,
			m_dREtaInterfaces,
			m_dREtaLevels,
			false);

		m_opDiffREdgeToREdge.InitializeInterfaceMethod(
			LinearColumnDiffFEM::InterpSource_Interfaces,
			m_nVerticalOrder,
			m_dREtaLevels,
			m_dREtaInterfaces,
			m_dREtaInterfaces,
			false);

		// Initialize second derivative operators
		m_opDiffDiffNodeToNode.Initialize(
			LinearColumnDiffDiffFEM::InterpSource_Levels,
			m_nVerticalOrder,
			m_dREtaLevels,
			m_dREtaInterfaces);

		m_opDiffDiffREdgeToREdge.Initialize(
			LinearColumnDiffDiffFEM::InterpSource_Interfaces,
			m_nVerticalOrder,
			m_dREtaLevels,
			m_dREtaInterfaces);

		// Initialize discontinuous penalty operator
		m_opPenaltyNodeToNode.Initialize(
			m_nVerticalOrder,
			m_dREtaLevels,
			m_dREtaInterfaces);

		// Initialize integral operator
		m_opIntNodeToNode.InitializeNodeToNodeInterfaceMethod(
			m_nVerticalOrder,
			m_dREtaLevels,
			m_dREtaInterfaces);
	}
}

///////////////////////////////////////////////////////////////////////////////

void GridGLL::InitializeVerticalCoordinate() {

	// Initialize variable locations
	Grid::InitializeVerticalCoordinate();

	// Initialize 
	if (m_model.GetEquationSet().GetDimensionality() == 2) {
		return;
	}

	// For finite volume discretization, backup vertical order
	int nNodesPerFiniteElement = m_nVerticalOrder;
	if (m_eVerticalDiscretization == VerticalDiscretization_FiniteVolume) {
		nNodesPerFiniteElement = 1;
	}

	// Initialize the vertical coordinate (INT staggering)
	if (m_eVerticalStaggering == VerticalStaggering_Interfaces) {

		// Check number of levels
		if ((m_nRElements - 1) % (nNodesPerFiniteElement - 1) != 0) {
			_EXCEPTIONT("(vertorder - 1) must divide (levels - 1) equally");
		}

		// Get Gauss-Lobatto points
		DataArray1D<double> dGL;
		DataArray1D<double> dWL;

		GaussLobattoQuadrature::GetPoints(
			nNodesPerFiniteElement, 0.0, 1.0, dGL, dWL);

		// Number of finite elements
		int nFiniteElements =
			(m_nRElements - 1) / (nNodesPerFiniteElement - 1);

		double dAvgDeltaElement =
			1.0 / static_cast<double>(nFiniteElements);

		// Uniform stretching
		if (m_pVerticalStretchF == NULL) {
			for (int k = 0; k < m_nRElements; k++) {
				double dA =
					static_cast<double>(k / (nNodesPerFiniteElement - 1));

				int kx = k % (nNodesPerFiniteElement - 1);

				m_dREtaLevels[k] = (dGL[kx] + dA) * dAvgDeltaElement;
				m_dREtaLevelsNormArea[k] = dWL[kx] * dAvgDeltaElement;
			}
			for (int a = 1; a < nFiniteElements; a++) {
				m_dREtaLevelsNormArea[a * (nNodesPerFiniteElement-1)] *= 2.0;
			}

		// Non-uniform stretching
		} else {

			double dStretchREta = 0.0;
			double dDxStretchREta = 0.0;

			double dREta0 = 0.0;
			double dREta1 = 0.0;

			EvaluateVerticalStretchF(
				dAvgDeltaElement,
				dREta1,
				dDxStretchREta);

			m_dREtaLevelsNormArea.Zero();

			for (int a = 0; a < nFiniteElements; a++) {

				double dDeltaElement = dREta1 - dREta0;
				for (int k = 0; k < nNodesPerFiniteElement; k++) {
					int kx = a * (nNodesPerFiniteElement - 1) + k;
					m_dREtaLevels[kx] = dREta0 + dGL[k] * dDeltaElement;
					m_dREtaLevelsNormArea[kx] += dWL[k] * dDeltaElement;
				}

				// Update element bounds
				dREta0 = dREta1;

				EvaluateVerticalStretchF(
					static_cast<double>(a+2) * dAvgDeltaElement,
					dREta1,
					dDxStretchREta);
			}
		}

	// Initialize the vertical coordinate (LEV / LOR / CPH staggering)
	} else {

		// Check number of levels
		if (m_nRElements % nNodesPerFiniteElement != 0) {
			_EXCEPTIONT("vertorder must divide levels equally");
		}

		// Get Gauss points
		DataArray1D<double> dG;
		DataArray1D<double> dW;

		GaussQuadrature::GetPoints(
			nNodesPerFiniteElement, 0.0, 1.0, dG, dW);

		// Get Gauss-Lobatto points
		DataArray1D<double> dGL;
		DataArray1D<double> dWL;

		GaussLobattoQuadrature::GetPoints(
			nNodesPerFiniteElement+1, 0.0, 1.0, dGL, dWL);

		// Number of finite elements
		int nFiniteElements = m_nRElements / nNodesPerFiniteElement;
		double dAvgDeltaElement = 1.0 / static_cast<double>(nFiniteElements);

		// Uniform stretching
		if (m_pVerticalStretchF == NULL) {

			for (int k = 0; k < m_nRElements; k++) {
				double dA = static_cast<double>(k / nNodesPerFiniteElement);
				int kx = k % nNodesPerFiniteElement;

				m_dREtaLevels[k] = (dG[kx] + dA) * dAvgDeltaElement;
				m_dREtaLevelsNormArea[k] = dW[kx] * dAvgDeltaElement;
			}

			for (int k = 0; k <= m_nRElements; k++) {
				double dA = static_cast<double>(k / nNodesPerFiniteElement);
				int kx = k % nNodesPerFiniteElement;

				m_dREtaInterfaces[k] = (dGL[kx] + dA) * dAvgDeltaElement;
				m_dREtaInterfacesNormArea[k] = dWL[kx] * dAvgDeltaElement;
			}
			for (int a = 1; a < nFiniteElements; a++) {
				m_dREtaInterfacesNormArea[a * nNodesPerFiniteElement] *= 2.0;
			}

		// Non-uniform stretching
		} else {

			double dStretchREta = 0.0;
			double dDxStretchREta = 0.0;

			double dREta0 = 0.0;
			double dREta1 = 0.0;

			EvaluateVerticalStretchF(
				dAvgDeltaElement,
				dREta1,
				dDxStretchREta);

			m_dREtaLevelsNormArea.Zero();
			m_dREtaInterfacesNormArea.Zero();

			for (int a = 0; a < nFiniteElements; a++) {

				double dDeltaElement = dREta1 - dREta0;

				for (int k = 0; k < nNodesPerFiniteElement; k++) {
					int kx = a * nNodesPerFiniteElement + k;
					m_dREtaLevels[kx] = dREta0 + dG[k] * dDeltaElement;
					m_dREtaLevelsNormArea[kx] = dW[k] * dDeltaElement;
				}

				for (int k = 0; k <= nNodesPerFiniteElement; k++) {
					int kx = a * nNodesPerFiniteElement + k;
					m_dREtaInterfaces[kx] = dREta0 + dGL[k] * dDeltaElement;
					m_dREtaInterfacesNormArea[kx] += dWL[k] * dDeltaElement;
				}

				// Update element bounds
				dREta0 = dREta1;

				EvaluateVerticalStretchF(
					static_cast<double>(a+2) * dAvgDeltaElement,
					dREta1,
					dDxStretchREta);
			}

		}
	}

	// Copy to stretch arrays (deprecated)
	m_dREtaStretchLevels = m_dREtaLevels;
	m_dREtaStretchInterfaces = m_dREtaInterfaces;

}

///////////////////////////////////////////////////////////////////////////////

void GridGLL::EvaluateTestCase(
	const TestCase & test,
	const Time & time,
	int iDataIndex
) {
	// Call the standard implementation of EvaluateTestCase
	Grid::EvaluateTestCase(test, time, iDataIndex);

	// Use DSS to average topographic derivatives
	ApplyDSS(0, DataType_TopographyDeriv);
}

///////////////////////////////////////////////////////////////////////////////

void GridGLL::PostProcessSubstage(
	int iDataUpdate,
	DataType eDataType
) {
	// Block parallel exchanges
	if (m_fBlockParallelExchange) {
		return;
	}

	// Apply Direct Stiffness Summation
	ApplyDSS(iDataUpdate, eDataType);

}

///////////////////////////////////////////////////////////////////////////////

void GridGLL::ComputeVorticityDivergence(
	int iDataIndex
) {
	// Block parallel exchanges
	if (m_fBlockParallelExchange) {
		return;
	}

	// Compute vorticity on all grid patches
	Grid::ComputeVorticityDivergence(iDataIndex);

	// Apply DSS
	ApplyDSS(0, DataType_Vorticity);
	ApplyDSS(0, DataType_Divergence);

}

///////////////////////////////////////////////////////////////////////////////

double GridGLL::InterpolateNodeToREdge(
	const double * dDataNode,
	const double * dDataRefNode,
	int iRint,
	double dDataRefREdge,
	int nStride
) const {

	// Apply operator
	if (dDataRefNode == NULL) {
		return m_opInterpNodeToREdge.Apply(
			dDataNode,
			iRint,
			nStride);

	} else {
		return m_opInterpNodeToREdge.Apply(
			dDataNode,
			dDataRefNode,
			dDataRefREdge,
			iRint,
			nStride);
	} 
}

///////////////////////////////////////////////////////////////////////////////

double GridGLL::InterpolateREdgeToNode(
	const double * dDataREdge,
	const double * dDataRefREdge,
	int iRnode,
	double dDataRefNode,
	int nStride
) const {

	// Apply operator
	if (dDataRefREdge == NULL) {
		return m_opInterpREdgeToNode.Apply(
			dDataREdge,
			iRnode,
			nStride);

	} else {
		return m_opInterpREdgeToNode.Apply(
			dDataREdge,
			dDataRefREdge,
			dDataRefNode,
			iRnode,
			nStride);
	}
}

///////////////////////////////////////////////////////////////////////////////

double GridGLL::DifferentiateNodeToNode(
	const double * dDataNode,
	int iRnode,
	int nStride
) const {

	// Apply operator
	return m_opDiffNodeToNode.Apply(
		dDataNode,
		iRnode,
		nStride);
}

///////////////////////////////////////////////////////////////////////////////

double GridGLL::DifferentiateNodeToREdge(
	const double * dDataNode,
	int iRnode,
	int nStride
) const {

	// Apply operator
	return m_opDiffNodeToREdge.Apply(
		dDataNode,
		iRnode,
		nStride);
}

///////////////////////////////////////////////////////////////////////////////

double GridGLL::DifferentiateREdgeToNode(
	const double * dDataREdge,
	int iRnode,
	int nStride
) const {

	// Apply operator
	return m_opDiffREdgeToNode.Apply(
		dDataREdge,
		iRnode,
		nStride);
}

///////////////////////////////////////////////////////////////////////////////

double GridGLL::DifferentiateREdgeToREdge(
	const double * dDataREdge,
	int iRint,
	int nStride
) const {

	// Apply operator
	return m_opDiffREdgeToREdge.Apply(
		dDataREdge,
		iRint,
		nStride);
}

///////////////////////////////////////////////////////////////////////////////

void GridGLL::InterpolateNodeToREdge(
	const double * dDataNode,
	double * dDataREdge,
	bool fZeroBoundaries
) const {

	// Apply operator
	m_opInterpNodeToREdge.Apply(
		dDataNode,
		dDataREdge);

	// Override boundary values if zero
	if (fZeroBoundaries) {
		dDataREdge[0] = 0.0;
		dDataREdge[GetRElements()] = 0.0;
	}
}

///////////////////////////////////////////////////////////////////////////////

void GridGLL::InterpolateREdgeToNode(
	const double * dDataREdge,
	double * dDataNode
) const {

	// Apply operator
	m_opInterpREdgeToNode.Apply(
		dDataREdge,
		dDataNode);
}

///////////////////////////////////////////////////////////////////////////////

void GridGLL::DifferentiateNodeToNode(
	const double * dDataNode,
	double * dDiffNode,
	bool fZeroBoundaries
) const {

	// Apply operator
	if (fZeroBoundaries) {
		m_opDiffNodeToNodeZeroBoundaries.Apply(
			dDataNode,
			dDiffNode);

	} else {
		m_opDiffNodeToNode.Apply(
			dDataNode,
			dDiffNode);
	}
}

///////////////////////////////////////////////////////////////////////////////

void GridGLL::DifferentiateNodeToREdge(
	const double * dDataNode,
	double * dDiffREdge,
	bool fZeroBoundaries
) const {

	// Apply operator
	m_opDiffNodeToREdge.Apply(
		dDataNode,
		dDiffREdge);

	if (fZeroBoundaries) {
		_EXCEPTIONT("ZeroBoundaries broken - sorry!");
	}
}

///////////////////////////////////////////////////////////////////////////////

void GridGLL::DifferentiateREdgeToNode(
	const double * dDataREdge,
	double * dDiffNode
) const {

	// Apply operator
	m_opDiffREdgeToNode.Apply(
		dDataREdge,
		dDiffNode);
}

///////////////////////////////////////////////////////////////////////////////

void GridGLL::DifferentiateREdgeToREdge(
	const double * dDataREdge,
	double * dDiffREdge
) const {

	// Apply operator
	m_opDiffREdgeToREdge.Apply(
		dDataREdge,
		dDiffREdge);
}

///////////////////////////////////////////////////////////////////////////////

void GridGLL::DiffDiffNodeToNode(
	const double * dDataNode,
	double * dDiffDiffNode
) const {

	// Apply operator
	m_opDiffDiffNodeToNode.Apply(
		dDataNode,
		dDiffDiffNode);
}

///////////////////////////////////////////////////////////////////////////////

void GridGLL::DiffDiffREdgeToREdge(
	const double * dDataREdge,
	double * dDiffDiffREdge
) const {

	// Apply operator
	m_opDiffDiffREdgeToREdge.Apply(
		dDataREdge,
		dDiffDiffREdge);
}

///////////////////////////////////////////////////////////////////////////////

