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
#include "HorizontalDynamicsDG.h"

#include "Direction.h"
#include "FluxCorrectionFunction.h"
#include "PolynomialInterp.h"
#include "GaussLobattoQuadrature.h"

///////////////////////////////////////////////////////////////////////////////

GridGLL::GridGLL(
	Model & model,
	int nBaseResolutionA,
	int nBaseResolutionB,
	int nRefinementRatio,
	int nHorizontalOrder,
	int nVerticalOrder,
	int nRElements,
	VerticalStaggering eVerticalStaggering
) :
	// Call up the stack
	Grid::Grid(
		model,
		nBaseResolutionA,
		nBaseResolutionB,
		nRefinementRatio,
		nRElements,
		eVerticalStaggering),
	m_nHorizontalOrder(nHorizontalOrder),
	m_nVerticalOrder(nVerticalOrder)
{
}

///////////////////////////////////////////////////////////////////////////////

void GridGLL::Initialize() {

	// Call up the stack
	Grid::Initialize();

	// Initialize the vertical coordinate (INT staggering)
	if (m_eVerticalStaggering == VerticalStaggering_Interfaces) {
		double dDeltaElement =
			static_cast<double>(m_nVerticalOrder - 1)
			/ static_cast<double>(m_nRElements - 1);

		InitializeVerticalCoordinate(
			GridSpacingGaussLobatto(dDeltaElement, 0.0, m_nVerticalOrder)
		);

	// Initialize the vertical coordinate (LEV / LOR / CPH staggering)
	} else {
		double dDeltaElement =
			static_cast<double>(m_nVerticalOrder)
			/ static_cast<double>(m_nRElements);

		InitializeVerticalCoordinate(
			GridSpacingMixedGaussLobatto(dDeltaElement, 0.0, m_nVerticalOrder)
		);
	}

	// Quadrature points for Gauss and Gauss-Lobatto quadrature
	DataVector<double> dG;
	DataVector<double> dW;

	DataVector<double> dGL;
	DataVector<double> dWL;

	///////////////////////////////////////////////////////////////////////////
	// Get quadrature points for Gauss-Lobatto quadrature (horizontal)
	GaussLobattoQuadrature::GetPoints(m_nHorizontalOrder, 0.0, 1.0, dGL, dWL);

	// Derivatives of the 1D basis functions at each point on the reference
	// element [0, 1]
	m_dDxBasis1D.Initialize(m_nHorizontalOrder, m_nHorizontalOrder);
	m_dStiffness1D.Initialize(m_nHorizontalOrder, m_nHorizontalOrder);

	DataVector<double> dCoeffs;
	dCoeffs.Initialize(m_nHorizontalOrder);

	for (int i = 0; i < m_nHorizontalOrder; i++) {
		PolynomialInterp::DiffLagrangianPolynomialCoeffs(
			m_nHorizontalOrder, dGL, dCoeffs, dGL[i]);
/*
		if (m_nHorizontalOrder == 4) {
			if (i == 0) {
				dCoeffs[0] = -6.0;
				dCoeffs[1] =  8.090169943749478;
				dCoeffs[2] = -3.090169943749475;
				dCoeffs[3] =  1.0;
			} else if (i == 1) {
				dCoeffs[0] = -1.618033988749895;
				dCoeffs[1] =  0.0;
				dCoeffs[2] =  2.236067977499781;
				dCoeffs[3] = -0.618033988749895;
			} else if (i == 2) {
				dCoeffs[0] =  0.618033988749895;
				dCoeffs[1] = -2.236067977499781;
				dCoeffs[2] =  0.0;
				dCoeffs[3] =  1.618033988749895;
			} else {
				dCoeffs[0] = -1.0;
				dCoeffs[1] =  3.090169943749476;
				dCoeffs[2] = -8.090169943749475;
				dCoeffs[3] =  6.0;
			}
		}
*/
#pragma message "Verify that DxBasis1D adds up to 0"
		for (int m = 0; m < m_nHorizontalOrder; m++) {
			m_dDxBasis1D[m][i] = dCoeffs[m];

			m_dStiffness1D[m][i] = m_dDxBasis1D[m][i] * dWL[i] / dWL[m];
		}
	}

	// Store GLL weights
	m_dGLLWeights1D = dWL;

	// Get the derivatives of the flux correction function
	m_dFluxDeriv1D.Initialize(m_nHorizontalOrder);
	FluxCorrectionFunction::GetDerivatives(
		2, m_nHorizontalOrder, dGL, m_dFluxDeriv1D);

	// Initialize differentiation operators
	if (m_eVerticalStaggering == VerticalStaggering_Interfaces) {

		// Special differentiation operator
		m_opDiffNodeToNode.InitializeGLLNodes(
			m_nVerticalOrder,
			m_dREtaLevels,
			m_dREtaLevels);

		// Special differentiation operator
		m_opDiffDiffNodeToNode.InitializeGLLNodes(
			m_nVerticalOrder,
			m_dREtaLevels);

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
			false);

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
	}
/*
	FILE * fp1 = fopen("op1.txt", "w");
	const DataMatrix<double> & dCoeff1 = m_opDiffNodeToNode.GetCoeffs();
	for (int n = 0; n < dCoeff1.GetRows(); n++) {
		for (int m = 0; m < dCoeff1.GetColumns(); m++) {
			fprintf(fp1, "%1.16e\t", dCoeff1[n][m]);
		}
		fprintf(fp1, "\n");
	}
	fclose(fp1);

	FILE * fp2 = fopen("lev.txt", "w");
	for (int n = 0; n < m_dREtaStretchLevels.GetRows(); n++) {
		fprintf(fp2, "%1.16e\n", m_dREtaStretchLevels[n]);
	}
	fclose(fp2);
	_EXCEPTION();
*/
/*
	const DataVector<int> & ixBegin1 = m_opDiffNodeToNode.GetIxBegin();
	const DataVector<int> & ixEnd1 = m_opDiffNodeToNode.GetIxEnd();
	for (int n = 0; n < dCoeff1.GetRows(); n++) {
		fprintf(fp1, "%i\t%i\n", ixBegin1[n], ixEnd1[n]);
	}
	for (int n = 0; n < m_dREtaStretchLevels.GetRows(); n++) {
		fprintf(fp1, "%1.5e %1.5e\n", m_dREtaStretchInterfaces[n], m_dREtaStretchLevels[n]);
	}
*/
/*
	FILE * fp2 = fopen("op2.txt", "w");
	const DataMatrix<double> & dCoeff2 = m_opInterpREdgeToNode.GetCoeffs();
	for (int n = 0; n < dCoeff2.GetRows(); n++) {
		for (int m = 0; m < dCoeff2.GetColumns(); m++) {
			fprintf(fp2, "%1.5e\t", dCoeff2[n][m]);
		}
		fprintf(fp2, "\n");
	}
	const DataVector<int> & ixBegin2 = m_opInterpREdgeToNode.GetIxBegin();
	const DataVector<int> & ixEnd2 = m_opInterpREdgeToNode.GetIxEnd();
	for (int n = 0; n < dCoeff2.GetRows(); n++) {
		fprintf(fp2, "%i\t%i\n", ixBegin2[n], ixEnd2[n]);
	}
	fclose(fp2);
*/
/*
	///////////////////////////////////////////////////////////////////////////
	// Get quadrature points for Gauss quadrature (vertical)
	GaussQuadrature::GetPoints(m_nVerticalOrder, 0.0, 1.0, dG, dW);

	// Get quadrature points for Gauss-Lobatto quadrature (vertical)
	GaussLobattoQuadrature::GetPoints(m_nVerticalOrder+1, 0.0, 1.0, dGL, dWL);

	// Storage for variables on element interfaces (used for computing
	// derivatives from nodes)
	int nFiniteElements = GetRElements() / m_nVerticalOrder;
	if (GetRElements() % m_nVerticalOrder != 0) {
		_EXCEPTIONT("Logic error: Vertical order must divide RElements");
	}

	// Vertical elemental grid spacing
	double dElementDeltaXi =
		static_cast<double>(m_nVerticalOrder)
		/ static_cast<double>(GetRElements());

	// Interpolation coefficients from nodes to interfaces and vice versa
	m_dInterpNodeToREdge.Initialize(m_nVerticalOrder+1, m_nVerticalOrder);
	m_dInterpREdgeToNode.Initialize(m_nVerticalOrder, m_nVerticalOrder+1);

	for (int n = 0; n < m_nVerticalOrder+1; n++) {
		PolynomialInterp::LagrangianPolynomialCoeffs(
			m_nVerticalOrder,
			m_dREtaLevels,
			m_dInterpNodeToREdge[n],
			m_dREtaInterfaces[n]);
	}
	for (int n = 0; n < m_nVerticalOrder; n++) {
		PolynomialInterp::LagrangianPolynomialCoeffs(
			m_nVerticalOrder+1,
			m_dREtaInterfaces,
			m_dInterpREdgeToNode[n],
			m_dREtaLevels[n]);
	}

	// Differentiation coefficients
	m_dDiffREdgeToNode.Initialize(m_nVerticalOrder, m_nVerticalOrder+1);
	m_dDiffREdgeToREdge.Initialize(m_nVerticalOrder+1, m_nVerticalOrder+1);
	m_dDiffNodeToREdge.Initialize(m_nVerticalOrder+1, m_nVerticalOrder);
	m_dDiffNodeToNode.Initialize(m_nVerticalOrder, m_nVerticalOrder);

	// Compute differentiation coefficients
	for (int n = 0; n < m_nVerticalOrder; n++) {
		PolynomialInterp::DiffLagrangianPolynomialCoeffs(
			m_nVerticalOrder+1, dGL, m_dDiffREdgeToNode[n], dG[n]);
		PolynomialInterp::DiffLagrangianPolynomialCoeffs(
			m_nVerticalOrder, dG, m_dDiffNodeToNode[n], dG[n]);

		for (int m = 0; m <= m_nVerticalOrder; m++) {
			m_dDiffREdgeToNode[n][m] /= dElementDeltaXi;
		}
		for (int m = 0; m < m_nVerticalOrder; m++) {
			m_dDiffNodeToNode[n][m] /= dElementDeltaXi;
		}
	}
	for (int n = 0; n <= m_nVerticalOrder; n++) {
		PolynomialInterp::DiffLagrangianPolynomialCoeffs(
			m_nVerticalOrder, dG, m_dDiffNodeToREdge[n], dGL[n]);
		PolynomialInterp::DiffLagrangianPolynomialCoeffs(
			m_nVerticalOrder+1, dGL, m_dDiffREdgeToREdge[n], dGL[n]);

		for (int m = 0; m <= m_nVerticalOrder; m++) {
			m_dDiffREdgeToREdge[n][m] /= dElementDeltaXi;
		}
		for (int m = 0; m < m_nVerticalOrder; m++) {
			m_dDiffNodeToREdge[n][m] /= dElementDeltaXi;
		}
	}

	// Get derivatives of flux reconstruction function and scale to the
	// element [0, dElementDeltaXi]
	FluxCorrectionFunction::GetDerivatives(
		m_nReconstructionPolyType,
		m_nVerticalOrder+1, dG, m_dDiffReconsPolyNode);

	FluxCorrectionFunction::GetDerivatives(
		m_nReconstructionPolyType,
		m_nVerticalOrder+1, dGL, m_dDiffReconsPolyREdge);

	for (int n = 0; n < m_dDiffReconsPolyNode.GetRows(); n++) {
		m_dDiffReconsPolyNode[n] /= dElementDeltaXi;
	}
	for (int n = 0; n < m_dDiffReconsPolyREdge.GetRows(); n++) {
		m_dDiffReconsPolyREdge[n] /= dElementDeltaXi;
	}

	// Compute amalgamated differentiation coefficients
	m_dDiffNodeToREdgeAmal.Initialize(
		m_nVerticalOrder+1, 3*m_nVerticalOrder);

	m_dDiffNodeToREdgeLeft.Initialize(
		m_nVerticalOrder+1, 2*m_nVerticalOrder);

	m_dDiffNodeToREdgeRight.Initialize(
		m_nVerticalOrder+1, 2*m_nVerticalOrder);

	for (int n = 0; n <= m_nVerticalOrder; n++) {
	for (int m = 0; m < m_nVerticalOrder; m++) {
		m_dDiffNodeToREdgeAmal[n][m_nVerticalOrder + m] =
			m_dDiffNodeToREdge[n][m];
		m_dDiffNodeToREdgeLeft[n][m] =
			m_dDiffNodeToREdge[n][m];
		m_dDiffNodeToREdgeRight[n][m_nVerticalOrder + m] =
			m_dDiffNodeToREdge[n][m];
	}
	}

	// Overlay differentiation stencils
	for (int m = 0; m < 2 * m_nVerticalOrder; m++) {
		m_dDiffNodeToREdgeAmal[0][m] = 0.5 * (
			  m_dDiffNodeToREdgeAmal[0][m]
			+ m_dDiffNodeToREdgeAmal[m_nVerticalOrder][m_nVerticalOrder+m]);

		m_dDiffNodeToREdgeAmal[m_nVerticalOrder][m_nVerticalOrder+m] =
			m_dDiffNodeToREdgeAmal[0][m];
	}

	// Contributions due to interface values
	for (int n = 0; n <= m_nVerticalOrder; n++) {
	for (int m = 0; m < m_nVerticalOrder; m++) {
		// Contribution from element on the right
		m_dDiffNodeToREdgeAmal[n][m_nVerticalOrder + m] -=
			0.5 * m_dInterpNodeToREdge[m_nVerticalOrder][m]
			* m_dDiffReconsPolyREdge[n];

		m_dDiffNodeToREdgeAmal[n][2 * m_nVerticalOrder + m] +=
			0.5 * m_dInterpNodeToREdge[0][m]
			* m_dDiffReconsPolyREdge[n];

		m_dDiffNodeToREdgeLeft[n][m] -=
			0.5 * m_dInterpNodeToREdge[m_nVerticalOrder][m] * (
				m_dDiffReconsPolyREdge[n]
				+ m_dDiffReconsPolyREdge[m_nVerticalOrder - n]);

		m_dDiffNodeToREdgeLeft[n][m_nVerticalOrder + m] +=
			0.5 * m_dInterpNodeToREdge[0][m] * (
				m_dDiffReconsPolyREdge[n]
				+ m_dDiffReconsPolyREdge[m_nVerticalOrder - n]);

		// Contribution from element on the left
		m_dDiffNodeToREdgeAmal[n][m] +=
			- 0.5 * m_dInterpNodeToREdge[m_nVerticalOrder][m]
			* m_dDiffReconsPolyREdge[m_nVerticalOrder - n];

		m_dDiffNodeToREdgeAmal[n][m_nVerticalOrder + m] -=
			- 0.5* m_dInterpNodeToREdge[0][m]
			* m_dDiffReconsPolyREdge[m_nVerticalOrder - n];

		m_dDiffNodeToREdgeRight[n][m] +=
			- 0.5 * m_dInterpNodeToREdge[m_nVerticalOrder][m] * (
				m_dDiffReconsPolyREdge[m_nVerticalOrder - n]
				+ m_dDiffReconsPolyREdge[n]);

		m_dDiffNodeToREdgeRight[n][m_nVerticalOrder + m] -=
			- 0.5 * m_dInterpNodeToREdge[0][m] * (
				m_dDiffReconsPolyREdge[m_nVerticalOrder - n]
				+ m_dDiffReconsPolyREdge[n]);
	}
	}
*/
/*
	// Test
	double dValue[12];
	double dDiffValue[12];
	for (int k = 0; k < 12; k++) {
		dValue[k] = 15.0 * m_dREtaStretchLevels[k];
	}

	DifferentiateNodeToNode(
		dValue,
		dDiffValue,
		false);

	for (int k = 0; k < 12; k++) {
		printf("%i: %1.5e\n", k, dDiffValue[k]);
	}
	m_opDiffNodeToNode.Apply(dValue, dDiffValue);
	for (int k = 0; k < 12; k++) {
		printf("%i: %1.5e\n", k, dDiffValue[k]);
	}
	_EXCEPTION();
*/
}

///////////////////////////////////////////////////////////////////////////////

void GridGLL::EvaluateTestCase(
	const TestCase & test,
	const Time & time,
	int iDataIndex
) {
	// Call the standard implementation of EvaluateTestCase
	Grid::EvaluateTestCase(test, time, iDataIndex);
/*
	const GridPatch * pPatch3 = GetPatch(3);
	const GridPatch * pPatch4 = GetPatch(4);

	const GridData3D & dataTopographyDeriv3 = pPatch3->GetTopographyDeriv();
	const GridData3D & dataTopographyDeriv4 = pPatch4->GetTopographyDeriv();

	printf("\n\n");
	printf("3i %1.15e %1.15e\n",
		dataTopographyDeriv3[0][43][64],
		dataTopographyDeriv3[1][43][64]);
	printf("4i %1.15e %1.15e\n",
		dataTopographyDeriv4[0][1][65-43],
		dataTopographyDeriv4[1][1][65-43]);
*/
	// Use DSS to average topographic derivatives
	ApplyDSS(0, DataType_TopographyDeriv);
/*
	printf("3o %1.15e %1.15e\n",
		dataTopographyDeriv3[0][43][64],
		dataTopographyDeriv3[1][43][64]);
	printf("4o %1.15e %1.15e\n",
		dataTopographyDeriv4[0][1][65-43],
		dataTopographyDeriv4[1][1][65-43]);

	printf("3x %1.15e %1.15e\n",
		tan(pPatch3->GetPatchBox().GetANode(43)),
		tan(pPatch3->GetPatchBox().GetBNode(64)));
	printf("4x %1.15e %1.15e\n",
		tan(pPatch4->GetPatchBox().GetANode(1)),
		tan(pPatch4->GetPatchBox().GetBNode(65-43)));

	printf("\n\n");
	_EXCEPTION();
*/
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

	// Get the pointer to the HorizontalDynamics object
	const HorizontalDynamicsDG * pHorizontalDynamics =
		dynamic_cast<const HorizontalDynamicsDG *>(
			m_model.GetHorizontalDynamics());

	// If SpectralElement dynamics are used, apply direct stiffness summation
	if (pHorizontalDynamics == NULL) {
		ApplyDSS(iDataUpdate, eDataType);
	}

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

