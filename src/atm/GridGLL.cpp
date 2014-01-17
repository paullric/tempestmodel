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

#include "Direction.h"
#include "FluxReconstructionFunction.h"
#include "PolynomialInterp.h"
#include "GaussLobattoQuadrature.h"
#include "GaussQuadrature.h"

///////////////////////////////////////////////////////////////////////////////

GridGLL::GridGLL(
	Model & model,
	int nBaseResolutionA,
	int nBaseResolutionB,
	int nRefinementRatio,
	int nHorizontalOrder,
	int nVerticalOrder,
	int nRElements
) :
	// Call up the stack
	Grid::Grid(
		model,
		nBaseResolutionA,
		nBaseResolutionB,
		nRefinementRatio,
		nRElements),
	m_nHorizontalOrder(nHorizontalOrder),
	m_nVerticalOrder(nVerticalOrder),
	m_nReconstructionPolyType(2)
{
}

///////////////////////////////////////////////////////////////////////////////

void GridGLL::Initialize() {

	// Call up the stack
	Grid::Initialize();

	// Initialize the vertical coordinate
	double dDeltaElement =
		static_cast<double>(m_nVerticalOrder)
		/ static_cast<double>(m_nRElements);

	InitializeVerticalCoordinate(
		GridSpacingMixedGaussLobatto(dDeltaElement, 0.0, m_nVerticalOrder)
	);

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

#pragma message "Verify that DxBasis1D adds up to 0"
		for (int m = 0; m < m_nHorizontalOrder; m++) {
			m_dDxBasis1D[m][i] = dCoeffs[m];

			m_dStiffness1D[m][i] = m_dDxBasis1D[m][i] * dWL[i] / dWL[m];
		}
	}

	///////////////////////////////////////////////////////////////////////////
	// Get quadrature points for Gauss quadrature (vertical)
	GaussQuadrature::GetPoints(m_nVerticalOrder, 0.0, 1.0, dG, dW);

	// Get quadrature points for Gauss-Lobatto quadrature (vertical)
	GaussLobattoQuadrature::GetPoints(m_nVerticalOrder+1, 0.0, 1.0, dGL, dWL);

	// Vertical grid spacing
	double dElementDeltaXi =
		  GetREtaInterface(m_nVerticalOrder) - GetREtaInterface(0);

	// Storage for variables on element interfaces (used for computing
	// derivatives from nodes)
	int nFiniteElements = GetRElements() / m_nVerticalOrder;
	if (GetRElements() % m_nVerticalOrder != 0) {
		_EXCEPTIONT("Logic error: Vertical order must divide RElements");
	}

	m_dStateFEEdge.Initialize(nFiniteElements+1, 2);

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
	FluxReconstructionFunction::GetDerivatives(
		m_nReconstructionPolyType,
		m_nVerticalOrder+1, dG, m_dDiffReconsPolyNode);

	FluxReconstructionFunction::GetDerivatives(
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
}

///////////////////////////////////////////////////////////////////////////////

void GridGLL::ComputeVorticityDivergence(
	int iDataIndex
) {
	// Compute vorticity on all grid patches
	Grid::ComputeVorticityDivergence(iDataIndex);

	// Apply DSS
	ApplyDSS(0, DataType_Vorticity);
	ApplyDSS(0, DataType_Divergence);
}

///////////////////////////////////////////////////////////////////////////////

void GridGLL::InterpolateNodeToFEEdges(
	const double * dDataNode,
	bool fZeroBoundaries
) const {
	const int Left = 0;
	const int Right = 1;

	// Number of radial elements
	int nRElements = GetRElements();

	int nFiniteElements = nRElements / m_nVerticalOrder;

	// Interpolate to interfaces (left and right)
	m_dStateFEEdge.Zero();

	for (int k = 0; k < nRElements; k++) {
		int a = k / m_nVerticalOrder;
		int m = k % m_nVerticalOrder;

		// Apply node value to left side of interface
		m_dStateFEEdge[a+1][Left] +=
			m_dInterpNodeToREdge[m_nVerticalOrder][m] * dDataNode[k];

		// Apply node value to right side of interface
		m_dStateFEEdge[a][Right] +=
			m_dInterpNodeToREdge[0][m] * dDataNode[k];
	}

	// Calculate average interpolant
	for (int a = 1; a < nFiniteElements; a++) {
		double dAvg = 0.5 * (m_dStateFEEdge[a][0] + m_dStateFEEdge[a][1]);

		m_dStateFEEdge[a][Left] -= dAvg;
		m_dStateFEEdge[a][Right] -= dAvg;
	}

#pragma message "Understand why this works for free boundaries"
	// Ignore contributions due to upper and lower boundary
	// NOTE: This needs to be changed for W to enforce boundary conditions
	if (fZeroBoundaries) {

	} else if (nFiniteElements == 1) {
		m_dStateFEEdge[0][Right] = 0.0;
		m_dStateFEEdge[nFiniteElements][Left] = 0.0;

	} else {
		m_dStateFEEdge[0][Right] =
			- m_dStateFEEdge[1][Left];
		m_dStateFEEdge[nFiniteElements][Left] =
			- m_dStateFEEdge[nFiniteElements-1][Right];
	}
}

///////////////////////////////////////////////////////////////////////////////

void GridGLL::InterpolateNodeToREdge(
	const double * dDataNode,
	const double * dDataRefNode,
	double * dDataREdge,
	const double * dDataRefREdge,
	bool fZeroBoundaries
) const {
	// Number of radial elements
	int nRElements = GetRElements();

	int nFiniteElements = nRElements / m_nVerticalOrder;

	// Zero the memory
	memset(dDataREdge, 0, (nRElements+1) * sizeof(double));

	// Loop over all nodes and apply the value to all interfaces
	for (int k = 0; k < nRElements; k++) {
		int a = k / m_nVerticalOrder;
		int m = k % m_nVerticalOrder;
		int lBegin = a * m_nVerticalOrder;

		// Apply node value to interface
		for (int l = 0; l <= m_nVerticalOrder; l++) {
			dDataREdge[lBegin + l] +=
				m_dInterpNodeToREdge[l][m]
					* (dDataNode[k] - dDataRefNode[k]);
		}
	}

	// Scale interior element edges
	for (int a = 1; a < nFiniteElements; a++) {
		dDataREdge[a * m_nVerticalOrder] *= 0.5;
	}
	for (int k = 0; k <= nRElements; k++) {
		dDataREdge[k] += dDataRefREdge[k];
	}

	// Override boundary values if zero
	if (fZeroBoundaries) {
		dDataREdge[0] = 0.0;
		dDataREdge[nRElements] = 0.0;
	}
}

///////////////////////////////////////////////////////////////////////////////

void GridGLL::InterpolateREdgeToNode(
	const double * dDataREdge,
	const double * dDataRefREdge,
	double * dDataNode,
	const double * dDataRefNode
) const {
	// Number of radial elements
	int nRElements = GetRElements();

	// Zero the memory
	memset(dDataNode, 0, nRElements * sizeof(double));

	// Loop over all nodes
	for (int k = 0; k < nRElements; k++) {
		int a = k / m_nVerticalOrder;
		int m = k % m_nVerticalOrder;
		int lBegin = a * m_nVerticalOrder;

		// Apply interface values to nodes
		for (int l = 0; l <= m_nVerticalOrder; l++) {
			dDataNode[k] +=
				m_dInterpREdgeToNode[m][l] *
					(dDataREdge[lBegin + l] - dDataRefREdge[lBegin + l]);
		}
	}

	// Restore the reference state
	for (int k = 0; k < nRElements; k++) {
		dDataNode[k] += dDataRefNode[k];
	}
}

///////////////////////////////////////////////////////////////////////////////

void GridGLL::DifferentiateNodeToNode(
	const double * dDataNode,
	double * dDiffNode,
	bool fZeroBoundaries
) const {
	const int Left = 0;
	const int Right = 1;

	// Number of radial elements
	int nRElements = GetRElements();

	// Zero the output
	memset(dDiffNode, 0, nRElements * sizeof(double));

	// Interpolate nodes to finite-element edges
	InterpolateNodeToFEEdges(dDataNode, fZeroBoundaries);

	// Calculate derivative at each node
	for (int k = 0; k < nRElements; k++) {
		int a = k / m_nVerticalOrder;
		int m = k % m_nVerticalOrder;

		int lBegin = a * m_nVerticalOrder;

		// Calculate derivatives due to internal bits
		for (int l = 0; l < m_nVerticalOrder; l++) {
			dDiffNode[k] += m_dDiffNodeToNode[m][l] * dDataNode[lBegin+l];
		}

		// Calculate derivatives due to interfaces
		dDiffNode[k] -=
			m_dDiffReconsPolyNode[m]
			* m_dStateFEEdge[a+1][Left];
		dDiffNode[k] +=
			m_dDiffReconsPolyNode[m_nVerticalOrder - m - 1]
			* m_dStateFEEdge[a][Right];
	}
}

///////////////////////////////////////////////////////////////////////////////

void GridGLL::DifferentiateNodeToREdge(
	const double * dDataNode,
	double * dDiffREdge,
	bool fZeroBoundaries
) const {
	const int Left = 0;
	const int Right = 1;

	// Number of radial elements
	int nRElements = GetRElements();

	int nFiniteElements = nRElements / m_nVerticalOrder;

	// Zero the output
	memset(dDiffREdge, 0, (nRElements+1) * sizeof(double));

	// Handle the single element case
	if (nFiniteElements == 1) {
		for (int k = 0; k < m_nVerticalOrder; k++) {
		for (int l = 0; l <= m_nVerticalOrder; l++) {
			dDiffREdge[l] +=
				m_dDiffNodeToREdge[l][k]
				* dDataNode[k];
		}
		}

		return;
	}

	if (nFiniteElements == 2) {
		_EXCEPTIONT("UNIMPLEMENTED: Still working on this...");
	}

	if (fZeroBoundaries) {
		_EXCEPTIONT("UNIMPLEMENTED: Still working on transitioning this...");
	}

	// Interpolate nodes to finite-element edges
	InterpolateNodeToFEEdges(dDataNode, fZeroBoundaries);

	// Calculate derivatives at interfaces
	for (int k = 0; k < nRElements; k++) {
		int a = k / m_nVerticalOrder;
		int m = k % m_nVerticalOrder;

		int lBegin = a * m_nVerticalOrder;

		// Push value of each node onto interface derivatives
		for (int l = 0; l <= m_nVerticalOrder; l++) {
			dDiffREdge[lBegin+l] +=
				m_dDiffNodeToREdge[l][m]
				* dDataNode[k];
		}
	}

	// Calculate derivatives due to interfaces
	for (int a = 0; a < nFiniteElements; a++) {
		int lBegin = a * m_nVerticalOrder;

		for (int l = 0; l <= m_nVerticalOrder; l++) {
			dDiffREdge[lBegin+l] -=
				m_dDiffReconsPolyREdge[l]
				* m_dStateFEEdge[a+1][Left];
			dDiffREdge[lBegin+l] +=
				m_dDiffReconsPolyREdge[m_nVerticalOrder - l]
				* m_dStateFEEdge[a][Right];
		}
	}

	// Scale interior finite element edges
	for (int a = 1; a < nFiniteElements; a++) {
		dDiffREdge[a * m_nVerticalOrder] *= 0.5;
	}
/*
#pragma message "Doesn't work with ReconstructionFunctionType = 1"
	if (ParamFluxReconstructionType != 2) {
		_EXCEPTIONT("UNIMPLEMENTED");
	}

	// Compute interior derivatives
	{
		int kBegin = m_nVerticalOrder;
		int kLast = (nFiniteElements - 1) * m_nVerticalOrder;

		for (int k = kBegin; k <= kLast; k++) {
			int a = k / m_nVerticalOrder;
			int l = k % m_nVerticalOrder;
			if (k == kLast) {
				a--;
				l = m_nVerticalOrder;
			}

			int lPrev = (a-1) * m_nVerticalOrder;

			for (int m = 0; m < 3 * m_nVerticalOrder; m++) {
				dDiffREdge[k] +=
					m_dDiffNodeToREdgeAmal[l][m]
					* dDataNode[lPrev + m];
			}
		}
	}

	// Compute derivatives at left boundary
	for (int l = 0; l < m_nVerticalOrder; l++) {
	for (int m = 0; m < 2 * m_nVerticalOrder; m++) {
		dDiffREdge[l] +=
			m_dDiffNodeToREdgeLeft[l][m]
			* dDataNode[m];
	}
	}

	// Compute derivatives at right boundary
	{
		int lBegin = (nFiniteElements - 1) * m_nVerticalOrder;
		int lEnd = lBegin + m_nVerticalOrder;
		for (int l = 1; l <= m_nVerticalOrder; l++) {
		for (int m = 0; m < 2 * m_nVerticalOrder; m++) {
			dDiffREdge[lBegin + l] +=
				m_dDiffNodeToREdgeRight[l][m]
				* dDataNode[lBegin - m_nVerticalOrder + m];
		}
		}
	}
*/
}

///////////////////////////////////////////////////////////////////////////////

void GridGLL::DifferentiateREdgeToNode(
	const double * dDataREdge,
	double * dDiffNode
) const {
	// Number of radial elements
	int nRElements = GetRElements();

	int nFiniteElements = nRElements / m_nVerticalOrder;

	// Zero the data
	memset(dDiffNode, 0, nRElements * sizeof(double));

	// Loop through all nodes
	for (int k = 0; k < nRElements; k++) {

		// Differentiate from neighboring interface values
		int a = k / m_nVerticalOrder;
		int m = k % m_nVerticalOrder;
		int lBegin = a * m_nVerticalOrder;

		for (int l = 0; l <= m_nVerticalOrder; l++) {
			dDiffNode[k] +=
				  m_dDiffREdgeToNode[m][l]
				* dDataREdge[lBegin + l];
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void GridGLL::DifferentiateREdgeToREdge(
	const double * dDataREdge,
	double * dDiffREdge
) const {
	// Number of radial elements
	int nRElements = GetRElements();

	int nFiniteElements = nRElements / m_nVerticalOrder;

	// Zero the data
	memset(dDiffREdge, 0, (nRElements+1) * sizeof(double));

	// Apply all interfaces values to all interfaces within element
	for (int a = 0; a < nFiniteElements; a++) {
	for (int l = 0; l <= m_nVerticalOrder; l++) {

		int lBegin = a * m_nVerticalOrder;

		for (int m = 0; m <= m_nVerticalOrder; m++) {
			dDiffREdge[lBegin + l] +=
				  m_dDiffREdgeToREdge[l][m]
				* dDataREdge[lBegin + m];
		}
	}
	}

	// Halve interior element interface values
	for (int a = 1; a < nFiniteElements; a++) {
		dDiffREdge[a * m_nVerticalOrder] *= 0.5;
	}
}

///////////////////////////////////////////////////////////////////////////////

