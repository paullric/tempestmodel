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
#include "CubedSphereTrans.h"
#include "PolynomialInterp.h"
#include "GaussQuadrature.h"
#include "GaussLobattoQuadrature.h"

#include "Announce.h"
#include "MathHelper.h"

#include <cmath>

///////////////////////////////////////////////////////////////////////////////

GridGLL::GridGLL(
	const Model & model,
	int nBaseResolution,
	int nRefinementRatio,
	int nHorizontalOrder,
	int nVerticalOrder,
	int nRElements
) :
	// Call up the stack
	Grid::Grid(
		model,
		nBaseResolution,
		nBaseResolution,
		nRefinementRatio,
		nRElements),
	m_nHorizontalOrder(nHorizontalOrder),
	m_nVerticalOrder(nVerticalOrder)
{
	// Get quadrature points for Gauss-Lobatto quadrature
	DataVector<double> dG;
	DataVector<double> dW;

	GaussLobattoQuadrature::GetPoints(m_nHorizontalOrder, 0.0, 1.0, dG, dW);

	// Derivatives of the 1D basis functions at each point on the reference
	// element [0, 1]
	m_dDxBasis1D.Initialize(m_nHorizontalOrder, m_nHorizontalOrder);
	m_dStiffness1D.Initialize(m_nHorizontalOrder, m_nHorizontalOrder);

	DataVector<double> dCoeffs;
	dCoeffs.Initialize(m_nHorizontalOrder);

	for (int i = 0; i < m_nHorizontalOrder; i++) {
		PolynomialInterp::DiffLagrangianPolynomialCoeffs(
			m_nHorizontalOrder, dG, dCoeffs, dG[i]);

		for (int m = 0; m < m_nHorizontalOrder; m++) {
			m_dDxBasis1D[m][i] = dCoeffs[m];

			m_dStiffness1D[m][i] = m_dDxBasis1D[m][i] * dW[i] / dW[m];
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

