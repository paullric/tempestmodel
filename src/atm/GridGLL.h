///////////////////////////////////////////////////////////////////////////////
///
///	\file    GridCSGLL.h
///	\author  Paul Ullrich
///	\version February 25, 2013
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

#ifndef _GRIDGLL_H_
#define _GRIDGLL_H_

#include "Grid.h"

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		The CubedSphere grid with degrees of freedom at Gauss-Lobatto-Legendre
///		quadrature nodes.
///	</summary>
class GridGLL : public Grid {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	GridGLL(
		const Model & model,
		int nBaseResolution,
		int nRefinementRatio,
		int nHorizontalOrder,
		int nVerticalOrder,
		int nRElements
	);

public:
	///	<summary>
	///		Apply the direct stiffness summation (DSS) operation on the grid.
	///	</summary>
	virtual void ApplyDSS(
		int iDataUpdate,
		DataType eDataType = DataType_State
	) {
		_EXCEPTIONT("Unimplemented");
	}

	///	<summary>
	///		Compute vorticity on the grid.
	///	</summary>
	virtual void ComputeVorticityDivergence(
		int iDataIndex
	);

public:
	///	<summary>
	///		Get the order of accuracy of the method.
	///	</summary>
	int GetHorizontalOrder() const {
		return m_nHorizontalOrder;
	}

	///	<summary>
	///		Get the order of accuracy in the vertical.
	///	</summary>
	int GetVerticalOrder() const {
		return m_nVerticalOrder;
	}

	///	<summary>
	///		Get the derivatives of the basis functions at nodal points on the
	///		horizontal reference element.
	///	</summary>
	const DataMatrix<double> & GetDxBasis1D() const {
		return m_dDxBasis1D;
	}

	///	<summary>
	///		Get the stiffness matrix coefficients at nodal points on the
	///		horizontal reference element.
	///	</summary>
	const DataMatrix<double> & GetStiffness1D() const {
		return m_dStiffness1D;
	}

	///	<summary>
	///		Get the interpolation coefficients from levels to interfaces.
	///	</summary>
	const DataMatrix<double> & GetInterpNodeToREdge() const {
		return m_dInterpNodeToREdge;
	}

	///	<summary>
	///		Get the interpolation coefficients from interfaces to levels.
	///	</summary>
	const DataMatrix<double> & GetInterpREdgeToNode() const {
		return m_dInterpREdgeToNode;
	}

protected:
	///	<summary>
	///		Order of accuracy of the method (number of nodes per element).
	///	</summary>
	int m_nHorizontalOrder;

	///	<summary>
	///		Order of accuracy in the vertical.
	///	</summary>
	int m_nVerticalOrder;

	///	<summary>
	///		Derivatives of the basis functions at nodal points on the
	///		horizontal reference element.
	///	</summary>
	DataMatrix<double> m_dDxBasis1D;

	///	<summary>
	///		Stiffness matrix coefficients at nodal points on the
	///		horizontal reference element.
	///	</summary>
	DataMatrix<double> m_dStiffness1D;

	///	<summary>
	///		Interpolation coefficients from levels to interfaces.
	///	</summary>
	DataMatrix<double> m_dInterpNodeToREdge;

	///	<summary>
	///		Interpolation coefficients from interfaces to levels.
	///	</summary>
	DataMatrix<double> m_dInterpREdgeToNode;

};

///////////////////////////////////////////////////////////////////////////////

#endif

