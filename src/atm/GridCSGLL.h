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

#ifndef _GRIDCSGLL_H_
#define _GRIDCSGLL_H_

#include "Grid.h"

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		The CubedSphere grid with degrees of freedom at Gauss-Lobatto-Legendre
///		quadrature nodes.
///	</summary>
class GridCSGLL : public Grid {

public:
	///	<summary>
	///		Get the type of grid.
	///	</summary>
	virtual Grid::Type GetType() const {
		return Grid::CubedSphereGLLGrid;
	}

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	GridCSGLL(
		const Model & model,
		int nBaseResolution,
		int nRefinementRatio,
		int nOrder,
		int nVerticalOrder,
		int nRElements
	);

	///	<summary>
	///		Initialize grid patches.
	///	</summary>
	virtual void Initialize();

	///	<summary>
	///		Initialize the vertical coordinate.
	///	</summary>
	virtual void InitializeVerticalCoordinate(
		const GridSpacing & aGridSpacing
	);

public:
	///	<summary>
	///		Convert an array of coordinate variables to coordinates on the
	///		reference grid (RLL on the sphere)
	///	</summary>
	virtual void ConvertReferenceToABP(
		const DataVector<double> & dXReference,
		const DataVector<double> & dYReference,
		DataVector<double> & dAlpha,
		DataVector<double> & dBeta,
		DataVector<int> & iPanel
	) const;

public:
	///	<summary>
	///		Apply the DSS operation on the grid.
	///	</summary>
	void ApplyDSS(
		int iDataUpdate,
		DataType eDataType = DataType_State
	);

	///	<summary>
	///		Compute vorticity on the grid.
	///	</summary>
	virtual void ComputeVorticity(
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

