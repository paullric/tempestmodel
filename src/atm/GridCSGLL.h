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
#include "GridPatch.h"

///////////////////////////////////////////////////////////////////////////////

class GridCSGLL;

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A grid patch on the CubedSphere grid with degrees of freedom at
///		Gauss-Lobatto-Legendre quadrature nodes.
///	</summary>
class GridPatchCSGLL : public GridPatch {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	GridPatchCSGLL(
		GridCSGLL & grid,
		int ixPatch,
		const PatchBox & box,
		int nOrder,
		int nVerticalOrder = 1
	);

public:
	///	<summary>
	///		Initialize the patch data.
	///	</summary>
	virtual void InitializeDataLocal();

private:
	///	<summary>
	///		Initialize geometric source term coefficients and topographical
	///		data from a TestCase.
	///	</summary>
	void EvaluateGeometricTerms(
		const TestCase & test
	);

public:
	///	<summary>
	///		Initialize this grid from a TestCase.
	///	</summary>
	void EvaluateTestCase(
		const TestCase & test,
		double dTime = 0.0,
		int iDataInstance = 0
	);

public:
	///	<summary>
	///		Compute vorticity on the grid.
	///	</summary>
	virtual void ComputeVorticity(
		int iDataIndex
	);

public:
	///	<summary>
	///		Interpolate data vertically from Nodes to REdges.
	///	</summary>
	virtual void InterpolateNodeToREdge(
		int iVar,
		int iDataIndex
	);

	///	<summary>
	///		Interpolate data vertically from REdges to Nodes.
	///	</summary>
	virtual void InterpolateREdgeToNode(
		int iVar,
		int iDataIndex
	);

public:
	///	<summary>
	///		Linearly interpolate data to the specified point.
	///	</summary>
	virtual void InterpolateData(
		const DataVector<double> & dAlpha,
		const DataVector<double> & dBeta,
		const DataVector<int> & iPanel,
		DataType eDataType,
		DataMatrix3D<double> & dInterpData,
		bool fIncludeReferenceState = true,
		bool fConvertToPrimitive = true
	);

private:
	///	<summary>
	///		Order of accuracy of this patch.
	///	</summary>
	int m_nHorizontalOrder;

	///	<summary>
	///		Vertical order of accuracy of this patch.
	///	</summary>
	int m_nVerticalOrder;
};

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
	int GetOrder() const {
		return m_nHorizontalOrder;
	}

	///	<summary>
	///		Get the order of accuracy in the vertical.
	///	</summary>
	int GetVerticalOrder() const {
		return m_nVerticalOrder;
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

