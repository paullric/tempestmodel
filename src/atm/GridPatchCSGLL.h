///////////////////////////////////////////////////////////////////////////////
///
///	\file    GridPatchCSGLL.h
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

#ifndef _GRIDPATCHCSGLL_H_
#define _GRIDPATCHCSGLL_H_

#include "Grid.h"
#include "GridPatch.h"

///////////////////////////////////////////////////////////////////////////////

class GridCSGLL;
class GridData3D;

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
	///		Compute the radial component of the curl on the grid given two
	///		contravariant vector fields.
	///	</summary>
	virtual void ComputeCurlR(
		const GridData3D & dataUa,
		const GridData3D & dataUb,
		GridData3D & dataCurlUr
	) const;

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

#endif

