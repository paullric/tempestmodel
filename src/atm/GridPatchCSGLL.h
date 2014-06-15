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
#include "GridPatchGLL.h"

///////////////////////////////////////////////////////////////////////////////

class Time;
class GridCSGLL;
class GridData3D;

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A grid patch on the CubedSphere grid with degrees of freedom at
///		Gauss-Lobatto-Legendre quadrature nodes.
///	</summary>
class GridPatchCSGLL : public GridPatchGLL {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	GridPatchCSGLL(
		GridCSGLL & grid,
		int ixPatch,
		const PatchBox & box,
		int nHorizontalOrder,
		int nVerticalOrder = 1
	);

public:
	///	<summary>
	///		Initialize the patch data.
	///	</summary>
	virtual void InitializeDataLocal();

private:
	///	<summary>
	///		Initialize topographical data from a TestCase.
	///	</summary>
	void EvaluateTopography(
		const TestCase & test
	);

	///	<summary>
	///		Initialize geometric terms.
	///	</summary>
	virtual void EvaluateGeometricTerms();

	///	<summary>
	///		Initialize state and tracer data from a TestCase.  Also adjust
	///		geometric quantities that are dependent on the TestCase.
	///	</summary>
	virtual void EvaluateTestCase(
		const TestCase & test,
		const Time & time,
		int iDataInstance = 0
	);

	///	<summary>
	///		Initialize state and tracer data from a TestCase.
	///	</summary>
	virtual void EvaluateTestCase_StateOnly(
		const TestCase & test,
		const Time & time,
		int iDataInstance = 0
	);

public:
	///	<summary>
	///		Compute the radial component of the curl on the grid given two
	///		contravariant vector fields.
	///	</summary>
	virtual void ComputeCurlAndDiv(
		const GridData3D & dataUa,
		const GridData3D & dataUb
	) const;

	///	<summary>
	///		Compute the divergence and vorticity of the vector velocity
	///		on the grid.
	///	</summary>
	virtual void ComputeVorticityDivergence(
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
	///		Transform vectors received from other panels to this panel's
	///		coordinate system.
	///	</summary>
	virtual void TransformHaloVelocities(
		int iDataUpdate
	);

public:
	///	<summary>
	///		Linearly interpolate data to the specified point.
	///	</summary>
	virtual void InterpolateData(
		const DataVector<double> & dAlpha,
		const DataVector<double> & dBeta,
		const DataVector<int> & iPatch,
		DataType eDataType,
		DataLocation eDataLocation,
		bool fInterpAllVariables,
		DataMatrix3D<double> & dInterpData,
		bool fIncludeReferenceState = true,
		bool fConvertToPrimitive = true
	);

};

///////////////////////////////////////////////////////////////////////////////

#endif

