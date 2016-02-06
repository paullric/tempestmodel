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
#include "DataArray3D.h"

///////////////////////////////////////////////////////////////////////////////

class Time;
class GridCSGLL;

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
	virtual void InitializeDataLocal(
		bool fAllocateGeometric = true,
		bool fAllocateActiveState = true,
		bool fAllocateBufferState = true,
		bool fAllocateAuxiliary = true
	);

protected:
	///	<summary>
	///		Initialize coordinate data on patch.
	///	</summary>
	virtual void InitializeCoordinateData();

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
		const DataArray3D<double> & dataUa,
		const DataArray3D<double> & dataUb
	);

	///	<summary>
	///		Compute the divergence and vorticity of the vector velocity
	///		on the grid.
	///	</summary>
	virtual void ComputeVorticityDivergence(
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

	///	<summary>
	///		Transform derivatives of the topography from other panels
	///		to this panel's coordinate system.
	///	</summary>
	virtual void TransformTopographyDeriv();

public:
	///	<summary>
	///		Linearly interpolate data to the specified point.
	///	</summary>
	virtual void InterpolateData(
		const DataArray1D<double> & dAlpha,
		const DataArray1D<double> & dBeta,
		const DataArray1D<int> & iPatch,
		DataType eDataType,
		DataLocation eDataLocation,
		bool fInterpAllVariables,
		DataArray3D<double> & dInterpData,
		bool fIncludeReferenceState = true,
		bool fConvertToPrimitive = true
	);

private:
	///	<summary>
	///		Buffer space for computing contravariant velocities in
	///		ComputeCurlAndDiv().
	///	</summary>
	DataArray3D<double> m_dBufferConU;
};

///////////////////////////////////////////////////////////////////////////////

#endif

