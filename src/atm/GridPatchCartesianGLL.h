///////////////////////////////////////////////////////////////////////////////
///
///	\file    GridPatchCartesianGLL.h
///	\author  Paul Ullrich, Jorge Guerra
///	\version September 26, 2013
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

#ifndef _GRIDPATCHCARTESIANGLL_H_
#define _GRIDPATCHCARTESIANGLL_H_

#include "Grid.h"
#include "GridCartesianGLL.h"
#include "GridPatchGLL.h"
#include "DataArray3D.h"

///////////////////////////////////////////////////////////////////////////////

class Time;
class GridCartesianGLL;

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A grid patch on the CubedSphere grid with degrees of freedom at
///		Gauss-Lobatto-Legendre quadrature nodes.
///	</summary>
class GridPatchCartesianGLL : public GridPatchGLL {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	GridPatchCartesianGLL(
		GridCartesianGLL & grid,
		int ixPatch,
		const PatchBox & box,
		int nHorizontalOrder,
		int nVerticalOrder
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
	///		Initialize geometric source term coefficients and topographical
	///		data from a TestCase.
	///	</summary>
	virtual void EvaluateGeometricTerms();

public:
	///	<summary>
	///		Initialize this grid from a TestCase.
	///	</summary>
	void EvaluateTestCase(
		const TestCase & test,
		const Time & time,
		int iDataInstance = 0
	);

public:
	///	<summary>
	///		Apply boundary conditions to this patch.
	///	</summary>
	virtual void ApplyBoundaryConditions(
		int iDataUpdate,
		DataType eDataType,
		int iAlphaBCPatch
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
		DataType eDataType,
		const DataArray1D<double> & dREta,
		const DataArray1D<double> & dAlpha,
		const DataArray1D<double> & dBeta,
		const DataArray1D<int> & iPatch,
		DataArray3D<double> & dInterpData,
		DataLocation eOnlyVariablesAt = DataLocation_None,
		bool fIncludeReferenceState = true,
		bool fConvertToPrimitive = true
	);
/*
public:

	///	<summary>
	///		Get the maximum of the topography field
	///	</summary>
	double GetMaximumTopo() {
		return m_dTopoHeight;
	}

	///	<summary>
	///		Get the reference latitude for f-plane or b-plane models
	///	</summary>
	double GetTopoScaleHeight() {
		return m_dSL;
	}

public:
	///	<summary>
	///		Maximum height of a topography feature
	///	</summary>
	double m_dTopoHeight;

	///	<summary>
	///		Scale height for exponential decay of topography
	///	</summary>
	double m_dSL;
*/
};

///////////////////////////////////////////////////////////////////////////////

#endif

