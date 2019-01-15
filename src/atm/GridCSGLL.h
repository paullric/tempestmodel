///////////////////////////////////////////////////////////////////////////////
///
///	\file    GridCSGLL.h
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

#ifndef _GRIDCSGLL_H_
#define _GRIDCSGLL_H_

#include "GridGLL.h"

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		The CubedSphere grid with degrees of freedom at Gauss-Lobatto-Legendre
///		quadrature nodes.
///	</summary>
class GridCSGLL : public GridGLL {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	GridCSGLL(
		Model & model
	);

	///	<summary>
	///		Define the parameters for the Grid.
	///	</summary>
	virtual void DefineParameters();

	///	<summary>
	///		Set the parameters for the GridGLL.
	///	</summary>
	virtual void SetParameters(
		int nRElements,
		int nMaxPatchCount,
		int nBaseResolution,
		int nRefinementRatio,
		int nHorizontalOrder,
		int nVerticalOrder,
		VerticalDiscretization eVerticalDiscretization,
		VerticalStaggering eVerticalStaggering
	);

public:
	///	<summary>
	///		Build the default patch layout.
	///	</summary>
	virtual void ApplyDefaultPatchLayout(
		int nPatchCount
	);

	///	<summary>
	///		Return a pointer to a new GridPatchCSGLL.
	///	</summary>
	virtual GridPatch * NewPatch(
		int ixPatch
	);

	///	<summary>
	///		Get the bounds on the reference grid.
	///	</summary>
	virtual void GetReferenceGridBounds(
		double & dX0,
		double & dX1,
		double & dY0,
		double & dY1
	);

	///	<summary>
	///		Convert an array of coordinate variables to coordinates on the
	///		reference grid (RLL on the sphere)
	///	</summary>
	virtual void ConvertReferenceToPatchCoord(
		const DataArray1D<double> & dXReference,
		const DataArray1D<double> & dYReference,
		DataArray1D<double> & dAlpha,
		DataArray1D<double> & dBeta,
		DataArray1D<int> & iPanel
	) const;

	///	<summary>
	///		Get the patch and coordinate index for the specified node.
	///	</summary>
	virtual void GetPatchFromCoordinateIndex(
		int iRefinementLevel,
		const DataArray1D<int> & vecIxA,
		const DataArray1D<int> & vecIxB,
		const DataArray1D<int> & vecPanel,
		DataArray1D<int> & vecPatchIndex,
		int nVectorLength = (-1)
	);

	///	<summary>
	///		Get the relation between coordinate vectors across panel
	///		boundaries.
	///	</summary>
	virtual void GetOpposingDirection(
		int ixPanelSrc,
		int ixPanelDest,
		Direction dir,
		Direction & dirOpposing,
		bool & fSwitchParallel,
		bool & fSwitchPerpendicular
	) const;

public:
	///	<summary>
	///		Apply the direct stiffness summation (DSS) operation on the grid.
	///	</summary>
	virtual void ApplyDSS(
		int iDataUpdate,
		DataType eDataType = DataType_State
	);
};

///////////////////////////////////////////////////////////////////////////////

#endif

