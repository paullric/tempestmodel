///////////////////////////////////////////////////////////////////////////////
///
///	\file    GridCartesianGLL.h
///	\author  Paul Ullrich, Jorge Guerra
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

#ifndef _GRIDCARTESIANGLL_H_
#define _GRIDCARTESIANGLL_H_

#include "GridGLL.h"

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		The CubedSphere grid with degrees of freedom at Gauss-Lobatto-Legendre
///		quadrature nodes.
///	</summary>
class GridCartesianGLL : public GridGLL {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	GridCartesianGLL(
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
		int nABaseResolution,
		int nBBaseResolution,
		int nRefinementRatio,
		int nHorizontalOrder,
		int nVerticalOrder,
		double dGDim[],
		double dRefLat,
		int iLatBC[],
		bool fCartesianXZ,
		VerticalDiscretization eVerticalDiscretization,
		VerticalStaggering eVerticalStaggering
	);

public:
	///	<summary
	///		Get the bounds on the reference grid.
	///	</summary>
	virtual void GetReferenceGridBounds(
		double & dX0,
		double & dX1,
		double & dY0,
		double & dY1
	);
	
	///	<summary>
	///		Build the default patch layout.
	///	</summary>
	virtual void ApplyDefaultPatchLayout(
		int nPatchCount
	);

	///	<summary>
	///		Return a pointer to a new GridPatchCartesianGLL.
	///	</summary>
	virtual GridPatch * NewPatch(
		int ixPatch
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
	///		Apply the boundary conditions.
	///	</summary>
	void ApplyBoundaryConditions(
		int iDataUpdate,
		DataType eDataType
	);

	///	<summary>
	///		Apply the direct stiffness summation (DSS) operation on the grid.
	///	</summary>
	virtual void ApplyDSS(
		int iDataUpdate,
		DataType eDataType = DataType_State
	);

public:
	///	<summary>
	///		Get the reference latitude for f-plane or b-plane models
	///	</summary>
	double GetReferenceLatitude() const {
		return m_dRefLat;
	}

	///	<summary>
	///		Get the minimum X dimension
	///	</summary>
	double GetMinimumX() const {
		return m_dGDim[0];
	}

	///	<summary>
	///		Get the maximum X dimension
	///	</summary>
	double GetMaximumX() const {
		return m_dGDim[1];
	}

	///	<summary>
	///		Get the minimum Y dimension
	///	</summary>
	double GetMinimumY() const {
		return m_dGDim[2];
	}

	///	<summary>
	///		Get the maximum Y dimension
	///	</summary>
	double GetMaximumY() const {
		return m_dGDim[3];
	}

	///	<summary>
	///		Get the maximum topography height
	///	</summary>
	double GetMaxTopoHeight() const {
		return m_dTopoHeight;
	}

	///	<summary>
	///		Get the reference latitude for f-plane or b-plane models
	///	</summary>
	double GetTopoScaleHeight() {
		return m_dSL;
	}

	///	<summary>
	///		Get the 2D Cartesian flag
	///	</summary>
	//bool GetIsCartesianXZ() {
	//	return m_fCartesianXZ;
	//}

public:
	///	<summary>
	///		2D Cartesian flag
	///	</summary>
	//bool m_fCartesianXZ;

	///	<summary>
	///		Maximum topography height
	///	</summary>
	double m_dTopoHeight;

	///	<summary>
	///		Scale height for exponential decay of topography
	///	</summary>
	double m_dSL;

	///	<summary>
	///		Dimension of the grid - private to cartesian grids.
	///	</summary>
	DataArray1D<double> m_dGDim;

	///	<summary>
	///		Referece latitude (for beta plane cases)
	///	</summary>
	DataStruct<double> m_dRefLat;
};

///////////////////////////////////////////////////////////////////////////////

#endif
