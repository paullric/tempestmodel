///////////////////////////////////////////////////////////////////////////////
///
///	\file    GridCartesianGLL.h
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

#ifndef _GRIDCARTESIANGLL_H_
#define _GRIDCARTESIANGLL_H_

#include "GridGLL.h"

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		The CubedSphere grid with degrees of freedom at Gauss-Lobatto-Legendre
///		quadrature nodes.
///	</summary>
class GridCartesianGLL : public GridGLL {
    
private:
    ///	<summary>
	///		Dimension of the grid - private to cartesian grids.
	///	</summary>
    double m_dGDim[];

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	GridCartesianGLL(
		const Model & model,
		int nBaseResolution,
		int nRefinementRatio,
		int nHorizontalOrder,
		int nVerticalOrder,
		int nRElements,
        double nGDim[]
	);

	///	<summary>
	///		Initialize grid patches.
	///	</summary>
	virtual void Initialize();
    
public:
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
	virtual void ConvertReferenceToPatchCoord(
		const DataVector<double> & dXReference,
		const DataVector<double> & dYReference,
		DataVector<double> & dAlpha,
		DataVector<double> & dBeta,
		DataVector<int> & iPanel
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
