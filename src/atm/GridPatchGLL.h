///////////////////////////////////////////////////////////////////////////////
///
///	\file    GridPatchGLL.h
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

#ifndef _GRIDPATCHGLL_H_
#define _GRIDPATCHGLL_H_

#include "GridPatch.h"
#include "GridGLL.h"

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A grid patch with degrees of freedom at Gauss-Lobatto-Legendre
///		quadrature nodes.
///	</summary>
class GridPatchGLL : public GridPatch {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	GridPatchGLL(
		Grid & grid,
		int ixPatch,
		const PatchBox & box,
		int nHorizontalOrder,
		int nVerticalOrder = 1
	);

protected:
	///	<summary>
	///		Initialize coordinate data on patch.
	///	</summary>
	virtual void InitializeCoordinateData();

public:
	///	<summary>
	///		Interpolate data vertically from nodes to interfaces.
	///	</summary>
	virtual void InterpolateNodeToREdge(
		int iVar,
		int iDataIndex
	);

	///	<summary>
	///		Interpolate data vertically from interfaces to nodes.
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
	) {
	}

public:
	///	<summary>
	///		Get the number of finite elements in the alpha direction.
	///	</summary>
	inline int GetElementCountA() const {
		return m_nElementCountA;
	}

	///	<summary>
	///		Get the number of finite elements in the beta direction.
	///	</summary>
	inline int GetElementCountB() const {
		return m_nElementCountB;
	}

	///	<summary>
	///		Get the element grid spacing in the alpha direction.
	///	</summary>
	inline double GetElementDeltaA() const {
		return m_dElementDeltaA;
	}

	///	<summary>
	///		Get the element grid spacing in the beta direction.
	///	</summary>
	inline double GetElementDeltaB() const {
		return m_dElementDeltaB;
	}

protected:
	///	<summary>
	///		Order of accuracy of this patch.
	///	</summary>
	int m_nHorizontalOrder;

	///	<summary>
	///		Vertical order of accuracy of this patch.
	///	</summary>
	int m_nVerticalOrder;

	///	<summary>
	///		Number of finite elements in the alpha direction.
	///	</summary>
	int m_nElementCountA;

	///	<summary>
	///		Number of finite elements in the beta direction.
	///	</summary>
	int m_nElementCountB;

	///	<summary>
	///		Element grid spacing in the alpha direction.
	///	</summary>
	double m_dElementDeltaA;

	///	<summary>
	///		Element grid spacing in the beta direction.
	///	</summary>
	double m_dElementDeltaB;
};

///////////////////////////////////////////////////////////////////////////////

#endif

