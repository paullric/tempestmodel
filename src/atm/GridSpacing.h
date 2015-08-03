///////////////////////////////////////////////////////////////////////////////
///
///	\file    GridSpacing.h
///	\author  Paul Ullrich
///	\version March 8, 2013
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

#ifndef _GRIDSPACING_H_
#define _GRIDSPACING_H_

#include "DataArray1D.h"

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Interface for spacing of nodal points on the grid.
///	</summary>
class GridSpacing  {
protected:
	///	<summary>
	///		Constructor.
	///	</summary>
	GridSpacing(
		double dDeltaElement,
		double dZeroCoord
	);

	///	<summary>
	///		Virtual destructor.
	///	</summary>
	virtual ~GridSpacing()
	{ }

public:
	///	<summary>
	///		Check if this grid spacing supports the given total number of
	///		nodes along this coordinate axis.
	///	</summary>
	virtual bool DoesNodeCountAgree(int nNodes) const {
		return true;
	}

	///	<summary>
	///		Get the nodal coordinate with the given index.
	///	</summary>
	virtual double GetNode(int ix) const = 0;

	///	<summary>
	///		Get the edge coordinate with the given index.
	///	</summary>
	virtual double GetEdge(int ix) const = 0;

	///	<summary>
	///		Get the normalized area of the given node.
	///	</summary>
	virtual double GetNodeNormArea(int ix) const = 0;

	///	<summary>
	///		Get the normalized area of the given edge.
	///	</summary>
	virtual double GetEdgeNormArea(int ix) const = 0;

public:
	///	<summary>
	///		Get the element width.
	///	</summary>
	double GetDeltaElement() const {
		return m_dDeltaElement;
	}

	///	<summary>
	///		Get the zero point of this coordinate.
	///	</summary>
	double GetZeroCoord() const {
		return m_dZeroCoord;
	}

protected:
	///	<summary>
	///		Grid spacing of each element.
	///	</summary>
	double m_dDeltaElement;

	///	<summary>
	///		Coordinate of zeroth element.
	///	</summary>
	double m_dZeroCoord;
};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Uniform spacing of nodal points on the grid.
///	</summary>
class GridSpacingUniform : public GridSpacing {
public:
	///	<summary>
	///		Constructor.
	///	</summary>
	GridSpacingUniform(
		double dDeltaElement,
		double dZeroCoord
	);

public:
	///	<summary>
	///		Get the nodal coordinate with the given index.
	///	</summary>
	virtual double GetNode(int ix) const;

	///	<summary>
	///		Get the edge coordinate with the given index.
	///	</summary>
	virtual double GetEdge(int ix) const;

	///	<summary>
	///		Get the normalized area of the given node.
	///	</summary>
	virtual double GetNodeNormArea(int ix) const;

	///	<summary>
	///		Get the normalized area of the given edge.
	///	</summary>
	virtual double GetEdgeNormArea(int ix) const;
};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Non-uniform Gauss-Lobatto spacing of nodal points on the grid.  Both
///		nodal points and edges are at Gauss-Lobatto points.
///	</summary>
class GridSpacingGaussLobatto : public GridSpacing {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	GridSpacingGaussLobatto(
		double dDeltaElement,
		double dZeroCoord,
		int nOrder
	);

public:
	///	<summary>
	///		Check if this grid spacing supports the given total number of
	///		nodes along this coordinate axis.
	///	</summary>
	virtual bool DoesNodeCountAgree(int nNodes) const {
		if (((nNodes - 1) % (m_nOrder - 1)) == 0) {
			return true;
		}
		return false;
	}

	///	<summary>
	///		Get the nodal coordinate with the given index.
	///	</summary>
	virtual double GetNode(int ix) const;

	///	<summary>
	///		Get the edge coordinate with the given index.
	///	</summary>
	virtual double GetEdge(int ix) const;

	///	<summary>
	///		Get the normalized area of the given node.
	///	</summary>
	virtual double GetNodeNormArea(int ix) const;

	///	<summary>
	///		Get the normalized area of the given edge.
	///	</summary>
	virtual double GetEdgeNormArea(int ix) const;

protected:
	///	<summary>
	///		Order of accuracy stored in each cell.
	///	</summary>
	int m_nOrder;

	///	<summary>
	///		Gauss-Lobatto points on the reference element.
	///	</summary>
	DataArray1D<double> m_dG;

	///	<summary>
	///		Gauss-Lobatto weights on the reference element.
	///	</summary>
	DataArray1D<double> m_dW;
};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Non-uniform Gauss-Lobatto spacing of nodal points on the grid.  Both
///		nodal points and edges are at Gauss-Lobatto points.
///	</summary>
class GridSpacingGaussLobattoRepeated : public GridSpacing {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	GridSpacingGaussLobattoRepeated(
		double dDeltaElement,
		double dZeroCoord,
		int nOrder
	);

public:
	///	<summary>
	///		Check if this grid spacing supports the given total number of
	///		nodes along this coordinate axis.
	///	</summary>
	virtual bool DoesNodeCountAgree(int nNodes) const {
		if ((nNodes % m_nOrder) == 0) {
			return true;
		}
		return false;
	}

	///	<summary>
	///		Get the nodal coordinate with the given index.
	///	</summary>
	virtual double GetNode(int ix) const;

	///	<summary>
	///		Get the edge coordinate with the given index.
	///	</summary>
	virtual double GetEdge(int ix) const;

	///	<summary>
	///		Get the normalized area of the given node.
	///	</summary>
	virtual double GetNodeNormArea(int ix) const;

	///	<summary>
	///		Get the normalized area of the given edge.
	///	</summary>
	virtual double GetEdgeNormArea(int ix) const;

protected:
	///	<summary>
	///		Order of accuracy stored in each cell.
	///	</summary>
	int m_nOrder;

	///	<summary>
	///		Gauss-Lobatto points on the reference element.
	///	</summary>
	DataArray1D<double> m_dG;

	///	<summary>
	///		Gauss-Lobatto weights on the reference element.
	///	</summary>
	DataArray1D<double> m_dW;
};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Non-uniform mixed Gauss-Lobatto spacing of nodal points on the grid.
///		Nodal points are at Gauss nodes and edges are at Gauss-Lobatto points.
///	</summary>
class GridSpacingMixedGaussLobatto : public GridSpacing {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	GridSpacingMixedGaussLobatto(
		double dDeltaElement,
		double dZeroCoord,
		int nOrder
	);

public:
	///	<summary>
	///		Check if this grid spacing supports the given total number of
	///		nodes along this coordinate axis.
	///	</summary>
	virtual bool DoesNodeCountAgree(int nNodes) const {
		if ((nNodes % m_nOrder) == 0) {
			return true;
		}
		return false;
	}

	///	<summary>
	///		Get the nodal coordinate with the given index.
	///	</summary>
	virtual double GetNode(int ix) const;

	///	<summary>
	///		Get the edge coordinate with the given index.
	///	</summary>
	virtual double GetEdge(int ix) const;

	///	<summary>
	///		Get the normalized area of the given node.
	///	</summary>
	virtual double GetNodeNormArea(int ix) const;

	///	<summary>
	///		Get the normalized area of the given edge.
	///	</summary>
	virtual double GetEdgeNormArea(int ix) const;

protected:
	///	<summary>
	///		Order of accuracy stored in each cell.
	///	</summary>
	int m_nOrder;

	///	<summary>
	///		Gauss points on the reference element.
	///	</summary>
	DataArray1D<double> m_dG;

	///	<summary>
	///		Gauss-Lobatto points on the reference element.
	///	</summary>
	DataArray1D<double> m_dGL;

	///	<summary>
	///		Gauss weights on the reference element.
	///	</summary>
	DataArray1D<double> m_dW;

	///	<summary>
	///		Gauss-Lobatto weights on the reference element.
	///	</summary>
	DataArray1D<double> m_dWL;

};

///////////////////////////////////////////////////////////////////////////////

#endif
