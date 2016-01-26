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

#ifndef _GRIDGLL_H_
#define _GRIDGLL_H_

#include "Grid.h"
#include "LinearColumnOperatorFEM.h"

#include "DataArray1D.h"
#include "DataArray2D.h"

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		The CubedSphere grid with degrees of freedom at Gauss-Lobatto-Legendre
///		quadrature nodes.
///	</summary>
class GridGLL : public Grid {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	GridGLL(
		Model & model
	);

	///	<summary>
	///		Define the parameters for the GridGLL.
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
		VerticalStaggering eVerticalStaggering
			= VerticalStaggering_CharneyPhillips
	);

	///	<summary>
	///		Initializer.
	///	</summary>
	virtual void Initialize();

	///	<summary>
	///		Initialize the vertical coordinate.
	///	</summary>
	virtual void InitializeVerticalCoordinate();

	///	<summary>
	///		Initialize topography height/derivatives, state and tracer data
	///		from a TestCase.
	///	</summary>
	virtual void EvaluateTestCase(
		const TestCase & test,
		const Time & time,
		int iDataIndex = 0
	);

public:
	///	<summary>
	///		Perform post-processing of variables on the grid after each
	///		TimeStep substage.
	///	</summary>
	virtual void PostProcessSubstage(
		int iDataUpdate,
		DataType eDataType = DataType_State
	);

	///	<summary>
	///		Apply the direct stiffness summation (DSS) operation on the grid.
	///	</summary>
	virtual void ApplyDSS(
		int iDataUpdate,
		DataType eDataType = DataType_State
	) {
		_EXCEPTIONT("Unimplemented");
	}

	///	<summary>
	///		Compute vorticity on the grid.
	///	</summary>
	virtual void ComputeVorticityDivergence(
		int iDataIndex
	);

public:
	///	<summary>
	///		Interpolate one column of data from nodes to the given interface.
	///	</summary>
	double InterpolateNodeToREdge(
		const double * dDataNode,
		const double * dDataRefNode,
		int iRint,
		double dDataRefREdge,
		int nStride = 1
	) const;

	///	<summary>
	///		Interpolate one column of data from interfaces to the given node.
	///	</summary>
	double InterpolateREdgeToNode(
		const double * dDataREdge,
		const double * dDataRefREdge,
		int iRnode,
		double dDataRefNode,
		int nStride = 1
	) const;

	///	<summary>
	///		Differentiate a variable from nodes to the given node.
	///	</summary>
	double DifferentiateNodeToNode(
		const double * dDataNode,
		int iRnode,
		int nStride = 1
	) const;

	///	<summary>
	///		Differentiate a variable from nodes to the given interface.
	///	</summary>
	double DifferentiateNodeToREdge(
		const double * dDataNode,
		int iRnode,
		int nStride = 1
	) const;

	///	<summary>
	///		Differentiate a variable from interfaces to the given node.
	///	</summary>
	double DifferentiateREdgeToNode(
		const double * dDataREdge,
		int iRnode,
		int nStride = 1
	) const;

	///	<summary>
	///		Differentiate a variable from interfaces to the given interface.
	///	</summary>
	double DifferentiateREdgeToREdge(
		const double * dDataREdge,
		int iRint,
		int nStride = 1
	) const;

public:
	///	<summary>
	///		Interpolate one column of data from nodes to interfaces.
	///	</summary>
	void InterpolateNodeToREdge(
		const double * dDataNode,
		double * dDataREdge,
		bool fZeroBoundaries = false
	) const;

	///	<summary>
	///		Interpolate one column of data from interfaces to nodes.
	///	</summary>
	void InterpolateREdgeToNode(
		const double * dDataREdge,
		double * dDataNode
	) const;

	///	<summary>
	///		Differentiate a variable from nodes to nodes.
	///	</summary>
	void DifferentiateNodeToNode(
		const double * dDataNode,
		double * dDiffNode,
		bool fZeroBoundaries = false
	) const;

	///	<summary>
	///		Differentiate a variable from nodes to interfaces.
	///	</summary>
	void DifferentiateNodeToREdge(
		const double * dDataNode,
		double * dDiffREdge,
		bool fZeroBoundaries = false
	) const;

	///	<summary>
	///		Differentiate a variable from interfaces to nodes.
	///	</summary>
	void DifferentiateREdgeToNode(
		const double * dDataREdge,
		double * dDiffNode
	) const;

	///	<summary>
	///		Differentiate a variable from interfaces to interfaces.
	///	</summary>
	void DifferentiateREdgeToREdge(
		const double * dDataREdge,
		double * dDiffREdge
	) const;

	///	<summary>
	///		Differentiate twice a variable from nodes to nodes.
	///	</summary>
	void DiffDiffNodeToNode(
		const double * dDataNode,
		double * dDiffDiffNode
	) const;

	///	<summary>
	///		Differentiate twice a variable from interfaces to interfaces.
	///	</summary>
	void DiffDiffREdgeToREdge(
		const double * dDataREdge,
		double * dDiffDiffREdge
	) const;

public:
	///	<summary>
	///		Get the order of accuracy of the method.
	///	</summary>
	int GetHorizontalOrder() const {
		return m_nHorizontalOrder;
	}

	///	<summary>
	///		Get the order of accuracy in the vertical.
	///	</summary>
	int GetVerticalOrder() const {
		return m_nVerticalOrder;
	}

	///	<summary>
	///		Get the derivatives of the basis functions at nodal points on the
	///		1D reference element.
	///	</summary>
	const DataArray2D<double> & GetDxBasis1D() const {
		return m_dDxBasis1D;
	}

	///	<summary>
	///		Get the derivatives of the flux reconstruction function, which
	///		is used for DiscontinuousGalerkin dynamics.
	///	</summary>
	const DataArray1D<double> & GetFluxDeriv1D() const {
		return m_dFluxDeriv1D;
	}

	///	<summary>
	///		Get the stiffness matrix coefficients at nodal points on the
	///		1D reference element.
	///	</summary>
	const DataArray2D<double> & GetStiffness1D() const {
		return m_dStiffness1D;
	}

	///	<summary>
	///		Get the GLL weights at nodal points on the 1D reference element.
	///	</summary>
	const DataArray1D<double> & GetGLLWeights1D() const {
		return m_dGLLWeights1D;
	}

public:
	///	<summary>
	///		Get the interpolation operator from levels to interfaces.
	///	</summary>
	const LinearColumnInterpFEM & GetOpInterpNodeToREdge() const {
		return m_opInterpNodeToREdge;
	}

	///	<summary>
	///		Get the interpolation operator from interfaces to levels.
	///	</summary>
	const LinearColumnInterpFEM & GetOpInterpREdgeToNode() const {
		return m_opInterpREdgeToNode;
	}

	///	<summary>
	///		Get the differentiation operator from levels to levels.
	///	</summary>
	const LinearColumnDiffFEM & GetOpDiffNodeToNode() const {
		return m_opDiffNodeToNode;
	}

	///	<summary>
	///		Get the differentiation operator from levels to interfaces.
	///	</summary>
	const LinearColumnDiffFEM & GetOpDiffNodeToREdge() const {
		return m_opDiffNodeToREdge;
	}

	///	<summary>
	///		Get the differentiation operator from interfacs to levels.
	///	</summary>
	const LinearColumnDiffFEM & GetOpDiffREdgeToNode() const {
		return m_opDiffREdgeToNode;
	}

	///	<summary>
	///		Get the differentiation operator from interfacs to interfaces.
	///	</summary>
	const LinearColumnDiffFEM & GetOpDiffREdgeToREdge() const {
		return m_opDiffREdgeToREdge;
	}

	///	<summary>
	///		Get the second derivative operator from levels to levels.
	///	</summary>
	const LinearColumnDiffDiffFEM & GetOpDiffDiffNodeToNode() const {
		return m_opDiffDiffNodeToNode;
	}

	///	<summary>
	///		Get the second derivative operator from interfaces to interfaces.
	///	</summary>
	const LinearColumnDiffDiffFEM & GetOpDiffDiffREdgeToREdge() const {
		return m_opDiffDiffREdgeToREdge;
	}

	///	<summary>
	///		Get the second derivative operator from levels to levels.
	///	</summary>
	const LinearColumnDiscPenaltyFEM & GetOpPenaltyNodeToNode() const {
		return m_opPenaltyNodeToNode;
	}

protected:
	///	<summary>
	///		Order of accuracy of the method (number of nodes per element).
	///	</summary>
	DataStruct<int> m_nHorizontalOrder;

	///	<summary>
	///		Order of accuracy in the vertical.
	///	</summary>
	DataStruct<int> m_nVerticalOrder;

protected:
	///	<summary>
	///		Derivatives of the basis functions at nodal points on the
	///		horizontal reference element.
	///	</summary>
	DataArray2D<double> m_dDxBasis1D;

	///	<summary>
	///		Stiffness matrix coefficients at nodal points on the
	///		horizontal reference element.
	///	</summary>
	DataArray2D<double> m_dStiffness1D;

	///	<summary>
	///		Pointwise GLL weights at nodal points on the 1D
	///		reference element.
	///	</summary>
	DataArray1D<double> m_dGLLWeights1D;

	///	<summary>
	///		Derivatives of the flux reconstruction function (used by
	///		discontinuous Galerkin dynamics).
	///	</summary>
	DataArray1D<double> m_dFluxDeriv1D;

	///	<summary>
	///		Interpolation operator from levels to interfaces.
	///	</summary>
	LinearColumnInterpFEM m_opInterpNodeToREdge;

	///	<summary>
	///		Interpolation operator from interfaces to levels.
	///	</summary>
	LinearColumnInterpFEM m_opInterpREdgeToNode;

	///	<summary>
	///		Differentiation operator from levels to levels.
	///	</summary>
	LinearColumnDiffFEM m_opDiffNodeToNode;

	///	<summary>
	///		Differentiation operator from levels to levels with zero boundaries.
	///	</summary>
	LinearColumnDiffFEM m_opDiffNodeToNodeZeroBoundaries;

	///	<summary>
	///		Differentiation operator from levels to interfaces.
	///	</summary>
	LinearColumnDiffFEM m_opDiffNodeToREdge;

	///	<summary>
	///		Differentiation operator from interfaces to levels.
	///	</summary>
	LinearColumnDiffFEM m_opDiffREdgeToNode;

	///	<summary>
	///		Differentiation operator from interfaces to interfaces.
	///	</summary>
	LinearColumnDiffFEM m_opDiffREdgeToREdge;

	///	<summary>
	///		Second derivative operator from levels to levels.
	///	</summary>
	LinearColumnDiffDiffFEM m_opDiffDiffNodeToNode;

	///	<summary>
	///		Second derivative operator from interfaces to interfaces.
	///	</summary>
	LinearColumnDiffDiffFEM m_opDiffDiffREdgeToREdge;

	///	<summary>
	///		Discontinuous penalty operator from levels to levels.
	///	</summary>
	LinearColumnDiscPenaltyFEM m_opPenaltyNodeToNode;

	///	<summary>
	///		Integral operator from levels to levels.
	///	</summary>
	LinearColumnIntFEM m_opIntNodeToNode;
};

///////////////////////////////////////////////////////////////////////////////

#endif

