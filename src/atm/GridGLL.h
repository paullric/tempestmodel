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

#include "DataVector.h"
#include "DataMatrix.h"

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
		Model & model,
		int nBaseResolutionA,
		int nBaseResolutionB,
		int nRefinementRatio,
		int nHorizontalOrder,
		int nVerticalOrder,
		int nRElements,
		VerticalStaggering eVerticalStaggering
			= VerticalStaggering_CharneyPhillips
	);

	///	<summary>
	///		Initializer.
	///	</summary>
	void Initialize();

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

private:
	///	<summary>
	///		Compute the central difference penalty (difference of left and
	///		right state from average of left and right state) at interfaces.
	///	</summary>
	void InterpolateCentralDiffPenalty(
		const double * dDataNode,
		bool fZeroBoundaries
	) const;

	///	<summary>
	///		Compute interpolated values of a variable array at edges and
	///		store in m_dStateFEEdge.
	///	</summary>
	void InterpolateNodeToFEEdges(
		const double * dDataNode,
		bool fZeroBoundaries
	) const;

public:
	///	<summary>
	///		Interpolate one column of data from nodes to interfaces.
	///	</summary>
	void InterpolateNodeToREdge(
		const double * dDataNode,
		const double * dDataRefNode,
		double * dDataREdge,
		const double * dDataRefREdge,
		bool fZeroBoundaries = false
	) const;

	///	<summary>
	///		Interpolate one column of data from interfaces to nodes.
	///	</summary>
	void InterpolateREdgeToNode(
		const double * dDataREdge,
		const double * dDataRefREdge,
		double * dDataNode,
		const double * dDataRefNode
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
	///		Differentiate a variable from interfaces to interfaces.
	///	</summary>
	void DifferentiateNodeToREdge(
		const double * dDataNode,
		double * dDiffREdge,
		bool fZeroBoundaries = false
	) const;

	///	<summary>
	///		Differentiate a variable from interfaces to interfaces.
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
	///		Apply discontinuous penalization to nodes.
	///	</summary>
	void CalculateDiscontinuousPenalty(
		const double * dWaveSpeedREdge,
		const double * dDataNode,
		double * dDataUpdate,
		bool fZeroBoundaries = false
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
	const DataMatrix<double> & GetDxBasis1D() const {
		return m_dDxBasis1D;
	}

	///	<summary>
	///		Get the derivatives of the flux reconstruction function, which
	///		is used for DiscontinuousGalerkin dynamics.
	///	</summary>
	const DataVector<double> & GetFluxDeriv1D() const {
		return m_dFluxDeriv1D;
	}

	///	<summary>
	///		Get the stiffness matrix coefficients at nodal points on the
	///		1D reference element.
	///	</summary>
	const DataMatrix<double> & GetStiffness1D() const {
		return m_dStiffness1D;
	}

	///	<summary>
	///		Get the GLL weights at nodal points on the 1D reference element.
	///	</summary>
	const DataVector<double> & GetGLLWeights1D() const {
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

public:
	///	<summary>
	///		Get the interpolation coefficients from levels to interfaces.
	///	</summary>
	const DataMatrix<double> & GetInterpNodeToREdge() const {
		return m_dInterpNodeToREdge;
	}

	///	<summary>
	///		Get the interpolation coefficients from interfaces to levels.
	///	</summary>
	const DataMatrix<double> & GetInterpREdgeToNode() const {
		return m_dInterpREdgeToNode;
	}

	///	<summary>
	///		Get the differentiation coefficients from interfaces to nodes.
	///	</summary>
	const DataMatrix<double> & GetDiffREdgeToNode() const {
		return m_dDiffREdgeToNode;
	}

	///	<summary>
	///		Get the differentiation coefficients from interfaces to interfaces.
	///	</summary>
	const DataMatrix<double> & GetDiffREdgeToREdge() const {
		return m_dDiffREdgeToREdge;
	}

	///	<summary>
	///		Get the differentiation coefficients from nodes to interfaces.
	///	</summary>
	const DataMatrix<double> & GetDiffNodeToREdge() const {
		return m_dDiffNodeToREdge;
	}

	///	<summary>
	///		Get the differentiation coefficients from nodes to nodes.
	///	</summary>
	const DataMatrix<double> & GetDiffNodeToNode() const {
		return m_dDiffNodeToNode;
	}

	///	<summary>
	///		Get amalgamated differentiation coefficients from nodes to
	///		interfaces.
	///	</summary>
	const DataMatrix<double> & GetDiffNodeToREdgeAmal() const {
		return m_dDiffNodeToREdgeAmal;
	}

	///	<summary>
	///		Get amalgamated differentiation coefficients from nodes to
	///		interfaces at left edge.
	///	</summary>
	const DataMatrix<double> & GetDiffNodeToREdgeLeft() const {
		return m_dDiffNodeToREdgeLeft;
	}

	///	<summary>
	///		Get amalgamated differentiation coefficients from nodes to
	///		interfaces at right edge.
	///	</summary>
	const DataMatrix<double> & GetDiffNodeToREdgeRight() const {
		return m_dDiffNodeToREdgeRight;
	}

	///	<summary>
	///		Get the type of reconstruction polynomial.
	///	</summary>
	int GetReconstructionPolyType() const {
		return m_nReconstructionPolyType;
	}

	///	<summary>
	///		Get the derivatives of reconstruction polynomial on nodes.
	///	</summary>
	const DataVector<double> & GetDiffReconsPolyNode() const {
		return m_dDiffReconsPolyNode;
	}

	///	<summary>
	///		Get the derivatives of reconstruction polynomial on interfaces.
	///	</summary>
	const DataVector<double> & GetDiffReconsPolyREdge() const {
		return m_dDiffReconsPolyREdge;
	}

private:
	///	<summary>
	///		Variable evaluated on finite element interfaces (left,
	///		right and average components)
	///	</summary>
	mutable DataMatrix<double> m_dStateFEEdge;

protected:
	///	<summary>
	///		Order of accuracy of the method (number of nodes per element).
	///	</summary>
	int m_nHorizontalOrder;

	///	<summary>
	///		Order of accuracy in the vertical.
	///	</summary>
	int m_nVerticalOrder;

	///	<summary>
	///		Derivatives of the basis functions at nodal points on the
	///		horizontal reference element.
	///	</summary>
	DataMatrix<double> m_dDxBasis1D;

	///	<summary>
	///		Stiffness matrix coefficients at nodal points on the
	///		horizontal reference element.
	///	</summary>
	DataMatrix<double> m_dStiffness1D;

	///	<summary>
	///		Pointwise GLL weights at nodal points on the 1D
	///		reference element.
	///	</summary>
	DataVector<double> m_dGLLWeights1D;

	///	<summary>
	///		Derivatives of the flux reconstruction function (used by
	///		discontinuous Galerkin dynamics).
	///	</summary>
	DataVector<double> m_dFluxDeriv1D;

	///	<summary>
	///		Interpolation coefficients from levels to interfaces.
	///	</summary>
	DataMatrix<double> m_dInterpNodeToREdge;

	///	<summary>
	///		Interpolation coefficients from interfaces to levels.
	///	</summary>
	DataMatrix<double> m_dInterpREdgeToNode;

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
	///		Differentiation coefficients from interfaces to nodes.
	///	</summary>
	DataMatrix<double> m_dDiffREdgeToNode;

	///	<summary>
	///		Differentiation coefficients from interfaces to interfaces.
	///	</summary>
	DataMatrix<double> m_dDiffREdgeToREdge;

	///	<summary>
	///		Differentiation coefficients from nodes to interfaces.
	///	</summary>
	DataMatrix<double> m_dDiffNodeToREdge;

	///	<summary>
	///		Differentiation coefficients from nodes to nodes.
	///	</summary>
	DataMatrix<double> m_dDiffNodeToNode;

	///	<summary>
	///		Amalgamated differentiation coefficients from nodes to interfaces;
	///		includes both interior derivative terms and derivatives from
	///		reconstruction polynomial.
	///	</summary>
	DataMatrix<double> m_dDiffNodeToREdgeAmal;

	///	<summary>
	///		Amalgamated differentiation coefficients from nodes to interfaces
	///		at left edge with extrapolated boundary conditions.
	///	</summary>
	DataMatrix<double> m_dDiffNodeToREdgeLeft;

	///	<summary>
	///		Amalgamated differentiation coefficients from nodes to interfaces
	///		at right edge with extrapolated boundary conditions.
	///	</summary>
	DataMatrix<double> m_dDiffNodeToREdgeRight;

	///	<summary>
	///		Type of reconstruction polynomial to use.
	///	</summary>
	int m_nReconstructionPolyType;

	///	<summary>
	///		Derivatives of reconstruction polynomial on nodes.
	///	</summary>
	DataVector<double> m_dDiffReconsPolyNode;

	///	<summary>
	///		Derivatives of reconstruction polynomial on interfaces.
	///	</summary>
	DataVector<double> m_dDiffReconsPolyREdge;

};

///////////////////////////////////////////////////////////////////////////////

#endif

