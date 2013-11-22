///////////////////////////////////////////////////////////////////////////////
///
///	\file    VerticalDynamicsFEM.h
///	\author  Paul Ullrich
///	\version May 20, 2013
///
///	<remarks>
///		Copyright 2000-2013 Paul Ullrich
///
///		This file is distributed as part of the Tempest source code package.
///		Permission is granted to use, copy, modify and distribute this
///		source code and its documentation under the terms of the GNU General
///		Public License.  This software is provided "as is" without express
///		or implied warranty.
///	</remarks>

#ifndef _VERTICALDYNAMICSFEM_H_
#define _VERTICALDYNAMICSFEM_H_

///////////////////////////////////////////////////////////////////////////////

#include "Defines.h"
#include "VerticalDynamics.h"
#include "JacobianFreeNewtonKrylov.h"
#include "DataVector.h"

#ifdef USE_PETSC
#include <petscsnes.h>
#endif

class Time;

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Finite-element based atmospheric vertical dynamics.
///	</summary>
class VerticalDynamicsFEM :
	public VerticalDynamics,
	public JacobianFreeNewtonKrylov
{

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	VerticalDynamicsFEM(
		Model & model,
		int nHorizontalOrder,
		int nVerticalOrder
	);

	///	<summary>
	///		Destructor.
	///	</summary>
	~VerticalDynamicsFEM();

public:
	///	<summary>
	///		Initializer.
	///	</summary>
	virtual void Initialize();

protected:
	///	<summary>
	///		Interpolate one column of data from nodes to interfaces.
	///	</summary>
	void InterpolateNodeToREdge(
		const double * dDataNode,
		const double * dDataRefNode,
		double * dDataREdge,
		const double * dDataRefREdge
	);

	///	<summary>
	///		Interpolate one column of data from interfaces to nodes.
	///	</summary>
	void InterpolateREdgeToNode(
		const double * dDataREdge,
		const double * dDataRefREdge,
		double * dDataNode,
		const double * dDataRefNode
	);

	///	<summary>
	///		Interpolate a state variable from nodes to finite-element edges.
	///	</summary>
	void InterpolateNodeToFEEdges(
		const double * dDataNode,
		bool fZeroBoundaries
	);

	///	<summary>
	///		Differentiate a state variable from nodes to nodes.
	///	</summary>
	void DifferentiateNodeToNode(
		const double * dDataNode,
		double * dDiffNode,
		bool fZeroBoundaries = false
	);

	///	<summary>
	///		Differentiate a state varaible from interfaces to interfaces.
	///	</summary>
	void DifferentiateNodeToREdge(
		const double * dDataNode,
		double * dDiffREdge,
		bool fZeroBoundaries = false
	);

	///	<summary>
	///		Differentiate a state varaible from interfaces to interfaces.
	///	</summary>
	void DifferentiateREdgeToNode(
		const double * dDataREdge,
		double * dDiffNode
	);

	///	<summary>
	///		Differentiate a state varaible from interfaces to interfaces.
	///	</summary>
	void DifferentiateREdgeToREdge(
		const double * dDataREdge,
		double * dDiffREdge
	);

public:
	///	<summary>
	///		Advance explicit terms of the vertical column one substep.
	///	</summary>
	virtual void StepExplicit(
		int iDataInitial,
		int iDataUpdate,
		const Time & time,
		double dDeltaT
	);

	///	<summary>
	///		Build the Jacobian matrix.
	///	</summary>
	void BuildJacobian();

	///	<summary>
	///		Advance implicit terms of the vertical column one substep.
	///	</summary>
	virtual void StepImplicit(
		int iDataInitial,
		int iDataUpdate,
		const Time & time,
		double dDeltaT
	);

protected:
	///	<summary>
	///		Evaluate the RHS when all variables are on model levels.
	///	</summary>
	void EvaluateDG(
		const double * dX,
		double * dF
	);

	///	<summary>
	///		Evaluate the RHS when all variables are on model interfaces.
	///	</summary>
	void EvaluateSE(
		const double * dX,
		double * dF
	);

	///	<summary>
	///		Evaluate the RHS when rho is on model levels.
	///	</summary>
	void EvaluateMixed(
		const double * dX,
		double * dF
	);

public:
	///	<summary>
	///		Evaluate the RHS (used by JacobianFreeNewtonKrylov)
	///	</summary>
	void Evaluate(
		const double * dX,
		double * dF
	);

protected:
	///	<summary>
	///		Horizontal order of accuracy of the method.
	///	</summary>
	int m_nHorizontalOrder;

	///	<summary>
	///		Vertical order of accuacy of the method.
	///	</summary>
	int m_nVerticalOrder;

	///	<summary>
	///		Timestep size.
	///	</summary>
	double m_dDeltaT;

	///	<summary>
	///		Domain height.
	///	</summary>
	double m_dDomainHeight;

	///	<summary>
	///		Number of degrees of freedom in vertical solution vector.
	///	</summary>
	int m_nColumnStateSize;

	///	<summary>
	///		First index of theta in the column state vector.
	///	</summary>
	int m_ixTBegin;

	///	<summary>
	///		First index of W in the column state vector.
	///	</summary>
	int m_ixWBegin;

	///	<summary>
	///		First index of rho in the column state vector.
	///	</summary>
	int m_ixRBegin;

	///	<summary>
	///		State variable column.
	///	</summary>
	DataVector<double> m_dColumnState;

	///	<summary>
	///		Reference state column on nodes.
	///	</summary>
	DataMatrix<double> m_dStateRefNode;

	///	<summary>
	///		Reference state column on interfaces.
	///	</summary>
	DataMatrix<double> m_dStateRefREdge;

	///	<summary>
	///		State vector on model levels, used by StepImplicit.
	///	</summary>
	DataMatrix<double> m_dStateNode;

	///	<summary>
	///		State vector on model interfaces, used by StepImplicit.
	///	</summary>
	DataMatrix<double> m_dStateREdge;

	///	<summary>
	///		Auxiliary state, used by StepImplicit
	///	</summary>
	DataVector<double> m_dStateAux;

	///	<summary>
	///		Derivative of auxiliary state, used by StepImplicit
	///	</summary>
	DataVector<double> m_dStateAuxDiff;

	///	<summary>
	///		State variable evaluated on finite element interfaces (left,
	///		right and average components)
	///	</summary>
	DataMatrix<double> m_dStateFEEdge;

	///	<summary>
	///		Mass flux on model levels.
	///	</summary>
	DataVector<double> m_dMassFluxNode;

	///	<summary>
	///		Mass flux on model interfaces.
	///	</summary>
	DataVector<double> m_dMassFluxREdge;

	///	<summary>
	///		Derivatives of mass flux on model levels.
	///	</summary>
	DataVector<double> m_dDiffMassFluxNode;

	///	<summary>
	///		Derivatives of mass flux on model interfaces.
	///	</summary>
	DataVector<double> m_dDiffMassFluxREdge;

	///	<summary>
	///		Exner pressure perturbation at model levels.
	///	</summary>
	DataVector<double> m_dExnerPertNode;

	///	<summary>
	///		Reference Exner pressure at model levels.
	///	</summary>
	DataVector<double> m_dExnerRefNode;

	///	<summary>
	///		Exner pressure perturbation at model interfaces.
	///	</summary>
	DataVector<double> m_dExnerPertREdge;

	///	<summary>
	///		Reference Exner pressure at model interfaces.
	///	</summary>
	DataVector<double> m_dExnerRefREdge;

	///	<summary>
	///		Derivative of Exner pressure perturbation at model levels.
	///	</summary>
	DataVector<double> m_dDiffExnerPertNode;

	///	<summary>
	///		Derivative of reference Exner pressure at model levels.
	///	</summary>
	DataVector<double> m_dDiffExnerRefNode;

	///	<summary>
	///		Derivative of Exner pressure perturbation at model interfaces.
	///	</summary>
	DataVector<double> m_dDiffExnerPertREdge;

	///	<summary>
	///		Derivative of reference Exner pressure at model interfaces.
	///	</summary>
	DataVector<double> m_dDiffExnerRefREdge;

	///	<summary>
	///		Solution vector from the implicit solve.
	///	</summary>
	DataVector<double> m_dSoln;

	///	<summary>
	///		Interpolation coefficients from nodes to interfaces.
	///	</summary>
	const DataMatrix<double> * m_pInterpNodeToREdge;

	///	<summary>
	///		Interpolation coefficients from interfaces to nodes.
	///	</summary>
	const DataMatrix<double> * m_pInterpREdgeToNode;

	///	<summary>
	///		Differentiation coefficients from extended levels (levels plus
	///		the element endpoints) to levels.
	///	</summary>
	DataMatrix<double> m_dDiffExtToLevel;

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
	///		Derivatives of reconstruction polynomial on nodes.
	///	</summary>
	DataVector<double> m_dDiffReconsPolyNode;

	///	<summary>
	///		Derivatives of reconstruction polynomial on interfaces.
	///	</summary>
	DataVector<double> m_dDiffReconsPolyREdge;

#ifdef USE_PETSC
private:
	///	<summary>
	///		PetSc solver context
	///	</summary>
	SNES m_snes;

	///	<summary>
	///		Jacobian matrix.
	///	</summary>
	Mat m_matJ;

	///	<summary>
	///		PetSc vector for storing the evaluated solution.
	///	</summary>
	Vec m_vecX;

	///	<summary>
	///		PetSc vector for storing the function value.
	///	</summary>
	Vec m_vecR;

#endif
};

///////////////////////////////////////////////////////////////////////////////

#ifdef USE_PETSC
PetscErrorCode VerticalDynamicsFEM_FormFunction(
	SNES snes,
	Vec x,
	Vec f,
	void * pDyn
);
#endif

///////////////////////////////////////////////////////////////////////////////

#endif

