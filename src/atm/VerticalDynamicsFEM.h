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
#include "DataMatrix.h"
#include "GridData3D.h"
#include "GridData4D.h"

#ifdef USE_JFNK_PETSC
#include <petscsnes.h>
#endif

//#ifdef USE_JACOBIAN_DEBUG
#include "Model.h"
//#endif

class GridPatch;
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
		int nVerticalOrder,
		int nHypervisOrder = 0,
		bool fFullyExplicit = false,
		bool fUseReferenceState = true,
		bool fExnerPressureOnLevels = true,
		bool fMassFluxOnLevels = false
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
	///		Component indices into the F vector.
	///	</summary>
	typedef int FComp;
	static const FComp FPIx = 0;
	static const FComp FWIx = 1;
	static const FComp FRIx = 2;

	static const FComp FTot = 3;

	///	<summary>
	///		Get the index of the component and level of the F vector.
	///	</summary>
	inline int VecFIx(FComp c, int k) {
#if defined(USE_JACOBIAN_DEBUG)
		int nREdges = m_model.GetGrid()->GetRElements() + 1;
		return (nREdges * c + k);
#else
		return (FTot*k + c);
#endif
	}

	///	<summary>
	///		Get the index of for the component / level pair of the F Jacobian.
	///		Note:  Is composed using Fortran ordering.
	///	</summary>
	inline int MatFIx(FComp c0, int k0, FComp c1, int k1) {
#if defined(USE_JACOBIAN_DEBUG)
		int nREdges = m_model.GetGrid()->GetRElements() + 1;
		return (m_nColumnStateSize * (nREdges * c0 + k0) + (nREdges * c1 + k1));
#elif defined(USE_JACOBIAN_GENERAL)
		return (m_nColumnStateSize * (FTot*k0 + c0) + (FTot*k1 + c1));
#elif defined(USE_JACOBIAN_DIAGONAL)
		return (2 * m_nJacobianFKL + (FTot*k1 + c1) - (FTot*k0 + c0))
			+ m_nColumnStateSize * (FTot*k0 + c0);
#else
		_EXCEPTION();
#endif
	}

protected:
	///	<summary>
	///		Second derivative of a variable from interfaces to interfaces.
	///	</summary>
	void DiffDiffREdgeToREdge(
		const double * dDataREdge,
		double * dDiffDiffREdge
	);

public:
	///	<summary>
	///		Setup the reference column.
	///	</summary>
	void SetupReferenceColumn(
		GridPatch * pPatch,
		int iA,
		int iB,
		const GridData4D & dataRefNode,
		const GridData4D & dataInitialNode,
		const GridData4D & dataRefREdge,
		const GridData4D & dataInitialREdge
/*
		const GridData3D & dataExnerNode,
		const GridData3D & dataDiffExnerNode,
		const GridData3D & dataExnerREdge,
		const GridData3D & dataDiffExnerREdge
*/
	);
 
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
	void BootstrapJacobian();

	///	<summary>
	///		Advance implicit terms of the vertical column one substep.
	///	</summary>
	virtual void StepImplicit(
		int iDataInitial,
		int iDataUpdate,
		const Time & time,
		double dDeltaT
	);

public:
	///	<summary>
	///		Prepare interpolated and differentiated column data.
	///	</summary>
	void PrepareColumn(
		const double * dX
	);

	///	<summary>
	///		Evaluate the zero equations.
	///	</summary>
	void BuildF(
		const double * dX,
		double * dF
	);

	///	<summary>
	///		Build the Jacobian matrix.
	///	</summary>
	void BuildJacobianF(
		const double * dX,
		double * dDG
	);

	///	<summary>
	///		Prepare the column then evaluate the zero equations
	///		(used by JacobianFreeNewtonKrylov)
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
	///		Execute vertical solve as fully explicit.
	///	</summary>
	bool m_fFullyExplicit;

	///	<summary>
	///		Use the background reference profile when computing vertical
	///		derivatives.
	///	</summary>
	bool m_fUseReferenceState;

	///	<summary>
	///		Flag indicating that Exner pressure should be stored on model
	///		levels.
	///	</summary>
	bool m_fExnerPressureOnLevels;

	///	<summary>
	///		Flag indicating that mass flux should be stored on model
	///		levels.
	///	</summary>
	bool m_fMassFluxOnLevels;

	///	<summary>
	///		Hypervis coefficient.
	///	</summary>
	double m_dHypervisCoeff;

	///	<summary>
	///		Order of hyperdiffusion to apply (must be even).
	///	</summary>
	int m_nHypervisOrder;

	///	<summary>
	///		Timestep size.
	///	</summary>
	double m_dDeltaT;

	///	<summary>
	///		Pointer to active patch.
	///	</summary>
	GridPatch * m_pPatch;

	///	<summary>
	///		Active alpha index on m_pPatch.
	///	</summary>
	int m_iA;

	///	<summary>
	///		Active beta index on m_pPatch.
	///	</summary>
	int m_iB;

	///	<summary>
	///		Number of degrees of freedom in vertical solution vector.
	///	</summary>
	int m_nColumnStateSize;

protected:
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
	///		Auxiliary state, used by StepExplicit.
	///	</summary>
	DataVector<double> m_dStateAux;

	///	<summary>
	///		Derivative of auxiliary state, used by StepExplicit.
	///	</summary>
	DataVector<double> m_dStateAuxDiff;

	///	<summary>
	///		Velocity across xi surfaces (xi_dot) at nodes.
	///	</summary>
	DataVector<double> m_dXiDotNode;

	///	<summary>
	///		Velocity across xi surfaces (xi_dot) at interfaces.
	///	</summary>
	DataVector<double> m_dXiDotREdge;

	///	<summary>
	///		Auxiliary storage for derivative of theta.
	///	</summary>
	DataVector<double> m_dDiffP;

	///	<summary>
	///		Auxiliary storage for higher derivatives of theta.
	///	</summary>
	DataVector<double> m_dDiffDiffP;

	///	<summary>
	///		Horizontal Kinetic energy on model levels.
	///	</summary>
	DataVector<double> m_dHorizKineticEnergyNode;

	///	<summary>
	///		Kinetic energy on model levels.
	///	</summary>
	DataVector<double> m_dKineticEnergyNode;

	///	<summary>
	///		Derivatives of kinetic energy on model levels.
	///	</summary>
	DataVector<double> m_dDiffKineticEnergyNode;

	///	<summary>
	///		Derivatives of kinetic energy on model levels.
	///	</summary>
	DataVector<double> m_dDiffKineticEnergyREdge;

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
	///		Pressure flux on model levels.
	///	</summary>
	DataVector<double> m_dPressureFluxNode;

	///	<summary>
	///		Pressure flux on model interfaces.
	///	</summary>
	DataVector<double> m_dPressureFluxREdge;

	///	<summary>
	///		Derivatives of pressure flux on model levels.
	///	</summary>
	DataVector<double> m_dDiffPressureFluxNode;

	///	<summary>
	///		Derivatives of pressure flux on model interfaces.
	///	</summary>
	DataVector<double> m_dDiffPressureFluxREdge;

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
	///		Second differentiation coefficients from nodes to nodes.
	///	</summary>
	DataMatrix<double> m_dDiffDiffNodeToNode;

	///	<summary>
	///		Hyperviscosity coefficients from nodes to nodes.
	///	</summary>
	DataMatrix<double> m_dHypervisNodeToNode;

	///	<summary>
	///		Second differentiation coefficients from edges to edges.
	///	</summary>
	DataMatrix<double> m_dDiffDiffREdgeToREdge;

	///	<summary>
	///		Hyperviscosity coefficients from edges to edges.
	///	</summary>
	DataMatrix<double> m_dHypervisREdgeToREdge;

#ifdef USE_JFNK_PETSC
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
#ifdef USE_DIRECTSOLVE_APPROXJ
private:
	///	<summary>
	///		Jacobian matrix used in direct solve.
	///	</summary>
	DataMatrix<double> m_matJacobianF;

	///	<summary>
	///		Pivot matrix used in direct solve.
	///	</summary>
	DataVector<int> m_vecIPiv;
#endif
#ifdef USE_DIRECTSOLVE
private:
	///	<summary>
	///		Jacobian matrix used in direct solve.
	///	</summary>
	DataMatrix<double> m_matJacobianF;

	///	<summary>
	///		Pivot matrix used in direct solve.
	///	</summary>
	DataVector<int> m_vecIPiv;
#endif
#ifdef USE_JACOBIAN_DIAGONAL
private:
	///	<summary>
	///		Number of sub-diagonals in Jacobian.
	///	</summary>
	int m_nJacobianFKL;

	///	<summary>
	///		Number of super-diagonals in Jacobian.
	///	</summary>
	int m_nJacobianFKU;
#endif
};

///////////////////////////////////////////////////////////////////////////////

#ifdef USE_JFNK_PETSC
PetscErrorCode VerticalDynamicsFEM_FormFunction(
	SNES snes,
	Vec x,
	Vec f,
	void * pDyn
);
#endif

///////////////////////////////////////////////////////////////////////////////

#endif

