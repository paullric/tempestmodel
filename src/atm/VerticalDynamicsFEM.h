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
#include "DataArray1D.h"
#include "DataArray2D.h"
#include "DataArray3D.h"
#include "DataArray4D.h"

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
		bool fForceMassFluxOnLevels = false
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
	///		Convert from standard state indices to component indices.
	///	</summary>
	FComp FIxFromCIx(int iC) {
		return (iC - 2);
	}

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

	///	<summary>
	///		Get the index of for the component / level pair of the F Jacobian.
	///		Note:  Is composed using Fortran ordering.
	///	</summary>
	inline int TracerMatFIx(int k0, int k1) {
#if defined(USE_JACOBIAN_DEBUG) || defined(USE_JACOBIAN_GENERAL)
		return (m_nRElements * k0 + k1);
#elif defined(USE_JACOBIAN_DIAGONAL)
		int nOffDiagonals = 2 * (2 * m_nVerticalOrder - 1);
		return (nOffDiagonals + k1 - k0) + m_nRElements * k0;
#else
		_EXCEPTION();
#endif
	}

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
	///		Advance an explicit update of implicit terms one time
	///	</summary>
	virtual void StepImplicitTermsExplicitly(
		int iDataInitial,
		int iDataUpdate,
		const Time & time,
		double dDeltaT
	);

	///	<summary>
	///		Force a full explicit update one time only
	///	</summary>
	virtual void ForceStepExplicit(
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
	///		Set up the reference column.  This function is called once for
	///		each column prior to the solve.
	///	</summary>
	void SetupReferenceColumn(
		GridPatch * pPatch,
		int iA,
		int iB,
		const DataArray4D<double> & dataRefNode,
		const DataArray4D<double> & dataInitialNode,
		const DataArray4D<double> & dataRefREdge,
		const DataArray4D<double> & dataInitialREdge
	);

	///	<summary>
	///		Prepare interpolated and differentiated column data.  This
	///		function contains shared initialization between BuildF and
	///		BuildJacobianF.
	///	</summary>
	void PrepareColumn(
		const double * dX
	);

	///	<summary>
	///		Evaluate the zero equations for the implicit solve.
	///	</summary>
	void BuildF(
		const double * dX,
		double * dF
	);

	///	<summary>
	///		Build the Jacobian matrix associated with the zero equations.
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
	///		Update tracers in the vertical.
	///	</summary>
	void UpdateColumnTracers(
		double dDeltaT,
		const DataArray4D<double> & dataInitialNode,
		const DataArray4D<double> & dataUpdateNode,
		const DataArray4D<double> & dataInitialREdge,
		const DataArray4D<double> & dataUpdateREdge,
		const DataArray4D<double> & dataInitialTracer,
		DataArray4D<double> & dataUpdateTracer
	);

public:
	///	<summary>
	///		Apply a positive definite filter to tracers in each column.
	///	</summary>
	virtual void FilterNegativeTracers(
		int iDataUpdate
	);

public:
	///	<summary>
	///		Execute vertical solve as fully explicit.
	///	</summary>
	bool IsFullyExplicit() const {
		return m_fFullyExplicit;
	}

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
	///		Flag indicating that mass flux should be calculated on model
	///		levels.
	///	</summary>
	bool m_fForceMassFluxOnLevels;

	///	<summary>
	///		Hypervis coefficient.
	///	</summary>
	double m_dHypervisCoeff;

	///	<summary>
	///		Variables to upwind.
	///	</summary>
	DataArray1D<bool> m_fUpwind;

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
	///		Number of radial elements in the vertical column.
	///	</summary>
	int m_nRElements;

	///	<summary>
	///		Number of degrees of freedom in vertical solution vector.
	///	</summary>
	int m_nColumnStateSize;

protected:
	///	<summary>
	///		Finite element upwinding weights.
	///	</summary>
	DataArray1D<double> m_dUpwindWeights;

	///	<summary>
	///		State variable column.
	///	</summary>
	DataArray1D<double> m_dColumnState;

	///	<summary>
	///		Reference state column on nodes.
	///	</summary>
	DataArray2D<double> m_dStateRefNode;

	///	<summary>
	///		Reference state column on interfaces.
	///	</summary>
	DataArray2D<double> m_dStateRefREdge;

	///	<summary>
	///		State vector on model levels, used by StepImplicit.
	///	</summary>
	DataArray2D<double> m_dStateNode;

	///	<summary>
	///		State vector on model interfaces, used by StepImplicit.
	///	</summary>
	DataArray2D<double> m_dStateREdge;

	///	<summary>
	///		Auxiliary state data.
	///	</summary>
	DataArray1D<double> m_dStateAux;

	///	<summary>
	///		Derivative of auxiliary state data.
	///	</summary>
	DataArray1D<double> m_dStateAuxDiff;

	///	<summary>
	///		Velocity across xi surfaces (xi_dot) at nodes.
	///	</summary>
	DataArray1D<double> m_dXiDotNode;

	///	<summary>
	///		Velocity across xi surfaces (xi_dot) at interfaces.
	///	</summary>
	DataArray1D<double> m_dXiDotREdge;

	///	<summary>
	///		Auxiliary storage for derivative of alpha velocity.
	///	</summary>
	DataArray1D<double> m_dDiffUa;

	///	<summary>
	///		Auxiliary storage for derivative of beta velocity.
	///	</summary>
	DataArray1D<double> m_dDiffUb;

	///	<summary>
	///		Auxiliary storage for derivative of pressure on nodes.
	///	</summary>
	DataArray1D<double> m_dDiffPNode;

	///	<summary>
	///		Auxiliary storage for derivative of pressure on interfaces.
	///	</summary>
	DataArray1D<double> m_dDiffPREdge;

	///	<summary>
	///		Auxiliary storage for second derivatives of the state
	///	</summary>
	DataArray2D<double> m_dDiffDiffState;

	///	<summary>
	///		Auxiliary storage for higher derivatives of the state
	///	</summary>
	DataArray2D<double> m_dHyperDiffState;

	///	<summary>
	///		Auxiliary storage for derivative of theta on nodes.
	///	</summary>
	DataArray1D<double> m_dDiffThetaNode;

	///	<summary>
	///		Auxiliary storage for derivative of theta on interfaces.
	///	</summary>
	DataArray1D<double> m_dDiffThetaREdge;

	///	<summary>
	///		Auxiliary storage for derivative of vertical velocity on nodes.
	///	</summary>
	DataArray1D<double> m_dDiffWNode;

	///	<summary>
	///		Auxiliary storage for derivative of vertical velocity on interfaces.
	///	</summary>
	DataArray1D<double> m_dDiffWREdge;

	///	<summary>
	///		Horizontal Kinetic energy on model levels.
	///	</summary>
	DataArray1D<double> m_dHorizKineticEnergyNode;

	///	<summary>
	///		Kinetic energy on model levels.
	///	</summary>
	DataArray1D<double> m_dKineticEnergyNode;

	///	<summary>
	///		Derivatives of kinetic energy on model levels.
	///	</summary>
	DataArray1D<double> m_dDiffKineticEnergyNode;

	///	<summary>
	///		Derivatives of kinetic energy on model levels.
	///	</summary>
	DataArray1D<double> m_dDiffKineticEnergyREdge;

	///	<summary>
	///		Mass flux on model levels.
	///	</summary>
	DataArray1D<double> m_dMassFluxNode;

	///	<summary>
	///		Mass flux on model interfaces.
	///	</summary>
	DataArray1D<double> m_dMassFluxREdge;

	///	<summary>
	///		Derivatives of mass flux on model levels.
	///	</summary>
	DataArray1D<double> m_dDiffMassFluxNode;

	///	<summary>
	///		Derivatives of mass flux on model interfaces.
	///	</summary>
	DataArray1D<double> m_dDiffMassFluxREdge;

	///	<summary>
	///		Pressure flux on model levels.
	///	</summary>
	DataArray1D<double> m_dPressureFluxNode;

	///	<summary>
	///		Pressure flux on model interfaces.
	///	</summary>
	DataArray1D<double> m_dPressureFluxREdge;

	///	<summary>
	///		Derivatives of pressure flux on model levels.
	///	</summary>
	DataArray1D<double> m_dDiffPressureFluxNode;

	///	<summary>
	///		Derivatives of pressure flux on model interfaces.
	///	</summary>
	DataArray1D<double> m_dDiffPressureFluxREdge;

	///	<summary>
	///		Exner pressure perturbation at model levels.
	///	</summary>
	DataArray1D<double> m_dExnerNode;

	///	<summary>
	///		Reference Exner pressure at model levels.
	///	</summary>
	DataArray1D<double> m_dExnerRefNode;

	///	<summary>
	///		Exner pressure perturbation at model interfaces.
	///	</summary>
	DataArray1D<double> m_dExnerREdge;

	///	<summary>
	///		Reference Exner pressure at model interfaces.
	///	</summary>
	DataArray1D<double> m_dExnerRefREdge;

	///	<summary>
	///		Derivative of Exner pressure perturbation at model levels.
	///	</summary>
	DataArray1D<double> m_dDiffExnerPertNode;

	///	<summary>
	///		Derivative of reference Exner pressure at model levels.
	///	</summary>
	DataArray1D<double> m_dDiffExnerRefNode;

	///	<summary>
	///		Derivative of Exner pressure perturbation at model interfaces.
	///	</summary>
	DataArray1D<double> m_dDiffExnerPertREdge;

	///	<summary>
	///		Derivative of reference Exner pressure at model interfaces.
	///	</summary>
	DataArray1D<double> m_dDiffExnerRefREdge;

	///	<summary>
	///		Tracer density on model levels.
	///	</summary>
	DataArray1D<double> m_dTracerDensityNode;

	///	<summary>
	///		Tracer density on model interfaces.
	///	</summary>
	DataArray1D<double> m_dTracerDensityREdge;

	///	<summary>
	///		Initial density on model levels.
	///	</summary>
	DataArray1D<double> m_dInitialDensityNode;

	///	<summary>
	///		Initial density on model interfaces.
	///	</summary>
	DataArray1D<double> m_dInitialDensityREdge;

	///	<summary>
	///		Updated density on model levels.
	///	</summary>
	DataArray1D<double> m_dUpdateDensityNode;

	///	<summary>
	///		Updated density on model interfaces.
	///	</summary>
	DataArray1D<double> m_dUpdateDensityREdge;

	///	<summary>
	///		Solution vector from the implicit solve.
	///	</summary>
	DataArray1D<double> m_dSoln;

private:
	///	<summary>
	///		Jacobian in the column on Nodes.
	///	</summary>
	DataArray1D<double> m_dColumnJacobianNode;

	///	<summary>
	///		Jacobian in the column on REdges.
	///	</summary>
	DataArray1D<double> m_dColumnJacobianREdge;

	///	<summary>
	///		Element area on Nodes.
	///	</summary>
	DataArray1D<double> m_dColumnElementArea;

	///	<summary>
	///		Inverse Jacobian in the column on Nodes.
	///	</summary>
	DataArray1D<double> m_dColumnInvJacobianNode;

	///	<summary>
	///		Inverse Jacobian in the column on REdges.
	///	</summary>
	DataArray1D<double> m_dColumnInvJacobianREdge;

	///	<summary>
	///		Vertical derivative transform in the column on Nodes.
	///	</summary>
	DataArray2D<double> m_dColumnDerivRNode;

	///	<summary>
	///		Vertical derivative transform in the column on REdges.
	///	</summary>
	DataArray2D<double> m_dColumnDerivRREdge;

	///	<summary>
	///		Contravariant metric (alpha component) in the column on Nodes.
	///	</summary>
	DataArray2D<double> m_dColumnContraMetricA;

	///	<summary>
	///		Contravariant metric (beta component) in the column on Nodes.
	///	</summary>
	DataArray2D<double> m_dColumnContraMetricB;

	///	<summary>
	///		Contravariant metric (xi component) in the column on Nodes.
	///	</summary>
	DataArray2D<double> m_dColumnContraMetricXi;

	///	<summary>
	///		Contravariant metric (alpha component) in the column on REdges.
	///	</summary>
	DataArray2D<double> m_dColumnContraMetricAREdge;

	///	<summary>
	///		Contravariant metric (beta component) in the column on REdges.
	///	</summary>
	DataArray2D<double> m_dColumnContraMetricBREdge;

	///	<summary>
	///		Contravariant metric (xi component) in the column on REdges.
	///	</summary>
	DataArray2D<double> m_dColumnContraMetricXiREdge;

private:
	///	<summary>
	///		Flux vector for tracer advection.
	///	</summary>
	DataArray1D<double> m_vecTracersF;

	///	<summary>
	///		LU decomposition of Jacobian matrix used for tracer advection.
	///	</summary>
	DataArray2D<double> m_matTracersLUDF;

	///	<summary>
	///		Pivot matrix used for updating tracers.
	///	</summary>
	DataArray1D<int> m_vecTracersIPiv;

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
#if defined(USE_DIRECTSOLVE_APPROXJ) \
 || defined(USE_DIRECTSOLVE)
private:
	///	<summary>
	///		Jacobian matrix used in direct solve.
	///	</summary>
	DataArray2D<double> m_matJacobianF;

	///	<summary>
	///		Pivot matrix used in direct solve.
	///	</summary>
	DataArray1D<int> m_vecIPiv;
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

