///////////////////////////////////////////////////////////////////////////////
///
///	\file    HorizontalDynamicsFEM.h
///	\author  Paul Ullrich
///	\version June 18, 2013
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

#ifndef _HORIZONTALDYNAMICSFEM_H_
#define _HORIZONTALDYNAMICSFEM_H_

#include "HorizontalDynamics.h"
#include "DataArray1D.h"
#include "DataArray2D.h"
#include "DataArray3D.h"
#include "DataArray4D.h"

///////////////////////////////////////////////////////////////////////////////

class Time;
class GridPatchGLL;

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Finite-element method (FEM) based atmospheric horizontal dynamics.
///	</summary>
class HorizontalDynamicsFEM : public HorizontalDynamics {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	HorizontalDynamicsFEM(
		Model & model,
		int nHorizontalOrder,
		int nHypervisOrder,
		double dNuScalar,
		double dNuDiv,
		double dNuVort,
		double dInstepNuDiv
	);

	///	<summary>
	///		Initializer.
	///	</summary>
	virtual void Initialize();

public:
	///	<summary>
	///		Get the number of halo elements needed by the model.
	///	</summary>
	virtual int GetHaloElements() const {
		return 1;
	}

public:
	///	<summary>
	///		Apply a positive definite filter to all tracers.
	///	</summary>
	void FilterNegativeTracers(
		int iDataUpdate
	);

	///	<summary>
	///		Perform one Forward Euler step for the interior terms of the
	///		shallow water equations.
	///	</summary>
	void StepShallowWater(
		int iDataInitial,
		int iDataUpdate,
		const Time & time,
		double dDeltaT
	);

	///	<summary>
	///		Perform one horizontal Forward Euler step for the interior terms of
	///		the non-hydrostatic primitive equations.
	///	</summary>
	void StepNonhydrostaticPrimitive(
		int iDataInitial,
		int iDataUpdate,
		const Time & time,
		double dDeltaT
	);

public:
	///	<summary>
	///		Perform one horizontal Forward Euler step.
	///	</summary>
	virtual void StepExplicit(
		int iDataInitial,
		int iDataUpdate,
		const Time & time,
		double dDeltaT
	);

protected:
	///	<summary>
	///		Apply the scalar Laplacian operator.
	///	</summary>
	void ApplyScalarHyperdiffusion(
		int iDataInitial,
		int iDataUpdate,
		double dDeltaT,
		double dNu,
		bool fScaleNuLocally,
		bool fDiffuseMass = true
	);

	///	<summary>
	///		Apply the vector Laplacian operator.
	///	</summary>
	void ApplyVectorHyperdiffusion(
		int iDataInitial,
		int iDataUpdate,
		double dDeltaT,
		double dNuDiff,
		double dNuVort,
		bool fScaleNuLocally
	);

	///	<summary>
	///		Apply Rayleigh damping.
	///	</summary>
	void ApplyRayleighFriction(
		int iDataUpdate,
		double dDeltaT
	);

public:
	///	<summary>
	///		Number of sub-steps after sub-cycling.
	///	</summary>
	virtual int GetSubStepAfterSubCycleCount();

	///	<summary>
	///		Filters, fixers and diffusion.
	///	</summar>
	virtual int SubStepAfterSubCycle(
		int iDataInitial,
		int iDataUpdate,
		int iDataWorking,
		const Time & time,
		double dDeltaT,
		int iSubStep
	);

	///	<summary>
	///		Apply hyperdiffusion.
	///	</summar>
	virtual void StepAfterSubCycle(
		int iDataInitial,
		int iDataUpdate,
		int iDataWorking,
		const Time & time,
		double dDeltaT
	);

protected:
	///	<summary>
	///		Spatial order of accuracy.
	///	</summary>
	int m_nHorizontalOrder;

	///	<summary>
	///		Nodal alpha mass fluxes.
	///	</summary>
	DataArray2D<double> m_dAlphaMassFlux;

	///	<summary>
	///		Nodal beta mass fluxes.
	///	</summary>
	DataArray2D<double> m_dBetaMassFlux;

	///	<summary>
	///		Nodal alpha pressure fluxes.
	///	</summary>
	DataArray2D<double> m_dAlphaPressureFlux;

	///	<summary>
	///		Nodal beta pressure fluxes.
	///	</summary>
	DataArray2D<double> m_dBetaPressureFlux;

	///	<summary>
	///		Nodal alpha tracer fluxes.
	///	</summary>
	DataArray3D<double> m_dAlphaTracerFlux;

	///	<summary>
	///		Nodal beta tracer fluxes.
	///	</summary>
	DataArray3D<double> m_dBetaTracerFlux;

	///	<summary>
	///		Auxiliary data within an element (on nodes).
	///	</summary>
	DataArray4D<double> m_dAuxDataNode;

	///	<summary>
	///		Auxiliary data within an element (on edges).
	///	</summary>
	DataArray4D<double> m_dAuxDataREdge;

	///	<summary>
	///		Divergence within an element (on nodes).
	///	</summary>
	DataArray3D<double> m_dDivergence;

protected:
	///	<summary>
	///		Nodal pointwise gradient of Jacobian in alpha direction (buffer).
	///	</summary>
	DataArray2D<double> m_dJGradientA;

	///	<summary>
	///		Nodal pointwise gradient of Jacobian in beta direction (buffer).
	///	</summary>
	DataArray2D<double> m_dJGradientB;

protected:
	///	<summary>
	///		Viscosity / hyperviscosity order.
	///	</summary>
	int m_nHyperviscosityOrder;

	///	<summary>
	///		Scalar hyperviscosity coefficient (at 1 degree resolution).
	///	</summary>
	double m_dNuScalar;

	///	<summary>
	///		Divergent hyperviscosity coefficient (at 1 degree resolution).
	///	</summary>
	double m_dNuDiv;

	///	<summary>
	///		Vortical hyperviscosity coefficient (at 1 degree resolution).
	///	</summary>
	double m_dNuVort;

	///	<summary>
	///		Instep divergent viscosity coefficient.
	///	</summary>
	double m_dInstepNuDiv;
};

///////////////////////////////////////////////////////////////////////////////

#endif

