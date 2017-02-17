///////////////////////////////////////////////////////////////////////////////
///
///	\file    SplitExplicitDynamics.h
///	\author  Paul Ullrich
///	\version February 6, 2017
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

#ifndef _SPLITEXPLICITDYNAMICS_H_
#define _SPLITEXPLICITDYNAMICS_H_

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
class SplitExplicitDynamics : public HorizontalDynamics {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	SplitExplicitDynamics(
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
	///		Additional data instances needed for acoustic loop.
	///	</summary>
	virtual int GetComponentDataInstances() const {
		return 4;
	}

public:
	///	<summary>
	///		Apply a positive definite filter to all tracers.
	///	</summary>
	void FilterNegativeTracers(
		int iDataUpdate
	);

	///	<summary>
	///		Calculate the explicit tendencies as part of the update.
	///	</summary>
	void CalculateTendencies(
		int iDataInitial,
		int iDataTendencies,
		double dDeltaT
	);

	///	<summary>
	///		Perform the updates from the first acoustic loop.
	///	</summary>
	void FirstAcousticLoop(
		int iDataInitial,
		int iDataTendencies,
		int iDataAcoustic2,
		double dDeltaT
	);

	///	<summary>
	///		Perform the updates from the acoustic loop.
	///	</summary>
	void PerformAcousticLoop(
		int iDataInitial,
		int iDataTendencies,
		int iDataAcoustic0,
		int iDataAcoustic1,
		int iDataAcoustic2,
		double dDeltaT
	);

public:
	///	<summary>
	///		Perform one Forward Euler step.
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
		int iComponent = (-1),
		bool fRemoveRefState = false
	);

	///	<summary>
	///		Apply the vector Laplacian operator.
	///	</summary>
	void ApplyVectorHyperdiffusion(
		int iDataInitial,
		int iDataWorking,
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
	///		Off-centering pressure parameter.
	///	</summary>
	double m_dBd;

	///	<summary>
	///		Off-centering implicit update parameter.
	///	</summary>
	double m_dBs;

	///	<summary>
	///		Column-wise A terms for vertical implicit solve.
	///	</summary>
	DataArray3D<double> m_dA;

	///	<summary>
	///		Column-wise B terms for vertical implicit solve.
	///	</summary>
	DataArray3D<double> m_dB;

	///	<summary>
	///		Column-wise C terms for vertical implicit solve.
	///	</summary>
	DataArray3D<double> m_dC;

	///	<summary>
	///		Column-wise D terms for vertical implicit solve.
	///	</summary>
	DataArray3D<double> m_dD;

	///	<summary>
	///		 2D Kinetic energy.
	///	</summary>
	DataArray3D<double> m_dK2;

	///	<summary>
	///		 Vertical alpha momentum flux on interfaces.
	///	</summary>
	DataArray3D<double> m_dSDotREdge;

	///	<summary>
	///		 Vertical alpha momentum flux on interfaces.
	///	</summary>
	DataArray3D<double> m_dSDotUaREdge;

	///	<summary>
	///		 Vertical beta momentum flux on interfaces.
	///	</summary>
	DataArray3D<double> m_dSDotUbREdge;

	///	<summary>
	///		Vertical theta flux on interfaces.
	///	</summary>
	DataArray3D<double> m_dSDotThetaREdge;

	///	<summary>
	///		Vertical vertical momentum flux on levels.
	///	</summary>
	DataArray3D<double> m_dSDotWNode;

	///	<summary>
	///		 alpha contravariant velocity.
	///	</summary>
	DataArray3D<double> m_d2DConUa;

	///	<summary>
	///		 beta contravariant velocity.
	///	</summary>
	DataArray3D<double> m_d2DConUb;

	///	<summary>
	///		 alpha covariant velocity.
	///	</summary>
	DataArray3D<double> m_d2DCovUa;

	///	<summary>
	///		 beta covariant velocity.
	///	</summary>
	DataArray3D<double> m_d2DCovUb;

	///	<summary>
	///		Nodal alpha mass fluxes.
	///	</summary>
	DataArray3D<double> m_dAlphaMassFlux;

	///	<summary>
	///		Nodal beta mass fluxes.
	///	</summary>
	DataArray3D<double> m_dBetaMassFlux;

	///	<summary>
	///		Vertical mass flux on interfaces (acoustic loop 1).
	///	</summary>
	DataArray3D<double> m_dZMassFluxREdge1;

	///	<summary>
	///		Vertical mass flux on interfaces (acoustic loop 2).
	///	</summary>
	DataArray3D<double> m_dZMassFluxREdge2;

	///	<summary>
	///		Vertical momentum flux in the alpha direction.
	///	</summary>
	DataArray3D<double> m_dAlphaVerticalMomentumFluxREdge;

	///	<summary>
	///		Vertical momentum flux in the beta direction.
	///	</summary>
	DataArray3D<double> m_dBetaVerticalMomentumFluxREdge;

	///	<summary>
	///		Nodal alpha pressure fluxes.
	///	</summary>
	DataArray3D<double> m_dAlphaPressureFlux;

	///	<summary>
	///		Nodal beta pressure fluxes.
	///	</summary>
	DataArray3D<double> m_dBetaPressureFlux;

	///	<summary>
	///		Vertical pressure flux on interfaces (acoustic loop 1).
	///	</summary>
	DataArray3D<double> m_dZPressureFluxREdge1;

	///	<summary>
	///		Vertical pressure flux on interfaces (acoustic loop 2).
	///	</summary>
	DataArray3D<double> m_dZPressureFluxREdge2;

	///	<summary>
	///		Nodal mass update.
	///	</summary>
	DataArray3D<double> m_dNodalMassUpdate;

	///	<summary>
	///		Nodal pressure update.
	///	</summary>
	DataArray3D<double> m_dNodalPressureUpdate;

	///	<summary>
	///		Acoustic pressure term used in acoustic loop.
	///	</summary>
	DataArray3D<double> m_dAcousticPressureTerm;

	///	<summary>
	///		Output from the tridiagonal solve
	///	</summary>
	DataArray2D<int> m_nInfo;

	///	<summary>
	///		Buffer state.
	///	</summary>
	DataArray2D<double> m_dBufferState;

	///	<summary>
	///		Buffer J gradient, used in hyperviscosity calculation.
	///	</summary>
	DataArray2D<double> m_dJGradientA;

	///	<summary>
	///		Buffer J gradient, used in hyperviscosity calculation.
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
