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
#include "DataVector.h"
#include "DataMatrix.h"
#include "DataMatrix3D.h"

///////////////////////////////////////////////////////////////////////////////

class Time;
class GridData3D;
class GridData4D;

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
		bool fNoHyperdiffusion = false
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
		bool fScaleNuLocally
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
	///		Nodal alpha fluxes.
	///	</summary>
	DataMatrix<double> m_dAlphaFlux;

	///	<summary>
	///		Nodal beta fluxes.
	///	</summary>
	DataMatrix<double> m_dBetaFlux;

	///	<summary>
	///		Nodal pointwise pressures.
	///	</summary>
	DataMatrix<double> m_dPressure;

	///	<summary>
	///		Zero vector of length RElements+1.
	///	</summary>
	DataVector<double> m_dZeroColumn;

	///	<summary>
	///		Nodal pointwise pressures in a column.
	///	</summary>
	DataVector<double> m_dColumnPressure;

	///	<summary>
	///		Nodal pointwise pressure alpha derivatives in a column.
	///	</summary>
	DataVector<double> m_dColumnDaPressure;

	///	<summary>
	///		Nodal pointwise pressure beta derivatives in a column.
	///	</summary>
	DataVector<double> m_dColumnDbPressure;

	///	<summary>
	///		Nodal pointwise pressure xi derivatives in a column.
	///	</summary>
	DataVector<double> m_dColumnDxPressure;

	///	<summary>
	///		Nodal pointwise pressure alpha derivatives in a column.
	///	</summary>
	DataVector<double> m_dColumnDaPressureREdge;

	///	<summary>
	///		Nodal pointwise pressure beta derivatives in a column.
	///	</summary>
	DataVector<double> m_dColumnDbPressureREdge;

protected:
	///	<summary>
	///		Flag indicating whether or not hyperdiffusion should be used.
	///	</summary>
	bool m_fNoHyperdiffusion;

	///	<summary>
	///		Nodal pointwise gradient.
	///	</summary>
	DataMatrix3D<double> m_dGradient;

	///	<summary>
	///		Nodal pointwise gradient of Jacobian in alpha direction (buffer).
	///	</summary>
	DataMatrix<double> m_dJGradientA;

	///	<summary>
	///		Nodal pointwise gradient of Jacobian in beta direction (buffer).
	///	</summary>
	DataMatrix<double> m_dJGradientB;

protected:
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
};

///////////////////////////////////////////////////////////////////////////////

#endif

