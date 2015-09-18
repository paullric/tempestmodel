///////////////////////////////////////////////////////////////////////////////
///
///	\file    HorizontalDynamicsDG.h
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

#ifndef _HORIZONTALDYNAMICSDG_H_
#define _HORIZONTALDYNAMICSDG_H_

#include "HorizontalDynamicsFEM.h"
#include "DataArray1D.h"
#include "DataArray2D.h"
#include "DataArray3D.h"
#include "DataArray4D.h"

///////////////////////////////////////////////////////////////////////////////

class Time;

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Discontinuous Galerkin (DG) based atmospheric horizontal dynamics.
///	</summary>
class HorizontalDynamicsDG : public HorizontalDynamicsFEM {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	HorizontalDynamicsDG(
		Model & model,
		int nHorizontalOrder,
		int nHyperviscosityOrder,
		double dNuScalar,
		double dNuDiv,
		double dNuVort,
		double dInstepNuDiv
	);

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
	///		Perform one Forward Euler step for the element fluxes of the
	///		shallow water equations.
	///	</summary>
	void ElementFluxesShallowWater(
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
	///		Initialize the application of hyperdiffusion to the boundary.
	///	</summary>
	void InitializeApplyHyperdiffusionToBoundary(
		int iDataInitial
	);

	///	<summary>
	///		Apply the scalar Laplacian operator across element boundaries.
	///	</summary>
	void ApplyScalarHyperdiffusionToBoundary(
		int iDataState,
		int iDataUpdate,
		double dDeltaT,
		double dNu,
		bool fScaleNuLocally
	);

	///	<summary>
	///		Apply the vector Laplacian operator across element boundaries.
	///	</summary>
	void ApplyVectorHyperdiffusionToBoundary(
		int iDataState,
		int iDataUpdate,
		double dDeltaT,
		double dNuDiv,
		double dNuVort,
		bool fScaleNuLocally
	);

	///	<summary>
	///		Finalize the application of hyperdiffusion to the boundary.
	///	</summary>
	void FinalizeApplyHyperdiffusionToBoundary(
		int iDataState,
		int iDataUpdate
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
};

///////////////////////////////////////////////////////////////////////////////

#endif

