///////////////////////////////////////////////////////////////////////////////
///
///	\file    DCMIPPhysics.cpp
///	\author  Paul Ullrich
///	\version May 9, 2016
///
///	<remarks>
///		Copyright 2000-2016 Paul Ullrich
///
///		This file is distributed as part of the Tempest source code package.
///		Permission is granted to use, copy, modify and distribute this
///		source code and its documentation under the terms of the GNU General
///		Public License.  This software is provided "as is" without express
///		or implied warranty.
///	</remarks>

#ifndef _DCMIPPHYSICS_H_
#define _DCMIPPHYSICS_H_

#include "Defines.h"
#include "DataArray1D.h"
#include "DataArray2D.h"
#include "DataArray3D.h"
#include "DataArray4D.h"

#include "WorkflowProcess.h"

///////////////////////////////////////////////////////////////////////////////

class GridPatch;
class Time;

///	<summary>
///		DCMIPPhysics type physics forcing.
///	</summary>
class DCMIPPhysics : public WorkflowProcess {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	DCMIPPhysics(
		Model & model,
		const Time & timeFrequency);
	
public:
  	///	<summary>
	///		Initializer.
	///	</summary>
	virtual void Initialize(
		const Time & timeStart
	);
	
	///	<summary>
	///		Apply DCMIPPhysics physics.
	///	</summary>
	virtual void Perform(
		const Time & time
	);
	
protected:
	///	<summary>
	///		Column water vapor mixing ratio.
	///	</summary>
	DataArray1D<double> m_dQv;

	///	<summary>
	///		Column cloud water mixing ratio.
	///	</summary>
	DataArray1D<double> m_dQc;

	///	<summary>
	///		Column rain water mixing ratio.
	///	</summary>
	DataArray1D<double> m_dQr;

	///	<summary>
	///		Column dry density.
	///	</summary>
	DataArray1D<double> m_dRho;

	///	<summary>
	///		Column model levels.
	///	</summary>
	DataArray1D<double> m_dZc;

	///	<summary>
	///		Column model interfaces.
	///	</summary>
	DataArray1D<double> m_dZi;

	///	<summary>
	///		Pressure on model levels.
	///	</summary>
	DataArray1D<double> m_dPmid;

	///	<summary>
	///		Pressure on model interfaces.
	///	</summary>
	DataArray1D<double> m_dPint;

	///	<summary>
	///		Pressure difference (layer thickness) over model levels.
	///	</summary>
	DataArray1D<double> m_dPdel;

	///	<summary>
	///		Reciprocal of pressure difference (layer thickness) over
	///		model levels.
	///	</summary>
	DataArray1D<double> m_dRPdel;

	///	<summary>
	///		Column temperature.
	///	</summary>
	DataArray1D<double> m_dT;

	///	<summary>
	///		Column zonal velocity.
	///	</summary>
	DataArray1D<double> m_dU;

	///	<summary>
	///		Column meridional velocity.
	///	</summary>
	DataArray1D<double> m_dV;

	///	<summary>
	///		Column virtual potential temperature on levels.
	///	</summary>
	DataArray1D<double> m_dThetaVNode;

	///	<summary>
	///		Column potential temperature.
	///	</summary>
	DataArray1D<double> m_dTheta;

	///	<summary>
	///		Column virtual potential temperature on model levels.
	///	</summary>
	DataArray1D<double> m_dTvNode;

	///	<summary>
	///		Column virtual potential temperature on model interfaces.
	///	</summary>
	DataArray1D<double> m_dTvREdge;

	///	<summary>
	///		Jacobian coefficients for implicit solve in the boundary layer.
	///	</summary>
	DataArray2D<double> m_dBLJacobian;

	///	<summary>
	///		Altitude of each model level.
	///	</summary>
	DataArray1D<double> m_dZNode;

	///	<summary>
	///		Altitude of each model interface.
	///	</summary>
	DataArray1D<double> m_dZREdge;

	///	<summary>
	///		Column eddy fluxes on model levels.
	///	</summary>
	DataArray2D<double> m_dEddyStateNode;

	///	<summary>
	///		Column eddy fluxes on model interfaces.
	///	</summary>
	DataArray2D<double> m_dEddyStateREdge;

	///	<summary>
	///		Column eddy fluxes on model levels.
	///	</summary>
	DataArray2D<double> m_dEddyTracerNode;

	///	<summary>
	///		Column eddy fluxes on model interfaces.
	///	</summary>
	DataArray2D<double> m_dEddyTracerREdge;

protected:
	///	<summary>
	///		Column turbulent mixing strength of velocity on interfaces.
	///	</summary>
	DataArray1D<double> m_dKm;

	///	<summary>
	///		Column turbulent mixing strength of scalars on interfaces.
	///	</summary>
	DataArray1D<double> m_dKE;

	///	<summary>
	///		Column E coefficients for boundary layer parameterization.
	///	</summary>
	DataArray1D<double> m_dEm;

	///	<summary>
	///		Column E coefficients for boundary layer parameterization.
	///	</summary>
	DataArray1D<double> m_dEE;

	///	<summary>
	///		Column F coefficients for boundary layer parameterization.
	///	</summary>
	DataArray2D<double> m_dF;

	///	<summary>
	///		Column density on interfaces.
	///	</summary>
	DataArray1D<double> m_dRhoREdge;

	///	<summary>
	///		Column potential temperature on interfaces.
	///	</summary>
	DataArray1D<double> m_dThetaREdge;

};

///////////////////////////////////////////////////////////////////////////////

#endif
