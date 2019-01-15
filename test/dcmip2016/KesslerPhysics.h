///////////////////////////////////////////////////////////////////////////////
///
///	\file    KesslerPhysics.cpp
///	\author  Antonin Verlet-Banide
///	\version May 3, 2016
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

#ifndef _KESSLERPHYSICS_H_
#define _KESSLERPHYSICS_H_

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
///		KesslerPhysics type physics forcing.
///	</summary>
class KesslerPhysics : public WorkflowProcess {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	KesslerPhysics(
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
	///		Apply KesslerPhysics physics.
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
	///		Column Exner function.
	///	</summary>
	DataArray1D<double> m_dPk;

	///	<summary>
	///		Column potential temperature.
	///	</summary>
	DataArray1D<double> m_dTheta;

	///	<summary>
	///		Column virtual potential temperature on model levels.
	///	</summary>
	DataArray1D<double> m_dThetaVNode;

	///	<summary>
	///		Column virtual potential temperature on model interfaces.
	///	</summary>
	DataArray1D<double> m_dThetaVREdge;

};

///////////////////////////////////////////////////////////////////////////////

#endif
