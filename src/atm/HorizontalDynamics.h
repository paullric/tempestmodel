///////////////////////////////////////////////////////////////////////////////
///
///	\file    HorizontalDynamics.h
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

#ifndef _HORIZONTALDYNAMICS_H_
#define _HORIZONTALDYNAMICS_H_

#include "Exception.h"

///////////////////////////////////////////////////////////////////////////////

class Time;
class Model;

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Interface for atmospheric horizontal dynamics.
///	</summary>
class HorizontalDynamics {

protected:
	///	<summary>
	///		Constructor.  Do not allow initialization of standalone
	///		HorizontalDynamics.
	///	</summary>
	HorizontalDynamics(
		Model & model
	) :
		m_model(model)
	{ }

public:
	///	<summary>
	///		Virtual destructor.
	///	</summary>
	virtual ~HorizontalDynamics() { }

public:
	///	<summary>
	///		Get the number of halo elements needed by the model.
	///	</summary>
	virtual int GetHaloElements() const {
		return 1;
	}

public:
	///	<summary>
	///		Initializer.  Called prior to execution.
	///	</summary>
	virtual void Initialize() { }

public:
	///	<summary>
	///		Perform one explicit time step.
	///	</summary>
	virtual void StepExplicit(
		int iDataInitial,
		int iDataUpdate,
		const Time & time,
		double dDeltaT
	) {
	}

	///	<summary>
	///		Perform one implicit time step.
	///	</summary>
	virtual void StepImplicit(
		int iDataInitial,
		int iDataUpdate,
		const Time & time,
		double dDeltaT
	) {
	}

	///	<summary>
	///		Number of sub-steps after sub-cycling.
	///	</summary>
	virtual int GetSubStepAfterSubCycleCount() {
		return 0;
	}

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
	) {
		_EXCEPTIONT("Not implemented");
	}

	///	<summary>
	///		Perform one time step after all sub-cycles are complete.
	///	</summary>
	virtual void StepAfterSubCycle(
		int iDataInitial,
		int iDataUpdate,
		int iDataWorking,
		const Time & time,
		double dDeltaT
	) {
	}

protected:
	///	<summary>
	///		Reference to the model.
	///	</summary>
	Model & m_model;
};

///////////////////////////////////////////////////////////////////////////////

#endif

