///////////////////////////////////////////////////////////////////////////////
///
///	\file    VerticalDynamics.h
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

#ifndef _VERTICALDYNAMICS_H_
#define _VERTICALDYNAMICS_H_

///////////////////////////////////////////////////////////////////////////////

class Time;
class Model;

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Interface for vertical dynamics in the atmospheric model.
///	</summary>
class VerticalDynamics {

protected:
	///	<summary>
	///		Constructor.  Do not allow initialization of standalone
	///		VerticalDynamics.
	///	</summary>
	VerticalDynamics(
		Model & model
	) :
		m_model(model)
	{ }

public:
	///	<summary>
	///		Virtual destructor.
	///	</summary>
	virtual ~VerticalDynamics() { }

public:
	///	<summary>
	///		Number of auxiliary data objects needed.
	///	</summary>
	virtual int GetAuxDataCount() const {
		return 2;
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
	///		Force a full explicit update one time only
	///	</summary>
	virtual void ForceStepExplicit(
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

	///	<summary>
	///		Filter negative tracers in all columns.
	///	</summary>
	virtual void FilterNegativeTracers(
		int iDataUpdate
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

