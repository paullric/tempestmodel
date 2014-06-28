///////////////////////////////////////////////////////////////////////////////
///
///	\file    WorkflowProcess.h
///	\author  Paul Ullrich
///	\version August 14, 2013
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

#ifndef _WORKFLOWPROCESS_H_
#define _WORKFLOWPROCESS_H_

#include "TimeObj.h"

#include <vector>

///////////////////////////////////////////////////////////////////////////////

class Model;

///////////////////////////////////////////////////////////////////////////////

class WorkflowProcess {

protected:
	///	<summary>
	///		Constructor.  Do not allow initialization of standalone
	///		WorkflowProcess.
	///	</summary>
	WorkflowProcess(
		Model & model,
		const Time & timeFrequency
	);

public:
	///	<summary>
	///		Virtual destructor.
	///	</summary>
	virtual ~WorkflowProcess() { }

public:
	///	<summary>
	///		Initializer.  Called prior to timestep loop.
	///	</summary>
	virtual void Initialize(
		const Time & timeStart
	);

public:
	///	<summary>
	///		Determine if this WorkflowProcess is ready to be activated.
	///	</summary>
	virtual bool IsReady(
		const Time & time
	);

public:
	///	<summary>
	///		Perform a task.
	///	</summary>
	virtual void Perform(
		const Time & time
	);

protected:
	///	<summary>
	///		Reference to the model.
	///	</summary>
	Model & m_model;

	///	<summary>
	///		Frequency of activation of this process.
	///	</summary>
	Time m_timeFrequency;

	///	<summary>
	///		Time when next activation is required.
	///	</summary>
	Time m_timeNextPerform;
};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Vector of pointers to WorkflowProcesses
///	</summary>
typedef std::vector<WorkflowProcess *> WorkflowProcessVector;

///////////////////////////////////////////////////////////////////////////////

#endif

