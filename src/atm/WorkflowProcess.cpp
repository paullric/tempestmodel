///////////////////////////////////////////////////////////////////////////////
///
///	\file    WorkflowProcess.cpp
///	\author  Paul Ullrich
///	\version June 27, 2010
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

#include "WorkflowProcess.h"

///////////////////////////////////////////////////////////////////////////////

WorkflowProcess::WorkflowProcess(
	Model & model,
	const Time & timeFrequency
) :
	m_model(model),
	m_timeFrequency(timeFrequency),
	m_timeNextPerform()
{ }

///////////////////////////////////////////////////////////////////////////////

void WorkflowProcess::Initialize(
	const Time & timeStart
) {
	m_timeNextPerform = timeStart;
	m_timeNextPerform += m_timeFrequency;
}

///////////////////////////////////////////////////////////////////////////////

bool WorkflowProcess::IsReady(
	const Time & time
) {
	// Only output initial and final results if no output time is set
	if (m_timeFrequency.IsZero()) {
		return false;
	}

	// Check current time against next output time
	if (time < m_timeNextPerform) {
		return false;
	}

	return true;
}

///////////////////////////////////////////////////////////////////////////////

void WorkflowProcess::Perform(
	const Time & time
) {
	// Update the next performance time
	m_timeNextPerform += m_timeFrequency;
}

///////////////////////////////////////////////////////////////////////////////

