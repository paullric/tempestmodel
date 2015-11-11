///////////////////////////////////////////////////////////////////////////////
///
///	\file    Model.h
///	\author  Paul Ullrich
///	\version February 25, 2013
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

#ifndef _MODEL_H_
#define _MODEL_H_

#include "TimeObj.h"
#include "EquationSet.h"
#include "PhysicalConstants.h"
#include "GridStaggering.h"
#include "Grid.h"
#include "TestCase.h"
#include "TimestepScheme.h"
#include "HorizontalDynamics.h"
#include "VerticalDynamics.h"
#include "OutputManager.h"
#include "WorkflowProcess.h"

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Relevant parameters for initialization and execution of the model.
///	</summary>
class ModelParameters {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	ModelParameters() :
		m_strRestartFile(""),
		m_timeDeltaT(),
		m_timeStart(),
		m_timeEnd()
	{ }

public:
	///	<summary>
	///		Name of the restart file to use.
	///	</summary>
	std::string m_strRestartFile;

	///	<summary>
	///		Time step size.
	///	</summary>
	Time m_timeDeltaT;

	///	<summary>
	///		Start time of the simulation.
	///	</summary>
	Time m_timeStart;

	///	<summary>
	///		End time of the simulation.
	///	</summary>
	Time m_timeEnd;
};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Master class for atmospheric model.
///	</summary>
class Model {

public:
	///	<summary>
	///		Constructor that specifies an EquationSet::Type
	///	</summary>
	Model(EquationSet::Type eEquationSetType);

	///	<summary>
	///		Constructor that specifies an EquationSet
	///	</summary>
	Model(const EquationSet & eqn);

public:
	///	<summary>
	///		Virtual destructor.
	///	</summary>
	virtual ~Model();

public:
	///	<summary>
	///		Set the model parameters.
	///	</summary>
	void SetParameters(const ModelParameters & param);

	///	<summary>
	///		Set the Grid from a pointer.  Model assumes ownership of the
	///		pointer once it is assigned.
	///	</summary>
	void SetGrid(
		Grid * pGrid,
		int nPatchCount = (-1),
		bool fInitializeConnectivity = true
	);

	///	<summary>
	///		Set the TimestepScheme from a pointer.  Model assumes ownership
	///		of the pointer once it is assigned.
	///	</summary>
	void SetTimestepScheme(TimestepScheme * pTimestepScheme);

	///	<summary>
	///		Set the HorizontalDynamics from a pointer.  Model assumes
	///		ownership of the pointer once it is assigned.
	///	</summary>
	void SetHorizontalDynamics(HorizontalDynamics * pHorizontalDynamics);

	///	<summary>
	///		Set the VerticalDynamics from a pointer.  Model assumes
	///		ownership of the pointer once it is assigned.
	///	</summary>
	void SetVerticalDynamics(VerticalDynamics * pVerticalDynamics);

	///	<summary>
	///		Set the TestCase from a pointer.  Model assumes ownership
	///		of the pointer once it is assigned.
	///	</summary>
	void SetTestCase(TestCase * pTestCase);

	///	<summary>
	///		Attach an OutputManager to this model.  Model assumes ownership
	///		of the pointer once it is assigned.
	///	</summary>
	void AttachOutputManager(OutputManager * pOutMan);

	///	<summary>
	///		Attach a WorkflowProcess to this model.  Model assumes ownership
	///		of the pointer once it is assigned.
	///	</summary>
	void AttachWorkflowProcess(WorkflowProcess * pWorkflowProcess);

protected:
	///	<summary>
	///		Evaluate the state from a Restart file.
	///	</summary>
	void EvaluateStateFromRestartFile();

public:
	///	<summary>
	///		Get the number of halo elements needed by the model.
	///	</summary>
	virtual int GetHaloElements() const {
		if (m_pHorizontalDynamics == NULL) {
			_EXCEPTIONT("HorizontalDynamics not initialized");
		}

		return m_pHorizontalDynamics->GetHaloElements();
	}

	///	<summary>
	///		Get the number of component data instances.
	///	</summary>
	int GetComponentDataInstances() const {
		if (m_pTimestepScheme == NULL) {
			_EXCEPTIONT("TimestepScheme not initialized");
		}

		return m_pTimestepScheme->GetComponentDataInstances();
	}

	///	<summary>
	///		Get the number of tracer data instances.
	///	</summary>
	int GetTracerDataInstances() const {
		if (m_pTimestepScheme == NULL) {
			_EXCEPTIONT("TimestepScheme not initialized");
		}

		return m_pTimestepScheme->GetTracerDataInstances();
	}

public:
	///	<summary>
	///		Perform one substep.
	///	</summary>
	void SubStep(
		bool fFirstStep,
		bool fLastStep,
		int iSubStep
	);

public:
	///	<summary>
	///		Begin the model.
	///	</summary>
	virtual void Go();

public:
	///	<summary>
	///		Compute error norms.
	///	</summary>
	virtual void ComputeErrorNorms();

public:
	///	<summary>
	///		Get a pointer to the grid.
	///	</summary>
	Grid * GetGrid() {
		return m_pGrid;
	}

	///	<summary>
	///		Get a pointer to the grid.
	///	</summary>
	const Grid * GetGrid() const {
		return m_pGrid;
	}

	///	<summary>
	///		Get a pointer to the TimestepScheme.
	///	</summary>
	TimestepScheme * GetTimestepScheme() const {
		return m_pTimestepScheme;
	}

	///	<summary>
	///		Get a pointer to the HorizontalDynamics.
	///	</summary>
	const HorizontalDynamics * GetHorizontalDynamics() const {
		return m_pHorizontalDynamics;
	}

	///	<summary>
	///		Get a pointer to the HorizontalDynamics.
	///	</summary>
	HorizontalDynamics * GetHorizontalDynamics() {
		return m_pHorizontalDynamics;
	}

	///	<summary>
	///		Get a pointer to the VerticalDynamics.
	///	</summary>
	const VerticalDynamics * GetVerticalDynamics() const {
		return m_pVerticalDynamics;
	}

	///	<summary>
	///		Get a pointer to the VerticalDynamics.
	///	</summary>
	VerticalDynamics * GetVerticalDynamics() {
		return m_pVerticalDynamics;
	}

	///	<summary>
	///		Get the physical constants for this model.
	///	</summary>
	const PhysicalConstants & GetPhysicalConstants() const {
		return m_phys;
	}

	///	<summary>
	///		Get the physical constants for this model.
	///	</summary>
	PhysicalConstants & GetPhysicalConstants() {
		return m_phys;
	}

	///	<summary>
	///		Get the equation set being solved.
	///	</summary>
	const EquationSet & GetEquationSet() const {
		return m_eqn;
	}

public:
	///	<summary>
	///		Get the time step size.
	///	</summary>
	const Time & GetDeltaT() const {
		return m_param.m_timeDeltaT;
	}

	///	<summary>
	///		Get the start time of the simulation.
	///	</summary>
	const Time & GetStartTime() const {
		return m_param.m_timeStart;
	}

	///	<summary>
	///		Set the start time of the simulation.
	///	</summary>
	void SetStartTime(const Time & timeStart) {
		m_param.m_timeStart = timeStart;
	}

protected:
	///	<summary>
	///		Model parameters.
	///	</summary>
	ModelParameters m_param;

	///	<summary>
	///		Pointer to grid
	///	</summary>
	Grid * m_pGrid;

	///	<summary>
	///		Pointer to timestepping scheme.
	///	</summary>
	TimestepScheme * m_pTimestepScheme;

	///	<summary>
	///		Pointer to horizontal dynamics.
	///	</summary>
	HorizontalDynamics * m_pHorizontalDynamics;

	///	<summary>
	///		Pointer to vertical dynamics.
	///	</summary>
	VerticalDynamics * m_pVerticalDynamics;

	///	<summary>
	///		Vector of WorkflowProcesses.
	///	</summary>
	WorkflowProcessVector m_vecWorkflowProcess;

	///	<summary>
	///		Vector of OutputManagers.
	///	</summary>
	OutputManagerVector m_vecOutMan;

	///	<summary>
	///		Pointer to test case.
	///	</summary>
	TestCase * m_pTestCase;

	///	<summary>
	///		Physical constants for model.
	///	</summary>
	PhysicalConstants m_phys;

	///	<summary>
	///		Equation set to solve.
	///	</summary>
	EquationSet m_eqn;

	///	<summary>
	///		Type of staggering used for each variable.
	///	</summary>
	GridStaggering m_stag;

private:
	///	<summary>
	///		Current model time.
	///	</summary>
	Time m_time;
};

///////////////////////////////////////////////////////////////////////////////

#endif

