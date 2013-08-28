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

#include "EquationSet.h"
#include "PhysicalConstants.h"
#include "GridStaggering.h"
#include "TimestepScheme.h"
#include "HorizontalDynamics.h"
#include "VerticalDynamics.h"

///////////////////////////////////////////////////////////////////////////////

class Grid;
class TestCase;
class OutputManager;

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Relevant paramters for initialization and execution of the model.
///	</summary>
class ModelParameters {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	ModelParameters() :
		m_dDeltaT(0.0),
		m_dBeginTime(0.0),
		m_dEndTime(0.0)
	{ }

	///	<summary>
	///		Virtual destructor.
	///	</summary>
	virtual ~ModelParameters() {
	}
/*
public:
	///	<summary>
	///		Verify model parameters
	///	</summary>
	virtual void Verify()
*/
public:
	///	<summary>
	///		Time step size.
	///	</summary>
	double m_dDeltaT;

	///	<summary>
	///		Simulation begin time.
	///	</summary>
	double m_dBeginTime;

	///	<summary>
	///		Simulation end time.
	///	</summary>
	double m_dEndTime;
};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Master class for atmospheric model.
///	</summary>
class Model {

public:
	///	<summary>
	///		Vector of pointers to output managers.
	///	</summary>
	typedef std::vector<OutputManager*> OutputManagerVector;

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	Model(
		EquationSet::Type eEquationSetType
	);

public:
	///	<summary>
	///		Virtual destructor.
	///	</summary>
	virtual ~Model() { }

public:
	///	<summary>
	///		Set the model parameters.
	///	</summary>
	void SetParameters(ModelParameters * pParam);

	///	<summary>
	///		Set the grid.
	///	</summary>
	void SetGrid(Grid * pGrid);

	///	<summary>
	///		Set the timestep scheme.
	///	</summary>
	void SetTimestepScheme(TimestepScheme * pTimestepScheme);

	///	<summary>
	///		Set the horizontal dynamics.
	///	</summary>
	void SetHorizontalDynamics(HorizontalDynamics * pHorizontalDynamics);

	///	<summary>
	///		Set the vertical dynamics.
	///	</summary>
	void SetVerticalDynamics(VerticalDynamics * pVerticalDynamics);

	///	<summary>
	///		Set the test case.
	///	</summary>
	void SetTestCase(TestCase * pTestCase);

	///	<summary>
	///		Attach an output manager to this model.
	///	</summary>
	void AttachOutputManager(OutputManager * pOutMan);

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

	///	<summary>
	///		Get the number of auxiliary data objects required by
	///		HorizontalDynamics.
	///	</summary>
	int GetHorizontalDynamicsAuxDataCount() const {
		if (m_pHorizontalDynamics == NULL) {
			_EXCEPTIONT("HorizontalDynamics not initialized");
		}

		return 0;
	}

	///	<summary>
	///		Get the number of auxiliary data objects required by
	///		VerticalDynamics.
	///	</summary>
	int GetVerticalDynamicsAuxDataCount() const {
		if (m_pVerticalDynamics == NULL) {
			_EXCEPTIONT("VerticalDynamics not initialized");
		}

		return m_pVerticalDynamics->GetAuxDataCount();
	}

protected:
	///	<summary>
	///		Perform one time step.
	///	</summary>
	virtual void Step(
		bool fFirstStep,
		bool fLastStep,
		double dTime,
		double dDeltaT
	) {
	}

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
	HorizontalDynamics * GetHorizontalDynamics() {
		return m_pHorizontalDynamics;
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
	///		Get the equation set being solved.
	///	</summary>
	const EquationSet & GetEquationSet() const {
		return m_eqn;
	}

protected:
	///	<summary>
	///		Pointer to model parameters.
	///	</summary>
	ModelParameters * m_pParam;

	///	<summary>
	///		Pointer to grid
	///	</summary>
	Grid * m_pGrid;

	///	<summary>
	///		Pointer to test case.
	///	</summary>
	TestCase * m_pTestCase;

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
	///		Vector of output managers.
	///	</summary>
	OutputManagerVector m_vecOutMan;

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
	double m_dTime;
};

///////////////////////////////////////////////////////////////////////////////

#endif

