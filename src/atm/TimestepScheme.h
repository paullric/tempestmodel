///////////////////////////////////////////////////////////////////////////////
///
///	\file    TimestepScheme.h
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

#ifndef _TIMESTEPSCHEME_H_
#define _TIMESTEPSCHEME_H_

#include "Exception.h"

///////////////////////////////////////////////////////////////////////////////

class Model;
class Time;

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Interface for time stepping in the atmospheric model.
///	</summary>
class TimestepScheme {

protected:
	///	<summary>
	///		Constructor.  Do not allow initialization of standalone
	///		TimestepScheme.
	///	</summary>
	TimestepScheme(
		Model & model
	) :
		m_model(model)
	{ }

public:
	///	<summary>
	///		Virtual destructor.
	///	</summary>
	virtual ~TimestepScheme() { }

public:
	///	<summary>
	///		Get the number of component data instances.
	///	</summary>
	virtual int GetComponentDataInstances() const = 0;

	///	<summary>
	///		Get the number of tracer data instances.
	///	</summary>
	virtual int GetTracerDataInstances() const = 0;

	///	<summary>
	///		Get the number of substeps.
	///	</summary>
	virtual int GetSubStepCount() const {
		return 0;
	}

public:
	///	<summary>
	///		Mixed method part.
	///	</summary>
	enum MixedMethodPart {
		DiscontinuousPart,
		ContinuousPart
	};

	///	<summary>
	///		Get the maximum stable Courant number for the explicit part of the
	///		Timesteps scheme.
	///	</summary>
	virtual double GetMaximumStableCourantNumber(
		MixedMethodPart eMixedMethodPart,
		int nOrder
	) const {
		return (0.0);
	}

public:
	///	<summary>
	///		Initializer.  Called prior to model execution.
	///	</summary>
	virtual void Initialize() { }

public:
	///	<summary>
	///		Perform one time substep.
	///	</summary>
	virtual int SubStep(
		bool fFirstStep,
		bool fLastStep,
		const Time & time,
		double dDeltaT,
		int iSubStep
	) {
		_EXCEPTIONT("Not implemented");
	}

	///	<summary>
	///		Perform one time step.
	///	</summary>
	virtual void Step(
		bool fFirstStep,
		bool fLastStep,
		const Time & time,
		double dDeltaT
	) = 0;

protected:
	///	<summary>
	///		Reference to the model.
	///	</summary>
	Model & m_model;
};

///////////////////////////////////////////////////////////////////////////////

#endif

