///////////////////////////////////////////////////////////////////////////////
///
///	\file    TimestepSchemeGARK2.h
///	\author  Jorge Guerra
///	\version January 29, 2016
///
///	<remarks>
///		Copyright 2000-2010 Paul Ullrich, Jorge Guerra
///
///		This file is distributed as part of the Tempest source code package.
///		Permission is granted to use, copy, modify and distribute this
///		source code and its documentation under the terms of the GNU General
///		Public License.  This software is provided "as is" without express
///		or implied warranty.
///	</remarks>

#ifndef _TIMESTEPSCHEMEGARK2_H_
#define _TIMESTEPSCHEMEGARK2_H_

#include "TimestepScheme.h"
#include "Exception.h"
#include "DataArray1D.h"

///////////////////////////////////////////////////////////////////////////////

class Model;
class Time;

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Adaptive Fourth-Order Runge-Kutta time stepping.
///	</summary>
class TimestepSchemeGARK2 : public TimestepScheme {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	TimestepSchemeGARK2(
		Model & model
	);

public:
	///	<summary>
	///		Get the number of component data instances.
	///	</summary>
	virtual int GetComponentDataInstances() const {
		return 5;
	}

	///	<summary>
	///		Get the number of tracer data instances.
	///	</summary>
	virtual int GetTracerDataInstances() const {
		return 5;
	}

protected:
	///	<summary>
	///		Perform one time step.
	///	</summary>
	virtual void Step(
		bool fFirstStep,
		bool fLastStep,
		const Time & time,
		double dDeltaT
	);

private:
	///	<summary>
	///		GARK2 parameter gamma
	///	</summary>
	static const double m_dgamma;
	///	<summary>
	///		GARK2 parameter alpha
	///	</summary>
	static const double m_dalpha;


	///	<summary>
	///		Coefficients for the explicit GARK2.
	///	</summary>
	static const double m_dExpCf[2][2];

	///	<summary>
	///		Coefficients for the explicit GARK2.
	///	</summary>
	static const double m_dImpCf[2][2];

	///	<summary>
	///		Coefficients for the explicit to implicit GARK2.
	///	</summary>
	static const double m_dEICf[2][2];

	///	<summary>
	///		Coefficients for the implicit to explicit GARK2.
	///	</summary>
	static const double m_dIECf[2][2];

	///	<summary>
	///		Y evaluation at the 2nd substage
	///	</summary>
	DataArray1D<double> m_du2fCombo;

	///	<summary>
	///		Z evaluation at the 2nd substage
	///	</summary>
	DataArray1D<double> m_du3fCombo;
};

///////////////////////////////////////////////////////////////////////////////

#endif


