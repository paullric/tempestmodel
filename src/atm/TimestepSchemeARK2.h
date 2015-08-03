///////////////////////////////////////////////////////////////////////////////
///
///	\file    TimestepSchemeARK2.h
///	\author  Paul Ullrich
///	\version April 22, 2014
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

#ifndef _TIMESTEPSCHEMEARK2_H_
#define _TIMESTEPSCHEMEARK2_H_

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
class TimestepSchemeARK2 : public TimestepScheme {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	TimestepSchemeARK2(
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
	///		ARK2 parameter gamma
	///	</summary>
	static const double m_dgamma;
	///	<summary>
	///		ARK2 parameter delta
	///	</summary>
	static const double m_ddelta;

    ///	<summary>
	///		Coefficients for the time increment ARK2.
	///	</summary>
	static const double m_dTimeCf[2];

    ///	<summary>
	///		Coefficients for the explicit ARK2.
	///	</summary>
	static const double m_dExpCf[3][3];

	///	<summary>
	///		Coefficients for the explicit ARK2.
	///	</summary>
	static const double m_dImpCf[2][2];

	///	<summary>
	///		Coefficients for the implicit ARK2.
	///	</summary>
	static const double m_dBCoeff1[2];

	///	<summary>
	///		Coefficients for the explicit ARK2.
	///	</summary>
	static const double m_dBCoeff2[3];

	///	<summary>
	///		K1 vector combination (Implicit)
	///	</summary>
	DataArray1D<double> m_dK1Combo;

	///	<summary>
	///		Kh1 vector combination (Explicit)
	///	</summary>
	DataArray1D<double> m_dKh1Combo;

	///	<summary>
	///		Explicit evaluation at the 2nd substage
	///	</summary>
	DataArray1D<double> m_du2fCombo;
};

///////////////////////////////////////////////////////////////////////////////

#endif


