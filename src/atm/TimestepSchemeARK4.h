///////////////////////////////////////////////////////////////////////////////
///
///	\file    TimestepSchemeARK4.h
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

#ifndef _TIMESTEPSCHEMEARK4_H_
#define _TIMESTEPSCHEMEARK4_H_

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
class TimestepSchemeARK4 : public TimestepScheme {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	TimestepSchemeARK4(
		Model & model
	);

public:
	///	<summary>
	///		Get the number of component data instances.
	///	</summary>
	virtual int GetComponentDataInstances() const {
		return 14;
	}

	///	<summary>
	///		Get the number of tracer data instances.
	///	</summary>
	virtual int GetTracerDataInstances() const {
		return 14;
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
	///		ARK4 parameter gamma
	///	</summary>
	static const double m_dgamma;
	///	<summary>
	///		ARK4 parameter delta
	///	</summary>
	static const double m_ddelta;

    ///	<summary>
	///		Coefficients for the time increment ARK4.
	///	</summary>
	static const double m_dTimeCf[6];

    ///	<summary>
	///		Coefficients for the explicit ARK4.
	///	</summary>
	static const double m_dExpCf[7][7];

	///	<summary>
	///		Coefficients for the explicit ARK4.
	///	</summary>
	static const double m_dImpCf[7][7];

	///	<summary>
	///		K0 vector combination (Implicit)
	///	</summary>
	DataArray1D<double> m_dK0Combo;

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
	DataArray1D<double> m_du1fCombo;

	///	<summary>
	///		Explicit evaluation at the 2nd substage
	///	</summary>
	DataArray1D<double> m_du2fCombo;

	///	<summary>
	///		K2 vector combination (Implicit)
	///	</summary>
	DataArray1D<double> m_dK2Combo;

	///	<summary>
	///		Kh2 vector combination (Explicit)
	///	</summary>
	DataArray1D<double> m_dKh2Combo;

	///	<summary>
	///		Explicit evaluation at the 3rd substage
	///	</summary>
	DataArray1D<double> m_du3fCombo;

	///	<summary>
	///		K3 vector combination (Implicit)
	///	</summary>
	DataArray1D<double> m_dK3Combo;

	///	<summary>
	///		Kh3 vector combination (Explicit)
	///	</summary>
	DataArray1D<double> m_dKh3Combo;

	///	<summary>
	///		Explicit evaluation at the 4th substage
	///	</summary>
	DataArray1D<double> m_du4fCombo;

	///	<summary>
	///		K4 vector combination (Implicit)
	///	</summary>
	DataArray1D<double> m_dK4Combo;

	///	<summary>
	///		Kh4 vector combination (Explicit)
	///	</summary>
	DataArray1D<double> m_dKh4Combo;

	///	<summary>
	///		Explicit evaluation at the 5th substage
	///	</summary>
	DataArray1D<double> m_du5fCombo;

	///	<summary>
	///		K5 vector combination (Implicit)
	///	</summary>
	DataArray1D<double> m_dK5Combo;

	///	<summary>
	///		Kh5 vector combination (Explicit)
	///	</summary>
	DataArray1D<double> m_dKh5Combo;

	///	<summary>
	///		Explicit evaluation at the 6th substage
	///	</summary>
	DataArray1D<double> m_du6fCombo;
};

///////////////////////////////////////////////////////////////////////////////

#endif


