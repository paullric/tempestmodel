///////////////////////////////////////////////////////////////////////////////
///
///	\file    TimestepSchemeARK343.h
///	\author  Jorge Guerra
///	\version January 27, 2016
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

#ifndef _TIMESTEPSCHEMEARK343_H_
#define _TIMESTEPSCHEMEARK343_H_

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
class TimestepSchemeARK343 : public TimestepScheme {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	TimestepSchemeARK343(
		Model & model
	);

public:
	///	<summary>
	///		Get the number of component data instances.
	///	</summary>
	virtual int GetComponentDataInstances() const {
		return 9;
	}

	///	<summary>
	///		Get the number of tracer data instances.
	///	</summary>
	virtual int GetTracerDataInstances() const {
		return 9;
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
	///		ARK343 parameter gamma
	///	</summary>
	static const double m_dgamma;

        ///	<summary>
	///		ARK343 parameter b1
	///	</summary>
	static const double m_db1;

        ///	<summary>
	///		ARK343 parameter b2
	///	</summary>
	static const double m_db2;

        ///	<summary>
	///		Coefficients for the time increment ARK343.
	///	</summary>
	static const double m_dbCf[3];

        ///	<summary>
	///		Coefficients for the explicit ARK343.
	///	</summary>
	static const double m_dExpCf[4][4];

	///	<summary>
	///		Coefficients for the explicit ARK343.
	///	</summary>
	static const double m_dImpCf[3][3];

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
};

///////////////////////////////////////////////////////////////////////////////

#endif


