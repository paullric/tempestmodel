///////////////////////////////////////////////////////////////////////////////
///
///	\file    TimestepSchemeHARK332.h
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

#ifndef _TIMESTEPSCHEMEHARK332_H_
#define _TIMESTEPSCHEMEHARK332_H_

#include "TimestepScheme.h"
#include "Exception.h"
#include "DataArray1D.h"

///////////////////////////////////////////////////////////////////////////////

class Model;
class Time;

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Additive SSP IMEX RK2 by (Higueras, 2012)
///	</summary>
class TimestepSchemeHARK332 : public TimestepScheme {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	TimestepSchemeHARK332(
		Model & model
	);

public:
	///	<summary>
	///		Get the number of component data instances.
	///	</summary>
	virtual int GetComponentDataInstances() const {
		return 8;
	}

	///	<summary>
	///		Get the number of tracer data instances.
	///	</summary>
	virtual int GetTracerDataInstances() const {
		return 8;
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
	///		Coefficients for final weights ARK2.
	///	</summary>
	static const double m_dbCf[3];

        ///	<summary>
	///		Coefficients for the explicit ARK2.
	///	</summary>
	static const double m_dExpCf[3][3];

	///	<summary>
	///		Coefficients for the explicit ARK2.
	///	</summary>
	static const double m_dImpCf[3][3];

	///	<summary>
	///		K1 vector combination (Implicit)
	///	</summary>
	DataArray1D<double> m_dKv1Combo;

	///	<summary>
	///		K2 vector combination (Implicit)
	///	</summary>
	DataArray1D<double> m_dKv2Combo;

	///	<summary>
	///		K3 vector combination (Implicit)
	///	</summary>
	DataArray1D<double> m_dKv3Combo;

	///	<summary>
	///		Kh1 vector combination (Explicit)
	///	</summary>
	DataArray1D<double> m_dKh1Combo;

	///	<summary>
	///		Kh2 vector combination (Explicit)
	///	</summary>
	DataArray1D<double> m_dKh2Combo;

	///	<summary>
	///		Explicit evaluation at the 2nd stage
	///	</summary>
	DataArray1D<double> m_du2Combo;

	///	<summary>
	///		Explicit evaluation at the 2nd stage
	///	</summary>
	DataArray1D<double> m_du3Combo;

	///	<summary>
	///		Explicit evaluation at the final stage
	///	</summary>
	DataArray1D<double> m_dufCombo;
};

///////////////////////////////////////////////////////////////////////////////

#endif


