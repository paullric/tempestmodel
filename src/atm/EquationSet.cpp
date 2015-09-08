///////////////////////////////////////////////////////////////////////////////
///
///	\file    EquationSet.cpp
///	\author  Paul Ullrich
///	\version May 1, 2015
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

#include "Defines.h"
#include "EquationSet.h"
#include "PhysicalConstants.h"

///////////////////////////////////////////////////////////////////////////////

EquationSet::EquationSet(
	Type eEquationSetType
) :
	m_eEquationSetType(eEquationSetType),
	m_nTracers(0)
{
	// Advection equations
	if (eEquationSetType == AdvectionEquations) {
		m_nDimensionality = 3;

		m_nComponents = 0;

	// Shallow water equations
	} else if (eEquationSetType == ShallowWaterEquations) {
		m_nDimensionality = 2;

		m_nComponents = 3;

		m_strComponentShortNames.push_back("U");
		m_strComponentShortNames.push_back("V");
		m_strComponentShortNames.push_back("H");

		m_strComponentFullNames.push_back("Alpha velocity");
		m_strComponentFullNames.push_back("Beta velocity");
		m_strComponentFullNames.push_back("Free surface height");

	// Primitive nonhydrostatic equations
	} else if (eEquationSetType == PrimitiveNonhydrostaticEquations) {
		m_nDimensionality = 3;

		m_nComponents = 5;

		m_strComponentShortNames.push_back("U");
		m_strComponentShortNames.push_back("V");

		m_strComponentFullNames.push_back("Alpha velocity");
		m_strComponentFullNames.push_back("Beta velocity");

#ifdef FORMULATION_PRESSURE
		m_strComponentShortNames.push_back("P");
		m_strComponentFullNames.push_back("Pressure");
#endif
#if defined(FORMULATION_THETA) || defined(FORMULATION_THETA_FLUX)
		m_strComponentShortNames.push_back("Theta");
		m_strComponentFullNames.push_back("Potential Temperature");
#endif
#if defined(FORMULATION_RHOTHETA_PI) || defined(FORMULATION_RHOTHETA_P)
		m_strComponentShortNames.push_back("RhoTheta");
		m_strComponentFullNames.push_back("Potential Temperature Density");
#endif

		m_strComponentShortNames.push_back("W");
		m_strComponentFullNames.push_back("Vertical velocity");

		m_strComponentShortNames.push_back("Rho");
		m_strComponentFullNames.push_back("Density");

	// Primitive nonhydrostatic equations (pressure surfaces)
	} else if (eEquationSetType == PrimitiveNonhydrostaticEquationsMassCoord) {
		m_nDimensionality = 3;

		m_nComponents = 6;

		m_strComponentShortNames.push_back("U");
		m_strComponentShortNames.push_back("V");

		m_strComponentFullNames.push_back("Alpha velocity");
		m_strComponentFullNames.push_back("Beta velocity");

		m_strComponentShortNames.push_back("Theta");
		m_strComponentFullNames.push_back("Potential Temperature");

		m_strComponentShortNames.push_back("W");
		m_strComponentFullNames.push_back("Vertical velocity");

		m_strComponentShortNames.push_back("Pressure");
		m_strComponentFullNames.push_back("Pressure");

		m_strComponentShortNames.push_back("Column Mass");
		m_strComponentFullNames.push_back("Column Mass");

	// Invalid equation set
	} else {
		_EXCEPTIONT("Invalid equation set.");
	}
}

///////////////////////////////////////////////////////////////////////////////

void EquationSet::ConvertComponents(
	const PhysicalConstants & phys,
	DataArray1D<double> & dState,
	DataArray1D<double> & dTracer
) const {

	// Indices of EquationSet variables
	const int UIx = 0;
	const int VIx = 1;
	const int HIx = 2;
	const int PIx = 2;
	const int WIx = 3;
	const int RIx = 4;

	// Shallow water equations
	if (m_eEquationSetType == ShallowWaterEquations) {
		if (dState.GetRows() != 3) {
			_EXCEPTIONT("Invalid state vector length");
		}
	}

	// Primitive non-hydrostatic equations
	if (m_eEquationSetType == PrimitiveNonhydrostaticEquations) {
		if (dState.GetRows() != 5) {
			_EXCEPTIONT("Invalid state vector length");
		}
#ifdef FORMULATION_PRESSURE
		dState[PIx] = phys.PressureFromRhoTheta(dState[PIx] * dState[RIx]);
#endif
#if defined(FORMULATION_THETA) || defined(FORMULATION_THETA_FLUX)
#endif
#if defined(FORMULATION_RHOTHETA_PI) || defined(FORMULATION_RHOTHETA_P)
		dState[PIx] *= dState[RIx];
#endif
	}
}

///////////////////////////////////////////////////////////////////////////////

