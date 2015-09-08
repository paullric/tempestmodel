///////////////////////////////////////////////////////////////////////////////
///
///	\file    EquationSet.h
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

#ifndef _EQUATIONSET_H_
#define _EQUATIONSET_H_

#include "Exception.h"
#include "DataArray1D.h"

#include <vector>
#include <string>

///////////////////////////////////////////////////////////////////////////////

class PhysicalConstants;

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Meta data describing the equation set being solved.
///	</summary>
class EquationSet {

public:
	///	<summary>
	///		Type of equations described by this equation set.
	///	</summary>
	typedef int Type;

	///	<summary>
	///		Advection equations.
	///	</summary>
	static const Type AdvectionEquations = 0;

	///	<summary>
	///		Shallow-water equations.
	///	</summary>
	static const Type ShallowWaterEquations = 1;

	///	<summary>
	///		Primitive nonhydrostatic equations.
	///	</summary>
	static const Type PrimitiveNonhydrostaticEquations = 2;

	///	<summary>
	///		Primitive nonhydrostatic equations (mass coordinate)
	///	</summary>
	static const Type PrimitiveNonhydrostaticEquationsMassCoord = 3;

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	EquationSet(Type eEquationSetType);

public:
	///	<summary>
	///		Get the name of the equation set.
	///	</summary>
	std::string GetName() const {
		if (m_eEquationSetType == AdvectionEquations) {
			return std::string("advection");
		} else if (m_eEquationSetType == ShallowWaterEquations) {
			return std::string("shallowwater");
		} else if (m_eEquationSetType == PrimitiveNonhydrostaticEquations) {
			return std::string("primitivenonhydro");
		} else {
			_EXCEPTIONT("Invalid equation set.");
		}
	}

public:
	///	<summary>
	///		Insert a new tracer variable.
	///	</summary>
	void InsertTracer(
		std::string strTracerShortName,
		std::string strTracerFullName
	) {
		m_nTracers++;
		m_strTracerShortNames.push_back(strTracerShortName);
		m_strTracerFullNames.push_back(strTracerFullName);
	}

public:
	///	<summary>
	///		Get the equation set type.
	///	</summary>
	inline Type GetType() const {
		return m_eEquationSetType;
	}

	///	<summary>
	///		Get the dimensionality of the problem.
	///	</summary>
	inline int GetDimensionality() const {
		return m_nDimensionality;
	}

	///	<summary>
	///		Get the number of components.
	///	</summary>
	inline int GetComponents() const {
		return m_nComponents;
	}

	///	<summary>
	///		Get the short name of the specified component.
	///	</summary>
	inline const std::string & GetComponentShortName(int ix) const {
		return m_strComponentShortNames[ix];
	}

	///	<summary>
	///		Get the full name of the specified component.
	///	</summary>
	inline const std::string & GetComponentFullName(int ix) const {
		return m_strComponentFullNames[ix];
	}

	///	<summary>
	///		Get the number of tracers.
	///	</summary>
	inline int GetTracers() const {
		return m_nTracers;
	}

	///	<summary>
	///		Get the short name of the specified tracer.
	///	</summary>
	inline const std::string & GetTracerShortName(int ix) const {
		return m_strTracerShortNames[ix];
	}

	///	<summary>
	///		Get the full name of the specified tracer.
	///	</summary>
	inline const std::string & GetTracerFullName(int ix) const {
		return m_strTracerFullNames[ix];
	}

public:
	///	<summary>
	///		Convert output from TestCase::EvaluatePointwiseState to
	///		state variables consistent with this EquationSet.
	///	</summary>
	void ConvertComponents(
		const PhysicalConstants & phys,
		DataArray1D<double> & dState,
		DataArray1D<double> & dTracer
	) const;

private:
	///	<summary>
	///		Type of equations described by this equation set.
	///	</summary>
	Type m_eEquationSetType;

	///	<summary>
	///		Dimensionality of the problem.
	///	</summary>
	int m_nDimensionality;

	///	<summary>
	///		Number of components per node.
	///	</summary>
	int m_nComponents;

	///	<summary>
	///		Short name of each component.
	///	</summary>
	std::vector<std::string> m_strComponentShortNames;

	///	<summary>
	///		Full name of each component.
	///	</summary>
	std::vector<std::string> m_strComponentFullNames;

	///	<summary>
	///		Number of tracers per node.
	///	</summary>
	int m_nTracers;

	///	<summary>
	///		Short name of each component.
	///	</summary>
	std::vector<std::string> m_strTracerShortNames;

	///	<summary>
	///		Full name of each component.
	///	</summary>
	std::vector<std::string> m_strTracerFullNames;

};

///////////////////////////////////////////////////////////////////////////////

#endif

