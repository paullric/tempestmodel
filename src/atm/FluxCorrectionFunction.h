///////////////////////////////////////////////////////////////////////////////
///
///	\file    FluxCorrectionFunction.h
///	\author  Paul Ullrich
///	\version September 11, 2013
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

#ifndef _FLUXCORRECTIONFUNCTION_H_
#define _FLUXCORRECTIONFUNCTION_H_

#include "DataArray1D.h"

///////////////////////////////////////////////////////////////////////////////

class FluxCorrectionFunction {

public:
	///	<summary>
	///		Compute the derivatives of the right flux correction function of
	///		the specified type at the specified nodes.
	///	</summary>
	///	<param name="iType">
	///		Type of flux correction function.
	///	</param>
	///	<param name="nOrder">
	///		Order of the corresponding discontinuous Galerkin method.
	///	</summary>
	///	<param name="dNodes">
	///		Set of nodes defined on the reference element [0,1].
	///	</param>
	///	<param name="dDeriv">
	///		Derivatives on the reference element [0,1].
	///	</param>
	static void GetDerivatives(
		int iType,
		int nOrder,
		const DataArray1D<double> & dNodes,
		DataArray1D<double> & dDeriv
	);
};

///////////////////////////////////////////////////////////////////////////////

#endif

