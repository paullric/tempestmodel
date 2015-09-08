///////////////////////////////////////////////////////////////////////////////
///
///	\file    LinearColumnOperatorFEM.h
///	\author  Paul Ullrich
///	\version July 29, 2014
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

#ifndef _LINEARCOLUMNOPERATORFEM_H_
#define _LINEARCOLUMNOPERATORFEM_H_

#include "LinearColumnOperator.h"

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Linear interpolation operator working on a column.
///	</summary>
class LinearColumnInterpFEM : public LinearColumnOperator {

public:
	///	<summary>
	///		Source of interpolation.
	///	</summary>
	enum InterpSource {
		InterpSource_Levels,
		InterpSource_Interfaces
	};

public:
	///	<summary>
	///		Default constructor.
	///	</summary>
	LinearColumnInterpFEM() :
		LinearColumnOperator()
	{ }

	///	<summary>
	///		Initialize the operator
	///	</summary>
	///	<param name="fContinuous">
	///		Set to 'true' if a continuous element formulation is used for
	///		the input column or set to 'false' if a discontinuous element
	///		formulation is used.
	///	</param>
	void Initialize(
		InterpSource eInterpSource,
		int nVerticalOrder,
		const DataArray1D<double> & dREtaNode,
		const DataArray1D<double> & dREtaREdge,
		const DataArray1D<double> & dREtaOut,
		bool fZeroBoundaries = false
	);
};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Linear differentiation operator working on a column.
///	</summary>
class LinearColumnDiffFEM : public LinearColumnOperator {

public:
	///	<summary>
	///		Source of interpolation.
	///	</summary>
	enum InterpSource {
		InterpSource_Levels,
		InterpSource_Interfaces
	};

public:
	///	<summary>
	///		Default constructor.
	///	</summary>
	LinearColumnDiffFEM() :
		LinearColumnOperator()
	{ }

	///	<summary>
	///		Initialize the operator by interpolating to interfaces and
	///		differentiating.
	///	</summary>
	void InitializeInterfaceMethod(
		InterpSource eInterpSource,
		int nVerticalOrder,
		const DataArray1D<double> & dREtaNode,
		const DataArray1D<double> & dREtaREdge,
		const DataArray1D<double> & dREtaOut,
		bool fZeroBoundaries
	);

	///	<summary>
	///		Initialize a differentiation operator using a flux reconstruction
	///		function based approach.
	///	</summary>
	void InitializeFluxCorrectionMethod(
		InterpSource eInterpSource,
		int nVerticalOrder,
		const DataArray1D<double> & dREtaNode,
		const DataArray1D<double> & dREtaREdge,
		const DataArray1D<double> & dREtaOut,
		bool fZeroBoundaries
	);

	///	<summary>
	///		Initialize a differentiation operator applied to GLL nodes.
	///	</summary>
	void InitializeGLLNodes(
		int nVerticalOrder,
		const DataArray1D<double> & dREtaNode,
		const DataArray1D<double> & dREtaOut
	);
};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Linear second derivative operator working on a column.
///	</summary>
class LinearColumnDiffDiffFEM : public LinearColumnOperator {

public:
	///	<summary>
	///		Source of interpolation.
	///	</summary>
	enum InterpSource {
		InterpSource_Levels,
		InterpSource_Interfaces
	};

public:
	///	<summary>
	///		Default constructor.
	///	</summary>
	LinearColumnDiffDiffFEM() :
		LinearColumnOperator()
	{ }

	///	<summary>
	///		Initialize the operator
	///	</summary>
	void Initialize(
		InterpSource eInterpSource,
		int nVerticalOrder,
		const DataArray1D<double> & dREtaNode,
		const DataArray1D<double> & dREtaREdge
	);

	///	<summary>
	///		Initialize a differentiation operator applied to GLL nodes.
	///	</summary>
	void InitializeGLLNodes(
		int nVerticalOrder,
		const DataArray1D<double> & dREtaNode
	);
};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Linear integration operator working on a column.
///	</summary>
class LinearColumnIntFEM : public LinearColumnOperator {

public:
	///	<summary>
	///		Source of interpolation.
	///	</summary>
	enum InterpSource {
		InterpSource_Levels,
		InterpSource_Interfaces
	};

public:
	///	<summary>
	///		Default constructor.
	///	</summary>
	LinearColumnIntFEM() :
		LinearColumnOperator()
	{ }

	///	<summary>
	///		Initialize the operator
	///	</summary>
	void InitializeNodeToNodeInterfaceMethod(
		int nVerticalOrder,
		const DataArray1D<double> & dREtaNode,
		const DataArray1D<double> & dREtaREdge
	);
};

///////////////////////////////////////////////////////////////////////////////

#endif
