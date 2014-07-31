///////////////////////////////////////////////////////////////////////////////
///
///	\file    LinearColumnOperator.h
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

#ifndef _LINEARCOLUMNOPERATOR_H_
#define _LINEARCOLUMNOPERATOR_H_

#include "DataVector.h"
#include "DataMatrix.h"

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A linear operator applied to state variables in a column.
///	</summary>
class LinearColumnOperator {

public:
	///	<summary>
	///		Default constructor.
	///	</summary>
	LinearColumnOperator();

	///	<summary>
	///		Constructor.
	///	</summary>
	LinearColumnOperator(
		int nRElementsIn,
		int nRElementsOut
	);

	///	<summary>
	///		Initialize the operator.
	///	</summary>
	void Initialize(
		int nRElementsIn,
		int nRElementsOut
	);

public:
	///	<summary>
	///		Apply the operator to a column.
	///	</summary>
	void Apply(
		const double * dColumnIn,
		double * dColumnOut
	) const;

	///	<summary>
	///		Apply the operator to a column.  Stride indicates the sparsity
	///		of the column array.
	///	</summary>
	void Apply(
		const double * dColumnIn,
		double * dColumnOut,
		int nStrideIn,
		int nStrideOut
	) const;

	///	<summary>
	///		Apply the operator to a column with the reference state.
	///	</summary>
	void Apply(
		const double * dColumnIn,
		const double * dColumnRefIn,
		double * dColumnOut,
		const double * dColumnRefOut
	) const;

	///	<summary>
	///		Apply the operator to a column with the reference state.  Stride
	///		indicates the sparsity of the column array.
	///	</summary>
	void Apply(
		const double * dColumnIn,
		const double * dColumnRefIn,
		double * dColumnOut,
		const double * dColumnRefOut,
		int nStrideIn,
		int nStrideOut
	) const;

public:
	///	<summary>
	///		Get the matrix of linear operator coefficients.
	///	</summary>
	DataMatrix<double> & GetCoeffs() {
		return m_dCoeff;
	}

	///	<summary>
	///		Get the matrix of linear operator coefficients.
	///	</summary>
	const DataMatrix<double> & GetCoeffs() const {
		return m_dCoeff;
	}

	///	<summary>
	///		Get the vector of begin indices.
	///	</summary>
	DataVector<int> & GetIxBegin() {
		return m_iBegin;
	}

	///	<summary>
	///		Get the vector of begin indices.
	///	</summary>
	const DataVector<int> & GetIxBegin() const {
		return m_iBegin;
	}

	///	<summary>
	///		Get the vector of end indices.
	///	</summary>
	DataVector<int> & GetIxEnd() {
		return m_iEnd;
	}

	///	<summary>
	///		Get the vector of end indices.
	///	</summary>
	const DataVector<int> & GetIxEnd() const {
		return m_iEnd;
	}

protected:
	///	<summary>
	///		Matrix of linear operator coefficients
	///	</summary>
	DataMatrix<double> m_dCoeff;

	///	<summary>
	///		Vector of begin indices.
	///	</summary>
	DataVector<int> m_iBegin;

	///	<summary>
	///		Vector of end indices.
	///	</summary>
	DataVector<int> m_iEnd;

};

///////////////////////////////////////////////////////////////////////////////

#endif

