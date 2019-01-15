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

#include "DataArray1D.h"
#include "DataArray2D.h"

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

	///	<summary>
	///		Initialize an identity operator.
	///	</summary>
	void InitializeIdentity(
		int nRElements
	);

	///	<summary>
	///		Compose this operator with another LinearColumnOperator.  The
	///		resulting operator will be (this) o (op)
	///	</summary>
	void ComposeWith(
		const LinearColumnOperator & op
	);

	///	<summary>
	///		Output debug output to file.
	///	</summary>
	void DebugOutput(
		const DataArray1D<double> * pREtaNode = NULL,
		const DataArray1D<double> * pREtaREdge = NULL,
		const std::string & strTag = "",
		bool fExcept = true
	);

public:
	///	<summary>
	///		Apply the column operator to a column, but only produce output on
	///		one level.  Stride indicates the sparsity of the dColumnIn array.
	///	</summary>
	inline double Apply(
		const double * dColumnIn,
		int iRout
	) const {
#ifdef DEBUG
		if (!m_fInitialized) {
			_EXCEPTIONT("Attempting to Apply uninitialized LinearColumnOperator");
		}
#endif

		double dOut = 0.0;
		for (int l = m_iBegin[iRout]; l < m_iEnd[iRout]; l++) {
			dOut += m_dCoeff[iRout][l] * dColumnIn[l];
		}
		return dOut;
	}

	///	<summary>
	///		Apply the column operator to a column, but only produce output on
	///		one level.  Stride indicates the sparsity of the dColumnIn array.
	///	</summary>
	inline double Apply(
		const double * dColumnIn,
		int iRout,
		int nStride
	) const {
#ifdef DEBUG
		if (!m_fInitialized) {
			_EXCEPTIONT("Attempting to Apply uninitialized LinearColumnOperator");
		}
#endif

		double dOut = 0.0;
		for (int l = m_iBegin[iRout]; l < m_iEnd[iRout]; l++) {
			int lx = l * nStride;

			dOut += m_dCoeff[iRout][l] * dColumnIn[lx];
		}
		return dOut;
	}

	///	<summary>
	///		Apply the column operator to a column, but only produce output on
	///		one level.  Stride indicates the sparsity of the dColumnIn array.
	///	</summary>
	inline double Apply(
		const double * dColumnIn,
		const double * dColumnRefIn,
		double dColumnRefOut,
		int iRout,
		int nStride
	) const {
#ifdef DEBUG
		if (!m_fInitialized) {
			_EXCEPTIONT("Attempting to Apply uninitialized LinearColumnOperator");
		}
#endif

		double dOut = dColumnRefOut;
		for (int l = m_iBegin[iRout]; l < m_iEnd[iRout]; l++) {
			int lx = l * nStride;

			dOut += m_dCoeff[iRout][l] * (dColumnIn[lx] - dColumnRefIn[lx]);
		}
		return dOut;
	}

public:
	///	<summary>
	///		Apply the operator to a column.
	///	</summary>
	inline void Apply(
		const double * dColumnIn,
		double * dColumnOut
	) const {
#ifdef DEBUG
		if (!m_fInitialized) {
			_EXCEPTIONT("Attempting to Apply uninitialized LinearColumnOperator");
		}
#endif

		const int nRElementsOut = m_dCoeff.GetRows();

		for (int k = 0; k < nRElementsOut; k++) {
			dColumnOut[k] = 0.0;

			for (int l = m_iBegin[k]; l < m_iEnd[k]; l++) {
				dColumnOut[k] += m_dCoeff[k][l] * dColumnIn[l];
			}
		}
	}

	///	<summary>
	///		Apply the operator to a column.  Stride indicates the sparsity
	///		of the column array.
	///	</summary>
	inline void Apply(
		const double * dColumnIn,
		double * dColumnOut,
		int nStrideIn,
		int nStrideOut
	) const {
#ifdef DEBUG
		if (!m_fInitialized) {
			_EXCEPTIONT("Attempting to Apply uninitialized LinearColumnOperator");
		}
#endif

		const int nRElementsOut = m_dCoeff.GetRows();

		for (int k = 0; k < nRElementsOut; k++) {
			int kx = k * nStrideOut;

			dColumnOut[kx] = 0.0;

			for (int l = m_iBegin[k]; l < m_iEnd[k]; l++) {
				int lx = l * nStrideIn;

				dColumnOut[kx] += m_dCoeff[k][l] * dColumnIn[lx];
			}
		}
	}

	///	<summary>
	///		Apply the operator to a column with the reference state.
	///	</summary>
	inline void ApplyWithRef(
		const double * dColumnIn,
		const double * dColumnRefIn,
		double * dColumnOut,
		const double * dColumnRefOut
	) const {
#ifdef DEBUG
		if (!m_fInitialized) {
			_EXCEPTIONT("Attempting to Apply uninitialized LinearColumnOperator");
		}
#endif

		const int nRElementsOut = m_dCoeff.GetRows();

		for (int k = 0; k < nRElementsOut; k++) {
			dColumnOut[k] = dColumnRefOut[k];

			for (int l = m_iBegin[k]; l < m_iEnd[k]; l++) {
				dColumnOut[k] +=
					m_dCoeff[k][l] * (dColumnIn[l] - dColumnRefIn[l]);
			}
		}
	}

	///	<summary>
	///		Apply the operator to a column with the reference state.  Stride
	///		indicates the sparsity of the column array.
	///	</summary>
	inline void ApplyWithRef(
		const double * dColumnIn,
		const double * dColumnRefIn,
		double * dColumnOut,
		const double * dColumnRefOut,
		int nStrideIn,
		int nStrideOut
	) const {
#ifdef DEBUG
		if (!m_fInitialized) {
			_EXCEPTIONT("Attempting to Apply uninitialized LinearColumnOperator");
		}
#endif

		const int nRElementsOut = m_dCoeff.GetRows();

		for (int k = 0; k < nRElementsOut; k++) {
			int kx = k * nStrideOut;

			dColumnOut[kx] = dColumnRefOut[kx];

			for (int l = m_iBegin[k]; l < m_iEnd[k]; l++) {
				int lx = l * nStrideIn;

				dColumnOut[kx] +=
					m_dCoeff[k][l] * (dColumnIn[lx] - dColumnRefIn[lx]);
			}
		}
	}

public:
	///	<summary>
	///		Get the matrix of linear operator coefficients.
	///	</summary>
	DataArray2D<double> & GetCoeffs() {
		return m_dCoeff;
	}

	///	<summary>
	///		Get the matrix of linear operator coefficients.
	///	</summary>
	const DataArray2D<double> & GetCoeffs() const {
		return m_dCoeff;
	}

	///	<summary>
	///		Get the vector of begin indices.
	///	</summary>
	DataArray1D<int> & GetIxBegin() {
		return m_iBegin;
	}

	///	<summary>
	///		Get the vector of begin indices.
	///	</summary>
	const DataArray1D<int> & GetIxBegin() const {
		return m_iBegin;
	}

	///	<summary>
	///		Get the vector of end indices.
	///	</summary>
	DataArray1D<int> & GetIxEnd() {
		return m_iEnd;
	}

	///	<summary>
	///		Get the vector of end indices.
	///	</summary>
	const DataArray1D<int> & GetIxEnd() const {
		return m_iEnd;
	}

protected:
	///	<summary>
	///		A flag indicating this LinearColumnOperator has been initialized.
	///	</summary>
	bool m_fInitialized;

	///	<summary>
	///		Matrix of linear operator coefficients
	///	</summary>
	DataArray2D<double> m_dCoeff;

	///	<summary>
	///		Vector of begin indices.
	///	</summary>
	DataArray1D<int> m_iBegin;

	///	<summary>
	///		Vector of end indices.
	///	</summary>
	DataArray1D<int> m_iEnd;

};

///////////////////////////////////////////////////////////////////////////////

#endif

