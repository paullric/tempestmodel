///////////////////////////////////////////////////////////////////////////////
///
///	\file    DifferentiableFunctor.h
///	\author  Paul Ullrich
///	\version July 26, 2010
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

#ifndef _DIFFERENTIABLEFUNCTOR_H_
#define _DIFFERENTIABLEFUNCTOR_H_

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		This class provides functions which allow a functor to be
///		numerically differentiated.
///	</summary>
class DifferentiableFunctor {

	public:
		///	<summary>
		///		Constructor.
		///	</summary>
		DifferentiableFunctor(
			double dEpsilon = 1e-5
		) :
			m_dEpsilon(dEpsilon)
		{ }

		///	<summary>
		///		Virtual destructor.
		///	</summary>
		virtual ~DifferentiableFunctor()
		{ }

	public:
		///	<summary>
		///		Function evaluation.
		///	</summary>
		virtual double operator()(double dX) const = 0;

	public:
		///	<summary>
		///		Calculate the first derivative of the functor, numerically.
		///	</summary>
		double Diff(double dX) const {
			return
				((*this)(dX + m_dEpsilon) - (*this)(dX - m_dEpsilon)) /
				(2.0 * m_dEpsilon);
		}

		///	<summary>
		///		Calculate the second derivative of this functor, numerically.
		///	</summary>
		double Diff2(double dX) const {
			return
				((*this)(dX + m_dEpsilon)
					- 2.0 * (*this)(dX)
					+ (*this)(dX - m_dEpsilon)
				) / (m_dEpsilon * m_dEpsilon);
		}

	private:
		///	<summary>
		///		A small parameter for taking numerical derivatives.
		///	</summary>
		double m_dEpsilon;
};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		This class provides functions which allow a 2D functor to be
///		numerically differentiated.
///	</summary>
class DifferentiableFunctor2D {

	public:
		///	<summary>
		///		Constructor.
		///	</summary>
		DifferentiableFunctor2D(
			double dEpsilon = 1e-5
		) :
			m_dEpsilon(dEpsilon)
		{ }

		///	<summary>
		///		Virtual destructor.
		///	</summary>
		virtual ~DifferentiableFunctor2D()
		{ }

	public:
		///	<summary>
		///		Function evaluation.
		///	</summary>
		virtual double operator()(double dX, double dY) const = 0;

	public:
		///	<summary>
		///		Calculate the first derivative of the functor, numerically.
		///	</summary>
		double DiffX(double dX, double dY) const {
			return
				((*this)(dX + m_dEpsilon, dY) - (*this)(dX - m_dEpsilon, dY)) /
				(2.0 * m_dEpsilon);
		}

		///	<summary>
		///		Calculate the first derivative of the functor, numerically.
		///	</summary>
		double DiffY(double dX, double dY) const {
			return
				((*this)(dX, dY + m_dEpsilon) - (*this)(dX, dY - m_dEpsilon)) /
				(2.0 * m_dEpsilon);
		}

		///	<summary>
		///		Calculate the second derivative of this functor, numerically.
		///	</summary>
		double DiffX2(double dX, double dY) const {
			return
				((*this)(dX + m_dEpsilon, dY)
					- 2.0 * (*this)(dX, dY)
					+ (*this)(dX - m_dEpsilon, dY)
				) / (m_dEpsilon * m_dEpsilon);
		}

		///	<summary>
		///		Calculate the second derivative of this functor, numerically.
		///	</summary>
		double DiffY2(double dX, double dY) const {
			return
				((*this)(dX, dY + m_dEpsilon)
					- 2.0 * (*this)(dX, dY)
					+ (*this)(dX, dY - m_dEpsilon)
				) / (m_dEpsilon * m_dEpsilon);
		}

	private:
		///	<summary>
		///		A small parameter for taking numerical derivatives.
		///	</summary>
		double m_dEpsilon;
};

///////////////////////////////////////////////////////////////////////////////

#endif

