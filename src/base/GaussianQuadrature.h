///////////////////////////////////////////////////////////////////////////////
///
///	\file    GaussianQuadrature.h
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

#ifndef _GAUSSIANQUADRATURE_H_
#define _GAUSSIANQUADRATURE_H_

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		GaussianQuadrature provides functionality for computing accurate
///		quadrature integrals.
///	</summary>
class GaussianQuadrature {

public:
	///	<summary>
	///		Get the number of Gauss points and weight needed for the specified
	///		order of accuracy.
	///	</summary>
	static int GetGaussPointCount(
		int nOrder
	);

	///	<summary>
	///		Get the Gauss points and weights for Gaussian quadrature of the
	///		given order of accuracy.
	///	</summary>
	static void GetGaussPointsWeights(
		int nOrder,
		double * dG,
		double * dW
	);
};

///////////////////////////////////////////////////////////////////////////////

template <class Functor>
class GaussianQuadrature1 {

public:
inline static double Integrate(
	const Functor &F,
	double dA,
	double dB
) {
	return (dB - dA) * F(0.5 * (dA + dB));
}
};

///////////////////////////////////////////////////////////////////////////////

template <class Functor>
class GaussianQuadrature2 {

public:
inline static double Integrate(
	const Functor &F,
	double dA,
	double dB
) {
	static const double dX = 0.57735026918962576451;

	// Remap to [-1,1]
	double dDiff = 0.5 * (dB - dA);
	double dAve = 0.5 * (dA + dB);

	// Calculate weighted sum
	return 0.5 * (dB - dA) * (F(-dDiff * dX + dAve) + F(dDiff * dX + dAve));
}
};

///////////////////////////////////////////////////////////////////////////////

template <class Functor>
class GaussianQuadrature3 {

public:
inline static double Integrate(
	const Functor &F,
	double dA,
	double dB
) {
	static const double dX  = 0.77459666924148337704;
	static const double dW1 = 0.55555555555555555556;
	static const double dW2 = 0.88888888888888888889;

	// Remap to [-1,1]
	double dDiff = 0.5 * (dB - dA);
	double dAve = 0.5 * (dA + dB);

	// Calculate weighted sum
	return
		0.5 * (dB - dA) * (
			dW1 * F(-dDiff * dX + dAve) +
			dW2 * F(dAve) +
			dW1 * F(dDiff * dX + dAve));
}
};

///////////////////////////////////////////////////////////////////////////////

template <class Functor>
class GaussianQuadrature4 {

public:
inline static double Integrate(
	const Functor &F,
	double dA,
	double dB
) {
	static const double dX1 = 0.86113631159405257524;
	static const double dX2 = 0.33998104358485626481;
	static const double dW1 = 0.34785484513745385737;
	static const double dW2 = 0.65214515486254614263;

	// Remap to [-1,1]
	double dDiff = 0.5 * (dB - dA);
	double dAve = 0.5 * (dA + dB);

	// Calculate weighted sum
	return
		0.5 * (dB - dA) * (
			dW1 * F(-dDiff * dX1 + dAve) +
			dW2 * F(-dDiff * dX2 + dAve) +
			dW2 * F(dDiff * dX2 + dAve) +
			dW1 * F(dDiff * dX1 + dAve));
}
};

///////////////////////////////////////////////////////////////////////////////

template <class Functor>
class GaussianQuadrature5 {

public:
inline static double Integrate(
	const Functor &F,
	double dA,
	double dB
) {
	static const double dX1 = 0.90617984593866399282;
	static const double dX2 = 0.53846931010568309105;
	static const double dW1 = 0.23692688505618908752;
	static const double dW2 = 0.47862867049936646804;
	static const double dW3 = 0.56888888888888888889;

	// Remap to [-1,1]
	double dDiff = 0.5 * (dB - dA);
	double dAve = 0.5 * (dA + dB);

	// Calculate weighted sum
	return
		0.5 * (dB - dA) * (
			dW1 * F(-dDiff * dX1 + dAve) +
			dW2 * F(-dDiff * dX2 + dAve) +
			dW3 * F(dAve) +
			dW2 * F(dDiff * dX2 + dAve) +
			dW1 * F(dDiff * dX1 + dAve));
}
};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A wrapper for a functor of two variables so that it may be called
///		with only one variable.
///	</summary>
template <class Functor>
class TwoParameterFunctorWrapper {
	public:
		///	<summary>
		///		Constructor.
		///	</summary>
		TwoParameterFunctorWrapper(const Functor &F, double dY = 0.0) :
			m_F(F),
			m_dY(dY)
		{ }

	public:
		///	<summary>
		///		Primary functor operator.
		///	</summary>
		double operator()(double dX) const {
			return m_F(dX, m_dY);
		}

	public:
		///	<summary>
		///		Set the Y value for all future function calls.
		///	</summary>
		void SetY(double dY) {
			m_dY = dY;
		}

	private:
		///	<summary>
		///		Y value for all function calls.
		///	</sumamry>
		double m_dY;

		///	<summary>
		///		Reference to a function of two parameters.
		///	</summary>
		const Functor &m_F;
};

///////////////////////////////////////////////////////////////////////////////

template <class Functor>
class GaussianQuadrature2D1 {

public:
inline static double Integrate(
	const Functor &F,
	double dX1,
	double dX2,
	double dY1,
	double dY2
) {
	return (dX2 - dX1) * (dY2 - dY1) * F(0.5 * (dX1 + dX2), 0.5 * (dY1 + dY2));
}
};

///////////////////////////////////////////////////////////////////////////////

template <class Functor>
class GaussianQuadrature2D2 {

public:
inline static double Integrate(
	const Functor &F,
	double dX1,
	double dX2,
	double dY1,
	double dY2
) {
	static const double dP = 0.57735026918962576451;

	double dDiff = 0.5 * (dY2 - dY1);
	double dAve  = 0.5 * (dY2 + dY1);

	double dTotal = 0.0;

	TwoParameterFunctorWrapper<Functor> wF(F);

	wF.SetY(-dP * dDiff + dAve);
	dTotal += GaussianQuadrature2< TwoParameterFunctorWrapper<Functor> >
	          ::Integrate(wF, dX1, dX2);

	wF.SetY(dP * dDiff + dAve);
	dTotal += GaussianQuadrature2< TwoParameterFunctorWrapper<Functor> >
	          ::Integrate(wF, dX1, dX2);

	return 0.5 * (dY2 - dY1) * dTotal;
}
};

///////////////////////////////////////////////////////////////////////////////

template <class Functor>
class GaussianQuadrature2D3 {

public:
inline static double Integrate(
	const Functor &F,
	double dX1,
	double dX2,
	double dY1,
	double dY2
) {
	static const double dP  = 0.77459666924148337704;
	static const double dW1 = 0.55555555555555555556;
	static const double dW2 = 0.88888888888888888889;

	double dDiff = 0.5 * (dY2 - dY1);
	double dAve  = 0.5 * (dY2 + dY1);

	double dTotal = 0.0;

	TwoParameterFunctorWrapper<Functor> wF(F);

	wF.SetY(-dP * dDiff + dAve);
	dTotal += dW1 * GaussianQuadrature3< TwoParameterFunctorWrapper<Functor> >
	          ::Integrate(wF, dX1, dX2);

	wF.SetY(dAve);
	dTotal += dW2 * GaussianQuadrature3< TwoParameterFunctorWrapper<Functor> >
	          ::Integrate(wF, dX1, dX2);

	wF.SetY(dP * dDiff + dAve);
	dTotal += dW1 * GaussianQuadrature3< TwoParameterFunctorWrapper<Functor> >
	          ::Integrate(wF, dX1, dX2);

	return 0.5 * (dY2 - dY1) * dTotal;
}
};

///////////////////////////////////////////////////////////////////////////////

template <class Functor>
class GaussianQuadrature2D4 {

public:
inline static double Integrate(
	const Functor &F,
	double dX1,
	double dX2,
	double dY1,
	double dY2
) {
	static const double dP1 = 0.86113631159405257524;
	static const double dP2 = 0.33998104358485626481;
	static const double dW1 = 0.34785484513745385737;
	static const double dW2 = 0.65214515486254614263;

	double dDiff = 0.5 * (dY2 - dY1);
	double dAve  = 0.5 * (dY2 + dY1);

	double dTotal = 0.0;

	TwoParameterFunctorWrapper<Functor> wF(F);

	wF.SetY(-dP1 * dDiff + dAve);
	dTotal += dW1 * GaussianQuadrature4< TwoParameterFunctorWrapper<Functor> >
	          ::Integrate(wF, dX1, dX2);

	wF.SetY(-dP2 * dDiff + dAve);
	dTotal += dW2 * GaussianQuadrature4< TwoParameterFunctorWrapper<Functor> >
	          ::Integrate(wF, dX1, dX2);

	wF.SetY(dP2 * dDiff + dAve);
	dTotal += dW2 * GaussianQuadrature4< TwoParameterFunctorWrapper<Functor> >
	          ::Integrate(wF, dX1, dX2);

	wF.SetY(dP1 * dDiff + dAve);
	dTotal += dW1 * GaussianQuadrature4< TwoParameterFunctorWrapper<Functor> >
	          ::Integrate(wF, dX1, dX2);

	return 0.5 * (dY2 - dY1) * dTotal;
}
};

///////////////////////////////////////////////////////////////////////////////

template <class Functor>
class GaussianQuadrature2D5 {

public:
inline static double Integrate(
	const Functor &F,
	double dX1,
	double dX2,
	double dY1,
	double dY2
) {
	static const double dP1 = 0.90617984593866399282;
	static const double dP2 = 0.53846931010568309105;
	static const double dW1 = 0.23692688505618908752;
	static const double dW2 = 0.47862867049936646804;
	static const double dW3 = 0.56888888888888888889;

	double dDiff = 0.5 * (dY2 - dY1);
	double dAve  = 0.5 * (dY2 + dY1);

	double dTotal = 0.0;

	TwoParameterFunctorWrapper<Functor> wF(F);

	wF.SetY(-dP1 * dDiff + dAve);
	dTotal += dW1 * GaussianQuadrature5< TwoParameterFunctorWrapper<Functor> >
	          ::Integrate(wF, dX1, dX2);

	wF.SetY(-dP2 * dDiff + dAve);
	dTotal += dW2 * GaussianQuadrature5< TwoParameterFunctorWrapper<Functor> >
	          ::Integrate(wF, dX1, dX2);

	wF.SetY(dAve);
	dTotal += dW3 * GaussianQuadrature5< TwoParameterFunctorWrapper<Functor> >
	          ::Integrate(wF, dX1, dX2);

	wF.SetY(dP2 * dDiff + dAve);
	dTotal += dW2 * GaussianQuadrature5< TwoParameterFunctorWrapper<Functor> >
	          ::Integrate(wF, dX1, dX2);

	wF.SetY(dP1 * dDiff + dAve);
	dTotal += dW1 * GaussianQuadrature5< TwoParameterFunctorWrapper<Functor> >
	          ::Integrate(wF, dX1, dX2);

	return 0.5 * (dY2 - dY1) * dTotal;
}
};

///////////////////////////////////////////////////////////////////////////////

#endif

