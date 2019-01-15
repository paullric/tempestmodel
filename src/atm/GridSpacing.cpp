///////////////////////////////////////////////////////////////////////////////
///
///	\file    GridSpacing.h
///	\author  Paul Ullrich
///	\version March 8, 2013
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

#include "GridSpacing.h"
#include "GaussQuadrature.h"
#include "GaussLobattoQuadrature.h"

///////////////////////////////////////////////////////////////////////////////
// GridSpacing
///////////////////////////////////////////////////////////////////////////////

GridSpacing::GridSpacing(
	double dDeltaElement,
	double dZeroCoord
) :
	m_dDeltaElement(dDeltaElement),
	m_dZeroCoord(dZeroCoord)
{ }

///////////////////////////////////////////////////////////////////////////////
// GridSpacingUniform
///////////////////////////////////////////////////////////////////////////////

GridSpacingUniform::GridSpacingUniform(
	double dDeltaElement,
	double dZeroCoord
) :
	GridSpacing(dDeltaElement, dZeroCoord)
{ }

///////////////////////////////////////////////////////////////////////////////

double GridSpacingUniform::GetNode(int ix) const {
	return (m_dZeroCoord
		+ (static_cast<double>(ix) + 0.5) * m_dDeltaElement);
}

///////////////////////////////////////////////////////////////////////////////

double GridSpacingUniform::GetEdge(int ix) const {
	return (m_dZeroCoord + static_cast<double>(ix) * m_dDeltaElement);
}

///////////////////////////////////////////////////////////////////////////////

double GridSpacingUniform::GetNodeNormArea(int ix) const {
	return (m_dDeltaElement);
}

///////////////////////////////////////////////////////////////////////////////

double GridSpacingUniform::GetEdgeNormArea(int ix) const {
	return (m_dDeltaElement);
}

///////////////////////////////////////////////////////////////////////////////
// GridSpacingGaussLobatto
///////////////////////////////////////////////////////////////////////////////

GridSpacingGaussLobatto::GridSpacingGaussLobatto(
	double dDeltaElement,
	double dZeroCoord,
	int nOrder
) :
	GridSpacing(dDeltaElement, dZeroCoord),
	m_nOrder(nOrder)
{
	// Check order
	if (m_nOrder < 2) {
		_EXCEPTION1("Invalid order of accuracy %i (Order < 2)", m_nOrder);
	}

	// Obtain GLL nodes
	GaussLobattoQuadrature::GetPoints(
		nOrder, 0.0, dDeltaElement, m_dG, m_dW);
}

///////////////////////////////////////////////////////////////////////////////

double GridSpacingGaussLobatto::GetNode(int ix) const {
	int ixElement    = (ix / (m_nOrder - 1));
	int ixSubElement = (ix % (m_nOrder - 1));

	// Case of a negative index
	if (ix < 0) {
		ixElement = ixElement - 1;
		ixSubElement = m_nOrder + ixSubElement;
	}

	if (ixSubElement >= m_nOrder) {
		_EXCEPTIONT("Logic error");
	}

	return (m_dZeroCoord
		+ m_dDeltaElement * static_cast<double>(ixElement)
		+ m_dG[ixSubElement]);
}

///////////////////////////////////////////////////////////////////////////////

double GridSpacingGaussLobatto::GetEdge(int ix) const {
	int ixElement    = (ix / (m_nOrder - 1));
	int ixSubElement = (ix % (m_nOrder - 1));

	if (ix < 0) {
		ixElement = ixElement - 1;
		ixSubElement = m_nOrder + ixSubElement;
	}

	if (ixSubElement >= m_nOrder) {
		_EXCEPTIONT("Logic error");
	}

	return (m_dZeroCoord
		+ m_dDeltaElement * static_cast<double>(ixElement)
		+ m_dG[ixSubElement]);
}

///////////////////////////////////////////////////////////////////////////////

double GridSpacingGaussLobatto::GetNodeNormArea(int ix) const {
	int ixSubElement = (ix % (m_nOrder - 1));

	if (ixSubElement == 0) {
		return (2.0 * m_dW[0]);
	} else {
		return (m_dW[ixSubElement]);
	}
}

///////////////////////////////////////////////////////////////////////////////

double GridSpacingGaussLobatto::GetEdgeNormArea(int ix) const {
	int ixSubElement = (ix % (m_nOrder - 1));

	if (ixSubElement == 0) {
		return (2.0 * m_dW[0]);
	} else {
		return (m_dW[ixSubElement]);
	}
}

///////////////////////////////////////////////////////////////////////////////
// GridSpacingGaussLobattoRepeated
///////////////////////////////////////////////////////////////////////////////

GridSpacingGaussLobattoRepeated::GridSpacingGaussLobattoRepeated(
	double dDeltaElement,
	double dZeroCoord,
	int nOrder
) :
	GridSpacing(dDeltaElement, dZeroCoord),
	m_nOrder(nOrder)
{
	// Check order
	if (m_nOrder < 2) {
		_EXCEPTION1("Invalid order of accuracy %i (Order < 2)", m_nOrder);
	}

	// Obtain GLL nodes
	GaussLobattoQuadrature::GetPoints(
		nOrder, 0.0, dDeltaElement, m_dG, m_dW);
}

///////////////////////////////////////////////////////////////////////////////

double GridSpacingGaussLobattoRepeated::GetNode(int ix) const {
	int ixElement    = (ix / m_nOrder);
	int ixSubElement = (ix % m_nOrder);

	// Case of a negative index
	if (ix < 0) {
		ixElement = ixElement - 1;
		ixSubElement = m_nOrder + ixSubElement;
	}

	if (ixSubElement >= m_nOrder) {
		_EXCEPTIONT("Logic error");
	}

	return (m_dZeroCoord
		+ m_dDeltaElement * static_cast<double>(ixElement)
		+ m_dG[ixSubElement]);
}

///////////////////////////////////////////////////////////////////////////////

double GridSpacingGaussLobattoRepeated::GetEdge(int ix) const {
	int ixElement    = (ix / m_nOrder);
	int ixSubElement = (ix % m_nOrder);

	if (ix < 0) {
		ixElement = ixElement - 1;
		ixSubElement = m_nOrder + ixSubElement;
	}

	if (ixSubElement >= m_nOrder) {
		_EXCEPTIONT("Logic error");
	}

	return (m_dZeroCoord
		+ m_dDeltaElement * static_cast<double>(ixElement)
		+ m_dG[ixSubElement]);
}

///////////////////////////////////////////////////////////////////////////////

double GridSpacingGaussLobattoRepeated::GetNodeNormArea(int ix) const {
	int ixSubElement = (ix % m_nOrder);

	return (m_dW[ixSubElement]);
}

///////////////////////////////////////////////////////////////////////////////

double GridSpacingGaussLobattoRepeated::GetEdgeNormArea(int ix) const {
	int ixSubElement = (ix % m_nOrder);

	return (m_dW[ixSubElement]);
}

///////////////////////////////////////////////////////////////////////////////
// GridSpacingMixedGaussLobatto
///////////////////////////////////////////////////////////////////////////////

GridSpacingMixedGaussLobatto::GridSpacingMixedGaussLobatto(
	double dDeltaElement,
	double dZeroCoord,
	int nOrder
) :
	GridSpacing(dDeltaElement, dZeroCoord),
	m_nOrder(nOrder)
{
	// Check order
	if (m_nOrder < 1) {
		_EXCEPTION1("Invalid order of accuracy %i (Order < 1)", m_nOrder);
	}

	// Obtain GLL nodes
	GaussLobattoQuadrature::GetPoints(
		nOrder+1, 0.0, dDeltaElement, m_dGL, m_dWL);

	// Obtain GL nodes
	GaussQuadrature::GetPoints(
		nOrder, 0.0, dDeltaElement, m_dG, m_dW);
}

///////////////////////////////////////////////////////////////////////////////

double GridSpacingMixedGaussLobatto::GetNode(int ix) const {
	int ixElement    = (ix / m_nOrder);
	int ixSubElement = (ix % m_nOrder);

	// Case of a negative index
	if (ix < 0) {
		ixElement = ixElement - 1;
		ixSubElement = m_nOrder + ixSubElement;
	}

	if (ixSubElement >= m_nOrder) {
		_EXCEPTIONT("Logic error");
	}

	return (m_dZeroCoord
		+ m_dDeltaElement * static_cast<double>(ixElement)
		+ m_dG[ixSubElement]);
}

///////////////////////////////////////////////////////////////////////////////

double GridSpacingMixedGaussLobatto::GetEdge(int ix) const {
	int ixElement    = (ix / m_nOrder);
	int ixSubElement = (ix % m_nOrder);

	// Case of a negative index
	if (ix < 0) {
		ixElement = ixElement - 1;
		ixSubElement = m_nOrder + ixSubElement;
	}

	if (ixSubElement >= m_nOrder) {
		_EXCEPTIONT("Logic error");
	}

	return (m_dZeroCoord
		+ m_dDeltaElement * static_cast<double>(ixElement)
		+ m_dGL[ixSubElement]);
}

///////////////////////////////////////////////////////////////////////////////

double GridSpacingMixedGaussLobatto::GetNodeNormArea(int ix) const {
	int ixSubElement = (ix % m_nOrder);

	return (m_dW[ixSubElement]);
}

///////////////////////////////////////////////////////////////////////////////

double GridSpacingMixedGaussLobatto::GetEdgeNormArea(int ix) const {
	int ixSubElement = (ix % m_nOrder);

	if (ixSubElement == 0) {
		return (2.0 * m_dWL[ixSubElement]);
	} else {
		return (m_dWL[ixSubElement]);
	}
}

///////////////////////////////////////////////////////////////////////////////

