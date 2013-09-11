///////////////////////////////////////////////////////////////////////////////
///
///	\file    GridPatchCSGLL.cpp
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

#include "GridPatchCSGLL.h"
#include "GridCSGLL.h"
#include "Model.h"
#include "TestCase.h"
#include "GridSpacing.h"

#include "Direction.h"
#include "CubedSphereTrans.h"
#include "PolynomialInterp.h"
#include "GaussQuadrature.h"
#include "GaussLobattoQuadrature.h"

#include "Announce.h"
#include "MathHelper.h"

#include <cmath>

///////////////////////////////////////////////////////////////////////////////

GridPatchCSGLL::GridPatchCSGLL(
	GridCSGLL & grid,
	int ixPatch,
	const PatchBox & box,
	int nOrder,
	int nVerticalOrder
) :
	GridPatch(grid, ixPatch, box),
	m_nHorizontalOrder(nOrder),
	m_nVerticalOrder(nVerticalOrder)
{
	// Verify that box boundaries are aligned with elements
	if (((box.GetAGlobalInteriorBegin() % nOrder) != 0) ||
		((box.GetAGlobalInteriorEnd()   % nOrder) != 0) ||
		((box.GetBGlobalInteriorBegin() % nOrder) != 0) ||
		((box.GetBGlobalInteriorEnd()   % nOrder) != 0)
	) {
		_EXCEPTION4(
			"CSGLL grid patch must be aligned with elements: %i %i %i %i",
			box.GetAGlobalInteriorBegin(),
			box.GetAGlobalInteriorEnd(),
			box.GetBGlobalInteriorBegin(),
			box.GetBGlobalInteriorEnd());
	}

	// Get panels in each coordinate direction
	if (grid.GetABaseResolution() != grid.GetBBaseResolution()) {
		_EXCEPTIONT("Invalid grid; CubedSphere grids must be square");
	}

	int ixDest;
	int jxDest;
	bool fSwitchAB;
	bool fSwitchPar;
	bool fSwitchPerp;

	m_ixNeighborPanel.resize(8);

	// Towards the right
	CubedSphereTrans::RelativeCoord(
		m_nHorizontalOrder * grid.GetABaseResolution(),
		box.GetPanel(),
		box.GetAGlobalInteriorEnd(),
		box.GetBGlobalInteriorBegin(),
		m_ixNeighborPanel[(int)(Direction_Right)],
		ixDest, jxDest,
		fSwitchAB, fSwitchPar, fSwitchPerp);

	// Towards the top
	CubedSphereTrans::RelativeCoord(
		m_nHorizontalOrder * grid.GetABaseResolution(),
		box.GetPanel(),
		box.GetAGlobalInteriorBegin(),
		box.GetBGlobalInteriorEnd(),
		m_ixNeighborPanel[(int)(Direction_Top)],
		ixDest, jxDest,
		fSwitchAB, fSwitchPar, fSwitchPerp);

	// Towards the left
	CubedSphereTrans::RelativeCoord(
		m_nHorizontalOrder * grid.GetABaseResolution(),
		box.GetPanel(),
		box.GetAGlobalInteriorBegin()-1,
		box.GetBGlobalInteriorBegin(),
		m_ixNeighborPanel[(int)(Direction_Left)],
		ixDest, jxDest,
		fSwitchAB, fSwitchPar, fSwitchPerp);

	// Towards the bottom
	CubedSphereTrans::RelativeCoord(
		m_nHorizontalOrder * grid.GetABaseResolution(),
		box.GetPanel(),
		box.GetAGlobalInteriorBegin(),
		box.GetBGlobalInteriorBegin()-1,
		m_ixNeighborPanel[(int)(Direction_Bottom)],
		ixDest, jxDest,
		fSwitchAB, fSwitchPar, fSwitchPerp);

	// Towards the top-right
	CubedSphereTrans::RelativeCoord(
		m_nHorizontalOrder * grid.GetABaseResolution(),
		box.GetPanel(),
		box.GetAGlobalInteriorEnd(),
		box.GetBGlobalInteriorEnd(),
		m_ixNeighborPanel[(int)(Direction_TopRight)],
		ixDest, jxDest,
		fSwitchAB, fSwitchPar, fSwitchPerp);

	// Towards the top-left
	CubedSphereTrans::RelativeCoord(
		m_nHorizontalOrder * grid.GetABaseResolution(),
		box.GetPanel(),
		box.GetAGlobalInteriorBegin()-1,
		box.GetBGlobalInteriorEnd(),
		m_ixNeighborPanel[(int)(Direction_TopLeft)],
		ixDest, jxDest,
		fSwitchAB, fSwitchPar, fSwitchPerp);

	// Towards the bottom-left
	CubedSphereTrans::RelativeCoord(
		m_nHorizontalOrder * grid.GetABaseResolution(),
		box.GetPanel(),
		box.GetAGlobalInteriorBegin()-1,
		box.GetBGlobalInteriorBegin()-1,
		m_ixNeighborPanel[(int)(Direction_BottomLeft)],
		ixDest, jxDest,
		fSwitchAB, fSwitchPar, fSwitchPerp);

	// Towards the bottom-right
	CubedSphereTrans::RelativeCoord(
		m_nHorizontalOrder * grid.GetABaseResolution(),
		box.GetPanel(),
		box.GetAGlobalInteriorEnd(),
		box.GetBGlobalInteriorBegin()-1,
		m_ixNeighborPanel[(int)(Direction_BottomRight)],
		ixDest, jxDest,
		fSwitchAB, fSwitchPar, fSwitchPerp);
}

///////////////////////////////////////////////////////////////////////////////

void GridPatchCSGLL::InitializeDataLocal() {

	// Allocate data
	GridPatch::InitializeDataLocal();

	// Initialize the longitude and latitude at each node
	for (int i = 0; i < m_box.GetATotalWidth(); i++) {
	for (int j = 0; j < m_box.GetBTotalWidth(); j++) {
		CubedSphereTrans::RLLFromABP(
			m_box.GetANode(i),
			m_box.GetBNode(j),
			m_box.GetPanel(),
			m_dataLon[i][j],
			m_dataLat[i][j]);
	}
	}
/*
	// Physical constants
	const PhysicalConstants & phys = m_grid.GetModel().GetPhysicalConstants();

	// Range of patch
	int ixPatchEnd = m_grid.GetABaseResolution()
		* IntPow(m_grid.GetRefinementRatio(), m_box.GetRefinementLevel());

	// Metric terms at each node; note that these metrics do not include
	// forcings due to topography.
	for (int k = 0; k < m_grid.GetRElements(); k++) {
	for (int i = 0; i < m_box.GetATotalWidth(); i++) {
	for (int j = 0; j < m_box.GetBTotalWidth(); j++) {
		double dX = tan(m_box.GetANode(i));
		double dY = tan(m_box.GetBNode(j));

		double dDelta2 = (1.0 + dX * dX + dY * dY);
		double dDelta  = sqrt(dDelta2);

		// Find index within sub-element
		int ix = (i - m_box.GetHaloElements()) % m_nHorizontalOrder;
		int jx = (j - m_box.GetHaloElements()) % m_nHorizontalOrder;

		if (ix < 0) {
			ix = ix + m_nHorizontalOrder;
		}
		if (jx < 0) {
			jx = jx + m_nHorizontalOrder;
		}

		int kx = k % m_nVerticalOrder;

		// Calculate pointwise Jacobian
		m_dataJacobian[k][i][j] =
			(1.0 + dX * dX) * (1.0 + dY * dY) / (dDelta * dDelta * dDelta);

		m_dataJacobian[k][i][j] *=
			  phys.GetEarthRadius()
			* phys.GetEarthRadius()
			* m_grid.GetZtop();

		// Element area associated with each GLL node
		m_dataElementArea[k][i][j] =
			m_dataJacobian[k][i][j]
			* dWL[ix] * dElementDeltaA
			* dWL[jx] * dElementDeltaB
			* dWNode[kx] * dElementDeltaXi;

		// Contravariant metric components at each node
		m_dataContraMetricA[k][i][j][0] = dDelta2 / (1.0 + dX * dX);
		m_dataContraMetricA[k][i][j][1] =
			dDelta2 * dX * dY / (1.0 + dX * dX) / (1.0 + dY * dY);
		m_dataContraMetricA[k][i][j][2] = 0.0;

		m_dataContraMetricB[k][i][j][0] =
			dDelta2 * dX * dY / (1.0 + dX * dX) / (1.0 + dY * dY);
		m_dataContraMetricB[k][i][j][1] = dDelta2 / (1.0 + dY * dY);
		m_dataContraMetricB[k][i][j][2] = 0.0;

		m_dataContraMetricA[k][i][j][0] /=
			phys.GetEarthRadius() * phys.GetEarthRadius();
		m_dataContraMetricA[k][i][j][1] /=
			phys.GetEarthRadius() * phys.GetEarthRadius();
		m_dataContraMetricA[k][i][j][2] /=
			phys.GetEarthRadius() * phys.GetEarthRadius();

		m_dataContraMetricB[k][i][j][0] /=
			phys.GetEarthRadius() * phys.GetEarthRadius();
		m_dataContraMetricB[k][i][j][1] /=
			phys.GetEarthRadius() * phys.GetEarthRadius();
		m_dataContraMetricB[k][i][j][2] /=
			phys.GetEarthRadius() * phys.GetEarthRadius();

		// Christoffel symbol components at each node
		// (off-diagonal element are doubled due to symmetry)
		m_dataChristoffelA[k][i][j][0] = 2.0 * dX * dY * dY / dDelta2;
		m_dataChristoffelA[k][i][j][1] = - 2.0 * dY * (1.0 + dY * dY) / dDelta2;
		m_dataChristoffelA[k][i][j][2] = 0.0;
		m_dataChristoffelB[k][i][j][0] = 0.0;
		m_dataChristoffelB[k][i][j][1] = - 2.0 * dX * (1.0 + dX * dX) / dDelta2;
		m_dataChristoffelB[k][i][j][2] = 2.0 * dX * dX * dY / dDelta2;
	}
	}
	}

	// Metric terms at each interface
	for (int k = 0; k <= m_grid.GetRElements(); k++) {
	for (int i = 0; i < m_box.GetATotalWidth(); i++) {
	for (int j = 0; j < m_box.GetBTotalWidth(); j++) {

		double dX = tan(m_box.GetANode(i));
		double dY = tan(m_box.GetBNode(j));

		double dDelta2 = (1.0 + dX * dX + dY * dY);
		double dDelta  = sqrt(dDelta2);

		// Find index within sub-element
		int ix = (i - m_box.GetHaloElements()) % m_nHorizontalOrder;
		int jx = (j - m_box.GetHaloElements()) % m_nHorizontalOrder;

		if (ix < 0) {
			ix = ix + m_nHorizontalOrder;
		}
		if (jx < 0) {
			jx = jx + m_nHorizontalOrder;
		}

		int kx = k % m_nVerticalOrder;

		// Calculate pointwise Jacobian
		double dJacobian =
			(1.0 + dX * dX) * (1.0 + dY * dY) / (dDelta * dDelta * dDelta);

		dJacobian *=
			  phys.GetEarthRadius()
			* phys.GetEarthRadius()
			* m_grid.GetZtop();

		// Element area associated with each vertical interface node
		m_dataElementAreaREdge[k][i][j] =
			dJacobian
			* dWL[ix] * dElementDeltaA
			* dWL[jx] * dElementDeltaB
			* dWREdge[kx] * dElementDeltaXi;

	}
	}
	}
*/
}

///////////////////////////////////////////////////////////////////////////////

void GridPatchCSGLL::EvaluateGeometricTerms(
	const TestCase & test
) {
	// Physical constants
	const PhysicalConstants & phys = m_grid.GetModel().GetPhysicalConstants();

	// Obtain Gauss Lobatto quadrature nodes and weights
	DataVector<double> dGL;
	DataVector<double> dWL;

	GaussLobattoQuadrature::GetPoints(m_nHorizontalOrder, 0.0, 1.0, dGL, dWL);

	// Obtain Gaussian quadrature nodes and weights in the vertical
	DataVector<double> dGNode;
	DataVector<double> dWNode;

	GaussQuadrature::GetPoints(
		m_nVerticalOrder, 0.0, 1.0, dGNode, dWNode);

	// Obtain Gauss Lobatto quadrature nodes and weights in the vertical
	DataVector<double> dGREdge;
	DataVector<double> dWREdge;

	GaussLobattoQuadrature::GetPoints(
		m_nVerticalOrder+1, 0.0, 1.0, dGREdge, dWREdge);

	// Element spacing at this refinement level
	double dElementDeltaA =
		  m_box.GetAEdge(m_box.GetHaloElements() + m_nHorizontalOrder)
		- m_box.GetAEdge(m_box.GetHaloElements());

	double dElementDeltaB =
		  m_box.GetBEdge(m_box.GetHaloElements() + m_nHorizontalOrder)
		- m_box.GetBEdge(m_box.GetHaloElements());

#pragma message "Implement rectangular grid elements"
	if (dElementDeltaA != dElementDeltaB) {
		_EXCEPTIONT("Not implemented.");
	}

	// Vertical grid spacing
	double dElementDeltaXi =
		  m_grid.GetREtaInterface(m_nVerticalOrder)
		- m_grid.GetREtaInterface(0);

	// Initialize metric and Christoffel symbols in terrain-following coords
	for (int i = 0; i < m_box.GetATotalWidth(); i++) {
	for (int j = 0; j < m_box.GetBTotalWidth(); j++) {

		// Find index within sub-element
		int ix = (i - m_box.GetHaloElements()) % m_nHorizontalOrder;
		int jx = (j - m_box.GetHaloElements()) % m_nHorizontalOrder;

		if (ix < 0) {
			ix = ix + m_nHorizontalOrder;
		}
		if (jx < 0) {
			jx = jx + m_nHorizontalOrder;
		}

		// Evaluate derivatives of topography
		double dEpsilon = 1.0e-5;
		double dTopography[3][3];
		for (int s = 0; s < 3; s++) {
		for (int t = 0; t < 3; t++) {
			double dAlpha =
				m_box.GetANode(i) + static_cast<double>(s-1) * dEpsilon;
			double dBeta =
				m_box.GetBNode(j) + static_cast<double>(t-1) * dEpsilon;

			double dLon;
			double dLat;
			CubedSphereTrans::RLLFromABP(
				dAlpha, dBeta, m_box.GetPanel(),
				dLon, dLat);

			dTopography[s][t] =
				test.EvaluateTopography(phys, dLon, dLat);
		}
		}

		// Gnomonic coordinates
		double dX = tan(m_box.GetANode(i));
		double dY = tan(m_box.GetBNode(j));
		double dDelta2 = (1.0 + dX * dX + dY * dY);
		double dDelta = sqrt(dDelta2);

		// Topography height and its derivatives
		double dZs = dTopography[1][1];
		double dDaZs =
			(dTopography[2][1] - dTopography[0][1]) / (2.0 * dEpsilon);
		double dDbZs =
			(dTopography[1][2] - dTopography[1][0]) / (2.0 * dEpsilon);
		double dDaaZs =
			(dTopography[2][1] - 2.0 * dZs + dTopography[0][1])
				/ (dEpsilon * dEpsilon);
		double dDbbZs =
			(dTopography[1][2] - 2.0 * dZs + dTopography[1][0])
				/ (dEpsilon * dEpsilon);
		double dDabZs =
			( dTopography[2][2] - dTopography[2][0]
			- dTopography[0][2] + dTopography[0][0])
				/ (dEpsilon * dEpsilon);

		// Vertical coordinate transform and its derivatives
		for (int k = 0; k < m_grid.GetRElements(); k++) {

			// Sub-element index
			int kx = k % m_nVerticalOrder;

			// Gal-Chen and Somerville (1975) linear terrain-following coord
			double dDaR = (1.0 - m_grid.GetREtaLevel(k)) * dDaZs;
			double dDbR = (1.0 - m_grid.GetREtaLevel(k)) * dDbZs;
			double dDaaR = (1.0 - m_grid.GetREtaLevel(k)) * dDaaZs;
			double dDabR = (1.0 - m_grid.GetREtaLevel(k)) * dDabZs;
			double dDbbR = (1.0 - m_grid.GetREtaLevel(k)) * dDbbZs;

			double dDxR = m_grid.GetZtop() - dZs;
			double dDaxR = - dDaZs;
			double dDbxR = - dDbZs;
			double dDxxR = 0.0;

			// Calculate pointwise Jacobian
			m_dataJacobian[k][i][j] =
				(1.0 + dX * dX) * (1.0 + dY * dY) / (dDelta * dDelta * dDelta);

#pragma message "Why does the Jacobian not include the DxR term? -- its inclusion shouldn't mess up the Coriolis force"
			m_dataJacobian[k][i][j] *=
				  phys.GetEarthRadius()
				* phys.GetEarthRadius();

			// Element area associated with each model level GLL node
			m_dataElementArea[k][i][j] =
				m_dataJacobian[k][i][j]
				* dWL[ix] * dElementDeltaA
				* dWL[jx] * dElementDeltaB
				* dWNode[kx] * dElementDeltaXi * dDxR;

			// Contravariant metric components
			double dContraMetricScale = 
				dDelta2 / (1.0 + dX * dX) / (1.0 + dY * dY)
				/ (phys.GetEarthRadius() * phys.GetEarthRadius());

			m_dataContraMetricA[k][i][j][0] =
				dContraMetricScale * (1.0 + dY * dY);
			m_dataContraMetricA[k][i][j][1] =
				dContraMetricScale * dX * dY;
			m_dataContraMetricA[k][i][j][2] = 0.0;

			m_dataContraMetricB[k][i][j][0] =
				dContraMetricScale * dX * dY;
			m_dataContraMetricB[k][i][j][1] =
				dContraMetricScale * (1.0 + dX * dX);
			m_dataContraMetricB[k][i][j][2] = 0.0;

			m_dataContraMetricXi[k][i][j][0] = 
				- dContraMetricScale / dDxR
					* ((1.0 + dY * dY) * dDaR + dX * dY * dDbR);
			m_dataContraMetricXi[k][i][j][1] =
				- dContraMetricScale / dDxR
					* (dX * dY * dDaR + (1.0 + dX * dX) * dDbR);
			m_dataContraMetricXi[k][i][j][2] =
				1.0 / (dDxR * dDxR)
					* (1.0 + dContraMetricScale
						* (  (1.0 + dY * dY) * dDaR * dDaR
						   + 2.0 * dX * dY * dDaR * dDbR
						   + (1.0 + dX * dX) * dDbR * dDbR));

			// Christoffel symbol components at each node
			// (off-diagonal element are doubled due to symmetry)
			m_dataChristoffelA[k][i][j][0] =
				2.0 * dX * dY * dY / dDelta2;
			m_dataChristoffelA[k][i][j][1] =
				- 2.0 * dY * (1.0 + dY * dY) / dDelta2;
			m_dataChristoffelA[k][i][j][2] = 0.0;
			m_dataChristoffelB[k][i][j][0] = 0.0;
			m_dataChristoffelB[k][i][j][1] =
				- 2.0 * dX * (1.0 + dX * dX) / dDelta2;
			m_dataChristoffelB[k][i][j][2] =
				2.0 * dX * dX * dY / dDelta2;

			// Vertical Christoffel symbol components
			m_dataChristoffelXi[k][i][j][0] =
				- 2.0 * dX * dY * dY / dDelta2 * dDaR + dDaaR;

			m_dataChristoffelXi[k][i][j][1] =
				+ 2.0 * dY * (1.0 + dY * dY) / dDelta2 * dDaR
				+ 2.0 * dX * (1.0 + dX * dX) / dDelta2 * dDbR
				+ 2.0 * dDabR;

			m_dataChristoffelXi[k][i][j][2] = 2.0 * dDaxR;

			m_dataChristoffelXi[k][i][j][3] =
				- 2.0 * dX * dX * dY / dDelta2 * dDbR + dDbbR;

			m_dataChristoffelXi[k][i][j][4] = 2.0 * dDbxR;

			m_dataChristoffelXi[k][i][j][5] = dDxxR;

			m_dataChristoffelXi[k][i][j][0] /= dDxR;
			m_dataChristoffelXi[k][i][j][1] /= dDxR;
			m_dataChristoffelXi[k][i][j][2] /= dDxR;
			m_dataChristoffelXi[k][i][j][3] /= dDxR;
			m_dataChristoffelXi[k][i][j][4] /= dDxR;
			m_dataChristoffelXi[k][i][j][5] /= dDxR;
		}

		// Metric terms at vertical interfaces
		for (int k = 0; k <= m_grid.GetRElements(); k++) {
			int kx = k % m_nVerticalOrder;

			// Gal-Chen and Somerville (1975) linear terrain-following coord
			double dDxR = m_grid.GetZtop() - dZs;

			// Calculate pointwise Jacobian
			double dJacobian =
				(1.0 + dX * dX) * (1.0 + dY * dY) / (dDelta * dDelta * dDelta);

			dJacobian *=
				  phys.GetEarthRadius()
				* phys.GetEarthRadius();

			// Element area associated with each model interface GLL node
			m_dataElementAreaREdge[k][i][j] =
				dJacobian
				* dWL[ix] * dElementDeltaA
				* dWL[jx] * dElementDeltaB
				* dWREdge[kx] * dElementDeltaXi * dDxR;

			if ((k != 0) && (k != m_grid.GetRElements()) && (kx == 0)) {
				m_dataElementAreaREdge[k][i][j] *= 2.0;
			}
		}
	}
	}
}

///////////////////////////////////////////////////////////////////////////////

void GridPatchCSGLL::EvaluateTestCase(
	const TestCase & test,
	double dTime,
	int iDataIndex
) {
	// Initialize the data at each node
	if (m_datavecStateNode.size() == 0) {
		_EXCEPTIONT("InitializeData must be called before InitialConditions");
	}
	if (iDataIndex >= m_datavecStateNode.size()) {
		_EXCEPTIONT("Invalid iDataIndex (out of range)");
	}

	// Check dimensionality
	if ((m_grid.GetModel().GetEquationSet().GetDimensionality() == 2) &&
		(m_nVerticalOrder != 1)
	) {
		_EXCEPTIONT("VerticalOrder / Dimensionality mismatch:\n"
			"For 2D problems vertical order must be 1.");
	}

	// Evaluate geometric terms
	EvaluateGeometricTerms(test);

	// Physical constants
	const PhysicalConstants & phys = m_grid.GetModel().GetPhysicalConstants();

	// Initialize the topography at each node
	for (int i = 0; i < m_box.GetATotalWidth(); i++) {
	for (int j = 0; j < m_box.GetBTotalWidth(); j++) {
		m_dataTopography[i][j] =
			test.EvaluateTopography(
				phys,
				m_dataLon[i][j],
				m_dataLat[i][j]);

		if (m_dataTopography[i][j] >= m_grid.GetZtop()) {
			_EXCEPTIONT("TestCase topography exceeds model top.");
		}

		// Gal-Chen and Sommerville vertical coordinate
		for (int k = 0; k < m_grid.GetRElements(); k++) {
			m_dataZLevels[k][i][j] =
				m_dataTopography[i][j]
					+ m_grid.GetREtaLevel(k)
						* (m_grid.GetZtop() - m_dataTopography[i][j]);
		}
		for (int k = 0; k <= m_grid.GetRElements(); k++) {
			m_dataZInterfaces[k][i][j] =
				m_dataTopography[i][j]
					+ m_grid.GetREtaInterface(k)
						* (m_grid.GetZtop() - m_dataTopography[i][j]);
		}
	}
	}

	// Buffer vector for storing pointwise states
	int nComponents = m_grid.GetModel().GetEquationSet().GetComponents();
	int nTracers = m_grid.GetModel().GetEquationSet().GetTracers();

	DataVector<double> dPointwiseState;
	dPointwiseState.Initialize(nComponents);

	DataVector<double> dPointwiseRefState;
	dPointwiseRefState.Initialize(nComponents);

	DataVector<double> dPointwiseTracers;
	if (m_datavecTracers.size() > 0) {
		dPointwiseTracers.Initialize(nTracers);
	}

	// Evaluate the state on model levels
	for (int k = 0; k < m_grid.GetRElements(); k++) {
	for (int i = 0; i < m_box.GetATotalWidth(); i++) {
	for (int j = 0; j < m_box.GetBTotalWidth(); j++) {

		// Evaluate pointwise state
		test.EvaluatePointwiseState(
			m_grid.GetModel().GetPhysicalConstants(),
			dTime,
			m_dataZLevels[k][i][j],
			m_dataLon[i][j],
			m_dataLat[i][j],
			dPointwiseState,
			dPointwiseTracers);

		for (int c = 0; c < dPointwiseState.GetRows(); c++) {
			m_datavecStateNode[iDataIndex][c][k][i][j] = dPointwiseState[c];
		}

		// Transform state velocities
		double dUlon;
		double dUlat;

		dUlon = m_datavecStateNode[iDataIndex][0][k][i][j];
		dUlat = m_datavecStateNode[iDataIndex][1][k][i][j];

		dUlon /= phys.GetEarthRadius();
		dUlat /= phys.GetEarthRadius();

		CubedSphereTrans::VecTransABPFromRLL(
			tan(m_box.GetANode(i)),
			tan(m_box.GetBNode(j)),
			m_box.GetPanel(),
			dUlon, dUlat,
			m_datavecStateNode[iDataIndex][0][k][i][j],
			m_datavecStateNode[iDataIndex][1][k][i][j]);

		// Evaluate reference state
		if (m_grid.HasReferenceState()) {
			test.EvaluateReferenceState(
				m_grid.GetModel().GetPhysicalConstants(),
				m_dataZLevels[k][i][j],
				m_dataLon[i][j],
				m_dataLat[i][j],
				dPointwiseRefState);

			for (int c = 0; c < dPointwiseState.GetRows(); c++) {
				m_dataRefStateNode[c][k][i][j] = dPointwiseRefState[c];
			}

			// Transform reference velocities
			dUlon = m_dataRefStateNode[0][k][i][j];
			dUlat = m_dataRefStateNode[1][k][i][j];

			dUlon /= phys.GetEarthRadius();
			dUlat /= phys.GetEarthRadius();

			CubedSphereTrans::VecTransABPFromRLL(
				tan(m_box.GetANode(i)),
				tan(m_box.GetBNode(j)),
				m_box.GetPanel(),
				dUlon, dUlat,
				m_dataRefStateNode[0][k][i][j],
				m_dataRefStateNode[1][k][i][j]);
		}

		// Evaluate tracers
		for (int c = 0; c < dPointwiseTracers.GetRows(); c++) {
			m_datavecTracers[iDataIndex][c][k][i][j] = dPointwiseTracers[c];
		}
	}
	}
	}

	// Evaluate the state on model interfaces
	for (int k = 0; k <= m_grid.GetRElements(); k++) {
	for (int i = 0; i < m_box.GetATotalWidth(); i++) {
	for (int j = 0; j < m_box.GetBTotalWidth(); j++) {

		// Evaluate pointwise state
		test.EvaluatePointwiseState(
			m_grid.GetModel().GetPhysicalConstants(),
			dTime,
			m_dataZInterfaces[k][i][j],
			m_dataLon[i][j],
			m_dataLat[i][j],
			dPointwiseState,
			dPointwiseTracers);

		for (int c = 0; c < dPointwiseState.GetRows(); c++) {
			m_datavecStateREdge[iDataIndex][c][k][i][j] = dPointwiseState[c];
		}

		// Transform state velocities
		double dUlon;
		double dUlat;

		dUlon = m_datavecStateREdge[iDataIndex][0][k][i][j];
		dUlat = m_datavecStateREdge[iDataIndex][1][k][i][j];

		dUlon /= phys.GetEarthRadius();
		dUlat /= phys.GetEarthRadius();

		CubedSphereTrans::VecTransABPFromRLL(
			tan(m_box.GetANode(i)),
			tan(m_box.GetBNode(j)),
			m_box.GetPanel(),
			dUlon, dUlat,
			m_datavecStateREdge[iDataIndex][0][k][i][j],
			m_datavecStateREdge[iDataIndex][1][k][i][j]);

		if (m_grid.HasReferenceState()) {
			test.EvaluateReferenceState(
				m_grid.GetModel().GetPhysicalConstants(),
				m_dataZInterfaces[k][i][j],
				m_dataLon[i][j],
				m_dataLat[i][j],
				dPointwiseRefState);

			for (int c = 0; c < dPointwiseState.GetRows(); c++) {
				m_dataRefStateREdge[c][k][i][j] = dPointwiseRefState[c];
			}

			// Transform reference velocities
			dUlon = m_dataRefStateREdge[0][k][i][j];
			dUlat = m_dataRefStateREdge[1][k][i][j];

			dUlon /= phys.GetEarthRadius();
			dUlat /= phys.GetEarthRadius();

			CubedSphereTrans::VecTransABPFromRLL(
				tan(m_box.GetANode(i)),
				tan(m_box.GetBNode(j)),
				m_box.GetPanel(),
				dUlon, dUlat,
				m_dataRefStateREdge[0][k][i][j],
				m_dataRefStateREdge[1][k][i][j]);
		}
	}
	}
	}
}

///////////////////////////////////////////////////////////////////////////////

void GridPatchCSGLL::ComputeCurlAndDiv(
	const GridData3D & dataUa,
	const GridData3D & dataUb
) const {
	// Parent grid
	const GridCSGLL & gridCSGLL = dynamic_cast<const GridCSGLL &>(m_grid);

	// Compute derivatives of the field
	const DataMatrix<double> & dDxBasis1D = gridCSGLL.GetDxBasis1D();

	// Element spacing
	double dElementDeltaA =
		  m_box.GetAEdge(m_box.GetHaloElements() + m_nHorizontalOrder)
		- m_box.GetAEdge(m_box.GetHaloElements());

	// Number of finite elements in each direction
	int nAFiniteElements = m_box.GetAInteriorWidth() / m_nHorizontalOrder;
	int nBFiniteElements = m_box.GetBInteriorWidth() / m_nHorizontalOrder;

	// Loop over all elements in the box
	for (int k = 0; k < gridCSGLL.GetRElements(); k++) {
	for (int a = 0; a < nAFiniteElements; a++) {
	for (int b = 0; b < nBFiniteElements; b++) {

		// Index of lower-left corner node
		int iA = a * m_nHorizontalOrder + m_box.GetHaloElements();
		int iB = b * m_nHorizontalOrder + m_box.GetHaloElements();

		for (int i = 0; i < m_nHorizontalOrder; i++) {
		for (int j = 0; j < m_nHorizontalOrder; j++) {

			// Pointwise field values
			double dUa = dataUa[k][iA+i][iB+j];
			double dUb = dataUb[k][iA+i][iB+j];

			// Compute derivatives at each node
			double dDaUa = 0.0;
			double dDaUb = 0.0;
			double dDbUa = 0.0;
			double dDbUb = 0.0;

			for (int s = 0; s < m_nHorizontalOrder; s++) {
				dDaUa += dataUa[k][iA+s][iB+j] * dDxBasis1D[s][i];
				dDaUb += dataUb[k][iA+s][iB+j] * dDxBasis1D[s][i];
				dDbUa += dataUa[k][iA+i][iB+s] * dDxBasis1D[s][j];
				dDbUb += dataUb[k][iA+i][iB+s] * dDxBasis1D[s][j];
			}

			dDaUa /= dElementDeltaA;
			dDaUb /= dElementDeltaA;
			dDbUa /= dElementDeltaA;
			dDbUb /= dElementDeltaA;

			// Compute covariant derivatives at node
			double dCovDaUa = dDaUa
				+ m_dataChristoffelA[k][iA+i][iB+j][0] * dUa
				+ m_dataChristoffelA[k][iA+i][iB+j][1] * 0.5 * dUb;

			double dCovDaUb = dDaUb
				+ m_dataChristoffelB[k][iA+i][iB+j][0] * dUa
				+ m_dataChristoffelB[k][iA+i][iB+j][1] * 0.5 * dUb;

			double dCovDbUa = dDbUa
				+ m_dataChristoffelA[k][iA+i][iB+j][1] * 0.5 * dUa
				+ m_dataChristoffelA[k][iA+i][iB+j][2] * dUb;

			double dCovDbUb = dDbUb
				+ m_dataChristoffelB[k][iA+i][iB+j][1] * 0.5 * dUa
				+ m_dataChristoffelB[k][iA+i][iB+j][2] * dUb;

			// Compute curl at node
			m_dataVorticity[k][iA+i][iB+j] = m_dataJacobian[k][iA+i][iB+j] * (
				+ m_dataContraMetricA[k][iA+i][iB+j][0] * dCovDaUb
				+ m_dataContraMetricA[k][iA+i][iB+j][1] * dCovDbUb
				- m_dataContraMetricB[k][iA+i][iB+j][0] * dCovDaUa
				- m_dataContraMetricB[k][iA+i][iB+j][1] * dCovDbUa);

			// Compute the divergence at node
			m_dataDivergence[k][iA+i][iB+j] = dCovDaUa + dCovDbUb;
		}
		}
	}
	}
	}
}

///////////////////////////////////////////////////////////////////////////////

void GridPatchCSGLL::ComputeVorticityDivergence(
	int iDataIndex
) {
	// Physical constants
	const PhysicalConstants & phys = m_grid.GetModel().GetPhysicalConstants();

	// Working data
	const GridData4D & dataState = GetDataState(iDataIndex, DataLocation_Node);

	if (dataState.GetComponents() < 2) {
		_EXCEPTIONT(
			"Insufficient components for vorticity calculation");
	}

	// Get the alpha and beta components of vorticity
	GridData3D dataUa;
	GridData3D dataUb;

	dataState.GetAsGridData3D(0, dataUa);
	dataState.GetAsGridData3D(1, dataUb);

	// Compute the radial component of the curl of the velocity field
	ComputeCurlAndDiv(dataUa, dataUb);

/*
	// Lagrangian differentiation coefficients element [0,1]
	DataVector<double> dG;
	DataVector<double> dW;

	GaussLobattoQuadrature::GetPoints(m_nHorizontalOrder, dG, dW);

	DataMatrix<double> dDiffCoeff;
	dDiffCoeff.Initialize(m_nHorizontalOrder, m_nHorizontalOrder);

	for (int i = 0; i < m_nHorizontalOrder; i++) {
		PolynomialInterp::DiffLagrangianPolynomialCoeffs(
			m_nHorizontalOrder, dG, dDiffCoeff[i], dG[i]);
	}

	// Number of finite elements in each direction
	int nAFiniteElements = m_box.GetAInteriorWidth() / m_nHorizontalOrder;
	int nBFiniteElements = m_box.GetBInteriorWidth() / m_nHorizontalOrder;

	// Loop over all elements in the box
	for (int k = 0; k < dataState.GetRElements(); k++) {
	for (int a = 0; a < nAFiniteElements; a++) {
	for (int b = 0; b < nBFiniteElements; b++) {
	
		for (int i = 0; i < m_nHorizontalOrder; i++) {
		for (int j = 0; j < m_nHorizontalOrder; j++) {

			// Sub-element index
			int iElementA = m_box.GetAInteriorBegin() + a * m_nHorizontalOrder;
			int iElementB = m_box.GetBInteriorBegin() + b * m_nHorizontalOrder;

			int iA = iElementA + i;
			int iB = iElementB + j;

			// Calculate pointwise derivatives
			double dDaUa = 0.0;
			double dDbUa = 0.0;
			double dDaUb = 0.0;
			double dDbUb = 0.0;

			for (int s = 0; s < m_nHorizontalOrder; s++) {
				double dPtX = tan(m_box.GetAEdge(iElementA + s));
				double dPtY = tan(m_box.GetBEdge(iB));

				double dPtDelta2 = (1.0 + dPtX * dPtX + dPtY * dPtY);

				// Conversion factors to the unit basis
				double dStretchA =
					phys.GetEarthRadius()
					* (1.0 + dPtX * dPtX)
					* sqrt(1.0 + dPtY * dPtY)
					/ dPtDelta2;

				double dStretchB =
					phys.GetEarthRadius()
					* sqrt(1.0 + dPtX * dPtX)
					* (1.0 + dPtY * dPtY)
					/ dPtDelta2;

				dDaUa +=
					  dDiffCoeff[i][s]
					* dataState[0][k][iElementA + s][iB]
					* dStretchA;

				dDaUb +=
					  dDiffCoeff[i][s]
					* dataState[1][k][iElementA + s][iB]
					* dStretchB;
			}

			for (int s = 0; s < m_nHorizontalOrder; s++) {
				double dPtX = tan(m_box.GetAEdge(iA));
				double dPtY = tan(m_box.GetBEdge(iElementB + s));

				double dPtDelta2 = (1.0 + dPtX * dPtX + dPtY * dPtY);

				// Conversion factors to the unit basis
				double dStretchA =
					phys.GetEarthRadius()
					* (1.0 + dPtX * dPtX)
					* sqrt(1.0 + dPtY * dPtY)
					/ dPtDelta2;

				double dStretchB =
					phys.GetEarthRadius()
					* sqrt(1.0 + dPtX * dPtX)
					* (1.0 + dPtY * dPtY)
					/ dPtDelta2;

				dDbUa +=
					  dDiffCoeff[j][s]
					* dataState[0][k][iA][iElementB + s]
					* dStretchA;

				dDbUb +=
					  dDiffCoeff[j][s]
					* dataState[1][k][iA][iElementB + s]
					* dStretchB;
			}

			// Compute vorticity
			double dAlpha = m_box.GetANode(iA);
			double dBeta  = m_box.GetBNode(iB);

			double dX = tan(dAlpha);
			double dY = tan(dBeta);

			double dDelta = sqrt(1.0 + dX * dX + dY * dY);
			double dC = 1.0 / sqrt(1.0 + dX * dX);
			double dD = 1.0 / sqrt(1.0 + dY * dY);

			dataVorticity[k][iA][iB] = dDelta / phys.GetEarthRadius() * (
				dX * dY * dC * dD * (dD * dDbUb - dC * dDaUa)
				- dD * dDbUa + dC * dDaUb);
		}
		}
	}
	}
	}
*/
}

///////////////////////////////////////////////////////////////////////////////

void GridPatchCSGLL::InterpolateNodeToREdge(
	int iVar,
	int iDataIndex
) {

	// Working data
	GridData4D & dataNode  = GetDataState(iDataIndex, DataLocation_Node);
	GridData4D & dataREdge = GetDataState(iDataIndex, DataLocation_REdge);

	// Parent grid, containing the vertical remapping information
	GridCSGLL * pCSGLLGrid = dynamic_cast<GridCSGLL*>(&m_grid);
	if (pCSGLLGrid == NULL) {
		_EXCEPTIONT("Logic error");
	}

	const DataMatrix<double> & dInterp = pCSGLLGrid->GetInterpNodeToREdge();

	int nVerticalOrder = pCSGLLGrid->GetVerticalOrder();
	if (dataNode.GetRElements() % nVerticalOrder != 0) {
		_EXCEPTIONT("Logic error");
	}

	// Loop over all elements in the box
	for (int i = m_box.GetAInteriorBegin(); i < m_box.GetAInteriorEnd(); i++) {
	for (int j = m_box.GetBInteriorBegin(); j < m_box.GetBInteriorEnd(); j++) {
		int nFiniteElements = dataNode.GetRElements() / nVerticalOrder;

		for (int k = 0; k < dataREdge.GetRElements(); k++) {
			dataREdge[iVar][k][i][j] = 0.0;
		}

		// Loop over all nodes
		for (int k = 0; k < dataNode.GetRElements(); k++) {
			int a = k / nVerticalOrder;
			int m = k % nVerticalOrder;
			int lBegin = a * nVerticalOrder;

			// Apply node value to interface
			for (int l = 0; l <= nVerticalOrder; l++) {
				dataREdge[iVar][lBegin + l][i][j] +=
					dInterp[l][m] * dataNode[iVar][k][i][j];
			}
		}

		// Halve interior element interface values
		for (int a = 1; a < nFiniteElements; a++) {
			dataREdge[iVar][a * nVerticalOrder][i][j] *= 0.5;
		}
	}
	}
}

///////////////////////////////////////////////////////////////////////////////

void GridPatchCSGLL::InterpolateREdgeToNode(
	int iVar,
	int iDataIndex
) {

	// Working data
	GridData4D & dataREdge = GetDataState(iDataIndex, DataLocation_REdge);
	GridData4D & dataNode  = GetDataState(iDataIndex, DataLocation_Node);

	// Parent grid, containing the vertical remapping information
	GridCSGLL * pCSGLLGrid = dynamic_cast<GridCSGLL*>(&m_grid);
	if (pCSGLLGrid == NULL) {
		_EXCEPTIONT("Logic error");
	}

	const DataMatrix<double> & dInterp = pCSGLLGrid->GetInterpREdgeToNode();

	int nVerticalOrder = pCSGLLGrid->GetVerticalOrder();
	if (dataNode.GetRElements() % nVerticalOrder != 0) {
		_EXCEPTIONT("Logic error");
	}

	// Loop over all elements in the box
	for (int i = m_box.GetAInteriorBegin(); i < m_box.GetAInteriorEnd(); i++) {
	for (int j = m_box.GetBInteriorBegin(); j < m_box.GetBInteriorEnd(); j++) {

		for (int k = 0; k < dataNode.GetRElements(); k++) {
			dataNode[iVar][k][i][j] = 0.0;
		}

		// Loop over all nodes
		for (int k = 0; k < dataNode.GetRElements(); k++) {
			int a = k / nVerticalOrder;
			int m = k % nVerticalOrder;
			int lBegin = a * nVerticalOrder;

			// Apply interface values to nodes
			for (int l = 0; l <= nVerticalOrder; l++) {
				dataNode[iVar][k][i][j] +=
					dInterp[m][l] * dataREdge[iVar][lBegin + l][i][j];
			}
		}
	}
	}
}

///////////////////////////////////////////////////////////////////////////////

void GridPatchCSGLL::InterpolateData(
	const DataVector<double> & dAlpha,
	const DataVector<double> & dBeta,
	const DataVector<int> & iPanel,
	DataType eDataType,
	DataMatrix3D<double> & dInterpData,
	bool fIncludeReferenceState,
	bool fConvertToPrimitive
) {
	if (dAlpha.GetRows() != dBeta.GetRows()) {
		_EXCEPTIONT("Point vectors must have equivalent length.");
	}

	// Vector for storage interpolated points
	DataVector<double> dAInterpCoeffs;
	dAInterpCoeffs.Initialize(m_nHorizontalOrder);

	DataVector<double> dBInterpCoeffs;
	dBInterpCoeffs.Initialize(m_nHorizontalOrder);

	DataVector<double> dADiffCoeffs;
	dADiffCoeffs.Initialize(m_nHorizontalOrder);

	DataVector<double> dBDiffCoeffs;
	dBDiffCoeffs.Initialize(m_nHorizontalOrder);

	DataVector<double> dAInterpPt;
	dAInterpPt.Initialize(m_nHorizontalOrder);

	// Element-wise grid spacing
	double dDeltaA =
		+ m_box.GetAEdge(m_box.GetHaloElements() + m_nHorizontalOrder)
		- m_box.GetAEdge(m_box.GetHaloElements());

	double dDeltaB =
		+ m_box.GetBEdge(m_box.GetHaloElements() + m_nHorizontalOrder)
		- m_box.GetBEdge(m_box.GetHaloElements());

	// Physical constants
	const PhysicalConstants & phys = m_grid.GetModel().GetPhysicalConstants();

	// Loop throught all points
	for (int i = 0; i < dAlpha.GetRows(); i++) {

		// Element index
		if ((iPanel[i] != m_box.GetPanel()) ||
			(dAlpha[i] < m_box.GetAEdge(m_box.GetAInteriorBegin())) ||
			(dAlpha[i] > m_box.GetAEdge(m_box.GetAInteriorEnd())) ||
			(dBeta[i] < m_box.GetBEdge(m_box.GetBInteriorBegin())) ||
			(dBeta[i] > m_box.GetBEdge(m_box.GetBInteriorEnd()))
		) {
			continue;
		}

		int iA =
			(dAlpha[i] - m_box.GetAEdge(m_box.GetAInteriorBegin())) / dDeltaA;

		int iB =
			(dBeta[i] - m_box.GetBEdge(m_box.GetBInteriorBegin())) / dDeltaB;

		// Bound the index within the element
		if (iA < 0) {
			iA = 0;
		}
		if (iA >= (m_box.GetAInteriorWidth() / m_nHorizontalOrder)) {
			iA = m_box.GetAInteriorWidth() / m_nHorizontalOrder - 1;
		}
		if (iB < 0) {
			iB = 0;
		}
		if (iB >= (m_box.GetBInteriorWidth() / m_nHorizontalOrder)) {
			iB = m_box.GetBInteriorWidth() / m_nHorizontalOrder - 1;
		}

		iA = m_box.GetHaloElements() + iA * m_nHorizontalOrder;
		iB = m_box.GetHaloElements() + iB * m_nHorizontalOrder;

		// Compute interpolation coefficients
		PolynomialInterp::LagrangianPolynomialCoeffs(
			m_nHorizontalOrder,
			&(m_box.GetAEdges()[iA]),
			dAInterpCoeffs,
			dAlpha[i]);

		PolynomialInterp::LagrangianPolynomialCoeffs(
			m_nHorizontalOrder,
			&(m_box.GetBEdges()[iB]),
			dBInterpCoeffs,
			dBeta[i]);

		// Perform interpolation on all variables
		int nComponents;
		if (eDataType == DataType_State) {
			nComponents = m_datavecStateNode[0].GetComponents();
		} else if (eDataType == DataType_Tracers) {
			nComponents = m_datavecTracers[0].GetComponents();
		} else if (eDataType == DataType_Vorticity) {
			nComponents = 1;
		} else if (eDataType == DataType_Divergence) {
			nComponents = 1;
		} else {
			_EXCEPTIONT("Invalid DataType");
		}

		for (int c = 0; c < nComponents; c++) {

			const double *** pData;
			if (eDataType == DataType_State) {
				pData = (const double ***)(m_datavecStateNode[0][c]);
			} else if (eDataType == DataType_Tracers) {
				pData = (const double ***)(m_datavecTracers[0][c]);
			} else if (eDataType == DataType_Vorticity) {
				pData = (const double ***)(double ***)(m_dataVorticity);
			} else if (eDataType == DataType_Divergence) {
				pData = (const double ***)(double ***)(m_dataDivergence);
			}

			// Perform interpolation on all levels
			for (int k = 0; k < m_grid.GetRElements(); k++) {

				dInterpData[c][k][i] = 0.0;

				for (int m = 0; m < m_nHorizontalOrder; m++) {
				for (int n = 0; n < m_nHorizontalOrder; n++) {
					dInterpData[c][k][i] +=
						  dAInterpCoeffs[m]
						* dBInterpCoeffs[n]
						* pData[k][iA+m][iB+n];
				}
				}

				// Do not include the reference state
				if ((eDataType == DataType_State) &&
					(!fIncludeReferenceState)
				) {
					for (int m = 0; m < m_nHorizontalOrder; m++) {
					for (int n = 0; n < m_nHorizontalOrder; n++) {
						dInterpData[c][k][i] -=
							  dAInterpCoeffs[m]
							* dBInterpCoeffs[n]
							* m_dataRefStateNode[c][k][iA+m][iB+n];
					}
					}
				}
			}
		}

		// Convert to primitive variables
		if ((eDataType == DataType_State) && (fConvertToPrimitive)) {
			for (int k = 0; k < m_grid.GetRElements(); k++) {
				double dUalpha = phys.GetEarthRadius()
					* dInterpData[0][k][i];
				double dUbeta = phys.GetEarthRadius()
					* dInterpData[1][k][i];

				CubedSphereTrans::VecTransRLLFromABP(
					tan(dAlpha[i]),
					tan(dBeta[i]),
					iPanel[i],
					dUalpha,
					dUbeta,
					dInterpData[0][k][i],
					dInterpData[1][k][i]);
			}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

