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

#include "Defines.h"
#include "GridPatchCSGLL.h"
#include "GridCSGLL.h"
#include "Model.h"
#include "TestCase.h"
#include "GridSpacing.h"
#include "VerticalStretch.h"

#include "Direction.h"
#include "CubedSphereTrans.h"
#include "PolynomialInterp.h"
#include "GaussLobattoQuadrature.h"

#include "Announce.h"
#include "MathHelper.h"

#include <cmath>

///////////////////////////////////////////////////////////////////////////////

GridPatchCSGLL::GridPatchCSGLL(
	GridCSGLL & grid,
	int ixPatch,
	const PatchBox & box,
	int nHorizontalOrder,
	int nVerticalOrder
) :
	GridPatchGLL(
		grid,
		ixPatch,
		box,
		nHorizontalOrder,
		nVerticalOrder)
{

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

	// Allocate buffer space for ComputeCurlAndDiv()
	m_dBufferConU.Allocate(
		2,
		m_nHorizontalOrder,
		m_nHorizontalOrder,
		grid.GetRElements());
}

///////////////////////////////////////////////////////////////////////////////

void GridPatchCSGLL::InitializeDataLocal(
	bool fAllocateGeometric,
	bool fAllocateActiveState,
	bool fAllocateBufferState,
	bool fAllocateAuxiliary
) {
	// Allocate data
	GridPatch::InitializeDataLocal(
		fAllocateGeometric,
		fAllocateActiveState,
		fAllocateBufferState,
		fAllocateAuxiliary
	);

	// Initialize coordinate data
	if (fAllocateGeometric) {
		InitializeCoordinateData();
	} else {
		//Announce("WARNING: Geometric data not initialized");
	}
}

///////////////////////////////////////////////////////////////////////////////

void GridPatchCSGLL::InitializeCoordinateData() {

	double dDeltaA = 0.5 * M_PI / m_grid.GetABaseResolution();

	GridSpacingGaussLobattoRepeated
		glspacing(dDeltaA, -0.25 * M_PI, m_nHorizontalOrder);

	for (int i = m_box.GetAGlobalBegin(); i < m_box.GetAGlobalEnd(); i++) {
		m_dANode[i - m_box.GetAGlobalBegin()] = glspacing.GetNode(i);
	}
	for (int i = m_box.GetAGlobalBegin(); i <= m_box.GetAGlobalEnd(); i++) {
		m_dAEdge[i - m_box.GetAGlobalBegin()] = glspacing.GetEdge(i);
	}

	for (int j = m_box.GetBGlobalBegin(); j < m_box.GetBGlobalEnd(); j++) {
		m_dBNode[j - m_box.GetBGlobalBegin()] = glspacing.GetNode(j);
	}
	for (int j = m_box.GetBGlobalBegin(); j <= m_box.GetBGlobalEnd(); j++) {
		m_dBEdge[j - m_box.GetBGlobalBegin()] = glspacing.GetEdge(j);
	}

	for (int i = 0; i < m_box.GetATotalWidth(); i++) {
	for (int j = 0; j < m_box.GetBTotalWidth(); j++) {
		CubedSphereTrans::RLLFromABP(
			m_dANode[i],
			m_dBNode[j],
			m_box.GetPanel(),
			m_dataLon(i,j),
			m_dataLat(i,j));
	}
	}

	// Compute gnomonic coordinate equivalents of alpha and beta
	m_dXNode.Allocate(m_box.GetATotalWidth());
	for (int i = 0; i < m_box.GetATotalWidth(); i++) {
		m_dXNode[i] = tan(m_dANode[i]);
	}

	m_dYNode.Allocate(m_box.GetBTotalWidth());
	for (int j = 0; j < m_box.GetBTotalWidth(); j++) {
		m_dYNode[j] = tan(m_dBNode[j]);
	}

	GridPatchGLL::InitializeCoordinateData();
}

///////////////////////////////////////////////////////////////////////////////

void GridPatchCSGLL::EvaluateTopography(
	const TestCase & test
) {
	const PhysicalConstants & phys = m_grid.GetModel().GetPhysicalConstants();

	// Compute values of topography
	for (int i = 0; i < m_box.GetATotalWidth(); i++) {
	for (int j = 0; j < m_box.GetBTotalWidth(); j++) {

		double dLon;
		double dLat;

		CubedSphereTrans::RLLFromABP(
			m_dANode[i],
			m_dBNode[j],
			m_box.GetPanel(),
			dLon,
			dLat);

		m_dataTopography(i,j) = test.EvaluateTopography(phys, dLon, dLat);
	}
	}

	// Get derivatves from basis
	GridCSGLL & gridCSGLL = dynamic_cast<GridCSGLL &>(m_grid);

	const DataArray2D<double> & dDxBasis1D = gridCSGLL.GetDxBasis1D();

	// Compute derivatives of topography
	for (int a = 0; a < GetElementCountA(); a++) {
	for (int b = 0; b < GetElementCountB(); b++) {

		for (int i = 0; i < m_nHorizontalOrder; i++) {
		for (int j = 0; j < m_nHorizontalOrder; j++) {

			// Nodal points
			int iElementA = m_box.GetAInteriorBegin() + a * m_nHorizontalOrder;
			int iElementB = m_box.GetBInteriorBegin() + b * m_nHorizontalOrder;

			int iA = iElementA + i;
			int iB = iElementB + j;

			// Topography height and its derivatives
			double dZs = m_dataTopography(iA,iB);

			double dDaZs = 0.0;
			double dDbZs = 0.0;

			for (int s = 0; s < m_nHorizontalOrder; s++) {
				dDaZs += dDxBasis1D(s,i)
					* m_dataTopography(iElementA+s,iB);
				dDbZs += dDxBasis1D(s,j)
					* m_dataTopography(iA,iElementB+s);
			}

			dDaZs /= GetElementDeltaA();
			dDbZs /= GetElementDeltaB();

			m_dataTopographyDeriv(iA,iB,0) = dDaZs;
			m_dataTopographyDeriv(iA,iB,1) = dDbZs;
		}
		}
	}
	}
}

///////////////////////////////////////////////////////////////////////////////

void GridPatchCSGLL::EvaluateGeometricTerms() {

	// Physical constants
	const PhysicalConstants & phys = m_grid.GetModel().GetPhysicalConstants();

	// 2D equation set
	bool fIs2DEquationSet = false;
	if (m_grid.GetModel().GetEquationSet().GetDimensionality() == 2) {
		fIs2DEquationSet = true;
	}

	if ((fIs2DEquationSet) && (m_grid.GetZtop() != 1.0)) {
		_EXCEPTIONT("Ztop must be 1.0 for 2D equation sets");
	}

	// Obtain Gauss Lobatto quadrature nodes and weights
	DataArray1D<double> dGL;
	DataArray1D<double> dWL;

	GaussLobattoQuadrature::GetPoints(m_nHorizontalOrder, 0.0, 1.0, dGL, dWL);

	// Obtain normalized areas in the vertical
	const DataArray1D<double> & dWNode =
		m_grid.GetREtaLevelsNormArea();
	const DataArray1D<double> & dWREdge =
		m_grid.GetREtaInterfacesNormArea();

	// Verify that normalized areas are correct
	double dWNodeSum = 0.0;
	for (int k = 0; k < dWNode.GetRows(); k++) {
		dWNodeSum += dWNode[k];
	}
	if (fabs(dWNodeSum - 1.0) > 1.0e-13) {
		_EXCEPTION1("Error in normalized areas (%1.15e)", dWNodeSum);
	}

	if (m_grid.GetVerticalStaggering() !=
	    Grid::VerticalStaggering_Interfaces
	) {
		double dWREdgeSum = 0.0;
		for (int k = 0; k < dWREdge.GetRows(); k++) {
			dWREdgeSum += dWREdge[k];
		}
		if (fabs(dWREdgeSum - 1.0) > 1.0e-13) {
			_EXCEPTION1("Error in normalized areas (%1.15e)", dWREdgeSum);
		}
	}

	// Derivatives of basis functions
	GridCSGLL & gridCSGLL = dynamic_cast<GridCSGLL &>(m_grid);

	const DataArray2D<double> & dDxBasis1D = gridCSGLL.GetDxBasis1D();

	// Initialize the Coriolis force at each node
	for (int i = 0; i < m_box.GetATotalWidth(); i++) {
	for (int j = 0; j < m_box.GetBTotalWidth(); j++) {
		m_dataCoriolisF(i,j) =
			2.0 * phys.GetOmega() * sin(m_dataLat(i,j));
	}
	}

	// Initialize metric in terrain-following coords
	for (int a = 0; a < GetElementCountA(); a++) {
	for (int b = 0; b < GetElementCountB(); b++) {

	for (int i = 0; i < m_nHorizontalOrder; i++) {
	for (int j = 0; j < m_nHorizontalOrder; j++) {

		// Nodal points
		int iElementA = m_box.GetAInteriorBegin() + a * m_nHorizontalOrder;
		int iElementB = m_box.GetBInteriorBegin() + b * m_nHorizontalOrder;

		int iA = iElementA + i;
		int iB = iElementB + j;

		// Gnomonic coordinates
		double dX = m_dXNode[iA];
		double dY = m_dYNode[iB];
		double dDelta2 = (1.0 + dX * dX + dY * dY);
		double dDelta = sqrt(dDelta2);

		// Topography height and its derivatives
		double dZs = m_dataTopography(iA,iB);
		double dDaZs = m_dataTopographyDeriv(iA,iB,0);
		double dDbZs = m_dataTopographyDeriv(iA,iB,1);

		// 2D equations
		if (fIs2DEquationSet) {
			dZs = 0.0;
			dDaZs = 0.0;
			dDbZs = 0.0;
		}

		// Initialize 2D Jacobian
		m_dataJacobian2D(iA,iB) =
			(1.0 + dX * dX) * (1.0 + dY * dY) / (dDelta * dDelta * dDelta);

		m_dataJacobian2D(iA,iB) *=
			  phys.GetEarthRadius()
			* phys.GetEarthRadius();

		// Initialize 2D contravariant metric
		double dContraMetricScale = 
			dDelta2 / (1.0 + dX * dX) / (1.0 + dY * dY)
			/ (phys.GetEarthRadius() * phys.GetEarthRadius());

		m_dataContraMetric2DA(iA,iB,0) =
			dContraMetricScale * (1.0 + dY * dY);
		m_dataContraMetric2DA(iA,iB,1) =
			dContraMetricScale * dX * dY;

		m_dataContraMetric2DB(iA,iB,0) =
			dContraMetricScale * dX * dY;
		m_dataContraMetric2DB(iA,iB,1) =
			dContraMetricScale * (1.0 + dX * dX);

		// Initialize 2D covariant metric
		double dCovMetricScale =
			phys.GetEarthRadius() * phys.GetEarthRadius()
			* (1.0 + dX * dX) * (1.0 + dY * dY)
			/ (dDelta2 * dDelta2);

		m_dataCovMetric2DA(iA,iB,0) =
			dCovMetricScale * (1.0 + dX * dX);
		m_dataCovMetric2DA(iA,iB,1) =
			dCovMetricScale * (- dX * dY);

		m_dataCovMetric2DB(iA,iB,0) =
			dCovMetricScale * (- dX * dY);
		m_dataCovMetric2DB(iA,iB,1) =
			dCovMetricScale * (1.0 + dY * dY);

		// Vertical coordinate transform and its derivatives
		for (int k = 0; k < m_grid.GetRElements(); k++) {

			// Gal-Chen and Somerville (1975) linear terrain-following coord
			double dREta = m_grid.GetREtaLevel(k);
/*
			double dREtaStretch;
			double dDxREtaStretch;
			m_grid.EvaluateVerticalStretchF(
				dREta, dREtaStretch, dDxREtaStretch);

			double dZ = dZs + (m_grid.GetZtop() - dZs) * dREtaStretch;
			double dDaR = (1.0 - dREtaStretch) * dDaZs;
			double dDbR = (1.0 - dREtaStretch) * dDbZs;
			double dDxR = (m_grid.GetZtop() - dZs) * dDxREtaStretch;
*/

			double dZ = dZs + (m_grid.GetZtop() - dZs) * dREta;
			double dDaR = (1.0 - dREta) * dDaZs;
			double dDbR = (1.0 - dREta) * dDbZs;
			double dDxR = (m_grid.GetZtop() - dZs);

			// Calculate pointwise Jacobian
			m_dataJacobian(iA,iB,k) =
				dDxR * m_dataJacobian2D(iA,iB);

			// Element area associated with each model level GLL node
			m_dataElementAreaNode(iA,iB,k) =
				m_dataJacobian(iA,iB,k)
				* dWL[i] * GetElementDeltaA()
				* dWL[j] * GetElementDeltaB()
				* dWNode[k];

			// Contravariant metric components
			m_dataContraMetricA(iA,iB,k,0) =
				m_dataContraMetric2DA(iA,iB,0);
			m_dataContraMetricA(iA,iB,k,1) =
				m_dataContraMetric2DA(iA,iB,1);
			m_dataContraMetricA(iA,iB,k,2) =
				- dContraMetricScale / dDxR * (
					(1.0 + dY * dY) * dDaR + dX * dY * dDbR);

			m_dataContraMetricB(iA,iB,k,0) =
				m_dataContraMetric2DB(iA,iB,0);
			m_dataContraMetricB(iA,iB,k,1) =
				m_dataContraMetric2DB(iA,iB,1);
			m_dataContraMetricB(iA,iB,k,2) =
				- dContraMetricScale / dDxR * (
					dX * dY * dDaR + (1.0 + dX * dX) * dDbR);

			m_dataContraMetricXi(iA,iB,k,0) =
				m_dataContraMetricA(iA,iB,k,2);
			m_dataContraMetricXi(iA,iB,k,1) =
				m_dataContraMetricB(iA,iB,k,2);
			m_dataContraMetricXi(iA,iB,k,2) =
				  1.0 / (dDxR * dDxR)
				- 1.0 / dDxR * (
					  m_dataContraMetricXi(iA,iB,k,0) * dDaR
					+ m_dataContraMetricXi(iA,iB,k,1) * dDbR);

			// Derivatives of the vertical coordinate transform
			m_dataDerivRNode(iA,iB,k,0) = dDaR;
			m_dataDerivRNode(iA,iB,k,1) = dDbR;
			m_dataDerivRNode(iA,iB,k,2) = dDxR;
		}

		// Metric terms at vertical interfaces
		for (int k = 0; k <= m_grid.GetRElements(); k++) {

			// Gal-Chen and Somerville (1975) linear terrain-following coord
			double dREta = m_grid.GetREtaInterface(k);
/*
			double dREtaStretch;
			double dDxREtaStretch;
			m_grid.EvaluateVerticalStretchF(
				dREta, dREtaStretch, dDxREtaStretch);

			double dZ = dZs + (m_grid.GetZtop() - dZs) * dREtaStretch;

			double dDaR = (1.0 - dREtaStretch) * dDaZs;
			double dDbR = (1.0 - dREtaStretch) * dDbZs;
			double dDxR = (m_grid.GetZtop() - dZs) * dDxREtaStretch;
*/
			double dZ = dZs + (m_grid.GetZtop() - dZs) * dREta;
			double dDaR = (1.0 - dREta) * dDaZs;
			double dDbR = (1.0 - dREta) * dDbZs;
			double dDxR = (m_grid.GetZtop() - dZs);

			// Calculate pointwise Jacobian
			m_dataJacobianREdge(iA,iB,k) =
				(1.0 + dX * dX) * (1.0 + dY * dY) / (dDelta * dDelta * dDelta);

			m_dataJacobianREdge(iA,iB,k) *=
				dDxR
				* phys.GetEarthRadius()
				* phys.GetEarthRadius();

			// Element area associated with each model interface GLL node
			m_dataElementAreaREdge(iA,iB,k) =
				m_dataJacobianREdge(iA,iB,k)
				* dWL[i] * GetElementDeltaA()
				* dWL[j] * GetElementDeltaB()
				* dWREdge[k];

			// Contravariant metric (alpha)
			m_dataContraMetricAREdge(iA,iB,k,0) =
				m_dataContraMetric2DA(iA,iB,0);
			m_dataContraMetricAREdge(iA,iB,k,1) =
				m_dataContraMetric2DA(iA,iB,1);
			m_dataContraMetricAREdge(iA,iB,k,2) =
				- dContraMetricScale / dDxR * (
					(1.0 + dY * dY) * dDaR + dX * dY * dDbR);

			// Contravariant metric (beta)
			m_dataContraMetricBREdge(iA,iB,k,0) =
				m_dataContraMetric2DB(iA,iB,0);
			m_dataContraMetricBREdge(iA,iB,k,1) =
				m_dataContraMetric2DB(iA,iB,1);
			m_dataContraMetricBREdge(iA,iB,k,2) =
				- dContraMetricScale / dDxR * (
					dX * dY * dDaR + (1.0 + dX * dX) * dDbR);

			// Contravariant metric (xi)
			m_dataContraMetricXiREdge(iA,iB,k,0) =
				- dContraMetricScale / dDxR * (
					(1.0 + dY * dY) * dDaR + dX * dY * dDbR);

			m_dataContraMetricXiREdge(iA,iB,k,1) =
				- dContraMetricScale / dDxR * (
					dX * dY * dDaR + (1.0 + dX * dX) * dDbR);

			m_dataContraMetricXiREdge(iA,iB,k,2) =
				  1.0 / (dDxR * dDxR)
				- 1.0 / dDxR * (
					  m_dataContraMetricXiREdge(iA,iB,k,0) * dDaR
					+ m_dataContraMetricXiREdge(iA,iB,k,1) * dDbR);

			// Derivatives of the vertical coordinate transform
			m_dataDerivRREdge(iA,iB,k,0) = dDaR;
			m_dataDerivRREdge(iA,iB,k,1) = dDbR;
			m_dataDerivRREdge(iA,iB,k,2) = dDxR;
		}
	}
	}

	}
	}
}

///////////////////////////////////////////////////////////////////////////////

void GridPatchCSGLL::EvaluateTestCase(
	const TestCase & test,
	const Time & time,
	int iDataIndex
) {
	static const int UIx = 0;
	static const int VIx = 1;

	// Initialize the data at each node
	if (m_datavecStateNode.size() == 0) {
		_EXCEPTIONT("InitializeData must be called before InitialConditions");
	}
	if (iDataIndex >= m_datavecStateNode.size()) {
		_EXCEPTIONT("Invalid iDataIndex (out of range)");
	}

	// 2D equation set
	bool fIs2DEquationSet = false;
	if (m_grid.GetModel().GetEquationSet().GetDimensionality() == 2) {
		fIs2DEquationSet = true;
	}

	// Check dimensionality
	if (fIs2DEquationSet && (m_nVerticalOrder != 1)) {
		_EXCEPTIONT("VerticalOrder / Dimensionality mismatch:\n"
			"For 2D problems vertical order must be 1.");
	}

	// Evaluate topography
	EvaluateTopography(test);

	// Physical constants
	const PhysicalConstants & phys = m_grid.GetModel().GetPhysicalConstants();

	// Data state
	DataArray4D<double> & dataStateNode =
		m_datavecStateNode[iDataIndex];
	DataArray4D<double> & dataStateREdge =
		m_datavecStateREdge[iDataIndex];
	DataArray4D<double> & dataTracers =
		m_datavecTracers[iDataIndex];

	// Initialize the vertical height in each node
	if (fIs2DEquationSet) {
		for (int i = 0; i < m_box.GetATotalWidth(); i++) {
		for (int j = 0; j < m_box.GetBTotalWidth(); j++) {
			m_dataZLevels(i,j,0) = 0.0;
			m_dataZInterfaces(i,j,0) = 0.0;
			m_dataZInterfaces(i,j,1) = 1.0;
		}
		}

	} else {
		for (int i = 0; i < m_box.GetATotalWidth(); i++) {
		for (int j = 0; j < m_box.GetBTotalWidth(); j++) {

			// Gal-Chen and Sommerville (1975) vertical coordinate
			for (int k = 0; k < m_grid.GetRElements(); k++) {
				double dREta = m_grid.GetREtaLevel(k);
/*
				double dREtaStretch;
				double dDxREtaStretch;
				m_grid.EvaluateVerticalStretchF(
					dREta, dREtaStretch, dDxREtaStretch);
*/
				m_dataZLevels(i,j,k) =
					m_dataTopography(i,j)
						+ dREta * (m_grid.GetZtop() - m_dataTopography(i,j));
			}
			for (int k = 0; k <= m_grid.GetRElements(); k++) {
				double dREta = m_grid.GetREtaInterface(k);
/*
				double dREtaStretch;
				double dDxREtaStretch;
				m_grid.EvaluateVerticalStretchF(
					dREta, dREtaStretch, dDxREtaStretch);
*/
				m_dataZInterfaces(i,j,k) =
					m_dataTopography(i,j)
						+ dREta * (m_grid.GetZtop() - m_dataTopography(i,j));
			}
		}
		}
	}

	// Initialize the Rayleigh friction strength at each node
	if (test.HasRayleighFriction()) {
		for (int i = 0; i < m_box.GetATotalWidth(); i++) {
		for (int j = 0; j < m_box.GetBTotalWidth(); j++) {
			for (int k = 0; k < m_grid.GetRElements(); k++) {
				m_dataRayleighStrengthNode(i,j,k) =
					test.EvaluateRayleighStrength(
						m_dataZLevels(i,j,k),
						m_dataLon(i,j),
						m_dataLat(i,j));
			}
			for (int k = 0; k <= m_grid.GetRElements(); k++) {
				m_dataRayleighStrengthREdge(i,j,k) =
					test.EvaluateRayleighStrength(
						m_dataZInterfaces(i,j,k),
						m_dataLon(i,j),
						m_dataLat(i,j));
			}
		}
		}
	}

	// Buffer vector for storing pointwise states
	const EquationSet & eqns = m_grid.GetModel().GetEquationSet();

	int nComponents = m_grid.GetModel().GetEquationSet().GetComponents();
	int nTracers = m_grid.GetModel().GetEquationSet().GetTracers();

	DataArray1D<double> dPointwiseState(nComponents);
	DataArray1D<double> dPointwiseRefState(nComponents);
	DataArray1D<double> dPointwiseTracers;
	DataArray1D<double> dPointwiseRefTracers;

	if (m_datavecTracers.size() > 0) {
		if (nTracers > 0) {
			dPointwiseTracers.Allocate(nTracers);
			dPointwiseRefTracers.Allocate(nTracers);
		}
	}

	// Evaluate the state on model levels
	for (int i = 0; i < m_box.GetATotalWidth(); i++) {
	for (int j = 0; j < m_box.GetBTotalWidth(); j++) {
	for (int k = 0; k < m_grid.GetRElements(); k++) {

		// Evaluate pointwise state
		test.EvaluatePointwiseState(
			phys,
			time,
			m_dataZLevels(i,j,k),
			m_dataLon(i,j),
			m_dataLat(i,j),
			dPointwiseState,
			dPointwiseTracers);

		eqns.ConvertComponents(
			phys, dPointwiseState, dPointwiseTracers);

		for (int c = 0; c < dPointwiseState.GetRows(); c++) {
			dataStateNode(c,i,j,k) =
				dPointwiseState[c];
		}

		// Transform state velocities
		double dUlon;
		double dUlat;

		dUlon = dataStateNode(UIx,i,j,k);
		dUlat = dataStateNode(VIx,i,j,k);

#if defined(PROGNOSTIC_CONTRAVARIANT_MOMENTA)
		dUlon /= phys.GetEarthRadius();
		dUlat /= phys.GetEarthRadius();

		CubedSphereTrans::VecTransABPFromRLL(
			m_dXNode[i],
			m_dYNode[j],
			m_box.GetPanel(),
			dUlon, dUlat,
			dataStateNode(UIx,i,j,k),
			dataStateNode(VIx,i,j,k));
#else
		dUlon *= phys.GetEarthRadius();
		dUlat *= phys.GetEarthRadius();

		CubedSphereTrans::CoVecTransABPFromRLL(
			m_dXNode[i],
			m_dYNode[j],
			m_box.GetPanel(),
			dUlon, dUlat,
			dataStateNode(UIx,i,j,k),
			dataStateNode(VIx,i,j,k));
#endif

		// Evaluate reference state
		if (m_grid.HasReferenceState()) {
			test.EvaluateReferenceState(
				m_grid.GetModel().GetPhysicalConstants(),
				m_dataZLevels(i,j,k),
				m_dataLon(i,j),
				m_dataLat(i,j),
				dPointwiseRefState,
				dPointwiseRefTracers);

			eqns.ConvertComponents(
				phys, dPointwiseRefState, dPointwiseRefTracers);

			for (int c = 0; c < dPointwiseRefState.GetRows(); c++) {
				m_dataRefStateNode(c,i,j,k) =
					dPointwiseRefState[c];
			}

			for (int c = 0; c < dPointwiseRefTracers.GetRows(); c++) {
				m_dataRefTracers(c,i,j,k) =
					dPointwiseRefTracers[c];
			}

			// Transform reference velocities
			dUlon = m_dataRefStateNode(UIx,i,j,k);
			dUlat = m_dataRefStateNode(VIx,i,j,k);

#if defined(PROGNOSTIC_CONTRAVARIANT_MOMENTA)
			dUlon /= phys.GetEarthRadius();
			dUlat /= phys.GetEarthRadius();

			CubedSphereTrans::VecTransABPFromRLL(
				m_dXNode[i],
				m_dYNode[j],
				m_box.GetPanel(),
				dUlon, dUlat,
				m_dataRefStateNode(UIx,i,j,k),
				m_dataRefStateNode(VIx,i,j,k));
#else
			dUlon *= phys.GetEarthRadius();
			dUlat *= phys.GetEarthRadius();

			CubedSphereTrans::CoVecTransABPFromRLL(
				m_dXNode[i],
				m_dYNode[j],
				m_box.GetPanel(),
				dUlon, dUlat,
				m_dataRefStateNode(UIx,i,j,k),
				m_dataRefStateNode(VIx,i,j,k));
#endif
		}

		// Evaluate tracers
		for (int c = 0; c < dPointwiseTracers.GetRows(); c++) {
			dataTracers(c,i,j,k) =
				dPointwiseTracers[c];
		}
	}
	}
	}

	// Evaluate the state on model interfaces
	for (int i = 0; i < m_box.GetATotalWidth(); i++) {
	for (int j = 0; j < m_box.GetBTotalWidth(); j++) {
	for (int k = 0; k <= m_grid.GetRElements(); k++) {

		// Evaluate pointwise state
		test.EvaluatePointwiseState(
			m_grid.GetModel().GetPhysicalConstants(),
			time,
			m_dataZInterfaces(i,j,k),
			m_dataLon(i,j),
			m_dataLat(i,j),
			dPointwiseState,
			dPointwiseTracers);

		eqns.ConvertComponents(
			phys, dPointwiseState, dPointwiseTracers);

		for (int c = 0; c < dPointwiseState.GetRows(); c++) {
			dataStateREdge(c,i,j,k) =
				dPointwiseState[c];
		}

		// Transform state velocities
		double dUlon;
		double dUlat;

		dUlon = dataStateREdge(UIx,i,j,k);
		dUlat = dataStateREdge(VIx,i,j,k);

#if defined(PROGNOSTIC_CONTRAVARIANT_MOMENTA)
		dUlon /= phys.GetEarthRadius();
		dUlat /= phys.GetEarthRadius();

		CubedSphereTrans::VecTransABPFromRLL(
			m_dXNode[i],
			m_dYNode[j],
			m_box.GetPanel(),
			dUlon, dUlat,
			dataStateREdge(UIx,i,j,k),
			dataStateREdge(VIx,i,j,k));
#else
		dUlon *= phys.GetEarthRadius();
		dUlat *= phys.GetEarthRadius();

		CubedSphereTrans::CoVecTransABPFromRLL(
			m_dXNode[i],
			m_dYNode[j],
			m_box.GetPanel(),
			dUlon, dUlat,
			dataStateREdge(UIx,i,j,k),
			dataStateREdge(VIx,i,j,k));
#endif

		if (m_grid.HasReferenceState()) {
			test.EvaluateReferenceState(
				m_grid.GetModel().GetPhysicalConstants(),
				m_dataZInterfaces(i,j,k),
				m_dataLon(i,j),
				m_dataLat(i,j),
				dPointwiseRefState,
				dPointwiseRefTracers);

			eqns.ConvertComponents(
				phys, dPointwiseRefState, dPointwiseRefTracers);

			for (int c = 0; c < dPointwiseState.GetRows(); c++) {
				m_dataRefStateREdge(c,i,j,k) =
					dPointwiseRefState[c];
			}

			// Transform reference velocities
			dUlon = m_dataRefStateREdge(UIx,i,j,k);
			dUlat = m_dataRefStateREdge(VIx,i,j,k);

#if defined(PROGNOSTIC_CONTRAVARIANT_MOMENTA)
			dUlon /= phys.GetEarthRadius();
			dUlat /= phys.GetEarthRadius();

			CubedSphereTrans::VecTransABPFromRLL(
				m_dXNode[i],
				m_dYNode[j],
				m_box.GetPanel(),
				dUlon, dUlat,
				m_dataRefStateREdge(UIx,i,j,k),
				m_dataRefStateREdge(VIx,i,j,k));
#else
			dUlon *= phys.GetEarthRadius();
			dUlat *= phys.GetEarthRadius();

			CubedSphereTrans::CoVecTransABPFromRLL(
				m_dXNode[i],
				m_dYNode[j],
				m_box.GetPanel(),
				dUlon, dUlat,
				m_dataRefStateREdge(UIx,i,j,k),
				m_dataRefStateREdge(VIx,i,j,k));
#endif
		}
	}
	}
	}
}

///////////////////////////////////////////////////////////////////////////////

void GridPatchCSGLL::EvaluateTestCase_StateOnly(
	const TestCase & test,
	const Time & time,
	int iDataIndex
) {
	_EXCEPTIONT("Unmaintained function");
}

///////////////////////////////////////////////////////////////////////////////

void GridPatchCSGLL::ComputeCurlAndDiv(
	const DataArray3D<double> & dataUa,
	const DataArray3D<double> & dataUb,
	const DataArray3D<double> & dataRho
) {
	const int UIx = 0;
	const int VIx = 1;

	// Parent grid
	const GridCSGLL & gridCSGLL = dynamic_cast<const GridCSGLL &>(m_grid);

#if defined(FIXED_HORIZONTAL_ORDER)
	const int nHorizontalOrder = FIXED_HORIZONTAL_ORDER;
	if (nHorizontalOrder != gridCSGLL.GetHorizontalOrder()) {
		_EXCEPTIONT("Command line order must match FIXED_HORIZONTAL_ORDER");
	}
#else
	const int nHorizontalOrder = gridCSGLL.GetHorizontalOrder();
#endif

#if defined(FIXED_RELEMENTS)
	const int nRElements = FIXED_RELEMENTS;
	if (nRElements != gridCSGLL.GetRElements()) {
		_EXCEPTIONT("Command line levels must match FIXED_RELEMENTS");
	}
#else
	const int nRElements = gridCSGLL.GetRElements();
#endif

	// Get derivatives of the basis functions
	const DataArray2D<double> & dDxBasis1D = gridCSGLL.GetDxBasis1D();

	// Get derivatives of the flux reconstruction function
	const DataArray1D<double> & dFluxDeriv1D = gridCSGLL.GetFluxDeriv1D();

	// Number of finite elements in each direction
	const int nAFiniteElements = m_box.GetAInteriorWidth() / m_nHorizontalOrder;
	const int nBFiniteElements = m_box.GetBInteriorWidth() / m_nHorizontalOrder;

	// Inverse grid spacings
	const double dInvElementDeltaA = 1.0 / GetElementDeltaA();
	const double dInvElementDeltaB = 1.0 / GetElementDeltaB();

	// Loop over all elements in the box
	for (int a = 0; a < nAFiniteElements; a++) {
	for (int b = 0; b < nBFiniteElements; b++) {

		// Index of lower-left corner node
		const int iElementA = a * m_nHorizontalOrder + m_box.GetHaloElements();
		const int iElementB = b * m_nHorizontalOrder + m_box.GetHaloElements();

		// Calculate contravariant velocities in element
		for (int i = 0; i < m_nHorizontalOrder; i++) {
		for (int j = 0; j < m_nHorizontalOrder; j++) {
#pragma simd
		for (int k = 0; k < nRElements; k++) {

			const int iA = iElementA + i;
			const int iB = iElementB + j;

#if defined(PROGNOSTIC_CONTRAVARIANT_MOMENTA)
			const double dInvRho = 1.0 / dataRho(iA,iB,k);

			m_dBufferConU(0,i,j,k) = dInvRho * (
				  m_dataCovMetric2DA(iA,iB,0)
					* dataUa(iA,iB,k)
				+ m_dataCovMetric2DA(iA,iB,1)
					* dataUb(iA,iB,k));

			m_dBufferConU(1,i,j,k) = dInvRho * (
				  m_dataCovMetric2DB(iA,iB,0)
					* dataUa(iA,iB,k)
				+ m_dataCovMetric2DB(iA,iB,1)
					* dataUb(iA,iB,k));
#else
			m_dBufferConU(0,i,j,k) =
				+ m_dataContraMetric2DA(iA,iB,0)
					* dataUa(iA,iB,k)
				+ m_dataContraMetric2DA(iA,iB,1)
					* dataUb(iA,iB,k);

			m_dBufferConU(1,i,j,k) =
				+ m_dataContraMetric2DB(iA,iB,0)
					* dataUa(iA,iB,k)
				+ m_dataContraMetric2DB(iA,iB,1)
					* dataUb(iA,iB,k);
#endif
		}
		}
		}

		// Calculate curl and div
		for (int i = 0; i < m_nHorizontalOrder; i++) {
		for (int j = 0; j < m_nHorizontalOrder; j++) {
#pragma simd
		for (int k = 0; k < nRElements; k++) {

			const int iA = iElementA + i;
			const int iB = iElementB + j;

			// Pointwise field values
			double dUa = dataUa(iA,iB,k);
			double dUb = dataUb(iA,iB,k);

			// Compute derivatives at each node
			double dDaUb = 0.0;
			double dDbUa = 0.0;

			double dDaJUa = 0.0;
			double dDbJUb = 0.0;

#if defined(PROGNOSTIC_CONTRAVARIANT_MOMENTA)

			for (int s = 0; s < m_nHorizontalOrder; s++) {
				dDaUb +=
					m_dBufferConU(VIx,s,j,k)
					* dDxBasis1D(s,i);

				dDbUa +=
					m_dBufferConU(UIx,i,s,k)
					* dDxBasis1D(s,j);

				dDaJUa +=
					m_dataJacobian2D(iElementA+s,iB)
					* dataUa(iElementA+s,iB,k)
					/ dataRho(iElementA+s,iB,k)
					* dDxBasis1D(s,i);

				dDbJUb +=
					m_dataJacobian2D(iA,iElementB+s)
					* dataUb(iA,iElementB+s,k)
					/ dataRho(iA,iElementB+s,k)
					* dDxBasis1D(s,j);
			}
#else
			for (int s = 0; s < m_nHorizontalOrder; s++) {
				dDaUb +=
					dataUb(iElementA+s,iB,k)
					* dDxBasis1D(s,i);

				dDbUa +=
					dataUa(iA,iElementB+s,k)
					* dDxBasis1D(s,j);

				dDaJUa +=
					m_dataJacobian2D(iElementA+s,iB)
					* m_dBufferConU(UIx,s,j,k)
					* dDxBasis1D(s,i);

				dDbJUb +=
					m_dataJacobian2D(iA,iElementB+s)
					* m_dBufferConU(VIx,i,s,k)
					* dDxBasis1D(s,j);
			}
#endif

			dDaUb *= dInvElementDeltaA;
			dDbUa *= dInvElementDeltaB;

			dDaJUa *= dInvElementDeltaA;
			dDbJUb *= dInvElementDeltaB;

			const double dInvJacobian2D = 1.0 / m_dataJacobian2D(iA,iB);

			m_dataDivergence(iA,iB,k) =
				(dDaJUa + dDbJUb) * dInvJacobian2D;
			m_dataVorticity(iA,iB,k) =
				(dDaUb - dDbUa) * dInvJacobian2D;
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
	const int UIx = 0;
	const int VIx = 1;
	const int HIx = 2;
	const int RIx = 4;

	// Physical constants
	const PhysicalConstants & phys = m_grid.GetModel().GetPhysicalConstants();

	// Working data
	DataArray4D<double> & dataState =
		GetDataState(iDataIndex, DataLocation_Node);

	if (dataState.GetSize(0) < 2) {
		_EXCEPTIONT(
			"Insufficient components for vorticity calculation");
	}

	// Get the alpha and beta components of vorticity
	DataArray3D<double> dataUa;
	dataUa.SetSize(
		dataState.GetSize(1),
		dataState.GetSize(2),
		dataState.GetSize(3));

	DataArray3D<double> dataUb;
	dataUb.SetSize(
		dataState.GetSize(1),
		dataState.GetSize(2),
		dataState.GetSize(3));

	DataArray3D<double> dataRho;
	dataRho.SetSize(
		dataState.GetSize(1),
		dataState.GetSize(2),
		dataState.GetSize(3));

	dataUa.AttachToData(&(dataState(UIx,0,0,0)));
	dataUb.AttachToData(&(dataState(VIx,0,0,0)));

	if (dataState.GetSize(0) == 3) {
		dataRho.AttachToData(&(dataState(HIx,0,0,0)));
	} else if (dataState.GetSize(0) > 4) {
		dataRho.AttachToData(&(dataState(RIx,0,0,0)));
	} else {
		_EXCEPTIONT("Invalid equation set");
	}

	// Compute the radial component of the curl of the velocity field
	ComputeCurlAndDiv(dataUa, dataUb, dataRho);
}

///////////////////////////////////////////////////////////////////////////////

void GridPatchCSGLL::InterpolateData(
	DataType eDataType,
	const DataArray1D<double> & dREta,
	const DataArray1D<double> & dAlpha,
	const DataArray1D<double> & dBeta,
	const DataArray1D<int> & iPatch,
	DataArray3D<double> & dInterpData,
	DataLocation eOnlyVariablesAt,
	bool fIncludeReferenceState,
	bool fConvertToPrimitive
) {
	if ((dAlpha.GetRows() != dBeta.GetRows()) ||
		(dAlpha.GetRows() != iPatch.GetRows())
	) {
		_EXCEPTIONT("Point vectors must have equivalent length.");
	}

	// Vector for storage interpolated points
	DataArray1D<double> dAInterpCoeffs(m_nHorizontalOrder);
	DataArray1D<double> dBInterpCoeffs(m_nHorizontalOrder);
	DataArray1D<double> dADiffCoeffs(m_nHorizontalOrder);
	DataArray1D<double> dBDiffCoeffs(m_nHorizontalOrder);
	DataArray1D<double> dAInterpPt(m_nHorizontalOrder);

	// Physical constants
	const PhysicalConstants & phys = m_grid.GetModel().GetPhysicalConstants();

	// Perform interpolation on all variables
	int nComponents = 0;
	int nRElements = m_grid.GetRElements();

	// Discretization type
	Grid::VerticalDiscretization eVerticalDiscType =
		m_grid.GetVerticalDiscretization();

	// State Data: Perform interpolation on all variables
	if (eDataType == DataType_State) {
		nComponents = m_datavecStateNode[0].GetSize(0);
		nRElements = m_grid.GetRElements() + 1;

	// Tracer Data: Perform interpolation on all variables
	} else if (eDataType == DataType_Tracers) {
		nComponents = m_datavecTracers[0].GetSize(0);

	// Topography Data
	} else if (eDataType == DataType_Topography) {
		nComponents = 1;
		nRElements = 1;

	// Vorticity Data
	} else if (eDataType == DataType_Vorticity) {
		nComponents = 1;

	// Divergence Data
	} else if (eDataType == DataType_Divergence) {
		nComponents = 1;

	// Temperature Data
	} else if (eDataType == DataType_Temperature) {
		nComponents = 1;

	// Surface Pressure Data
	} else if (eDataType == DataType_SurfacePressure) {
		nComponents = 1;
		nRElements = 1;

	// 2D User Data
	} else if (eDataType == DataType_Auxiliary2D) {
		nComponents = m_dataUserData2D.GetSize(0);
		nRElements = 1;

	} else {
		_EXCEPTIONT("Invalid DataType");
	}

	// Buffer storage in column
	DataArray1D<double> dColumnDataOut(dREta.GetRows());

	// Loop through all components
	for (int c = 0; c < nComponents; c++) {

		DataLocation eDataLocation = DataLocation_Node;

		if (eDataType == DataType_State) {
			eDataLocation = m_grid.GetVarLocation(c);

			// Exclude variables not at the specified DataLocation
			if ((eOnlyVariablesAt != DataLocation_None) &&
			    (eOnlyVariablesAt != eDataLocation)
			) {
				continue;
			}

			// Adjust RElements depending on state data location
			if (eDataLocation == DataLocation_Node) {
				nRElements = m_grid.GetRElements();
			} else if (eDataLocation == DataLocation_REdge) {
				nRElements = m_grid.GetRElements() + 1;
			} else {
				_EXCEPTIONT("Invalid DataLocation");
			}
		}

		// Vertical interpolation operator
		LinearColumnInterpFEM opInterp;

		if (nRElements != 1) {

			// Finite element interpolation
			if (eVerticalDiscType ==
				Grid::VerticalDiscretization_FiniteElement
			) {
				if (eDataLocation == DataLocation_Node) {
					opInterp.Initialize(
						LinearColumnInterpFEM::InterpSource_Levels,
						m_nVerticalOrder,
						m_grid.GetREtaLevels(),
						m_grid.GetREtaInterfaces(),
						dREta);

				} else if (eDataLocation == DataLocation_REdge) {
					opInterp.Initialize(
						LinearColumnInterpFEM::InterpSource_Interfaces,
						m_nVerticalOrder,
						m_grid.GetREtaLevels(),
						m_grid.GetREtaInterfaces(),
						dREta);

				} else {
					_EXCEPTIONT("Invalid DataLocation");
				}

			// Finite volume interpolation
			} else if (
				eVerticalDiscType ==
				Grid::VerticalDiscretization_FiniteVolume
			) {
				if (eDataLocation == DataLocation_Node) {
					opInterp.Initialize(
						LinearColumnInterpFEM::InterpSource_Levels,
						1,
						m_grid.GetREtaLevels(),
						m_grid.GetREtaInterfaces(),
						dREta);

				} else if (eDataLocation == DataLocation_REdge) {
					opInterp.Initialize(
						LinearColumnInterpFEM::InterpSource_Interfaces,
						1,
						m_grid.GetREtaLevels(),
						m_grid.GetREtaInterfaces(),
						dREta);

				} else {
					_EXCEPTIONT("Invalid DataLocation");
				}

			// Invalid vertical discretization type
			} else {
				_EXCEPTIONT("Invalid VerticalDiscretization");
			}

		} else {
			opInterp.InitializeIdentity(1);
		}

		// Buffer storage in column
		DataArray1D<double> dColumnData(nRElements);

		// Get a pointer to the 3D data structure
		DataArray3D<double> pData;
		DataArray3D<double> pDataRef;

		pData.SetSize(
			m_box.GetATotalWidth(),
			m_box.GetBTotalWidth(),
			nRElements);

		pDataRef.SetSize(
			m_box.GetATotalWidth(),
			m_box.GetBTotalWidth(),
			nRElements);

		if (eDataType == DataType_State) {
			if (eDataLocation == DataLocation_Node) {
				DataArray4D<double> & dataStateNode = m_datavecStateNode[0];
				pData.AttachToData(&(dataStateNode(c,0,0,0)));
				pDataRef.AttachToData(&(m_dataRefStateNode(c,0,0,0)));

			} else if (eDataLocation == DataLocation_REdge) {
				DataArray4D<double> & dataStateREdge = m_datavecStateREdge[0];
				pData.AttachToData(&(dataStateREdge(c,0,0,0)));
				pDataRef.AttachToData(&(m_dataRefStateREdge(c,0,0,0)));

			} else {
				_EXCEPTIONT("Invalid DataLocation");
			}

		} else if (eDataType == DataType_Tracers) {
			DataArray4D<double> & dataStateTracers = m_datavecTracers[0];
			pData.AttachToData(&(dataStateTracers(c,0,0,0)));

		} else if (eDataType == DataType_Topography) {
			pData.AttachToData(&(m_dataTopography(0,0)));

		} else if (eDataType == DataType_Vorticity) {
			pData.AttachToData(&(m_dataVorticity(0,0,0)));

		} else if (eDataType == DataType_Divergence) {
			pData.AttachToData(&(m_dataDivergence(0,0,0)));

		} else if (eDataType == DataType_Temperature) {
			pData.AttachToData(&(m_dataTemperature(0,0,0)));

		} else if (eDataType == DataType_SurfacePressure) {
			pData.AttachToData(&(m_dataSurfacePressure(0,0)));

		} else if (eDataType == DataType_Auxiliary2D) {
			pData.AttachToData(&(m_dataUserData2D(c,0,0)));
		}

		// Loop throught all points
		for (int i = 0; i < dAlpha.GetRows(); i++) {

			// Element index
			if (iPatch[i] != GetPatchIndex()) {
				continue;
			}

			// Verify point lies within domain of patch
			const double Eps = 1.0e-10;
			if ((dAlpha[i] < m_dAEdge[m_box.GetAInteriorBegin()] - Eps) ||
				(dAlpha[i] > m_dAEdge[m_box.GetAInteriorEnd()] + Eps) ||
				(dBeta[i] < m_dBEdge[m_box.GetBInteriorBegin()] - Eps) ||
				(dBeta[i] > m_dBEdge[m_box.GetBInteriorEnd()] + Eps)
			) {
				_EXCEPTIONT("Point out of range");
			}

			// Determine finite element index
			int iA =
				(dAlpha[i] - m_dAEdge[m_box.GetAInteriorBegin()])
					/ GetElementDeltaA();

			int iB =
				(dBeta[i] - m_dBEdge[m_box.GetBInteriorBegin()])
					/ GetElementDeltaB();

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
				&(m_dAEdge[iA]),
				dAInterpCoeffs,
				dAlpha[i]);

			PolynomialInterp::LagrangianPolynomialCoeffs(
				m_nHorizontalOrder,
				&(m_dBEdge[iB]),
				dBInterpCoeffs,
				dBeta[i]);

			// Perform interpolation on all levels
			for (int k = 0; k < nRElements; k++) {

				dColumnData[k] = 0.0;

				// Rescale vertical velocity
				const int WIx = 3;
				if ((c == WIx) && (fConvertToPrimitive)) {
					if (m_grid.GetVarLocation(WIx) == DataLocation_REdge) {
						for (int m = 0; m < m_nHorizontalOrder; m++) {
						for (int n = 0; n < m_nHorizontalOrder; n++) {

#if defined(PROGNOSTIC_CONTRAVARIANT_MOMENTA)
							dColumnData[k] +=
								  dAInterpCoeffs[m]
								* dBInterpCoeffs[n]
								* pData(iA+m,iB+n,k);
#else
							dColumnData[k] +=
								  dAInterpCoeffs[m]
								* dBInterpCoeffs[n]
								* pData(iA+m,iB+n,k)
								/ m_dataDerivRREdge(iA,iB,k,2);
#endif
						}
						}

					} else {
						for (int m = 0; m < m_nHorizontalOrder; m++) {
						for (int n = 0; n < m_nHorizontalOrder; n++) {
							dColumnData[k] +=
								  dAInterpCoeffs[m]
								* dBInterpCoeffs[n]
								* pData(iA+m,iB+n,k)
								/ m_dataDerivRNode(iA,iB,k,2);
						}
						}
					}

				} else {
					for (int m = 0; m < m_nHorizontalOrder; m++) {
					for (int n = 0; n < m_nHorizontalOrder; n++) {
						dColumnData[k] +=
							  dAInterpCoeffs[m]
							* dBInterpCoeffs[n]
							* pData(iA+m,iB+n,k);
					}
					}
				}

				// Do not include the reference state
				if ((eDataType == DataType_State) &&
					(!fIncludeReferenceState)
				) {
					for (int m = 0; m < m_nHorizontalOrder; m++) {
					for (int n = 0; n < m_nHorizontalOrder; n++) {
						dColumnData[k] -=
							  dAInterpCoeffs[m]
							* dBInterpCoeffs[n]
							* pDataRef(iA+m,iB+n,k);
					}
					}
				}
			}

			// Interpolate vertically
			opInterp.Apply(
				&(dColumnData[0]),
				&(dColumnDataOut[0]));

			// Store data
			for (int k = 0; k < dREta.GetRows(); k++) {
				dInterpData(c,k,i) =
					dColumnDataOut[k];
			}

		}
	}

	// Convert to primitive variables
	if ((eDataType == DataType_State) && (fConvertToPrimitive)) {

		const int UIx = 0;
		const int VIx = 1;
		const int WIx = 3;
		int RIx = 4;
		if (dInterpData.GetRows() < 5) {
			if (dInterpData.GetRows() == 3) {
				RIx = 2;
			} else {
				_EXCEPTIONT("Invalid EquationSet");
			}
		}

		for (int i = 0; i < dAlpha.GetRows(); i++) { 
			if (iPatch[i] != GetPatchIndex()) {
				continue;
			}

			for (int k = 0; k < dREta.GetRows(); k++) {
#if defined(PROGNOSTIC_CONTRAVARIANT_MOMENTA)
				double dUalpha =
					dInterpData(UIx,k,i)
					/ dInterpData(RIx,k,i)
					* phys.GetEarthRadius();
				double dUbeta =
					dInterpData(VIx,k,i)
					/ dInterpData(RIx,k,i)
					* phys.GetEarthRadius();

				CubedSphereTrans::VecTransRLLFromABP(
					tan(dAlpha[i]),
					tan(dBeta[i]),
					GetPatchBox().GetPanel(),
					dUalpha,
					dUbeta,
					dInterpData(UIx,k,i),
					dInterpData(VIx,k,i));
#else
				double dUalpha =
					dInterpData(UIx,k,i) / phys.GetEarthRadius();
				double dUbeta =
					dInterpData(VIx,k,i) / phys.GetEarthRadius();

				CubedSphereTrans::CoVecTransRLLFromABP(
					tan(dAlpha[i]),
					tan(dBeta[i]),
					GetPatchBox().GetPanel(),
					dUalpha,
					dUbeta,
					dInterpData(UIx,k,i),
					dInterpData(VIx,k,i));
#endif
			}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void GridPatchCSGLL::TransformHaloVelocities(
	int iDataUpdate
) {
	// Indices of velocities
	static const int UIx = 0;
	static const int VIx = 1;

	// Velocity data
	DataArray4D<double> * pDataVelocity =
		&(GetDataState(iDataUpdate, m_grid.GetVarLocation(UIx)));

	if (pDataVelocity->GetSize(0) < 2) {
		_EXCEPTIONT("Invalid number of components.");
	}

	// Panels in each coordinate direction
	const int ixRightPanel  = GetNeighborPanel(Direction_Right);
	const int ixTopPanel    = GetNeighborPanel(Direction_Top);
	const int ixLeftPanel   = GetNeighborPanel(Direction_Left);
	const int ixBottomPanel = GetNeighborPanel(Direction_Bottom);

	const int ixTopRightPanel    = GetNeighborPanel(Direction_TopRight);
	const int ixTopLeftPanel     = GetNeighborPanel(Direction_TopLeft);
	const int ixBottomLeftPanel  = GetNeighborPanel(Direction_BottomLeft);
	const int ixBottomRightPanel = GetNeighborPanel(Direction_BottomRight);

	// Post-process velocities across right edge
	if (ixRightPanel != m_box.GetPanel()) {
		const int jBegin = m_box.GetBInteriorBegin()-1;
		const int jEnd = m_box.GetBInteriorEnd()+1;

		int i = m_box.GetAInteriorEnd();
		for (int k = 0; k < pDataVelocity->GetSize(3); k++) {
		for (int j = jBegin; j < jEnd; j++) {
#if defined(PROGNOSTIC_CONTRAVARIANT_MOMENTA)
			CubedSphereTrans::VecPanelTrans(
				ixRightPanel,
				m_box.GetPanel(),
				(*pDataVelocity)(UIx,i,j,k),
				(*pDataVelocity)(VIx,i,j,k),
				m_dXNode[i],
				m_dYNode[j]);
#else
			CubedSphereTrans::CoVecPanelTrans(
				ixRightPanel,
				m_box.GetPanel(),
				(*pDataVelocity)(UIx,i,j,k),
				(*pDataVelocity)(VIx,i,j,k),
				m_dXNode[i],
				m_dYNode[j]);
#endif
		}
		}
	}

	// Post-process velocities across top edge
	if (ixTopPanel != m_box.GetPanel()) {
		const int iBegin = m_box.GetAInteriorBegin()-1;
		const int iEnd = m_box.GetAInteriorEnd()+1;

		int j = m_box.GetBInteriorEnd();
		for (int k = 0; k < pDataVelocity->GetSize(3); k++) {
		for (int i = iBegin; i < iEnd; i++) {
#if defined(PROGNOSTIC_CONTRAVARIANT_MOMENTA)
			CubedSphereTrans::VecPanelTrans(
				ixTopPanel,
				m_box.GetPanel(),
				(*pDataVelocity)(UIx,i,j,k),
				(*pDataVelocity)(VIx,i,j,k),
				m_dXNode[i],
				m_dYNode[j]);
#else
			CubedSphereTrans::CoVecPanelTrans(
				ixTopPanel,
				m_box.GetPanel(),
				(*pDataVelocity)(UIx,i,j,k),
				(*pDataVelocity)(VIx,i,j,k),
				m_dXNode[i],
				m_dYNode[j]);
#endif
		}
		}
	}

	// Post-process velocities across left edge
	if (ixLeftPanel != m_box.GetPanel()) {
		const int jBegin = m_box.GetBInteriorBegin()-1;
		const int jEnd = m_box.GetBInteriorEnd()+1;

		int i = m_box.GetAInteriorBegin()-1;
		for (int k = 0; k < pDataVelocity->GetSize(3); k++) {
		for (int j = jBegin; j < jEnd; j++) {
#if defined(PROGNOSTIC_CONTRAVARIANT_MOMENTA)
			CubedSphereTrans::VecPanelTrans(
				ixLeftPanel,
				m_box.GetPanel(),
				(*pDataVelocity)(UIx,i,j,k),
				(*pDataVelocity)(VIx,i,j,k),
				m_dXNode[i],
				m_dYNode[j]);
#else
			CubedSphereTrans::CoVecPanelTrans(
				ixLeftPanel,
				m_box.GetPanel(),
				(*pDataVelocity)(UIx,i,j,k),
				(*pDataVelocity)(VIx,i,j,k),
				m_dXNode[i],
				m_dYNode[j]);
#endif
		}
		}
	}

	// Post-process velocities across bottom edge
	if (ixBottomPanel != m_box.GetPanel()) {
		const int iBegin = m_box.GetAInteriorBegin()-1;
		const int iEnd = m_box.GetAInteriorEnd()+1;

		int j = m_box.GetBInteriorBegin()-1;
		for (int k = 0; k < pDataVelocity->GetSize(3); k++) {
		for (int i = iBegin; i < iEnd; i++) {
#if defined(PROGNOSTIC_CONTRAVARIANT_MOMENTA)
			CubedSphereTrans::VecPanelTrans(
				ixBottomPanel,
				m_box.GetPanel(),
				(*pDataVelocity)(UIx,i,j,k),
				(*pDataVelocity)(VIx,i,j,k),
				m_dXNode[i],
				m_dYNode[j]);
#else
			CubedSphereTrans::CoVecPanelTrans(
				ixBottomPanel,
				m_box.GetPanel(),
				(*pDataVelocity)(UIx,i,j,k),
				(*pDataVelocity)(VIx,i,j,k),
				m_dXNode[i],
				m_dYNode[j]);
#endif
		}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void GridPatchCSGLL::TransformTopographyDeriv() {

	// Panels in each coordinate direction
	const int ixRightPanel  = GetNeighborPanel(Direction_Right);
	const int ixTopPanel    = GetNeighborPanel(Direction_Top);
	const int ixLeftPanel   = GetNeighborPanel(Direction_Left);
	const int ixBottomPanel = GetNeighborPanel(Direction_Bottom);

	const int ixTopRightPanel    = GetNeighborPanel(Direction_TopRight);
	const int ixTopLeftPanel     = GetNeighborPanel(Direction_TopLeft);
	const int ixBottomLeftPanel  = GetNeighborPanel(Direction_BottomLeft);
	const int ixBottomRightPanel = GetNeighborPanel(Direction_BottomRight);

	// Post-process velocities across right edge
	if (ixRightPanel != m_box.GetPanel()) {
		const int jBegin = m_box.GetBInteriorBegin()-1;
		const int jEnd = m_box.GetBInteriorEnd()+1;

		const int i = m_box.GetAInteriorEnd();
		for (int j = jBegin; j < jEnd; j++) {
			CubedSphereTrans::CoVecPanelTrans(
				ixRightPanel,
				m_box.GetPanel(),
				m_dataTopographyDeriv(i,j,0),
				m_dataTopographyDeriv(i,j,1),
				m_dXNode[i],
				m_dYNode[j]);
		}
	}

	// Post-process velocities across top edge
	if (ixTopPanel != m_box.GetPanel()) {
		const int iBegin = m_box.GetAInteriorBegin()-1;
		const int iEnd = m_box.GetAInteriorEnd()+1;

		int j = m_box.GetBInteriorEnd();
		for (int i = iBegin; i < iEnd; i++) {
			CubedSphereTrans::CoVecPanelTrans(
				ixTopPanel,
				m_box.GetPanel(),
				m_dataTopographyDeriv(i,j,0),
				m_dataTopographyDeriv(i,j,1),
				m_dXNode[i],
				m_dYNode[j]);
		}
	}

	// Post-process velocities across left edge
	if (ixLeftPanel != m_box.GetPanel()) {
		const int jBegin = m_box.GetBInteriorBegin()-1;
		const int jEnd = m_box.GetBInteriorEnd()+1;

		int i = m_box.GetAInteriorBegin()-1;
		for (int j = jBegin; j < jEnd; j++) {
			CubedSphereTrans::CoVecPanelTrans(
				ixLeftPanel,
				m_box.GetPanel(),
				m_dataTopographyDeriv(i,j,0),
				m_dataTopographyDeriv(i,j,1),
				m_dXNode[i],
				m_dYNode[j]);
		}
	}

	// Post-process velocities across bottom edge
	if (ixBottomPanel != m_box.GetPanel()) {
		const int iBegin = m_box.GetAInteriorBegin()-1;
		const int iEnd = m_box.GetAInteriorEnd()+1;

		int j = m_box.GetBInteriorBegin()-1;
		for (int i = iBegin; i < iEnd; i++) {
			CubedSphereTrans::CoVecPanelTrans(
				ixBottomPanel,
				m_box.GetPanel(),
				m_dataTopographyDeriv(i,j,0),
				m_dataTopographyDeriv(i,j,1),
				m_dXNode[i],
				m_dYNode[j]);
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

