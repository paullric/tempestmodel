
///////////////////////////////////////////////////////////////////////////////
///
///	\file    GridPatchCartesianGLL.cpp
///	\author  Paul Ullrich, Jorge Guerra
///	\version September 26, 2013
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

#include "GridPatchCartesianGLL.h"
#include "GridCartesianGLL.h"
#include "Model.h"
#include "TestCase.h"
#include "GridSpacing.h"
#include "VerticalStretch.h"
#include "Defines.h"

#include "Direction.h"
#include "CubedSphereTrans.h"
#include "PolynomialInterp.h"
#include "GaussLobattoQuadrature.h"

#include "Announce.h"
#include "MathHelper.h"

#include <cmath>
#include <algorithm>

///////////////////////////////////////////////////////////////////////////////

GridPatchCartesianGLL::GridPatchCartesianGLL(
	GridCartesianGLL & grid,
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

}

///////////////////////////////////////////////////////////////////////////////

void GridPatchCartesianGLL::InitializeDataLocal(
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

	// Physical constants
	const PhysicalConstants & phys = m_grid.GetModel().GetPhysicalConstants();

	// Initialize coordinate data
	if (fAllocateGeometric) {
		InitializeCoordinateData();
	} else {
		//Announce("WARNING: Geometric data not initialized");
	}
}

///////////////////////////////////////////////////////////////////////////////

void GridPatchCartesianGLL::InitializeCoordinateData() {

	// Get the Cartesian grid
	GridCartesianGLL & gridCartesianGLL =
		dynamic_cast<GridCartesianGLL &>(m_grid);

	double dElementDeltaA = (gridCartesianGLL.GetMaximumX() -
							 gridCartesianGLL.GetMinimumX())
		/ static_cast<double>(m_grid.GetABaseResolution());
	double dElementDeltaB = (gridCartesianGLL.GetMaximumY() - 
							 gridCartesianGLL.GetMinimumY())
		/ static_cast<double>(m_grid.GetBBaseResolution());

	
	double X0 = gridCartesianGLL.GetMinimumX();
	double Y0 = gridCartesianGLL.GetMinimumY();
	GridSpacingGaussLobattoRepeated
		glspacingA(dElementDeltaA, X0, m_nHorizontalOrder);
	GridSpacingGaussLobattoRepeated
		glspacingB(dElementDeltaB, Y0, m_nHorizontalOrder);

	for (int i = m_box.GetAGlobalBegin(); i < m_box.GetAGlobalEnd(); i++) {
		m_dANode[i - m_box.GetAGlobalBegin()] = glspacingA.GetNode(i);
	}
	for (int i = m_box.GetAGlobalBegin(); i <= m_box.GetAGlobalEnd(); i++) {
		m_dAEdge[i - m_box.GetAGlobalBegin()] = glspacingA.GetEdge(i);
	}
	for (int j = m_box.GetBGlobalBegin(); j < m_box.GetBGlobalEnd(); j++) {
		m_dBNode[j - m_box.GetBGlobalBegin()] = glspacingB.GetNode(j);
	}
	for (int j = m_box.GetBGlobalBegin(); j <= m_box.GetBGlobalEnd(); j++) {
		m_dBEdge[j - m_box.GetBGlobalBegin()] = glspacingB.GetEdge(j);
	}

	for (int i = 0; i < m_box.GetATotalWidth(); i++) {
	for (int j = 0; j < m_box.GetBTotalWidth(); j++) {
		m_dataLon[i][j] = m_dANode[i];
		m_dataLat[i][j] = m_dBNode[j];
	}
	}
	GridPatchGLL::InitializeCoordinateData();
}

///////////////////////////////////////////////////////////////////////////////

void GridPatchCartesianGLL::EvaluateTopography(
	const TestCase & test
) {
	const PhysicalConstants & phys = m_grid.GetModel().GetPhysicalConstants();
	
	// Get the cartesian grid
	GridCartesianGLL & gridCartesianGLL =
		dynamic_cast<GridCartesianGLL &>(m_grid);

	// Compute values of topography and find the top
	for (int i = 0; i < m_box.GetATotalWidth(); i++) {
	for (int j = 0; j < m_box.GetBTotalWidth(); j++) {

		double dX = m_dANode[i];
		double dY = m_dBNode[j];

		m_dataTopography[i][j] = test.EvaluateTopography(phys, dX, dY);

		if (m_dataTopography[i][j] >= m_grid.GetZtop()) {
			_EXCEPTIONT("TestCase topography exceeds model top.");
		}

	}
	}

	// Get derivatves from basis
	const DataArray2D<double> & dDxBasis1D =
		gridCartesianGLL.GetDxBasis1D();

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
			double dZs = m_dataTopography[iA][iB];

			double dDaZs = 0.0;
			double dDbZs = 0.0;

			for (int s = 0; s < m_nHorizontalOrder; s++) {
				dDaZs += dDxBasis1D[s][i] * m_dataTopography[iElementA+s][iB];
				dDbZs += dDxBasis1D[s][j] * m_dataTopography[iA][iElementB+s];
			}

			dDaZs /= GetElementDeltaA();
			dDbZs /= GetElementDeltaB();

			m_dataTopographyDeriv[0][iA][iB] = dDaZs;
			m_dataTopographyDeriv[1][iA][iB] = dDbZs;
		}
		}
	}
	}
}

///////////////////////////////////////////////////////////////////////////////

void GridPatchCartesianGLL::EvaluateGeometricTerms() {

	// Physical constants
	const PhysicalConstants & phys = m_grid.GetModel().GetPhysicalConstants();

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
	GridCartesianGLL & gridCartesianGLL =
		dynamic_cast<GridCartesianGLL &>(m_grid);

	const DataArray2D<double> & dDxBasis1D = gridCartesianGLL.GetDxBasis1D();

	// Initialize the Coriolis force at each node
	bool fCartesianXZ = gridCartesianGLL.GetIsCartesianXZ();
	double dRefLat = gridCartesianGLL.GetReferenceLatitude();
	double dy0 = 0.5 * fabs(gridCartesianGLL.GetMaximumY() - 
							gridCartesianGLL.GetMinimumY());
	double dfp = 2.0 * phys.GetOmega() * sin(dRefLat);
	double dbetap = 2.0 * phys.GetOmega() * cos(dRefLat) / 
					phys.GetEarthRadius();
	for (int i = 0; i < m_box.GetATotalWidth(); i++) {
	for (int j = 0; j < m_box.GetBTotalWidth(); j++) {
		if (!fCartesianXZ) {
			// Coriolis force by beta approximation
			m_dataCoriolisF[i][j] = dfp + dbetap * (m_dataLat[i][j] - dy0);
			//m_dataCoriolisF[i][j] = dfp;
		}
		else {
			m_dataCoriolisF[i][j] = 0.0;
		}
	}
	}

	// Initialize metric and Christoffel symbols in terrain-following coords
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
			double dZs = m_dataTopography[iA][iB];
			double dDaZs = m_dataTopographyDeriv[0][iA][iB];
			double dDbZs = m_dataTopographyDeriv[1][iA][iB];

			// Initialize 2D Jacobian
			m_dataJacobian2D[iA][iB] = 1.0;

			// Initialize 2D contravariant metric
			m_dataContraMetric2DA[iA][iB][0] = 1.0;
			m_dataContraMetric2DA[iA][iB][1] = 0.0;

			m_dataContraMetric2DB[iA][iB][0] = 0.0;
			m_dataContraMetric2DB[iA][iB][1] = 1.0;

			// Initialize 2D covariant metric
			m_dataCovMetric2DA[iA][iB][0] = 1.0;
			m_dataCovMetric2DA[iA][iB][1] = 0.0;

			m_dataCovMetric2DB[iA][iB][0] = 0.0;
			m_dataCovMetric2DB[iA][iB][1] = 1.0;

			// Metric terms at vertical levels
			for (int k = 0; k < m_grid.GetRElements(); k++) {

				// Gal-Chen and Somerville (1975) terrain following coord
				// 2nd order polynomial decay terrain following coord
				// 4th order polynomial decay terrain following coord
				// Schar Exponential Decay terrain following coord
				double dREta = m_grid.GetREtaLevel(k);

				double dZ = dZs + (m_grid.GetZtop() - dZs) * dREta;
				double dDaZ = (1.0 - dREta) * dDaZs;
				double dDbZ = (1.0 - dREta) * dDbZs;
				double dDxZ = (m_grid.GetZtop() - dZs);
				//double dDXZs = -hc * exp(-pow(m_dataLon[iA][iB] / ac,2.0))
				//* (M_PI / lm * sin(2 * M_PI * m_dataLon[iA][iB] / lm)
				//+ 2.0 * m_dataLon[iA][iB] / (ac * ac) 
				//* pow(cos(M_PI * m_dataLon[iA][iB] / lm),2.0));
				//double dDaX = sqrt(1.0 + dDXZs * dDXZs);
				//double dDaZ = (1.0 - dREta) * dDXZs * dDaX / (1.0 + dDXZs * dDXZs);

/*
				double dZ = m_grid.GetZtop() * dREta + 
							std::pow((1.0 - dREta), 2.0) * dZs;
				double dDaZ = std::pow((1.0 - dREta), 2.0) * dDaZs;
				double dDbZ = std::pow((1.0 - dREta), 2.0) * dDbZs;
				double dDxZ = m_grid.GetZtop() - 
							2.0 * (1.0 - dREta) * dZs;
*/
/*
				double dZ = m_grid.GetZtop() * dREta + 
							std::pow((1.0 - dREta), 4.0) * dZs;
				double dDaZ = std::pow((1.0 - dREta), 4.0) * dDaZs;
				double dDbZ = std::pow((1.0 - dREta), 4.0) * dDbZs;
				double dDxZ = m_grid.GetZtop() - 
							4.0 * std::pow((1.0 - dREta), 3.0) * dZs;
*/
/*
				double dZ = m_grid.GetZtop() * dREta 
					+ std::pow(std::cos(0.5 * M_PI * dREta), 4.0) * dZs;
				double dDaZ = std::pow(std::cos(0.5 * M_PI * dREta), 4.0)
					* dDaZs;
				double dDbZ = std::pow(std::cos(0.5 * M_PI * dREta), 4.0)
					* dDbZs;
				double dDxZ = m_grid.GetZtop() - 2.0 * M_PI
					* std::sin(0.5 * M_PI * dREta)
					* std::pow(std::cos(0.5 * M_PI * dREta), 3.0) * dZs;
*/
/*
				double dZ = m_grid.GetZtop() * dREta 
					+ std::pow(std::cos(0.5 * M_PI * dREta), 8.0) * dZs;
				double dDaZ = std::pow(std::cos(0.5 * M_PI * dREta), 8.0)
					* dDaZs;
				double dDbZ = std::pow(std::cos(0.5 * M_PI * dREta), 8.0)
					* dDbZs;
				double dDxZ = m_grid.GetZtop() - 4.0 * M_PI
					* std::sin(0.5 * M_PI * dREta)
					* std::pow(std::cos(0.5 * M_PI * dREta), 7.0) * dZs;
*/
/*
				double dScaleH = gridCartesianGLL.m_dSL;
				double dH = m_grid.GetZtop();
				double dbZ = sinh(dH * (1.0 - dREta) / dScaleH)
							/ sinh(dH / dScaleH);
				double dZ = m_grid.GetZtop() * dREta + dZs * dbZ;
				double dDaZ = dbZ * dDaZs;
				double dDbZ = dbZ * dDbZs;
				double dDxZ = dH - dZs * dH * 
					cosh(dH * (1.0 - dREta) / dScaleH) /
					(dScaleH * sinh(dH / dScaleH));
*/
/*
				double power = 10.0;
				double botRate = 1.0;
				double dZ = m_grid.GetZtop() * dREta 
					+ (1.0 - botRate * dREta)
					* std::pow(std::cos(0.5 * M_PI * dREta), power) * dZs;
				double dDaZ = (1.0 - botRate * dREta)
					* std::pow(std::cos(0.5 * M_PI * dREta), power)
					* dDaZs;
				double dDbZ = (1.0 - botRate * dREta)
					* std::pow(std::cos(0.5 * M_PI * dREta), power)
					* dDbZs;
				double dDxZ = m_grid.GetZtop() - 0.5 * power * M_PI
					* (1.0 - botRate * dREta) * std::sin(0.5 * M_PI * dREta)
					* std::pow(std::cos(0.5 * M_PI * dREta), power - 1.0) * dZs
					- botRate * std::pow(std::cos(0.5 * M_PI * dREta), power) * dZs;
*/
//printf("%.10E %.10E %.10E %.10E %.10E %.10E \n",m_dataLon[iA][iB],m_dataLat[iA][iB],dZ,dDaZ,dDbZ,dDxZ);

				// Calculate pointwise Jacobian
				m_dataJacobian[k][iA][iB] =
					dDxZ * m_dataJacobian2D[iA][iB];

				// Element area associated with each model level GLL node
				m_dataElementArea[k][iA][iB] =
					m_dataJacobian[k][iA][iB]
					* dWL[i] * GetElementDeltaA()
					* dWL[j] * GetElementDeltaB()
					* dWNode[k];

				// Contravariant metric components
				m_dataContraMetricA[k][iA][iB][0] =
					m_dataContraMetric2DA[iA][iB][0];
				m_dataContraMetricA[k][iA][iB][1] =
					m_dataContraMetric2DA[iA][iB][1];
				m_dataContraMetricA[k][iA][iB][2] =
					- dDaZ / dDxZ;

				m_dataContraMetricB[k][iA][iB][0] =
					m_dataContraMetric2DB[iA][iB][0];
				m_dataContraMetricB[k][iA][iB][1] =
					m_dataContraMetric2DB[iA][iB][1];
				m_dataContraMetricB[k][iA][iB][2] =
					- dDbZ / dDxZ;

				m_dataContraMetricXi[k][iA][iB][0] =
					m_dataContraMetricA[k][iA][iB][2];
				m_dataContraMetricXi[k][iA][iB][1] =
					m_dataContraMetricB[k][iA][iB][2];
				m_dataContraMetricXi[k][iA][iB][2] =
					(1.0 + dDaZ * dDaZ + dDbZ * dDbZ) / (dDxZ * dDxZ);

				// Covariant metric components
				m_dataCovMetricA[k][iA][iB][0] =
					m_dataCovMetric2DA[iA][iB][0] + dDaZ * dDaZ;
				m_dataCovMetricA[k][iA][iB][1] =
					m_dataCovMetric2DA[iA][iB][1] + dDaZ * dDbZ;
				m_dataCovMetricA[k][iA][iB][2] =
					dDaZ * dDxZ;

				m_dataCovMetricB[k][iA][iB][0] =
					m_dataCovMetric2DB[iA][iB][0] + dDbZ * dDaZ;
				m_dataCovMetricB[k][iA][iB][1] =
					m_dataCovMetric2DB[iA][iB][1] + dDbZ * dDbZ;
				m_dataCovMetricB[k][iA][iB][2] =
					dDbZ * dDxZ;

				m_dataCovMetricXi[k][iA][iB][0] =
					dDaZ * dDxZ;
				m_dataCovMetricXi[k][iA][iB][1] =
					dDbZ * dDxZ;
				m_dataCovMetricXi[k][iA][iB][2] =
					dDxZ * dDxZ;

				// Derivatives of the vertical coordinate transform
				m_dataDerivRNode[k][iA][iB][0] = dDaZ;
				m_dataDerivRNode[k][iA][iB][1] = dDbZ;
				m_dataDerivRNode[k][iA][iB][2] = dDxZ;
			}

			// Metric terms at vertical interfaces
			for (int k = 0; k <= m_grid.GetRElements(); k++) {

				// Gal-Chen and Somerville (1975) terrain following coord
				// 2nd order polynomial decay terrain following coord
				// 4th order polynomial decay terrain following coord
				// Schar Exponential decay terrain following coord
				double dREta = m_grid.GetREtaInterface(k);
				double dZ = dZs + (m_grid.GetZtop() - dZs) * dREta;
				double dDaZ = (1.0 - dREta) * dDaZs;
				double dDbZ = (1.0 - dREta) * dDbZs;
				double dDxZ = (m_grid.GetZtop() - dZs);
				//double dDXZs = -hc * exp(-pow(m_dataLon[iA][iB] / ac,2.0))
				//* (M_PI / lm * sin(2 * M_PI * m_dataLon[iA][iB] / lm)
				//+ 2.0 * m_dataLon[iA][iB] / (ac * ac) 
				//* pow(cos(M_PI * m_dataLon[iA][iB] / lm),2.0));
				//double dDaX = sqrt(1.0 + dDXZs * dDXZs);
				//double dDaZ = (1.0 - dREta) * dDXZs * dDaX / (1.0 + dDXZs * dDXZs);
/*
				double dZ = m_grid.GetZtop() * dREta + 
							std::pow((1.0 - dREta), 2.0) * dZs;
				double dDaZ = std::pow((1.0 - dREta), 2.0) * dDaZs;
				double dDbZ = std::pow((1.0 - dREta), 2.0) * dDbZs;
				double dDxZ = m_grid.GetZtop() - 
							2.0 * (1.0 - dREta) * dZs;
*/
/*
				double dZ = m_grid.GetZtop() * dREta + 
							std::pow((1.0 - dREta), 4.0) * dZs;
				double dDaZ = std::pow((1.0 - dREta), 4.0) * dDaZs;
				double dDbZ = std::pow((1.0 - dREta), 4.0) * dDbZs;
				double dDxZ = m_grid.GetZtop() - 
							4.0 * std::pow((1.0 - dREta), 3.0) * dZs;
*/
/*
				double dZ = m_grid.GetZtop() * dREta 
					+ std::pow(std::cos(0.5 * M_PI * dREta), 4.0) * dZs;
				double dDaZ = std::pow(std::cos(0.5 * M_PI * dREta), 4.0)
					* dDaZs;
				double dDbZ = std::pow(std::cos(0.5 * M_PI * dREta), 4.0)
					* dDbZs;
				double dDxZ = m_grid.GetZtop() - 2.0 * M_PI
					* std::sin(0.5 * M_PI * dREta)
					* std::pow(std::cos(0.5 * M_PI * dREta), 3.0) * dZs;
*/
/*
				double dZ = m_grid.GetZtop() * dREta 
					+ std::pow(std::cos(0.5 * M_PI * dREta), 8.0) * dZs;
				double dDaZ = std::pow(std::cos(0.5 * M_PI * dREta), 8.0)
					* dDaZs;
				double dDbZ = std::pow(std::cos(0.5 * M_PI * dREta), 8.0)
					* dDbZs;
				double dDxZ = m_grid.GetZtop() - 2.0 * M_PI
					* std::sin(0.5 * M_PI * dREta)
					* std::pow(std::cos(0.5 * M_PI * dREta), 7.0) * dZs;
*/
/*
				double dScaleH = gridCartesianGLL.m_dSL;
				double dH = m_grid.GetZtop();
				double dbZ = sinh(dH * (1.0 - dREta) / dScaleH)
							/ sinh(dH / dScaleH);
				double dZ = m_grid.GetZtop() * dREta + dZs * dbZ;
				double dDaZ = dbZ * dDaZs;
				double dDbZ = dbZ * dDbZs;
				double dDxZ = dH - dZs * dH * 
					cosh(dH * (1.0 - dREta) / dScaleH) /
					(dScaleH * sinh(dH / dScaleH));
*/
/*
				double power = 10.0;
				double botRate = 1.0;
				double dZ = m_grid.GetZtop() * dREta 
					+ (1.0 - botRate * dREta)
					* std::pow(std::cos(0.5 * M_PI * dREta), power) * dZs;
				double dDaZ = (1.0 - botRate * dREta)
					* std::pow(std::cos(0.5 * M_PI * dREta), power)
					* dDaZs;
				double dDbZ = (1.0 - botRate * dREta)
					* std::pow(std::cos(0.5 * M_PI * dREta), power)
					* dDbZs;
				double dDxZ = m_grid.GetZtop() - 0.5 * power * M_PI
					* (1.0 - botRate * dREta) * std::sin(0.5 * M_PI * dREta)
					* std::pow(std::cos(0.5 * M_PI * dREta), power - 1.0) * dZs
					- botRate * std::pow(std::cos(0.5 * M_PI * dREta), power) * dZs;
*/
//printf("%.16E %.16E %.16E %.16E %.16E %.16E \n",m_dataLon[iA][iB],m_dataLat[iA][iB],dZ,dDaZ,dDbZ,dDxZ);

				// Calculate pointwise Jacobian
				m_dataJacobianREdge[k][iA][iB] =
					dDxZ * m_dataJacobian2D[iA][iB];

				// Element area associated with each model interface GLL node
				m_dataElementAreaREdge[k][iA][iB] =
					m_dataJacobianREdge[k][iA][iB]
					* dWL[i] * GetElementDeltaA()
					* dWL[j] * GetElementDeltaB()
					* dWREdge[k];

				// Components of the contravariant metric
				m_dataContraMetricAREdge[k][iA][iB][0] =
					m_dataContraMetric2DA[iA][iB][0];
				m_dataContraMetricAREdge[k][iA][iB][1] =
					m_dataContraMetric2DA[iA][iB][1];
				m_dataContraMetricAREdge[k][iA][iB][2] =
					- dDaZ / dDxZ;

				m_dataContraMetricBREdge[k][iA][iB][0] =
					m_dataContraMetric2DB[iA][iB][0];
				m_dataContraMetricBREdge[k][iA][iB][1] =
					m_dataContraMetric2DB[iA][iB][1];
				m_dataContraMetricBREdge[k][iA][iB][2] =
					- dDbZ / dDxZ;

				m_dataContraMetricXiREdge[k][iA][iB][0] =
					- dDaZ / dDxZ;
				m_dataContraMetricXiREdge[k][iA][iB][1] =
					- dDbZ / dDxZ;
				m_dataContraMetricXiREdge[k][iA][iB][2] =
					(1.0 + dDaZ * dDaZ + dDbZ * dDbZ) / (dDxZ * dDxZ);

				// Derivatives of the vertical coordinate transform
				m_dataDerivRREdge[k][iA][iB][0] = dDaZ;
				m_dataDerivRREdge[k][iA][iB][1] = dDbZ;
				m_dataDerivRREdge[k][iA][iB][2] = dDxZ;
			}
		}
		}
	}
	}
}

///////////////////////////////////////////////////////////////////////////////

void GridPatchCartesianGLL::EvaluateTestCase(
	const TestCase & test,
	const Time & time,
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

	// Evaluate topography
	EvaluateTopography(test);

	// Physical constants
	const PhysicalConstants & phys = m_grid.GetModel().GetPhysicalConstants();

	GridCartesianGLL & gridCartesianGLL =
		dynamic_cast<GridCartesianGLL &>(m_grid);

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

/*
		// 2nd order decay vertical coordinate
		for (int k = 0; k < m_grid.GetRElements(); k++) {
			m_dataZLevels[k][i][j] =
				m_grid.GetREtaLevel(k) * m_grid.GetZtop()
				+ m_dataTopography[i][j]
				* std::pow(1.0 - m_grid.GetREtaLevel(k), 2.0);
		}
		for (int k = 0; k <= m_grid.GetRElements(); k++) {
			m_dataZInterfaces[k][i][j] =
				m_grid.GetREtaInterface(k) * m_grid.GetZtop()
				+ m_dataTopography[i][j]
				* std::pow(1.0 - m_grid.GetREtaInterface(k), 2.0);
		}
*/
/*
		// 4th order decay vertical coordinate
		for (int k = 0; k < m_grid.GetRElements(); k++) {
			m_dataZLevels[k][i][j] =
				m_grid.GetREtaLevel(k) * m_grid.GetZtop()
				+ m_dataTopography[i][j]
				* std::pow(1.0 - m_grid.GetREtaLevel(k), 4.0);
		}
		for (int k = 0; k <= m_grid.GetRElements(); k++) {
			m_dataZInterfaces[k][i][j] =
				m_grid.GetREtaInterface(k) * m_grid.GetZtop()
				+ m_dataTopography[i][j]
				* std::pow(1.0 - m_grid.GetREtaInterface(k), 4.0);
		}
*/
/*
		// 4th order cosine decay vertical coordinate
		for (int k = 0; k < m_grid.GetRElements(); k++) {
			m_dataZLevels[k][i][j] =
				m_grid.GetREtaLevel(k) * m_grid.GetZtop()
				+ m_dataTopography[i][j]
				* std::pow(std::cos(0.5 * M_PI * m_grid.GetREtaLevel(k)), 4.0);
		}
		for (int k = 0; k <= m_grid.GetRElements(); k++) {
			m_dataZInterfaces[k][i][j] =
				m_grid.GetREtaInterface(k) * m_grid.GetZtop()
				+ m_dataTopography[i][j]
				* std::pow(std::cos(0.5 * M_PI * m_grid.GetREtaInterface(k)), 4.0);
		}
*/
/*
		// 8th order cosine decay vertical coordinate
		for (int k = 0; k < m_grid.GetRElements(); k++) {
			m_dataZLevels[k][i][j] =
				m_grid.GetREtaLevel(k) * m_grid.GetZtop()
				+ m_dataTopography[i][j]
				* std::pow(std::cos(0.5 * M_PI * m_grid.GetREtaLevel(k)), 8.0);
		}
		for (int k = 0; k <= m_grid.GetRElements(); k++) {
			m_dataZInterfaces[k][i][j] =
				m_grid.GetREtaInterface(k) * m_grid.GetZtop()
				+ m_dataTopography[i][j]
				* std::pow(std::cos(0.5 * M_PI * m_grid.GetREtaInterface(k)), 8.0);
		}
*/
/*
		// Schar Exponential Decay vertical coordinate
		for (int k = 0; k < m_grid.GetRElements(); k++) {
			m_dataZLevels[k][i][j] = m_grid.GetZtop() * m_grid.GetREtaLevel(k) + 
			m_dataTopography[i][j] * sinh(m_grid.GetZtop() * 
			(1.0 - m_grid.GetREtaLevel(k)) / gridCartesianGLL.m_dSL) / 
			sinh(m_grid.GetZtop() / gridCartesianGLL.m_dSL);
		}
		for (int k = 0; k <= m_grid.GetRElements(); k++) {
			m_dataZInterfaces[k][i][j] = m_grid.GetZtop() * m_grid.GetREtaInterface(k) + 
			m_dataTopography[i][j] * sinh(m_grid.GetZtop() * 
			(1.0 - m_grid.GetREtaInterface(k)) / gridCartesianGLL.m_dSL) / 
			sinh(m_grid.GetZtop() / gridCartesianGLL.m_dSL);
		}
*/
/*
		// 10th order cosine decay with finite bottom slope vertical coordinate
		double power = 10.0;
		double botRate = 1.0;
		for (int k = 0; k < m_grid.GetRElements(); k++) {
			m_dataZLevels[k][i][j] =
				m_grid.GetREtaLevel(k) * m_grid.GetZtop()
				+ m_dataTopography[i][j] 
				* (1.0 - botRate * m_grid.GetREtaLevel(k))
				* std::pow(std::cos(0.5 * M_PI * m_grid.GetREtaLevel(k)), power);
		}
		for (int k = 0; k <= m_grid.GetRElements(); k++) {
			m_dataZInterfaces[k][i][j] =
				m_grid.GetREtaInterface(k) * m_grid.GetZtop()
				+ m_dataTopography[i][j] 
				* (1.0 - botRate * m_grid.GetREtaInterface(k))
				* std::pow(std::cos(0.5 * M_PI * m_grid.GetREtaInterface(k)), power);
		}
*/
	}
	}

	// Initialize the Rayleigh friction strength at each node IN COMPUTATIONAL GRID
	if (test.HasRayleighFriction()) {
		for (int i = 0; i < m_box.GetATotalWidth(); i++) {
		for (int j = 0; j < m_box.GetBTotalWidth(); j++) {
			for (int k = 0; k < m_grid.GetRElements(); k++) {
				m_dataRayleighStrengthNode[k][i][j] =
					test.EvaluateRayleighStrength(
						m_dataZLevels[k][i][j],
						m_dataLon[i][j],
						m_dataLat[i][j]);
			}
			for (int k = 0; k < m_grid.GetRElements(); k++) {
				m_dataRayleighStrengthREdge[k][i][j] =
					test.EvaluateRayleighStrength(
						m_dataZInterfaces[k][i][j],
						m_dataLon[i][j],
						m_dataLat[i][j]);
			}
		}
		}
	}

	// Buffer vector for storing pointwise states
	const EquationSet & eqns = m_grid.GetModel().GetEquationSet();

	int nComponents = eqns.GetComponents();
	int nTracers = eqns.GetTracers();

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
	for (int k = 0; k < m_grid.GetRElements(); k++) {
	for (int i = 0; i < m_box.GetATotalWidth(); i++) {
	for (int j = 0; j < m_box.GetBTotalWidth(); j++) {

		// Evaluate pointwise state
		test.EvaluatePointwiseState(
			m_grid.GetModel().GetPhysicalConstants(),
			time,
			m_dataZLevels[k][i][j],
			m_dataLon[i][j],
			m_dataLat[i][j],
			dPointwiseState,
			dPointwiseTracers);

		eqns.ConvertComponents(
			phys, dPointwiseState, dPointwiseTracers);

		for (int c = 0; c < dPointwiseState.GetRows(); c++) {
			m_datavecStateNode[iDataIndex][c][k][i][j] = dPointwiseState[c];
		}

		// Evaluate reference state
		if (m_grid.HasReferenceState()) {
			test.EvaluateReferenceState(
				m_grid.GetModel().GetPhysicalConstants(),
				m_dataZLevels[k][i][j],
				m_dataLon[i][j],
				m_dataLat[i][j],
				dPointwiseRefState);

			eqns.ConvertComponents(
				phys, dPointwiseRefState, dPointwiseRefTracers);

			for (int c = 0; c < dPointwiseState.GetRows(); c++) {
				m_dataRefStateNode[c][k][i][j] = dPointwiseRefState[c];
			}

#if defined(TRACER_REFERENCE_STATE)
			for (int c = 0; c < dPointwiseRefTracers.GetRows(); c++) {
				m_dataRefTracers[c][k][i][j] =
					dPointwiseRefTracers[c];
			}
#endif
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
			phys,
			time,
			m_dataZInterfaces[k][i][j],
			m_dataLon[i][j],
			m_dataLat[i][j],
			dPointwiseState,
			dPointwiseTracers);

		eqns.ConvertComponents(phys, dPointwiseState, dPointwiseTracers);

		for (int c = 0; c < dPointwiseState.GetRows(); c++) {
			m_datavecStateREdge[iDataIndex][c][k][i][j] = dPointwiseState[c];
		}

		if (m_grid.HasReferenceState()) {
			test.EvaluateReferenceState(
				phys,
				m_dataZInterfaces[k][i][j],
				m_dataLon[i][j],
				m_dataLat[i][j],
				dPointwiseRefState);

			DataArray1D<double> dPointwiseRefTracers;
			eqns.ConvertComponents(
				phys, dPointwiseRefState, dPointwiseRefTracers);

			for (int c = 0; c < dPointwiseState.GetRows(); c++) {
				m_dataRefStateREdge[c][k][i][j] = dPointwiseRefState[c];
			}
		}
	}
	}
	}
}

///////////////////////////////////////////////////////////////////////////////

void GridPatchCartesianGLL::ApplyBoundaryConditions(
	int iDataIndex,
	DataType eDataType,
	int iAlphaBCPatch
) {
	// Indices of EquationSet variables
	const int UIx = 0;
	const int VIx = 1;
	const int TIx = 2;
	const int WIx = 3;
	const int RIx = 4;

	double ua_hat = 0.0;
	double ub_hat = 0.0;
	double ux_hat = 0.0;
	double gaa = 0.0; double gab = 0.0; double gax = 0.0;
	double gbb = 0.0; double gba = 0.0; double gbx = 0.0;

//std::cout << m_box.GetAGlobalBegin() << "  " << m_box.GetAGlobalEnd() << std::endl;

	// Impose boundary conditions (everything on levels)
	if (m_grid.GetVerticalStaggering() ==
		Grid::VerticalStaggering_Levels
	) {
		_EXCEPTIONT("Not implemented: BC for all variables on levels!");

	// Impose boundary conditions (everything on interfaces)
	} else if (m_grid.GetVerticalStaggering() ==
		Grid::VerticalStaggering_Interfaces
	) {
		// Impose boundary conditions along right edge
		Grid::BoundaryCondition eBoundaryRight =
			m_grid.GetBoundaryCondition(Direction_Right);
		if (eBoundaryRight != Grid::BoundaryCondition_Periodic) {
			int i = m_box.GetATotalWidth() - 1;

			for (int k = 0; k < m_grid.GetRElements(); k++) {
			for (int j = 0; j < m_box.GetBTotalWidth(); j++) {
				if (eBoundaryRight != Grid::BoundaryCondition_NoSlip) {
					m_datavecStateNode[iDataIndex][UIx][k][i][j] =
						- m_datavecStateNode[iDataIndex][UIx][k][i-1][j];

					m_datavecStateNode[iDataIndex][VIx][k][i][j] =
						- m_datavecStateNode[iDataIndex][VIx][k][i-1][j];

					m_datavecStateNode[iDataIndex][WIx][k][i][j] =
						- m_datavecStateNode[iDataIndex][WIx][k][i-1][j];
				} else if (eBoundaryRight != Grid::BoundaryCondition_NoFlux) {
					// DSS the local boundary u_beta and u_xi
					ub_hat = 0.5 * 
						m_datavecStateNode[iDataIndex][VIx][k][i][j] +
						m_datavecStateNode[iDataIndex][VIx][k][i-1][j];
					ux_hat = 0.5 * 
						m_datavecStateNode[iDataIndex][WIx][k][i][j] +
						m_datavecStateNode[iDataIndex][WIx][k][i-1][j];
					ux_hat *= m_dataDerivRNode[k][i-1][j][2];
					// Get the local alpha contravariant metric components
					gaa = m_dataContraMetricA[k][i-1][j][0];
					gab = m_dataContraMetricA[k][i-1][j][1];
					gax = m_dataContraMetricA[k][i-1][j][2];
					// Compute the convariant boundary velocity u^hat
					ua_hat = -1.0 / gaa * (gab * ub_hat + gax * ux_hat);
					// Compute the pre-DSS halo velocity that imposes no-flux
					m_datavecStateNode[iDataIndex][UIx][k][i][j] =
						2.0 * ua_hat
						- m_datavecStateNode[iDataIndex][UIx][k][i-1][j];
				} else {
					_EXCEPTIONT("+X NoSlip or NoFlux only: ApplyBoundaryConditions");
				}
			}
			}
		}

		// Impose boundary conditions along top edge
		Grid::BoundaryCondition eBoundaryTop =
			m_grid.GetBoundaryCondition(Direction_Top);
		if (eBoundaryTop != Grid::BoundaryCondition_Periodic) {
			int j = m_box.GetBTotalWidth() - 1;

			for (int k = 0; k < m_grid.GetRElements(); k++) {
			for (int i = 1; i < m_box.GetATotalWidth()-1; i++) {
				if (eBoundaryTop != Grid::BoundaryCondition_NoSlip) {
					m_datavecStateNode[iDataIndex][UIx][k][i][j] =
						- m_datavecStateNode[iDataIndex][UIx][k][i][j-1];

					m_datavecStateNode[iDataIndex][VIx][k][i][j] =
						- m_datavecStateNode[iDataIndex][VIx][k][i][j-1];

					m_datavecStateNode[iDataIndex][WIx][k][i][j] =
						- m_datavecStateNode[iDataIndex][WIx][k][i][j-1];
				} else if (eBoundaryTop != Grid::BoundaryCondition_NoFlux) {
					// DSS the local boundary u_alpha and u_xi
					ua_hat = 0.5 * 
						m_datavecStateNode[iDataIndex][UIx][k][i][j] +
						m_datavecStateNode[iDataIndex][UIx][k][i][j-1];
					ux_hat = 0.5 * 
						m_datavecStateNode[iDataIndex][WIx][k][i][j] +
						m_datavecStateNode[iDataIndex][WIx][k][i][j-1];
					ux_hat *= m_dataDerivRNode[k][i][j-1][2];
					// Get the local beta contravariant metric components
					gba = m_dataContraMetricB[k][i][j-1][0];
					gbb = m_dataContraMetricB[k][i][j-1][1];
					gbx = m_dataContraMetricB[k][i][j-1][2];
					//std::cout << gbb << "\n";
					// Compute the convariant boundary velocity u^hat
					ub_hat = -1.0 / gbb * (gba * ua_hat + gbx * ux_hat);
					// Compute the pre-DSS halo velocity that imposes no-flux
					m_datavecStateNode[iDataIndex][VIx][k][i][j] =
						2.0 * ub_hat
						- m_datavecStateNode[iDataIndex][VIx][k][i][j-1];
				} else {
					_EXCEPTIONT("+Y NoSlip or NoFlux only: ApplyBoundaryConditions");
				}
			}
			}
		}

		// Impose boundary conditions along left edge
		Grid::BoundaryCondition eBoundaryLeft =
			m_grid.GetBoundaryCondition(Direction_Left);
		if (eBoundaryLeft != Grid::BoundaryCondition_Periodic) {
			int i = 0;

			for (int k = 0; k < m_grid.GetRElements(); k++) {
			for (int j = 1; j < m_box.GetBTotalWidth()-1; j++) {
				if (eBoundaryLeft == Grid::BoundaryCondition_NoSlip) {
					m_datavecStateNode[iDataIndex][UIx][k][i][j] =
						- m_datavecStateNode[iDataIndex][UIx][k][i+1][j];

					m_datavecStateNode[iDataIndex][VIx][k][i][j] =
						- m_datavecStateNode[iDataIndex][VIx][k][i+1][j];

					m_datavecStateNode[iDataIndex][WIx][k][i][j] =
						- m_datavecStateNode[iDataIndex][WIx][k][i+1][j];
				} else if (eBoundaryLeft == Grid::BoundaryCondition_NoFlux) {
					// DSS the local boundary u_beta and u_xi				
					ub_hat = 0.5 * 
						m_datavecStateNode[iDataIndex][VIx][k][i][j] +
						m_datavecStateNode[iDataIndex][VIx][k][i+1][j];
					ux_hat = 0.5 * 
						m_datavecStateNode[iDataIndex][WIx][k][i][j] +
						m_datavecStateNode[iDataIndex][WIx][k][i+1][j];
					ux_hat *= m_dataDerivRNode[k][i+1][j][2];
					// Get the local alpha contravariant metric components
					gaa = m_dataContraMetricA[k][i+1][j][0];
					gab = m_dataContraMetricA[k][i+1][j][1];
					gax = m_dataContraMetricA[k][i+1][j][2];
					// Compute the convariant boundary velocity u^hat
					ua_hat = -1.0 / gaa * (gab * ub_hat + gax * ux_hat);
					// Compute the pre-DSS halo velocity that imposes no-flux
					m_datavecStateNode[iDataIndex][UIx][k][i][j] =
						2.0 * ua_hat
						- m_datavecStateNode[iDataIndex][UIx][k][i+1][j];
				} else {
					_EXCEPTIONT("-X NoSlip or NoFlux only: ApplyBoundaryConditions");
				}
			}
			}
		}

		// Impose boundary conditions along bottom edge
		Grid::BoundaryCondition eBoundaryBottom =
			m_grid.GetBoundaryCondition(Direction_Bottom);
		if (eBoundaryBottom != Grid::BoundaryCondition_Periodic) {
			int j = 0;

			for (int k = 0; k < m_grid.GetRElements(); k++) {
			for (int i = 1; i < m_box.GetATotalWidth()-1; i++) {
				if (eBoundaryBottom == Grid::BoundaryCondition_NoSlip) {
					m_datavecStateNode[iDataIndex][UIx][k][i][j] =
						- m_datavecStateNode[iDataIndex][UIx][k][i][j+1];

					m_datavecStateNode[iDataIndex][VIx][k][i][j] =
						- m_datavecStateNode[iDataIndex][VIx][k][i][j+1];

					m_datavecStateNode[iDataIndex][WIx][k][i][j] =
						- m_datavecStateNode[iDataIndex][WIx][k][i][j+1];
				} else if (eBoundaryBottom == Grid::BoundaryCondition_NoFlux) {
					// DSS the local boundary u_alpha and u_xi				
					ua_hat = 0.5 * 
						m_datavecStateNode[iDataIndex][UIx][k][i][j] +
						m_datavecStateNode[iDataIndex][UIx][k][i][j+1];
					ux_hat = 0.5 * 
						m_datavecStateNode[iDataIndex][WIx][k][i][j] +
						m_datavecStateNode[iDataIndex][WIx][k][i][j+1];
					ux_hat *= m_dataDerivRNode[k][i][j+1][2];
					// Get the local beta contravariant metric components
					gba = m_dataContraMetricB[k][i][j+1][0];
					gbb = m_dataContraMetricB[k][i][j+1][1];
					gbx = m_dataContraMetricB[k][i][j+1][2];
					// Compute the convariant boundary velocity u^hat
					ub_hat = -1.0 / gbb * (gba * ua_hat + gbx * ux_hat);
					// Compute the pre-DSS halo velocity that imposes no-flux
					m_datavecStateNode[iDataIndex][VIx][k][i][j] =
						2.0 * ub_hat
						- m_datavecStateNode[iDataIndex][VIx][k][i][j+1];
				} else {
					_EXCEPTIONT("-Y NoSlip or NoFlux only: ApplyBoundaryConditions");
				}
			}
			}
		}

	// Impose boundary conditions (vertical velocity on interfaces)
	} else {
		// Impose boundary conditions along right (alpha+) edge
		Grid::BoundaryCondition eBoundaryRight =
			m_grid.GetBoundaryCondition(Direction_Right);
		if (eBoundaryRight != Grid::BoundaryCondition_Periodic) {
			int i = m_box.GetATotalWidth() - 1;
			
			for (int j = 1; j < m_box.GetBTotalWidth()-1; j++) {
			for (int k = 0; k < m_grid.GetRElements(); k++) {
				if (eBoundaryRight == Grid::BoundaryCondition_NoSlip) {
					m_datavecStateNode[iDataIndex][UIx][k][i][j] =
						- m_datavecStateNode[iDataIndex][UIx][k][i-1][j];

					m_datavecStateNode[iDataIndex][VIx][k][i][j] =
						- m_datavecStateNode[iDataIndex][VIx][k][i-1][j];

					m_datavecStateNode[iDataIndex][WIx][k][i][j] =
						- m_datavecStateNode[iDataIndex][WIx][k][i-1][j];
				} else if (eBoundaryRight == Grid::BoundaryCondition_NoFlux) {
					// DSS the local boundary u_beta and u_xi				
					ub_hat = 0.5 * 
						m_datavecStateNode[iDataIndex][VIx][k][i][j] +
						m_datavecStateNode[iDataIndex][VIx][k][i-1][j];
					ux_hat = 0.5 * 
						m_datavecStateNode[iDataIndex][WIx][k][i][j] +
						m_datavecStateNode[iDataIndex][WIx][k][i-1][j];
					ux_hat *= m_dataDerivRNode[k][i-1][j][2];
					// Get the local alpha contravariant metric components
					gaa = m_dataContraMetricA[k][i-1][j][0];
					gab = m_dataContraMetricA[k][i-1][j][1];
					gax = m_dataContraMetricA[k][i-1][j][2];
					// Compute the convariant boundary velocity u^hat
					ua_hat = -1.0 / gaa * (gab * ub_hat + gax * ux_hat);
					// Compute the pre-DSS halo velocity that imposes no-flux
					m_datavecStateNode[iDataIndex][UIx][k][i][j] =
						2.0 * ua_hat
						- m_datavecStateNode[iDataIndex][UIx][k][i-1][j];
				} else {
					_EXCEPTIONT("+X NoSlip or NoFlux only: ApplyBoundaryConditions");
				}
			}

			for (int k = 0; k <= m_grid.GetRElements(); k++) {
				if (eBoundaryRight == Grid::BoundaryCondition_NoSlip) {
					m_datavecStateREdge[iDataIndex][UIx][k][i][j] =
						- m_datavecStateREdge[iDataIndex][UIx][k][i-1][j];

					m_datavecStateREdge[iDataIndex][VIx][k][i][j] =
						- m_datavecStateREdge[iDataIndex][VIx][k][i-1][j];

					m_datavecStateREdge[iDataIndex][WIx][k][i][j] =
						- m_datavecStateREdge[iDataIndex][WIx][k][i-1][j];
				} else if (eBoundaryRight == Grid::BoundaryCondition_NoFlux) {
					// DSS the local boundary u_beta and u_xi				
					ub_hat = 0.5 * 
						m_datavecStateREdge[iDataIndex][VIx][k][i][j] +
						m_datavecStateREdge[iDataIndex][VIx][k][i-1][j];
					ux_hat = 0.5 * 
						m_datavecStateREdge[iDataIndex][WIx][k][i][j] +
						m_datavecStateREdge[iDataIndex][WIx][k][i-1][j];
					ux_hat *= m_dataDerivRREdge[k][i-1][j][2];
					// Get the local alpha contravariant metric components
					gaa = m_dataContraMetricAREdge[k][i-1][j][0];
					gab = m_dataContraMetricAREdge[k][i-1][j][1];
					gax = m_dataContraMetricAREdge[k][i-1][j][2];
					// Compute the convariant boundary velocity u^hat
					ua_hat = -1.0 / gaa * (gab * ub_hat + gax * ux_hat);
					// Compute the pre-DSS halo velocity that imposes no-flux
					m_datavecStateREdge[iDataIndex][UIx][k][i][j] =
						2.0 * ua_hat
						- m_datavecStateREdge[iDataIndex][UIx][k][i-1][j];
				} else {
					_EXCEPTIONT("+X NoSlip or NoFlux only: ApplyBoundaryConditions");
				}
			}
			}
		}

		// Impose boundary conditions along top (beta+) edge
		Grid::BoundaryCondition eBoundaryTop =
			m_grid.GetBoundaryCondition(Direction_Top);
		if (eBoundaryTop != Grid::BoundaryCondition_Periodic) {
			int j = m_box.GetBTotalWidth() - 1;

			for (int i = 1; i < m_box.GetATotalWidth()-1; i++) {
			for (int k = 0; k < m_grid.GetRElements(); k++) {
				if (eBoundaryTop == Grid::BoundaryCondition_NoSlip) {
					m_datavecStateNode[iDataIndex][UIx][k][i][j] =
						- m_datavecStateNode[iDataIndex][UIx][k][i][j-1];

					m_datavecStateNode[iDataIndex][VIx][k][i][j] =
						- m_datavecStateNode[iDataIndex][VIx][k][i][j-1];

					m_datavecStateNode[iDataIndex][WIx][k][i][j] =
						- m_datavecStateNode[iDataIndex][WIx][k][i][j-1];
				} else if (eBoundaryTop == Grid::BoundaryCondition_NoFlux) {
					// DSS the local boundary u_alpha and u_xi				
					ua_hat = 0.5 * 
						m_datavecStateNode[iDataIndex][UIx][k][i][j] +
						m_datavecStateNode[iDataIndex][UIx][k][i][j-1];
					ux_hat = 0.5 * 
						m_datavecStateNode[iDataIndex][WIx][k][i][j] +
						m_datavecStateNode[iDataIndex][WIx][k][i][j-1];
					ux_hat *= m_dataDerivRNode[k][i][j-1][2];
					// Get the local beta contravariant metric components
					gba = m_dataContraMetricB[k][i][j-1][0];
					gbb = m_dataContraMetricB[k][i][j-1][1];
					gbx = m_dataContraMetricB[k][i][j-1][2];
					//std::cout << gbb << "\n";
					// Compute the convariant boundary velocity u^hat
					ub_hat = -1.0 / gbb * (gba * ua_hat + gbx * ux_hat);
					// Compute the pre-DSS halo velocity that imposes no-flux
					m_datavecStateNode[iDataIndex][VIx][k][i][j] =
						2.0 * ub_hat
						- m_datavecStateNode[iDataIndex][VIx][k][i][j-1];
				} else {
					_EXCEPTIONT("+Y NoSlip or NoFlux only: ApplyBoundaryConditions");
				}
			}

			for (int k = 0; k <= m_grid.GetRElements(); k++) {
				if (eBoundaryTop == Grid::BoundaryCondition_NoSlip) {
					m_datavecStateREdge[iDataIndex][UIx][k][i][j] =
						- m_datavecStateREdge[iDataIndex][UIx][k][i][j-1];

					m_datavecStateREdge[iDataIndex][VIx][k][i][j] =
						- m_datavecStateREdge[iDataIndex][VIx][k][i][j-1];

					m_datavecStateREdge[iDataIndex][WIx][k][i][j] =
						- m_datavecStateREdge[iDataIndex][WIx][k][i][j-1];
				} else if (eBoundaryTop == Grid::BoundaryCondition_NoFlux) {
					// DSS the local boundary u_alpha and u_xi				
					ua_hat = 0.5 * 
						m_datavecStateREdge[iDataIndex][UIx][k][i][j] +
						m_datavecStateREdge[iDataIndex][UIx][k][i][j-1];
					ux_hat = 0.5 * 
						m_datavecStateREdge[iDataIndex][WIx][k][i][j] +
						m_datavecStateREdge[iDataIndex][WIx][k][i][j-1];
					ux_hat *= m_dataDerivRREdge[k][i][j-1][2];
					// Get the local beta contravariant metric components
					gba = m_dataContraMetricBREdge[k][i][j-1][0];
					gbb = m_dataContraMetricBREdge[k][i][j-1][1];
					gbx = m_dataContraMetricBREdge[k][i][j-1][2];
					// Compute the convariant boundary velocity u^hat
					ub_hat = -1.0 / gbb * (gba * ua_hat + gbx * ux_hat);
					// Compute the pre-DSS halo velocity that imposes no-flux
					m_datavecStateREdge[iDataIndex][VIx][k][i][j] =
						2.0 * ub_hat
						- m_datavecStateREdge[iDataIndex][VIx][k][i][j-1];
				} else {
					_EXCEPTIONT("+Y NoSlip or NoFlux only: ApplyBoundaryConditions");
				}
			}
			}
		}

		// Impose boundary conditions along left (alpha-) edge
		Grid::BoundaryCondition eBoundaryLeft =
			m_grid.GetBoundaryCondition(Direction_Left);
		if (eBoundaryLeft != Grid::BoundaryCondition_Periodic) {
			int i = 0;
			
			for (int j = 1; j < m_box.GetBTotalWidth()-1; j++) {
			for (int k = 0; k < m_grid.GetRElements(); k++) {
				if (eBoundaryLeft == Grid::BoundaryCondition_NoSlip) {
					m_datavecStateNode[iDataIndex][UIx][k][i][j] =
						- m_datavecStateNode[iDataIndex][UIx][k][i+1][j];

					m_datavecStateNode[iDataIndex][VIx][k][i][j] =
						- m_datavecStateNode[iDataIndex][VIx][k][i+1][j];

					m_datavecStateNode[iDataIndex][WIx][k][i][j] =
						- m_datavecStateNode[iDataIndex][WIx][k][i+1][j];
				} else if (eBoundaryLeft == Grid::BoundaryCondition_NoFlux) {
					// DSS the local boundary u_beta and u_xi				
					ub_hat = 0.5 * 
						m_datavecStateNode[iDataIndex][VIx][k][i][j] +
						m_datavecStateNode[iDataIndex][VIx][k][i+1][j];
					ux_hat = 0.5 * 
						m_datavecStateNode[iDataIndex][WIx][k][i][j] +
						m_datavecStateNode[iDataIndex][WIx][k][i+1][j];
					ux_hat *= m_dataDerivRNode[k][i+1][j][2];
					// Get the local alpha contravariant metric components
					gaa = m_dataContraMetricA[k][i+1][j][0];
					gab = m_dataContraMetricA[k][i+1][j][1];
					gax = m_dataContraMetricA[k][i+1][j][2];
					// Compute the convariant boundary velocity u^hat
					ua_hat = -1.0 / gaa * (gab * ub_hat + gax * ux_hat);
					// Compute the pre-DSS halo velocity that imposes no-flux
					m_datavecStateNode[iDataIndex][UIx][k][i][j] =
						2.0 * ua_hat
						- m_datavecStateNode[iDataIndex][UIx][k][i+1][j];
				} else {
					_EXCEPTIONT("-X NoSlip or NoFlux only: ApplyBoundaryConditions");
				}
			}

			for (int k = 0; k <= m_grid.GetRElements(); k++) {
				if (eBoundaryLeft == Grid::BoundaryCondition_NoSlip) {
					m_datavecStateREdge[iDataIndex][UIx][k][i][j] =
						- m_datavecStateREdge[iDataIndex][UIx][k][i+1][j];

					m_datavecStateREdge[iDataIndex][VIx][k][i][j] =
						- m_datavecStateREdge[iDataIndex][VIx][k][i+1][j];

					m_datavecStateREdge[iDataIndex][WIx][k][i][j] =
						- m_datavecStateREdge[iDataIndex][WIx][k][i+1][j];
				} else if (eBoundaryLeft == Grid::BoundaryCondition_NoFlux) {
					// DSS the local boundary u_beta and u_xi				
					ub_hat = 0.5 * 
						m_datavecStateREdge[iDataIndex][VIx][k][i][j] +
						m_datavecStateREdge[iDataIndex][VIx][k][i+1][j];
					ux_hat = 0.5 * 
						m_datavecStateREdge[iDataIndex][WIx][k][i][j] +
						m_datavecStateREdge[iDataIndex][WIx][k][i+1][j];
					ux_hat *= m_dataDerivRREdge[k][i+1][j][2];
					// Get the local alpha contravariant metric components
					gaa = m_dataContraMetricAREdge[k][i+1][j][0];
					gab = m_dataContraMetricAREdge[k][i+1][j][1];
					gax = m_dataContraMetricAREdge[k][i+1][j][2];
					// Compute the convariant boundary velocity u^hat
					ua_hat = -1.0 / gaa * (gab * ub_hat + gax * ux_hat);
					// Compute the pre-DSS halo velocity that imposes no-flux
					m_datavecStateREdge[iDataIndex][UIx][k][i][j] =
						2.0 * ua_hat
						- m_datavecStateREdge[iDataIndex][UIx][k][i+1][j];
				} else {
					_EXCEPTIONT("-X NoSlip or NoFlux only: ApplyBoundaryConditions");
				}
			}
			}
		}

		// Impose boundary conditions along bottom (beta-) edge
		Grid::BoundaryCondition eBoundaryBottom =
			m_grid.GetBoundaryCondition(Direction_Bottom);
		if (eBoundaryBottom != Grid::BoundaryCondition_Periodic) {
			int j = 0;
			
			for (int i = 1; i < m_box.GetATotalWidth()-1; i++) {
			for (int k = 0; k < m_grid.GetRElements(); k++) {
				if (eBoundaryBottom == Grid::BoundaryCondition_NoSlip) {
					m_datavecStateNode[iDataIndex][UIx][k][i][j] =
						- m_datavecStateNode[iDataIndex][UIx][k][i][j+1];

					m_datavecStateNode[iDataIndex][VIx][k][i][j] =
						- m_datavecStateNode[iDataIndex][VIx][k][i][j+1];

					m_datavecStateNode[iDataIndex][WIx][k][i][j] =
						- m_datavecStateNode[iDataIndex][WIx][k][i][j+1];
				} else if (eBoundaryBottom == Grid::BoundaryCondition_NoFlux) {
					// DSS the local boundary u_alpha and u_xi				
					ua_hat = 0.5 * 
						m_datavecStateNode[iDataIndex][UIx][k][i][j] +
						m_datavecStateNode[iDataIndex][UIx][k][i][j+1];
					ux_hat = 0.5 * 
						m_datavecStateNode[iDataIndex][WIx][k][i][j] +
						m_datavecStateNode[iDataIndex][WIx][k][i][j+1];
					ux_hat *= m_dataDerivRNode[k][i][j+1][2];
					// Get the local beta contravariant metric components
					gba = m_dataContraMetricB[k][i][j+1][0];
					gbb = m_dataContraMetricB[k][i][j+1][1];
					gbx = m_dataContraMetricB[k][i][j+1][2];
					// Compute the convariant boundary velocity u^hat
					ub_hat = -1.0 / gbb * (gba * ua_hat + gbx * ux_hat);
					// Compute the pre-DSS halo velocity that imposes no-flux
					m_datavecStateNode[iDataIndex][VIx][k][i][j] =
						2.0 * ub_hat
						- m_datavecStateNode[iDataIndex][VIx][k][i][j+1];
				} else {
					_EXCEPTIONT("-Y NoSlip or NoFlux only: ApplyBoundaryConditions");
				}
			}

			for (int k = 0; k <= m_grid.GetRElements(); k++) {
				if (eBoundaryBottom == Grid::BoundaryCondition_NoSlip) {
					m_datavecStateREdge[iDataIndex][UIx][k][i][j] =
						- m_datavecStateREdge[iDataIndex][UIx][k][i][j+1];

					m_datavecStateREdge[iDataIndex][VIx][k][i][j] =
						- m_datavecStateREdge[iDataIndex][VIx][k][i][j+1];

					m_datavecStateREdge[iDataIndex][WIx][k][i][j] =
						- m_datavecStateREdge[iDataIndex][WIx][k][i][j+1];
				} else if (eBoundaryBottom == Grid::BoundaryCondition_NoFlux) {
					// DSS the local boundary u_alpha and u_xi				
					ua_hat = 0.5 * 
						m_datavecStateREdge[iDataIndex][UIx][k][i][j] +
						m_datavecStateREdge[iDataIndex][UIx][k][i][j+1];
					ux_hat = 0.5 * 
						m_datavecStateREdge[iDataIndex][WIx][k][i][j] +
						m_datavecStateREdge[iDataIndex][WIx][k][i][j+1];
					ux_hat *= m_dataDerivRREdge[k][i][j+1][2];
					// Get the local beta contravariant metric components
					gba = m_dataContraMetricBREdge[k][i][j+1][0];
					gbb = m_dataContraMetricBREdge[k][i][j+1][1];
					gbx = m_dataContraMetricBREdge[k][i][j+1][2];
					// Compute the convariant boundary velocity u^hat
					ub_hat = -1.0 / gbb * (gba * ua_hat + gbx * ux_hat);
					// Compute the pre-DSS halo velocity that imposes no-flux
					m_datavecStateREdge[iDataIndex][VIx][k][i][j] =
						2.0 * ub_hat
						- m_datavecStateREdge[iDataIndex][VIx][k][i][j+1];
				} else {
					_EXCEPTIONT("-Y NoSlip or NoFlux only: ApplyBoundaryConditions");
				}
			}
			}
		}
#pragma message "Interpolate boundary conditions for W back to interfaces somehow"
	}
}

///////////////////////////////////////////////////////////////////////////////

void GridPatchCartesianGLL::ComputeCurlAndDiv(
	const DataArray3D<double> & dataUa,
	const DataArray3D<double> & dataUb
) {
	// Parent grid
	const GridCartesianGLL & gridCSGLL =
		dynamic_cast<const GridCartesianGLL &>(m_grid);

	// Compute derivatives of the field
	const DataArray2D<double> & dDxBasis1D = gridCSGLL.GetDxBasis1D();

	// Number of finite elements in each direction
	int nAFiniteElements = m_box.GetAInteriorWidth() / m_nHorizontalOrder;
	int nBFiniteElements = m_box.GetBInteriorWidth() / m_nHorizontalOrder;

	// Contravariant velocity within an element
	DataArray2D<double> dConUa(m_nHorizontalOrder, m_nHorizontalOrder);
	DataArray2D<double> dConUb(m_nHorizontalOrder, m_nHorizontalOrder);

	// Loop over all elements in the box
	for (int k = 0; k < gridCSGLL.GetRElements(); k++) {
	for (int a = 0; a < nAFiniteElements; a++) {
	for (int b = 0; b < nBFiniteElements; b++) {

		// Index of lower-left corner node
		int iA = a * m_nHorizontalOrder + m_box.GetHaloElements();
		int iB = b * m_nHorizontalOrder + m_box.GetHaloElements();

		// Calculate contravariant velocity at each node within the element
		for (int i = 0; i < m_nHorizontalOrder; i++) {
		for (int j = 0; j < m_nHorizontalOrder; j++) {
			dConUa[i][j] =
				  m_dataContraMetric2DA[iA+i][iB+j][0]
				  	* dataUa[k][iA+i][iB+j]
				+ m_dataContraMetric2DA[iA+i][iB+j][1]
					* dataUb[k][iA+i][iB+j];

			dConUb[i][j] =
				  m_dataContraMetric2DB[iA+i][iB+j][0]
				  	* dataUa[k][iA+i][iB+j]
				+ m_dataContraMetric2DB[iA+i][iB+j][1]
					* dataUb[k][iA+i][iB+j];
		}
		}

		// Calculate divergance and curl
		for (int i = 0; i < m_nHorizontalOrder; i++) {
		for (int j = 0; j < m_nHorizontalOrder; j++) {

			// Compute derivatives at each node
			double dDaJUa = 0.0;
			double dDbJUb = 0.0;

			double dCovDaUb = 0.0;
			double dCovDbUa = 0.0;

			for (int s = 0; s < m_nHorizontalOrder; s++) {
				dDaJUa += dDxBasis1D[s][i]
					* m_dataJacobian2D[iA+s][iB+j]
					* dConUa[s][j];

				dDbJUb += dDxBasis1D[s][j]
					* m_dataJacobian2D[iA+i][iB+s]
					* dConUb[i][s];

				dCovDaUb += dDxBasis1D[s][i] * dataUb[k][iA+s][iB+j];
					//( m_dataCovMetric2DB[iA+s][iB+j][0] * dataUa[k][iA+s][iB+j]
					//+ m_dataCovMetric2DB[iA+s][iB+j][1] * dataUb[k][iA+s][iB+j]);
				dCovDbUa += dDxBasis1D[s][j] * dataUa[k][iA+i][iB+s];
					//( m_dataCovMetric2DA[iA+i][iB+s][0] * dataUa[k][iA+i][iB+s]
					//+ m_dataCovMetric2DA[iA+i][iB+s][1] * dataUb[k][iA+i][iB+s]);
			}

			dDaJUa /= GetElementDeltaA();
			dDbJUb /= GetElementDeltaB();

			dCovDaUb /= GetElementDeltaA();
			dCovDbUa /= GetElementDeltaB();

			// Radial vorticity (vertical component)
			m_dataVorticity[k][iA+i][iB+j] =
				(dCovDaUb - dCovDbUa) / m_dataJacobian2D[iA+i][iB+j];

			// Compute the divergence at node
			m_dataDivergence[k][iA+i][iB+j] =
				(dDaJUa + dDbJUb) / m_dataJacobian2D[iA+i][iB+j];

/*
			// Compute covariant derivatives at node
			double dCovDaUa = dDaUa
				+ m_dataChristoffelA[iA+i][iB+j][0] * dUa
				+ m_dataChristoffelA[iA+i][iB+j][1] * 0.5 * dUb;

			double dCovDaUb = dDaUb
				+ m_dataChristoffelB[iA+i][iB+j][0] * dUa
				+ m_dataChristoffelB[iA+i][iB+j][1] * 0.5 * dUb;

			double dCovDbUa = dDbUa
				+ m_dataChristoffelA[iA+i][iB+j][1] * 0.5 * dUa
				+ m_dataChristoffelA[iA+i][iB+j][2] * dUb;

			double dCovDbUb = dDbUb
				+ m_dataChristoffelB[iA+i][iB+j][1] * 0.5 * dUa
				+ m_dataChristoffelB[iA+i][iB+j][2] * dUb;

			// Compute curl at node
			m_dataVorticity[k][iA+i][iB+j] = m_dataJacobian2D[iA+i][iB+j] * (
				+ m_dataContraMetric2DA[iA+i][iB+j][0] * dDaUb
				+ m_dataContraMetric2DA[iA+i][iB+j][1] * dDbUb
				- m_dataContraMetric2DB[iA+i][iB+j][0] * dDaUa
				- m_dataContraMetric2DB[iA+i][iB+j][1] * dDbUa);

			// Compute the divergence at node
			m_dataDivergence[k][iA+i][iB+j] = dDaUa + dDbUb;
*/
		}
		}
	}
	}
	}
}

///////////////////////////////////////////////////////////////////////////////

void GridPatchCartesianGLL::ComputeVorticityDivergence(
	int iDataIndex
) {
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

	dataUa.AttachToData(&(dataState[0][0][0][0]));
	dataUb.AttachToData(&(dataState[1][0][0][0]));

	// Compute the radial component of the curl of the velocity field
	ComputeCurlAndDiv(dataUa, dataUb);
}

///////////////////////////////////////////////////////////////////////////////

void GridPatchCartesianGLL::InterpolateData(
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
	if (dAlpha.GetRows() != dBeta.GetRows()) {
		_EXCEPTIONT("Point vectors must have equivalent length.");
	}

	// Vector for storage interpolated points
	DataArray1D<double> dAInterpCoeffs(m_nHorizontalOrder);
	DataArray1D<double> dBInterpCoeffs(m_nHorizontalOrder);
	DataArray1D<double> dADiffCoeffs(m_nHorizontalOrder);
	DataArray1D<double> dBDiffCoeffs(m_nHorizontalOrder);
	DataArray1D<double> dAInterpPt(m_nHorizontalOrder);

	// Element-wise grid spacing
	double dDeltaA = GetElementDeltaA();
	double dDeltaB = GetElementDeltaB();

	// Physical constants
	const PhysicalConstants & phys = m_grid.GetModel().GetPhysicalConstants();

	// Perform interpolation on all variables
	int nComponents;
	int nRElements = m_grid.GetRElements();

	// State Data: Perform interpolation on all variables
	if (eDataType == DataType_State) {
		nComponents = m_datavecStateNode[0].GetSize(0);

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

	} else {
		_EXCEPTIONT("Invalid DataType");
	}

	// Check dInterpData
	if (dInterpData.GetRows() != nComponents) {
		_EXCEPTIONT("Invalid size in InterpData (0)");
	}
	if (dInterpData.GetColumns() != dREta.GetRows()) {
		_EXCEPTIONT("Invalid size in InterpData (1)");
	}
	if (dInterpData.GetSubColumns() != dAlpha.GetRows()) {
		_EXCEPTIONT("Invalid size in InterpData (2)");
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

		} else {
			opInterp.InitializeIdentity(1);
		}

		// Buffer storage in column
		DataArray1D<double> dColumnData(nRElements);

		// Get a pointer to the 3D data structure
		DataArray3D<double> pData;
		DataArray3D<double> pDataRef;

		pData.SetSize(
			nRElements,
			m_box.GetATotalWidth(),
			m_box.GetBTotalWidth());

		pDataRef.SetSize(
			nRElements,
			m_box.GetATotalWidth(),
			m_box.GetBTotalWidth());

		if (eDataType == DataType_State) {
			if (eDataLocation == DataLocation_Node) {
				pData.AttachToData(&(m_datavecStateNode[0][c][0][0][0]));
				pDataRef.AttachToData(&(m_dataRefStateNode[c][0][0][0]));
			} else if (eDataLocation == DataLocation_REdge) {
				pData.AttachToData(&(m_datavecStateREdge[0][c][0][0][0]));
				pDataRef.AttachToData(&(m_dataRefStateREdge[c][0][0][0]));
			} else {
				_EXCEPTIONT("Invalid DataLocation");
			}

		} else if (eDataType == DataType_Tracers) {
			pData.AttachToData(&(m_datavecTracers[0][c][0][0][0]));

		} else if (eDataType == DataType_Topography) {
			pData.AttachToData(&(m_dataTopography[0][0]));

		} else if (eDataType == DataType_Vorticity) {
			pData.AttachToData(&(m_dataVorticity[0][0][0]));

		} else if (eDataType == DataType_Divergence) {
			pData.AttachToData(&(m_dataDivergence[0][0][0]));

		} else if (eDataType == DataType_Temperature) {
			pData.AttachToData(&(m_dataTemperature[0][0][0]));

		} else {
			_EXCEPTIONT("Invalid DataType");
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
				(dAlpha[i] - m_dAEdge[m_box.GetAInteriorBegin()]) / dDeltaA;

			int iB =
				(dBeta[i] - m_dBEdge[m_box.GetBInteriorBegin()]) / dDeltaB;

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

				for (int m = 0; m < m_nHorizontalOrder; m++) {
				for (int n = 0; n < m_nHorizontalOrder; n++) {
					dColumnData[k] +=
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
						dColumnData[k] -=
							  dAInterpCoeffs[m]
							* dBInterpCoeffs[n]
							* pDataRef[k][iA+m][iB+n];
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
				dInterpData[c][k][i] = dColumnDataOut[k];
			}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void GridPatchCartesianGLL::TransformHaloVelocities(
	int iDataUpdate
) {
	// Transform not necessary on Cartesian grid
}

///////////////////////////////////////////////////////////////////////////////

void GridPatchCartesianGLL::TransformTopographyDeriv() {

	// Transform not necessary on Cartesian grid
}

///////////////////////////////////////////////////////////////////////////////

