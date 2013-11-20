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

#include "Direction.h"
#include "CubedSphereTrans.h"
#include "PolynomialInterp.h"
#include "GaussQuadrature.h"
#include "GaussLobattoQuadrature.h"

#include "Announce.h"
#include "MathHelper.h"

#include <cmath>

///////////////////////////////////////////////////////////////////////////////

GridPatchCartesianGLL::GridPatchCartesianGLL(
	GridCartesianGLL & grid,
	int ixPatch,
	const PatchBox & box,
	int nHorizontalOrder,
	int nVerticalOrder,
    double dGDim[]
) :
	GridPatchGLL(
		grid,
		ixPatch,
		box,
		nHorizontalOrder,
		nVerticalOrder)
{

	// Bring in the grid dimensions as a member variable
    // Bring through the grid dimensions
    m_dGDim[0] = dGDim[0]; m_dGDim[1] = dGDim[1];
    m_dGDim[2] = dGDim[2]; m_dGDim[3] = dGDim[3];
    m_dGDim[4] = dGDim[4]; m_dGDim[5] = dGDim[5];
    
    m_ixNeighborPanel.resize(8);
    m_ixNeighborPanel[0] = 0.; m_ixNeighborPanel[1] = 0.;
    m_ixNeighborPanel[2] = 0.; m_ixNeighborPanel[3] = 0.;
    m_ixNeighborPanel[4] = 0.; m_ixNeighborPanel[5] = 0.;
    m_ixNeighborPanel[6] = 0.; m_ixNeighborPanel[7] = 0.;
}

///////////////////////////////////////////////////////////////////////////////

void GridPatchCartesianGLL::InitializeDataLocal() {

	// Allocate data
	GridPatch::InitializeDataLocal();

	// Physical constants
	const PhysicalConstants & phys = m_grid.GetModel().GetPhysicalConstants();
    
    // Get the maximum value of the ABP coordinates (to scale dataLon and Lat)
    double maxA = m_box.GetANode(m_box.GetATotalWidth()-1);
    double maxB = m_box.GetBNode(m_box.GetBTotalWidth()-1);
    double XL = (m_dGDim[1] - m_dGDim[0]);
    double YL = (m_dGDim[3] - m_dGDim[2]);
    
    // Initialize the longitude and latitude at each node
	for (int i = 0; i < m_box.GetATotalWidth(); i++) {
	for (int j = 0; j < m_box.GetBTotalWidth(); j++) {
        // TODO: Initialize the force with user input
        // (latitude input for reference force)
        m_dataLat[i][j] = 0.5 * XL *
                        m_box.GetANode(i) / maxA + 0.5 * XL;
        m_dataLon[i][j] = 0.5 * YL *
                        m_box.GetBNode(j) / maxB;
		m_dataCoriolisF[i][j] = m_dataCoriolisF[i][j] +
                                phys.GetOmega() * m_box.GetBNode(j);
        
        //std::cout << m_dataLat[i][j] << " "
                  //<< m_dataLon[i][j] << "\n";
                  //<< m_dataZLevels[k][i][j] << "\n";

	}
	}
}

///////////////////////////////////////////////////////////////////////////////

void GridPatchCartesianGLL::EvaluateTopography(
                                        const TestCase & test
                                        ) {
	const PhysicalConstants & phys = m_grid.GetModel().GetPhysicalConstants();
    
	for (int i = 0; i < m_box.GetATotalWidth(); i++) {
        for (int j = 0; j < m_box.GetBTotalWidth(); j++) {
            double dLon;
            double dLat;
            
            dLon = m_box.GetANode(i);
            dLat = m_box.GetBNode(j);
            
            m_dataTopography[i][j] = test.EvaluateTopography(phys, dLon, dLat);
            
            if (m_dataTopography[i][j] >= m_grid.GetZtop()) {
                _EXCEPTIONT("TestCase topography exceeds model top.");
            }
        }
	}
}

///////////////////////////////////////////////////////////////////////////////

void GridPatchCartesianGLL::EvaluateGeometricTerms() {
    
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
    
	// Derivatives of basis functions
	GridCartesianGLL & gridCartesianGLL =
    dynamic_cast<GridCartesianGLL &>(m_grid);
    
	const DataMatrix<double> & dDxBasis1D = gridCartesianGLL.GetDxBasis1D();
    
	// Element spacing at this refinement level
	double dElementDeltaA =
    m_box.GetAEdge(m_box.GetHaloElements() + m_nHorizontalOrder)
    - m_box.GetAEdge(m_box.GetHaloElements());
    
	double dElementDeltaB =
    m_box.GetBEdge(m_box.GetHaloElements() + m_nHorizontalOrder)
    - m_box.GetBEdge(m_box.GetHaloElements());
    
#pragma message "Implement rectangular grid elements"
	if (fabs(dElementDeltaA - dElementDeltaB) > 1.0e-12) {
		_EXCEPTIONT("Not implemented.");
	}
    
	// Vertical grid spacing
	double dElementDeltaXi =
    m_grid.GetREtaInterface(m_nVerticalOrder)
    - m_grid.GetREtaInterface(0);
    
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
                    
                    // Gnomonic coordinates
                    double dX = tan(m_box.GetANode(iA));
                    double dY = tan(m_box.GetBNode(iB));
                    double dDelta2 = (1.0 + dX * dX + dY * dY);
                    double dDelta = sqrt(dDelta2);
                    
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
                    
                    double dDaaZs = 0.0;
                    double dDabZs = 0.0;
                    double dDbbZs = 0.0;

                    #pragma message "FIX second derivatives"
                    // Initialize 2D Jacobian
                    m_dataJacobian2D[iA][iB] = 1.0;
                    
                    // Initialize 2D contravariant metric
                    double dContraMetricScale = 1.0;
                    
                    m_dataContraMetric2DA[iA][iB][0] = dContraMetricScale;
                    m_dataContraMetric2DA[iA][iB][1] = 0.0;
                    
                    m_dataContraMetric2DB[iA][iB][0] = 0.0;
                    m_dataContraMetric2DB[iA][iB][1] = dContraMetricScale;
                    
                    // Christoffel symbol components at each node
                    // (off-diagonal element are doubled due to symmetry)
                    m_dataChristoffelA[i][j][0] = 0.0;
                    m_dataChristoffelA[i][j][1] = 0.0;
                    m_dataChristoffelA[i][j][2] = 0.0;
                    m_dataChristoffelB[i][j][0] = 0.0;
                    m_dataChristoffelB[i][j][1] = 0.0;
                    m_dataChristoffelB[i][j][2] = 0.0;
                    
                    // Vertical coordinate transform and its derivatives
                    for (int k = 0; k < m_grid.GetRElements(); k++) {
                        
                        // Sub-element index
                        int kx = k % m_nVerticalOrder;
                        
                        // Z coordinate
                        double dREta = m_grid.GetREtaLevel(k);
                        double dZ = dZs + (m_grid.GetZtop() - dZs) * dREta;
                        //std::cout << dZ << "\n";
                        
                        // Gal-Chen and Somerville (1975) linear terrain
                        // following coord
                        double dDaZ = (1.0 - m_grid.GetREtaLevel(k)) * dDaZs;
                        double dDbZ = (1.0 - m_grid.GetREtaLevel(k)) * dDbZs;
                        double dDaaZ = (1.0 - m_grid.GetREtaLevel(k)) * dDaaZs;
                        double dDabZ = (1.0 - m_grid.GetREtaLevel(k)) * dDabZs;
                        double dDbbZ = (1.0 - m_grid.GetREtaLevel(k)) * dDbbZs;
                        
                        double dDxZ = m_grid.GetZtop() - dZs;
                        double dDaxZ = - dDaZs;
                        double dDbxZ = - dDbZs;
                        double dDxxZ = 0.0;
                        
                        // Calculate pointwise Jacobian
                        m_dataJacobian[k][i][j] = dDxZ;
                        
                        // Element area associated with each model level GLL node
                        m_dataElementArea[k][i][j] =
                        m_dataJacobian[k][i][j]
                        * dWL[i] * dElementDeltaA
                        * dWL[j] * dElementDeltaB
                        * dWNode[kx] * dElementDeltaXi * dDxZ;
                        
                        // Contravariant metric components
                        double dContraMetricScale = 1.0;
                        
                        m_dataContraMetricA[k][iA][iB][0] =
                            m_dataContraMetric2DA[iA][iB][0];
                        m_dataContraMetricA[k][iA][iB][1] =
                            m_dataContraMetric2DA[iA][iB][1];
                        m_dataContraMetricA[k][iA][iB][2] =
                            - dContraMetricScale * dDxZ / dDaZ;
                        
                        m_dataContraMetricB[k][iA][iB][0] =
                            m_dataContraMetric2DB[iA][iB][0];
                        m_dataContraMetricB[k][iA][iB][1] =
                            m_dataContraMetric2DB[iA][iB][1];
                        m_dataContraMetricB[k][iA][iB][2] =
                            - dContraMetricScale * dDxZ / dDbZ;
                        
                        m_dataContraMetricXi[k][i][j][0] =
                            - dContraMetricScale * dDxZ / dDaZ;
                        m_dataContraMetricXi[k][i][j][1] =
                            - dContraMetricScale * dDxZ / dDbZ;
                        m_dataContraMetricXi[k][i][j][2] = (1.0 + dDaZ * dDaZ + dDbZ * dDbZ) / (dDxZ * dDxZ);
                        
                        // Vertical Christoffel symbol components
                        m_dataChristoffelXi[k][i][j][0] = dDaaZ;
                        m_dataChristoffelXi[k][i][j][1] = dDabZ;
                        m_dataChristoffelXi[k][i][j][2] = dDaxZ;
                        m_dataChristoffelXi[k][i][j][3] = dDbbZ;
                        m_dataChristoffelXi[k][i][j][4] = dDbxZ;
                        m_dataChristoffelXi[k][i][j][5] = dDxxZ;
                        
                        m_dataChristoffelXi[k][i][j][0] /= dDxZ;
                        m_dataChristoffelXi[k][i][j][1] /= dDxZ;
                        m_dataChristoffelXi[k][i][j][2] /= dDxZ;
                        m_dataChristoffelXi[k][i][j][3] /= dDxZ;
                        m_dataChristoffelXi[k][i][j][4] /= dDxZ;
                        m_dataChristoffelXi[k][i][j][5] /= dDxZ;
                        
                        // Orthonormalization coefficients
                        m_dataOrthonormNode[k][iA][iB][0] =
                        - dDaZs * (m_grid.GetZtop() - dZ)
                        / (m_grid.GetZtop() - dZs)
                        / (m_grid.GetZtop() - dZs);
                        
                        m_dataOrthonormNode[k][iA][iB][1] =
                        - dDbZs * (m_grid.GetZtop() - dZ)
                        / (m_grid.GetZtop() - dZs)
                        / (m_grid.GetZtop() - dZs);
                        
                        m_dataOrthonormNode[k][iA][iB][2] =
                        m_grid.GetZtop() / (m_grid.GetZtop() - dZs);
                    }
                    
                    // Metric terms at vertical interfaces
                    for (int k = 0; k <= m_grid.GetRElements(); k++) {
                        int kx = k % m_nVerticalOrder;
                        
                        // Z coordinate
                        double dREta = m_grid.GetREtaInterface(k);
                        double dZ = dZs + (m_grid.GetZtop() - dZs) * dREta;
                        
                        // Gal-Chen and Somerville (1975) linear terrain-following coord
                        double dDxZ = m_grid.GetZtop() - dZs;
                        
                        // Calculate pointwise Jacobian
                        double dJacobian = dDxZ;
                        
                        // Element area associated with each model interface GLL node
                        m_dataElementAreaREdge[k][iA][iB] =
                        dJacobian
                        * dWL[i] * dElementDeltaA
                        * dWL[j] * dElementDeltaB
                        * dWREdge[kx] * dElementDeltaXi;
                        
                        if ((k != 0) && (k != m_grid.GetRElements()) && (kx == 0)) {
                            m_dataElementAreaREdge[k][iA][iB] *= 2.0;
                        }
                        
                        // Orthonormalization coefficients
                        m_dataOrthonormREdge[k][iA][iB][0] =
                        - dDaZs * (m_grid.GetZtop() - dZ)
                        / (m_grid.GetZtop() - dZs)
                        / (m_grid.GetZtop() - dZs);
                        
                        m_dataOrthonormREdge[k][iA][iB][1] =
                        - dDbZs * (m_grid.GetZtop() - dZ)
                        / (m_grid.GetZtop() - dZs)
                        / (m_grid.GetZtop() - dZs);
                        
                        m_dataOrthonormREdge[k][iA][iB][2] =
                        m_grid.GetZtop() / (m_grid.GetZtop() - dZs);
                    }
                }
            }
            
        }
	}
}

///////////////////////////////////////////////////////////////////////////////

void GridPatchCartesianGLL::EvaluateTestCase(
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
    
    // Evaluate topography
	EvaluateTopography(test);

	// Evaluate geometric terms
	EvaluateGeometricTerms();

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
		}
	}
	}
	}
}

///////////////////////////////////////////////////////////////////////////////

void GridPatchCartesianGLL::ComputeCurlAndDiv(
	const GridData3D & dataUa,
	const GridData3D & dataUb
) const {
	// Parent grid
	const GridCartesianGLL & gridCSGLL = dynamic_cast<const GridCartesianGLL &>(m_grid);

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
				+ m_dataContraMetric2DA[iA+i][iB+j][0] * dCovDaUb
				+ m_dataContraMetric2DA[iA+i][iB+j][1] * dCovDbUb
				- m_dataContraMetric2DB[iA+i][iB+j][0] * dCovDaUa
				- m_dataContraMetric2DB[iA+i][iB+j][1] * dCovDbUa);

			// Compute the divergence at node
			m_dataDivergence[k][iA+i][iB+j] = dCovDaUa + dCovDbUb;
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
}

///////////////////////////////////////////////////////////////////////////////

void GridPatchCartesianGLL::InterpolateNodeToREdge(
	int iVar,
	int iDataIndex
) {

	// Working data
	GridData4D & dataNode  = GetDataState(iDataIndex, DataLocation_Node);
	GridData4D & dataREdge = GetDataState(iDataIndex, DataLocation_REdge);

	// Parent grid, containing the vertical remapping information
	GridCartesianGLL * pCSGLLGrid = dynamic_cast<GridCartesianGLL*>(&m_grid);
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

void GridPatchCartesianGLL::InterpolateREdgeToNode(
	int iVar,
	int iDataIndex
) {

	// Working data
	GridData4D & dataREdge = GetDataState(iDataIndex, DataLocation_REdge);
	GridData4D & dataNode  = GetDataState(iDataIndex, DataLocation_Node);

	// Parent grid, containing the vertical remapping information
	GridCartesianGLL * pCSGLLGrid = dynamic_cast<GridCartesianGLL*>(&m_grid);
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

void GridPatchCartesianGLL::InterpolateData(
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
	}
}

///////////////////////////////////////////////////////////////////////////////
// TODO: Transform not necessary in cartesian grid (delete method)
void GridPatchCartesianGLL::TransformHaloVelocities(
	int iDataUpdate
) {
	// Indices of velocities
	const int UIx = 0;
	const int VIx = 1;

	// Velocity data
	GridData4D * pDataVelocity =
		&(GetDataState(iDataUpdate, m_grid.GetVarLocation(UIx)));

	if (pDataVelocity->GetComponents() < 2) {
		_EXCEPTIONT("Invalid number of components.");
	}
}

///////////////////////////////////////////////////////////////////////////////

