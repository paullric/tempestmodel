///////////////////////////////////////////////////////////////////////////////
///
///	\file    GridCartesianGLL.cpp
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

#include "GridCartesianGLL.h"
#include "GridPatchCartesianGLL.h"
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
#include <sstream>

///////////////////////////////////////////////////////////////////////////////

GridCartesianGLL::GridCartesianGLL(
	const Model & model,
	int nBaseResolution,
	int nRefinementRatio,
	int nHorizontalOrder,
	int nVerticalOrder,
	int nRElements,
    double dGDim[]
) :
    // Bring through the grid dimensions
    //m_dGDim({0.0, 0.0, 0.0, 0.0, 0.0, 0.0}),

	// Call up the stack
	GridGLL::GridGLL(
		model,
		nBaseResolution,
		nRefinementRatio,
		nHorizontalOrder,
		nVerticalOrder,
		nRElements)
{
	// Set the reference length scale
	m_dReferenceLength = 10000.0;
    
    // Bring through the grid dimensions
    m_dGDim[0] = dGDim[0]; m_dGDim[1] = dGDim[1];
    m_dGDim[2] = dGDim[2]; m_dGDim[3] = dGDim[3];
    m_dGDim[4] = dGDim[4]; m_dGDim[5] = dGDim[5];
}

///////////////////////////////////////////////////////////////////////////////

void GridCartesianGLL::Initialize() {

	// Call up the stack
	Grid::Initialize();

	// Initialize the vertical coordinate
	double dDeltaElement =
		static_cast<double>(m_nVerticalOrder)
		/ static_cast<double>(m_nRElements);

	InitializeVerticalCoordinate(
		GridSpacingMixedGaussLobatto(dDeltaElement, 0.0, m_nVerticalOrder)
	);
    
    // Interpolation coefficients from nodes to interfaces and vice versa
	m_dInterpNodeToREdge.Initialize(m_nVerticalOrder+1, m_nVerticalOrder);
	m_dInterpREdgeToNode.Initialize(m_nVerticalOrder, m_nVerticalOrder+1);
    
	for (int n = 0; n < m_nVerticalOrder+1; n++) {
		PolynomialInterp::LagrangianPolynomialCoeffs(
            m_nVerticalOrder,
            m_dREtaLevels,
            m_dInterpNodeToREdge[n],
            m_dREtaInterfaces[n]);
	}
	for (int n = 0; n < m_nVerticalOrder; n++) {
		PolynomialInterp::LagrangianPolynomialCoeffs(
            m_nVerticalOrder+1,
            m_dREtaInterfaces,
            m_dInterpREdgeToNode[n],
            m_dREtaLevels[n]);
	}
    
	// Distribute patches to processors
	Grid::DistributePatches();
    
	// Set up connectivity
	Grid::InitializeConnectivity();

}

///////////////////////////////////////////////////////////////////////////////

void GridCartesianGLL::AddDefaultPatches() {

	// Determine number of usable processors
	int nCommSize;
	MPI_Comm_size(MPI_COMM_WORLD, &nCommSize);

	int nProcsPerDirection = Max((int)ISqrt(nCommSize), 1);

	int nProcsPerPanel = nProcsPerDirection * nProcsPerDirection;

	int nDistributedPatches = nProcsPerPanel;

	if (nDistributedPatches < nCommSize) {
		Announce("WARNING: Patch / thread mismatch: "
			"%i threads will be unutilized",
			nCommSize - nDistributedPatches);
	}

	// Determine arrangement of elements on processors
	if (GetABaseResolution() % nProcsPerDirection != 0) {
		_EXCEPTIONT("\n(UNIMPLEMENTED) Currently elements must be "
			"equally divided among processors.");
	}

	int nElementsPerDirection = GetABaseResolution() / nProcsPerDirection;
	DataVector<int> iBoxBegin;
	iBoxBegin.Initialize(nProcsPerDirection + 1);

	iBoxBegin[0] = 0;
	for (int n = 1; n < nProcsPerDirection; n++) {
		iBoxBegin[n] = n * nElementsPerDirection;
	}
	iBoxBegin[nProcsPerDirection] = GetABaseResolution();

	// Create master patch for each panel
	std::vector<GridPatch *> pPatches;
	pPatches.reserve(nDistributedPatches);

	// Single panel 0 implementation (Cartesian Grid)
	// Rectangular alpha-wise patches that span all of the beta direction
	// (as many as there are processors available)
    //std::cout << "nProcsPerDirection" << nProcsPerDirection << "\n";
    int n = 0;
	for (int i = 0; i < nProcsPerDirection; i++) {
	int j = 0;
		double dDeltaA = 0.5 * M_PI / GetABaseResolution();

		GridSpacingGaussLobattoRepeated
			glspacing(dDeltaA, -0.25 * M_PI, m_nHorizontalOrder);

        // Patch strips that span beta
		PatchBox boxMaster(
			n, 0, m_model.GetHaloElements(),
			m_nHorizontalOrder * iBoxBegin[i],
			m_nHorizontalOrder * iBoxBegin[i+1],
			0.0,
			m_nHorizontalOrder * nElementsPerDirection,
			glspacing,
			glspacing);

		int ixPatch = n * nProcsPerPanel + i * nProcsPerDirection + j;
        //std::cout << "ixPatch" << ixPatch << "\n";

		pPatches.push_back(
			AddPatch(
				new GridPatchCartesianGLL(
					(*this),
					ixPatch,
					boxMaster,
					m_nHorizontalOrder,
					m_nVerticalOrder)));
	}

	if (pPatches.size() != nDistributedPatches) {
		_EXCEPTIONT("Logic error");
	}
}

///////////////////////////////////////////////////////////////////////////////

void GridCartesianGLL::InitializeVerticalCoordinate(
	const GridSpacing & aGridSpacing
) {
	// Call to Grid
	Grid::InitializeVerticalCoordinate(aGridSpacing);
}

///////////////////////////////////////////////////////////////////////////////

void GridCartesianGLL::GetReferenceGridBounds(
    double & dX0,
    double & dX1,
    double & dY0,
    double & dY1
) {
	dX0 = m_dGDim[0];
	dX1 = m_dGDim[1];
	dY0 = m_dGDim[2];
	dY1 = m_dGDim[3];
}

///////////////////////////////////////////////////////////////////////////////

void GridCartesianGLL::ConvertReferenceToPatchCoord(
	const DataVector<double> & dXReference,
	const DataVector<double> & dYReference,
	DataVector<double> & dAlpha,
	DataVector<double> & dBeta,
	DataVector<int> & iPatch
) const {
    
    // No conversion needed for cartesian grid but left the dimension check
	if ((dXReference.GetRows() != dYReference.GetRows()) ||
		(dXReference.GetRows() != dAlpha.GetRows()) ||
		(dXReference.GetRows() != dBeta.GetRows()) ||
		(dXReference.GetRows() != iPatch.GetRows())
        ) {
		_EXCEPTIONT("Dimension mismatch: All arrays must have same length");
	}
}

///////////////////////////////////////////////////////////////////////////////

void GridCartesianGLL::GetPatchFromCoordinateIndex(
    int iRefinementLevel,
    const DataVector<int> & vecIxA,
    const DataVector<int> & vecIxB,
    const DataVector<int> & vecPanel,
    DataVector<int> & vecPatchIndex,
    int nVectorLength
) {
	// Set vector length
	if (nVectorLength == (-1)) {
		nVectorLength = vecIxA.GetRows();
	}
    
	// Check arguments
	if ((vecIxA.GetRows() < nVectorLength) ||
		(vecIxB.GetRows() < nVectorLength) ||
		(vecPanel.GetRows() < nVectorLength)
        ) {
		_EXCEPTIONT("Argument vector length mismatch");
	}
	if (iRefinementLevel < 0) {
		_EXCEPTIONT("Refinement level must be positive");
	}
    
	// Calculate local resolution
	int nLocalResolution =
    m_nHorizontalOrder * GetABaseResolution(iRefinementLevel);
    
	// Loop through all entries
	int iLastPatch = GridPatch::InvalidIndex;
	for (int i = 0; i < nVectorLength; i++) {
        
		int iA = vecIxA[i];
		int iB = vecIxB[i];
		int iP = vecPanel[i];
        
		// Transform to alternative panel if necessary
		if ((vecIxA[i] < 0) || (vecIxA[i] >= nLocalResolution) ||
			(vecIxB[i] < 0) || (vecIxB[i] >= nLocalResolution)
            ) {
			bool fSwitchAB;
			bool fSwitchPar;
			bool fSwitchPerp;
            
			// Coordinate out of bounds
			if ((iA == (-1)) && (iB == (-1))) {
				vecPatchIndex[i] = GridPatch::InvalidIndex;
				continue;
			}
		}
        
		// Check the last patch searched
		if (iLastPatch != GridPatch::InvalidIndex) {
			const GridPatch * pPatch = GetPatch(iLastPatch);
            
			const PatchBox & box = pPatch->GetPatchBox();
            
			if (box.ContainsGlobalPoint(iP, iA, iB)) {
				vecPatchIndex[i] = pPatch->GetPatchIndex();
				continue;
			}
		}
        
		// Check all other patches
		int n;
		for (n = 0; n < GetPatchCount(); n++) {
			const GridPatch * pPatch = GetPatch(n);
            
			const PatchBox & box = pPatch->GetPatchBox();
            
			if (box.ContainsGlobalPoint(iP, iA, iB)) {
				vecPatchIndex[i] = pPatch->GetPatchIndex();
				iLastPatch = n;
				break;
			}
		}
		if (n == GetPatchCount()) {
			vecPatchIndex[i] = GridPatch::InvalidIndex;
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void GridCartesianGLL::GetOpposingDirection(
    int ixPanelSrc,
    int ixPanelDest,
    Direction dir,
    Direction & dirOpposing,
    bool & fSwitchParallel
) const {
	if ((ixPanelSrc < 0) || (ixPanelSrc > 5)) {
		_EXCEPTIONT("Invalid value for ixPanelSrc: Out of range");
	}
	if ((ixPanelDest < 0) || (ixPanelDest > 5)) {
		_EXCEPTIONT("Invalid value for ixPanelDest: Out of range");
	}
    
	// Get the opposing direction for Cartesian panels Source = Destination
    if (ixPanelDest != ixPanelSrc) {
        _EXCEPTIONT("ERROR: Soure and Destination panels must be equal!");
    }
    
	if (dir == Direction_Right) {
        dirOpposing = Direction_Left;
    } else if (dir == Direction_Top) {
        dirOpposing = Direction_Bottom;
    } else if (dir == Direction_Left) {
        dirOpposing = Direction_Right;
    } else if (dir == Direction_Bottom) {
        dirOpposing = Direction_Top;
    } else if (dir == Direction_TopRight) {
        dirOpposing = Direction_BottomLeft;
    } else if (dir == Direction_TopLeft) {
        dirOpposing = Direction_BottomRight;
    } else if (dir == Direction_BottomLeft) {
        dirOpposing = Direction_TopRight;
    } else if (dir == Direction_BottomRight) {
        dirOpposing = Direction_TopLeft;
    }
    
	/*/ Check panel indices to determine parallel switch of vectors
	fSwitchParallel = false;
	if ((ixPanelSrc == 1) && (ixPanelDest == 5)) {
		fSwitchParallel = true;
        
	} else if ((ixPanelSrc == 2) && (ixPanelDest > 3)) {
		fSwitchParallel = true;
        
	} else if ((ixPanelSrc == 3) && (ixPanelDest == 4)) {
		fSwitchParallel = true;
        
	} else if (
               (ixPanelSrc == 4) &&
               ((ixPanelDest == 2) || (ixPanelDest == 3))
               ) {
		fSwitchParallel = true;
        
	} else if (
               (ixPanelSrc == 5) &&
               ((ixPanelDest == 2) || (ixPanelDest == 1))
               ) {
		fSwitchParallel = true;
	}
    */
}

///////////////////////////////////////////////////////////////////////////////

void GridCartesianGLL::ApplyDSS(
	int iDataUpdate,
	DataType eDataType
) {
	// Exchange data between nodes
	Exchange(eDataType, iDataUpdate);

	// Post-process velocities across panel edges and
	// perform direct stiffness summation (DSS)
	for (int n = 0; n < GetActivePatchCount(); n++) {
        //std::cout << "DSS patch loop count: " << n << "\n";
		GridPatchCartesianGLL * pPatch =
			dynamic_cast<GridPatchCartesianGLL*>(GetActivePatch(n));

		const PatchBox & box = pPatch->GetPatchBox();
        
        // Patch-specific quantities
		int nElementCountA = pPatch->GetElementCountA();
		int nElementCountB = pPatch->GetElementCountB();
        
		// Apply panel transforms to velocity data
		if (eDataType == DataType_State) {
			pPatch->TransformHaloVelocities(iDataUpdate);
		}

		// Loop through all components associated with this DataType
		int nComponents;
		if (eDataType == DataType_State) {
			nComponents = m_model.GetEquationSet().GetComponents();
		} else if (eDataType == DataType_Tracers) {
			nComponents = m_model.GetEquationSet().GetTracers();
		} else if (eDataType == DataType_Vorticity) {
			nComponents = 1;
		} else if (eDataType == DataType_Divergence) {
			nComponents = 1;
		} else {
			_EXCEPTIONT("Invalid DataType");
		}

		// Perform Direct Stiffness Summation (DSS)
		for (int c = 0; c < nComponents; c++) {

			// Obtain the array of working data
			int nRElements = GetRElements();
			double *** pDataUpdate;
			if (eDataType == DataType_State) {
				pDataUpdate =
					pPatch->GetDataState(iDataUpdate, GetVarLocation(c))[c];

				if (GetVarLocation(c) == DataLocation_REdge) {
					nRElements++;
				}

			} else if (eDataType == DataType_Tracers) {
				pDataUpdate =
					pPatch->GetDataTracers(iDataUpdate)[c];
			} else if (eDataType == DataType_Vorticity) {
				pDataUpdate = pPatch->GetDataVorticity();
			} else if (eDataType == DataType_Divergence) {
				pDataUpdate = pPatch->GetDataDivergence();
			}
            
            /* BAD IMPLEMENTATION OF BOUNDARY CONDITIONS - JEG
            if (n == 0) {
            //std::cout << nComponents << "\n";
            }
            std::cout << pDataUpdate[0][0][0] << "\n";
            
            // Apply horizontal direction BC's (PERIODIC ONLY FOR NOW)
            for (int k = 0; k < nRElements; k++) {
                // Loop in beta over first and last alpha edge. Apply X BC
                int aFirst = 0;
                int aLast= nAElements * m_nHorizontalOrder;
                for (int j = 0; j <= nBElements * m_nHorizontalOrder; j++) {
                    // Periodic condition in X
                    pDataUpdate[k][aFirst][j] = pDataUpdate[k][aLast][j];
                }
                // Loop in alpha over first and last beta edge. Apply Y BC
                int bFirst = 0;
                int bLast = nBElements * m_nHorizontalOrder;
                for (int i = 0; i <= nAElements * m_nHorizontalOrder; i++) {
                    // Periodic condition in Y
                    pDataUpdate[k][i][bFirst] = pDataUpdate[k][i][bLast];
                }
            }
        
            // Apply the vertical direction BC (NO FLUX TOP AND BOTTOM FOR NOW)
            int kFirst = 0;
            int kLast = nRElements - 1;
            //std::cout << nRElements << "\n";
            
            for (int i = 0; i < nAElements * m_nHorizontalOrder; i++) {
                for (int j = 0; j < nBElements * m_nHorizontalOrder; j++) {
                    pDataUpdate[kFirst][i][j] = pDataUpdate[kFirst+1][i][j];
                    pDataUpdate[kLast][i][j] = pDataUpdate[kLast-1][i][j];
                }
            }
            */
            
            //std::cout << "Entering DSS averaging loop!";
            // Averaging DSS across patch boundaries
            for (int k = 0; k < nRElements; k++) {
                
				// Average in the alpha direction
				for (int a = 0; a <= nElementCountA; a++) {
					int iA = a * m_nHorizontalOrder + box.GetHaloElements();
                    
					// Averaging done at the corners of the panel
					int jBegin = box.GetBInteriorBegin()-1;
					int jEnd = box.GetBInteriorEnd()+1;
                    
					// Perform averaging across edge of patch
					for (int j = jBegin; j < jEnd; j++) {
						pDataUpdate[k][iA][j] = 0.5 * (
                        + pDataUpdate[k][iA  ][j]
                        + pDataUpdate[k][iA-1][j]);
                        
						pDataUpdate[k][iA-1][j] = pDataUpdate[k][iA][j];
					}
				}
                
				// Average in the beta direction
				for (int b = 0; b <= nElementCountB; b++) {
					int iB = b * m_nHorizontalOrder + box.GetHaloElements();
                    
					// Averaging done at the corners of the panel
					int iBegin = box.GetAInteriorBegin()-1;
					int iEnd = box.GetAInteriorEnd()+1;
                    
					for (int i = iBegin; i < iEnd; i++) {
						pDataUpdate[k][i][iB] = 0.5 * (
                        + pDataUpdate[k][i][iB  ]
                        + pDataUpdate[k][i][iB-1]);
                        
						pDataUpdate[k][i][iB-1] = pDataUpdate[k][i][iB];
					}
				}
            }
		}
	}
}

///////////////////////////////////////////////////////////////////////////////
