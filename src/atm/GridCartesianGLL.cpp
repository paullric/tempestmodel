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

/////////////////////////////////////////////////////////////////////////////

GridCartesianGLL::GridCartesianGLL(
	Model & model,
	int nBaseResolutionA,
	int nBaseResolutionB,
	int nRefinementRatio,
	int nHorizontalOrder,
	int nVerticalOrder,
	int nRElements,
	double dGDim[],
	double dRefLat,
	VerticalStaggering eVerticalStaggering
) :
	// Call up the stack
	GridGLL::GridGLL(
		model,
		nBaseResolutionA,
		nBaseResolutionB,
		nRefinementRatio,
		nHorizontalOrder,
		nVerticalOrder,
		nRElements,
		eVerticalStaggering)
{
	// Set the reference length scale (110km)
	m_dReferenceLength = 110000.0;

	// Bring in the reference latitude (if any) for large regions where the
	// beta plane approximation is necessary in the equations
	m_dRefLat = dRefLat;

	// Bring through the grid dimensions
	m_dGDim[0] = dGDim[0]; m_dGDim[1] = dGDim[1];
	m_dGDim[2] = dGDim[2]; m_dGDim[3] = dGDim[3];
	m_dGDim[4] = dGDim[4]; m_dGDim[5] = dGDim[5];

	m_dZtop = dGDim[5];
}

///////////////////////////////////////////////////////////////////////////////

void GridCartesianGLL::Initialize() {

	// Call up the stack
	GridGLL::Initialize();

	// Distribute patches to processors
	Grid::DistributePatches();

	// Set up connectivity
	Grid::InitializeConnectivity();

}

///////////////////////////////////////////////////////////////////////////////

void GridCartesianGLL::AddDefaultPatches() {

	// Verify no Patches have been previously added
	if (GetPatchCount() != 0) {
		_EXCEPTIONT("AddDefaultPatches() must be called on an empty Grid");
	}

	// Determine number of usable processors
	int nCommSize;
	MPI_Comm_size(MPI_COMM_WORLD, &nCommSize);

	int nProcsPerDirection = nCommSize;

	int nDistributedPatches = nProcsPerDirection;

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

	// Patch grid spacing
	double dDeltaA = (m_dGDim[1] - m_dGDim[0])
		/ static_cast<double>(GetABaseResolution());
	double dDeltaB = (m_dGDim[3] - m_dGDim[2])
		/ static_cast<double>(GetBBaseResolution());

	GridSpacingGaussLobattoRepeated
		glspacingAlpha(dDeltaA, m_dGDim[0], m_nHorizontalOrder);
	GridSpacingGaussLobattoRepeated
		glspacingBeta(dDeltaB, m_dGDim[2], m_nHorizontalOrder);

	// Single panel 0 implementation (Cartesian Grid)
	// Rectangular alpha-wise patches that span all of the beta direction
	// (as many as there are processors available)
	int n = 0;
	for (int i = 0; i < nProcsPerDirection; i++) {

		// Patch strips that span beta
		PatchBox boxMaster(
			n, 0, m_model.GetHaloElements(),
			m_nHorizontalOrder * iBoxBegin[i],
			m_nHorizontalOrder * iBoxBegin[i+1],
			0,
			m_nHorizontalOrder * GetBBaseResolution(),
			glspacingAlpha,
			glspacingBeta);

		Grid::AddPatch(
			new GridPatchCartesianGLL(
				(*this),
				i,
				boxMaster,
				m_nHorizontalOrder,
				m_nVerticalOrder,
				m_dGDim,
				m_dRefLat));
	}
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

	// Loop over all coordinates
	for (int i = 0; i < dXReference.GetRows(); i++) {

		// Reference and computational coordinates are the same
		dAlpha[i] = dXReference[i];
		dBeta[i]  = dYReference[i];

		// Loop over all patches
		int n = 0;

		for (; n < GetPatchCount(); n++) {
			const GridPatch * pPatch = GetPatch(n);
			const PatchBox & box = pPatch->GetPatchBox();

			if ((dAlpha[i] >= box.GetAEdge(box.GetAInteriorBegin())) &&
				(dAlpha[i] <= box.GetAEdge(box.GetAInteriorEnd())) &&
				(dBeta[i] >= box.GetBEdge(box.GetBInteriorBegin())) &&
				(dBeta[i] <= box.GetBEdge(box.GetBInteriorEnd()))
			) {
				iPatch[i] = pPatch->GetPatchIndex();
				break;
			}
		}

		if (n == GetPatchCount()) {
			_EXCEPTION4("Unable to find associated patch for node:\n"
				"(%1.5e, %1.5e) : (%1.5e, %1.5e)",
				dXReference[i], dYReference[i],
				dAlpha[i], dBeta[i]);
		}
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
	int nLocalResolutionA =
		m_nHorizontalOrder * GetABaseResolution(iRefinementLevel);
	int nLocalResolutionB =
		m_nHorizontalOrder * GetBBaseResolution(iRefinementLevel);

	// Loop through all entries
	int iLastPatch = GridPatch::InvalidIndex;
	for (int i = 0; i < nVectorLength; i++) {

		int iA = vecIxA[i];
		int iB = vecIxB[i];
		int iP = vecPanel[i];

		// Wrap global indices
		if (iA < 0) {
			iA += nLocalResolutionA;
		}
		if (iA >= nLocalResolutionA) {
			iA -= nLocalResolutionA;
		}
		if (iB < 0) {
			iB += nLocalResolutionB;
		}
		if (iB >= nLocalResolutionB) {
			iB -= nLocalResolutionB;
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

			_EXCEPTIONT("(LOGIC ERROR) Invalid global coordinate");
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void GridCartesianGLL::GetOpposingDirection(
	int ixPanelSrc,
	int ixPanelDest,
	Direction dir,
	Direction & dirOpposing,
	bool & fSwitchParallel,
	bool & fSwitchPerpendicular
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

	// Do not switch directions across this edge
	fSwitchParallel = false;
	fSwitchPerpendicular = false;
}

///////////////////////////////////////////////////////////////////////////////

void GridCartesianGLL::ApplyBoundaryConditions(
	int iDataUpdate,
	DataType eDataType
) {
	_EXCEPTION();
	for (int n = 0; n < GetActivePatchCount(); n++) {
		GridPatchCartesianGLL * pPatch =
			dynamic_cast<GridPatchCartesianGLL*>(GetActivePatch(n));

		pPatch->ApplyBoundaryConditions(iDataUpdate, eDataType);
	}
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
