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

#include "Announce.h"
#include "MathHelper.h"

#include <cmath>
#include <sstream>

/////////////////////////////////////////////////////////////////////////////

GridCartesianGLL::GridCartesianGLL(
	Model & model
) :
	GridGLL::GridGLL(model)
{ }

///////////////////////////////////////////////////////////////////////////////

void GridCartesianGLL::DefineParameters() {

	if (m_dcGridParameters.IsAttached()) {
		_EXCEPTIONT("Attempting to recall SetParameters");
	}

	// Set the dimension of the Cartesian domain
	m_dGDim.SetSize(6);

	// Add parameters for GridCartesianGLL
	m_dcGridParameters.PushDataChunk(&m_dGDim);
	m_dcGridParameters.PushDataChunk(&m_dRefLat);
	//m_dcGridParameters.PushDataChunk(&m_fCartesianXZ);

	// Call up the stack
	GridGLL::DefineParameters();
}

///////////////////////////////////////////////////////////////////////////////

void GridCartesianGLL::SetParameters(
	int nRElements,
	int nMaxPatchCount,
	int nABaseResolution,
	int nBBaseResolution,
	int nRefinementRatio,
	int nHorizontalOrder,
	int nVerticalOrder,
	double dGDim[],
	double dRefLat,
	int iLatBC[],
	bool fCartesianXZ,
	VerticalDiscretization eVerticalDiscretization,
	VerticalStaggering eVerticalStaggering
) {

	// Call up the stack
	GridGLL::SetParameters(
		nRElements,
		nMaxPatchCount,
		nABaseResolution,
		nBBaseResolution,
		nRefinementRatio,
		nHorizontalOrder,
		nVerticalOrder,
		eVerticalDiscretization,
		eVerticalStaggering
	);

	// Set the flag indicating a 2D XZ domain (turn off Ub update and Coriolis)
	m_fCartesianXZ = fCartesianXZ;

	// Bring through the grid dimensions
	m_dGDim[0] = dGDim[0]; m_dGDim[1] = dGDim[1];
	m_dGDim[2] = dGDim[2]; m_dGDim[3] = dGDim[3];
	m_dGDim[4] = dGDim[4]; m_dGDim[5] = dGDim[5];

	m_dZtop = dGDim[5];

	// Bring in the reference latitude (if any) for large regions where the
	// beta plane approximation is necessary in the equations
	m_dRefLat = dRefLat;
	
	// Set the boundary conditions for this run
	m_eBoundaryCondition[Direction_Right] = iLatBC[0];
	m_eBoundaryCondition[Direction_Top] = iLatBC[1];
	m_eBoundaryCondition[Direction_Left] = iLatBC[2];
	m_eBoundaryCondition[Direction_Bottom] = iLatBC[3];

	// Initialize maximum topographic height (8000 meters by default)
	m_dTopoHeight = 8000.0;
	// Initialize the scale height for the decay of topography features (SLEVE)
	m_dSL = m_dZtop / 5.0;
}

///////////////////////////////////////////////////////////////////////////////

void GridCartesianGLL::ApplyDefaultPatchLayout(
	int nPatchCount
) {

	// Verify patch count is positive
	if (nPatchCount < 1) {
		_EXCEPTIONT("nPatchCount must be a positive integer");
	}

	// Verify no Patches have been previously added
	if (m_nInitializedPatchBoxes != 0) {
		_EXCEPTIONT("ApplyDefaultPatchLayout() must be called on an empty Grid");
	}

	// Determine number of usable processors
	int nProcsPerDirection = nPatchCount;

	int nDistributedPatches = nProcsPerDirection;

	// Determine arrangement of elements on processors
	// SINGLE RECTANGULAR PATCH
	if (GetABaseResolution() * GetBBaseResolution() % nProcsPerDirection != 0) {
		_EXCEPTIONT("\n(UNIMPLEMENTED) Currently Alpha-Beta elements must be "
			"equally divided among processors.");
	}

	bool fCartXZ = GetIsCartesianXZ();
	int nProcsADirection = 0;
	int nProcsBDirection = 0;

	// Handle special case of 2 or 3 processors (laptop testing)
	if (fCartXZ || (nDistributedPatches < 4)) {
		nProcsADirection = nProcsPerDirection;
		nProcsBDirection = 1;
	}
/*
	// Allow 3D cases to split the Beta elements once (EVEN NUMBER OF PROCESSORS)
	else if (!fCartXZ && (nProcsPerDirection % 2 == 0)) {
		# pragma message "3D Cartesian processor layout..."
		nProcsADirection = (int) (nProcsPerDirection / 2);
		nProcsBDirection = 2;
	}
	else if (!fCartXZ && 
			 (nProcsADirection * nProcsBDirection != nDistributedPatches)) {
		_EXCEPTIONT("Even number of processors required"
			"for 3D Cartesian cases.");
	}
*/
	// Default to alpha wise strips
	else {
		nProcsADirection = nProcsPerDirection;
		nProcsBDirection = 1;
	}

	//nProcsADirection = nProcsPerDirection;
	//nProcsBDirection = 1;

	int nElementsPerDirectionA = GetABaseResolution() / nProcsADirection;
	DataArray1D<int> iBoxBeginA(nProcsADirection + 1);

	int nElementsPerDirectionB = GetBBaseResolution() / nProcsBDirection;
	DataArray1D<int> iBoxBeginB(nProcsBDirection + 1);

	// Set the patch array for Alpha
	iBoxBeginA[0] = 0;
	for (int n = 1; n < nProcsADirection; n++) {
		iBoxBeginA[n] = n * nElementsPerDirectionA;
	}
	iBoxBeginA[nProcsADirection] = GetABaseResolution();

	// Set the patch array for Beta
	iBoxBeginB[0] = 0;
	for (int n = 1; n < nProcsBDirection; n++) {
		iBoxBeginB[n] = n * nElementsPerDirectionB;
	}
	iBoxBeginB[nProcsBDirection] = GetBBaseResolution();

/*
	// Single panel 0 implementation (Cartesian Grid)
	// Rectangular alpha-wise patches that span all of the beta direction
	// (as many as there are processors available)
	for (int i = 0; i < nProcsPerDirection; i++) {

		// Patch strips that span beta
		m_aPatchBoxes[i] = PatchBox(
			0, 0, m_model.GetHaloElements(),
			m_nHorizontalOrder * iBoxBegin[i],
			m_nHorizontalOrder * iBoxBegin[i+1],
			0,
			m_nHorizontalOrder * GetBBaseResolution());
	}
*/
	// Create master patch for each panel
	int ixPatch = 0;

	for (int i = 0; i < nProcsADirection; i++) {
	for (int j = 0; j < nProcsBDirection; j++) {

		m_aPatchBoxes[ixPatch] = PatchBox(
			0, 0, m_model.GetHaloElements(),
			m_nHorizontalOrder * iBoxBeginA[i],
			m_nHorizontalOrder * iBoxBeginA[i+1],
			m_nHorizontalOrder * iBoxBeginB[j],
			m_nHorizontalOrder * iBoxBeginB[j+1]);

		ixPatch++;
	}
	}

	m_nInitializedPatchBoxes = ixPatch;
}

///////////////////////////////////////////////////////////////////////////////

GridPatch * GridCartesianGLL::NewPatch(
	int ixPatch
) {
	const PatchBox & box = GetPatchBox(ixPatch);

	return (
		new GridPatchCartesianGLL(
			(*this),
			ixPatch,
			box,
			m_nHorizontalOrder,
			m_nVerticalOrder));
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
	const DataArray1D<double> & dXReference,
	const DataArray1D<double> & dYReference,
	DataArray1D<double> & dAlpha,
	DataArray1D<double> & dBeta,
	DataArray1D<int> & iPatch
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
			const PatchBox & box = m_aPatchBoxes[n];

			double dElementDeltaA = (m_dGDim[1] - m_dGDim[0])
				/ static_cast<double>(GetABaseResolution());
			double dElementDeltaB = (m_dGDim[3] - m_dGDim[2])
				/ static_cast<double>(GetBBaseResolution());

			int iAElementInteriorBegin =
				box.GetAGlobalInteriorBegin() / m_nHorizontalOrder;
			int iAElementInteriorEnd =
				box.GetAGlobalInteriorEnd() / m_nHorizontalOrder;

			int iBElementInteriorBegin =
				box.GetBGlobalInteriorBegin() / m_nHorizontalOrder;
			int iBElementInteriorEnd =
				box.GetBGlobalInteriorEnd() / m_nHorizontalOrder;

			if ((box.GetAGlobalInteriorBegin() % m_nHorizontalOrder != 0) ||
				(box.GetAGlobalInteriorEnd()   % m_nHorizontalOrder != 0) ||
				(box.GetBGlobalInteriorBegin() % m_nHorizontalOrder != 0) ||
				(box.GetBGlobalInteriorEnd()   % m_nHorizontalOrder != 0)
			) {
				_EXCEPTIONT("Elements must be aligned with HorizontalOrder");
			}

			double dAInteriorBegin =
				m_dGDim[0] + iAElementInteriorBegin * dElementDeltaA;
			double dAInteriorEnd =
				m_dGDim[0] + iAElementInteriorEnd * dElementDeltaA;
			double dBInteriorBegin =
				m_dGDim[2] + iBElementInteriorBegin * dElementDeltaB;
			double dBInteriorEnd =
				m_dGDim[2] + iBElementInteriorEnd * dElementDeltaB;

			if ((dAlpha[i] >= dAInteriorBegin) &&
				(dAlpha[i] <= dAInteriorEnd) &&
				(dBeta[i] >= dBInteriorBegin) &&
				(dBeta[i] <= dBInteriorEnd)
			) {
				iPatch[i] = n;
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
	const DataArray1D<int> & vecIxA,
	const DataArray1D<int> & vecIxB,
	const DataArray1D<int> & vecPanel,
	DataArray1D<int> & vecPatchIndex,
	int nVectorLength
) {
	// Set vector length
	if (nVectorLength == (-1)) {
		nVectorLength = vecIxA.GetRows();
	}

//std::cout << GetPatchCount() << '\n';

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
		int iP = 0;

		// Wrap global indices
		if (iA < 0) {
			BoundaryCondition eLeftBoundary =
				m_eBoundaryCondition[Direction_Left];

			if (eLeftBoundary == BoundaryCondition_Periodic) {
				iA += nLocalResolutionA;
			} else {
				vecPatchIndex[i] = GridPatch::InvalidIndex;
				continue;
			}
		}
		if (iA >= nLocalResolutionA) {
			BoundaryCondition eRightBoundary =
				m_eBoundaryCondition[Direction_Right];

			if (eRightBoundary == BoundaryCondition_Periodic) {
				iA -= nLocalResolutionA;
			} else {
				vecPatchIndex[i] = GridPatch::InvalidIndex;
				continue;
			}
		}
		if (iB < 0) {
			BoundaryCondition eBottomBoundary =
				m_eBoundaryCondition[Direction_Bottom];

			if (eBottomBoundary == BoundaryCondition_Periodic) {
				iB += nLocalResolutionB;
			} else {
				vecPatchIndex[i] = GridPatch::InvalidIndex;
				continue;
			}
		}
		if (iB >= nLocalResolutionB) {
			BoundaryCondition eTopBoundary =
				m_eBoundaryCondition[Direction_Top];

			if (eTopBoundary == BoundaryCondition_Periodic) {
				iB -= nLocalResolutionB;
			} else {
				vecPatchIndex[i] = GridPatch::InvalidIndex;
				continue;
			}
		}

		// Check the last patch searched
		if (iLastPatch != GridPatch::InvalidIndex) {
			const PatchBox & box = m_aPatchBoxes[iLastPatch];

			if (box.ContainsGlobalPoint(iP, iA, iB)) {
				vecPatchIndex[i] = iLastPatch;
				continue;
			}
		}

		// Check all other patches
		int n;
		for (n = 0; n < GetPatchCount(); n++) {
			const PatchBox & box = m_aPatchBoxes[n];

			if (box.ContainsGlobalPoint(iP, iA, iB)) {
				vecPatchIndex[i] = n;
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
		if (eDataType == DataType_TopographyDeriv) {
			pPatch->TransformTopographyDeriv();
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
		} else if (eDataType == DataType_TopographyDeriv) {
			nComponents = 2;
		} else {
			_EXCEPTIONT("Invalid DataType");
		}

		// Apply BC only to state DSS
		if (eDataType == DataType_State) {
			pPatch->ApplyBoundaryConditions(iDataUpdate, DataType_State, n);
		}

		// Perform Direct Stiffness Summation (DSS)
		for (int c = 0; c < nComponents; c++) {

			// Obtain the array of working data
			int nRElements = GetRElements();

			DataArray3D<double> pDataUpdate;

			if ((eDataType == DataType_State) &&
				(GetVarLocation(c) == DataLocation_REdge)
			) {
				nRElements++;
			}
			if (eDataType == DataType_TopographyDeriv) {
				nRElements = 2;
			}

			pDataUpdate.SetSize(
				box.GetATotalWidth(),
				box.GetBTotalWidth(),
				nRElements);

			// State data
			if (eDataType == DataType_State) {
				DataArray4D<double> & dState =
					pPatch->GetDataState(iDataUpdate, GetVarLocation(c));

				pDataUpdate.AttachToData(&(dState[c][0][0][0]));

			// Tracer data
			} else if (eDataType == DataType_Tracers) {
				DataArray4D<double> & dTracers =
					pPatch->GetDataTracers(iDataUpdate);

				pDataUpdate.AttachToData(&(dTracers[c][0][0][0]));

			// Vorticity data
			} else if (eDataType == DataType_Vorticity) {
				DataArray3D<double> & dVorticity =
					pPatch->GetDataVorticity();

				pDataUpdate.AttachToData(&(dVorticity[0][0][0]));

			// Divergence data
			} else if (eDataType == DataType_Divergence) {
				DataArray3D<double> & dDivergence =
					pPatch->GetDataDivergence();

				pDataUpdate.AttachToData(&(dDivergence[0][0][0]));

			// Topographic derivative data
			} else if (eDataType == DataType_TopographyDeriv) {
				DataArray3D<double> & dTopographyDeriv =
					pPatch->GetTopographyDeriv();

				pDataUpdate.AttachToData(&(dTopographyDeriv[0][0][0]));
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
						pDataUpdate[iA][j][k] = 0.5 * (
							+ pDataUpdate[iA  ][j][k]
							+ pDataUpdate[iA-1][j][k]);

						pDataUpdate[iA-1][j][k] = pDataUpdate[iA][j][k];
					}
				}

				// Average in the beta direction
				for (int b = 0; b <= nElementCountB; b++) {
					int iB = b * m_nHorizontalOrder + box.GetHaloElements();

					// Averaging done at the corners of the panel
					int iBegin = box.GetAInteriorBegin()-1;
					int iEnd = box.GetAInteriorEnd()+1;

					for (int i = iBegin; i < iEnd; i++) {
						pDataUpdate[i][iB][k] = 0.5 * (
							+ pDataUpdate[i][iB  ][k]
							+ pDataUpdate[i][iB-1][k]);

						pDataUpdate[i][iB-1][k] = pDataUpdate[i][iB][k];
					}
				}
			}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////
