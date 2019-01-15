///////////////////////////////////////////////////////////////////////////////
///
///	\file    GridCSGLL.cpp
///	\author  Paul Ullrich
///	\version September 19, 2013
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

#include "GridCSGLL.h"
#include "GridPatchCSGLL.h"
#include "Model.h"
#include "TestCase.h"
#include "GridSpacing.h"

#include "Direction.h"
#include "CubedSphereTrans.h"
#include "PolynomialInterp.h"

#include "Announce.h"
#include "MathHelper.h"

#include <cmath>

///////////////////////////////////////////////////////////////////////////////

GridCSGLL::GridCSGLL(
	Model & model
) :
	GridGLL::GridGLL(model)
{ }

///////////////////////////////////////////////////////////////////////////////

void GridCSGLL::DefineParameters() {

	if (m_dcGridParameters.IsAttached()) {
		_EXCEPTIONT("Attempting to recall DefineParameters");
	}

	// Call up the stack
	GridGLL::DefineParameters();
}

///////////////////////////////////////////////////////////////////////////////

void GridCSGLL::SetParameters(
	int nRElements,
	int nMaxPatchCount,
	int nBaseResolution,
	int nRefinementRatio,
	int nHorizontalOrder,
	int nVerticalOrder,
	VerticalDiscretization eVerticalDiscretization,
	VerticalStaggering eVerticalStaggering
) {
	if (!m_dcGridParameters.IsAttached()) {
		_EXCEPTIONT("DefineParameters() must be called before SetParameters()");
	}

	// At least 6 patches needed
	if (nMaxPatchCount < 6) {
		_EXCEPTIONT("nMinPatchCount must be >= 6");
	}

	// Call up the stack
	GridGLL::SetParameters(
		nRElements,
		nMaxPatchCount,
		nBaseResolution,
		nBaseResolution,
		nRefinementRatio,
		nHorizontalOrder,
		nVerticalOrder,
		eVerticalDiscretization,
		eVerticalStaggering
	);

	// Set the reference length scale
	m_dReferenceLength = 0.5 * M_PI / 30.0;
}

///////////////////////////////////////////////////////////////////////////////

void GridCSGLL::ApplyDefaultPatchLayout(
	int nPatchCount
) {

	// Verify no Patches have been previously added
	if (m_nInitializedPatchBoxes != 0) {
		_EXCEPTIONT("ApplyDefaultPatchLayout() must be called on an empty Grid");
	}

	int nProcsPerDirection = Max((int)ISqrt(nPatchCount / 6), 1);

	int nProcsPerPanel = nProcsPerDirection * nProcsPerDirection;

	int nDistributedPatches = 6 * nProcsPerPanel;

	if (nDistributedPatches < nPatchCount) {
		Announce("WARNING: Patch / thread mismatch: "
			"%i threads will be unutilized",
			nPatchCount - nDistributedPatches);
	}

	// Determine arrangement of elements on processors
	if (GetABaseResolution() % nProcsPerDirection != 0) {
		_EXCEPTIONT("\n(UNIMPLEMENTED) Currently elements must be "
			"equally divided among processors.");
	}

	int nElementsPerDirection = GetABaseResolution() / nProcsPerDirection;
	DataArray1D<int> iBoxBegin(nProcsPerDirection + 1);

	iBoxBegin[0] = 0;
	for (int n = 1; n < nProcsPerDirection; n++) {
		iBoxBegin[n] = n * nElementsPerDirection;
	}
	iBoxBegin[nProcsPerDirection] = GetABaseResolution();

	// Create master patch for each panel
	int ixPatch = 0;

	for (int n = 0; n < 6; n++) {
	for (int i = 0; i < nProcsPerDirection; i++) {
	for (int j = 0; j < nProcsPerDirection; j++) {

		m_aPatchBoxes[ixPatch] = PatchBox(
			n, 0, m_model.GetHaloElements(),
			m_nHorizontalOrder * iBoxBegin[i],
			m_nHorizontalOrder * iBoxBegin[i+1],
			m_nHorizontalOrder * iBoxBegin[j],
			m_nHorizontalOrder * iBoxBegin[j+1]);

		ixPatch++;
	}
	}
	}

	m_nInitializedPatchBoxes = ixPatch;
}

///////////////////////////////////////////////////////////////////////////////

GridPatch * GridCSGLL::NewPatch(
	int ixPatch
) {
	const PatchBox & box = GetPatchBox(ixPatch);

	return (
		new GridPatchCSGLL(
			(*this),
			ixPatch,
			box,
			m_nHorizontalOrder,
			m_nVerticalOrder));
}

///////////////////////////////////////////////////////////////////////////////

void GridCSGLL::GetReferenceGridBounds(
	double & dX0,
	double & dX1,
	double & dY0,
	double & dY1
) {
	dX0 = - 180.0;
	dX1 = + 180.0;
	dY0 = -  90.0;
	dY1 = +  90.0;
}

///////////////////////////////////////////////////////////////////////////////

void GridCSGLL::ConvertReferenceToPatchCoord(
	const DataArray1D<double> & dXReference,
	const DataArray1D<double> & dYReference,
	DataArray1D<double> & dAlpha,
	DataArray1D<double> & dBeta,
	DataArray1D<int> & iPatch
) const {
	if ((dXReference.GetRows() != dYReference.GetRows()) ||
		(dXReference.GetRows() != dAlpha.GetRows()) ||
		(dXReference.GetRows() != dBeta.GetRows()) ||
		(dXReference.GetRows() != iPatch.GetRows())
	) {
		_EXCEPTIONT("Dimension mismatch: All arrays must have same length");
	}

	// Loop over all coordinates
	for (int i = 0; i < dXReference.GetRows(); i++) {

		int iPanel;

		CubedSphereTrans::ABPFromRLL(
			dXReference[i] * M_PI / 180.0,
			dYReference[i] * M_PI / 180.0,
			dAlpha[i],
			dBeta[i],
			iPanel);

		// Loop over all patches
		int n = 0;

		for (; n < GetPatchCount(); n++) {
			const PatchBox & box = m_aPatchBoxes[n];

			if (iPanel != box.GetPanel()) {
				continue;
			}

			double dElementDeltaA =
				0.5 * M_PI / static_cast<double>(GetABaseResolution());
			double dElementDeltaB =
				0.5 * M_PI / static_cast<double>(GetBBaseResolution());

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
				-0.25 * M_PI + iAElementInteriorBegin * dElementDeltaA;
			double dAInteriorEnd =
				-0.25 * M_PI + iAElementInteriorEnd * dElementDeltaA;
			double dBInteriorBegin =
				-0.25 * M_PI + iBElementInteriorBegin * dElementDeltaB;
			double dBInteriorEnd =
				-0.25 * M_PI + iBElementInteriorEnd * dElementDeltaB;

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
			_EXCEPTION5("Unable to find associated patch for node:\n"
				"(%1.5e, %1.5e) : (%1.5e, %1.5e, %i)",
				dXReference[i], dYReference[i],
				dAlpha[i], dBeta[i], iPanel);
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void GridCSGLL::GetPatchFromCoordinateIndex(
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
		if ((vecIxA[i] < 0) ||
			(vecIxA[i] >= nLocalResolution) ||
			(vecIxB[i] < 0) ||
			(vecIxB[i] >= nLocalResolution)
		) {
			bool fSwitchAB;
			bool fSwitchPar;
			bool fSwitchPerp;

			CubedSphereTrans::RelativeCoord(
				nLocalResolution,
				vecPanel[i],
				vecIxA[i],
				vecIxB[i],
				iP, iA, iB,
				fSwitchAB,
				fSwitchPar,
				fSwitchPerp);

			// Coordinate out of bounds
			if ((iA == (-1)) && (iB == (-1))) {
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
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void GridCSGLL::GetOpposingDirection(
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

	// Get the opposing direction
	dirOpposing =
		CubedSphereTrans::OpposingDirection(
			ixPanelSrc, ixPanelDest, dir);

	// Check panel indices to determine parallel switch of vectors
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

	// Check panel indices to determine perpendicular switch of vectors
	fSwitchPerpendicular = false;
	if ((ixPanelSrc == 1) && (ixPanelDest == 4)) {
		fSwitchPerpendicular = true;

	} else if ((ixPanelSrc == 2) && (ixPanelDest > 3)) {
		fSwitchPerpendicular = true;

	} else if ((ixPanelSrc == 3) && (ixPanelDest == 5)) {
		fSwitchPerpendicular = true;

	} else if (
		(ixPanelSrc == 4) &&
		((ixPanelDest == 1) || (ixPanelDest == 2))
	) {
		fSwitchPerpendicular = true;

	} else if (
		(ixPanelSrc == 5) &&
		((ixPanelDest == 3) || (ixPanelDest == 2))
	) {
		fSwitchPerpendicular = true;
	}
}

///////////////////////////////////////////////////////////////////////////////

void GridCSGLL::ApplyDSS(
	int iDataUpdate,
	DataType eDataType
) {
	// Exchange data between nodes
	Exchange(eDataType, iDataUpdate);

	// Post-process velocities across panel edges and
	// perform direct stiffness summation (DSS)
	for (int n = 0; n < GetActivePatchCount(); n++) {
		GridPatchCSGLL * pPatch =
			dynamic_cast<GridPatchCSGLL*>(GetActivePatch(n));

		const PatchBox & box = pPatch->GetPatchBox();
	
	// Patch-specific quantities
		const int nElementCountA = pPatch->GetElementCountA();
		const int nElementCountB = pPatch->GetElementCountB();
	
		// Apply panel transforms to velocity data
		if (eDataType == DataType_State) {
			pPatch->TransformHaloVelocities(iDataUpdate);
		}
		if (eDataType == DataType_TopographyDeriv) {
			pPatch->TransformTopographyDeriv();

			const DataArray3D<double> & dataTopographyDeriv =
				pPatch->GetTopographyDeriv();
		}

		// Panels in each coordinate direction
		const int ixRightPanel =
			pPatch->GetNeighborPanel(Direction_Right);
		const int ixTopPanel =
			pPatch->GetNeighborPanel(Direction_Top);
		const int ixLeftPanel =
			pPatch->GetNeighborPanel(Direction_Left);
		const int ixBottomPanel =
			pPatch->GetNeighborPanel(Direction_Bottom);

		const int ixTopRightPanel =
			pPatch->GetNeighborPanel(Direction_TopRight);
		const int ixTopLeftPanel =
			pPatch->GetNeighborPanel(Direction_TopLeft);
		const int ixBottomLeftPanel =
			pPatch->GetNeighborPanel(Direction_BottomLeft);
		const int ixBottomRightPanel =
			pPatch->GetNeighborPanel(Direction_BottomRight);

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
			nComponents = 1;
		} else {
			_EXCEPTIONT("Invalid DataType");
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

				pDataUpdate.AttachToData(&(dState(c,0,0,0)));

			// Tracer data
			} else if (eDataType == DataType_Tracers) {
				DataArray4D<double> & dTracers =
					pPatch->GetDataTracers(iDataUpdate);

				pDataUpdate.AttachToData(&(dTracers(c,0,0,0)));

			// Vorticity data
			} else if (eDataType == DataType_Vorticity) {
				DataArray3D<double> & dVorticity =
					pPatch->GetDataVorticity();

				pDataUpdate.AttachToData(&(dVorticity(0,0,0)));

			// Divergence data
			} else if (eDataType == DataType_Divergence) {
				DataArray3D<double> & dDivergence =
					pPatch->GetDataDivergence();

				pDataUpdate.AttachToData(&(dDivergence(0,0,0)));

			// Topographic derivative data
			} else if (eDataType == DataType_TopographyDeriv) {
				DataArray3D<double> & dTopographyDeriv =
					pPatch->GetTopographyDeriv();

				pDataUpdate.AttachToData(&(dTopographyDeriv(0,0,0)));
			}

			/// start of loop nest 1

			// Average in the alpha direction (patch interior)
			{
				const int jBegin = box.GetBInteriorBegin()-1;
				const int jEnd = box.GetBInteriorEnd()+1;

				// Perform averaging across edge
				for (int a = 1; a < nElementCountA; a++) {
					const int iA = a * m_nHorizontalOrder + box.GetHaloElements();
				    for (int j = jBegin; j < jEnd; j++) {
#pragma simd
						for (int k = 0; k < nRElements; k++) {
						    pDataUpdate(iA,j,k) = 0.5 * (
								pDataUpdate(iA,j,k)
								+ pDataUpdate(iA-1,j,k));

						    pDataUpdate(iA-1,j,k) =
								pDataUpdate(iA,j,k);
						}
				    }
				}
			}

			// Average in the alpha direction (left edge of patch)
			{
				const int iA = box.GetHaloElements();

				// Do not average across cubed-sphere corners
				int jBegin = box.GetBInteriorBegin()-1;
				int jEnd = box.GetBInteriorEnd()+1;

				if (ixTopLeftPanel  == InvalidPanel) {
					jEnd -= 2;
				}
				if (ixBottomLeftPanel == InvalidPanel) {
					jBegin += 2;
				}

				// Perform averaging across edge
				for (int j = jBegin; j < jEnd; j++) {
#pragma simd
				    for (int k = 0; k < nRElements; k++) {
						pDataUpdate(iA,j,k) = 0.5 * (
							pDataUpdate(iA,j,k)
							+ pDataUpdate(iA-1,j,k));

						pDataUpdate(iA-1,j,k) =
							pDataUpdate(iA,j,k);
					}
				}
			}

			// Average in the alpha direction (right edge of patch)
			{
				const int iA =
					nElementCountA * m_nHorizontalOrder
					+ box.GetHaloElements();

				// Do not average across cubed-sphere corners
				int jBegin = box.GetBInteriorBegin()-1;
				int jEnd = box.GetBInteriorEnd()+1;

				if (ixTopRightPanel == InvalidPanel) {
					jEnd -= 2;
				}
				if (ixBottomRightPanel == InvalidPanel) {
					jBegin += 2;
				}

			    // Perform averaging across edge
				for (int j = jBegin; j < jEnd; j++) {
#pragma simd
				    for (int k = 0; k < nRElements; k++) {
						pDataUpdate(iA,j,k) = 0.5 * (
							pDataUpdate(iA,j,k)
							+ pDataUpdate(iA-1,j,k));

						pDataUpdate(iA-1,j,k) =
							pDataUpdate(iA,j,k);
					}
				}
			}
			/// end of loop 1

			/// start of loop nest 2

			// Average in the beta direction (patch interior)
			{
				const int iBegin = box.GetAInteriorBegin()-1;
				const int iEnd = box.GetAInteriorEnd()+1;
				for (int b = 1; b < nElementCountB; b++) {
					int iB = b * m_nHorizontalOrder + box.GetHaloElements();

					for (int i = iBegin; i < iEnd; i++) {
#pragma simd
						for (int k = 0; k < nRElements; k++) {
							pDataUpdate(i,iB,k) = 0.5 * (
								pDataUpdate(i,iB,k)
								+ pDataUpdate(i,iB-1,k));

							pDataUpdate(i,iB-1,k) =
								pDataUpdate(i,iB,k);
						}
					}
				}
			}

			// Average in the beta direction (bottom edge of patch)
			{
				const int iB = box.GetHaloElements();

				int iBegin = box.GetAInteriorBegin()-1;
				int iEnd = box.GetAInteriorEnd()+1;

				if (ixBottomLeftPanel == InvalidPanel) {
					iBegin += 2;
				}
				if (ixBottomRightPanel == InvalidPanel) {
					iEnd -= 2;
				}

				for (int i = iBegin; i < iEnd; i++) {
#pragma simd
					for (int k = 0; k < nRElements; k++) {
						pDataUpdate(i,iB,k) = 0.5 * (
							pDataUpdate(i,iB,k)
							+ pDataUpdate(i,iB-1,k));
						pDataUpdate(i,iB-1,k) =
							pDataUpdate(i,iB,k);
					}
				}
			}

			// Average in the beta direction (top edge of patch)
			{
				const int iB =
					nElementCountB * m_nHorizontalOrder
					+ box.GetHaloElements();

				int iBegin = box.GetAInteriorBegin()-1;
				int iEnd = box.GetAInteriorEnd()+1;

				if (ixTopLeftPanel == InvalidPanel) {
					iBegin += 2;
				}
				if (ixTopRightPanel == InvalidPanel) {
					iEnd -= 2;
				}

				for (int i = iBegin; i < iEnd; i++) {
#pragma simd
					for (int k = 0; k < nRElements; k++) {
						pDataUpdate(i,iB,k) = 0.5 * (
							pDataUpdate(i,iB,k)
							+ pDataUpdate(i,iB-1,k));

						pDataUpdate(i,iB-1,k) =
							pDataUpdate(i,iB,k);
					}
				}
			}

			/// end of loop 2

			/// start of loop nest 3

			// Handle patch corners that are coincident with panel corners.
			// These are also nodes of connectivity 3.
			if (ixTopRightPanel == InvalidPanel) {
				const int iA = box.GetAInteriorEnd()-1;
				const int iB = box.GetBInteriorEnd()-1;

#pragma simd
				for (int k = 0; k < nRElements; k++) {
					pDataUpdate(iA,iB,k) = (1.0/3.0) * (
						+ pDataUpdate(iA  ,iB  ,k)
						+ pDataUpdate(iA+1,iB  ,k)
						+ pDataUpdate(iA  ,iB+1,k));
				}
			}

			if (ixTopLeftPanel == InvalidPanel) {
				const int iA = box.GetAInteriorBegin();
				const int iB = box.GetBInteriorEnd()-1;

#pragma simd
				for (int k = 0; k < nRElements; k++) {
					pDataUpdate(iA,iB,k) = (1.0/3.0) * (
						+ pDataUpdate(iA  ,iB  ,k)
						+ pDataUpdate(iA-1,iB  ,k)
						+ pDataUpdate(iA  ,iB+1,k));
				}
			}

			if (ixBottomLeftPanel == InvalidPanel) {
				const int iA = box.GetAInteriorBegin();
				const int iB = box.GetBInteriorBegin();

#pragma simd
				for (int k = 0; k < nRElements; k++) {
					pDataUpdate(iA,iB,k) = (1.0/3.0) * (
						+ pDataUpdate(iA  ,iB  ,k)
						+ pDataUpdate(iA-1,iB  ,k)
						+ pDataUpdate(iA  ,iB-1,k));
				}
			}

			if (ixBottomRightPanel == InvalidPanel) {
				const int iA = box.GetAInteriorEnd()-1;
				const int iB = box.GetBInteriorBegin();

#pragma simd
				for (int k = 0; k < nRElements; k++) {
					pDataUpdate(iA,iB,k) = (1.0/3.0) * (
						+ pDataUpdate(iA  ,iB  ,k)
						+ pDataUpdate(iA+1,iB  ,k)
						+ pDataUpdate(iA  ,iB-1,k));
				}
			}
			/// end of loop 3
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

