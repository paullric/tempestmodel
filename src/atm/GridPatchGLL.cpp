///////////////////////////////////////////////////////////////////////////////
///
///	\file    GridPatchGLL.cpp
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

#include "GridPatchGLL.h"
#include "Grid.h"
#include "Model.h"
#include "EquationSet.h"
#include "Defines.h"
#include "DataArray1D.h"

///////////////////////////////////////////////////////////////////////////////

GridPatchGLL::GridPatchGLL(
	Grid & grid,
	int ixPatch,
	const PatchBox & box,
	int nHorizontalOrder,
	int nVerticalOrder
) :
	GridPatch(grid, ixPatch, box),
	m_nHorizontalOrder(nHorizontalOrder),
	m_nVerticalOrder(nVerticalOrder)
{
	// Verify that box boundaries are aligned with elements
	if (((box.GetAGlobalInteriorBegin() % nHorizontalOrder) != 0) ||
		((box.GetAGlobalInteriorEnd()   % nHorizontalOrder) != 0) ||
		((box.GetBGlobalInteriorBegin() % nHorizontalOrder) != 0) ||
		((box.GetBGlobalInteriorEnd()   % nHorizontalOrder) != 0)
	) {
		_EXCEPTION4(
			"GLL grid patch must be aligned with elements: %i %i %i %i",
			box.GetAGlobalInteriorBegin(),
			box.GetAGlobalInteriorEnd(),
			box.GetBGlobalInteriorBegin(),
			box.GetBGlobalInteriorEnd());
	}

	// Get the number of finite elements in each coordinate direction
	m_nElementCountA =
		(box.GetAInteriorEnd() - box.GetAInteriorBegin()) / nHorizontalOrder;
	m_nElementCountB =
		(box.GetBInteriorEnd() - box.GetBInteriorBegin()) / nHorizontalOrder;

	if ((box.GetAInteriorWidth() % m_nHorizontalOrder) != 0) {
		_EXCEPTIONT("Logic error: Invalid PatchBox alpha spacing");
	}
	if ((box.GetBInteriorWidth() % m_nHorizontalOrder) != 0) {
		_EXCEPTIONT("Logic error: Invalid PatchBox beta spacing");
	}
}

///////////////////////////////////////////////////////////////////////////////

void GridPatchGLL::InitializeCoordinateData() {

	m_dElementDeltaA =
		  m_dAEdge[m_box.GetHaloElements() + m_nHorizontalOrder]
		- m_dAEdge[m_box.GetHaloElements()];

	m_dElementDeltaB =
		  m_dBEdge[m_box.GetHaloElements() + m_nHorizontalOrder]
		- m_dBEdge[m_box.GetHaloElements()];
}

///////////////////////////////////////////////////////////////////////////////

void GridPatchGLL::InterpolateNodeToREdge(
	int iVar,
	int iDataIndex
) {

	// Working data
	DataArray4D<double> & dataNode =
		GetDataState(iDataIndex, DataLocation_Node);
	DataArray4D<double> & dataREdge =
		GetDataState(iDataIndex, DataLocation_REdge);

	// Parent grid, containing the vertical remapping information
	GridGLL * pGLLGrid = dynamic_cast<GridGLL*>(&m_grid);
	if (pGLLGrid == NULL) {
		_EXCEPTIONT("Logic error");
	}

	// Loop over all elements in the box
	int nStride = dataNode.GetSize(2) * dataNode.GetSize(3);

	const LinearColumnInterpFEM & opInterpNodeToREdge =
		pGLLGrid->GetOpInterpNodeToREdge();

	for (int i = m_box.GetAInteriorBegin(); i < m_box.GetAInteriorEnd(); i++) {
	for (int j = m_box.GetBInteriorBegin(); j < m_box.GetBInteriorEnd(); j++) {

		opInterpNodeToREdge.Apply(
			&(dataNode[iVar][0][i][j]),
			&(dataREdge[iVar][0][i][j]),
			nStride,
			nStride);
	}
	}
}

///////////////////////////////////////////////////////////////////////////////

void GridPatchGLL::InterpolateREdgeToNode(
	int iVar,
	int iDataIndex
) {

	// Working data
	DataArray4D<double> & dataREdge =
		GetDataState(iDataIndex, DataLocation_REdge);
	DataArray4D<double> & dataNode =
		GetDataState(iDataIndex, DataLocation_Node);

	// Parent grid, containing the vertical remapping information
	GridGLL * pGLLGrid = dynamic_cast<GridGLL*>(&m_grid);
	if (pGLLGrid == NULL) {
		_EXCEPTIONT("Logic error");
	}

	// Loop over all elements in the box
	int nStride = dataNode.GetSize(2) * dataNode.GetSize(3);

	const LinearColumnInterpFEM & opInterpREdgeToNode =
		pGLLGrid->GetOpInterpREdgeToNode();

	// Loop over all elements in the box
	for (int i = m_box.GetAInteriorBegin(); i < m_box.GetAInteriorEnd(); i++) {
	for (int j = m_box.GetBInteriorBegin(); j < m_box.GetBInteriorEnd(); j++) {

		opInterpREdgeToNode.Apply(
			&(dataREdge[iVar][0][i][j]),
			&(dataNode[iVar][0][i][j]),
			nStride,
			nStride);
	}
	}
}

///////////////////////////////////////////////////////////////////////////////

void GridPatchGLL::ComputeRichardson(
	int iDataIndex,
	DataLocation loc
) {
	const PhysicalConstants & phys = m_grid.GetModel().GetPhysicalConstants();

	// Number of radial elements
	int nRElements = m_grid.GetRElements();

	// Get metric quantities
	const DataArray4D<double> & dDerivRNode =
		GetDerivRNode();
	const DataArray4D<double> & dDerivRREdge =
		GetDerivRREdge();

	// Indices of EquationSet variables
	const int UIx = 0;
	const int VIx = 1;
	const int PIx = 2;
	const int WIx = 3;
	const int RIx = 4;

// Column array of density and density gradient
	DataArray1D<double> dDensityNode;
	DataArray1D<double> dDensityREdge;
	DataArray1D<double> dDiffDensityNode;
	DataArray1D<double> dDiffDensityREdge;
	dDensityNode.Allocate(nRElements);
	dDiffDensityNode.Allocate(nRElements);
	dDensityREdge.Allocate(nRElements+1);
	dDiffDensityREdge.Allocate(nRElements+1);
// Column array of U_x and shear in the zonal wind
	DataArray1D<double> dUXNode;
	DataArray1D<double> dUXREdge;
	DataArray1D<double> dDiffUXNode;
	DataArray1D<double> dDiffUXREdge;
	dUXNode.Allocate(nRElements);
	dDiffUXNode.Allocate(nRElements);
	dUXREdge.Allocate(nRElements+1);
	dDiffUXREdge.Allocate(nRElements+1);
// Column array of V_y and shear in the meridional wind
	DataArray1D<double> dVYNode;
	DataArray1D<double> dVYREdge;
	DataArray1D<double> dDiffVYNode;
	DataArray1D<double> dDiffVYREdge;
	dVYNode.Allocate(nRElements);
	dDiffVYNode.Allocate(nRElements);
	dVYREdge.Allocate(nRElements+1);
	dDiffVYREdge.Allocate(nRElements+1);

	// Under this configuration, set fluxes at boundaries to zero
	bool fZeroBoundaries =
		(m_grid.GetVarLocation(WIx) == DataLocation_REdge);

	if (m_grid.GetModel().GetEquationSet().GetComponents() < 5) {
		_EXCEPTIONT("Invalid EquationSet.");
	}

	// Parent grid, containing the vertical remapping information
	GridGLL * pGLLGrid = dynamic_cast<GridGLL*>(&m_grid);
	if (pGLLGrid == NULL) {
		_EXCEPTIONT("Logic error");
	}

	int k;
	int i;
	int j;

	// Calculate Richardson number on nodes
	if (loc == DataLocation_Node) {
		if (m_grid.GetVarLocation(PIx) == DataLocation_REdge) {
			InterpolateREdgeToNode(PIx, iDataIndex);
		}
		if (m_grid.GetVarLocation(RIx) == DataLocation_REdge) {
			InterpolateREdgeToNode(RIx, iDataIndex);
		}
		if (m_grid.GetVarLocation(WIx) == DataLocation_REdge) {
			InterpolateREdgeToNode(WIx, iDataIndex);
		}

		const DataArray4D<double> & dataNode = m_datavecStateNode[iDataIndex];

		for (i = m_box.GetAInteriorBegin(); i < m_box.GetAInteriorEnd(); i++) {
		for (j = m_box.GetBInteriorBegin(); j < m_box.GetBInteriorEnd(); j++) {

			for (k = 0; k < nRElements; k++) {
				dDensityNode[k] = dataNode[RIx][k][i][j];

				// Convert U_alpha and V_beta to X and Y (Lon, Lat)
				dUXNode[k] = dataNode[UIx][k][i][j] - 
							dDerivRNode[k][i][j][0] * dataNode[WIx][k][i][j];

				dVYNode[k] = dataNode[VIx][k][i][j] - 
							dDerivRNode[k][i][j][1] * dataNode[WIx][k][i][j];
			}

			pGLLGrid->DifferentiateNodeToNode(
			dDensityNode,
			dDiffDensityNode,
			fZeroBoundaries);

			pGLLGrid->DifferentiateNodeToNode(
			dUXNode,
			dDiffUXNode,
			fZeroBoundaries);

			pGLLGrid->DifferentiateNodeToNode(
			dVYNode,
			dDiffVYNode,
			fZeroBoundaries);

			for (k = 0; k < nRElements; k++) {
				//[k][i][j] = dDiffDensityNode[k];
				//m_dataRichardson[k][i][j] = dDiffUXNode[k];
				m_dataRichardson[k][i][j] = phys.GetG() / dDensityNode[k] * 
				dDiffDensityNode[k] / 
				((dDiffUXNode[k] * dDiffUXNode[k]) + 
				 (dDiffVYNode[k] * dDiffVYNode[k]));

//printf("%i %.16E \n",k,m_dataRichardson[k][i][j]);
			}
		}
		}
	}

	// Calculate Richardson on interfaces
	if (loc == DataLocation_REdge) {
		//_EXCEPTIONT("Richardson number not implemented on interfaces");

		if (m_grid.GetVarLocation(PIx) == DataLocation_Node) {
			InterpolateNodeToREdge(PIx, iDataIndex);
		}
		if (m_grid.GetVarLocation(RIx) == DataLocation_Node) {
			InterpolateNodeToREdge(RIx, iDataIndex);
		}
		if (m_grid.GetVarLocation(WIx) == DataLocation_Node) {
			InterpolateNodeToREdge(WIx, iDataIndex);
		}

		const DataArray4D<double> & dataREdge = m_datavecStateREdge[iDataIndex];

		for (i = m_box.GetAInteriorBegin(); i < m_box.GetAInteriorEnd(); i++) {
		for (j = m_box.GetBInteriorBegin(); j < m_box.GetBInteriorEnd(); j++) {

			for (k = 0; k <= nRElements; k++) {
				dDensityREdge[k] = dataREdge[RIx][k][i][j];

				// Convert U_alpha and V_beta to X and Y (Lon, Lat)
				dUXREdge[k] = dataREdge[UIx][k][i][j] - 
							dDerivRREdge[k][i][j][0] * dataREdge[WIx][k][i][j];

				dVYREdge[k] = dataREdge[VIx][k][i][j] - 
							dDerivRREdge[k][i][j][1] * dataREdge[WIx][k][i][j];
			}

			pGLLGrid->DifferentiateREdgeToREdge(
			dDensityREdge,
			dDiffDensityREdge);

			pGLLGrid->DifferentiateREdgeToREdge(
			dUXREdge,
			dDiffUXREdge);

			pGLLGrid->DifferentiateREdgeToREdge(
			dVYREdge,
			dDiffVYREdge);

			for (k = 0; k <= nRElements; k++) {
				//m_dataRichardson[k][i][j] = dDiffUXREdge[k] * dDiffUXREdge[k];
				//m_dataRichardson[k][i][j] = dDiffDensityREdge[k] / dDensityREdge[k];
				m_dataRichardson[k][i][j] = phys.GetG() * 
				(dDiffDensityREdge[k] / dDensityREdge[k]) / 
				((dDiffUXREdge[k] * dDiffUXREdge[k]) + 
				 (dDiffVYREdge[k] * dDiffVYREdge[k]));
			}
		}
		}
	}
}

