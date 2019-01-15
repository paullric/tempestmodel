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
	const LinearColumnInterpFEM & opInterpNodeToREdge =
		pGLLGrid->GetOpInterpNodeToREdge();

	for (int i = m_box.GetAInteriorBegin(); i < m_box.GetAInteriorEnd(); i++) {
	for (int j = m_box.GetBInteriorBegin(); j < m_box.GetBInteriorEnd(); j++) {

		opInterpNodeToREdge.Apply(
			&(dataNode[iVar][i][j][0]),
			&(dataREdge[iVar][i][j][0]));
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
	const LinearColumnInterpFEM & opInterpREdgeToNode =
		pGLLGrid->GetOpInterpREdgeToNode();

	// Loop over all elements in the box
	for (int i = m_box.GetAInteriorBegin(); i < m_box.GetAInteriorEnd(); i++) {
	for (int j = m_box.GetBInteriorBegin(); j < m_box.GetBInteriorEnd(); j++) {

		opInterpREdgeToNode.Apply(
			&(dataREdge[iVar][i][j][0]),
			&(dataNode[iVar][i][j][0]));
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

	// Indices of EquationSet variables
	const int UIx = 0;
	const int VIx = 1;
	const int PIx = 2;
	const int WIx = 3;
	const int RIx = 4;

	// Column array of density and density gradient
	DataArray1D<double> dDensityNode;
	DataArray1D<double> dDiffDensityNode;
	dDensityNode.Allocate(nRElements);
	dDiffDensityNode.Allocate(nRElements);
	// Column array of U_x and shear in the zonal wind
	DataArray1D<double> dUXNode;
	DataArray1D<double> dDiffUXNode;
	dUXNode.Allocate(nRElements);
	dDiffUXNode.Allocate(nRElements);
	// Column array of V_y and shear in the meridional wind
	DataArray1D<double> dVYNode;
	DataArray1D<double> dDiffVYNode;
	dVYNode.Allocate(nRElements);
	dDiffVYNode.Allocate(nRElements);

	if (m_grid.GetModel().GetEquationSet().GetComponents() < 5) {
		_EXCEPTIONT("Invalid EquationSet.");
	}

	// Parent grid, containing the vertical remapping information
	GridGLL * pGLLGrid = dynamic_cast<GridGLL*>(&m_grid);
	if (pGLLGrid == NULL) {
		_EXCEPTIONT("Logic error");
	}

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

		for (int i = m_box.GetAInteriorBegin(); i < m_box.GetAInteriorEnd(); i++) {
		for (int j = m_box.GetBInteriorBegin(); j < m_box.GetBInteriorEnd(); j++) {

			for (int k = 0; k < nRElements; k++) {
				dDensityNode[k] = dataNode[RIx][i][j][k];

				// Convert U_alpha and V_beta to X and Y (Lon, Lat)
				dUXNode[k] = dataNode(UIx,i,j,k)
						- dDerivRNode(i,j,k,0)
						* (dataNode(WIx,i,j,k)
						/ dDerivRNode(i,j,k,2));

				dVYNode[k] = dataNode(VIx,i,j,k)
						- dDerivRNode(i,j,k,1)
						* (dataNode(WIx,i,j,k)
						/ dDerivRNode(i,j,k,2));
			}

			pGLLGrid->DifferentiateNodeToNode(
			dDensityNode,
			dDiffDensityNode);

			pGLLGrid->DifferentiateNodeToNode(
			dUXNode,
			dDiffUXNode);

			pGLLGrid->DifferentiateNodeToNode(
			dVYNode,
			dDiffVYNode);

			for (int k = 0; k < nRElements; k++) {
				m_dataRichardson[i][j][k] =
					- phys.GetG()
					/ (dDensityNode[k] * dDerivRNode(i,j,k,3))
					* dDiffDensityNode[k]
					/ ((dDiffUXNode[k] * dDiffUXNode[k]) +
					   (dDiffVYNode[k] * dDiffVYNode[k]));
				//
				if (m_dataRichardson[i][j][k] >= 1.0E6) {
					m_dataRichardson[i][j][k] = 1.0E6;
				} else if (m_dataRichardson[i][j][k] < 0.0) {
					m_dataRichardson[i][j][k] = 0.0;
				}
				//
			}
		}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void GridPatchGLL::ComputeConvectiveGrad(
	int iDataIndex,
	DataLocation loc
) {
	const PhysicalConstants & phys = m_grid.GetModel().GetPhysicalConstants();

	// Number of radial elements
	int nRElements = m_grid.GetRElements();

	// Get metric quantities
	const DataArray4D<double> & dDerivRNode =
		GetDerivRNode();

	// Indices of EquationSet variables
	const int UIx = 0;
	const int VIx = 1;
	const int PIx = 2;
	const int WIx = 3;
	const int RIx = 4;

	// Column array of theta and theta gradient
	DataArray1D<double> dThetaNode;
	DataArray1D<double> dDiffThetaNode;
	dThetaNode.Allocate(nRElements);
	dDiffThetaNode.Allocate(nRElements);
	DataArray1D<double> dTemperatureNode;
	dTemperatureNode.Allocate(nRElements);

	if (m_grid.GetModel().GetEquationSet().GetComponents() < 5) {
		_EXCEPTIONT("Invalid EquationSet.");
	}

	// Parent grid, containing the vertical remapping information
	GridGLL * pGLLGrid = dynamic_cast<GridGLL*>(&m_grid);
	if (pGLLGrid == NULL) {
		_EXCEPTIONT("Logic error");
	}

	// Calculate vertical derivative of potential temp. on nodes
	if (loc == DataLocation_Node) {
		if (m_grid.GetVarLocation(PIx) == DataLocation_REdge) {
			InterpolateREdgeToNode(PIx, iDataIndex);
		}

		const DataArray4D<double> & dataNode = m_datavecStateNode[iDataIndex];

		for (int i = m_box.GetAInteriorBegin(); i < m_box.GetAInteriorEnd(); i++) {
		for (int j = m_box.GetBInteriorBegin(); j < m_box.GetBInteriorEnd(); j++) {

			for (int k = 0; k < nRElements; k++) {
				dThetaNode[k] = 1.0;
#ifdef FORMULATION_PRESSURE
				const double dPressure = dataNode(PIx,i,j,k);
				dTemperatureNode[k] =
					dPressure / (dataNode(RIx,i,j,k) * phys.GetR());

				dThetaNode[k] = dTemperature[k] * pow(phys.GetP0()
					/ dPressure, phys.GetR / phys.GetCp);
#endif
#if defined(FORMULATION_RHOTHETA_PI) || defined(FORMULATION_RHOTHETA_P)
				const double dPressure = phys.PressureFromRhoTheta(dataNode(PIx,i,j,k));
				dTemperatureNode[k] = dPressure
						/ (dataNode(RIx,i,j,k) * phys.GetR());
				dThetaNode[k] = dataNode(PIx,i,j,k)
						/ dataNode(RIx,i,j,k);
#endif
#if defined(FORMULATION_THETA) || defined(FORMULATION_THETA_FLUX)
				const double dPressure =
					phys.PressureFromRhoTheta(
						dataNode(RIx,i,j,k)
						* dataNode(PIx,i,j,k));
				dTemperatureNode[k] =
					dPressure / (dataNode(RIx,i,j,k) * phys.GetR());

				dThetaNode[k] = dataNode(PIx,i,j,k);
#endif
			}

			pGLLGrid->DifferentiateNodeToNode(
			dThetaNode,
			dDiffThetaNode);

			for (int k = 0; k < nRElements; k++) {
				m_dataConvective[i][j][k] =
					(dTemperatureNode[k] / dThetaNode[k])
					* dDiffThetaNode[k]
					* dDerivRNode(i,j,k,3);
			}
		}
		}
	}
}
///////////////////////////////////////////////////////////////////////////////
