///////////////////////////////////////////////////////////////////////////////
///
///	\file    DataContainerTest.cpp
///	\author  Paul Ullrich, Jorge Guerra
///	\version December 18, 2013
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

#include "Tempest.h"
#include "GridPatchCartesianGLL.h"

#include <mpi.h>

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

	// Initialize MPI
	MPI_Init(&argc, &argv);

try {
	std::cout << "Initializing Model and Grid ... " << std::endl;

	// Model parameters
	const int nResolutionX = 10;
	const int nResolutionY = 10;
	const int nHorizontalOrder = 4;
	const int nVerticalOrder = 4;
	const int nLevels = 40;

	double dGDim[6];
	dGDim[0] = 0.0;
	dGDim[1] = 1000.0;
	dGDim[2] = -1000.0;
	dGDim[3] = 1000.0;
	dGDim[4] = 0.0;
	dGDim[5] = 1000.0;

	const double dRefLat = 0.0;
	const double dTopoHeight = 0.0;

	const Grid::VerticalStaggering eVerticalStaggering =
		Grid::VerticalStaggering_Levels;

	// Setup the Model
	Model model(EquationSet::PrimitiveNonhydrostaticEquations);

	// Set the parameters
	ModelParameters param;

	model.SetParameters(param);

	model.SetTimestepScheme(new TimestepSchemeStrang(model));

	model.SetHorizontalDynamics(
		new HorizontalDynamicsStub(model));

	model.SetVerticalDynamics(
		new VerticalDynamicsStub(model));

	// Set the model grid (one patch Cartesian grid)
	GridCartesianGLL * pGrid = 
		new GridCartesianGLL(
			model,
			nResolutionX,
			nResolutionY,
			4,
			nHorizontalOrder,
			nVerticalOrder,
			nLevels,
			dGDim,
			dRefLat,
			dTopoHeight,
			eVerticalStaggering);

	//////////////////////////////////////////////////////////////////
	// BEGIN MAIN PROGRAM BLOCK

	// Define a PatchBox
	double dDeltaA = (dGDim[1] - dGDim[0])
		/ static_cast<double>(nResolutionX);
	double dDeltaB = (dGDim[3] - dGDim[2])
		/ static_cast<double>(nResolutionY);

	GridSpacingGaussLobattoRepeated
		glspacingAlpha(dDeltaA, dGDim[0], nHorizontalOrder);
	GridSpacingGaussLobattoRepeated
		glspacingBeta(dDeltaB, dGDim[2], nHorizontalOrder);

	PatchBox box(
		0, 0, model.GetHaloElements(),
		0, nHorizontalOrder * nResolutionX,
		0, nHorizontalOrder * nResolutionY,
		glspacingAlpha,
		glspacingBeta);

	// Create a GridPatch
	GridPatch * pPatchFirst =
		new GridPatchCartesianGLL(
			(*pGrid),
			0,
			box,
			nHorizontalOrder,
			nVerticalOrder,
			dGDim,
			dRefLat,
			dTopoHeight);

	// Create a second GridPatch with the same box
	GridPatch * pPatchSecond =
		new GridPatchCartesianGLL(
			(*pGrid),
			0,
			box,
			nHorizontalOrder,
			nVerticalOrder,
			dGDim,
			dRefLat,
			dTopoHeight);

	// Build the DataContainer object for the GridPatch objects
	pPatchFirst->InitializeDataLocal(false, false, false, false);
	pPatchSecond->InitializeDataLocal(false, false, false, false);

	// Get DataContainers associated with GridPatch
	DataContainer & dataGeometric = pPatchFirst->GetDataContainerGeometric();
	DataContainer & dataActiveState = pPatchFirst->GetDataContainerActiveState();
	DataContainer & dataBufferState = pPatchFirst->GetDataContainerBufferState();
	DataContainer & dataAuxiliary = pPatchFirst->GetDataContainerAuxiliary();

	// Output the size requirements (in bytes) of each DataContainer
	std::cout << "GridPatch.Geometric Size:   " << dataGeometric.GetTotalByteSize() << " bytes" << std::endl;
	std::cout << "GridPatch.ActiveState Size: " << dataActiveState.GetTotalByteSize() << " bytes" << std::endl;
	std::cout << "GridPatch.BufferState Size: " << dataBufferState.GetTotalByteSize() << " bytes" << std::endl;
	std::cout << "GridPatch.Auxiliary Size:   " << dataAuxiliary.GetTotalByteSize() << " bytes" << std::endl;

	// Allocate data
	std::cout << "Allocating data ... " << std::endl;
	unsigned char * pDataGeometric = (unsigned char*)malloc(dataGeometric.GetTotalByteSize());
	unsigned char * pDataActiveState = (unsigned char*)malloc(dataActiveState.GetTotalByteSize());
	unsigned char * pDataBufferState = (unsigned char*)malloc(dataBufferState.GetTotalByteSize());
	unsigned char * pDataAuxiliary = (unsigned char*)malloc(dataAuxiliary.GetTotalByteSize());

	// Initialize data to zero
	memset(pDataGeometric, 0, dataGeometric.GetTotalByteSize());
	memset(pDataActiveState, 0, dataActiveState.GetTotalByteSize());
	memset(pDataBufferState, 0, dataBufferState.GetTotalByteSize());
	memset(pDataAuxiliary, 0, dataAuxiliary.GetTotalByteSize());

	// Attach data to DataContainers
	std::cout << "Attaching data to DataContainers ... " << std::endl;
	dataGeometric.AttachTo(pDataGeometric);
	dataActiveState.AttachTo(pDataActiveState);
	dataBufferState.AttachTo(pDataBufferState);
	dataAuxiliary.AttachTo(pDataAuxiliary);

	// Update data in GridPatch
	std::cout << "Updating GridPatch data ... " << std::endl;
	DataArray4D<double> & dataStateNode = pPatchFirst->GetDataState(0);
	for (int i = box.GetAInteriorBegin(); i < box.GetAInteriorEnd(); i++) {
	for (int j = box.GetBInteriorBegin(); j < box.GetAInteriorEnd(); j++) {
		dataStateNode[0][0][i][j] = dataStateNode[0][0][i][j] + 1.0;
	}
	}

	// Detach data from GridPatch
	std::cout << "Detaching data from pPatchFirst ... " << std::endl;
	dataGeometric.Detach();
	dataActiveState.Detach();
	dataBufferState.Detach();
	dataAuxiliary.Detach();

	// Attach data to second GridPatch
	std::cout << "Attaching data to pPatchSecond ... " << std::endl;
	pPatchSecond->GetDataContainerActiveState().AttachTo(pDataActiveState);

	// Output 
	std::cout << "Checking GridPatch data ... " << std::endl;
	DataArray4D<double> & dataStateNodeSecond = pPatchSecond->GetDataState(0);
	std::cout << "Test (should be 1.0): " << dataStateNodeSecond[0][0][1][1] << std::endl;

	if (dataStateNodeSecond[0][0][1][1] == 1.0) {
		std::cout << "Test Passed" << std::endl;
	} else {
		std::cout << "Test Failed" << std::endl;
	}

	// END MAIN PROGRAM BLOCK
	//////////////////////////////////////////////////////////////////

} catch(Exception & e) {
	std::cout << e.ToString() << std::endl;
}

	// Deinitialize Tempest
	MPI_Finalize();
}

///////////////////////////////////////////////////////////////////////////////

