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

static const int ParamTotalPatchCount = 1;

///////////////////////////////////////////////////////////////////////////////

void InitializeModelAndGrid(
	Model * pModel
) {
	// Model parameters
	const int nResolutionX = 10;
	const int nResolutionY = 10;
	const int nHorizontalOrder = 4;
	const int nVerticalOrder = 4;
	const int nLevels = 40;
	const int nMaxPatchCount = ParamTotalPatchCount;

	double dGDim[6];
	dGDim[0] = 0.0;
	dGDim[1] = 1000.0;
	dGDim[2] = -1000.0;
	dGDim[3] = 1000.0;
	dGDim[4] = 0.0;
	dGDim[5] = 1000.0;

	const double dRefLat = 0.0;

	const Grid::VerticalStaggering eVerticalStaggering =
		Grid::VerticalStaggering_Levels;

	// Set the parameters
	ModelParameters param;

	pModel->SetParameters(param);

	pModel->SetTimestepScheme(new TimestepSchemeStrang(*pModel));

	pModel->SetHorizontalDynamics(
		new HorizontalDynamicsStub(*pModel));

	pModel->SetVerticalDynamics(
		new VerticalDynamicsStub(*pModel));

	// Set the model grid (one patch Cartesian grid)
	GridCartesianGLL * pGrid = 
		new GridCartesianGLL(
			*pModel,
			nMaxPatchCount,
			nResolutionX,
			nResolutionY,
			4,
			nHorizontalOrder,
			nVerticalOrder,
			nLevels,
			dGDim,
			dRefLat,
			eVerticalStaggering);

	pModel->SetGrid(pGrid);
}

///////////////////////////////////////////////////////////////////////////////

void InitializeTask(
	int ixPatch,
	unsigned char ** pGeometricData,
	unsigned char ** pActiveStateData,
	unsigned char ** pBufferStateData,
	unsigned char ** pAuxiliaryData
) {
	// Setup the Model
	Model model(EquationSet::PrimitiveNonhydrostaticEquations);

	InitializeModelAndGrid(&model);

	Grid * pGrid = model.GetGrid();

	// Set the data for the GridPatch
	GridPatch * pPatch = pGrid->GetPatch(ixPatch);

	_EXCEPTION();

	pPatch->InitializeDataLocal(false, false, false, false);

	*pGeometricData = (unsigned char *)
		malloc(pPatch->GetDataContainerGeometric().GetTotalByteSize());
	*pActiveStateData = (unsigned char *)
		malloc(pPatch->GetDataContainerActiveState().GetTotalByteSize());
	*pBufferStateData = (unsigned char *)
		malloc(pPatch->GetDataContainerBufferState().GetTotalByteSize());
	*pAuxiliaryData = (unsigned char *)
		malloc(pPatch->GetDataContainerAuxiliary().GetTotalByteSize());

	// TODO: Create a new ExecuteTask with Patch ixPatch, Step 0, SubStep 0
}

///////////////////////////////////////////////////////////////////////////////

void ExecuteTask(
	int ixPatch,
	int ixStep,
	int ixSubStep,
	unsigned char * pGeometricData,
	unsigned char * pActiveStateData,
	unsigned char * pBufferStateData,
	unsigned char * pAuxiliaryData
) {
	const int nHorizontalOrder = 4;
	const int nVerticalOrder = 4;

	// Setup the Model
	Model model(EquationSet::PrimitiveNonhydrostaticEquations);

	InitializeModelAndGrid(&model);

	Grid * pGrid = model.GetGrid();

	// Set the data for the GridPatch
	GridPatch * pPatch = pGrid->GetPatch(ixPatch);

	pPatch->InitializeDataLocal(false, false, false, false);

	pPatch->GetDataContainerGeometric().AttachTo(pGeometricData);
	pPatch->GetDataContainerActiveState().AttachTo(pActiveStateData);
	pPatch->GetDataContainerBufferState().AttachTo(pBufferStateData);
	pPatch->GetDataContainerAuxiliary().AttachTo(pAuxiliaryData);

	// Perform sub-step
	model.SubStep((ixStep == 0), false, ixSubStep);

	// Pack data from Patch into exchange buffers
	// UNIMPLEMENTED

	// TODO: Create a new ExecuteTask with Patch ixPatch, Step ixStep,
	// SubStep ixSubStep+1
}

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

	// Initialize MPI
	MPI_Init(&argc, &argv);

try {

	// TODO: These should be in DataBlocks
	unsigned char * pGeometricData[ParamTotalPatchCount];
	unsigned char * pActiveStateData[ParamTotalPatchCount];
	unsigned char * pBufferStateData[ParamTotalPatchCount];
	unsigned char * pAuxiliaryData[ParamTotalPatchCount];

	//////////////////////////////////////////////////////////////////
	// BEGIN INITIALIZATION BLOCK

	std::cout << "BEGIN Initialization" << std::endl;
	for (int ixPatch = 0; ixPatch < ParamTotalPatchCount; ixPatch++) {
		InitializeTask(0,
			&(pGeometricData[ixPatch]),
			&(pActiveStateData[ixPatch]),
			&(pBufferStateData[ixPatch]),
			&(pAuxiliaryData[ixPatch]));
	}
	std::cout << "DONE Initialization" << std::endl;

	// END INITIALIZATION BLOCK
	//////////////////////////////////////////////////////////////////

	//////////////////////////////////////////////////////////////////
	// BEGIN MAIN PROGRAM BLOCK

	std::cout << "BEGIN Execution" << std::endl;

	for (int ixPatch = 0; ixPatch < ParamTotalPatchCount; ixPatch++) {
	for (int ixStep = 0; ixStep < 1; ixStep++) {
	for (int ixSubStep = 0; ixSubStep < 1; ixSubStep++) {
		ExecuteTask(0,
			ixStep,
			ixSubStep,
			pGeometricData[ixPatch],
			pActiveStateData[ixPatch],
			pBufferStateData[ixPatch],
			pAuxiliaryData[ixPatch]);
	}
	}
	}
	std::cout << "DONE Execution" << std::endl;

	// END MAIN PROGRAM BLOCK
	//////////////////////////////////////////////////////////////////

} catch(Exception & e) {
	std::cout << e.ToString() << std::endl;
}

	// Deinitialize Tempest
	MPI_Finalize();
}

///////////////////////////////////////////////////////////////////////////////

