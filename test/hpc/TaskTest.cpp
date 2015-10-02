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

static const int ParamHorizontalOrder = 4;
static const int ParamVerticalOrder = 4;

static const int ParamTotalPatchCount = 1;

static const int ParamTotalSteps = 1;

///////////////////////////////////////////////////////////////////////////////

void InitializeModelAndGrid(
	Model * pModel
) {
	// Model parameters
	const int nResolutionX = 10;
	const int nResolutionY = 10;
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

	param.m_timeDeltaT.FromFormattedString("1s");

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
			ParamHorizontalOrder,
			ParamVerticalOrder,
			nLevels,
			dGDim,
			dRefLat,
			eVerticalStaggering);

	pModel->SetGrid(pGrid, false);
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
	printf("CURR:  GridPatch %i performing step %i.%i\n",
		ixPatch, ixStep, ixSubStep);
/*
	model.SubStep(
		(ixStep == 0),
		(ixStep == ParamTotalSteps-1),
		ixSubStep);
*/
	// TODO: Pack data from Patch into exchange buffers and tag as complete

	// TODO: Create a new ExecuteTask on Patch ixPatch, ixNextStep, ixNextSubStep
	int nSubStepsPerStep = model.GetTimestepScheme()->GetSubStepCount();

	int ixNextSubStep;
	int ixNextStep;

	if (ixSubStep == nSubStepsPerStep-1) {
		ixNextSubStep = 0;
		ixNextStep = ixStep + 1;
	} else {
		ixNextSubStep = ixSubStep + 1;
		ixNextStep = ixStep;
	}

	printf("NEXT:  GridPatch %i to perform step %i.%i (new EDT)\n",
		ixPatch, ixNextStep, ixNextSubStep);
}

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

	// Initialize MPI
	MPI_Init(&argc, &argv);

try {

	// TODO: Create DataBlocks for all GridPatch data
	unsigned char * pGeometricData[ParamTotalPatchCount];
	unsigned char * pActiveStateData[ParamTotalPatchCount];
	unsigned char * pBufferStateData[ParamTotalPatchCount];
	unsigned char * pAuxiliaryData[ParamTotalPatchCount];

	//////////////////////////////////////////////////////////////////
	// BEGIN INITIALIZATION BLOCK

	std::cout << "BEGIN Initialization" << std::endl;

	for (int ixPatch = 0; ixPatch < ParamTotalPatchCount; ixPatch++) {
		InitializeTask(
			ixPatch,
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
	for (int ixStep = 0; ixStep < ParamTotalSteps; ixStep++) {
	for (int ixSubStep = 0; ixSubStep < 1; ixSubStep++) {
		ExecuteTask(
			ixPatch,
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

