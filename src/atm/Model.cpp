///////////////////////////////////////////////////////////////////////////////
///
///	\file    Model.h
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

#include "Model.h"

#include "Grid.h"
#include "TestCase.h"
#include "OutputManager.h"
#include "OutputManagerComposite.h"

#include "FunctionTimer.h"
#include "Announce.h"
#include "MemoryTools.h"

#include <cfloat>

///////////////////////////////////////////////////////////////////////////////

Model::Model(
	EquationSet::Type eEquationSetType
) :
	m_fGridFromRestartFile(false),
	m_pGrid(NULL),
	m_pTimestepScheme(NULL),
	m_pHorizontalDynamics(NULL),
	m_pVerticalDynamics(NULL),
	m_pTestCase(NULL),
	m_eqn(eEquationSetType),
	m_time()
{
}

///////////////////////////////////////////////////////////////////////////////

Model::Model(
	const EquationSet & eqn
) :
	m_fGridFromRestartFile(false),
	m_pGrid(NULL),
	m_pTimestepScheme(NULL),
	m_pHorizontalDynamics(NULL),
	m_pVerticalDynamics(NULL),
	m_pTestCase(NULL),
	m_eqn(eqn),
	m_metaUserData(),
	m_time()
{
}

///////////////////////////////////////////////////////////////////////////////

Model::Model(
	const EquationSet & eqn,
	const UserDataMeta & metaUserData
) :
	m_fGridFromRestartFile(false),
	m_pGrid(NULL),
	m_pTimestepScheme(NULL),
	m_pHorizontalDynamics(NULL),
	m_pVerticalDynamics(NULL),
	m_pTestCase(NULL),
	m_eqn(eqn),
	m_metaUserData(metaUserData),
	m_time()
{
}

///////////////////////////////////////////////////////////////////////////////

Model::~Model() {
	if (m_pGrid != NULL) {
		delete m_pGrid;
	}
	if (m_pTimestepScheme != NULL) {
		delete m_pTimestepScheme;
	}
	if (m_pHorizontalDynamics != NULL) {
		delete m_pHorizontalDynamics;
	}
	if (m_pVerticalDynamics != NULL) {
		delete m_pVerticalDynamics;
	}

	for (int n = 0; n < m_vecOutMan.size(); n++) {
		delete m_vecOutMan[n];
	}

	if (m_pTestCase != NULL) {
		delete m_pTestCase;
	}
}

///////////////////////////////////////////////////////////////////////////////

void Model::SetGrid(
	Grid * pGrid,
	int nPatchCount,
	bool fInitializeConnectivity
) {
	if (pGrid == NULL) {
		_EXCEPTIONT("Invalid Grid (NULL)");
	}
	if (m_pGrid != NULL) {
		_EXCEPTIONT("Grid already specified");
	}

	// Attach the grid
	m_pGrid = pGrid;

	// Set up patches
#ifdef TEMPEST_MPIOMP
	if (nPatchCount == (-1)) {
		MPI_Comm_size(MPI_COMM_WORLD, &nPatchCount);
	}
#else
	if (nPatchCount == (-1)) {
		_EXCEPTIONT("Unimplemented: PatchCount must be specified");
	}
#endif

	m_pGrid->ApplyDefaultPatchLayout(nPatchCount);

	// Initialize the grid
	m_pGrid->Initialize();

	if (fInitializeConnectivity) {
		m_pGrid->DistributePatches();
		m_pGrid->InitializeExchangeBuffersFromActivePatches();
		m_pGrid->InitializeConnectivity();
	}
}

///////////////////////////////////////////////////////////////////////////////

void Model::SetGridFromRestartFile(
	Grid * pGrid,
	const std::string & strRestartFile
) {
	if (pGrid == NULL) {
		_EXCEPTIONT("Invalid Grid (NULL)");
	}
	if (m_pGrid != NULL) {
		_EXCEPTIONT("Grid already specified");
	}

	// Attach the grid
	m_pGrid = pGrid;

	// Load the Grid data from the file
	OutputManagerComposite ompComposite(*m_pGrid, Time(), "", "", "");

	m_time = ompComposite.Input(strRestartFile);

	m_timeStart = m_time;

	// Initialize the grid
	m_pGrid->Initialize();

	// Initialize exchange buffers and connectivity
	m_pGrid->InitializeExchangeBuffersFromActivePatches();
	m_pGrid->InitializeConnectivity();

	// Set flag indicating Grid was initialized from restart file
	m_fGridFromRestartFile = true;
}

///////////////////////////////////////////////////////////////////////////////

void Model::SetTimestepScheme(TimestepScheme * pTimestepScheme) {
	if (pTimestepScheme == NULL) {
		_EXCEPTIONT("Invalid TimestepScheme (NULL)");
	}
	if (m_pTimestepScheme != NULL) {
		_EXCEPTIONT("TimestepScheme already specified");
	}

	m_pTimestepScheme = pTimestepScheme;
}

///////////////////////////////////////////////////////////////////////////////

void Model::SetHorizontalDynamics(HorizontalDynamics * pHorizontalDynamics) {
	if (pHorizontalDynamics == NULL) {
		_EXCEPTIONT("Invalid HorizontalDynamics (NULL)");
	}
	if (m_pHorizontalDynamics != NULL) {
		_EXCEPTIONT("HorizontalDynamics already specified");
	}

	m_pHorizontalDynamics = pHorizontalDynamics;
}

///////////////////////////////////////////////////////////////////////////////

void Model::SetVerticalDynamics(VerticalDynamics * pVerticalDynamics) {
	if (pVerticalDynamics == NULL) {
		_EXCEPTIONT("Invalid VerticalDynamics (NULL)");
	}
	if (m_pVerticalDynamics != NULL) {
		_EXCEPTIONT("VerticalDynamics already specified");
	}

	m_pVerticalDynamics = pVerticalDynamics;
}

///////////////////////////////////////////////////////////////////////////////

void Model::SetTestCase(
	TestCase * pTestCase
) {
	if (m_pGrid == NULL) {
		_EXCEPTIONT(
			"A grid must be specified before attaching a TestCase.");
	}
	if (pTestCase == NULL) {
		_EXCEPTIONT("Invalid TestCase (NULL)");
	}

	// Attach the test case
	m_pTestCase = pTestCase;

	// Evaluate physical constants and data from TestCase
	if (!m_fGridFromRestartFile) {

		// Evaluate physical constants
		m_pTestCase->EvaluatePhysicalConstants(m_phys);

		// Initialize the topography and data
		m_pGrid->EvaluateTestCase(*pTestCase, m_timeStart);
	}
}

///////////////////////////////////////////////////////////////////////////////

void Model::AttachOutputManager(OutputManager * pOutMan) {
	if (m_pGrid == NULL) {
		_EXCEPTIONT(
			"A grid must be specified before attaching an OutputManager.");
	}

	// Attach output manager
	m_vecOutMan.push_back(pOutMan);
}

///////////////////////////////////////////////////////////////////////////////

void Model::AttachWorkflowProcess(WorkflowProcess * pWorkflowProcess) {
	if (m_pGrid == NULL) {
		_EXCEPTIONT(
			"A grid must be specified before attaching an WorkflowProcess.");
	}

	// Attach output manager
	m_vecWorkflowProcess.push_back(pWorkflowProcess);
}

///////////////////////////////////////////////////////////////////////////////

void Model::SubStep(
	bool fFirstStep,
	bool fLastStep,
	int iSubStep
) {
	// Next time step
	double dDeltaT;

	Time timeNext = m_time;
	timeNext += m_timeDeltaT;

	if (timeNext >= m_timeEnd) {
		dDeltaT = m_timeEnd - m_time;
		fLastStep = true;

	} else {
		dDeltaT = timeNext - m_time;
	}

	// Substep
	m_pTimestepScheme->SubStep(
		fFirstStep,
		fLastStep,
		m_time,
		dDeltaT,
		iSubStep);
}

///////////////////////////////////////////////////////////////////////////////

void Model::Go() {

	// Check pointers
	if (m_pGrid == NULL) {
		_EXCEPTIONT("Grid not specified.");
	}
	if (m_pTimestepScheme == NULL) {
		_EXCEPTIONT("TimestepScheme not specified.");
	}
	if (m_pHorizontalDynamics == NULL) {
		_EXCEPTIONT("HorizontalDynamics not specified.");
	}
	if (m_pVerticalDynamics == NULL) {
		_EXCEPTIONT("VerticalDynamics not specified.");
	}
	if (m_vecOutMan.size() == 0) {
		Announce("WARNING: No OutputManager specified.");
	}
	if (m_pTestCase == NULL) {
		_EXCEPTIONT("TestCase not specified.");
	}

	// Check time step
	if (m_timeDeltaT.IsZero()) {
		_EXCEPTIONT("Dynamic timestepping not implemented. "
		            "DeltaT must be non-zero.");
	}

	// Evaluate geometric terms in the grid
	// NOTE: This needs to be called after EvaluateTestCase, since it relies
	// on information about topographic derivatives.
	m_pGrid->EvaluateGeometricTerms();

	// Apply boundary conditions
	m_pGrid->ApplyBoundaryConditions();

	// Initialize all components
	m_pTimestepScheme->Initialize();
	m_pHorizontalDynamics->Initialize();
	m_pVerticalDynamics->Initialize();

	// Set the current time
	m_time = m_timeStart;

	// Check time
	if (m_time >= m_timeEnd) {
		Announce("Warning: Simulation start time (%s)\n"
			"  equals or exceeds end time (%s)",
			m_time.ToString().c_str(),
			m_timeEnd.ToString().c_str());
		return;
	}

	// Initial output
	for (int om = 0; om < m_vecOutMan.size(); om++) {
/* COMMENT IN FOR MASS, ENERGY, AND MOMENTUM OUTPUTS
		if (om == 0) {
			Announce("%s %1.15e %1.15e %1.15e",
			"Energy:",
			m_pGrid->ComputeTotalEnergy(0),
			m_pGrid->ComputeTotalPotentialEnstrophy(0),
			m_pGrid->ComputeTotalVerticalMomentum(0));
		}
*/
		m_vecOutMan[om]->InitialOutput(m_time);
	}

	// Initialize WorkflowProcesses
	for (int wfp = 0; wfp < m_vecWorkflowProcess.size(); wfp++) {
		m_vecWorkflowProcess[wfp]->Initialize(m_time);
	}

	// First time step
	bool fFirstStep = true;

	// Reset the communication timer
	FunctionTimer::ResetGroupTimeRecord("Communicate");

	// Loop
	for(int iStep = 0;; iStep++) {

		//PrintMemoryLine();

		FunctionTimer timerLoop("Loop");

		// Last time step
		bool fLastStep = false;

		// Next time step
		double dDeltaT;

		Time timeNext = m_time;
		timeNext += m_timeDeltaT;

		if (timeNext >= m_timeEnd) {
			dDeltaT = m_timeEnd - m_time;
			fLastStep = true;

		} else {
			dDeltaT = timeNext - m_time;
		}

		// Perform one time step
		Announce("Step %s", m_time.ToString().c_str());
		m_pTimestepScheme->Step(fFirstStep, fLastStep, m_time, dDeltaT);
/*
		// Energy and enstrophy
		{
			if (m_eqn.GetDimensionality() == 3) {
				if (m_pGrid->GetVerticalStaggering() ==
				    Grid::VerticalStaggering_Lorenz
				) {
					m_pGrid->InterpolateREdgeToNode(3, 0);
				}

				if (m_pGrid->GetVerticalStaggering() ==
				    Grid::VerticalStaggering_CharneyPhillips
				) {
					m_pGrid->InterpolateREdgeToNode(2, 0);
					m_pGrid->InterpolateREdgeToNode(3, 0);
				}
			}

			Announce("%1.15e %1.15e",
				m_pGrid->ComputeTotalEnergy(0),
				m_pGrid->ComputeTotalPotentialEnstrophy(0));

		}
*/
/*
		// L2 errors of the height field
		{
			DataArray1D<double> dSums;
			DataArray1D<double> dNorms;

			m_pGrid->EvaluateTestCase_StateOnly(*m_pTestCase, m_time, 2);

			m_pGrid->CopyData(2, 3, DataType_State);

			DataArray1D<double> dDifference;
			dDifference.Initialize(3);
			dDifference[0] = -1.0;
			dDifference[1] =  0.0;
			dDifference[2] = +1.0;

			m_pGrid->LinearCombineData(dDifference, 2, DataType_State);

			m_pGrid->Checksum(DataType_State, dSums, 3, ChecksumType_L2);
			m_pGrid->Checksum(DataType_State, dNorms, 2, ChecksumType_L2);

			Announce("%1.10e", dNorms[2] / dSums[2]);
		}
*/
		// Update the timer
		if (timeNext >= m_timeEnd) {
			m_time = m_timeEnd;
		} else {
			m_time = timeNext;
		}

		// Check for WorkflowProcesses
		for (int wfp = 0; wfp < m_vecWorkflowProcess.size(); wfp++) {
			if (m_vecWorkflowProcess[wfp]->IsReady(m_time)) {
				m_vecWorkflowProcess[wfp]->Perform(m_time);
			}
		}

		// Check for output
		for (int om = 0; om < m_vecOutMan.size(); om++) {
			if (fLastStep) {
/* COMMENT IN FOR MASS, ENERGY, AND MOMENTUM OUTPUTS
				if (om == 0) {
						Announce("%s %1.15e %1.15e %1.15e",
						"Energy:",
						m_pGrid->ComputeTotalEnergy(0),
						m_pGrid->ComputeTotalPotentialEnstrophy(0),
						m_pGrid->ComputeTotalVerticalMomentum(0));
				}
*/
				m_vecOutMan[om]->FinalOutput(m_time);

			} else if (m_vecOutMan[om]->IsOutputNeeded(m_time)) {
/* COMMENT IN FOR MASS, ENERGY, AND MOMENTUM OUTPUTS
				if (om == 0) {
						Announce("%s %1.15e %1.15e %1.15e",
						"Energy:",
						m_pGrid->ComputeTotalEnergy(0),
						m_pGrid->ComputeTotalPotentialEnstrophy(0),
						m_pGrid->ComputeTotalVerticalMomentum(0));
				}
*/
				m_vecOutMan[om]->ManageOutput(m_time);
			}
		}

		// Exit on last step
		if (fLastStep) {
			break;
		}

		// No longer first time step
		fFirstStep = false;
	}

#if defined(TEMPEST_MPIOMP)
	{
		int nCommSize;
		MPI_Comm_size(MPI_COMM_WORLD, &nCommSize);

		long lGlobalTimeLoop[3];
		long lTimeLoop =
			FunctionTimer::GetAverageGroupTime("Loop");

#if defined(PROGNOSTIC_CONTRAVARIANT_MOMENTA)
		long lGlobalTimeTend[3];
		long lTimeTend =
			FunctionTimer::GetAverageGroupTime(
				"CalculateTendencies");

		long lGlobalTimeAcLo[3];
		long lTimeAcLo =
			FunctionTimer::GetAverageGroupTime(
				"AcousticLoop");
#else
		long lGlobalTimeSNHP[3];
		long lTimeHSNP =
			FunctionTimer::GetAverageGroupTime(
				"HorizontalStepNonhydrostaticPrimitive");

		long lGlobalTimeVSEx[3];
		long lTimeVSEx =
			FunctionTimer::GetAverageGroupTime(
				"VerticalStepExplicit");

		long lGlobalTimeVSIm[3];
		long lTimeVSIm =
			FunctionTimer::GetAverageGroupTime(
				"VerticalStepImplicit");
#endif
		long lGlobalTimeSaSc[3];
		long lTimeSaSc =
			FunctionTimer::GetAverageGroupTime(
				"StepAfterSubCycle");

		long lGlobalTimeComm[3];
		long lTimeComm =
			FunctionTimer::GetAverageGroupTime(
				"Communicate");

		MPI_Reduce(&lTimeLoop, &lGlobalTimeLoop[0],
			1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(&lTimeLoop, &lGlobalTimeLoop[1],
			1, MPI_LONG, MPI_MIN, 0, MPI_COMM_WORLD);
		MPI_Reduce(&lTimeLoop, &lGlobalTimeLoop[2],
			1, MPI_LONG, MPI_MAX, 0, MPI_COMM_WORLD);

		lGlobalTimeLoop[0] /= static_cast<long>(nCommSize);

#if defined(PROGNOSTIC_CONTRAVARIANT_MOMENTA)
		MPI_Reduce(&lTimeTend, &lGlobalTimeTend[0],
			1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(&lTimeTend, &lGlobalTimeTend[1],
			1, MPI_LONG, MPI_MIN, 0, MPI_COMM_WORLD);
		MPI_Reduce(&lTimeTend, &lGlobalTimeTend[2],
			1, MPI_LONG, MPI_MAX, 0, MPI_COMM_WORLD);

		lGlobalTimeTend[0] /= static_cast<long>(nCommSize);

		MPI_Reduce(&lTimeAcLo, &lGlobalTimeAcLo[0],
			1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(&lTimeAcLo, &lGlobalTimeAcLo[1],
			1, MPI_LONG, MPI_MIN, 0, MPI_COMM_WORLD);
		MPI_Reduce(&lTimeAcLo, &lGlobalTimeAcLo[2],
			1, MPI_LONG, MPI_MAX, 0, MPI_COMM_WORLD);

		lGlobalTimeAcLo[0] /= static_cast<long>(nCommSize);

#else
		MPI_Reduce(&lTimeHSNP, &lGlobalTimeSNHP[0],
			1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(&lTimeHSNP, &lGlobalTimeSNHP[1],
			1, MPI_LONG, MPI_MIN, 0, MPI_COMM_WORLD);
		MPI_Reduce(&lTimeHSNP, &lGlobalTimeSNHP[2],
			1, MPI_LONG, MPI_MAX, 0, MPI_COMM_WORLD);

		lGlobalTimeSNHP[0] /= static_cast<long>(nCommSize);

		MPI_Reduce(&lTimeVSEx, &lGlobalTimeVSEx[0],
			1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(&lTimeVSEx, &lGlobalTimeVSEx[1],
			1, MPI_LONG, MPI_MIN, 0, MPI_COMM_WORLD);
		MPI_Reduce(&lTimeVSEx, &lGlobalTimeVSEx[2],
			1, MPI_LONG, MPI_MAX, 0, MPI_COMM_WORLD);

		lGlobalTimeVSEx[0] /= static_cast<long>(nCommSize);

		MPI_Reduce(&lTimeVSIm, &lGlobalTimeVSIm[0],
			1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(&lTimeVSIm, &lGlobalTimeVSIm[1],
			1, MPI_LONG, MPI_MIN, 0, MPI_COMM_WORLD);
		MPI_Reduce(&lTimeVSIm, &lGlobalTimeVSIm[2],
			1, MPI_LONG, MPI_MAX, 0, MPI_COMM_WORLD);

		lGlobalTimeVSIm[0] /= static_cast<long>(nCommSize);
#endif

		MPI_Reduce(&lTimeSaSc, &lGlobalTimeSaSc[0],
			1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(&lTimeSaSc, &lGlobalTimeSaSc[1],
			1, MPI_LONG, MPI_MIN, 0, MPI_COMM_WORLD);
		MPI_Reduce(&lTimeSaSc, &lGlobalTimeSaSc[2],
			1, MPI_LONG, MPI_MAX, 0, MPI_COMM_WORLD);

		lGlobalTimeSaSc[0] /= static_cast<long>(nCommSize);

		MPI_Reduce(&lTimeComm, &lGlobalTimeComm[0],
			1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(&lTimeComm, &lGlobalTimeComm[1],
			1, MPI_LONG, MPI_MIN, 0, MPI_COMM_WORLD);
		MPI_Reduce(&lTimeComm, &lGlobalTimeComm[2],
			1, MPI_LONG, MPI_MAX, 0, MPI_COMM_WORLD);

		lGlobalTimeComm[0] /= static_cast<long>(nCommSize);

		int nEntriesLoop =
			FunctionTimer::GetNumberOfEntries("Loop");
		Announce("Time [Loop]: %li [%li, %li] (%i)",
			lGlobalTimeLoop[0], lGlobalTimeLoop[1], lGlobalTimeLoop[2], nEntriesLoop);

#if defined(PROGNOSTIC_CONTRAVARIANT_MOMENTA)
		int nEntriesTend =
			FunctionTimer::GetNumberOfEntries(
				"CalculateTendencies");
		Announce("Time [Tend]: %li [%li, %li] (%i)",
			lGlobalTimeTend[0], lGlobalTimeTend[1], lGlobalTimeTend[2], nEntriesTend);

		int nEntriesAcLo =
			FunctionTimer::GetNumberOfEntries(
				"AcousticLoop");
		Announce("Time [AcLo]: %li [%li, %li] (%i)",
			lGlobalTimeAcLo[0], lGlobalTimeAcLo[1], lGlobalTimeAcLo[2], nEntriesAcLo);

#else
		int nEntriesSNHP =
			FunctionTimer::GetNumberOfEntries(
				"HorizontalStepNonhydrostaticPrimitive");
		Announce("Time [SNHP]: %li [%li, %li] (%i)",
			lGlobalTimeSNHP[0], lGlobalTimeSNHP[1], lGlobalTimeSNHP[2], nEntriesSNHP);

		int nEntriesVSEx =
			FunctionTimer::GetNumberOfEntries(
				"VerticalStepExplicit");
		Announce("Time [VSEx]: %li [%li, %li] (%i)",
			lGlobalTimeVSEx[0], lGlobalTimeVSEx[1], lGlobalTimeVSEx[2], nEntriesVSEx);

		int nEntriesVSIm =
			FunctionTimer::GetNumberOfEntries(
				"VerticalStepImplicit");
		Announce("Time [VSIm]: %li [%li, %li] (%i)",
			lGlobalTimeVSIm[0], lGlobalTimeVSIm[1], lGlobalTimeVSIm[2], nEntriesVSIm);
#endif

		int nEntriesSaSc =
			FunctionTimer::GetNumberOfEntries(
				"StepAfterSubCycle");
		Announce("Time [SaSc]: %li [%li, %li] (%i)",
			lGlobalTimeSaSc[0], lGlobalTimeSaSc[1], lGlobalTimeSaSc[2], nEntriesSaSc);

		int nEntriesComm =
			FunctionTimer::GetNumberOfEntries(
				"Communicate");
		Announce("Time [Comm]: %li [%li, %li] (%i)",
			lGlobalTimeComm[0], lGlobalTimeComm[1], lGlobalTimeComm[2], nEntriesComm);
	}
#endif
}

///////////////////////////////////////////////////////////////////////////////

void Model::ComputeErrorNorms() {
	if (m_pTestCase == NULL) {
		Announce("Error: No TestCase specified; cannot compute error norms.");
		return;
	}
	if (m_pGrid == NULL) {
		_EXCEPTIONT("Model Grid not specified.");
	}

	// Local rank
	int nRank = 0;
#if defined(TEMPEST_MPIOMP)
	MPI_Comm_rank(MPI_COMM_WORLD, &nRank);
#endif

	// Initialize test case reference data
	m_pGrid->EvaluateTestCase(*m_pTestCase, m_time, 1);

	// Construct the reference state
	m_pGrid->CopyData(1, 2, DataType_State);
	m_pGrid->AddReferenceState(2);

	// Remove current state
	DataArray1D<double> dDifference;
	dDifference.Allocate(2);
	dDifference[0] = -1.0;
	dDifference[1] = +1.0;

	m_pGrid->LinearCombineData(dDifference, 1, DataType_State);

	// Compute error norms
	DataArray2D<double> dErrorSums(m_eqn.GetComponents(), 3);
	DataArray2D<double> dErrorNorms(m_eqn.GetComponents(), 3);

	DataArray1D<double> dSums;
	DataArray1D<double> dNorms;

	// L1 errors
	m_pGrid->Checksum(DataType_State, dSums, 2, ChecksumType_L1);
	m_pGrid->Checksum(DataType_State, dNorms, 1, ChecksumType_L1);

	if (nRank == 0) {
		for (int c = 0; c < m_eqn.GetComponents(); c++) {
			dErrorSums[c][0] = dSums[c];
			dErrorNorms[c][0] = dNorms[c];
		}
	}

	// L2 errors
	m_pGrid->Checksum(DataType_State, dSums, 2, ChecksumType_L2);
	m_pGrid->Checksum(DataType_State, dNorms, 1, ChecksumType_L2);

	if (nRank == 0) {
		for (int c = 0; c < m_eqn.GetComponents(); c++) {
			dErrorSums[c][1] = dSums[c];
			dErrorNorms[c][1] = dNorms[c];
		}
	}

	// Linf errors
	m_pGrid->Checksum(DataType_State, dSums, 2, ChecksumType_Linf);
	m_pGrid->Checksum(DataType_State, dNorms, 1, ChecksumType_Linf);

	if (nRank == 0) {
		for (int c = 0; c < m_eqn.GetComponents(); c++) {
			dErrorSums[c][2] = dSums[c];
			dErrorNorms[c][2] = dNorms[c];
		}
	}

	// Output
	if (nRank == 0) {
		printf("Variable  L1 Error     L2 Error     Linf Error\n");
		printf("--------  -----------  -----------  -----------\n");
		for (int c = 0; c < m_eqn.GetComponents(); c++) {
			double dL1Error = (dErrorSums[c][0] != 0.0)?
				(dErrorNorms[c][0] / dErrorSums[c][0]):(dErrorNorms[c][0]);
			double dL2Error = (dErrorSums[c][1] != 0.0)?
				(dErrorNorms[c][1] / dErrorSums[c][1]):(dErrorNorms[c][1]);
			double dLinfError = (dErrorSums[c][2] != 0.0)?
				(dErrorNorms[c][2] / dErrorSums[c][2]):(dErrorNorms[c][2]);

			printf("%6s    %1.5e  %1.5e  %1.5e\n",
				m_eqn.GetComponentShortName(c).c_str(),
				dL1Error, dL2Error, dLinfError);
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

