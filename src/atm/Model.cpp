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

#include <cfloat>

///////////////////////////////////////////////////////////////////////////////

Model::Model(
	EquationSet::Type eEquationSetType
) :
	m_pParam(NULL),
	m_pGrid(NULL),
	m_pTimestepScheme(NULL),
	m_pHorizontalDynamics(NULL),
	m_pVerticalDynamics(NULL),
	m_pTestCase(NULL),
	m_eqn(eEquationSetType),
	m_dTime(0.0)
{
	// Initialize staggering from equation set
	m_stag.Initialize(m_eqn);
}

///////////////////////////////////////////////////////////////////////////////

void Model::SetParameters(ModelParameters * pParam) {
	m_pParam = pParam;
}

///////////////////////////////////////////////////////////////////////////////

void Model::SetGrid(Grid * pGrid) {
	if (pGrid == NULL) {
		_EXCEPTIONT("Invalid Grid (NULL)");
	}
	if (m_pGrid != NULL) {
		_EXCEPTIONT("Grid already specified");
	}

	// Attach the grid
	m_pGrid = pGrid;

	// Initialize the grid
	m_pGrid->Initialize();
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

void Model::SetTestCase(TestCase * pTestCase) {
	if (m_pGrid == NULL) {
		_EXCEPTIONT(
			"A grid must be specified before attaching a TestCase.");
	}
	if (pTestCase == NULL) {
		_EXCEPTIONT("Invalid TestCase (NULL)");
	}

	// Attach the test case
	m_pTestCase = pTestCase;

	// Evaluate physical constants
	m_pTestCase->EvaluatePhysicalConstants(m_phys);

	// Initialize the topography and data
	m_pGrid->EvaluateTestCase(*pTestCase);
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

void Model::Go() {

	// Check pointers
	if (m_pParam == NULL) {
		_EXCEPTIONT("ModelParameters not specified.");
	}
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
	if (m_pTestCase == NULL) {
		_EXCEPTIONT("TestCase not specified.");
	}
	if (m_vecOutMan.size() == 0) {
		printf("WARNING: No OutputManager specified.");
	}

	// Check time step
	if (m_pParam->m_dDeltaT == 0.0) {
		_EXCEPTIONT("Dynamic timestepping not implemented. "
		            "DeltaT must be non-zero.");
	}

	// Initialize all components
	m_pTimestepScheme->Initialize();
	m_pHorizontalDynamics->Initialize();
	m_pVerticalDynamics->Initialize();

	// Current time
	m_dTime = m_pParam->m_dBeginTime;

	// Check time
	if (m_dTime / m_pParam->m_dEndTime >= 1.0 - DBL_EPSILON) {
		return;
	}

	// Initial output
	for (int om = 0; om < m_vecOutMan.size(); om++) {
		m_vecOutMan[om]->InitialOutput(m_dTime);
	}

	// First time step
	bool fFirstStep = true;

	// Loop
	for(;;) {

		// Last time step
		bool fLastStep = false;

		// Next time step
		double dTimeNext = m_dTime + m_pParam->m_dDeltaT;

		// Time step size
		double dDeltaT = m_pParam->m_dDeltaT;

		// Perform a semi-timestep if necessary to align timescales
		if (dTimeNext / m_pParam->m_dEndTime >= 1.0 - DBL_EPSILON) {
			dDeltaT = m_pParam->m_dEndTime - m_dTime;
			fLastStep = true;
		}

		// Perform one time step
		m_pTimestepScheme->Step(fFirstStep, fLastStep, m_dTime, dDeltaT);

		// Update the timer
		m_dTime += dDeltaT;

		if ((fLastStep) ||
			(m_dTime / m_pParam->m_dEndTime >= (1.0 - DBL_EPSILON))
		) {
			break;
		}

		// Check for output
		for (int om = 0; om < m_vecOutMan.size(); om++) {
			if (m_vecOutMan[om]->IsOutputNeeded(m_dTime)) {
				m_vecOutMan[om]->ManageOutput(m_dTime);
			}
		}

		// First time step
		fFirstStep = false;
	}

	// Final output
	for (int om = 0; om < m_vecOutMan.size(); om++) {
		m_vecOutMan[om]->FinalOutput(m_dTime);
	}
}

///////////////////////////////////////////////////////////////////////////////

void Model::ComputeErrorNorms() {
	if (m_pGrid == NULL) {
		_EXCEPTIONT("Model Grid not specified.");
	}
	if (m_pTestCase == NULL) {
		_EXCEPTIONT("TestCase not specified.");
	}

	// Local rank
	int nRank;
	MPI_Comm_rank(MPI_COMM_WORLD, &nRank);

	// Initialize test case reference data
	m_pGrid->EvaluateTestCase(*m_pTestCase, m_dTime, 1);

	// Construct the reference state
	m_pGrid->CopyData(1, 2, DataType_State);
	m_pGrid->AddReferenceState(2);

	// Remove current state
	DataVector<double> dDifference;
	dDifference.Initialize(2);
	dDifference[0] = -1.0;
	dDifference[1] = +1.0;

	m_pGrid->LinearCombineData(dDifference, 1, DataType_State);

	// Compute error norms
	DataMatrix<double> dErrorSums;
	dErrorSums.Initialize(m_eqn.GetComponents(), 3);

	DataMatrix<double> dErrorNorms;
	dErrorNorms.Initialize(m_eqn.GetComponents(), 3);

	DataVector<double> dSums;
	DataVector<double> dNorms;

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

