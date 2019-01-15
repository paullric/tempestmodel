///////////////////////////////////////////////////////////////////////////////
///
///	\file    TerminatorPhysics.cpp
///	\author  Paul Ullrich
///	\version June 1, 2016
///
///	<remarks>
///		Copyright 2000-2016 Paul Ullrich
///
///		This file is distributed as part of the Tempest source code package.
///		Permission is granted to use, copy, modify and distribute this
///		source code and its documentation under the terms of the GNU General
///		Public License.  This software is provided "as is" without express
///		or implied warranty.
///	</remarks>

#include "TerminatorPhysics.h"

#include "Model.h"
#include "GridGLL.h"

#include "Announce.h" 

///////////////////////////////////////////////////////////////////////////////

extern "C" {
	void tendency_Terminator(
		double * lat,
		double * lon,
		double * cl,
		double * cl2,
		double * dt,
		double * cl_f,
		double * cl2_f);
}

///////////////////////////////////////////////////////////////////////////////

TerminatorPhysics::TerminatorPhysics(
	Model & model,
	const Time & timeFrequency
) :
WorkflowProcess(
	model,
	timeFrequency)
{ }

//////////////////////////////////////////////////////////////////////////////

void TerminatorPhysics::Initialize(
	const Time & timeStart
) {
}

	
///////////////////////////////////////////////////////////////////////////////
	
void TerminatorPhysics::Perform(
	const Time & time
) {
	// Indices of EquationSet variables
	const int UIx = 0;
	const int VIx = 1;
	const int TIx = 2;
	const int WIx = 3;
	const int RIx = 4;

	// Get DeltaT
	double dDeltaT = m_timeFrequency.GetSeconds();

	// Get a copy of the GLL grid
	GridGLL * pGridGLL = dynamic_cast<GridGLL*>(m_model.GetGrid());
	if (pGridGLL == NULL) {
		_EXCEPTIONT("Not implemented for a general Grid");
	}

	// Number of radial elements
	int nRElements = pGridGLL->GetRElements();

	// Perform local update
	for (int n = 0; n < pGridGLL->GetActivePatchCount(); n++) {
		GridPatch * pPatch = pGridGLL->GetActivePatch(n);

		const PatchBox & box = pPatch->GetPatchBox();

		// Grid data
		DataArray4D<double> & dataNode =
			pPatch->GetDataState(0, DataLocation_Node);

		DataArray4D<double> & dataTracer =
			pPatch->GetDataTracers(0);

		const DataArray2D<double> & dataLatitude =
			pPatch->GetLatitude();

		const DataArray2D<double> & dataLongitude =
			pPatch->GetLongitude();

		if (dataTracer.GetSize(0) <= 4) {
			_EXCEPTIONT("Insufficient entries in tracer data");
		}


		// Loop over all nodes in GridPatch
		for (int k = 0; k < nRElements; k++) {
		for (int i = box.GetAInteriorBegin(); i < box.GetAInteriorEnd(); i++) {
		for (int j = box.GetBInteriorBegin(); j < box.GetBInteriorEnd(); j++) {

			double dRho = dataNode[RIx][k][i][j];

			double dLatDeg = dataLatitude[i][j] * 180.0 / M_PI;
			double dLonDeg = dataLongitude[i][j] * 180.0 / M_PI;
			
			double dQCl = dataTracer[3][k][i][j] / dRho;
			double dQCl2 = dataTracer[4][k][i][j] / dRho;

			double dTendencyQCl = 0.0;
			double dTendencyQCl2 = 0.0;

			tendency_Terminator(
				&dLatDeg,
				&dLonDeg,
				&dQCl,
				&dQCl2,
				&dDeltaT,
				&dTendencyQCl,
				&dTendencyQCl2);

			dataTracer[3][k][i][j] += dDeltaT * dRho * dTendencyQCl;
			dataTracer[4][k][i][j] += dDeltaT * dRho * dTendencyQCl2;
		}
		}
		}
	}

	// Call up the stack to update performance time
	WorkflowProcess::Perform(time);
}

