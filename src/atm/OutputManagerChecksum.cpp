///////////////////////////////////////////////////////////////////////////////
///
///	\file    OutputManagerChecksum.cpp
///	\author  Paul Ullrich
///	\version May 15, 2013
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

#include "OutputManagerChecksum.h"

#include "Model.h"
#include "Grid.h"

#include "Announce.h"

#include <mpi.h>

///////////////////////////////////////////////////////////////////////////////

OutputManagerChecksum::OutputManagerChecksum(
	Grid & grid,
	const Time & timeOutputFrequency
) :
	OutputManager(
		grid,
		timeOutputFrequency,
		"",
		"",
		-1)
{
}

///////////////////////////////////////////////////////////////////////////////

void OutputManagerChecksum::Output(
	const Time & time
) {
	// Get processor rank
	int nRank;
	MPI_Comm_rank(MPI_COMM_WORLD, &nRank);

	// Equation set
	const EquationSet & eqn = m_grid.GetModel().GetEquationSet();

	// Compute checksums
	DataArray1D<double> dChecksum;

	m_grid.Checksum(DataType_State, dChecksum);
	if (nRank == 0) {
		for (int c = 0; c < eqn.GetComponents(); c++) {
			Announce("..Checksum (%s): %1.15e",
				eqn.GetComponentShortName(c).c_str(), dChecksum[c]);
		}
	}

	m_grid.Checksum(DataType_Tracers, dChecksum);
	if (nRank == 0) {
		for (int c = 0; c < eqn.GetTracers(); c++) {
			Announce("..Checksum (%s): %1.15e",
				eqn.GetTracerShortName(c).c_str(), dChecksum[c]);
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

