///////////////////////////////////////////////////////////////////////////////
///
///	\file    ConsolidationStatus.h
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

#include "ConsolidationStatus.h"
#include "Grid.h"

///////////////////////////////////////////////////////////////////////////////

ConsolidationStatus::ConsolidationStatus(
	const Grid & grid,
	const std::vector<DataTypeLocationPair> & vecDataTypes
) :
	m_nCurrentSendRequest(0)
{

	// Get process id
	int nRank;
	MPI_Comm_rank(MPI_COMM_WORLD, &nRank);

	// Built map of data types and indices
	for (int m = 0; m < vecDataTypes.size(); m++) {
		m_mapDataTypes.insert(
			ConsolidationIndexMapPair(vecDataTypes[m], m));
	}

	// Built map of send requests
	m_vecSendRequests.resize(
		grid.GetActivePatchCount() * m_mapDataTypes.size());

	// Non-root process stores no receive statuses
	if (nRank != 0) {
		return;
	}

	// Number of received data packets per data type
	m_vecReceiveStatusCount.resize(vecDataTypes.size());

	// Reset the number of received data packets
	for (int m = 0; m < m_vecReceiveStatusCount.size(); m++) {
		m_vecReceiveStatusCount[m] = 0;
	}

	// Status of each DataTypeLocationPair / patch pair
	m_vecReceiveStatus.resize(grid.GetPatchCount());
	for (int n = 0; n < grid.GetPatchCount(); n++) {
		m_vecReceiveStatus[n].resize(vecDataTypes.size());
		for (int m = 0; m < m_vecReceiveStatus[n].size(); m++) {
			m_vecReceiveStatus[n][m] = false;
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

bool ConsolidationStatus::Contains(
	DataType eDataType,
	DataLocation eDataLocation
) const {
	DataTypeLocationPair key(eDataType, eDataLocation);

	if (m_mapDataTypes.find(key) != m_mapDataTypes.end()) {
		return true;
	}
	return false;
}

///////////////////////////////////////////////////////////////////////////////

bool ConsolidationStatus::Done() const {
	for (int m = 0; m < m_vecReceiveStatusCount.size(); m++) {
		if (m_vecReceiveStatusCount[m] != m_vecReceiveStatus.size()) {
			return false;
		}
	}
	return true;
}

///////////////////////////////////////////////////////////////////////////////

void ConsolidationStatus::SetReceiveStatus(
	int ixPatch,
	DataType eDataType,
	DataLocation eDataLocation
) {
	DataTypeLocationPair key(eDataType, eDataLocation);

	ConsolidationIndexMapIterator iter = m_mapDataTypes.find(key);
	if (iter == m_mapDataTypes.end()) {
		_EXCEPTIONT("Invalid DataType");
	}
	if ((ixPatch < 0) || (ixPatch >= m_vecReceiveStatus.size())) {
		_EXCEPTIONT("Patch index out of range.");
	}
	if (m_vecReceiveStatus[ixPatch][iter->second]) {
		_EXCEPTIONT("Receive status already set");
	}
	m_vecReceiveStatusCount[iter->second]++;
	m_vecReceiveStatus[ixPatch][iter->second] = true;
}

///////////////////////////////////////////////////////////////////////////////

MPI_Request * ConsolidationStatus::GetNextSendRequest() {
	if (m_nCurrentSendRequest >= m_vecSendRequests.size()) {
		_EXCEPTIONT("Too many MPI_Isend requests required.");
	}
	MPI_Request * req = &(m_vecSendRequests[m_nCurrentSendRequest]);
	m_nCurrentSendRequest++;
	return req;
}

///////////////////////////////////////////////////////////////////////////////

