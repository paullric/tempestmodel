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

#ifndef _CONSOLIDATIONSTATUS_H_
#define _CONSOLIDATIONSTATUS_H_

#include "DataType.h"

#include "mpi.h"

#include <vector>
#include <map>

///////////////////////////////////////////////////////////////////////////////

class Grid;

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		An array for holding the current consolidation status, which is used
///		for MPI reduce of data to the root processor.
///	</summary>
class ConsolidationStatus {
	friend class Grid;

public:
	typedef std::map<DataType, int>          DataTypeIndexMap;
	typedef DataTypeIndexMap::const_iterator DataTypeIndexMapIterator;
	typedef DataTypeIndexMap::value_type     DataTypeIndexMapPair;

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	ConsolidationStatus(
		const Grid & grid,
		const std::vector<DataType> & vecDataTypes
	);

public:
	///	<summary>
	///		Generate a MPI_TAG.
	///	</summary>
	static int GenerateTag(
		int ixPatch,
		DataType eDataType
	) {
		return (ixPatch << 4) + static_cast<int>(eDataType);
	}

	///	<summary>
	///		Parse the MPI_TAG for DataType and patch index.
	///	</summary>
	static void ParseTag(
		int nTag,
		int & ixPatch,
		DataType & eDataType
	) {
		eDataType = static_cast<DataType>(nTag & 0xF);
		ixPatch = (nTag >> 4);
	}

public:
	///	<summary>
	///		Return the number of DataTypes in this structure.
	///	</summary>
	int GetDataTypeCount() const {
		return m_mapDataTypes.size();
	}

	///	<summary>
	///		Determine if the specified DataType is to be consolidated.
	///	</summary>
	bool Contains(DataType eDataType) const;

	///	<summary>
	///		Determine if all receives have completed.
	///	</summary>
	bool Done() const;

protected:
	/// <summary>
	///		Set the receive status for the given DataType and patch.
	///	</summary>
	void SetReceiveStatus(
		int ixPatch,
		DataType eDataType
	);

	///	<summary>
	///		Get the specified send request.
	///	</summary>
	MPI_Request * GetNextSendRequest();

private:
	///	<summary>
	///		Data types stored in this data structure.
	///	</summary>
	DataTypeIndexMap m_mapDataTypes;

	///	<summary>
	///		Total number of receives set.
	///	</summary>
	std::vector<int> m_vecReceiveStatusCount;

	///	<summary>
	///		Vector of receive status.
	///	</summary>
	std::vector< std::vector<bool> > m_vecReceiveStatus;

	///	<summary>
	///		Current send request.
	///	</summary>
	int m_nCurrentSendRequest;

	///	<summary>
	///		Vector of send requests.
	///	</summary>
	std::vector<MPI_Request> m_vecSendRequests;
};

///////////////////////////////////////////////////////////////////////////////

#endif

