///////////////////////////////////////////////////////////////////////////////
///
///	\file    UserDataMeta.h
///	\author  Paul Ullrich
///	\version May 23, 2016
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

#ifndef _USERDATAMETA_H_
#define _USERDATAMETA_H_

#include "Exception.h"
#include "DataArray1D.h"

#include <vector>
#include <string>

///////////////////////////////////////////////////////////////////////////////

class PhysicalConstants;

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Meta data describing additional auxiliary user data.
///	</summary>
class UserDataMeta {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	UserDataMeta() :
		m_nUserData2DItemCount(0)
	{ }

public:
	///	<summary>
	///		Insert a new 2D user data item.
	///	</summary>
	void InsertDataItem2D(
		std::string strUserData2DItemName
	) {
		m_nUserData2DItemCount++;
		m_strUserData2DItemNames.push_back(strUserData2DItemName);
	}

public:
	///	<summary>
	///		Get the name of the specified user data item.
	///	</summary>
	inline const std::string & GetUserData2DItemName(int ix) const {
		return m_strUserData2DItemNames[ix];
	}

	///	<summary>
	///		Get the number of user data items.
	///	</summary>
	inline int GetUserData2DItemCount() const {
		return m_nUserData2DItemCount;
	}

private:
	///	<summary>
	///		Number of components per node.
	///	</summary>
	int m_nUserData2DItemCount;

	///	<summary>
	///		Short name of each component.
	///	</summary>
	std::vector<std::string> m_strUserData2DItemNames;
};

///////////////////////////////////////////////////////////////////////////////

#endif

