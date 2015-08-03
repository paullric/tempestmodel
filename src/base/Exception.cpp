///////////////////////////////////////////////////////////////////////////////
///
///	\file    Exception.cpp
///	\author  Paul Ullrich
///	\version July 26, 2010
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

#include "Exception.h"

///////////////////////////////////////////////////////////////////////////////

Exception::Exception(
	const char * szFile,
	unsigned int uiLine
) :
	m_strText("General exception"),
	m_strFile(szFile),
	m_uiLine(uiLine)
{ }

///////////////////////////////////////////////////////////////////////////////

Exception::Exception(
	const char * szFile,
	unsigned int uiLine,
	const char * szText,
	...
) :
	m_strFile(szFile),
	m_uiLine(uiLine)
{
	char szBuffer[ExceptionBufferSize];

	va_list arguments;

	// Initialize the argument list
	va_start(arguments, szText);

	// Write to string
	vsprintf(szBuffer, szText, arguments);

	m_strText = szBuffer;

	// Cleans up the argument list
	va_end(arguments);
}

///////////////////////////////////////////////////////////////////////////////

std::string Exception::ToString() const {
	std::string strReturn;

	char szBuffer[128];

	// Preamble
	sprintf(szBuffer, "EXCEPTION (");
	strReturn.append(szBuffer);

	// File name
	strReturn.append(m_strFile);

	// Line number
	sprintf(szBuffer, ", Line %u) ", m_uiLine);
	strReturn.append(szBuffer);

	// Text
	strReturn.append(m_strText);

	return strReturn;
}

///////////////////////////////////////////////////////////////////////////////

