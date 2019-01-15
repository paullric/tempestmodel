///////////////////////////////////////////////////////////////////////////////
///
///	\file    StringHelper.h
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

#ifndef _STRINGHELPER_H_
#define _STRINGHELPER_H_

#include <string.h>

///	<summary>
///		This class exposes some functionality which can be used to perform
///		certain special operations on strings.
///	</summary>
class StringHelper {

///////////////////////////////////////////////////////////////////////////////

private:
StringHelper() { }

public:

///////////////////////////////////////////////////////////////////////////////

inline static void ToLower(char * szString) {
	while (*szString) {
		(*szString) = tolower(*szString);
		szString++;
	}
}

///////////////////////////////////////////////////////////////////////////////

inline static void ToUpper(char * szString) {
	while (*szString) {
		(*szString) = toupper(*szString);
		szString++;
	}
}

///////////////////////////////////////////////////////////////////////////////

inline static const char * IgnoreWhitespace(const char * szString) {
	while (
		((*szString) == ' ') ||
		((*szString) == '\t') ||
		((*szString) == '\r') ||
		((*szString) == '\n')
	) {
		szString++;
	}

	return szString;
}

///////////////////////////////////////////////////////////////////////////////

inline static char * IgnoreWhitespace(char *szString) {
	while (
		((*szString) == ' ') ||
		((*szString) == '\t') ||
		((*szString) == '\r') ||
		((*szString) == '\n')
	) {
		szString++;
	}

	return szString;
}

///////////////////////////////////////////////////////////////////////////////

inline static const char * FindWhitespace(const char *szString) {
	while (
		((*szString) != ' ') &&
		((*szString) != '\t') &&
		((*szString) != '\r') &&
		((*szString) != '\n') &&
		((*szString) != '\0')
	) {
		szString++;
	}

	return szString;
}

///////////////////////////////////////////////////////////////////////////////

inline static char * FindWhitespace(char *szString) {
	while (
		((*szString) != ' ') &&
		((*szString) != '\t') &&
		((*szString) != '\r') &&
		((*szString) != '\n') &&
		((*szString) != '\0')
	) {
		szString++;
	}

	return szString;
}

///////////////////////////////////////////////////////////////////////////////

static int StringToHex(const char *szString, int nLength = 0) {
	int nValue;
	int nI;

	if (nLength == 0) {
		nLength = strlen(szString);
	}

	nValue = 0;
	for (nI = nLength-1; nI >= 0; nI--) {
		nValue = nValue << 4;
		if ((szString[nI] >= '0') && (szString[nI] <= '9')) {
			nValue = nValue + szString[nI] - '0';
		} else if ((szString[nI] >= 'A') && (szString[nI] <= 'F')) {
			nValue = nValue + 10 + szString[nI] - 'A';
		} else if ((szString[nI] >= 'a') && (szString[nI] <= 'f')) {
			nValue = nValue + 10 + szString[nI] - 'a';
		} else {
			return -1;
		}
	}

	return nValue;
}

///////////////////////////////////////////////////////////////////////////////

};

#endif
