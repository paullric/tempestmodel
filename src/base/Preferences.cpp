///////////////////////////////////////////////////////////////////////////////
///
///	\file    Preferences.cpp
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

#include "Preferences.h"
#include "Exception.h"
#include "STLStringHelper.h"

#include <iostream>
#include <fstream>
#include <cstring>
#include <cstdlib>
#include <cstdio>

////////////////////////////////////////////////////////////////////////////////

void Preferences::ParsePreferences(
	const char *szFilename
) {
	// Empty any existing preferences
	m_mapPreferences.clear();

	// Open the file
	std::fstream file(szFilename, std::ios_base::in);

	// Verify the file was properly opened
	if (!file.is_open()) {
		_EXCEPTION1("Error opening file %s", szFilename);
	}

	char szBuffer[256];

	// Read the file line by line
	for (;;) {
		file.getline(szBuffer, 256);
		if (file.eof()) {
			break;
		}

		// Ignore comment characters
		if (strlen(szBuffer) == 0) {
			continue;
		}
		if (szBuffer[0] == '#') {
			continue;
		}

		// Identify the location of the equal sign
		char* szName = szBuffer;
		char* szEqual = strchr(szBuffer, '=');
		if (NULL == szEqual) {
			_EXCEPTIONT("Malformed preferences file.");
		}

		(*szEqual) = static_cast<char>(NULL);
		szEqual++;

		// Buffer counter
		size_t i;

		// Remove whitespace from the name string
		size_t nLength = strlen(szName);
		for (i = nLength-1;; i--) {
			if (szName[i] == ' ') {
				szName[i] = static_cast<char>(NULL);
			} else {
				break;
			}

			// No name - throw exception
			if (i == 0) {
				_EXCEPTION1("Invalid name in Preference: %s\n", szName[i]);
			}
		}

		for (i = 0; i < nLength; i++) {
			if (szName[i] == ' ') {
				szName[i] = static_cast<char>(NULL);
				szName++;
				i--;
			} else {
				break;
			}
		}

		// Remove whitespace from the value string
		nLength = strlen(szEqual);
		for (i = nLength-1;; i--) {
			if (szEqual[i] == ' ') {
				szEqual[i] = static_cast<char>(NULL);
			} else {
				break;
			}

			// No value, set to NULL
			if (i == 0) {
				break;
			}
		}

		for (i = 0; i < nLength; i++) {
			if (szEqual[i] == ' ') {
				szEqual[i] = static_cast<char>(NULL);
				szEqual++;
				i--;
			} else {
				break;
			}
		}

		// Insert this string into the map
		m_mapPreferences.insert(
			PairType(std::string(szName), std::string(szEqual))
		);
	}
}

////////////////////////////////////////////////////////////////////////////////

const std::string *
Preferences::GetPreference_NoThrow(
	const char *szName
) const {
	ConstIterator iterIndex = m_mapPreferences.find(std::string(szName));

	if (iterIndex == m_mapPreferences.end()) {
		return NULL;
	} else {
		return &(iterIndex->second);
	}
}

////////////////////////////////////////////////////////////////////////////////

std::string
Preferences::GetPreferenceAsString_NoThrow(
	const char* szName,
	int *iError
) const {
	ConstIterator iterIndex = m_mapPreferences.find(std::string(szName));

	if (iterIndex == m_mapPreferences.end()) {
		if (iError != NULL) {
			(*iError) = 1;
		}
		return "";
	} else {
		if (iError != NULL) {
			(*iError) = 0;
		}
		return iterIndex->second;
	}
}

////////////////////////////////////////////////////////////////////////////////

std::string
Preferences::GetPreferenceAsString_NoCase_NoThrow(
	const char *szName,
	int *iError
) const {
	ConstIterator iterIndex = m_mapPreferences.find(std::string(szName));

	if (iterIndex == m_mapPreferences.end()) {
		if (iError != NULL) {
			(*iError) = 1;
		}
		return "";
	} else {
		if (iError != NULL) {
			(*iError) = 0;
		}
		std::string strOut = iterIndex->second;
		STLStringHelper::ToLower(strOut);
		return strOut;
	}
}

////////////////////////////////////////////////////////////////////////////////

const char *
Preferences::GetPreferenceAsCharArray_NoThrow(
	const char *szName
) const {
	const std::string *pstrPref = GetPreference_NoThrow(szName);
	if (pstrPref == NULL) {
		return NULL;
	} else {
		return pstrPref->c_str();
	}
}

////////////////////////////////////////////////////////////////////////////////

double Preferences::GetPreferenceAsDouble_NoThrow(
	const char *szName,
	int *iError
) const {
	const std::string *pstrPref = GetPreference_NoThrow(szName);
	if (pstrPref == NULL) {
		if (iError != NULL) {
			(*iError) = 1;
		}
		return 0.0;
	} else {
		if (iError != NULL) {
			(*iError) = 0;
		}
		return atof(GetPreference(szName)->c_str());
	}
}

////////////////////////////////////////////////////////////////////////////////

int Preferences::GetPreferenceAsInt_NoThrow(
	const char *szName,
	int *iError
) const {
	const std::string *pstrPref = GetPreference_NoThrow(szName);
	if (pstrPref == NULL) {
		if (iError != NULL) {
			(*iError) = 1;
		}
		return 0;
	} else {
		if (iError != NULL) {
			(*iError) = 0;
		}
		return atoi(GetPreference(szName)->c_str());
	}
}

////////////////////////////////////////////////////////////////////////////////

unsigned int Preferences::GetPreferenceAsUInt_NoThrow(
	const char *szName,
	int *iError
) const {
	const std::string* pstrPref = GetPreference_NoThrow(szName);
	if (pstrPref == NULL) {
		if (iError != NULL) {
			(*iError) = 1;
		}
		return 0;
	} else {
		if (iError != NULL) {
			(*iError) = 0;
		}
		return static_cast<unsigned int>(
			atoi(GetPreference(szName)->c_str())
		);
	}
}

////////////////////////////////////////////////////////////////////////////////

const std::string *
Preferences::GetPreference(
	const char *szName
) const {
	ConstIterator iterIndex = m_mapPreferences.find(std::string(szName));

	if (iterIndex == m_mapPreferences.end()) {
		_EXCEPTION1("Invalid preference:  %s not found.", szName);
	}

	return &(iterIndex->second);
}

////////////////////////////////////////////////////////////////////////////////

std::string
Preferences::GetPreferenceAsString(
	const char *szName
) const {
	return *(GetPreference(szName));
}

////////////////////////////////////////////////////////////////////////////////

std::string
Preferences::GetPreferenceAsString_NoCase(
	const char *szName
) const {
	std::string strOut = *(GetPreference(szName));
	STLStringHelper::ToLower(strOut);
	return strOut;
}

////////////////////////////////////////////////////////////////////////////////

const char *
Preferences::GetPreferenceAsCharArray(
	const char *szName
) const {
	return GetPreference(szName)->c_str();
}

////////////////////////////////////////////////////////////////////////////////

double Preferences::GetPreferenceAsDouble(
	const char *szName
) const {
	return atof(GetPreference(szName)->c_str());
}

////////////////////////////////////////////////////////////////////////////////

int Preferences::GetPreferenceAsInt(
	const char *szName
) const {
	return atoi(GetPreference(szName)->c_str());
}

////////////////////////////////////////////////////////////////////////////////

unsigned int Preferences::GetPreferenceAsUInt(
	const char *szName
) const {
	return static_cast<unsigned int>(
		atoi(GetPreference(szName)->c_str())
	);
}

////////////////////////////////////////////////////////////////////////////////
