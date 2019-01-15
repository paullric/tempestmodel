///////////////////////////////////////////////////////////////////////////////
///
///	\file    Preferences.h
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

#ifndef _PREFERENCES_H_
#define _PREFERENCES_H_

////////////////////////////////////////////////////////////////////////////////

#include <map>
#include <string>

////////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A class for handing preferences read from a file.
///	</summary>
///	<remarks>
///		Known issues:  These functions are vulnerable to buffer overflow, but
///		for scientific purposes this should not be problematic.
///	</remarks>
class Preferences {

	public:
		///	<summary>
		///		A dictionary of name-value pairs used for property lookup.
		///	</summary>
		typedef std::map<std::string, std::string> PreferencesMap;

		typedef PreferencesMap::iterator           Iterator;

		typedef PreferencesMap::const_iterator     ConstIterator;

		///	<summary>
		///		A name-value pair.
		///	</summary>
		typedef std::pair<std::string, std::string> PairType;

	public:
		///	<summary>
		///		Constructor.
		///	</summary>
		Preferences()
		{ }

		///	<summary>
		///		Constructor that specifies the preferences file.
		///	</summary>
		Preferences(const char *szFilename) {
			ParsePreferences(szFilename);
		}

	public:
		///	<summary>
		///		Parse the preferences file given by <c>szFilename</c> and
		///		load the contents into this Preferences object.
		///	</summary>
		void ParsePreferences(const char *szFilename);

	private:
		///	<summary>
		///		Get the specified preference as a pointer to a std::string.
		///		If the preference is not found, a NULL pointer will be
		///		returned.
		///	</summary>
		const std::string * GetPreference_NoThrow(
			const char *szName
		) const;

	public:
		///	<summary>
		///		Get the specified preference as a std::string.
		///	</summary>
		std::string GetPreferenceAsString_NoThrow(
			const char *szName,
			int *iError = NULL
		) const;

		///	<summary>
		///		Get the specified preference as a std::string where all
		///		characters appear in lower case.
		///	</summary>
		std::string GetPreferenceAsString_NoCase_NoThrow(
			const char *szName,
			int *iError = NULL
		) const;

		///	<summary>
		///		Get the specified preference as a pointer to a char array.
		///		If the preference is not found, a NULL pointer will be
		///		returned.
		///	</summary>
		const char * GetPreferenceAsCharArray_NoThrow(
			const char *szName
		) const;

		///	<summary>
		///		Get the specified preference as a double.  If the preference is
		///		not found a zero will be returned and the variable pointed to
		///		by iError, if it exists, will be set to 1 (further, it will
		///		be set to 0 if the preference is found).
		///	</summary>
		double GetPreferenceAsDouble_NoThrow(
			const char *szName,
			int *iError = NULL
		) const;

		///	<summary>
		///		Get the specified preference as an Int.  If the preference is
		///		not found a zero will be returned and the variable pointed to
		///		by iError, if it exists, will be set to 1 (further, it will
		///		be set to 0 if the preference is found).
		///	</summary>
		int GetPreferenceAsInt_NoThrow(
			const char *szName,
			int *iError = NULL
		) const;

		///	<summary>
		///		Get the specified preference as a uInt.  If the preference is
		///		not found a zero will be returned and the variable pointed to
		///		by iError, if it exists, will be set to 1 (further, it will
		///		be set to 0 if the preference is found).
		///	</summary>
		unsigned int GetPreferenceAsUInt_NoThrow(
			const char *szName,
			int *iError = NULL
		) const;

	private:
		///	<summary>
		///		Get the specified preference as a pointer to a string.  If the
		///		preference cannot be found, an exception is thrown.
		///	</summary>
		const std::string * GetPreference(
			const char *szName
		) const;

	public:
		///	<summary>
		///		Get the specified preference as a string.  If the preference
		///		cannot be found, an exception is thrown.
		///	</summary>
		std::string GetPreferenceAsString(
			const char *szName
		) const;

		///	<summary>
		///		Get the specified preference as a string where all characters
		///		are in lower case.  If the preference cannot be found, an
		///		exception is thrown.
		///	</summary>
		std::string GetPreferenceAsString_NoCase(
			const char *szName
		) const;

		///	<summary>
		///		Get the specified preference as a character array.  If the
		///		preference cannot be found, an exception is thrown.
		///	</summary>
		const char * GetPreferenceAsCharArray(
			const char *szName
		) const;

		///	<summary>
		///		Get the specified preference as a double.  If the
		///		preference cannot be found, an exception is thrown.
		///	</summary>
		double GetPreferenceAsDouble(
			const char *szName
		) const;

		///	<summary>
		///		Get the specified preference as an Int.  If the
		///		preference cannot be found, an exception is thrown.
		///	</summary>
		int GetPreferenceAsInt(
			const char *szName
		) const;

		///	<summary>
		///		Get the specified preference as a uInt.  If the preference
		///		cannot be found, an exception is thrown.
		///	</summary>
		unsigned int GetPreferenceAsUInt(
			const char *szName
		) const;

	private:
		///	<summary>
		///		The dictionary of name-value pairs used for property lookup.
		///	</summary>
		PreferencesMap m_mapPreferences;

};

////////////////////////////////////////////////////////////////////////////////

#endif

