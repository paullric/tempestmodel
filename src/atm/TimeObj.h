///////////////////////////////////////////////////////////////////////////////
///
///	\file    TimeObj.h
///	\author  Paul Ullrich
///	\version October 31, 2013
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

#ifndef _TIMEOBJ_H_
#define _TIMEOBJ_H_

#include "Exception.h"

#include <string>
#include <cstdlib>
#include <cmath>

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A class for storing a time as year, month, day and seconds.
///	</summary>
class Time {

public:
	///	<summary>
	///		Length of time (in seconds) which is considered "concurrent"
	///	</summary>
	static const double SecondsEpsilon;

public:
	///	<summary>
	///		Type of calendar.
	///	</summary>
	enum CalendarType {
		CalendarNone,
		CalendarNoLeap,
	};

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	Time() :
		m_iYear(0),
		m_iMonth(0),
		m_iDay(0),
		m_dSeconds(0),
		m_eCalendarType(CalendarNoLeap)
	{ }

	///	<summary>
	///		Constructor.
	///	</summary>
	Time(
		int iYear,
		int iMonth,
		int iDay,
		double dSeconds
	) :
		m_iYear(iYear),
		m_iMonth(iMonth),
		m_iDay(iDay),
		m_dSeconds(dSeconds),
		m_eCalendarType(CalendarNoLeap)
	{
		NormalizeTime();
	}

	///	<summary>
	///		Constructor from std::string.
	///	</summary>
	Time(const std::string & strTime) :
		m_iYear(0),
		m_iMonth(0),
		m_iDay(0),
		m_dSeconds(0),
		m_eCalendarType(CalendarNoLeap)
	{
		FromShortString(strTime);
	}

public:
	///	<summary>
	///		Determine if two Times occur on the same date.
	///	</summary>
	inline bool IsSameDate(const Time & time) {
		if (m_eCalendarType != time.m_eCalendarType) {
			return false;
		}

		if ((m_iYear  == time.m_iYear) &&
		    (m_iMonth == time.m_iMonth) &&
		    (m_iDay   == time.m_iDay)
		) {
			return true;
		}
		return false;
	}

	///	<summary>
	///		Equality between Times.
	///	</summary>
	bool operator==(const Time & time) const;

	///	<summary>
	///		Less-than between Times.
	///	</summary>
	bool operator<(const Time & time) const;

	///	<summary>
	///		Greater-than comparator between Times.
	///	</summary>
	bool operator>(const Time & time) const;

	///	<summary>
	///		Less-than-or-equal comparator between Times.
	///	</summary>
	bool operator<=(const Time & time) const {
		return !((*this) > time);
	}

	///	<summary>
	///		Greater-than-or-equal comparator between Times.
	///	</summary>
	bool operator>=(const Time & time) const {
		return !((*this) < time);
	}

protected:
	///	<summary>
	///		Verify that the Time is in accordance with the calendar.
	///	</summary>
	void VerifyTime();

	///	<summary>
	///		Normalize the Time, in accordance with the Calendar.
	///	</summary>
	void NormalizeTime();

public:
	///	<summary>
	///		Add a number of seconds to the Time.
	///	</summary>
	inline void AddSeconds(double dSeconds) {
		m_dSeconds += dSeconds;

		NormalizeTime();
	}

	///	<summary>
	///		Add a number of seconds to the Time (operator version).
	///	</summary>
	inline void operator+=(double dSeconds) {
		AddSeconds(dSeconds);
	}

	///	<summary>
	///		Add a number of days to the Time.
	///	</summary>
	inline void AddDays(int nDays) {
		m_iDay += nDays;

		NormalizeTime();
	}

	///	<summary>
	///		Add a number of months to the Time.
	///	</summary>
	inline void AddMonths(int nMonths) {
		m_iMonth += nMonths;

		NormalizeTime();
	}

	///	<summary>
	///		Add a number of years to the Time.
	///	</summary>
	inline void AddYears(int nYears) {
		m_iYear += nYears;

		NormalizeTime();
	}

public:
	///	<summary>
	///		Add a number of seconds to the Time to produce a new Time value.
	///	</summary>
	Time operator+(double dSeconds) const;

	///	<summary>
	///		Determine the number of seconds between two Times.
	///	</summary>
	double operator-(const Time & time) const;

public:
	///	<summary>
	///		Get the year.
	///	</summary>
	inline int GetYear() const {
		return m_iYear;
	}

	///	<summary>
	///		Get the month.
	///	</summary>
	inline int GetMonth() const {
		return m_iMonth;
	}

	///	<summary>
	///		Get the day.
	///	</summary>
	inline int GetDay() const {
		return m_iDay;
	}

	///	<summary>
	///		Get the number of seconds in the day.
	///	</summary>
	inline double GetSeconds() const {
		return m_dSeconds;
	}

	///	<summary>
	///		Get the CalendarType.
	///	</summary>
	inline CalendarType GetCalendarType() const {
		return m_eCalendarType;
	}

	///	<summary>
	///		Set the year.
	///	</summary>
	inline void SetYear(int iYear) {
		m_iYear = iYear;
	}

	///	<summary>
	///		Set the month.
	///	</summary>
	inline void SetMonth(int iMonth) {
		m_iMonth = iMonth;
	}

	///	<summary>
	///		Set the day.
	///	</summary>
	inline void SetDay(int iDay) {
		m_iDay = iDay;
	}

	///	<summary>
	///		Set the number of seconds in the day.
	///	</summary>
	inline void SetSeconds(double dSeconds) {
		m_dSeconds = dSeconds;
	}

public:
	///	<summary>
	///		Get the Time as a string, only including the date: "YYYY-MM-DD"
	///	</summary>
	std::string ToDateString() const;

	///	<summary>
	///		Get the Time as a short string: "YYYY-MM-DD-SSSSS"
	///	</summary>
	std::string ToShortString() const;

	///	<summary>
	///		Get the Time as a full-length string showing hours, days
	///		and seconds: "YYYY-MM-DD hh:mm:ss"
	///	</summary>
	std::string ToString() const;

	///	<summary>
	///		Set the Time using a short string: "YYYY-MM-DD-SSSSS"
	///	</summary>
	void FromShortString(const std::string & strShortTime);

	///	<summary>
	///		Get the name of the calendar.
	///	</summary>
	std::string GetCalendarName() const;

private:
	///	<summary>
	///		The year.
	///	</summary>
	int m_iYear;

	///	<summary>
	///		The month.
	///	</summary>
	int m_iMonth;

	///	<summary>
	///		The day.
	///	</summary>
	int m_iDay;

	///	<summary>
	///		Number of seconds in the day.
	///	</summary>
	double m_dSeconds;

	///	<summary>
	///		Calendar type.
	///	</summary>
	CalendarType m_eCalendarType;
};

///////////////////////////////////////////////////////////////////////////////

#endif

