///////////////////////////////////////////////////////////////////////////////
///
///	\file    TimeObj.cpp
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

#include "TimeObj.h"
#include "Exception.h"

///////////////////////////////////////////////////////////////////////////////

const double Time::SecondsEpsilon = 1.0e-12;

///////////////////////////////////////////////////////////////////////////////

bool Time::operator==(const Time & time) const {
	if ((m_eCalendarType != time.m_eCalendarType) ||
	    (m_iYear  != time.m_iYear) ||
	    (m_iMonth != time.m_iMonth) ||
	    (m_iDay   != time.m_iDay)
	) {
		return false;
	}

	if (fabs(m_dSeconds - time.m_dSeconds) > SecondsEpsilon) {
		return false;
	}

	return true;
}

///////////////////////////////////////////////////////////////////////////////

bool Time::operator<(const Time & time) const {
	if (m_eCalendarType != time.m_eCalendarType) {
		_EXCEPTIONT(
			"Cannot compare Time objects with different calendars");
	}

	if (m_iYear < time.m_iYear) {
		return true;
	} else if (m_iYear > time.m_iYear) {
		return false;
	}

	if (m_iMonth < time.m_iMonth) {
		return true;
	} else if (m_iMonth > time.m_iMonth) {
		return false;
	}

	if (m_iDay < time.m_iDay) {
		return true;
	} else if (m_iDay > time.m_iDay) {
		return false;
	}

	return (m_dSeconds < time.m_dSeconds - SecondsEpsilon);
}

///////////////////////////////////////////////////////////////////////////////

bool Time::operator>(const Time & time) const {
	if (m_eCalendarType != time.m_eCalendarType) {
		_EXCEPTIONT(
			"Cannot compare Time objects with different calendars");
	}

	if (m_iYear > time.m_iYear) {
		return true;
	} else if (m_iYear < time.m_iYear) {
		return false;
	}

	if (m_iMonth > time.m_iMonth) {
		return true;
	} else if (m_iMonth < time.m_iMonth) {
		return false;
	}

	if (m_iDay > time.m_iDay) {
		return true;
	} else if (m_iDay < time.m_iDay) {
		return false;
	}

	return (m_dSeconds > time.m_dSeconds + SecondsEpsilon);
}

///////////////////////////////////////////////////////////////////////////////

void Time::VerifyTime() {

	// Calendar with no leap years
	if (m_eCalendarType == CalendarNoLeap) {

		const int nDaysPerMonth[]
			= {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};

		if ((m_iYear < 0) || (m_iYear > 9999)) {
			_EXCEPTIONT("Year out of range");
		}
		if ((m_iMonth < 0) || (m_iMonth > 11)) {
			_EXCEPTIONT("Month out of range");
		}
		if ((m_iDay < 0) || (m_iDay > nDaysPerMonth[m_iMonth])) {
			_EXCEPTIONT("Day out of range");
		}
		if ((m_dSeconds < 0.0) || (m_dSeconds > 86400.0)) {
			_EXCEPTIONT("Seconds out of range");
		}

	// Operation not permitted on this CalendarType
	} else {
		_EXCEPTIONT("Invalid CalendarType");
	}
}

///////////////////////////////////////////////////////////////////////////////

void Time::NormalizeTime() {

	// Calendar with no leap years
	if (m_eCalendarType == CalendarNoLeap) {

		const int nDaysPerMonth[]
			= {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};

		// Add days
		int nAddedDays = 0;
		if (m_dSeconds >= 86400.0) {
			nAddedDays =
				static_cast<int>(m_dSeconds / 86400.0);

		} else if (m_dSeconds < 0.0) {
			nAddedDays =
				static_cast<int>(- (86400.0 - m_dSeconds) / 86400.0);
		}

		m_dSeconds -= static_cast<double>(nAddedDays) * 86400.0;
		m_iDay += nAddedDays;

		// Add years
		int nAddedYears = 0;
		if (m_iMonth >= 12) {
			nAddedYears = m_iMonth / 12;
		} else if (m_iMonth < 0) {
			nAddedYears = - (12 - m_iMonth) / 12;
		}
		m_iMonth -= nAddedYears * 12;
		m_iYear += nAddedYears;

		if ((m_iMonth < 0) || (m_iMonth >= 12)) {
			_EXCEPTIONT("Logic error");
		}

		// Subtract months
		while (m_iDay < 0) {
			m_iMonth--;
			if (m_iMonth < 0) {
				m_iMonth = 11;
				m_iYear--;
			}
			m_iDay += nDaysPerMonth[m_iMonth];
		}

		// Add months
		while (m_iDay > nDaysPerMonth[m_iMonth]) {
			m_iDay -= nDaysPerMonth[m_iMonth];
			m_iMonth++;
			if (m_iMonth > 11) {
				m_iMonth = 0;
				m_iYear++;
			}
		}

		// Check that the result is ok
		if ((m_iMonth < 0) || (m_iMonth >= 12) ||
		    (m_iDay < 0) || (m_iDay >= nDaysPerMonth[m_iMonth]) ||
		    (m_dSeconds < 0.0) || (m_dSeconds >= 86400.0)
		) {
			_EXCEPTION4("Logic error: %i %i %i %1.15e",
				m_iYear, m_iMonth, m_iDay, m_dSeconds);
		}

	// Operation not permitted on this CalendarType
	} else {
		_EXCEPTIONT("Invalid CalendarType");
	}
}

///////////////////////////////////////////////////////////////////////////////

Time Time::operator+(double dSeconds) const {
	Time timeNew = (*this);
	timeNew += dSeconds;
	return timeNew;
}

///////////////////////////////////////////////////////////////////////////////

double Time::operator-(const Time & time) const {

	double dDeltaSeconds = 0;

	// Calendar with no leap years
	if (m_eCalendarType == CalendarNoLeap) {

		const int nDaysPerMonth[]
			= {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};

		dDeltaSeconds =
			+ static_cast<double>(m_iYear - time.m_iYear) * 31536000.0
			+ (m_dSeconds - time.m_dSeconds);

		// Remove intervening months
		if (m_iMonth > time.m_iMonth) {
			dDeltaSeconds +=
				+ 86400.0 * static_cast<double>(
					nDaysPerMonth[time.m_iMonth] - time.m_iDay - 1)
				+ 86400.0 * static_cast<double>(m_iDay);

			for (int i = time.m_iMonth+1; i < m_iMonth; i++) {
				dDeltaSeconds +=
					static_cast<double>(nDaysPerMonth[i]) * 86400.0;
			}

		// Remove intervening months
		} else if (m_iMonth < time.m_iMonth) {
			dDeltaSeconds -=
				+ 86400.0 * static_cast<double>(
					nDaysPerMonth[m_iMonth] - m_iDay - 1)
				+ 86400.0 * static_cast<double>(time.m_iDay);

			for (int i = m_iMonth+1; i < time.m_iMonth; i++) {
				dDeltaSeconds -=
					86400.0 * static_cast<double>(nDaysPerMonth[i]);
			}

		// Same month, just take the difference in days
		} else {
			dDeltaSeconds +=
				86400.0 * static_cast<double>(m_iDay - time.m_iDay);
		}

	// Operation not permitted on this CalendarType
	} else {
		_EXCEPTIONT("Invalid CalendarType");
	}
	
	return dDeltaSeconds;
}

///////////////////////////////////////////////////////////////////////////////

std::string Time::ToDateString() const {
	char szBuffer[100];

	sprintf(szBuffer, "%04i-%02i-%02i",
		m_iYear, m_iMonth+1, m_iDay+1);

	return std::string(szBuffer);
}

///////////////////////////////////////////////////////////////////////////////

std::string Time::ToShortString() const {
	char szBuffer[100];

	sprintf(szBuffer, "%04i-%02i-%02i-%05i",
		m_iYear, m_iMonth+1, m_iDay+1, static_cast<int>(m_dSeconds));

	return std::string(szBuffer);
}

///////////////////////////////////////////////////////////////////////////////

std::string Time::ToString() const {
	char szBuffer[100];

	sprintf(szBuffer, "%04i-%02i-%02i %02i:%02i:%02i",
		m_iYear, m_iMonth+1, m_iDay+1,
		static_cast<int>(m_dSeconds / 3600.0),
		static_cast<int>(fmod(m_dSeconds, 3600.0) / 60.0),
		static_cast<int>(fmod(m_dSeconds, 60.0)));

	return std::string(szBuffer);
}

///////////////////////////////////////////////////////////////////////////////

void Time::FromShortString(const std::string & strShortTime) {
	if (strShortTime.length() != 16) {
		_EXCEPTIONT("Invalid Time ShortString");
	}
	m_iYear = atoi(strShortTime.c_str());
	m_iMonth = atoi(strShortTime.c_str()+5) - 1;
	m_iDay = atoi(strShortTime.c_str()+8) - 1;
	m_dSeconds = atof(strShortTime.c_str()+11);

	VerifyTime();
}

///////////////////////////////////////////////////////////////////////////////

std::string Time::GetCalendarName() const {
	if (m_eCalendarType == CalendarNone) {
		return std::string("none");
	} else if (m_eCalendarType == CalendarNoLeap) {
		return std::string("noleap");
	} else {
		_EXCEPTIONT("Invalid CalendarType");
	}
}

///////////////////////////////////////////////////////////////////////////////

