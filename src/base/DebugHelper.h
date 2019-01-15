///////////////////////////////////////////////////////////////////////////////
///
///	\file    DebugHelper.h
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

#ifndef _DEBUGHELPER_H_
#define _DEBUGHELPER_H_

///////////////////////////////////////////////////////////////////////////////

#include "Exception.h"

///////////////////////////////////////////////////////////////////////////////

#ifdef DEBUG

#include <cstdio>

#define _DEBUGT(x) {fprintf(stderr, x);}
#define _ASSERT(x) {if (!(x)) {_EXCEPTIONT("ASSERT FAIL");}}

#else

#define _DEBUGT(x)
#define _ASSERT(x)

#endif

///////////////////////////////////////////////////////////////////////////////

#endif

