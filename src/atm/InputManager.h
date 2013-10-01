///////////////////////////////////////////////////////////////////////////////
///
///	\file    InputManager.h
///	\author  Paul Ullrich
///	\version September 30, 2013
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

#ifndef _INPUTMANAGER_H_
#define _INPUTMANAGER_H_

///////////////////////////////////////////////////////////////////////////////

class Grid;

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Interface class for handling initial condition input.
///	</summary>
class InputManager {

public:
	///	<summary>
	///		Initialize the grid data from a file.
	///	</summary>
	virtual void Input(
		Grid & grid
	) const = 0;
};

///////////////////////////////////////////////////////////////////////////////

#endif
