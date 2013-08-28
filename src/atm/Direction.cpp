///////////////////////////////////////////////////////////////////////////////
///
///	\file    Direction.cpp
///	\author  Paul Ullrich
///	\version August 8, 2013
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

#include "Direction.h"

#include "Exception.h"

///////////////////////////////////////////////////////////////////////////////

Direction DirectionOpposite(
	Direction dir
) {
	if (dir == Direction_Unreachable) {
		return Direction_Unreachable;

	} else if (dir == Direction_Middle) {
		return Direction_Middle;

	} else if (dir == Direction_Right) {
		return Direction_Left;

	} else if (dir == Direction_Top) {
		return Direction_Bottom;

	} else if (dir == Direction_Left) {
		return Direction_Right;

	} else if (dir == Direction_Bottom) {
		return Direction_Top;

	} else if (dir == Direction_TopRight) {
		return Direction_BottomLeft;

	} else if (dir == Direction_TopLeft) {
		return Direction_BottomRight;

	} else if (dir == Direction_BottomLeft) {
		return Direction_TopRight;

	} else if (dir == Direction_BottomRight) {
		return Direction_TopLeft;

	} else {
		_EXCEPTIONT("Invalid direction");
	}
}

///////////////////////////////////////////////////////////////////////////////

void DirectionIncrement(
	Direction dir,
	int & i,
	int & j
) {
	if (dir == Direction_Unreachable) {
		_EXCEPTIONT("Invalid direction");

	} else if (dir == Direction_Middle) {

	} else if (dir == Direction_Right) {
		i++;

	} else if (dir == Direction_Top) {
		j++;

	} else if (dir == Direction_Left) {
		i--;

	} else if (dir == Direction_Bottom) {
		j--;

	} else if (dir == Direction_TopRight) {
		i++;
		j++;

	} else if (dir == Direction_TopLeft) {
		i--;
		j++;

	} else if (dir == Direction_BottomLeft) {
		i--;
		j--;

	} else if (dir == Direction_BottomRight) {
		i++;
		j--;

	} else {
		_EXCEPTIONT("Invalid direction");
	}
}

///////////////////////////////////////////////////////////////////////////////

