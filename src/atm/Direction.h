///////////////////////////////////////////////////////////////////////////////
///
///	\file    Direction.h
///	\author  Paul Ullrich
///	\version March 24, 2013
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

#ifndef _DIRECTION_H_
#define _DIRECTION_H_

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Index of an invalid panel.
///	</summary>
static const int InvalidPanel = (-1);

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Directions that can be exchanged between processors.
///	</summary>
enum Direction {
	Direction_Unreachable = (-2),
	Direction_Middle = (-1),

	Direction_Right  = 0,
	Direction_Top    = 1,
	Direction_Left   = 2,
	Direction_Bottom = 3,

	Direction_TopRight    = 4,
	Direction_TopLeft     = 5,
	Direction_BottomLeft  = 6,
	Direction_BottomRight = 7,
};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Return the direction opposite to the specified direction.
///	</summary>
Direction DirectionOpposite(
	Direction dir
);

///	<summary>
///		Increment the coordinates by one in the specified direction.
///	</summary>
void DirectionIncrement(
	Direction dir,
	int & i,
	int & j
);

///////////////////////////////////////////////////////////////////////////////

#endif
