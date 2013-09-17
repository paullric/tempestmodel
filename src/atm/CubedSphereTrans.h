///////////////////////////////////////////////////////////////////////////////
///
///	\file    CubedSphereTrans.h
///	\author  Paul Ullrich
///	\version August 11, 2010
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

#ifndef _CUBEDSPHERETRANS_H_
#define _CUBEDSPHERETRANS_H_

#include "Exception.h"
#include "Direction.h"

#include <cmath>
#include <iostream>

////////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Coordinate transformation functions on the CubedSphere grid.
///	</summary>
class CubedSphereTrans {

private:
	///	<summary>
	///		No initialization for objects of this class.
	///	</summary>
	CubedSphereTrans()
	{ }

public:
	///	<summary>
	///		Determine the Cartesian (XYZ) coordinates of a point on a sphere
	///		given its Gnomonic (XYP) coordinate.
	///	</summary>
	///	<parameters>
	///		XX       - Gnomonic X coordinate
	///		YY       - Gnomonic Y coordinate
	///		np       - Panel coordinate (0-5, 4 = north, 5 = south)
	///		xx (OUT) - Calculated X coordinate
	///		yy (OUT) - Calculated Y coordinate
	///		zz (OUT) - Calculated Z coordinate
	///	</parameters>
	static void XYZFromXYP(
		double dXX,
		double dYY,
		int np,
		double & xx,
		double & yy,
		double & zz
	);

	///	<summary>
	///		Determine the (X, Y, idest) coordinates of a source point on 
	///		panel isource.
	///	</summary>
	///	<parameters>
	///		dX_in           - Gnomonic X coordinate on source panel
	///		dY_in           - Gnomonic Y coordinate on source panel
	///		isource         - Source panel (0-5, 4 = north, 5 = south)
	///		idest           - Destination panel (0-5, 4 = north, 5 = south)
	///		dX_out (OUT)    - Gnomonic X coordinate on destination panel
	///		dY_out (OUT)    - Gnomonic Y coordinate on destination panel
	///	</parameters>
	static void XYPFromXYP(
		double dX_in,
		double dY_in,
		int isource,
		int idest,
		double & dX_out,
		double & dY_out
	);

	///	<summary>
	///		Determine the (lat, lon) coordinates of a point in gnomonic XYP
	///		coordinates.
	///	</summary>
	///	<parameters>
	///		dX        - Gnomonic X coordinate
	///		dY        - Gnomonic Y coordinate
	///		nP        - Panel coordinate (0-5, 4 = north, 5 = south)
	///		lat (OUT) - Latitude
	///		lon (OUT) - Longitude
	///	</parameters>
	static void RLLFromXYP(
		double dX,
		double dY,
		int nP,
		double & lon,
		double & lat
	);

	///	<summary>
	///		Determine the (lat, lon) coordinates of a point in equiangular ABP 
	///		coordinates.
	///	</summary>
	///	<parameters>
	///		dA        - Equiangular alpha coordinate
	///		dB        - Equiangular beta coordinate
	///		nP        - Panel coordinate (0-5, 4 = north, 5 = south)
	///		lat (OUT) - Latitude
	///		lon (OUT) - Longitude
	///	</parameters>
	static void RLLFromABP(
		double dA,
		double dB,
		int nP,
		double & lon,
		double & lat
	);

	///	<summary>
	///		Determine the gnomonic XYP coordinates of a point in (lat, lon)
	///		coordinates.
	///	</summary>
	///	<parameters>
	///		lat      - Latitude
	///		lon      - Longitude
	///		dX (OUT) - Gnomonic X coordinate
	///		dY (OUT) - Gnomonic Y coordinate
	///		nP (OUT) - Panel coordinate (0-5, 4 = north, 5 = south)
	///	</parameters>
	static void XYPFromRLL(
		double lon,
		double lat,
		double & dX,
		double & dY,
		int & nP
	);

	///	<summary>
	///		Determine the equiangular ABP coordinates of a point in (lat, lon)
	///		coordinates.
	///	</summary>
	///	<parameters>
	///		lat      - Latitude
	///		lon      - Longitude
	///		dA (OUT) - Equiangular alpha coordinate
	///		dB (OUT) - Equiangular beta coordinate
	///		nP (OUT) - Panel coordinate (0-5, 4 = north, 5 = south)
	///	</parameters>
	static void ABPFromRLL(
		double lon,
		double lat,
		double & dA,
		double & dB,
		int & nP
	);

	///	<summary>
	///		Translate a vector in spherical coordinates to a vector in
	///		equiangular coordinates.  The components of the vector field in
	///		spherical coordinates are expected in the unit basis, whereas the
	///		components in equiangular coordinates are given in the geometric
	///		basis.
	///	</summary>
	///	<parameters>
	///		dX            - Gnomonic X coordinate of translation point
	///		dY            - Gnomonic Y coordinate of translation point
	///		nP            - Gnomonic panel coordinate (0-5, 4 = north,
	///		                5 = south) of translation point
	///		dUlon         - Longitudinal component of vector field
	///		dUlat         - Latitudinal component of vector field
	///		dUalpha (OUT) - Alpha component of vector field
	///		dUbeta (OUT)  - Beta component of vector field
	///	</parameters>
	static void VecTransABPFromRLL(
		double dX,
		double dY,
		int nP,
		double dUlon,
		double dUlat,
		double &dUalpha,
		double &dUbeta
	);

	///	<summary>
	///		Translate a vector in equiangular coordinates to a vector in
	///		spherical coordinates.  The components of the vector field in
	///		spherical coordinates are expected in the unit basis, whereas the
	///		components in equiangular coordinates are given in the geometric
	///		basis.
	///	</summary>
	///	<parameters>
	///		dX          - Gnomonic X coordinate of translation point
	///		dY          - Gnomonic Y coordinate of translation point
	///		nP          - Gnomonic panel coordinate (0-5, 4 = north, 5 = south)
	///		              of translation point
	///		dUalpha     - Alpha component of vector field
	///		dUbeta      - Beta component of vector field
	///		dUlon (OUT) - Longitudinal component of vector field (OUT)
	///		dUlat (OUT) - Latitudinal component of vector field (OUT)
	///	</parameters>
	static void VecTransRLLFromABP(
		double dX,
		double dY,
		int nP,
		double dUalpha,
		double dUbeta,
		double &dUlon,
		double &dUlat
	);

	///	<summary>
	///		Compute the angle between grid lines on a gnomonic or
	///		equiangular projection.
	///	</summary>
	///	<parameters>
	///		dX - Gnomonic X coordinate
	///		dY - Gnomonic Y coordinate
	///	</parameters>
	inline static double GnomonicGridAngle(
		double dX,
		double dY
	) {
		return acos(- dX * dY / (sqrt(1.0 + dX * dX) * sqrt(1.0 + dY * dY)));
	}

	///	<summary>
	///		Compute the area of a square region on the Gnomonic cubed sphere
	///		grid.
	///	</summary>
	///	<parameters>
	///		dX      - Gnomonic X coordinate of bottom-left corner
	///		dDeltaX - Delta X
	///		dY      - Gnomonic Y coordinate of the bottom-left corner
	///		dDeltaY - Delta Y
	///		dRadius - Radius of the sphere
	///	</parameters>
	inline static double GnomonicElementArea(
		double dX,
		double dDeltaX,
		double dY,
		double dDeltaY,
		double dRadius = 1.0
	) {
		return dRadius * dRadius *
		       (+ GnomonicGridAngle(dX          , dY   )
	    	    - GnomonicGridAngle(dX + dDeltaX, dY   )
				- GnomonicGridAngle(dX          , dY + dDeltaY)
				+ GnomonicGridAngle(dX + dDeltaX, dY + dDeltaY));
	}

	///	<summary>
	///		Determine the panel id in the given direction relative to
	///		a given panel.
	///	</summary>
	///	<parameters>
	///		iPanel            - Source panel
	///		dir               - Direction relative to source panel
	///		nRelPanel (OUT)   - Panel in the direction "dir"
	///		dirOpposing (OUT) - Direction on the opposing panel that returns
	///		                    to the source panel
	///	</parameters>
	inline static void RelativePanel(
		int iPanel,
		Direction dir,
		int & iRelPanel,
		Direction & dirOpposing
	) {

		// Equatorial panel 0
		if (iPanel == 0) {
			if (dir == Direction_Right) {
				iRelPanel = 1;
				dirOpposing = Direction_Left;

			} else if (dir == Direction_Top) {
				iRelPanel = 4;
				dirOpposing = Direction_Bottom;

			} else if (dir == Direction_Left) {
				iRelPanel = 3;
				dirOpposing = Direction_Right;

			} else if (dir == Direction_Bottom) {
				iRelPanel = 5;
				dirOpposing = Direction_Top;

			} else {
				_EXCEPTIONT("Invalid direction.");
			}

		// Equatorial panel 1
		} else if (iPanel == 1) {
			if (dir == Direction_Right) {
				iRelPanel = 2;
				dirOpposing = Direction_Left;

			} else if (dir == Direction_Top) {
				iRelPanel = 4;
				dirOpposing = Direction_Right;

			} else if (dir == Direction_Left) {
				iRelPanel = 0;
				dirOpposing = Direction_Right;

			} else if (dir == Direction_Bottom) {
				iRelPanel = 5;
				dirOpposing = Direction_Right;

			} else {
				_EXCEPTIONT("Invalid direction.");
			}

		// Equatorial panel 2
		} else if (iPanel == 2) {
			if (dir == Direction_Right) {
				iRelPanel = 3;
				dirOpposing = Direction_Left;

			} else if (dir == Direction_Top) {
				iRelPanel = 4;
				dirOpposing = Direction_Top;

			} else if (dir == Direction_Left) {
				iRelPanel = 1;
				dirOpposing = Direction_Right;

			} else if (dir == Direction_Bottom) {
				iRelPanel = 5;
				dirOpposing = Direction_Bottom;

			} else {
				_EXCEPTIONT("Invalid direction.");
			}

		// Equatorial panel 3
		} else if (iPanel == 3) {
			if (dir == Direction_Right) {
				iRelPanel = 0;
				dirOpposing = Direction_Left;

			} else if (dir == Direction_Top) {
				iRelPanel = 4;
				dirOpposing = Direction_Left;

			} else if (dir == Direction_Left) {
				iRelPanel = 2;
				dirOpposing = Direction_Right;

			} else if (dir == Direction_Bottom) {
				iRelPanel = 5;
				dirOpposing = Direction_Left;

			} else {
				_EXCEPTIONT("Invalid direction.");
			}

		// North polar panel
		} else if (iPanel == 4) {
			if (dir == Direction_Right) {
				iRelPanel = 1;
				dirOpposing = Direction_Top;

			} else if (dir == Direction_Top) {
				iRelPanel = 2;
				dirOpposing = Direction_Top;

			} else if (dir == Direction_Left) {
				iRelPanel = 3;
				dirOpposing = Direction_Top;

			} else if (dir == Direction_Bottom) {
				iRelPanel = 0;
				dirOpposing = Direction_Top;

			} else {
				_EXCEPTIONT("Invalid direction.");
			}

		// South polar panel
		} else if (iPanel == 5) {
			if (dir == Direction_Right) {
				iRelPanel = 1;
				dirOpposing = Direction_Bottom;

			} else if (dir == Direction_Top) {
				iRelPanel = 0;
				dirOpposing = Direction_Bottom;

			} else if (dir == Direction_Left) {
				iRelPanel = 3;
				dirOpposing = Direction_Bottom;

			} else if (dir == Direction_Bottom) {
				iRelPanel = 2;
				dirOpposing = Direction_Bottom;

			} else {
				_EXCEPTIONT("Invalid direction.");
			}

		} else {
			_EXCEPTIONT("Invalid panel.");
		}
	}

	///	<summary>
	///		Determine the panel id in the given direction relative to
	///		a given panel id.
	///	</summary>
	///	<parameters>
	///		Nc               - Resolution of the cubed sphere grid
	///		p_src            - Source panel
	///		ix_src           - Index in X direction of source element
	///		jx_src           - Index in Y direction of source element
	///		p_dest (OUT)     - Destination panel
	///		ix_dest (OUT)    - Index in the X direction of destination element
	///		jx_dest (OUT)    - Index in the Y direction of destination element
	///		switchAB (OUT)   - Flag indicating whether alpha and beta switch
	///		                   across the panel.
	///		switchPar (OUT)  - Flag indicating whether the parallel
	///		                   direction switches sign across the panel.
	///		switchPerp (OUT) - Flag indicating whether the perpendicular
	///		                   direction switches sign across the panel.
	///	</parameters>
	static void RelativeCoord(
		int Nc,
		int p_src,
		int ix_src,
		int jx_src,
		int & p_dest,
		int & ix_dest,
		int & jx_dest,
		bool & switchAB,
		bool & switchPar,
		bool & switchPerp
	) {
		// Internal points
		if ((jx_src >= 0) && (ix_src >= 0) && (jx_src < Nc) && (ix_src < Nc)) {
			jx_dest = jx_src;
			ix_dest = ix_src;
			p_dest = p_src;
			switchAB = false;
			switchPar = false;
			switchPerp = false;

		// Exterior points
		} else if (
			((ix_src < 0) || (ix_src >= Nc)) &&
			((jx_src < 0) || (jx_src >= Nc))
		) {
			jx_dest = (-1);
			ix_dest = (-1);
			p_dest = InvalidPanel;
			switchAB = false;
			switchPar = false;
			switchPerp = false;

		// Equatorial panel 0
		} else if (p_src == 0) {
			if (ix_src >= Nc) {
				p_dest = 1;
				jx_dest = jx_src;
				ix_dest = ix_src - Nc;
				switchAB = false;
				switchPar = false;
				switchPerp = false;

			} else if (jx_src >= Nc) {
				p_dest = 4;
				jx_dest = jx_src - Nc;
				ix_dest = ix_src;
				switchAB = false;
				switchPar = false;
				switchPerp = false;

			} else if (ix_src < 0) {
				p_dest = 3;
				jx_dest = jx_src;
				ix_dest = ix_src + Nc;
				switchAB = false;
				switchPar = false;
				switchPerp = false;

			} else if (jx_src < 0) {
				p_dest = 5;
				jx_dest = jx_src + Nc;
				ix_dest = ix_src;
				switchAB = false;
				switchPar = false;
				switchPerp = false;
			}

		// Equatorial panel 1
		} else if (p_src == 1) {
			if (ix_src >= Nc) {
				p_dest = 2;
				jx_dest = jx_src;
				ix_dest = ix_src - Nc;
				switchAB = false;
				switchPar = false;
				switchPerp = false;

			} else if (jx_src >= Nc) {
				p_dest = 4;
				jx_dest = ix_src;
				ix_dest = 2 * Nc - 1 - jx_src;
				switchAB = true;
				switchPar = false;
				switchPerp = true;

			} else if (ix_src < 0) {
				p_dest = 0;
				jx_dest = jx_src;
				ix_dest = ix_src + Nc;
				switchAB = false;
				switchPar = false;
				switchPerp = false;

			} else if (jx_src < 0) {
				p_dest = 5;
				jx_dest = Nc - 1 - ix_src;
				ix_dest = Nc + jx_src;
				switchAB = true;
				switchPar = true;
				switchPerp = false;
			}

		// Equatorial panel 2
		} else if (p_src == 2) {
			if (ix_src >= Nc) {
				p_dest = 3;
				jx_dest = jx_src;
				ix_dest = ix_src - Nc;
				switchAB = false;
				switchPar = false;
				switchPerp = false;

			} else if (jx_src >= Nc) {
				p_dest = 4;
				jx_dest = 2 * Nc - 1 - jx_src;
				ix_dest = Nc - ix_src - 1;
				switchAB = false;
				switchPar = true;
				switchPerp = true;

			} else if (ix_src < 0) {
				p_dest = 1;
				jx_dest = jx_src;
				ix_dest = ix_src + Nc;
				switchAB = false;
				switchPar = false;
				switchPerp = false;

			} else if (jx_src < 0) {
				p_dest = 5;
				jx_dest = - jx_src - 1;
				ix_dest = Nc - ix_src - 1;
				switchAB = false;
				switchPar = true;
				switchPerp = true;
			}

		// Equatorial panel 3
		} else if (p_src == 3) {
			if (ix_src >= Nc) {
				p_dest = 0;
				jx_dest = jx_src;
				ix_dest = ix_src - Nc;
				switchAB = false;
				switchPar = false;
				switchPerp = false;

			} else if (jx_src >= Nc) {
				p_dest = 4;
				jx_dest = Nc - ix_src - 1;
				ix_dest = jx_src - Nc;
				switchAB = true;
				switchPar = true;
				switchPerp = false;

			} else if (ix_src < 0) {
				p_dest = 2;
				jx_dest = jx_src;
				ix_dest = ix_src + Nc;
				switchAB = false;
				switchPar = false;
				switchPerp = false;

			} else if (jx_src < 0) {
				p_dest = 5;
				jx_dest = ix_src;
				ix_dest = - jx_src - 1;
				switchAB = true;
				switchPar = false;
				switchPerp = true;
			}

		// North polar panel
		} else if (p_src == 4) {
			if (ix_src >= Nc) {
				p_dest = 1;
				jx_dest = 2 * Nc - 1 - ix_src;
				ix_dest = jx_src;
				switchAB = true;
				switchPar = false;
				switchPerp = true;

			} else if (jx_src >= Nc) {
				p_dest = 2;
				jx_dest = 2 * Nc - 1 - jx_src;
				ix_dest = Nc - 1 - ix_src;
				switchAB = false;
				switchPar = true;
				switchPerp = true;

			} else if (ix_src < 0) {
				p_dest = 3;
				jx_dest = ix_src + Nc;
				ix_dest = Nc - 1 - jx_src;
				switchAB = true;
				switchPar = true;
				switchPerp = false;

			} else if (jx_src < 0) {
				p_dest = 0;
				jx_dest = jx_src + Nc;
				ix_dest = ix_src;
				switchAB = false;
				switchPar = false;
				switchPerp = false;
			}

		// South polar panel
		} else if (p_src == 5) {
			if (ix_src >= Nc) {
				p_dest = 1;
				jx_dest = ix_src - Nc;
				ix_dest = Nc - 1 - jx_src;
				switchAB = true;
				switchPar = true;
				switchPerp = false;

			} else if (jx_src >= Nc) {
				p_dest = 0;
				jx_dest = jx_src - Nc;
				ix_dest = ix_src;
				switchAB = false;
				switchPar = false;
				switchPerp = false;

			} else if (ix_src < 0) {
				p_dest = 3;
				jx_dest = - ix_src - 1;
				ix_dest = jx_src;
				switchAB = true;
				switchPar = false;
				switchPerp = true;

			} else if (jx_src < 0) {
				p_dest = 2;
				jx_dest = - jx_src - 1;
				ix_dest = Nc - 1 - ix_src;
				switchAB = false;
				switchPar = true;
				switchPerp = true;
			}
		}
	}

	///	<summary>
	///		Determine the panel in the given direction of the source panel.
	///	</summary>
	static inline int PanelFromDirection(
		int iPanel,
		Direction dir
	) {
		// Equatorial panels
		if (iPanel < 4) {
			if (dir == Direction_Middle) {
				return iPanel;
			} else if (dir == Direction_Right) {
				return (iPanel + 1) % 4;
			} else if (dir == Direction_Top) {
				return 4;
			} else if (dir == Direction_Left) {
				return (iPanel + 3) % 4;
			} else if (dir == Direction_Bottom) {
				return 5;
			} else {
				_EXCEPTIONT("Invalid direction.");
			}

		// North polar panel
		} else if (iPanel == 4) {
			if (dir == Direction_Middle) {
				return 4;
			} else if (dir == Direction_Right) {
				return 1;
			} else if (dir == Direction_Top) {
				return 2;
			} else if (dir == Direction_Left) {
				return 3;
			} else if (dir == Direction_Bottom) {
				return 0;
			} else {
				_EXCEPTIONT("Invalid direction.");
			}

		// South polar panel
		} else if (iPanel == 5) {
			if (dir == Direction_Middle) {
				return 5;
			} else if (dir == Direction_Right) {
				return 1;
			} else if (dir == Direction_Top) {
				return 0;
			} else if (dir == Direction_Left) {
				return 3;
			} else if (dir == Direction_Bottom) {
				return 2;
			} else {
				_EXCEPTIONT("Invalid direction.");
			}

		} else {
			_EXCEPTIONT("Invalid panel.");
		}
	}

	///	<summary>
	///		Determine the direction between two panels.
	///	</summary>
	///	<parameters>
	///		iPanelOrigin   - Origin panel
	///		iPanelDest     - Destination panel
	///	</parameters>
	static inline Direction DirectionFromPanel(
		int iPanelOrigin,
		int iPanelDest
	) {
		if (iPanelDest > 5) {
			_EXCEPTIONT("Invalid iPanelDest");
		}

		if (iPanelOrigin == 0) {
			if (iPanelDest == 0) {
				return Direction_Middle;
			} else if (iPanelDest == 1) {
				return Direction_Right;
			} else if (iPanelDest == 2) {
				return Direction_Unreachable;
			} else if (iPanelDest == 3) {
				return Direction_Left;
			} else if (iPanelDest == 4) {
				return Direction_Top;
			} else if (iPanelDest == 5) {
				return Direction_Bottom;
			}

		} else if (iPanelOrigin == 1) {
			if (iPanelDest == 0) {
				return Direction_Left;
			} else if (iPanelDest == 1) {
				return Direction_Middle;
			} else if (iPanelDest == 2) {
				return Direction_Right;
			} else if (iPanelDest == 3) {
				return Direction_Unreachable;
			} else if (iPanelDest == 4) {
				return Direction_Top;
			} else if (iPanelDest == 5) {
				return Direction_Bottom;
			}

		} else if (iPanelOrigin == 2) {
			if (iPanelDest == 0) {
				return Direction_Unreachable;
			} else if (iPanelDest == 1) {
				return Direction_Left;
			} else if (iPanelDest == 2) {
				return Direction_Middle;
			} else if (iPanelDest == 3) {
				return Direction_Right;
			} else if (iPanelDest == 4) {
				return Direction_Top;
			} else if (iPanelDest == 5) {
				return Direction_Bottom;
			}

		} else if (iPanelOrigin == 3) {
			if (iPanelDest == 0) {
				return Direction_Right;
			} else if (iPanelDest == 1) {
				return Direction_Unreachable;
			} else if (iPanelDest == 2) {
				return Direction_Left;
			} else if (iPanelDest == 3) {
				return Direction_Middle;
			} else if (iPanelDest == 4) {
				return Direction_Top;
			} else if (iPanelDest == 5) {
				return Direction_Bottom;
			}

		} else if (iPanelOrigin == 4) {
			if (iPanelDest == 0) {
				return Direction_Bottom;
			} else if (iPanelDest == 1) {
				return Direction_Right;
			} else if (iPanelDest == 2) {
				return Direction_Top;
			} else if (iPanelDest == 3) {
				return Direction_Left;
			} else if (iPanelDest == 4) {
				return Direction_Middle;
			} else if (iPanelDest == 5) {
				return Direction_Unreachable;
			}

		} else if (iPanelOrigin == 5) {
			if (iPanelDest == 0) {
				return Direction_Top;
			} else if (iPanelDest == 1) {
				return Direction_Right;
			} else if (iPanelDest == 2) {
				return Direction_Bottom;
			} else if (iPanelDest == 3) {
				return Direction_Left;
			} else if (iPanelDest == 4) {
				return Direction_Unreachable;
			} else if (iPanelDest == 5) {
				return Direction_Middle;
			}
		}

		_EXCEPTIONT("Invalid iPanelOrigin");
	}

	///	<summary>
	///		Determine the opposite direction when crossing a panel boundary.
	///	</summary>
	static inline Direction OpposingDirection(
		int iPanelSrc,
		int iPanelDest,
		Direction dir
	) {
		if (dir == Direction_Unreachable) {
			return Direction_Unreachable;
		}
		if (dir == Direction_Middle) {
			if (iPanelSrc != iPanelDest) {
				_EXCEPTIONT("Inconsistent parameters");
			}
			return Direction_Middle;
		}
		if (iPanelSrc == iPanelDest) {
			if (dir == Direction_Right) {
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
			}
		}
		if (iPanelSrc == 0) {
			if (dir == Direction_Right) {
				return Direction_Left;
			} else if (dir == Direction_Top) {
				return Direction_Bottom;
			} else if (dir == Direction_Left) {
				return Direction_Right;
			} else if (dir == Direction_Bottom) {
				return Direction_Top;

			} else if (dir == Direction_TopRight) {
				if ((iPanelDest != 1) && (iPanelDest != 4)) {
					_EXCEPTIONT("Inconsistent parameters");
				}
				return Direction_BottomLeft;

			} else if (dir == Direction_TopLeft) {
				if ((iPanelDest != 3) && (iPanelDest != 4)) {
					_EXCEPTIONT("Inconsistent parameters");
				}
				return Direction_BottomRight;

			} else if (dir == Direction_BottomLeft) {
				if ((iPanelDest != 3) && (iPanelDest != 5)) {
					_EXCEPTIONT("Inconsistent parameters");
				}
				return Direction_TopRight;

			} else if (dir == Direction_BottomRight) {
				if ((iPanelDest != 1) && (iPanelDest != 5)) {
					_EXCEPTIONT("Inconsistent parameters");
				}
				return Direction_TopLeft;
			}
		}
		if (iPanelSrc == 1) {
			if (dir == Direction_Right) {
				return Direction_Left;
			} else if (dir == Direction_Top) {
				return Direction_Right;
			} else if (dir == Direction_Left) {
				return Direction_Right;
			} else if (dir == Direction_Bottom) {
				return Direction_Right;

			} else if (dir == Direction_TopRight) {
				if (iPanelDest == 4) {
					return Direction_BottomRight;
				} else if (iPanelDest == 2) {
					return Direction_BottomLeft;
				} else {
					_EXCEPTIONT("Inconsistent parameters");
				}

			} else if (dir == Direction_TopLeft) {
				if (iPanelDest == 4) {
					return Direction_TopRight;
				} else if (iPanelDest == 0) {
					return Direction_BottomRight;
				} else {
					_EXCEPTIONT("Inconsistent parameters");
				}

			} else if (dir == Direction_BottomLeft) {
				if (iPanelDest == 5) {
					return Direction_BottomRight;
				} else if (iPanelDest == 0) {
					return Direction_TopRight;
				} else {
					_EXCEPTIONT("Inconsistent parameters");
				}

			} else if (dir == Direction_BottomRight) {
				if (iPanelDest == 5) {
					return Direction_TopRight;
				} else if (iPanelDest == 2) {
					return Direction_TopLeft;
				} else {
					_EXCEPTIONT("Inconsistent parameters");
				}
			}
		}
		if (iPanelSrc == 2) {
			if (dir == Direction_Right) {
				return Direction_Left;
			} else if (dir == Direction_Top) {
				return Direction_Top;
			} else if (dir == Direction_Left) {
				return Direction_Right;
			} else if (dir == Direction_Bottom) {
				return Direction_Bottom;

			} else if (dir == Direction_TopRight) {
				if (iPanelDest == 4) {
					return Direction_TopRight;
				} else if (iPanelDest == 3) {
					return Direction_BottomLeft;
				} else {
					_EXCEPTIONT("Inconsistent parameters");
				}

			} else if (dir == Direction_TopLeft) {
				if (iPanelDest == 4) {
					return Direction_TopLeft;
				} else if (iPanelDest == 1) {
					return Direction_BottomRight;
				} else {
					_EXCEPTIONT("Inconsistent parameters");
				}

			} else if (dir == Direction_BottomLeft) {
				if (iPanelDest == 5) {
					return Direction_BottomLeft;
				} else if (iPanelDest == 1) {
					return Direction_TopRight;
				} else {
					_EXCEPTIONT("Inconsistent parameters");
				}

			} else if (dir == Direction_BottomRight) {
				if (iPanelDest == 5) {
					return Direction_BottomRight;
				} else if (iPanelDest == 3) {
					return Direction_TopLeft;
				} else {
					_EXCEPTIONT("Inconsistent parameters");
				}
			}
		}
		if (iPanelSrc == 3) {
			if (dir == Direction_Right) {
				return Direction_Left;
			} else if (dir == Direction_Top) {
				return Direction_Left;
			} else if (dir == Direction_Left) {
				return Direction_Right;
			} else if (dir == Direction_Bottom) {
				return Direction_Left;

			} else if (dir == Direction_TopRight) {
				if (iPanelDest == 4) {
					return Direction_TopLeft;
				} else if (iPanelDest == 0) {
					return Direction_BottomLeft;
				} else {
					_EXCEPTIONT("Inconsistent parameters");
				}

			} else if (dir == Direction_TopLeft) {
				if (iPanelDest == 4) {
					return Direction_BottomLeft;
				} else if (iPanelDest == 2) {
					return Direction_BottomRight;
				} else {
					_EXCEPTIONT("Inconsistent parameters");
				}

			} else if (dir == Direction_BottomLeft) {
				if (iPanelDest == 5) {
					return Direction_TopLeft;
				} else if (iPanelDest == 2) {
					return Direction_TopRight;
				} else {
					_EXCEPTIONT("Inconsistent parameters");
				}

			} else if (dir == Direction_BottomRight) {
				if (iPanelDest == 5) {
					return Direction_BottomLeft;
				} else if (iPanelDest == 0) {
					return Direction_TopLeft;
				} else {
					_EXCEPTIONT("Inconsistent parameters");
				}
			}
		}
		if (iPanelSrc == 4) {
			if (dir == Direction_Right) {
				return Direction_Top;
			} else if (dir == Direction_Top) {
				return Direction_Top;
			} else if (dir == Direction_Left) {
				return Direction_Top;
			} else if (dir == Direction_Bottom) {
				return Direction_Top;

			} else if (dir == Direction_TopRight) {
				if (iPanelDest == 1) {
					return Direction_TopLeft;
				} else if (iPanelDest == 2) {
					return Direction_TopRight;
				} else {
					_EXCEPTIONT("Inconsistent parameters");
				}

			} else if (dir == Direction_TopLeft) {
				if (iPanelDest == 2) {
					return Direction_TopLeft;
				} else if (iPanelDest == 3) {
					return Direction_TopRight;
				} else {
					_EXCEPTIONT("Inconsistent parameters");
				}

			} else if (dir == Direction_BottomLeft) {
				if (iPanelDest == 3) {
					return Direction_TopLeft;
				} else if (iPanelDest == 0) {
					return Direction_TopRight;
				} else {
					_EXCEPTIONT("Inconsistent parameters");
				}

			} else if (dir == Direction_BottomRight) {
				if (iPanelDest == 0) {
					return Direction_TopLeft;
				} else if (iPanelDest == 1) {
					return Direction_TopRight;
				} else {
					_EXCEPTIONT("Inconsistent parameters");
				}
			}
		}
		if (iPanelSrc == 5) {
			if (dir == Direction_Right) {
				return Direction_Bottom;
			} else if (dir == Direction_Top) {
				return Direction_Bottom;
			} else if (dir == Direction_Left) {
				return Direction_Bottom;
			} else if (dir == Direction_Bottom) {
				return Direction_Bottom;

			} else if (dir == Direction_TopRight) {
				if (iPanelDest == 1) {
					return Direction_BottomRight;
				} else if (iPanelDest == 0) {
					return Direction_BottomLeft;
				} else {
					_EXCEPTIONT("Inconsistent parameters");
				}

			} else if (dir == Direction_TopLeft) {
				if (iPanelDest == 0) {
					return Direction_BottomRight;
				} else if (iPanelDest == 3) {
					return Direction_BottomLeft;
				} else {
					_EXCEPTIONT("Inconsistent parameters");
				}

			} else if (dir == Direction_BottomLeft) {
				if (iPanelDest == 3) {
					return Direction_BottomRight;
				} else if (iPanelDest == 2) {
					return Direction_BottomLeft;
				} else {
					_EXCEPTIONT("Inconsistent parameters");
				}

			} else if (dir == Direction_BottomRight) {
				if (iPanelDest == 2) {
					return Direction_BottomRight;
				} else if (iPanelDest == 1) {
					return Direction_BottomLeft;
				} else {
					_EXCEPTIONT("Inconsistent parameters");
				}
			}
		}

		_EXCEPTIONT("Inconsistent parameters");
	}

	///	<summary>
	///		Convert a ABP coordinate to a coordinate on the interior of a panel.
	///	</summary>
	static inline Direction DirectionFromXYP(
		int iPanel,
		double dX,
		double dY
	) {
		static const double TINY = 1e-14;

		Direction dir;

		// Equatorial panels
		if (iPanel < 4) {
			if (dY > 1.0 - TINY) {
				if (dX < -dY - TINY) {
					dir = Direction_Left;
				} else if (dX > dY + TINY) {
					dir = Direction_Right;
				} else {
					dir = Direction_Top;
				}

			} else if (dY < - 1.0 + TINY) {
				if (dX < dY - TINY) {
					dir = Direction_Left;
				} else if (dX > - dY + TINY) {
					dir = Direction_Right;
				} else {
					dir = Direction_Bottom;
				}

			} else if (dX < - 1.0 + TINY) {
				dir = Direction_Left;

			} else if (dX < 1.0 + TINY) {
				dir = Direction_Middle;

			} else {
				dir = Direction_Right;
			}

		// North polar panel
		} else if (iPanel == 4) {
			if (dY >= 1.0 + TINY) {
				if (dX > dY - TINY) {
					dir = Direction_Right;
				} else if (dX > - dY - TINY) {
					dir = Direction_Top;
				} else {
					dir = Direction_Left;
				}

			} else if (dY <= - 1.0 - TINY) {
				if (dX < dY + TINY) {
					dir = Direction_Left;
				} else if (dX < - dY + TINY) {
					dir = Direction_Bottom;
				} else {
					dir = Direction_Right;
				}

			} else if (dX <= - 1.0 - TINY) {
				dir = Direction_Left;

			} else if (dX >= 1.0 + TINY) {
				dir = Direction_Right;

			} else {
				dir = Direction_Middle;
			}

		// South polar panel
		} else {
			if (dY >= 1.0 + TINY) {
				if (dX < - dY + TINY) {
					dir = Direction_Left;
				} else if (dX < dY + TINY) {
					dir = Direction_Top;
				} else {
					dir = Direction_Right;
				}

			} else if (dY <= - 1.0 - TINY) {
				if (dX > - dY - TINY) {
					dir = Direction_Right;
				} else if (dX > dY - TINY) {
					dir = Direction_Bottom;
				} else {
					dir = Direction_Left;
				}

			} else if (dX <= - 1.0 - TINY) {
				dir = Direction_Left;

			} else if (dX >= 1.0 + TINY) {
				dir = Direction_Right;

			} else {
				dir = Direction_Middle;
			}
		}

		return dir;
	}

	///	<summary>
	///		Remap vectors from the given source panel to the given destination
	///		panel.  All calculations are performed in-place.
	///	</summary>
	///	<param name="p_src">
	///		Panel index of the source panel.
	///	</param>
	///	<param name="p_dest">
	///		Panel index of the destination panel.
	///	</param>
	///	<param name="dAlpha">
	///		Alpha component of the vector.
	///	</param>
	///	<param name="dBeta">
	///		Beta component of the vector.
	///	</param>
	///	<param name="dXdest">
	///		Gnomonic X coordinate (panel 2) of destination point.
	///	</param>
	///	<param name="dYdest">
	///		Gnomonic Y coordinate (panel 2) of destination point.
	///	</param>
	static void VecPanelTrans(
		int p_src,
		int p_dest,
		double &dAlpha,
		double &dBeta,
		double dXdest,
		double dYdest
	) {
		if ((p_dest < 4) && (p_src < 4)) {
			if ((p_dest + 7 - p_src) % 4 == 0) {
				VecToP2FromP1(dAlpha, dBeta, dXdest, dYdest);
			} else if ((p_dest + 5 - p_src) % 4 == 0) {
				VecToP1FromP2(dAlpha, dBeta, dXdest, dYdest);
			} else {
				_EXCEPTION();
			}

		} else if (p_src == 4) {
			if (p_dest == 0) {
				VecToP1FromP5(dAlpha, dBeta, dXdest, dYdest);
			} else if (p_dest == 1) {
				VecToP2FromP5(dAlpha, dBeta, dXdest, dYdest);
			} else if (p_dest == 2) {
				VecToP3FromP5(dAlpha, dBeta, dXdest, dYdest);
			} else if (p_dest == 3) {
				VecToP4FromP5(dAlpha, dBeta, dXdest, dYdest);
			} else {
				_EXCEPTION();
			}

		} else if (p_dest == 4) {
			if (p_src == 0) {
				VecToP5FromP1(dAlpha, dBeta, dXdest, dYdest);
			} else if (p_src == 1) {
				VecToP5FromP2(dAlpha, dBeta, dXdest, dYdest);
			} else if (p_src == 2) {
				VecToP5FromP3(dAlpha, dBeta, dXdest, dYdest);
			} else if (p_src == 3) {
				VecToP5FromP4(dAlpha, dBeta, dXdest, dYdest);
			} else {
				_EXCEPTION();
			}

		} else if (p_src == 5) {
			if (p_dest == 0) {
				VecToP1FromP6(dAlpha, dBeta, dXdest, dYdest);
			} else if (p_dest == 1) {
				VecToP2FromP6(dAlpha, dBeta, dXdest, dYdest);
			} else if (p_dest == 2) {
				VecToP3FromP6(dAlpha, dBeta, dXdest, dYdest);
			} else if (p_dest == 3) {
				VecToP4FromP6(dAlpha, dBeta, dXdest, dYdest);
			} else {
				_EXCEPTION();
			}

		} else if (p_dest == 5) {
			if (p_src == 0) {
				VecToP6FromP1(dAlpha, dBeta, dXdest, dYdest);
			} else if (p_src == 1) {
				VecToP6FromP2(dAlpha, dBeta, dXdest, dYdest);
			} else if (p_src == 2) {
				VecToP6FromP3(dAlpha, dBeta, dXdest, dYdest);
			} else if (p_src == 3) {
				VecToP6FromP4(dAlpha, dBeta, dXdest, dYdest);
			} else {
				_EXCEPTION();
			}

		} else {
			_EXCEPTION();
		}
	}

	///	<summary>
	///		Remap vectors to equatorial panel 2 from equatorial panel 1.  All
	///		calculations are performed in-place.
	///	</summary>
	///	<param name="dAlpha">
	///		Alpha component of the vector system.
	///	</param>
	///	<param name="dBeta">
	///		Beta component of the vector system.
	///	</param>
	///	<param name="dX2">
	///		Gnomonic X coordinate (panel 2) of point of rotation.
	///	</param>
	///	<param name="dY2">
	///		Gnomonic Y coordinate (panel 2) of point of rotation.
	///	</param>
	inline static void VecToP2FromP1(
		double &dAlpha,
		double &dBeta,
		double dX2,
		double dY2
	) {
		dBeta =
			dY2 * (1.0 + dX2 * dX2) / (dX2 * (1.0 + dY2 * dY2)) * dAlpha
			- (dX2 * dX2 + dY2 * dY2) / (dX2 * (1.0 + dY2 * dY2)) * dBeta;
	}

	///	<summary>
	///		As VecToP2FromP1 except remapping vectors to panel 1 from
	///		equatorial panel 2.
	///	</summary>
	inline static void VecToP1FromP2(
		double &dAlpha,
		double &dBeta,
		double dX1,
		double dY1
	) {
		dBeta =
			dY1 * (1.0 + dX1 * dX1) / (dX1 * (1.0 + dY1 * dY1)) * dAlpha
			+ (dX1 * dX1 + dY1 * dY1) / (dX1 * (1.0 + dY1 * dY1)) * dBeta;
	}

	///	<summary>
	///		As VecToP2FromP1 except remapping vectors to panel 5 (north pole)
	///		from equatorial panel 1.
	///	</summary>
	inline static void VecToP5FromP1(
		double &dAlpha,
		double &dBeta,
		double dX5,
		double dY5
	) {
		dAlpha =
			- (dX5 * dX5 + dY5 * dY5) / (dY5 * (dX5 * dX5 + 1.0)) * dAlpha
			+ dX5 * (1.0 + dY5 * dY5) / (dY5 * (dX5 * dX5 + 1.0)) * dBeta;
	}

	///	<summary>
	///		As VecToP2FromP1 except remapping vectors to equatorial panel 1
	///		from panel 5 (north pole).
	///	</summary>
	inline static void VecToP1FromP5(
		double &dAlpha,
		double &dBeta,
		double dX1,
		double dY1
	) {
		dAlpha =
			(dX1 * dX1 + dY1 * dY1) / (dY1 * (1.0 + dX1 * dX1)) * dAlpha
			+ dX1 * (1.0 + dY1 * dY1) / (dY1 * (1.0 + dX1 * dX1)) * dBeta;
	}

	///	<summary>
	///		As VecToP2FromP1 except remapping vectors to panel 5 (north pole)
	///		from equatorial panel 2.
	///	</summary>
	inline static void VecToP5FromP2(
		double &dAlpha,
		double &dBeta,
		double dX5,
		double dY5
	) {
		double dTemp = dBeta;

		dBeta =
			(dX5 * dX5 + dY5 * dY5) / (dX5 * (1.0 + dY5 * dY5)) * dAlpha
			- (dY5 * (dX5 * dX5 + 1.0)) / (dX5 * (1.0 + dY5 * dY5)) * dBeta;

		dAlpha = -dTemp;
	}

	///	<summary>
	///		As VecToP2FromP1 except remapping vectors to equatorial panel 2
	///		from panel 5 (north pole).
	///	</summary>
	inline static void VecToP2FromP5(
		double &dAlpha,
		double &dBeta,
		double dX2,
		double dY2
	) {
		double dTemp = dAlpha;

		dAlpha =
			- (dX2 * (1.0 + dY2 * dY2)) / (dY2 * (1.0 + dX2 * dX2)) * dAlpha
			+ (dX2 * dX2 + dY2 * dY2) / (dY2 * (1.0 + dX2 * dX2)) * dBeta;

		dBeta = -dTemp;
	}

	///	<summary>
	///		As VecToP2FromP1 except remapping vectors to panel 5 (north pole)
	///		from equatorial panel 3.
	///	</summary>
	inline static void VecToP5FromP3(
		double &dAlpha,
		double &dBeta,
		double dX5,
		double dY5
	) {
		dAlpha =
			- (dX5 * dX5 + dY5 * dY5) / (dY5 * (dX5 * dX5 + 1.0)) * dAlpha
			- dX5 * (1.0 + dY5 * dY5) / (dY5 * (dX5 * dX5 + 1.0)) * dBeta;

		dBeta = -dBeta;
	}

	///	<summary>
	///		As VecToP2FromP1 except remapping vectors to equatorial panel 3
	///		from panel 5 (north pole).
	///	</summary>
	inline static void VecToP3FromP5(
		double &dAlpha,
		double &dBeta,
		double dX3,
		double dY3
	) {
		dAlpha =
			- (dX3 * dX3 + dY3 * dY3) / (dY3 * (dX3 * dX3 + 1.0)) * dAlpha
			- dX3 * (dY3 * dY3 + 1.0) / (dY3 * (dX3 * dX3 + 1.0)) * dBeta;

		dBeta = -dBeta;
	}

	///	<summary>
	///		As VecToP2FromP1 except remapping vectors to panel 5 (north pole)
	///		from equatorial panel 4.
	///	</summary>
	inline static void VecToP5FromP4(
		double &dAlpha,
		double &dBeta,
		double dX5,
		double dY5
	) {
		double dTemp = dBeta;

		dBeta =
			(dX5 * dX5 + dY5 * dY5) / (dX5 * (1.0 + dY5 * dY5)) * dAlpha
			+ dY5 * (dX5 * dX5 + 1.0) / (dX5 * (1.0 + dY5 * dY5)) * dBeta;

		dAlpha = dTemp;
	}

	///	<summary>
	///		As VecToP2FromP1 except remapping vectors to equatorial panel 4
	///		from panel 5 (north pole).
	///	</summary>
	inline static void VecToP4FromP5(
		double &dAlpha,
		double &dBeta,
		double dX4,
		double dY4
	) {
		double dTemp = dAlpha;

		dAlpha =
			(dX4 * (1.0 + dY4 * dY4)) / (dY4 * (1.0 + dX4 * dX4)) * dAlpha
			- (dX4 * dX4 + dY4 * dY4) / (dY4 * (1.0 + dX4 * dX4)) * dBeta;

		dBeta = dTemp;
	}

	///	<summary>
	///		As VecToP2FromP1 except remapping vectors to panel 6 (south pole)
	///		from equatorial panel 1.
	///	</summary>
	inline static void VecToP6FromP1(
		double &dAlpha,
		double &dBeta,
		double dX6,
		double dY6
	) {
		dAlpha =
			(dX6 * dX6 + dY6 * dY6) / (dY6 * (1.0 + dX6 * dX6)) * dAlpha
			+ dX6 * (1.0 + dY6 * dY6) / (dY6 * (1.0 + dX6 * dX6)) * dBeta;
	}

	///	<summary>
	///		As VecToP2FromP1 except remapping vectors to equatorial panel 1
	///		from panel 6 (south pole).
	///	</summary>
	inline static void VecToP1FromP6(
		double &dAlpha,
		double &dBeta,
		double dX1,
		double dY1
	) {
		dAlpha =
			- (dX1 * dX1 + dY1 * dY1) / (dY1 * (1.0 + dX1 * dX1)) * dAlpha
			+ dX1 * (1.0 + dY1 * dY1) / (dY1 * (1.0 + dX1 * dX1)) * dBeta;
	}

	///	<summary>
	///		As VecToP2FromP1 except remapping vectors to panel 6 (south pole)
	///		from equatorial panel 2.
	///	</summary>
	inline static void VecToP6FromP2(
		double &dAlpha,
		double &dBeta,
		double dX6,
		double dY6
	) {
		double dTemp = dBeta;

		dBeta =
			- (dX6 * dX6 + dY6 * dY6) / (dX6 * (1.0 + dY6 * dY6)) * dAlpha
			+ (dY6 * (1.0 + dX6 * dX6)) / (dX6 * (1.0 + dY6 * dY6)) * dBeta;

		dAlpha = dTemp;
	}

	///	<summary>
	///		As VecToP2FromP1 except remapping vectors to equatorial panel 2
	///		from panel 6 (south pole).
	///	</summary>
	inline static void VecToP2FromP6(
		double &dAlpha,
		double &dBeta,
		double dX2,
		double dY2
	) {
		double dTemp = dAlpha;
		
		dAlpha =
			(dX2 * (1.0 + dY2 * dY2)) / (dY2 * (1.0 + dX2 * dX2)) * dAlpha
			+ (dX2 * dX2 + dY2 * dY2) / (dY2 * (1.0 + dX2 * dX2)) * dBeta;

		dBeta = dTemp;
	}

	///	<summary>
	///		As VecToP2FromP1 except remapping vectors to panel 6 (south pole)
	///		from equatorial panel 3.
	///	</summary>
	inline static void VecToP6FromP3(
		double &dAlpha,
		double &dBeta,
		double dX6,
		double dY6
	) {
		dAlpha =
			(dX6 * dX6 + dY6 * dY6) / (dY6 * (1.0 + dX6 * dX6)) * dAlpha
			- (dX6 * (1.0 + dY6 * dY6)) / (dY6 * (1.0 + dX6 * dX6)) * dBeta;

		dBeta = -dBeta;
	}

	///	<summary>
	///		As VecToP2FromP1 except remapping vectors to equatorial panel 3
	///		from panel 6 (south pole).
	///	</summary>
	inline static void VecToP3FromP6(
		double &dAlpha,
		double &dBeta,
		double dX3,
		double dY3
	) {
		dAlpha =
			(dX3 * dX3 + dY3 * dY3) / (dY3 * (dX3 * dX3 + 1.0)) * dAlpha
			- dX3 * (dY3 * dY3 + 1.0) / (dY3 * (dX3 * dX3 + 1.0)) * dBeta;

		dBeta = -dBeta;
	}

	///	<summary>
	///		As VecToP2FromP1 except remapping vectors to panel 6 (south pole)
	///		from equatorial panel 4.
	///	</summary>
	inline static void VecToP6FromP4(
		double &dAlpha,
		double &dBeta,
		double dX6,
		double dY6
	) {
		double dTemp = dBeta;

		dBeta =
			- (dX6 * dX6 + dY6 * dY6) / (dX6 * (1.0 + dY6 * dY6)) * dAlpha
			- (dY6 * (1.0 + dX6 * dX6)) / (dX6 * (1.0 + dY6 * dY6)) * dBeta;

		dAlpha = -dTemp;
	}

	///	<summary>
	///		As VecToP2FromP1 except remapping vectors to equatorial panel 4
	///		from panel 6 (south pole).
	///	</summary>
	inline static void VecToP4FromP6(
		double &dAlpha,
		double &dBeta,
		double dX4,
		double dY4
	) {
		double dTemp = dAlpha;

		dAlpha =
			- (dX4 * (1.0 + dY4 * dY4)) / (dY4 * (1.0 + dX4 * dX4)) * dAlpha
			- (dY4 * dY4 + dX4 * dX4) / (dY4 * (1.0 + dX4 * dX4)) * dBeta;

		dBeta = -dTemp;
	}
};

////////////////////////////////////////////////////////////////////////////////

#endif

