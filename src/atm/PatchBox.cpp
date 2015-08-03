///////////////////////////////////////////////////////////////////////////////
///
///	\file    GridPatch.h
///	\author  Paul Ullrich
///	\version March 8, 2013
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

#include "PatchBox.h"

///////////////////////////////////////////////////////////////////////////////

void PatchBox::InitializeCoordinates(
	GridSpacing & gridspacingA,
	GridSpacing & gridspacingB
) {
	// Coordinates in the alpha direction
	m_dANodes.Allocate(GetATotalWidth());
	m_dAEdges.Allocate(GetATotalWidth()+1);

	for (int i = m_ixABegin; i < m_ixAEnd; i++) {
		m_dANodes[i - m_ixABegin] = gridspacingA.GetNode(i);
	}
	for (int i = m_ixABegin; i <= m_ixAEnd; i++) {
		m_dAEdges[i - m_ixABegin] = gridspacingA.GetEdge(i);
	}

	// Coordinates in the beta direction
	m_dBNodes.Allocate(GetBTotalWidth());
	m_dBEdges.Allocate(GetBTotalWidth()+1);

	for (int j = m_ixBBegin; j < m_ixBEnd; j++) {
		m_dBNodes[j - m_ixBBegin] = gridspacingB.GetNode(j);
	}
	for (int j = m_ixBBegin; j <= m_ixBEnd; j++) {
		m_dBEdges[j - m_ixBBegin] = gridspacingB.GetEdge(j);
	}
}

///////////////////////////////////////////////////////////////////////////////

