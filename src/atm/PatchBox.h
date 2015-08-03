///////////////////////////////////////////////////////////////////////////////
///
///	\file    PatchBox.h
///	\author  Paul Ullrich
///	\version February 25, 2013
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

#ifndef _PATCHBOX_H_
#define _PATCHBOX_H_

#include "Direction.h"
#include "GridSpacing.h"

#include "DataArray1D.h"

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Computational and physical coordinates of a 2D patch.
///	</summary>
class PatchBox {
public:
	///	<summary>
	///		Constructor.
	///	</summary>
	PatchBox(
		int ixPanel,
		int iRefinementLevel,
		int nHaloElements,
		int ixAGlobalInteriorBegin,
		int ixAGlobalInteriorEnd,
		int ixBGlobalInteriorBegin,
		int ixBGlobalInteriorEnd,
		GridSpacing & gridspacingA,
		GridSpacing & gridspacingB
	) :
		m_ixPanel(ixPanel),
		m_iRefinementLevel(iRefinementLevel),
		m_nHaloElements(nHaloElements),
		m_ixABegin(ixAGlobalInteriorBegin - nHaloElements),
		m_ixAEnd(ixAGlobalInteriorEnd + nHaloElements),
		m_ixBBegin(ixBGlobalInteriorBegin - nHaloElements),
		m_ixBEnd(ixBGlobalInteriorEnd + nHaloElements)
	{
		// Initialize the coordinate arrays
		InitializeCoordinates(gridspacingA, gridspacingB);
	}

	///	<summary>
	///		Constructor with coordinate arrays.
	///	</summary>
	PatchBox(
		int ixPanel,
		int iRefinementLevel,
		int nHaloElements,
		int ixAGlobalInteriorBegin,
		int ixAGlobalInteriorEnd,
		int ixBGlobalInteriorBegin,
		int ixBGlobalInteriorEnd,
		const DataArray1D<double> & dANodes,
		const DataArray1D<double> & dBNodes,
		const DataArray1D<double> & dAEdges,
		const DataArray1D<double> & dBEdges
	) :
		m_ixPanel(ixPanel),
		m_iRefinementLevel(iRefinementLevel),
		m_nHaloElements(nHaloElements),
		m_ixABegin(ixAGlobalInteriorBegin - nHaloElements),
		m_ixAEnd(ixAGlobalInteriorEnd + nHaloElements),
		m_ixBBegin(ixBGlobalInteriorBegin - nHaloElements),
		m_ixBEnd(ixBGlobalInteriorEnd + nHaloElements),
		m_dANodes(dANodes),
		m_dBNodes(dBNodes),
		m_dAEdges(dAEdges),
		m_dBEdges(dBEdges)
	{
		// Check coordinate arrays?
	}

public:
	///	<summary>
	///		Initialize the coordinate arrays.
	///	</summary>
	void InitializeCoordinates(
		GridSpacing & gridspacingA,
		GridSpacing & gridspacingB
	);

public:
	///	<summary>
	///		Determine if the coordinate is within this PatchBox.
	///	</summary>
	inline bool ContainsGlobalPoint(
		int ixPanel,
		int ixGlobalA,
		int ixGlboalB
	) const {
		if (ixPanel != m_ixPanel) {
			return false;
		}
		if ((ixGlobalA >= GetAGlobalInteriorBegin()) &&
			(ixGlobalA <  GetAGlobalInteriorEnd()) &&
			(ixGlboalB >= GetBGlobalInteriorBegin()) &&
			(ixGlboalB <  GetBGlobalInteriorEnd())
		) {
			return true;
		}
		return false;
	}

public:
	///	<summary>
	///		Get the panel on which this box appears.
	///	</summary>
	inline int GetPanel() const {
		return m_ixPanel;
	}

	///	<summary>
	///		Get the level of refinement of this box.
	///	</summary>
	inline int GetRefinementLevel() const {
		return m_iRefinementLevel;
	}

	///	<summary>
	///		Get the number of halo elements of this box.
	///	</summary>
	inline int GetHaloElements() const {
		return m_nHaloElements;
	}

	///	<summary>
	///		Get the global begin index of the alpha coordinate on this patch.
	///	</summary>
	inline int GetAGlobalInteriorBegin() const {
		return (m_ixABegin + m_nHaloElements);
	}

	///	<summary>
	///		Get the global end index of the alpha coordinate on this patch.
	///	</summary>
	inline int GetAGlobalInteriorEnd() const {
		return (m_ixAEnd - m_nHaloElements);
	}

	///	<summary>
	///		Get the global begin index of the beta coordinate on this patch.
	///	</summary>
	inline int GetBGlobalInteriorBegin() const {
		return (m_ixBBegin + m_nHaloElements);
	}

	///	<summary>
	///		Get the global end index of the beta coordinate on this patch.
	///	</summary>
	inline int GetBGlobalInteriorEnd() const {
		return (m_ixBEnd - m_nHaloElements);
	}

	///	<summary>
	///		Get the total width of the patch box.
	///	</summary>
	inline int GetATotalWidth() const {
		return (m_ixAEnd - m_ixABegin);
	}

	///	<summary>
	///		Get the total width of the patch box.
	///	</summary>
	inline int GetBTotalWidth() const {
		return (m_ixBEnd - m_ixBBegin);
	}

	///	<summary>
	///		Get the total number of nodes in this patch.
	///	</summary>
	inline int GetTotalNodes() const {
		return (GetATotalWidth() * GetBTotalWidth());
	}

	///	<summary>
	///		Get the length of the interior perimeter of this patch.
	///	</summary>
	inline int GetInteriorPerimeter() const {
		return (2 * GetAInteriorWidth() + 2 * GetBInteriorWidth());
	}

public:
	///	<summary>
	///		Get the local interior begin index in the alpha direction.
	///	</summary>
	inline int GetAInteriorBegin() const {
		return m_nHaloElements;
	}

	///	<summary>
	///		Get the local interior end index in the alpha direction.
	///	</summary>
	inline int GetAInteriorEnd() const {
		return (m_ixAEnd - m_ixABegin - m_nHaloElements);
	}

	///	<summary>
	///		Get the local interior begin index in the beta direction.
	///	</summary>
	inline int GetBInteriorBegin() const {
		return m_nHaloElements;
	}

	///	<summary>
	///		Get the local interior end index in the beta direction.
	///	</summary>
	inline int GetBInteriorEnd() const {
		return (m_ixBEnd - m_ixBBegin - m_nHaloElements);
	}

	///	<summary>
	///		Get the interior width in the alpha direction.
	///	</summary>
	inline int GetAInteriorWidth() const {
		return (m_ixAEnd - m_ixABegin - 2 * m_nHaloElements);
	}

	///	<summary>
	///		Get the interior width in the beta direction.
	///	</summary>
	inline int GetBInteriorWidth() const {
		return (m_ixBEnd - m_ixBBegin - 2 * m_nHaloElements);
	}

public:
	///	<summary>
	///		Get the alpha node with specified local index.
	///	</summary>
	inline double GetANode(int ix) const {
		return m_dANodes[ix];
	}

	///	<summary>
	///		Get the alpha edge with specified local alpha edge index.
	///	</summary>
	inline double GetAEdge(int ix) const {
		return m_dAEdges[ix];
	}

	///	<summary>
	///		Get the array of alpha nodes.
	///	</summary>
	inline const DataArray1D<double> & GetANodes() const {
		return m_dANodes;
	}

	///	<summary>
	///		Get the array of alpha edges.
	///	</summary>
	inline const DataArray1D<double> & GetAEdges() const {
		return m_dAEdges;
	}

	///	<summary>
	///		Get the beta node with specified local index.
	///	</summary>
	inline double GetBNode(int ix) const {
		return m_dBNodes[ix];
	}

	///	<summary>
	///		Get the beta edge with specified local beta edge index.
	///	</summary>
	inline double GetBEdge(int ix) const {
		return m_dBEdges[ix];
	}

	///	<summary>
	///		Get the array of beta nodes.
	///	</summary>
	inline const DataArray1D<double> & GetBNodes() const {
		return m_dBNodes;
	}

	///	<summary>
	///		Get the array of beta edges.
	///	</summary>
	inline const DataArray1D<double> & GetBEdges() const {
		return m_dBEdges;
	}

private:
	///	<summary>
	///		Panel on which this box appears.
	///	</summary>
	int m_ixPanel;

	///	<summary>
	///		Level of refinement where this box is defined.
	///	</summary>
	int m_iRefinementLevel;

	///	<summary>
	///		Number of halo elements on this patch.
	///	</summary>
	int m_nHaloElements;

	///	<summary>
	///		Global indices of patch on panel.
	///	</summary>
	int m_ixABegin;
	int m_ixAEnd;

	int m_ixBBegin;
	int m_ixBEnd;

	///	<summary>
	///		Array of nodal values in the alpha direction.
	///	</summary>
	DataArray1D<double> m_dANodes;

	///	<summary>
	///		Array of edge values in the alpha direction.
	///	</summary>
	DataArray1D<double> m_dAEdges;

	///	<summary>
	///		Array of nodal values in the beta direction.
	///	</summary>
	DataArray1D<double> m_dBNodes;

	///	<summary>
	///		Array of edge values in the beta direction.
	///	</summary>
	DataArray1D<double> m_dBEdges;
};

///////////////////////////////////////////////////////////////////////////////

#endif

