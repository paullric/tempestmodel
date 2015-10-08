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

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Computational coordinates of a 2D patch.
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
		int ixBGlobalInteriorEnd
	) :
		m_ixPanel(ixPanel),
		m_iRefinementLevel(iRefinementLevel),
		m_nHaloElements(nHaloElements),
		m_ixABegin(ixAGlobalInteriorBegin - nHaloElements),
		m_ixAEnd(ixAGlobalInteriorEnd + nHaloElements),
		m_ixBBegin(ixBGlobalInteriorBegin - nHaloElements),
		m_ixBEnd(ixBGlobalInteriorEnd + nHaloElements)
	{ }

public:
	///	<summary>
	///		Determine if the coordinate is within this PatchBox.
	///	</summary>
	inline bool ContainsGlobalPoint(
		int ixPanel,
		int ixGlobalA,
		int ixGlobalB
	) const {
		if (ixPanel != m_ixPanel) {
			return false;
		}
		if ((ixGlobalA >= GetAGlobalInteriorBegin()) &&
			(ixGlobalA <  GetAGlobalInteriorEnd()) &&
			(ixGlobalB >= GetBGlobalInteriorBegin()) &&
			(ixGlobalB <  GetBGlobalInteriorEnd())
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
	inline int GetAGlobalBegin() const {
		return (m_ixABegin);
	}

	///	<summary>
	///		Get the global end index of the alpha coordinate on this patch.
	///	</summary>
	inline int GetAGlobalEnd() const {
		return (m_ixAEnd);
	}

	///	<summary>
	///		Get the global begin index of the beta coordinate on this patch.
	///	</summary>
	inline int GetBGlobalBegin() const {
		return (m_ixBBegin);
	}

	///	<summary>
	///		Get the global end index of the beta coordinate on this patch.
	///	</summary>
	inline int GetBGlobalEnd() const {
		return (m_ixBEnd);
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

	///	<summary>
	///		Get the total number of nodes in 2D on this GridPatch.
	///	</summary>
	inline int GetTotalNodeCount2D() const {
		return (GetATotalWidth() * GetBTotalWidth());
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
};

///////////////////////////////////////////////////////////////////////////////

#endif

