///////////////////////////////////////////////////////////////////////////////
///
///	\file    MatlabInput.h
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

#ifndef _MATLABINPUTFILE_H_
#define _MATLABINPUTFILE_H_

#include "DataVector.h"
#include "DataMatrix.h"
#include "DataMatrix3D.h"
#include "DataMatrix4D.h"

#include <map>

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A version 5 MATLAB data file.  This class allows DataVector or
///		DataMatrix objects to be read from MATLAB data files of this type.
///	</summary>
class MatlabInputFile {

private:
	///	<summary>
	///		MATLAB file format specific variables.
	///	</summary>
	static const unsigned int miINT8;
	static const unsigned int miUINT8;
	static const unsigned int miINT16;
	static const unsigned int miUINT16;
	static const unsigned int miINT32;
	static const unsigned int miUINT32;
	static const unsigned int miSINGLE;
	static const unsigned int miDOUBLE;
	static const unsigned int miINT64;
	static const unsigned int miUINT64;
	static const unsigned int miMATRIX;

	static const unsigned int mxDOUBLE_CLASS;

	static const unsigned int mfCOMPLEX;
	static const unsigned int mfGLOBAL;
	static const unsigned int mfLOGICAL;
	static const unsigned int mfNONE;

	static const unsigned int mv8;
	static const unsigned int mv12;

protected:
	///	<summary>
	///		A map between variable names and file positions.
	///	</summary>
	typedef std::map<std::string, long int> NamePositionMap;

	///	<summary>
	///		A pair object for the NamePositionMap.
	///	</summary>
	typedef std::pair<std::string, long int> NamePositionPair;

public:
	///	<summary>
	///		Default constructor.
	///	</summary>
	MatlabInputFile() :
		m_fp(NULL)
	{ }

	///	<summary>
	///		Constructor with specified filename.
	///	</summary>
	///	<param name="szFilename">
	///		A filename for the MATLAB output file.
	///	</param>
	MatlabInputFile(const char *szFilename) :
		m_fp(NULL)
	{
		Open(szFilename);
	}

	///	<summary>
	///		Destructor.
	///	</summary>
	~MatlabInputFile() {
		Close();
	}

public:
	///	<summary>
	///		Open the given MATLAB output file for writing.
	///	</summary>
	///	<param name="szFilename">
	///		A filename for the MATLAB output file.
	///	</param>
	void Open(const char *szFilename, bool fVerbose = false);

	///	<summary>
	///		Close the given MATLAB file.
	///	</summary>
	void Close();

protected:
	///	<summary>
	///		Load in a miMATRIX chunk header from a MATLAB formatted intput file.
	///	</summary>
	void ReadChunkHeader(
		unsigned int * uiArrayFlags,
		unsigned int * nDim,
		char * szArrayName,
		unsigned int nArrayName
	);

public:
	///	<summary>
	///		Determine the number of dimensions of the object with the
	///		given name.
	///	</summary>
	///	<param name="szName">
	///		A variable name used to identify the matrix object in the file.
	///	</param>
	unsigned int GetDataDimension(
		char *szName
	);

	///	<summary>
	///		Determine if the given matrix exists in the file.
	///	</summary>
	///	<param name="szName">
	///		A variable name used to identify the matrix object in the file.
	///	</param>
	bool HasMatrix(
		const char *szName
	);

public:
	///	<summary>
	///		Load in a DataVector from the MATLAB formatted input file.
	///	</summary>
	///	<param name="szName">
	///		A variable name used to identify <c>data</c> in the file.
	///	</param>
	///	<param name="data">
	///		The DataVector object.
	///	</param>
	void InputDataVector(
		const char * szName,
		DataVector<double> &data,
		bool fVerbose = false
	);

	///	<summary>
	///		Load in a DataMatrix from the MATLAB formatted input file.
	///	</summary>
	///	<param name="szName">
	///		A variable name used to identify <c>data</c> in the file.
	///	</param>
	///	<param name="data">
	///		The DataMatrix object.
	///	</param>
	void InputDataMatrix(
		const char * szName,
		DataMatrix<double> &data,
		bool fVerbose = false
	);

	///	<summary>
	///		Load in a DataMatrix3D from the MATLAB formatted input file.
	///	</summary>
	///	<param name="szName">
	///		A variable name used to identify <c>data</c> in the file.
	///	</param>
	///	<param name="data">
	///		The DataMatrix3D object.
	///	</param>
	void InputDataMatrix3D(
		const char * szName,
		DataMatrix3D<double> &data,
		bool fVerbose = false
	);

	///	<summary>
	///		Load in a DataMatrix4D from the MATLAB formatted input file.
	///	</summary>
	///	<param name="szName">
	///		A variable name used to identify <c>data</c> in the file.
	///	</param>
	///	<param name="data">
	///		The DataMatrix4D object.
	///	</param>
	void InputDataMatrix4D(
		const char * szName,
		DataMatrix4D<double> &data,
		bool fVerbose = false
	);

private:
	///	<summary>
	///		A file pointer to a MATLAB file being written.
	///	</summary>
	FILE * m_fp;

	///	<summary>
	///		Endian flag for this MATLAB file. (not used)
	///	</summary>
	bool m_fEndian;

	///	<summary>
	///		A map between variable names and file positions.
	///	</summary>
	NamePositionMap m_mapPos;
};

///////////////////////////////////////////////////////////////////////////////

#endif

