///////////////////////////////////////////////////////////////////////////////
///
///	\file    MatlabOutput.h
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

#ifndef _MATLABOUTPUTFILE_H_
#define _MATLABOUTPUTFILE_H_

#include "DataVector.h"
#include "DataMatrix.h"
#include "DataMatrix3D.h"
#include "DataMatrix4D.h"

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A version 5 MATLAB data file.  This class allows DataVector or
///		DataMatrix objects to be written to MATLAB data files of this type.
///	</summary>
class MatlabOutputFile {

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
	static const unsigned int mv16;

public:
	///	<summary>
	///		Default constructor.
	///	</summary>
	MatlabOutputFile() :
		m_fp(NULL)
	{ }

	///	<summary>
	///		Constructor with specified filename.
	///	</summary>
	///	<param name="szFilename">
	///		A filename for the MATLAB output file.
	///	</param>
	MatlabOutputFile(const char *szFilename) :
		m_fp(NULL)
	{
		Open(szFilename);
	}

	///	<summary>
	///		Destructor.
	///	</summary>
	~MatlabOutputFile() {
		Close();
	}

public:
	///	<summary>
	///		Open the given MATLAB output file for writing.
	///	</summary>
	///	<param name="szFilename">
	///		A filename for the MATLAB output file.
	///	</param>
	void Open(const char *szFilename);

	///	<summary>
	///		Close the given MATLAB file.
	///	</summary>
	void Close();

public:
	///	<summary>
	///		Append a MATLAB-formatted double to a file.
	///	</summary>
	///	<param name="data">
	///		The value to append to the file.
	///	</param>
	///	<param name="szName">
	///		A variable name used to identify <c>data</c> in the file.
	///	</param>
	void OutputDouble(
		const double data,
		const char * szName
	);

	///	<summary>
	///		Append a MATLAB-formatted DataVector to a file.
	///	</summary>
	///	<param name="data">
	///		The DataVector object.
	///	</param>
	///	<param name="szName">
	///		A variable name used to identify <c>data</c> in the file.
	///	</param>
	void OutputDataVector(
		const DataVector<double> &data,
		const char * szName
	);

	///	<summary>
	///		Append a MATLAB-formatted DataMatrix to a file.
	///	</summary>
	///	<param name="data">
	///		The DataMatrix object.
	///	</param>
	///	<param name="szName">
	///		A variable name used to identify <c>data</c> in the file.
	///	</param>
	void OutputDataMatrix(
		const DataMatrix<double> &data,
		const char * szName
	);

	///	<summary>
	///		Append a MATLAB-formatted DataMatrix3D to a file.
	///	</summary>
	///	<param name="data">
	///		The DataMatrix3D object.
	///	</param>
	///	<param name="szName">
	///		A variable name used to identify <c>data</c> in the file.
	///	</param>
	void OutputDataMatrix3D(
		const DataMatrix3D<double> &data,
		const char * szName
	);

	///	<summary>
	///		Append a MATLAB-formatted DataMatrix4D to a file.
	///	</summary>
	///	<param name="data">
	///		The DataMatrix4D object.
	///	</param>
	///	<param name="szName">
	///		A variable name used to identify <c>data</c> in the file.
	///	</param>
	void OutputDataMatrix4D(
		const DataMatrix4D<double> &data,
		const char * szName
	);

private:
	///	<summary>
	///		A file pointer to a MATLAB file being written.
	///	</summary>
	FILE * m_fp;
};

///////////////////////////////////////////////////////////////////////////////

#endif

