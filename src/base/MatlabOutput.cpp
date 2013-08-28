///////////////////////////////////////////////////////////////////////////////
///
///	\file    MatlabOutput.cpp
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

#include "MatlabOutput.h"

#include <cstring>

///////////////////////////////////////////////////////////////////////////////

const unsigned int MatlabOutputFile::miINT8      = 1;
const unsigned int MatlabOutputFile::miUINT8     = 2;
const unsigned int MatlabOutputFile::miINT16     = 3;
const unsigned int MatlabOutputFile::miUINT16    = 4;
const unsigned int MatlabOutputFile::miINT32     = 5;
const unsigned int MatlabOutputFile::miUINT32    = 6;
const unsigned int MatlabOutputFile::miSINGLE    = 7;
const unsigned int MatlabOutputFile::miDOUBLE    = 9;
const unsigned int MatlabOutputFile::miINT64     = 12;
const unsigned int MatlabOutputFile::miUINT64    = 13;
const unsigned int MatlabOutputFile::miMATRIX    = 14;

const unsigned int MatlabOutputFile::mxDOUBLE_CLASS = 6;

const unsigned int MatlabOutputFile::mfCOMPLEX = 0x80;
const unsigned int MatlabOutputFile::mfGLOBAL  = 0x40;
const unsigned int MatlabOutputFile::mfLOGICAL = 0x20;
const unsigned int MatlabOutputFile::mfNONE    = 0x00;

const unsigned int MatlabOutputFile::mv8      = 8;
const unsigned int MatlabOutputFile::mv12     = 12;
const unsigned int MatlabOutputFile::mv16     = 16;

///////////////////////////////////////////////////////////////////////////////

void MatlabOutputFile::Open(const char *szFilename) {

	// Check that this object does not correspond to an already-open
	// file pointer.
	if (m_fp != NULL) {
		_EXCEPTIONT("Attempting to open an already-open MATLAB file.");
	}

	// Open the file for writing
	m_fp = fopen(szFilename, "w");
	if (m_fp == NULL) {
		_EXCEPTION1("Unable to open file \'%s\' for writing.", szFilename);
	}

	// Header information
	const char * szHeader =
		"MATLAB 5.0 MAT-file, Platform: MACI, Created on: M"
		"on Jan  1 00:00:00 2000                           "
		"                        ";

	char szVersion[2];

	char szEndian[2];

	// Set the version number
	szVersion[0] = 0;
	szVersion[1] = 1;

	// Determine endian of the local system
	union _EndianUnion { unsigned int i; unsigned char bytes[4]; };

	_EndianUnion eu;

	eu.i = 0xFF000000;

	if (eu.bytes[0] == 0xFF) {
		szEndian[0] = 'M';
		szEndian[1] = 'I';
	} else {
		szEndian[0] = 'I';
		szEndian[1] = 'M';
	}

	fwrite(szHeader, sizeof(char), 124, m_fp);

	fwrite(szVersion, sizeof(char), 2, m_fp);

	fwrite(szEndian, sizeof(char), 2, m_fp);
}

///////////////////////////////////////////////////////////////////////////////

void MatlabOutputFile::Close() {
	if (m_fp == NULL) {
		return;
	}

	fclose(m_fp);
	m_fp = NULL;
}

///////////////////////////////////////////////////////////////////////////////

void MatlabOutputFile::OutputDouble(
	const double data,
	const char * szName
) {
	// Check that the file is open
	if (m_fp == NULL) {
		_EXCEPTIONT("Attempting to write to unopened MATLAB file.");
	}

	// Array of NULLs (useful for writing to file)
	static const unsigned char mNULL8[8] = "\0\0\0\0\0\0\0";

	// Buffer variables
	unsigned int i;
	unsigned int nValue;
	unsigned int nNameLength;

	// Check name to ensure it is given
	nNameLength = strlen(szName);

	if (nNameLength == 0) {
		_EXCEPTIONT("OutputDataVector received name string of length 0.");
	}

	// Calculate array name length
	nValue = (((nNameLength - 1) / 8) + 1) * 8;

	// Write matrix header
	fwrite(&miMATRIX, sizeof(int), 1, m_fp);

	nValue = nValue
	         + 16    // Array flags
	         + 16    // Array dimensions
			 + 8     // Array name header
	         + 8;    // Data header

	// Add length of included data
	nValue = nValue + sizeof(double);

	fwrite(&nValue, sizeof(int), 1, m_fp);

	// Array flags
	fwrite(&miUINT32, sizeof(int), 1, m_fp);
	fwrite(&mv8, sizeof(int), 1, m_fp);

	nValue = mxDOUBLE_CLASS | mfNONE;

	fwrite(&nValue, sizeof(int), 1, m_fp);

	fwrite(&mNULL8, sizeof(char), 4, m_fp);

	// Array dimensions
	fwrite(&miINT32, sizeof(int), 1, m_fp);
	fwrite(&mv8, sizeof(int), 1, m_fp);

	nValue = 1;
	fwrite(&nValue, sizeof(int), 1, m_fp);

	nValue = 1;
	fwrite(&nValue, sizeof(int), 1, m_fp);

	// Array name
	fwrite(&miINT8, sizeof(int), 1, m_fp);
	fwrite(&nNameLength, sizeof(int), 1, m_fp);

	fwrite(szName, sizeof(char), nNameLength, m_fp);
	fwrite(mNULL8, sizeof(char), 7 - ((nNameLength - 1) % 8), m_fp);

	// Values
	nValue = sizeof(double);

	fwrite(&miDOUBLE, sizeof(int), 1, m_fp);
	fwrite(&nValue, sizeof(int), 1, m_fp);

	fwrite(&data, sizeof(double), 1, m_fp);
}

///////////////////////////////////////////////////////////////////////////////

void MatlabOutputFile::OutputDataVector(
	const DataVector<double> &data,
	const char * szName
) {
	// Check that the file is open
	if (m_fp == NULL) {
		_EXCEPTIONT("Attempting to write to unopened MATLAB file.");
	}

	// Array of NULLs (useful for writing to file)
	static const unsigned char mNULL8[8] = "\0\0\0\0\0\0\0";

	// Buffer variables
	unsigned int i;
	unsigned int nValue;
	unsigned int nNameLength;

	// Check name to ensure it is given
	nNameLength = strlen(szName);

	if (nNameLength == 0) {
		_EXCEPTIONT("OutputDataVector received name string of length 0.");
	}

	// Calculate array name length
	nValue = (((nNameLength - 1) / 8) + 1) * 8;

	// Write matrix header
	fwrite(&miMATRIX, sizeof(int), 1, m_fp);

	nValue = nValue
	         + 16    // Array flags
	         + 16    // Array dimensions
			 + 8     // Array name header
	         + 8;    // Data header

	// Add length of included data
	nValue = nValue + sizeof(double) * data.GetRows();

	fwrite(&nValue, sizeof(int), 1, m_fp);

	// Array flags
	fwrite(&miUINT32, sizeof(int), 1, m_fp);
	fwrite(&mv8, sizeof(int), 1, m_fp);

	nValue = mxDOUBLE_CLASS | mfNONE;

	fwrite(&nValue, sizeof(int), 1, m_fp);

	fwrite(&mNULL8, sizeof(char), 4, m_fp);

	// Array dimensions
	fwrite(&miINT32, sizeof(int), 1, m_fp);
	fwrite(&mv8, sizeof(int), 1, m_fp);

	nValue = data.GetRows();
	fwrite(&nValue, sizeof(int), 1, m_fp);

	nValue = 1;
	fwrite(&nValue, sizeof(int), 1, m_fp);

	// Array name
	fwrite(&miINT8, sizeof(int), 1, m_fp);
	fwrite(&nNameLength, sizeof(int), 1, m_fp);

	fwrite(szName, sizeof(char), nNameLength, m_fp);
	fwrite(mNULL8, sizeof(char), 7 - ((nNameLength - 1) % 8), m_fp);

	// Values
	nValue = data.GetRows() * sizeof(double);

	fwrite(&miDOUBLE, sizeof(int), 1, m_fp);
	fwrite(&nValue, sizeof(int), 1, m_fp);

	for (i = 0; i < data.GetRows(); i++) {
		fwrite(&(data[i]), sizeof(double), 1, m_fp);
	}
}

///////////////////////////////////////////////////////////////////////////////

void MatlabOutputFile::OutputDataMatrix(
	const DataMatrix<double> &data,
	const char * szName
) {
	// Check that the file is open
	if (m_fp == NULL) {
		_EXCEPTIONT("Attempting to write to unopened MATLAB file.");
	}

	// Array of NULLs (useful for writing to file)
	static const unsigned char mNULL8[8] = "\0\0\0\0\0\0\0";

	// Buffer variables
	unsigned int i;
	unsigned int j;
	unsigned int nValue;
	unsigned int nNameLength;

	// Check name to ensure it is given
	nNameLength = strlen(szName);

	if (nNameLength == 0) {
		_EXCEPTIONT("OutputDataMatrix received name string of length 0.");
	}

	// Calculate array name length
	nValue = (((nNameLength - 1) / 8) + 1) * 8;

	// Write matrix header
	fwrite(&miMATRIX, sizeof(int), 1, m_fp);

	nValue = nValue
	         + 16    // Array flags
	         + 16    // Array dimensions
			 + 8     // Array name header
	         + 8;    // Data header

	// Add length of included data
	nValue = nValue + sizeof(double) * data.GetRows() * data.GetColumns();

	fwrite(&nValue, sizeof(int), 1, m_fp);

	// Array flags
	fwrite(&miUINT32, sizeof(int), 1, m_fp);
	fwrite(&mv8, sizeof(int), 1, m_fp);

	nValue = mxDOUBLE_CLASS | mfNONE;

	fwrite(&nValue, sizeof(int), 1, m_fp);

	fwrite(&mNULL8, sizeof(char), 4, m_fp);

	// Array dimensions
	fwrite(&miINT32, sizeof(int), 1, m_fp);
	fwrite(&mv8, sizeof(int), 1, m_fp);

	nValue = data.GetRows();
	fwrite(&nValue, sizeof(int), 1, m_fp);

	nValue = data.GetColumns();
	fwrite(&nValue, sizeof(int), 1, m_fp);

	// Array name
	fwrite(&miINT8, sizeof(int), 1, m_fp);
	fwrite(&nNameLength, sizeof(int), 1, m_fp);

	fwrite(szName, sizeof(char), nNameLength, m_fp);
	fwrite(mNULL8, sizeof(char), 7 - ((nNameLength - 1) % 8), m_fp);

	// Values
	nValue = data.GetColumns() * data.GetRows() * sizeof(double);

	fwrite(&miDOUBLE, sizeof(int), 1, m_fp);
	fwrite(&nValue, sizeof(int), 1, m_fp);

	for (j = 0; j < data.GetColumns(); j++) {
	for (i = 0; i < data.GetRows(); i++) {
		fwrite(&(data[i][j]), sizeof(double), 1, m_fp);
	}
	}
}

///////////////////////////////////////////////////////////////////////////////

void MatlabOutputFile::OutputDataMatrix3D(
	const DataMatrix3D<double> &data,
	const char * szName
) {
	// Check that the file is open
	if (m_fp == NULL) {
		_EXCEPTIONT("Attempting to write to unopened MATLAB file.");
	}

	// Array of NULLs (useful for writing to file)
	static const unsigned char mNULL8[8] = "\0\0\0\0\0\0\0";

	// Buffer variables
	unsigned int i;
	unsigned int j;
	unsigned int k;
	unsigned int nValue;
	unsigned int nNameLength;

	// Check name to ensure it is given
	nNameLength = strlen(szName);

	if (nNameLength == 0) {
		_EXCEPTIONT("OutputDataMatrix3D received name string of length 0.");
	}

	// Calculate array name length
	nValue = (((nNameLength - 1) / 8) + 1) * 8;

	// Write matrix header
	fwrite(&miMATRIX, sizeof(int), 1, m_fp);

	nValue = nValue
	         + 16    // Array flags
	         + 24    // Array dimensions
			 + 8     // Array name header
	         + 8;    // Data header

	// Add length of included data
	nValue += sizeof(double)
		* data.GetRows()
		* data.GetColumns()
		* data.GetSubColumns();

	fwrite(&nValue, sizeof(int), 1, m_fp);

	// Array flags
	fwrite(&miUINT32, sizeof(int), 1, m_fp);
	fwrite(&mv8, sizeof(int), 1, m_fp);

	nValue = mxDOUBLE_CLASS | mfNONE;

	fwrite(&nValue, sizeof(int), 1, m_fp);

	fwrite(&mNULL8, sizeof(char), 4, m_fp);

	// Array dimensions
	fwrite(&miINT32, sizeof(int), 1, m_fp);
	fwrite(&mv12, sizeof(int), 1, m_fp);

	nValue = data.GetRows();
	fwrite(&nValue, sizeof(int), 1, m_fp);

	nValue = data.GetColumns();
	fwrite(&nValue, sizeof(int), 1, m_fp);

	nValue = data.GetSubColumns();
	fwrite(&nValue, sizeof(int), 1, m_fp);

	fwrite(&mNULL8, sizeof(char), 4, m_fp);

	// Array name
	fwrite(&miINT8, sizeof(int), 1, m_fp);
	fwrite(&nNameLength, sizeof(int), 1, m_fp);

	fwrite(szName, sizeof(char), nNameLength, m_fp);
	fwrite(mNULL8, sizeof(char), 7 - ((nNameLength - 1) % 8), m_fp);

	// Values
	nValue = sizeof(double)
		* data.GetSubColumns()
		* data.GetColumns()
		* data.GetRows();

	fwrite(&miDOUBLE, sizeof(int), 1, m_fp);
	fwrite(&nValue, sizeof(int), 1, m_fp);

	for (k = 0; k < data.GetSubColumns(); k++) {
	for (j = 0; j < data.GetColumns(); j++) {
	for (i = 0; i < data.GetRows(); i++) {
		fwrite(&(data[i][j][k]), sizeof(double), 1, m_fp);
	}
	}
	}
}

///////////////////////////////////////////////////////////////////////////////

void MatlabOutputFile::OutputDataMatrix4D(
	const DataMatrix4D<double> &data,
	const char * szName
) {
	// Check that the file is open
	if (m_fp == NULL) {
		_EXCEPTIONT("Attempting to write to unopened MATLAB file.");
	}

	// Array of NULLs (useful for writing to file)
	static const unsigned char mNULL8[8] = "\0\0\0\0\0\0\0";

	// Buffer variables
	unsigned int i;
	unsigned int j;
	unsigned int k;
	unsigned int l;
	unsigned int nValue;
	unsigned int nNameLength;

	// Check name to ensure it is given
	nNameLength = strlen(szName);

	if (nNameLength == 0) {
		_EXCEPTIONT("OutputDataMatrix3D received name string of length 0.");
	}

	// Calculate array name length
	nValue = (((nNameLength - 1) / 8) + 1) * 8;

	// Write matrix header
	fwrite(&miMATRIX, sizeof(int), 1, m_fp);

	nValue = nValue
	         + 16    // Array flags
	         + 24    // Array dimensions
			 + 8     // Array name header
	         + 8;    // Data header

	// Add length of included data
	nValue += sizeof(double)
		* data.GetSize(0)
		* data.GetSize(1)
		* data.GetSize(2)
		* data.GetSize(3);

	fwrite(&nValue, sizeof(int), 1, m_fp);

	// Array flags
	fwrite(&miUINT32, sizeof(int), 1, m_fp);
	fwrite(&mv8, sizeof(int), 1, m_fp);

	nValue = mxDOUBLE_CLASS | mfNONE;

	fwrite(&nValue, sizeof(int), 1, m_fp);

	fwrite(&mNULL8, sizeof(char), 4, m_fp);

	// Array dimensions
	fwrite(&miINT32, sizeof(int), 1, m_fp);
	fwrite(&mv16, sizeof(int), 1, m_fp);

	nValue = data.GetSize(0);
	fwrite(&nValue, sizeof(int), 1, m_fp);

	nValue = data.GetSize(1);
	fwrite(&nValue, sizeof(int), 1, m_fp);

	nValue = data.GetSize(2);
	fwrite(&nValue, sizeof(int), 1, m_fp);

	nValue = data.GetSize(3);
	fwrite(&nValue, sizeof(int), 1, m_fp);

	// Array name
	fwrite(&miINT8, sizeof(int), 1, m_fp);
	fwrite(&nNameLength, sizeof(int), 1, m_fp);

	fwrite(szName, sizeof(char), nNameLength, m_fp);
	fwrite(mNULL8, sizeof(char), 7 - ((nNameLength - 1) % 8), m_fp);

	// Values
	nValue = sizeof(double)
		* data.GetSize(0)
		* data.GetSize(1)
		* data.GetSize(2)
		* data.GetSize(3);

	fwrite(&miDOUBLE, sizeof(int), 1, m_fp);
	fwrite(&nValue, sizeof(int), 1, m_fp);

	for (l = 0; l < data.GetSize(3); l++) {
	for (k = 0; k < data.GetSize(2); k++) {
	for (j = 0; j < data.GetSize(1); j++) {
	for (i = 0; i < data.GetSize(0); i++) {
		fwrite(&(data[i][j][k][l]), sizeof(double), 1, m_fp);
	}
	}
	}
	}
}

///////////////////////////////////////////////////////////////////////////////

