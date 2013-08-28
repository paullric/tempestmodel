///////////////////////////////////////////////////////////////////////////////
///
///	\file    MatlabInput.cpp
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

#include "MatlabInput.h"

#include <cstring>

///////////////////////////////////////////////////////////////////////////////

const unsigned int MatlabInputFile::miINT8      = 1;
const unsigned int MatlabInputFile::miUINT8     = 2;
const unsigned int MatlabInputFile::miINT16     = 3;
const unsigned int MatlabInputFile::miUINT16    = 4;
const unsigned int MatlabInputFile::miINT32     = 5;
const unsigned int MatlabInputFile::miUINT32    = 6;
const unsigned int MatlabInputFile::miSINGLE    = 7;
const unsigned int MatlabInputFile::miDOUBLE    = 9;
const unsigned int MatlabInputFile::miINT64     = 12;
const unsigned int MatlabInputFile::miUINT64    = 13;
const unsigned int MatlabInputFile::miMATRIX    = 14;

const unsigned int MatlabInputFile::mxDOUBLE_CLASS = 6;

const unsigned int MatlabInputFile::mfCOMPLEX = 0x80;
const unsigned int MatlabInputFile::mfGLOBAL  = 0x40;
const unsigned int MatlabInputFile::mfLOGICAL = 0x20;
const unsigned int MatlabInputFile::mfNONE    = 0x00;

const unsigned int MatlabInputFile::mv8      = 8;
const unsigned int MatlabInputFile::mv12     = 12;

///////////////////////////////////////////////////////////////////////////////

void MatlabInputFile::ReadChunkHeader(
	unsigned int *uiArrayFlags,
	unsigned int *nDim,
	char *szArrayName,
	unsigned int nArrayName
) {
	unsigned int i;

	unsigned int uiType;
	unsigned int uiSize;

	char szNull[8];

	// Read in array flags
	fread(&uiType, sizeof(int), 1, m_fp);
	fread(&uiSize, sizeof(int), 1, m_fp);

	if (uiType != miUINT32) {
		_EXCEPTIONT("Array flags must be type miUINT32.");
	}

	fread(uiArrayFlags, sizeof(int), 2, m_fp);

	// Read in array dimensions
	fread(&uiType, sizeof(int), 1, m_fp);
	fread(&uiSize, sizeof(int), 1, m_fp);

	if (((uiSize % 4) != 0) || ((uiSize / 4) > 4)) {
		_EXCEPTION1("Invalid array size in array dimension parameter: %i",
			uiSize);
	}

	for (i = 0; i < (uiSize / 4); i++) {
		fread(nDim + i, sizeof(int), 1, m_fp);
	}
	if ((uiSize % 8) != 0) {
		fread(nDim + i, sizeof(int), 1, m_fp);
	}
	for (i = (uiSize / 4); i < 4; i++) {
		nDim[i] = 0;
	}

	// Read in array name
	fread(&uiType, sizeof(int), 1, m_fp);
	fread(&uiSize, sizeof(int), 1, m_fp);

	if (uiType != miINT8) {
		_EXCEPTIONT("Array name must be type miINT8.");
	}
	if (uiSize > nArrayName) {
		_EXCEPTION1("Array name must be at most %i characters.", nArrayName);
	}

	fread(szArrayName, sizeof(char), uiSize, m_fp);
	fread(szNull, sizeof(char), 7 - ((uiSize - 1) % 8), m_fp);

	szArrayName[uiSize] = '\0';
}

///////////////////////////////////////////////////////////////////////////////

void MatlabInputFile::Open(const char *szFilename, bool fVerbose) {

	// Check that this object does not correspond to an already-open
	// file pointer.
	if (m_fp != NULL) {
		_EXCEPTIONT("Attempting to open an already-open MATLAB file.");
	}

	// Open the file for reading
	m_fp = fopen(szFilename, "r");
	if (m_fp == NULL) {
		_EXCEPTION1("Unable to open file \'%s\' for reading.", szFilename);
	}

	// Header information
	char szHeader[125];

	char szVersion[2];

	char szEndian[2];

	// Set the version number
	szVersion[0] = 0;
	szVersion[1] = 1;

	// Header and version number
	fread(szHeader, sizeof(char), 124, m_fp);

	fread(szVersion, sizeof(char), 2, m_fp);

	// Determine endian of the local system
	fread(szEndian, sizeof(char), 2, m_fp);

	if ((szEndian[0] != 'I') || (szEndian[1] != 'M')) {
		_EXCEPTION2("Invalid endian formatting found: %c:%c.",
			szEndian[0], szEndian[1]);
	}

	if ((szVersion[0] != 0) || (szVersion[1] != 1)) {
		_EXCEPTION2("Unknown version information found: %i:%i",
			szVersion[0], szVersion[1]);
	}

	// Buffer data
	char szArrayName[256];

	// Current position in file
	long int lPos;

	// Type of data being read in (currently only supports miMATRIX)
	unsigned int uiType;

	// Size of data in bytes being read in
	unsigned int uiSize;

	// Array flags information
	unsigned int uiArrayFlags[2];

	// Array dimensions
	unsigned int nDim[4];

	// Index the MATLAB file
	for (;;) {
		// Read in type and size of next data chunk
		fread(&uiType, sizeof(int), 1, m_fp);
		fread(&uiSize, sizeof(int), 1, m_fp);

		// Check for end of file
		if (feof(m_fp)) {
			break;
		}

		// Only support miMATRIX
		if (uiType != miMATRIX) {
			_EXCEPTIONT("Reader only supports MATLAB MATRIX type.");
		}

		// Store current position in file
		lPos = ftell(m_fp);

		// Read in the chunk header
		ReadChunkHeader(uiArrayFlags, nDim, szArrayName, 256); 

		// Values
		fread(&uiType, sizeof(int), 1, m_fp);
		fread(&uiSize, sizeof(int), 1, m_fp);

		if (uiType != miDOUBLE) {
			_EXCEPTIONT("Array values must be stored as miDOUBLE.");
		}

		fseek(m_fp, uiSize, SEEK_CUR);

		// Add the new object to the map
		NamePositionPair namepospair;
		namepospair.first = szArrayName;
		namepospair.second = lPos;

		m_mapPos.insert(namepospair);

		if (fVerbose) {
			std::cout << "Found matrix: \"" << namepospair.first << "\" ";
			std::cout << "with " << uiSize << " elements.";
			std::cout << std::endl;
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void MatlabInputFile::Close() {
	if (m_fp == NULL) {
		return;
	}

	fclose(m_fp);
	m_fp = NULL;
}

///////////////////////////////////////////////////////////////////////////////

unsigned int MatlabInputFile::GetDataDimension(
	char * szName
) {
	// Check that the file is open
	if (m_fp == NULL) {
		_EXCEPTIONT("Attempting to read from unopened MATLAB file.");
	}

	// Array name from the header
	char szArrayName[256];

	// Type of data being read in (currently only supports miMATRIX)
	unsigned int uiArrayFlags[2];

	// Buffer counter
	unsigned int i;

	// Array dimensions
	unsigned int nDim[4];

	// Find the position of this array in the file
	NamePositionMap::iterator iter;

	iter = m_mapPos.find(szName);

	if (iter == m_mapPos.end()) {
		_EXCEPTION1("Could not identify Matrix %s in MATLAB file for reading.",
			szName);
	}

	// Seek to the position in the file
	fseek(m_fp, iter->second, SEEK_SET);

	// Read in the chunk header
	ReadChunkHeader(uiArrayFlags, nDim, szArrayName, 256); 

	// Determine the dimension
	for (i = 0; i < 4; i++) {
		if (nDim[i] == 0) {
			return i;
		}
	}

	return 4;
}

///////////////////////////////////////////////////////////////////////////////

 bool MatlabInputFile::HasMatrix(
	const char * szName
) {

	// Find the position of this array in the file
	NamePositionMap::iterator iter;

	iter = m_mapPos.find(szName);

	if (iter == m_mapPos.end()) {
		return false;
	}

	return true;
}

///////////////////////////////////////////////////////////////////////////////

void MatlabInputFile::InputDataVector(
	const char * szName,
	DataVector<double> &data,
	bool fVerbose
) {
	// Check that the file is open
	if (m_fp == NULL) {
		_EXCEPTIONT("Attempting to read from unopened MATLAB file.");
	}

	// Array name from the header
	char szArrayName[256];

	// Type of data being read in (currently only supports miMATRIX)
	unsigned int uiArrayFlags[2];

	// Type and Size
	unsigned int uiType;
	unsigned int uiSize;

	// Buffer counter
	unsigned int i;

	// Array dimensions
	unsigned int nDim[4];

	// Find the position of this array in the file
	NamePositionMap::iterator iter;

	iter = m_mapPos.find(szName);

	if (iter == m_mapPos.end()) {
		_EXCEPTION1("Could not identify Matrix %s in MATLAB file for reading.",
			szName);
	}

	// Seek to the position in the file
	fseek(m_fp, iter->second, SEEK_SET);

	// Read in the chunk header
	ReadChunkHeader(uiArrayFlags, nDim, szArrayName, 256); 

	// Check array size
	if (nDim[1] != 1) {
		_EXCEPTION4(
			"Specified Matrix is not of type DataVector:\n"
			"Found dimensions: %i, %i, %i, %i",
			nDim[0], nDim[1], nDim[2], nDim[3]);
	}

	// Values
	fread(&uiType, sizeof(int), 1, m_fp);
	fread(&uiSize, sizeof(int), 1, m_fp);

	if (uiType != miDOUBLE) {
		_EXCEPTIONT("Array values must be stored as miDOUBLE.");
	}

	// Resize the DataVector object
	data.Initialize(nDim[0]);

	for (i = 0; i < nDim[0]; i++) {
		fread(data + i, sizeof(double), 1, m_fp);
	}

	if (fVerbose) {
		std::cout << "Successfully loaded matrix \"" << szName << "\"";
		std::cout << std::endl;
	}
}

///////////////////////////////////////////////////////////////////////////////

void MatlabInputFile::InputDataMatrix(
	const char * szName,
	DataMatrix<double> &data,
	bool fVerbose
) {
	// Check that the file is open
	if (m_fp == NULL) {
		_EXCEPTIONT("Attempting to read from unopened MATLAB file.");
	}

	// Array name from the header
	char szArrayName[256];

	// Type of data being read in (currently only supports miMATRIX)
	unsigned int uiArrayFlags[2];

	// Type and Size
	unsigned int uiType;
	unsigned int uiSize;

	// Buffer counter
	unsigned int i;
	unsigned int j;

	// Array dimensions
	unsigned int nDim[4];

	// Find the position of this array in the file
	NamePositionMap::iterator iter;

	iter = m_mapPos.find(szName);

	if (iter == m_mapPos.end()) {
		_EXCEPTION1("Could not identify Matrix %s in MATLAB file for reading.",
			szName);
	}

	// Seek to the position in the file
	fseek(m_fp, iter->second, SEEK_SET);

	// Read in the chunk header
	ReadChunkHeader(uiArrayFlags, nDim, szArrayName, 256); 

	// Check array size
	if (nDim[2] != 0) {
		_EXCEPTION4(
			"Specified Matrix is not of type DataMatrix:\n"
			"Found dimensions: %i, %i, %i, %i",
			nDim[0], nDim[1], nDim[2], nDim[3]);
	}

	// Values
	fread(&uiType, sizeof(int), 1, m_fp);
	fread(&uiSize, sizeof(int), 1, m_fp);

	if (uiType != miDOUBLE) {
		_EXCEPTIONT("Array values must be stored as miDOUBLE.");
	}

	// Resize the DataVector object
	data.Initialize(nDim[0], nDim[1]);

	for (j = 0; j < nDim[1]; j++) {
	for (i = 0; i < nDim[0]; i++) {
		fread(&(data[i][j]), sizeof(double), 1, m_fp);
	}
	}

	if (fVerbose) {
		std::cout << "Successfully loaded matrix \"" << szName << "\"";
		std::cout << std::endl;
	}

}

///////////////////////////////////////////////////////////////////////////////

void MatlabInputFile::InputDataMatrix3D(
	const char * szName,
	DataMatrix3D<double> &data,
	bool fVerbose
) {
	// Check that the file is open
	if (m_fp == NULL) {
		_EXCEPTIONT("Attempting to read from unopened MATLAB file.");
	}

	// Array name from the header
	char szArrayName[256];

	// Type of data being read in (currently only supports miMATRIX)
	unsigned int uiArrayFlags[2];

	// Type and Size
	unsigned int uiType;
	unsigned int uiSize;

	// Buffer counter
	unsigned int i;
	unsigned int j;
	unsigned int k;

	// Array dimensions
	unsigned int nDim[4];

	// Find the position of this array in the file
	NamePositionMap::iterator iter;

	iter = m_mapPos.find(szName);

	if (iter == m_mapPos.end()) {
		_EXCEPTION1("Could not identify Matrix %s in MATLAB file for reading.",
			szName);
	}

	// Seek to the position in the file
	fseek(m_fp, iter->second, SEEK_SET);

	// Read in the chunk header
	ReadChunkHeader(uiArrayFlags, nDim, szArrayName, 256); 

	// Check array size
	if (nDim[3] != 0) {
		_EXCEPTION4(
			"Specified Matrix is not of type DataMatrix3D:\n"
			"Found dimensions: %i, %i, %i, %i",
			nDim[0], nDim[1], nDim[2], nDim[3]);
	}

	// Values
	fread(&uiType, sizeof(int), 1, m_fp);
	fread(&uiSize, sizeof(int), 1, m_fp);

	if (uiType != miDOUBLE) {
		_EXCEPTIONT("Array values must be stored as miDOUBLE.");
	}

	// Resize the DataVector object
	data.Initialize(nDim[0], nDim[1], nDim[2]);

	for (k = 0; k < nDim[2]; k++) {
	for (j = 0; j < nDim[1]; j++) {
	for (i = 0; i < nDim[0]; i++) {
		fread(&(data[i][j][k]), sizeof(double), 1, m_fp);
	}
	}
	}
	
	if (fVerbose) {
		std::cout << "Successfully loaded matrix \"" << szName << "\"";
		std::cout << std::endl;
	}

}

///////////////////////////////////////////////////////////////////////////////

void MatlabInputFile::InputDataMatrix4D(
	const char * szName,
	DataMatrix4D<double> &data,
	bool fVerbose
) {
	// Check that the file is open
	if (m_fp == NULL) {
		_EXCEPTIONT("Attempting to read from unopened MATLAB file.");
	}

	// Array name from the header
	char szArrayName[256];

	// Type of data being read in (currently only supports miMATRIX)
	unsigned int uiArrayFlags[2];

	// Type and Size
	unsigned int uiType;
	unsigned int uiSize;

	// Buffer counter
	unsigned int i;
	unsigned int j;
	unsigned int k;
	unsigned int l;

	// Array dimensions
	unsigned int nDim[4];

	// Find the position of this array in the file
	NamePositionMap::iterator iter;

	iter = m_mapPos.find(szName);

	if (iter == m_mapPos.end()) {
		_EXCEPTION1("Could not identify Matrix %s in MATLAB file for reading.",
			szName);
	}

	// Seek to the position in the file
	fseek(m_fp, iter->second, SEEK_SET);

	// Read in the chunk header
	ReadChunkHeader(uiArrayFlags, nDim, szArrayName, 256); 

	// Check array size
	if (nDim[3] == 0) {
		_EXCEPTION4(
			"Specified Matrix is not of type DataMatrix4D:\n"
			"Found dimensions: %i, %i, %i, %i",
			nDim[0], nDim[1], nDim[2], nDim[3]);
	}

	// Values
	fread(&uiType, sizeof(int), 1, m_fp);
	fread(&uiSize, sizeof(int), 1, m_fp);

	if (uiType != miDOUBLE) {
		_EXCEPTIONT("Array values must be stored as miDOUBLE.");
	}

	// Resize the DataVector object
	data.Initialize(nDim[0], nDim[1], nDim[2], nDim[3]);

	for (l = 0; l < nDim[3]; l++) {
	for (k = 0; k < nDim[2]; k++) {
	for (j = 0; j < nDim[1]; j++) {
	for (i = 0; i < nDim[0]; i++) {
		fread(&(data[i][j][k][l]), sizeof(double), 1, m_fp);
	}
	}
	}
	}
	
	if (fVerbose) {
		std::cout << "Successfully loaded matrix \"" << szName << "\"";
		std::cout << std::endl;
	}

}

///////////////////////////////////////////////////////////////////////////////

