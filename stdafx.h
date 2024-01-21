// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#pragma once

#include "targetver.h"

#include <stdio.h>
#include <tchar.h>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <sstream>
#include <Windows.h>

#include "DataTypeDefinitions.h"
#include "Matrix\inc\Matrix.h"
#include "LeastSquares\inc\LeastSquares.h"
using namespace std;

#define NETWORK_NUM_OPTIONS ((ULONG) 19)
#define NETWORK_INPUTFILE_INDEX ((ULONG) 0)
#define NETWORK_P1X_INDEX ((ULONG) 1)
#define NETWORK_P1Y_INDEX ((ULONG) 2)
#define NETWORK_P2X_INDEX ((ULONG) 3)
#define NETWORK_P2Y_INDEX ((ULONG) 4)
#define NETWORK_P3X_INDEX ((ULONG) 5)
#define NETWORK_P3Y_INDEX ((ULONG) 6)
#define NETWORK_P4X_INDEX ((ULONG) 7)
#define NETWORK_P4Y_INDEX ((ULONG) 8)
#define NETWORK_P5X_INDEX ((ULONG) 9)
#define NETWORK_P5Y_INDEX ((ULONG) 10)
#define NETWORK_P6X_INDEX ((ULONG) 11)
#define NETWORK_P6Y_INDEX ((ULONG) 12)
#define NETWORK_FIXPOINTS_INDEX ((ULONG) 13)
#define NETWORK_RESA_INDEX ((ULONG) 14)
#define NETWORK_RESB_INDEX ((ULONG) 15)
#define NETWORK_RESC_INDEX ((ULONG) 16)
#define NETWORK_OUTPATH_INDEX ((ULONG) 17)
#define NETWORK_METHOD_INDEX ((ULONG) 18)

//----------------------------------------------------------------------------
// Assigns coordinates to the coordinate map
// pmpCoordinates_: Coordinate map
// sFixedPoints_: Comma separated list of points that have fixed coordinates
// pdConstants_: Assign constants to vector
//----------------------------------------------------------------------------
void AssignCoordinates(map<ULONG, PointCoordinate>* pmpCoordinates_, string sFixedPoints_, DOUBLE* pdConstants_);

//----------------------------------------------------------------------------
// Determines the points that have fixed coordinates
// sFixedPoints_: Comma-separated list of points that have fixed coordinates
// pvFixedPoints_: List of points having fixed coordinates
//----------------------------------------------------------------------------
void FixedPointsVector(const string& sFixedPoints_, vector<ULONG>* pvFixedPoints_);

//----------------------------------------------------------------------------
// Reads input file
// sInputFile_: Full path to input file
// pvMeasurements_: Measurement list
//----------------------------------------------------------------------------
void ReadInputFile(string sInputFile_, vector<DistanceMeasurement>* pvMeasurements_);

//----------------------------------------------------------------------------
// Parse data item, determine if it is numeric
// sData_: string to parse
//----------------------------------------------------------------------------
BOOLEANO bParseData(string sData_);

//----------------------------------------------------------------------------
// Output results to file
// pmpCoordinates_: Computed coordinates
// pclLS_: Object of type LeastSquares with state covariance and residuals
// sOutPath_: Output path for file output
//----------------------------------------------------------------------------
void Output(map<ULONG, PointCoordinate>* pmpCoordinates_, LeastSquares* pclLS_, string sOutPath_);

//----------------------------------------------------------------------------
// Output residuals to file
// pvMeasurements_: Measurements
// pclLS_: Object of type LeastSquares with state covariance and residuals
// sOutPath_: Output path for file output
//----------------------------------------------------------------------------
void Residuals(vector<DistanceMeasurement>* pvMeasurements_, LeastSquares* pclLS_, string sOutPath_);
