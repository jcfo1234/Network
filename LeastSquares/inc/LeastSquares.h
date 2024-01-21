//----------------------------------------------------------------------
// LSCommon.h
// Module that handles common least squares functions
//----------------------------------------------------------------------

#ifndef LEASTSQUARES_H
#define LEASTSQUARES_H

#include <stdio.h>
#include <tchar.h>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <sstream>

#include "..\..\DataTypeDefinitions.h"
#include "..\..\Matrix\inc\Matrix.h"
using namespace std;

#define LEASTSQUARES_WEIGHT_A_INDEX ((ULONG) 0)
#define LEASTSQUARES_WEIGHT_B_INDEX ((ULONG) 1)
#define LEASTSQUARES_WEIGHT_C_INDEX ((ULONG) 2)
#define LEASTSQUARES_ERROR_TOLERANCE ((DOUBLE) 1e-6)

class LeastSquares
{
public:
   // Default constructor
   LeastSquares() {};
   // destructor
   ~LeastSquares() {};
   // Class assignment operator
   LeastSquares& operator = (const LeastSquares &other);
   // Copy Constructor
   LeastSquares(const LeastSquares &other);
   // Get residuals
   CMatrix GetResiduals() { return this->clResiduals; };
   // Get state covariance matrix
   CMatrix GetStateCovariance() { return this->clCX; };
   // Set residuals
   void SetResiduals(CMatrix& other) { this->clResiduals = other; };
   // Set state covariance matrix
   void SetStateCovariance(CMatrix& other) { this->clCX = other; };
private:
   CMatrix clResiduals;
   CMatrix clCX;
};

//------------------------------------------------------------------------------
// Find the number of fixed points
// pmpCoordinates_: Current iteration coordinates
//------------------------------------------------------------------------------
void FixedPoints(map<ULONG, PointCoordinate>* pmpCoordinates_, vector<ULONG>* pvPointID_);

//------------------------------------------------------------------------------
// Compute misclosures from measurement vector
// pmpCoordinates_: Current iteration coordinates
// pvMeasurements_: Measurement vector
//------------------------------------------------------------------------------
CMatrix Misclosures(map<ULONG, PointCoordinate>* pmpCoordinates_, vector<DistanceMeasurement>* pvMeasurements_);

//------------------------------------------------------------------------------
// Compute linear measurement matrix A
// pmpCoordinates_: Current iteration coordinates
// pvMeasurements_: Measurement vector
//------------------------------------------------------------------------------
CMatrix MeasurementMatrix(map<ULONG, PointCoordinate>* pmpCoordinates_, vector<DistanceMeasurement>* pvMeasurements_);

//------------------------------------------------------------------------------
// Compute observation covariance matrix Cl
// pvMeasurements_: Measurement vector
//------------------------------------------------------------------------------
CMatrix ObservationCovariance(vector<DistanceMeasurement>* pvMeasurements_);

//------------------------------------------------------------------------------
// Compute least squares robust weight matrix
// pclA_: Measurements matrix
// pclObsCov_: Measurements observation covariance matrix
// pclMisc_: Misclosures
// pdConstants_: Hampel's constants
//------------------------------------------------------------------------------
CMatrix HampelWeight(CMatrix* pclObsCov_, CMatrix* pclA_, CMatrix* pclMisc_, DOUBLE* pdConstants_);

//------------------------------------------------------------------------------
// Compute least squares robust weight matrix
// pclA_: Measurements matrix
// pclObsCov_: Measurements observation covariance matrix
// pclMisc_: Misclosures
// pdConstants_: Hampel's constants
//------------------------------------------------------------------------------
CMatrix HuberWeight(CMatrix* pclObsCov_, CMatrix* pclA_, CMatrix* pclMisc_, DOUBLE* pdConstants_);

//------------------------------------------------------------------------------
// Compute least squares robust weight matrix
// pclA_: Measurements matrix
// pclObsCov_: Measurements observation covariance matrix
// pclMisc_: Misclosures
// pdConstants_: Hampel's constants
//------------------------------------------------------------------------------
CMatrix AndrewWeight(CMatrix* pclObsCov_, CMatrix* pclA_, CMatrix* pclMisc_, DOUBLE* pdConstants_);

//------------------------------------------------------------------------------
// Compute least squares robust weight matrix
// pclObsCov_: Measurements observation covariance matrix
//------------------------------------------------------------------------------
CMatrix LSWeight(CMatrix* pclObsCov_, CMatrix* pclA_ = NULL, CMatrix* pclMisc_=NULL, DOUBLE* pdConstants_=NULL);

//------------------------------------------------------------------------------
// Compute least squares robust state variance covariance
// pclWeightMatrix_: Weight matrix for robust least squares
// pclDesignMatrixA_: Measurements matrix
//------------------------------------------------------------------------------
CMatrix ComputeStateCovariance(CMatrix* pclWeightMatrix_, CMatrix* pclDesignMatrixA_);

//------------------------------------------------------------------------------
// Compute least squares robust state correction
// pclWeightMatrix_: Weight matrix for robust least squares
// pclDesignMatrixA_: Measurements matrix
// pclMisclosure_: Misclosure vector
//------------------------------------------------------------------------------
CMatrix ComputeStateCorrection(CMatrix* pclWeightMatrix_, CMatrix* pclDesignMatrixA_, CMatrix* pclMisclosure_);

//------------------------------------------------------------------------------
// Refine network coordinates through lest squares solution
// pmpCoordinates_: Current iteration coordinates
// pclDelta_: Correction vector for coordinates
//------------------------------------------------------------------------------
void Update(map<ULONG, PointCoordinate>* pmpCoordinates_, CMatrix* pclDelta_);

//------------------------------------------------------------------------------
// Process data
// pmpCoordinates_: Initial iteration coordinates
// pvMeasurements_: Measurement vector
// sMethod_: Selects the robust least squares method to apply
// pdConstants_: Cutoff values for robust least squares methods
//------------------------------------------------------------------------------
LeastSquares Process(map<ULONG, PointCoordinate>* pmpCoordinates_, vector<DistanceMeasurement>* pvMeasurements_, string sMethod_, DOUBLE* pdConstants_);

#endif