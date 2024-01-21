//----------------------------------------------------------------------
// LeastSquares.cpp
// Source file of LeastSquares.h
// Source code for Least Squares algorithm
//----------------------------------------------------------------------
#include "..\inc\LeastSquares.h"

// Function pointer to the appropriate data processing
typedef CMatrix(*RobustWeightMatrix)(CMatrix* pclObsCov_, CMatrix* pclA_, CMatrix* pclMisc_, DOUBLE* pdConstants_);

//------------------------------------------------------------------------------
LeastSquares::LeastSquares(const LeastSquares &other)
{
   this->clResiduals = other.clResiduals;
   this->clCX = other.clCX;
}

//------------------------------------------------------------------------------
LeastSquares& LeastSquares::operator=(const LeastSquares &other)
{
   this->clResiduals = other.clResiduals;
   this->clCX = other.clCX;
   // Copied class
   return *this;
}

//------------------------------------------------------------------------------
void FixedPoints(map<ULONG, PointCoordinate>* pmpCoordinates_, vector<ULONG>* pvPointID_)
{
   pvPointID_->clear();
   // Point coordinates iteration
   for (map<ULONG, PointCoordinate>::iterator pIterCoord = pmpCoordinates_->begin(); pIterCoord != pmpCoordinates_->end(); pIterCoord++)
   {
      // Save the point that is not fixed to a vector of points
      if (!pIterCoord->second.bFIXED)
      {
         pvPointID_->push_back(pIterCoord->first);
      }
   }
   sort(pvPointID_->begin(), pvPointID_->end());
}

//------------------------------------------------------------------------------
CMatrix Misclosures(map<ULONG, PointCoordinate>* pmpCoordinates_, vector<DistanceMeasurement>* pvMeasurements_)
{
   // Misclosure vector with dimensions equal to the measurements vector
   CMatrix clMisclosures("Misclosures", pvMeasurements_->size(), 1);
   // Distance pseudo-measurement
   DOUBLE dDistance = 0;
   DOUBLE dDistX = 0;
   DOUBLE dDistY = 0;
   // Index of misclosure vector
   INT iIndex = 0;

   // Compute misclosure vector
   for (vector<DistanceMeasurement>::iterator pIterMeas = pvMeasurements_->begin(); pIterMeas != pvMeasurements_->end(); pIterMeas++)
   {
      // Make sure points read from file exist in the database of points
      if (pmpCoordinates_->find(pIterMeas->ulPointA) != pmpCoordinates_->end() && pmpCoordinates_->find(pIterMeas->ulPointB) != pmpCoordinates_->end())
      {
         // X distance
         dDistX = pmpCoordinates_->at(pIterMeas->ulPointA).P_X - pmpCoordinates_->at(pIterMeas->ulPointB).P_X;
         // Y distance
         dDistY = pmpCoordinates_->at(pIterMeas->ulPointA).P_Y - pmpCoordinates_->at(pIterMeas->ulPointB).P_Y;
         // Distance
         dDistance = sqrt(pow(dDistX, 2) + pow(dDistY, 2));
         // Misclosure component
         clMisclosures.SetComponent(iIndex, 0, pIterMeas->dDistMeas - dDistance);
         iIndex = iIndex + 1;
      }
   }

   return clMisclosures;
}

//------------------------------------------------------------------------------
CMatrix MeasurementMatrix(map<ULONG, PointCoordinate>* pmpCoordinates_, vector<DistanceMeasurement>* pvMeasurements_)
{
   vector<ULONG> vPoints;
   FixedPoints(pmpCoordinates_, &vPoints);
   // Measurement matrix
   CMatrix clA("Measurement Matrix", pvMeasurements_->size(), 2 * vPoints.size());
   // Index of misclosure vector
   INT iIndex = 0;
   // Column Indexes
   INT iColXA = -1;
   INT iColYA = -1;
   INT iColXB = -1;
   INT iColYB = -1;
   // Distances
   DOUBLE dDistX = 0;
   DOUBLE dDistY = 0;
   DOUBLE dDistance = 0;

   // Iterate over measurement vector
   for (vector<DistanceMeasurement>::iterator pIterMeas = pvMeasurements_->begin(); pIterMeas != pvMeasurements_->end(); pIterMeas++)
   {
      // Make sure points read from file exist in the database of points
      if (pmpCoordinates_->find(pIterMeas->ulPointA) != pmpCoordinates_->end() && pmpCoordinates_->find(pIterMeas->ulPointB) != pmpCoordinates_->end())
      {
         // Compute the column position for a measurement reading
         if (find(vPoints.begin(), vPoints.end(), pIterMeas->ulPointA) != vPoints.end())
         {
            iColXA = 2 * (find(vPoints.begin(), vPoints.end(), pIterMeas->ulPointA) - vPoints.begin());
            iColYA = 2 * (find(vPoints.begin(), vPoints.end(), pIterMeas->ulPointA) - vPoints.begin()) + 1;
            // Compute the matrix entry for point A
            dDistX = pmpCoordinates_->at(pIterMeas->ulPointA).P_X - pmpCoordinates_->at(pIterMeas->ulPointB).P_X;
            dDistY = pmpCoordinates_->at(pIterMeas->ulPointA).P_Y - pmpCoordinates_->at(pIterMeas->ulPointB).P_Y;
            dDistance = sqrt(pow(dDistX, 2) + pow(dDistY, 2));
            clA.SetComponent(iIndex, iColXA, -dDistX / dDistance);
            clA.SetComponent(iIndex, iColYA, -dDistY / dDistance);
         }
         // Compute the column position for a measurement reading
         if (find(vPoints.begin(), vPoints.end(), pIterMeas->ulPointB) != vPoints.end())
         {
            iColXB = 2 * (find(vPoints.begin(), vPoints.end(), pIterMeas->ulPointB) - vPoints.begin());
            iColYB = 2 * (find(vPoints.begin(), vPoints.end(), pIterMeas->ulPointB) - vPoints.begin()) + 1;
            // Compute the matrix entry for point B
            dDistX = pmpCoordinates_->at(pIterMeas->ulPointB).P_X - pmpCoordinates_->at(pIterMeas->ulPointA).P_X;
            dDistY = pmpCoordinates_->at(pIterMeas->ulPointB).P_Y - pmpCoordinates_->at(pIterMeas->ulPointA).P_Y;
            dDistance = sqrt(pow(dDistX, 2) + pow(dDistY, 2));
            clA.SetComponent(iIndex, iColXB, -dDistX / dDistance);
            clA.SetComponent(iIndex, iColYB, -dDistY / dDistance);
         }
         iIndex = iIndex + 1;
      }
   }

   return clA;
}

//------------------------------------------------------------------------------
CMatrix ObservationCovariance(vector<DistanceMeasurement>* pvMeasurements_)
{
   INT iIndex = 0;
   CMatrix clCl("Observation Covariance Matrix", pvMeasurements_->size(), pvMeasurements_->size());
   // Set the diagonal components of the observation covariance matrix
   for (vector<DistanceMeasurement>::iterator pIterMeas = pvMeasurements_->begin(); pIterMeas != pvMeasurements_->end(); pIterMeas++)
   {
      iIndex = pIterMeas - pvMeasurements_->begin();
      clCl.SetComponent(iIndex, iIndex, pow(pIterMeas->dDistStdDev, 2));
   }

   return clCl;
}

//------------------------------------------------------------------------------
CMatrix HampelWeight(CMatrix* pclObsCov_, CMatrix* pclA_, CMatrix* pclMisc_, DOUBLE* pdConstants_)
{
   // Declare needed matrices
   CMatrix clW("Hampel Matrix", pclObsCov_->GetNumRows(), pclObsCov_->GetNumCols());
   CMatrix clI("Identity", pclObsCov_->GetNumRows(), pclObsCov_->GetNumCols());
   CMatrix clResiduals("Residuals", pclMisc_->GetNumRows(), 1);
   CMatrix clTemp("Temporal", pclObsCov_->GetNumRows(), pclObsCov_->GetNumCols());
   CMatrix clObsCovInv(pclObsCov_->NumericInverse2(1e-6));
   DOUBLE dEntry = 0;

   // Identity matrix
   clI.SetIdentity();
   // Compute residuals
   clTemp = *pclA_ * (pclA_->Transpose() * clObsCovInv * (*pclA_)).NumericInverse2(1e-6) * pclA_->Transpose() * clObsCovInv;
   clResiduals = (clI - clTemp) * (*pclMisc_);

   // Compute Hapel Matrix Weights
   for (INT i = 0; i < clResiduals.GetNumRows(); i++)
   {
      // Hampels first cutoff
      if ( abs(clResiduals.GetComponent(i, 0)) < pdConstants_[LEASTSQUARES_WEIGHT_A_INDEX] )
      {
         dEntry = sqrt(clObsCovInv.GetComponent(i, i));
      }
      // Hampels second cutoff a <= |v| < b
      else if (pdConstants_[LEASTSQUARES_WEIGHT_A_INDEX] <= abs(clResiduals.GetComponent(i, 0)) &&
               pdConstants_[LEASTSQUARES_WEIGHT_B_INDEX] > abs(clResiduals.GetComponent(i, 0)))
      {
         dEntry = sqrt(clObsCovInv.GetComponent(i, i)) * pdConstants_[LEASTSQUARES_WEIGHT_A_INDEX] / abs(clResiduals.GetComponent(i, 0));
      }
      // Hampels third cutoff b <= |v| < c
      else if (pdConstants_[LEASTSQUARES_WEIGHT_B_INDEX] <= abs(clResiduals.GetComponent(i, 0)) &&
               pdConstants_[LEASTSQUARES_WEIGHT_C_INDEX] > abs(clResiduals.GetComponent(i, 0)))
      {
         dEntry = sqrt(clObsCovInv.GetComponent(i, i)) * pdConstants_[LEASTSQUARES_WEIGHT_C_INDEX] / abs(clResiduals.GetComponent(i, 0)) - sqrt(clObsCovInv.GetComponent(i, i));
      }
      // Hampels fourth cutoff |v| >= c
      else
      {
         dEntry = 0;
      }
      clW.SetComponent(i, i, dEntry);
   }
   return clW * clW.Transpose();
}

//------------------------------------------------------------------------------
CMatrix HuberWeight(CMatrix* pclObsCov_, CMatrix* pclA_, CMatrix* pclMisc_, DOUBLE* pdConstants_)
{
   // Declare needed matrices
   CMatrix clW("Huber Matrix", pclObsCov_->GetNumRows(), pclObsCov_->GetNumCols());
   CMatrix clI("Identity", pclObsCov_->GetNumRows(), pclObsCov_->GetNumCols());
   CMatrix clResiduals("Residuals", pclMisc_->GetNumRows(), 1);
   CMatrix clTemp("Temporal", pclObsCov_->GetNumRows(), pclObsCov_->GetNumCols());
   CMatrix clObsCovInv(pclObsCov_->NumericInverse2(1e-6));
   DOUBLE dEntry = 0;

   // Identity matrix
   clI.SetIdentity();
   // Compute residuals
   clTemp = *pclA_ * (pclA_->Transpose() * clObsCovInv * (*pclA_)).NumericInverse2(1e-6) * pclA_->Transpose() * clObsCovInv;
   clResiduals = (clI - clTemp) * (*pclMisc_);

   // Compute Huber Matrix Weights
   for (INT i = 0; i < clResiduals.GetNumRows(); i++)
   {
      // Huber first cutoff
      if (abs(clResiduals.GetComponent(i, 0)) < pdConstants_[LEASTSQUARES_WEIGHT_C_INDEX])
      {
         dEntry = sqrt(clObsCovInv.GetComponent(i, i));
      }
      // Huber second cutoff c <= |v|
      else if (pdConstants_[LEASTSQUARES_WEIGHT_C_INDEX] <= abs(clResiduals.GetComponent(i, 0)))
      {
         dEntry = sqrt(clObsCovInv.GetComponent(i, i)) * pdConstants_[LEASTSQUARES_WEIGHT_C_INDEX] / abs(clResiduals.GetComponent(i, 0));
      }
      clW.SetComponent(i, i, dEntry);
   }
   return clW * clW.Transpose();
}

//------------------------------------------------------------------------------
CMatrix AndrewWeight(CMatrix* pclObsCov_, CMatrix* pclA_, CMatrix* pclMisc_, DOUBLE* pdConstants_)
{
   // Declare needed matrices
   CMatrix clW("Andrew Matrix", pclObsCov_->GetNumRows(), pclObsCov_->GetNumCols());
   CMatrix clI("Identity", pclObsCov_->GetNumRows(), pclObsCov_->GetNumCols());
   CMatrix clResiduals("Residuals", pclMisc_->GetNumRows(), 1);
   CMatrix clTemp("Temporal", pclObsCov_->GetNumRows(), pclObsCov_->GetNumCols());
   CMatrix clObsCovInv(pclObsCov_->NumericInverse2(1e-6));
   DOUBLE dEntry = 0;

   // Identity matrix
   clI.SetIdentity();
   // Compute residuals
   clTemp = *pclA_ * (pclA_->Transpose() * clObsCovInv * (*pclA_)).NumericInverse2() * pclA_->Transpose() * clObsCovInv;
   clResiduals = (clI - clTemp) * (*pclMisc_);

   // Compute Andrew Matrix Weights
   for (INT i = 0; i < clResiduals.GetNumRows(); i++)
   {
      // Andrew first cutoff
      if (abs(clResiduals.GetComponent(i, 0) * sqrt(clObsCovInv.GetComponent(i, i))) < DATATYPEDEFINITIONS_PI * pdConstants_[LEASTSQUARES_WEIGHT_C_INDEX])
      {
         dEntry = sin(sqrt(clObsCovInv.GetComponent(i, i)) * clResiduals.GetComponent(i, 0) / pdConstants_[LEASTSQUARES_WEIGHT_C_INDEX]);
         dEntry = dEntry / (sqrt(clObsCovInv.GetComponent(i, i)) * clResiduals.GetComponent(i, 0) / pdConstants_[LEASTSQUARES_WEIGHT_C_INDEX]);
      }
      // Andrew second cutoff c <= |v|
      else
      {
         dEntry = 0;
      }
      clW.SetComponent(i, i, dEntry);
   }
   return clW * clW.Transpose();
}

//------------------------------------------------------------------------------
CMatrix LSWeight(CMatrix* pclObsCov_, CMatrix* pclA_, CMatrix* pclMisc_, DOUBLE* pdConstants_)
{
   CMatrix clW(pclObsCov_->NumericInverse2(1e-6));
   return clW;
}

//------------------------------------------------------------------------------
CMatrix ComputeStateCovariance(CMatrix* pclWeightMatrix_, CMatrix* pclDesignMatrixA_)
{
   // Normal matrix (A^T * Cl^(-1) * A)
   CMatrix clN(pclDesignMatrixA_->Transpose() * (*pclWeightMatrix_) * (*pclDesignMatrixA_));
   // State covariance matrix (A^T * W^(-1) * A) ^ (-1)
   CMatrix clStateCov(clN.NumericInverse2());

   return clStateCov;
}

//------------------------------------------------------------------------------
CMatrix ComputeStateCorrection(CMatrix* pclWeightMatrix_, CMatrix* pclDesignMatrixA_, CMatrix* pclMisclosure_)
{
   // Zero matrix
   CMatrix clZero("Zero Matrix", pclDesignMatrixA_->GetNumCols(), 1);
   clZero.SetZero();
   // Normal matrix
   CMatrix clN(pclDesignMatrixA_->Transpose() * (*pclWeightMatrix_) * (*pclDesignMatrixA_));
   // Correction vector
   CMatrix clDelta(clN.NumericInverse2() * pclDesignMatrixA_->Transpose() * (*pclWeightMatrix_) * (*pclMisclosure_));
   clDelta = clZero - clDelta;
   return clDelta;
}

//------------------------------------------------------------------------------
void Update(map<ULONG, PointCoordinate>* pmpCoordinates_, CMatrix* pclDelta_)
{
   vector<ULONG> vPoints;
   FixedPoints(pmpCoordinates_, &vPoints);
   // Rows of correction vector
   INT iRowX = -1;
   INT iRowY = -1;
   DOUBLE dCorrectionX = 0;
   DOUBLE dCorrectionY = 0;

   // Point ID iteration
   for (vector<ULONG>::iterator pIterPoint = vPoints.begin(); pIterPoint != vPoints.end(); pIterPoint++)
   {
      iRowX = 2 * (pIterPoint - vPoints.begin());
      iRowY = 2 * (pIterPoint - vPoints.begin()) + 1;
      dCorrectionX = pclDelta_->GetComponent(iRowX, 0);
      dCorrectionY = pclDelta_->GetComponent(iRowY, 0);
      // Update point of expansion
      pmpCoordinates_->at(*pIterPoint).P_X = pmpCoordinates_->at(*pIterPoint).P_X + dCorrectionX;
      pmpCoordinates_->at(*pIterPoint).P_Y = pmpCoordinates_->at(*pIterPoint).P_Y + dCorrectionY;
   }
}

//------------------------------------------------------------------------------
LeastSquares Process(map<ULONG, PointCoordinate>* pmpCoordinates_, vector<DistanceMeasurement>* pvMeasurements_, string sMethod_, DOUBLE* pdConstants_)
{
   LeastSquares clLS;
   // Weight matrix selection
   map<string, RobustWeightMatrix> mpWeightFunc;
   mpWeightFunc["LEASTSQUARES"] = &LSWeight;
   mpWeightFunc["HAMPELWEIGHT"] = &HampelWeight;
   mpWeightFunc["HUBERWEIGHT"] = &HuberWeight;
   mpWeightFunc["ANDREWWEIGHT"] = &AndrewWeight;

   // Misclosures
   CMatrix clMisc(Misclosures(pmpCoordinates_, pvMeasurements_));
   // Geometry matrix
   CMatrix clA(MeasurementMatrix(pmpCoordinates_, pvMeasurements_));
   // Observation Covariance matrix
   CMatrix clCL(ObservationCovariance(pvMeasurements_));
   // Weight Matrix
   CMatrix clW(mpWeightFunc.at(sMethod_)(&clCL, &clA, &clMisc, pdConstants_));
   // State covariance matrix
   CMatrix clCX(ComputeStateCovariance(&clW, &clA));
   // Error correction vector
   CMatrix clDelta(ComputeStateCorrection(&clW, &clA, &clMisc));
   // Correction vector norm
   DOUBLE dErrorNorm2 = clDelta.VectorNorm2();
   // Verify sanity of norm
   if (dErrorNorm2 < 0)
   {
      throw 7;
   }
   // Residuals computation
   CMatrix clI("Identity", clCL.GetNumRows(), clCL.GetNumCols());
   CMatrix clResiduals("Residuals", clMisc.GetNumRows(), 1);
   CMatrix clTemp("Temporal", clCL.GetNumRows(), clCL.GetNumCols());
   CMatrix clObsCovInv(clCL.NumericInverse2(1e-6));

   while (dErrorNorm2 > clCX.MatrixNorm2() * LEASTSQUARES_ERROR_TOLERANCE)
   {
      // Update point of expansion
      Update(pmpCoordinates_, &clDelta);
      // Misclosures
      clMisc = Misclosures(pmpCoordinates_, pvMeasurements_);
      // Geometry Matrix
      clA = MeasurementMatrix(pmpCoordinates_, pvMeasurements_);
      // Weight Matrix
      clW = mpWeightFunc.at(sMethod_)(&clCL, &clA, &clMisc, pdConstants_);
      // State covariance matrix
      clCX = ComputeStateCovariance(&clW, &clA);
      // Error correction vector
      clDelta = ComputeStateCorrection(&clW, &clA, &clMisc);
      // Correction vector norm
      dErrorNorm2 = clDelta.VectorNorm2();
      // Verify sanity of norm
      if (dErrorNorm2 < 0)
      {
         throw 8;
      }
   }
   // State covariance matrix
   clLS.SetStateCovariance(clCX);
   // Residuals vector r = [I - A(A^T * Cl^-1 * A)^-1 * A^T * Cl^-1] * w
   clI.SetIdentity();
   clTemp = clA * (clA.Transpose() * clObsCovInv * clA).NumericInverse2(1e-6) * clA.Transpose() * clObsCovInv;
   clResiduals = (clI - clTemp) * clMisc;
   clLS.SetResiduals(clResiduals);

   return clLS;
}
