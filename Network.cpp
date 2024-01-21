// Network.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

// P1 Coordinates
DOUBLE P1_X = 2000.10;
DOUBLE P1_Y = 4999.90;
// P2 Coordinates
DOUBLE P2_X = 3000.10;
DOUBLE P2_Y = 3000.10;
// P3 Coordinates
DOUBLE P3_X = 5999.85;
DOUBLE P3_Y = 3999.90;
// P4 Coordinates
DOUBLE P4_X = 7999.90;
DOUBLE P4_Y = 6000.15;
// P5 Coordinates
DOUBLE P5_X = 1000.00;
DOUBLE P5_Y = 2000.00;
// P6 Coordinates
DOUBLE P6_X = 7000.00;
DOUBLE P6_Y = 1000.00;
// Influence function constants
DOUBLE da = 0.449;
DOUBLE db = 0.898;
DOUBLE dc = 1.347;
// Input file path
string sInputFile("C:\\Users\\Juanchelsea\\Desktop\\MASTERS\\ENGO629\\Labs\\Lab2\\Input1.txt");
string sOutPath("C:\\Users\\Juanchelsea\\Desktop\\MASTERS\\ENGO629\\Labs\\Lab2");
// Fix points
string sFixPoints("5,6");
// Solution method
string sMethod("LEASTSQUARES");

//-----------------------------------------------------------------------------------
// Control switch options
// Control switches are:
// /-INPUT             : Input measurement file
// /-P1X               : Change X coordinate of P1
// /-P1Y               : Change Y coordinate of P1
// /-P2X               : Change X coordinate of P2
// /-P2Y               : Change Y coordinate of P2
// /-P3X               : Change X coordinate of P3
// /-P3Y               : Change Y coordinate of P3
// /-P4X               : Change X coordinate of P4
// /-P4Y               : Change Y coordinate of P4
// /-P5X               : Change X coordinate of P5
// /-P5Y               : Change Y coordinate of P5
// /-P6X               : Change X coordinate of P6
// /-P6Y               : Change Y coordinate of P6
// /-FIXPOINT          : Change the fix points
// /-RESA              : Residual 'a' threshold
// /-RESB              : Residual 'b' threshold
// /-RESC              : Residual 'c' threshold
//-----------------------------------------------------------------------------------
const CHAR* acOptions[NETWORK_NUM_OPTIONS] = { "INPUT", "P1X", "P1Y", "P2X", "P2Y",
                                               "P3X", "P3Y", "P4X", "P4Y", "P5X",
                                               "P5Y", "P6X", "P6Y", "FIXPOINT", "RESA",
                                               "RESB", "RESC", "OUTPATH", "METHOD"};
vector<string> vOptions(acOptions, acOptions + NETWORK_NUM_OPTIONS);
map<string, void*> mpTheSwitchOptions;
// Error codes
map<INT, string> mpTheErrorCodes;

//-----------------------------------------------------------------------------------
INT main(INT iargc_, CHAR** ppcargv_)
{
   // Point coordinate map
   map<ULONG, PointCoordinate> mpCoordinates;
   // Point measurements
   vector<DistanceMeasurement> vMeasurements;
   // Constants for robust weight least squares methods
   DOUBLE adConstants[3] = {0};
   // Class that holds state covariance and residuals
   LeastSquares LS;
   // Initialize control switch options
   mpTheSwitchOptions[vOptions.at(NETWORK_INPUTFILE_INDEX)] = &sInputFile;
   mpTheSwitchOptions[vOptions.at(NETWORK_P1X_INDEX)] = &P1_X;
   mpTheSwitchOptions[vOptions.at(NETWORK_P1Y_INDEX)] = &P1_Y;
   mpTheSwitchOptions[vOptions.at(NETWORK_P2X_INDEX)] = &P2_X;
   mpTheSwitchOptions[vOptions.at(NETWORK_P2Y_INDEX)] = &P2_Y;
   mpTheSwitchOptions[vOptions.at(NETWORK_P3X_INDEX)] = &P3_X;
   mpTheSwitchOptions[vOptions.at(NETWORK_P3Y_INDEX)] = &P3_Y;
   mpTheSwitchOptions[vOptions.at(NETWORK_P4X_INDEX)] = &P4_X;
   mpTheSwitchOptions[vOptions.at(NETWORK_P4Y_INDEX)] = &P4_Y;
   mpTheSwitchOptions[vOptions.at(NETWORK_P5X_INDEX)] = &P5_X;
   mpTheSwitchOptions[vOptions.at(NETWORK_P5Y_INDEX)] = &P5_Y;
   mpTheSwitchOptions[vOptions.at(NETWORK_P6X_INDEX)] = &P6_X;
   mpTheSwitchOptions[vOptions.at(NETWORK_P6Y_INDEX)] = &P6_Y;
   mpTheSwitchOptions[vOptions.at(NETWORK_FIXPOINTS_INDEX)] = &sFixPoints;
   mpTheSwitchOptions[vOptions.at(NETWORK_RESA_INDEX)] = &da;
   mpTheSwitchOptions[vOptions.at(NETWORK_RESB_INDEX)] = &db;
   mpTheSwitchOptions[vOptions.at(NETWORK_RESC_INDEX)] = &dc;
   mpTheSwitchOptions[vOptions.at(NETWORK_OUTPATH_INDEX)] = &sOutPath;
   mpTheSwitchOptions[vOptions.at(NETWORK_METHOD_INDEX)] = &sMethod;
   // Error codes
   mpTheErrorCodes[1] = "Option used is invalid (main)";
   mpTheErrorCodes[2] = "Switch option used is invalid (main)";
   mpTheErrorCodes[3] = "Parameter out of range (FixedPointsVector)";
   mpTheErrorCodes[4] = "Argument is invalid (FixedPointsVector)";
   mpTheErrorCodes[5] = "Point does not exist (AssignCoordinates)";
   mpTheErrorCodes[6] = "Input file format is incorrect (ReadInputFile)";
   mpTheErrorCodes[7] = "Invalid vector dimensions for norm computation (Process)";
   mpTheErrorCodes[8] = "Invalid vector dimensions for norm computation (Process)";

   try
   {
      // Parse program command arguments
      for (INT i = 1; i < iargc_; i++)
      {
         // Transform one string to uppercase and keep original string
         string stOption(ppcargv_[i]);
         string stOptionOrig(ppcargv_[i]);
         BOOLEANO bOptionFound = FALSE;
         transform(stOption.begin(), stOption.end(), stOption.begin(), toupper);
         switch (stOption.at(0))
         {
            // Allowed switch delimiters
         case '/':
         case '-':
            for (vector<string>::iterator itOptions = vOptions.begin(); itOptions != vOptions.end(); itOptions++)
            {
               // Option found
               if (stOption.find(*itOptions) != string::npos)
               {
                  // File path assignment changes from default
                  if (*itOptions == vOptions.at(NETWORK_INPUTFILE_INDEX) ||
                      *itOptions == vOptions.at(NETWORK_FIXPOINTS_INDEX) ||
                      *itOptions == vOptions.at(NETWORK_OUTPATH_INDEX) ||
                      *itOptions == vOptions.at(NETWORK_METHOD_INDEX))
                  {
                     string stChangeValue = stOptionOrig.substr(stOption.find(*itOptions) + itOptions->length());
                     *static_cast<string*>(mpTheSwitchOptions.at(*itOptions)) = stChangeValue;
                  }
                  else
                  {
                     string stChangeValue = stOptionOrig.substr(stOption.find(*itOptions) + itOptions->length());
                     *((DOUBLE*)mpTheSwitchOptions.at(*itOptions)) = stod(stChangeValue);
                  }
                  bOptionFound = TRUE;
                  break;
               }
            }
            // Option not found
            if (!bOptionFound)
            {
               throw 1;
            }
            break;
            // Invalid switch delimiter
         default:
            throw 2;
         }
      }
      // Point coordinates
      AssignCoordinates(&mpCoordinates, sFixPoints, &adConstants[0]);
      // Get measurements
      ReadInputFile(sInputFile, &vMeasurements);
      // Compute results
      LS = Process(&mpCoordinates, &vMeasurements, sMethod, &adConstants[0]);
      // Output results
      Output(&mpCoordinates, &LS, sOutPath);
      // Output residuals
      Residuals(&vMeasurements, &LS, sOutPath);
   }
   // Catch errors discovered in the program
   catch (INT iError)
   {
      cout << mpTheErrorCodes.at(iError) << endl;
      return EXIT_FAILURE;
   }
   // Catch unknown errors discovered in sub-functions and dependencies
   catch (...)
   {
      cout << "Unknown Error" << endl;
      return EXIT_FAILURE;
   }
   return EXIT_SUCCESS;
}

//----------------------------------------------------------------------------
void AssignCoordinates(map<ULONG, PointCoordinate>* pmpCoordinates_, string sFixedPoints_, DOUBLE* pdConstants_)
{
   pdConstants_[0] = da;
   pdConstants_[1] = db;
   pdConstants_[2] = dc;
   // Fix point coordiantes vector
   vector<ULONG> vFixPoints;
   FixedPointsVector(sFixedPoints_, &vFixPoints);
   // Store point 1 coordinates
   pmpCoordinates_->operator[](1).P_X = P1_X;
   pmpCoordinates_->operator[](1).P_Y = P1_Y;
   pmpCoordinates_->operator[](1).bFIXED = FALSE;
   // Store point 2 coordinates
   pmpCoordinates_->operator[](2).P_X = P2_X;
   pmpCoordinates_->operator[](2).P_Y = P2_Y;
   pmpCoordinates_->operator[](2).bFIXED = FALSE;
   // Store point 3 coordinates
   pmpCoordinates_->operator[](3).P_X = P3_X;
   pmpCoordinates_->operator[](3).P_Y = P3_Y;
   pmpCoordinates_->operator[](3).bFIXED = FALSE;
   // Store point 4 coordinates
   pmpCoordinates_->operator[](4).P_X = P4_X;
   pmpCoordinates_->operator[](4).P_Y = P4_Y;
   pmpCoordinates_->operator[](4).bFIXED = FALSE;
   // Store point 5 coordinates
   pmpCoordinates_->operator[](5).P_X = P5_X;
   pmpCoordinates_->operator[](5).P_Y = P5_Y;
   pmpCoordinates_->operator[](5).bFIXED = FALSE;
   // Store point 6 coordinates
   pmpCoordinates_->operator[](6).P_X = P6_X;
   pmpCoordinates_->operator[](6).P_Y = P6_Y;
   pmpCoordinates_->operator[](6).bFIXED = FALSE;
   // Assign points that are fixed coordinates
   for (vector<ULONG>::iterator pIterFixed = vFixPoints.begin(); pIterFixed != vFixPoints.end(); pIterFixed++)
   {
      // Make sure the point is in the list of point coordiantes
      if (pmpCoordinates_->find(*pIterFixed) != pmpCoordinates_->end())
      {
         pmpCoordinates_->at(*pIterFixed).bFIXED = TRUE;
      }
      else
      {
         throw 5;
      }
   }
}

//----------------------------------------------------------------------------
void FixedPointsVector(const string& sFixedPoints_, vector<ULONG>* pvFixedPoints_)
{
   stringstream ss(sFixedPoints_);
   string sData;
   // Parse line
   while (getline(ss, sData, ','))
   {
      if (!sData.empty())
      {
         try
         {
            // Fix point assignment
            pvFixedPoints_->push_back(stoul(sData));
         }
         // Out of range data
         catch (std::out_of_range ex)
         {
            throw 3;
         }
         // Invalid argument
         catch (std::invalid_argument ex)
         {
            throw 4;
         }
      }
   }
}

//----------------------------------------------------------------------------
void ReadInputFile(string sInputFile_, vector<DistanceMeasurement>* pvMeasurements_)
{
   // File handler
   fstream fInputFile;
   // Open the files
   fInputFile.open(sInputFile_, ios::in | ios::binary);
   string sMeasLine;
   string sData;
   INT iIndex;
   // Reading input file
   while (getline(fInputFile, sMeasLine, '\n'))
   {
      stringstream ss(sMeasLine);
      iIndex = 0;
      // Structure data
      DistanceMeasurement stMeasurement;
      // Parse line
      while (getline(ss, sData, ','))
      {
         if ( !sData.empty() && bParseData(sData) )
         {
            switch (iIndex)
            {
               // Point A
            case 0:
               stMeasurement.ulPointA = stoul(sData);
               break;
               // Point B
            case 1:
               stMeasurement.ulPointB = stoul(sData);
               break;
               // Distance measurement
            case 2:
               stMeasurement.dDistMeas = stod(sData);
               break;
               // Measurement standard deviation
            case 3:
               stMeasurement.dDistStdDev = stod(sData);
               pvMeasurements_->push_back(stMeasurement);
               break;
            default:
               throw 6;
            }
            iIndex = iIndex + 1;
         }
      }
   }
   fInputFile.close();
}

//----------------------------------------------------------------------------
BOOLEANO bParseData(string sData_)
{
   INT iOffset = sData_.find_first_not_of("0123456789");
   BOOLEANO bNumeric = TRUE;
   // Determine if the string is numeric
   if (iOffset != string::npos)
   {
      if (iOffset == 0 && sData_.at(iOffset) != '+' && sData_.at(iOffset) != '-')
      {
         bNumeric = FALSE;
      }
      else if (sData_.at(iOffset) == '.' && sData_.substr(iOffset + 1).find('.') != string::npos)
      {
         bNumeric = FALSE;
      }
   }
   return bNumeric;
}

//----------------------------------------------------------------------------
void Output(map<ULONG, PointCoordinate>* pmpCoordinates_, LeastSquares* pclLS_, string sOutPath_)
{
   vector<ULONG> vPoints;
   INT iColX = -1;
   INT iColY = -1;
   string sOutFile(sOutPath_);
   sOutFile.append("\\Output.txt");
   CreateDirectory(sOutPath_.c_str(), NULL);
   // File handler
   fstream fOutputFile;
   FixedPoints(pmpCoordinates_, &vPoints);
   // Open the files
   fOutputFile.open(sOutFile, ios::out | ios::binary);
   // File header
   fOutputFile << "Point ID, X, Y, SigmaX, SigmaY" << endl;
   // Print computed point coordinates
   for (vector<ULONG>::iterator pIterPoint = vPoints.begin(); pIterPoint != vPoints.end(); pIterPoint++)
   {
      iColX = 2 * (pIterPoint - vPoints.begin());
      iColY = 2 * (pIterPoint - vPoints.begin()) + 1;
      fOutputFile << *pIterPoint << ",";
      fOutputFile << fixed << setprecision(6) << pmpCoordinates_->at(*pIterPoint).P_X << ",";
      fOutputFile << fixed << setprecision(6) << pmpCoordinates_->at(*pIterPoint).P_Y << ",";
      fOutputFile << fixed << setprecision(3) << sqrt(pclLS_->GetStateCovariance().GetComponent(iColX, iColX)) << ",";
      fOutputFile << fixed << setprecision(3) << sqrt(pclLS_->GetStateCovariance().GetComponent(iColY, iColY)) << endl;
   }
   // Print fixed point coordinates
   for (map<ULONG, PointCoordinate>::iterator pIterPoint = pmpCoordinates_->begin(); pIterPoint != pmpCoordinates_->end(); pIterPoint++)
   {
      if (pIterPoint->second.bFIXED)
      {
         fOutputFile << pIterPoint->first << ",";
         fOutputFile << fixed << setprecision(6) << pIterPoint->second.P_X << ",";
         fOutputFile << fixed << setprecision(6) << pIterPoint->second.P_Y << ",";
         fOutputFile << "0.000,";
         fOutputFile << "0.000" << endl;
      }
   }
}

//----------------------------------------------------------------------------
void Residuals(vector<DistanceMeasurement>* pvMeasurements_, LeastSquares* pclLS_, string sOutPath_)
{
   INT iIndex = 0;
   string sOutFile(sOutPath_);
   sOutFile.append("\\Residuals.txt");
   // File handler
   fstream fResidualFile;
   // Open the files
   fResidualFile.open(sOutFile, ios::out | ios::binary);
   // File header
   fResidualFile << "Point A, Point B, Measurement, Residuals" << endl;
   // Print computed residuals
   for (vector<DistanceMeasurement>::iterator pIterMeas = pvMeasurements_->begin(); pIterMeas != pvMeasurements_->end(); pIterMeas++)
   {
      iIndex = pIterMeas - pvMeasurements_->begin();
      fResidualFile << pIterMeas->ulPointA << ",";
      fResidualFile << pIterMeas->ulPointB << ",";
      fResidualFile << fixed << setprecision(6) << pIterMeas->dDistMeas << ",";
      fResidualFile << fixed << setprecision(6) << pclLS_->GetResiduals().GetComponent(iIndex, 0) << endl;
   }
}
