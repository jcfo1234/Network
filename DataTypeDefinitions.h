#ifndef DATATYPEDEFINITIONS_H
#define DATATYPEDEFINITIONS_H

// Data type definitions
typedef double DOUBLE;
typedef float FLOAT;
typedef int INT;
typedef unsigned int UINT;
typedef unsigned long ULONG;
typedef unsigned long long ULONGLONG;
typedef long LONG;
typedef long long LONGLONG;
typedef short SHORT;
typedef unsigned short USHORT;
typedef unsigned char UCHAR;
typedef char CHAR;
typedef bool BOOLEANO;

// Constant definitions
#define DATATYPEDEFINITIONS_PI ((DOUBLE) 3.14159265358979323846)
#ifndef TRUE
#define TRUE true
#endif
#ifndef FALSE
#define FALSE false
#endif

// Coordinates structure
typedef struct PointCoordinate
{
   DOUBLE P_X;
   DOUBLE P_Y;
   BOOLEANO bFIXED;
} PointCoordinate;

// Measurement structure
typedef struct DistanceMeasurement
{
   ULONG ulPointA;
   ULONG ulPointB;
   DOUBLE dDistMeas;
   DOUBLE dDistStdDev;
} DistanceMeasurement;

#endif