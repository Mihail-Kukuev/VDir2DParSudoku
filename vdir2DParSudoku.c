#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "mpi.h"
#include "task_def.h"

#define MAX(a,b) ((a)>(b)?(a):(b))

//#define WITH_PROCESS_LOG

void swap_pointers_d(double** a, double** b);

void swap_pointers(double*** a, double*** b);

/* Example of 4 processes grid:
 * x2 axis
 * ^________________
 * | 3 | 0 | 1 | 2 |
 * | 2 | 3 | 0 | 1 |
 * | 1 | 2 | 3 | 0 |
 * |_0_|_1_|_2_|_3_|-> x1 axis
 */

int main(int argc, char **argv )
{
    int i, i1, i2, j, i1_gl, i2_gl;
    int n1 = 0;
    int n2 = 0;
    int j0 = 0;
    double h1 = 0.0;
    double h2 = 0.0;
    double tau = 0.0;
    double h1Square = 0.0;      // h1*h1
    double h2Square = 0.0;      // h2*h2
    double h1SquareTau = 0.0;   // 2*h1*h1/tau
    double h2SquareTau = 0.0;   // 2*h2*h2/tau
    double h1_2Square = 0.0;    // (h1*h1)/(h2*h2)
    double h2_1Square = 0.0;    // (h2*h2)/(h1*h1)
    double x1;
    double x2;
    double x1start = 0.0;
    double x2start = 0.0;
    double t;
    double tPrev;
    double tPrevSub = 0.0;
    double *pX1start = NULL;
    double *pX2start = NULL;

    double dAValue = 0.0;
    double dBValue = 0.0;
    double dCValue = 0.0;
    double dFValue = 0.0;
    double kappa1;
    double eta1;
    double kappa2;
    double eta2;
    double *dAlphaBeta = NULL;
    double *pABStartPos = NULL;
    double *pABLine = NULL;
    double *pABInLinePos = NULL;
    double *pABLastLine = NULL;
    double dCoefA1_R;
    double dCoefA1_L;
    double dCoefA2_R;
    double dCoefA2_L;
    double dK;
    double dDiv;

    double *dY = NULL;
    double *yCurrent = NULL;
    double *yPrev = NULL;
    double *pYPrevStartPos = NULL;
    double *pYPrevLine = NULL;
    double *pYCurrentStartPos = NULL;
    double *pYCurrentLine = NULL;
    double *pYLastLine = NULL;
    double* pYLastCol = NULL;
    double *p = NULL;
    int iYBufferSize = 0;
    double eps = 0.0;
    FILE *F = NULL;

    int iMyRank = 0;
    int iPCount = 0;
    int iPrevRank = 0;
    int iNextRank = 0;
    char cStr[10];
    int iCurTile1_Size = 0;
    int iCurTile2_Size = 0;
    int iBufferABTileSize = 0;
    int *iTile2_Size = NULL; // Distribution of mesh line counts per processes to find solution
    // in direction x1. Mesh line fixed by value of coordinate i2.
    int *iTile1_Size = NULL; // Distribution of mesh line counts per processes to find solution in direction x2
    // Mesh line fixed by value of coordinate i1.
    int tag = 0;
    MPI_Status status;
    MPI_Request recvRequest = MPI_REQUEST_NULL;
    MPI_Request sendRequest = MPI_REQUEST_NULL;
    FILE *FL = NULL;

    clock_t cp1, cp2;
    double time1, time2;


    // Initialization of MPI
    MPI_Init ( &argc, &argv );
    MPI_Comm_rank ( MPI_COMM_WORLD, &iMyRank );
    MPI_Comm_size ( MPI_COMM_WORLD, &iPCount );

    iPrevRank = (iMyRank - 1 + iPCount) % iPCount;
    iNextRank = (iMyRank + 1) % iPCount;

    if ( iMyRank == 0 )
    {
        printf ( "\nSolving the heat equation by variable directions method\n" );
        printf ( "Parallel mode, version 1.0\n" );
        printf ( "Process count: %d\n", iPCount );
    }

#ifdef WITH_PROCESS_LOG
    // Open log-file for current process
 sprintf ( cStr, "log_%d", iMyRank );
 FL = fopen ( ( const char* ) cStr, "w" );
 fprintf ( FL, "Process: %d\nProcess count: %d\n", iMyRank, iPCount );
#endif


    // Load mesh parameters n1, n2, j0 from command line
    n1 = atoi ( argv[1] );
    n2 = atoi ( argv[2] );
    j0 = atoi ( argv[3] );

    // Setting of discretization steps and its derivatives
    h1 = ( ( double ) L1 ) / n1;
    h2 = ( ( double ) L2 ) / n2;
    tau = ( ( double ) T ) / j0;
    h1Square = h1 * h1;
    h1SquareTau = 2 * h1Square / tau;
    h2Square = h2 * h2;
    h2SquareTau = 2 * h2Square / tau;
    h1_2Square = h1Square / h2Square;
    h2_1Square = h2Square / h1Square;

    if ( iMyRank == 0 )
    {
        printf ( "Mesh parameters:\n  n1=%d  n2=%d  j0=%d\n  h1=%e  h2=%e  tau=%e\n",
                 n1, n2, j0, h1, h2, tau );
    }

#ifdef WITH_PROCESS_LOG
    fprintf ( FL, "Mesh parameters:\n  n1=%d  n2=%d  j0=%d\n  h1=%e  h2=%e  tau=%e\n",
           n1, n2, j0, h1, h2, tau );
#endif


    // Calculate the count of mesh lines for each process in both direction x1 and x2.
    // Tile size is count of such mesh lines without lines with boundary values.
    iTile2_Size = ( int* ) calloc ( iPCount, sizeof ( int ) );
    iTile1_Size = ( int* ) calloc ( iPCount, sizeof ( int ) );
    for ( i = 0; i < iPCount; i++ )
    {
        iTile2_Size[i] = ( n2 - 1 ) / iPCount + ( i < ( n2 - 1 ) % iPCount );
        iTile1_Size[i] = ( n1 - 1 ) / iPCount + ( i < ( n1 - 1 ) % iPCount );
    }

    // Calculate start values of variables x1 and x2 for each block in "sudoku" grid, and start positions of solutions
    pX1start = ( double* ) calloc ( iPCount, sizeof ( double ) );
    pX2start = ( double* ) calloc ( iPCount, sizeof ( double ) );
    for ( i = 0, x1start = h1, x2start = h2; i < iPCount; i++ )
    {
        pX1start[i] = x1start;
        for (j = 0; j < iTile1_Size[i]; j++)
            x1start += h1;
        pX2start[i] = x2start;
        for (j = 0; j < iTile2_Size[i]; j++)
            x2start += h2;
    }


#ifdef WITH_PROCESS_LOG
    fprintf ( FL, "Count of processed mesh lines: %d (x1), %d (x2)\n",
              iTile2_Size[iMyRank], iTile1_Size[iMyRank] );
#endif

    // For each tile we should have 2 extra columns of solutions, one from previous tile and one from next.
    // They also will be used for boundary values.
    int iX1_YLineSize = n1 - 1 + 2 * iPCount;

    int iMaxTile1_Size = iTile1_Size[0];
    int iMaxTile2_Size = iTile2_Size[0];

    // Memory allocation to store mesh data and coefficients of linear equations
    iYBufferSize = ( iMaxTile2_Size + 2 ) * iX1_YLineSize;
    yCurrent = ( double* ) calloc ( iYBufferSize, sizeof ( double ) );
    yPrev = ( double* ) calloc ( iYBufferSize, sizeof ( double ) );

    iBufferABTileSize = 2 * (iMaxTile1_Size + 1) * (iMaxTile2_Size + 1);
    dAlphaBeta = ( double* ) calloc ( iBufferABTileSize * iPCount, sizeof ( double ) );


    // Pointers to start position of solution in tile
    double** pCurYTileStart = ( double** ) calloc ( iPCount, sizeof ( double* ) );
    double** pPrevYTileStart = ( double** ) calloc ( iPCount, sizeof ( double* ) );
    for ( i = 0; i < iPCount; i++ )
    {
        pCurYTileStart[i] = (i != 0) ? pCurYTileStart[i - 1] + iTile1_Size[i - 1] + 2: yCurrent;
        pPrevYTileStart[i] = (i != 0) ? pPrevYTileStart[i - 1] + iTile1_Size[i - 1] + 2: yPrev;
    }

    // Custom types for sending/receiving coefficients and Y values
    MPI_Datatype *AB_COLUMN_TYPE = ( MPI_Datatype* ) malloc ( iPCount * sizeof ( MPI_Datatype ) );
    MPI_Datatype *Y_COLUMN_TYPE = ( MPI_Datatype* ) malloc ( iPCount * sizeof ( MPI_Datatype ) );
    for ( i1_gl = 0; i1_gl < iPCount; i1_gl++ )
    {
        i2_gl = (iMyRank + iPCount - i1_gl) % iPCount;

        MPI_Type_vector(iTile2_Size[i2_gl], 2, (iTile1_Size[i1_gl] + 1) << 1, MPI_DOUBLE, AB_COLUMN_TYPE + i1_gl);
        MPI_Type_commit(AB_COLUMN_TYPE + i1_gl);

        MPI_Type_vector(iTile2_Size[i2_gl], 1, iX1_YLineSize, MPI_DOUBLE, &Y_COLUMN_TYPE[i1_gl]);
        MPI_Type_commit(&Y_COLUMN_TYPE[i1_gl]);
    }


    // Initialization mesh function on 0-layer and save it as yPrev array
    for (i1_gl = 0; i1_gl < iPCount; i1_gl++)
    {
        pYCurrentStartPos = pCurYTileStart[i1_gl];
        i2_gl = (iMyRank + iPCount - i1_gl) % iPCount;

        for (i2 = 0, x2 = pX2start[i2_gl] - h2, pYCurrentLine = pYCurrentStartPos;
             i2 <= iTile2_Size[i2_gl] + 1;
             i2++, x2 += h2, pYCurrentLine += iX1_YLineSize)
        {
            for (i1 = 0, x1 = pX1start[i1_gl] - h1;
                 i1 <= iTile1_Size[i1_gl] + 1;
                 i1++, x1 += h1, p++)
            {
                pYCurrentLine[i1] = u(x1, x2, 0.0);
            }
        }
    }


    // Main cycle
    MPI_Barrier ( MPI_COMM_WORLD );
    cp1 = clock ();
    time1 = MPI_Wtime ();
    tag = 0;

    for ( j = 1, tPrev = 0.0;  j <= j0;  j++, tPrev += tau )
    {
        //---------------------------------------------------------------//
        //   Calculate mesh function on j-1/2 sublayer in direction x1   //
        //---------------------------------------------------------------//
        t = tPrev + tau / 2;

        swap_pointers_d(&yCurrent, &yPrev);
        swap_pointers(&pCurYTileStart, &pPrevYTileStart);
        memset ( yCurrent, 0, iYBufferSize * sizeof ( double ) );

        for (i1_gl = 0, pABStartPos = dAlphaBeta;
             i1_gl < iPCount;
             i1_gl++, pABStartPos += iBufferABTileSize)
        {
            iCurTile1_Size = iTile1_Size[i1_gl];
            pYPrevStartPos = pPrevYTileStart[i1_gl] + iX1_YLineSize;

            i2_gl = (iMyRank + iPCount - i1_gl) % iPCount;
            iCurTile2_Size = iTile2_Size[i2_gl];

            if (i1_gl != 0)
            {
                // Process receives bound values of coefficients
                // alpha and beta from previous process for current tile.
                MPI_Recv(pABStartPos, 1, AB_COLUMN_TYPE[i1_gl], iPrevRank, tag, MPI_COMM_WORLD, &status);
            }
            else {
                // Process calculates all coefficients alpha(1,i2) and beta(1,i2)
                for (i2 = 1, x2 = pX2start[i2_gl], pABInLinePos = pABStartPos;
                     i2 <= iCurTile2_Size;
                     i2++, x2 += h2, pABInLinePos += (iCurTile1_Size + 1) << 1)
                {
                    dCoefA1_L = 0.5 * ( k2 ( 0.0, x2, tPrev ) + k2 ( 0.0, x2 - h2, tPrev ) );
                    dCoefA1_R = 0.5 * ( k2 ( 0.0, x2 + h2, tPrev ) + k2 ( 0.0, x2, tPrev ) );
                    dCoefA2_L = 0.5 * ( k2 ( 0.0, x2, tPrev + tau ) + k2 ( 0.0, x2 - h2, tPrev + tau ) );
                    dCoefA2_R = 0.5 * ( k2 ( 0.0, x2 + h2, tPrev + tau ) + k2 ( 0.0, x2, tPrev + tau ) );
                    eta1 = 0.5 * ( m1 ( x2, tPrev ) + m1 ( x2, tPrev + tau ) ) -
                           ( ( dCoefA2_L * m1 ( x2 - h2, tPrev + tau ) - ( dCoefA2_L + dCoefA2_R ) * m1 ( x2, tPrev + tau ) +
                               dCoefA2_R * m1 ( x2 + h2, tPrev + tau ) ) -
                             ( dCoefA1_L * m1 ( x2 - h2, tPrev ) - ( dCoefA1_L + dCoefA1_R ) * m1 ( x2, tPrev ) +
                               dCoefA1_R * m1 ( x2 + h2, tPrev ) ) ) *
                           0.5 / h2SquareTau;
                    kappa1 = 0.0;
                    pABInLinePos[0] = kappa1;
                    pABInLinePos[1] = eta1;
                }
            }

            // Calculate coefficients alpha[] and beta[] in current tile
            for (i2 = 1, x2 = pX2start[i2_gl], pABLine = pABStartPos, pYPrevLine = pYPrevStartPos;
                 i2 <= iCurTile2_Size;
                 i2++, x2 += h2, pABLine += (iCurTile1_Size + 1) << 1, pYPrevLine += iX1_YLineSize)
            {
                // Scan nodes in line of tile and consiquently calculate coefficients alpha and beta
                // associated with nodes
                x1start = pX1start[i1_gl];

                dCoefA1_L = 0.5 * (k1(x1start, x2, t) + k1(x1start - h1, x2, t));
                for (i1 = 1, x1 = pX1start[i1_gl], pABInLinePos = pABLine;
                     i1 <= iCurTile1_Size;
                     i1++, x1 += h1, pABInLinePos += 2)
                {
                    // Calculate a1(x1+h1,x2,t), a2(x1,x2,t-tau/2) and a2(x1,x2+h2,t-tau/2).
                    // Coefficient a1(x1,x2,t) calculated on previous i1-1 iteration of cycle
                    dCoefA1_R = 0.5 * ( k1 ( x1 + h1, x2, t ) + k1 ( x1, x2, t ) );
                    dK = k2 ( x1, x2, tPrev );
                    dCoefA2_L = 0.5 * ( dK + k2 ( x1, x2 - h2, tPrev ) );
                    dCoefA2_R = 0.5 * ( k2 ( x1, x2 + h2, tPrev ) + dK );

                    // Calculate coefficients A, B, C, F in linear equation
                    dAValue = dCoefA1_L;
                    dCValue = dCoefA1_L + dCoefA1_R + h1SquareTau;
                    dBValue = dCoefA1_R;
                    dFValue = h1_2Square *
                              ( dCoefA2_L * pYPrevLine[i1 - iX1_YLineSize] -
                                ( dCoefA2_L + dCoefA2_R ) * pYPrevLine[i1] +
                                dCoefA2_R * pYPrevLine[i1 + iX1_YLineSize] ) +
                              h1SquareTau * pYPrevLine[i1] + h1Square * f ( x1, x2, tPrev );
                    dCoefA1_L = dCoefA1_R;

                    // Calculate coefficients alpha(x1,x2+h2) and beta(x1,x2+h2)
                    dDiv = dCValue - pABInLinePos[0] * dAValue;
                    pABInLinePos[2] = dBValue / dDiv;
                    pABInLinePos[3] = (dAValue * pABInLinePos[1] + dFValue) / dDiv;
                }
            }

            if (i1_gl != iPCount - 1)
            {
                // Process sends bound values of coefficients alpha and beta to next process
                MPI_Wait(&sendRequest, &status);
                MPI_Isend(pABStartPos + (iCurTile1_Size << 1), 1, AB_COLUMN_TYPE[i1_gl],
                          iNextRank, tag, MPI_COMM_WORLD, &sendRequest);
            }
        }

        MPI_Wait(&sendRequest, &status);

        ////////////////////////////////////////////////////////////////////////////////////
        // Second (reverse) stage of progonka - calculation of solution.

        for (i1_gl = iPCount - 1, pABStartPos = dAlphaBeta + iBufferABTileSize * (iPCount - 1);
             i1_gl >= 0;
             i1_gl--, pABStartPos -= iBufferABTileSize)
        {
            iCurTile1_Size = iTile1_Size[i1_gl];
            pYCurrentStartPos = pCurYTileStart[i1_gl] + iX1_YLineSize;
            pYLastCol = pYCurrentStartPos + iCurTile1_Size + 1;

            i2_gl = (iMyRank + iPCount - i1_gl) % iPCount;
            iCurTile2_Size = iTile2_Size[i2_gl];

            if (i1_gl != iPCount - 1)
            {
                MPI_Recv(pYLastCol, 1, Y_COLUMN_TYPE[i1_gl], iNextRank, tag, MPI_COMM_WORLD, &status );
            }
            else {
                for (i2 = 0, x2 = pX2start[i2_gl], pABInLinePos = pABStartPos + (iCurTile1_Size << 1);
                     i2 < iCurTile2_Size;
                     i2++, x2 += h2, pABInLinePos += (iCurTile1_Size + 1) << 1)
                {
                    dCoefA1_L = 0.5 * ( k2 ( L1, x2, tPrev ) + k2 ( L1, x2 - h2, tPrev ) );
                    dCoefA1_R = 0.5 * ( k2 ( L1, x2 + h2, tPrev ) + k2 ( L1, x2, tPrev ) );
                    dCoefA2_L = 0.5 * ( k2 ( L1, x2, tPrev + tau ) + k2 ( L1, x2 - h2, tPrev + tau ) );
                    dCoefA2_R = 0.5 * ( k2 ( L1, x2 + h2, tPrev + tau ) + k2 ( L1, x2, tPrev + tau ) );
                    eta2 = 0.5 * ( m2 ( x2, tPrev ) + m2 ( x2, tPrev + tau ) ) -
                           ( ( dCoefA2_L * m2 ( x2 - h2, tPrev + tau ) - ( dCoefA2_L + dCoefA2_R ) * m2 ( x2, tPrev + tau ) +
                               dCoefA2_R * m2 ( x2 + h2, tPrev + tau ) ) -
                             ( dCoefA1_L * m2 ( x2 - h2, tPrev ) - ( dCoefA1_L + dCoefA1_R ) * m2 ( x2, tPrev ) +
                               dCoefA1_R * m2 ( x2 + h2, tPrev ) ) ) *
                           0.5 / h2SquareTau;
                    kappa2 = 0.0;
                    pYLastCol[i2 * iX1_YLineSize] = (eta2 + kappa2 * pABInLinePos[1]) / (1 - kappa2 * pABInLinePos[0]);
                }
            }

            // Calculate solution in current tile.
            // Scan lines of tile (line is a set of tile nodes with fixed i2)
            for (i2 = 0, pABLine = pABStartPos + (iCurTile1_Size << 1), pYCurrentLine = pYCurrentStartPos;
                 i2 < iCurTile2_Size;
                 i2++, pABLine += (iCurTile1_Size + 1) << 1, pYCurrentLine += iX1_YLineSize)
            {
                // Scan nodes in line and consiquently calculate solution associated with nodes
                for (i1 = iCurTile1_Size + 1, pABInLinePos = pABLine;
                     i1 > 0;
                     i1--, pABInLinePos -= 2)
                {
                    pYCurrentLine[i1 - 1] = pABInLinePos[0] * pYCurrentLine[i1] + pABInLinePos[1];
                }
            }

            // Process sends bound values of solution to previous process.
            if (i1_gl != 0)
            {
                MPI_Wait(&sendRequest, &status);
                MPI_Isend(pYCurrentStartPos + 1, 1, Y_COLUMN_TYPE[i1_gl], iPrevRank, tag, MPI_COMM_WORLD, &sendRequest);
            }
        }

        MPI_Wait(&sendRequest, &status);

        //--------------------------------------------------------//
        //   Calculate mesh function on j layer in direction x2   //
        //--------------------------------------------------------//

        tPrevSub = t;
        t = tPrev + tau;

        // Shift solution obtained on j-1 step
        swap_pointers_d(&yCurrent, &yPrev);
        swap_pointers(&pCurYTileStart, &pPrevYTileStart);

        ///////////////////////////////////////////////////////////////////////////////////////
        // First (forward) stage of progonka - calculation of coefficients alpha[] and beta[]

        for (i2_gl = 0; i2_gl < iPCount; i2_gl++)
        {
            iCurTile2_Size = iTile2_Size[i2_gl];

            i1_gl = (iMyRank + iPCount - i2_gl) % iPCount;
            iCurTile1_Size = iTile1_Size[i1_gl];
            pYPrevStartPos = pPrevYTileStart[i1_gl] + iX1_YLineSize;
            pABStartPos = dAlphaBeta + i1_gl * iBufferABTileSize;

            if (i2_gl != 0)
            {
                // Process receives bound values of coefficients alpha and beta from previous process.
                MPI_Recv(pABStartPos, iCurTile1_Size << 1, MPI_DOUBLE, iPrevRank, tag, MPI_COMM_WORLD, &status);
            }
            else {
                // Process calculates all coefficients alpha(i1,1) and beta(i1,1)
                for (i1 = 1, x1 = pX1start[i1_gl], pABInLinePos = pABStartPos;
                     i1 <= iCurTile1_Size;
                     i1++, x1 += h1, pABInLinePos += 2)
                {
                    kappa1 = 0.0;
                    eta1 = m3(x1, t);
                    pABInLinePos[0] = kappa1;
                    pABInLinePos[1] = eta1;
                }
            }

            // Calculate coefficients alpha[] and beta[] in current tile
            for (i2 = 1, x2 = pX2start[i2_gl], pABLine = pABStartPos, pYPrevLine = pYPrevStartPos;
                 i2 <= iCurTile2_Size;
                 i2++, x2 += h2, pABLine += iCurTile1_Size << 1, pYPrevLine += iX1_YLineSize)
            {
                // Scan nodes in line of tile and consiquently calculate coefficients alpha and beta
                // associated with nodes
                x1start = pX1start[i1_gl];

                dCoefA1_L = 0.5 * (k1(x1start, x2, tPrevSub) + k1(x1start - h1, x2, tPrevSub));
                for (i1 = 1, x1 = x1start, pABInLinePos = pABLine;
                     i1 <= iCurTile1_Size;
                     i1++, x1 += h1, pABInLinePos += 2)
                {
                    // Mesh node is fixed by current value of (x1,x2).
                    // Calculate a1(x1+h1,x2,t), a2(x1,x2,t-tau/2) and a2(x1,x2+h2,t-tau/2).
                    // Coefficient a1(x1,x2,t) calculated on previous i1-1 iteration of cycle
                    dCoefA1_R = 0.5 * (k1(x1 + h1, x2, tPrevSub) + k1(x1, x2, tPrevSub));
                    dK = k2(x1, x2, t);
                    dCoefA2_L = 0.5 * (dK + k2(x1, x2 - h2, t));
                    dCoefA2_R = 0.5 * (k2(x1, x2 + h2, t) + dK);

                    // Calculate coefficients A, B, C, F in linear equation.
                    dAValue = dCoefA2_L;
                    dCValue = dCoefA2_L + dCoefA2_R + h2SquareTau;
                    dBValue = dCoefA2_R;
                    dFValue = h2_1Square *
                              (dCoefA1_L * pYPrevLine[i1 - 1] -
                               (dCoefA1_L + dCoefA1_R) * pYPrevLine[i1] +
                               dCoefA1_R * pYPrevLine[i1 + 1]) +
                              h2SquareTau * pYPrevLine[i1] + h2Square * f(x1, x2, tPrev);
                    dCoefA1_L = dCoefA1_R;

                    // Calculate coefficients alpha(x1,x2+h2) and beta(x1,x2+h2)
                    dDiv = dCValue - pABInLinePos[0] * dAValue;
                    pABInLinePos[iCurTile1_Size << 1] = dBValue / dDiv;
                    pABInLinePos[(iCurTile1_Size << 1) + 1] = (dAValue * pABInLinePos[1] + dFValue) / dDiv;
                }
            }

            // Process sends bound values of coefficients alpha and beta to next process.
            if (i2_gl != iPCount - 1)
            {
                MPI_Wait(&sendRequest, &status);
                MPI_Isend(pABLine, iCurTile1_Size << 1, MPI_DOUBLE, iNextRank, tag, MPI_COMM_WORLD, &sendRequest);
            }
        }

        MPI_Wait(&sendRequest, &status);

        ////////////////////////////////////////////////////////////////////////////////////
        // Second (reverse) stage of progonka - calculation of solution.

        for (i2_gl = iPCount - 1; i2_gl >= 0; i2_gl--)
        {
            iCurTile2_Size = iTile2_Size[i2_gl];

            i1_gl = (iMyRank + iPCount - i2_gl) % iPCount;
            iCurTile1_Size = iTile1_Size[i1_gl];
            pYCurrentStartPos = pCurYTileStart[i1_gl];
            pABStartPos = dAlphaBeta + i1_gl * iBufferABTileSize;

            pABLastLine = pABStartPos + iCurTile2_Size * (iCurTile1_Size << 1);
            pYLastLine = pYCurrentStartPos + (iCurTile2_Size + 1) * iX1_YLineSize;

            if (i2_gl != iPCount - 1)
            {
                MPI_Recv(pYLastLine + 1, iCurTile1_Size, MPI_DOUBLE, iNextRank, tag, MPI_COMM_WORLD, &status);
            }
            else {
                for (i1 = 1, x1 = pX1start[i1_gl], pABInLinePos = pABLastLine;
                     i1 <= iCurTile1_Size;
                     i1++, x1 += h1, pABInLinePos += 2)
                {
                    kappa2 = 0.0;
                    eta2 = m4(x1, t);
                    pYLastLine[i1] = (eta2 + kappa2 * pABInLinePos[1]) / (1 - kappa2 * pABInLinePos[0]);
                }
            }

            // Calculate solution in current tile.
            // Scan lines of tile (the line is a set of tile nodes with fixed i2)
            for (i2 = iCurTile2_Size, pABLine = pABLastLine, pYCurrentLine = pYLastLine;
                 i2 >= 0;
                 i2--, pABLine -= iCurTile1_Size << 1, pYCurrentLine -= iX1_YLineSize)
            {
                // Scan nodes in line and consiquently calculate solution associated with nodes
                for (i1 = 1, pABInLinePos = pABLine;
                     i1 <= iCurTile1_Size;
                     i1++, pABInLinePos += 2)
                {
                    pYCurrentLine[i1 - iX1_YLineSize] = pABInLinePos[0] * pYCurrentLine[i1] + pABInLinePos[1];
                }
            }

            // Process sends bound values of solution to previous process.
            if (i2_gl != 0)
            {
                MPI_Wait(&sendRequest, &status);
                MPI_Isend(pYCurrentLine + iX1_YLineSize + 1, iCurTile1_Size, MPI_DOUBLE,
                         iPrevRank, tag, MPI_COMM_WORLD, &sendRequest);
            }

            // Initialize boundary values at x1 = 0 and x1 = L1
            if (i1_gl == 0 || i1_gl == iPCount - 1)
            {
                p = yCurrent;
                for (i2 = 0, x2 = pX2start[i2_gl] - h2;
                     i2 <= iCurTile2_Size + 1;
                     i2++, x2 += h2, p += iX1_YLineSize)
                {
                    if (i1_gl == 0)
                        p[0] = m1(x2, t);
                    else
                        p[iX1_YLineSize - 1] = m2(x2, t);
                }
            }
        }

        MPI_Wait(&sendRequest, &status);
        tag++;
    }

    MPI_Barrier ( MPI_COMM_WORLD );
    cp2 = clock ();
    time2 = MPI_Wtime ();

    // Gather mesh solution onto zero process

    if ( iMyRank == 0 )
    {
        dY = (double *) calloc((n1 + 1) * (n2 + 1), sizeof(double));
    }

    double *dYStartPos;
    int iLinesCount, iLineSize;
    for (i2_gl = 0, dYStartPos = dY;
         i2_gl < iPCount;
         i2_gl++, dYStartPos += iLinesCount * (n1 + 1))
    {
        iLinesCount = iTile2_Size[i2_gl];
        if ( i2_gl == iPCount - 1 || i2_gl == 0 )
            iLinesCount++;

        for (i1_gl = 0, p = dYStartPos;
             i1_gl < iPCount;
             i1_gl++, p += iLineSize)
        {
            pYCurrentStartPos = pCurYTileStart[i1_gl];
            iCurTile1_Size = iTile1_Size[i1_gl];
            iLineSize = iCurTile1_Size;
            if ( i1_gl == iPCount - 1 || i1_gl == 0 )
                iLineSize++;

            double* pYCurrentSendStartPos = pYCurrentStartPos;
            if (i2_gl != 0)
                pYCurrentSendStartPos += iX1_YLineSize;
            if (i1_gl != 0)
                pYCurrentSendStartPos++;

            int iCurProc = (i1_gl + i2_gl) % iPCount;

            //send-receive iLinesCount of lines-solutions with length = iLineLength, placed with offset = iX1_YLineSize
            if ( iMyRank == 0 )
            {
                MPI_Datatype blocktype;
                MPI_Type_vector(iLinesCount, iLineSize, n1 + 1, MPI_DOUBLE, &blocktype);
                MPI_Type_commit(&blocktype);
                MPI_Irecv(p, 1, blocktype, iCurProc, tag, MPI_COMM_WORLD, &recvRequest);
            }

            if (iMyRank == iCurProc)
            {
                MPI_Datatype blocktype;
                MPI_Type_vector(iLinesCount, iLineSize, iX1_YLineSize, MPI_DOUBLE, &blocktype);
                MPI_Type_commit(&blocktype);
                MPI_Send(pYCurrentSendStartPos, 1, blocktype, 0, tag, MPI_COMM_WORLD);
            }

            if (iMyRank == 0)
                MPI_Wait(&recvRequest, &status);
        }
    }


// Save mesh solution into file
#ifndef WITHOUT_RES_SAVE
    if ( iMyRank == 0 )
    {
        F = fopen ( "res.parSudoku", "w" );
        p = dY;
        for ( i2 = 0;  i2 <= n2;  i2++ )
        {
            for ( i1 = 0;  i1 <= n1;  i1++, p++ )  fprintf ( F, "%e ", *p );
            fprintf ( F, "\n" );
        }

        fclose ( F );
    }
#endif


    // Calculate error of computation
    if ( iMyRank == 0 )
    {
        p = dY;
        eps = 0.0;
        for ( i2 = 0, x2 = 0.0;  i2 <= n2;  i2++, x2 += h2 )
        {
            for ( i1 = 0, x1 = 0.0;  i1 <= n1;  i1++, x1 += h1, p++ )
            {
                eps = MAX(eps,fabs ( *p - u ( x1, x2, T ) ));
            }
        }
    }


    if ( iMyRank == 0 )
    {
        printf ( "Time of execution (real/process): %.1fsec / %.1fsec\n",
                 time2 - time1, ( ( double ) ( cp2 - cp1 ) ) / CLOCKS_PER_SEC );
        printf ( "Error of computation: %e\n\n", eps );
    }


    if ( dY != NULL ) free ( dY );
    if ( pCurYTileStart != NULL ) free ( pCurYTileStart );
    if ( pPrevYTileStart != NULL ) free ( pPrevYTileStart );
    if ( AB_COLUMN_TYPE != NULL ) free ( AB_COLUMN_TYPE );
    if ( Y_COLUMN_TYPE != NULL ) free ( Y_COLUMN_TYPE );
    if ( yCurrent != NULL ) free ( yCurrent );
    if ( yPrev != NULL ) free ( yPrev );
    if ( dAlphaBeta != NULL ) free ( dAlphaBeta );
    if ( iTile2_Size != NULL ) free ( iTile2_Size );
    if ( iTile1_Size != NULL ) free ( iTile1_Size );
    if ( pX1start != NULL ) free ( pX1start );
    if ( pX2start != NULL ) free ( pX2start );
    if ( FL != NULL ) fclose ( FL );

    MPI_Finalize ();

    return 0;
}

inline void swap_pointers(double*** a, double*** b)
{
    double** p = *b;
    *b = *a;
    *a = p;
}

inline void swap_pointers_d(double** a, double** b)
{
    double* p = *b;
    *b = *a;
    *a = p;
}
