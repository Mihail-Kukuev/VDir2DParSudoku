#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "mpi.h"
#include "task_def.h"

#define MAX(a,b) ((a)>(b)?(a):(b))

/* Example of processes grid:
 * x2 axis
 * ^________________
 * | 3 | 0 | 1 | 2 |
 * | 2 | 3 | 0 | 1 |
 * | 1 | 2 | 3 | 0 |
 * |_0_|_1_|_2_|_3_|-> x1 axis
 */

int main(int argc, char **argv )
{
    int i, i1, i2, j;
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
    double *pX1start = NULL;
    double *pX2start = NULL;
    double t;
    double tPrev;
    double tPrevSub = 0.0;

    double dAValue = 0.0;
    double dBValue = 0.0;
    double dCValue = 0.0;
    double dFValue = 0.0;
    double kappa1;
    double eta1;
    double kappa2;
    double eta2;
    double *dAlphaBetaX2 = NULL;
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
    double *p = NULL;
    int iYBufferSize = 0;
    double eps = 0.0;
    FILE *F = NULL;
    FILE *FE = NULL;

    int iMyRank = 0;
    int iPCount = 0;
    int iPrevRank = 0;
    int iNextRank = 0;
    char cStr[10];
    int *iPLineCountX1 = NULL; // Distribution of mesh line counts per processes to find solution
    // in direction x1. Mesh line fixed by value of coordinate i2.
    int *iPLineCountX2 = NULL; // Distribution of mesh line counts per processes to find solution in direction x2
    // Mesh line fixed by value of coordinate i1.
    int tag = 0;
    MPI_Status status;
    FILE *FL = NULL;

    clock_t cp1, cp2;
    double time1, time2;
    int iError = 0;


    // Initialization of MPI
    MPI_Init ( &argc, &argv );
    MPI_Comm_rank ( MPI_COMM_WORLD, &iMyRank );
    MPI_Comm_size ( MPI_COMM_WORLD, &iPCount );

    iPrevRank = (iMyRank - 1 + iPCount) % iPCount;
    iNextRank = (iMyRank + 1) % iPCount;

    if ( iMyRank == 0 )
    {
        printf ( "\nSolving the heat equation by variable directions method\n" );
        printf ( "Parallel mode, version 2.2b\n" );
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


    // Calculate the count of mesh lines for each process in both direction x1 and x2
    iPLineCountX1 = ( int* ) calloc ( iPCount, sizeof ( int ) );
    iPLineCountX2 = ( int* ) calloc ( iPCount, sizeof ( int ) );
    for ( i = 0; i < iPCount; i++ )
    {
        iPLineCountX1[i] = ( n2 - 1 ) / iPCount + ( i < ( n2 - 1 ) % iPCount );
        iPLineCountX2[i] = ( n1 - 1 ) / iPCount + ( i < ( n1 - 1 ) % iPCount );
    }

    // Calculate start values of variables x1 and x2 for each block in "sudoku" grid
    pX1start = ( double* ) calloc ( iPCount, sizeof ( int ) );
    pX2start = ( double* ) calloc ( iPCount, sizeof ( int ) );
    for ( i = 0, x1start = h1, x2start = h2; i < iPCount; i++ )
    {
        pX1start[i] = x1start;
        pX2start[i] = x2start;
        for ( j = 0; j < iPLineCountX1[i]; j++ )
        {
            x1start += h1;
            x2start += h2;
        }
    }


#ifdef WITH_PROCESS_LOG
    fprintf ( FL, "Count of processed mesh lines in direction x1: %d \n", iPLineCountX1[iMyRank] );
#endif

    int i1_gl, i2_gl;
    int iCurTile1_Size = 0;
    int iCurTile2_Size = 0;
    int iBufferABTileSize = 0;
    double *pYNextStartPos = NULL;
    double *pYLastLine = NULL;
    double *pYPrevCol = NULL;
    double *pYCurrentCol = NULL;
    double *pYLastCol = NULL;

    int iNextTile2_Size = 0;
    int i1_gl_next = 0;

    int iMaxTile1_Size = iPLineCountX2[0];
    int iMaxTile2_Size = iPLineCountX1[0];

    // For each tile we should have a column to receive solutions from previous process in reverse stage.
    // For comfort we have column on each boundary of tile,
    // but both of them will be used only for first tile in X1 direction.
    int iX1_YLineSize = n1 - 1 + 2 * iPCount;

    // Memory allocation to store mesh data and coefficients of linear equations
    iYBufferSize = ( iMaxTile2_Size + 2 ) * iX1_YLineSize;
    yCurrent = ( double* ) calloc ( iYBufferSize, sizeof ( double ) );
    yPrev = ( double* ) calloc ( iYBufferSize, sizeof ( double ) );
    dAlphaBetaX2 = ( double* ) calloc ( MAX( 2 * ( n1 - 1 ) * ( iMaxTile2_Size + 1 ),
                                             2 * ( n2 - 1 ) * ( iMaxTile1_Size + 1 ) ),
                                        sizeof ( double ) );
    iBufferABTileSize = 2 * (iMaxTile2_Size + 1) * (iMaxTile1_Size + 1);


    // Initialization of mesh function on 0-layer
    for (i1_gl = 0, pYCurrentStartPos = yCurrent;
         i1_gl < iPCount;
         i1_gl++, pYCurrentStartPos += iPLineCountX2[i1_gl] + 2)
    {
        int x2index = (iMyRank - i1_gl + iPCount) % iPCount;

        for (i2 = 0, x2 = pX2start[x2index] - h2, pYCurrentLine = pYCurrentStartPos;
             i2 <= iPLineCountX1[x2index] + 1;
             i2++, x2 += h2, pYCurrentLine += iX1_YLineSize)
        {
            for (i1 = 0, x1 = pX1start[i1_gl] - h1, p = pYCurrentLine;
                 i1 <= iPLineCountX2[i1_gl] + 1;
                 i1++, x1 += h1, p++)
            {
                *p = u(x1, x2, 0.0);
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

        swap_pointers(yCurrent, yPrev);

        for (i1_gl = 0, pABStartPos = dAlphaBetaX2, pYPrevStartPos = yCurrent + 1;
             i1_gl < iPCount;
             i1_gl++, pABStartPos += iBufferABTileSize, pYPrevStartPos += iCurTile1_Size + 2 )
        {
            iCurTile1_Size = iPLineCountX2[i1_gl];

            for (i2_gl = 0; i2_gl < iPCount; i2_gl++)
            {
                iCurTile2_Size = iPLineCountX1[i2_gl];

                if (iMyRank == (i1_gl + i2_gl) % iPCount)
                {
                    // If i1_gl is not 0, process receives bound values of coefficients
                    // alpha and beta from process iPrevRank for current tile.
                    // If i1_gl is 0, processes calculates all coefficients alpha(1,i2) and beta(1,i2)

                    if (i1_gl != 0)
                    {
                        MPI_Recv(dAlphaBetaX2, iMaxTile2_Size << 1, MPI_DOUBLE, iPrevRank, tag, MPI_COMM_WORLD,
                                 &status);
                    }
                    else
                        {
                        for (i2 = 1, x2 = pX2start[i2_gl], pABInLinePos = pABStartPos;
                             i2 <= iCurTile2_Size;
                             i2++, x2 += h2, pABInLinePos += 2) {

                            dCoefA1_L = 0.5 * ( k2 ( 0.0, x2, tPrev ) + k2 ( 0.0, x2 - h2, tPrev ) );
                            dCoefA1_R = 0.5 * ( k2 ( 0.0, x2 + h2, tPrev ) + k2 ( 0.0, x2, tPrev ) );
                            dCoefA2_L = 0.5 * ( k2 ( 0.0, x2, tPrev + tau ) + k2 ( 0.0, x2 - h2, tPrev + tau ) );
                            dCoefA2_R = 0.5 * ( k2 ( 0.0, x2 + h2, tPrev + tau ) + k2 ( 0.0, x2, tPrev + tau ) );
                            kappa1 = 0.0;
                            eta1 = 0.5 * ( m1 ( x2, tPrev ) + m1 ( x2, tPrev + tau ) ) -
                                   ( ( dCoefA2_L * m1 ( x2 - h2, tPrev + tau ) - ( dCoefA2_L + dCoefA2_R ) * m1 ( x2, tPrev + tau ) +
                                       dCoefA2_R * m1 ( x2 + h2, tPrev + tau ) ) -
                                     ( dCoefA1_L * m1 ( x2 - h2, tPrev ) - ( dCoefA1_L + dCoefA1_R ) * m1 ( x2, tPrev ) +
                                       dCoefA1_R * m1 ( x2 + h2, tPrev ) ) ) *
                                   0.5 / h2SquareTau;
                            pABInLinePos[i2] = kappa1;
                            pABInLinePos[i2 + 1] = eta1;
                        }
                    }

                    // Calculate coefficients alpha[] and beta[] in current tile
                    for (i1 = 1, x1 = pX1start[i1_gl], pABLine = pABStartPos, pYPrevCol = pYPrevStartPos;
                         i1 <= iCurTile1_Size;
                         i1++, x1 += h1, pABLine += iCurTile2_Size << 1, pYPrevCol += 1)
                    {
                        // Scan nodes in column of tile and consiquently calculate coefficients alpha and beta
                        // associated with nodes
                        x1start = pX1start[i1_gl];

                        dCoefA1_L = 0.5 * (k1(x1start, x2, t) + k1(x1start - h1, x2, t));
                        for (i2 = 1, x2 = pX2start[i2_gl], pABInLinePos = pABLine;
                             i2 <= iCurTile2_Size;
                             i2++, x2 += h2, pABInLinePos += 2)
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
                                      ( dCoefA2_L * pYPrevCol[i1 - iX1_YLineSize] -
                                       ( dCoefA2_L + dCoefA2_R ) * pYPrevCol[i1] +
                                       dCoefA2_R * pYPrevCol[i1 + iX1_YLineSize] ) +
                                      h1SquareTau * pYPrevCol[i1] + h1Square * f ( x1, x2, tPrev );
                            dCoefA1_L = dCoefA1_R;

                            // Calculate coefficients alpha(x1,x2+h2) and beta(x1,x2+h2)
                            dDiv = dCValue - pABInLinePos[0] * dAValue;
                            pABInLinePos[iCurTile2_Size << 1] = dBValue / dDiv;
                            pABInLinePos[(iCurTile2_Size << 1) + 1] = (dAValue * pABInLinePos[1] + dFValue) / dDiv;
                        }
                    }

                    // If i1_gl is not last tile, process sends bound values of coefficients alpha and beta
                    // to process iNextRank.
                    if (i1_gl != iPCount - 1) {
                        MPI_Send(pABLine, iCurTile1_Size << 1, MPI_DOUBLE, iNextRank, tag, MPI_COMM_WORLD);
                    }

                }
            }
        }


        ////////////////////////////////////////////////////////////////////////////////////
        // Second (reverse) stage of progonka - calculation of solution.

        for (i1_gl = iPCount - 1,
             pABStartPos = dAlphaBetaX2 + iBufferABTileSize * (iPCount - 1),
             pYCurrentStartPos = pYPrevStartPos - 1;
             i1_gl >= 0;
             i1_gl--, pABStartPos -= iBufferABTileSize, pYCurrentStartPos -= iCurTile1_Size + 2 )
        {
            iCurTile1_Size = iPLineCountX2[i1_gl];

            for (i2_gl = 0; i2_gl < iPCount; i2_gl++)
            {
                iCurTile2_Size = iPLineCountX1[i2_gl];

                if (iMyRank == (i1_gl + i2_gl) % iPCount)
                {
                    pABLastLine = pABStartPos + iCurTile1_Size * (iCurTile2_Size << 1);
                    pYLastCol = pYCurrentStartPos;

                    if (i1_gl != iPCount - 1)
                    {
                        MPI_Datatype coltype;
                        MPI_Type_vector( iCurTile2_Size, 1, iX1_YLineSize, MPI_DOUBLE, &coltype);
                        MPI_Type_commit(&coltype);
                        MPI_Recv( pYCurrentStartPos - iCurTile1_Size - 2, 1, coltype, iNextRank, tag, MPI_COMM_WORLD, &status );
                    }
                    else
                    {
                        for (i2 = 0, x2 = pX2start[i2_gl], pABInLinePos = pABLastLine;
                             i2 < iCurTile2_Size;
                             i2++, x2 += h2, pABInLinePos += 2)
                        {
                            kappa2 = 0.0;
                            eta2 = m4(x1, t);
                            pYLastCol[i2 * iX1_YLineSize] = (eta2 + kappa2 * pABInLinePos[1]) / (1 - kappa2 * pABInLinePos[0]);
                        }
                    }

                    // Calculate solution in current tile.
                    // Scan columns of tile (the column is a set of tile nodes with fixed i1)
                    for (i1 = iCurTile1_Size, pABLine = pABLastLine,
                         pYCurrentCol = pYLastCol;
                         i1 >= 0;
                         i1--, pABLine -= iCurTile2_Size << 1, pYCurrentCol -= 1)
                    {
                        // Scan nodes in column and consiquently calculate solution associated with nodes
                        for (i2 = 0, pABInLinePos = pABLine;
                             i2 < iCurTile2_Size;
                             i2++, pABInLinePos += 2)
                        {
                            pYCurrentCol[i2 * iX1_YLineSize - 1] = pABInLinePos[0] * pYCurrentCol[i2 * iX1_YLineSize] + pABInLinePos[1];
                        }
                    }

                    // If i1_gl is not zero, process sends bound values of solution to process iPrevRank.
                    if (i1_gl != 0)
                    {
                        MPI_Datatype coltype;
                        MPI_Type_vector( iCurTile2_Size, 1, iX1_YLineSize, MPI_DOUBLE, &coltype);
                        MPI_Type_commit(&coltype);
                        MPI_Send(pYCurrentCol + 1, 1, coltype, iPrevRank, tag, MPI_COMM_WORLD);
                    }

                    // Initialize boundary values at x2 = 0 and x2 = L2
                    if (i2_gl == 0 || i2_gl == iPCount - 1)
                    {
                        p = pYCurrentCol + 1;
                        for (i1 = 0, x1 = pX1start[i1_gl] - h1;
                             i1 <= iPLineCountX2[i1_gl] + 1;
                             i1++, x1 += h1, p += 1) {
                            if (i1_gl == 0) {
                                p[0] = m3(x2, t);
                            }
                            else {
                                p[(iCurTile2_Size + 1) * iX1_YLineSize] = m4(x2, t);
                            }
                        }
                    }
                }
            }
        }

        tag++;

        //--------------------------------------------------------//
        //   Calculate mesh function on j layer in direction x2   //
        //--------------------------------------------------------//

        tPrevSub = t;
        t = tPrev + tau;

        // Shift solution obtained on j-1 step
        swap_pointers(yCurrent, yPrev);

        ///////////////////////////////////////////////////////////////////////////////////////
        // First (forward) stage of progonka - calculation of coefficients alpha[] and beta[]

        for (i2_gl = 0; i2_gl < iPCount; i2_gl++)
        {
            iCurTile2_Size = iPLineCountX1[i2_gl];
            pYPrevStartPos = yPrev + iX1_YLineSize + 1;

            for (i1_gl = 0, pABStartPos = dAlphaBetaX2, pYPrevStartPos = yCurrent + iX1_YLineSize + 1;
                 i1_gl < iPCount;
                 i1_gl++, pABStartPos += iBufferABTileSize, pYPrevStartPos += iCurTile1_Size + 2)
            {
                iCurTile1_Size = iPLineCountX2[i1_gl];

                if (iMyRank == (i1_gl + i2_gl) % iPCount)
                {
                    // If i2_gl is not 0, process receives bound values of coefficients
                    // alpha and beta from process iPrevRank for current tile.
                    // If i2_gl is 0, processes calculates all coefficients alpha(i1,1) and beta(i1,1)

                    if (i2_gl != 0) {
                        MPI_Recv(pABStartPos, iMaxTile1_Size << 1, MPI_DOUBLE, iPrevRank, tag, MPI_COMM_WORLD, &status);
                    }
                    else {
                        for (i1 = 1, x1 = pX1start[i1_gl], pABInLinePos = pABStartPos;
                             i1 <= iCurTile1_Size;
                             i1++, x1 += h1, pABInLinePos += 2) {
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

                    // If i2_gl is not last tile, process sends bound values of coefficients alpha and beta
                    // to process iNextRank.
                    if (i2_gl != iPCount - 1) {
                        MPI_Send(pABLine, iCurTile1_Size << 1, MPI_DOUBLE, iNextRank, tag, MPI_COMM_WORLD);
                    }
                }
            }
        }


        ////////////////////////////////////////////////////////////////////////////////////
        // Second (reverse) stage of progonka - calculation of solution.

        for (i2_gl = iPCount - 1; i2_gl >= 0; i2_gl--)
        {
            iCurTile2_Size = iPLineCountX1[i2_gl];

            for (i1_gl = 0, pABStartPos = dAlphaBetaX2, pYCurrentStartPos = yCurrent + 1;
                 i1_gl < iPCount;
                 i1_gl++, pABStartPos += iBufferABTileSize, pYCurrentStartPos += iCurTile1_Size + 2)
            {
                iCurTile1_Size = iPLineCountX2[i1_gl];

                if (iMyRank == (i1_gl + i2_gl) % iPCount)
                {
                    pABLastLine = pABStartPos + iCurTile2_Size * (iCurTile1_Size << 1);
                    pYLastLine = pYCurrentStartPos + (iCurTile2_Size + 1) * iX1_YLineSize;

                    if (i2_gl != iPCount - 1) {
                        MPI_Recv(pYCurrentStartPos + (iCurTile2_Size + 1) * iX1_YLineSize + 1,
                                  iMaxTile1_Size, MPI_DOUBLE, iNextRank, tag, MPI_COMM_WORLD, &status);
                    }
                    else {
                        for (i1 = 0, x1 = pX1start[i1_gl], pABInLinePos = pABLastLine;
                             i1 < iCurTile1_Size;
                             i1++, x1 += h1, pABInLinePos += 2)
                        {
                            kappa2 = 0.0;
                            eta2 = m4(x1, t);
                            pYLastLine[i1] = (eta2 + kappa2 * pABInLinePos[1]) / (1 - kappa2 * pABInLinePos[0]);
                        }
                    }

                    // Calculate solution in current tile.
                    // Scan lines of tile (the line is a set of tile nodes with fixed i2)
                    for (i2 = iCurTile2_Size, pABLine = pABLastLine,
                         pYCurrentLine = pYLastLine;
                         i2 >= 0;
                         i2--, pABLine -= iCurTile1_Size << 1, pYCurrentLine -= iX1_YLineSize)
                    {
                        // Scan nodes in line and consiquently calculate solution associated with nodes
                        for (i1 = 0, pABInLinePos = pABLine;
                             i1 < iCurTile1_Size;
                             i1++, pABInLinePos += 2)
                        {
                            pYCurrentLine[i1 - iX1_YLineSize] = pABInLinePos[0] * pYCurrentLine[i1] + pABInLinePos[1];
                        }
                    }

                    // If i2_gl is not zero, process sends bound values of solution to process iPrevRank.
                    // All processes send their data in blocking mode
                    if (i2_gl != 0)
                    {
                        MPI_Send(pYCurrentLine + iX1_YLineSize, iCurTile1_Size, MPI_DOUBLE,
                                  iPrevRank, tag, MPI_COMM_WORLD);
                    }

                    // Initialize boundary values at x1 = 0 and x1 = L1
                    if (i1_gl == 0 || i1_gl == iPCount - 1)
                    {
                        p = yCurrent;
                        for (i2 = 0, x2 = pX2start[i2_gl] - h2;
                             i2 <= iPLineCountX1[i2_gl] + 1;
                             i2++, x2 += h2, p += iX1_YLineSize) {
                            if (i1_gl == 0) {
                                p[0] = m1(x2, t);
                            }
                            else {
                                p[iX1_YLineSize - 1] = m2(x2, t);
                            }
                        }
                    }
                }
            }
        }
        tag++;
    }

    MPI_Barrier ( MPI_COMM_WORLD );
    cp2 = clock ();
    time2 = MPI_Wtime ();


    // Gather mesh solution onto zero process

    if ( iMyRank == 0 ) {
        dY = (double *) calloc(iX1_YLineSize * (n2 + 1), sizeof(double));
    }
    p = dY;

    double *dYStartPos;
    int iLinesCount, iLineSize;
    for (i2_gl = 0, dYStartPos = dY;
         i2_gl < iPCount;
         i2_gl++, dYStartPos += iLinesCount * iX1_YLineSize)
    {
        iLinesCount = iPLineCountX1[i2_gl];
        if ( i2_gl == iPCount - 1 || i2_gl == 0 )
            iLinesCount++;

        for (i1_gl = 0, pYCurrentStartPos = yCurrent, p = dYStartPos;
             i1_gl < iPCount;
             i1_gl++, pYCurrentStartPos += iPLineCountX2[i1_gl] + 2, p += iLineSize)
        {
            iLineSize = iPLineCountX2[i1_gl];
            int iCurProc = (i1_gl + i2_gl) % iPCount;

            MPI_Datatype blocktype;
            MPI_Type_vector( iLinesCount, iLineSize, iX1_YLineSize, MPI_DOUBLE, &blocktype);
            MPI_Type_commit(&blocktype);

            //send-receive iLinesCount of lines-solutions with length = iLineLength, placed with offset = iX1_YLineSize
            if (iMyRank == iCurProc)
            {
                if (i2_gl != 0)
                    MPI_Send(pYCurrentStartPos + iX1_YLineSize, 1, blocktype, 0, MPI_ANY_TAG, MPI_COMM_WORLD);
                else
                    MPI_Send(pYCurrentStartPos, 1, blocktype, 0, MPI_ANY_TAG, MPI_COMM_WORLD);
            }

            if ( iMyRank == 0 )
            {
                MPI_Recv ( p, 1, blocktype, iCurProc, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
            }
        }
    }


// Save mesh solution into file
#ifndef WITHOUT_RES_SAVE
    if ( iMyRank == 0 )
    {
        F = fopen ( "res.parSudoku", "w" );
        p = dY;
        for ( i2 = 0, x2 = 0.0;  i2 <= n2;  i2++, x2 += h2 )
        {
            for ( i1 = 0, x1 = 0.0;  i1 <= n1;  i1++, x1 += h1, p++ )  fprintf ( F, "%e ", *p );
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
    if ( yCurrent != NULL ) free ( yCurrent );
    if ( yPrev != NULL ) free ( yPrev );
    if ( dAlphaBetaX2 != NULL ) free ( dAlphaBetaX2 );
    if ( iPLineCountX1 != NULL ) free ( iPLineCountX1 );
    if ( FL != NULL ) fclose ( FL );

    MPI_Finalize ();

    return 0;
}

void swap_pointers(double** yCurrent, double** yPrev) {
    double* p = *yPrev;
    *yPrev = *yCurrent;
    *yCurrent = p;
}
