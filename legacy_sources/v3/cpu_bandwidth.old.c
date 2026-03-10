

/*----------------------------------------------------------------------------*/
/*                                                                            */
/*                               LPlib V4.00                                  */
/*                                                                            */
/*----------------------------------------------------------------------------*/
/*                                                                            */
/*   Description:       benchmarking code to evaluate memory bandwidth        */
/*                      with various access patterns: direct, indirect,       */
/*                      ragged tables and vectorized ragged tables            */
/*                      direct and indirect memory reads                      */
/*   Author:            Loic MARECHAL                                         */
/*   Creation date:     mar 20 2025                                           */
/*   Last modification: mar 28 2025                                           */
/*                                                                            */
/*----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------*/
/* Includes                                                                   */
/*----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "lplib3.h"
#include "libmeshb7.h"


/*----------------------------------------------------------------------------*/
/* Defines                                                                    */
/*----------------------------------------------------------------------------*/

typedef float v4f   __attribute__ ((vector_size ( 16)));
typedef int   v4i   __attribute__ ((vector_size ( 16)));
typedef float v32f  __attribute__ ((vector_size (128)));
typedef int   v32i  __attribute__ ((vector_size (128)));
typedef int   v64i  __attribute__ ((vector_size (256)));


/*----------------------------------------------------------------------------*/
/* Structure's prototype                                                      */
/*----------------------------------------------------------------------------*/

typedef struct
{
   int      TetTyp, VerTyp, NmbVer, NmbTet, *VerBal, *VerDeg, *VerAdr;
   int      VerTypD32, VerTypD64, BalD32Idx, BalD64Idx, NmbVerD32, NmbVerD64;
   v4i      *TetVer;
   v32i     *BalD32;
   v64i     *BalD64;
   float    *VerScaSrc, *VerScaDst, *TetScaSrc, *TetScaDst;
   v4f      *VerCrd;
   v32f     *VerVecSrc, *VerVecDst, *TetVecSrc, *TetVecDst;
   int64_t  ParIdx;
}MshSct;


/*----------------------------------------------------------------------------*/
/* Macro instructions                                                         */
/*----------------------------------------------------------------------------*/

#define MIN(a,b)        ((a) < (b) ? (a) : (b))
#define MAX(a,b)        ((a) > (b) ? (a) : (b))


/*----------------------------------------------------------------------------*/
/* Loop over vertices, read and write data in a direct way                    */
/*----------------------------------------------------------------------------*/

void DirMemVerSca(int BegIdx, int EndIdx, int PthIdx, MshSct *msh)
{
   for(int i=BegIdx;i<=EndIdx;i++)
      msh->VerScaDst[i] = msh->VerScaSrc[i];
}

void DirMemVerVec(int BegIdx, int EndIdx, int PthIdx, MshSct *msh)
{
   for(int i=BegIdx;i<=EndIdx;i++)
      msh->VerVecDst[i] = msh->VerVecSrc[i];
}


/*----------------------------------------------------------------------------*/
/* Loop over tets, read and write data in a direct way                        */
/*----------------------------------------------------------------------------*/

void DirMemTetSca(int BegIdx, int EndIdx, int PthIdx, MshSct *msh)
{
   for(int i=BegIdx;i<=EndIdx;i++)
      msh->TetScaDst[i] = msh->TetScaSrc[i];
}

void DirMemTetVec(int BegIdx, int EndIdx, int PthIdx, MshSct *msh)
{
   for(int i=BegIdx;i<=EndIdx;i++)
      msh->TetVecDst[i] = msh->TetVecSrc[i];
}


/*----------------------------------------------------------------------------*/
/* Indirect memory access from tets to vertices data                          */
/*----------------------------------------------------------------------------*/

void IndMemSca(int BegIdx, int EndIdx, int PthIdx, MshSct *msh)
{
   for(int i=BegIdx;i<=EndIdx;i++)
   {
      msh->TetScaDst[i] = msh->VerScaSrc[ msh->TetVer[i][0] ]
                        + msh->VerScaSrc[ msh->TetVer[i][1] ]
                        + msh->VerScaSrc[ msh->TetVer[i][2] ]
                        + msh->VerScaSrc[ msh->TetVer[i][3] ];
   }
}

void IndMemVec(int BegIdx, int EndIdx, int PthIdx, MshSct *msh)
{
   for(int i=BegIdx;i<=EndIdx;i++)
   {
      msh->TetVecDst[i] = msh->VerVecSrc[ msh->TetVer[i][0] ]
                        + msh->VerVecSrc[ msh->TetVer[i][1] ]
                        + msh->VerVecSrc[ msh->TetVer[i][2] ]
                        + msh->VerVecSrc[ msh->TetVer[i][3] ];
   }
}


/*----------------------------------------------------------------------------*/
/* Indirect ragged memory access from vertices to tets data (ball of point)   */
/*----------------------------------------------------------------------------*/

void RagMemSca(int BegIdx, int EndIdx, int PthIdx, MshSct *msh)
{
   int i, j, adr, deg;
   float dat;

   for(i=BegIdx;i<=EndIdx;i++)
   {
      adr = msh->VerAdr[i];
      deg = msh->VerDeg[i];
      dat = 0.;

      for(j=0;j<deg;j++)
         dat += msh->TetScaSrc[ msh->VerBal[ adr + j ] ];

      msh->VerScaDst[i] = dat;
   }
}

void RagMemVec(int BegIdx, int EndIdx, int PthIdx, MshSct *msh)
{
   int i, j, adr, deg;
   v32f dat;

   for(i=BegIdx;i<=EndIdx;i++)
   {
      adr = msh->VerAdr[i];
      deg = msh->VerDeg[i];
      dat = (v32f){0.};

      for(j=0;j<deg;j++)
         dat += msh->TetVecSrc[ msh->VerBal[ adr + j ] ];

      msh->VerVecDst[i] = dat;
   }
}


/*----------------------------------------------------------------------------*/
/* Vectorized Indirect ragged memory access from vertices to tets data        */
/*----------------------------------------------------------------------------*/

void RagMemD32Sca(int BegIdx, int EndIdx, int PthIdx, MshSct *msh)
{
   int i, j;
   float dat;
   v32i bal;

   // Handle only vertices whose degree is <= 32
   for(i=BegIdx;i<=EndIdx;i++)
   {
      dat = 0.;
      bal = msh->BalD32[i];

      for(j=0;j<32;j++)
         if(bal[j])
            dat += msh->TetScaSrc[ bal[j] ];

      msh->VerScaDst[i] = dat;
   }
}

void RagMemD32Vec(int BegIdx, int EndIdx, int PthIdx, MshSct *msh)
{
   int i, j;
   v32f dat;
   v32i bal;

   // Handle only vertices whose degree is <= 32
   for(i=BegIdx;i<=EndIdx;i++)
   {
      dat = (v32f){0.};
      bal = msh->BalD32[i];

      for(j=0;j<32;j++)
         if(bal[j])
            dat += msh->TetVecSrc[ bal[j] ];

      msh->VerVecDst[i] = dat;
   }
}


/*----------------------------------------------------------------------------*/
/* Vectorized Indirect ragged memory access from vertices to tets data        */
/*----------------------------------------------------------------------------*/

void RagMemD64Sca(int BegIdx, int EndIdx, int PthIdx, MshSct *msh)
{
   int i, j;
   float dat;
   v64i bal;

   // Handle only vertices whose degree is > 32
   for(i=BegIdx;i<=EndIdx;i++)
   {
      dat = 0.;
      bal = msh->BalD64[i];

      for(j=0;j<32;j++)
         dat += msh->TetScaSrc[ bal[j] ];

      for(j=32;j<64;j++)
         if(bal[j])
            dat += msh->TetScaSrc[ bal[j] ];

      msh->VerScaDst[ msh->BalD32Idx + i ] = dat;
   }
}

void RagMemD64Vec(int BegIdx, int EndIdx, int PthIdx, MshSct *msh)
{
   int i, j;
   v32f dat;
   v64i bal;

   // Handle only vertices whose degree is > 32
   for(i=BegIdx;i<=EndIdx;i++)
   {
      dat = (v32f){0.};
      bal = msh->BalD64[i];

      for(j=0;j<32;j++)
         dat += msh->TetVecSrc[ bal[j] ];

      for(j=32;j<64;j++)
         if(bal[j])
            dat += msh->TetVecSrc[ bal[j] ];

      msh->VerVecDst[ msh->BalD32Idx + i ] = dat;
   }
}


/*----------------------------------------------------------------------------*/
/* Setup the LPlib and launch direct and indirect memory access loops         */
/*----------------------------------------------------------------------------*/

int main(int ArgCnt, char **ArgVec)
{
   int      i, j, k, ver, dim, ref, NmbCpu, NmbItr;
   int64_t  InpMsh, NmbBal = 0;
   double   tim;
   char     *MshNam;
   MshSct   msh = {0};


   // -----
   // SETUP
   // -----

   // Read the number of threads to launch from the command line argument
   if(ArgCnt != 4)
   {
      puts("\ncpu_bandwidth   NmbThreads   NmbIterations   MeshFile\n");
      puts("In order to fully evaluate all kinds of memory access performances you should generate");
      puts("three different numberings from the same test mesh and feed them to cpu_bandwith:\n");
      puts(" hilbert -in MyTestMesh -out RandomMesh -scheme 2");
      puts(" hilbert -in MyTestMesh -out HilbertMesh");
      puts(" hilbert -in MyTestMesh -out SlicedMesh -gmlib generic\n");
      exit(0);
   }
   else
   {
      NmbCpu = atoi(*++ArgVec);
      NmbItr = atoi(*++ArgVec);
      MshNam = *++ArgVec;

      if(!NmbCpu)
         NmbCpu = GetNumberOfCores();
   }


   // --------------------------
   // I/O AND MEMORY ALLOCATIONS
   // --------------------------

   // Open the input mesh
   if(!(InpMsh = GmfOpenMesh(MshNam, GmfRead, &ver, &dim)))
      return(1);

   msh.NmbVer = GmfStatKwd(InpMsh, GmfVertices);
   msh.NmbTet = GmfStatKwd(InpMsh, GmfTetrahedra);

   // Allocate the memory
   msh.TetVer = malloc( (int64_t)(msh.NmbTet+1) * sizeof(v4i) );
   msh.TetScaSrc = malloc( (int64_t)(msh.NmbTet+1) * sizeof(float) );
   msh.TetScaDst = malloc( (int64_t)(msh.NmbTet+1) * sizeof(float) );
   msh.TetVecSrc = malloc( (int64_t)(msh.NmbTet+1) * sizeof(v32f) );
   msh.TetVecDst = malloc( (int64_t)(msh.NmbTet+1) * sizeof(v32f) );
   msh.VerCrd = malloc( (int64_t)(msh.NmbVer+1) * sizeof(v4f) );
   msh.VerScaSrc = malloc( (int64_t)(msh.NmbVer+1) * sizeof(float) );
   msh.VerScaDst = malloc( (int64_t)(msh.NmbVer+1) * sizeof(float) );
   msh.VerVecSrc = malloc( (int64_t)(msh.NmbVer+1) * sizeof(v32f) );
   msh.VerVecDst = malloc( (int64_t)(msh.NmbVer+1) * sizeof(v32f) );

   if(!msh.TetVer || !msh.VerCrd)
   {
      puts("Failed to allocate memory");
      exit(1);
   }

   printf("\nAllocated %.2f GB\n", (float)(msh.NmbTet + msh.NmbVer) * 70. * sizeof(float) / 1E9);

   // Read the vertices
   GmfGetBlock(InpMsh, GmfVertices, 1, msh.NmbVer, 0, NULL, NULL,
               GmfFloatVec, 4, &msh.VerCrd[1], &msh.VerCrd[ msh.NmbVer ],
               GmfInt, &ref, &ref);

   // Read the tets
   GmfGetBlock(InpMsh, GmfTetrahedra, 1, msh.NmbTet, 0, NULL, NULL,
               GmfIntVec, 4, &msh.TetVer[1], &msh.TetVer[ msh.NmbTet ],
               GmfInt, &ref, &ref);

   GmfCloseMesh(InpMsh);

   printf(  "Imported %d vertices and %d tets from the mesh file\n\n",
            msh.NmbVer, msh.NmbTet );


   // ----------
   // INIT LPLIB
   // ----------

   // Initialize the LPlib and setup the data types
   msh.ParIdx = InitParallel(NmbCpu);
   msh.TetTyp = NewType(msh.ParIdx, msh.NmbTet);
   msh.VerTyp = NewType(msh.ParIdx, msh.NmbVer);

   if(!msh.ParIdx || !msh.TetTyp || !msh.VerTyp)
   {
      puts("Error while initializing the LPlib or datatypes.");
      exit(1);
   }


   // ------------------------
   // BUILD RAGGED BALLS TABLE
   // ------------------------

   // Build the balls table
   msh.VerDeg = calloc(msh.NmbVer + 1, sizeof(int));
   msh.VerAdr = malloc((int64_t)(msh.NmbVer + 1) * sizeof(int));

   if(!msh.VerDeg || !msh.VerAdr )
   {
      puts("Failed to allocate memory");
      exit(1);
   }

   for(i=1;i<=msh.NmbTet;i++)
      for(j=0;j<4;j++)
         msh.VerDeg[ msh.TetVer[i][j] ]++;

   for(i=1;i<=msh.NmbVer;i++)
   {
      msh.VerAdr[i] = NmbBal;
      NmbBal += msh.VerDeg[i];
      msh.VerDeg[i] = 0;
   }

   msh.VerBal = malloc(NmbBal * sizeof(int));

   if(!msh.VerBal)
   {
      puts("Failed to allocate memory");
      exit(1);
   }

   for(i=1;i<=msh.NmbTet;i++)
      for(j=0;j<4;j++)
      {
         k = msh.TetVer[i][j];
         msh.VerBal[ msh.VerAdr[k] + msh.VerDeg[k] ] = i;
         msh.VerDeg[k]++;
      }


   // --------------------------------
   // BUILD THE VECTORIZED BALLS TABLE
   // --------------------------------

   // Build the vectorized balls table
   for(i=1;i<=msh.NmbVer;i++)
      if(msh.VerDeg[i] > 32)
      {
         msh.BalD32Idx = i;
         msh.BalD64Idx = msh.NmbVer;
         break;
      }

   msh.NmbVerD32 = msh.BalD32Idx;
   msh.NmbVerD64 = msh.BalD64Idx - msh.BalD32Idx;

   msh.BalD32 = calloc(msh.NmbVerD32 + 1, sizeof(v32i));
   msh.BalD64 = calloc(msh.NmbVerD64 + 1, sizeof(v64i));

   if(!msh.BalD32 || !msh.BalD64 )
   {
      puts("Failed to allocate memory");
      exit(1);
   }

   msh.VerTypD32 = NewType(msh.ParIdx, msh.NmbVerD32);
   msh.VerTypD64 = NewType(msh.ParIdx, msh.NmbVerD64);

   // Copy the balls with a degree <= 32
   for(i=1;i<=msh.BalD32Idx;i++)
      for(j=0;j<msh.VerDeg[i];j++)
         msh.BalD32[i][j] = msh.VerBal[ msh.VerAdr[i] + j ];

   // Copy the balls with a degree > 32 and cap them to 64
   for(i=msh.BalD32Idx + 1; i<=msh.BalD64Idx; i++)
      for(j=0;j<MIN(64,msh.VerDeg[i]);j++)
         msh.BalD64[ i - msh.BalD32Idx ][j] = msh.VerBal[ msh.VerAdr[i] + j ];


   // -----------------------------
   // RUN DIRECT ACCESS MEMORY TEST
   // -----------------------------

   // Perform parallel direct memory access loops on vertices and tets
   tim = GetWallClock();

   for(i=1;i<=NmbItr;i++)
   {
      LaunchParallel(msh.ParIdx, msh.VerTyp, 0, (void *)DirMemVerSca, (void *)&msh);
      LaunchParallel(msh.ParIdx, msh.TetTyp, 0, (void *)DirMemTetSca, (void *)&msh);
   }

   tim = GetWallClock() - tim;

   printf("Direct reads            : %d steps, run time = %7.3fs, bandwidth = %6.1f GB/s\n",
            NmbItr, tim, (NmbItr * (8LL * msh.NmbVer + 8LL * msh.NmbTet)) / (tim * 1E9) );


   // -------------------------------
   // RUN INDIRECT ACCESS MEMORY TEST
   // -------------------------------

   // Perform parallel indirect memory access loops on tets
   tim = GetWallClock();

   for(i=1;i<=NmbItr;i++)
      LaunchParallel(msh.ParIdx, msh.TetTyp, 0, (void *)IndMemSca, (void *)&msh);

   tim = GetWallClock() - tim;

   printf("Indirect reads          : %d steps, run time = %7.3fs, bandwidth = %6.1f GB/s (unique reads = %6.1f GB/s)\n",
            NmbItr, tim, (NmbItr * 32LL * msh.NmbTet) / (tim * 1E9),
            (NmbItr * (4LL * msh.NmbVer + 20LL * msh.NmbTet)) / (tim * 1E9) );


   // --------------------------------------
   // RUN INDIRECT RAGGED ACCESS MEMORY TEST
   // --------------------------------------

   // Perform parallel ragged indirect memory access loops on vertices
   tim = GetWallClock();

   for(i=1;i<=NmbItr;i++)
      LaunchParallel(msh.ParIdx, msh.VerTyp, 0, (void *)RagMemSca, (void *)&msh);

   tim = GetWallClock() - tim;

   printf("Ragged reads            : %d steps, run time = %7.3fs, bandwidth = %6.1f GB/s (unique reads = %6.1f GB/s)\n",
            NmbItr, tim, (NmbItr * (8LL * NmbBal + 12LL * msh.NmbVer)) / (tim * 1E9),
            (NmbItr * (12LL * msh.NmbVer + 4LL * msh.NmbTet + 4LL * NmbBal)) / (tim * 1E9) );


   // -------------------------------------------------
   // RUN INDIRECT VECTORIZED RAGGED ACCESS MEMORY TEST
   // -------------------------------------------------

   // Perform parallel vectorized ragged indirect memory access loops on vertices
   tim = GetWallClock();

   for(i=1;i<=NmbItr;i++)
   {
      LaunchParallel(msh.ParIdx, msh.VerTypD32, 0, (void *)RagMemD32Sca, (void *)&msh);
      LaunchParallel(msh.ParIdx, msh.VerTypD64, 0, (void *)RagMemD64Sca, (void *)&msh);
   }

   tim = GetWallClock() - tim;

   printf("Vectorized ragged reads : %d steps, run time = %7.3fs, bandwidth = %6.1f GB/s (unique reads = %6.1f GB/s)\n\n",
            NmbItr, tim, (NmbItr * (4LL * NmbBal + 136LL * msh.NmbVerD32 + 264LL * msh.NmbVerD64)) / (tim * 1E9),
            (NmbItr * (136LL * msh.NmbVerD32 + 264LL * msh.NmbVerD64 + 4LL * msh.NmbTet)) / (tim * 1E9) );


   // -----------------------------
   // RUN DIRECT ACCESS MEMORY TEST
   // -----------------------------

   // Perform parallel direct memory access loops on vertices and tets
   tim = GetWallClock();

   for(i=1;i<=NmbItr;i++)
   {
      LaunchParallel(msh.ParIdx, msh.VerTyp, 0, (void *)DirMemVerVec, (void *)&msh);
      LaunchParallel(msh.ParIdx, msh.TetTyp, 0, (void *)DirMemTetVec, (void *)&msh);
   }

   tim = GetWallClock() - tim;

   printf("Direct reads            : %d steps, run time = %7.3fs, bandwidth = %6.1f GB/s\n",
            NmbItr, tim, (NmbItr * (256LL * msh.NmbVer + 256LL * msh.NmbTet)) / (tim * 1E9) );


   // -------------------------------
   // RUN INDIRECT ACCESS MEMORY TEST
   // -------------------------------

   // Perform parallel indirect memory access loops on tets
   tim = GetWallClock();

   for(i=1;i<=NmbItr;i++)
      LaunchParallel(msh.ParIdx, msh.TetTyp, 0, (void *)IndMemVec, (void *)&msh);

   tim = GetWallClock() - tim;

   printf("Indirect reads          : %d steps, run time = %7.3fs, bandwidth = %6.1f GB/s (unique reads = %6.1f GB/s)\n",
            NmbItr, tim, (NmbItr * 656LL * msh.NmbTet) / (tim * 1E9),
            (NmbItr * (128LL * msh.NmbVer + 144LL * msh.NmbTet)) / (tim * 1E9) );


   // --------------------------------------
   // RUN INDIRECT RAGGED ACCESS MEMORY TEST
   // --------------------------------------

   // Perform parallel ragged indirect memory access loops on vertices
   tim = GetWallClock();

   for(i=1;i<=NmbItr;i++)
      LaunchParallel(msh.ParIdx, msh.VerTyp, 0, (void *)RagMemVec, (void *)&msh);

   tim = GetWallClock() - tim;

   printf("Ragged reads            : %d steps, run time = %7.3fs, bandwidth = %6.1f GB/s (unique reads = %6.1f GB/s)\n",
            NmbItr, tim, (NmbItr * (132LL * NmbBal + 140LL * msh.NmbVer)) / (tim * 1E9),
            (NmbItr * (128LL * msh.NmbVer + 16LL * msh.NmbTet + 4LL * NmbBal)) / (tim * 1E9) );


   // -------------------------------------------------
   // RUN INDIRECT VECTORIZED RAGGED ACCESS MEMORY TEST
   // -------------------------------------------------

   // Perform parallel vectorized ragged indirect memory access loops on vertices
   tim = GetWallClock();

   for(i=1;i<=NmbItr;i++)
   {
      LaunchParallel(msh.ParIdx, msh.VerTypD32, 0, (void *)RagMemD32Vec, (void *)&msh);
      LaunchParallel(msh.ParIdx, msh.VerTypD64, 0, (void *)RagMemD64Vec, (void *)&msh);
   }

   tim = GetWallClock() - tim;

   printf("Vectorized ragged reads : %d steps, run time = %7.3fs, bandwidth = %6.1f GB/s (unique reads = %6.1f GB/s)\n\n",
            NmbItr, tim, (NmbItr * (128LL * NmbBal + 260LL * msh.NmbVerD32 + 388LL * msh.NmbVerD64)) / (tim * 1E9),
            (NmbItr * (160LL * msh.NmbVerD32 + 284LL * msh.NmbVerD64 + 16LL * msh.NmbTet)) / (tim * 1E9) );


   // -------
   // CLEANUP
   // -------

   StopParallel(msh.ParIdx);

   free(msh.TetVer);
   free(msh.VerCrd);
   free(msh.VerDeg);
   free(msh.VerAdr);
   free(msh.VerBal);
   free(msh.BalD32);
   free(msh.BalD64);

   return(0);
}
