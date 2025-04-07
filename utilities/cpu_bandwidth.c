

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
/*   Last modification: apr 07 2025                                           */
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

typedef float v4f   __attribute__ ((vector_size (16)));
typedef int   v4i   __attribute__ ((vector_size (16)));


/*----------------------------------------------------------------------------*/
/* Structure's prototype                                                      */
/*----------------------------------------------------------------------------*/

typedef struct
{
   int      TetTyp, VerTyp, NmbVer, NmbTet, *VerBal, *VerDeg, *VerAdr;
   v4i      *TetVer, *TetInt;
   v4f      *VerCrd, *VerDat, *TetDat;
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

void DirMemVer(int BegIdx, int EndIdx, int PthIdx, MshSct *msh)
{
   for(int i=BegIdx;i<=EndIdx;i++)
      msh->VerDat[i] = msh->VerCrd[i];
}


/*----------------------------------------------------------------------------*/
/* Loop over tets, read and write data in a direct way                        */
/*----------------------------------------------------------------------------*/

void DirMemTet(int BegIdx, int EndIdx, int PthIdx, MshSct *msh)
{
   for(int i=BegIdx;i<=EndIdx;i++)
      msh->TetInt[i] = msh->TetVer[i];
}


/*----------------------------------------------------------------------------*/
/* Indirect memory access from tets to vertices data                          */
/*----------------------------------------------------------------------------*/

void IndMem(int BegIdx, int EndIdx, int PthIdx, MshSct *msh)
{
   for(int i=BegIdx;i<=EndIdx;i++)
   {
      msh->TetDat[i] = msh->VerDat[ msh->TetVer[i][0] ]
                     + msh->VerDat[ msh->TetVer[i][1] ]
                     + msh->VerDat[ msh->TetVer[i][2] ]
                     + msh->VerDat[ msh->TetVer[i][3] ];
   }
}


/*----------------------------------------------------------------------------*/
/* Indirect ragged memory access from vertices to tets data (ball of point)   */
/*----------------------------------------------------------------------------*/

void RagMem(int BegIdx, int EndIdx, int PthIdx, MshSct *msh)
{
   int i, j, adr, deg;
   v4f dat;

   for(i=BegIdx;i<=EndIdx;i++)
   {
      adr = msh->VerAdr[i];
      deg = msh->VerDeg[i];
      dat = (v4f){0.,0.,0.,0.};

      for(j=0;j<deg;j++)
         dat += msh->TetDat[ msh->VerBal[ adr + j ] ];

      msh->VerDat[i] = dat;
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
      puts("two different numberings from the same test mesh and feed them to cpu_bandwith:\n");
      puts(" hilbert -in MyTestMesh -out RandomMesh -scheme 2");
      puts(" hilbert -in MyTestMesh -out HilbertMesh\n");
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
   msh.TetDat = malloc( (int64_t)(msh.NmbTet+1) * sizeof(v4f) );
   msh.VerCrd = malloc( (int64_t)(msh.NmbVer+1) * sizeof(v4f) );
   msh.VerDat = malloc( (int64_t)(msh.NmbVer+1) * sizeof(v4f) );
   msh.TetInt = (v4i *)msh.TetDat;

   if(!msh.TetVer || !msh.TetDat || !msh.VerCrd || !msh.VerDat)
   {
      puts("Failed to allocate memory");
      exit(1);
   }

   printf("\nAllocated %.2f GB\n", (float)(msh.NmbTet + msh.NmbVer) * 2. * sizeof(v4f) / 1E9);

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


   // -----------------------------
   // RUN DIRECT ACCESS MEMORY TEST
   // -----------------------------

   // Perform parallel direct memory access loops on vertices and tets
   tim = GetWallClock();

   for(i=1;i<=NmbItr;i++)
   {
      LaunchParallel(msh.ParIdx, msh.VerTyp, 0, (void *)DirMemVer, (void *)&msh);
      LaunchParallel(msh.ParIdx, msh.TetTyp, 0, (void *)DirMemTet, (void *)&msh);
   }

   tim = GetWallClock() - tim;

   printf("Direct reads   : %d steps, run time = %7.3fs, bandwidth = %6.1f GB/s\n",
            NmbItr, tim, (NmbItr * (32LL * msh.NmbVer + 32LL * msh.NmbTet)) / (tim * 1E9) );


   // -------------------------------
   // RUN INDIRECT ACCESS MEMORY TEST
   // -------------------------------

   // Perform parallel indirect memory access loops on tets
   tim = GetWallClock();

   for(i=1;i<=NmbItr;i++)
      LaunchParallel(msh.ParIdx, msh.TetTyp, 0, (void *)IndMem, (void *)&msh);

   tim = GetWallClock() - tim;

   printf("Indirect reads : %d steps, run time = %7.3fs, bandwidth = %6.1f GB/s (unique reads = %6.1f GB/s)\n",
            NmbItr, tim, (NmbItr * 96LL * msh.NmbTet) / (tim * 1E9),
            (NmbItr * (16LL * msh.NmbVer + 32LL * msh.NmbTet)) / (tim * 1E9) );


   // --------------------------------------
   // RUN INDIRECT RAGGED ACCESS MEMORY TEST
   // --------------------------------------

   // Perform parallel ragged indirect memory access loops on vertices
   tim = GetWallClock();

   for(i=1;i<=NmbItr;i++)
      LaunchParallel(msh.ParIdx, msh.VerTyp, 0, (void *)RagMem, (void *)&msh);

   tim = GetWallClock() - tim;

   printf("Ragged reads   : %d steps, run time = %7.3fs, bandwidth = %6.1f GB/s (unique reads = %6.1f GB/s)\n",
            NmbItr, tim, (NmbItr * (20LL * NmbBal + 24LL * msh.NmbVer)) / (tim * 1E9),
            (NmbItr * (24LL * msh.NmbVer + 16LL * msh.NmbTet + 4LL * NmbBal)) / (tim * 1E9) );


   // -------
   // CLEANUP
   // -------

   StopParallel(msh.ParIdx);

   free(msh.TetVer);
   free(msh.TetDat);
   free(msh.VerCrd);
   free(msh.VerDat);
   free(msh.VerDeg);
   free(msh.VerAdr);
   free(msh.VerBal);

   return(0);
}
