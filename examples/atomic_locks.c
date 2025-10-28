

/*----------------------------------------------------------------------------*/
/*                                                                            */
/*                      DEPENDNCY-LOOPS VERSUS ATMIC LOCKS                    */
/*                                                                            */
/*----------------------------------------------------------------------------*/
/*                                                                            */
/*   Description:       perform indirect memory writes and compare the runtime*/
/*                      between sequential, dependency-loop and atomic locks  */
/*   Author:            Loic MARECHAL                                         */
/*   Creation date:     oct 27 2025                                           */
/*   Last modification: oct 27 2025                                           */
/*                                                                            */
/*----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------*/
/* Includes                                                                   */
/*----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <stdint.h>
#include "lplib4.h"
#include "libmeshb8.h"


/*----------------------------------------------------------------------------*/
/* Structure's prototype                                                      */
/*----------------------------------------------------------------------------*/

typedef struct
{
   int      TetTyp, VerTyp, NmbVer, NmbTet, (*TetVer)[4], *VerDeg;
   int64_t  ParIdx;
   double   *VerTem, *TetTem;
}MshSct;


/*----------------------------------------------------------------------------*/
/* Procedure called in parallel adding the vertices'value                     */
/*  to the tetrahedron's value (INDIRECT READ)                                */
/* Arguments are imposed by the library:                                      */
/* BegIdx: begining loop index                                                */
/* EndIdx: ending loop index                                                  */
/* PthIdx: index of the thread executing this instance                        */
/* arg   : pointer to a user's structure containing every                     */
/*        arguments needed by the loop                                        */
/*----------------------------------------------------------------------------*/

void TetTem(int BegIdx, int EndIdx, int PthIdx, MshSct *msh)
{
   int i, j;
   double t;

   for(i=BegIdx;i<=EndIdx;i++)
   {
      t = 0;

      for(j=0;j<4;j++)
         t += msh->VerTem[ msh->TetVer[i][j] ];

      msh->TetTem[i] = t / 4.;
   }
}


/*----------------------------------------------------------------------------*/
/* Procedure called in parallel adding a tetrahedron's                        */
/*  temperature to the vertices' values (INDIRECT WRITE)                      */
/* Arguments are imposed by the library:                                      */
/* BegIdx: begining loop index                                                */
/* EndIdx: ending loop index                                                  */
/* PthIdx: index of the thread executing this instance                        */
/* arg   : pointer to a user's structure containing every                     */
/*        arguments needed by the loop                                        */
/*----------------------------------------------------------------------------*/

void VerTem(int BegIdx, int EndIdx, int PthIdx, MshSct *msh)
{
   int i, j;

   for(i=BegIdx;i<=EndIdx;i++)
      for(j=0;j<4;j++)
         msh->VerTem[ msh->TetVer[i][j] ] +=
            msh->TetTem[i] / (4. * msh->VerDeg[ msh->TetVer[i][j] ]);
}

void LokVer(int BegIdx, int EndIdx, int PthIdx, MshSct *msh)
{
   int i, j, VerIdx;

   for(i=BegIdx;i<=EndIdx;i++)
      for(j=0;j<4;j++)
      {
         VerIdx = msh->TetVer[i][j];
         AtomicLock(msh->ParIdx, msh->VerTyp, VerIdx);
         msh->VerTem[ VerIdx ] += msh->TetTem[i] / (4. * msh->VerDeg[ VerIdx ]);
         AtomicUnlock(msh->ParIdx, msh->VerTyp, VerIdx);
      }
}


/*----------------------------------------------------------------------------*/
/* Setup LPlib and launch 100 temperature soomthing steps                     */
/*----------------------------------------------------------------------------*/

int main(int ArgCnt, char **ArgVec)
{
   int      i, j, NmbCpu = GetNumberOfCores(), ParMod, ver, dim, ref, NmbItr = 100;
   int64_t  InpMsh;
   float    sta[2], acc=0;
   double   tim=0, res;
   MshSct   msh;

   // Read the number of threads to launch from the command line argument
   if(ArgCnt != 3)
   {
      puts("atomic_locks nbThreads mode");
      puts("choose mode from   1:serial,  2:dependency-loop,  3:atomic-locks");
      exit(0);
   }

   NmbCpu = atoi(*++ArgVec);
   ParMod = atoi(*++ArgVec);

   if(ParMod < 1 || ParMod > 3)
      exit(1);

   // Open the input mesh
   if(!(InpMsh = GmfOpenMesh("../sample_meshes/tet.meshb", GmfRead, &ver, &dim)))
      return(1);

   puts("");
   printf("Input mesh: idx = %lld, version = %d, dimension = %d\n",
            InpMsh, ver, dim);

   msh.NmbVer = (int)GmfStatKwd(InpMsh, GmfVertices);
   printf("Input mesh: nmb vertices = %d\n", msh.NmbVer);

   msh.NmbTet = (int)GmfStatKwd(InpMsh, GmfTetrahedra);
   printf("Input mesh: nmb tets = %d\n", msh.NmbTet);

   // Allocate the memory
   msh.TetVer = malloc( (int64_t)(msh.NmbTet+1) * 4 * sizeof(int) );
   msh.TetTem = malloc( (int64_t)(msh.NmbTet+1) * sizeof(double) );
   msh.VerTem = malloc( (int64_t)(msh.NmbVer+1) * sizeof(double) );
   msh.VerDeg = calloc( (int64_t)(msh.NmbVer+1), sizeof(int) );
   puts("alloc done");

   // Read the tets
   GmfGetBlock(InpMsh, GmfTetrahedra, 1, msh.NmbTet, 0, NULL, NULL,
               GmfIntVec, 4, msh.TetVer[1], msh.TetVer[ msh.NmbTet ],
               GmfInt, &ref, &ref);

   GmfCloseMesh(InpMsh);
   puts("read done");



   // Initialize the vertices' temperature with some crap values
   for(i=1;i<=msh.NmbVer;i++)
      msh.VerTem[i] = 1.;

   for(i=1;i<=msh.NmbTet;i++)
      for(j=0;j<4;j++)
         msh.VerDeg[ msh.TetVer[i][j] ]++;

   // Initialize the LPlib and setup the data types
   if(!(msh.ParIdx = InitParallel(NmbCpu)))
   {
      puts("Error initializing the LPlib.");
      exit(1);
   }

   if(!(msh.TetTyp = NewType(msh.ParIdx, msh.NmbTet)))
   {
      puts("Error while creating tetrahedra data type.");
      exit(1);
   }

   if(!(msh.VerTyp = NewType(msh.ParIdx, msh.NmbVer)))
   {
      puts("Error while creating vertices data type.");
      exit(1);
   }

   puts("");
   printf("TetTyp = %d, VerTyp = %d, NmbCpu = %d\n",
         msh.TetTyp, msh.VerTyp, NmbCpu);

   if(ParMod == 1)
   {
      tim = GetWallClock();

      for(i=1;i<=NmbItr;i++)
      {
         TetTem(1, msh.NmbTet, 0, &msh);
         VerTem(1, msh.NmbTet, 0, &msh);
      }

      tim = GetWallClock() - tim;

      printf("%d steps, sequential running time = %gs\n", NmbItr, tim);
   }
   else if(ParMod == 2)
   {
      // Setup dependencies between tets and vertices 
      // and compute vertices' degree on the flight
      BeginDependency(msh.ParIdx, msh.TetTyp, msh.VerTyp);

      for(i=1;i<=msh.NmbTet;i++)
         for(j=0;j<4;j++)
            AddDependency(msh.ParIdx, i, msh.TetVer[i][j]);

      EndDependency(msh.ParIdx, sta);

      printf("dependencies stats: average = %g %%, maximum = %g %%\n",
               sta[0], sta[1]);

      // Perform temperature smoothing steps
      tim = GetWallClock();

      puts("");
      for(i=1;i<=NmbItr;i++)
      {
         acc += LaunchParallel(msh.ParIdx, msh.TetTyp,          0, TetTem, &msh);
         acc += LaunchParallel(msh.ParIdx, msh.TetTyp, msh.VerTyp, VerTem, &msh);
      }

      tim = GetWallClock() - tim;

      printf("%d steps, average concurency = %g, // running time = %gs\n",
               NmbItr, acc / (2 * NmbItr), tim);
   }
   else if(ParMod == 3)
   {
      if(!AllocAtomicLocks(msh.ParIdx, msh.VerTyp))
      {
         puts("Falied to allocate the vertices' atomic locks");
         exit(1);
      }

      // Perform temperature smoothing steps
      tim = GetWallClock();

      puts("");
      for(i=1;i<=NmbItr;i++)
      {
         acc += LaunchParallel(msh.ParIdx, msh.TetTyp, 0, TetTem, &msh);
         acc += LaunchParallel(msh.ParIdx, msh.TetTyp, 0, LokVer, &msh);
      }

      tim = GetWallClock() - tim;

      printf("%d steps, average concurency = %g, // running time = %gs\n",
               NmbItr, acc / (2 * NmbItr), tim);

      FreeAtomicLocks(msh.ParIdx, msh.VerTyp);
   }

   for(i=1;i<=msh.NmbVer;i++)
      res += msh.VerTem[i];

   printf("Residual = %f\n", res);

   // Stop and free everything
   StopParallel(msh.ParIdx);

   free(msh.TetVer);
   free(msh.TetTem);
   free(msh.VerTem);
   free(msh.VerDeg);

   return(0);
}
