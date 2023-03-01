

/*----------------------------------------------------------------------------*/
/*                                                                            */
/*          ADAPT BLOCKS SIZES FOR INDIRECT MEMORY WRITES USING LPlib         */
/*                                                                            */
/*----------------------------------------------------------------------------*/
/*                                                                            */
/*   Description:       adaptive blocking direct and indirect memory writes   */
/*   Author:            Loic MARECHAL                                         */
/*   Creation date:     feb 28 2023                                           */
/*   Last modification: mar 01 2023                                           */
/*                                                                            */
/*----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------*/
/* Includes                                                                   */
/*----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <stdint.h>
#include "lplib3.h"
#include "libmeshb7.h"


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

      msh->TetTem[i] = t;
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
            msh->TetTem[i] / msh->VerDeg[ msh->TetVer[i][j] ];
}


/*----------------------------------------------------------------------------*/
/* Setup LPlib and launch 100 temperature soomthing steps                     */
/*----------------------------------------------------------------------------*/

int main(int ArgCnt, char **ArgVec)
{
   int      i, j, res, NmbCpu = GetNumberOfCores(), ver, dim, ref, NmbItr=10;
   int64_t  InpMsh;
   float    sta[2], acc=0;
   double   tim=0, OldTim;
   MshSct   msh;

   // Read the number of threads to launch from the command line argument
   if(ArgCnt > 1)
      NmbCpu = atoi(*++ArgVec);


   // ************
   // MESH READING
   // ************

   // Open the input mesh
   if(!(InpMsh = GmfOpenMesh("../sample_meshes/tet.meshb", GmfRead, &ver, &dim)))
      return(1);

   // Print mesh stats
   puts("-----------------------------------------------------------");
   printf(" Input mesh: idx = %8d, version = %1d, dimension = %1d\n",
            (int)InpMsh, ver, dim);

   msh.NmbVer = (int)GmfStatKwd(InpMsh, GmfVertices);
   printf(" Input mesh: nmb vertices  = %8d\n", msh.NmbVer);

   msh.NmbTet = (int)GmfStatKwd(InpMsh, GmfTetrahedra);
   printf(" Input mesh: nmb tetrahedra = %8d\n", msh.NmbTet);
   puts("-----------------------------------------------------------\n");

   // Allocate the memory
   msh.TetVer = malloc( (msh.NmbTet+1) * 4 * sizeof(int) );
   msh.TetTem = malloc( (msh.NmbTet+1) * sizeof(double) );
   msh.VerTem = malloc( (msh.NmbVer+1) * sizeof(double) );
   msh.VerDeg = calloc( (msh.NmbVer+1), sizeof(int) );

   // Read the tets
   GmfGetBlock(InpMsh, GmfTetrahedra, 1, msh.NmbTet, 0, NULL, NULL,
               GmfIntVec, 4, msh.TetVer[1], msh.TetVer[ msh.NmbTet ],
               GmfInt, &ref, &ref);

   GmfCloseMesh(InpMsh);


   // ********************
   // MULTITHREADING SETUP
   // ********************

   puts("-----------------------------------------------------------");

   // Initialize the LPlib and setup the data types
   if(!(msh.ParIdx = InitParallel(NmbCpu)))
   {
      puts("Error initializing the LPlib.");
      exit(1);
   }

   // Disable the block sorting so we can modify the number of blocks
   SetExtendedAttributes(msh.ParIdx, DisableBlockSorting);

   // Set the default number of blocks with a very high value
   SetExtendedAttributes(msh.ParIdx, SetSmallBlock, 1024);
   SetExtendedAttributes(msh.ParIdx, SetDependencyBlock, 1024);

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

   printf(" TetTyp = %d, VerTyp = %d, NmbCpu = %d\n",
         msh.TetTyp, msh.VerTyp, NmbCpu);

   // Setup dependencies between tets and vertices 
   // and compute vertices' degree on the flight
   BeginDependency(msh.ParIdx, msh.TetTyp, msh.VerTyp);

   for(i=1;i<=msh.NmbTet;i++)
      for(j=0;j<4;j++)
      {
         AddDependency(msh.ParIdx, i, msh.TetVer[i][j]);
         msh.VerDeg[ msh.TetVer[i][j] ]++;
      }

   EndDependency(msh.ParIdx, sta);

   printf(" Average dependencies = %3.2f%%, max dependencies = %3.2f%%\n",
            sta[0], sta[1]);

   puts("-----------------------------------------------------------\n");


   // *************
   // REFERENCE RUN
   // *************

   // Initialize the vertices' temperature with some crap values
   for(i=1;i<=msh.NmbVer;i++)
#ifdef WIN32
      msh.VerTem[i] = rand();
#else
      msh.VerTem[i] = random();
#endif

   // Perform temperature smoothing steps
   tim = GetWallClock();

   for(i=1;i<=NmbItr;i++)
   {
      LaunchParallel(msh.ParIdx, msh.TetTyp, 0, (void *)TetTem, (void *)&msh);
      acc += LaunchParallel(msh.ParIdx, msh.TetTyp, msh.VerTyp, (void *)VerTem, (void *)&msh);
   }

   tim = GetWallClock() - tim;

   puts("-----------------------------------------------------------");
   printf(" Init run: %d steps, concurency = %3.2f, time = %3.2fs\n",
            NmbItr, 1. + .5 * (acc / NmbItr), tim);
   puts("-----------------------------------------------------------\n");


   // ********************************
   // ADAPT THE NUMBER OF SMALL BLOCKS
   // ********************************

   // Perform the dependency loop only
   tim = GetWallClock();
   acc = 0;
   for(i=1;i<=NmbItr;i++)
      acc += LaunchParallel(  msh.ParIdx, msh.TetTyp, msh.VerTyp,
                              (void *)VerTem, (void *)&msh );
   tim = GetWallClock() - tim;
   puts("-----------------------------------------------------------");
   puts(" Halve small blocks as long as it speeds things up");
   printf(" reference run:   concurency: %3.2f, time: %2.2fs\n", acc / NmbItr, tim);

   // Halve the number of small blocks as long as the run time goes down
   do
   {
      OldTim = tim;

      if(!(res = HalveSmallBlocks(msh.ParIdx, msh.TetTyp, msh.VerTyp)))
         break;

      // Perform temperature smoothing steps
      tim = GetWallClock();
      acc = 0;
      for(i=1;i<=NmbItr;i++)
         acc += LaunchParallel(  msh.ParIdx, msh.TetTyp, msh.VerTyp,
                                 (void *)VerTem, (void *)&msh );
      tim = GetWallClock() - tim;

      printf(" blocks: %8d concurency: %3.2f, time: %2.2fs\n",
               res, acc / NmbItr, tim );
   }while(tim < OldTim);

   puts("-----------------------------------------------------------\n");

   puts("-----------------------------------------------------------");
   puts(" Halve dependency blocks as long as it speeds things up");
   do
   {
      OldTim = tim;

      if(!(res = HalveDependencyBlocks(msh.ParIdx, msh.TetTyp, msh.VerTyp)))
         break;

      // Perform temperature smoothing steps
      tim = GetWallClock();
      acc = 0;
      for(i=1;i<=NmbItr;i++)
         acc += LaunchParallel(  msh.ParIdx, msh.TetTyp, msh.VerTyp,
                                 (void *)VerTem, (void *)&msh );
      tim = GetWallClock() - tim;

      printf(" blocks: %8d concurency: %3.2f, time: %2.2fs\n",
               res, acc / NmbItr, tim );
   }while(tim < OldTim);

   puts("-----------------------------------------------------------\n");


   // *******************
   // FINAL OPTIMIZED RUN
   // *******************

   // Perform temperature smoothing steps
   tim = GetWallClock();
   acc = 0;

   for(i=1;i<=NmbItr;i++)
   {
      LaunchParallel(msh.ParIdx, msh.TetTyp, 0, (void *)TetTem, (void *)&msh);
      acc += LaunchParallel(msh.ParIdx, msh.TetTyp, msh.VerTyp, (void *)VerTem, (void *)&msh);
   }

   tim = GetWallClock() - tim;

   puts("-----------------------------------------------------------");
   printf(" Final run: %d steps, concurency = %3.2f, time = %3.2fs\n",
            NmbItr, 1. + .5 * (acc / NmbItr), tim);
   puts("-----------------------------------------------------------\n");

   // Stop and free everything
   StopParallel(msh.ParIdx);

   free(msh.TetVer);
   free(msh.TetTem);
   free(msh.VerTem);
   free(msh.VerDeg);

   return(0);
}
