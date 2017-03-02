

/*----------------------------------------------------------------------------*/
/*                                                                            */
/*                PARALLEL INDIRECT MEMORY WRITES USING LPlib                 */
/*                                                                            */
/*----------------------------------------------------------------------------*/
/*                                                                            */
/*   Description:       direct and indirect memory writes                     */
/*   Author:            Loic MARECHAL                                         */
/*   Creation date:     jan 15 2015                                           */
/*   Last modification: feb 01 2017                                           */
/*                                                                            */
/*----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------*/
/* Includes                                                                   */
/*----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include "lplib3.h"
#include "libmeshb7.h"


/*----------------------------------------------------------------------------*/
/* Structure's prototype                                                      */
/*----------------------------------------------------------------------------*/

typedef struct
{
   int TetTyp, VerTyp, NmbVer, NmbTet, (*TetVer)[4], *VerDeg;
   long long ParIdx;
   double *VerTem, *TetTem;
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
         msh->VerTem[ msh->TetVer[i][j] ] += \
            msh->TetTem[i] / msh->VerDeg[ msh->TetVer[i][j] ];
}


/*----------------------------------------------------------------------------*/
/* Setup LPlib and launch 100 temperature soomthing steps                     */
/*----------------------------------------------------------------------------*/

int main(int ArgCnt, char **ArgVec)
{
   int i, j, NmbCpu=GetNumberOfCores(), ver, dim, ref, NmbItr=100;
   long long InpMsh;
   float sta[2], acc=0;
   double tim=0;
   MshSct msh;

   // Read the number of threads to launch from the command line argument
   if(ArgCnt > 1)
      NmbCpu = atoi(*++ArgVec);

   // Open the input mesh
   if(!(InpMsh = GmfOpenMesh("../sample_meshes/tet.meshb", GmfRead, &ver, &dim)))
      return(1);

   puts("");
   printf("Input mesh : idx = %lld, version = %d, dimension = %d\n", \
            InpMsh, ver, dim);

   msh.NmbVer = GmfStatKwd(InpMsh, GmfVertices);
   printf("Input mesh : nmb vertices = %d\n", msh.NmbVer);

   msh.NmbTet = GmfStatKwd(InpMsh, GmfTetrahedra);
   printf("Input mesh : nmb tets = %d\n", msh.NmbTet);

   // Allocate the memory
   msh.TetVer = malloc( (msh.NmbTet+1) * 4 * sizeof(int) );
   msh.TetTem = malloc( (msh.NmbTet+1) * sizeof(double) );
   msh.VerTem = malloc( (msh.NmbVer+1) * sizeof(double) );
   msh.VerDeg = calloc( (msh.NmbVer+1), sizeof(int) );

   // Read the mesh
   GmfGotoKwd(InpMsh, GmfTetrahedra);

   for(i=1;i<=msh.NmbTet;i++)
      GmfGetLin(  InpMsh, GmfTetrahedra, &msh.TetVer[i][0], &msh.TetVer[i][1], \
                  &msh.TetVer[i][2], &msh.TetVer[i][3], &ref );

   GmfCloseMesh(InpMsh);

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
   printf("TetTyp = %d, VerTyp = %d, NmbCpu = %d\n", \
         msh.TetTyp, msh.VerTyp, NmbCpu);

   /* Setup dependencies between tets and vertices 
      and compute vertices' degree on the flight */

   BeginDependency(msh.ParIdx, msh.TetTyp, msh.VerTyp);

   for(i=1;i<=msh.NmbTet;i++)
      for(j=0;j<4;j++)
      {
         AddDependency(msh.ParIdx, i, msh.TetVer[i][j]);
         msh.VerDeg[ msh.TetVer[i][j] ]++;
      }

   EndDependency(msh.ParIdx, sta);

   printf("dependencies stats : average = %g %%, maximum = %g %%\n", \
            sta[0], sta[1]);

   // Initialize the vertices' temperature with some crap values
   for(i=1;i<=msh.NmbVer;i++)
      msh.VerTem[i] = random();

   // Perform temperature smoothing steps
   tim = GetWallClock();

   puts("");
   for(i=1;i<=NmbItr;i++)
   {
      if(!(acc += LaunchParallel(msh.ParIdx, msh.TetTyp, msh.VerTyp, \
                                 (void *)TetTem, (void *)&msh)))
      {
         puts("Error while launching the parallel loop TetTem.");
         exit(1);
      }

      if(!(acc += LaunchParallel(msh.ParIdx, msh.TetTyp, 0, \
                                 (void *)VerTem, (void *)&msh)))
      {
         puts("Error while launching the parallel loop VerTem.");
         exit(1);
      }
   }

   tim = GetWallClock() - tim;

   printf(" %d steps, average concurency = %g, // running time = %gs\n", \
         NmbItr, acc / (2 * NmbItr), tim);

   // Stop and free everything
   StopParallel(msh.ParIdx);

   free(msh.TetVer);
   free(msh.TetTem);
   free(msh.VerTem);
   free(msh.VerDeg);

   return(0);
}
