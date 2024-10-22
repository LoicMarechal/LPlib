

/*----------------------------------------------------------------------------*/
/*                                                                            */
/*                PARALLEL COLORED GRAINS SCHEME USING LPlib                  */
/*                                                                            */
/*----------------------------------------------------------------------------*/
/*                                                                            */
/*   Description:       direct and indirect memory writes                     */
/*   Author:            Loic MARECHAL                                         */
/*   Creation date:     sep 09 2024                                           */
/*   Last modification: oct 22 2024                                           */
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
/* Macros                                                                     */
/*----------------------------------------------------------------------------*/

#define MIN(a,b)  ((a) < (b) ? (a) : (b))
#define MAX(a,b)  ((a) > (b) ? (a) : (b))


/*----------------------------------------------------------------------------*/
/* Structure's prototype                                                      */
/*----------------------------------------------------------------------------*/

typedef struct
{
   int      TetTyp, VerTyp, NmbVer, NmbTet, *VerDeg, (*TetVer)[4];
   int      NmbCol, NmbGrn, (*ColPar)[2], (*VerGrnPar)[2], (*TetGrnPar)[2];
   int64_t  ParIdx;
}MshSct;


/*----------------------------------------------------------------------------*/
/* Procedure called in parallel for each grain                                */
/*----------------------------------------------------------------------------*/

void ColGrnPar(int BegIdx, int EndIdx, int GrnIdx, MshSct *msh)
{
   int i, j;

   //printf("  grain %d, idx %d -> %d\n", GrnIdx, BegIdx, EndIdx);

   for(i=BegIdx;i<=EndIdx;i++)
      for(j=0;j<4;j++)
         msh->VerDeg[ msh->TetVer[i][j] ]++;
}


/*----------------------------------------------------------------------------*/
/* Setup LPlib and launch 100 temperature soomthing steps                     */
/*----------------------------------------------------------------------------*/

int main(int ArgCnt, char **ArgVec)
{
   int      i, ref, NmbCpu = 0, ver, dim;
   int64_t  InpMsh;
   float    sta[2], acc = 0;
   double   tim = 0;
   MshSct   msh;

   // Read the number of threads to launch from the command line argument
   if(ArgCnt > 1)
      NmbCpu = atoi(*++ArgVec);

   // Open the input mesh
   if(!(InpMsh = GmfOpenMesh("../sample_meshes/colorgrains.meshb", GmfRead, &ver, &dim)))
      return(1);

   puts("");
   printf("Input mesh: idx = %lld, version = %d, dimension = %d\n",
            InpMsh, ver, dim);

   msh.NmbVer = (int)GmfStatKwd(InpMsh, GmfVertices);
   msh.NmbTet = (int)GmfStatKwd(InpMsh, GmfTetrahedra);
   msh.NmbCol = (int)GmfStatKwd(InpMsh, GmfColorPartitions);
   msh.NmbGrn = (int)GmfStatKwd(InpMsh, GmfVertexGrainPartitions);

   printf("Input mesh: nmb vertices = %d\n", msh.NmbVer);
   printf("Input mesh: nmb colors   = %d\n", msh.NmbCol);
   printf("Input mesh: nmb grains   = %d\n", msh.NmbGrn);
   printf("Input mesh: nmb tets     = %d\n", msh.NmbTet);

   // Allocate the memory
   msh.TetVer = malloc( (msh.NmbTet+1) * 4 * sizeof(int) );
   msh.ColPar = malloc( (msh.NmbCol+1) * 2 * sizeof(int) );
   msh.VerGrnPar = malloc( (msh.NmbGrn+1) * 2 * sizeof(int) );
   msh.TetGrnPar = malloc( (msh.NmbGrn+1) * 2 * sizeof(int) );
   msh.VerDeg = malloc( (msh.NmbVer+1) * sizeof(int) );

   // Read the tets
   puts("read tet");
   GmfGetBlock(InpMsh, GmfTetrahedra, 1, msh.NmbTet, 0, NULL, NULL,
               GmfIntVec, 4, msh.TetVer[1], msh.TetVer[ msh.NmbTet ],
               GmfInt, &ref, &ref);

   puts("read col");
   GmfGetBlock(InpMsh, GmfColorPartitions, 1, msh.NmbCol, 0, NULL, NULL,
               GmfIntVec, 2, msh.ColPar[1], msh.ColPar[ msh.NmbCol ]);

   puts("read grain");
   GmfGetBlock(InpMsh, GmfVertexGrainPartitions, 1, msh.NmbGrn, 0, NULL, NULL,
               GmfIntVec, 2, msh.VerGrnPar[1], msh.VerGrnPar[ msh.NmbGrn ]);

   GmfGetBlock(InpMsh, GmfTetrahedronGrainPartitions, 1, msh.NmbGrn, 0, NULL, NULL,
               GmfIntVec, 2, msh.TetGrnPar[1], msh.TetGrnPar[ msh.NmbGrn ]);

   puts("done");
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
   printf("TetTyp = %d, VerTyp = %d, NmbCpu = %d\n",
         msh.TetTyp, msh.VerTyp, NmbCpu);

   printf("\nNmb grains = %d\n", msh.NmbGrn);
   /*
   puts("Vertex grains:\n");
   for(i=1;i<=msh.NmbGrn;i++)
      printf("grain %d: ver %d -> %d\n", i, msh.VerGrnPar[i][0], msh.VerGrnPar[i][1]);

   puts("Tet grains:\n");
   for(i=1;i<=msh.NmbGrn;i++)
      printf("grain %d: tet %d -> %d\n", i, msh.TetGrnPar[i][0], msh.TetGrnPar[i][1]);
*/
   printf("\nNmb colors = %d\n", msh.NmbCol);
   /*
   for(i=1;i<=msh.NmbCol;i++)
      printf("color %d: grain %2d -> %2d\n", i, msh.ColPar[i][0], msh.ColPar[i][1]);
*/

   SetColorGrains(msh.ParIdx, msh.VerTyp, msh.NmbCol, (int *)msh.ColPar,
                  msh.NmbGrn, (int *)msh.VerGrnPar);

   SetColorGrains(msh.ParIdx, msh.TetTyp, msh.NmbCol, (int *)msh.ColPar,
                  msh.NmbGrn, (int *)msh.TetGrnPar);

   for(i=1;i<=1000;i++)
      LaunchParallel(msh.ParIdx, msh.TetTyp, ColorGrainScheduling, ColGrnPar, &msh);

   for(i=1;i<=10;i++)
      printf("Ver %d, deg %d\n", i, msh.VerDeg[i]);

   // Stop and free everything
   StopParallel(msh.ParIdx);

   free(msh.TetVer);
   free(msh.ColPar);
   free(msh.VerGrnPar);
   free(msh.TetGrnPar);

   return(0);
}
