

/*----------------------------------------------------------------------------*/
/*                                                                            */
/*                PARALLEL COLORED GRAINS SCHEME USING LPlib                  */
/*                                                                            */
/*----------------------------------------------------------------------------*/
/*                                                                            */
/*   Description:       direct and indirect memory writes                     */
/*   Author:            Loic MARECHAL                                         */
/*   Creation date:     sep 09 2024                                           */
/*   Last modification: oct 31 2024                                           */
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
#include "lplib3_helpers.h"
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
   int      TetTyp, EdgTyp, VerTyp, NmbVer, NmbEdg, NmbTet, NmbCol, NmbGrn;
   int      *VerDeg,  (*EdgTab)[2], (*TetTab)[4];
   int      (*ColPar)[2], (*VerGrnPar)[2], (*TetGrnPar)[2];
   int64_t  ParIdx;
}MshSct;


/*----------------------------------------------------------------------------*/
/* Procedure called in parallel for each grain                                */
/*----------------------------------------------------------------------------*/

void ColGrnPar(int BegIdx, int EndIdx, int GrnIdx, MshSct *msh)
{
   int i, j;

   // Loop over tets, loop over their vertices and increment the degree
   for(i=BegIdx;i<=EndIdx;i++)
      for(j=0;j<4;j++)
         msh->VerDeg[ msh->TetTab[i][j] ]++;
}


/*----------------------------------------------------------------------------*/
/* Launch a 1000 loop over tets to increment vertex degree on colored grains  */
/*----------------------------------------------------------------------------*/

int main(int ArgCnt, char **ArgVec)
{
   int      i, ref, NmbCpu = 0, ver, dim, ret, *EdgCol = NULL, *EdgGrn = NULL;
   int64_t  InpMsh;
   float    sta[2], acc = 0;
   double   tim = 0;
   MshSct   msh;


   // --------------------------------------
   // INITIALIZATION, OPENING AND ALLOCATION
   // --------------------------------------

   // Read the number of threads to launch from the command line argument
   if(ArgCnt > 1)
      NmbCpu = atoi(*++ArgVec);

   // Open the input mesh
   if(!(InpMsh = GmfOpenMesh("../sample_meshes/colorgrains.meshb", GmfRead, &ver, &dim)))
   {
      puts("Cannot open the file ../sample_meshes/colorgrains.meshb");
      return(1);
   }

   puts("");
   printf("Input mesh: idx = %lld, version = %d, dimension = %d\n",
            InpMsh, ver, dim);

   // Query the mesh file
   msh.NmbVer = (int)GmfStatKwd(InpMsh, GmfVertices);
   msh.NmbTet = (int)GmfStatKwd(InpMsh, GmfTetrahedra);
   msh.NmbCol = (int)GmfStatKwd(InpMsh, GmfColorPartitions);
   msh.NmbGrn = (int)GmfStatKwd(InpMsh, GmfVertexGrainPartitions);

   if(!msh.NmbVer || !msh.NmbTet || !msh.NmbCol || !msh.NmbGrn)
   {
      puts("Unsupported mesh file");
      return(2);
   }

   // Allocate the memory
   msh.TetTab    = malloc( (msh.NmbTet+1) * 4 * sizeof(int) );
   msh.ColPar    = malloc( (msh.NmbCol+1) * 2 * sizeof(int) );
   msh.VerGrnPar = malloc( (msh.NmbGrn+1) * 2 * sizeof(int) );
   msh.TetGrnPar = malloc( (msh.NmbGrn+1) * 2 * sizeof(int) );
   msh.VerDeg    = malloc( (msh.NmbVer+1)     * sizeof(int) );

   if(!msh.TetTab || !msh.ColPar || !msh.VerGrnPar || !msh.TetGrnPar | !msh.VerDeg)
   {
      puts("Failed to allocate memory");
      return(3);
   }


   // -----------------------------------
   // FILE READING AND DATA PREPROCESSING
   // -----------------------------------

   // Read the tets
   GmfGetBlock(InpMsh, GmfTetrahedra, 1, msh.NmbTet, 0, NULL, NULL,
               GmfIntVec, 4, msh.TetTab[1], msh.TetTab[ msh.NmbTet ],
               GmfInt, &ref, &ref);

   GmfGetBlock(InpMsh, GmfColorPartitions, 1, msh.NmbCol, 0, NULL, NULL,
               GmfIntVec, 2, msh.ColPar[1], msh.ColPar[ msh.NmbCol ]);

   GmfGetBlock(InpMsh, GmfVertexGrainPartitions, 1, msh.NmbGrn, 0, NULL, NULL,
               GmfIntVec, 2, msh.VerGrnPar[1], msh.VerGrnPar[ msh.NmbGrn ]);

   GmfGetBlock(InpMsh, GmfTetrahedronGrainPartitions, 1, msh.NmbGrn, 0, NULL, NULL,
               GmfIntVec, 2, msh.TetGrnPar[1], msh.TetGrnPar[ msh.NmbGrn ]);

   GmfCloseMesh(InpMsh);

   // Extract internal edges
/*   msh.NmbEdg = ParallelBuildEdges( msh.NmbTet, LplTet,
                                    (int *)msh.TetTab, (int **)&msh.EdgTab );

   if(!msh.NmbEdg)
   {
      puts("Failed to extract internal edges");
      exit(4);
   }
*/
   // Assign a color and a grain to each edge

   // And renumber the edges against the color, grain and hilbert code

   printf("Input mesh: nmb vertices = %d\n", msh.NmbVer);
   printf("Input mesh: nmb colors   = %d\n", msh.NmbCol);
   printf("Input mesh: nmb grains   = %d\n", msh.NmbGrn);
   //printf("Input mesh: nmb edges    = %d\n", msh.NmbEdg);
   printf("Input mesh: nmb tets     = %d\n", msh.NmbTet);


   // ----------------------------
   // LPLIB AN PARALLEL DATA SETUP
   // ----------------------------

   // Initialize the LPlib and setup the data types
   if(!(msh.ParIdx = InitParallel(NmbCpu)))
   {
      puts("Error initializing the LPlib.");
      exit(5);
   }

   if(!(msh.TetTyp = NewType(msh.ParIdx, msh.NmbTet)))
   {
      puts("Error while creating tetrahedra data type.");
      exit(6);
   }

/*   if(!(msh.EdgTyp = NewType(msh.ParIdx, msh.NmbEdg)))
   {
      puts("Error while creating edges data type.");
      exit(7);
   }
*/
   if(!(msh.VerTyp = NewType(msh.ParIdx, msh.NmbVer)))
   {
      puts("Error while creating vertices data type.");
      exit(8);
   }

   puts("");
   printf("TetTyp = %d, VerTyp = %d, NmbCpu = %d\n",
         msh.TetTyp, msh.VerTyp, NmbCpu);

   // Send vertices color and grain to the LPlib
   SetColorGrains(msh.ParIdx, msh.VerTyp, msh.NmbCol, (int *)msh.ColPar,
                  msh.NmbGrn, (int *)msh.VerGrnPar);

   // Send tets color and grain to the LPlib
   SetColorGrains(msh.ParIdx, msh.TetTyp, msh.NmbCol, (int *)msh.ColPar,
                  msh.NmbGrn, (int *)msh.TetGrnPar);

/*   ret = SetElementsColorGrain(  msh.ParIdx,  msh.VerTyp,  msh.EdgTyp,
                                 2, (int *)msh.EdgTab, &EdgCol, &EdgGrn );

   if(ret || !EdgCol || !EdgGrn)
   {
      printf(  "SetElementsColorGrain failed with code %d, EdgCol=%p, EdgGrn=%p\n",
               ret, EdgCol, EdgGrn );
      exit(9);
   }*/


   // -------------------
   // MAIN PARALLEL LOOPS
   // -------------------

   // Loop over tets and acces vertices
   for(i=1;i<=1000;i++)
   {
      ret = LaunchColorGrains(msh.ParIdx, msh.TetTyp, ColGrnPar, &msh);

      if(ret)
      {
         printf("LaunchColorGrains exited with error code %d\n", ret);
         exit(9);
      }
   }

   for(i=1;i<=10;i++)
      printf("Ver %d, deg %d\n", i, msh.VerDeg[i]);


   // ------------------------
   // STOP AND FREE EVERYTHING
   // ------------------------

   StopParallel(msh.ParIdx);

   free(msh.TetTab);
   free(msh.ColPar);
   free(msh.VerGrnPar);
   free(msh.TetGrnPar);

   return(0);
}
