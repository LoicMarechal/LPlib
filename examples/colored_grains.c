

/*----------------------------------------------------------------------------*/
/*                                                                            */
/*                PARALLEL COLORED GRAINS SCHEME USING LPLIB                  */
/*                                                                            */
/*----------------------------------------------------------------------------*/
/*                                                                            */
/*   Description:       handle indirect memory writes with colors and grains  */
/*   Author:            Loic MARECHAL                                         */
/*   Creation date:     sep 09 2024                                           */
/*   Last modification: dec 10 2024                                           */
/*                                                                            */
/*----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------*/
/* Includes                                                                   */
/*----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <libmeshb7.h>
#include "lplib3.h"
#include "lplib3_helpers.h"


/*----------------------------------------------------------------------------*/
/* Macros                                                                     */
/*----------------------------------------------------------------------------*/

#define MIN(a,b)  ((a) < (b) ? (a) : (b))
#define MAX(a,b)  ((a) > (b) ? (a) : (b))
#define NMBITR 10


/*----------------------------------------------------------------------------*/
/* Structure's prototype                                                      */
/*----------------------------------------------------------------------------*/

typedef struct
{
   int      TetTyp, EdgTyp, VerTyp, NmbCol, NmbGrn;
   int      *VerDeg, (*EdgTab)[2], (*TetTab)[4];
   int      (*ColPar)[2], (*VerGrnPar)[2], (*TetGrnPar)[2];
   int64_t  NmbVer, NmbEdg, NmbTet;
   double   *VerSol;
   int64_t  ParIdx;
}MshSct;


/*----------------------------------------------------------------------------*/
/* Procedure called in parallel for each grain                                */
/* increment each edge's vertex degree                                        */
/*----------------------------------------------------------------------------*/

void EdgPar(int BegIdx, int EndIdx, int GrnIdx, MshSct *msh)
{
   int i, j, idx;
   double sol;

   // Loop over the edges, loop over their vertices and increment the degree
   for(i=BegIdx;i<=EndIdx;i++)
      for(j=0;j<2;j++)
      {
         idx = msh->EdgTab[i][j];
         msh->VerDeg[ idx ]++;
         sol = msh->VerSol[ idx ];
         sol = tan(atan(exp(log(sqrt(sol * sol))))) + 1.;
         msh->VerSol[ idx ] = sol;
      }
}


/*----------------------------------------------------------------------------*/
/* Procedure called in parallel for each grain                                */
/* increment each tet's vertex degree                                         */
/*----------------------------------------------------------------------------*/

void TetPar(int BegIdx, int EndIdx, int GrnIdx, MshSct *msh)
{
   int i, j, idx;
   double sol;

   // Loop over the tets, loop over their vertices and increment the degree
   for(i=BegIdx;i<=EndIdx;i++)
      for(j=0;j<4;j++)
      {
         idx = msh->TetTab[i][j];
         msh->VerDeg[ idx ]++;
         sol = msh->VerSol[ idx ];
         sol = tan(atan(exp(log(sqrt(sol * sol))))) + 1.;
         msh->VerSol[ idx ] += sol;
      }
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/

void ClrTet(int BegIdx, int EndIdx, int GrnIdx, MshSct *msh)
{
   int i, j;

   for(i=BegIdx;i<=EndIdx;i++)
      for(j=0;j<4;j++)
         msh->TetTab[i][j] = 0;
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/

void ClrVer(int BegIdx, int EndIdx, int GrnIdx, MshSct *msh)
{
   int i;

   for(i=BegIdx;i<=EndIdx;i++)
      msh->VerDeg[i] = 0;

   for(i=BegIdx;i<=EndIdx;i++)
      msh->VerSol[i] = 0.;
}


/*----------------------------------------------------------------------------*/
/* Comparison of two items for the qsort                                      */
/*----------------------------------------------------------------------------*/

int CmpEdg(const void *a, const void *b)
{
   int *pa = (int *)a, *pb = (int *)b;

   if(pa[0] > pb[0])
      return(1);
   else if(pa[0] < pb[0])
      return(-1);
   else if(pa[1] > pb[1])
      return(1);
   else if(pa[1] < pb[1])
      return(-1);
   else
      return(0);
}


/*----------------------------------------------------------------------------*/
/* Launch NMBITR loops over elements to increment their vertex degree         */
/* and perform some calulation with colored grains or dynamic scheduling      */
/*----------------------------------------------------------------------------*/

int main(int ArgCnt, char **ArgVec)
{
   int      i, j, ref, NmbCpu = 0, ver, dim, ret;
   int      *EdgCol = NULL, *EdgGrn = NULL;
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
   msh.TetTab    = malloc( (msh.NmbTet + 1) * 4 * sizeof(int) );
   msh.ColPar    = malloc( (msh.NmbCol + 1) * 2 * sizeof(int) );
   msh.VerGrnPar = malloc( (msh.NmbGrn + 1) * 2 * sizeof(int) );
   msh.TetGrnPar = malloc( (msh.NmbGrn + 1) * 2 * sizeof(int) );
   msh.VerDeg    = malloc( (msh.NmbVer + 1)     * sizeof(int) );
   msh.VerSol    = malloc( (msh.NmbVer + 1)     * sizeof(double) );

   if(!msh.TetTab || !msh.ColPar || !msh.VerGrnPar || !msh.TetGrnPar | !msh.VerDeg)
   {
      puts("Failed to allocate memory");
      return(3);
   }


   // --------------------------------
   // READ PARTITIONS AND TOUCH MEMORY
   // --------------------------------

   GmfGetBlock(InpMsh, GmfColorPartitions, 1, msh.NmbCol, 0, NULL, NULL,
               GmfIntVec, 2, msh.ColPar[1], msh.ColPar[ msh.NmbCol ]);

   GmfGetBlock(InpMsh, GmfVertexGrainPartitions, 1, msh.NmbGrn, 0, NULL, NULL,
               GmfIntVec, 2, msh.VerGrnPar[1], msh.VerGrnPar[ msh.NmbGrn ]);

   GmfGetBlock(InpMsh, GmfTetrahedronGrainPartitions, 1, msh.NmbGrn, 0, NULL, NULL,
               GmfIntVec, 2, msh.TetGrnPar[1], msh.TetGrnPar[ msh.NmbGrn ]);


   // -----------------------------
   // LPLIB AND PARALLEL DATA SETUP
   // -----------------------------

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


   LaunchColorGrains(msh.ParIdx, msh.VerTyp, ClrVer, &msh);
   LaunchColorGrains(msh.ParIdx, msh.TetTyp, ClrTet, &msh);


   // Read the tets
   GmfGetBlock(InpMsh, GmfTetrahedra, 1, msh.NmbTet, 0, NULL, NULL,
               GmfIntVec, 4, msh.TetTab[1], msh.TetTab[ msh.NmbTet ],
               GmfInt, &ref, &ref);

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
   // Sort the edges against their color, grain and hilbert number
   //qsort(msh.EdgTab[1], msh.NmbEdg, 2 * sizeof(int), CmpEdg);

   printf("Input mesh: nmb vertices = %lld\n", msh.NmbVer);
   printf("Input mesh: nmb colors   = %d\n", msh.NmbCol);
   printf("Input mesh: nmb grains   = %d\n", msh.NmbGrn);
   //printf("Input mesh: nmb edges    = %lld\n", msh.NmbEdg);
   printf("Input mesh: nmb tets     = %lld\n", msh.NmbTet);
   /*
   if(!(msh.EdgTyp = NewType(msh.ParIdx, msh.NmbEdg)))
   {
      puts("Error while creating edges data type.");
      exit(7);
   }

   // Build the edges colored grains partitions
   ret = SetElementsColorGrain(  msh.ParIdx,  msh.VerTyp,  msh.EdgTyp,
                                 2, (int *)msh.EdgTab );

   if(ret)
   {
      printf("SetElementsColorGrain failed with code %d\n", ret);
      exit(9);
   }

*/
   // ---------------------------------
   // MAIN COLORED GRAINS LOOP ON EDGES
   // ---------------------------------

/*   puts("\nColored grains scheduling on edges:");
   tim = GetWallClock();
*/
   // Loop over edges and access vertices
/*   for(i=1;i<=NMBITR;i++)
   {
      ret = LaunchColorGrains(msh.ParIdx, msh.EdgTyp, EdgPar, &msh);

      if(ret)
      {
         printf("LaunchColorGrains exited with error code %d\n", ret);
         exit(10);
      }
   }

   printf("run time = %g\n", GetWallClock() - tim);
*/

   // -----------------------------
   // MAIN DEPENDENCY LOOP ON EDGES
   // -----------------------------

/*   puts("\nDependency loop on edges:");
   BeginDependency(msh.ParIdx, msh.EdgTyp, msh.VerTyp);

   for(i=1;i<=msh.NmbEdg;i++)
      for(j=0;j<2;j++)
         AddDependency(msh.ParIdx, i, msh.EdgTab[i][j]);

   EndDependency(msh.ParIdx, sta);

   tim = GetWallClock();
*/
   // Loop over tets and acces vertices
/*   for(i=1;i<=NMBITR;i++)
      LaunchParallel(msh.ParIdx, msh.EdgTyp, msh.VerTyp, EdgPar, &msh);
   
   printf("Run time = %g\n", GetWallClock() - tim);
*/

   // --------------------------------
   // MAIN COLORED GRAINS LOOP ON TETS
   // --------------------------------

   puts("\nColored grains scheduling on tets:");
   tim = GetWallClock();

   // Loop over edges and access vertices
   for(i=1;i<=NMBITR;i++)
   {
      ret = LaunchColorGrains(msh.ParIdx, msh.TetTyp, TetPar, &msh);

      if(ret)
      {
         printf("LaunchColorGrains exited with error code %d\n", ret);
         exit(10);
      }
   }

   printf("run time = %g\n", GetWallClock() - tim);


   // ----------------------------
   // MAIN DEPENDENCY LOOP ON TETS
   // ----------------------------
   /*
   puts("\nDependency loop on tets:");
   BeginDependency(msh.ParIdx, msh.TetTyp, msh.VerTyp);

   for(i=1;i<=msh.NmbTet;i++)
      for(j=0;j<4;j++)
         AddDependency(msh.ParIdx, i, msh.TetTab[i][j]);

   EndDependency(msh.ParIdx, sta);

   tim = GetWallClock();
*/
   // Loop over tets and acces vertices
/*   for(i=1;i<=NMBITR;i++)
      LaunchParallel(msh.ParIdx, msh.TetTyp, msh.VerTyp, TetPar, &msh);
   
   printf("Run time = %g\n", GetWallClock() - tim);
*/

   // ------------------------
   // STOP AND FREE EVERYTHING
   // ------------------------

   StopParallel(msh.ParIdx);

   free(msh.VerDeg);
   free(msh.VerSol);
   //free(msh.EdgTab);
   free(msh.TetTab);
   free(msh.ColPar);
   free(msh.VerGrnPar);
   free(msh.TetGrnPar);

   return(0);
}
