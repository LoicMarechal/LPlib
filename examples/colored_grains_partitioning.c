

/*----------------------------------------------------------------------------*/
/*                                                                            */
/*                PARALLEL COLORED GRAINS SCHEME USING LPLIB                  */
/*                                                                            */
/*----------------------------------------------------------------------------*/
/*                                                                            */
/*   Description:       handle indirect memory writes with colors and grains  */
/*   Author:            Loic MARECHAL                                         */
/*   Creation date:     sep 09 2024                                           */
/*   Last modification: aug 20 2025                                           */
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
#include <libmeshb8.h>
#include "lplib4.h"
#include "lplib4_helpers.h"


/*----------------------------------------------------------------------------*/
/* Macros                                                                     */
/*----------------------------------------------------------------------------*/

#define MIN(a,b)  ((a) < (b) ? (a) : (b))
#define MAX(a,b)  ((a) > (b) ? (a) : (b))
#define NMBITR 100


/*----------------------------------------------------------------------------*/
/* Structure's prototype                                                      */
/*----------------------------------------------------------------------------*/

typedef struct
{
   int      TetTyp, EdgTyp, VerTyp, NmbCol, NmbGrn, *VerDeg;
   itg      (*EdgTab)[2], (*TetTab)[4], *TetRef;
   int      (*ColPar)[2], (*VerGrnPar)[2], (*TetGrnPar)[2];
   int64_t  NmbVer, NmbEdg, NmbTet;
   double   *VerSol, (*VerCrd)[3];
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
         msh->VerSol[ idx ] = sqrt(sol * sol + 1.) + 1.;
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
         msh->VerSol[ idx ] = sqrt(sol * sol + 1.) + 1.;
      }
}


void SetTetRef(int BegIdx, int EndIdx, int GrnIdx, MshSct *msh)
{
   int i;

   for(i=BegIdx;i<=EndIdx;i++)
      msh->TetRef[i] = GrnIdx;
}

/*----------------------------------------------------------------------------*/
/* Touch tet's memory: for ccNUMA experiments only                            */
/*----------------------------------------------------------------------------*/

void ClrTet(int BegIdx, int EndIdx, int GrnIdx, MshSct *msh)
{
   int i, j;

   for(i=BegIdx;i<=EndIdx;i++)
      for(j=0;j<4;j++)
         msh->TetTab[i][j] = 0;
}


/*----------------------------------------------------------------------------*/
/* Reset vertex degree and solution                                           */
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
/* Comparison between two items for the qsort                                 */
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
   int         i, j, ref, NmbCpu = 0, NmbGrn, ver, dim, ret;
   int         *EdgCol = NULL, *EdgGrn = NULL;
   int64_t     InpMsh, OutMsh, DegTot = 0;
   float       sta[2], acc = 0;
   double      tim = 0;
   MshSct      msh;
   LplSct      *RenNfo;


   // --------------------------------------
   // INITIALIZATION, OPENING AND ALLOCATION
   // --------------------------------------

   // Read the number of threads to launch from the command line argument
   if(ArgCnt == 3)
   {
      NmbCpu = atoi(*++ArgVec);
      NmbGrn = atoi(*++ArgVec);
   }
   else
   {
      puts("colored_grains_partitioning   nbThreads   nbGrains");
      exit(0);
   }

   // Open the input mesh
   if(!(InpMsh = GmfOpenMesh("../sample_meshes/tet.meshb", GmfRead, &ver, &dim)))
   {
      puts("Cannot open the file ../sample_meshes/tet.meshb");
      return(1);
   }

   puts("");
   printf("Input mesh: idx = %lld, version = %d, dimension = %d\n",
            InpMsh, ver, dim);

   // Query the mesh file
   msh.NmbVer = (int)GmfStatKwd(InpMsh, GmfVertices);
   msh.NmbTet = (int)GmfStatKwd(InpMsh, GmfTetrahedra);

   if(!msh.NmbVer || !msh.NmbTet)
   {
      puts("Unsupported mesh file");
      return(2);
   }

   // Allocate the memory
   puts("Allocate Mesh");
   msh.TetTab = malloc( (msh.NmbTet + 1) * 4 * sizeof(itg) );
   msh.VerCrd = malloc( (msh.NmbVer + 1) * 3 * sizeof(double) );
   msh.VerDeg = calloc( (msh.NmbVer + 1)     , sizeof(int) );
   msh.VerSol = calloc( (msh.NmbVer + 1)     , sizeof(double) );

   if(!msh.TetTab || !msh.VerSol || !msh.VerDeg || !msh.VerCrd)
   {
      puts("Failed to allocate memory");
      return(3);
   }


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

   // Read the vertices
   puts("Read vertices");

   GmfGetBlock(InpMsh, GmfVertices, 1, msh.NmbVer, 0, NULL, NULL,
               GmfDoubleVec, 3, msh.VerCrd[1], msh.VerCrd[ msh.NmbVer ],
               GmfInt, &ref, &ref);

   // Read the tets
   puts("Read tets");

#ifdef INT64
   GmfGetBlock(InpMsh, GmfTetrahedra, 1, msh.NmbTet, 0, NULL, NULL,
               GmfLongVec, 4, msh.TetTab[1], msh.TetTab[ msh.NmbTet ],
               GmfLong, &ref, &ref);
#else
   GmfGetBlock(InpMsh, GmfTetrahedra, 1, msh.NmbTet, 0, NULL, NULL,
               GmfIntVec, 4, msh.TetTab[1], msh.TetTab[ msh.NmbTet ],
               GmfInt, &ref, &ref);
#endif

   GmfCloseMesh(InpMsh);

   // Extract internal edges
   puts("Build edges");
   msh.NmbEdg = ParallelBuildEdges( msh.NmbTet, LplTet,
                                    (itg *)msh.TetTab, (itg **)&msh.EdgTab );

   if(!msh.NmbEdg)
   {
      puts("Failed to extract internal edges");
      exit(4);
   }

   printf("Nmb edges extracted: %lld\n",  msh.NmbEdg);

   // Sort the edges against their color, grain and hilbert number
   puts("Sort edges");
   qsort(msh.EdgTab[1], msh.NmbEdg, 2 * sizeof(int), CmpEdg);

   printf("Input mesh: nmb vertices = %lld\n", msh.NmbVer);
   printf("Input mesh: nmb tets     = %lld\n", msh.NmbTet);

   if(!(msh.EdgTyp = NewType(msh.ParIdx, msh.NmbEdg)))
   {
      puts("Error while creating edges data type.");
      exit(7);
   }

   printf("VerTyp = %d, EdgTyp = %d, TetTyp = %d, NmbCpu = %d\n",
         msh.VerTyp, msh.EdgTyp, msh.TetTyp, NmbCpu);

   if(!(OutMsh = GmfOpenMesh("/tmp/ini.meshb", GmfWrite, 2, 3)))
   {
      puts("Cannot open the file /tmp/ini.meshb");
      return(1);
   }

   GmfSetKwd(OutMsh, GmfVertices, msh.NmbVer);
   GmfSetBlock(OutMsh, GmfVertices, 1, msh.NmbVer, 0, NULL, NULL,
               GmfDoubleVec, 3, msh.VerCrd[1], msh.VerCrd[ msh.NmbVer ],
               GmfInt, &ref, &ref);

   GmfSetKwd(OutMsh, GmfTetrahedra, msh.NmbTet);
   GmfSetBlock(OutMsh, GmfTetrahedra, 1, msh.NmbTet, 0, NULL, NULL,
               GmfIntVec, 4, msh.TetTab[1], msh.TetTab[ msh.NmbTet ],
               GmfInt, &ref, &ref);

   GmfCloseMesh(OutMsh);

   RenNfo = MeshRenumbering(  msh.ParIdx, NmbGrn, LplHilbert, 0, 3,
                              LplVer, msh.VerTyp, msh.NmbVer, msh.VerCrd, NULL,
                              LplEdg, msh.EdgTyp, msh.NmbEdg, msh.EdgTab, NULL,
                              LplTet, msh.TetTyp, msh.NmbTet, msh.TetTab, NULL,
                              LplMax );


   // -----------------------------
   // MAIN DEPENDENCY LOOP ON EDGES
   // -----------------------------

/*   puts("Set dependency edges - >vertices");
   BeginDependency(msh.ParIdx, msh.EdgTyp, msh.VerTyp);

   for(i=1;i<=msh.NmbEdg;i++)
      for(j=0;j<2;j++)
         AddDependency(msh.ParIdx, i, msh.EdgTab[i][j]);

   EndDependency(msh.ParIdx, sta);
   printf("collisions: %g / %g\n", sta[0], sta[1]);
*/

   // -----------------------------
   // MAIN DEPENDENCY LOOP ON EDGES
   // -----------------------------

/*   puts("\nDependency loop on edges:");
   tim = GetWallClock();

   for(i=1;i<=NMBITR;i++)
      LaunchParallel(msh.ParIdx, msh.EdgTyp, msh.VerTyp, EdgPar, &msh);

   printf("Run time = %g\n", GetWallClock() - tim);

   for(i=1;i<=msh.NmbVer;i++)
   {
      DegTot += msh.VerDeg[i];
      msh.VerDeg[i] = 0;
   }

   printf("Vertex total degree = %lld\n", DegTot);
*/


   msh.TetRef = malloc( (msh.NmbTet + 1) * sizeof(int) );
   LaunchColorGrains(msh.ParIdx, LplTet, SetTetRef, &msh);

   if(!(OutMsh = GmfOpenMesh("/tmp/col.meshb", GmfWrite, 2, 3)))
   {
      puts("Cannot open the file /tmp/col.meshb");
      return(1);
   }

   GmfSetKwd(OutMsh, GmfVertices, msh.NmbVer);
   GmfSetBlock(OutMsh, GmfVertices, 1, msh.NmbVer, 0, NULL, NULL,
               GmfDoubleVec, 3, msh.VerCrd[1], msh.VerCrd[ msh.NmbVer ],
               GmfInt, &ref, &ref);

   GmfSetKwd(OutMsh, GmfTetrahedra, msh.NmbTet);
   GmfSetBlock(OutMsh, GmfTetrahedra, 1, msh.NmbTet, 0, NULL, NULL,
               GmfIntVec, 4, msh.TetTab[1], msh.TetTab[ msh.NmbTet ],
               GmfInt, &msh.TetRef[1], &msh.TetRef[ msh.NmbTet ]);

   GmfCloseMesh(OutMsh);
   free(msh.TetRef);
   //exit(0);

   // ---------------------------------
   // MAIN COLORED GRAINS LOOP ON EDGES
   // ---------------------------------

   puts("\nColored grains scheduling on edges:");
   tim = GetWallClock();

   // Loop over edges and access vertices
   for(i=1;i<=NMBITR;i++)
      LaunchColorGrains(msh.ParIdx, LplEdg, EdgPar, &msh);

   printf("Run time = %g\n\n", GetWallClock() - tim);

   DegTot = 0;

   for(i=1;i<=msh.NmbVer;i++)
   {
      DegTot += msh.VerDeg[i];
      msh.VerDeg[i] = 0;
   }

   printf("Vertex total degree = %lld\n", DegTot);


   // --------------------------------
   // MAIN COLORED GRAINS LOOP ON TETS
   // --------------------------------

/*   puts("\nColored grains scheduling on tets:");
   tim = GetWallClock();
  */    
   // Loop over tets and access vertices
/*   for(i=1;i<=NMBITR;i++)
      LaunchColorGrains(msh.ParIdx, LplTet, TetPar, &msh);

   printf("Run time = %g\n", GetWallClock() - tim);

   DegTot = 0;

   for(i=1;i<=msh.NmbVer;i++)
      DegTot += msh.VerDeg[i];

   printf("Vertex total degree = %lld\n", DegTot);
*/

   // ------------------------
   // STOP AND FREE EVERYTHING
   // ------------------------

   StopParallel(msh.ParIdx);

   free(msh.VerDeg);
   free(msh.VerCrd);
   free(msh.VerSol);
   free(msh.EdgTab);
   free(msh.TetTab);

   return(0);
}
