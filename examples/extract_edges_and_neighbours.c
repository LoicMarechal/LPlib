

/*----------------------------------------------------------------------------*/
/*                                                                            */
/*    PARALLEL NEIGHBOURS AND EDGES EXTRACTIONUSING THE LPLIB HELPERS         */
/*                                                                            */
/*----------------------------------------------------------------------------*/
/*                                                                            */
/*   Description:       extract inner edges and face neighbours               */
/*                      from a volume only tetrahedral mesh                   */
/*   Author:            Loic MARECHAL                                         */
/*   Creation date:     mar 30 2026                                           */
/*   Last modification: mar 30 2026                                           */
/*                                                                            */
/*----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------*/
/* Includes                                                                   */
/*----------------------------------------------------------------------------*/

#ifdef _WIN32
#define _CRT_SECURE_NO_WARNINGS
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libmeshb8.h>
#include <lplib4.h>
#include <lplib4_helpers.h>


/*----------------------------------------------------------------------------*/
/* Read the volume, extract the surface and write the mesh                    */
/*----------------------------------------------------------------------------*/

int main(int ArgCnt, char **ArgVec)
{
   char     *PtrArg, *TmpStr, InpNam[1000], OutNam[1000], *VoyTab;
   int      NmbEdg, NmbTet, *TetTab, *NgbTab, *EdgTab;
   int      i, j, ref, ver, dim, NmbCpu = 0;
   int64_t  InpMsh;
   double   timer;


   // 
   // Command line parsing
   // 

   if(ArgCnt == 1)
   {
      puts("\nextract_edges_and_neighbours v1.00 March 30 2026   Loic MARECHAL / INRIA");
      puts(" Usage      : extract_edges_and_neighbours -in volume_mesh");
      puts(" -in name   : name of the input tetrahedral-only mesh");
      puts(" -nproc n   : n is the number of threads to be launched (default = all available threads)\n");
      exit(0);
   }

   for(i=2;i<=ArgCnt;i++)
   {
      PtrArg = *++ArgVec;

      if(!strcmp(PtrArg,"-in"))
      {
         TmpStr = *++ArgVec;
         ArgCnt--;
         strcpy(InpNam, TmpStr);

         if(!strstr(InpNam, ".mesh"))
            strcat(InpNam, ".meshb");

         continue;
      }

      if(!strcmp(PtrArg,"-nproc"))
      {
         NmbCpu = atoi(*++ArgVec);
         ArgCnt--;
         continue;
      }
   }

   if(!strlen(InpNam))
   {
      puts("No input mesh provided");
      exit(1);
   }


   //
   // Mesh reading
   //
   printf("\nReading mesh        : ");
   timer = GetWallClock();

   // Check mesh format
   if(!(InpMsh = GmfOpenMesh(InpNam, GmfRead, &ver, &dim)))
   {
      printf("Cannot open mesh %s\n", InpNam);
      exit(1);
   }

   // Check dimension and number of tets
   NmbTet = GmfStatKwd(InpMsh, GmfTetrahedra);

   if(dim != 3 || !NmbTet)
   {
      puts("Can only handle 3D tet-meshes\n");
      exit(1);
   }

   // Get stats and allocate tables
   if(!(TetTab = malloc((size_t)(NmbTet+1) * 4 * sizeof(int))))
   {
      puts("Failed to allocate tets memory");
      exit(1);
   }

   // Read the tetrahedra
   GmfGetBlock(InpMsh, GmfTetrahedra, 1, NmbTet, 0, NULL, NULL,
               GmfIntVec, 4, &TetTab[1*4], &TetTab[ (size_t)NmbTet*4 ],
               GmfInt, &ref, &ref);

   GmfCloseMesh(InpMsh);

   printf("%g s\n", GetWallClock() - timer);
   printf("Input mesh: version = %d, tets = %d\n", ver, NmbTet);


   //
   // Set tets' neighbours and unique inner edges
   //

   NgbTab = malloc( (size_t)(NmbTet + 1) * 4 * sizeof(int));
   VoyTab = malloc( (size_t)(NmbTet + 1) * 4 * sizeof(char));

   if(!NgbTab || !VoyTab)
   {
      puts("Failed to allocate neighbours memory");
      exit(1);
   }

   timer = GetWallClock();
   ParallelNeighbours(NmbCpu, NmbTet, LplTet, TetTab, NgbTab, VoyTab);
   printf("ParallelNeighbours = %f s\n",GetWallClock() - timer);

   timer = GetWallClock();
   NmbEdg = ParallelBuildEdges(NmbCpu, NmbTet, LplTet, TetTab, &EdgTab);
   printf("ParallelBuildEdges = %f s, %d edges extracted\n",GetWallClock() - timer, NmbEdg);
}
