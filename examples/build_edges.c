

/*----------------------------------------------------------------------------*/
/*                                                                            */
/*                 PARALLEL EDGE LIST BUILDING USING LPlib                    */
/*                                                                            */
/*----------------------------------------------------------------------------*/
/*                                                                            */
/*   Description:       build the list of unique edges                        */
/*                      from a volumic tetrahedral mesh                       */
/*   Author:            Loic MARECHAL                                         */
/*   Creation date:     feb 13 2015                                           */
/*   Last modification: mar 06 2024                                           */
/*                                                                            */
/*----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------*/
/* Includes                                                                   */
/*----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <string.h>
#include <inttypes.h>
#include <assert.h>
#include "libmeshb7.h"
#include "lplib3.h"


/*----------------------------------------------------------------------------*/
/* Defines                                                                    */
/*----------------------------------------------------------------------------*/

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MaxEdg 1000

#ifdef INT64
#define GmfItg GmfLong
#define GmfItgVec GmfLongVec
#else
#define GmfItg GmfInt
#define GmfItgVec GmfIntVec
#endif


/*----------------------------------------------------------------------------*/
/* Structures                                                                 */
/*----------------------------------------------------------------------------*/

typedef struct
{
   double   crd[3];
   int      ref;
}VerSct;

typedef struct
{
   itg idx[2];
}EdgSct;

typedef struct
{
   itg idx[3];
   int ref;
}TriSct;

typedef struct
{
   itg idx[4];
}TetSct;

typedef struct
{
   itg MinIdx, MaxIdx, NexBuc;
}HshSct;

typedef struct
{
   itg      NmbVer, NmbEdg, NmbTri, NmbTet;
   int      MshVer, LibIdx;
   VerSct   *ver;
   EdgSct   *edg;
   TriSct   *tri;
   TetSct   *tet;
}MshSct;

typedef struct
{
   itg      beg, end, HshSiz, ColPos, NmbEdg, EdgAdr;
   int      NmbCpu;
   HshSct   *HshTab;
   MshSct   *msh;
}ParSct;


/*----------------------------------------------------------------------------*/
/* Prototypes of local procedures                                             */
/*----------------------------------------------------------------------------*/

void SetEdgSer (MshSct *);
void SetEdgPar (MshSct *, int);
void ParEdg1   (itg, itg, int, ParSct *);
void ParEdg2   (itg, itg, int, ParSct *);
void ScaMsh    (char *, MshSct *);
void RecMsh    (char *, MshSct *);
void GetTim    (double *);
void PrtVerDeg (MshSct *);


/*----------------------------------------------------------------------------*/
/* Global tables                                                              */
/*----------------------------------------------------------------------------*/

const int tvpe[6][2] = { {0,1}, {1,2}, {2,0}, {3,0}, {3,1}, {3,2} };


/*----------------------------------------------------------------------------*/
/* Read, build the edge list and write the mesh                               */
/*----------------------------------------------------------------------------*/

int main(int ArgCnt, char **ArgVec)
{
   char *PtrArg, *TmpStr, InpNam[1000], OutNam[1000];
   int i, NmbCpu = 0, VecMod = 0;
   MshSct msh;

   // Command line parsing
   memset(&msh, 0, sizeof(MshSct));

   if(ArgCnt == 1)
   {
      puts("\nbuild_edges v1.01 feb 28 2024   Loic MARECHAL / INRIA");
      puts(" Usage       : build_edges -in volume_mesh -out edge_mesh");
      puts(" -in  name   : name of the input tetrahedral-only mesh");
      puts(" -out name   : name of the output mesh that will contain tets, edges and vertices");
      puts(" -vector     : print statistics on vertex connectivity and vector padding");
      puts(" -serial     : use the serial optimized version (different from -nproc 1)");
      puts(" -nproc n    : n is the number of threads (default = all available threads)\n");
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

      if(!strcmp(PtrArg,"-out"))
      {
         TmpStr = *++ArgVec;
         ArgCnt--;
         strcpy(OutNam, TmpStr);

         if(!strstr(OutNam, ".mesh"))
            strcat(OutNam, ".meshb");

         continue;
      }

      if(!strcmp(PtrArg,"-vector"))
      {
         VecMod = 1;
         continue;
      }

      if(!strcmp(PtrArg,"-serial"))
      {
         NmbCpu = -1;
         continue;
      }

      if(!strcmp(PtrArg,"-nproc"))
      {
         NmbCpu = atoi(*++ArgVec);
         NmbCpu = MAX(NmbCpu, 1);
         NmbCpu = MIN(NmbCpu, MaxPth);
         ArgCnt--;
         continue;
      }
   }

   if(!strlen(InpNam))
   {
      puts("No input mesh provided");
      exit(1);
   }

   if(!strlen(OutNam))
   {
      puts("No output name provided");
      exit(1);
   }

   // Mesh reading
   ScaMsh(InpNam, &msh);

   // Launch the parallel neeighbour procedure
   if(NmbCpu == -1)
      SetEdgSer(&msh);
   else
      SetEdgPar(&msh, NmbCpu);

   // Print statistics on vertex connectivity
   if(VecMod)
      PrtVerDeg(&msh);

   // Mesh writing
   RecMsh(OutNam, &msh);

   free(msh.ver);
   free(msh.edg);
   free(msh.tri);
   free(msh.tet);
}


/*----------------------------------------------------------------------------*/
/* Build the list of unique edges sequentialy                                 */
/*----------------------------------------------------------------------------*/

void SetEdgSer(MshSct *msh)
{
   itg i, j, idx0, idx1, key, MinIdx, MaxIdx, siz = msh->NmbTet, col = siz;
   double timer = 0.;
   TetSct *tet;
   HshSct *hsh, *buc;

   // Allocate the hash table to store all edges
   hsh = malloc(6 * siz * sizeof(HshSct));
   assert(hsh);

   // Start timer
   printf("Build edges sequentialy     : ");
   GetTim(&timer);

   // Clear the hash table's direct entries,
   // there is no need to clear the collision entries
   memset(hsh, 0, siz * sizeof(HshSct));

   // Loop over each tet and each tet's edges
   for(i=1;i<=msh->NmbTet;i++)
   {
      tet = &msh->tet[i];

      for(j=0;j<6;j++)
      {
         // Compute the hashing key from the edge's vertex indices
         idx0 = tet->idx[ tvpe[j][0] ];
         idx1 = tet->idx[ tvpe[j][1] ];

         if(idx0 < idx1)
         {
            MinIdx = idx0;
            MaxIdx = idx1;
         }
         else
         {
            MinIdx = idx1;
            MaxIdx = idx0;
         }

         key = (3 * MinIdx + 5 * MaxIdx) % siz;

         // If the bucket is empty, store the edge
         if(!hsh[ key ].MinIdx)
         {
            hsh[ key ].MinIdx = MinIdx;
            hsh[ key ].MaxIdx = MaxIdx;
            msh->NmbEdg++;
            continue;
         }

         // Otherwise, search through the linked list
         do
         {
            // If the same edge is found in the hash table, do nothing
            if( (hsh[ key ].MinIdx == MinIdx) && (hsh[ key ].MaxIdx == MaxIdx) )
               break;

            // If not, allocate a new bucket from the overflow table
            // and link it to the main entry
            if(hsh[ key ].NexBuc)
               key = hsh[ key ].NexBuc;
            else
            {
               hsh[ key ].NexBuc = col;
               key = col++;
               hsh[ key ].MinIdx = MinIdx;
               hsh[ key ].MaxIdx = MaxIdx;
               msh->NmbEdg++;
               break;
            }
         }while(1);
      }
   }

   // Now allocate and edge table
   msh->edg = malloc((msh->NmbEdg+1) * sizeof(EdgSct));
   assert(msh->edg);
   msh->NmbEdg = 0;

   // Loop over the hash table's direct entry and store the edges
   for(i=0;i<msh->NmbTet;i++)
   {
      key = i;

      // Follow the link if this bucket has some collisions
      do
      {
         buc = &hsh[ key ];

         if(buc->MinIdx)
         {
            msh->NmbEdg++;
            msh->edg[ msh->NmbEdg ].idx[0] = buc->MinIdx;
            msh->edg[ msh->NmbEdg ].idx[1] = buc->MaxIdx;
         }
      }while((key = buc->NexBuc));
   }

   // Free the hash table and stop the timer
   free(hsh);
   GetTim(&timer);
   printf("%g s\n", timer);
   printf("Unique edges found          : %lld\n", (int64_t)msh->NmbEdg);
}


/*----------------------------------------------------------------------------*/
/* Build edges in parallel: head procedure                                    */
/*----------------------------------------------------------------------------*/

void SetEdgPar(MshSct *msh, int NmbCpu)
{
   itg i, HshSiz, IncSiz, adr = 0;
   int64_t LibIdx;
   int TetTyp;
   float sta[2];
   double timer = 0.;
   HshSct *HshTab;
   ParSct par[ MaxPth ];

   // Setup LPlib and datatypes
   if(!NmbCpu)
      NmbCpu = GetNumberOfCores();

   printf("Build edges with %3d threads: ", NmbCpu);

   GetTim(&timer);
   LibIdx = InitParallel(NmbCpu);
   TetTyp = NewType(LibIdx, msh->NmbTet);

   // Setup parallel parameters
   IncSiz = (msh->NmbTet / NmbCpu) / NmbCpu;
   HshSiz =  IncSiz * NmbCpu;

   for(i=0;i<NmbCpu;i++)
   {
      par[i].beg = i * IncSiz;
      par[i].end = (i + 1) * IncSiz;
      par[i].HshSiz = HshSiz;
      par[i].ColPos = HshSiz;
      par[i].msh = msh;
      par[i].NmbCpu = NmbCpu;
      par[i].EdgAdr = 0;
   }

   // Each thread builds a local edge table
   LaunchParallel(LibIdx, TetTyp, 0, (void *)ParEdg1, (void *)par);

   // Count the number of unique edges in each thread
   LaunchParallel(LibIdx, TetTyp, 0, (void *)ParEdg2, (void *)par);

   // Allocate the global edge table and give a slice of it to each threads
   msh->NmbEdg = 0;

   for(i=0;i<NmbCpu;i++)
   {
      par[i].EdgAdr = msh->NmbEdg + 1;
      msh->NmbEdg += par[i].NmbEdg;
   }

   msh->edg = malloc((msh->NmbEdg+1) * sizeof(EdgSct));
   assert(msh->edg);

   // Now each threads counts and stores the unique edges
   LaunchParallel(LibIdx, TetTyp, 0, (void *)ParEdg2, (void *)par);

   // Free the local hash tables
   for(i=0;i<NmbCpu;i++)
      free(par[i].HshTab);

   StopParallel(LibIdx);

   GetTim(&timer);
   printf("%g s\n", timer);

   printf("Unique edges found          : %lld\n", (int64_t)msh->NmbEdg);
}


/*----------------------------------------------------------------------------*/
/* Build thread subdomain's edges                                             */
/*----------------------------------------------------------------------------*/

void ParEdg1(itg BegIdx, itg EndIdx, int PthIdx, ParSct *par)
{
   itg i, j, key, idx0, idx1, MinIdx, MaxIdx;
   itg siz = par[ PthIdx ].HshSiz, col = par[ PthIdx ].ColPos;
   HshSct *hsh;
   MshSct *msh = par[ PthIdx ].msh;
   TetSct *tet;

   // allocate a thread local hash table
   hsh = par[ PthIdx ].HshTab = malloc(6 * siz * sizeof(HshSct));
   assert(hsh);

   // Clear the hash table's direct entries,
   // there is no need to clear the collision entries
   memset(hsh, 0, siz * sizeof(HshSct));

   // Loop over each tet and each tet's edges
   for(i=BegIdx; i<=EndIdx; i++)
   {
      tet = &msh->tet[i];

      for(j=0;j<6;j++)
      {
         // Compute the hashing key from the edge's vertices indices
         idx0 = tet->idx[ tvpe[j][0] ];
         idx1 = tet->idx[ tvpe[j][1] ];

         if(idx0 < idx1)
         {
            MinIdx = idx0;
            MaxIdx = idx1;
         }
         else
         {
            MinIdx = idx1;
            MaxIdx = idx0;
         }

         key = (3 * MinIdx + 5 * MaxIdx) % siz;

         // If the bucket is empty, store the edge
         if(!hsh[ key ].MinIdx)
         {
            hsh[ key ].MinIdx = MinIdx;
            hsh[ key ].MaxIdx = MaxIdx;
            continue;
         }

         // Otherwise, search through the linked list
         do
         {
            // If the same edge is found in the hash table, do nothing
            if( (hsh[ key ].MinIdx == MinIdx) && (hsh[ key ].MaxIdx == MaxIdx) )
               break;

            // If not, allocate a new bucket from the overflow table
            // and link it to the main entry */
            if(hsh[ key ].NexBuc)
               key = hsh[ key ].NexBuc;
            else
            {
               hsh[ key ].NexBuc = col;
               key = col++;
               hsh[ key ].MinIdx = MinIdx;
               hsh[ key ].MaxIdx = MaxIdx;
               break;
            }
         }while(1);
      }
   }
}


/*----------------------------------------------------------------------------*/
/* Setup the missing links between tets that cross subdomains                 */
/*----------------------------------------------------------------------------*/

void ParEdg2(itg BegIdx, itg EndIdx, int PthIdx, ParSct *par)
{
   itg i, key, PthNmbEdg = 0, edg[ MaxEdg ][2];
   itg siz = par[ PthIdx ].HshSiz, NmbCpu = par[ PthIdx ].NmbCpu;
   int NmbEdg, flg, j, k;
   HshSct *hsh = par[ PthIdx ].HshTab, *buc;
   MshSct *msh = par[ PthIdx ].msh;
   TetSct *tet;

   // Loop over the hash table direct entries following an interleaved stencil
   for(i=par[ PthIdx ].beg; i<par[ PthIdx ].end; i++)
   {
      NmbEdg = 0;

      // Loop over every entries with the same key
      // among all threads' local hash tables
      for(j=0;j<NmbCpu;j++)
      {
         key = i;

         // In case of collision, follow the links
         do
         {
            buc = &par[j].HshTab[ key ];

            if(buc->MinIdx)
            {
               // Since edges from different local hash tables may be the same, 
               // they compared again to avoid duplicates
               flg = 0;

               for(k=0;k<NmbEdg;k++)
                  if( (buc->MinIdx == edg[k][0]) && (buc->MaxIdx == edg[k][1]) )
                  {
                     flg= 1;
                     break;
                  }

               // If this edge does not belong to the list, add it to the end
               if(!flg)
               {
                  edg[ NmbEdg ][0] = buc->MinIdx;
                  edg[ NmbEdg ][1] = buc->MaxIdx;
                  NmbEdg++;

                  if(NmbEdg >= MaxEdg)
                  {
                     puts("Too many local edges, increase MaxEdg value.");
                     exit(1);
                  }
               }
            }
         }while((key = buc->NexBuc));
      }

      // On the second run, add the list of local edges to the global ones
      if(par[ PthIdx ].EdgAdr)
         for(j=0;j<NmbEdg;j++)
         {
            msh->edg[ par[ PthIdx ].EdgAdr + PthNmbEdg + j ].idx[0] = edg[j][0];
            msh->edg[ par[ PthIdx ].EdgAdr + PthNmbEdg + j ].idx[1] = edg[j][1];
         }

      PthNmbEdg += NmbEdg;
   }

   par[ PthIdx ].NmbEdg = PthNmbEdg;
}


/*----------------------------------------------------------------------------*/
/* Wall clock timer                                                           */
/*----------------------------------------------------------------------------*/

void GetTim(double *timer)
{
   *timer = GetWallClock() - *timer;
}


/*----------------------------------------------------------------------------*/
/* Read mesh                                                                  */
/*----------------------------------------------------------------------------*/

void ScaMsh(char *InpNam, MshSct *msh)
{
   int dim, ref;
   int64_t InpMsh;
   float flt[3];
   double timer = 0.;

   printf("\nRead mesh                   : ");
   GetTim(&timer);

   // Check mesh format
   if(!(InpMsh = GmfOpenMesh(InpNam, GmfRead, &msh->MshVer, &dim)))
   {
      printf("Cannot open mesh %s\n", InpNam);
      exit(1);
   }

   if(dim != 3)
   {
      puts("Can only handle 3D meshes");
      exit(1);
   }

   // Get stats and allocate tables
   if(!(msh->NmbVer = GmfStatKwd(InpMsh, GmfVertices)))
   {
      puts("No vertices found");
      exit(1);
   }

   if(!(msh->NmbTri = GmfStatKwd(InpMsh, GmfTriangles)))
   {
      puts("No triangles found");
      exit(1);
   }

   if(!(msh->NmbTet = GmfStatKwd(InpMsh, GmfTetrahedra)))
   {
      puts("No tetrahedra found");
      exit(1);
   }

   msh->ver = malloc((msh->NmbVer+1) * sizeof(VerSct));
   assert(msh->ver);
   msh->tri = malloc((msh->NmbTri+1) * sizeof(TriSct));
   assert(msh->tri);
   msh->tet = malloc((msh->NmbTet+1) * sizeof(TetSct));
   assert(msh->tet);

   // Read the vertices
   GmfGetBlock(InpMsh, GmfVertices, 1, msh->NmbVer, 0, NULL, NULL,
               GmfDoubleVec, 3, msh->ver[1].crd, msh->ver[ msh->NmbVer ].crd,
               GmfInt, &msh->ver[1].ref, &msh->ver[ msh->NmbVer ].ref);
   
   // Read the triangles
   GmfGetBlock(InpMsh, GmfTriangles, 1, msh->NmbTri, 0, NULL, NULL,
               GmfItgVec, 3, msh->tri[1].idx, msh->tri[ msh->NmbTri ].idx,
               GmfInt, &msh->tri[1].ref, &msh->tri[ msh->NmbTri ].ref);
   
   // Read the tetrahedra
   GmfGetBlock(InpMsh, GmfTetrahedra, 1, msh->NmbTet, 0, NULL, NULL,
               GmfItgVec, 4, msh->tet[1].idx, msh->tet[ msh->NmbTet ].idx,
               GmfInt, &ref, &ref);

   GmfCloseMesh(InpMsh);

   GetTim(&timer);
   printf("%g s\n", timer);
   printf("Input mesh                  : version = %d, %lld vertices, %lld triangles, %lld tets\n",
         msh->MshVer, (int64_t)msh->NmbVer, (int64_t)msh->NmbTri, (int64_t)msh->NmbTet);
}


/*----------------------------------------------------------------------------*/
/* Write mesh                                                                 */
/*----------------------------------------------------------------------------*/

void RecMsh(char *OutNam, MshSct *msh)
{
   int ref=0;
   int64_t OutMsh;
   double timer = 0.;

   printf("Write mesh                  : ");
   GetTim(&timer);

   // Create the output mesh
   if(!msh->NmbVer || !msh->NmbEdg || !msh->NmbTet \
   || !(OutMsh = GmfOpenMesh(OutNam, GmfWrite, msh->MshVer, 3)) )
   {
      printf("Cannot create mesh %s\n", OutNam);
      exit(1);
   }

   // Save the vertices from the input mesh
   GmfSetKwd(OutMsh, GmfVertices, msh->NmbVer);
   GmfSetBlock(OutMsh, GmfVertices, 1, msh->NmbVer, 0, NULL, NULL,
               GmfDoubleVec, 3, msh->ver[1].crd, msh->ver[ msh->NmbVer ].crd,
               GmfInt, &msh->ver[1].ref, &msh->ver[ msh->NmbVer ].ref);

   // Save the extracted edges
   GmfSetKwd(OutMsh, GmfEdges, msh->NmbEdg);
   GmfSetBlock(OutMsh, GmfEdges,  1, msh->NmbEdg, 0, NULL, NULL,
               GmfItgVec, 2, msh->edg[1].idx, &msh->edg[ msh->NmbEdg ].idx,
               GmfInt, &ref, &ref);

   // Save the triangles from the input mesh
   GmfSetKwd(OutMsh, GmfTriangles, msh->NmbTri);
   GmfSetBlock(OutMsh, GmfTriangles, 1, msh->NmbTri, 0, NULL, NULL,
               GmfItgVec, 3, msh->tri[1].idx, msh->tri[ msh->NmbTri ].idx,
               GmfInt, &msh->tri[1].ref, &msh->tri[ msh->NmbTri ].ref);

   // Save the tetrahedra from the input mesh
   GmfSetKwd(OutMsh, GmfTetrahedra, msh->NmbTet);
   GmfSetBlock(OutMsh, GmfTetrahedra, 1, msh->NmbTet, 0, NULL, NULL,
               GmfItgVec, 4, msh->tet[1].idx, msh->tet[ msh->NmbTet ].idx,
               GmfInt, &ref, &ref);

   GmfCloseMesh(OutMsh);

   GetTim(&timer);
   printf("%g s\n\n", timer);
}


/*----------------------------------------------------------------------------*/
/* Print statistics about vertex connectivity                                 */
/*----------------------------------------------------------------------------*/

void PrtVerDeg(MshSct *msh)
{
   int i, j, *DegTab = calloc(msh->NmbVer+1, sizeof(int));
   int64_t DegTot=0, VecTot=0, deg16=0, deg32=0, deg64=0, deg128=0, deg256=0, ovf=0;

   assert(DegTab);

   for(i=1;i<=msh->NmbEdg;i++)
      for(j=0;j<2;j++)
         DegTab[ msh->edg[i].idx[j] ]++;

   for(i=1;i<=msh->NmbVer;i++)
   {
      DegTot += DegTab[i];

      if(DegTab[i] <= 16)
         deg16++;
      else if(DegTab[i] <= 32)
         deg32++;
      else if(DegTab[i] <= 64)
         deg64++;
      else if(DegTab[i] <= 128)
         deg128++;
      else if(DegTab[i] <= 256)
         deg256++;
      else
         ovf++;
   }

   free(DegTab);

   VecTot = 16*deg16 + 32*deg32 + 64*deg64 + 128*deg128 + 256*deg256;

   puts("");
   puts  ("vector |  %age  | number");
   puts  ("----------------------------");
   printf("   16  | %6.2f | %10lld\n", (float)(100 * deg16 ) / msh->NmbVer, deg16);
   printf("   32  | %6.2f | %10lld\n", (float)(100 * deg32 ) / msh->NmbVer, deg32);
   printf("   64  | %6.2f | %10lld\n", (float)(100 * deg64 ) / msh->NmbVer, deg64);
   printf("  128  | %6.2f | %10lld\n", (float)(100 * deg128) / msh->NmbVer, deg128);
   printf("  256  | %6.2f | %10lld\n", (float)(100 * deg256) / msh->NmbVer, deg256);
   printf("  OUT  | %6.2f | %10lld\n", (float)(100 * ovf   ) / msh->NmbVer, ovf);

   puts("");
   printf("vector filling : %3.2f%%\n", (float)(100*DegTot)/VecTot);
   printf("real non-zero  : %10lld\n", DegTot);
   printf("vector non-zero: %10lld\n", VecTot);
   puts("");
}
