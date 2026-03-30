

/*----------------------------------------------------------------------------*/
/*                                                                            */
/*                               LPlib Helpers V1.2                           */
/*                                                                            */
/*----------------------------------------------------------------------------*/
/*                                                                            */
/* Description:         lplib's helper functions' headers                     */
/* Author:              Loic MARECHAL                                         */
/* Creation date:       may 16 2024                                           */
/* Last modification:   mar 30 2026                                           */
/*                                                                            */
/*----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------*/
/* Includes                                                                   */
/*----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdarg.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <libmeshb8.h>
#include "lplib4.h"
#include "lplib4_helpers.h"


/*----------------------------------------------------------------------------*/
/* Defintion of macro commands and constants                                  */
/*----------------------------------------------------------------------------*/

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define POW(a)    ((a) * (a))
#define MAXEDG    1000


/*----------------------------------------------------------------------------*/
/* Structures                                                                 */
/*----------------------------------------------------------------------------*/

typedef struct
{
   itg MinIdx, MaxIdx, NexBuc;
}HshSct;

typedef struct
{
   int      tet;
   char     voy, min, mid, max;
   size_t   nex;
}HshTriSct;


typedef struct
{
   char        *FlgTab, *VoyTab;
   itg         beg, end, NmbEdg, EdgAdr, *EleTab, (*EdgTab)[2], *NgbTab;
   int         NmbCpu;
   int64_t     HshSiz, HshPos, HshMsk, ColPos;
   HshTriSct   *tab;
   HshSct      *HshTab;
}ParSct;


/*----------------------------------------------------------------------------*/
/* Prototypes of local procedures                                             */
/*----------------------------------------------------------------------------*/

static void ParEdg1(itg, itg, int, ParSct *);
static void ParEdg2(itg, itg, int, ParSct *);
static void ParNgb1(int, int, int, ParSct *);
static void ParNgb2(int, int, int, ParSct *);


/*----------------------------------------------------------------------------*/
/* Build edges in parallel: head procedure                                    */
/*----------------------------------------------------------------------------*/

itg ParallelBuildEdges( int NmbCpu, itg NmbEle, int EleTyp,
                        itg *EleTab, itg **UsrEdg )
{
   itg      i, HshSiz, IncSiz, NmbEdg = 0, (*EdgTab)[2];
   int64_t  LibIdx;
   int      TetTyp;
   ParSct   par[ MaxPth ];

   // As for now only tets are supported
   if(EleTyp != LplTet)
      return(0);

   // Setup LPlib and datatypes
   LibIdx = InitParallel(NmbCpu);
   TetTyp = NewType(LibIdx, NmbEle);

   // Setup parallel parameters
   IncSiz = (NmbEle / NmbCpu) / NmbCpu;
   HshSiz =  IncSiz * NmbCpu;

   for(i=0;i<NmbCpu;i++)
   {
      par[i].beg = i * IncSiz;
      par[i].end = (i + 1) * IncSiz;
      par[i].HshSiz = HshSiz;
      par[i].ColPos = HshSiz;
      par[i].EleTab = EleTab;
      par[i].NmbCpu = NmbCpu;
      par[i].EdgAdr = 0;
   }

   // Each thread builds a local edge table
   LaunchParallel(LibIdx, TetTyp, 0, (void *)ParEdg1, (void *)par);

   // Count the number of unique edges in each thread
   LaunchParallel(LibIdx, TetTyp, 0, (void *)ParEdg2, (void *)par);

   // Allocate the global edge table and give a slice of it to each threads

   for(i=0;i<NmbCpu;i++)
   {
      par[i].EdgAdr = NmbEdg + 1;
      NmbEdg += par[i].NmbEdg;
   }

   EdgTab = malloc((size_t)(NmbEdg+1) * 2 * sizeof(itg));
   assert(EdgTab);

   for(i=0;i<NmbCpu;i++)
      par[i].EdgTab = EdgTab;

   // Now each threads counts and stores the unique edges
   LaunchParallel(LibIdx, TetTyp, 0, (void *)ParEdg2, (void *)par);

   // Free the local hash tables
   for(i=0;i<NmbCpu;i++)
      free(par[i].HshTab);

   StopParallel(LibIdx);

   *UsrEdg = (itg *)EdgTab;

   return(NmbEdg);
}


/*----------------------------------------------------------------------------*/
/* Build thread subdomain's edges                                             */
/*----------------------------------------------------------------------------*/

static void ParEdg1(itg BegIdx, itg EndIdx, int PthIdx, ParSct *par)
{
   itg i, j, key, idx0, idx1, MinIdx, MaxIdx;
   itg siz = par[ PthIdx ].HshSiz, col = par[ PthIdx ].ColPos;
   HshSct *hsh;
   itg *tet, *EleTab = par[ PthIdx ].EleTab;
   const int tvpe[6][2] = { {0,1}, {1,2}, {2,0}, {3,0}, {3,1}, {3,2} };

   // allocate a thread local hash table
   hsh = par[ PthIdx ].HshTab = malloc(7LL * (size_t)siz * sizeof(HshSct));
   assert(hsh);

   // Clear the hash table's direct entries,
   // there is no need to clear the collision entries
   memset(hsh, 0, (size_t)siz * sizeof(HshSct));

   // Loop over each tet and each tet's edges
   for(i=BegIdx; i<=EndIdx; i++)
   {
      tet = &EleTab[ i * 4 ];

      for(j=0;j<6;j++)
      {
         // Compute the hashing key from the edge's vertices indices
         idx0 = tet[ tvpe[j][0] ];
         idx1 = tet[ tvpe[j][1] ];

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

         key = (3 * (size_t)MinIdx + 5 * (size_t)MaxIdx) % siz;

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
               hsh[ key ].NexBuc = 0;
               break;
            }
         }while(1);
      }
   }
}


/*----------------------------------------------------------------------------*/
/* Setup the missing links between tets that cross subdomains                 */
/*----------------------------------------------------------------------------*/

static void ParEdg2(itg BegIdx, itg EndIdx, int PthIdx, ParSct *par)
{
   itg i, key, PthNmbEdg = 0, edg[ MAXEDG ][2];
   itg NmbCpu = par[ PthIdx ].NmbCpu;
   int NmbEdg, flg, j, k;
   HshSct *buc;
   itg (*EdgTab)[2] = par[ PthIdx ].EdgTab;

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

                  if(NmbEdg >= MAXEDG)
                  {
                     puts("Too many local edges, increase MAXEDG value.");
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
            EdgTab[ par[ PthIdx ].EdgAdr + PthNmbEdg + j ][0] = edg[j][0];
            EdgTab[ par[ PthIdx ].EdgAdr + PthNmbEdg + j ][1] = edg[j][1];
         }

      PthNmbEdg += NmbEdg;
   }

   par[ PthIdx ].NmbEdg = PthNmbEdg;
}


/*----------------------------------------------------------------------------*/
/* Parallel neighbours between tets                                           */
/*----------------------------------------------------------------------------*/

int ParallelNeighbours( int NmbCpu, itg NmbEle, int EleTyp,
                        itg *EleTab, itg *UsrNgb, char *UsrVoy )
{
   char     *FlgTab;
   int      i, j, k, dec, NgbIdx, NmbTyp, TetTyp;
   itg      *NgbTab;
   int      tvpf[4][3] = { {1,2,3}, {2,0,3}, {3,0,1}, {0,2,1} };
   size_t   HshSiz, MshSiz;
   int64_t  LibIdx;
   ParSct   par[ MaxPth ];

   // As for now only tets are supported
   if(EleTyp != LplTet)
      return(0);

   // Setup LPlib and datatypes
   if(!(LibIdx = InitParallel(NmbCpu)))
      return(0);

   if(!(TetTyp = NewType(LibIdx, NmbEle)))
      return(0);

   // Allocate a hash table and overflow buffer
   dec =    ceil(log(1. + 2. * NmbEle / NmbCpu) / log(2));
   HshSiz = 1LL << dec;
   MshSiz = NmbEle / NmbCpu;

   if(!(FlgTab = calloc((NmbEle+1), sizeof(char))))
      return(0);

   ParallelMemClear(LibIdx, UsrNgb, (size_t)(NmbEle + 1) * 4 * sizeof(itg));

   if(UsrVoy)
      ParallelMemClear(LibIdx, UsrVoy, (size_t)(NmbEle + 1) * 4 * sizeof(char));

   // Setup parallel parameters
   for(i=0;i<NmbCpu;i++)
   {
      par[i].beg = i * MshSiz + 1;
      par[i].end = (i+1) * MshSiz;
      par[i].HshSiz = HshSiz;
      par[i].ColPos = HshSiz;
      par[i].HshMsk = HshSiz - 1;
      par[i].EleTab = EleTab;
      par[i].FlgTab = FlgTab;
      par[i].NmbCpu = NmbCpu;
      par[i].NgbTab = UsrNgb;
      par[i].VoyTab = UsrVoy;
   }

   par[ NmbCpu-1 ].end = NmbEle;

   /* Launch parallel loops: the first one build local neighbours
      among each subdomains and the second one build neighbourhood
      information between cross block elements */

   LaunchParallel(LibIdx, TetTyp, 0, (void *)ParNgb1, (void *)par);

   if(NmbCpu > 1)
      LaunchParallel(LibIdx, TetTyp, 0, (void *)ParNgb2, (void *)par);

   free(FlgTab);

   for(i=0;i<NmbCpu;i++)
      free(par[i].tab);

   return(1);
}


/*----------------------------------------------------------------------------*/
/* Set links between tets from this local subdomain                           */
/*----------------------------------------------------------------------------*/

static void ParNgb1(int BegIdx, int EndIdx, int c, ParSct *par)
{
   char        *FlgTab = par[c].FlgTab, *VoyTab = par[c].VoyTab;
   int         i, j, k;
   uint32_t    min, mid, max;
   uint64_t    key;
   itg         *EleTab = par[c].EleTab, *NgbTab = par[c].NgbTab, *tet, *ngb;
   HshTriSct   *tab = par[c].tab = calloc(5LL * par[c].HshSiz, sizeof(HshTriSct));

   // Allocate a local hash table and loop over the local elements
   for(i=par[c].beg; i<=par[c].end; i++)
   {
      tet = &EleTab[ (size_t)i * 4LL ];

      for(j=0;j<4;j++)
      {
         // Compute the hashing key from the face's vertices indices
         min = max = (j+1)%4;

         for(k=0;k<4;k++)
            if(k != j)
            {
               if(tet[k] < tet[ min ])
                  min = k;
               else if(tet[k] > tet[ max ])
                  max = k;
            }

         mid = 6 - min - max - j;
         key = (31LL * (size_t)tet[ min ]
               + 7LL * (size_t)tet[ mid ]
               + 3LL * (size_t)tet[ max ]) & par[c].HshMsk;

         // If the bucket is empty, store the face
         if(!tab[ key ].tet)
         {
            tab[ key ].tet = i;
            tab[ key ].voy = j;
            tab[ key ].min = min;
            tab[ key ].mid = mid;
            tab[ key ].max = max;
            continue;
         }

         // Otherwise, search through the linked list
         do
         {
            ngb = &EleTab[ (size_t)tab[ key ].tet * 4LL ];

            // If the same face is found in the hash table, 
            // setup a link between both tetrahedra
            if( (ngb[ (int)tab[ key ].min ] == tet[ min ])
            &&  (ngb[ (int)tab[ key ].mid ] == tet[ mid ])
            &&  (ngb[ (int)tab[ key ].max ] == tet[ max ]) )
            {
               NgbTab[ (size_t)i * 4LL + (size_t)j ] = tab[ key ].tet;
               FlgTab[i]++;
               NgbTab[ (size_t)tab[ key ].tet * 4LL + (size_t)tab[ key ].voy ] = i;
               FlgTab[ tab[ key ].tet ]++;

               if(VoyTab)
               {
                  VoyTab[ (size_t)i * 4LL + (size_t)j ] = tab[ key ].voy;
                  VoyTab[ (size_t)tab[ key ].tet * 4LL + (size_t)tab[ key ].voy ] = j;
               }

               break;
            }

            // If not, allocate a new bucket from the overflow table
            // and link it to the main entry
            if(tab[ key ].nex)
               key = tab[ key ].nex;
            else
            {
               tab[ key ].nex = par[c].ColPos;
               key = par[c].ColPos++;
               tab[ key ].tet = i;
               tab[ key ].voy = j;
               tab[ key ].min = min;
               tab[ key ].mid = mid;
               tab[ key ].max = max;
               break;
            }
         }while(1);
      }
   }
}


/*----------------------------------------------------------------------------*/
/* Setup the missing links between tets that cross subdomains                 */
/*----------------------------------------------------------------------------*/

static void ParNgb2(int BegIdx, int EndIdx, int c, ParSct *par)
{
   char        *VoyTab = par[c].VoyTab;
   int         n, i, j, k, key, BasKey, flg;
   uint32_t    min, mid, max;
   itg         *tet, *ngb, *EleTab = par[c].EleTab, *NgbTab = par[c].NgbTab;
   HshTriSct   *tab;

   for(i=par[c].beg; i<=par[c].end; i++)
   {
      // If a tetrahedron has already 4 links,
      // there is no need to find a missing ones
      if(par[c].FlgTab[i] == 4)
         continue;

      tet = &EleTab[ (size_t)i * 4 ];

      for(j=0;j<4;j++)
      {
         // If there is no neighbour through this face,
         // try to find on among other subdomains local hash tables
         if(NgbTab[ (size_t)i * 4 + j ])
            continue;

         min = max = (j+1)%4;

         for(k=0;k<4;k++)
            if(k != j)
            {
               if(tet   [k] < tet[ min ])
                  min = k;
               else if(tet[k] > tet[ max ])
                  max = k;
            }

         mid = 6 - min - max - j;
         flg = 0;
         BasKey = (31 * (size_t)tet[ min ] + 7 * (size_t)tet[ mid ] + 3 * (size_t)tet[ max ]) & par[c].HshMsk;

         for(n=0; n<par[c].NmbCpu; n++)
         {
            if(n == c)
               continue;

            tab = par[n].tab;
            key = BasKey;

            do
            {
               ngb = &EleTab[ (size_t)tab[ key ].tet * 4 ];

               if( (ngb[ (int)tab[ key ].min ] == tet[ min ])
               &&  (ngb[ (int)tab[ key ].mid ] == tet[ mid ])
               &&  (ngb[ (int)tab[ key ].max ] == tet[ max ]) )
               {
                  NgbTab[ (size_t)i * 4 + j ] = tab[ key ].tet;

                  if(VoyTab)
                     VoyTab[ (size_t)i * 4 + j ] = tab[ key ].voy;

                  flg = 1;
                  break;
               }

               if(tab[ key ].nex)
                  key = tab[ key ].nex;
               else
                  break;
            }while(1);

            if(flg)
               break;
         }
      }
   }
}
