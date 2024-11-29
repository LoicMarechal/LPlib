

/*----------------------------------------------------------------------------*/
/*                                                                            */
/*                               LPlib Helpers V0.1                           */
/*                                                                            */
/*----------------------------------------------------------------------------*/
/*                                                                            */
/* Description:         lplib's helper functions' headers                     */
/* Author:              Loic MARECHAL                                         */
/* Creation date:       may 16 2024                                           */
/* Last modification:   oct 31 2024                                           */
/*                                                                            */
/*----------------------------------------------------------------------------*/


#ifndef LPLIB3_HELPERS_H
#define LPLIB3_HELPERS_H


/*----------------------------------------------------------------------------*/
/* Includes                                                                   */
/*----------------------------------------------------------------------------*/

#include <assert.h>
#include <string.h>
#include <stdio.h>
#include "lplib3.h"
#include "lplib3_helpers.h"


/*----------------------------------------------------------------------------*/
/* Defines                                                                    */
/*----------------------------------------------------------------------------*/

#ifdef INT64
#define itg int64_t
#else
#define itg int32_t
#endif

#ifdef REAL32
#define fpn float
#else
#define fpn double
#endif


/*----------------------------------------------------------------------------*/
/* Defintion of macro commands                                                */
/*----------------------------------------------------------------------------*/

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define POW(a) ((a) * (a))
#define MAXEDG 1000


/*----------------------------------------------------------------------------*/
/* Structures                                                                 */
/*----------------------------------------------------------------------------*/

typedef struct
{
   itg MinIdx, MaxIdx, NexBuc;
}HshSct;

typedef struct
{
   itg      beg, end, HshSiz, ColPos, NmbEdg, EdgAdr, *EleTab, (*EdgTab)[2];
   int      NmbCpu;
   HshSct   *HshTab;
}ParSct;


/*----------------------------------------------------------------------------*/
/* Prototypes of local procedures                                             */
/*----------------------------------------------------------------------------*/

void ParEdg1(itg, itg, int, ParSct *);
void ParEdg2(itg, itg, int, ParSct *);


/*----------------------------------------------------------------------------*/
/* Global tables                                                              */
/*----------------------------------------------------------------------------*/

const int tvpe[6][2] = { {0,1}, {1,2}, {2,0}, {3,0}, {3,1}, {3,2} };


/*----------------------------------------------------------------------------*/
/* Build edges in parallel: head procedure                                    */
/*----------------------------------------------------------------------------*/

itg ParallelBuildEdges(itg NmbEle, int EleTyp, itg *EleTab, itg **UsrEdg)
{
   itg      i, HshSiz, IncSiz, adr = 0, NmbEdg = 0, (*EdgTab)[2];
   int64_t  LibIdx;
   int      TetTyp, NmbCpu;
   HshSct   *HshTab;
   ParSct   par[ MaxPth ];

   // As for now only tets are supported
   //if(EleTyp != LplTet)
     //return(0);

   // Setup LPlib and datatypes
   NmbCpu = GetNumberOfCores();
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

   EdgTab = malloc((NmbEdg+1) * 2 * sizeof(itg));
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

void ParEdg1(itg BegIdx, itg EndIdx, int PthIdx, ParSct *par)
{
   itg i, j, key, idx0, idx1, MinIdx, MaxIdx;
   itg siz = par[ PthIdx ].HshSiz, col = par[ PthIdx ].ColPos;
   HshSct *hsh;
   itg *tet, *EleTab = par[ PthIdx ].EleTab;

   // allocate a thread local hash table
   hsh = par[ PthIdx ].HshTab = malloc(7 * siz * sizeof(HshSct));
   assert(hsh);

   // Clear the hash table's direct entries,
   // there is no need to clear the collision entries
   memset(hsh, 0, siz * sizeof(HshSct));

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

void ParEdg2(itg BegIdx, itg EndIdx, int PthIdx, ParSct *par)
{
   itg i, key, PthNmbEdg = 0, edg[ MAXEDG ][2];
   itg siz = par[ PthIdx ].HshSiz, NmbCpu = par[ PthIdx ].NmbCpu;
   int NmbEdg, flg, j, k;
   HshSct *hsh = par[ PthIdx ].HshTab, *buc;
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

#endif
