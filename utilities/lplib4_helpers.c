

/*----------------------------------------------------------------------------*/
/*                                                                            */
/*                               LPlib Helpers V1.0                           */
/*                                                                            */
/*----------------------------------------------------------------------------*/
/*                                                                            */
/* Description:         lplib's helper functions' headers                     */
/* Author:              Loic MARECHAL                                         */
/* Creation date:       may 16 2024                                           */
/* Last modification:   jul 25 2025                                           */
/*                                                                            */
/*----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------*/
/* Includes                                                                   */
/*----------------------------------------------------------------------------*/

#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <stdarg.h>

#ifdef WITH_METIS
#include <metis.h>
#include <math.h>
#endif

#include "lplib4.h"
#include "lplib4_helpers.h"


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

#ifdef WITH_METIS
typedef struct
{
   idx_t nvtxs, nedges, ncon, nparts, *xadj, *adjncy, *VerDeg, objval, *part;
   idx_t *adjwgt, options[ METIS_NOPTIONS ];
}MtsSct;
#endif


/*----------------------------------------------------------------------------*/
/* Prototypes of local procedures                                             */
/*----------------------------------------------------------------------------*/

void           ParEdg1        (itg, itg, int, ParSct *);
void           ParEdg2        (itg, itg, int, ParSct *);
static void    SetBnbBox      (RenSct *);
static double *SetMidCrd      (RenSct *, int);

#ifdef WITH_METIS
static void    SetBal         (LplMshSct *);
static int     GetBal         (LplMshSct *, int, int *, int *);
static void    GetMtsRef      (LplMshSct *, MtsSct *);
static void    BuildMetisGraph(LplMshSct *, MtsSct *);
#endif


/*----------------------------------------------------------------------------*/
/* Global tables                                                              */
/*----------------------------------------------------------------------------*/

static const int tvpe[6][2] = { {0,1}, {1,2}, {2,0}, {3,0}, {3,1}, {3,2} };
static const int EleSiz[ LplMax ] = {0,2,3,4,4,5,6,8};


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


/*----------------------------------------------------------------------------*/
/* Read, renumber through a Hilbert SFC and write the mesh                    */
/*----------------------------------------------------------------------------*/

RenSct *MeshRenumbering(int dim, int NmbVer, double *CrdTab, ...)
{
   int      i, j, t, siz, NmbCpu = 0, *TmpEle;
   int64_t  LibParIdx;
   double   *EleCrd, *TmpCrd;
   RenSct   *ren;
   va_list  VarArg;

   // Check mandatory inputs
   if( (dim != 3) || (NmbVer < 1) || !CrdTab)
      return(NULL);

   // Allocate and setup the renumbering structure
   if(!(ren = calloc(1, sizeof(RenSct))))
      return(NULL);

   ren->dim = dim;
   ren->NmbVer = NmbVer;
   ren->CrdTab = CrdTab;
   SetBnbBox(ren);

   // Decode the variable list of arguments
   va_start(VarArg, CrdTab);

   do
   {
      t = va_arg(VarArg, int);

      if( (t > 0) && (t < LplMax) )
      {
         ren->NmbEle[t] = va_arg(VarArg, int);
         ren->EleTab[t] = va_arg(VarArg, int *);
      }
   }while(t);

   va_end(VarArg);

   // Launch the required number of threads
   if(!(LibParIdx = InitParallel(NmbCpu)))
      return(NULL);


   // --------------------
   // Vertices renumbering
   // --------------------

   if(!(ren->RenTab[0] = malloc((ren->NmbVer+1) * 2 * sizeof(int64_t))))
      return(NULL);

   if(!HilbertRenumbering(LibParIdx, NmbVer, ren->box, (double (*)[3])ren->CrdTab, ren->RenTab[0]))
      return(NULL);


   // --------------------
   // Elements renumbering
   // --------------------

   for(t=LplEdg;t<LplMax;t++)
   {
      if(!ren->NmbEle[t])
         continue;

      if(!(ren->RenTab[t] = malloc( (ren->NmbEle[t] + 1) * 2 * sizeof(int64_t) )))
         return(NULL);

      if(!(EleCrd = SetMidCrd(ren, t)))
         return(NULL);

      if(!HilbertRenumbering(LibParIdx, ren->NmbEle[t], ren->box, (double (*)[3])EleCrd, ren->RenTab[t]))
         return(NULL);

      free(EleCrd);
      siz = EleSiz[t];

      if(!(TmpEle = malloc( (ren->NmbEle[t] + 1) * siz * sizeof(int) )))
         return(NULL);

      for(i=1;i<=ren->NmbEle[t];i++)
         for(j=0;j<siz;j++)
            TmpEle[ i * siz + j ] = ren->RenTab[0][ ren->EleTab[t][ ren->RenTab[t][i][1] * siz + j ] ][0];

      memcpy(ren->EleTab[t], TmpEle, (ren->NmbEle[t] + 1) * siz * sizeof(int));
      free(TmpEle);
   }

   // Move coordinates

   if(!(TmpCrd = malloc( (ren->NmbVer+1) * 3 * sizeof(double) )))
      return(NULL);

   for(i=1;i<=ren->NmbVer;i++)
      for(j=0;j<3;j++)
         TmpCrd[ i*3 + j ] = ren->CrdTab[ ren->RenTab[0][i][1] * 3 + j ];

   memcpy(ren->CrdTab, TmpCrd, (ren->NmbVer+1) * 3 * sizeof(double));
   free(TmpCrd);

   StopParallel(LibParIdx);

   return(ren);
}


/*----------------------------------------------------------------------------*/
/* Free all elements' numbering tables and the global structure itself        */
/*----------------------------------------------------------------------------*/

void FreeNumberingStruct(RenSct *ren)
{
   int t;

   for(t=0;t<LplMax;t++)
      if(ren->RenTab[t])
         free(ren->RenTab[t]);

   free(ren);
}


/*----------------------------------------------------------------------------*/
/* Evaluate elements' numbering ranging from 0 (perfect) to 1 (random)        */
/*----------------------------------------------------------------------------*/

double EvaluateRenumbering(int EleTyp, int NmbEle, int *EleTab)
{
   int      i, MinVer = EleTab[3], MaxVer = EleTab[3];
   double   dlt, NumQal = 0.;

   for(i=3;i<NmbEle * EleSiz[ EleTyp ];i++)
   {
      if(EleTab[i] != EleTab[ i+1 ])
      {
         dlt = (double)EleTab[i] - (double)EleTab[ i+1 ];
         NumQal += POW(dlt);
      }

      MinVer = MIN(MinVer, EleTab[i]);
      MaxVer = MAX(MaxVer, EleTab[i]);
   }

   dlt = (double)MaxVer - (double)MinVer;
   NumQal /= (NmbEle * EleSiz[ EleTyp ] * POW(dlt));

   return(NumQal);
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/

int RestoreNumbering(RenSct *ren, int NmbVer, double *CrdTab, ... )
{
   return(0);
}


/*----------------------------------------------------------------------------*/
/* Get an entity's index in the old numbering from its current one            */
/*----------------------------------------------------------------------------*/

int GetOldIndex(RenSct *ren, int typ, int idx)
{
   return(ren->RenTab[ typ ][ idx ][0]);
}


/*----------------------------------------------------------------------------*/
/* Get an entity's index in the current numbering from its old one            */
/*----------------------------------------------------------------------------*/

int GetNewIndex(RenSct *ren, int typ, int idx)
{
   return(ren->RenTab[ typ ][ idx ][1]);
}


/*----------------------------------------------------------------------------*/
/* Compute a mesh's bounding box                                              */
/*----------------------------------------------------------------------------*/

static void SetBnbBox(RenSct *ren)
{
   int i, j;

   ren->box[0] = ren->box[3] = ren->CrdTab[3];
   ren->box[1] = ren->box[4] = ren->CrdTab[4];
   ren->box[2] = ren->box[5] = ren->CrdTab[5];

   for(i=1;i<=ren->NmbVer;i++)
      for(j=0;j<3;j++)
      {
         ren->box[j  ] = MIN(ren->box[j  ], ren->CrdTab[ i*3 + j ]);
         ren->box[j+3] = MAX(ren->box[j+3], ren->CrdTab[ i*3 + j ]);
      }
}


/*----------------------------------------------------------------------------*/
/* Compute the barycenter of any kind of element                              */
/*----------------------------------------------------------------------------*/

static double *SetMidCrd(RenSct *ren, int typ)
{
   int i, j, k, siz = EleSiz[ typ ];
   double *crd;

   if(!(crd = malloc( (ren->NmbEle[ typ ] + 1) * 3 * sizeof(double) )))
      return(NULL);

   for(i=1;i<=ren->NmbEle[ typ ];i++)
   {
      for(j=0;j<3;j++)
         crd[ i * 3 + j ] = 0.;

      for(j=0;j<siz;j++)
         for(k=0;k<3;k++)
            crd[ i * 3 + k ] += ren->CrdTab[ ren->EleTab[ typ ][ i * siz + j ] * 3 + k ];

      for(j=0;j<3;j++)
         crd[ i * 3 + j ] /= siz;
   }

   return(crd);
}


#ifdef WITH_METIS

int MetisPartitioning(LplMshSct *msh, int NmbPar)
{
   MtsSct mts;

   SetBal(msh);
   BuildMetisGraph(msh, &mts);
   mts.nparts = NmbPar;

   if(METIS_PartGraphKway( &mts.nvtxs, &mts.ncon, mts.xadj, mts.adjncy,
                           NULL, NULL, NULL, &mts.nparts, NULL, NULL,
                           mts.options, &mts.objval, mts.part ) != METIS_OK)
   {
      return(0);
   }

   GetMtsRef(msh, &mts);

   return(1);
}


/*----------------------------------------------------------------------------*/
/* Extract a metis graph from a tet mesh                                      */
/*----------------------------------------------------------------------------*/

static void BuildMetisGraph(LplMshSct *msh, MtsSct *mts)
{
   int i, j, BalTab[1000], WgtTab[1000], TotDeg = 0;

   METIS_SetDefaultOptions(mts->options);

   mts->options[ METIS_OPTION_CTYPE ]     = METIS_CTYPE_RM;
   mts->options[ METIS_OPTION_CTYPE ] = 0;
   mts->options[ METIS_OPTION_CONTIG ]    = 1;
   mts->options[ METIS_OPTION_OBJTYPE ]   = METIS_OBJTYPE_VOL;
   //mts->options[ METIS_OPTION_DBGLVL ]    = METIS_DBG_INFO;
   mts->options[ METIS_OPTION_IPTYPE ]    = METIS_IPTYPE_RANDOM;

   for(i=1;i<=msh->NmbVer;i++)
      TotDeg += GetBal(msh, i, BalTab, WgtTab);

   mts->nvtxs  = msh->NmbVer;
   mts->nedges = TotDeg;
   mts->ncon   = 1;
   mts->VerDeg = malloc(msh->NmbVer * sizeof(int));
   mts->xadj   = malloc( (msh->NmbVer + 1) * sizeof(idx_t));
   mts->part   = malloc(msh->NmbVer * sizeof(idx_t));
   mts->adjncy = malloc(TotDeg      * sizeof(idx_t));
   mts->adjwgt = malloc(TotDeg      * sizeof(idx_t));

   if(!mts->xadj || !mts->VerDeg || ! mts->adjncy || !mts->part)
      exit(1);

   TotDeg = 0;

   for(i=0;i<msh->NmbVer;i++)
   {
      mts->xadj[i] = TotDeg;
      mts->VerDeg[i] = GetBal(msh, i+1, &mts->adjncy[ TotDeg ], &mts->adjwgt[ TotDeg ]);
      TotDeg += mts->VerDeg[i];
   }

   mts->xadj[ msh->NmbVer ] = TotDeg;
}


/*----------------------------------------------------------------------------*/
/* Extract the ball of unique vertices from the ball of tets                  */
/*----------------------------------------------------------------------------*/

static int GetBal(LplMshSct *msh, int VerIdx, int *VerTab, int *WgtTab)
{
   int i, j, k, flg, NmbVer = 0, *TetVer, wgt;
   double siz;

   for(i=0;i<msh->VerDeg[ VerIdx ];i++)
   {
      TetVer = msh->TetTab[ msh->BalTab[ msh->VerBal[ VerIdx ] + i ][0] ];

      for(j=0;j<4;j++)
      {
         if(TetVer[j] == VerIdx)
            continue;

         flg = 0;

         for(k=0;k<NmbVer;k++)
            if(VerTab[k] == TetVer[j])
            {
               flg = 1;
               break;
            }

         if(!flg)
            VerTab[ NmbVer++ ] = TetVer[j];
      }
   }

   for(i=0;i<NmbVer;i++)
   {
      siz = POW(msh->CrdTab[ VerIdx ][0] - msh->CrdTab[ VerTab[i] ][0])
          + POW(msh->CrdTab[ VerIdx ][1] - msh->CrdTab[ VerTab[i] ][1])
          + POW(msh->CrdTab[ VerIdx ][2] - msh->CrdTab[ VerTab[i] ][2]);

      siz = sqrt(siz);
      wgt = siz / msh->MinSiz;
      //wgt = msh->MaxSiz / siz;

      VerTab[i]--;
      WgtTab[i] = wgt;
   }

   return(NmbVer);
}


/*----------------------------------------------------------------------------*/
/* Setup the balls table                                                      */
/*----------------------------------------------------------------------------*/

static void SetBal(LplMshSct *msh)
{
   int i, j, VerIdx, TabSiz = 0;

   msh->VerDeg = calloc(msh->NmbVer + 1, sizeof(int));
   msh->VerBal = calloc(msh->NmbVer + 1, sizeof(int));

   // Set vertices degree
   for(i=1;i<=msh->NmbTet;i++)
      for(j=0;j<4;j++)
         msh->VerDeg[ msh->TetTab[i][j] ]++;

   for(i=1;i<=msh->NmbVer;i++)
      TabSiz += msh->VerDeg[i];

   // Allocate the global balls table and give each vertex a pointer to its own table
   msh->BalTab = malloc(TabSiz * 2 * sizeof(int));
   TabSiz = 0;

   if(!msh->VerDeg || !msh->VerBal || !msh->BalTab)
   {
      puts("Failed to allocate memory");
      exit(1);
   }

   for(i=1;i<=msh->NmbVer;i++)
   {
      msh->VerBal[i] = TabSiz;
      TabSiz += msh->VerDeg[i];
      msh->VerDeg[i] = 0;
   }

   // Fill the ball tables with elements type and index
   for(i=1;i<=msh->NmbTet;i++)
      for(j=0;j<4;j++)
      {
         VerIdx = msh->TetTab[i][j];
         msh->BalTab[ msh->VerBal[ VerIdx ] + msh->VerDeg[ VerIdx ] ][0] = i;
         msh->BalTab[ msh->VerBal[ VerIdx ] + msh->VerDeg[ VerIdx ] ][1] = j;
         msh->VerDeg[ VerIdx ]++;
      }
}


/*----------------------------------------------------------------------------*/
/* Set mesh vertices and tets' ref from the Metis partitions                  */
/*----------------------------------------------------------------------------*/

static void GetMtsRef(LplMshSct *msh, MtsSct *mts)
{
   int i, j, ref;

   for(i=0;i<msh->NmbVer;i++)
      msh->RefTab[i+1] = mts->part[i];

   for(i=1;i<=msh->NmbTet;i++)
   {
      ref = msh->RefTab[ msh->TetTab[i][0] ];

      for(j=1;j<4;j++)
         if(msh->RefTab[ msh->TetTab[i][j] ] < ref)
            ref = msh->RefTab[ msh->TetTab[i][j] ];

      msh->TetTab[i][4] = ref;
   }
}
#endif
