

/*----------------------------------------------------------------------------*/
/*                                                                            */
/*                               LPlib Helpers V1.1                           */
/*                                                                            */
/*----------------------------------------------------------------------------*/
/*                                                                            */
/* Description:         lplib's helper functions' headers                     */
/* Author:              Loic MARECHAL                                         */
/* Creation date:       may 16 2024                                           */
/* Last modification:   aug 14 2025                                           */
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

typedef itg int1d;
typedef itg int2d[2];
typedef itg int11d[11];
typedef itg int56d[56];


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


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/

int Wlf_ComparVer1d(const void *v1, const void *v2)
{
  int *r1 = ( int *)v1;
  int *r2 = ( int *)v2;

	int rank1 = (*r1);
	int rank2 = (*r2);

	if ( rank1 > rank2 )
		return 1;
	else
		return -1;
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/

int1d Wlf2_CheckIfVertexIsInTheList(int1d iVer, int1d nbv, int1d *verLst)
{
  int1d k; 
  
  for (k=0; k<nbv; k++) {
    if ( verLst[k] == iVer ) 
      return k;
  }
  
  return -1;
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/

void Wlf_AddHeapListCol(int1d iVer, int1d *NbrHepLst, int1d *HepLst, int1d *Pos, int2d *ColVer2Ver)
{
  int1d son, fat, jVer;

  if ( Pos[iVer] > 0 )
    return;
  
  (*NbrHepLst)++;
  son = *NbrHepLst;
  //--- reorder the list
  while ( son > 1 ) {
    fat  = son/2;  // son>>1;  
    jVer = HepLst[fat];   
    if ( (ColVer2Ver[iVer][0] > ColVer2Ver[jVer][0]) || (ColVer2Ver[iVer][0] == ColVer2Ver[jVer][0] && ColVer2Ver[iVer][1] > ColVer2Ver[jVer][1]) ) {   // switch father and son
      HepLst[son] = jVer;
      Pos[jVer]   = son;
      son         = fat;
    }
    else {    // set vertex at is position
      HepLst[son] = iVer;
      Pos[iVer]   = son;
      return;
    } 
  }

  if ( son != 1 ) {
    fprintf(stderr, "  ## ERROR Red2_AddHeapList: Add heap list\n");
    exit(1);
  }
  
  HepLst[son] = iVer;
  Pos[iVer]   = son;

  return;
}





/*

  Remove the top vertex from the list
   - First, remove the top vertex
   - Second, position the last vertex at the top of the list
   - Third, traverse the list from top downwards to reposition the first element 

*/
int1d Wlf_RemoveHeapListCol(int1d *NbrHepLst, int1d *HepLst, int1d *Pos, int2d *ColVer2Ver)
{
  int1d  iVer, jVer, lVer, son1, son2, fat;

  //--- Remove the top vertex
  iVer      = HepLst[1];
  Pos[iVer] = -1;
  
  //--- Put the last vertex at the top of the list
  jVer = HepLst[*NbrHepLst];
  fat  = 1;
  (*NbrHepLst)--;

  //--- Traverse the list from top downwards to reposition the first element
  while ( (*NbrHepLst) > 2*fat ) {
    son1 = 2*fat;
    son2 = son1 + 1;
    if ( (ColVer2Ver[HepLst[son1]][0] > ColVer2Ver[HepLst[son2]][0]) || (ColVer2Ver[HepLst[son1]][0] == ColVer2Ver[HepLst[son2]][0] && ColVer2Ver[HepLst[son1]][1] > ColVer2Ver[HepLst[son2]][1])) {
      if ( (ColVer2Ver[HepLst[son1]][0] > ColVer2Ver[jVer][0]) || (ColVer2Ver[HepLst[son1]][0] == ColVer2Ver[jVer][0] && ColVer2Ver[HepLst[son1]][1] > ColVer2Ver[jVer][1])){
        lVer        = HepLst[son1];
        HepLst[fat] = lVer;
        Pos[lVer]   = fat;
        fat         = son1;     
      }
      else {   // list is updated
        HepLst[fat] = jVer;
        Pos[jVer]   = fat;
        return iVer;
      }
    }
    else {
      if ( (ColVer2Ver[HepLst[son2]][0] > ColVer2Ver[jVer][0]) || (ColVer2Ver[HepLst[son2]][0] == ColVer2Ver[jVer][0] && ColVer2Ver[HepLst[son2]][1] > ColVer2Ver[jVer][1])){    // switch jVer and son2 (Peut-etre fat and son2?)
        lVer        = HepLst[son2];
        HepLst[fat] = lVer;
        Pos[lVer]   = fat;
        fat         = son2;
      }
      else {   // list is updated
        HepLst[fat] = jVer;
        Pos[jVer]   = fat;
        return iVer;
      }
    }
  }
  
  //--- I am at the bottom of the list: v take this position
  HepLst[fat] = jVer;
  Pos[jVer]   = fat;

  return iVer;
}




/*

  Update the position of a vertex in the heap list if its distance is modified

*/
void Wlf_UpdateHeapListCol(int1d iVer, int1d *HepLst, int1d *Pos, int2d *ColVer2Ver)
{
  int1d  jVer, son, fat;
  son = Pos[iVer];
   
  while ( son > 1 ) {
    fat  = son/2;
    jVer = HepLst[fat];
    
    if ( (ColVer2Ver[iVer][0] > ColVer2Ver[jVer][0]) || (ColVer2Ver[iVer][0] == ColVer2Ver[jVer][0] && ColVer2Ver[iVer][1] > ColVer2Ver[jVer][1])) {   // switch father and son
      HepLst[son] = jVer;
      Pos[jVer]   = son;
      son         = fat;
    }
    else {    // set vertex at is position
      HepLst[son] = iVer;
      Pos[iVer]   = son;
      return;
    } 
  }

  //--- I am at the top of the list ie son = 1
  if ( son != 1 ) {
    fprintf(stderr, "  ## ERROR Wlf_UpdateHeapListCol: Update heap list son = %d \n", son);
    exit(1);
  }

  HepLst[son] = iVer;
  Pos[iVer]   = son;

  return;    
}


void Wlf_ComputeNeighborsRank2_3d(int *PtrVer2Ver, int *Ver2Ver, int1d iVer, int1d **VerLst, int1d *NbrEdgRk2)
{
  int1d   jVer, iVoi, iVfr, nbv, nbr, iLay, cptLow, cptUpp, beg, end, cpt, nbrLow, nbrUpp, idxLow, idxUpp, nbrVoiRk2;
  int1d   i, ive, idx, flg, sizLst, ref, nbvRkLow, nbvRkUpp, cptLay =0, NbrLay;
    int1d  *verLst;
	
  //--- Init
  sizLst = 512;
  verLst = (int1d *)malloc(sizeof(int1d)*(sizLst));
  memset(verLst, 0, sizeof(int1d)*512);
    
  nbv    = 0;
  NbrLay = 1;

  verLst[0] = iVer;  // Put the first to not treat it the vertices of Rk 2 BUT don't put it in the final list
  
  nbv++;

  cpt = 0;
  
  //--- Vertices neighbors of rank 1: for each vertex get its ball of elements  
  nbrVoiRk2 = 0;

  jVer = iVer;
  nbr  = PtrVer2Ver[jVer+1] - PtrVer2Ver[jVer];
  
  for (i=0; i<nbr; i++) {
    idx  = PtrVer2Ver[jVer] + i;
    iVoi = Ver2Ver[idx];

    
    //- add new vertex neighbor
    verLst[nbv] = iVoi;
    nbv++;

    cpt++;
   
    if ( nbv >= sizLst ) {
      sizLst   *= 2;
      verLst = (int1d *)realloc(verLst, sizeof(int1d)*(sizLst));
    }      
  }

  beg = 0;
  end = nbv;
  //--- Vertices neighbors of rank 2: for each vertex get its ball of elements
  
  for(iLay = 0 ; iLay<NbrLay ; iLay++ ) {
    for (ive=beg; ive<end; ++ive) {  
      jVer = verLst[ive];
      nbr  = PtrVer2Ver[jVer+1] - PtrVer2Ver[jVer];
    
      for (i=0; i<nbr; i++) {
        idx  = PtrVer2Ver[jVer] + i;
        iVoi = Ver2Ver[idx];

        flg = Wlf2_CheckIfVertexIsInTheList(iVoi, nbv, verLst);

        if ( flg >= 0 ) 
          continue;
        
          //- add new vertex neighbor
        verLst[nbv] = iVoi;
        nbv++;
        cptLay++; // count the numver of vertices addes in that layer for computing the neigboors oh the same layer for the enxt pass

        cpt++;
     
        if ( nbv >= sizLst ) {
          sizLst *= 2;
          verLst = (int1d *)realloc(verLst, sizeof(int1d)*(sizLst));
            
        }      
      }
    }
    beg = end;
    end += cptLay;
    cptLay=0;
  }

  *NbrEdgRk2 = nbv-1;

  VerLst[iVer] = (int1d*)malloc(sizeof(int1d)*cpt);  
  for (i=1; i<nbv; i++){   
    VerLst[iVer][i-1] = verLst[i];

  }
  free(verLst);
  verLst = NULL;
      
  return;


}


void Wlf_InitializeWolfNscDataBaseNeighboorRank2(int NmbVer, int NmbEdg, int *EdgTab, int *PtrVer2Ver, int *Ver2Ver)
{
  int iluLvl;
  int iVer, iEdg;
  int1d delete, periodic;
  int2d nbvRk1;
  int1d  *NbrVoiRk2 = (int1d*)malloc(sizeof(int1d)*(NmbVer+1));
   int1d **LstVoiRk2 = (int1d**)malloc(sizeof(int1d*)*(NmbVer+1));

    NbrVoiRk2[0] = 0;

    for (iVer=1 ; iVer<=NmbVer ; iVer ++) {
      Wlf_ComputeNeighborsRank2_3d(PtrVer2Ver, Ver2Ver, iVer, LstVoiRk2, &(NbrVoiRk2[iVer]));
    }
}

void Wlf_ColoringPartitionImplicit3dNew(int NmbPth, int nbrPar, int NmbVer, int *PtrVer2Ver, int *Ver2Ver, int *NbrVoiRk2, int56d *LstVoiRk2,int *VidPar)
{
  int1d       NbrCol, i, ColUpp, ColLow, NbrPar = 64, addCol =0, flag1, flag2, flag=0;
  int1d       iVer, iCol, jCol, iPar, jPar,  iVoi, idx, idxPar, jdxPar, tgtPar , jVer, delete, periodic, nbrVoi, degMaxCol=0, degMax=0;
  int1d       idxCol, jdxCol, NbrHepLst, NbrParTgt, IteBck, test, colRef;
  int1d       *cntPar =NULL, *colPar=NULL, *cntColPar=NULL;
  int1d       *cntParSrt=NULL, *cntParColVoi=NULL, *cntVoiNoCol=NULL, *LstCol=NULL, *LstColUpp=NULL, *LstColLow=NULL ;
  int1d       NbrTyp = 2;
  double      *cntCol = NULL , *tagDbl = NULL;
  int1d   *HepLst, *Pos ,*bufInt, *GphColPar; 
  int2d   *ColVer2Ver;
  int56d  CntCol;
  int11d  BackUpCol;
  NbrParTgt = NmbPth;


  tagDbl    = (double *)calloc(NmbVer+1, sizeof(double));

  for (i=0;i<30;i++) {
    CntCol[i]    = 0;
  }

  bufInt     = (int1d *)malloc((nbrPar+1)*sizeof(int1d));
  HepLst     = (int1d *)malloc((nbrPar+1)*sizeof(int1d));
  Pos        = (int1d *)malloc((nbrPar+1)*sizeof(int1d));
  ColVer2Ver = (int2d *)malloc((nbrPar+1)*sizeof(int2d));
  cntPar             = (int1d *)calloc(nbrPar+1, sizeof(int1d));
  cntColPar          = (int1d *)calloc(nbrPar+1, sizeof(int1d));
  colPar             = (int1d *)calloc((nbrPar+1)*1000, sizeof(int1d));
  GphColPar        = (int1d *)calloc((nbrPar+1),sizeof(int1d));
  cntParSrt          = (int1d *)calloc(2*(nbrPar+1), sizeof(int1d));  // use to order the partition according to their increasing degree, second int stock the idPar after ordering

  // Treat the Sub domains like Vertex in the point implicit then give the color to each vertex of each partition

    for(iVer=1; iVer<=NmbVer; iVer++) {
      nbrVoi = NbrVoiRk2[iVer];
      idxPar = VidPar[iVer];

      for ( iVoi = 0 ; iVoi < nbrVoi ; iVoi++) {
        jVer = LstVoiRk2[iVer][iVoi];
        jdxPar = VidPar[jVer];

        if ( idxPar==jdxPar ) continue;

        else {
          flag1 = Wlf2_CheckIfVertexIsInTheList(jdxPar, cntPar[idxPar], &colPar[1000*idxPar]); // flag1 = -1 if jdxPar is not in &colPar[1000*idxPar] ()

          if ( cntPar[idxPar]>1000 || cntPar[jdxPar]>1000 ) {
            printf("MESH IS TO COARSE TO BE PARTITIONNED AND COLORED PROPERLY \n");
            exit(1);
          }

          if ( flag1==-1) {
            colPar[1000*idxPar+cntPar[idxPar]]=jdxPar;
            cntPar[idxPar]++;
          }

          else{
            continue;
          }
        }
      }
    }
  cntParColVoi            = (int1d *)calloc(nbrPar+1, sizeof(int1d));
  cntVoiNoCol             = (int1d *)calloc(nbrPar+1, sizeof(int1d));
  LstCol                  = (int1d *)calloc(100, sizeof(int1d));

  NbrHepLst=0;
  idx=0;
  NbrCol=0;

  for (iPar = 1 ; iPar<=nbrPar; iPar++) {
    nbrVoi = cntPar[iPar];
    ColVer2Ver[iPar][0] = 0;   // Nbr Voisins coloriés
    ColVer2Ver[iPar][1] = nbrVoi; // Personne n'est colorié
    Pos[iPar]=0;

    if (nbrVoi > idx ) {
      flag = iPar;
      idx = nbrVoi;
    }
  }

  Wlf_AddHeapListCol(flag, &NbrHepLst, HepLst, Pos, ColVer2Ver);

  for (iPar = 1 ; iPar<=nbrPar; iPar++) {
    iCol = 0;
    jPar = HepLst[1];

    Wlf_RemoveHeapListCol(&NbrHepLst, HepLst, Pos, ColVer2Ver);

    nbrVoi = cntPar[jPar];

    if (nbrVoi > 500 ){
      printf("Maillage de merde \n");
      exit(1);
    }

    for ( iVoi = 0 ; iVoi < nbrVoi ; iVoi++) {
      idxPar  = colPar[1000*jPar + iVoi];
      ColVer2Ver[idxPar][0]++;
      // ColVer2Ver[kVer][0]++;
      ColVer2Ver[idxPar][1]--;

      if (ColVer2Ver[idxPar][1] <0) {
        printf("ON A COLORIÉ TROP DE VOISINS A idxPar = %d\n", idxPar);
        exit(1);
      }
      // printf(" \t iVoi = %d kVer = %d",iVoi, kVer);
      if (ColVer2Ver[idxPar][0] == 1 ){
        Wlf_AddHeapListCol(idxPar, &NbrHepLst, HepLst, Pos, ColVer2Ver);
      }
      if (ColVer2Ver[idxPar][0] > 1 && GphColPar[idxPar] == 0) {
        Wlf_UpdateHeapListCol(idxPar, HepLst, Pos, ColVer2Ver);
      }

      LstCol[iVoi] = GphColPar[idxPar];


    }

    int1d nbrMin = nbrPar * 2;
    for (i = 1; i<30; i++){
      if ( (CntCol[i] != 0 )|| i<=10) {
      // if ( i<=9) {
        iVoi = 0;
        jCol = i;
        
        while ( iVoi < nbrVoi ) {
          if ( LstCol[iVoi] == jCol) // On skip la boucle car 'jCol est déja pris
            goto nexCol;
          //else if ( LstCol[iVoi] < jCol && LstCol[iVoi+1] > jCol ) // On a trouvé que la couleur 'flag' c'est bon
          //  iCol = jCol;
          //else
          iVoi++;
        }

        // this color is available
        if ( iCol == 0 ) {
          iCol = jCol;
          nbrMin = CntCol[i];
        } 
        else {
          if ( CntCol[i] < nbrMin ) {
            iCol   = jCol;
            nbrMin = CntCol[i];
          }
        }
      }

      nexCol:
      continue;
    }



    if (iCol == 0) { // On a fait tous les voisins et toutes les couleurs consécutives sont données
      iCol = LstCol[nbrVoi-1] +  1;
    }

    GphColPar[jPar] = iCol;

    CntCol[iCol]++;

    // if (iCol == 16 || iCol == 2) printf(" \t %%%%%%%% Color %d for iVer %d in total %d \n ", iCol , jVer ,CntCol[iCol]);

    if (iCol > NbrCol){
      NbrCol = iCol;
    }
  }

  //---- Check up of evrything
  flag = 0;
  colRef = 0;
  for( iCol =1 ; iCol<=NbrCol ; iCol++){
    BackUpCol[iCol] = CntCol[iCol] - NbrParTgt;
    if (BackUpCol[iCol] != 0) {
      flag++; //Coutn for the number of problem (must be even)
      printf("\t Color %d has a wrong number of partition %d (Target = %d) \t", iCol, CntCol[iCol], NbrParTgt);      
    }
  }

  // printf("flag = %d so we have %d partitions to transfer \n", flag, flag/2);
  if (flag == 0) IteBck = 0;
  else IteBck = 5;
  
  while ( flag > 0 && IteBck > 0 ) {
    printf("flag = %d IteBck = %d \n", flag, IteBck);
    idxCol =0;
    jdxCol = 0;
    for( iCol =1 ; iCol<=NbrCol ; iCol++) {
      //--- find idxCol such that idxCol has to less and jdxCol to much partitions for a transfer
      if ( BackUpCol[iCol] < 0 && idxCol == 0 && iCol != colRef ) idxCol = iCol; 

      if ( BackUpCol[iCol] > 0 && jdxCol == 0) jdxCol = iCol; 
  
    }
    //--if we find a couple of color for transfer lets try
    if ( idxCol != 0 && jdxCol != 0 ) {  
      colRef =  idxCol; //-- to prevent from treating the same color again and again
      //-- Looking for a partition for the transfer
      for (iPar = 1 ; iPar<=nbrPar; iPar++) { 
        test = 0; 
        nbrVoi = cntPar[iPar];
        iCol = GphColPar[iPar];
        

        if ( iCol == jdxCol ) { //-- iPar is a partition that is in color with to many partition : it is candidate for transfer, test for jdxCol
          // printf("Voisins de %d : ",iPar);
          for ( iVoi = 0 ; iVoi < nbrVoi ; iVoi++){
            jPar  = colPar[1000*(iPar)+iVoi];
            jCol = GphColPar[jPar];
            // printf("\t%% %d col : %d", jPar, jCol);
            if ( jCol == idxCol ) 
              test = 1;
          }
          // printf("\n");
          if ( test == 0 ) { // Test is good iPar is a champion !!!
            CntCol[idxCol]++;
            CntCol[jdxCol]--;
            GphColPar[iPar] = idxCol;

            BackUpCol[idxCol]++;
            BackUpCol[jdxCol]--;

            flag=flag-2;
            goto end;
            // BAD IDEA !!!
            // if (BackUpCol[idxCol]== 0  && BackUpCol[jdxCol]== 0 ) // A transfer managed to balance 2 colors
            //   flag=flag-2;
          } 
            
        }
      }
    }
    end:
    IteBck--;
  }

  if ( flag !=0 ) printf("Balance is not good\n");

  return;
}



#endif
