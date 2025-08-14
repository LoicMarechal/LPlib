

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
#include <float.h>
#include <math.h>

#ifdef WITH_METIS
#include <metis.h>
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
#define MAXITR    21
#define MAXELE    14
#define HILMOD    0
#define OCTMOD    1
#define RNDMOD    2
enum {TypEdg, TypTri, TypQad, TypTet, TypPyr, TypPri, TypHex};
enum {HilMod=0, OctMod, RndMod, IniMod, TopMod};


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

typedef struct
{
   uint64_t cod;
   double   crd[3];
   int      idx, ref, col, grn, deg;
}VerSct;

typedef struct
{
   uint64_t cod;
   int      *idx, siz, col, grn, ref, gid;
}EleSct;

typedef struct
{
   int      NmbVer, *Old2New, MshVer, dim, mod, TypIdx, VerTyp, GmlMod;
   int      NmbEle[ MAXELE ], *IdxTab[ MAXELE ], EleTyp[ MAXELE ];
   int      MaxDeg[ MAXELE ], DegVec[ MAXELE ], HghDeg, OvrDeg;
   int      ColGrnFlg, ColGrnMod, NmbCol, NmbGrn;
   int      (*ColPar)[3], (*VerGrnPar)[ MAXELE ][4], (*EleGrnPar)[ MAXELE ][2];
   int      ColBit, GrnBit, DegBit, RefBit, VerHilBit, FacHilBit, VolHilBit;
   int      ColLft, GrnLft, DegLft, RefLft, VerHilRgt, FacHilRgt, VolHilRgt;
   uint64_t ColMsk, GrnMsk, DegMsk, RefMsk;
   double   box[6];
   VerSct   *ver;
   EleSct   *ele[ MAXELE ];
}MshSct;

#ifdef WITH_METIS
typedef struct
{
   idx_t nvtxs, nedges, ncon, nparts, *xadj, *adjncy, *VerDeg, objval, *part;
   idx_t *adjwgt, options[ METIS_NOPTIONS ];
}MtsSct;
#endif


/*----------------------------------------------------------------------------*/
/* Global variables                                                           */
/*----------------------------------------------------------------------------*/

// For each kind of element: number of node, number of faces, GMF keyword
int EleTab[ MAXELE ][3] = {
   { 2, 2, 0},
   { 3, 3, 0},
   { 4, 4, 0},
   { 4, 4, 0},
   { 5, 5, 0},
   { 6, 5, 0},
   { 8, 6, 0},
   { 3, 2, 0},
   { 6, 3, 0},
   { 9, 4, 0},
   {10, 4, 0},
   {14, 5, 0},
   {18, 5, 0},
   {27, 6, 0} };

// For each kind of element: GMF keyword strings
char *EleNam[ MAXELE ] = {
   "Edges           ",
   "Triangles       ",
   "Quadrilaterals  ",
   "Tetrahedra      ",
   "Pyramids        ",
   "Prisms          ",
   "Hexahedra       ",
   "EdgesP2         ",
   "TrianglesP2     ",
   "QuadrilateralsQ2",
   "TetrahedraP2    ",
   "PyramidsP2      ",
   "PrismsP2        ",
   "HexahedraQ2     " };

   // For each kind of element: low vertex degree and highest vertex degree
   int MaxDeg[ MAXELE ][2] = {
      { 2,   8},
      { 8,  32},
      { 4,  16},
      {32, 128},
      {16,  64},
      {16,  64},
      { 8,  32},
      { 2,   8},
      { 8,  32},
      { 4,  16},
      {32, 128},
      {16,  64},
      {16,  64},
      { 8,  32} };


/*----------------------------------------------------------------------------*/
/* Prototypes of local procedures                                             */
/*----------------------------------------------------------------------------*/

void           ParEdg1        (itg, itg, int, ParSct *);
void           ParEdg2        (itg, itg, int, ParSct *);
static void    SetBnbBox      (RenSct *);
static double *SetMidCrd1      (RenSct *, int);
void SetMidCrd2(int , int *, MshSct *, double *);

int      CmpFnc   (const void *, const void *);
void     RenVer   (int, int, int, MshSct *);
void     RenEle   (int, int, int, MshSct *);
void     SwpMem   (MshSct *, int);
void     SetVerDeg(MshSct *);
void     SetMatSlc(MshSct *);
int      SetDeg   (MshSct *, int *);
void     ChkColGrn(MshSct *);
uint64_t HilCod   (double *, double *, int, int);
int      CmpFnc   (const void *, const void *);

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

      if(!(EleCrd = SetMidCrd1(ren, t)))
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

static double *SetMidCrd1(RenSct *ren, int typ)
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

int Renum()
{
   char     *PtrArg, *TmpStr, *BndTab, InpNam[1000], OutNam[1000];
   char     MtsNodNam[1000], MtsEleNam[1000];
   char     *SchStr[4] = {"Hilbert", "Z-curve", "random", "initial"};
   int      i, j, t, NmbCpu=0, StaFlg=0, BndFlg=0, *IdxTab, NmbSrf, NmbVol;
   int      NmbGrn, CurGrn, NmbCol, CurCol, grn;
   int64_t  LibParIdx;
   double   timer = 0.;
   MshSct   msh;
   VerSct   *OldVer, *NewVer, *VolVer;

   memset(&msh, 0, sizeof(MshSct));

   // If the color renumbering is set but no colors or grains are present in the
   // input file, the mode is disabled and set to default hilbert renumbering
   if(msh.ColGrnMod && !msh.ColGrnFlg)
   {
      msh.ColGrnMod = 0;
      puts("Could not find colors and grains information: switching back to default renumbering");
   }


   // --------------------------------------------
   // Compute initial block dependenies statistics
   // --------------------------------------------

   LibParIdx = InitParallel(NmbCpu);


   // ---------------------------------------------------------
   // Set the size of bit fields needed by color and grain data
   // ---------------------------------------------------------

   if(msh.ColGrnMod)
   {
      msh.NmbCol = msh.NmbGrn = 0;

      for(i=1;i<=msh.NmbVer;i++)
      {
         msh.NmbCol = MAX(msh.NmbCol, msh.ver[i].col);
         msh.NmbGrn = MAX(msh.NmbGrn, msh.ver[i].grn);
      }

      printf("Found %d colors and %d grains in the input file\n", msh.NmbCol, msh.NmbGrn);

      for(t=0;t<MAXELE;t++)
      {
         for(i=1;i<=msh.NmbEle[t];i++)
         {
            msh.ele[t][i].col = msh.ver[ msh.ele[t][i].idx[0] ].col;
            msh.ele[t][i].grn = msh.ver[ msh.ele[t][i].idx[0] ].grn;
         }
      }

      ChkColGrn(&msh);

      // Set the bit field size, maks and left shift to store the color value
      msh.ColBit = ceil(log(msh.NmbCol) / log(2));
      msh.ColMsk = (1ULL << msh.ColBit) - 1ULL;
      msh.ColLft = 64 - msh.ColBit;

      // Set the bit field size, maks and left shift to store the grain value
      msh.GrnBit = ceil(log(msh.NmbGrn) / log(2));
      msh.GrnMsk = (1ULL << msh.GrnBit) - 1ULL;
      msh.GrnLft = 64 - msh.ColBit - msh.GrnBit;
   }


   // ----------------------------------------------------------------------
   // Prepare references and degree related data for the GMlib modified sort
   // ----------------------------------------------------------------------

   if(msh.GmlMod == 1)
   {
      SetVerDeg(&msh);

      // Only one bit is need to encode the high or low degree information
      msh.DegBit = 1;
      msh.DegMsk = 1ULL;
      msh.DegLft = 64 - msh.ColBit - msh.GrnBit - msh.DegBit;

      // Use the 8 leftmost bits to store face ref as boundary condition tag
      msh.RefBit = 8;
      msh.RefMsk = (1ULL << msh.RefBit) - 1ULL;
      msh.RefLft = 64 - msh.RefBit;
   }
   else if(msh.GmlMod == 2)
   {
      SetMatSlc(&msh);

      // three bits are needed to encode the 5 possible vertex degrees
      msh.DegBit = 3;
      msh.DegMsk = 7ULL;
      msh.RefMsk = (1ULL << msh.DegBit) - 1ULL;
      msh.DegLft = 64 - msh.ColBit - msh.GrnBit - msh.DegBit;
   }


   // ------------------------------------------------------
   // Finaly, set the Hilbert bit field size and right shift
   // ------------------------------------------------------

   // The hilbert code is right shifted with the number of bit
   // used by all other sorting keys
   msh.VerHilBit = 64 - msh.ColBit - msh.GrnBit - msh.DegBit;
   msh.VerHilRgt = msh.ColBit + msh.GrnBit + msh.DegBit;

   msh.FacHilBit = 64 - msh.ColBit - msh.GrnBit - msh.RefBit;
   msh.FacHilRgt = msh.ColBit + msh.GrnBit + msh.RefBit;

   msh.VolHilBit = 64 - msh.ColBit - msh.GrnBit;
   msh.VolHilRgt = msh.ColBit + msh.GrnBit;


   // --------------------
   // Vertices renumbering
   // --------------------

   printf("Sorting keys table: number of bit per key for each dimension of mesh entities\n");
   printf(" Entity | rank4 (color) | rank3 (grain) | rank2 (degree or ref) | rank1 (%s)\n",
            SchStr[ msh.mod ]);
   printf(" Vertex |      %2d       |      %2d       |           %2d          |      %2d\n",
            msh.ColBit, msh.GrnBit, msh.DegBit, msh.VerHilBit);

   printf(" Face   |      %2d       |      %2d       |           %2d          |      %2d\n",
            msh.ColBit, msh.GrnBit, msh.RefBit, msh.FacHilBit);

   printf(" Volume |      %2d       |      %2d       |           %2d          |      %2d\n",
            msh.ColBit, msh.GrnBit, 0, msh.VolHilBit);

   msh.VerTyp = NewType(LibParIdx, msh.NmbVer);
   LaunchParallel(LibParIdx, msh.VerTyp, 0, (void *)RenVer, (void *)&msh);
   ParallelQsort(LibParIdx, &msh.ver[1], msh.NmbVer, sizeof(VerSct), CmpFnc);

   msh.Old2New = malloc( (size_t)(msh.NmbVer+1) * sizeof(int) );
   assert(msh.Old2New);

   for(i=1;i<=msh.NmbVer;i++)
      msh.Old2New[ msh.ver[i].idx ] = i;



   // --------------------
   // Elements renumbering
   // --------------------

   for(t=0;t<MAXELE;t++)
      if(msh.NmbEle[t])
      {
         msh.EleTyp[t] = NewType(LibParIdx, msh.NmbEle[t]);
         msh.TypIdx = t;
         LaunchParallel(LibParIdx, msh.EleTyp[t], 0, (void *)RenEle, (void *)&msh);
         ParallelQsort(LibParIdx, &msh.ele[t][1], msh.NmbEle[t], sizeof(EleSct), CmpFnc);
         SwpMem(&msh, t);
      }


   // -----------------------------------------------------------------
   // Color-grain mode: setup global colouring and each entities grains
   // -----------------------------------------------------------------

   if(msh.ColGrnMod)
   {
      CurCol = msh.ver[1].col;
      NmbCol = 1;
      CurGrn = msh.ver[1].grn;
      NmbGrn = 1;

      for(i=1;i<=msh.NmbVer;i++)
      {
         if(msh.ver[i].col == CurCol)
            msh.ver[i].col = NmbCol;
         else
         {
            CurCol = msh.ver[i].col;
            msh.ver[i].col = ++NmbCol;
         }

         if(msh.ver[i].grn == CurGrn)
            msh.ver[i].grn = NmbGrn;
         else
         {
            CurGrn = msh.ver[i].grn;
            msh.ver[i].grn = ++NmbGrn;
         }
      }

      // Derive element's colors and grains from the vertices
      for(t=0;t<MAXELE;t++)
      {
         for(i=1;i<=msh.NmbEle[t];i++)
         {
            msh.ele[t][i].col = msh.ver[ msh.ele[t][i].idx[0] ].col;
            msh.ele[t][i].grn = msh.ver[ msh.ele[t][i].idx[0] ].grn;
         }
      }

      // Allocate a single colour table for the whole mesh as colours
      // are based on vertices only
      msh.ColPar = malloc( (msh.NmbCol + 1) * 3 * sizeof(int) );
      assert(msh.ColPar);

      // Each kind of entity needs a dedicated grain table to store
      // the begin and ending indices. Some grains may be empty
      msh.VerGrnPar = malloc( (msh.NmbGrn + 1) * 4 * sizeof(int) );
      assert(msh.VerGrnPar);
      msh.EleGrnPar = malloc( (msh.NmbGrn + 1) * MAXELE * 2 * sizeof(int) );
      assert(msh.EleGrnPar);

      // Setup vertex colours and grains partitions
      CurGrn = msh.ver[1].grn;
      NmbGrn = 1;
      msh.VerGrnPar[ NmbGrn ][0][0] = 1;
      msh.VerGrnPar[ NmbGrn ][0][2] = msh.ver[1].col;
      msh.VerGrnPar[ NmbGrn ][0][3] = msh.ver[1].grn;

      CurCol = msh.ver[1].col;
      NmbCol = 1;
      msh.ColPar[ NmbCol ][0] = 1;
      msh.ColPar[ NmbCol ][2] = NmbGrn;

      for(i=1;i<msh.NmbVer;i++)
      {
         if(msh.ver[i].grn != CurGrn)
         {
            msh.VerGrnPar[ NmbGrn ][0][1] = i - 1;
            NmbGrn++;
            msh.VerGrnPar[ NmbGrn ][0][0] = i;
            CurGrn = msh.ver[i].grn;
            msh.VerGrnPar[ NmbGrn ][0][2] = msh.ver[i].col;
            msh.VerGrnPar[ NmbGrn ][0][3] = msh.ver[i].grn;

            if(msh.ver[i].col != CurCol)
            {
               msh.ColPar[ NmbCol ][1] = NmbGrn - 1;
               NmbCol++;
               msh.ColPar[ NmbCol ][0] = NmbGrn;
               CurCol = msh.ver[i].col;
               msh.ColPar[ NmbCol ][2] = msh.ver[i].col;
            }
         }
      }

      msh.VerGrnPar[ NmbGrn ][0][1] = msh.NmbVer;
      msh.ColPar[ NmbCol ][1] = NmbGrn;

      for(i=1;i<=NmbGrn;i++)
         printf(  "vertex grain %3d (%3d/%3d): %8d -> %8d, size: %8d\n",
                  i, msh.VerGrnPar[i][0][2], msh.VerGrnPar[i][0][3],
                  msh.VerGrnPar[i][0][0], msh.VerGrnPar[i][0][1],
                  msh.VerGrnPar[i][0][1] - msh.VerGrnPar[i][0][0] + 1);

      for(i=1;i<=NmbCol;i++)
         printf(  "vertex color %3d (%3d): %8d -> %8d, size: %8d\n",
                  i, msh.ColPar[i][2], msh.ColPar[i][0], msh.ColPar[i][1],
                  msh.ColPar[i][1] - msh.ColPar[i][0] + 1);

      for(t=0;t<MAXELE;t++)
      {
         if(!msh.NmbEle[t])
            continue;

         for(i=1;i<=NmbGrn;i++)
            msh.EleGrnPar[i][t][0] = msh.EleGrnPar[i][t][1] = 0;

         printf("ele %d = %d items\n",t,msh.NmbEle[t]);
         for(i=1;i<=msh.NmbEle[t];i++)
         {
            grn = msh.ele[t][i].grn;

            if(!msh.EleGrnPar[ grn ][t][0])
            {
               msh.EleGrnPar[ grn ][t][0] = i;
               msh.EleGrnPar[ grn ][t][1] = i;
            }
            else
            {
               msh.EleGrnPar[ grn ][t][0] = MIN(msh.EleGrnPar[ grn ][t][0], i);
               msh.EleGrnPar[ grn ][t][1] = MAX(msh.EleGrnPar[ grn ][t][1], i);
            }
         }
         puts("done");

         for(i=1;i<=NmbGrn;i++)
            printf(  "%s grain %3d: %8d -> %8d, size: %8d\n",
                     EleNam[t], i,
                     msh.EleGrnPar[i][t][0], msh.EleGrnPar[i][t][1],
                     msh.EleGrnPar[i][t][1] - msh.EleGrnPar[i][t][0] + 1 );
      }

      ChkColGrn(&msh);
   }

   StopParallel(LibParIdx);

   return(0);
}


/*----------------------------------------------------------------------------*/
/* Compute the bounding box                                                   */
/*----------------------------------------------------------------------------*/

static void SetBndBox(int64_t BegIdx, int64_t EndIdx, MshSct *msh)
{
   int64_t  i;
   int      j;

   for(i=BegIdx; i<=EndIdx; i++)
   {
      msh->ver[i].idx = (int)i;

      if(msh->dim == 2)
         msh->ver[i].crd[2] = 0.;

      if(i==1)
         for(j=0;j<3;j++)
            msh->box[j] = msh->box[j+3] = msh->ver[i].crd[j];
      else
         for(j=0;j<3;j++)
            if(msh->ver[i].crd[j] < msh->box[j])
               msh->box[j] = msh->ver[i].crd[j];
            else if(msh->ver[i].crd[j] > msh->box[j+3])
               msh->box[j+3] = msh->ver[i].crd[j];
   }
}


/*----------------------------------------------------------------------------*/
/* Comparison of two items for the qsort                                      */
/*----------------------------------------------------------------------------*/

int CmpFnc(const void *a, const void *b)
{
   uint64_t *pa = (uint64_t *)a, *pb = (uint64_t *)b;

   if(*pa > *pb)
      return(1);
   else
      return(-1);
}


/*----------------------------------------------------------------------------*/
/* Compute the hilbert code from 3d coordinates                               */
/*----------------------------------------------------------------------------*/

uint64_t HilCod(double crd[3], double box[6], int itr, int mod)
{
   uint64_t IntCrd[3], m=1ULL<<63, cod;
   int      j, b, GeoWrd, NewWrd, BitTab[3] = {1,2,4};
   double   TmpCrd[3];
   int      rot[8], GeoCod[8]={0,3,7,4,1,2,6,5}; // Hilbert
   int      OctCod[8] = {5,4,7,6,1,0,3,2}; // octree or Z curve
   int      HilCod[8][8] = {
            {0,7,6,1,2,5,4,3}, {0,3,4,7,6,5,2,1},
            {0,3,4,7,6,5,2,1}, {2,3,0,1,6,7,4,5},
            {2,3,0,1,6,7,4,5}, {6,5,2,1,0,3,4,7},
            {6,5,2,1,0,3,4,7}, {4,3,2,5,6,1,0,7} };

   if(mod == RNDMOD)
      return(rand());

   // Convert double precision coordinates to integers
   for(j=0;j<3;j++)
   {
      TmpCrd[j] = (crd[j] - box[j]) * box[j+3];
      IntCrd[j] = (uint64_t)TmpCrd[j];
   }

   // Binary hilbert renumbering loop
   cod = 0;

   for(j=0;j<8;j++)
      rot[j] = GeoCod[j];

   for(b=0;b<itr;b++)
   {
      GeoWrd = 0;

      for(j=0;j<3;j++)
      {
         if(IntCrd[j] & m)
            GeoWrd |= BitTab[j];

         IntCrd[j] = IntCrd[j]<<1;
      }

      if(mod == OCTMOD)
      {
         NewWrd = OctCod[ GeoWrd ];
         cod = cod<<3 | NewWrd;
      }
      else
      {
         NewWrd = rot[ GeoWrd ];
         cod = cod<<3 | NewWrd;

         for(j=0;j<8;j++)
            rot[j] = HilCod[ NewWrd ][ rot[j] ];
      }
   }

   return(cod);
}


/*----------------------------------------------------------------------------*/
/* Parallel loop renumbering vertices                                         */
/*----------------------------------------------------------------------------*/

void RenVer(int BegIdx, int EndIdx, int PthIdx, MshSct *msh)
{
   int i;
   uint64_t ColCod = 0, GrnCod = 0, DegCod = 0, SchCod = 0;

   // Set the vertex code with the hilbert code, shifted by 11 bits in the case of
   // GMlib mode to make room for the high degree bit and the 10 bit ref hash tag
   for(i=BegIdx; i<=EndIdx; i++)
   {
      if(msh->ColGrnMod)
      {
         ColCod = (msh->ver[i].col & msh->ColMsk) << msh->ColLft;
         GrnCod = (msh->ver[i].grn & msh->GrnMsk) << msh->GrnLft;
      }

      if(msh->GmlMod)
         DegCod = (msh->ver[i].deg & msh->DegMsk) << msh->DegLft;

      if(msh->mod == IniMod)
         SchCod = i;
      else
         SchCod = HilCod(msh->ver[i].crd, msh->box, MAXITR, msh->mod);

      SchCod = SchCod >> msh->VerHilRgt;

      msh->ver[i].cod = ColCod | GrnCod | DegCod | SchCod;
   }
}


/*----------------------------------------------------------------------------*/
/* Compute the barycenter of any kind of element                              */
/*----------------------------------------------------------------------------*/

void SetMidCrd2(int NmbVer, int *IdxTab, MshSct *msh, double *crd)
{
   int         i, j, MaxEdg;
   const int   tvpe[6][2] = { {0,1}, {1,2}, {2,0}, {3,0}, {3,1}, {3,2} };
   double      *TetCrd[4], len, MinLen, MaxLen;

   // Default treatment is to compute the element's barycenter
   for(j=0;j<3;j++)
      crd[j] = 0.;

   for(i=0;i<NmbVer;i++)
      for(j=0;j<3;j++)
         crd[j] += msh->ver[ IdxTab[i] ].crd[j];

   for(j=0;j<3;j++)
      crd[j] /= NmbVer;

   // Special handling of tetrahedra: if the element is anisotropic,
   // use the longest edge midpoint instedad of the tet's barycenter
   // as it better represents an elongated element
   if(msh->TypIdx == TypTet)
   {
      MinLen = FLT_MAX;
      MaxEdg = -1;

      for(i=0;i<4;i++)
         TetCrd[i] = msh->ver[ IdxTab[i] ].crd;

      // Compute each edge length and search for the shortest and longest ones
      for(i=0;i<6;i++)
      {
         len = POW(TetCrd[ tvpe[i][0] ][0] - TetCrd[ tvpe[i][1] ][0])
             + POW(TetCrd[ tvpe[i][0] ][1] - TetCrd[ tvpe[i][1] ][1])
             + POW(TetCrd[ tvpe[i][0] ][2] - TetCrd[ tvpe[i][1] ][2]);

         MinLen = MIN(MinLen, len);

         if(MaxEdg == -1 || len > MaxLen)
         {
            MaxLen = len;
            MaxEdg = i;
         }
      }

      // If the longest edge is more than three times longer than the shortest one
      // the tet is anisotropic so its barycenter coordinates are replaced by
      // the longest edge's center point
      if(MaxEdg != -1 && MaxLen > POW(3) * MinLen)
         for(j=0;j<3;j++)
            crd[j] = (TetCrd[ tvpe[ MaxEdg ][0] ][j] + TetCrd[ tvpe[ MaxEdg ][1] ][j]) / 2.;
   }
}


/*----------------------------------------------------------------------------*/
/* Set the old to new index convertion table                                  */
/*----------------------------------------------------------------------------*/

void SetNewIdx(int NmbVer, int *IdxTab, int *Old2New)
{
   int i;

   for(i=0;i<NmbVer;i++)
      IdxTab[i] = Old2New[ IdxTab[i] ];
}


/*----------------------------------------------------------------------------*/
/* Parallel loop renumbering any kind of element                              */
/*----------------------------------------------------------------------------*/

void RenEle(int BegIdx, int EndIdx, int PthIdx, MshSct *msh)
{
   char *BytPtr;
   int      i, ref;
   double   crd[3];
   uint64_t ColCod = 0, GrnCod = 0, RefCod = 0, SchCod = 0;

   // Set the elements SFC code from their barycenter
   for(i=BegIdx; i<=EndIdx; i++)
   {
      SetNewIdx(EleTab[ msh->TypIdx ][0], msh->ele[ msh->TypIdx ][i].idx, msh->Old2New);

      if(msh->ColGrnMod)
      {
         ColCod = (msh->ele[ msh->TypIdx ][i].col & msh->ColMsk) << msh->ColLft;
         GrnCod = (msh->ele[ msh->TypIdx ][i].grn & msh->GrnMsk) << msh->GrnLft;
      }

      if(msh->GmlMod && ((msh->TypIdx == TypTri || (msh->TypIdx == TypQad))))
      {
         BytPtr = (char *)&msh->ele[ msh->TypIdx ][i].ref;
         ref = BytPtr[0] + BytPtr[1] + BytPtr[2] + BytPtr[3];
         RefCod = (ref & msh->RefMsk) << msh->RefLft;
      }
      else
         RefCod = 0;

      if(msh->mod == IniMod)
         SchCod = i;
      else
      {
         SetMidCrd2(EleTab[ msh->TypIdx ][0], msh->ele[ msh->TypIdx ][i].idx, msh, crd);
         SchCod = HilCod(crd, msh->box, MAXITR, msh->mod);
      }

      if((msh->TypIdx == TypTri) || (msh->TypIdx == TypQad))
         SchCod = SchCod >> msh->FacHilRgt;
      else if((msh->TypIdx >= TypTet) && (msh->TypIdx <= TypHex))
         SchCod = SchCod >> msh->VolHilRgt;
      else
         SchCod = 0;

      msh->ele[ msh->TypIdx ][i].cod = ColCod | GrnCod | RefCod | SchCod;
   }
}


/*----------------------------------------------------------------------------*/
/* Compute and print elements / vertices dependencies stats                   */
/*----------------------------------------------------------------------------*/

void PrtSta(MshSct *msh, int64_t LibParIdx)
{
   int      i, t;
   float    sta[2];

   for(t=0;t<MAXELE;t++)
      if(msh->EleTyp[t])
      {
         BeginDependency(LibParIdx, msh->EleTyp[t], msh->VerTyp);

         for(i=1;i<=msh->NmbEle[t];i++)
            AddDependencyFast(LibParIdx, 1, &i, EleTab[t][0], msh->ele[t][i].idx);

         EndDependency(LibParIdx, sta);
         printf(" %s : %3.2f%% / %3.2f%%\n", EleNam[t], sta[0],sta[1]);
      }

   puts("");
}


/*----------------------------------------------------------------------------*/
/* Allocate a new index buffer and copy the old indices to the new table      */
/*----------------------------------------------------------------------------*/

void SwpMem(MshSct *msh, int typ)
{
   int i, *IdxTab, siz = EleTab[ typ ][0];

   // Allocate the new table
   IdxTab = malloc( (size_t)(msh->NmbEle[ typ ] + 1) * siz * sizeof(int));
   assert(IdxTab);

   // Copy the old indices to the new table
   // and save the new pointer to each elements index address
   for(i=1;i<=msh->NmbEle[ typ ];i++)
   {
      memcpy(&IdxTab[ i * siz ], msh->ele[ typ ][i].idx, siz * sizeof(int));
      msh->ele[ typ ][i].idx = &IdxTab[ i * siz ];
   }

   // Free the old table and set the element type pointer with the new one
   free(msh->IdxTab[ typ ]);
   msh->IdxTab[ typ ] = IdxTab;
}


/*----------------------------------------------------------------------------*/
/* Set the code 10 upper bit with a hash key based on element's ref           */
/*----------------------------------------------------------------------------*/

void SetVerDeg(MshSct *msh)
{
   int i, j, t, (*DegTab)[ MAXELE ];
   uint64_t cod;
   VerSct *ver;
   EleSct *ele;

   // Allocate a vertex degree table with one scalar per kind of element
   DegTab = calloc( (size_t)(msh->NmbVer + 1), MAXELE * sizeof(int));
   assert(DegTab);

   // Add each element's vertices to the degree associated to its kind
   for(t=0;t<MAXELE;t++)
      for(i=1;i<=msh->NmbEle[t];i++)
         for(j=0;j<EleTab[t][0];j++)
            DegTab[ msh->ele[t][i].idx[j] ][t]++;

   // Set the high or low degree flag
   for(i=1;i<=msh->NmbVer;i++)
   {
      ver = &msh->ver[i];

      if( (DegTab[i][1] > MaxDeg[1][0])
      ||  (DegTab[i][2] > MaxDeg[2][0])
      ||  (DegTab[i][3] > MaxDeg[3][0])
      ||  (DegTab[i][6] > MaxDeg[6][0]) )
      {
         ver->deg = 1;
         msh->HghDeg++;

         for(j=0;j<MAXELE;j++)
            msh->MaxDeg[j] = MAX(msh->MaxDeg[j], DegTab[i][j]);
      }
      else
         ver->deg = 0;

      // Count the number of over connected vertices
      if( (DegTab[i][1] > MaxDeg[1][1])
      ||  (DegTab[i][2] > MaxDeg[2][1])
      ||  (DegTab[i][3] > MaxDeg[3][1])
      ||  (DegTab[i][6] > MaxDeg[6][1]) )
      {
         msh->OvrDeg++;
      }
   }

   // Set the right vector size for each kind of element ball
   for(j=0;j<MAXELE;j++)
      if(msh->MaxDeg[j])
         msh->DegVec[j] = pow(2., ceil(log2(msh->MaxDeg[j])));

   // Print statistics about vertex connectivity
   printf(  "High-connected vertices      : %3.6f%%\n",
            (100. * (float)msh->HghDeg) / (float)msh->NmbVer );

   printf(  "Over-connected vertices      : %3.6f%%\n",
            (100. * (float)msh->OvrDeg) / (float)msh->NmbVer );

   for(j=0;j<MAXELE;j++)
      if(msh->MaxDeg[j])
         printf(  "Ball of %s     : max deg = %3d, vec size = %3d\n",
                  EleNam[j], msh->MaxDeg[j],  msh->DegVec[j] );

   puts("");
}


/*----------------------------------------------------------------------------*/
/* Set the code 10 upper bit with a hash key based on element's ref           */
/*----------------------------------------------------------------------------*/

void SetMatSlc(MshSct *msh)
{
   char     *BytPtr;
   int      i, NmbEdg, *DegTab, VecCnt[6] = {0};
   uint64_t cod, DegTot = 0, VecTot = 0;
   VerSct   *ver;
   EleSct   *ele;

   // Allocate a degree table with one scalar per kind of element
   DegTab = calloc( (size_t)(msh->NmbVer + 1), sizeof(int));
   assert(DegTab);

   NmbEdg = SetDeg(msh, DegTab);
   printf("Unique edges extracted       : %d\n", NmbEdg);

   // Set the vertex degree with its slice vector size code: 1 -> 5
   for(i=1;i<=msh->NmbVer;i++)
   {
      ver = &msh->ver[i];

      if(DegTab[i] <= 16)
         cod = 1;
      else if(DegTab[i] <= 32)
         cod = 2;
      else if(DegTab[i] <= 64)
         cod = 3;
      else if(DegTab[i] <= 128)
         cod = 4;
      else
         cod = 5;

      ver->deg = cod;
      VecCnt[ cod ]++;
      DegTot += DegTab[i];
   }

   VecTot = 16 * VecCnt[1] + 32 * VecCnt[2] + 64 * VecCnt[3]
         + 128 * VecCnt[4] + 256 * VecCnt[5];

   puts("");
   puts  ("vector |  %age  | number");
   puts  ("----------------------------");
   printf("   16  | %6.2f | %10d\n", (float)(100 * VecCnt[1]) / msh->NmbVer, VecCnt[1]);
   printf("   32  | %6.2f | %10d\n", (float)(100 * VecCnt[2]) / msh->NmbVer, VecCnt[2]);
   printf("   64  | %6.2f | %10d\n", (float)(100 * VecCnt[3]) / msh->NmbVer, VecCnt[3]);
   printf("  128  | %6.2f | %10d\n", (float)(100 * VecCnt[4]) / msh->NmbVer, VecCnt[4]);
   printf("  256  | %6.2f | %10d\n", (float)(100 * VecCnt[5]) / msh->NmbVer, VecCnt[5]);

   puts("");
   printf("vector filling : %3.2f%%\n", (float)(100 * DegTot) / VecTot);
   printf("real non-zero  : %lld\n", DegTot);
   printf("vector non-zero: %lld\n", VecTot);
   puts("");
}


/*----------------------------------------------------------------------------*/
/* Build the list of unique edges sequentialy                                 */
/*----------------------------------------------------------------------------*/

int SetDeg(MshSct *msh, int *DegTab)
{
   int i, j, idx0, idx1, key, MinIdx, MaxIdx, siz, col, *tet, NmbEdg = 0;
   HshSct *hsh, *buc;

   // Allocate the hash table to store all edges
   col = siz = msh->NmbEle[ TypTet ];
   hsh = malloc( 7LL * (size_t)siz * sizeof(HshSct));
   assert(hsh);

   // Clear the hash table's direct entries,
   // there is no need to clear the collision entries
   memset(hsh, 0, siz * sizeof(HshSct));

   // Loop over each tet and each tet's edges
   for(i=1;i<=siz;i++)
   {
      tet = msh->ele[ TypTet ][i].idx;

      for(j=0;j<6;j++)
      {
         // Compute the hashing key from the edge's vertex indices
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
            NmbEdg++;
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
               hsh[ key ].NexBuc = 0;
               NmbEdg++;
               break;
            }
         }while(1);
      }
   }

   for(i=0;i<siz;i++)
   {
      buc = &hsh[i];

      if(buc->MinIdx)
      {
         DegTab[ buc->MinIdx ]++;
         DegTab[ buc->MaxIdx ]++;
      }
   }

   buc = &hsh[ siz ];

   while(buc->MinIdx)
   {
      DegTab[ buc->MinIdx ]++;
      DegTab[ buc->MaxIdx ]++;
      buc++;
   }

   free(hsh);

   return(NmbEdg);
}


/*----------------------------------------------------------------------------*/
/* Look for possible vertex pinch between different grains from the same color*/
/*----------------------------------------------------------------------------*/

void ChkColGrn(MshSct *msh)
{
   int i, j, k, VerIdx, NmbCol = 0;
   int *VerColTab = calloc(msh->NmbVer + 1, sizeof(int));
   int *VerGrnTab = calloc(msh->NmbVer + 1, sizeof(int));

   printf("Check color-grain consistency: ");

   for(i=1;i<=msh->NmbCol;i++)
   {
      for(j=1;j<=msh->NmbEle[ TypTet ];j++)
      {
         if(msh->ele[ TypTet ][j].col != i)
            continue;

         for(k=0;k<4;k++)
         {
            VerIdx = msh->ele[ TypTet ][j].idx[k];

            if(i > VerColTab[ VerIdx ])
            {
               VerColTab[ VerIdx ] = i;
               VerGrnTab[ VerIdx ] = msh->ele[ TypTet ][j].grn;
            }
            else if( (VerColTab[ VerIdx ] == i)
                  && (VerGrnTab[ VerIdx ] != msh->ele[ TypTet ][j].grn) )
            {
               printf("vertex %d / tet %d / grain %d / color %d, collide with grain %d\n",
                        VerIdx, j, msh->ele[ TypTet ][j].grn, i, VerGrnTab[ VerIdx ]);
               NmbCol++;
            }
         }
      }
   }

   if(!NmbCol)
      puts("OK");
   else
      printf("%d collisions\n", NmbCol);

   free(VerColTab);
   free(VerGrnTab);
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
