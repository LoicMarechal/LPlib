

/*----------------------------------------------------------------------------*/
/*                                                                            */
/*                               LPlib V4.13                                  */
/*                                                                            */
/*----------------------------------------------------------------------------*/
/*                                                                            */
/*   Description:       Handles threads, scheduling & dependencies            */
/*   Author:            Loic MARECHAL                                         */
/*   Creation date:     feb 25 2008                                           */
/*   Last modification: oct 13 2025                                           */
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
#include <stdarg.h>
#include <assert.h>
#include <float.h>
#include <math.h>
#include <errno.h>
#include <limits.h>
#include "lplib4.h"

#ifdef _WIN32
#include <windows.h>
#include <sys/timeb.h>
#else
#include <unistd.h>
#include <sys/time.h>
#endif

#if defined(_WIN32) && !defined(NATIVE_WINDOWS)
#include "winpthreads.h"
#else
#include <pthread.h>
#endif

#ifdef WITH_LIBMEMBLOCKS
#include <libmemblocks1.h>
#endif

#ifdef WITH_METIS
#include <metis.h>
#endif


/*----------------------------------------------------------------------------*/
/* Defintion of macro commands and constants                                  */
/*----------------------------------------------------------------------------*/

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define POW(a)    ((a) * (a))

#define MaxLibPar 10
#define MaxTyp    100
#define DefSmlBlk 64
#define DefDepBlk 256
#define MaxTotPip 65536
#define MaxPipDep 100
#define MaxHsh    10
#define HshBit    16
#define MaxVarArg 20
#define MaxF77Arg 20
#define WrkPerGrp 8
#define BigMemSiz 100000000ULL
#define MAXEDG    1000
#define MAXITR    21
#define HILMOD    0
#define OCTMOD    1
#define RNDMOD    2

enum {HilMod=0, OctMod, RndMod, IniMod, TopMod};
enum ParCmd {  RunBigWrk, RunSmlWrk, RunDetWrk, RunColWrk,
               ClrMem, CpyMem, RunGrnWrk, EndPth };


/*----------------------------------------------------------------------------*/
/* Structures' prototypes                                                     */
/*----------------------------------------------------------------------------*/

typedef struct WrkSct
{
   itg               BegIdx, EndIdx, ItlTab[ MaxPth ][2];
   int               NmbDep, *DepWrdTab, GrpIdx, rnd, GrnIdx;
   double            RunTim;
   struct WrkSct     *pre, *nex;
}WrkSct;

typedef struct GrpSct
{
   int               idx, NmbSmlWrk[ MaxPth ];
   WrkSct            *SmlWrkTab[ MaxPth ][ WrkPerGrp ];
   struct GrpSct     *nex;
}GrpSct;

typedef struct
{
   itg               NmbLin, MaxNmbLin;
   int               NmbSmlWrk, SmlWrkSiz, DepWrkSiz, NmbGrp;
   int               NmbDepWrd, *DepWrdMat, *RunDepTab;
   WrkSct            *SmlWrkTab, *BigWrkTab;
   GrpSct            *NexGrp;
}TypSct;

typedef struct
{
   int               idx, NmbDetWrk;
   char              *ClrAdr, *DstAdr, *SrcAdr;
   size_t            StkSiz, CpyMemSiz, ClrMemSiz;
   void *            *UsrStk;
   WrkSct            *wrk, **DetWrkTab;
   pthread_mutex_t   mtx;
   pthread_cond_t    cnd;
   pthread_t         pth;
   pthread_attr_t    atr;
   struct ParSct     *par;
}PthSct;

typedef struct PipSct
{
   int               idx, NmbVarArg, NmbDep, DepTab[ MaxPipDep ], BegIdx, EndIdx, GrnIdx;
   void              *prc, *arg, *VarArgTab[ MaxVarArg ];
   size_t            StkSiz;
   void              *UsrStk;
   pthread_attr_t    atr;
   pthread_t         pth;
   struct ParSct     *par;
}PipSct;

typedef struct
{
   int               BegIdx, EndIdx, GrnIdx;
   void              *prc, *arg;
   pthread_t         pth;
   struct ParSct     *par;
}GrnSct;

typedef struct ParSct
{
   int               NmbCpu, WrkCpt, NmbPip, PenPip, RunPip, NmbTyp, DynSch;
   int               req, cmd, *PipWrd, SizMul, NmbVarArg, NmbDep;
   int               WrkSizSrt, NmbItlBlk, ItlBlkSiz, BufMax, BufCpt;
   int               NmbSmlBlk, NmbDepBlk, NmbColGrn, GrnNxt, GrnDon, clk;
   int               NmbGrnWrk, GrnWrkSiz, DepGrnSiz, NmbCol, NmbGrn, *GrnMat;
   int               (*GrnTab[ LplMax ])[2], (*ColTab)[2], CurCol;
   int               NmbDepWrd, *RunDepTab, *ColCpt, *GrnCol;
   int               NmbGrnWrd, *GrnWrdMat, *RunGrnTab, TypIdx[ LplMax ];
   size_t            StkSiz;
   void              *lmb, *VarArgTab[ MaxVarArg ];
   float             sta[2];
   void              (*prc)(itg, itg, int, void *), *arg;
   pthread_cond_t    ParCnd, PipCnd;
   pthread_mutex_t   ParMtx, PipMtx;
   pthread_t         PipPth;
   PthSct            *PthTab;
   TypSct            *TypTab, *CurTyp, *DepTyp, *typ1, *typ2;
   WrkSct            *NexWrk, *BufWrk[ MaxPth / 4 ], *GrnWrkTab;
}ParSct;

typedef struct
{
   uint64_t          (*idx)[2];
   double            box[6], (*crd)[3], (*crd2)[2];
}ArgSct;

typedef struct
{
   void              *base;
   size_t            nel, width;
   int               (*compar)(const void *, const void *);
}PipArgSct;

typedef struct
{
   itg MinIdx, MaxIdx, NexBuc;
}HshSct;

#ifdef WITH_METIS
typedef struct
{
   idx_t             nvtxs, nedges, ncon, nparts, *xadj, *adjncy, *VerDeg, objval, *part;
   idx_t            *adjwgt, options[ METIS_NOPTIONS ];
}MtsSct;

typedef itg          int1d;
typedef itg          int2d[2];
typedef itg          int11d[11];
typedef itg          int56d[56];

#endif


/*----------------------------------------------------------------------------*/
/* Global tables                                                              */
/*----------------------------------------------------------------------------*/

// Tetrahderon's table of edges
static const int tvpe[6][2] = { {0,1}, {1,2}, {2,0}, {3,0}, {3,1}, {3,2} };

// number of nodes per element kind
static const int EleSiz[ LplMax ] = {0,2,3,4,4,5,6,8,3,6,9,10,14,18,27};

// For each kind of element: GMF keyword strings
static const char *EleNam[ LplMax ] = {
   "Vertices        ",
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
static const int MaxDeg[ LplMax ][2] = {
   { 0,   0},
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
/* Private procedures' prototypes                                             */
/*----------------------------------------------------------------------------*/

static int        SetBit         (int *, int);
static int        GetBit         (int *, int);
static void       ClrBit         (int *, int);
static int        AndWrd         (int, int *, int *);
static void       AddWrd         (int, int *, int *);
static void       SubWrd         (int, int *, int *);
static void       ClrWrd         (int, int *);
static void       CpyWrd         (int, int *, int *);
int               CmpWrk         (const void *, const void *);
static void      *PipHdl         (void *);
static void      *PthHdl         (void *);
static WrkSct    *NexWrk         (ParSct *, int);
static WrkSct    *NexGrn         (ParSct *, int);
void              PipSrt         (PipArgSct *);
static void       CalVarArgPip   (PipSct *, void *);
static void       CalVarArgPrc   (itg, itg, int, ParSct *);
static int64_t    IniPar         (int, size_t, void *);
static void       SetItlBlk      (ParSct *, TypSct *);
static int        SetGrp         (ParSct *, TypSct *);
static void      *LPL_malloc     (void *, int64_t);
static void      *LPL_calloc     (void *, int64_t, int64_t);
static void       LPL_free       (void *, void *);
static void       UpdBlkSiz      (ParSct *, TypSct *);
static void       SetBndBox      (LplSct *);
static void       SetMidCrd      (int , int *, LplSct *, double *);
static int        CmpFnc         (const void *, const void *);
static void       RenVer         (int, int, int, LplSct *);
static void       RenEle         (int, int, int, LplSct *);
static void       SetVerDeg      (LplSct *);
static void       SetMatSlc      (LplSct *);
static int        SetDeg         (LplSct *, int *);
static uint64_t   GetHilCod      (double *, double *, int, int);

#ifdef WITH_METIS
static void       ChkColGrn      (LplSct *);
static int        SetColors      (int64_t, int, int (*)[2], int);
static int        SetGrains      (int64_t, int, int, int (*)[2]);
static void       SetGrnDep      (ParSct *, LplSct *);
static void       SetV2VBal      (LplSct *);
static void       SetRk2Bal      (LplSct *);
static void       BuildMetisGraph   (LplSct *, MtsSct *);
static int        MetisPartitioning (LplSct *, int);
static int       *ColorPartImplicit (int, int, int, uint64_t *, int *,
                                     uint64_t *, int *,int *, int);
#endif


/*----------------------------------------------------------------------------*/
/* Init structures, scheduler and launch threads                              */
/*----------------------------------------------------------------------------*/

int64_t InitParallel(int NmbCpu)
{
   return(IniPar(NmbCpu, 0, NULL));
}

int64_t InitParallelAttr(int NmbCpu, size_t StkSiz, void *lmb)
{
#ifdef PTHREAD_STACK_MIN
   if(StkSiz < PTHREAD_STACK_MIN)
      StkSiz = PTHREAD_STACK_MIN;
#endif
   return(IniPar(NmbCpu, StkSiz, lmb));
}

static int64_t IniPar(int NmbCpu, size_t StkSiz, void *lmb)
{
   int i;
   int64_t ParIdx;
   ParSct *par;
   PthSct *pth;

   // Check the number of requested cpus
   if(NmbCpu < 1)
   {
      NmbCpu = GetNumberOfCores();

      if(NmbCpu < 1)
         return(0);
   }

   if(NmbCpu > MaxPth)
      NmbCpu = MaxPth;

   // Allocate and build main parallel structure
   if(!(par = LPL_calloc(lmb, 1, sizeof(ParSct))))
      return(0);

   // Pass along a potential libMemBlocks structure
   par->lmb = lmb;

   if(!(par->PthTab = LPL_calloc(par->lmb, NmbCpu, sizeof(PthSct))))
      return(0);

   if(!(par->TypTab = LPL_calloc(par->lmb, (MaxTyp + 1), sizeof(TypSct))))
      return(0);

   if(!(par->PipWrd = LPL_calloc(par->lmb, MaxTotPip/32, sizeof(int))))
      return(0);

   par->NmbCpu = NmbCpu;
   par->WrkCpt = par->NmbPip = par->PenPip = par->RunPip = 0;
   par->SizMul = 2;
   par->StkSiz = StkSiz;
   par->NmbItlBlk = 1;
   par->WrkSizSrt = 1;
   par->DynSch = 1;
   par->NmbSmlBlk = DefSmlBlk;
   par->NmbDepBlk = DefDepBlk;

   // Set the size of WP buffer
   if(NmbCpu >= 4)
      par->BufMax = NmbCpu / 4;
   else
      par->BufMax = 1;

   pthread_mutex_init(&par->ParMtx, NULL);
   pthread_mutex_init(&par->PipMtx, NULL);
   pthread_cond_init(&par->ParCnd, NULL);
   pthread_cond_init(&par->PipCnd, NULL);

   // Launch pthreads
   for(i=0;i<par->NmbCpu;i++)
   {
      pth = &par->PthTab[i];
      pth->idx = i;
      pth->par = par;
      pthread_mutex_init(&pth->mtx, NULL);
      pthread_cond_init(&pth->cnd, NULL);

      if(StkSiz)
      {
         pthread_attr_init(&pth->atr);
         pth->StkSiz = StkSiz;
         pth->UsrStk = LPL_malloc(par->lmb, pth->StkSiz);
#ifdef _WIN32
         pthread_attr_setstackaddr(&pth->atr, pth->UsrStk);
         pthread_attr_setstacksize(&pth->atr, pth->StkSiz);
#else
         pthread_attr_setstack(&pth->atr, pth->UsrStk, pth->StkSiz);
#endif
         pthread_create(&pth->pth, &pth->atr, PthHdl, (void *)pth);
      }
      else
      {
         pth->StkSiz = 0;
         pth->UsrStk = NULL;
         pthread_create(&pth->pth, NULL, PthHdl, (void *)pth);
      }
   }

   // Wait for all threads to be up and wainting
   pthread_mutex_lock(&par->ParMtx);

   while(par->WrkCpt < par->NmbCpu)
      pthread_cond_wait(&par->ParCnd, &par->ParMtx);

   pthread_mutex_unlock(&par->ParMtx);

   ParIdx = (int64_t)par;

   return(ParIdx);
}


/*----------------------------------------------------------------------------*/
/* Stop all threads and free memories                                         */
/*----------------------------------------------------------------------------*/

void StopParallel(int64_t ParIdx)
{
   int i;
   PthSct *pth;
   ParSct *par = (ParSct *)ParIdx;

   // Get and check lib parallel instance
   if(!ParIdx)
      return;

   // Send stop to all threads
   pthread_mutex_lock(&par->ParMtx);
   par->cmd = EndPth;
   pthread_mutex_unlock(&par->ParMtx);

   // Wait for all threads to complete
   for(i=0;i<par->NmbCpu;i++)
   {
      pth = &par->PthTab[i];
      pthread_mutex_lock(&pth->mtx);
      pthread_cond_signal(&pth->cnd);
      pthread_mutex_unlock(&pth->mtx);
      pthread_join(pth->pth, NULL);

      if(pth->UsrStk)
         LPL_free(par->lmb, pth->UsrStk);
   }

   pthread_mutex_destroy(&par->ParMtx);
   pthread_cond_destroy(&par->ParCnd);

   WaitPipeline(ParIdx);

   pthread_mutex_destroy(&par->PipMtx);
   pthread_cond_destroy(&par->PipCnd);

   // Free memories
   for(i=1;i<=MaxTyp;i++)
      if(par->TypTab[i].NmbLin)
         FreeType(ParIdx, i);

   LPL_free(par->lmb, par->PthTab);
   LPL_free(par->lmb, par->TypTab);
   LPL_free(par->lmb, par->PipWrd);
   free(par);
}


/*----------------------------------------------------------------------------*/
/* Try to detect the number of system's threads                               */
/*----------------------------------------------------------------------------*/

int GetNumberOfCores()
{
#ifdef _WIN32
   SYSTEM_INFO info;
   GetSystemInfo(&info);
   return(info.dwNumberOfProcessors);
#else
#ifdef _SC_NPROCESSORS_ONLN
   return((itg)sysconf(_SC_NPROCESSORS_ONLN));
#else
   return(1);
#endif
#endif
}


/*----------------------------------------------------------------------------*/
/* Set extra attributes to fine tune the blocks processing order              */
/*----------------------------------------------------------------------------*/

int SetExtendedAttributes(int64_t ParIdx, ...)
{
   int NmbArg = 0, ArgCod, ArgVal;
   ParSct *par = (ParSct *)ParIdx;
   va_list ArgLst;

   // Attributes cannot be modified the a prallel loop is running
   if(par->typ1)
      return(0);

   // Read the argument list that must be terminated by a zero
   va_start(ArgLst, ParIdx);

   ArgCod = va_arg(ArgLst, int);

   switch(ArgCod)
   {
      // Set the number of interleaved blocks in independant loops
      case SetInterleavingFactor :
      {
         ArgVal = va_arg(ArgLst, int);

         if(ArgVal > 0)
         {
            par->NmbItlBlk = ArgVal;
            par->ItlBlkSiz = 0;
            NmbArg++;
         }
      }break;

      // Set the interleaved blocks size in independant loops
      case SetInterleavingSize :
      {
         ArgVal = va_arg(ArgLst, int);

         if(ArgVal > 0)
         {
            par->NmbItlBlk = 0;
            par->ItlBlkSiz = ArgVal;
            NmbArg++;
         }
      }break;

      // Desable blocks interleaving (default)
      case DisableInterleaving :
      {
         par->NmbItlBlk = 1;
         par->ItlBlkSiz = 0;
         NmbArg++;
      }break;

      // Sort depedancy-loop WP through their number of dependencies (default)
      case EnableBlockSorting :
      {
         par->WrkSizSrt = 1;
         NmbArg++;
      }break;

      // Disable WP sorting: it lowers concurrency but enhances cache reuse
      case DisableBlockSorting :
      {
         par->WrkSizSrt = 0;
         NmbArg++;
      }break;

      // Static scheduling makes the library deterministic
      case StaticScheduling :
      {
         // WP sorting is useless in this mode so it is disabled
         par->WrkSizSrt = par->DynSch = 0;
         NmbArg++;
      }break;

      case SetSmallBlock :
      {
         ArgVal = va_arg(ArgLst, int);

         if(ArgVal > 0)
         {
            par->NmbSmlBlk = ArgVal;
            NmbArg++;
         }
      }break;

      case SetDependencyBlock :
      {
         ArgVal = va_arg(ArgLst, int);

         if(ArgVal > 0)
         {
            par->NmbDepBlk = ArgVal;
            NmbArg++;
         }
      }break;

      case EnableAdaptiveSizing :
      {
         par->clk = 1;
         NmbArg++;
      }break;

      case DisableAdaptiveSizing :
      {
         par->clk = 0;
         NmbArg++;
      }break;
   }

   va_end(ArgLst);

   return(NmbArg);
}


/*----------------------------------------------------------------------------*/
/* Launch the loop prc on typ1 element depending on typ2                      */
/*----------------------------------------------------------------------------*/

float LaunchParallel(int64_t ParIdx, int TypIdx1, int TypIdx2,
                     void *prc, void *PtrArg )
{
   int      i;
   float    acc = 0.;
   PthSct   *pth;
   ParSct   *par = (ParSct *)ParIdx;
   TypSct   *typ1, *typ2 = NULL;
   GrpSct   *grp;

   // Get and check lib parallel instance
   if(!ParIdx)
      return(-1.);

   // Check bounds
   if( (TypIdx1 < 1) || (TypIdx1 > MaxTyp)
   ||  (TypIdx2 > MaxTyp) || (TypIdx1 == TypIdx2) )
   {
      return(-1.);
   }

   typ1 =  &par->TypTab[ TypIdx1 ];

   // Launch small WP with static scheduling
   if( (TypIdx2 > 0) && !par->DynSch )
   {
      grp = typ1->NexGrp;

      do
      {
         // Lock acces to global parameters
         pthread_mutex_lock(&par->ParMtx);

         par->cmd = RunDetWrk;
         par->prc = (void (*)(itg, itg, int, void *))prc;
         par->arg = PtrArg;
         par->typ1 = typ1;
         par->typ2 = NULL;
         par->WrkCpt = 0;

         for(i=0;i<par->NmbCpu;i++)
         {
            pth = &par->PthTab[i];
            pth->NmbDetWrk = grp->NmbSmlWrk[i];
            pth->DetWrkTab = grp->SmlWrkTab[i];
            acc += (float)grp->NmbSmlWrk[i];
         }

         for(i=0;i<par->NmbCpu;i++)
         {
            pth = &par->PthTab[i];
            pthread_mutex_lock(&pth->mtx);
            pthread_cond_signal(&pth->cnd);
            pthread_mutex_unlock(&pth->mtx);
         }

         pthread_cond_wait(&par->ParCnd, &par->ParMtx);

         pthread_mutex_unlock(&par->ParMtx);
         grp = grp->nex;
      }while(grp);

      acc /= (float)(par->NmbSmlBlk * typ1->NmbGrp) / (float)WrkPerGrp;
   }
   else if( (TypIdx2 > 0) && par->DynSch )
   {
      // Launch small WP with dynamic scheduling

      // Lock acces to global parameters
      pthread_mutex_lock(&par->ParMtx);

      par->cmd = RunSmlWrk;
      par->prc = (void (*)(itg, itg, int, void *))prc;
      par->arg = PtrArg;
      par->typ1 = typ1;
      par->typ2 = typ2 = &par->TypTab[ TypIdx2 ];
      par->NexWrk = typ1->SmlWrkTab;
      par->BufCpt = 0;
      par->WrkCpt = 0;
      par->sta[0] = par->sta[1] = 0.;
      par->req = 0;
      par->NmbDep = 0;

      // Clear running wp
      for(i=0;i<par->NmbCpu;i++)
         par->PthTab[i].wrk = NULL;

      ClrWrd(typ1->NmbDepWrd, typ1->RunDepTab);

      // Build a linked list of wp
      for(i=0;i<par->typ1->NmbSmlWrk;i++)
      {
         typ1->SmlWrkTab[i].pre = &typ1->SmlWrkTab[ i-1 ];
         typ1->SmlWrkTab[i].nex = &typ1->SmlWrkTab[ i+1 ];
      }

      typ1->SmlWrkTab[0].pre = typ1->SmlWrkTab[ typ1->NmbSmlWrk - 1 ].nex = NULL;

      // Main loop: wake up threads and wait for completion or blocked threads
      do
      {
         // Search for some idle threads
         par->req = 0;

         for(i=0;i<par->NmbCpu;i++)
         {
            pth = &par->PthTab[i];

            if(pth->wrk)
               continue;

            if(!(pth->wrk = NexWrk(par, i)))
            {
               par->req = 1;
               break;
            }

            // Wake up the thread and provide it with a WP list
            pthread_mutex_lock(&pth->mtx);
            pthread_cond_signal(&pth->cnd);
            pthread_mutex_unlock(&pth->mtx);
         }

         // If every WP are done : exit the parallel loop
         if(par->WrkCpt == typ1->NmbSmlWrk)
            break;

         // Otherwise, wait for a blocked thread
         pthread_cond_wait(&par->ParCnd, &par->ParMtx);
      }while(1);

      pthread_mutex_unlock(&par->ParMtx);

      // Compute the average concurrency factor
      acc = par->sta[0] ? (par->sta[1] / par->sta[0]) : 0;
   }
   else if(!TypIdx2)
   {
      // Launch big WP with static scheduling

      // Lock acces to global parameters
      pthread_mutex_lock(&par->ParMtx);

      par->cmd = RunBigWrk;
      par->prc = (void (*)(itg, itg, int, void *))prc;
      par->arg = PtrArg;
      par->typ1 = typ1;
      par->typ2 = NULL;
      par->WrkCpt = 0;

      for(i=0;i<par->NmbCpu;i++)
      {
         pth = &par->PthTab[i];
         pth->wrk = &typ1->BigWrkTab[i];
      }

      // Update block interleaving according to the current attributes
      if( (par->NmbItlBlk != 1) || par->ItlBlkSiz)
         SetItlBlk(par, typ1);

      for(i=0;i<par->NmbCpu;i++)
      {
         pth = &par->PthTab[i];
         pthread_mutex_lock(&pth->mtx);
         pthread_cond_signal(&pth->cnd);
         pthread_mutex_unlock(&pth->mtx);
      }

      pthread_cond_wait(&par->ParCnd, &par->ParMtx);

      pthread_mutex_unlock(&par->ParMtx);

      // Arbitrary set the average concurrency factor
      acc = (float)par->NmbCpu;

      if(par->clk)
         UpdBlkSiz(par, typ1);
   }
   else
      return(-1.);

   // Clear the main datatyp loop to indicate that no LaunchParallel is running
   par->typ1 = 0;

   // Return the concurrency factor
   return(acc);
}


/*----------------------------------------------------------------------------*/
/* Launch a parallel procudure with variable arguments.                       */
/* Arguments are passed as pointer to void.                                   */
/*----------------------------------------------------------------------------*/

float LaunchParallelMultiArg( int64_t ParIdx, int TypIdx1, int TypIdx2,
                              void *prc, int NmbArg, ... )
{
   int i;
   float acc;
   va_list ArgLst;
   ParSct *par = (ParSct *)ParIdx;

   if(NmbArg > 20)
      return(-1.);

   par->NmbVarArg = NmbArg;
   va_start(ArgLst, NmbArg);

   for(i=0;i<NmbArg;i++)
      par->VarArgTab[i] = va_arg(ArgLst, void *);

   va_end(ArgLst);

   acc = LaunchParallel(ParIdx, TypIdx1, TypIdx2, prc, NULL);
   par->NmbVarArg = 0;

   return(acc);
}


/*----------------------------------------------------------------------------*/
/* Pthread handler, waits for job, does it, then signal end                   */
/*----------------------------------------------------------------------------*/

static void *PthHdl(void *ptr)
{
   itg i, beg, end;
   PthSct *pth = (PthSct *)ptr;
   ParSct *par = pth->par;

   // Tell the scheduler if all threads are ready
   pthread_mutex_lock(&par->ParMtx);
   par->WrkCpt++;
   pthread_cond_signal(&par->ParCnd);
   pthread_mutex_lock(&pth->mtx);
   pthread_mutex_unlock(&par->ParMtx);

   // Enter main loop until StopParallel is send
   do
   {
      // Wait for a wake-up signal from the main loop
      pthread_cond_wait(&pth->cnd, &pth->mtx);

      // Update stats
      par->sta[0]++;

      for(i=0;i<par->NmbCpu;i++)
         if(par->PthTab[i].wrk)
            par->sta[1]++;

      switch(par->cmd)
      {
         // Call user's procedure with big WP
         case RunBigWrk :
         {
            // Loop over the interleaved blocks
            for(i=0;i<par->NmbItlBlk;i++)
            {
               beg = pth->wrk->ItlTab[i][0];
               end = pth->wrk->ItlTab[i][1];

               if(!beg || !end || (end < beg))
                  continue;

               if(par->clk)
                  pth->wrk->RunTim = GetWallClock();

               // Launch a single big wp and signal completion to the scheduler
               if(par->NmbVarArg)
                  CalVarArgPrc(beg, end, pth->idx, par);
               else
                  par->prc(beg, end, pth->idx, par->arg);

               if(par->clk)
                  pth->wrk->RunTim = GetWallClock() - pth->wrk->RunTim;
            }

            pthread_mutex_lock(&par->ParMtx);
            par->WrkCpt++;

            if(par->WrkCpt >= par->NmbCpu)
               pthread_cond_signal(&par->ParCnd);

            pthread_mutex_unlock(&par->ParMtx);
         }break;

         // Call user's procedure with small WP using dynamic scheduling
         case RunSmlWrk :
         {
            do
            {
               // Run the WP
               if(par->NmbVarArg)
                  CalVarArgPrc(pth->wrk->BegIdx, pth->wrk->EndIdx, pth->idx, par);
               else
                  par->prc(pth->wrk->BegIdx, pth->wrk->EndIdx, pth->idx, par->arg);

               // Locked acces to global parameters: 
               // update WP count, tag WP done and signal the main loop
               pthread_mutex_lock(&par->ParMtx);

               par->WrkCpt++;

               if(!(pth->wrk = NexWrk(par, pth->idx)))
               {
                  par->req = 1;
                  pthread_cond_signal(&par->ParCnd);
                  pthread_mutex_unlock(&par->ParMtx);
                  break;
               }

               if(par->req)
                  pthread_cond_signal(&par->ParCnd);

               pthread_mutex_unlock(&par->ParMtx);
            }while(1);
         }break;

         // Call user's procedure with small WP using dynamic scheduling
         case RunGrnWrk :
         {
            do
            {
               // Run the WP
               if(par->NmbVarArg)
                  CalVarArgPrc(pth->wrk->BegIdx, pth->wrk->EndIdx, pth->wrk->GrnIdx, par);
               else
                  par->prc(pth->wrk->BegIdx, pth->wrk->EndIdx, pth->wrk->GrnIdx, par->arg);

               // Locked acces to global parameters: 
               // update WP count, tag WP done and signal the main loop
               pthread_mutex_lock(&par->ParMtx);

               par->WrkCpt++;

               if(!(pth->wrk = NexGrn(par, pth->idx)))
               {
                  par->req = 1;
                  pthread_cond_signal(&par->ParCnd);
                  pthread_mutex_unlock(&par->ParMtx);
                  break;
               }

               if(par->req)
                  pthread_cond_signal(&par->ParCnd);

               pthread_mutex_unlock(&par->ParMtx);
            }while(1);
         }break;

         // Call user's procedure with small WP using static scheduling
         case RunDetWrk :
         {
            // Loop over the groups' WP
            for(i=0;i<pth->NmbDetWrk;i++)
            {
               beg = pth->DetWrkTab[i]->BegIdx;
               end = pth->DetWrkTab[i]->EndIdx;

               if(par->NmbVarArg)
                  CalVarArgPrc(beg, end, pth->idx, par);
               else
                  par->prc(beg, end, pth->idx, par->arg);
            }

            pthread_mutex_lock(&par->ParMtx);
            par->WrkCpt++;

            if(par->WrkCpt >= par->NmbCpu)
               pthread_cond_signal(&par->ParCnd);

            pthread_mutex_unlock(&par->ParMtx);
         }break;

         case ClrMem :
         {
            // Clear memory and signal completion to the scheduler
            memset(pth->ClrAdr, 0, pth->ClrMemSiz);

            pthread_mutex_lock(&par->ParMtx);
            par->WrkCpt++;
            pthread_cond_signal(&par->ParCnd);
            pthread_mutex_unlock(&par->ParMtx);
         }break;

         case CpyMem :
         {
            // Copy memory and signal completion to the scheduler
            memcpy(pth->DstAdr, pth->SrcAdr, pth->CpyMemSiz);

            pthread_mutex_lock(&par->ParMtx);
            par->WrkCpt++;
            pthread_cond_signal(&par->ParCnd);
            pthread_mutex_unlock(&par->ParMtx);
         }break;

         case EndPth :
         {
            // Destroy the thread mutex and condition and call for join
            pthread_mutex_unlock(&pth->mtx);
            pthread_mutex_destroy(&pth->mtx);
            pthread_cond_destroy(&pth->cnd);
            return(NULL);
         }break;
      }
   }while(1);

   return(NULL);
}


/*----------------------------------------------------------------------------*/
/* Get the next WP to be computed                                             */
/*----------------------------------------------------------------------------*/

static WrkSct *NexWrk(ParSct *par, int PthIdx)
{
   PthSct *pth = &par->PthTab[ PthIdx ];
   WrkSct *wrk;

   // Remove previous work's tags
   if(pth->wrk)
      SubWrd(par->typ1->NmbDepWrd, pth->wrk->DepWrdTab, par->typ1->RunDepTab);

   // If the wp's buffer is empty search for some new compatible wp to fill in
   if(!par->BufCpt)
   {
      wrk = par->NexWrk;

      while(wrk)
      {
         // Check for dependencies
         if(!AndWrd(par->typ1->NmbDepWrd, wrk->DepWrdTab, par->typ1->RunDepTab))
         {
            par->BufWrk[ par->BufCpt++ ] = wrk;

            // Unlink wp
            if(wrk->pre)
               wrk->pre->nex = wrk->nex;
            else
               par->NexWrk = wrk->nex;

            if(wrk->nex)
               wrk->nex->pre = wrk->pre;

            // Add new work's tags
            AddWrd(par->typ1->NmbDepWrd, wrk->DepWrdTab, par->typ1->RunDepTab);

            if(par->BufCpt == par->BufMax)
               break;
         }

         wrk = wrk->nex;
      }
   }

   // Return the next available wp in buffer and unlink it from the todo list
   return(par->BufCpt ? par->BufWrk[ --par->BufCpt ] : NULL);
}


/*----------------------------------------------------------------------------*/
/* Get the next WP to be computed                                             */
/*----------------------------------------------------------------------------*/

static WrkSct *NexGrn(ParSct *par, int PthIdx)
{
   PthSct *pth = &par->PthTab[ PthIdx ];
   WrkSct *wrk;

   // Remove previous work's tags
   if(pth->wrk)
   {
      ClrBit(par->RunDepTab, pth->wrk->GrnIdx);
      par->ColCpt[ par->GrnCol[ pth->wrk->GrnIdx ] ]--;

      if(!par->ColCpt[ par->GrnCol[ pth->wrk->GrnIdx ] ])
         par->CurCol++;
   }

   wrk = par->NexWrk;

   while(wrk && (par->GrnCol[ wrk->GrnIdx ] <= par->CurCol + 1) )
   {
      // Check for dependencies
      if((wrk->GrnIdx <= par->ColTab[ par->CurCol ][1])
      || ( par->DynSch && (par->CurCol < par->NmbCol)
         && (wrk->GrnIdx <= par->ColTab[ par->CurCol+1 ][1])
         && !AndWrd(par->NmbDepWrd, wrk->DepWrdTab, par->RunDepTab) ) )
      {
         // Unlink wp
         if(wrk->pre)
            wrk->pre->nex = wrk->nex;
         else
            par->NexWrk = wrk->nex;

         if(wrk->nex)
            wrk->nex->pre = wrk->pre;

         // Add new work's tags
         SetBit(par->RunDepTab, wrk->GrnIdx);

         return(wrk);
      }

      wrk = wrk->nex;
   }

   return(NULL);
}


/*----------------------------------------------------------------------------*/
/* Allocate a new kind of elements and set work-packages                      */
/*----------------------------------------------------------------------------*/

int NewType(int64_t ParIdx, itg NmbLin)
{
   int      TypIdx = 0;
   itg      i, idx, BigWrkSiz, NmbBigWrk;
   TypSct   *typ;
   ParSct   *par = (ParSct *)ParIdx;

   // Get and check lib parallel instance
   if(!ParIdx)
      return(0);

   if(NmbLin <= 0)
      return(0);

   // Search for a free type structure
   for(i=1;i<=MaxTyp;i++)
      if(!par->TypTab[i].NmbLin)
      {
         TypIdx = i;
         break;
      }

   if(!TypIdx)
      return(0);

   typ = &par->TypTab[ TypIdx ];
   typ->NmbLin = NmbLin;
   typ->MaxNmbLin = NmbLin * par->SizMul;
   typ->NexGrp = NULL;

   // Compute the size of small work-packages
   if(NmbLin >= par->NmbSmlBlk * par->NmbCpu)
   {
      typ->SmlWrkSiz = NmbLin / (par->NmbSmlBlk * par->NmbCpu);
      typ->NmbSmlWrk = NmbLin / typ->SmlWrkSiz;

      if(NmbLin != typ->NmbSmlWrk * typ->SmlWrkSiz)
         typ->NmbSmlWrk++;
   }
   else
   {
      typ->SmlWrkSiz = NmbLin;
      typ->NmbSmlWrk = 1;
   }

   if(!(typ->SmlWrkTab = LPL_calloc(par->lmb, typ->NmbSmlWrk * par->SizMul , sizeof(WrkSct))))
      return(0);

   // Set small work-packages
   idx = 0;

   for(i=0;i<typ->NmbSmlWrk;i++)
   {
      typ->SmlWrkTab[i].BegIdx = idx + 1;
      typ->SmlWrkTab[i].EndIdx = idx + typ->SmlWrkSiz;
      idx += typ->SmlWrkSiz;
   }

   typ->SmlWrkTab[ typ->NmbSmlWrk - 1 ].EndIdx = NmbLin;

   // Compute the size of big work-packages
   if(!(typ->BigWrkTab = LPL_calloc(par->lmb, par->NmbCpu * par->SizMul , sizeof(WrkSct))))
      return(0);

   // Compute the size of big work-packages
	if(NmbLin >= par->NmbCpu)
	{
		BigWrkSiz = NmbLin / par->NmbCpu;
		NmbBigWrk = par->NmbCpu;
	}
	else
	{
		BigWrkSiz = NmbLin;
		NmbBigWrk = 1;
	}

	// Set big work-packages
	idx = 0;

	for(i=0;i<NmbBigWrk;i++)
	{
		typ->BigWrkTab[i].ItlTab[0][0] = idx + 1;
		typ->BigWrkTab[i].ItlTab[0][1] = idx + BigWrkSiz;
		idx += BigWrkSiz;
	}

	typ->BigWrkTab[ NmbBigWrk - 1 ].ItlTab[0][1] = NmbLin;

   return(TypIdx);
}


/*----------------------------------------------------------------------------*/
/* Adaptive big block resizing based on run time                              */
/*----------------------------------------------------------------------------*/

static void UpdBlkSiz(ParSct *par, TypSct *typ)
{
   int      i, CurBlk = 0;
   itg      CurLin = 0, NewBlk[ MaxPth ][2];
   double   sca, AvgTim, RemTim, TotTim = 0.;

   for(i=0;i<par->NmbCpu;i++)
      TotTim += typ->BigWrkTab[i].RunTim;

   AvgTim = TotTim / par->NmbCpu;

   for(i=0;i<par->NmbCpu;i++)
   {
      NewBlk[i][0] = CurLin + 1;
      TotTim = 0;

      do
      {
         RemTim = (typ->BigWrkTab[ CurBlk ].RunTim
                * (typ->BigWrkTab[ CurBlk ].ItlTab[0][1] - CurLin))
                / (typ->BigWrkTab[ CurBlk ].ItlTab[0][1]
                - typ->BigWrkTab[ CurBlk ].ItlTab[0][0]);

         if(TotTim + RemTim < AvgTim)
         {
            CurLin = typ->BigWrkTab[ CurBlk ].ItlTab[0][1];
            TotTim += RemTim;
            CurBlk++;
         }
         else
         {
            sca = (AvgTim - TotTim) / typ->BigWrkTab[ CurBlk ].RunTim;
            CurLin += sca * (typ->BigWrkTab[ CurBlk ].ItlTab[0][1]
                           - typ->BigWrkTab[ CurBlk ].ItlTab[0][0]);
            break;
         }
      }while(TotTim < AvgTim);

      NewBlk[i][1] = CurLin;
   }

   NewBlk[ par->NmbCpu - 1 ][1] = typ->NmbLin;

	for(i=0;i<par->NmbCpu;i++)
   {
      typ->BigWrkTab[i].ItlTab[0][0] = NewBlk[i][0];
      typ->BigWrkTab[i].ItlTab[0][1] = NewBlk[i][1];
   }
}


/*----------------------------------------------------------------------------*/
/* Resize a data type up to a two fold increase                               */
/*----------------------------------------------------------------------------*/

int ResizeType(int64_t ParIdx, int TypIdx, itg NmbLin)
{
   itg      i, idx, BigWrkSiz, NmbBigWrk;
   TypSct   *typ;
   ParSct   *par = (ParSct *)ParIdx;

   // Get and check lib parallel instance
   if(!ParIdx)
      return(0);

   // Check bounds and free mem
   if( (TypIdx < 1) || (TypIdx > MaxTyp) )
      return(0);

   typ = &par->TypTab[ TypIdx ];

   if( (NmbLin < typ->NmbLin) || (NmbLin > typ->MaxNmbLin) )
      return(0);

   // Set small work-packages
   i = typ->NmbSmlWrk;
   idx = typ->NmbLin;
   typ->NmbLin = NmbLin;

   while(idx < NmbLin)
   {
      typ->SmlWrkTab[i].BegIdx = idx + 1;
      typ->SmlWrkTab[i].EndIdx = idx + typ->SmlWrkSiz;
      i++;
      idx += typ->SmlWrkSiz;
      typ->NmbSmlWrk++;
   }

   typ->SmlWrkTab[ typ->NmbSmlWrk - 1 ].EndIdx = NmbLin;

   // Compute the size of big work-packages
	if(NmbLin >= par->NmbCpu)
	{
		BigWrkSiz = NmbLin / par->NmbCpu;
		NmbBigWrk = par->NmbCpu;
	}
	else
	{
		BigWrkSiz = NmbLin;
		NmbBigWrk = 1;
	}

	// Set big work-packages
	idx = 0;

	for(i=0;i<NmbBigWrk;i++)
	{
		typ->BigWrkTab[i].ItlTab[0][0] = idx + 1;
		typ->BigWrkTab[i].ItlTab[0][1] = idx + BigWrkSiz;
		idx += BigWrkSiz;
	}

	typ->BigWrkTab[ NmbBigWrk - 1 ].ItlTab[0][1] = NmbLin;

   return(TypIdx);
}


/*----------------------------------------------------------------------------*/
/* Set big work blocks interleaving                                           */
/*----------------------------------------------------------------------------*/

static void SetItlBlk(ParSct *par, TypSct *typ)
{
   itg      i, j, BegIdx, EndIdx, CpuIdx = 0, PagIdx[ MaxPth ] = {0};
   double   ItlSiz, ItlIdx = 0.;

   // Set big WP interleaved indices and block sizes if requested
   if(par->NmbItlBlk)
      ItlSiz = (double)typ->NmbLin / (double)(par->NmbItlBlk * par->NmbCpu);
   else if(par->ItlBlkSiz)
   {
      par->NmbItlBlk = typ->NmbLin / (par->ItlBlkSiz * par->NmbCpu);
      ItlSiz = (double)par->ItlBlkSiz;
   }

   // In case the block or lines numbers is too small, deactivate interleaving
   if(!par->NmbItlBlk || !ItlSiz)
   {
      par->NmbItlBlk = 1;
      ItlSiz = (double)typ->NmbLin / (double)(par->NmbCpu);
   }

   for(i=0;i<par->NmbCpu;i++)
   {
      for(j=0;j<par->NmbItlBlk;j++)
      {
         BegIdx = (itg)(ItlIdx + 1.);
         EndIdx = (itg)(ItlIdx + ItlSiz);
         ItlIdx += ItlSiz;

         if(BegIdx <= EndIdx)
         {
            typ->BigWrkTab[ CpuIdx ].ItlTab[ PagIdx[ CpuIdx ] ][0] = BegIdx;
            typ->BigWrkTab[ CpuIdx ].ItlTab[ PagIdx[ CpuIdx ] ][1] = EndIdx;
            PagIdx[ CpuIdx ]++;
            CpuIdx++;
            CpuIdx = CpuIdx % par->NmbCpu;
         }
      }
   }

   typ->BigWrkTab[ par->NmbCpu - 1 ].ItlTab[ par->NmbItlBlk - 1 ][1] = typ->NmbLin;
}


/*----------------------------------------------------------------------------*/
/* Add this kind of element to the free-list                                  */
/*----------------------------------------------------------------------------*/

void FreeType(int64_t ParIdx, int TypIdx)
{
   TypSct *typ;
   ParSct *par = (ParSct *)ParIdx;
   GrpSct *grp, *NexGrp;

   // Get and check lib parallel instance
   if(!ParIdx)
      return;

   // Check bounds and free mem
   if( (TypIdx < 1) || (TypIdx > MaxTyp) )
      return;

   typ = &par->TypTab[ TypIdx ];

   if(typ->SmlWrkTab)
      LPL_free(par->lmb, typ->SmlWrkTab);

   if(typ->BigWrkTab)
      LPL_free(par->lmb, typ->BigWrkTab);

   if(typ->RunDepTab)
      LPL_free(par->lmb, typ->RunDepTab);

   if(typ->DepWrdMat)
      LPL_free(par->lmb, typ->DepWrdMat);

   NexGrp = typ->NexGrp;

   while((grp = NexGrp))
   {
      NexGrp = grp->nex;
      LPL_free(par->lmb, grp);
   }

   memset(typ, 0, sizeof(TypSct));
}


/*----------------------------------------------------------------------------*/
/* Allocate a dependency matrix linking both types                            */
/*----------------------------------------------------------------------------*/

int BeginDependency(int64_t ParIdx, int TypIdx1, int TypIdx2)
{
   int i;
   TypSct *typ1, *typ2;
   ParSct *par = (ParSct *)ParIdx;

   // Get and check lib parallel instance
   if(!ParIdx)
      return(0);

   // Check bounds
   par->CurTyp = typ1 = &par->TypTab[ TypIdx1 ];
   par->DepTyp = typ2 = &par->TypTab[ TypIdx2 ];

   if( (TypIdx1 < 1) || (TypIdx1 > MaxTyp) || (TypIdx2 < 1)
   ||  (TypIdx2 > MaxTyp) || (typ1 == typ2) || !typ1->NmbLin || !typ2->NmbLin)
   {
      return(0);
   }

   // Compute dependency table's size
   if( (typ2->NmbLin >= par->NmbDepBlk * par->NmbCpu)
   &&  (typ2->NmbLin >= typ1->DepWrkSiz * 32) )
   {
      typ1->DepWrkSiz = typ2->NmbLin / (par->NmbDepBlk * par->NmbCpu);
      typ1->NmbDepWrd = typ2->NmbLin / (typ1->DepWrkSiz * 32);

      if(typ2->NmbLin != typ1->NmbDepWrd * typ1->DepWrkSiz * 32)
         typ1->NmbDepWrd++;
   }
   else
   {
      typ1->DepWrkSiz = typ2->NmbLin;
      typ1->NmbDepWrd = 1;
   }

   // Allocate a global dependency table
   if(!(typ1->DepWrdMat =
      LPL_calloc(par->lmb, typ1->NmbSmlWrk * typ1->NmbDepWrd * par->SizMul, sizeof(int))))
   {
      return(0);
   }

   // Then spread sub-tables among WP
   for(i=0;i<typ1->NmbSmlWrk;i++)
   {
      typ1->SmlWrkTab[i].NmbDep = 0;
      typ1->SmlWrkTab[i].DepWrdTab =
         &typ1->DepWrdMat[ i * typ1->NmbDepWrd * par->SizMul ];
   }

   // Allocate a running tags table
   if(!(typ1->RunDepTab = LPL_calloc(par->lmb, typ1->NmbDepWrd * par->SizMul, sizeof(int))))
      return(0);

   return(typ1->NmbDepWrd);
}


/*----------------------------------------------------------------------------*/
/* Type1 element idx1 depends on type2 element idx2                           */
/*----------------------------------------------------------------------------*/

int AddDependency(int64_t ParIdx, itg idx1, itg idx2)
{
   WrkSct *wrk;
   ParSct *par = (ParSct *)ParIdx;

   // Get and check lib parallel instance
   if( !par || !par->CurTyp || !par->CurTyp->SmlWrkSiz
   ||  !par->CurTyp->DepWrkSiz || !par->DepTyp || (idx1 < 1)
   ||  (idx1 > par->CurTyp->NmbLin) || (idx2 < 1)
   ||  (idx2 > par->DepTyp->NmbLin) )
   {
      return(0);
   }

   // Set and count dependency bit
   wrk = &par->CurTyp->SmlWrkTab[ (idx1-1) / par->CurTyp->SmlWrkSiz ];

   if(!SetBit(wrk->DepWrdTab, (idx2-1) / par->CurTyp->DepWrkSiz ))
      wrk->NmbDep++;

   return(wrk->NmbDep);
}


/*----------------------------------------------------------------------------*/
/* Set all to all dependencies without any pre-checking                       */
/*----------------------------------------------------------------------------*/

void AddDependencyFast( int64_t ParIdx, int NmbTyp1, itg *TabIdx1,
                        int NmbTyp2, itg *TabIdx2 )
{
   int i, j;
   ParSct *par = (ParSct *)ParIdx;
   WrkSct *wrk;

   for(i=0;i<NmbTyp1;i++)
   {
      wrk = &par->CurTyp->SmlWrkTab[ (TabIdx1[i] - 1) / par->CurTyp->SmlWrkSiz ];

      for(j=0;j<NmbTyp2;j++)
         if( !SetBit(wrk->DepWrdTab, (TabIdx2[j] - 1) / par->CurTyp->DepWrkSiz ) )
            wrk->NmbDep++;
   }
}


/*----------------------------------------------------------------------------*/
/* Type1 element idx1 depends on type2 element idx2                           */
/*----------------------------------------------------------------------------*/

int UpdateDependency(int64_t ParIdx, int TypIdx1, int TypIdx2,
                     itg idx1, itg idx2 )
{
   WrkSct *wrk;
   ParSct *par = (ParSct *)ParIdx;
   TypSct *typ1, *typ2;

   // Get and check lib parallel instance
   if(!ParIdx)
      return(0);

   // Check bounds
   typ1 = &par->TypTab[ TypIdx1 ];
   typ2 = &par->TypTab[ TypIdx2 ];

   if( (TypIdx1 < 1) || (TypIdx1 > MaxTyp) || (TypIdx2 < 1)
   ||  (TypIdx2 > MaxTyp) || (typ1 == typ2) || !typ1->NmbLin || !typ2->NmbLin
   ||  !typ1->SmlWrkSiz || !typ1->DepWrkSiz || (idx1 < 1)
   ||  (idx1 > typ1->MaxNmbLin) || (idx2 < 1) || (idx2 > typ2->MaxNmbLin) )
   {
      return(0);
   }

   // Set and count dependency bit
   wrk = &typ1->SmlWrkTab[ (idx1-1) / typ1->SmlWrkSiz ];

   if(!SetBit(wrk->DepWrdTab, (idx2-1) / typ1->DepWrkSiz ))
      wrk->NmbDep++;

   return(wrk->NmbDep);
}


/*----------------------------------------------------------------------------*/
/* Same as above, without any pre-checking of data                            */
/*----------------------------------------------------------------------------*/

void UpdateDependencyFast( int64_t ParIdx,  int TypIdx1, int NmbTyp1,
                           itg *TabIdx1, int TypIdx2, int NmbTyp2,
                           itg *TabIdx2 )
{
   int i, j;
   ParSct *par = (ParSct *)ParIdx;
   TypSct *typ1 = &par->TypTab[ TypIdx1 ];
   WrkSct *wrk;
   (void)(TypIdx2);

   for(i=0;i<NmbTyp1;i++)
   {
      wrk = &typ1->SmlWrkTab[ (TabIdx1[i] - 1) / typ1->SmlWrkSiz ];

      for(j=0;j<NmbTyp2;j++)
         if( !SetBit(wrk->DepWrdTab, (TabIdx2[j] - 1) / typ1->DepWrkSiz ) )
            wrk->NmbDep++;
   }
}


/*----------------------------------------------------------------------------*/
/* Sort wp depending on their number of dependencies                          */
/*----------------------------------------------------------------------------*/

int EndDependency(int64_t ParIdx, float DepSta[2])
{
   int      i, NmbDepBit, TotNmbDep = 0;
   ParSct   *par = (ParSct *)ParIdx;
   TypSct   *typ1, *typ2;

   // Get and check lib parallel instance
   if( !ParIdx || !DepSta )
      return(0);

   // Compute average number of collisions
   DepSta[1] = 0.;
   typ1 = par->CurTyp;
   typ2 = par->DepTyp;

   if(!typ1 || !typ2 || !typ1->DepWrkSiz)
      return(0);

   for(i=0;i<typ1->NmbSmlWrk;i++)
   {
      TotNmbDep += typ1->SmlWrkTab[i].NmbDep;
      typ1->SmlWrkTab[i].rnd = rand();

      if(typ1->SmlWrkTab[i].NmbDep > DepSta[1])
         DepSta[1] = (float)typ1->SmlWrkTab[i].NmbDep;
   }

   if(!TotNmbDep)
      return(0);

   DepSta[0] = (float)TotNmbDep;

   // Compute stats
   if(typ2->NmbLin >= typ1->DepWrkSiz)
   {
      NmbDepBit = typ2->NmbLin / typ1->DepWrkSiz;

      if(typ2->NmbLin - NmbDepBit * typ1->DepWrkSiz)
         NmbDepBit++;
   }
   else
      NmbDepBit = 1;

   if(!NmbDepBit)
      return(0);

   DepSta[0] = 100 * DepSta[0] / (typ1->NmbSmlWrk * NmbDepBit);
   DepSta[1] = 100 * DepSta[1] / NmbDepBit;

   // Sort WP from highest collision number to the lowest
   if(par->WrkSizSrt && par->DynSch)
      qsort(typ1->SmlWrkTab, typ1->NmbSmlWrk, sizeof(WrkSct), CmpWrk);

   // If the dynamic scheduling is disabled, set static WP
   if(!par->DynSch && !SetGrp(par, typ1))
      return(0);

   return(1);
}


/*----------------------------------------------------------------------------*/
/* Sort wp depending on their number of dependencies                          */
/*----------------------------------------------------------------------------*/

void GetDependencyStats(int64_t ParIdx, int TypIdx1,
                        int TypIdx2, float DepSta[2])
{
   int i, NmbDepBit, TotNmbDep = 0;
   ParSct *par = (ParSct *)ParIdx;
   TypSct *typ1, *typ2;

   // Get and check lib parallel instance
   if(!ParIdx)
      return;

   // Check bounds
   typ1 = &par->TypTab[ TypIdx1 ];
   typ2 = &par->TypTab[ TypIdx2 ];

   if( (TypIdx1 < 1) || (TypIdx1 > MaxTyp) || (TypIdx2 < 1)
   ||  (TypIdx2 > MaxTyp) || (typ1 == typ2) || !typ1->NmbLin
   ||  !typ2->NmbLin )
   {
      return;
   }

   // Compute average number of collisions
   DepSta[1] = 0.;

   for(i=0;i<typ1->NmbSmlWrk;i++)
   {
      TotNmbDep += typ1->SmlWrkTab[i].NmbDep;

      if(typ1->SmlWrkTab[i].NmbDep > DepSta[1])
         DepSta[1] = (float)typ1->SmlWrkTab[i].NmbDep;
   }

   DepSta[0] = (float)TotNmbDep;

   // Compute stats
   if(typ2->NmbLin >= typ1->DepWrkSiz)
   {
      NmbDepBit = typ2->NmbLin / typ1->DepWrkSiz;

      if(typ2->NmbLin - NmbDepBit * typ1->DepWrkSiz)
         NmbDepBit++;
   }
   else
      NmbDepBit = 1;

   DepSta[0] = 100 * DepSta[0] / (typ1->NmbSmlWrk * NmbDepBit);
   DepSta[1] = 100 * DepSta[1] / NmbDepBit;
}


/*----------------------------------------------------------------------------*/
/* Halve the number of small blocks by compining pairs of consecutive blocks  */
/*----------------------------------------------------------------------------*/

int HalveSmallBlocks(int64_t ParIdx, int TypIdx1, int TypIdx2)
{
   int i, j;
   ParSct *par = (ParSct *)ParIdx;
   WrkSct *EvnWrk, *OddWrk, *NewWrk;
   TypSct *typ1, *typ2;

   // Get and check lib parallel instance
   if(!ParIdx)
      return(0);

   // Check bounds
   typ1 = &par->TypTab[ TypIdx1 ];
   typ2 = &par->TypTab[ TypIdx2 ];

   if( (TypIdx1 < 1) || (TypIdx1 > MaxTyp) || (TypIdx2 < 1)
   ||  (TypIdx2 > MaxTyp) || (typ1 == typ2) || !typ1->NmbLin
   ||  !typ2->NmbLin )
   {
      return(0);
   }

   // Do not halve the number of blocks if there is only one left
   // nor if the small blocks have been sorted
   if(typ1->NmbSmlWrk < 2 || par->WrkSizSrt)
      return(0);

   // For each new block, compute the logical OR between two consecutive old blocks
   // The new data is copied on top of former one
   for(i=0;i<typ1->NmbSmlWrk;i+=2)
   {
      EvnWrk = &typ1->SmlWrkTab[ i     ];
      OddWrk = &typ1->SmlWrkTab[ i + 1 ];
      NewWrk = &typ1->SmlWrkTab[ i / 2 ];
      NewWrk->BegIdx = EvnWrk->BegIdx;
      NewWrk->EndIdx = OddWrk->EndIdx;

      for(j=0;j<typ1->NmbDepWrd;j++)
         if(i+1 < typ1->NmbSmlWrk)
            NewWrk->DepWrdTab[j] = EvnWrk->DepWrdTab[j] | OddWrk->DepWrdTab[j];
         else
            NewWrk->DepWrdTab[j] = EvnWrk->DepWrdTab[j];
   }

   // Halve the number of blocks and add one if the number was odd
   typ1->SmlWrkSiz *= typ1->NmbSmlWrk;

   if(typ1->NmbSmlWrk & 1)
      typ1->NmbSmlWrk = typ1->NmbSmlWrk / 2 + 1;
   else
      typ1->NmbSmlWrk /= 2;

   typ1->SmlWrkSiz /= typ1->NmbSmlWrk;

   return(typ1->NmbSmlWrk);
}


/*----------------------------------------------------------------------------*/
/* Halve the number of dependency words by ORing consecutive pairs of bits    */
/*----------------------------------------------------------------------------*/

int HalveDependencyBlocks(int64_t ParIdx, int TypIdx1, int TypIdx2)
{
   int i, j;
   WrkSct *wrk;
   ParSct *par = (ParSct *)ParIdx;
   TypSct *typ1, *typ2;

   // Get and check lib parallel instance
   if(!ParIdx)
      return(0);

   // Check bounds
   typ1 = &par->TypTab[ TypIdx1 ];
   typ2 = &par->TypTab[ TypIdx2 ];

   if( (TypIdx1 < 1) || (TypIdx1 > MaxTyp) || (TypIdx2 < 1)
   ||  (TypIdx2 > MaxTyp) || (typ1 == typ2) || !typ1->NmbLin
   ||  !typ2->NmbLin )
   {
      return(0);
   }

   // Do not halve the number of blocks if there is only one left
   if(typ1->NmbDepWrd < 2)
      return(0);

   for(i=0;i<typ1->NmbSmlWrk;i++)
   {
      wrk = &typ1->SmlWrkTab[i];

      for(j=0;j<typ1->NmbDepWrd;j+=2)
         if(GetBit(wrk->DepWrdTab, j) || GetBit(wrk->DepWrdTab, j+1))
            SetBit(wrk->DepWrdTab, j/2);
         else
            ClrBit(wrk->DepWrdTab, j/2);
   }

   typ1->DepWrkSiz *= typ1->NmbDepWrd;

   if(typ1->NmbDepWrd & 1)
      typ1->NmbDepWrd = typ1->NmbDepWrd / 2 + 1;
   else
      typ1->NmbDepWrd /= 2;

   typ1->DepWrkSiz /= typ1->NmbDepWrd;

   return(typ1->NmbDepWrd * 32);
}


/*----------------------------------------------------------------------------*/
/* Return a block of elements' begin and ending indices                       */
/*----------------------------------------------------------------------------*/

int GetBigBlkNfo(int64_t ParIdx, int TypIdx, int BlkIdx, int *BegIdx, int *EndIdx)
{
   ParSct *par = (ParSct *)ParIdx;
   TypSct *typ;

   // Get and check lib parallel instance, type and block indices
   if(!ParIdx || TypIdx < 1 || TypIdx >= MaxTyp || BlkIdx < 0 || BlkIdx >= par->NmbCpu)
      return(0);

   // Check if thi type has been allocated
   typ = &par->TypTab[ TypIdx ];

   if(!typ->NmbLin)
      return(0);

   // Set begin and end indices
   *BegIdx = typ->BigWrkTab[ BlkIdx ].ItlTab[0][0];
   *EndIdx = typ->BigWrkTab[ BlkIdx ].ItlTab[0][1];

   return(1);
}


/*----------------------------------------------------------------------------*/
/* Return the block ID containing the given element ID                        */
/*----------------------------------------------------------------------------*/

int GetBlkIdx(int64_t ParIdx, int typ, int idx)
{
   ParSct *par = (ParSct *)ParIdx;

   // Get and check lib parallel instance
   if(!ParIdx)
      return(-1);

   return( (idx - 1) / par->TypTab[ typ ].SmlWrkSiz );
}


/*----------------------------------------------------------------------------*/
/* Check whether two blocks are dependant or not from each others             */
/*----------------------------------------------------------------------------*/

int ChkBlkDep(int64_t ParIdx, int TypIdx, int blk1, int blk2)
{
   ParSct *par = (ParSct *)ParIdx;
   TypSct *typ;

   // Get and check lib parallel instance
   if(!ParIdx)
      return(-1);

   typ = &par->TypTab[ TypIdx ];

   return(AndWrd( typ->NmbDepWrd,
                  typ->SmlWrkTab[ blk1 ].DepWrdTab,
                  typ->SmlWrkTab[ blk2 ].DepWrdTab ));
}


/*----------------------------------------------------------------------------*/
/* Test and set a bit in a multibyte word                                     */
/*----------------------------------------------------------------------------*/

static int SetBit(int *tab, int idx)
{
   int res = ( tab[ idx >> 5 ] & (1UL << (idx & 31)) );
   tab[ idx >> 5 ] |= 1UL << (idx & 31);
   return(res);
}


/*----------------------------------------------------------------------------*/
/* Test a bit in a multibyte word                                             */
/*----------------------------------------------------------------------------*/

static int GetBit(int *tab, int idx)
{
   return( tab[ idx >> 5 ] & (1UL << (idx & 31)) );
}


/*----------------------------------------------------------------------------*/
/* Clear a bit in a multibyte word                                            */
/*----------------------------------------------------------------------------*/

static void ClrBit(int *tab, int idx)
{
   tab[ idx >> 5 ] &= ~(1UL << (idx & 31));
}


/*----------------------------------------------------------------------------*/
/* Check wether two WP share common resources -> locked                       */
/*----------------------------------------------------------------------------*/

static int AndWrd(int NmbWrd, int *wrd1, int *wrd2)
{
   int i;

   for(i=0;i<NmbWrd;i++)
      if(wrd1[i] & wrd2[i])
         return(1);

   return(0);
}


/*----------------------------------------------------------------------------*/
/* Logical OR between two multibyte words                                     */
/*----------------------------------------------------------------------------*/

static void AddWrd(int NmbWrd, int *SrcWrd, int *DstWrd)
{
   int i;

   for(i=0;i<NmbWrd;i++)
      DstWrd[i] |= SrcWrd[i];
}


/*----------------------------------------------------------------------------*/
/* Exclusive OR between two multibyte words                                   */
/*----------------------------------------------------------------------------*/

static void SubWrd(int NmbWrd, int *SrcWrd, int *DstWrd)
{
   int i;

   for(i=0;i<NmbWrd;i++)
      DstWrd[i] &= ~SrcWrd[i];
}


/*----------------------------------------------------------------------------*/
/* Clear a multibyte word                                                     */
/*----------------------------------------------------------------------------*/

static void ClrWrd(int NmbWrd, int *wrd)
{
   int i;

   for(i=0;i<NmbWrd;i++)
      wrd[i] = 0;
}


/*----------------------------------------------------------------------------*/
/* Copy a multibyte word into another one                                     */
/*----------------------------------------------------------------------------*/

static void CpyWrd(int NmbWrd, int *SrcWrd, int *DstWrd)
{
   int i;

   for(i=0;i<NmbWrd;i++)
      DstWrd[i] = SrcWrd[i];
}


/*----------------------------------------------------------------------------*/
/* Compare two workpackages number of bits                                    */
/*----------------------------------------------------------------------------*/

int CmpWrk(const void *ptr1, const void *ptr2)
{
   WrkSct *w1, *w2;

   w1 = (WrkSct *)ptr1;
   w2 = (WrkSct *)ptr2;

   if(w1->NmbDep > w2->NmbDep)
      return(-1);
   else if(w1->NmbDep < w2->NmbDep)
      return(1);
   else
      return(0);
}


/*----------------------------------------------------------------------------*/
/* Generate static scheduling groups of small WP for each thread              */
/*----------------------------------------------------------------------------*/

static int SetGrp(ParSct *par, TypSct *typ)
{
   int      i, NmbSmlWrk, *GrpWrd, *AllWrd, *TstWrd;
   int      IncFlg, siz = typ->NmbDepWrd;
   GrpSct   *grp;
   WrkSct   *wrk, *NexWrk;

   // Initialize the static group list
   NmbSmlWrk = typ->NmbSmlWrk;
   typ->NmbGrp = 0;
   typ->NexGrp = NULL;

   // Allocate a dependency word to contain all threads
   if(!(GrpWrd = LPL_malloc(par->lmb, par->NmbCpu * siz * sizeof(int))))
      return(0);

   // Allocate a dependency word to concatenate all thread words
   if(!(AllWrd = LPL_malloc(par->lmb, siz * sizeof(int))))
      return(0);

   // Allocate a working dependency word for testings
   if(!(TstWrd = LPL_malloc(par->lmb, siz * sizeof(int))))
      return(0);

   // Link all WP together to make a free list
   NexWrk = &typ->SmlWrkTab[0];

   for(i=0;i<typ->NmbSmlWrk;i++)
   {
      typ->SmlWrkTab[i].pre = &typ->SmlWrkTab[ i - 1 ];
      typ->SmlWrkTab[i].nex = &typ->SmlWrkTab[ i + 1 ];
   }

   typ->SmlWrkTab[0].pre = NULL;
   typ->SmlWrkTab[ typ->NmbSmlWrk - 1 ].nex = NULL;

   // Kepp on creating new groups as long as there are free WP
   do
   {
      // Allocate a new group and link it
      if(!(grp = LPL_calloc(par->lmb, 1, sizeof(GrpSct))))
         return(0);

      grp->nex = typ->NexGrp;
      typ->NexGrp = grp;
      grp->idx = ++typ->NmbGrp;

      // Reset all the dependencies words
      memset(GrpWrd, 0, par->NmbCpu * siz * sizeof(int));
      memset(AllWrd, 0, siz * sizeof(int));

      // Keep on adding WP to each groups' threads until they reach WrkPerGrp
      do
      {
         // Could we add some WP to this group during this iteration ?
         IncFlg = 0;

         // Loop over the threads and try adding them as many WP as possible
         for(i=0;i<par->NmbCpu;i++)
         {
            // Do not add WP is this threads has reached its target number
            if(!NexWrk || (grp->NmbSmlWrk[i] >= WrkPerGrp))
               continue;

            // Initialize the WP pointer with the first
            // available one from the linked list
            wrk = NexWrk;

            // Then loop over the linked list of available WP and add them to
            // the curent group-thread list as long as they do not interfere
            // with the combined dependency word of other threads
            do
            {
               // To avoid making a AND with all threads' words, we remode (XOR)
               // the tested thread' word from the combined one and test (AND)
               // it against this WP's word
               CpyWrd(siz, AllWrd, TstWrd);
               SubWrd(siz, &GrpWrd[ i * siz ], TstWrd);

               if(!AndWrd(siz, wrk->DepWrdTab, TstWrd))
               {
                  // If this WP is compatible, add its word to the thread
                  // dependency word and to the combined one
                  AddWrd(siz, wrk->DepWrdTab, AllWrd);
                  AddWrd(siz, wrk->DepWrdTab, &GrpWrd[ i * siz ]);

                  // Add this WP to the list, decrease the number of available
                  // WPs and set the flag to indicate that we found some work
                  grp->SmlWrkTab[i][ grp->NmbSmlWrk[i] ] = wrk;
                  grp->NmbSmlWrk[i]++;
                  NmbSmlWrk--;
                  IncFlg = 1;

                  // Unlink the WP from the available WP linked list
                  if(wrk->pre)
                     wrk->pre->nex = wrk->nex;
                  else
                     NexWrk = wrk->nex;

                  if(wrk->nex)
                     wrk->nex->pre = wrk->pre;
               }

               // Stop adding WP to this group-thread if its table is full
               if(grp->NmbSmlWrk[i] >= WrkPerGrp)
                  break;

               wrk = wrk->nex;
            }while(wrk && NmbSmlWrk);
         }
      }while(IncFlg && NmbSmlWrk);
   }while(NmbSmlWrk);

   LPL_free(par->lmb, GrpWrd);
   LPL_free(par->lmb, AllWrd);
   LPL_free(par->lmb, TstWrd);

   return(1);
}


/*----------------------------------------------------------------------------*/
/* Launch the loop prc on typ1 element depending on typ2                      */
/*----------------------------------------------------------------------------*/

float LaunchColorGrains(int64_t ParIdx, int typ, void *prc, void *PtrArg)
{
   int      i;
   float    acc = 0.;
   PthSct   *pth;
   ParSct   *par = (ParSct *)ParIdx;

   // Get and check lib parallel instance
   if(!ParIdx)
      return(-1.);

   // Lock acces to global parameters
   pthread_mutex_lock(&par->ParMtx);

   par->cmd = RunGrnWrk;
   par->prc = (void (*)(itg, itg, int, void *))prc;
   par->arg = PtrArg;
   par->typ1 = NULL;
   par->typ2 = NULL;
   par->NexWrk = par->GrnWrkTab;
   par->WrkCpt = 0;
   par->sta[0] = par->sta[1] = 0.;
   par->req = 0;
   par->NmbDep = 0;
   par->CurCol = 1;

   // Clear running wp
   for(i=0;i<par->NmbCpu;i++)
      par->PthTab[i].wrk = NULL;

   ClrWrd(par->NmbGrnWrd, par->RunGrnTab);

   // Build a linked list of wp
   for(i=0;i<par->NmbGrnWrk;i++)
   {
      par->GrnWrkTab[i].pre = &par->GrnWrkTab[ i-1 ];
      par->GrnWrkTab[i].nex = &par->GrnWrkTab[ i+1 ];
      par->GrnWrkTab[i].GrnIdx = i + 1;
      par->GrnWrkTab[i].BegIdx = par->GrnTab[ typ ][i+1][0];
      par->GrnWrkTab[i].EndIdx = par->GrnTab[ typ ][i+1][1];
   }

   par->GrnWrkTab[0].pre = par->GrnWrkTab[ par->NmbGrnWrk - 1 ].nex = NULL;

   for(i=1;i<=par->NmbCol;i++)
      par->ColCpt[i] = par->ColTab[i][1] - par->ColTab[i][0] + 1;


   // Main loop: wake up threads and wait for completion or blocked threads
   do
   {
      // Search for some idle threads
      par->req = 0;

      for(i=0;i<par->NmbCpu;i++)
      {
         pth = &par->PthTab[i];

         if(pth->wrk)
            continue;

         if(!(pth->wrk = NexGrn(par, i)))
         {
            par->req = 1;
            break;
         }

         // Wake up the thread and provide it with a WP list
         pthread_mutex_lock(&pth->mtx);
         pthread_cond_signal(&pth->cnd);
         pthread_mutex_unlock(&pth->mtx);
      }

      // If every WP are done : exit the parallel loop
      if(par->WrkCpt == par->NmbGrnWrk)
         break;

      // Otherwise, wait for a blocked thread
      pthread_cond_wait(&par->ParCnd, &par->ParMtx);
   }while(1);

   pthread_mutex_unlock(&par->ParMtx);

   // Compute the average concurrency factor
   acc = par->sta[0] ? (par->sta[1] / par->sta[0]) : 0;

   // Clear the main datatyp loop to indicate that no LaunchParallel is running
   par->typ1 = 0;

   // Return the concurrency factor
   return(acc);
}


/*----------------------------------------------------------------------------*/
/* Tell the pthread handler to clear a chunmk of memory in parallel           */
/*----------------------------------------------------------------------------*/

int ParallelMemClear(int64_t ParIdx, void *PtrArg, size_t siz)
{
   int      i;
   size_t   StdSiz, EndSiz;
   PthSct   *pth;
   ParSct   *par = (ParSct *)ParIdx;
   char     *tab = (char *)PtrArg;

   // Get and check lib parallel instance, adresse and size
   if(!ParIdx || !tab || (siz < (size_t)par->NmbCpu) )
      return(0);

   // If the memory chunk is too small, clear it sequentially
   if(siz < BigMemSiz || par->NmbCpu == 1)
   {
      memset(PtrArg, 0, siz);
      return(1);
   }

   // Lock acces to global parameters
   pthread_mutex_lock(&par->ParMtx);

   par->cmd = ClrMem;
   par->WrkCpt = 0;

   // Share the size between threads and adjust the value of the final one
   StdSiz = siz / par->NmbCpu;
   EndSiz = siz - StdSiz * (par->NmbCpu - 1);

   // Spread the buffer among each thread and wake then up
   for(i=0;i<par->NmbCpu;i++)
   {
      pth = &par->PthTab[i];
      pth->ClrAdr = &tab[ i * StdSiz ];

      // The last thread processes the remaining size
      if(i < par->NmbCpu - 1)
         pth->ClrMemSiz = StdSiz;
      else
         pth->ClrMemSiz = EndSiz;

      pthread_mutex_lock(&pth->mtx);
      pthread_cond_signal(&pth->cnd);
      pthread_mutex_unlock(&pth->mtx);
   }

   // Wait for each thread to complete
   while(par->WrkCpt < par->NmbCpu)
      pthread_cond_wait(&par->ParCnd, &par->ParMtx);

   pthread_mutex_unlock(&par->ParMtx);

   return(1);
}


/*----------------------------------------------------------------------------*/
/* Tell the pthread handler to copy a chunmk of memory in parallel            */
/*----------------------------------------------------------------------------*/

int ParallelMemCopy(int64_t ParIdx, void *PtrDst, void *PtrSrc, size_t siz)
{
   int      i;
   PthSct   *pth;
   ParSct   *par = (ParSct *)ParIdx;
   size_t   StdSiz, EndSiz;
   char     *DstTab = (char *)PtrDst;
   char     *SrcTab = (char *)PtrSrc;

   // Get and check lib parallel instance, adresse and size
   if(!ParIdx || !PtrDst || !PtrSrc || (siz < (size_t)par->NmbCpu) )
      return(0);

   // If the memory chunk is too small, clear it sequentially
   if(siz < BigMemSiz || par->NmbCpu == 1)
   {
      memcpy(PtrDst, PtrSrc, siz);
      return(1);
   }

   // Lock acces to global parameters
   pthread_mutex_lock(&par->ParMtx);

   par->cmd = CpyMem;
   par->WrkCpt = 0;

   // Share the size between threads and adjust the value of the final one
   StdSiz = siz / par->NmbCpu;
   EndSiz = siz - StdSiz * (par->NmbCpu - 1);

   // Spread the buffer among each thread and wake then up
   for(i=0;i<par->NmbCpu;i++)
   {
      pth = &par->PthTab[i];
      pth->DstAdr = &DstTab[ i * StdSiz ];
      pth->SrcAdr = &SrcTab[ i * StdSiz ];

      // The last thread processes the remaining size
      if(i < par->NmbCpu - 1)
         pth->CpyMemSiz = StdSiz;
      else
         pth->CpyMemSiz = EndSiz;

      pthread_mutex_lock(&pth->mtx);
      pthread_cond_signal(&pth->cnd);
      pthread_mutex_unlock(&pth->mtx);
   }

   // Wait for each thread to complete
   while(par->WrkCpt < par->NmbCpu)
      pthread_cond_wait(&par->ParCnd, &par->ParMtx);

   pthread_mutex_unlock(&par->ParMtx);

   return(1);
}


/*----------------------------------------------------------------------------*/
/* Wait for a condition, launch and detach a user procedure                   */
/*----------------------------------------------------------------------------*/

int LaunchPipeline(  int64_t ParIdx, void *prc,
                     void *PtrArg, int NmbDep, int *DepTab )
{
   int i;
   PipSct *NewPip=NULL;
   ParSct *par = (ParSct *)ParIdx;

   // Get and check lib parallel instance and the number of pipes and dependencies
   if( !ParIdx || (NmbDep > MaxPipDep) || (par->NmbPip >= MaxTotPip) )
      return(0);

   // Allocate and setup a new pipe
   if(!(NewPip = LPL_calloc(par->lmb, 1, sizeof(PipSct))))
      return(0);

   NewPip->prc = prc;
   NewPip->arg = PtrArg;
   NewPip->par = par;
   NewPip->NmbDep = NmbDep;

   for(i=0;i<NmbDep;i++)
      NewPip->DepTab[i] = DepTab[i];

   // In case variable arguments where passed through the common LPlib structure
   // Copy them in the pipeline's own arguments table
   if(par->NmbVarArg)
   {
      NewPip->NmbVarArg = par->NmbVarArg;

      for(i=0;i<par->NmbVarArg;i++)
         NewPip->VarArgTab[i] = par->VarArgTab[i];
   }

   // Lock pipe mutex, increment pipe counter
   // and launch the pipe regardless dependencies
   pthread_mutex_lock(&par->PipMtx);
   NewPip->idx = ++par->NmbPip;
   par->PenPip++;

   if(par->StkSiz)
   {
      pthread_attr_init(&NewPip->atr);
      NewPip->StkSiz = par->StkSiz;
      NewPip->UsrStk = LPL_malloc(par->lmb, par->StkSiz);
#ifdef _WIN32
      pthread_attr_setstackaddr(&NewPip->atr, NewPip->UsrStk);
      pthread_attr_setstacksize(&NewPip->atr, NewPip->StkSiz);
#else
      pthread_attr_setstack(&NewPip->atr, NewPip->UsrStk, NewPip->StkSiz);
#endif
      pthread_create(&NewPip->pth, &NewPip->atr, PipHdl, (void *)NewPip);
   }
   else
      pthread_create(&NewPip->pth, NULL, PipHdl, (void *)NewPip);

   // Make the thread unjoinable to work around a memory leak in some systems
   pthread_detach(NewPip->pth);
   pthread_mutex_unlock(&par->PipMtx);

   return(NewPip->idx);
}


/*----------------------------------------------------------------------------*/
/* Launch a pipeline procedure with variable numbers of arguments             */
/*----------------------------------------------------------------------------*/

int LaunchPipelineMultiArg(int64_t ParIdx, int NmbDep, int *DepTab, 
                           void *prc, int NmbArg, ... )
{
   int i, ret;
   va_list ArgLst;
   ParSct *par = (ParSct *)ParIdx;

   if(NmbArg > MaxVarArg)
      return(-1);

   // Get the variable arguments and store then in the common LPlib structure
   // This is a temporary storage as the will be copied to the pipeline's own
   // arguments table after it has been allocated by LaunchPipeline()
   par->NmbVarArg = NmbArg;
   va_start(ArgLst, NmbArg);

   for(i=0;i<NmbArg;i++)
      par->VarArgTab[i] = va_arg(ArgLst, void *);

   va_end(ArgLst);

   ret = LaunchPipeline(ParIdx, prc, NULL, NmbDep, DepTab);

   // Clear the temporary arguments table
   par->NmbVarArg = 0;

   return(ret);
}


/*----------------------------------------------------------------------------*/
/* Thread handler launching and waitint for user's procedure                  */
/*----------------------------------------------------------------------------*/

static void *PipHdl(void *ptr)
{
   int RunFlg=0, i;
   PipSct *pip = (PipSct *)ptr;
   ParSct *par = pip->par;
   void (*prc)(void *), (*prcgrn)(int, int, int, void *);

   if(pip->GrnIdx)
   {
      prcgrn = (void (*)(int, int, int, void *))pip->prc;
      par->RunPip++;

      prcgrn(pip->BegIdx, pip->EndIdx, pip->GrnIdx, pip->arg);

      pthread_mutex_lock(&par->PipMtx);
      par->PenPip--;
      par->RunPip--;
      LPL_free(par->lmb, pip);
      pthread_mutex_unlock(&par->PipMtx);
   }
   else
   {
      // Wait for conditions to be met
      do
      {
         pthread_mutex_lock(&par->PipMtx);

         if(par->RunPip < par->NmbCpu)
         {
            RunFlg = 1;

            for(i=0;i<pip->NmbDep;i++)
               if(!GetBit(par->PipWrd, pip->DepTab[i]))
               {
                  RunFlg = 0;
                  break;
               }
         }

         if(!RunFlg)
         {
            pthread_mutex_unlock(&par->PipMtx);
#ifdef _WIN32
            Sleep(1);
#else
            usleep(1000);
#endif
         }
      }while(!RunFlg);

   // Execute the user's procedure and set the flag to 2 (done)
      prc = (void (*)(void *))pip->prc;
      par->RunPip++;

      pthread_mutex_unlock(&par->PipMtx);

      if(pip->NmbVarArg)
         CalVarArgPip(pip, pip->prc);
      else
         prc(pip->arg);

      pthread_mutex_lock(&par->PipMtx);

      SetBit(par->PipWrd, pip->idx);

      par->PenPip--;
      par->RunPip--;
      LPL_free(par->lmb, pip);
      pthread_mutex_unlock(&par->PipMtx);
   }

   return(NULL);
}


/*----------------------------------------------------------------------------*/
/* Wait for all pipelined procedures to complete                              */
/*----------------------------------------------------------------------------*/

void WaitPipeline(int64_t ParIdx)
{
   int PenPip;
   ParSct *par = (ParSct *)ParIdx;

   // Get and check lib parallel instance
   if(!ParIdx)
      return;

   do
   {
      pthread_mutex_lock(&par->PipMtx);
      PenPip = par->PenPip;
      pthread_mutex_unlock(&par->PipMtx);
#ifdef _WIN32
      Sleep(1);
#else
      //usleep(1000);
      usleep(10);
#endif
   }while(PenPip);
}


/*----------------------------------------------------------------------------*/
/* Wait for all pipelined procedures to complete                              */
/*----------------------------------------------------------------------------*/

void GetLplibInformation(int64_t ParIdx, int *NmbCpu, int *NmbTyp)
{
   int i;
   ParSct *par = (ParSct *)ParIdx;

   // Get and check lib parallel instance
   if(!ParIdx)
      return;

   *NmbCpu = par->NmbCpu;
   *NmbTyp = 0;

   // Count the number of used types
   for(i=1;i<=MaxTyp;i++)
      if(par->TypTab[i].NmbLin)
         (*NmbTyp)++;
}


/*----------------------------------------------------------------------------*/
/* Call qsort in parallel then call merge sort in serial                      */
/*----------------------------------------------------------------------------*/

void ParallelQsort(  int64_t ParIdx, void *base, size_t nel, size_t width,
                     int (*compar)(const void *, const void *) )
{
   (void)(ParIdx);
#ifdef __MACH__
   psort(base, nel, width, compar);
#else
   qsort(base, nel, width, compar);
#endif
}


/*----------------------------------------------------------------------------*/
/* Compute the hilbert code from 3d coordinates                               */
/*----------------------------------------------------------------------------*/

static void RenPrc(itg BegIdx, itg EndIdx, int PthIdx, ArgSct *arg)
{
   uint64_t IntCrd[3], m=1ULL<<63, cod;
   itg i;
   int j, b, GeoWrd, NewWrd, BitTab[3] = {1,2,4};
   double dbl;
   int rot[8], GeoCod[8]={0,3,7,4,1,2,6,5};  // Z curve = {5,4,7,6,1,0,3,2}
   int HilCod[8][8] = {
      {0,7,6,1,2,5,4,3}, {0,3,4,7,6,5,2,1},
      {0,3,4,7,6,5,2,1}, {2,3,0,1,6,7,4,5},
      {2,3,0,1,6,7,4,5}, {6,5,2,1,0,3,4,7},
      {6,5,2,1,0,3,4,7}, {4,3,2,5,6,1,0,7} };
   (void)(PthIdx);

   for(i=BegIdx; i<=EndIdx; i++)
   {
      // Convert double precision coordinates to integers
      for(j=0;j<3;j++)
      {
         dbl = (arg->crd[i][j] - arg->box[j]) * arg->box[j+3];
         IntCrd[j] = (uint64_t)dbl;
      }

      // Binary hilbert renumbering loop
      cod = 0;

      for(j=0;j<8;j++)
         rot[j] = GeoCod[j];

      for(b=0;b<21;b++)
      {
         GeoWrd = 0;

         for(j=0;j<3;j++)
         {
            if(IntCrd[j] & m)
               GeoWrd |= BitTab[j];

            IntCrd[j] = IntCrd[j]<<1;
         }

         NewWrd = rot[ GeoWrd ];
         cod = cod<<3 | NewWrd;

         for(j=0;j<8;j++)
            rot[j] = HilCod[ NewWrd ][ rot[j] ];
      }

      arg->idx[i][0] = cod;
      arg->idx[i][1] = i;
   }
}


/*----------------------------------------------------------------------------*/
/* Comparison of two items for the qsort                                      */
/*----------------------------------------------------------------------------*/

int CmpPrc(const void *a, const void *b)
{
   uint64_t *pa = (uint64_t *)a, *pb = (uint64_t *)b;
   return(*pa > *pb ? 1 : -1);
}

void PipSrt(PipArgSct *arg)
{
   qsort(arg->base, arg->nel, arg->width, arg->compar);
}


/*----------------------------------------------------------------------------*/
/* Renumber a set of coordinates through a Hilbert SFC                        */
/*----------------------------------------------------------------------------*/

int HilbertRenumbering( int64_t ParIdx, itg NmbLin, double box[6],
                        double (*crd)[3], uint64_t (*idx)[2] )
{
   int      i, NewTyp;
   double   len = pow(2,64);
   ArgSct   arg;

   // Get and check lib parallel instance
   if(!ParIdx)
     return(0);

   // Setup the bounding box and a data type,
   // then give a Hilbert code to each entries
   NewTyp = NewType(ParIdx, NmbLin);
   arg.crd = crd;
   arg.idx = idx;
   arg.box[0] = box[0];
   arg.box[1] = box[1];
   arg.box[2] = box[2];
   arg.box[3] = len / (box[3] - box[0]);
   arg.box[4] = len / (box[4] - box[1]);
   arg.box[5] = len / (box[5] - box[2]);


   if(NmbLin < 10000)
      RenPrc(1, NmbLin, 0, (void *)&arg);
   else
      LaunchParallel(ParIdx, NewTyp, 0, (void *)RenPrc, (void *)&arg);

   qsort(&idx[1][0], NmbLin, 2 * sizeof(int64_t), CmpPrc);

   for(i=1;i<=NmbLin;i++)
      idx[ idx[i][1] ][0] = i;

   return(1);
}


/*----------------------------------------------------------------------------*/
/* Compute the hilbert code from 2d coordinates                               */
/*----------------------------------------------------------------------------*/

static void RenPrc2D(itg BegIdx, itg EndIdx, int PthIdx, ArgSct *arg)
{
   uint64_t IntCrd[2], m=1ULL<<62, cod;
   itg i;
   int j, b, GeoWrd, NewWrd, BitTab[2] = {1,2};
   double dbl;
   int rot[4], GeoCod[4]={1,2,0,3};
   int HilCod[4][4] = {{0,3,2,1}, {0,1,2,3}, {0,1,2,3}, {2,1,0,3}};
   (void)(PthIdx);

   for(i=BegIdx; i<=EndIdx; i++)
   {
      // Convert double precision coordinates to integers
      for(j=0;j<2;j++)
      {
         dbl = (arg->crd2[i][j] - arg->box[j]) * arg->box[j+2];
         IntCrd[j] = (uint64_t)dbl;
      }

      // Binary hilbert renumbering loop
      cod = 0;

      for(j=0;j<4;j++)
         rot[j] = GeoCod[j];

      for(b=0;b<31;b++)
      {
         GeoWrd = 0;

         for(j=0;j<2;j++)
         {
            if(IntCrd[j] & m)
               GeoWrd |= BitTab[j];

            IntCrd[j] = IntCrd[j]<<1;
         }

         NewWrd = rot[ GeoWrd ];
         cod = cod<<2 | NewWrd;

         for(j=0;j<4;j++)
            rot[j] = HilCod[ NewWrd ][ rot[j] ];
      }

      arg->idx[i][0] = cod;
      arg->idx[i][1] = i;
   }
}


/*----------------------------------------------------------------------------*/
/* Renumber a set of 2D coordinates through a Hilbert SFC                     */
/*----------------------------------------------------------------------------*/

int HilbertRenumbering2D(  int64_t ParIdx, itg NmbLin, double box[4],
                           double (*crd)[2], uint64_t (*idx)[2] )
{
   itg i;
   int NewTyp;
   double len = pow(2,62);
   ArgSct arg;

   // Get and check lib parallel instance
   if(!ParIdx)
      return(0);

   NewTyp = NewType(ParIdx, NmbLin);
   arg.crd2 = crd;
   arg.idx = idx;
   arg.box[0] = box[0];
   arg.box[1] = box[1];
   arg.box[2] = len / (box[2] - box[0]);
   arg.box[3] = len / (box[3] - box[1]);

   LaunchParallel(ParIdx, NewTyp, 0, (void *)RenPrc2D, (void *)&arg);
   ParallelQsort(ParIdx, &idx[1][0], NmbLin, 2 * sizeof(int64_t), CmpPrc);

   for(i=1;i<=NmbLin;i++)
      idx[ idx[i][1] ][0] = i;

   return(1);
}


/*----------------------------------------------------------------------------*/
/* Starts or stops the given timer                                            */
/*----------------------------------------------------------------------------*/

double GetWallClock()
{
#ifdef _WIN32
   struct __timeb64 tb;
   _ftime64(&tb);
   return((double)tb.time + (double)tb.millitm/1000.);
#else
   struct timeval tp;
   gettimeofday(&tp, NULL);
   return(tp.tv_sec + tp.tv_usec / 1000000.);
#endif
}


/*----------------------------------------------------------------------------*/
/* Encapsulate the selection between libMemBlock and regular libc malloc      */
/*----------------------------------------------------------------------------*/

static void *LPL_malloc(void *lmb, int64_t siz)
{
#ifdef WITH_LIBMEMBLOCKS
   return(LmbAlcPag((LmbSct *)lmb, siz, 0, LMB_ALLOC_DONT_CLEAR));
#else
  return(malloc(siz));
#endif
}


/*----------------------------------------------------------------------------*/
/* Encapsulate the selection between libMemBlock and regular libc calloc      */
/*----------------------------------------------------------------------------*/

static void *LPL_calloc(void *lmb, int64_t itm, int64_t siz)
{
#ifdef WITH_LIBMEMBLOCKS
   return(LmbAlcPag((LmbSct *)lmb, itm * siz, 0, LMB_ALLOC_AND_CLEAR));
#else
   return(calloc(itm, siz));
#endif
}


/*----------------------------------------------------------------------------*/
/* Encapsulate the selection between libMemBlock and regular libc free        */
/*----------------------------------------------------------------------------*/

static void LPL_free(void *lmb, void *adr)
{
#ifdef WITH_LIBMEMBLOCKS
   LmbRlsPag((LmbSct *)lmb, adr);
#else
   free(adr);
#endif
}


/*----------------------------------------------------------------------------*/
/* Duplication macros                                                        */
/*----------------------------------------------------------------------------*/

#define DUP(s,n) DUP ## n (s)
#define DUP1(s) s
#define DUP2(s) DUP1(s),s
#define DUP3(s) DUP2(s),s
#define DUP4(s) DUP3(s),s
#define DUP5(s) DUP4(s),s
#define DUP6(s) DUP5(s),s
#define DUP7(s) DUP6(s),s
#define DUP8(s) DUP7(s),s
#define DUP9(s) DUP8(s),s
#define DUP10(s) DUP9(s),s
#define DUP11(s) DUP10(s),s
#define DUP12(s) DUP11(s),s
#define DUP13(s) DUP12(s),s
#define DUP14(s) DUP13(s),s
#define DUP15(s) DUP14(s),s
#define DUP16(s) DUP15(s),s
#define DUP17(s) DUP16(s),s
#define DUP18(s) DUP17(s),s
#define DUP19(s) DUP18(s),s
#define DUP20(s) DUP19(s),s


#define ARG(a,n) ARG ## n (a)
#define ARG1(a) a[0]
#define ARG2(a) ARG1(a),a[1]
#define ARG3(a) ARG2(a),a[2]
#define ARG4(a) ARG3(a),a[3]
#define ARG5(a) ARG4(a),a[4]
#define ARG6(a) ARG5(a),a[5]
#define ARG7(a) ARG6(a),a[6]
#define ARG8(a) ARG7(a),a[7]
#define ARG9(a) ARG8(a),a[8]
#define ARG10(a) ARG9(a),a[9]
#define ARG11(a) ARG10(a),a[10]
#define ARG12(a) ARG11(a),a[11]
#define ARG13(a) ARG12(a),a[12]
#define ARG14(a) ARG13(a),a[13]
#define ARG15(a) ARG14(a),a[14]
#define ARG16(a) ARG15(a),a[15]
#define ARG17(a) ARG16(a),a[16]
#define ARG18(a) ARG17(a),a[17]
#define ARG19(a) ARG18(a),a[18]
#define ARG20(a) ARG19(a),a[19]


/*----------------------------------------------------------------------------*/
/* Call a C thread with 1 to 20 arguments                                     */
/*----------------------------------------------------------------------------*/

static void CalVarArgPrc(itg BegIdx, itg EndIdx, int PthIdx, ParSct *par)
{
   switch(par->NmbVarArg)
   {
      case 1 :
      {
         void (*prc1)(itg, itg, int, DUP(void *, 1)) =
            (void (*)(itg, itg, int, DUP(void *, 1)))par->prc;
         prc1(BegIdx, EndIdx, PthIdx, ARG(par->VarArgTab, 1));
      }break;

      case 2 :
      {
         void (*prc1)(itg, itg, int, DUP(void *, 2)) =
            (void (*)(itg, itg, int, DUP(void *, 2)))par->prc;
         prc1(BegIdx, EndIdx, PthIdx, ARG(par->VarArgTab, 2));
      }break;

      case 3 :
      {
         void (*prc1)(itg, itg, int, DUP(void *, 3)) =
            (void (*)(itg, itg, int, DUP(void *, 3)))par->prc;
         prc1(BegIdx, EndIdx, PthIdx, ARG(par->VarArgTab, 3));
      }break;

      case 4 :
      {
         void (*prc1)(itg, itg, int, DUP(void *, 4)) =
            (void (*)(itg, itg, int, DUP(void *, 4)))par->prc;
         prc1(BegIdx, EndIdx, PthIdx, ARG(par->VarArgTab, 4));
      }break;

      case 5 :
      {
         void (*prc1)(itg, itg, int, DUP(void *, 5)) =
            (void (*)(itg, itg, int, DUP(void *, 5)))par->prc;
         prc1(BegIdx, EndIdx, PthIdx, ARG(par->VarArgTab, 5));
      }break;

      case 6 :
      {
         void (*prc1)(itg, itg, int, DUP(void *, 6)) =
            (void (*)(itg, itg, int, DUP(void *, 6)))par->prc;
         prc1(BegIdx, EndIdx, PthIdx, ARG(par->VarArgTab, 6));
      }break;

      case 7 :
      {
         void (*prc1)(itg, itg, int, DUP(void *, 7)) =
            (void (*)(itg, itg, int, DUP(void *, 7)))par->prc;
         prc1(BegIdx, EndIdx, PthIdx, ARG(par->VarArgTab, 7));
      }break;

      case 8 :
      {
         void (*prc1)(itg, itg, int, DUP(void *, 8)) =
            (void (*)(itg, itg, int, DUP(void *, 8)))par->prc;
         prc1(BegIdx, EndIdx, PthIdx, ARG(par->VarArgTab, 8));
      }break;

      case 9 :
      {
         void (*prc1)(itg, itg, int, DUP(void *, 9)) =
            (void (*)(itg, itg, int, DUP(void *, 9)))par->prc;
         prc1(BegIdx, EndIdx, PthIdx, ARG(par->VarArgTab, 9));
      }break;

      case 10 :
      {
         void (*prc1)(itg, itg, int, DUP(void *, 10)) =
            (void (*)(itg, itg, int, DUP(void *, 10)))par->prc;
         prc1(BegIdx, EndIdx, PthIdx, ARG(par->VarArgTab, 10));
      }break;

      case 11 :
      {
         void (*prc1)(itg, itg, int, DUP(void *, 11)) =
            (void (*)(itg, itg, int, DUP(void *, 11)))par->prc;
         prc1(BegIdx, EndIdx, PthIdx, ARG(par->VarArgTab, 11));
      }break;

      case 12 :
      {
         void (*prc1)(itg, itg, int, DUP(void *, 12)) =
            (void (*)(itg, itg, int, DUP(void *, 12)))par->prc;
         prc1(BegIdx, EndIdx, PthIdx, ARG(par->VarArgTab, 12));
      }break;

      case 13 :
      {
         void (*prc1)(itg, itg, int, DUP(void *, 13)) =
            (void (*)(itg, itg, int, DUP(void *, 13)))par->prc;
         prc1(BegIdx, EndIdx, PthIdx, ARG(par->VarArgTab, 13));
      }break;

      case 14 :
      {
         void (*prc1)(itg, itg, int, DUP(void *, 14)) =
            (void (*)(itg, itg, int, DUP(void *, 14)))par->prc;
         prc1(BegIdx, EndIdx, PthIdx, ARG(par->VarArgTab, 14));
      }break;

      case 15 :
      {
         void (*prc1)(itg, itg, int, DUP(void *, 15)) =
            (void (*)(itg, itg, int, DUP(void *, 15)))par->prc;
         prc1(BegIdx, EndIdx, PthIdx, ARG(par->VarArgTab, 15));
      }break;

      case 16 :
      {
         void (*prc1)(itg, itg, int, DUP(void *, 16)) =
            (void (*)(itg, itg, int, DUP(void *, 16)))par->prc;
         prc1(BegIdx, EndIdx, PthIdx, ARG(par->VarArgTab, 16));
      }break;

      case 17 :
      {
         void (*prc1)(itg, itg, int, DUP(void *, 17)) =
            (void (*)(itg, itg, int, DUP(void *, 17)))par->prc;
         prc1(BegIdx, EndIdx, PthIdx, ARG(par->VarArgTab, 17));
      }break;

      case 18 :
      {
         void (*prc1)(itg, itg, int, DUP(void *, 18)) =
            (void (*)(itg, itg, int, DUP(void *, 18)))par->prc;
         prc1(BegIdx, EndIdx, PthIdx, ARG(par->VarArgTab, 18));
      }break;

      case 19 :
      {
         void (*prc1)(itg, itg, int, DUP(void *, 19)) =
            (void (*)(itg, itg, int, DUP(void *, 19)))par->prc;
         prc1(BegIdx, EndIdx, PthIdx, ARG(par->VarArgTab, 19));
      }break;

      case 20 :
      {
         void (*prc1)(itg, itg, int, DUP(void *, 20)) =
            (void (*)(itg, itg, int, DUP(void *, 20)))par->prc;
         prc1(BegIdx, EndIdx, PthIdx, ARG(par->VarArgTab, 20));
      }break;
   }
}


/*----------------------------------------------------------------------------*/
/* Call a C pipeline with 1 to 20 arguments                                   */
/*----------------------------------------------------------------------------*/

static void CalVarArgPip(PipSct *pip, void *prc)
{
   switch(pip->NmbVarArg)
   {
      case 1 :
      {
         void (*prc1)(DUP(void *, 1)) = (void (*)(DUP(void *, 1)))prc;
         prc1(ARG(pip->VarArgTab, 1));
      }break;

      case 2 :
      {
         void (*prc1)(DUP(void *, 2)) = (void (*)(DUP(void *, 2)))prc;
         prc1(ARG(pip->VarArgTab, 2));
      }break;

      case 3 :
      {
         void (*prc1)(DUP(void *, 3)) = (void (*)(DUP(void *, 3)))prc;
         prc1(ARG(pip->VarArgTab, 3));
      }break;

      case 4 :
      {
         void (*prc1)(DUP(void *, 4)) = (void (*)(DUP(void *, 4)))prc;
         prc1(ARG(pip->VarArgTab, 4));
      }break;

      case 5 :
      {
         void (*prc1)(DUP(void *, 5)) = (void (*)(DUP(void *, 5)))prc;
         prc1(ARG(pip->VarArgTab, 5));
      }break;

      case 6 :
      {
         void (*prc1)(DUP(void *, 6)) = (void (*)(DUP(void *, 6)))prc;
         prc1(ARG(pip->VarArgTab, 6));
      }break;

      case 7 :
      {
         void (*prc1)(DUP(void *, 7)) = (void (*)(DUP(void *, 7)))prc;
         prc1(ARG(pip->VarArgTab, 7));
      }break;

      case 8 :
      {
         void (*prc1)(DUP(void *, 8)) = (void (*)(DUP(void *, 8)))prc;
         prc1(ARG(pip->VarArgTab, 8));
      }break;

      case 9 :
      {
         void (*prc1)(DUP(void *, 9)) = (void (*)(DUP(void *, 9)))prc;
         prc1(ARG(pip->VarArgTab, 9));
      }break;

      case 10 :
      {
         void (*prc1)(DUP(void *, 10)) = (void (*)(DUP(void *, 10)))prc;
         prc1(ARG(pip->VarArgTab, 10));
      }break;

      case 11 :
      {
         void (*prc1)(DUP(void *, 11)) = (void (*)(DUP(void *, 11)))prc;
         prc1(ARG(pip->VarArgTab, 11));
      }break;

      case 12 :
      {
         void (*prc1)(DUP(void *, 12)) = (void (*)(DUP(void *, 12)))prc;
         prc1(ARG(pip->VarArgTab, 12));
      }break;

      case 13 :
      {
         void (*prc1)(DUP(void *, 13)) = (void (*)(DUP(void *, 13)))prc;
         prc1(ARG(pip->VarArgTab, 13));
      }break;

      case 14 :
      {
         void (*prc1)(DUP(void *, 14)) = (void (*)(DUP(void *, 14)))prc;
         prc1(ARG(pip->VarArgTab, 14));
      }break;

      case 15 :
      {
         void (*prc1)(DUP(void *, 15)) = (void (*)(DUP(void *, 15)))prc;
         prc1(ARG(pip->VarArgTab, 15));
      }break;

      case 16 :
      {
         void (*prc1)(DUP(void *, 16)) = (void (*)(DUP(void *, 16)))prc;
         prc1(ARG(pip->VarArgTab, 16));
      }break;

      case 17 :
      {
         void (*prc1)(DUP(void *, 17)) = (void (*)(DUP(void *, 17)))prc;
         prc1(ARG(pip->VarArgTab, 17));
      }break;

      case 18 :
      {
         void (*prc1)(DUP(void *, 18)) = (void (*)(DUP(void *, 18)))prc;
         prc1(ARG(pip->VarArgTab, 18));
      }break;

      case 19 :
      {
         void (*prc1)(DUP(void *, 19)) = (void (*)(DUP(void *, 19)))prc;
         prc1(ARG(pip->VarArgTab, 19));
      }break;

      case 20 :
      {
         void (*prc1)(DUP(void *, 20)) = (void (*)(DUP(void *, 20)))prc;
         prc1(ARG(pip->VarArgTab, 20));
      }break;
   }
}


/*----------------------------------------------------------------------------*/
/* Read, renumber through a Hilbert SFC and write the mesh                    */
/*----------------------------------------------------------------------------*/

LplSct *MeshRenumbering(int64_t ParIdx, int NmbGrn,
                        int RenTyp, int GmlMod, int dim, ...)
{
   int      i, j, t, siz, NmbCpu = 10, *TmpEle;
   int      RenEleTyp[ LplMax ], RenVerTyp;
   int64_t  LibParIdx;
   double   *TmpCrd;
   LplSct   *msh;
   va_list  VarArg;
   char     *SchStr[4] = {"Hilbert", "Z-curve", "random", "initial"};
#ifdef WITH_METIS
   int      grn, *NewGrn, *ColTab, CurCol, CurGrn;
   int      NmbCol, (*ColPar)[3], (*VerGrnPar)[4];
   ParSct   *par = (ParSct *)ParIdx;
#endif


   // ---------------------------------------------------------------
   // Allocate the main data structure, parse and check the arguments
   // ---------------------------------------------------------------

   // Check mandatory inputs
   if(dim != 3)
      return(NULL);

   // Allocate and setup the renumbering structure
   if(!(msh = calloc(1, sizeof(LplSct))))
      return(NULL);

   // Decode the variable list of arguments
   va_start(VarArg, dim);

   do
   {
      t = va_arg(VarArg, int);

      if(t == LplVer)
      {
         msh->VerTyp = va_arg(VarArg, int);
         msh->NmbVer = va_arg(VarArg, int);
         msh->CrdTab = va_arg(VarArg, double *);
         msh->VerRef = va_arg(VarArg, int *);
      }
      else if( (t >= LplEdg) && (t < LplMax) )
      {
         msh->EleTyp[t] = va_arg(VarArg, int);
         msh->NmbEle[t] = va_arg(VarArg, int);
         msh->EleTab[t] = va_arg(VarArg, int *);
         msh->EleRef[t] = va_arg(VarArg, int *);
      }
   }while(t != LplMax);

   va_end(VarArg);

   // We need at least some vertices
   if(!msh->NmbVer || !msh->CrdTab)
   {
      free(msh);
      return(NULL);
   }


   // -------------------------------------------------------------------------
   // Call metis partitioning, the colouring, then setup colors and grains data
   // -------------------------------------------------------------------------

#ifdef WITH_METIS
   if(NmbGrn)
   {
      // Activate the color-grain mode
      msh->NmbGrn = NmbGrn;
      msh->ColGrnMod = 1;

      // Allocate a color and a grain table for each mesh entities
      msh->VerGrn = calloc(msh->NmbVer + 1, sizeof(int));
      assert(msh->VerGrn);

      msh->VerCol = calloc(msh->NmbVer + 1, sizeof(int));
      assert(msh->VerCol);

      for(t=LplEdg; t<LplMax; t++)
      {
         if(!msh->NmbEle[t])
            continue;

         msh->EleCol[t] = calloc(msh->NmbEle[t] + 1, sizeof(int));
         assert(msh->EleCol[t]);

         msh->EleGrn[t] = calloc(msh->NmbEle[t] + 1, sizeof(int));
         assert(msh->EleGrn[t]);
      }

      // Allocate and setup rank 1 vertex balls
      SetV2VBal(msh);

      // Partition the mesh into grains with Metis
      MetisPartitioning(msh, NmbGrn);

      // Allocate and setup rank 2 vertex balls
      SetRk2Bal(msh);

      // Give a color to each grain
      ColTab = ColorPartImplicit(NmbCpu, NmbGrn, msh->NmbVer,
                                 msh->AdrBalRk1, msh->LstBalRk1,
                                 msh->AdrBalRk2, msh->LstBalRk2,
                                 msh->VerGrn, msh->VrbLvl);

      // Propagate the color-grain information to elements
      for(i=1;i<=msh->NmbVer;i++)
         msh->VerCol[i] = ColTab[ msh->VerGrn[i] ];

      for(t=LplEdg; t<LplMax; t++)
      {
         if(!msh->NmbEle[t])
            continue;

         for(i=1;i<=msh->NmbEle[t];i++)
         {
            msh->EleCol[t][i] = msh->VerCol[ msh->EleTab[t][ i * EleSiz[t] ] ];
            msh->EleGrn[t][i] = msh->VerGrn[ msh->EleTab[t][ i * EleSiz[t] ] ];
         }
      }

      // Look for the number a color and grains and define the required
      // number of bits needed to store these information
      msh->NmbCol = 0;

      for(i=1;i<=NmbGrn;i++)
         msh->NmbCol = MAX(msh->NmbCol, ColTab[i]);

      // Check that there are no 2nd order pinchments
      ChkColGrn(msh);

      // Set the bit field size, maks and left shift to store the color value
      msh->ColBit = ceil(log(msh->NmbCol) / log(2));
      msh->ColMsk = (1ULL << msh->ColBit) - 1ULL;
      msh->ColLft = 64 - msh->ColBit;

      // Set the bit field size, maks and left shift to store the grain value
      msh->GrnBit = ceil(log(msh->NmbGrn) / log(2));
      msh->GrnMsk = (1ULL << msh->GrnBit) - 1ULL;
      msh->GrnLft = 64 - msh->ColBit - msh->GrnBit;
   }
#endif


   // ----------------------------------------------------------------------
   // Prepare references and degree related data for the GMlib modified sort
   // ----------------------------------------------------------------------

   if(GmlMod == 1)
   {
      SetVerDeg(msh);

      // Only one bit is need to encode the high or low degree information
      msh->DegBit = 1;
      msh->DegMsk = 1ULL;
      msh->DegLft = 64 - msh->ColBit - msh->GrnBit - msh->DegBit;

      // Use the 8 leftmost bits to store face ref as boundary condition tag
      msh->RefBit = 8;
      msh->RefMsk = (1ULL << msh->RefBit) - 1ULL;
      msh->RefLft = 64 - msh->RefBit;
   }
   else if(GmlMod == 2)
   {
      SetMatSlc(msh);

      // three bits are needed to encode the 5 possible vertex degrees
      msh->DegBit = 3;
      msh->DegMsk = 7ULL;
      msh->RefMsk = (1ULL << msh->DegBit) - 1ULL;
      msh->DegLft = 64 - msh->ColBit - msh->GrnBit - msh->DegBit;
   }


   // ------------------------------------------------------
   // Finaly, set the Hilbert bit field size and right shift
   // ------------------------------------------------------

   // The hilbert code is right shifted with the number of bit
   // used by all other sorting keys
   msh->VerHilBit = 64 - msh->ColBit - msh->GrnBit - msh->DegBit;
   msh->VerHilRgt = msh->ColBit + msh->GrnBit + msh->DegBit;

   msh->FacHilBit = 64 - msh->ColBit - msh->GrnBit - msh->RefBit;
   msh->FacHilRgt = msh->ColBit + msh->GrnBit + msh->RefBit;

   msh->VolHilBit = 64 - msh->ColBit - msh->GrnBit;
   msh->VolHilRgt = msh->ColBit + msh->GrnBit;


   // --------------------
   // Vertices renumbering
   // --------------------

   if(msh->VrbLvl >= 1)
   {
      printf("Sorting keys table: number of bit per key for each dimension of mesh entities\n");
      printf(" Entity | rank4 (color) | rank3 (grain) | rank2 (degree or ref) | rank1 (%s)\n",
               SchStr[ msh->mod ]);
      printf(" Vertex |      %2d       |      %2d       |           %2d          |      %2d\n",
               msh->ColBit, msh->GrnBit, msh->DegBit, msh->VerHilBit);

      printf(" Face   |      %2d       |      %2d       |           %2d          |      %2d\n",
               msh->ColBit, msh->GrnBit, msh->RefBit, msh->FacHilBit);

      printf(" Volume |      %2d       |      %2d       |           %2d          |      %2d\n",
               msh->ColBit, msh->GrnBit, 0, msh->VolHilBit);
   }

   msh->VerCod = malloc( (msh->NmbVer + 1) * 2 * sizeof(int64_t) );
   assert(msh->VerCod);

   LibParIdx = InitParallel(NmbCpu);

   RenVerTyp = NewType(LibParIdx, msh->NmbVer);

   SetBndBox(msh);

   LaunchParallel(LibParIdx, RenVerTyp, 0, (void *)RenVer, (void *)msh);
   ParallelQsort(LibParIdx, msh->VerCod[1], msh->NmbVer, 2 * sizeof(int64_t), CmpFnc);

   msh->Old2New = malloc( (size_t)(msh->NmbVer+1) * sizeof(int) );
   assert(msh->Old2New);

   for(i=1;i<=msh->NmbVer;i++)
      msh->Old2New[ msh->VerCod[i][0] ] = i;

   TmpCrd = malloc( (msh->NmbVer + 1) * 3 * sizeof(double) );
   assert(TmpCrd);

   for(i=1;i<=msh->NmbVer;i++)
      for(j=0;j<3;j++)
         TmpCrd[ i * 3 + j ] = msh->CrdTab[ msh->VerCod[i][0] * 3 + j ];

   memcpy(msh->CrdTab, TmpCrd, (msh->NmbVer+1) * 3 * sizeof(double));
   free(TmpCrd);

#ifdef WITH_METIS
   if(NmbGrn)
   {
      // As the vertices renumbered, their colors and grains tables must be rebuilt
      NewGrn = malloc( (msh->NmbVer + 1) * sizeof(int) );
      assert(NewGrn);

      for(i=1;i<=msh->NmbVer;i++)
         NewGrn[i] = msh->VerGrn[ msh->VerCod[i][0] ];

      free(msh->VerGrn);
      msh->VerGrn = NewGrn;

      for(i=1;i<=msh->NmbVer;i++)
         msh->VerCol[i] = ColTab[ msh->VerGrn[i] ];
   }
#endif


   // --------------------
   // Elements renumbering
   // --------------------

   for(t=LplEdg; t<LplMax; t++)
   {
      if(!msh->NmbEle[t])
         continue;

      // Renumber all kinds of elements
      RenEleTyp[t] = NewType(LibParIdx, msh->NmbEle[t]);
      msh->TypIdx = t;
      msh->EleCod[t] = malloc( (msh->NmbEle[t] + 1) * 2 * sizeof(int64_t) );
      assert(msh->EleCod[t]);
      LaunchParallel(LibParIdx, RenEleTyp[t], 0, (void *)RenEle, (void *)msh);
      ParallelQsort(LibParIdx, msh->EleCod[t][1], msh->NmbEle[t], 2 * sizeof(int64_t), CmpFnc);

      siz = EleSiz[t];

      TmpEle = malloc( (msh->NmbEle[t] + 1) * siz * sizeof(int) );
      assert(TmpEle);

      for(i=1;i<=msh->NmbEle[t];i++)
         for(j=0;j<siz;j++)
            TmpEle[ i * siz + j ] = msh->EleTab[t][ msh->EleCod[t][i][0] * siz + j ];

      memcpy(msh->EleTab[t], TmpEle, (msh->NmbEle[t] + 1) * siz * sizeof(int));
      free(TmpEle);

#ifdef WITH_METIS
      if(NmbGrn)
      {
         // And rebuild their color and grain tables
         for(i=1;i<=msh->NmbEle[t];i++)
         {
            msh->EleCol[t][i] = msh->VerCol[ msh->EleTab[t][ i * siz ] ];
            msh->EleGrn[t][i] = msh->VerGrn[ msh->EleTab[t][ i * siz ] ];
         }
      }
#endif
   }

   StopParallel(LibParIdx);


   // -----------------------------------------------------------------
   // Color-grain mode: setup global colouring and each entities grains
   // -----------------------------------------------------------------

#ifdef WITH_METIS
   if(NmbGrn)
   {
      CurCol = msh->VerCol[1];
      NmbCol = 1;
      CurGrn = msh->VerGrn[1];
      NmbGrn = 1;

      // Build the vertex color-grains partitions
      for(i=1;i<=msh->NmbVer;i++)
      {
         if(msh->VerCol[i] == CurCol)
            msh->VerCol[i] = NmbCol;
         else
         {
            CurCol = msh->VerCol[i];
            msh->VerCol[i] = ++NmbCol;
         }

         if(msh->VerGrn[i] == CurGrn)
            msh->VerGrn[i] = NmbGrn;
         else
         {
            CurGrn = msh->VerGrn[i];
            msh->VerGrn[i] = ++NmbGrn;
         }
      }

      // Derive element's colors and grains from the vertices
      for(t=LplEdg; t<LplMax; t++)
      {
         for(i=1;i<=msh->NmbEle[t];i++)
         {
            msh->EleCol[t][i] = msh->VerCol[ msh->EleTab[t][ i * EleSiz[t] ] ];
            msh->EleGrn[t][i] = msh->VerGrn[ msh->EleTab[t][ i * EleSiz[t] ] ];
         }
      }

      // Allocate a single colour table for the whole mesh as colours
      // are based on vertices only
      ColPar = malloc( (msh->NmbCol + 1) * 3 * sizeof(int) );
      assert(ColPar);

      // Each kind of entity needs a dedicated grain table to store
      // the begin and ending indices. Some grains may be empty
      VerGrnPar = malloc( (msh->NmbGrn + 1) * 4 * sizeof(int) );
      assert(VerGrnPar);

      for(t=LplEdg; t<LplMax; t++)
         if(msh->NmbEle[t])
         {
            msh->EleGrnPar[t] = malloc( (msh->NmbGrn + 1) * 2 * sizeof(int) );
            assert(msh->EleGrnPar[t]);
         }

      // Setup vertex colours and grains partitions
      CurGrn = msh->VerGrn[1];
      NmbGrn = 1;
      VerGrnPar[ NmbGrn ][0] = 1;
      VerGrnPar[ NmbGrn ][2] = msh->VerCol[1];
      VerGrnPar[ NmbGrn ][3] = msh->VerGrn[1];

      CurCol = msh->VerCol[1];
      NmbCol = 1;
      ColPar[ NmbCol ][0] = 1;
      ColPar[ NmbCol ][2] = NmbGrn;

      for(i=1;i<msh->NmbVer;i++)
      {
         if(msh->VerGrn[i] != CurGrn)
         {
            VerGrnPar[ NmbGrn ][1] = i - 1;
            NmbGrn++;
            VerGrnPar[ NmbGrn ][0] = i;
            CurGrn = msh->VerGrn[i];
            VerGrnPar[ NmbGrn ][2] = msh->VerCol[i];
            VerGrnPar[ NmbGrn ][3] = msh->VerGrn[i];

            if(msh->VerCol[i] != CurCol)
            {
               ColPar[ NmbCol ][1] = NmbGrn - 1;
               NmbCol++;
               ColPar[ NmbCol ][0] = NmbGrn;
               CurCol = msh->VerCol[i];
               ColPar[ NmbCol ][2] = msh->VerCol[i];
            }
         }
      }

      VerGrnPar[ NmbGrn ][1] = msh->NmbVer;
      ColPar[ NmbCol ][1] = NmbGrn;

      if(msh->VrbLvl >= 1)
      {
         for(i=1;i<=NmbCol;i++)
            printf(  "vertex color %2d (%2d): %6d -> %6d, size: %6d\n",
                     i, ColPar[i][2], ColPar[i][0], ColPar[i][1],
                     ColPar[i][1] - ColPar[i][0] + 1);

         for(i=1;i<=NmbGrn;i++)
            printf(  "vertex grain %6d (%6d/%6d): %10d -> %10d, size: %10d\n",
                     i, VerGrnPar[i][2], VerGrnPar[i][3],
                     VerGrnPar[i][0], VerGrnPar[i][1],
                     VerGrnPar[i][1] - VerGrnPar[i][0] + 1);
      }

      msh->ColPar = malloc((NmbCol + 1) * 2 * sizeof(int));
      assert(msh->ColPar);

      msh->VerGrnPar = malloc((NmbGrn + 1) * 2 * sizeof(int));
      assert(msh->VerGrnPar);

      for(i=1;i<=NmbGrn;i++)
         for(j=0;j<2;j++)
            msh->VerGrnPar[i][j] = VerGrnPar[i][j];

      for(i=1;i<=NmbCol;i++)
         for(j=0;j<2;j++)
            msh->ColPar[i][j] = ColPar[i][j];

      free(ColPar);
      free(VerGrnPar);

      // Derive elements' grain and color from the vertices
      for(t=LplEdg; t<LplMax; t++)
      {
         if(!msh->NmbEle[t])
            continue;

         for(i=0;i<=NmbGrn;i++)
            msh->EleGrnPar[t][i][0] = msh->EleGrnPar[t][i][1] = 0;

         for(i=1;i<=msh->NmbEle[t];i++)
         {
            grn = msh->EleGrn[t][i];

            if(!msh->EleGrnPar[t][ grn ][0])
            {
               msh->EleGrnPar[t][ grn ][0] = i;
               msh->EleGrnPar[t][ grn ][1] = i;
            }
            else
            {
               msh->EleGrnPar[t][ grn ][0] = MIN(msh->EleGrnPar[t][ grn ][0], i);
               msh->EleGrnPar[t][ grn ][1] = MAX(msh->EleGrnPar[t][ grn ][1], i);
            }
         }
      }

      // Check the colouring validity again after the renumbering
      ChkColGrn(msh);

      // Send colors to grains information to the LPlib
      SetColors(ParIdx, msh->NmbCol, msh->ColPar, msh->NmbGrn);

      // Send grains of vertices to the LPlib
      SetGrains(ParIdx, LplVer, msh->VerTyp, msh->VerGrnPar);

      // Send grains of edges to the LPlib
      SetGrains(ParIdx, LplEdg, msh->EleTyp[ LplEdg ], msh->EleGrnPar[ LplEdg ]);

      // Send grains of tets to the LPlib
      SetGrains(ParIdx, LplTet, msh->EleTyp[ LplTet ], msh->EleGrnPar[ LplTet ]);

      // Rebuild the rank 1 & 2 balls of vertices after the renumbering
      free(msh->VerDeg);
      free(msh->AdrBalRk1);
      free(msh->LstBalRk1);
      SetV2VBal(msh);

      free(msh->AdrBalRk2);
      free(msh->LstBalRk2);
      SetRk2Bal(msh);

      // Build the dependencies between grains for the dynamic scheduling
      SetGrnDep(par, msh);

      // Allocate and setup the grain wroks for use by the Lplib's scheduler
      par->GrnWrkTab = calloc(par->NmbGrn , sizeof(WrkSct));
      assert(par->GrnWrkTab);

      par->NmbGrnWrk = msh->NmbGrn;
      par->NmbCol = msh->NmbCol;

      par->ColCpt = calloc(par->NmbCol + 1, sizeof(int));
      assert(par->ColCpt);

      par->RunDepTab = calloc(par->NmbGrn / 32 + 1 , sizeof(int));
      assert(par->RunDepTab);

      for(i=0;i<par->NmbGrn;i++)
         par->GrnWrkTab[i].DepWrdTab = &par->GrnMat[ (i+1) * par->NmbDepWrd ];
   }

   // Free all working tables
   free(msh->VerGrn);
   free(msh->VerCol);
   free(msh->VerCod);
   free(msh->Old2New);

   for(t=LplEdg; t<LplMax; t++)
      if(msh->NmbEle[t])
      {
         free(msh->EleCol[t]);
         free(msh->EleGrn[t]);
         free(msh->EleCod[t]);
         free(msh->EleGrnPar[t]);
         
      }
#endif

   // Retun the Lpl mesh datastructure that keeps track of the previous numbering
   return(msh);
}


/*----------------------------------------------------------------------------*/
/* Free all elements' numbering tables and the global structure itself        */
/*----------------------------------------------------------------------------*/

void FreeNumberingStruct(LplSct *ren)
{
   int t;

   for(t=0; t<LplMax; t++)
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

int RestoreNumbering(LplSct *ren, int NmbVer, double *CrdTab, ... )
{
   return(0);
}


/*----------------------------------------------------------------------------*/
/* Get an entity's index in the old numbering from its current one            */
/*----------------------------------------------------------------------------*/

int GetOldIndex(LplSct *ren, int t, int idx)
{
   return(ren->RenTab[t][ idx ][0]);
}


/*----------------------------------------------------------------------------*/
/* Get an entity's index in the current numbering from its old one            */
/*----------------------------------------------------------------------------*/

int GetNewIndex(LplSct *ren, int t, int idx)
{
   return(ren->RenTab[t][ idx ][1]);
}


/*----------------------------------------------------------------------------*/
/* Compute a mesh's bounding box                                              */
/*----------------------------------------------------------------------------*/

static void SetBndBox(LplSct *msh)
{
   int i, j;

   msh->box[0] = msh->box[3] = msh->CrdTab[3];
   msh->box[1] = msh->box[4] = msh->CrdTab[4];
   msh->box[2] = msh->box[5] = msh->CrdTab[5];

   for(i=1;i<=msh->NmbVer;i++)
      for(j=0;j<3;j++)
      {
         msh->box[j  ] = MIN(msh->box[j  ], msh->CrdTab[ i*3 + j ]);
         msh->box[j+3] = MAX(msh->box[j+3], msh->CrdTab[ i*3 + j ]);
      }

   // normalize the bounding box to map the geometry on a 64 bit cube
   for(j=0;j<3;j++)
      msh->box[j+3] = pow(2,64) / (msh->box[j+3] - msh->box[j]);
}


/*----------------------------------------------------------------------------*/
/* Comparison of two items for the qsort                                      */
/*----------------------------------------------------------------------------*/

static int CmpFnc(const void *a, const void *b)
{
   uint64_t *pa = (uint64_t *)a, *pb = (uint64_t *)b;

   if(pa[1] > pb[1])
      return(1);
   else
      return(-1);
}


/*----------------------------------------------------------------------------*/
/* Compute the hilbert code from 3d coordinates                               */
/*----------------------------------------------------------------------------*/

static uint64_t GetHilCod(double crd[3], double box[6], int itr, int mod)
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

static void RenVer(int BegIdx, int EndIdx, int PthIdx, LplSct *msh)
{
   int i;
   uint64_t ColCod = 0, GrnCod = 0, DegCod = 0, SchCod = 0;

   // Set the vertex code with the hilbert code, shifted by 11 bits in the case of
   // GMlib mode to make room for the high degree bit and the 10 bit ref hash tag
   for(i=BegIdx; i<=EndIdx; i++)
   {
      if(msh->ColGrnMod)
      {
         ColCod = (msh->VerCol[i] & msh->ColMsk) << msh->ColLft;
         GrnCod = (msh->VerGrn[i] & msh->GrnMsk) << msh->GrnLft;
      }

      if(msh->GmlMod)
         DegCod = (msh->VerDeg[i] & msh->DegMsk) << msh->DegLft;

      if(msh->mod == IniMod)
         SchCod = i;
      else
         SchCod = GetHilCod(&msh->CrdTab[ i*3 ], msh->box, MAXITR, msh->mod);

      SchCod = SchCod >> msh->VerHilRgt;

      msh->VerCod[i][0] = i;
      msh->VerCod[i][1] = ColCod | GrnCod | DegCod | SchCod;
   }
}


/*----------------------------------------------------------------------------*/
/* Compute the barycenter of any kind of element                              */
/*----------------------------------------------------------------------------*/

static void SetMidCrd(int NmbVer, int *EleTab, LplSct *msh, double *crd)
{
   int         i, j, MaxEdg;
   double      *TetCrd[4], len, MinLen, MaxLen;

   // Default treatment is to compute the element's barycenter
   for(j=0;j<3;j++)
      crd[j] = 0.;

   for(i=0;i<NmbVer;i++)
      for(j=0;j<3;j++)
         crd[j] += msh->CrdTab[ EleTab[i] * 3 + j ];

   for(j=0;j<3;j++)
      crd[j] /= NmbVer;


   return;
   // Special handling of tetrahedra: if the element is anisotropic,
   // use the longest edge midpoint instedad of the tet's barycenter
   // as it better represents an elongated element
   if(msh->TypIdx == LplTet)
   {
      MinLen = FLT_MAX;
      MaxEdg = -1;

      for(i=0;i<4;i++)
         TetCrd[i] = &msh->CrdTab[ EleTab[ i * NmbVer ] * 3 ];

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
/* Parallel loop renumbering any kind of element                              */
/*----------------------------------------------------------------------------*/

static void RenEle(int BegIdx, int EndIdx, int PthIdx, LplSct *msh)
{
   char     *BytPtr;
   int      i, j, ref, t = msh->TypIdx, siz = EleSiz[t];
   double   crd[3];
   uint64_t ColCod = 0, GrnCod = 0, RefCod = 0, SchCod = 0;

   // Set the elements SFC code from their barycenter
   for(i=BegIdx; i<=EndIdx; i++)
   {
      for(j=0;j<siz;j++)
         msh->EleTab[t][ i * siz + j ] =
            msh->Old2New[ msh->EleTab[t][ i * siz + j ] ];

      if(msh->ColGrnMod)
      {
         ColCod = (msh->EleCol[t][i] & msh->ColMsk) << msh->ColLft;
         GrnCod = (msh->EleGrn[t][i] & msh->GrnMsk) << msh->GrnLft;
      }

      if(msh->GmlMod && ((t == LplTri || (t == LplQad))))
      {
         BytPtr = (char *)&msh->EleRef[t][i];
         ref = BytPtr[0] + BytPtr[1] + BytPtr[2] + BytPtr[3];
         RefCod = (ref & msh->RefMsk) << msh->RefLft;
      }
      else
         RefCod = 0;

      if(msh->mod == IniMod)
         SchCod = i;
      else
      {
         SetMidCrd(siz, &msh->EleTab[t][ i * siz ], msh, crd);
         SchCod = GetHilCod(crd, msh->box, MAXITR, msh->mod);
      }

      if((t == LplTri) || (t == LplQad))
         SchCod = SchCod >> msh->FacHilRgt;
      else if((t >= LplTet) && (t <= LplHex))
         SchCod = SchCod >> msh->VolHilRgt;
      else
         SchCod = 0;

      msh->EleCod[t][i][0] = i;
      msh->EleCod[t][i][1] = ColCod | GrnCod | RefCod | SchCod;
   }
}


/*----------------------------------------------------------------------------*/
/* Set the code 10 upper bit with a hash key based on element's ref           */
/*----------------------------------------------------------------------------*/

static void SetVerDeg(LplSct *msh)
{
   int i, j, t, (*DegTab)[ LplMax ];
 
   // Allocate a vertex degree table with one scalar per kind of element
   DegTab = calloc( (size_t)(msh->NmbVer + 1), LplMax * sizeof(int));
   assert(DegTab);

   // Add each element's vertices to the degree associated to its kind
   for(t=0; t<LplMax; t++)
      for(i=1;i<=msh->NmbEle[t];i++)
         for(j=0;j<EleSiz[t];j++)
            DegTab[ msh->EleTab[t][i * EleSiz[t] + j ] ][t]++;

   // Set the high or low degree flag
   for(i=1;i<=msh->NmbVer;i++)
   {
      if( (DegTab[i][1] > MaxDeg[1][0])
      ||  (DegTab[i][2] > MaxDeg[2][0])
      ||  (DegTab[i][3] > MaxDeg[3][0])
      ||  (DegTab[i][6] > MaxDeg[6][0]) )
      {
         msh->VerDeg[i] = 1;
         msh->HghDeg++;

         for(j=0;j<LplMax;j++)
            msh->MaxDeg[j] = MAX(msh->MaxDeg[j], DegTab[i][j]);
      }
      else
         msh->VerDeg[i] = 0;

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
   for(j=0;j<LplMax;j++)
      if(msh->MaxDeg[j])
         msh->DegVec[j] = pow(2., ceil(log2(msh->MaxDeg[j])));

   // Print statistics about vertex connectivity
   printf(  "High-connected vertices      : %3.6f%%\n",
            (100. * (float)msh->HghDeg) / (float)msh->NmbVer );

   printf(  "Over-connected vertices      : %3.6f%%\n",
            (100. * (float)msh->OvrDeg) / (float)msh->NmbVer );

   for(j=0;j<LplMax;j++)
      if(msh->MaxDeg[j])
         printf(  "Ball of %s     : max deg = %3d, vec size = %3d\n",
                  EleNam[j], msh->MaxDeg[j],  msh->DegVec[j] );

   puts("");
}


/*----------------------------------------------------------------------------*/
/* Set the code 10 upper bit with a hash key based on element's ref           */
/*----------------------------------------------------------------------------*/

static void SetMatSlc(LplSct *msh)
{
   int      i, NmbEdg, *DegTab, VecCnt[6] = {0};
   uint64_t cod, DegTot = 0, VecTot = 0;

   // Allocate a degree table with one scalar per kind of element
   DegTab = calloc( (size_t)(msh->NmbVer + 1), sizeof(int));
   assert(DegTab);

   NmbEdg = SetDeg(msh, DegTab);
   printf("Unique edges extracted       : %d\n", NmbEdg);

   // Set the vertex degree with its slice vector size code: 1 -> 5
   for(i=1;i<=msh->NmbVer;i++)
   {
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

      msh->VerDeg[i] = cod;
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

static int SetDeg(LplSct *msh, int *DegTab)
{
   int i, j, idx0, idx1, key, MinIdx, MaxIdx, siz, col, *tet, NmbEdg = 0;
   HshSct *hsh, *buc;

   // Allocate the hash table to store all edges
   col = siz = msh->NmbEle[ LplTet ];
   hsh = malloc( 7LL * (size_t)siz * sizeof(HshSct));
   assert(hsh);

   // Clear the hash table's direct entries,
   // there is no need to clear the collision entries
   memset(hsh, 0, siz * sizeof(HshSct));

   // Loop over each tet and each tet's edges
   for(i=1;i<=siz;i++)
   {
      tet = &msh->EleTab[ LplTet ][i * EleSiz[ LplTet ] ];

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


#ifdef WITH_METIS


/*----------------------------------------------------------------------------*/
/* Provide the LPlib with color and grain index for each entities             */
/*----------------------------------------------------------------------------*/

static int SetColors(int64_t ParIdx, int NmbCol, int (*ColTab)[2], int NmbGrn)
{
   ParSct *par = (ParSct *)ParIdx;

   // Allocate and copy memory to store colors information
   par->NmbCol = NmbCol;
   par->NmbGrn = NmbGrn;
   par->ColTab = LPL_malloc(par->lmb, (NmbCol+1) * 2 * sizeof(int));

   if(!par->ColTab)
      return(1);

   memcpy(par->ColTab, ColTab, (NmbCol+1) * 2 * sizeof(int));

   return(0);
}


/*----------------------------------------------------------------------------*/
/* Provide the LPlib with color and grain index for each entities             */
/*----------------------------------------------------------------------------*/

static int SetGrains(int64_t ParIdx, int MshTyp, int TypIdx, int (*GrnTab)[2])
{
   ParSct *par = (ParSct *)ParIdx;

   par->TypIdx[ MshTyp ] = TypIdx;

   // Allocate and copy memory to store grains information
   par->GrnTab[ MshTyp ] = LPL_malloc(par->lmb, (par->NmbGrn+1) * 2 * sizeof(int));

   if(!par->GrnTab[ MshTyp ])
      return(1);

   memcpy(par->GrnTab[ MshTyp ], GrnTab, (par->NmbGrn+1) * 2 * sizeof(int));

   return(0);
}


/*----------------------------------------------------------------------------*/
/* Build the grain depedencies table to enable colorgrains dynamic scheduling */
/*----------------------------------------------------------------------------*/

static void SetGrnDep(ParSct *par, LplSct *msh)
{
   int i, j, *VerGrn, v1, v2, g1, g2, c1, c2;

   // Build the grains forward compatibility matrix
   par->NmbDepWrd = (msh->NmbGrn + 1) / 32 + 1;
   par->GrnMat = LPL_calloc(  par->lmb, msh->NmbGrn + 1,
                              par->NmbDepWrd * sizeof(int) );

   // Set all next grain entries to 0
   par->GrnCol = LPL_calloc(par->lmb, msh->NmbGrn + 1, sizeof(int) );

   // Set the identity
   for(i=1;i<=msh->NmbCol;i++)
      for(j=msh->ColPar[i][0]; j<=msh->ColPar[i][1]; j++)
         par->GrnCol[j] = i;

   VerGrn = malloc( (msh->NmbVer + 1) * sizeof(int) );
   assert(VerGrn);

   // Set grain bits with dependencies
   for(i=1;i<=msh->NmbGrn;i++)
      for(j=msh->VerGrnPar[i][0]; j<=msh->VerGrnPar[i][1]; j++)
         VerGrn[j] = i;

   for(i=1;i<=msh->NmbVer;i++)
   {
      v1 = i;
      g1 = VerGrn[ v1 ];
      c1 = par->GrnCol[ g1 ];

      for(j=msh->AdrBalRk2[i]; j<msh->AdrBalRk2[i+1]; j++)
      {
         v2 = msh->LstBalRk2[j];
         g2 = VerGrn[ v2 ];
         c2 = par->GrnCol[ g2 ];

         if(c2 == c1 - 1)
            SetBit(&par->GrnMat[ g1 * par->NmbDepWrd ], g2);
      }
   }

   free(VerGrn);
}


/*----------------------------------------------------------------------------*/
/* Look for possible vertex pinch between different grains from the same color*/
/*----------------------------------------------------------------------------*/

static void ChkColGrn(LplSct *msh)
{
   int i, j, k, VerIdx, NmbCol = 0, *VerColTab, *VerGrnTab;

   VerColTab = calloc(msh->NmbVer + 1, sizeof(int));
   assert(VerColTab);

   VerGrnTab = calloc(msh->NmbVer + 1, sizeof(int));
   assert(VerGrnTab);

   if(msh->VrbLvl >= 1)
      printf("Check color-grain consistency: ");

   for(i=1;i<=msh->NmbCol;i++)
   {
      for(j=1;j<=msh->NmbEle[ LplTet ];j++)
      {
         if(msh->EleCol[ LplTet ][j] != i)
            continue;

         for(k=0;k<4;k++)
         {
            VerIdx = msh->EleTab[ LplTet ][ j * EleSiz[ LplTet ] + k ];

            if(i > VerColTab[ VerIdx ])
            {
               VerColTab[ VerIdx ] = i;
               VerGrnTab[ VerIdx ] = msh->EleGrn[ LplTet ][j];
            }
            else if( (VerColTab[ VerIdx ] == i)
                  && (VerGrnTab[ VerIdx ] != msh->EleGrn[ LplTet ][j]) )
            {
               if(msh->VrbLvl >= 1)
               {
                  printf("vertex %d / tet %d / grain %d / color %d, collide with grain %d\n",
                           VerIdx, j, msh->EleGrn[ LplTet ][j], i, VerGrnTab[ VerIdx ]);
               }

               NmbCol++;
            }
         }
      }
   }

   if(msh->VrbLvl >= 1)
   {
      if(!NmbCol)
         puts("OK");
      else
         printf("%d collisions\n", NmbCol);
   }

   free(VerColTab);
   free(VerGrnTab);
}


/*----------------------------------------------------------------------------*/
/* Build the dual graph from a mesh and call Metis to get the grains          */
/*----------------------------------------------------------------------------*/

static int MetisPartitioning(LplSct *msh, int NmbPar)
{
   int i;
   MtsSct mts;

   BuildMetisGraph(msh, &mts);
   mts.nparts = NmbPar;

   if(METIS_PartGraphKway( &mts.nvtxs, &mts.ncon, mts.xadj, mts.adjncy,
                           NULL, NULL, NULL, &mts.nparts, NULL, NULL,
                           mts.options, &mts.objval, mts.part ) != METIS_OK)
   {
      return(0);
   }

   for(i=0;i<msh->NmbVer;i++)
      msh->VerGrn[i+1] = mts.part[i]+1;

   return(1);
}


/*----------------------------------------------------------------------------*/
/* Extract a metis graph from a tet mesh                                      */
/*----------------------------------------------------------------------------*/

static void BuildMetisGraph(LplSct *msh, MtsSct *mts)
{
   int      i, j;
   int64_t  TotDeg = 0;

   METIS_SetDefaultOptions(mts->options);

   mts->options[ METIS_OPTION_CTYPE ]     = METIS_CTYPE_RM;
   mts->options[ METIS_OPTION_CTYPE ]     = 0;
   mts->options[ METIS_OPTION_CONTIG ]    = 1;
   mts->options[ METIS_OPTION_OBJTYPE ]   = METIS_OBJTYPE_VOL;
   mts->options[ METIS_OPTION_IPTYPE ]    = METIS_IPTYPE_RANDOM;

   mts->nvtxs  = msh->NmbVer;
   mts->nedges = msh->TotDeg;
   mts->ncon   = 1;
   mts->VerDeg = malloc( msh->NmbVer       * sizeof(idx_t) );
   mts->xadj   = malloc( (msh->NmbVer + 1) * sizeof(idx_t) );
   mts->part   = malloc( msh->NmbVer       * sizeof(idx_t) );
   mts->adjncy = malloc( msh->TotDeg       * sizeof(idx_t) );

   if(!mts->xadj || !mts->VerDeg || ! mts->adjncy || !mts->part)
      exit(1);

   for(i=0;i<msh->NmbVer;i++)
   {
      mts->xadj[i] = TotDeg;
      mts->VerDeg[i] = msh->AdrBalRk1[i+2] - msh->AdrBalRk1[i+1];

      for(j=msh->AdrBalRk1[i+1]; j<msh->AdrBalRk1[i+2]; j++)
         mts->adjncy[j] = msh->LstBalRk1[j] - 1;

      TotDeg += mts->VerDeg[i];
   }

   mts->xadj[ msh->NmbVer ] = TotDeg;
}


/*----------------------------------------------------------------------------*/
/* Setup the vertex to vertex balls table compatible with colorgrains         */
/*----------------------------------------------------------------------------*/

static void SetV2VBal(LplSct *msh)
{
   int      i, j, a, b;
   int64_t  TabSiz = 0;

   msh->AdrBalRk1 = calloc(msh->NmbVer + 2, sizeof(int64_t));
   assert(msh->AdrBalRk1);

   msh->VerDeg = calloc(msh->NmbVer + 1, sizeof(int));
   assert(msh->VerDeg);

   // Set vertices degree
   for(i=1;i<=msh->NmbEle[ LplEdg ];i++)
      for(j=0;j<2;j++)
         msh->VerDeg[ msh->EleTab[ LplEdg ][ i * 2 + j ] ]++;

   for(i=1;i<=msh->NmbVer;i++)
   {
      msh->AdrBalRk1[i] = TabSiz;
      TabSiz += msh->VerDeg[i];
      msh->VerDeg[i] = 0;
   }

   msh->AdrBalRk1[ msh->NmbVer + 1 ] = TabSiz;

   // Allocate the global balls table and give each vertex a pointer to its own table
   msh->LstBalRk1 = malloc(TabSiz * sizeof(int));
   assert(msh->LstBalRk1);

   msh->TotDeg = TabSiz;

   if(!msh->AdrBalRk1 || !msh->LstBalRk1)
   {
      puts("Failed to allocate memory");
      exit(1);
   }

   // Fill the ball tables with elements type and index
   for(i=1;i<=msh->NmbEle[ LplEdg ];i++)
   {
      a = msh->EleTab[ LplEdg ][ i * 2     ];
      b = msh->EleTab[ LplEdg ][ i * 2 + 1 ];
      msh->LstBalRk1[ msh->AdrBalRk1[a] + msh->VerDeg[a]++ ] = b;
      msh->LstBalRk1[ msh->AdrBalRk1[b] + msh->VerDeg[b]++ ] = a;
   }
}


/*----------------------------------------------------------------------------*/
/* Set vertex to vertex rank 2 balls                                          */
/*----------------------------------------------------------------------------*/

static void SetRk2Bal(LplSct *msh)
{
   int *MrkTab;
   uint64_t i, j, k, TotDeg = 0, siz = 100, TabSiz = siz * msh->NmbVer;

   // Allocate a mark table and an address table
   MrkTab = calloc(msh->NmbVer + 1, sizeof(int));
   assert(MrkTab);

   msh->AdrBalRk2 = malloc((msh->NmbVer + 2) * sizeof(uint64_t));
   assert(msh->AdrBalRk2);

   msh->LstBalRk2 = malloc(TabSiz * sizeof(int));
   assert(msh->LstBalRk2);

   // Count the total rank 2 vertex degree and fill the address table
   for(i=1;i<=msh->NmbVer;i++)
   {
      msh->AdrBalRk2[i] = TotDeg;

      for(j=msh->AdrBalRk1[i]; j<msh->AdrBalRk1[ i+1 ]; j++)
      {
         for(k=msh->AdrBalRk1[ msh->LstBalRk1[j] ]; k<msh->AdrBalRk1[ msh->LstBalRk1[j] + 1 ]; k++)
         {
            if( (MrkTab[ msh->LstBalRk1[k] ] != i) && (msh->LstBalRk1[k] != i) )
            {
               if(TotDeg >= TabSiz)
               {
                  TabSiz *= 2;
                  msh->LstBalRk2 = realloc(msh->LstBalRk2, TabSiz * sizeof(int));
                  assert(msh->LstBalRk2);
                  puts("WARNING: rank 2 ball table realloc");
               }

               MrkTab[ msh->LstBalRk1[k] ] = i;
               msh->LstBalRk2[ TotDeg ] = msh->LstBalRk1[k];
               TotDeg++;
            }
         }
      }
   }

   // Set boundary elements
   msh->AdrBalRk2[0] = 0;
   msh->AdrBalRk2[ msh->NmbVer+1 ] = TotDeg;

   free(MrkTab);
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/

static int1d Wlf2_CheckIfVertexIsInTheList(int1d iVer, int1d nbv, int1d *verLst)
{
   int1d k;

   for (k = 0; k < nbv; k++)
   {
      if (verLst[k] == iVer) return k;
   }

   return -1;
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/

static void AddHeapListCol(int1d iVer, int1d *NbrHepLst, int1d *HepLst, int1d *Pos, int2d *ColVer2Ver)
{
   int1d son, fat, jVer;

   if (Pos[iVer] > 0) return;

   (*NbrHepLst)++;
   son = *NbrHepLst;
   //--- reorder the list
   while (son > 1)
   {
      fat  = son / 2;   // son>>1;
      jVer = HepLst[fat];
      if ((ColVer2Ver[iVer][0] > ColVer2Ver[jVer][0]) ||
          (ColVer2Ver[iVer][0] == ColVer2Ver[jVer][0] && ColVer2Ver[iVer][1] > ColVer2Ver[jVer][1]))
      {   // switch father and son
         HepLst[son] = jVer;
         Pos[jVer]   = son;
         son         = fat;
      }
      else
      {   // set vertex at is position
         HepLst[son] = iVer;
         Pos[iVer]   = son;
         return;
      }
   }

   if (son != 1)
   {
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
static int1d RemoveHeapListCol(int1d *NbrHepLst, int1d *HepLst, int1d *Pos, int2d *ColVer2Ver)
{
   int1d iVer, jVer, lVer, son1, son2, fat;

   //--- Remove the top vertex
   iVer      = HepLst[1];
   Pos[iVer] = -1;

   //--- Put the last vertex at the top of the list
   jVer = HepLst[*NbrHepLst];
   fat  = 1;
   (*NbrHepLst)--;

   //--- Traverse the list from top downwards to reposition the first element
   while ((*NbrHepLst) > 2 * fat)
   {
      son1 = 2 * fat;
      son2 = son1 + 1;
      if ((ColVer2Ver[HepLst[son1]][0] > ColVer2Ver[HepLst[son2]][0]) ||
          (ColVer2Ver[HepLst[son1]][0] == ColVer2Ver[HepLst[son2]][0] &&
           ColVer2Ver[HepLst[son1]][1] > ColVer2Ver[HepLst[son2]][1]))
      {
         if ((ColVer2Ver[HepLst[son1]][0] > ColVer2Ver[jVer][0]) ||
             (ColVer2Ver[HepLst[son1]][0] == ColVer2Ver[jVer][0] && ColVer2Ver[HepLst[son1]][1] > ColVer2Ver[jVer][1]))
         {
            lVer        = HepLst[son1];
            HepLst[fat] = lVer;
            Pos[lVer]   = fat;
            fat         = son1;
         }
         else
         {   // list is updated
            HepLst[fat] = jVer;
            Pos[jVer]   = fat;
            return iVer;
         }
      }
      else
      {
         if ((ColVer2Ver[HepLst[son2]][0] > ColVer2Ver[jVer][0]) ||
             (ColVer2Ver[HepLst[son2]][0] == ColVer2Ver[jVer][0] && ColVer2Ver[HepLst[son2]][1] > ColVer2Ver[jVer][1]))
         {   // switch jVer and son2 (Peut-etre fat and son2?)
            lVer        = HepLst[son2];
            HepLst[fat] = lVer;
            Pos[lVer]   = fat;
            fat         = son2;
         }
         else
         {   // list is updated
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
static void UpdateHeapListCol(int1d iVer, int1d *HepLst, int1d *Pos, int2d *ColVer2Ver)
{
   int1d jVer, son, fat;
   son = Pos[iVer];

   while (son > 1)
   {
      fat  = son / 2;
      jVer = HepLst[fat];

      if ((ColVer2Ver[iVer][0] > ColVer2Ver[jVer][0]) ||
          (ColVer2Ver[iVer][0] == ColVer2Ver[jVer][0] && ColVer2Ver[iVer][1] > ColVer2Ver[jVer][1]))
      {   // switch father and son
         HepLst[son] = jVer;
         Pos[jVer]   = son;
         son         = fat;
      }
      else
      {   // set vertex at is position
         HepLst[son] = iVer;
         Pos[iVer]   = son;
         return;
      }
   }

   //--- I am at the top of the list ie son = 1
   if (son != 1)
   {
      fprintf(stderr, "  ## ERROR UpdateHeapListCol: Update heap list son = %d \n", son);
      exit(1);
   }

   HepLst[son] = iVer;
   Pos[iVer]   = son;

   return;
}


static int *ColorPartImplicit(int NmbPth, int nbrPar, int NmbVer,
                              uint64_t *AdrBalRk1, int *LstBalRk1,
                              uint64_t *AdrBalRk2, int *LstBalRk2,
                              int *VidPar, int VrbLvl)
{
   int1d   NbrCol, i, flag1, flag = 0;
   int1d   iVer, iCol, jCol, iPar, jPar, iVoi, idx, idxPar, jdxPar;
   int1d   jVer, nbrVoi;
   int1d   idxCol, jdxCol, NbrHepLst, NbrParTgt, IteBck, test, colRef;
   int1d  *cntPar = NULL, *colPar = NULL;
   int1d  *LstCol = NULL;
   int1d  *HepLst, *Pos, *GphColPar;
   int2d  *ColVer2Ver;
   int56d  CntCol;
   int11d  BackUpCol;
   NbrParTgt = NmbPth;

   for (i = 0; i < 30; i++)
   {
      CntCol[i] = 0;
   }

   HepLst     = (int1d *)malloc((nbrPar + 1) * sizeof(int1d));
   Pos        = (int1d *)malloc((nbrPar + 1) * sizeof(int1d));
   ColVer2Ver = (int2d *)malloc((nbrPar + 1) * sizeof(int2d));
   cntPar     = (int1d *)calloc(nbrPar + 1, sizeof(int1d));
   colPar     = (int1d *)calloc((nbrPar + 1) * 1000, sizeof(int1d));
   GphColPar  = (int1d *)calloc((nbrPar + 1), sizeof(int1d));
   // use to order the partition according to their increasing degree, second int stock the idPar after ordering

   // Treat the Sub domains like Vertex in the point implicit then give the color to each vertex of each partition

   for (iVer = 1; iVer <= NmbVer; iVer++)
   {
      nbrVoi = AdrBalRk2[iVer+1] - AdrBalRk2[iVer];
      idxPar = VidPar[iVer];

      for (iVoi =  AdrBalRk2[iVer]; iVoi < AdrBalRk2[iVer+1]; iVoi++)
      {
         jVer   = LstBalRk2[iVoi];
         jdxPar = VidPar[jVer];

         if (idxPar == jdxPar)
            continue;

         else
         {
            flag1 = Wlf2_CheckIfVertexIsInTheList(jdxPar, cntPar[idxPar], &colPar[1000 * idxPar]);
            // flag1 = -1 if jdxPar is not in &colPar[1000*idxPar] ()

            if (cntPar[idxPar] > 1000 || cntPar[jdxPar] > 1000)
            {
               printf("MESH IS TO COARSE TO BE PARTITIONNED AND COLORED PROPERLY \n");
               exit(1);
            }

            if (flag1 == -1)
            {
               colPar[1000 * idxPar + cntPar[idxPar]] = jdxPar;
               cntPar[idxPar]++;
            }

            else
            {
               continue;
            }
         }
      }
   }

   LstCol       = (int1d *)calloc(100, sizeof(int1d));

   NbrHepLst = 0;
   idx       = 0;
   NbrCol    = 0;

   for (iPar = 1; iPar <= nbrPar; iPar++)
   {
      nbrVoi              = cntPar[iPar];
      ColVer2Ver[iPar][0] = 0;        // Nbr Voisins coloriés
      ColVer2Ver[iPar][1] = nbrVoi;   // Personne n'est colorié
      Pos[iPar]           = 0;

      if (nbrVoi > idx)
      {
         flag = iPar;
         idx  = nbrVoi;
      }
   }

   AddHeapListCol(flag, &NbrHepLst, HepLst, Pos, ColVer2Ver);

   for (iPar = 1; iPar <= nbrPar; iPar++)
   {
      iCol = 0;
      jPar = HepLst[1];

      RemoveHeapListCol(&NbrHepLst, HepLst, Pos, ColVer2Ver);

      nbrVoi = cntPar[jPar];

      if (nbrVoi > 500)
      {
         printf("Maillage de merde \n");
         exit(1);
      }

      for (iVoi = 0; iVoi < nbrVoi; iVoi++)
      {
         idxPar = colPar[1000 * jPar + iVoi];
         ColVer2Ver[idxPar][0]++;
         // ColVer2Ver[kVer][0]++;
         ColVer2Ver[idxPar][1]--;

         if (ColVer2Ver[idxPar][1] < 0)
         {
            printf("ON A COLORIÉ TROP DE VOISINS A idxPar = %d\n", idxPar);
            exit(1);
         }
         // printf(" \t iVoi = %d kVer = %d",iVoi, kVer);
         if (ColVer2Ver[idxPar][0] == 1)
         {
            AddHeapListCol(idxPar, &NbrHepLst, HepLst, Pos, ColVer2Ver);
         }
         if (ColVer2Ver[idxPar][0] > 1 && GphColPar[idxPar] == 0)
         {
            UpdateHeapListCol(idxPar, HepLst, Pos, ColVer2Ver);
         }

         LstCol[iVoi] = GphColPar[idxPar];
      }

      int1d nbrMin = nbrPar * 2;
      for (i = 1; i < 30; i++)
      {
         if ((CntCol[i] != 0) || i <= 10)
         {
            // if ( i<=9) {
            iVoi = 0;
            jCol = i;

            while (iVoi < nbrVoi)
            {
               if (LstCol[iVoi] == jCol)   // On skip la boucle car 'jCol est déja pris
                  goto nexCol;
               // else if ( LstCol[iVoi] < jCol && LstCol[iVoi+1] > jCol ) // On a trouvé que la couleur 'flag' c'est
               // bon
               //   iCol = jCol;
               // else
               iVoi++;
            }

            // this color is available
            if (iCol == 0)
            {
               iCol   = jCol;
               nbrMin = CntCol[i];
            }
            else
            {
               if (CntCol[i] < nbrMin)
               {
                  iCol   = jCol;
                  nbrMin = CntCol[i];
               }
            }
         }

      nexCol:
         continue;
      }


      if (iCol == 0)
      {   // On a fait tous les voisins et toutes les couleurs consécutives sont données
         iCol = LstCol[nbrVoi - 1] + 1;
      }

      GphColPar[jPar] = iCol;

      CntCol[iCol]++;

      if (iCol > NbrCol)
      {
         NbrCol = iCol;
      }
   }

   //---- Check up of evrything
   flag   = 0;
   colRef = 0;
   for (iCol = 1; iCol <= NbrCol; iCol++)
   {
      BackUpCol[iCol] = CntCol[iCol] - NbrParTgt;
      if (BackUpCol[iCol] != 0)
      {
         flag++;   // Coutn for the number of problem (must be even)
         if(VrbLvl >= 1)
            printf("  Color %d has a wrong number of partition %d (Target = %d)\n", iCol, CntCol[iCol], NbrParTgt);
      }
   }

   // printf("flag = %d so we have %d partitions to transfer \n", flag, flag/2);
   if (flag == 0)
      IteBck = 0;
   else
      IteBck = 5;

   while (flag > 0 && IteBck > 0)
   {
      idxCol = 0;
      jdxCol = 0;
      for (iCol = 1; iCol <= NbrCol; iCol++)
      {
         //--- find idxCol such that idxCol has to less and jdxCol to much partitions for a transfer
         if (BackUpCol[iCol] < 0 && idxCol == 0 && iCol != colRef) idxCol = iCol;

         if (BackUpCol[iCol] > 0 && jdxCol == 0) jdxCol = iCol;
      }
      //--if we find a couple of color for transfer lets try
      if (idxCol != 0 && jdxCol != 0)
      {
         colRef = idxCol;   //-- to prevent from treating the same color again and again
         //-- Looking for a partition for the transfer
         for (iPar = 1; iPar <= nbrPar; iPar++)
         {
            test   = 0;
            nbrVoi = cntPar[iPar];
            iCol   = GphColPar[iPar];


            if (iCol == jdxCol)
            {   //-- iPar is a partition that is in color with to many partition:
                // it is candidate for transfer, test for jdxCol
               // printf("Voisins de %d : ",iPar);
               for (iVoi = 0; iVoi < nbrVoi; iVoi++)
               {
                  jPar = colPar[1000 * (iPar) + iVoi];
                  jCol = GphColPar[jPar];
                  // printf("\t%% %d col : %d", jPar, jCol);
                  if (jCol == idxCol) test = 1;
               }
               // printf("\n");
               if (test == 0)
               {   // Test is good iPar is a champion !!!
                  CntCol[idxCol]++;
                  CntCol[jdxCol]--;
                  GphColPar[iPar] = idxCol;

                  BackUpCol[idxCol]++;
                  BackUpCol[jdxCol]--;

                  flag = flag - 2;
                  goto end;
               }
            }
         }
      }
   end:
      IteBck--;
   }

   if (VrbLvl >= 1 && flag != 0) printf("Balance is not good\n");

   return (GphColPar);
}

#endif
