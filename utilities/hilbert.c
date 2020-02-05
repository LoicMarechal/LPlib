

/*----------------------------------------------------------------------------*/
/*                                                                            */
/*                               HILBERT V 2.01                               */
/*                                                                            */
/*----------------------------------------------------------------------------*/
/*                                                                            */
/* Description:         renumber .meshb files                                 */
/* Author:              Loic MARECHAL                                         */
/* Creation date:       mar 11 2010                                           */
/* Last modification:   feb 05 2020                                           */
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
#include "libmeshb7.h"
#include "lplib3.h"


/*----------------------------------------------------------------------------*/
/* Defines                                                                    */
/*----------------------------------------------------------------------------*/

#define MAXITR    21
#define MAXELE    14
#define HILMOD    0
#define OCTMOD    1
#define RNDMOD    2
#define MIN(a,b)  ((a) < (b) ? (a) : (b))
#define MAX(a,b)  ((a) > (b) ? (a) : (b))


/*----------------------------------------------------------------------------*/
/* Structures                                                                 */
/*----------------------------------------------------------------------------*/

typedef struct
{
   uint64_t cod;
   double   crd[3];
   int      idx, ref;
}VerSct;

typedef struct
{
   uint64_t cod;
   int      *idx, siz, ref;
}EleSct;

typedef struct
{
   int      NmbVer, *Old2New, MshVer, dim, mod, TypIdx, VerTyp;
   int      NmbEle[ MAXELE ], *IdxTab[ MAXELE ], EleTyp[ MAXELE ];
   double   box[6];
   VerSct   *ver;
   EleSct   *ele[ MAXELE ];
}MshSct;


/*----------------------------------------------------------------------------*/
/* Global variables                                                           */
/*----------------------------------------------------------------------------*/

int EleTab[ MAXELE ][2] = {
   { 2, GmfEdges},
   { 3, GmfTriangles},
   { 4, GmfQuadrilaterals},
   { 4, GmfTetrahedra},
   { 5, GmfPyramids},
   { 6, GmfPrisms},
   { 8, GmfHexahedra},
   { 3, GmfEdgesP2},
   { 6, GmfTrianglesP2},
   { 9, GmfQuadrilateralsQ2},
   {10, GmfTetrahedraP2},
   {14, GmfPyramidsP2},
   {18, GmfPrismsP2},
   {27, GmfHexahedraQ2} };

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


/*----------------------------------------------------------------------------*/
/* Prototypes of local procedures                                             */
/*----------------------------------------------------------------------------*/

void     ScaMsh(char *, MshSct *);
void     RecMsh(char *, MshSct *);
uint64_t HilCod(double *, double *, int, int);
void     GetTim(double *);
int      CmpFnc(const void *, const void *);
void     RenVer(int, int, int, MshSct *);
void     RenEle(int, int, int, MshSct *);
void     PrtSta(MshSct *, int64_t);
void     SwpMem(MshSct *, int);


/*----------------------------------------------------------------------------*/
/* Read, renumber through a Hilbert SFC and write the mesh                    */
/*----------------------------------------------------------------------------*/

int main(int ArgCnt, char **ArgVec)
{
   char     *PtrArg, *TmpStr, InpNam[1000], OutNam[1000];
   int      i, j, t, NmbCpu = 0, StaFlg=0;
   int64_t  LibParIdx;
   float    flt[3], sta[2];
   double   timer=0;
   MshSct   msh;

   // Command line parsing
   memset(&msh, 0, sizeof(MshSct));

   if(ArgCnt == 1)
   {
      puts("\nHILBERT v2.01 february 05 2020   Loic MARECHAL / INRIA");
      puts(" Usage       : hilbert -in input_mesh -out renumbered_mesh");
      puts(" -in name    : name of the input mesh");
      puts(" -out name   : name of the output renumbered mesh");
      puts(" -stats      : print element blocks dependencies stats before and after renumbering");
      puts(" -scheme n   : renumbering scheme: 0 = Hilbert,  1 = Z curve,  2 = random,  (default = 0)");
      puts(" -nproc n    : n is the number of threads to be launched (default = all available threads)\n");
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

      if(!strcmp(PtrArg,"-stats"))
      {
         StaFlg = 1;
         continue;
      }

      if(!strcmp(PtrArg,"-scheme"))
      {
         msh.mod = atoi(*++ArgVec);
         msh.mod = MAX(msh.mod, 0);
         msh.mod = MIN(msh.mod, 2);
         ArgCnt--;
         continue;
      }

      if(!strcmp(PtrArg,"-nproc"))
      {
         NmbCpu = atoi(*++ArgVec);
         NmbCpu = MAX(NmbCpu, 1);
         NmbCpu = MIN(NmbCpu, 128);
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
   printf("\nReading mesh          : ");
   timer = 0.;
   GetTim(&timer);
   ScaMsh(InpNam, &msh);
   GetTim(&timer);
   printf("%g s\n", timer);

   printf("\nInput mesh : version %d\n", msh.MshVer);
   printf(" %d Vertices\n", msh.NmbVer);

   for(t=0;t<MAXELE;t++)
      if(msh.NmbEle[t])
         printf(" %d %s", msh.NmbEle[t], EleNam[t]);

   puts("");

   // Compute initial stats
   LibParIdx = InitParallel(NmbCpu);

   if(StaFlg)
   {
      msh.VerTyp = NewType(LibParIdx, msh.NmbVer);

      for(t=0;t<MAXELE;t++)
         if(msh.NmbEle[t])
            msh.EleTyp[t] = NewType(LibParIdx, msh.NmbEle[t]);

      puts("\nDependencies before renumbering (average / MAX) :");
      PrtSta(&msh, LibParIdx);
   }

   // Vertices renumbering
   printf("Renumbering vertices         : ");
   timer = 0.;
   GetTim(&timer);
   msh.VerTyp = NewType(LibParIdx, msh.NmbVer);
   LaunchParallel(LibParIdx, msh.VerTyp, 0, (void *)RenVer, (void *)&msh);
   ParallelQsort(LibParIdx, &msh.ver[1], msh.NmbVer, sizeof(VerSct), CmpFnc);

   msh.Old2New = malloc( (msh.NmbVer+1) * sizeof(int) );

   for(i=1;i<=msh.NmbVer;i++)
      msh.Old2New[ msh.ver[i].idx ] = i;

   GetTim(&timer);
   printf("%g s\n", timer);

   // Elements renumbering
   for(t=0;t<MAXELE;t++)
      if(msh.NmbEle[t])
      {
         printf("Renumbering %s : ", EleNam[t]);
         timer = 0.;
         GetTim(&timer);
         msh.EleTyp[t] = NewType(LibParIdx, msh.NmbEle[t]);
         msh.TypIdx = t;
         LaunchParallel(LibParIdx, msh.EleTyp[t], 0, (void *)RenEle, (void *)&msh);
         ParallelQsort(LibParIdx, &msh.ele[t][1], msh.NmbEle[t], sizeof(EleSct), CmpFnc);
         SwpMem(&msh, t);
         GetTim(&timer);
         printf("%g s\n", timer);
      }

   // Compute dependencies
   if(StaFlg)
   {
      puts("\nDependencies after renumbering (average / MAX) :");
      PrtSta(&msh, LibParIdx);
   }

   StopParallel(LibParIdx);

   // Mesh writing
   printf("Writing mesh          : ");
   timer = 0.;
   GetTim(&timer);
   RecMsh(OutNam, &msh);
   GetTim(&timer);
   printf("%g s\n\n", timer);

   // Release memory
   if(msh.ver)
      free(msh.ver);

   for(t=0;t<MAXELE;t++)
      if(msh.ele[t])
         free(msh.ele[t]);

   if(msh.Old2New)
      free(msh.Old2New);
}


/*----------------------------------------------------------------------------*/
/* Compute the bounding box                                                   */
/*----------------------------------------------------------------------------*/

static void SetBndBox(int64_t BegIdx, int64_t EndIdx, MshSct *msh)
{
   int i, j;

   for(i=BegIdx; i<=EndIdx; i++)
   {
      msh->ver[i].idx = i;

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

   for(j=0;j<3;j++)
      msh->box[j+3] = pow(2,64) / (msh->box[j+3] - msh->box[j]);
}


/*----------------------------------------------------------------------------*/
/* Open, allocate and read the mesh fields                                    */
/*----------------------------------------------------------------------------*/

void ScaMsh(char *InpNam, MshSct *msh)
{
   int      i, j, n, t, *PtrIdx;
   int64_t  InpMsh;
   float    flt[3];

   // Check mesh format
   if(!(InpMsh = GmfOpenMesh(InpNam, GmfRead, &msh->MshVer, &msh->dim)))
   {
      printf("Cannot open mesh %s\n", InpNam);
      exit(1);
   }

   // Get stats and allocate tables
   if((msh->NmbVer = GmfStatKwd(InpMsh, GmfVertices)))
      msh->ver = malloc((msh->NmbVer+1) * sizeof(VerSct));
   else
   {
      puts("Cannot renumber a mesh without vertices");
      exit(1);
   }

   // Read the vertices
   if(msh->NmbVer)
   {
      GmfGetBlock(InpMsh, GmfVertices, 1, msh->NmbVer, 0, NULL, SetBndBox, msh,
                  GmfDoubleVec, msh->dim, msh->ver[1].crd,  msh->ver[ msh->NmbVer ].crd,
                  GmfInt,                &msh->ver[1].ref, &msh->ver[ msh->NmbVer ].ref);
   }

   // Allocate and read elements
   for(t=0;t<MAXELE;t++)
   {
      n = msh->NmbEle[t] = GmfStatKwd(InpMsh, EleTab[t][1]);

      if(!n)
         continue;

      msh->ele[t]    = malloc( (n+1) * sizeof(EleSct) );
      msh->IdxTab[t] = malloc( (n+1) * EleTab[t][0] * sizeof(int));
      PtrIdx         = msh->IdxTab[t];

      for(i=1;i<=n;i++)
         msh->ele[t][i].idx = &PtrIdx[ i * EleTab[t][0] ];

      GmfGetBlock(InpMsh,    EleTab[t][1], 1, n, 0, NULL, NULL,
                  GmfIntVec, EleTab[t][0], msh->ele[t][1].idx,  msh->ele[t][n].idx,
                  GmfInt,                 &msh->ele[t][1].ref, &msh->ele[t][n].ref);
   }

   GmfCloseMesh(InpMsh);
}


/*----------------------------------------------------------------------------*/
/* Write mesh                                                                 */
/*----------------------------------------------------------------------------*/

void RecMsh(char *OutNam, MshSct *msh)
{
   int      i, t, n;
   int64_t  OutMsh;

   if(!(OutMsh = GmfOpenMesh(OutNam, GmfWrite, msh->MshVer, msh->dim)))
   {
      printf("Cannot create mesh %s\n", OutNam);
      exit(1);
   }

   if(msh->NmbVer)
   {
      GmfSetKwd(OutMsh, GmfVertices, msh->NmbVer);
      GmfSetBlock(OutMsh, GmfVertices, 1, msh->NmbVer, 0, NULL, NULL,
                  GmfDoubleVec, msh->dim, msh->ver[1].crd,  msh->ver[ msh->NmbVer ].crd,
                  GmfInt,                &msh->ver[1].ref, &msh->ver[ msh->NmbVer ].ref);
   }

   for(t=0;t<MAXELE;t++)
   {
      if(!(n = msh->NmbEle[t]))
         continue;

      GmfSetKwd  (OutMsh,    EleTab[t][1], n);
      GmfSetBlock(OutMsh,    EleTab[t][1], 1, n, 0, NULL, NULL,
                  GmfIntVec, EleTab[t][0], msh->ele[t][1].idx,  msh->ele[t][n].idx,
                  GmfInt,                 &msh->ele[t][1].ref, &msh->ele[t][n].ref);
   }

   GmfCloseMesh(OutMsh);
}


/*----------------------------------------------------------------------------*/
/* Compute the hilbert code from 3d coordinates                               */
/*----------------------------------------------------------------------------*/

uint64_t HilCod(double crd[3], double box[6], int itr, int mod)
{
   uint64_t IntCrd[3], m=1LL<<63, cod;
   int      i, j, b, GeoWrd, NewWrd, BitTab[3] = {1,2,4};
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
      IntCrd[j] = TmpCrd[j];
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
/* Wall clock timer                                                           */
/*----------------------------------------------------------------------------*/

void GetTim(double *timer)
{
   *timer = GetWallClock() - *timer;
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
/* Parallel loop renumbering vertices                                         */
/*----------------------------------------------------------------------------*/

void RenVer(int BegIdx, int EndIdx, int PthIdx, MshSct *msh)
{
   int i;

   for(i=BegIdx; i<=EndIdx; i++)
      msh->ver[i].cod = HilCod(msh->ver[i].crd, msh->box, MAXITR, msh->mod);
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
/* Compute the barycenter of any kind of element                              */
/*----------------------------------------------------------------------------*/

void SetMidCrd(int NmbVer, int *IdxTab, MshSct *msh, double *crd)
{
   int i, j;

   for(j=0;j<3;j++)
      crd[j] = 0.;

   for(i=0;i<NmbVer;i++)
      for(j=0;j<3;j++)
         crd[j] += msh->ver[ IdxTab[i] ].crd[j];

   for(j=0;j<3;j++)
      crd[j] /= NmbVer;

}

/*----------------------------------------------------------------------------*/
/* Parallel loop renumbering any kind of element                              */
/*----------------------------------------------------------------------------*/

void RenEle(int BegIdx, int EndIdx, int PthIdx, MshSct *msh)
{
   int      i;
   double   crd[3];

    for(i=BegIdx; i<=EndIdx; i++)
    {
       SetNewIdx(EleTab[ msh->TypIdx ][0], msh->ele[ msh->TypIdx ][i].idx, msh->Old2New);
       SetMidCrd(EleTab[ msh->TypIdx ][0], msh->ele[ msh->TypIdx ][i].idx, msh, crd);
       msh->ele[ msh->TypIdx ][i].cod = HilCod(crd, msh->box, MAXITR, msh->mod);
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
   int i, *IdxTab, *PtrIdx;

   // Allocate the new table
   IdxTab = malloc( (msh->NmbEle[ typ ] + 1) * EleTab[ typ ][0] * sizeof(int));
   PtrIdx = IdxTab;

   // Copy the old indices to the new table
   // and save the new pointer to each elements index adress
   for(i=1;i<=msh->NmbEle[ typ ];i++)
   {
      memcpy(PtrIdx, msh->ele[ typ ][i].idx, EleTab[ typ ][0] * sizeof(int));
      msh->ele[ typ ][i].idx = PtrIdx;
      PtrIdx += EleTab[ typ ][0];
   }

   // Free the old table and set the element type pointer with the new one
   free(msh->IdxTab[ typ ]);
   msh->IdxTab[ typ ] = IdxTab;
}
