

/*----------------------------------------------------------------------------*/
/*                                                                            */
/*                               HILBERT V 2.00                               */
/*                                                                            */
/*----------------------------------------------------------------------------*/
/*                                                                            */
/* Description:         renumber .meshb files                                 */
/* Author:              Loic MARECHAL                                         */
/* Creation date:       mar 11 2010                                           */
/* Last modification:   jul 25 2017                                           */
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

#define MAXITR 21
#define MAXELE 14
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define HILMOD 0
#define OCTMOD 1
#define RNDMOD 2


/*----------------------------------------------------------------------------*/
/* Structures                                                                 */
/*----------------------------------------------------------------------------*/

typedef struct
{
   uint64_t cod;
   double crd[3];
   int idx, ref;
}VerSct;

typedef struct
{
   uint64_t cod;
   int *idx, siz, ref;
}EleSct;

typedef struct
{
   int NmbVer, *Old2New, MshVer, dim, mod, TypIdx, VerTyp;
   int NmbEle[ MAXELE ], *IdxTab[ MAXELE ], EleTyp[ MAXELE ];
   double box[6];
   VerSct *ver;
   EleSct *ele[ MAXELE ];
}MshSct;


/*----------------------------------------------------------------------------*/
/* Global variables                                                           */
/*----------------------------------------------------------------------------*/

int EleTab[ MAXELE ][3] = {
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
uint64_t hilbert(double *, double *, int, int);
void     GetTim(double *);
int      CmpFnc(const void *, const void *);
void     RenVer(int, int, int, MshSct *);
void     RenEle(int, int, int, MshSct *);
void     PrtSta(MshSct *, int64_t);


/*----------------------------------------------------------------------------*/
/* Read, renumber through a Hilbert SFC and write the mesh                    */
/*----------------------------------------------------------------------------*/

int main(int ArgCnt, char **ArgVec)
{
   char *PtrArg, *TmpStr, InpNam[1000], OutNam[1000];
   int i, j, t, NmbCpu = 0, StaFlg=0;
   int64_t LibParIdx;
   float flt[3], sta[2];
   double timer=0;
   MshSct msh;

   // Command line parsing
   memset(&msh, 0, sizeof(MshSct));

   if(ArgCnt == 1)
   {
      puts("\nHILBERT v1.12 march 21 2018   Loic MARECHAL / INRIA");
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
/* Read mesh                                                                  */
/*----------------------------------------------------------------------------*/

void ScaMsh(char *InpNam, MshSct *msh)
{
   int i, j, t, *PtrIdx;
   int64_t InpMsh;
   float flt[3];

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

   for(t=0;t<MAXELE;t++)
      if((msh->NmbEle[t] = GmfStatKwd(InpMsh, EleTab[t][1])))
      {
         msh->ele[t] = malloc((msh->NmbEle[t]+1) * sizeof(EleSct));
         msh->IdxTab[t] = malloc((msh->NmbEle[t]+1) * EleTab[t][0] * sizeof(int));
         PtrIdx = msh->IdxTab[t];

         for(i=1;i<=msh->NmbEle[t];i++)
         {
            msh->ele[t][i].idx = PtrIdx;
            PtrIdx += EleTab[t][0];
         }
      }

   // Read fields
   if(msh->NmbVer)
   {
      if(msh->dim == 2)
      {
         GmfGetBlock(InpMsh, GmfVertices, 1, msh->NmbVer, 0, NULL, NULL, \
                     GmfDouble, &msh->ver[1].crd[0], &msh->ver[ msh->NmbVer ].crd[0], \
                     GmfDouble, &msh->ver[1].crd[1], &msh->ver[ msh->NmbVer ].crd[1], \
                     GmfInt, &msh->ver[1].ref,       &msh->ver[ msh->NmbVer ].ref);
      }
      else
      {
         GmfGetBlock(InpMsh, GmfVertices, 1, msh->NmbVer, 0, NULL, NULL, \
                     GmfDouble, &msh->ver[1].crd[0], &msh->ver[ msh->NmbVer ].crd[0], \
                     GmfDouble, &msh->ver[1].crd[1], &msh->ver[ msh->NmbVer ].crd[1], \
                     GmfDouble, &msh->ver[1].crd[2], &msh->ver[ msh->NmbVer ].crd[2], \
                     GmfInt, &msh->ver[1].ref,       &msh->ver[ msh->NmbVer ].ref);
      }

      for(i=1;i<=msh->NmbVer;i++)
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

   for(t=0;t<MAXELE;t++)
   {
      if(!msh->NmbEle[t])
         continue;

      switch(EleTab[t][0])
      {
         case 2 :
         {
            GmfGetBlock(InpMsh, EleTab[t][1], 1, msh->NmbEle[t], 0, NULL, NULL,
                        GmfInt, &msh->ele[t][1].idx[0], &msh->ele[t][ msh->NmbEle[t] ].idx[0],
                        GmfInt, &msh->ele[t][1].idx[1], &msh->ele[t][ msh->NmbEle[t] ].idx[1],
                        GmfInt, &msh->ele[t][1].ref,    &msh->ele[t][ msh->NmbEle[t] ].ref);
         }break;

         case 3 :
         {
            GmfGetBlock(InpMsh, EleTab[t][1], 1, msh->NmbEle[t], 0, NULL, NULL,
                        GmfInt, &msh->ele[t][1].idx[0], &msh->ele[t][ msh->NmbEle[t] ].idx[0],
                        GmfInt, &msh->ele[t][1].idx[1], &msh->ele[t][ msh->NmbEle[t] ].idx[1],
                        GmfInt, &msh->ele[t][1].idx[2], &msh->ele[t][ msh->NmbEle[t] ].idx[2],
                        GmfInt, &msh->ele[t][1].ref,    &msh->ele[t][ msh->NmbEle[t] ].ref);
         }break;

         case 4 :
         {
            GmfGetBlock(InpMsh, EleTab[t][1], 1, msh->NmbEle[t], 0, NULL, NULL,
                        GmfInt, &msh->ele[t][1].idx[0], &msh->ele[t][ msh->NmbEle[t] ].idx[0],
                        GmfInt, &msh->ele[t][1].idx[1], &msh->ele[t][ msh->NmbEle[t] ].idx[1],
                        GmfInt, &msh->ele[t][1].idx[2], &msh->ele[t][ msh->NmbEle[t] ].idx[2],
                        GmfInt, &msh->ele[t][1].idx[3], &msh->ele[t][ msh->NmbEle[t] ].idx[3],
                        GmfInt, &msh->ele[t][1].ref,    &msh->ele[t][ msh->NmbEle[t] ].ref);
         }break;

         case 5 :
         {
            GmfGetBlock(InpMsh, EleTab[t][1], 1, msh->NmbEle[t], 0, NULL, NULL,
                        GmfInt, &msh->ele[t][1].idx[0], &msh->ele[t][ msh->NmbEle[t] ].idx[0],
                        GmfInt, &msh->ele[t][1].idx[1], &msh->ele[t][ msh->NmbEle[t] ].idx[1],
                        GmfInt, &msh->ele[t][1].idx[2], &msh->ele[t][ msh->NmbEle[t] ].idx[2],
                        GmfInt, &msh->ele[t][1].idx[3], &msh->ele[t][ msh->NmbEle[t] ].idx[3],
                        GmfInt, &msh->ele[t][1].idx[4], &msh->ele[t][ msh->NmbEle[t] ].idx[4],
                        GmfInt, &msh->ele[t][1].ref,    &msh->ele[t][ msh->NmbEle[t] ].ref);
         }break;

         case 6 :
         {
            GmfGetBlock(InpMsh, EleTab[t][1], 1, msh->NmbEle[t], 0, NULL, NULL,
                        GmfInt, &msh->ele[t][1].idx[0], &msh->ele[t][ msh->NmbEle[t] ].idx[0],
                        GmfInt, &msh->ele[t][1].idx[1], &msh->ele[t][ msh->NmbEle[t] ].idx[1],
                        GmfInt, &msh->ele[t][1].idx[2], &msh->ele[t][ msh->NmbEle[t] ].idx[2],
                        GmfInt, &msh->ele[t][1].idx[3], &msh->ele[t][ msh->NmbEle[t] ].idx[3],
                        GmfInt, &msh->ele[t][1].idx[4], &msh->ele[t][ msh->NmbEle[t] ].idx[4],
                        GmfInt, &msh->ele[t][1].idx[5], &msh->ele[t][ msh->NmbEle[t] ].idx[5],
                        GmfInt, &msh->ele[t][1].ref,    &msh->ele[t][ msh->NmbEle[t] ].ref);
         }break;

         case 8 :
         {
            GmfGetBlock(InpMsh, EleTab[t][1], 1, msh->NmbEle[t], 0, NULL, NULL,
                        GmfInt, &msh->ele[t][1].idx[0], &msh->ele[t][ msh->NmbEle[t] ].idx[0],
                        GmfInt, &msh->ele[t][1].idx[1], &msh->ele[t][ msh->NmbEle[t] ].idx[1],
                        GmfInt, &msh->ele[t][1].idx[2], &msh->ele[t][ msh->NmbEle[t] ].idx[2],
                        GmfInt, &msh->ele[t][1].idx[3], &msh->ele[t][ msh->NmbEle[t] ].idx[3],
                        GmfInt, &msh->ele[t][1].idx[4], &msh->ele[t][ msh->NmbEle[t] ].idx[4],
                        GmfInt, &msh->ele[t][1].idx[5], &msh->ele[t][ msh->NmbEle[t] ].idx[5],
                        GmfInt, &msh->ele[t][1].idx[6], &msh->ele[t][ msh->NmbEle[t] ].idx[6],
                        GmfInt, &msh->ele[t][1].idx[7], &msh->ele[t][ msh->NmbEle[t] ].idx[7],
                        GmfInt, &msh->ele[t][1].ref,    &msh->ele[t][ msh->NmbEle[t] ].ref);
         }break;

         case 9 :
         {
            GmfGetBlock(InpMsh, EleTab[t][1], 1, msh->NmbEle[t], 0, NULL, NULL,
                        GmfInt, &msh->ele[t][1].idx[0], &msh->ele[t][ msh->NmbEle[t] ].idx[0],
                        GmfInt, &msh->ele[t][1].idx[1], &msh->ele[t][ msh->NmbEle[t] ].idx[1],
                        GmfInt, &msh->ele[t][1].idx[2], &msh->ele[t][ msh->NmbEle[t] ].idx[2],
                        GmfInt, &msh->ele[t][1].idx[3], &msh->ele[t][ msh->NmbEle[t] ].idx[3],
                        GmfInt, &msh->ele[t][1].idx[4], &msh->ele[t][ msh->NmbEle[t] ].idx[4],
                        GmfInt, &msh->ele[t][1].idx[5], &msh->ele[t][ msh->NmbEle[t] ].idx[5],
                        GmfInt, &msh->ele[t][1].idx[6], &msh->ele[t][ msh->NmbEle[t] ].idx[6],
                        GmfInt, &msh->ele[t][1].idx[7], &msh->ele[t][ msh->NmbEle[t] ].idx[7],
                        GmfInt, &msh->ele[t][1].idx[8], &msh->ele[t][ msh->NmbEle[t] ].idx[8],
                        GmfInt, &msh->ele[t][1].ref,    &msh->ele[t][ msh->NmbEle[t] ].ref);
         }break;

         case 10 :
         {
            GmfGetBlock(InpMsh, EleTab[t][1], 1, msh->NmbEle[t], 0, NULL, NULL,
                        GmfInt, &msh->ele[t][1].idx[0], &msh->ele[t][ msh->NmbEle[t] ].idx[0],
                        GmfInt, &msh->ele[t][1].idx[1], &msh->ele[t][ msh->NmbEle[t] ].idx[1],
                        GmfInt, &msh->ele[t][1].idx[2], &msh->ele[t][ msh->NmbEle[t] ].idx[2],
                        GmfInt, &msh->ele[t][1].idx[3], &msh->ele[t][ msh->NmbEle[t] ].idx[3],
                        GmfInt, &msh->ele[t][1].idx[4], &msh->ele[t][ msh->NmbEle[t] ].idx[4],
                        GmfInt, &msh->ele[t][1].idx[5], &msh->ele[t][ msh->NmbEle[t] ].idx[5],
                        GmfInt, &msh->ele[t][1].idx[6], &msh->ele[t][ msh->NmbEle[t] ].idx[6],
                        GmfInt, &msh->ele[t][1].idx[7], &msh->ele[t][ msh->NmbEle[t] ].idx[7],
                        GmfInt, &msh->ele[t][1].idx[8], &msh->ele[t][ msh->NmbEle[t] ].idx[8],
                        GmfInt, &msh->ele[t][1].idx[9], &msh->ele[t][ msh->NmbEle[t] ].idx[9],
                        GmfInt, &msh->ele[t][1].ref,   &msh->ele[t][ msh->NmbEle[t] ].ref);
         }break;

         case 14 :
         {
            GmfGetBlock(InpMsh, EleTab[t][1], 1, msh->NmbEle[t], 0, NULL, NULL,
                        GmfInt, &msh->ele[t][1].idx[ 0], &msh->ele[t][ msh->NmbEle[t] ].idx[ 0],
                        GmfInt, &msh->ele[t][1].idx[ 1], &msh->ele[t][ msh->NmbEle[t] ].idx[ 1],
                        GmfInt, &msh->ele[t][1].idx[ 2], &msh->ele[t][ msh->NmbEle[t] ].idx[ 2],
                        GmfInt, &msh->ele[t][1].idx[ 3], &msh->ele[t][ msh->NmbEle[t] ].idx[ 3],
                        GmfInt, &msh->ele[t][1].idx[ 4], &msh->ele[t][ msh->NmbEle[t] ].idx[ 4],
                        GmfInt, &msh->ele[t][1].idx[ 5], &msh->ele[t][ msh->NmbEle[t] ].idx[ 5],
                        GmfInt, &msh->ele[t][1].idx[ 6], &msh->ele[t][ msh->NmbEle[t] ].idx[ 6],
                        GmfInt, &msh->ele[t][1].idx[ 7], &msh->ele[t][ msh->NmbEle[t] ].idx[ 7],
                        GmfInt, &msh->ele[t][1].idx[ 8], &msh->ele[t][ msh->NmbEle[t] ].idx[ 8],
                        GmfInt, &msh->ele[t][1].idx[ 9], &msh->ele[t][ msh->NmbEle[t] ].idx[ 9],
                        GmfInt, &msh->ele[t][1].idx[10], &msh->ele[t][ msh->NmbEle[t] ].idx[10],
                        GmfInt, &msh->ele[t][1].idx[11], &msh->ele[t][ msh->NmbEle[t] ].idx[11],
                        GmfInt, &msh->ele[t][1].idx[12], &msh->ele[t][ msh->NmbEle[t] ].idx[12],
                        GmfInt, &msh->ele[t][1].idx[13], &msh->ele[t][ msh->NmbEle[t] ].idx[13],
                        GmfInt, &msh->ele[t][1].ref,    &msh->ele[t][ msh->NmbEle[t] ].ref);
         }break;

         case 18 :
         {
            GmfGetBlock(InpMsh, EleTab[t][1], 1, msh->NmbEle[t], 0, NULL, NULL,
                        GmfInt, &msh->ele[t][1].idx[ 0], &msh->ele[t][ msh->NmbEle[t] ].idx[ 0],
                        GmfInt, &msh->ele[t][1].idx[ 1], &msh->ele[t][ msh->NmbEle[t] ].idx[ 1],
                        GmfInt, &msh->ele[t][1].idx[ 2], &msh->ele[t][ msh->NmbEle[t] ].idx[ 2],
                        GmfInt, &msh->ele[t][1].idx[ 3], &msh->ele[t][ msh->NmbEle[t] ].idx[ 3],
                        GmfInt, &msh->ele[t][1].idx[ 4], &msh->ele[t][ msh->NmbEle[t] ].idx[ 4],
                        GmfInt, &msh->ele[t][1].idx[ 5], &msh->ele[t][ msh->NmbEle[t] ].idx[ 5],
                        GmfInt, &msh->ele[t][1].idx[ 6], &msh->ele[t][ msh->NmbEle[t] ].idx[ 6],
                        GmfInt, &msh->ele[t][1].idx[ 7], &msh->ele[t][ msh->NmbEle[t] ].idx[ 7],
                        GmfInt, &msh->ele[t][1].idx[ 8], &msh->ele[t][ msh->NmbEle[t] ].idx[ 8],
                        GmfInt, &msh->ele[t][1].idx[ 9], &msh->ele[t][ msh->NmbEle[t] ].idx[ 9],
                        GmfInt, &msh->ele[t][1].idx[10], &msh->ele[t][ msh->NmbEle[t] ].idx[10],
                        GmfInt, &msh->ele[t][1].idx[11], &msh->ele[t][ msh->NmbEle[t] ].idx[11],
                        GmfInt, &msh->ele[t][1].idx[12], &msh->ele[t][ msh->NmbEle[t] ].idx[12],
                        GmfInt, &msh->ele[t][1].idx[13], &msh->ele[t][ msh->NmbEle[t] ].idx[13],
                        GmfInt, &msh->ele[t][1].idx[14], &msh->ele[t][ msh->NmbEle[t] ].idx[14],
                        GmfInt, &msh->ele[t][1].idx[15], &msh->ele[t][ msh->NmbEle[t] ].idx[15],
                        GmfInt, &msh->ele[t][1].idx[16], &msh->ele[t][ msh->NmbEle[t] ].idx[16],
                        GmfInt, &msh->ele[t][1].idx[17], &msh->ele[t][ msh->NmbEle[t] ].idx[17],
                        GmfInt, &msh->ele[t][1].ref,    &msh->ele[t][ msh->NmbEle[t] ].ref);
         }break;

         case 27 :
         {
            GmfGetBlock(InpMsh, EleTab[t][1], 1, msh->NmbEle[t], 0, NULL, NULL,
                        GmfInt, &msh->ele[t][1].idx[ 0], &msh->ele[t][ msh->NmbEle[t] ].idx[ 0],
                        GmfInt, &msh->ele[t][1].idx[ 1], &msh->ele[t][ msh->NmbEle[t] ].idx[ 1],
                        GmfInt, &msh->ele[t][1].idx[ 2], &msh->ele[t][ msh->NmbEle[t] ].idx[ 2],
                        GmfInt, &msh->ele[t][1].idx[ 3], &msh->ele[t][ msh->NmbEle[t] ].idx[ 3],
                        GmfInt, &msh->ele[t][1].idx[ 4], &msh->ele[t][ msh->NmbEle[t] ].idx[ 4],
                        GmfInt, &msh->ele[t][1].idx[ 5], &msh->ele[t][ msh->NmbEle[t] ].idx[ 5],
                        GmfInt, &msh->ele[t][1].idx[ 6], &msh->ele[t][ msh->NmbEle[t] ].idx[ 6],
                        GmfInt, &msh->ele[t][1].idx[ 7], &msh->ele[t][ msh->NmbEle[t] ].idx[ 7],
                        GmfInt, &msh->ele[t][1].idx[ 8], &msh->ele[t][ msh->NmbEle[t] ].idx[ 8],
                        GmfInt, &msh->ele[t][1].idx[ 9], &msh->ele[t][ msh->NmbEle[t] ].idx[ 9],
                        GmfInt, &msh->ele[t][1].idx[10], &msh->ele[t][ msh->NmbEle[t] ].idx[10],
                        GmfInt, &msh->ele[t][1].idx[11], &msh->ele[t][ msh->NmbEle[t] ].idx[11],
                        GmfInt, &msh->ele[t][1].idx[12], &msh->ele[t][ msh->NmbEle[t] ].idx[12],
                        GmfInt, &msh->ele[t][1].idx[13], &msh->ele[t][ msh->NmbEle[t] ].idx[13],
                        GmfInt, &msh->ele[t][1].idx[14], &msh->ele[t][ msh->NmbEle[t] ].idx[14],
                        GmfInt, &msh->ele[t][1].idx[15], &msh->ele[t][ msh->NmbEle[t] ].idx[15],
                        GmfInt, &msh->ele[t][1].idx[16], &msh->ele[t][ msh->NmbEle[t] ].idx[16],
                        GmfInt, &msh->ele[t][1].idx[17], &msh->ele[t][ msh->NmbEle[t] ].idx[17],
                        GmfInt, &msh->ele[t][1].idx[18], &msh->ele[t][ msh->NmbEle[t] ].idx[18],
                        GmfInt, &msh->ele[t][1].idx[19], &msh->ele[t][ msh->NmbEle[t] ].idx[19],
                        GmfInt, &msh->ele[t][1].idx[20], &msh->ele[t][ msh->NmbEle[t] ].idx[20],
                        GmfInt, &msh->ele[t][1].idx[21], &msh->ele[t][ msh->NmbEle[t] ].idx[21],
                        GmfInt, &msh->ele[t][1].idx[22], &msh->ele[t][ msh->NmbEle[t] ].idx[22],
                        GmfInt, &msh->ele[t][1].idx[23], &msh->ele[t][ msh->NmbEle[t] ].idx[23],
                        GmfInt, &msh->ele[t][1].idx[24], &msh->ele[t][ msh->NmbEle[t] ].idx[24],
                        GmfInt, &msh->ele[t][1].idx[25], &msh->ele[t][ msh->NmbEle[t] ].idx[25],
                        GmfInt, &msh->ele[t][1].idx[26], &msh->ele[t][ msh->NmbEle[t] ].idx[26],
                        GmfInt, &msh->ele[t][1].ref,    &msh->ele[t][ msh->NmbEle[t] ].ref);
         }break;
      }
   }

   GmfCloseMesh(InpMsh);
}


/*----------------------------------------------------------------------------*/
/* Write mesh                                                                 */
/*----------------------------------------------------------------------------*/

void RecMsh(char *OutNam, MshSct *msh)
{
   int i, t;
   int64_t OutMsh;

   if(!(OutMsh = GmfOpenMesh(OutNam, GmfWrite, msh->MshVer, msh->dim)))
   {
      printf("Cannot create mesh %s\n", OutNam);
      exit(1);
   }

   if(msh->NmbVer)
   {
      GmfSetKwd(OutMsh, GmfVertices, msh->NmbVer);

      if(msh->dim == 2)
      {
         GmfSetBlock(OutMsh, GmfVertices, 1, msh->NmbVer, 0, NULL, NULL, \
                     GmfDouble, &msh->ver[1].crd[0], &msh->ver[ msh->NmbVer ].crd[0], \
                     GmfDouble, &msh->ver[1].crd[1], &msh->ver[ msh->NmbVer ].crd[1], \
                     GmfInt, &msh->ver[1].ref,       &msh->ver[ msh->NmbVer ].ref);
      }
      else
      {
         GmfSetBlock(OutMsh, GmfVertices, 1, msh->NmbVer, 0, NULL, NULL, \
                     GmfDouble, &msh->ver[1].crd[0], &msh->ver[ msh->NmbVer ].crd[0], \
                     GmfDouble, &msh->ver[1].crd[1], &msh->ver[ msh->NmbVer ].crd[1], \
                     GmfDouble, &msh->ver[1].crd[2], &msh->ver[ msh->NmbVer ].crd[2], \
                     GmfInt, &msh->ver[1].ref,       &msh->ver[ msh->NmbVer ].ref);
      }
   }

   for(t=0;t<MAXELE;t++)
   {
      if(!msh->NmbEle[t])
         continue;

      GmfSetKwd(OutMsh, EleTab[t][1], msh->NmbEle[t]);

      switch(EleTab[t][0])
      {
         case 2 :
         {
            for(i=1;i<=msh->NmbEle[t];i++)
               GmfSetLin(  OutMsh, EleTab[t][1],  msh->ele[t][i].idx[0],
                           msh->ele[t][i].idx[1], msh->ele[t][i].ref);
         }break;

         case 3 :
         {
            for(i=1;i<=msh->NmbEle[t];i++)
               GmfSetLin(  OutMsh, EleTab[t][1],  msh->ele[t][i].idx[0],
                           msh->ele[t][i].idx[1], msh->ele[t][i].idx[2],
                           msh->ele[t][i].ref);
         }break;

         case 4 :
         {
            for(i=1;i<=msh->NmbEle[t];i++)
               GmfSetLin(  OutMsh, EleTab[t][1],  msh->ele[t][i].idx[0],
                           msh->ele[t][i].idx[1], msh->ele[t][i].idx[2],
                           msh->ele[t][i].idx[3], msh->ele[t][i].ref);
         }break;

         case 5 :
         {
            for(i=1;i<=msh->NmbEle[t];i++)
               GmfSetLin(  OutMsh, EleTab[t][1],  msh->ele[t][i].idx[0],
                           msh->ele[t][i].idx[1], msh->ele[t][i].idx[2],
                           msh->ele[t][i].idx[3], msh->ele[t][i].idx[4],
                           msh->ele[t][i].ref);
         }break;

         case 6 :
         {
            for(i=1;i<=msh->NmbEle[t];i++)
               GmfSetLin(  OutMsh, EleTab[t][1],  msh->ele[t][i].idx[0],
                           msh->ele[t][i].idx[1], msh->ele[t][i].idx[2],
                           msh->ele[t][i].idx[3], msh->ele[t][i].idx[4],
                           msh->ele[t][i].idx[5], msh->ele[t][i].ref);
         }break;

         case 8 :
         {
            for(i=1;i<=msh->NmbEle[t];i++)
               GmfSetLin(  OutMsh, EleTab[t][1],  msh->ele[t][i].idx[0],
                           msh->ele[t][i].idx[1], msh->ele[t][i].idx[2],
                           msh->ele[t][i].idx[3], msh->ele[t][i].idx[4],
                           msh->ele[t][i].idx[5], msh->ele[t][i].idx[6],
                           msh->ele[t][i].idx[7], msh->ele[t][i].ref);
         }break;

         case 9 :
         {
            for(i=1;i<=msh->NmbEle[t];i++)
               GmfSetLin(  OutMsh, EleTab[t][1],  msh->ele[t][i].idx[0],
                           msh->ele[t][i].idx[1], msh->ele[t][i].idx[2],
                           msh->ele[t][i].idx[3], msh->ele[t][i].idx[4],
                           msh->ele[t][i].idx[5], msh->ele[t][i].idx[6],
                           msh->ele[t][i].idx[7], msh->ele[t][i].idx[8],
                           msh->ele[t][i].ref);
         }break;

         case 10 :
         {
            for(i=1;i<=msh->NmbEle[t];i++)
               GmfSetLin(  OutMsh, EleTab[t][1],  msh->ele[t][i].idx[0],
                           msh->ele[t][i].idx[1], msh->ele[t][i].idx[2],
                           msh->ele[t][i].idx[3], msh->ele[t][i].idx[4],
                           msh->ele[t][i].idx[5], msh->ele[t][i].idx[6],
                           msh->ele[t][i].idx[7], msh->ele[t][i].idx[8],
                           msh->ele[t][i].idx[9], msh->ele[t][i].ref);
         }break;

         case 14 :
         {
            for(i=1;i<=msh->NmbEle[t];i++)
               GmfSetLin(  OutMsh, EleTab[t][1],   msh->ele[t][i].idx[ 0],
                           msh->ele[t][i].idx[ 1], msh->ele[t][i].idx[ 2],
                           msh->ele[t][i].idx[ 3], msh->ele[t][i].idx[ 4],
                           msh->ele[t][i].idx[ 5], msh->ele[t][i].idx[ 6],
                           msh->ele[t][i].idx[ 7], msh->ele[t][i].idx[ 8],
                           msh->ele[t][i].idx[ 9], msh->ele[t][i].idx[10],
                           msh->ele[t][i].idx[11], msh->ele[t][i].idx[12],
                           msh->ele[t][i].idx[13], msh->ele[t][i].ref);
         }break;

         case 18 :
         {
            for(i=1;i<=msh->NmbEle[t];i++)
               GmfSetLin(  OutMsh, EleTab[t][1],   msh->ele[t][i].idx[ 0],
                           msh->ele[t][i].idx[ 1], msh->ele[t][i].idx[ 2],
                           msh->ele[t][i].idx[ 3], msh->ele[t][i].idx[ 4],
                           msh->ele[t][i].idx[ 5], msh->ele[t][i].idx[ 6],
                           msh->ele[t][i].idx[ 7], msh->ele[t][i].idx[ 8],
                           msh->ele[t][i].idx[ 9], msh->ele[t][i].idx[10],
                           msh->ele[t][i].idx[11], msh->ele[t][i].idx[12],
                           msh->ele[t][i].idx[13], msh->ele[t][i].idx[14],
                           msh->ele[t][i].idx[15], msh->ele[t][i].idx[16],
                           msh->ele[t][i].idx[17], msh->ele[t][i].ref);
         }break;

         case 27 :
         {
            for(i=1;i<=msh->NmbEle[t];i++)
               GmfSetLin(  OutMsh, EleTab[t][1],   msh->ele[t][i].idx[ 0],
                           msh->ele[t][i].idx[ 1], msh->ele[t][i].idx[ 2],
                           msh->ele[t][i].idx[ 3], msh->ele[t][i].idx[ 4],
                           msh->ele[t][i].idx[ 5], msh->ele[t][i].idx[ 6],
                           msh->ele[t][i].idx[ 7], msh->ele[t][i].idx[ 8],
                           msh->ele[t][i].idx[ 9], msh->ele[t][i].idx[10],
                           msh->ele[t][i].idx[11], msh->ele[t][i].idx[12],
                           msh->ele[t][i].idx[13], msh->ele[t][i].idx[14],
                           msh->ele[t][i].idx[15], msh->ele[t][i].idx[16],
                           msh->ele[t][i].idx[17], msh->ele[t][i].idx[18],
                           msh->ele[t][i].idx[19], msh->ele[t][i].idx[20],
                           msh->ele[t][i].idx[21], msh->ele[t][i].idx[22],
                           msh->ele[t][i].idx[23], msh->ele[t][i].idx[24],
                           msh->ele[t][i].idx[25], msh->ele[t][i].idx[26],
                           msh->ele[t][i].ref);
         }break;
      }
   }

   GmfCloseMesh(OutMsh);
}


/*----------------------------------------------------------------------------*/
/* Compute the hilbert code from 3d coordinates                               */
/*----------------------------------------------------------------------------*/

uint64_t hilbert(double crd[3], double box[6], int itr, int mod)
{
   uint64_t IntCrd[3], m=1LL<<63, cod;
   int i, j, b, GeoWrd, NewWrd, BitTab[3] = {1,2,4};
   double TmpCrd[3];
   int rot[8], GeoCod[8]={0,3,7,4,1,2,6,5}; // Hilbert
   int OctCod[8] = {5,4,7,6,1,0,3,2}; // octree or Z curve
   int HilCod[8][8] = {
      {0,7,6,1,2,5,4,3}, {0,3,4,7,6,5,2,1}, {0,3,4,7,6,5,2,1}, {2,3,0,1,6,7,4,5},
      {2,3,0,1,6,7,4,5}, {6,5,2,1,0,3,4,7}, {6,5,2,1,0,3,4,7}, {4,3,2,5,6,1,0,7}};

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
      msh->ver[i].cod = hilbert(msh->ver[i].crd, msh->box, MAXITR, msh->mod);
}


void SetNewIdx(int NmbVer, int *IdxTab, int *Old2New)
{
   int i;

   for(i=0;i<NmbVer;i++)
      IdxTab[i] = Old2New[ IdxTab[i] ];
}

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
/* Parallel loop renumbering tets                                             */
/*----------------------------------------------------------------------------*/

void RenEle(int BegIdx, int EndIdx, int PthIdx, MshSct *msh)
{
   int i;
   double crd[3];

    for(i=BegIdx; i<=EndIdx; i++)
    {
       SetNewIdx(EleTab[ msh->TypIdx ][0], msh->ele[ msh->TypIdx ][i].idx, msh->Old2New);
       SetMidCrd(EleTab[ msh->TypIdx ][0], msh->ele[ msh->TypIdx ][i].idx, msh, crd);
       msh->ele[ msh->TypIdx ][i].cod = hilbert(crd, msh->box, MAXITR, msh->mod);
    }
}


/*----------------------------------------------------------------------------*/
/* Compute and print elements / vertices dependencies stats                   */
/*----------------------------------------------------------------------------*/

void PrtSta(MshSct *msh, int64_t LibParIdx)
{
   int i, t;
   float  sta[2];

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
