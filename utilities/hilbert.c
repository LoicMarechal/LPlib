

/*----------------------------------------------------------------------------*/
/*                                                                            */
/*                               HILBERT V 2.23                               */
/*                                                                            */
/*----------------------------------------------------------------------------*/
/*                                                                            */
/* Description:         renumber .meshb files                                 */
/* Author:              Loic MARECHAL                                         */
/* Creation date:       mar 11 2010                                           */
/* Last modification:   mar 31 2022                                           */
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
#include <float.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include "libmeshb7.h"
#include "libmeshb7_helpers.h"
#include "lplib3.h"


/*----------------------------------------------------------------------------*/
/* Defines                                                                    */
/*----------------------------------------------------------------------------*/

#define MAXITR    21
#define MAXELE    14
#define MAXDOM    1024
#define HILMOD    0
#define OCTMOD    1
#define RNDMOD    2

enum {TypEdg, TypTri, TypQad, TypTet, TypPyr, TypPri, TypHex};


/*----------------------------------------------------------------------------*/
/* Macros                                                                     */
/*----------------------------------------------------------------------------*/

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
#ifdef WITH_GMLIB
   int      dom, gid, lid;
#endif
}VerSct;

typedef struct
{
   uint64_t cod;
   int      *idx, siz, ref, gid;
#ifdef WITH_GMLIB
   int      dom, lid;
#endif
}EleSct;

typedef struct PtrBuc
{
   int      idx[3], ele;
   char     typ, voy;
   struct   PtrBuc *nex;
}BucSct;

typedef struct
{
   int      NmbVer, *Old2New, MshVer, dim, mod, TypIdx, VerTyp, GmlMod;
   int      NmbEle[ MAXELE ], *IdxTab[ MAXELE ], EleTyp[ MAXELE ];
   int      MaxDeg[ MAXELE ], DegVec[ MAXELE ], HghDeg, OvrDeg;
   int      MtsEleTyp, MtsNmbDom;
   double   box[6];
   VerSct   *ver;
   EleSct   *ele[ MAXELE ];
}MshSct;


/*----------------------------------------------------------------------------*/
/* Global variables                                                           */
/*----------------------------------------------------------------------------*/

// For each kind of element: number of node, number of faces, GMF keyword
int EleTab[ MAXELE ][3] = {
   { 2, 2, GmfEdges},
   { 3, 3, GmfTriangles},
   { 4, 4, GmfQuadrilaterals},
   { 4, 4, GmfTetrahedra},
   { 5, 5, GmfPyramids},
   { 6, 5, GmfPrisms},
   { 8, 6, GmfHexahedra},
   { 3, 2, GmfEdgesP2},
   { 6, 3, GmfTrianglesP2},
   { 9, 4, GmfQuadrilateralsQ2},
   {10, 4, GmfTetrahedraP2},
   {14, 5, GmfPyramidsP2},
   {18, 5, GmfPrismsP2},
   {27, 6, GmfHexahedraQ2} };

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

#ifdef WITH_GMLIB
// For each kind of element: low vertex degree and highest vertex degree
int MaxDeg[ MAXELE ][2] = {
   { 2,   8},
   { 8,  32},
   { 4,  16},
   {28, 128},
   {16,  64},
   {16,  64},
   { 8,  32},
   { 2,   8},
   { 8,  32},
   { 4,  16},
   {28, 128},
   {16,  64},
   {16,  64},
   { 8,  32} };
#endif

// For each kind of element: give each face number of nodes
int FacDeg[7][6] = { 
   {0,0,0,0,0,0}, {3,0,0,0,0,0}, {4,0,0,0,0,0},
   {3,3,3,3,0,0}, {3,3,3,3,4,0}, {3,3,4,4,4,0}, {4,4,4,4,4,4} };

// For each kind of element: give eahc face list of nodes
int EleFac[7][6][4] = {
   { {0,0,0,0}, {0,0,0,0}, {0,0,0,0}, {0,0,0,0}, {0,0,0,0}, {0,0,0,0} },
   { {0,1,2,0}, {0,0,0,0}, {0,0,0,0}, {0,0,0,0}, {0,0,0,0}, {0,0,0,0} },
   { {0,1,2,3}, {0,0,0,0}, {0,0,0,0}, {0,0,0,0}, {0,0,0,0}, {0,0,0,0} },
   { {1,2,3,0}, {2,0,3,0}, {3,0,1,0}, {0,2,1,0}, {0,0,0,0}, {0,0,0,0} },
   { {0,1,4,0}, {1,2,4,0}, {2,3,4,0}, {3,0,4,0}, {3,2,1,0}, {0,0,0,0} },
   { {0,2,1,0}, {3,4,5,0}, {0,1,4,3}, {1,2,5,4}, {3,5,2,0}, {0,0,0,0} },
   { {0,4,7,3}, {1,2,6,5}, {0,1,5,4}, {3,7,6,2}, {0,3,2,1}, {4,5,6,7} } };


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
char    *SetNgb(MshSct *);
void     SrtFac(EleSct *, int, int, int [3]);
#ifdef WITH_GMLIB
void     GmlSrt(MshSct *);
void     ScaMts(char *, char *, MshSct *);
#endif


/*----------------------------------------------------------------------------*/
/* Read, renumber through a Hilbert SFC and write the mesh                    */
/*----------------------------------------------------------------------------*/

int main(int ArgCnt, char **ArgVec)
{
   char     *PtrArg, *TmpStr, *BndTab, InpNam[1000], OutNam[1000];
   int      i, t, NmbCpu = 0, StaFlg = 0, BndFlg = 0, *IdxTab, NmbSrf, NmbVol;
   int64_t  LibParIdx;
   double   timer = 0.;
   MshSct   msh;
   VerSct   *OldVer, *NewVer, *VolVer;
#ifdef WITH_GMLIB
   int      j;
   char     MtsNodNam[1000], MtsEleNam[1000];
#endif


   // --------------------
   // Command line parsing
   // --------------------

   memset(&msh, 0, sizeof(MshSct));

   if(ArgCnt == 1)
   {
      puts("\nHILBERT v2.23 march 31 2022   Loic MARECHAL / INRIA");
      puts(" Usage       : hilbert -in input_mesh -out renumbered_mesh");
      puts(" -in name    : input mesh(b) name");
      puts(" -out name   : output renumbered mesh(b)");
      puts(" -stats      : print element blocks dependencies stats before and after renumbering");
      puts(" -fixbnd     : do not renumber boundary nodes");
      puts(" -scheme n   : renumbering scheme: 0 = Hilbert,  1 = Z curve,  2 = random,  (default = 0)");
      puts(" -nproc n    : n is the number of threads to be launched (default = all available threads)");
#ifdef WITH_GMLIB
      puts(" -gmlib      : special sorting to make the mesh fit for the GMlib");
      puts(" -ndom n     : n: number of domains generated by Metis");
      puts(" -npart pf   : pf: Metis node partition filename");
      puts(" -epart pf t : pf: Metis element partition filename, t: element type");
#endif
      puts("");
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

      if(!strcmp(PtrArg,"-fixbnd"))
      {
         BndFlg = 1;
         continue;
      }

#ifdef WITH_GMLIB
      if(!strcmp(PtrArg,"-ndom"))
      {
         msh.MtsNmbDom = atoi(*++ArgVec);
         msh.MtsNmbDom = MIN(msh.MtsNmbDom, MAXDOM);
         ArgCnt--;
         continue;
      }

      if(!strcmp(PtrArg,"-npart"))
      {
         TmpStr = *++ArgVec;
         ArgCnt--;
         strcpy(MtsNodNam, TmpStr);
         continue;
      }

      if(!strcmp(PtrArg,"-epart"))
      {
         TmpStr = *++ArgVec;
         ArgCnt--;
         strcpy(MtsEleNam, TmpStr);
         msh.MtsEleTyp = atoi(*++ArgVec);
         ArgCnt--;
         continue;
      }

      if(!strcmp(PtrArg,"-gmlib"))
      {
         msh.GmlMod = 1;
         continue;
      }
#endif

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


   // ------------
   // Mesh reading
   // ------------

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


   // --------------------------------------------
   // Compute initial block dependenies statistics
   // --------------------------------------------

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


   // ----------------------------------------------------------------------
   // Prepare references and degree related data for the GMlib modified sort
   // ----------------------------------------------------------------------

#ifdef WITH_GMLIB
   if(msh.GmlMod)
   {
      if(msh.MtsNmbDom)
         ScaMts(MtsNodNam, MtsEleNam, &msh);

      GmlSrt(&msh);

      if(StaFlg)
      {
         printf(  "High-connected vertices      : %3.6f%%\n",
                  (100. * (float)msh.HghDeg) / (float)msh.NmbVer );

         printf(  "Over-connected vertices      : %3.6f%%\n",
                  (100. * (float)msh.OvrDeg) / (float)msh.NmbVer );

         for(j=0;j<MAXELE;j++)
            if(msh.MaxDeg[j])
               printf(  "Ball of %s     : max deg = %3d, vec size = %3d\n",
                        EleNam[j], msh.MaxDeg[j],  msh.DegVec[j] );

         puts("");
      }
   }
#endif


   // --------------------
   // Vertices renumbering
   // --------------------

   printf("Renumbering vertices         : ");
   timer = 0.;
   GetTim(&timer);

   if(BndFlg)
   {
      // Set the neighbours between elements and set boundary vertices tags
      BndTab = SetNgb(&msh);

      // Extract the inner volume vertices in a separate table
      IdxTab = malloc( (msh.NmbVer + 1) * sizeof(int) );
      assert(IdxTab);
      NmbVol = 0;

      for(i=1;i<=msh.NmbVer;i++)
         IdxTab[i] = BndTab[i] ? 0 : ++NmbVol;

      VolVer = malloc( (NmbVol + 1) * sizeof(VerSct) );
      assert(VolVer);

      for(i=1;i<=msh.NmbVer;i++)
         if(!BndTab[i])
            memcpy(&VolVer[ IdxTab[i] ], &msh.ver[i], sizeof(VerSct));

      OldVer = msh.ver;
      msh.ver = VolVer;
   
      // Renumber the inner volume vertices only
      msh.VerTyp = NewType(LibParIdx, NmbVol);
      LaunchParallel(LibParIdx, msh.VerTyp, 0, (void *)RenVer, (void *)&msh);
      ParallelQsort(LibParIdx, &msh.ver[1], NmbVol, sizeof(VerSct), CmpFnc);

      // Join the surface vertices with the renumbered volume ones
      NewVer = malloc( (msh.NmbVer + 1) * sizeof(VerSct) );
      assert(NewVer);

      for(i=1;i<=msh.NmbVer;i++)
         if(BndTab[i])
            memcpy(&NewVer[i], &OldVer[i], sizeof(VerSct));
         else
            memcpy(&NewVer[i], &VolVer[ IdxTab[i] ], sizeof(VerSct));

      msh.ver = NewVer;

      // Free everything
      free(VolVer);
      free(OldVer);
      free(IdxTab);
      free(BndTab);
   }
   else
   {
      msh.VerTyp = NewType(LibParIdx, msh.NmbVer);
      LaunchParallel(LibParIdx, msh.VerTyp, 0, (void *)RenVer, (void *)&msh);
      ParallelQsort(LibParIdx, &msh.ver[1], msh.NmbVer, sizeof(VerSct), CmpFnc);
   }

   msh.Old2New = malloc( (msh.NmbVer+1) * sizeof(int) );

   for(i=1;i<=msh.NmbVer;i++)
      msh.Old2New[ msh.ver[i].idx ] = i;

   GetTim(&timer);
   printf("%g s\n", timer);


   // --------------------
   // Elements renumbering
   // --------------------

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


   // --------------------------------------
   // Compute dependencies after renumbering
   // --------------------------------------
   if(StaFlg)
   {
      puts("\nDependencies after renumbering (average / MAX) :");
      PrtSta(&msh, LibParIdx);
   }

   StopParallel(LibParIdx);


   // ------------
   // Mesh writing
   // ------------

   printf("Writing mesh          : ");
   timer = 0.;
   GetTim(&timer);
   RecMsh(OutNam, &msh);
   GetTim(&timer);
   printf("%g s\n\n", timer);


   // --------------
   // Release memory
   // --------------

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
   int64_t  i;
   int      j;

   for(i=BegIdx; i<=EndIdx; i++)
   {
      msh->ver[i].idx = (itg)i;

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
/* Open, allocate and read the mesh fields                                    */
/*----------------------------------------------------------------------------*/

void ScaMsh(char *InpNam, MshSct *msh)
{
   int      i, j, n, t, *PtrIdx;
   int64_t  InpMsh;

   // Check mesh format
   if(!(InpMsh = GmfOpenMesh(InpNam, GmfRead, &msh->MshVer, &msh->dim)))
   {
      printf("Cannot open mesh %s\n", InpNam);
      exit(1);
   }

   // Get stats and allocate tables
   if((msh->NmbVer = (int)GmfStatKwd(InpMsh, GmfVertices)))
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

      // normalize the bounding box to map the geometry on a 64 bit cube
      for(j=0;j<3;j++)
         msh->box[j+3] = pow(2,64) / (msh->box[j+3] - msh->box[j]);
   }

   // Allocate and read elements
   for(t=0;t<MAXELE;t++)
   {
      n = msh->NmbEle[t] = (int)GmfStatKwd(InpMsh, EleTab[t][2]);

      if(!n)
         continue;

      msh->ele[t]    = malloc( (n+1) * sizeof(EleSct) );
      msh->IdxTab[t] = malloc( (n+1) * EleTab[t][0] * sizeof(int));
      PtrIdx         = msh->IdxTab[t];

      for(i=1;i<=n;i++)
         msh->ele[t][i].idx = &PtrIdx[ i * EleTab[t][0] ];

      GmfGetBlock(InpMsh,    EleTab[t][2], 1, n, 0, NULL, NULL,
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
   int      t, n;
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

      GmfSetKwd  (OutMsh,    EleTab[t][2], n);
      GmfSetBlock(OutMsh,    EleTab[t][2], 1, n, 0, NULL, NULL,
                  GmfIntVec, EleTab[t][0], msh->ele[t][1].idx,  msh->ele[t][n].idx,
                  GmfInt,                 &msh->ele[t][1].ref, &msh->ele[t][n].ref);
   }

#ifdef WITH_GMLIB
   if(msh->GmlMod)
   {
      // Save the number of domains
      GmfSetKwd(OutMsh, GmfDomains, 1);
      GmfSetLin(OutMsh, GmfDomains, 0, msh->MtsNmbDom);

      // Save the vertices Global ID / Domain / Local ID
      GmfSetKwd(OutMsh, GmfVerticesGID, msh->NmbVer);
      GmfSetBlock(OutMsh, GmfVerticesGID, 1, msh->NmbVer, 0, NULL, NULL,
                  GmfInt, &msh->ver[1].gid, &msh->ver[ msh->NmbVer ].gid,
                  GmfInt, &msh->ver[1].dom, &msh->ver[ msh->NmbVer ].dom,
                  GmfInt, &msh->ver[1].lid, &msh->ver[ msh->NmbVer ].lid);

      // Save the tets Global ID / Domain / Local ID
      GmfSetKwd(OutMsh, GmfTetrahedraGID, msh->NmbEle[3]);
      GmfSetBlock(OutMsh, GmfTetrahedraGID, 1, msh->NmbEle[3], 0, NULL, NULL,
                  GmfInt, &msh->ele[3][1].gid, &msh->ele[3][ msh->NmbEle[3] ].gid,
                  GmfInt, &msh->ele[3][1].dom, &msh->ele[3][ msh->NmbEle[3] ].dom,
                  GmfInt, &msh->ele[3][1].lid, &msh->ele[3][ msh->NmbEle[3] ].lid);
   }

#endif

   GmfCloseMesh(OutMsh);
}


/*----------------------------------------------------------------------------*/
/* Compute the hilbert code from 3d coordinates                               */
/*----------------------------------------------------------------------------*/

uint64_t HilCod(double crd[3], double box[6], int itr, int mod)
{
   uint64_t IntCrd[3], m=1LL<<63, cod;
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

   // Set the vertex code with the hilbert code, shifted by 11 bits in the case of
   // GMlib mode to make room for the high degree bit and the 10 bit ref hash tag
#ifdef WITH_GMLIB
   if(msh->GmlMod)
      for(i=BegIdx; i<=EndIdx; i++)
         msh->ver[i].cod |= HilCod(msh->ver[i].crd, msh->box, MAXITR, msh->mod) >> 11;
   else
#endif
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

   // Set the elements code with their barycenter hilbert code, shifted
   // by 10 bits in GMlib mode to make room for the 10 bits ref hash tag
#ifdef WITH_GMLIB
   if(msh->GmlMod)
      for(i=BegIdx; i<=EndIdx; i++)
      {
         SetNewIdx(EleTab[ msh->TypIdx ][0], msh->ele[ msh->TypIdx ][i].idx, msh->Old2New);
         SetMidCrd(EleTab[ msh->TypIdx ][0], msh->ele[ msh->TypIdx ][i].idx, msh, crd);
         msh->ele[ msh->TypIdx ][i].cod |= HilCod(crd, msh->box, MAXITR, msh->mod) >> 10;
      }
   else
#endif
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
   int i, *IdxTab, siz = EleTab[ typ ][0];

   // Allocate the new table
   IdxTab = malloc( (msh->NmbEle[ typ ] + 1) * siz * sizeof(int));

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
/* Setup packed neighbours + voyeurs information                              */
/*----------------------------------------------------------------------------*/

char *SetNgb(MshSct *msh)
{
   int      i, j, k, t, idx[3], key, NmbCol = 0, siz, gid = 1, (*ngb)[6];
   char     (*voy)[6], *BndTab;
   BucSct   *hsh, *col, *buc, *cur = NULL;
   EleSct   *ele;

   // Compute the total number of faces
   siz = 4 * msh->NmbEle[ TypTet ]
       + 5 * msh->NmbEle[ TypPyr ]
       + 5 * msh->NmbEle[ TypPri ]
       + 6 * msh->NmbEle[ TypHex ];

   // Allocate but do not clear a colision table for each face
   col = malloc(siz * sizeof(BucSct));
   assert(col);

   // Allocate and clear the base hash table with half the number of faces
   siz /= 2;
   hsh = calloc(siz, sizeof(BucSct));
   assert(hsh);

   // Set each element's global ID
   for(t=TypTet; t<=TypHex; t++)
      for(i=1;i<=msh->NmbEle[t];i++)
         msh->ele[t][i].gid = gid++;

   // Allocate the neighbours and a voyeurs tables
   ngb = calloc(gid + 1, 6 * sizeof(int));
   assert(ngb);
   voy = calloc(gid + 1, 6 * sizeof(char));
   assert(voy);

   // Scan each element type
   for(t=TypTet; t<=TypHex; t++)
   {
      // Scan each element
      for(i=1;i<=msh->NmbEle[t];i++)
      {
         ele = &msh->ele[t][i];

         // Scan each face
         for(j=0;j<EleTab[t][1];j++)
         {
            // Get the face's caracteristic node indices
            SrtFac(ele, t, j, idx);
            key = (3 * idx[0] + 5 * idx[1] + 7 * idx[2]) % siz;
            buc = &hsh[ key ];

            // If the bucket is empty, store the face
            if(!buc->ele)
            {
               for(k=0;k<3;k++)
                  buc->idx[k] = idx[k];

               buc->ele = ele->gid;
               buc->voy = j;
               buc->typ = t;
            }
            else
            {
               // Search for the face through the linked list
               do
               {
                  // If found, setup the neighborhood information
                  if((idx[0] == buc->idx[0])
                  && (idx[1] == buc->idx[1])
                  && (idx[2] == buc->idx[2]))
                  {
                     ngb[ ele->gid ][j] = buc->ele;
                     ngb[ buc->ele ][ buc->voy ] = ele->gid;
                     break;
                  }

                  cur = buc;
               }while((buc = buc->nex));

               // If not, add the face to the collision buffer
               if(!buc)
               {
                  buc = cur->nex = &col[ NmbCol++ ];
                  buc->ele = ele->gid;
                  buc->voy = j;
                  buc->typ = t;
                  buc->nex = NULL;

                  for(k=0;k<3;k++)
                     buc->idx[k] = idx[k];
               }
            }
         }
      }
   }

   free(hsh);
   free(col);

   // Allocate and fill a boundary vertex flag table and return it
   BndTab = calloc(msh->NmbVer + 1, sizeof(char));
   assert(BndTab);

   for(t=TypTet; t<=TypHex; t++)
      for(i=1;i<=msh->NmbEle[t];i++)
         for(j=0;j<EleTab[t][1];j++)
            if(!ngb[ msh->ele[t][i].gid ][j])
               for(k=0;k<FacDeg[t][j];k++)
                  BndTab[ msh->ele[t][i].idx[ EleFac[t][j][k] ] ] = 1;

   free(ngb);
   free(voy);

   return(BndTab);
}


/*----------------------------------------------------------------------------*/
/* Fill idx table with the sorted face vertices                               */
/*----------------------------------------------------------------------------*/

void SrtFac(EleSct *ele, int typ, int fac, int nod[3])
{
   int a, b, c, i, MinPos;

   // Triangle face: return the three sorted indices
   if(FacDeg[ typ ][ fac ] == 3)
   {
      a = ele->idx[ EleFac[ typ ][ fac ][0] ];
      b = ele->idx[ EleFac[ typ ][ fac ][1] ];
      c = ele->idx[ EleFac[ typ ][ fac ][2] ];

      if(a < b)
         if(b < c)
         {
            nod[0] = a;
            nod[1] = b;
            nod[2] = c;
         }
         else
            if(a < c)
            {
               nod[0] = a;
               nod[1] = c;
               nod[2] = b;
            }
            else
            {
               nod[0] = c;
               nod[1] = a;
               nod[2] = b;
            }
      else
         if(a < c)
         {
            nod[0] = b;
            nod[1] = a;
            nod[2] = c;
         }
         else
            if(b < c)
            {
               nod[0] = b;
               nod[1] = c;
               nod[2] = a;
            }
            else
            {
               nod[0] = c;
               nod[1] = b;
               nod[2] = a;
            }
   }
   else if(FacDeg[ typ ][ fac ] == 4)
   {
      // Quad face: return the minimum index vertex and its diagonal opposite
      MinPos = 0;

      for(i=1;i<4;i++)
         if(ele->idx[ EleFac[ typ ][ fac ][i] ] < ele->idx[ EleFac[ typ ][ fac ][ MinPos ] ])
            MinPos = i;

      nod[0] = ele->idx[ EleFac[ typ ][ fac ][ MinPos ] ];
      nod[1] = ele->idx[ EleFac[ typ ][ fac ][ (MinPos + 2) % 4 ] ];
      nod[2] = 0;
   }
}


#ifdef WITH_GMLIB

/*----------------------------------------------------------------------------*/
/* Read the nodes and elements domain from the Metis files                    */
/*----------------------------------------------------------------------------*/

void ScaMts(char *MtsNodNam, char *MtsEleNam, MshSct *msh)
{
   int      i, ref, dom, VerCpt[ MAXDOM ], TetCpt[ MAXDOM ];
   FILE     *mts;
   VerSct   *ver;
   EleSct   *tet;

   // Read the Metis vertex domain and build the local idx / domain information
   printf("Reading node partitions from %s\n", MtsNodNam);
   if(!(mts = fopen(MtsNodNam, "r")))
      return;

   for(i=0;i<msh->MtsNmbDom;i++)
      VerCpt[i] = 0;

   // Deduce the local ID from the Global ID and the domain ID
   for(i=1;i<=msh->NmbVer;i++)
   {
      ver = &msh->ver[i];
      fscanf(mts, "%d", &dom);
      dom = MIN(dom, msh->MtsNmbDom - 1);
      dom = MAX(dom, 0);
      VerCpt[ dom ]++;
      ver->gid = i;
      ver->dom = dom;
      ver->lid = VerCpt[ dom ];
   }

   fclose(mts);

   printf("Reading element %d partitions from %s\n", msh->MtsEleTyp, MtsNodNam);
   if(!(mts = fopen(MtsEleNam, "r")))
      return;

   if( (msh->MtsEleTyp != GmfTetrahedra) || !msh->NmbEle[3] )
      return;

   for(i=1;i<=msh->MtsNmbDom;i++)
      TetCpt[i] = 0;

   // Deduce the local ID from the Global ID and the domain ID
   for(i=1;i<=msh->NmbEle[3];i++)
   {
      tet = &msh->ele[3][i];
      fscanf(mts, "%d", &dom);
      dom = MIN(dom, msh->MtsNmbDom - 1);
      dom = MAX(dom, 0);
      TetCpt[ dom ]++;
      tet->gid = i;
      tet->dom = dom;
      tet->lid = TetCpt[ dom ];
   }

   fclose(mts);
}


/*----------------------------------------------------------------------------*/
/* Set the code 10 upper bit with a hash key based on element's ref           */
/*----------------------------------------------------------------------------*/

void GmlSrt(MshSct *msh)
{
   char *BytPtr;
   int i, j, t, (*DegTab)[ MAXELE ];
   uint64_t cod;
   VerSct *ver;
   EleSct *ele;

   // Allocate a degree table with one scalar per kind of element
   DegTab = calloc( msh->NmbVer + 1, MAXELE * sizeof(int));

   // Add each element's vertices to the degree associated to its kind
   // And compute each elements hash tag based on their ref
   for(t=0;t<MAXELE;t++)
   {
      if(!msh->EleTyp[t])
         continue;

      for(i=1;i<=msh->NmbEle[t];i++)
      {
         // Compute the element ref hash tag as bits 54 to 63
         ele = &msh->ele[t][i];
         BytPtr = (char *)&ele->ref;
         cod = BytPtr[0] + BytPtr[1] + BytPtr[2] + BytPtr[3];
         ele->cod = cod << 54;

         // Increment the vertex degrees
         for(j=0;j<EleTab[t][0];j++)
            DegTab[ ele->idx[j] ][t]++;
      }
   }

   // Compute the same ref hash tag for each vertices as bits 53 to 62
   // and set a special high degree flag as leftmost bit (63rd)
   for(i=1;i<=msh->NmbVer;i++)
   {
      ver = &msh->ver[i];
      BytPtr = (char *)&ver->ref;
      cod = BytPtr[0] + BytPtr[1] + BytPtr[2] + BytPtr[3];
      ver->cod = cod << 53;

      if( (DegTab[i][1] > MaxDeg[1][0])
      ||  (DegTab[i][2] > MaxDeg[2][0])
      ||  (DegTab[i][3] > MaxDeg[3][0])
      ||  (DegTab[i][6] > MaxDeg[6][0]) )
      {
         ver->cod |= 1L<<63;
         msh->HghDeg++;

         for(j=0;j<MAXELE;j++)
            msh->MaxDeg[j] = MAX(msh->MaxDeg[j], DegTab[i][j]);
      }

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
}
#endif
