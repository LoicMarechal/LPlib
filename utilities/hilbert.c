

/*----------------------------------------------------------------------------*/
/*                                                                            */
/*                               HILBERT V3.11                                */
/*                                                                            */
/*----------------------------------------------------------------------------*/
/*                                                                            */
/* Description:         renumber .mesh(b) files                               */
/* Author:              Loic MARECHAL                                         */
/* Creation date:       mar 11 2010                                           */
/* Last modification:   dec 05 2024                                           */
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
#include <libmeshb7.h>
#include <libmeshb7_helpers.h>
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
enum {HilMod=0, OctMod, RndMod, IniMod, TopMod};


/*----------------------------------------------------------------------------*/
/* Macros                                                                     */
/*----------------------------------------------------------------------------*/

#define MIN(a,b)  ((a) < (b) ? (a) : (b))
#define MAX(a,b)  ((a) > (b) ? (a) : (b))
#define POW(a)    ((a) * (a))


/*----------------------------------------------------------------------------*/
/* Structures                                                                 */
/*----------------------------------------------------------------------------*/

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

typedef struct PtrBuc
{
   int      idx[3], ele;
   char     typ, voy;
   struct   PtrBuc *nex;
}BucSct;

typedef struct
{
   int MinIdx, MaxIdx, NexBuc;
}HshSct;

typedef struct
{
   int      NmbVer, *Old2New, MshVer, dim, mod, TypIdx, VerTyp, GmlMod;
   int      NmbEle[ MAXELE ], *IdxTab[ MAXELE ], EleTyp[ MAXELE ];
   int      MaxDeg[ MAXELE ], DegVec[ MAXELE ], HghDeg, OvrDeg, ColGrnMsh;
   int      ColGrnMod, NmbGrnPar, NmbTypGrnPar[ MAXELE ], (*GrnPar)[ MAXELE ][4];
   int      NmbColPar, NmbTypColPar[ MAXELE ], (*ColPar)[ MAXELE ][3];
   int      ColBit, GrnBit, DegBit, RefBit, VerHilBit, FacHilBit, VolHilBit;
   int      ColLft, GrnLft, DegLft, RefLft, VerHilRgt, FacHilRgt, VolHilRgt;
   uint64_t ColMsk, GrnMsk, DegMsk, RefMsk;
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

const int tvpe[6][2] = { {0,1}, {1,2}, {2,0}, {3,0}, {3,1}, {3,2} };


/*----------------------------------------------------------------------------*/
/* Prototypes of local procedures                                             */
/*----------------------------------------------------------------------------*/

void     ScaMsh   (char *, MshSct *);
void     RecMsh   (char *, MshSct *);
uint64_t HilCod   (double *, double *, int, int);
uint64_t IntHilCod(int *, int, int);
void     GetTim   (double *);
int      CmpFnc   (const void *, const void *);
void     RenVer   (int, int, int, MshSct *);
void     RenEle   (int, int, int, MshSct *);
void     PrtSta   (MshSct *, int64_t);
void     SwpMem   (MshSct *, int);
char    *SetNgb   (MshSct *);
void     SrtFac   (EleSct *, int, int, int [3]);
void     SetVerDeg(MshSct *);
void     SetMatSlc(MshSct *);
void     ScaMts   (char *, char *, MshSct *);
int      SetDeg   (MshSct *, int *);


/*----------------------------------------------------------------------------*/
/* Read, renumber through a Hilbert SFC and write the mesh                    */
/*----------------------------------------------------------------------------*/

int main(int ArgCnt, char **ArgVec)
{
   char     *PtrArg, *TmpStr, *BndTab, InpNam[1000], OutNam[1000];
   char     MtsNodNam[1000], MtsEleNam[1000];
   char     *SchStr[4] = {"Hilbert", "Z-curve", "random", "initial"};
   int      i, j, t, NmbCpu=0, StaFlg=0, BndFlg=0, *IdxTab, NmbSrf, NmbVol;
   int      NmbGrn, CurGrn, NmbCol, CurCol;
   int64_t  LibParIdx;
   double   timer = 0.;
   MshSct   msh;
   VerSct   *OldVer, *NewVer, *VolVer;


   // --------------------
   // Command line parsing
   // --------------------

   memset(&msh, 0, sizeof(MshSct));

   if(ArgCnt == 1)
   {
      puts("\nHILBERT v3.10 december 04 2024   Loic MARECHAL / INRIA\n");
      puts(" Usage         : hilbert -in input_mesh -out renumbered_mesh");
      puts("   -in name    : input mesh(b) name");
      puts("   -out name   : output renumbered mesh(b)");
      puts("   -stats      : print element blocks dependencies stats before and after renumbering");
      puts("   -fixbnd     : do not renumber boundary nodes");
      puts("   -nproc n    : n is the number of threads to be launched (default = all threads)\n");
      puts(" Sorting       : optional arguments to control the entities sorting (see bottom explanation)");
      puts("   -colors     : use color as rank4 and grain as rank3 if such fields are present in the input file");
      puts("   -gmlib type : special vertex sorting to make the mesh fit for the GMlib depending on type");
      puts("                 generic: set rank2 with high/low degree for vertices and reference for faces");
      puts("                 matrix : set rank2 with the matrix slice size for each vertex");
      puts("   -scheme s   : set rank1 with a value computed by a renumbering scheme");
      puts("                 0: geometrical Hilbert (default)");
      puts("                 1: Z curve (octree like numbering)");
      puts("                 2: random");
      puts("                 3: no sort (preserve initial numbering");
      puts("                 4: geometrical Hilbert for vertices and topological Hilbert for elements\n");
      puts(" All entities are sorted against four keys ranging from rank 4 (highest) to 1 (lowest)");
      puts(" rank4: color, rank3: grain, rank2: vertex degree or face ref, rank1: local scheme");
      puts(" all ranks are optional and cam be controled by the above arguments\n");
      exit(0);
   }

   for(i=2;i<=ArgCnt;i++)
   {
      PtrArg = *++ArgVec;

      // Get input file name and add the .meshb extension if unspecified
      if(!strcmp(PtrArg,"-in"))
      {
         TmpStr = *++ArgVec;
         ArgCnt--;
         strcpy(InpNam, TmpStr);

         if(!strstr(InpNam, ".mesh"))
            strcat(InpNam, ".meshb");

         continue;
      }

      // Get output file name and add the .meshb extension if unspecified
      if(!strcmp(PtrArg,"-out"))
      {
         TmpStr = *++ArgVec;
         ArgCnt--;
         strcpy(OutNam, TmpStr);

         if(!strstr(OutNam, ".mesh"))
            strcat(OutNam, ".meshb");

         continue;
      }

      // Print statistics about blocks collisions before and after renumbering
      if(!strcmp(PtrArg,"-stats"))
      {
         StaFlg = 1;
         continue;
      }

      // Do not renumber boundary nodes and elements
      if(!strcmp(PtrArg,"-fixbnd"))
      {
         BndFlg = 1;
         continue;
      }

      // Use the color and grain information stored at vertices if present in the file
      // and use those two fields as primary and secondary keys to the sorting step
      if(!strcmp(PtrArg,"-colors"))
      {
          msh.ColGrnMod = 1;
         continue;
      }

      // Compute the vertex degrees and sort them against it
      // ball mode is for general purpose high/low degree vertex sorting
      // matrix mode slices the matrix in 5 slabs of 16,32,64,128,256 size
      if(!strcmp(PtrArg,"-gmlib"))
      {
         TmpStr = *++ArgVec;
         ArgCnt--;

         if(!strcmp(TmpStr,"generic"))
            msh.GmlMod = 1;
         else if(!strcmp(TmpStr,"matrix"))
            msh.GmlMod = 2;
            
         continue;
      }

      // Select the renumbering scheme: default if Hibert
      // but Z-curve and random can be specified
      if(!strcmp(PtrArg,"-scheme"))
      {
         msh.mod = atoi(*++ArgVec);
         msh.mod = MAX(msh.mod, 0);
         msh.mod = MIN(msh.mod, 4);
         ArgCnt--;
         continue;
      }

      // Set the number of threads used by the renumbering step
      // 0 uses all available system's cores
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


   // ------------
   // Mesh reading
   // ------------

   printf("\nReading mesh                 : ");
   timer = 0.;
   GetTim(&timer);
   ScaMsh(InpNam, &msh);
   GetTim(&timer);
   printf("%g s\n", timer);

   printf("Input mesh version           : %d\n", msh.MshVer);
   printf("Vertices                     : %d\n", msh.NmbVer);

   for(t=0;t<MAXELE;t++)
      if(msh.NmbEle[t])
         printf("%s             : %d\n", EleNam[t], msh.NmbEle[t]);

   puts("");

   // If the color renumbering is set but no colors or grains are present in the
   // input file, the mode is disabled and set to default hilbert renumbering
   if(msh.ColGrnMod && !msh.ColGrnMsh)
   {
      msh.ColGrnMod = 0;
      puts("Could not find colors and grains information: switching back to default renumbering");
   }


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


   // ---------------------------------------------------------
   // Set the size of bit fields needed by color and grain data
   // ---------------------------------------------------------

   if(msh.ColGrnMod)
   {
      msh.NmbColPar = msh.NmbGrnPar = 0;

      for(i=1;i<=msh.NmbVer;i++)
      {
         msh.NmbColPar = MAX(msh.NmbColPar, msh.ver[i].col);
         msh.NmbGrnPar = MAX(msh.NmbGrnPar, msh.ver[i].grn);
      }

      printf("Found %d colors and %d grains in the input file\n", msh.NmbColPar, msh.NmbGrnPar);

      // Derive element's colors and grains from the vertices
      for(t=0;t<MAXELE;t++)
      {
         for(i=1;i<=msh.NmbEle[t];i++)
         {
            msh.ele[t][i].col = msh.ver[ msh.ele[t][i].idx[0] ].col;
            msh.ele[t][i].grn = msh.ver[ msh.ele[t][i].idx[0] ].grn;
         }
      }

      // Set the bit field size, maks and left shift to store the color value
      msh.ColBit = ceil(log(msh.NmbColPar) / log(2));
      msh.ColMsk = (1ULL << msh.ColBit) - 1ULL;
      msh.ColLft = 64 - msh.ColBit;

      // Set the bit field size, maks and left shift to store the grain value
      msh.GrnBit = ceil(log(msh.NmbGrnPar) / log(2));
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

   printf("\nRenumbering vertices         : ");
   timer = 0.;
   GetTim(&timer);

   if(BndFlg)
   {
      // Set the neighbours between elements and set boundary vertices tags
      BndTab = SetNgb(&msh);

      // Extract the inner volume vertices in a separate table
      IdxTab = malloc( (size_t)(msh.NmbVer + 1) * sizeof(int) );
      assert(IdxTab);
      NmbVol = 0;

      for(i=1;i<=msh.NmbVer;i++)
         IdxTab[i] = BndTab[i] ? 0 : ++NmbVol;

      VolVer = malloc( (size_t)(NmbVol + 1) * sizeof(VerSct) );
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
      NewVer = malloc( (size_t)(msh.NmbVer + 1) * sizeof(VerSct) );
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

   msh.Old2New = malloc( (size_t)(msh.NmbVer+1) * sizeof(int) );
   assert(msh.Old2New);

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

   //

   if(msh.ColGrnMod)
   {
      msh.ColPar = malloc( (msh.NmbColPar + 1) * MAXELE * 3 * sizeof(int) );
      assert(msh.ColPar);

      msh.GrnPar = malloc( (msh.NmbGrnPar + 1) * MAXELE * 4 * sizeof(int) );
      assert(msh.GrnPar);

      CurGrn = msh.ver[1].grn;
      NmbGrn = 1;
      msh.GrnPar[ NmbGrn ][0][0] = 1;
      msh.GrnPar[ NmbGrn ][0][2] = msh.ver[1].col;
      msh.GrnPar[ NmbGrn ][0][3] = msh.ver[1].grn;

      CurCol = msh.ver[1].col;
      NmbCol = 1;
      msh.ColPar[ NmbCol ][0][0] = 1;
      msh.ColPar[ NmbCol ][0][2] = NmbGrn;

      for(i=1;i<msh.NmbVer;i++)
      {
         if(msh.ver[i].grn != CurGrn)
         {
            msh.GrnPar[ NmbGrn ][0][1] = i - 1;
            NmbGrn++;
            msh.GrnPar[ NmbGrn ][0][0] = i;
            CurGrn = msh.ver[i].grn;
            msh.GrnPar[ NmbGrn ][0][2] = msh.ver[i].col;
            msh.GrnPar[ NmbGrn ][0][3] = msh.ver[i].grn;

            if(msh.ver[i].col != CurCol)
            {
               msh.ColPar[ NmbCol ][0][1] = NmbGrn - 1;
               NmbCol++;
               msh.ColPar[ NmbCol ][0][0] = NmbGrn;
               CurCol = msh.ver[i].col;
               msh.ColPar[ NmbCol ][0][2] = msh.ver[i].col;
            }
         }
      }

      msh.GrnPar[ NmbGrn ][0][1] = msh.NmbVer;
      msh.ColPar[ NmbCol ][0][1] = NmbGrn;
      msh.NmbTypGrnPar[0] = NmbGrn;
      msh.NmbTypColPar[0] = NmbCol;

      for(i=1;i<=NmbGrn;i++)
         printf("vertex grain %3d (%3d/%3d): %8d -> %8d, size: %8d\n",
                  i, msh.GrnPar[i][0][2], msh.GrnPar[i][0][3],
                  msh.GrnPar[i][0][0], msh.GrnPar[i][0][1],
                  msh.GrnPar[i][0][1] - msh.GrnPar[i][0][0] + 1);

      for(i=1;i<=NmbCol;i++)
         printf("vertex color %3d (%3d): %8d -> %8d, size: %8d\n",
                  i, msh.ColPar[i][0][2],
                  msh.ColPar[i][0][0], msh.ColPar[i][0][1],
                  msh.ColPar[i][0][1] - msh.ColPar[i][0][0] + 1);

      for(t=0;t<MAXELE;t++)
      {
         if(!msh.NmbEle[t])
            continue;

         CurGrn = msh.ele[t][1].grn;
         NmbGrn = 1;
         msh.GrnPar[ NmbGrn ][t][0] = 1;
         msh.GrnPar[ NmbGrn ][t][2] = msh.ele[t][1].col;
         msh.GrnPar[ NmbGrn ][t][3] = msh.ele[t][1].grn;

         CurCol = msh.ele[t][1].col;
         NmbCol = 1;
         msh.ColPar[ NmbCol ][0][0] = 1;
         msh.ColPar[ NmbCol ][0][2] = NmbGrn;

         for(i=1;i<msh.NmbEle[t];i++)
         {
            if(msh.ele[t][i].grn != CurGrn)
            {
               msh.GrnPar[ NmbGrn ][t][1] = i - 1;
               NmbGrn++;
               msh.GrnPar[ NmbGrn ][t][0] = i;
               CurGrn = msh.ele[t][i].grn;
               msh.GrnPar[ NmbGrn ][t][2] = msh.ele[t][i].col;
               msh.GrnPar[ NmbGrn ][t][3] = msh.ele[t][i].grn;

               if(msh.ele[t][i].col != CurCol)
               {
                  msh.ColPar[ NmbCol ][t][1] = NmbGrn - 1;
                  NmbCol++;
                  msh.ColPar[ NmbCol ][t][0] = NmbGrn;
                  CurCol = msh.ele[t][i].col;
                  msh.ColPar[ NmbCol ][t][2] = msh.ele[t][i].col;
               }
            }
         }

         msh.GrnPar[ NmbGrn ][t][1] = msh.NmbEle[t];
         msh.ColPar[ NmbCol ][t][1] = NmbGrn;
         msh.NmbTypGrnPar[t] = NmbGrn;
         msh.NmbTypColPar[t] = NmbCol;
         /*
         for(i=1;i<=NmbGrn;i++)
            printf(  "%s grain %3d (%3d/%3d): %8d -> %8d, size: %8d\n",
                     EleNam[t], i,
                     msh.GrnPar[i][t][2], msh.GrnPar[i][t][3],
                     msh.GrnPar[i][t][0], msh.GrnPar[i][t][1],
                     msh.GrnPar[i][t][1] - msh.GrnPar[i][t][0] + 1 );

        for(i=1;i<=NmbCol;i++)
           printf("%s color %3d (%3d): %8d -> %8d, size: %8d\n",
                    EleNam[t], i, msh.ColPar[i][t][2],
                    msh.ColPar[i][t][0], msh.ColPar[i][t][1],
                    msh.ColPar[i][t][1] - msh.ColPar[i][t][0] + 1);*/
      }
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

   printf("Writing mesh                 : ");
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
   {
      msh->ver = malloc( (size_t)(msh->NmbVer+1) * sizeof(VerSct));
      assert(msh->ver);
   }
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

      if((GmfStatKwd(InpMsh, GmfVerticesColour) == msh->NmbVer)
      && (GmfStatKwd(InpMsh, GmfVerticesGrain)  == msh->NmbVer))
      {
         msh->ColGrnMsh = 1;
         GmfGetBlock(InpMsh, GmfVerticesColour, 1, msh->NmbVer, 0, NULL, NULL,
                     GmfInt, &msh->ver[1].col, &msh->ver[ msh->NmbVer ].col);

         GmfGetBlock(InpMsh, GmfVerticesGrain, 1, msh->NmbVer, 0, NULL, NULL,
                     GmfInt, &msh->ver[1].grn, &msh->ver[ msh->NmbVer ].grn);
      }
   }

   // Allocate and read elements
   for(t=0;t<MAXELE;t++)
   {
      n = msh->NmbEle[t] = (int)GmfStatKwd(InpMsh, EleTab[t][2]);

      if(!n)
         continue;

      msh->ele[t]    = malloc( (size_t)(n+1) * sizeof(EleSct) );
      assert(msh->ele[t]);

      msh->IdxTab[t] = malloc( (size_t)(n+1) * (size_t)EleTab[t][0] * sizeof(int));
      assert(msh->IdxTab[t]);

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

   if(msh->ColGrnMod)
   {
      GmfSetKwd(OutMsh, GmfVertexGrainPartitions, msh->NmbTypGrnPar[0]);
      GmfSetBlock(OutMsh, GmfVertexGrainPartitions, 1, msh->NmbTypGrnPar[0], 0, NULL, NULL,
                  GmfIntVec, 2, msh->GrnPar[1][0], msh->GrnPar[ msh->NmbTypGrnPar[0] ][0]);

      if(msh->NmbEle[ TypTet ])
      {
         GmfSetKwd(OutMsh, GmfTetrahedronGrainPartitions, msh->NmbTypGrnPar[ TypTet ]);
         GmfSetBlock(OutMsh, GmfTetrahedronGrainPartitions, 1,
                     msh->NmbTypGrnPar[ TypTet ], 0, NULL, NULL,
                     GmfIntVec, 2, msh->GrnPar[1][ TypTet ],
                     msh->GrnPar[ msh->NmbTypGrnPar[ TypTet ] ][ TypTet ]);
      }

      GmfSetKwd(OutMsh, GmfColorPartitions, msh->NmbColPar);
      GmfSetBlock(OutMsh, GmfColorPartitions, 1, msh->NmbColPar, 0, NULL, NULL,
                  GmfIntVec, 2, msh->ColPar[1], msh->ColPar[ msh->NmbColPar ]);
   }

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
/* Compute the topological Hilbert code from a pair of node indices           */
/*----------------------------------------------------------------------------*/

uint64_t IntHilCod(int *EleCrd, int NmbVer, int dim)
{
   uint64_t IntCrd[2], m=1LL<<63, cod, min = 1LL<<63, max = 0;
   int      lft = 63 - log2(NmbVer);
   int      j, b, GeoWrd, NewWrd, BitTab[3] = {1,2,4};
   int      rot[8], GeoCod[8]={0,3,7,4,1,2,6,5};
   int      HilCod[8][8] = {
            {0,7,6,1,2,5,4,3}, {0,3,4,7,6,5,2,1},
            {0,3,4,7,6,5,2,1}, {2,3,0,1,6,7,4,5},
            {2,3,0,1,6,7,4,5}, {6,5,2,1,0,3,4,7},
            {6,5,2,1,0,3,4,7}, {4,3,2,5,6,1,0,7} };

   // Binary hilbert renumbering loop
   cod = 0;

   // Initialize IntCrd[2] with the min and max element's node indices
   for(j=0;j<dim;j++)
   {
      min = MIN(min, EleCrd[j]);
      max = MAX(max, EleCrd[j]);
   }

   IntCrd[0] = min << 32;
   IntCrd[1] = max << 32;

   for(j=0;j<8;j++)
      rot[j] = GeoCod[j];

   // Interleave the pair of 32-bit coordinates into a single 64-bit Hilbert code
   for(b=0;b<32;b++)
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

      for(j=0;j<8;j++)
         rot[j] = HilCod[ NewWrd ][ rot[j] ];
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
      else if(msh->mod == TopMod)
         SchCod = IntHilCod(msh->ele[ msh->TypIdx ][i].idx, msh->NmbVer, EleTab[ msh->TypIdx ][0]);
      else
      {
         SetMidCrd(EleTab[ msh->TypIdx ][0], msh->ele[ msh->TypIdx ][i].idx, msh, crd);
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
   col = malloc( (size_t)siz * sizeof(BucSct));
   assert(col);

   // Allocate and clear the base hash table with half the number of faces
   siz /= 2;
   hsh = calloc( (size_t)siz, sizeof(BucSct));
   assert(hsh);

   // Set each element's global ID
   for(t=TypTet; t<=TypHex; t++)
      for(i=1;i<=msh->NmbEle[t];i++)
         msh->ele[t][i].gid = gid++;

   // Allocate the neighbours and a voyeurs tables
   ngb = calloc( (size_t)(gid + 1), 6 * sizeof(int));
   assert(ngb);
   voy = calloc( (size_t)(gid + 1), 6 * sizeof(char));
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
   BndTab = calloc( (size_t)(msh->NmbVer + 1), sizeof(char));
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
