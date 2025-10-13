#ifndef _LPLIB_H
#define _LPLIB_H


/*----------------------------------------------------------------------------*/
/*                                                                            */
/*                               LPlib V4.12                                  */
/*                                                                            */
/*----------------------------------------------------------------------------*/
/*                                                                            */
/*   Description:       Handles threads, scheduling, pipelines & dependencies */
/*   Author:            Loic MARECHAL                                         */
/*   Creation date:     feb 25 2008                                           */
/*   Last modification: oct 03 2025                                           */
/*                                                                            */
/*----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------*/
/* Types definitions                                                          */
/*----------------------------------------------------------------------------*/

#include <stdint.h>
#include <stdlib.h>
#include <stddef.h>

#ifdef INT64
#define itg int64_t
#else
#define itg int32_t
#endif

#ifdef __cplusplus
extern "C" {
#endif


/*----------------------------------------------------------------------------*/
/* Defines                                                                    */
/*----------------------------------------------------------------------------*/

#define LplGpu 1

enum LplTyp {  LplVer, LplEdg, LplTri, LplQad, LplTet, LplPyr, LplPri, LplHex,
               LplEdgP2, LplTriP2, LplQadQ2, LplTetP2, LplPyrP2, LplPriP2,
               LplHexQ2, LplMax };

enum RenTyp {LplNoRenum, LplHilbert, LplZcurve};


/*----------------------------------------------------------------------------*/
/* Structures                                                                 */
/*----------------------------------------------------------------------------*/

typedef struct
{
   int      *Old2New, MshVer, dim, mod, TypIdx, VerTyp, VrbLvl;
   int      EleTyp[ LplMax ], GmlMod, *EleTab[ LplMax ];
   int      MaxDeg[ LplMax ], DegVec[ LplMax ], HghDeg, OvrDeg;
   int      ColGrnFlg, ColGrnMod, NmbCol, NmbGrn, *VerGrn, *VerCol;
   int      (*ColPar)[2], (*VerGrnPar)[2], (*EleGrnPar[ LplMax ])[2];
   int      ColBit, GrnBit, DegBit, RefBit, VerHilBit, FacHilBit, VolHilBit;
   int      ColLft, GrnLft, DegLft, RefLft, VerHilRgt, FacHilRgt, VolHilRgt;
   int      *VerDeg, *VerBal, *LstBalRk1, *RefTab;
   int      (*BalTab)[2], *EleCol[ LplMax ], *EleGrn[ LplMax ];
   int      *LstBalRk2, *VerRef, *EleRef[ LplMax ];
   uint64_t NmbVer, NmbEle[ LplMax ];
   uint64_t ColMsk, GrnMsk, DegMsk, RefMsk, *AdrBalRk1, *AdrBalRk2, TotDeg;
   uint64_t (*RenTab[ LplMax ])[2], (*VerCod)[2], (*EleCod[ LplMax ])[2];
   double   box[6];
   double   MinSiz, MaxSiz, *CrdTab;
}LplSct;


/*----------------------------------------------------------------------------*/
/* User available procedures' prototypes                                      */
/*----------------------------------------------------------------------------*/

int      AddDependency           (int64_t, itg, itg);
void     AddDependencyFast       (int64_t, int, itg *, int, itg *);
int      BeginDependency         (int64_t, int, int);
int      EndDependency           (int64_t, float [2]);
void     FreeType                (int64_t, int);
void     GetDependencyStats      (int64_t, int, int, float [2]);
void     GetLplibInformation     (int64_t, int *, int *);
int      GetNumberOfCores        ();
double   GetWallClock            ();
int      HilbertRenumbering      (int64_t, itg, double [6],
                                  double (*)[3], uint64_t (*)[2]);
int      HilbertRenumbering2D    (int64_t, itg, double [4],
                                  double (*)[2], uint64_t (*)[2]);
int64_t  InitParallel            (int);
int64_t  InitParallelAttr        (int, size_t, void *);
float    LaunchParallel          (int64_t, int, int, void *, void *);
float    LaunchParallelMultiArg  (int64_t, int, int, void *, int, ...);
int      LaunchPipeline          (int64_t, void *, void *, int, int *);
int      LaunchPipelineMultiArg  (int64_t, int, int *, void *prc, int, ...);
int      NewType                 (int64_t, itg);
int      ParallelMemClear        (int64_t, void *, size_t);
int      ParallelMemCopy         (int64_t, void *, void *, size_t);
void     ParallelQsort           (int64_t, void *, size_t, size_t, 
                                  int (*)(const void *, const void *));
int      ResizeType              (int64_t, int, itg);
void     StopParallel            (int64_t);
int      UpdateDependency        (int64_t, int, int, itg, itg);
void     UpdateDependencyFast    (int64_t, int, int, itg *, int, int, itg *);
void     WaitPipeline            (int64_t);
int      GetBigBlkNfo            (int64_t, int, int, int *, int *);
int      GetBlkIdx               (int64_t, int, int);
int      ChkBlkDep               (int64_t, int, int, int);
int      SetExtendedAttributes   (int64_t , ...);
int      HalveSmallBlocks        (int64_t, int, int);
int      HalveDependencyBlocks   (int64_t, int, int);
float    LaunchColorGrains       (int64_t, int, void *, void *);
LplSct  *MeshRenumbering         (int64_t, int, int, int, int, ...);
void     FreeNumberingStruct     (LplSct *);
double   EvaluateRenumbering     (int, int, int *);
int      RestoreNumbering        (LplSct *, int, double *, ...);
int      GetOldIndex             (LplSct *, int, int);
int      GetNewIndex             (LplSct *, int, int);

#ifdef __cplusplus
} // end extern "C"
#endif


/*----------------------------------------------------------------------------*/
/* Public defines                                                             */
/*----------------------------------------------------------------------------*/

#define MaxPth 256

enum ArgAtr {
   SetInterleavingFactor = 1,
   SetInterleavingSize,
   DisableInterleaving,
   EnableBlockSorting,
   DisableBlockSorting,
   StaticScheduling,
   SetSmallBlock,
   SetDependencyBlock,
   EnableAdaptiveSizing,
   DisableAdaptiveSizing
};


#endif  //-- define _LPLIB_H
