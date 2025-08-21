

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


#ifndef LPLIB4_HELPERS_H
#define LPLIB4_HELPERS_H


/*----------------------------------------------------------------------------*/
/* Defines                                                                    */
/*----------------------------------------------------------------------------*/

#include <stdint.h>

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

#define LplGpu 1
typedef itg int56d[56];

enum LplTyp {LplVer, LplEdg, LplTri, LplQad, LplTet, LplPyr, LplPri, LplHex, LplMax};
enum RenTyp {LplNoRenum, LplHilbert, LplZcurve};


/*----------------------------------------------------------------------------*/
/* Prototypes of public structures                                            */
/*----------------------------------------------------------------------------*/

typedef struct
{
   int      dim, NmbVer, NmbEle[ LplMax ], *EleTab[ LplMax ];
   uint64_t (*RenTab[ LplMax ])[2];
   double   *CrdTab, box[6];
}LplRenSct;

#ifdef WITH_METIS
typedef struct
{
   int      NmbVer, NmbEdg, NmbTet;
   int      *VerDeg, *VerBal, *EdgTab, *RefTab, (*TetTab)[4];
   int      (*BalTab)[2], *TetRef, *VerDegRk2;
   int      **VerBalRk2;
   double   MinSiz, MaxSiz, (*CrdTab)[3];
}LplMshSct;
#endif


/*----------------------------------------------------------------------------*/
/* Prototypes of public procedures                                            */
/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#endif

itg      ParallelBuildEdges   (itg, int, itg *, itg **);
LplRenSct  *MeshRenumbering      (int, int, int, int, ...);
void     FreeNumberingStruct  (LplRenSct *);
double   EvaluateRenumbering  (int, int, int *);
int      RestoreNumbering     (LplRenSct *, int, double *, ...);
int      GetOldIndex          (LplRenSct *, int, int);
int      GetNewIndex          (LplRenSct *, int, int);

#ifdef WITH_METIS
int      MetisPartitioning    (LplMshSct *, int);
#endif
   
#ifdef __cplusplus
} // end extern "C"
#endif

#endif
