

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

enum LplTyp {LplVer, LplEdg, LplTri, LplQad, LplTet, LplPyr, LplPri, LplHex, LplMax};


/*----------------------------------------------------------------------------*/
/* Prototypes of public structures                                            */
/*----------------------------------------------------------------------------*/

typedef struct
{
   int      dim, NmbVer, NmbEle[ LplMax ], *EleTab[ LplMax ];
   uint64_t (*RenTab[ LplMax ])[2];
   double   *CrdTab, box[6];
}RenSct;


/*----------------------------------------------------------------------------*/
/* Prototypes of public procedures                                            */
/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#endif

itg      ParallelBuildEdges   (itg, int, itg *, itg **);
RenSct  *MeshRenumbering      (int, int, double *, ...);
void     FreeNumberingStruct  (RenSct *);
double   EvaluateRenumbering  (int, int, int *);
int      RestoreNumbering     (RenSct *, int, double *, ...);
int      GetOldIndex          (RenSct *, int, int);
int      GetNewIndex          (RenSct *, int, int);

#ifdef __cplusplus
} // end extern "C"
#endif

#endif
