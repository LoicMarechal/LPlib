

/*----------------------------------------------------------------------------*/
/*                                                                            */
/*                               LPlib Helpers V0.1                           */
/*                                                                            */
/*----------------------------------------------------------------------------*/
/*                                                                            */
/* Description:         lplib's helper functions' headers                     */
/* Author:              Loic MARECHAL                                         */
/* Creation date:       may 16 2024                                           */
/* Last modification:   oct 31 2024                                           */
/*                                                                            */
/*----------------------------------------------------------------------------*/


#ifndef LPLIB3_HELPERS_H
#define LPLIB3_HELPERS_H


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

enum LplTyp {LplVer, LplEdg, LplTri, LplQad, LplTet, LplPyr, LplPri, LplHex};


/*----------------------------------------------------------------------------*/
/* Prototypes of public structures                                            */
/*----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------*/
/* Prototypes of public procedures                                            */
/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#endif

itg ParallelBuildEdges(itg, int, itg *, itg **);

#ifdef __cplusplus
} // end extern "C"
#endif

#endif
