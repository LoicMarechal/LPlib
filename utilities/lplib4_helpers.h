

/*----------------------------------------------------------------------------*/
/*                                                                            */
/*                               LPlib Helpers V1.2                           */
/*                                                                            */
/*----------------------------------------------------------------------------*/
/*                                                                            */
/* Description:         lplib's helper functions' headers                     */
/* Author:              Loic MARECHAL                                         */
/* Creation date:       may 16 2024                                           */
/* Last modification:   mar 30 2026                                           */
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


/*----------------------------------------------------------------------------*/
/* Prototypes of public procedures                                            */
/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#endif

itg ParallelBuildEdges(int, itg, int, itg *, itg **);
int ParallelNeighbours(int, itg, int, itg *, itg *, char *);
   
#ifdef __cplusplus
} // end extern "C"
#endif

#endif
