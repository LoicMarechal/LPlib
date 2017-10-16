

/*----------------------------------------------------------------------------*/
/*                                                                            */
/*                               LPlib V3.53                                  */
/*                                                                            */
/*----------------------------------------------------------------------------*/
/*                                                                            */
/*   Description:       Handles threads, scheduling                           */
/*                      & dependencies                                        */
/*   Author:            Loic MARECHAL                                         */
/*   Creation date:     feb 25 2008                                           */
/*   Last modification: sep 07 2017                                           */
/*                                                                            */
/*----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------*/
/* Types definitions                                                          */
/*----------------------------------------------------------------------------*/

#ifdef INT64
#define LplInt int64_t
#else
#define LplInt int
#endif

#include <stdint.h>


/*----------------------------------------------------------------------------*/
/* User available procedures' prototypes                                      */
/*----------------------------------------------------------------------------*/

int      AddDependency           (int64_t, LplInt, LplInt);
void     AddDependencyFast       (int64_t, int, LplInt *, int, LplInt *);
int      BeginDependency         (int64_t, int, int);
int      EndDependency           (int64_t, float [2]);
void     FreeType                (int64_t, int);
void     GetDependencyStats      (int64_t, int, int, float [2]);
void     GetLplibInformation     (int64_t, int *, int *);
int      GetNumberOfCores        ();
double   GetWallClock            ();
int      HilbertRenumbering      (int64_t, LplInt, double [6], double (*)[3], uint64_t (*)[2]);
int      HilbertRenumbering2D    (int64_t, LplInt, double [4], double (*)[2], uint64_t (*)[2]);
int64_t  InitParallel            (int);
float    LaunchParallel          (int64_t, int, int, void *, void *);
float    LaunchParallelMultiArg  (int64_t, int, int, void *, int, ...);
int      LaunchPipeline          (int64_t, void *, void *, int, int *);
int      LaunchPipelineMultiArg  (int64_t, int, int *, void *prc, int, ...);
int      NewType                 (int64_t, LplInt);
int      ParallelMemClear        (int64_t, void *, size_t);
void     ParallelQsort           (int64_t, void *, size_t, size_t, int (*)(const void *, const void *));
int      ResizeType              (int64_t, int, LplInt);
void     StopParallel            (int64_t);
int      UpdateDependency        (int64_t, int, int, LplInt, LplInt);
void     UpdateDependencyFast    (int64_t, int, int, LplInt *, int, int, LplInt *);
void     WaitPipeline            (int64_t);


/*----------------------------------------------------------------------------*/
/* Public defines                                                             */
/*----------------------------------------------------------------------------*/

#define MaxPth 128
