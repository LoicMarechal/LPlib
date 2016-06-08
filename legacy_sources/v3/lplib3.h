

/*----------------------------------------------------------*/
/*															*/
/*						LPLIB	V3.51						*/
/*															*/
/*----------------------------------------------------------*/
/*															*/
/*	Description:		Handles threads, scheduling			*/
/*						& dependencies						*/
/*	Author:				Loic MARECHAL						*/
/*	Creation date:		feb 25 2008							*/
/*	Last modification:	dec 04 2015							*/
/*															*/
/*----------------------------------------------------------*/


#ifdef INT64
#define LplInt long long
#else
#define LplInt int
#endif


/*----------------------------------------------------------*/
/* User available procedures' prototypes					*/
/*----------------------------------------------------------*/

long long InitParallel(LplInt);
void StopParallel(long long);
LplInt NewType(long long, LplInt);
LplInt ResizeType(long long, LplInt, LplInt);
void FreeType(long long, LplInt);
LplInt BeginDependency(long long, LplInt, LplInt);
LplInt AddDependency(long long, LplInt, LplInt);
void AddDependencyFast(long long, LplInt, LplInt *, LplInt, LplInt *);
LplInt UpdateDependency(long long, LplInt, LplInt, LplInt, LplInt);
void UpdateDependencyFast(long long, LplInt, LplInt, LplInt *, LplInt, LplInt, LplInt *);
LplInt EndDependency(long long, float [2]);
void GetDependencyStats(long long, LplInt, LplInt, float [2]);
float LaunchParallel(long long, LplInt, LplInt, void *, void *);
LplInt LaunchPipeline(long long, void *, void *, LplInt, LplInt *);
void WaitPipeline(long long);
LplInt ParallelMemClear(long long , void *, size_t);
void ParallelQsort(long long, void *, size_t, size_t, int (*)(const void *, const void *));
LplInt HilbertRenumbering(long long, LplInt, double [6], double (*)[3], unsigned long long (*)[2]);
LplInt HilbertRenumbering2D(long long, LplInt, double [4], double (*)[2], unsigned long long (*)[2]);
void GetLplibInformation(long long, LplInt *, LplInt *);
LplInt GetNumberOfCores();
double GetWallClock();


/*----------------------------------------------------------*/
/* Public defines											*/
/*----------------------------------------------------------*/

#define MaxPth 128
