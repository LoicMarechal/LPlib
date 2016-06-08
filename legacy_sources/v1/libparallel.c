

/*----------------------------------------------------------*/
/*															*/
/*					LIBPARALLEL	V1.0						*/
/*															*/
/*----------------------------------------------------------*/
/*															*/
/*	Description:		Handles threads, scheduling			*/
/*						& dependencies						*/
/*	Author:				Loic MARECHAL						*/
/*	Creation date:		feb 25 2008							*/
/*	Last modification:	feb 27 2008							*/
/*															*/
/*----------------------------------------------------------*/


#ifndef SERIAL

/*----------------------------------------------------------*/
/* Includes													*/
/*----------------------------------------------------------*/

#include <unistd.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "libparallel.h"


/*----------------------------------------------------------*/
/* Defines													*/
/*----------------------------------------------------------*/

enum ParSta {ParIdl, ParEnd, ParLok, ParRun};

#define MaxTyp 100
#define WrkSiz 1024
#define DepSiz 1024
#define MaxPth 64
#define IdlTim 1000


/*----------------------------------------------------------*/
/* Structures' prototypes									*/
/*----------------------------------------------------------*/

typedef struct
{
	int BegIdx, EndIdx, flg, NmbBit, *dep;
}WrkSct;

typedef struct
{
	int NmbLin, NmbWrk, NmbDep, *DepTab;
	WrkSct *WrkTab;
}TypSct;

typedef struct
{
	int idx, sta;
	WrkSct *wrk;
	pthread_t pth;
}PthSct;

typedef struct
{
	int NmbCpu, WrkCpt;
	float sta[2];
	void (*prc)(int, int, int, void *), *arg;
	pthread_mutex_t *mtx;
	PthSct PthTab[ MaxPth ];
	TypSct TypTab[ MaxTyp+1 ], *CurTyp, *typ1, *typ2;
}ParSct;


/*----------------------------------------------------------*/
/* Private procedures' prototypes							*/
/*----------------------------------------------------------*/

static void *PthSch(void *);
static int WrkLok(WrkSct *, WrkSct *);
static void SetBit(int *, int);
static int GetBit(int *, int);
static int AndBit(int *, int *);
int CmpWrk(const void *, const void *);


/*----------------------------------------------------------*/
/* Global variables											*/
/*----------------------------------------------------------*/

ParSct *par;
pthread_attr_t PthAtr;
pthread_mutex_t mtx=PTHREAD_MUTEX_INITIALIZER;


/*----------------------------------------------------------*/
/* Init structures, scheduler and launch threads			*/
/*----------------------------------------------------------*/

int InitParallel(int NmbCpu)
{
	int i;

	/* Alloc and init Pthreads */

	pthread_attr_init(&PthAtr);
	pthread_attr_setdetachstate(&PthAtr, PTHREAD_CREATE_JOINABLE);

	/* Allocate and build main parallel structure */

	if(!(par = calloc(1, sizeof(ParSct))))
		return(0);

	par->NmbCpu = NmbCpu;
	par->mtx = &mtx;

	/* Launch Pthreads */

	for(i=0;i<NmbCpu;i++)
	{
		par->PthTab[i].idx = i;
		par->PthTab[i].sta = ParIdl;
		par->PthTab[i].wrk = NULL;
		memset(&par->PthTab[i].pth, 0, sizeof(pthread_t));
		pthread_create(&par->PthTab[i].pth, &PthAtr, PthSch, (void *)&par->PthTab[i]);
	}

	return(1);
}


/*----------------------------------------------------------*/
/* Stop all threads and free memories						*/
/*----------------------------------------------------------*/

void StopParallel(void)
{
	int i, sta;
	pthread_mutex_lock(par->mtx);

	/* Stop threads */

	for(i=0;i<par->NmbCpu;i++)
		par->PthTab[i].sta = ParEnd;

	pthread_mutex_unlock(par->mtx);

	for(i=0;i<par->NmbCpu;i++)
		pthread_join(par->PthTab[i].pth, (void **)&sta);

	/* Free memories */

	for(i=0;i<MaxTyp;i++)
	{
		if(par->TypTab[i].WrkTab)
			free(par->TypTab[i].WrkTab);

		if(par->TypTab[i].DepTab)
			free(par->TypTab[i].DepTab);
	}

	free(par);
	pthread_attr_destroy(&PthAtr);
}


/*----------------------------------------------------------*/
/* Set work-packages for this kind of element				*/
/*----------------------------------------------------------*/

int InitType(int TypIdx, int NmbLin)
{
	int i, res, idx=0;
	TypSct *typ;

	/* Check bounds */

	if( (TypIdx < 1) || (TypIdx > MaxTyp) || (NmbLin < 1) )
		return(0);

	typ = &par->TypTab[ TypIdx ];

	/* Free previously allocated memory in case of type resizing */

	if(typ->WrkTab)
		free(typ->WrkTab);

	if(typ->DepTab)
		free(typ->DepTab);

	typ->NmbLin = NmbLin;
	typ->NmbWrk = NmbLin / WrkSiz;

	if(res = NmbLin - typ->NmbWrk * WrkSiz)
		typ->NmbWrk++;

	if(!(typ->WrkTab = calloc(typ->NmbWrk , sizeof(WrkSct))))
		return(0);

	/* Set work-packages */

	for(i=0;i<typ->NmbWrk;i++)
	{
		typ->WrkTab[i].BegIdx = idx + 1;
		typ->WrkTab[i].EndIdx = idx + WrkSiz;
		idx += WrkSiz;
	}

	if(res)
		typ->WrkTab[ typ->NmbWrk - 1 ].EndIdx = NmbLin;

	return(1);
}


/*----------------------------------------------------------*/
/* Allocate a dependency matrix linking both types			*/
/*----------------------------------------------------------*/

int BeginDependency(int TypIdx1, int TypIdx2)
{
	int i, res;
	TypSct *typ1 = &par->TypTab[ TypIdx1 ], *typ2 = &par->TypTab[ TypIdx2 ];

	/* Check bounds */

	if( (TypIdx1 < 1) || (TypIdx1 > MaxTyp) || (TypIdx2 < 1) || (TypIdx2 > MaxTyp) || (typ1 == typ2) \
	|| !typ1->NmbLin || !typ2->NmbLin)
	{
		return(0);
	}

	par->CurTyp = typ1;
	typ1->NmbDep = typ2->NmbLin / (DepSiz * 32);

	if(res = typ2->NmbLin - typ1->NmbDep * DepSiz * 32)
		typ1->NmbDep++;

	/* Allocate a global dependency table */

	if(!(typ1->DepTab = calloc(typ1->NmbWrk * typ1->NmbDep , sizeof(int))))
		return(0);

	/* Then spread sub-tables among WP */

	for(i=0;i<typ1->NmbWrk;i++)
	{
		typ1->WrkTab[i].NmbBit = 0;
		typ1->WrkTab[i].dep = &typ1->DepTab[ i * typ1->NmbDep ];
	}

	return(1);
}


/*----------------------------------------------------------*/
/* Element idx1 of type1 depends on element idx2 of type2	*/
/*----------------------------------------------------------*/

void AddDependency(int idx1, int idx2)
{
	WrkSct *wrk;

	/* Set and count dependency bit */

	wrk = &par->CurTyp->WrkTab[ (idx1-1) / WrkSiz ];

	if(!GetBit(wrk->dep, (idx2-1) / DepSiz ))
	{
		SetBit(wrk->dep, (idx2-1) / DepSiz );
		wrk->NmbBit++;
	}
}


/*----------------------------------------------------------*/
/* Sort wp depending on their number of dependencies		*/
/*----------------------------------------------------------*/

void EndDependency(float DepSta[2])
{
	int i;

	/* Compute average and max number of collisions */

	DepSta[0] = DepSta[1] = 0.;

	for(i=0;i<par->CurTyp->NmbWrk;i++)
	{
		DepSta[0] += par->CurTyp->WrkTab[i].NmbBit;

		if(par->CurTyp->WrkTab[i].NmbBit > DepSta[1])
			DepSta[1] = par->CurTyp->WrkTab[i].NmbBit;
	}

	DepSta[0] = 100 * DepSta[0] / (par->CurTyp->NmbWrk * par->CurTyp->NmbDep * 32);
	DepSta[1] = 100 * DepSta[1] / (par->CurTyp->NmbDep * 32);

	/* Sort WP from highest collision number to the lowest */

	qsort(par->CurTyp->WrkTab, par->CurTyp->NmbWrk, sizeof(WrkSct), CmpWrk);

	par->CurTyp = NULL;
}


/*----------------------------------------------------------*/
/* Launch the loop prc on typ1 element depending on typ2	*/
/*----------------------------------------------------------*/

float LaunchParallel(int typ1, int typ2, void *prc, void *arg)
{
	int i, sta;

	/* Check bounds */

	if( (typ1 < 1) || (typ1 > MaxTyp) || (typ2 < 0) || (typ2 > MaxTyp) || (typ1 == typ2) )
		return(0.);

	/* Protect launch command within a mutex */

	pthread_mutex_lock(par->mtx);

	par->prc = (void (*)(int, int, int, void *))prc;
	par->arg = arg;
	par->typ1 = &par->TypTab[ typ1 ];
	par->WrkCpt = 0;
	par->sta[0] = par->sta[1] = 0.;

	if(typ2)
		par->typ2 = &par->TypTab[ typ2 ];
	else
		par->typ2 = NULL;

	for(i=0;i<par->NmbCpu;i++)
		par->PthTab[i].sta = ParLok;

	for(i=0;i<par->typ1->NmbWrk;i++)
		par->typ1->WrkTab[i].flg = 0;

	pthread_mutex_unlock(par->mtx);

	/* Now wait for every threads to get idle */

	do
	{
		pthread_mutex_lock(par->mtx);
		sta = ParIdl;

		for(i=0;i<par->NmbCpu;i++)
			if(par->PthTab[i].sta != ParIdl)
				sta = ParRun;

		pthread_mutex_unlock(par->mtx);
		usleep(IdlTim);
	}while(sta == ParRun);

	/* Return the average speedup */

	return(par->sta[1] / par->sta[0]);
}


/*----------------------------------------------------------*/
/* the scheduler set threads states : stop, idle, wait, run */
/*----------------------------------------------------------*/

static void *PthSch(void *arg)
{
	int i, j, sta;
	WrkSct *wrk;
	PthSct *pth = (PthSct *)arg;

	do
	{
		pthread_mutex_lock(par->mtx);
		sta = pth->sta;

		if(sta == ParLok)
		{
			wrk = NULL;

			for(i=0;i<par->typ1->NmbWrk;i++)
			{
				wrk = &par->typ1->WrkTab[i];

				if(wrk->flg)
				{
					wrk = NULL;
					continue;
				}

				if(par->typ2)
					for(j=0;j<par->NmbCpu;j++)
					{
						if(j == pth->idx)
							continue;
            
						if(par->PthTab[j].wrk && AndBit(par->PthTab[j].wrk->dep, wrk->dep))
						{
							wrk = NULL;
							break;
						}
					}

				if(wrk)
					break;
			}

			if(wrk)
			{
				sta = pth->sta = ParRun;
				pth->wrk = wrk;
				wrk->flg = 1;
				par->WrkCpt++;
				par->sta[0]++;

				for(i=0;i<par->NmbCpu;i++)
					if(par->PthTab[i].wrk)
						par->sta[1]++;
			}
			else if(par->WrkCpt == par->typ1->NmbWrk)
				sta = pth->sta = ParIdl;
		}

		pthread_mutex_unlock(par->mtx);

		if( (sta == ParIdl) || (sta == ParLok) )
			usleep(IdlTim);
		else if(sta == ParRun)
		{
			par->prc(wrk->BegIdx, wrk->EndIdx, pth->idx, par->arg);
			pth->sta = ParLok;
			pth->wrk = NULL;
		}
	}while(sta != ParEnd);

	return(NULL);
}


/*----------------------------------------------------------*/
/* Set/get a bit in a multibyte number						*/
/*----------------------------------------------------------*/

static void SetBit(int *tab, int idx)
{
	int adr = idx >> 5;
	tab[ adr ] |= 1 << (idx - (adr << 5));
}

static int GetBit(int *tab, int idx)
{
	int adr = idx >> 5;
	return( (tab[ adr ] & (1 << (idx - (adr << 5)))) );
}


/*----------------------------------------------------------*/
/* Check wether two WP share common resources -> locked		*/
/*----------------------------------------------------------*/

static int AndBit(int *t1, int *t2)
{
	int i;

	for(i=0;i<par->typ1->NmbDep;i++)
		if(t1[i] & t2[i])
			return(1);

	return(0);
}


/*----------------------------------------------------------*/
/* Compare two workpackages number of bits					*/
/*----------------------------------------------------------*/

int CmpWrk(const void *ptr1, const void *ptr2)
{
	WrkSct *w1, *w2;

	w1 = (WrkSct *)ptr1;
	w2 = (WrkSct *)ptr2;

	if(w1->NmbBit > w2->NmbBit)
		return(-1);
	else if(w1->NmbBit < w2->NmbBit)
		return(1);
	else
		return(0);
}

#else

#include "libparallel.h"
#define MaxTyp 100
int TypTab[ MaxTyp ];

int InitParallel(int NmbCpu)
{
	int i;

	for(i=0;i<MaxTyp;i++)
		TypTab[i] = 0;

	return(1);
}


void StopParallel(void)
{
}


int InitType(int TypIdx, int NmbLin)
{
	if( (TypIdx < 1) || (TypIdx > MaxTyp) || (NmbLin < 1) )
		return(0);

	TypTab[ TypIdx ] = NmbLin;

	return(1);
}


int BeginDependency(int TypIdx1, int TypIdx2)
{
	return(1);
}


void AddDependency(int idx1, int idx2)
{
}


void EndDependency(float DepSta[2])
{
}


float LaunchParallel(int typ1, int typ2, void *prc, void *arg)
{
	int i, sta;
	void (*UsrPrc)(int, int, int, void *) = (void (*)(int, int, int, void *))prc;

	/* Check bounds */

	if( (typ1 < 1) || (typ1 > MaxTyp) || (typ2 < 0) || (typ2 > MaxTyp) || (typ1 == typ2) )
		return(0.);

	UsrPrc(1, TypTab[ typ1 ], 0, arg);

	return(1.);
}

#endif
