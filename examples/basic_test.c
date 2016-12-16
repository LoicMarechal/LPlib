

/*----------------------------------------------------------*/
/*															*/
/*				BASIC PARALLEL TEST USING LPLib3			*/
/*															*/
/*----------------------------------------------------------*/
/*															*/
/*	Description:		compute v3 = f(v1) + g(v2)			*/
/*	Author:				Loic MARECHAL						*/
/*	Creation date:		apr 22 2008							*/
/*	Last modification:	jan 20 2015							*/
/*															*/
/*----------------------------------------------------------*/


/*----------------------------------------------------------*/
/* Includes													*/
/*----------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <sys/time.h>
#include "lplib3.h"


/*----------------------------------------------------------*/
/* Defines													*/
/*----------------------------------------------------------*/

#define size 100000000


/*----------------------------------------------------------*/
/* Structures' prototypes									*/
/*----------------------------------------------------------*/

typedef struct
{
	double *vec1, *vec2, *vec3;
}AddArgSct;


/*----------------------------------------------------------*/
/* Procedure called in parallel computing v3 = f(v1) + g(v2)*/
/* Arguments are imposed by the library :					*/
/* BegIdx : begining loop index								*/
/* EndIdx : ending loop index								*/
/* PthIdx : index of the thread executing this instance		*/
/* arg    : pointer to a user's structure containing every	*/
/*          arguments needed by the loop					*/
/*----------------------------------------------------------*/

void AddVec(int BegIdx, int EndIdx, int PthIdx, AddArgSct *arg)
{
	int i;
	double *vec1, *vec2, *vec3;

	/* Copy arguments in local variables (it may help avoiding some false memory conflicts) */
	vec1 = arg->vec1;
	vec2 = arg->vec2;
	vec3 = arg->vec3;

	/* Compute v3 = f(v1) + g(v2) within the given range */
	for(i=BegIdx; i<=EndIdx; i++)
		vec3[i] = exp(log(acos(cos(vec1[i]+1.)))) + exp(log(acos(cos(vec2[i]+2.))));
}


/*----------------------------------------------------------*/
/* The main procedure reads the number of threads to launch	*/
/*----------------------------------------------------------*/

int main(int ArgCnt, char **ArgVec)
{
	int i, NmbCpu, TypIdx;
    long long LibParIdx;
	float acc;
	double *vec1, *vec2, *vec3, tim=0.;
	AddArgSct arg;
	struct timeval tp;

	/* Read the command line arguments */
	if(ArgCnt > 1)
		NmbCpu = atoi(*++ArgVec);

	/* Alloocate and initialize all three vectors */
	vec1 = malloc(size * sizeof(double));
	vec2 = malloc(size * sizeof(double));
	vec3 = malloc(size * sizeof(double));

	for(i=0;i<size;i++)
	{
		vec1[i] = i;
		vec2[i] = i*2;
		vec3[i] = 0;
	}

	/* Build the arguments structure: it stores only pointers to the three vectors */
	arg.vec1 = vec1;
	arg.vec2 = vec2;
	arg.vec3 = vec3;

	/* Initialize the LPLib3: if the return value is 0, something went wrong, otherwise,
		an instanciation index is returned. This index should be provided to any further
		library function */
	if(!(LibParIdx = InitParallel(NmbCpu)))
	{
		puts("Error initializing the LPLib3.");
		exit(1);
	}

	/* Setup a new data type: provide the instanciation number and the size of the vector.
		The procedure return 0 for a failure and a data type index otherwise. This new index
		must be provided to launch any parallel loop on this vector */
	if(!(TypIdx = NewType(LibParIdx, size)))
	{
		puts("Error while creating new data type.");
		exit(1);
	}

	gettimeofday(&tp, NULL);
	tim = tp.tv_sec + tp.tv_usec / 1000000.;

	/* Launch the parallel loop computing v3 = f(v1) + g(v2) with the following arguments :
		-library instanciation index
		-base data type index
		-dependency index (0 since there is no dependencies in this example)
		-pointer to the procedure to be launched in parallel
		-pointer to arguments structure needed by the loop
		The return value is 0 is something went wrong or an acceleration factor otherwise */
	if(!(acc = LaunchParallel(LibParIdx, TypIdx, 0, (void *)AddVec, (void *)&arg)))
	{
		puts("Error while launching the parallel loop AddVec.");
		exit(1);
	}

	gettimeofday(&tp, NULL);
	tim = tp.tv_sec + tp.tv_usec / 1000000. - tim;

	printf("Theoretical speedup for loop AddVec = %g, wall clock = %g s\n", acc, tim);

	/* Stop LPLib3 and free everything */
	StopParallel(LibParIdx);

	free(vec1);
	free(vec2);
	free(vec3);

	return(0);
}
