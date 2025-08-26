

/*----------------------------------------------------------------------------*/
/*                                                                            */
/*             MEM CLEAR AND COPY PARALLEL TEST USING LPLib4                  */
/*                                                                            */
/*----------------------------------------------------------------------------*/
/*                                                                            */
/*   Description:       clear v1(i) and copy v2(i) to v3(i)                   */
/*   Author:            Loic MARECHAL                                         */
/*   Creation date:     aug 26 2025                                           */
/*   Last modification: aug 26 2025                                           */
/*                                                                            */
/*----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------*/
/* Includes                                                                   */
/*----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include "lplib4.h"


/*----------------------------------------------------------------------------*/
/* Defines                                                                    */
/*----------------------------------------------------------------------------*/

#define size 1000000000ULL


/*----------------------------------------------------------------------------*/
/* The main procedure reads the number of threads to launch                   */
/*----------------------------------------------------------------------------*/

int main(int ArgCnt, char **ArgVec)
{
   int    i, NmbCpu=0, TypIdx;
   int64_t LibParIdx;
   float   acc;
   double  *vec1, *vec2, *vec3, tim=0.;

   // Read the command line arguments
   if(ArgCnt > 1)
      NmbCpu = atoi(*++ArgVec);

   // Alloocate and initialize all three vectors
   vec1 = malloc((size_t)(size+1) * sizeof(double));
   vec2 = malloc((size_t)(size+1) * sizeof(double));
   vec3 = malloc((size_t)(size+1) * sizeof(double));

   if(!vec1 || !vec2 || !vec3)
   {
      puts("malloc failed");
      exit(1);
   }

   for(i=0;i<size;i++)
   {
      vec1[i] = i;
      vec2[i] = i*2;
   }

   /* Initialize the LPLib: if the return value is 0, something went wrong,
      otherwise, an instanciation index is returned. 
      This index should be provided to any further library function */

   if(!(LibParIdx = InitParallel(NmbCpu)))
   {
      puts("Error initializing the LPLib4.");
      exit(1);
   }

   // System's sequential memcpy
   tim = GetWallClock();
   memcpy(vec3, vec2, size * sizeof(double));
   tim = GetWallClock() - tim;
   printf("sequential memcpy wall clock = %g s\n", tim);

   // LPlib's parallel call to memcpy
   tim = GetWallClock();
   ParallelMemCopy(LibParIdx, vec3, vec2, size * sizeof(double));
   tim = GetWallClock() - tim;
   printf("parallel copy wall clock     = %g s\n", tim);


   // System's sequential memset
   tim = GetWallClock();
   memset(vec1, 0, size * sizeof(double));
   tim = GetWallClock() - tim;
   printf("sequential memset wall clock = %g s\n", tim);

   // LPlib's parallel call to memset
   tim = GetWallClock();
   ParallelMemClear(LibParIdx, vec1, size * sizeof(double));
   tim = GetWallClock() - tim;
   printf("parallel clear wall clock    = %g s\n", tim);


   /* Stop LPLib and free everything */
   StopParallel(LibParIdx);

   free(vec1);
   free(vec2);
   free(vec3);

   return(0);
}
