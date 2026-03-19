

/*----------------------------------------------------------------------------*/
/*                                                                            */
/*          COMPARE QSORT, PSORT AGAINST LPLIB'S RADIX SORT                   */
/*                                                                            */
/*----------------------------------------------------------------------------*/
/*                                                                            */
/*    Description:         Sort an array of integers with three algorithms    */
/*    Author:              Loic MARECHAL                                      */
/*    Creation date:       mar 10 2026                                        */
/*    Last modification:   mar 19 2026                                        */
/*                                                                            */
/*----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------*/
/* Includes                                                                   */
/*----------------------------------------------------------------------------*/

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include "lplib4.h"


/*----------------------------------------------------------------------------*/
/* Defines                                                                    */
/*----------------------------------------------------------------------------*/

#define nel 112000000LL


/*----------------------------------------------------------------------------*/
/* QSORT comparison function                                                  */
/*----------------------------------------------------------------------------*/

int compar(const void *a, const void *b)
{
   uint32_t *pa = (uint32_t *)a, *pb = (uint32_t *)b;
   return(*pa > *pb ? 1 : -1);
}


/*----------------------------------------------------------------------------*/
/* Compare all three sorts run time                                           */
/*----------------------------------------------------------------------------*/

int main()
{
   int64_t ParIdx;
   uint32_t i, (*a)[2];
   double tim;

   if(!(a = malloc(nel * 2 * sizeof(uint32_t))))
   {
      printf("Failed to allocate %lld bytes\n", nel);
      exit(1);
   }

   if(!(ParIdx = InitParallel(0)))
   {
      puts("Failed to initialize the LPlib");
      exit(1);
   }


   //-------------
   // SERIAL QSORT
   //-------------

   puts("\n\nSerial QSORT\n");

   // Initialize the table with a random number and an index
   puts("Before sorting");
   for(i=0;i<nel;i++)
   {
      a[i][0] = random();
      a[i][1] = i;
   }

   // Print the first and last lines before the sorting
   for(i=0;i<5;i++)
      printf("line %13d : %13u | %13u\n", i, a[i][0], a[i][1]);

   for(i=nel-5;i<nel;i++)
      printf("line %13d : %13u | %13u\n", i, a[i][0], a[i][1]);

   // Perform a sequential qsort
   tim = GetWallClock();
   qsort(a, nel, 2 * sizeof(uint32_t), compar);
   printf("\ntime for serial qsort = %g s\n\n", GetWallClock() - tim);

   // Print the first and last lines after the sorting
   puts("After sorting");
   for(i=0;i<5;i++)
      printf("line %13d : %13u | %13u\n", i, a[i][0], a[i][1]);

   for(i=nel-5;i<nel;i++)
      printf("line %13d : %13u | %13u\n", i, a[i][0], a[i][1]);


#ifdef __MACH__

   //---------------
   // PARALLEL QSORT
   //---------------

   puts("\n\nParallel QSORT");
   for(i=0;i<5;i++)
      printf("line %13d : %13u | %13u\n", i, a[i][0], a[i][1]);

   for(i=nel-5;i<nel;i++)
      printf("line %13d : %13u | %13u\n", i, a[i][0], a[i][1]);

   // Perform a parallel qsort
   tim = GetWallClock();
   psort(a, nel, 2 * sizeof(uint32_t), compar);
   printf("\ntime for parallel qsort = %g s\n\n", GetWallClock() - tim);

   // Print the first and last lines after the sorting
   for(i=0;i<5;i++)
      printf("line %13d : %13u | %13u\n", i, a[i][0], a[i][1]);

   for(i=nel-5;i<nel;i++)
      printf("line %13d : %13u | %13u\n", i, a[i][0], a[i][1]);

#endif


   //------------------
   // SERIAL RADIX SORT
   //------------------

   puts("\n\nSerial RADIX sort");
   for(i=0;i<5;i++)
      printf("line %13d : %13u | %13u\n", i, a[i][0], a[i][1]);

   for(i=nel-5;i<nel;i++)
      printf("line %13d : %13u | %13u\n", i, a[i][0], a[i][1]);

   // Perform a parallel qsort
   tim = GetWallClock();
   RadixSort(ParIdx, nel, a, LplInt, NULL, 0);
   printf("\ntime for radix sort = %g s\n\n", GetWallClock() - tim);

   // Print the first and last lines after the sorting
   for(i=0;i<5;i++)
      printf("line %13d : %13u | %13u\n", i, a[i][0], a[i][1]);

   for(i=nel-5;i<nel;i++)
      printf("line %13d : %13u | %13u\n", i, a[i][0], a[i][1]);

   return(0);
}
