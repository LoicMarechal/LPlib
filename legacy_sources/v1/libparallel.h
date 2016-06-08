

/*----------------------------------------------------------*/
/*															*/
/*							LIBPARAL						*/
/*															*/
/*----------------------------------------------------------*/
/*															*/
/*	Description:		Handles scheduling & threads launch	*/
/*	Author:				Loic MARECHAL						*/
/*	Creation date:		feb 25 2008							*/
/*	Last modification:	feb 27 2008							*/
/*															*/
/*----------------------------------------------------------*/


/*----------------------------------------------------------*/
/* User available procedures' prototypes					*/
/*----------------------------------------------------------*/

int InitParallel(int);
void StopParallel(void);
int InitType(int, int);
int BeginDependency(int, int);
void AddDependency(int, int);
void EndDependency(float [2]);
float LaunchParallel(int, int, void *, void *);
