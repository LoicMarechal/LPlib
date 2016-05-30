

/*----------------------------------------------------------*/
/*															*/
/*		PARALLEL EDGE LIST BUILDING USING LP3				*/
/*															*/
/*----------------------------------------------------------*/
/*															*/
/*	Description:		build the list of unique edges		*/
/*                      from a volumic tetrahedral mesh		*/
/*	Author:				Loic MARECHAL						*/
/*	Creation date:		feb 13 2015							*/
/*	Last modification:	feb 18 2015							*/
/*															*/
/*----------------------------------------------------------*/


/*----------------------------------------------------------*/
/* Includes													*/
/*----------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>
#include "libmesh6.h"
#include "lplib3.h"


/*----------------------------------------------------------*/
/* Defines													*/
/*----------------------------------------------------------*/

#define min(a,b) ((a) < (b) ? (a) : (b))
#define max(a,b) ((a) > (b) ? (a) : (b))
#define MaxEdg 1000

#ifdef i8
#define lng long long
#define LngTyp GmfLong
#else
#define lng int
#define LngTyp GmfInt
#endif


/*----------------------------------------------------------*/
/* Structures												*/
/*----------------------------------------------------------*/

typedef struct
{
	double crd[3];
	int ref;
}VerSct;

typedef struct
{
	lng idx[2];
}EdgSct;

typedef struct
{
	lng idx[3];
	int ref;
}TriSct;

typedef struct
{
	lng idx[4];
}TetSct;

typedef struct
{
	lng MinIdx, MaxIdx, NexBuc;
}HshSct;

typedef struct
{
	lng NmbVer, NmbEdg, NmbTri, NmbTet;
	int MshVer, LibIdx;
	VerSct *ver;
	EdgSct *edg;
	TriSct *tri;
	TetSct *tet;
}MshSct;

typedef struct
{
	lng beg, end, HshSiz, ColPos, NmbEdg, EdgAdr;
	int NmbCpu;
	HshSct *HshTab;
	MshSct *msh;
}ParSct;


/*----------------------------------------------------------*/
/* Prototypes of local procedures							*/
/*----------------------------------------------------------*/

void SetEdgSer(MshSct *);
void SetEdgPar(MshSct *, int);
void ParEdg1(lng, lng, int, ParSct *);
void ParEdg2(lng, lng, int, ParSct *);
void ScaMsh(char *, MshSct *);
void RecMsh(char *, MshSct *);
void GetTim(double *);


/*----------------------------------------------------------*/
/* Global tables											*/
/*----------------------------------------------------------*/

const int tvpe[6][2] = { {0,1}, {1,2}, {2,0}, {3,0}, {3,1}, {3,2} };


/*----------------------------------------------*/
/* Read, build the edge list and write the mesh	*/
/*----------------------------------------------*/

int main(int ArgCnt, char **ArgVec)
{
	char *PtrArg, *TmpStr, InpNam[1000], OutNam[1000];
	int i, NmbCpu = 0;
	MshSct msh;

	/* Command line parsing */

	memset(&msh, 0, sizeof(MshSct));

	if(ArgCnt == 1)
	{
		puts("\ntetrahedra_edges v1.00 feb 13 2015   Loic MARECHAL / INRIA");
		puts(" Usage       : tetrahedra_edges -in volume_mesh -out edge_mesh");
		puts(" -in name    : name of the input tetrahedral-only mesh");
		puts(" -out name   : name of the output mesh that will contain tets, edges and vertices");
		puts(" -serial     : use the serial optimized version (different from -nproc 1)");
		puts(" -nproc n    : n is the number of threads to be launched (default = all available threads)\n");
		exit(0);
	}

	for(i=2;i<=ArgCnt;i++)
	{
		PtrArg = *++ArgVec;

		if(!strcmp(PtrArg,"-in"))
		{
			TmpStr = *++ArgVec;
			ArgCnt--;
			strcpy(InpNam, TmpStr);

			if(!strstr(InpNam, ".mesh"))
				strcat(InpNam, ".meshb");

			continue;
		}

		if(!strcmp(PtrArg,"-out"))
		{
			TmpStr = *++ArgVec;
			ArgCnt--;
			strcpy(OutNam, TmpStr);

			if(!strstr(OutNam, ".mesh"))
				strcat(OutNam, ".meshb");

			continue;
		}

		if(!strcmp(PtrArg,"-serial"))
		{
			NmbCpu = -1;
			continue;
		}

		if(!strcmp(PtrArg,"-nproc"))
		{
			NmbCpu = atoi(*++ArgVec);
			NmbCpu = max(NmbCpu, 1);
			NmbCpu = min(NmbCpu, 128);
			ArgCnt--;
			continue;
		}
	}

	if(!strlen(InpNam))
	{
		puts("No input mesh provided");
		exit(1);
	}

	if(!strlen(OutNam))
	{
		puts("No output name provided");
		exit(1);
	}

	// Mesh reading
	ScaMsh(InpNam, &msh);

	// Launch the parallel neeighbour procedure

	if(NmbCpu == -1)
		SetEdgSer(&msh);
	else
		SetEdgPar(&msh, NmbCpu);

	// Mesh writing
	//RecMsh(OutNam, &msh);

	//FreMem(msh.ver, (msh.NmbVer+1) * sizeof(VerSct));
	FreMem(msh.edg, (msh.NmbEdg+1) * sizeof(EdgSct));
	//FreMem(msh.tri, (msh.NmbTri+1) * sizeof(TriSct));
	FreMem(msh.tet, (msh.NmbTet+1) * sizeof(TetSct));
}


/*----------------------------------------------------------*/
/* Build the list of unique edges sequentialy				*/
/*----------------------------------------------------------*/

void SetEdgSer(MshSct *msh)
{
	lng i, j, idx0, idx1, key, MinIdx, MaxIdx, siz = msh->NmbTet, col = siz;
	double timer = 0.;
	TetSct *tet;
	HshSct *hsh = AloMem(6 * siz * sizeof(HshSct), GloMem), *buc;

	/* Start timer */
	printf("Build edges list sequentialy : ");
	GetTim(&timer);

	/* Clear the hash table's direct entries, there is no need to clear the collision entries */
	memset(hsh, 0, siz * sizeof(HshSct));

	/* Loop over each tet and each tet's edges */
	for(i=1;i<=msh->NmbTet;i++)
	{
		tet = &msh->tet[i];

		for(j=0;j<6;j++)
		{
			/* Compute the hashing key from the edge's vertex indices */

			idx0 = tet->idx[ tvpe[j][0] ];
			idx1 = tet->idx[ tvpe[j][1] ];

			if(idx0 < idx1)
			{
				MinIdx = idx0;
				MaxIdx = idx1;
			}
			else
			{
				MinIdx = idx1;
				MaxIdx = idx0;
			}

			key = (3 * MinIdx + 5 * MaxIdx) % siz;

			/* If the bucket is empty, store the edge */

			if(!hsh[ key ].MinIdx)
			{
				hsh[ key ].MinIdx = MinIdx;
				hsh[ key ].MaxIdx = MaxIdx;
				msh->NmbEdg++;
				continue;
			}

			/* Otherwise, search through the linked list */

			do
			{
				/* If the same edge is found in the hash table, do nothing */

				if( (hsh[ key ].MinIdx == MinIdx) && (hsh[ key ].MaxIdx == MaxIdx) )
					break;

				/* If not, allocate a new bucket from the overflow table and link it to the main entry */

				if(hsh[ key ].NexBuc)
					key = hsh[ key ].NexBuc;
				else
				{
					hsh[ key ].NexBuc = col;
					key = col++;
					hsh[ key ].MinIdx = MinIdx;
					hsh[ key ].MaxIdx = MaxIdx;
					msh->NmbEdg++;
					break;
				}
			}while(1);
		}
	}

	/* Now allocate and edge table */

	msh->edg = AloMem( (msh->NmbEdg+1) * sizeof(EdgSct), GloMem);
	msh->NmbEdg = 0;

	/* Loop over the hash table's direct entry and store the edges */
	for(i=0;i<msh->NmbTet;i++)
	{
		key = i;

		/* Follow the link if this bucket has some collisions */
		do
		{
			buc = &hsh[ key ];

			if(buc->MinIdx)
			{
				msh->NmbEdg++;
				msh->edg[ msh->NmbEdg ].idx[0] = buc->MinIdx;
				msh->edg[ msh->NmbEdg ].idx[1] = buc->MaxIdx;
			}
		}while((key = buc->NexBuc));
	}

	/* Free the hash table and stop the timer */

	FreMem(hsh, 6 * siz * sizeof(HshSct));
	GetTim(&timer);
	printf("%g s\n", timer);
	printf("%zd unique edges found\n", msh->NmbEdg);
}


/*----------------------------------------------------------*/
/* Build edges in parallel: head procedure					*/
/*----------------------------------------------------------*/

void SetEdgPar(MshSct *msh, int NmbCpu)
{
	lng i, HshSiz, IncSiz, adr = 0;
	int LibIdx, TetTyp;
	float sta[2];
	double timer = 0.;
	HshSct *HshTab;
	ParSct par[ MaxPth ];

	/* Setup LP3 lib and datatypes */

	if(!NmbCpu)
		NmbCpu = GetNumberOfCores();

	printf("Build edges list with %d threads     : ", NmbCpu);

	GetTim(&timer);

	LibIdx = InitParallel(NmbCpu);
	TetTyp = NewType(LibIdx, msh->NmbTet);

	/* Setup parallel parameters */

	IncSiz = (msh->NmbTet / NmbCpu) / NmbCpu;
	HshSiz =  IncSiz * NmbCpu;

	for(i=0;i<NmbCpu;i++)
	{
		par[i].beg = i * IncSiz;
		par[i].end = (i + 1) * IncSiz;
		par[i].HshSiz = HshSiz;
		par[i].ColPos = HshSiz;
		par[i].msh = msh;
		par[i].NmbCpu = NmbCpu;
		par[i].EdgAdr = 0;
	}

	/* Each thread builds a local edge table */
	LaunchParallel(LibIdx, TetTyp, 0, (void *)ParEdg1, (void *)par);

	/* Count the number of unique edges in each thread */
	LaunchParallel(LibIdx, TetTyp, 0, (void *)ParEdg2, (void *)par);

	/* Allocate the global edge table and give a slice of it to each threads */

	msh->NmbEdg = 0;

	for(i=0;i<NmbCpu;i++)
	{
		par[i].EdgAdr = msh->NmbEdg + 1;
		msh->NmbEdg += par[i].NmbEdg;
	}

	msh->edg = AloMem( (msh->NmbEdg+1) * sizeof(EdgSct), GloMem);

	/* Now each threads counts and stores the unique edges */
	LaunchParallel(LibIdx, TetTyp, 0, (void *)ParEdg2, (void *)par);

	/* Free the local hash tables */
	for(i=0;i<NmbCpu;i++)
		FreMem(par[i].HshTab, 6 * HshSiz * sizeof(HshSct));

	StopParallel(LibIdx);

	GetTim(&timer);
	printf("%g s\n", timer);

	printf("%zd unique edges found\n", msh->NmbEdg);
}


/*----------------------------------------------------------*/
/* Build thread subdomain's edges							*/
/*----------------------------------------------------------*/

void ParEdg1(lng BegIdx, lng EndIdx, int PthIdx, ParSct *par)
{
	lng i, j, key, idx0, idx1, MinIdx, MaxIdx;
	lng siz = par[ PthIdx ].HshSiz, col = par[ PthIdx ].ColPos;
	HshSct *hsh = par[ PthIdx ].HshTab = AloMem( 6 * siz * sizeof(HshSct), LocMem );
	MshSct *msh = par[ PthIdx ].msh;
	TetSct *tet;

	/* Clear the hash table's direct entries, there is no need to clear the collision entries */
	memset(hsh, 0, siz * sizeof(HshSct));

	/* Loop over each tet and each tet's edges */
	for(i=BegIdx; i<=EndIdx; i++)
	{
		tet = &msh->tet[i];

		for(j=0;j<6;j++)
		{
			/* Compute the hashing key from the edge's vertices indices */

			idx0 = tet->idx[ tvpe[j][0] ];
			idx1 = tet->idx[ tvpe[j][1] ];

			if(idx0 < idx1)
			{
				MinIdx = idx0;
				MaxIdx = idx1;
			}
			else
			{
				MinIdx = idx1;
				MaxIdx = idx0;
			}

			key = (3 * MinIdx + 5 * MaxIdx) % siz;

			/* If the bucket is empty, store the edge */

			if(!hsh[ key ].MinIdx)
			{
				hsh[ key ].MinIdx = MinIdx;
				hsh[ key ].MaxIdx = MaxIdx;
				continue;
			}

			/* Otherwise, search through the linked list */

			do
			{
				/* If the same edge is found in the hash table, do nothing */

				if( (hsh[ key ].MinIdx == MinIdx) && (hsh[ key ].MaxIdx == MaxIdx) )
					break;

				/* If not, allocate a new bucket from the overflow table and link it to the main entry */

				if(hsh[ key ].NexBuc)
					key = hsh[ key ].NexBuc;
				else
				{
					hsh[ key ].NexBuc = col;
					key = col++;
					hsh[ key ].MinIdx = MinIdx;
					hsh[ key ].MaxIdx = MaxIdx;
					break;
				}
			}while(1);
		}
	}
}


/*------------------------------------------------------------*/
/* Setup the missing links between tets that cross subdomains */
/*------------------------------------------------------------*/

void ParEdg2(lng BegIdx, lng EndIdx, int PthIdx, ParSct *par)
{
	lng i, key, PthNmbEdg = 0, edg[ MaxEdg ][2];
	lng siz = par[ PthIdx ].HshSiz, NmbCpu = par[ PthIdx ].NmbCpu;
	int NmbEdg, flg, j, k;
	HshSct *hsh = par[ PthIdx ].HshTab, *buc;
	MshSct *msh = par[ PthIdx ].msh;
	TetSct *tet;

	/* Loop over the hash table direct entries following an interleaved stencil */

	for(i=par[ PthIdx ].beg; i<par[ PthIdx ].end; i++)
	{
		NmbEdg = 0;

		/* Loop over every entries with the same key among all threads' local hash tables */

		for(j=0;j<NmbCpu;j++)
		{
			key = i;

			/* In case of collision, follow the links */

			do
			{
				buc = &par[j].HshTab[ key ];

				if(buc->MinIdx)
				{
					/* Since edges from different local hash tables may be the same, 
						they compared again to avoid duplicates */

					flg = 0;

					for(k=0;k<NmbEdg;k++)
						if( (buc->MinIdx == edg[k][0]) && (buc->MaxIdx == edg[k][1]) )
						{
							flg= 1;
							break;
						}

					/* If this edge does not belong to the list, add it to the end */

					if(!flg)
					{
						edg[ NmbEdg ][0] = buc->MinIdx;
						edg[ NmbEdg ][1] = buc->MaxIdx;
						NmbEdg++;

						if(NmbEdg >= MaxEdg)
						{
							puts("Too many local edges, increase MaxEdg value.");
							exit(1);
						}
					}
				}
			}while((key = buc->NexBuc));
		}

		/* On the second run, add the list of local edges to the global ones */

		if(par[ PthIdx ].EdgAdr)
			for(j=0;j<NmbEdg;j++)
			{
				msh->edg[ par[ PthIdx ].EdgAdr + PthNmbEdg + j ].idx[0] = edg[j][0];
				msh->edg[ par[ PthIdx ].EdgAdr + PthNmbEdg + j ].idx[1] = edg[j][1];
			}

		PthNmbEdg += NmbEdg;
	}

	par[ PthIdx ].NmbEdg = PthNmbEdg;
}


/*----------------------------------------------------------*/
/* Wall clock timer											*/
/*----------------------------------------------------------*/

void GetTim(double *timer)
{
	struct timeval tp;
	gettimeofday(&tp, NULL);
	*timer = tp.tv_sec + tp.tv_usec / 1000000. - *timer;
}


/*----------------------------------------------------------*/
/* Read mesh 												*/
/*----------------------------------------------------------*/

void ScaMsh(char *InpNam, MshSct *msh)
{
	int InpMsh, dim, ref;
	float flt[3];
	double timer = 0.;

	printf("\nReading mesh          : ");
	GetTim(&timer);

	/* Check mesh format */

	if(!(InpMsh = GmfOpenMesh(InpNam, GmfRead, &msh->MshVer, &dim)))
	{
		printf("Cannot open mesh %s\n", InpNam);
		exit(1);
	}

	if(dim != 3)
	{
		puts("Can only handle 3D meshes");
		exit(1);
	}

	/* Get stats and allocate tables */

/*	if(!(msh->NmbVer = GmfStatKwd(InpMsh, GmfVertices)))
	{
		puts("No vertices found");
		exit(1);
	}

	if(!(msh->NmbTri = GmfStatKwd(InpMsh, GmfTriangles)))
	{
		puts("No triangles found");
		exit(1);
	}
*/
	if(!(msh->NmbTet = GmfStatKwd(InpMsh, GmfTetrahedra)))
	{
		puts("No tetrahedra found");
		exit(1);
	}

/*	msh->ver = AloMem( (msh->NmbVer+1) * sizeof(VerSct), GloMem);
	msh->tri = AloMem( (msh->NmbTri+1) * sizeof(TriSct), GloMem);*/
	msh->tet = AloMem( (msh->NmbTet+1) * sizeof(TetSct), GloMem);

	/* Read vertices */

/*	GmfGotoKwd(InpMsh, GmfVertices);
	GmfGetBlock(InpMsh, GmfVertices, \
				GmfDouble, &msh->ver[1].crd[0], &msh->ver[2].crd[0], \
				GmfDouble, &msh->ver[1].crd[1], &msh->ver[2].crd[1], \
				GmfDouble, &msh->ver[1].crd[2], &msh->ver[2].crd[2], \
				GmfInt, &msh->ver[1].ref, &msh->ver[2].ref);
*/
	/* Read triangles */

/*	GmfGotoKwd(InpMsh, GmfTriangles);
	GmfGetBlock(InpMsh, GmfTriangles, \
				LngTyp, &msh->tri[1].idx[0], &msh->tri[2].idx[0], \
				LngTyp, &msh->tri[1].idx[1], &msh->tri[2].idx[1], \
				LngTyp, &msh->tri[1].idx[2], &msh->tri[2].idx[2], \
				GmfInt, &msh->tri[1].ref, &msh->tri[2].ref);
*/
	/* Read tetrahedra */

	GmfGotoKwd(InpMsh, GmfTetrahedra);
	GmfGetBlock(InpMsh, GmfTetrahedra, \
				LngTyp, &msh->tet[1].idx[0], &msh->tet[2].idx[0], \
				LngTyp, &msh->tet[1].idx[1], &msh->tet[2].idx[1], \
				LngTyp, &msh->tet[1].idx[2], &msh->tet[2].idx[2], \
				LngTyp, &msh->tet[1].idx[3], &msh->tet[2].idx[3], \
				GmfInt, &ref, &ref);

	GmfCloseMesh(InpMsh);

	GetTim(&timer);
	printf("%g s\n", timer);
	printf("\nInput mesh : version = %d, %zd vertices, %zd triangles, %zd tets\n", \
			msh->MshVer, msh->NmbVer, msh->NmbTri, msh->NmbTet);
}


/*----------------------------------------------------------*/
/* Write mesh 												*/
/*----------------------------------------------------------*/

void RecMsh(char *OutNam, MshSct *msh)
{
	int OutMsh, ref=0;
	double timer = 0.;

	printf("Writing mesh          : ");
	GetTim(&timer);

	/* Create the output mesh */

	if(	!msh->NmbVer || !msh->NmbEdg || !msh->NmbTet \
	||	!(OutMsh = GmfOpenMesh(OutNam, GmfWrite, msh->MshVer, 3)) )
	{
		printf("Cannot create mesh %s\n", OutNam);
		exit(1);
	}

	/* Save the vertices from the input mesh */

	GmfSetKwd(OutMsh, GmfVertices, msh->NmbVer);
	GmfSetBlock(OutMsh, GmfVertices, \
				GmfDouble, &msh->ver[1].crd[0], &msh->ver[2].crd[0], \
				GmfDouble, &msh->ver[1].crd[1], &msh->ver[2].crd[1], \
				GmfDouble, &msh->ver[1].crd[2], &msh->ver[2].crd[2], \
				GmfInt, &msh->ver[1].ref, &msh->ver[2].ref);

	/* Save the extracted edges */

	GmfSetKwd(OutMsh, GmfEdges, msh->NmbEdg);
	GmfSetBlock(OutMsh, GmfEdges, \
				LngTyp, &msh->edg[1].idx[0], &msh->edg[2].idx[0], \
				LngTyp, &msh->edg[1].idx[1], &msh->edg[2].idx[1], \
				GmfInt, &ref, &ref);

	/* Save the triangles from the input mesh */

	GmfSetKwd(OutMsh, GmfTriangles, msh->NmbTri);
	GmfSetBlock(OutMsh, GmfTriangles, \
				LngTyp, &msh->tri[1].idx[0], &msh->tri[2].idx[0], \
				LngTyp, &msh->tri[1].idx[1], &msh->tri[2].idx[1], \
				LngTyp, &msh->tri[1].idx[2], &msh->tri[2].idx[2], \
				GmfInt, &msh->tri[1].ref, &msh->tri[2].ref);

	/* Save the tetrahedra from the input mesh */

	GmfSetKwd(OutMsh, GmfTetrahedra, msh->NmbTet);
	GmfSetBlock(OutMsh, GmfTetrahedra, \
				LngTyp, &msh->tet[1].idx[0], &msh->tet[2].idx[0], \
				LngTyp, &msh->tet[1].idx[1], &msh->tet[2].idx[1], \
				LngTyp, &msh->tet[1].idx[2], &msh->tet[2].idx[2], \
				LngTyp, &msh->tet[1].idx[3], &msh->tet[2].idx[3], \
				GmfInt, &ref, &ref);

	GmfCloseMesh(OutMsh);

	GetTim(&timer);
	printf("%g s\n\n", timer);
}
