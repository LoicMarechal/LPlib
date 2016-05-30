

/*----------------------------------------------------------*/
/*															*/
/*			PARALLEL NEIGHBOURS USING LP3					*/
/*															*/
/*----------------------------------------------------------*/
/*															*/
/*	Description:		extract a triangulated surface mesh */
/*                      from a volumic tetrahedral mesh		*/
/*	Author:				Loic MARECHAL						*/
/*	Creation date:		mar 11 2010							*/
/*	Last modification:	jan 21 2015							*/
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


/*----------------------------------------------------------*/
/* Structures												*/
/*----------------------------------------------------------*/

typedef struct
{
	double crd[3];
	int idx, ref;
}VerSct;

typedef struct
{
	int idx[3], ref;
}TriSct;

typedef struct
{
	int idx[4], ref;
}TetSct;

typedef struct
{
	int tet, nex;
	char voy, min, mid, max;
}HshSct;

typedef struct
{
	int NmbVer, NmbTri, NmbTet, MshVer;
	VerSct *ver;
	TriSct *tri;
	TetSct *tet;
}MshSct;

typedef struct
{
	char *FlgTab;
	int beg, end, HshSiz, HshPos, ColPos, NmbCpu, (*NgbTab)[4];
	HshSct *tab;
	MshSct *msh;
}ParSct;


/*----------------------------------------------------------*/
/* Global variables											*/
/*----------------------------------------------------------*/

struct timeval tp;
int VerTyp, TriTyp, TetTyp;


/*----------------------------------------------------------*/
/* Prototypes of local procedures							*/
/*----------------------------------------------------------*/

static void ScaMsh(char *, MshSct *);
static void RecMsh(char *, MshSct *);
static void GetTim(double *);
static void SetNgb(MshSct *, int);
static void ParNgb1(int, int, int, ParSct *);
static void ParNgb2(int, int, int, ParSct *);


/*----------------------------------------------------------*/
/* Read the volume, extract the surface and write the mesh	*/
/*----------------------------------------------------------*/

int main(int ArgCnt, char **ArgVec)
{
	char *PtrArg, *TmpStr, InpNam[1000], OutNam[1000];
	int i, j, LibParIdx, NmbCpu = 0;
	double timer=0;
	TetSct *tet;
	MshSct msh;

	/* Command line parsing */

	memset(&msh, 0, sizeof(MshSct));

	if(ArgCnt == 1)
	{
		puts("\ntetrahedra_neighbours v1.02 jan 20 2015   Loic MARECHAL / INRIA");
		puts(" Usage       : tetrahedra_neighbours -in volume_mesh -out surface_mesh");
		puts(" -in name    : name of the input tetrahedral-only mesh");
		puts(" -out name   : name of the output surface mesh");
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
	printf("\nReading mesh          : ");
	timer = 0.;
	GetTim(&timer);
	ScaMsh(InpNam, &msh);
	GetTim(&timer);
	printf("%g s\n", timer);
	printf("\nInput mesh : version = %d, %d vertices, %d tets\n", msh.MshVer, msh.NmbVer, msh.NmbTet);

	// Setup LP3 lib and datatypes
	LibParIdx = InitParallel(NmbCpu);
	VerTyp = NewType(LibParIdx, msh.NmbVer);
	TetTyp = NewType(LibParIdx, msh.NmbTet);

	// Launch the parallel neeighbour procedure
	SetNgb(&msh, LibParIdx);

	// Mesh writing
	printf("Writing mesh          : ");
	timer = 0.;
	GetTim(&timer);
	RecMsh(OutNam, &msh);
	GetTim(&timer);
	printf("%g s\n\n", timer);

	/* Release memory and stop LP3 */

	if(msh.ver)
		free(msh.ver);

	if(msh.tri)
		free(msh.tri);

	if(msh.tet)
		free(msh.tet);

	StopParallel(LibParIdx);
}


/*----------------------------------------------------------*/
/* Wall clock timer											*/
/*----------------------------------------------------------*/

static void GetTim(double *timer)
{
	gettimeofday(&tp, NULL);
	*timer = tp.tv_sec + tp.tv_usec / 1000000. - *timer;
}


/*----------------------------------------------------------*/
/* Read mesh 												*/
/*----------------------------------------------------------*/

static void ScaMsh(char *InpNam, MshSct *msh)
{
	int i, j, InpMsh, dim;
	float flt[3];

	/* Check mesh format */

	if(!(InpMsh = GmfOpenMesh(InpNam, GmfRead, &msh->MshVer, &dim)))
	{
		printf("Cannot open mesh %s\n", InpNam);
		exit(1);
	}

	if(dim != 3)
	{
		puts("Can only handle 3D meshes\n");
		exit(1);
	}

	/* Get stats and allocate tables */

	if((msh->NmbVer = GmfStatKwd(InpMsh, GmfVertices)))
		msh->ver = malloc((msh->NmbVer+1) * sizeof(VerSct));
	else
	{
		puts("Cannot renumber a mesh without vertices");
		exit(1);
	}

	if((msh->NmbTet = GmfStatKwd(InpMsh, GmfTetrahedra)))
		msh->tet = malloc((msh->NmbTet+1) * sizeof(TetSct));

	/* Read vertices */

	if(msh->NmbVer)
	{
		GmfGotoKwd(InpMsh, GmfVertices);
		GmfGetBlock(InpMsh, GmfVertices, \
					GmfDouble, &msh->ver[1].crd[0], &msh->ver[2].crd[0], \
					GmfDouble, &msh->ver[1].crd[1], &msh->ver[2].crd[1], \
					GmfDouble, &msh->ver[1].crd[2], &msh->ver[2].crd[2], \
					GmfInt, &msh->ver[1].ref, &msh->ver[2].ref);
	}

	/* Read tetrahedra */

	if(msh->NmbTet)
	{
		GmfGotoKwd(InpMsh, GmfTetrahedra);
		GmfGetBlock(InpMsh, GmfTetrahedra, \
					GmfInt, &msh->tet[1].idx[0], &msh->tet[2].idx[0], \
					GmfInt, &msh->tet[1].idx[1], &msh->tet[2].idx[1], \
					GmfInt, &msh->tet[1].idx[2], &msh->tet[2].idx[2], \
					GmfInt, &msh->tet[1].idx[3], &msh->tet[2].idx[3], \
					GmfInt, &msh->tet[1].ref, &msh->tet[2].ref);
	}

	GmfCloseMesh(InpMsh);
}


/*----------------------------------------------------------*/
/* Write mesh 												*/
/*----------------------------------------------------------*/

static void RecMsh(char *OutNam, MshSct *msh)
{
	int i, OutMsh;

	/* Create the output mesh */

	if(!(OutMsh = GmfOpenMesh(OutNam, GmfWrite, msh->MshVer, 3)))
	{
		printf("Cannot create mesh %s\n", OutNam);
		exit(1);
	}

	/* Save the vertices from the input mesh */
	//puts("1");
	if(msh->NmbVer)
	{
		GmfSetKwd(OutMsh, GmfVertices, msh->NmbVer);
		GmfSetBlock(OutMsh, GmfVertices, \
					GmfDouble, &msh->ver[1].crd[0], &msh->ver[2].crd[0], \
					GmfDouble, &msh->ver[1].crd[1], &msh->ver[2].crd[1], \
					GmfDouble, &msh->ver[1].crd[2], &msh->ver[2].crd[2], \
					GmfInt, &msh->ver[1].ref, &msh->ver[2].ref);
	}

	/* Save the extracted triangles */

	//puts("2");
	if(msh->NmbTri)
	{
		GmfSetKwd(OutMsh, GmfTriangles, msh->NmbTri);
		GmfSetBlock(OutMsh, GmfTriangles, \
					GmfInt, &msh->tri[1].idx[0], &msh->tri[2].idx[0], \
					GmfInt, &msh->tri[1].idx[1], &msh->tri[2].idx[1], \
					GmfInt, &msh->tri[1].idx[2], &msh->tri[2].idx[2], \
					GmfInt, &msh->tri[1].ref, &msh->tri[2].ref);
	}

	/* Save the tetrahedra from the input mesh */

	//puts("3");
	if(msh->NmbTet)
	{
		GmfSetKwd(OutMsh, GmfTetrahedra, msh->NmbTet);
		GmfSetBlock(OutMsh, GmfTetrahedra, \
					GmfInt, &msh->tet[1].idx[0], &msh->tet[2].idx[0], \
					GmfInt, &msh->tet[1].idx[1], &msh->tet[2].idx[1], \
					GmfInt, &msh->tet[1].idx[2], &msh->tet[2].idx[2], \
					GmfInt, &msh->tet[1].idx[3], &msh->tet[2].idx[3], \
					GmfInt, &msh->tet[1].ref, &msh->tet[2].ref);
	}

	GmfCloseMesh(OutMsh);
}


/*----------------------------------------------------------*/
/* Parallel neighbours between tets							*/
/*----------------------------------------------------------*/

static void SetNgb(MshSct *msh, int LibParIdx)
{
	char *FlgTab;
	int i, j, k, NgbIdx, NmbCpu, NmbTyp, HshSiz, MshSiz, (*NgbTab)[4];
	int tvpf[4][3] = { {1,2,3}, {2,0,3}, {3,0,1}, {0,2,1} };
	double timer;
	ParSct par[ MaxPth ];

	printf("Tet neighbours        : ");
	timer = 0.;
	GetTim(&timer);

	/* Allocate a hash table and overflow buffer */

	GetLplibInformation(LibParIdx, &NmbCpu, &NmbTyp);
	HshSiz =  (msh->NmbTet * 2)  / NmbCpu;
	MshSiz =  msh->NmbTet / NmbCpu;
	FlgTab = calloc( (msh->NmbTet+1), sizeof(char) );
	NgbTab = calloc( (msh->NmbTet+1), 4 * sizeof(int) );

	/* Setup parallel parameters */

	for(i=0;i<NmbCpu;i++)
	{
		par[i].beg = i * MshSiz + 1;
		par[i].end = (i+1) * MshSiz;
		par[i].HshSiz = HshSiz;
		par[i].ColPos = HshSiz;
		par[i].msh = msh;
		par[i].FlgTab = FlgTab;
		par[i].NmbCpu = NmbCpu;
		par[i].NgbTab = NgbTab;
	}

	par[ NmbCpu-1 ].end = msh->NmbTet;

	/* Launch parallel loops: the first one build local neighbours among each subdomains 
		and the second one build neighbourhood information between cross block elements */

	LaunchParallel(LibParIdx, TetTyp, 0, (void *)ParNgb1, (void *)par);

	if(NmbCpu > 1)
		LaunchParallel(LibParIdx, TetTyp, 0, (void *)ParNgb2, (void *)par);

	GetTim(&timer);
	printf("%g s\n", timer);

	/* Compute the number of unique triangles */

	msh->NmbTri = 0;

	for(i=1;i<=msh->NmbTet;i++)
		for(j=0;j<4;j++)
			if( !(NgbIdx = NgbTab[i][j]) || ( (msh->tet[i].ref != msh->tet[ NgbIdx ].ref) && (i > NgbIdx) ) )
				msh->NmbTri++;

	/* Add a triangle table to the mesh */

	msh->tri = malloc( (msh->NmbTri+1) * sizeof(TriSct) );
	msh->NmbTri = 0;

	/* Then fill it with the boundary triangles */

	for(i=1;i<=msh->NmbTet;i++)
		for(j=0;j<4;j++)
			if( !(NgbIdx = NgbTab[i][j]) || ( (msh->tet[i].ref != msh->tet[ NgbIdx ].ref) && (i > NgbIdx) ) )
			{
				msh->NmbTri++;

				for(k=0;k<3;k++)
					msh->tri[ msh->NmbTri ].idx[k] = msh->tet[i].idx[ tvpf[j][k] ];

				if(!NgbIdx)
					msh->tri[ msh->NmbTri ].ref = 0;
				else
					msh->tri[ msh->NmbTri ].ref = 1;
			}

	free(NgbTab);
	free(FlgTab);

	for(i=0;i<NmbCpu;i++)
		free(par[i].tab);

	printf("%d boundary triangles extracted\n", msh->NmbTri);
}


/*----------------------------------------------------------*/
/* Set links between tets from this local subdomain			*/
/*----------------------------------------------------------*/

static void ParNgb1(int BegIdx, int EndIdx, int c, ParSct *par)
{
	char *FlgTab = par[c].FlgTab;
	int i, j, k, key, (*NgbTab)[4] = par[c].NgbTab;
	unsigned int min, mid, max;
	TetSct *tet, *ngb;
	MshSct *msh = par[c].msh;
	HshSct *tab = par[c].tab = calloc(2 * par[c].HshSiz, sizeof(HshSct));

	/* Allocate a local hash table and loop over the local elements */

	for(i=par[c].beg; i<=par[c].end; i++)
	{
		tet = &msh->tet[i];

		for(j=0;j<4;j++)
		{
			/* Compute the hashing key from the face's vertices indices */

			min = max = (j+1)%4;

			for(k=0;k<4;k++)
				if(k != j)
				{
					if(tet->idx[k] < tet->idx[ min ])
						min = k;
					else if(tet->idx[k] > tet->idx[ max ])
						max = k;
				}

			mid = 6 - min - max - j;
			key = (3 * tet->idx[ min ] + 5 * tet->idx[ mid ] + 7 * tet->idx[ max ]) % par[c].HshSiz;

			/* If the bucket is empty, store the face */

			if(!tab[ key ].tet)
			{
				tab[ key ].tet = i;
				tab[ key ].voy = j;
				tab[ key ].min = min;
				tab[ key ].mid = mid;
				tab[ key ].max = max;
				continue;
			}

			/* Otherwise, search through the linked list */

			do
			{
				ngb = &msh->tet[ tab[ key ].tet ];

				/* If the same face is found in the hash table, setup a link between both tetrahedra */

				if( (ngb->idx[ tab[ key ].min ] == tet->idx[ min ]) \
				&& 	(ngb->idx[ tab[ key ].mid ] == tet->idx[ mid ]) \
				&& 	(ngb->idx[ tab[ key ].max ] == tet->idx[ max ]) )
				{
					NgbTab[i][j] = tab[ key ].tet;
					FlgTab[i]++;
					NgbTab[ tab[ key ].tet ][ tab[ key ].voy ] = i;
					FlgTab[ tab[ key ].tet ]++;
					break;
				}

				/* If not, allocate a new bucket from the overflow table and link it to the main entry */

				if(tab[ key ].nex)
					key = tab[ key ].nex;
				else
				{
					tab[ key ].nex = par[c].ColPos;
					key = par[c].ColPos++;
					tab[ key ].tet = i;
					tab[ key ].voy = j;
					tab[ key ].min = min;
					tab[ key ].mid = mid;
					tab[ key ].max = max;
					break;
				}
			}while(1);
		}
	}
}


/*------------------------------------------------------------*/
/* Setup the missing links between tets that cross subdomains */
/*------------------------------------------------------------*/

static void ParNgb2(int BegIdx, int EndIdx, int c, ParSct *par)
{
	int n, i, j, k, key, BasKey, flg, (*NgbTab)[4] = par[c].NgbTab;
	unsigned int min, mid, max;
	TetSct *tet, *ngb;
	HshSct *tab;
	MshSct *msh = par[c].msh;

	for(i=par[c].beg; i<=par[c].end; i++)
	{
		/* If a tetrahedron has already 4 links, there is no need to find a missing ones */

		if(par[c].FlgTab[i] == 4)
			continue;

		tet = &msh->tet[i];

		for(j=0;j<4;j++)
		{
			/* If there is no neighbour through this face, try to find on among other subdomains local hash tables */

			if(NgbTab[i][j])
				continue;

			min = max = (j+1)%4;

			for(k=0;k<4;k++)
				if(k != j)
				{
					if(tet->idx[k] < tet->idx[ min ])
						min = k;
					else if(tet->idx[k] > tet->idx[ max ])
						max = k;
				}

			mid = 6 - min - max - j;
			flg = 0;
			BasKey = (3 * tet->idx[ min ] + 5 * tet->idx[ mid ] + 7 * tet->idx[ max ]) % par[c].HshSiz;
    
			for(n=0; n<par[c].NmbCpu; n++)
			{
				if(n == c)
					continue;

				tab = par[n].tab;
				key = BasKey;

				do
				{
					ngb = &msh->tet[ tab[ key ].tet ];
    
					if( (ngb->idx[ tab[ key ].min ] == tet->idx[ min ]) \
					&& 	(ngb->idx[ tab[ key ].mid ] == tet->idx[ mid ]) \
					&& 	(ngb->idx[ tab[ key ].max ] == tet->idx[ max ]) )
					{
						NgbTab[i][j] = tab[ key ].tet;
						flg = 1;
						break;
					}

					if(tab[ key ].nex)
						key = tab[ key ].nex;
					else
						break;
				}while(1);

				if(flg)
					break;
			}
		}
	}
}
