

/*----------------------------------------------------------*/
/*															*/
/*						HILBERT V 1.10						*/
/*															*/
/*----------------------------------------------------------*/
/*															*/
/*	Description:		renumber .meshb files				*/
/*	Author:				Loic MARECHAL						*/
/*	Creation date:		mar 11 2010							*/
/*	Last modification:	nov 09 2015							*/
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

#define MaxItr 21
#define min(a,b) ((a) < (b) ? (a) : (b))
#define max(a,b) ((a) > (b) ? (a) : (b))
#define BufSiz 250
#define MinVec 2
#define MaxVec 64


/*----------------------------------------------------------*/
/* Structures												*/
/*----------------------------------------------------------*/

typedef unsigned long long int u64;

typedef struct
{
	u64 cod;
	double crd[3];
	int idx, ref;
}VerSct;

typedef struct
{
	u64 cod;
	int idx[2], ref;
}EdgSct;

typedef struct
{
	u64 cod;
	int idx[3], ref;
}TriSct;

typedef struct
{
	u64 cod;
	int idx[4], ref;
}QadSct;

typedef struct
{
	u64 cod;
	int idx[4], ref;
}TetSct;

typedef struct
{
	u64 cod;
	int idx[5], ref;
}PyrSct;

typedef struct
{
	u64 cod;
	int idx[6], ref;
}PriSct;

typedef struct
{
	u64 cod;
	int idx[8], ref;
}HexSct;

typedef struct
{
	int tet, nex;
	char voy, min, mid, max;
}HshSct;

typedef struct
{
	int NmbVer, NmbEdg, NmbTri, NmbQad, NmbTet, NmbPyr, NmbPri, NmbHex, *Old2New, MshVer, EleTyp;
	double box[6];
	VerSct *ver;
	EdgSct *edg;
	TriSct *tri;
	QadSct *qad;
	TetSct *tet;
	PyrSct *pyr;
	PriSct *pri;
	HexSct *hex;
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
int VerTyp, EdgTyp, TriTyp, QadTyp, TetTyp, PyrTyp, PriTyp, HexTyp;


/*----------------------------------------------------------*/
/* Prototypes of local procedures							*/
/*----------------------------------------------------------*/

void ScaMsh(char *, MshSct *);
void RecMsh(char *, MshSct *);
u64 hilbert(double *, double *, int);
static void GetTim(double *);
int CmpFnc(const void *, const void *);
static void RenVer(int, int, int, MshSct *);
static void RenEle(int, int, int, MshSct *);
void PrtSta(MshSct *, int);
static void SetNgb(MshSct *, int);
static void ParNgb1(int, int, int, ParSct *);
static void ParNgb2(int, int, int, ParSct *);
void VecRen(MshSct *, int);


/*----------------------------------------------------------*/
/* Read, renumber through a Hilbert SFC and write the mesh	*/
/*----------------------------------------------------------*/

int main(int ArgCnt, char **ArgVec)
{
	char *PtrArg, *TmpStr, InpNam[1000], OutNam[1000];
	int i, j, LibParIdx, NmbCpu = 0, StaFlg=0, VecSiz=0;
	float flt[3], sta[2];
	double timer=0;
	TetSct *tet;
	MshSct msh;

	/* Command line parsing */

	memset(&msh, 0, sizeof(MshSct));

	if(ArgCnt == 1)
	{
		puts("\nHILBERT v1.10 nov 09 2015   Loic MARECHAL / INRIA");
		puts(" Usage       : hilbert -in input_mesh -out renumbered_mesh");
		puts(" -in name    : name of the input mesh");
		puts(" -out name   : name of the output renumbered mesh");
		puts(" -vec size   : optional vector safe post renumbering");
		puts(" -stats      : print element blocks dependencies stats before and after renumbering");
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

		if(!strcmp(PtrArg,"-vec"))
		{
			VecSiz = atoi(*++ArgVec);
			VecSiz = max(VecSiz, MinVec);
			VecSiz = min(VecSiz, MaxVec);
			ArgCnt--;
			continue;
		}

		if(!strcmp(PtrArg,"-stats"))
		{
			StaFlg = 1;
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

	/* Mesh reading */

	printf("\nReading mesh          : ");
	timer = 0.;
	GetTim(&timer);
	ScaMsh(InpNam, &msh);
	GetTim(&timer);
	printf("%g s\n", timer);

	printf("\nInput mesh : version = %d, %d vertices, %d edges, %d triangles, %d quads, %d tets, %d hexes\n",
	 		msh.MshVer, msh.NmbVer, msh.NmbEdg, msh.NmbTri, msh.NmbQad, msh.NmbTet, msh.NmbHex);

	/* Compute initial stats */

	LibParIdx = InitParallel(NmbCpu);

	if(StaFlg)
	{
		VerTyp = NewType(LibParIdx, msh.NmbVer);
		EdgTyp = NewType(LibParIdx, msh.NmbEdg);
		TriTyp = NewType(LibParIdx, msh.NmbTri);
		QadTyp = NewType(LibParIdx, msh.NmbQad);
		TetTyp = NewType(LibParIdx, msh.NmbTet);
		PyrTyp = NewType(LibParIdx, msh.NmbPyr);
		PriTyp = NewType(LibParIdx, msh.NmbPri);
		HexTyp = NewType(LibParIdx, msh.NmbHex);
		puts("\nDependencies before renumbering (average / max) :");
		PrtSta(&msh, LibParIdx);
	}

	/* Vertices renumbering */

	printf("Renumbering vertices  : ");
	timer = 0.;
	GetTim(&timer);
	VerTyp = NewType(LibParIdx, msh.NmbVer);
	LaunchParallel(LibParIdx, VerTyp, 0, (void *)RenVer, (void *)&msh);
	ParallelQsort(LibParIdx, &msh.ver[1], msh.NmbVer, sizeof(VerSct), CmpFnc);

	msh.Old2New = malloc( (msh.NmbVer+1) * sizeof(int) );

	for(i=1;i<=msh.NmbVer;i++)
		msh.Old2New[ msh.ver[i].idx ] = i;

	GetTim(&timer);
	printf("%g s\n", timer);

	/* Edges renumbering */

	if(msh.NmbEdg)
	{
		printf("Renumbering edges     : ");
		timer = 0.;
		GetTim(&timer);
		EdgTyp = NewType(LibParIdx, msh.NmbEdg);
		msh.EleTyp = GmfEdges;
		LaunchParallel(LibParIdx, EdgTyp, 0, (void *)RenEle, (void *)&msh);
		ParallelQsort(LibParIdx, &msh.edg[1], msh.NmbEdg, sizeof(EdgSct), CmpFnc);
		GetTim(&timer);
		printf("%g s\n", timer);
	}

	/* Triangles renumbering */

	if(msh.NmbTri)
	{
		printf("Renumbering triangles : ");
		timer = 0.;
		GetTim(&timer);
		TriTyp = NewType(LibParIdx, msh.NmbTri);
		msh.EleTyp = GmfTriangles;
		LaunchParallel(LibParIdx, TriTyp, 0, (void *)RenEle, (void *)&msh);
		ParallelQsort(LibParIdx, &msh.tri[1], msh.NmbTri, sizeof(TriSct), CmpFnc);
		GetTim(&timer);
		printf("%g s\n", timer);
	}

	/* Quads renumbering */

	if(msh.NmbQad)
	{
		printf("Renumbering quads     : ");
		timer = 0.;
		GetTim(&timer);
		QadTyp = NewType(LibParIdx, msh.NmbQad);
		msh.EleTyp = GmfQuadrilaterals;
		LaunchParallel(LibParIdx, QadTyp, 0, (void *)RenEle, (void *)&msh);
		ParallelQsort(LibParIdx, &msh.qad[1], msh.NmbQad, sizeof(QadSct), CmpFnc);
		GetTim(&timer);
		printf("%g s\n", timer);
	}

	/* Tets renumbering */

	if(msh.NmbTet)
	{
		printf("Renumbering tets      : ");
		timer = 0.;
		GetTim(&timer);
		TetTyp = NewType(LibParIdx, msh.NmbTet);
		msh.EleTyp = GmfTetrahedra;
		LaunchParallel(LibParIdx, TetTyp, 0, (void *)RenEle, (void *)&msh);
		ParallelQsort(LibParIdx, &msh.tet[1], msh.NmbTet, sizeof(TetSct), CmpFnc);
		GetTim(&timer);
		printf("%g s\n", timer);

		if(VecSiz)
		{
			VecRen(&msh, VecSiz);
			ParallelQsort(LibParIdx, &msh.tet[1], msh.NmbTet, sizeof(TetSct), CmpFnc);
		}
	}

	/* Pyramids renumbering */

	if(msh.NmbPyr)
	{
		printf("Renumbering pyramids  : ");
		timer = 0.;
		GetTim(&timer);
		PyrTyp = NewType(LibParIdx, msh.NmbPyr);
		msh.EleTyp = GmfPyramids;
		LaunchParallel(LibParIdx, PyrTyp, 0, (void *)RenEle, (void *)&msh);
		ParallelQsort(LibParIdx, &msh.pyr[1], msh.NmbPyr, sizeof(PyrSct), CmpFnc);
		GetTim(&timer);
		printf("%g s\n", timer);
	}

	/* Prisms renumbering */

	if(msh.NmbPri)
	{
		printf("Renumbering prisms    : ");
		timer = 0.;
		GetTim(&timer);
		PriTyp = NewType(LibParIdx, msh.NmbPri);
		msh.EleTyp = GmfPrisms;
		LaunchParallel(LibParIdx, PriTyp, 0, (void *)RenEle, (void *)&msh);
		ParallelQsort(LibParIdx, &msh.pri[1], msh.NmbPri, sizeof(PriSct), CmpFnc);
		GetTim(&timer);
		printf("%g s\n", timer);
	}

	/* Hexes renumbering */

	if(msh.NmbHex)
	{
		printf("Renumbering hexes     : ");
		timer = 0.;
		GetTim(&timer);
		HexTyp = NewType(LibParIdx, msh.NmbHex);
		msh.EleTyp = GmfHexahedra;
		LaunchParallel(LibParIdx, HexTyp, 0, (void *)RenEle, (void *)&msh);
		ParallelQsort(LibParIdx, &msh.hex[1], msh.NmbHex, sizeof(HexSct), CmpFnc);
		GetTim(&timer);
		printf("%g s\n", timer);
	}

	/* Compute dependencies */

	if(StaFlg)
	{
		puts("\nDependencies after renumbering (average / max) :");
		PrtSta(&msh, LibParIdx);
	}

//	SetNgb(&msh, LibParIdx);

	StopParallel(LibParIdx);

	/* Mesh writing */

	printf("Writing mesh          : ");
	timer = 0.;
	GetTim(&timer);
	RecMsh(OutNam, &msh);
	GetTim(&timer);
	printf("%g s\n\n", timer);

	/* Release memory */

	if(msh.ver)
		free(msh.ver);

	if(msh.edg)
		free(msh.edg);

	if(msh.tri)
		free(msh.tri);

	if(msh.qad)
		free(msh.qad);

	if(msh.tet)
		free(msh.tet);

	if(msh.pyr)
		free(msh.pyr);

	if(msh.pri)
		free(msh.pri);

	if(msh.hex)
		free(msh.hex);

	if(msh.Old2New)
		free(msh.Old2New);
}


/*----------------------------------------------------------*/
/* Read mesh 												*/
/*----------------------------------------------------------*/

void ScaMsh(char *InpNam, MshSct *msh)
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

	if((msh->NmbEdg = GmfStatKwd(InpMsh, GmfEdges)))
		msh->edg = malloc((msh->NmbEdg+1) * sizeof(EdgSct));

	if((msh->NmbTri = GmfStatKwd(InpMsh, GmfTriangles)))
		msh->tri = malloc((msh->NmbTri+1) * sizeof(TriSct));

	if((msh->NmbQad = GmfStatKwd(InpMsh, GmfQuadrilaterals)))
		msh->qad = malloc((msh->NmbQad+1) * sizeof(QadSct));

	if((msh->NmbTet = GmfStatKwd(InpMsh, GmfTetrahedra)))
		msh->tet = malloc((msh->NmbTet+1) * sizeof(TetSct));

	if((msh->NmbPyr = GmfStatKwd(InpMsh, GmfPyramids)))
		msh->pyr = malloc((msh->NmbPyr+1) * sizeof(PyrSct));

	if((msh->NmbPri = GmfStatKwd(InpMsh, GmfPrisms)))
		msh->pri = malloc((msh->NmbPri+1) * sizeof(PriSct));

	if((msh->NmbHex = GmfStatKwd(InpMsh, GmfHexahedra)))
		msh->hex = malloc((msh->NmbHex+1) * sizeof(HexSct));

	/* Read fields */

	if(msh->NmbVer)
	{
		GmfGotoKwd(InpMsh, GmfVertices);
		GmfGetBlock(InpMsh, GmfVertices, \
					GmfDouble, &msh->ver[1].crd[0], &msh->ver[2].crd[0], \
					GmfDouble, &msh->ver[1].crd[1], &msh->ver[2].crd[1], \
					GmfDouble, &msh->ver[1].crd[2], &msh->ver[2].crd[2], \
					GmfInt, &msh->ver[1].ref, &msh->ver[2].ref);

		for(i=1;i<=msh->NmbVer;i++)
		{
			msh->ver[i].idx = i;

			if(i==1)
				for(j=0;j<3;j++)
					msh->box[j] = msh->box[j+3] = msh->ver[i].crd[j];
			else
				for(j=0;j<3;j++)
					if(msh->ver[i].crd[j] < msh->box[j])
						msh->box[j] = msh->ver[i].crd[j];
					else if(msh->ver[i].crd[j] > msh->box[j+3])
						msh->box[j+3] = msh->ver[i].crd[j];
		}

		for(j=0;j<3;j++)
			msh->box[j+3] = pow(2,64) / (msh->box[j+3] - msh->box[j]);
	}

	if(msh->NmbEdg)
	{
		GmfGotoKwd(InpMsh, GmfEdges);
		GmfGetBlock(InpMsh, GmfEdges, \
					GmfInt, &msh->edg[1].idx[0], &msh->edg[2].idx[0], \
					GmfInt, &msh->edg[1].idx[1], &msh->edg[2].idx[1], \
					GmfInt, &msh->edg[1].ref, &msh->edg[2].ref);
	}

	if(msh->NmbTri)
	{
		GmfGotoKwd(InpMsh, GmfTriangles);
		GmfGetBlock(InpMsh, GmfTriangles, \
					GmfInt, &msh->tri[1].idx[0], &msh->tri[2].idx[0], \
					GmfInt, &msh->tri[1].idx[1], &msh->tri[2].idx[1], \
					GmfInt, &msh->tri[1].idx[2], &msh->tri[2].idx[2], \
					GmfInt, &msh->tri[1].ref, &msh->tri[2].ref);
	}

	if(msh->NmbQad)
	{
		GmfGotoKwd(InpMsh, GmfQuadrilaterals);
		GmfGetBlock(InpMsh, GmfQuadrilaterals, \
					GmfInt, &msh->qad[1].idx[0], &msh->qad[2].idx[0], \
					GmfInt, &msh->qad[1].idx[1], &msh->qad[2].idx[1], \
					GmfInt, &msh->qad[1].idx[2], &msh->qad[2].idx[2], \
					GmfInt, &msh->qad[1].idx[3], &msh->qad[2].idx[3], \
					GmfInt, &msh->qad[1].ref, &msh->qad[2].ref);
	}

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

	if(msh->NmbPyr)
	{
		GmfGotoKwd(InpMsh, GmfPyramids);
		GmfGetBlock(InpMsh, GmfPyramids, \
					GmfInt, &msh->pyr[1].idx[0], &msh->pyr[2].idx[0], \
					GmfInt, &msh->pyr[1].idx[1], &msh->pyr[2].idx[1], \
					GmfInt, &msh->pyr[1].idx[2], &msh->pyr[2].idx[2], \
					GmfInt, &msh->pyr[1].idx[3], &msh->pyr[2].idx[3], \
					GmfInt, &msh->pyr[1].idx[4], &msh->pyr[2].idx[4], \
					GmfInt, &msh->pyr[1].ref, &msh->pyr[2].ref);
	}

	if(msh->NmbPri)
	{
		GmfGotoKwd(InpMsh, GmfPrisms);
		GmfGetBlock(InpMsh, GmfPrisms, \
					GmfInt, &msh->pri[1].idx[0], &msh->pri[2].idx[0], \
					GmfInt, &msh->pri[1].idx[1], &msh->pri[2].idx[1], \
					GmfInt, &msh->pri[1].idx[2], &msh->pri[2].idx[2], \
					GmfInt, &msh->pri[1].idx[3], &msh->pri[2].idx[3], \
					GmfInt, &msh->pri[1].idx[4], &msh->pri[2].idx[4], \
					GmfInt, &msh->pri[1].idx[5], &msh->pri[2].idx[5], \
					GmfInt, &msh->pri[1].ref, &msh->pri[2].ref);
	}

	if(msh->NmbHex)
	{
		GmfGotoKwd(InpMsh, GmfHexahedra);
		GmfGetBlock(InpMsh, GmfHexahedra, \
					GmfInt, &msh->hex[1].idx[0], &msh->hex[2].idx[0], \
					GmfInt, &msh->hex[1].idx[1], &msh->hex[2].idx[1], \
					GmfInt, &msh->hex[1].idx[2], &msh->hex[2].idx[2], \
					GmfInt, &msh->hex[1].idx[3], &msh->hex[2].idx[3], \
					GmfInt, &msh->hex[1].idx[4], &msh->hex[2].idx[4], \
					GmfInt, &msh->hex[1].idx[5], &msh->hex[2].idx[5], \
					GmfInt, &msh->hex[1].idx[6], &msh->hex[2].idx[6], \
					GmfInt, &msh->hex[1].idx[7], &msh->hex[2].idx[7], \
					GmfInt, &msh->hex[1].ref, &msh->hex[2].ref);
	}

	GmfCloseMesh(InpMsh);
}


/*----------------------------------------------------------*/
/* Write mesh 												*/
/*----------------------------------------------------------*/

void RecMsh(char *OutNam, MshSct *msh)
{
	int i, OutMsh, NewNmbVer;

	if(!(OutMsh = GmfOpenMesh(OutNam, GmfWrite, msh->MshVer, 3)))
	{
		printf("Cannot create mesh %s\n", OutNam);
		exit(1);
	}

	if(msh->NmbVer)
	{
		GmfSetKwd(OutMsh, GmfVertices, msh->NmbVer);
		GmfSetBlock(OutMsh, GmfVertices, \
					GmfDouble, &msh->ver[1].crd[0], &msh->ver[2].crd[0], \
					GmfDouble, &msh->ver[1].crd[1], &msh->ver[2].crd[1], \
					GmfDouble, &msh->ver[1].crd[2], &msh->ver[2].crd[2], \
					GmfInt, &msh->ver[1].ref, &msh->ver[2].ref);
	}

	if(msh->NmbEdg)
	{
		GmfSetKwd(OutMsh, GmfEdges, msh->NmbEdg);
		GmfSetBlock(OutMsh, GmfEdges, \
					GmfInt, &msh->edg[1].idx[0], &msh->edg[2].idx[0], \
					GmfInt, &msh->edg[1].idx[1], &msh->edg[2].idx[1], \
					GmfInt, &msh->edg[1].ref, &msh->edg[2].ref);
	}

	if(msh->NmbTri)
	{
		GmfSetKwd(OutMsh, GmfTriangles, msh->NmbTri);
		GmfSetBlock(OutMsh, GmfTriangles, \
					GmfInt, &msh->tri[1].idx[0], &msh->tri[2].idx[0], \
					GmfInt, &msh->tri[1].idx[1], &msh->tri[2].idx[1], \
					GmfInt, &msh->tri[1].idx[2], &msh->tri[2].idx[2], \
					GmfInt, &msh->tri[1].ref, &msh->tri[2].ref);
	}

	if(msh->NmbQad)
	{
		GmfSetKwd(OutMsh, GmfQuadrilaterals, msh->NmbQad);
		GmfSetBlock(OutMsh, GmfQuadrilaterals, \
					GmfInt, &msh->qad[1].idx[0], &msh->qad[2].idx[0], \
					GmfInt, &msh->qad[1].idx[1], &msh->qad[2].idx[1], \
					GmfInt, &msh->qad[1].idx[2], &msh->qad[2].idx[2], \
					GmfInt, &msh->qad[1].idx[3], &msh->qad[2].idx[3], \
					GmfInt, &msh->qad[1].ref, &msh->qad[2].ref);
	}

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

	if(msh->NmbPyr)
	{
		GmfSetKwd(OutMsh, GmfPyramids, msh->NmbPyr);
		GmfSetBlock(OutMsh, GmfPyramids, \
					GmfInt, &msh->pyr[1].idx[0], &msh->pyr[2].idx[0], \
					GmfInt, &msh->pyr[1].idx[1], &msh->pyr[2].idx[1], \
					GmfInt, &msh->pyr[1].idx[2], &msh->pyr[2].idx[2], \
					GmfInt, &msh->pyr[1].idx[3], &msh->pyr[2].idx[3], \
					GmfInt, &msh->pyr[1].idx[4], &msh->pyr[2].idx[4], \
					GmfInt, &msh->pyr[1].ref, &msh->pyr[2].ref);
	}

	if(msh->NmbPri)
	{
		GmfSetKwd(OutMsh, GmfPrisms, msh->NmbPri);
		GmfSetBlock(OutMsh, GmfPrisms, \
					GmfInt, &msh->pri[1].idx[0], &msh->pri[2].idx[0], \
					GmfInt, &msh->pri[1].idx[1], &msh->pri[2].idx[1], \
					GmfInt, &msh->pri[1].idx[2], &msh->pri[2].idx[2], \
					GmfInt, &msh->pri[1].idx[3], &msh->pri[2].idx[3], \
					GmfInt, &msh->pri[1].idx[4], &msh->pri[2].idx[4], \
					GmfInt, &msh->pri[1].idx[5], &msh->pri[2].idx[5], \
					GmfInt, &msh->pri[1].ref, &msh->pri[2].ref);
	}

	if(msh->NmbHex)
	{
		GmfSetKwd(OutMsh, GmfHexahedra, msh->NmbHex);
		GmfSetBlock(OutMsh, GmfHexahedra, \
					GmfInt, &msh->hex[1].idx[0], &msh->hex[2].idx[0], \
					GmfInt, &msh->hex[1].idx[1], &msh->hex[2].idx[1], \
					GmfInt, &msh->hex[1].idx[2], &msh->hex[2].idx[2], \
					GmfInt, &msh->hex[1].idx[3], &msh->hex[2].idx[3], \
					GmfInt, &msh->hex[1].idx[4], &msh->hex[2].idx[4], \
					GmfInt, &msh->hex[1].idx[5], &msh->hex[2].idx[5], \
					GmfInt, &msh->hex[1].idx[6], &msh->hex[2].idx[6], \
					GmfInt, &msh->hex[1].idx[7], &msh->hex[2].idx[7], \
					GmfInt, &msh->hex[1].ref, &msh->hex[2].ref);
	}

	GmfCloseMesh(OutMsh);
}


/*----------------------------------------------------------*/
/* Compute the hilbert code from 3d coordinates				*/
/*----------------------------------------------------------*/

u64 hilbert(double crd[3], double box[6], int itr)
{
	u64 IntCrd[3], m=1LL<<63, cod;
	int i, j, b, GeoWrd, NewWrd, BitTab[3] = {1,2,4};
	double TmpCrd[3];
	int rot[8], GeoCod[8]={0,3,7,4,1,2,6,5};  /* Z curve = {5,4,7,6,1,0,3,2} */
	int HilCod[8][8] = {{0,7,6,1,2,5,4,3}, {0,3,4,7,6,5,2,1}, {0,3,4,7,6,5,2,1}, {2,3,0,1,6,7,4,5},\
						{2,3,0,1,6,7,4,5}, {6,5,2,1,0,3,4,7}, {6,5,2,1,0,3,4,7}, {4,3,2,5,6,1,0,7}};

	/* Convert double precision coordinates to integers */
    
	for(j=0;j<3;j++)
	{
		TmpCrd[j] = (crd[j] - box[j]) * box[j+3];
		IntCrd[j] = TmpCrd[j];
	}
    
	/* Binary hilbert renumbering loop */
    
	cod = 0;
    
	for(j=0;j<8;j++)
		rot[j] = GeoCod[j];
    
	for(b=0;b<itr;b++)
	{
		GeoWrd = 0;
    
		for(j=0;j<3;j++)
		{
			if(IntCrd[j] & m)
				GeoWrd |= BitTab[j];
    
			IntCrd[j] = IntCrd[j]<<1;
		}
    
		NewWrd = rot[ GeoWrd ];
    
		cod = cod<<3 | NewWrd;
    
		for(j=0;j<8;j++)
			rot[j] = HilCod[ NewWrd ][ rot[j] ];
    
	}

	return(cod);
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
/* Comparison of two items for the qsort					*/
/*----------------------------------------------------------*/

int CmpFnc(const void *a, const void *b)
{
	u64 *pa = (u64 *)a, *pb = (u64 *)b;

	if(*pa > *pb)
		return(1);
	else
		return(-1);
}


/*----------------------------------------------------------*/
/* Parallel loop renumbering vertices						*/
/*----------------------------------------------------------*/

static void RenVer(int BegIdx, int EndIdx, int PthIdx, MshSct *msh)
{
	int i;

	for(i=BegIdx; i<=EndIdx; i++)
		msh->ver[i].cod = hilbert(msh->ver[i].crd, msh->box, MaxItr);
}


static void SetNewIdx(int NmbVer, int *IdxTab, int *Old2New)
{
	int i;

	for(i=0;i<NmbVer;i++)
		IdxTab[i] = Old2New[ IdxTab[i] ];
}

static void SetMidCrd(int NmbVer, int *IdxTab, MshSct *msh, double *crd)
{
	int i, j;

	for(j=0;j<3;j++)
		crd[j] = 0.;

	for(i=0;i<NmbVer;i++)
		for(j=0;j<3;j++)
			crd[j] += msh->ver[ IdxTab[i] ].crd[j];

	for(j=0;j<3;j++)
		crd[j] /= NmbVer;

}

/*----------------------------------------------------------*/
/* Parallel loop renumbering tets							*/
/*----------------------------------------------------------*/

static void RenEle(int BegIdx, int EndIdx, int PthIdx, MshSct *msh)
{
	int i, j;
	double crd[3];

	switch(msh->EleTyp)
	{
		case GmfEdges :
		{
			for(i=BegIdx; i<=EndIdx; i++)
			{
				SetNewIdx(2, msh->edg[i].idx, msh->Old2New);
				SetMidCrd(2, msh->edg[i].idx, msh, crd);
				msh->edg[i].cod = hilbert(crd, msh->box, MaxItr);
			}
		}break;

		case GmfTriangles :
		{
			for(i=BegIdx; i<=EndIdx; i++)
			{
				SetNewIdx(3, msh->tri[i].idx, msh->Old2New);
				SetMidCrd(3, msh->tri[i].idx, msh, crd);
				msh->tri[i].cod = hilbert(crd, msh->box, MaxItr);
			}
		}break;

		case GmfQuadrilaterals :
		{
			for(i=BegIdx; i<=EndIdx; i++)
			{
				SetNewIdx(4, msh->qad[i].idx, msh->Old2New);
				SetMidCrd(4, msh->qad[i].idx, msh, crd);
				msh->qad[i].cod = hilbert(crd, msh->box, MaxItr);
			}
		}break;

		case GmfTetrahedra :
		{
			for(i=BegIdx; i<=EndIdx; i++)
			{
				SetNewIdx(4, msh->tet[i].idx, msh->Old2New);
				SetMidCrd(4, msh->tet[i].idx, msh, crd);
				msh->tet[i].cod = hilbert(crd, msh->box, MaxItr);
			}
		}break;

		case GmfPyramids :
		{
			for(i=BegIdx; i<=EndIdx; i++)
			{
				SetNewIdx(5, msh->pyr[i].idx, msh->Old2New);
				SetMidCrd(5, msh->pyr[i].idx, msh, crd);
				msh->pyr[i].cod = hilbert(crd, msh->box, MaxItr);
			}
		}break;

		case GmfPrisms :
		{
			for(i=BegIdx; i<=EndIdx; i++)
			{
				SetNewIdx(6, msh->pri[i].idx, msh->Old2New);
				SetMidCrd(6, msh->pri[i].idx, msh, crd);
				msh->pri[i].cod = hilbert(crd, msh->box, MaxItr);
			}
		}break;

		case GmfHexahedra :
		{
			for(i=BegIdx; i<=EndIdx; i++)
			{
				SetNewIdx(8, msh->hex[i].idx, msh->Old2New);
				SetMidCrd(8, msh->hex[i].idx, msh, crd);
				msh->hex[i].cod = hilbert(crd, msh->box, MaxItr);
			}
		}break;

	}
}


/*----------------------------------------------------------*/
/* Compute and print elements / vertices dependencies stats	*/
/*----------------------------------------------------------*/

void PrtSta(MshSct *msh, int LibParIdx)
{
	int i, j, k, l, NmbLin, AvgLin, NewLin, NewHit, NmbHit, buf[ BufSiz ][8], lin[8], idx;
	float  sta[2];

	if(msh->NmbEdg)
	{
		BeginDependency(LibParIdx, EdgTyp, VerTyp);

		for(i=1;i<=msh->NmbEdg;i++)
			AddDependencyFast(LibParIdx, 1, &i, 2, msh->edg[i].idx);

		EndDependency(LibParIdx, sta);
		printf(" edges     : %3.2f%% / %3.2f%%\n",sta[0],sta[1]);
	}

	if(msh->NmbTri)
	{
		BeginDependency(LibParIdx, TriTyp, VerTyp);

		for(i=1;i<=msh->NmbTri;i++)
			AddDependencyFast(LibParIdx, 1, &i, 3, msh->tri[i].idx);

		EndDependency(LibParIdx, sta);
		printf(" triangles : %3.2f%% / %3.2f%%\n",sta[0],sta[1]);
	}

	if(msh->NmbQad)
	{
		BeginDependency(LibParIdx, QadTyp, VerTyp);

		for(i=1;i<=msh->NmbQad;i++)
			AddDependencyFast(LibParIdx, 1, &i, 4, msh->qad[i].idx);

		EndDependency(LibParIdx, sta);
		printf(" quads     : %3.2f%% / %3.2f%%\n",sta[0],sta[1]);
	}

	if(msh->NmbTet)
	{
		BeginDependency(LibParIdx, TetTyp, VerTyp);

		for(i=1;i<=msh->NmbTet;i++)
			AddDependencyFast(LibParIdx, 1, &i, 4, msh->tet[i].idx);

		EndDependency(LibParIdx, sta);
		printf(" tets      : %3.2f%% / %3.2f%%\n", sta[0],sta[1]);
	}

	if(msh->NmbPyr)
	{
		BeginDependency(LibParIdx, PyrTyp, VerTyp);

		for(i=1;i<=msh->NmbPyr;i++)
			AddDependencyFast(LibParIdx, 1, &i, 5, msh->pyr[i].idx);

		EndDependency(LibParIdx, sta);
		printf(" pyramids  : %3.2f%% / %3.2f%%\n",sta[0],sta[1]);
	}

	if(msh->NmbPri)
	{
		BeginDependency(LibParIdx, PriTyp, VerTyp);

		for(i=1;i<=msh->NmbPri;i++)
			AddDependencyFast(LibParIdx, 1, &i, 6, msh->pri[i].idx);

		EndDependency(LibParIdx, sta);
		printf(" prisms    : %3.2f%% / %3.2f%%\n",sta[0],sta[1]);
	}

	if(msh->NmbHex)
	{
		BeginDependency(LibParIdx, HexTyp, VerTyp);

		for(i=1;i<=msh->NmbHex;i++)
			AddDependencyFast(LibParIdx, 1, &i, 8, msh->hex[i].idx);

		EndDependency(LibParIdx, sta);
		printf(" hexes     : %3.2f%% / %3.2f%%\n",sta[0],sta[1]);
	}

	puts("");
}


/*----------------------------------------------------------*/
/* Parallel neighbours between tets							*/
/*----------------------------------------------------------*/

static void SetNgb(MshSct *msh, int LibParIdx)
{
	char *FlgTab;
	int i, NmbCpu, NmbTyp, HshSiz, MshSiz, (*NgbTab)[4];
	double timer;
	HshSct *tab;

	printf("Tet neighbours        : ");
	timer = 0.;
	GetTim(&timer);

	GetLplibInformation(LibParIdx, &NmbCpu, &NmbTyp);
	HshSiz =  (msh->NmbTet * 2)  / NmbCpu;
	MshSiz =  msh->NmbTet / NmbCpu;
	FlgTab = calloc( (msh->NmbTet+1), sizeof(char) );
	NgbTab = calloc( (msh->NmbTet+1), 4 * sizeof(int) );
	ParSct par[ NmbCpu ];

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

	LaunchParallel(LibParIdx, TetTyp, 0, (void *)ParNgb1, (void *)par);

	if(NmbCpu > 1)
		LaunchParallel(LibParIdx, TetTyp, 0, (void *)ParNgb2, (void *)par);

	free(tab);
	free(NgbTab);
	free(FlgTab);

	GetTim(&timer);
	printf("%g s\n", timer);
}


/*----------------------------------------------------------*/
/* The first scatter pass setups local links				*/
/*----------------------------------------------------------*/

static void ParNgb1(int BegIdx, int EndIdx, int c, ParSct *par)
{
	char *FlgTab = par[c].FlgTab;
	int i, j, k, key, (*NgbTab)[4] = par[c].NgbTab;
	unsigned int min, mid, max;
	TetSct *tet, *ngb;
	MshSct *msh = par[c].msh;
	HshSct *tab = par[c].tab = calloc(2 * par[c].HshSiz, sizeof(HshSct));

	for(i=par[c].beg; i<=par[c].end; i++)
	{
		tet = &msh->tet[i];
    
		for(j=0;j<4;j++)
		{
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
    
			if(!tab[ key ].tet)
			{
				tab[ key ].tet = i;
				tab[ key ].voy = j;
				tab[ key ].min = min;
				tab[ key ].mid = mid;
				tab[ key ].max = max;
				continue;
			}
    
			do
			{
				ngb = &msh->tet[ tab[ key ].tet ];
    
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


/*----------------------------------------------------------*/
/* The second gather pass setups missing globak links		*/
/*----------------------------------------------------------*/

static void ParNgb2(int BegIdx, int EndIdx, int c, ParSct *par)
{
	int n, i, j, k, key, BasKey, flg, (*NgbTab)[4] = par[c].NgbTab;
	unsigned int min, mid, max;
	TetSct *tet, *ngb;
	HshSct *tab;
	MshSct *msh = par[c].msh;

	for(i=par[c].beg; i<=par[c].end; i++)
	{
		if(par[c].FlgTab[i] == 4)
			continue;

		tet = &msh->tet[i];

		for(j=0;j<4;j++)
		{
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


/*----------------------------------------------------------*/
/* Vector safe tet renumbering								*/
/*----------------------------------------------------------*/

void VecRen(MshSct *msh, int VecSiz)
{
	int i, j, k, TetBlk = 0, BlkSiz, VecEnd = 0;
	TetSct *tet, *tet2;

	for(i=1;i<=msh->NmbVer;i++)
		msh->ver[i].cod = 0;

	for(i=1;i<=msh->NmbTet;i++)
		msh->tet[i].cod = 0;

	for(i=1;i<=msh->NmbTet;i++)
	{
		tet = &msh->tet[i];

		if(tet->cod)
			continue;

		tet->cod = ++TetBlk;
		BlkSiz = 1;

		for(k=0;k<4;k++)
			msh->ver[ tet->idx[k] ].cod = TetBlk;

		for(j=i+1;j<=msh->NmbTet;j++)
		{
			tet2 = &msh->tet[j];

			if(tet2->cod \
			|| (msh->ver[ tet2->idx[0] ].cod == TetBlk) \
			|| (msh->ver[ tet2->idx[1] ].cod == TetBlk) \
			|| (msh->ver[ tet2->idx[2] ].cod == TetBlk) \
			|| (msh->ver[ tet2->idx[3] ].cod == TetBlk) )
			{
				continue;
			}

			tet2->cod = TetBlk;

			for(k=0;k<4;k++)
				msh->ver[ tet2->idx[k] ].cod = TetBlk;

			BlkSiz++;

			if(BlkSiz == VecSiz)
				break;
		}

		if(!VecEnd && (BlkSiz < VecSiz))
			VecEnd = msh->NmbTet - i;
	}

	printf("Vector renumbering    : %d vec%d and %d last tets partially vectorized\n", TetBlk, VecSiz, VecEnd);
}
