/*************************************************************************
	> File Name: myfun.h
	> Author:kuangxiong 
	> Mail: kuangxiong@lsec.cc.ac.cn
	> Created Time: 2016年10月21日 星期五 21时26分17秒
 ************************************************************************/

/* f_h: \bar{f}^n */
/* mean value of f in [t_(n - 1), t_n] */
/* use Simpson formula */
#include"math.h"
#include"stdio.h"
#include"mpi.h"
#include"../fun.h"
#include"phg.h"
///*************me**************************************/
FLOAT
phgQuadDofDotGradBasBas(ELEMENT *e, DOF *u, DOF *v, int m, int n, int order, DOF *TB1, DOF *TB2, DOF *TB3)
{
	int i;
	const FLOAT *g1, *g2, *gB1, *gB2, *gB3, *w, *lambda;
	FLOAT d, d0;
//	DOF *B1, *B2, *B3;
	QUAD *quad;
    GRID *g=u->g;
//	assert(u->dim == 3);

//		order = 0;
	quad = phgQuadGetQuad3D(order);

    gB1= phgQuadGetDofValues(e, TB1, quad);
	gB2= phgQuadGetDofValues(e, TB2, quad);
	gB3= phgQuadGetDofValues(e, TB3, quad);
	g1 = phgQuadGetBasisGradient(e, v, m, quad);
	g2 = phgQuadGetBasisValues(e, v, n, quad);
	d = 0.;
    lambda = quad->points;
	w = quad->weights;
	for (i = 0; i < quad->npoints; i++) 
	{
		d0 = 0.;
		d0 += (*(gB1++)) * (*(g1++)) * (*(g2));	/* B1 dphi_m/dx phi_n */
		d0 += (*(gB2++)) * (*(g1++)) * (*(g2));	/* B2 dphi_m/dy phi_n */
		d0 += (*(gB3++)) * (*(g1++)) * (*(g2));	/* B3 dphi_m/dz phi_n */
		g2++;
		d += d0 * (*(w++));
		lambda += Dim + 1;
	}

	return d * phgGeomGetVolume(u->g, e);
}
///******************************************************/
///*********************computing D11*******************/
//FLOAT 
//GetD11(DOF *u_h)
//{
//  ELEMENT *e;
//	DOF *grad_u = phgDofGradient(u_h, NULL, NULL, NULL);
//	FLOAT norm;
//#if USE_MPI
//	FLOAT a, b;
//#endif
//	norm = 0.0;
//	ForAllElements(u_h->g, e)
//   {  
//      norm +=phgQuadDofDotDof(e, grad_u, grad_u, QUAD_DEFAULT);
//   }
//#if USE_MPI
//	a = norm;
//	MPI_Allreduce(&a, &b, 1, PHG_MPI_FLOAT, PHG_SUM, u_h->g->comm);
//	norm = b;
//#endif	
//	phgDofFree(&grad_u);
//	return norm;
//}
//
void
build_rhs_Vec(FLOAT timestp, VEC *g1, DOF *B1, MAP *map, DOF *u_h)
{
  int i, j;
  GRID *g = B1->g;
  ELEMENT *e;
  ForAllElements(g, e)
  {
	 int N= DofGetNBas(B1, e);
	 FLOAT A[N][N], B[N];
	 INT I[N];
	 for(i = 0;i< N; i++)
		 I[i] = phgMapE2L(map, 0, e, i);
	 for(i = 0; i< N; i++)
	{
		/* right hand side: \int f * phi_i i*/
		phgQuadDofTimesBas(e, B1, u_h, i, QUAD_DEFAULT, &B[i]);
		B[i] = timestp * B[i];
	}
    phgVecAddEntries(g1, 0, N, I, &B);
  }
}

void
build_rhs_Mat(MAT *C, MAP *map, DOF *u_h)
{
  int i, j;
  GRID *g = u_h->g;
  ELEMENT *e;
  ForAllElements(g, e)
  {
  int N= DofGetNBas(u_h, e);
  FLOAT A[N][N];
  INT I[N];
  for(i = 0;i< N; i++)
  {
     I[i] = phgMapE2L(map, 0, e, i);
	 for(j=0; j< N; j++)
		 A[i][j] =  phgQuadBasDotBas(e, u_h, j, u_h, i, QUAD_DEFAULT);
  }
	 for(i = 0; i< N; i++)
		phgMatAddEntries(C, 1, I+i, N, I, A[i]); 
  }
}

void
build_stiff_Mat(FLOAT timestp, FLOAT a,  MAT *G1, MAP *map, DOF *u_h)
{  int i, j;
	GRID *g = u_h->g;
    ELEMENT *e;
   ForAllElements(g, e)
  {
	 int N= DofGetNBas(u_h, e);
	FLOAT A[N][N], rhs[N], B[N];
	FLOAT tmp;
	INT I[N];
	for(i = 0;i< N; i++)
	 {
		  I[i] = phgMapE2L(map, 0, e, i);
		 for(j=0; j< N; j++)
			 A[i][j] =
				a*phgQuadGradBasDotGradBas(e, u_h, j, u_h, i, QUAD_DEFAULT) * timestp
				+ phgQuadBasDotBas(e, u_h, j, u_h, i, QUAD_DEFAULT);
	 }
	 for(i = 0; i< N; i++)
		phgMatAddEntries(G1, 1, I+i, N, I, A[i]); 
  }
}
void
build_stiff_Mat1(FLOAT timestp, FLOAT a,  MAT *G1, MAP *map, DOF *B1, DOF *B2, DOF *B3, DOF *u_h)
{  int i, j;
	GRID *g = u_h->g;
    ELEMENT *e;
   phgPrintf("D0%f\n", D0);
   ForAllElements(g, e)
  {
	 int N= DofGetNBas(u_h, e);
	FLOAT A[N][N], rhs[N], B[N];
	FLOAT tmp;
	INT I[N];
	for(i = 0;i< N; i++)
	 {
		  I[i] = phgMapE2L(map, 0, e, i);
		 for(j=0; j< N; j++)
			 A[i][j] =
				a*phgQuadGradBasDotGradBas(e, u_h, j, u_h, i, QUAD_DEFAULT) * timestp
				+ phgQuadBasDotBas(e, u_h, j, u_h, i, QUAD_DEFAULT)
				+ phgQuadDofDotGradBasBas(e, u_h, u_h, j, i , 3, B1, B2, B3) * timestp;		
	 }
	 for(i = 0; i< N; i++)
		phgMatAddEntries(G1, 1, I+i, N, I, A[i]); 
  }
}
void
build_stiff_Mat2(FLOAT timestp, MAT *G1, MAP *map, DOF *B1, DOF *B2, DOF *B3, DOF *u_h)
{  int i, j;
	GRID *g = u_h->g;
    ELEMENT *e;

   ForAllElements(g, e)
  {
	 int N= DofGetNBas(u_h, e);
	FLOAT A[N][N], rhs[N], B[N];
	FLOAT tmp;
	INT I[N];
	for(i = 0;i< N; i++)
	 {
		  I[i] = phgMapE2L(map, 0, e, i);
		 for(j=0; j< N; j++)
			 A[i][j] =  phgQuadDofDotGradBasBas(e, u_h, u_h, j, i , 3, B1, B2, B3) * timestp;		
	 }
	 for(i = 0; i< N; i++)
		phgMatAddEntries(G1, 1, I+i, N, I, A[i]); 
  }
}

///* space error indicator*/
//FLOAT
//estimate_space_error(DOF *u_h, DOF *u_p, DOF *f_h, DOF *Btmp1, DOF *Btmp2, DOF *Btmp3, DOF *error)
//{
//    GRID *g = u_h->g;
//    ELEMENT *e;
//    DOF *tmp, *jump;
//    FLOAT eta, h, d;
//
//	int i;
//  
//    tmp = phgDofGradient(u_h, NULL, NULL, "tmp");
//    jump = phgQuadFaceJump(tmp, DOF_PROJ_DOT, NULL, QUAD_DEFAULT);
//    phgDofFree(&tmp);
//
//    eta = 1.0 / stime;
//    tmp = phgDofCopy(f_h, NULL, f_h->type, NULL);
//    phgDofAXPY(eta, u_p, &tmp);
//    eta *= -1.0;
//    phgDofAXPY(eta, u_h, &tmp);
//
//    ForAllElements(g, e){
//        eta = 0.0;
//        for (i = 0; i < NFace; i++) {
//            if (e->bound_type[i] & (DIRICHLET | NEUMANN))
//                continue;    /* boundary face */
//                h = phgGeomGetFaceDiameter(g, e, i);
//                eta +=  (*DofFaceData(jump, e->faces[i])) * h;
//        }
//        h = phgGeomGetDiameter(g, e);
//        eta += eta * 0.5  + h * h * phgQuadDofDotDof(e, tmp, tmp, QUAD_DEFAULT);
//       // eta += eta * 0.5  + h * h * phgQuadDofDotDof(e, tmp, tmp, QUAD_DEFAULT) ;
//
//        *DofElementData(error, e->index) = eta;
//    }
//    phgDofFree(&tmp);
//    phgDofFree(&jump);
//    return phgDofNormInftyVec(error);
//}
//
FLOAT
phgVecDotDofGrad(ELEMENT *e,DOF *B1, DOF *B2, DOF *B3, DOF *u_h)
{
	GRID *g = u_h->g;
	int i, j;
	QUAD *quad;
    FLOAT *gB1, *gB2, *gB3, d, d0, *g1, *w, *lambda;
    int  N= DofGetNBas(u_h, e);
	FLOAT data[N];
   
	d = 0.0;
	quad = phgQuadGetQuad3D(3);
    lambda = quad->points;
	w = quad->weights;
	phgDofGetElementDatas(u_h, e, data);
	for(i=0 ; i< N; i++)
	{
		gB1 = phgQuadGetDofValues(e, B1, quad);
		gB2 = phgQuadGetDofValues(e, B2, quad);
		gB3 = phgQuadGetDofValues(e, B3, quad);
		g1 =  phgQuadGetBasisGradient(e, u_h, i, quad);   
		for (j = 0; j < quad->npoints; j++) 
		{
			d0 = 0.;
			d0 += data[i]*(*(gB1++)) * (*(g1++)) ;	/* B1 dphi_m/dx phi_n */
			d0 += data[i]*(*(gB2++)) * (*(g1++)) ;	/* B2 dphi_m/dy phi_n */
			d0 += data[i]*(*(gB3++)) * (*(g1++)) ;	/* B3 dphi_m/dz phi_n */
			d += d0 * (*(w++));
			lambda += Dim + 1;
		}
	}

	return d * phgGeomGetVolume(g, e);
}

void
phgMatSvd(FLOAT *data, MAP *map, INT nT, FLOAT *U, FLOAT *S, INT  numprocs, int rank)
{  
	int context, M , N, nb, nprow, npcol, mVT, nVT, info, ONE, ZERO;
	int descU[9], descVT[9], descA[9], i, j, k, myrow, mycol, maxnlocal, lwork;
    FLOAT *work, wkopt, *adjustA, *VT;
	ONE = 1;
	ZERO = 0;
	MPI_Reduce(&map->nlocal, &maxnlocal, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
    
	MPI_Bcast(&maxnlocal, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
	M = numprocs * maxnlocal;

	N = nT;
    nb = 1;
	nprow = numprocs;
	npcol =1;
	Cblacs_pinfo(&rank, &numprocs);
	Cblacs_get(-1, 0, &context);
	Cblacs_gridinit(&context, "Row", nprow, npcol);
	Cblacs_gridinfo(context, &nprow, &npcol, &myrow, &mycol);	
	adjustA = (double *)malloc(maxnlocal *N * sizeof(double));
	for(i=0; i< maxnlocal; i++)
		for(j=0; j<N; j++)
		{	if(i< map->nlocal)
			adjustA[i+j*maxnlocal] = data[i+j*map->nlocal]; 
            if(i>= map->nlocal)
			adjustA[i+j*maxnlocal] = 0;
		}
/******************describe  matrix*****************************/
	descA[0]=1; descA[1]= 0; descA[2]= M; descA[3] = N; descA[4] =1;
	descA[5]=1; descA[6]= 0; descA[7]= 0;descA[8] = maxnlocal;

	descU[0]=1; descU[1]= 0; descU[2]= M; descU[3] = M; descU[4] =1;
	descU[5]=1; descU[6]= 0; descU[7]= 0;  descU[8] = maxnlocal;
    phgPrintf("/***********start SVD decomposition**********************/\n");
	lwork= -1;
	info = 0;
	pdgesvd_("V", "N", &M, &N, adjustA, &ONE, &ONE, descA, S, U, &ONE, &ONE, &descU, NULL, NULL, NULL, NULL, &wkopt, &lwork, &info);
	lwork=(int)wkopt;
	work = (double*)malloc(lwork*sizeof(double));
	pdgesvd_("V", "N", &M, &N, adjustA, &ONE, &ONE, descA, S, U, &ONE, &ONE, &descU, NULL, NULL, NULL, NULL, work, &lwork, &info);

}
int
phgGetPodNumbasis(FLOAT rate, FLOAT *S, int length)
{
	double globalsum, partsum;
	int i, j, n;
	if(rate > 1)
	{
	   phgPrintf("\n error: rate is biger than 1\n");
	   return 0;
	}
	globalsum = 0 ;
	partsum = 0;
	for(i=0 ; i< length; i++)
		globalsum +=S[i];
	for(i=0 ;i< length; i++)
	{
	   partsum = partsum + S[i];
	   if(partsum >= rate * globalsum)
	   {
	      n = i;
		  break;
	   }
	}
    return n+1;
}

//void
//phgFEMLoopTime(FLOAT ctime, DOF *u_h, DOF *u_p, MAT *C, MAT *G1, MAT *G2, VEC *g1, VEC *g2)
//{     
//	MAP *map = g1->map;
//	VEC *btmp, *btmp1, *b;
//	MAT *A;
//    SOLVER *solver;
//    
//    btmp = phgMapCreateVec(map, 1);
//    btmp1 = phgMapCreateVec(map, 1);
//
//    phgMapDofArraysToVec(map, 1, FALSE, &btmp, &u_p, NULL);
//    b = phgMapCreateVec(map, 1);
//	phgVecDisassemble(b);
//	phgMatVec(0, 1.0, C, btmp, 0.0, &b);
//		
//	phgVecAXPBY(-1.0, g1, 0.0 , &btmp1); 
//    phgVecAXPBY(-1.0*f_tmp(ctime), g2, 1.0, &btmp1);
//	phgVecAXPBY(1.0, btmp1, 1.0, &b);
//
//	A = phgMapCreateMat(map,map);
//	phgMatAXPBY(1.0, G1, 0, &A);
//	phgMatAXPBY(cos(ctime), G2, 1.0, &A);
//		
//    solver = phgSolverCreate(SOLVER_GMRES, u_h, NULL);
//	solver->mat = A;
//	solver->rhs = b;
//    phgSolverSolve(solver, TRUE, u_h, NULL);
//	phgSolverDestroy(&solver);
//}
//
void
phgPODDofToFEMDof(MAP *map, INT N0, FLOAT *U, INT maxnlocal, FLOAT *PODb, DOF *u_h)
{
	VEC *Vtf, *Vtmp;
    int i, j;	
	Vtmp = phgMapCreateVec(map, 1);
	Vtf = phgMapCreateVec(map, 1);
	for(j=0; j< N0; j++)
	{
	    memcpy(Vtmp->data, &U[j*maxnlocal], map->nlocal*sizeof(double));
        phgVecAXPBY(PODb[j], Vtmp, 1.0, &Vtf);
	}
	
	phgMapVecToDofArrays(map, Vtf, FALSE, &u_h,  NULL);
 }

/*****************have some problem*************************/
static void 
phgPODLoopTime(FLOAT ctime, INT N0, FLOAT *PODu_p, FLOAT *PODG1, FLOAT *PODG2, FLOAT *PODg1, FLOAT *PODg2, FLOAT *PODC, FLOAT *PODb)
{   int i, j, k, dzero, info, IPIV[N0], ONE;
	FLOAT PODmat[N0][N0], y[N0];

	ONE =1;
	dzero = 0.0;

	for(i=0; i<N0; i++)
         PODb[i] = -1.0*PODg1[i] - cos(ctime)* PODg2[i];

	for(i=0; i<N0; i++)
		for(j=0 ; j< N0; j++)
			PODmat[i][j] = PODG1[i+j*N0] + cos(ctime) * PODG2[i+j*N0];
//	for(i=0; i< N0; i++)
//		y[i] =0;
	dgemv_("T", &N0, &N0, &ONE, PODC, &N0, PODu_p, &ONE, &dzero, y, &ONE);
    for(i=0; i< N0; i++)
		PODb[i] += y[i];
     
    dgesv_(&N0, &ONE, PODmat, &N0, IPIV, PODb, &N0, &info);
}
