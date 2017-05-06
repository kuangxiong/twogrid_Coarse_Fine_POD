/*************************************************************************
	> File Name: APODfun.h
	> Author:kuangxiong 
	> Mail: kuangxiong@lsec.cc.ac.cn
	> Created Time: 2016年10月21日 星期五 21时26分17秒
 ************************************************************************/
#ifndef  __APODFUN_H_
#define __APODFUN_H_
#endif

#include<stdio.h>
#include"mpi.h"
#include"math.h"
#include"phg.h"

static FLOAT APODT0 = 4.0;
static INT CLT = 3; /* number of coarse grid dof for updata pod basis  */
static FLOAT T = 50;
static INT RefineN = 12;
static INT CosRefineN = 8;
static FLOAT threshold = 0.00005;
static FLOAT gamma1 = 0.9999;
static FLOAT gamma2 = 0.99999;
static FLOAT gamma3 = 0.999999;
static FLOAT stime1 = 0.05;
static FLOAT stime = 0.01;
static FLOAT D0 = 0.5;
static FLOAT crtime = 0.0;    /* current time. Modified by program */
static FLOAT ptime = 0.0;     /* previous time. */

static const char *fn="meshdata/testcube6.dat";

static void     /* rhs: f */
func_f(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
//   *value = -Cos(y)-Sin(y)*Cos(crtime);
//   *value = -Cos(y)-Sin(y)*Cos(crtime) + 2*exp(-(crtime-100)*(crtime -100)/0.05);
     *value = (1 + 3*crtime)*Cos(x)*Cos(y)*Cos(z);
//	   *value = y + cos(0.1*crtime) + 10*exp(-(crtime -50)*(crtime -50)/0.05);
//     *value = y*Cos(0.5*crtime) + 5*exp(-(crtime-200)*(crtime-200)/0.05);
//     *value = -2*sin(2*crtime)*(x+y+z)*0.5*exp(4*Cos(2*crtime));
//     *value = 0.1*crtime*Cos(crtime)*(x + y + z);
//     *value = y*0.1*crtime*Cos(x);
//       *value = y*Cos(0.5*crtime) + 5*exp(-(crtime-100)*(crtime-100)/0.05);
//     *value = Sin(crtime)*Sin(y) -Cos(y) - Sin(y)*Cos(crtime);
}
static void     /* rhs: f */
func_u(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
//   *value = -Cos(y)-Sin(y)*Cos(crtime);
//     *value = crtime*Cos(x)*Cos(y)*Cos(z);
}

static void 
func_B1(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
	*value = Cos(x);
//    *value = 0;
}

static void
func_B2(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
	*value = Cos(y) ;
//    *value = 0;
}

static void 
func_B3(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
	*value = Cos(z);
//    *value = 0;
}

static void     /* rhs: f */
func_C1(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
//  *value = Sin(x);
    *value = 0;
}
static void     /* rhs: f */
func_C2(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
//   *value = Sin(y);
    *value = 0;
}

static void     /* rhs: f */
func_C3(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
//   *value = Sin(z);
    *value = 0;
}

static void
func_u0(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
	*value = 0;
//   *value = -Cos(y)- Sin(y);
}
static FLOAT
eta_t(FLOAT t)
{   FLOAT z;
	z = cos(t);
	return z;
}

static FLOAT 
f_tmp(FLOAT t)
{  FLOAT z;
//  if(t <=500)
		z = cos(t);
//	else 
//		z = cos(t) + sin(t);
    return z;
}

///*********************computing D11*******************/
static FLOAT 
GetD11(DOF *u_h)
{
     ELEMENT *e;
	 DOF *grad_u = phgDofGradient(u_h, NULL, NULL, NULL);

	 FLOAT  norm;
#if USE_MPI
	 FLOAT a, b;
#endif

	 norm = 0;
     ForAllElements(u_h->g, e)
	 {
        norm += phgQuadDofDotDof(e, grad_u, grad_u, QUAD_DEFAULT);
	 }
#if USE_MPI
	 a = norm;
	 MPI_Allreduce(&a, &b, 1, PHG_MPI_FLOAT, PHG_SUM, u_h->g->comm);
	 norm = b;
#endif
	 phgDofFree(&grad_u);
	 
	 return norm;
}

static FLOAT
estimate_space_error3(FLOAT time, MAT *G1, MAT *G2, MAT *C, VEC *b, DOF *u_h, DOF *u_p)
{
    MAT *M;
	VEC *bb, *btmp,*btmp1,  *b1;
	MAP *map = b->map;
	GRID *g = u_h->g;
    ELEMENT *e;
    FLOAT eta, h, d;

	bb = phgMapCreateVec(map, 1);
	b1 = phgMapCreateVec(map, 1);
	btmp = phgMapCreateVec(map, 1);
	btmp1 = phgMapCreateVec(map, 1);
	M = phgMapCreateMat(map, map);
	phgMapDofArraysToVec(map, 1, FALSE, &btmp, &u_p, NULL);

	phgMatAXPBY(1.0, G1, 0.0, &M);
	phgMatAXPBY(1.0*Cos(time), G2, 1.0, &M);
    phgMatVec(0, 1.0, C, btmp, 0.0, &bb);
	phgVecAXPBY(1.0, b, 1.0, &bb);
    
	phgMapDofArraysToVec(map, 1, FALSE, &btmp1, &u_h, NULL);
	
	phgMatVec(0, 1.0, M, btmp1, 0.0, &b1);
	phgVecAXPBY(-1.0, b1, 1.0, &bb);

    phgMatDestroy(&M);

	//return phgVecNorm2(b, 0, NULL)/phgVecNorm2(btmp1, 0, NULL);
    return phgVecNorm2(bb, 0, NULL);
}
static FLOAT
estimate_space_error(FLOAT time, MAT *G1, MAT *G2, MAT *C, VEC *g1, VEC *g2, DOF *u_h, DOF *u_p)
{
    MAT *M;
	VEC *bb, *btmp,*btmp1,  *b1;
	MAP *map = g1->map;
	GRID *g = u_h->g;
    ELEMENT *e;
    FLOAT eta, h, d;

	bb = phgMapCreateVec(map, 1);
	b1 = phgMapCreateVec(map, 1);
	btmp = phgMapCreateVec(map, 1);
	btmp1 = phgMapCreateVec(map, 1);
	M = phgMapCreateMat(map, map);
	phgMapDofArraysToVec(map, 1, FALSE, &btmp, &u_p, NULL);

	phgMatAXPBY(1.0, G1, 0.0, &M);
	phgMatAXPBY(1.0*Cos(time), G2, 1.0, &M);
    phgMatVec(0, 1.0, C, btmp, 0.0, &bb);
	phgVecAXPBY(-1.0, g1, 1.0, &bb);
	phgVecAXPBY(-1.0*f_tmp(time), g2, 1.0, &bb);
    
	phgMapDofArraysToVec(map, 1, FALSE, &btmp1, &u_h, NULL);
	
	phgMatVec(0, 1.0, M, btmp1, 0.0, &b1);
	phgVecAXPBY(-1.0, b1, 1.0, &bb);

    phgMatDestroy(&M);

	//return phgVecNorm2(b, 0, NULL)/phgVecNorm2(btmp1, 0, NULL);
    return phgVecNorm2(bb, 0, NULL);
}
static FLOAT
Newestimate_space_error(DOF *f_h, DOF *u_h, DOF *u_p)
{
    DOF *error;
	FLOAT a;
	GRID *g = u_h->g;
    error = phgDofNew(g, DOF_P1, 1, "error", DofInterpolation);
    phgDofAXPBY(1.0, u_h, 0, &error);
    phgDofAXPBY(-1.0, u_p, 1.0, &error);
    phgDofAXPBY(1.0, f_h, -1.0/stime, &error);
    a = phgDofNormL2(error);
	phgDofFree(&error);
    return a;
}

static FLOAT
estimate_space_error1(DOF *u_h, DOF *u_p, DOF *f_h, DOF *error)
{
    GRID *g = u_h->g;
    ELEMENT *e;
    DOF *tmp, *jump;
    FLOAT eta, h;
    int i;
  
    tmp = phgDofGradient(u_h, NULL, NULL, "tmp");
    jump = phgQuadFaceJump(tmp, DOF_PROJ_DOT, NULL, QUAD_DEFAULT);
    phgDofFree(&tmp);

    eta = 1.0 / stime;
    tmp = phgDofCopy(f_h, NULL, f_h->type, NULL);
    phgDofAXPY(eta, u_p, &tmp);
    eta *= -1.0;
    phgDofAXPY(eta, u_h, &tmp);

    ForAllElements(g, e){
        eta = 0.0;
        for (i = 0; i < NFace; i++) {
            if (e->bound_type[i] & (DIRICHLET | NEUMANN))
                continue;    /* boundary face */
                h = phgGeomGetFaceDiameter(g, e, i);
                eta +=  (*DofFaceData(jump, e->faces[i])) * h;
        }
        h = phgGeomGetDiameter(g, e);
    //    eta += eta * 0.5  + h * h * phgQuadDofDotDof(e, tmp, tmp, QUAD_DEFAULT);
        eta =  h * h * phgQuadDofDotDof(e, tmp, tmp, QUAD_DEFAULT);

        *DofElementData(error, e->index) = eta;
    }
    phgDofFree(&tmp);
    phgDofFree(&jump);
    return phgDofNormInftyVec(error);
}

static FLOAT
estimate_time_error(DOF *u_h, DOF *u_p)
{
    DOF *tmp, *grad;
    FLOAT ss;

    tmp = phgDofCopy(u_h, NULL, DOF_DEFAULT, "tmp2");
    phgDofAXPY( -1.0, u_p, &tmp);
    grad = phgDofGradient(tmp, NULL, NULL, "grad");
    ss = phgDofNormL2(grad);

    phgDofFree(&tmp);
    phgDofFree(&grad);

    return ss * ss * (FLOAT)(1.L/3.L);
}

/* data oscillation */

static void
phgFEMLoopTime3(FLOAT ctime, DOF *u_h, DOF *u_p, MAT *C, MAT *G1, MAT *G2, VEC *bb)
{     
	MAP *map = bb->map;
	VEC *btmp, *b;
	MAT *A;
    SOLVER *solver;
    
    btmp = phgMapCreateVec(map, 1);

    phgMapDofArraysToVec(map, 1, FALSE, &btmp, &u_p, NULL);

    b = phgMapCreateVec(map, 1);
	phgVecDisassemble(b);
	phgMatVec(0, 1.0, C, btmp, 0.0, &b);
	phgVecAXPBY(1.0, bb, 1.0 , &b); 

	A = phgMapCreateMat(map,map);
	phgMatAXPBY(1.0, G1, 0, &A);
	phgMatAXPBY(eta_t(ctime), G2, 1.0, &A);
		
    solver = phgSolverCreate(SOLVER_GMRES, u_h, NULL);
	solver->mat = A;
	solver->rhs = b;
    phgSolverSolve(solver, TRUE, u_h, NULL);
	phgSolverDestroy(&solver);

}
static void
phgRefineFEMLoopTime1(GRID *g, INT TotalN,  DOF **uhs, INT rank, char *snD1)
{
    INT i, j, k, flag = 0;
	FLOAT err;
	MAT *AA, *GG1, *GG2, *CC;
	VEC *bb, *bbtmp, *bbtmp1;
	MAP *mmap;
	SOLVER *ssolver;
	ELEMENT *ee;
	DOF *uu_h, **uuhs, *ff_h;  /* numerical solution at time t_n */
	DOF *uu_p, *ggrad_u;  /* numerical solution at time t_{n-1} */
	DOF  *BB1, *BB2, *BB3, *CC1, *CC2, *CC3, *eerror;
    
	uuhs = phgAlloc(TotalN * sizeof(*uhs));
	uu_p = phgDofNew(g, DOF_P1, 1, "uu_p", DofInterpolation);
	ff_h = phgDofNew(g, DOF_P1, 1, "ff_h", DofInterpolation);
	phgDofSetDataByFunction(uu_p, func_u0);
	uu_h = phgDofNew(g, DOF_P1, 1, "uu_h", DofInterpolation);
	eerror = phgDofNew(g, DOF_P1, 1, "eerror", DofInterpolation);	
	uuhs[0] = phgDofCopy(uu_p, NULL, DOF_P1, "uu_p1");
	uuhs[0] ->userfunc = DofInterpolation;
	k = 1;
	BB1=phgDofNew(g, DOF_P1, 1, "BB1", func_B1);
	BB2=phgDofNew(g, DOF_P1, 1, "BB2", func_B2);
	BB3=phgDofNew(g, DOF_P1, 1, "BB3", func_B3);
	CC1=phgDofNew(g, DOF_P1, 1, "CC1", func_C1);
	CC2=phgDofNew(g, DOF_P1, 1, "CC2", func_C2);
	CC3=phgDofNew(g, DOF_P1, 1, "CC3", func_C3);

	mmap = phgMapCreate(uu_h, NULL);
	GG1 = phgMapCreateMat(mmap, mmap);
	GG2 = phgMapCreateMat(mmap, mmap);
	CC = phgMapCreateMat(mmap, mmap);

	bbtmp = phgMapCreateVec(mmap, 1);
	phgVecDisassemble(bbtmp);

	build_rhs_Mat(CC, mmap, uu_h);
	build_stiff_Mat1(D0, GG1, mmap, BB1, BB2, BB3, uu_h);
	build_stiff_Mat2(GG2, mmap, CC1, CC2, CC3, uu_h);
	crtime = 0;
	while(crtime < T - 1e-8)
	{
		/********************************************************************/	
		ptime = crtime;
		crtime += stime;
		phgPrintf("\n/********* start new time layer *************/\n");
		phgPrintf("current time layer: [%lf, %lf]\n", (double)ptime, (double)crtime);
		flag++;
		if (flag > 1)
		{   /* update u_p */
			phgDofFree(&uu_p);
			uu_p = phgDofCopy(uu_h, NULL, DOF_P1, "uu_p");
			uu_p->userfunc = DofInterpolation;
		}
		phgMapDofArraysToVec(mmap, 1, FALSE, &bbtmp, &uu_p, NULL);
		bb = phgMapCreateVec(mmap, 1);
		bbtmp1 = phgMapCreateVec(mmap, 1);
		phgVecDisassemble(bb);
		phgVecDisassemble(bbtmp1);

		phgMatVec(0, 1.0, CC, bbtmp, 0.0, &bb);	
		phgDofSetDataByFunction(ff_h, func_f);

		build_rhs_Vec(bbtmp1, ff_h, mmap, uu_h);
		phgVecAXPBY(1.0, bbtmp1, 1.0, &bb);

		AA = phgMapCreateMat(mmap, mmap);
		phgMatAXPBY(1.0, GG1, 0.0, &AA);
		phgMatAXPBY(eta_t(crtime), GG2, 1.0, &AA);

		ssolver = phgSolverCreate(SOLVER_GMRES, uu_h, NULL);
		ssolver->mat = AA;
		ssolver->rhs = bb;
		phgSolverSolve(ssolver, TRUE, uu_h, NULL);
		phgSolverDestroy(&ssolver);

		uuhs[k] = phgDofCopy(uu_h, NULL, DOF_P1, "uu_p1");
		uuhs[k] ->userfunc = DofInterpolation;
		k++;
		phgPrintf("\ntime step: %lg\n\n", (double)stime);
	}

	for(i=1; i< TotalN; i++)
	{    
		phgDofAXPBY(1.0, uhs[i], 0, &eerror);
		phgDofAXPBY(-1.0, uuhs[i], 1.0, &eerror);
		err = phgDofNormL2(eerror);//phgDofNormL2(uhs[i]);
		if(rank==0)
	    	fprintf(snD1, "%f\t%e\n", i*stime, err);
	}
	for(i=0 ;i< TotalN; i++)
		phgDofFree(uuhs+i);
	phgFree(uuhs);
	phgDofFree(&eerror);
	phgDofFree(&BB1);
	phgDofFree(&BB2);
	phgDofFree(&BB3);
	phgDofFree(&ff_h);

	phgDofFree(&CC1);
	phgDofFree(&CC2);
	phgDofFree(&CC3);
	phgDofFree(&uu_h);
	phgDofFree(&uu_p);
}

static int 
gettimeflag(GRID *Cosg, int *timeflag, FLOAT timestp, int TotalN1, int numprocs, INT rank)
{
	DOF **Cosuhs, **CosPODuhs, *Cosu_h, *Cosu_p, *error, *FEMu_p, *FEMu_h; 
    DOF *CosB1, *CosB2, *CosB3, *CosC1, *CosC2, *CosC3, *Cosf_h; 
	VEC *Cosb, *CosV, *Cosb1, *btmp, *CosV1; 
	MAT *CosG1, *CosG2, *CosC, *A;
	MAP *Cosmap;
    SOLVER *solver;
	int i, k, CosnCLT, CosM, CosN, j, maxnlocal, nT1, N0, N1, GetPOD =1;
	int ONE = 1, info, flag =0, FEMflag = 0;
	FLOAT *U, *S, *U1, *S1, *svdM;

   FEMu_p = phgDofNew(Cosg, DOF_P1, 1, "FEMu_p", DofInterpolation);
    FEMu_h = phgDofNew(Cosg, DOF_P1, 1, "FEMu_h", DofInterpolation);
	Cosu_p = phgDofNew(Cosg, DOF_P1, 1, "Cosu_p", DofInterpolation);
	phgDofSetDataByFunction(Cosu_p, func_u0);
    Cosf_h = phgDofNew(Cosg, DOF_P1, 1, "Cosf_h", DofInterpolation);
    error = phgDofNew(Cosg, DOF_P1, 1, "error", DofInterpolation);

	CosB1 = phgDofNew(Cosg, DOF_P1, 1, "B1", func_B1);
	CosB2 = phgDofNew(Cosg, DOF_P1, 1, "B2", func_B2);
	CosB3 = phgDofNew(Cosg, DOF_P1, 1, "B3", func_B3);
	CosC1 = phgDofNew(Cosg, DOF_P1, 1, "C1", func_C1);
	CosC2 = phgDofNew(Cosg, DOF_P1, 1, "C2", func_C2);
	CosC3 = phgDofNew(Cosg, DOF_P1, 1, "C3", func_C3);
    Cosu_h = phgDofNew(Cosg, DOF_P1, 1, "u_h", DofInterpolation);

	Cosmap = phgMapCreate(Cosu_h, NULL);

	CosG1 = phgMapCreateMat(Cosmap, Cosmap);
	CosG2 = phgMapCreateMat(Cosmap, Cosmap);
	CosC = phgMapCreateMat(Cosmap, Cosmap);

    btmp = phgMapCreateVec(Cosmap, 1);
	phgVecDisassemble(btmp);

	build_rhs_Mat(timestp, CosC, Cosmap, Cosu_h);
    build_stiff_Mat1(timestp, D0, CosG1, Cosmap, CosB1, CosB2, CosB3, Cosu_h);
    build_stiff_Mat2(timestp, CosG2, Cosmap, CosC1, CosC2, CosC3, Cosu_h);
	
	CosnCLT = CLT /timestp + 1;
	nT1 = APODT0 /timestp + 1;

	Cosuhs = phgAlloc(TotalN1* sizeof(*Cosuhs));

	Cosuhs[0] = phgDofCopy(Cosu_p, NULL, DOF_P1, NULL);
    Cosuhs[0]->userfunc = DofInterpolation;


    k = 0;   crtime = 0;
	while(crtime < APODT0 - 1e-8)
	{	
	/*********************************************************/
		crtime = crtime + timestp;
		phgPrintf("\n/********* start new time layer(Coarse FEM) *************/\n");
        phgPrintf("current time layer: [%lf, %lf]\n", (FLOAT)crtime - timestp, (FLOAT)crtime);
		Cosb = phgMapCreateVec(Cosmap, 1);
		phgVecDisassemble(Cosb);

		phgDofSetDataByFunction(Cosf_h, func_f);
		build_rhs_Vec(timestp, Cosb, Cosf_h, Cosmap, Cosu_h);
		phgFEMLoopTime3(crtime, Cosu_h, Cosuhs[k], CosC, CosG1, CosG2, Cosb);

        k++;
		Cosuhs[k] = phgDofCopy(Cosu_h, NULL, DOF_P1, NULL);
        Cosuhs[k]->userfunc = DofInterpolation;
	    phgDofFree(&FEMu_p);
		FEMu_p = phgDofCopy(Cosuhs[k], NULL, DOF_P1, NULL);
	}
/******************** POD *************************************************/
	CosV = phgMapCreateVec(Cosmap, nT1);
	phgMapDofArraysToVec(Cosmap, nT1, FALSE, &CosV, Cosuhs, NULL);
	CosV1 = phgMapCreateVec(Cosmap, CosnCLT);
/******************************SVD decomposition***********************/
	MPI_Reduce(&Cosmap->nlocal, &maxnlocal, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Bcast(&maxnlocal, 1, MPI_INT, 0, MPI_COMM_WORLD);
	CosM = numprocs * maxnlocal;
	CosN = nT1;

    U = (FLOAT *)malloc((maxnlocal * CosM) * sizeof(FLOAT));
    U1 = (FLOAT *)malloc((maxnlocal * CosM) * sizeof(FLOAT));
    S = (FLOAT *)malloc(CosN * sizeof(FLOAT));
    phgMatSvd(CosV->data, Cosmap, nT1, U, S, numprocs, rank);
	N0 = phgGetPodNumbasis(gamma1, S, CosnCLT);
//	for(i=0; i< CosnCLT; i++)
//		phgPrintf(" %f\n", S[i]);
	free(S);
/*********************************************************************/
    FLOAT *PODmat, *PODC, *PODu_p;
    VEC *PODu0, *VL, *VR, *Vtmp, *Vtf;
	FLOAT *PODG1, *PODG2, *PODg1, *PODg2, *PODb, *y;
    int  PODflag = 0, *IPIV;
	FLOAT  done, dzero, etaspace;
	
	done = 1.0;  dzero = 0.0;
    
	Vtmp = phgMapCreateVec(Cosmap, 1);
	PODu0 = phgMapCreateVec(Cosmap, 1);
    VL = phgMapCreateVec(Cosmap, 1);
	VR = phgMapCreateVec(Cosmap, 1);
/********************start  loop time****************/
    PODflag = 0;
//	crtime = APODT0;
    while(crtime < T - 1e-8)
   {
    if(phgDofNormL2(error) > threshold && PODflag > 0)
	  {	
		timeflag[k] = 1;
		k--;
		phgDofFree(Cosuhs + k + 1);
		crtime = crtime - timestp;
		ptime = crtime - stime;
		i = 0;
		/****************updata POD basis***************************/
		while(i < CosnCLT)
        {
		   if(i!= 0)	
		   {
			   phgDofFree(&FEMu_p);
		       FEMu_p = phgDofCopy(FEMu_h, NULL, DOF_P1, NULL);
		       FEMu_p->userfunc = DofInterpolation;
		   }

		   ptime = crtime;
		   crtime += timestp;
		   phgPrintf("\n/********* start new time layer *************/\n");
           phgPrintf("current time layer: [%lf, %lf]\n", (FLOAT)ptime, (FLOAT)crtime);

           Cosb = phgMapCreateVec(Cosmap, 1);
	       phgVecDisassemble(Cosb);
	       phgDofSetDataByFunction(Cosf_h, func_f);
		   build_rhs_Vec(timestp, Cosb, Cosf_h, Cosmap, Cosu_h);

		   phgFEMLoopTime3(crtime, FEMu_h, FEMu_p, CosC, CosG1, CosG2, Cosb);

	       phgFEMLoopTime3(crtime, Cosu_h, Cosuhs[k], CosC, CosG1, CosG2, Cosb);

		   k++;
		   Cosuhs[k] = phgDofCopy(Cosu_h, NULL, DOF_P1, "u_p1");
		   Cosuhs[k]->userfunc = DofInterpolation;
		   i++;
		   if(crtime > T - 1e-8)
		   {	  
			  phgPrintf("\n end of FEM :%f\n", crtime);
              return 0;		      
		   }
		}
		//**********************SVD1*************************/
		phgMapDofArraysToVec(Cosmap, CosnCLT, FALSE, &CosV1, &Cosuhs[k-CosnCLT+1], NULL);
		S1 = (FLOAT *)malloc(CosnCLT*sizeof(FLOAT));
        phgMatSvd(CosV1->data, Cosmap, CosnCLT, U1, S1, numprocs, rank);
		N1 = phgGetPodNumbasis(gamma2, S1, CosnCLT);
		phgPrintf("N1:%d\n", N1);
		/********************SVD2*************************/
		svdM =(FLOAT *) malloc(Cosmap->nlocal*(N0 + N1)*sizeof(FLOAT));
        for(i=0; i< N0 + N1; i++)
		{
			if(i< N0)
			memcpy(&svdM[i*Cosmap->nlocal], &U[i*maxnlocal], Cosmap->nlocal*sizeof(FLOAT));
	        if(i >= N0)
			memcpy(&svdM[i*Cosmap->nlocal], &U1[(i-N0)*maxnlocal], Cosmap->nlocal*sizeof(FLOAT));
		}
    	MPI_Barrier(MPI_COMM_WORLD);
		S = (FLOAT *)malloc((N0 + N1) * sizeof(FLOAT));
		phgMatSvd(svdM, Cosmap, N0+N1, U, S, numprocs, rank);
		N0 = phgGetPodNumbasis(gamma3, S, N0 + N1);
		/************************Free Following values******************************/
		free(S1);  free(S); free(svdM);
	    free(PODmat);  free(PODC); free(PODu_p);free(y);  free(IPIV);
		free(PODb);  free(PODG1); free(PODG2);

		PODflag = 0;
	   }
       if(PODflag ==0)
       {
	       PODmat = (FLOAT *)malloc(N0*N0*sizeof(FLOAT));
		   PODC = (FLOAT *)malloc(N0*N0 * sizeof(FLOAT));
		   PODu_p = (FLOAT *)malloc(N0 * sizeof(FLOAT));
    
		   PODG1 = (FLOAT *)malloc(N0*N0*sizeof(FLOAT));
		   PODG2 = (FLOAT *)malloc(N0*N0*sizeof(FLOAT));
		   PODb = (FLOAT *)malloc(N0 * sizeof(FLOAT));
		   y = (FLOAT *)malloc(N0 * sizeof(FLOAT));
		   IPIV = (int *)malloc(N0 * sizeof(int));
/************build POD linear system***********************/	
	       if(GetPOD ==1)
		   {
			   memcpy(PODu0->data, &CosV->data[(nT1-1)*Cosmap->nlocal], Cosmap->nlocal*sizeof(*PODu0->data));
		       GetPOD = 0;
		   }
		   else
			   memcpy(PODu0->data, &CosV1->data[(CosnCLT-1)*Cosmap->nlocal], Cosmap->nlocal*sizeof(*PODu0->data));

		   for(i=0; i< N0; i++)
	       { 
	           memcpy(VL->data, &U[i*maxnlocal], Cosmap->nlocal*sizeof(FLOAT));
		      for(j=0; j< N0; j++)
	  	    {
			   memcpy(VR->data, &U[j*maxnlocal], Cosmap->nlocal*sizeof(FLOAT));
			   phgMatVec(0, 1.0, CosG1, VR, 0.0, &Vtmp);
			   PODG1[j*N0 + i] = phgVecDot(Vtmp, 0, VL, 0, NULL);
			   phgMatVec(0, 1.0, CosG2, VR, 0.0, &Vtmp);
			   PODG2[j*N0 + i] = phgVecDot(Vtmp, 0, VL, 0, NULL);
			   phgMatVec(0, 1.0, CosC, VR, 0.0, &Vtmp);
			   PODC[j*N0 + i]  = phgVecDot(Vtmp, 0, VL, 0, NULL);
		     }
           PODu_p[i] = phgVecDot(VL, 0, PODu0, 0, NULL);
          }
    }

	ptime = crtime;
    crtime += timestp;
    phgPrintf("current time layer: [%lf, %lf]\n", (FLOAT)ptime, (FLOAT)crtime);

    if(PODflag > 0)
		memcpy(PODu_p, PODb, N0*sizeof(double));

	if(FEMflag > 0)
	{
        phgDofFree(&FEMu_p);
		FEMu_p = phgDofCopy(FEMu_h, NULL, DOF_P1, NULL);
		FEMu_p->userfunc = DofInterpolation;
	}

    Cosb = phgMapCreateVec(Cosmap, 1);
	phgVecDisassemble(Cosb);
	phgDofSetDataByFunction(Cosf_h, func_f);
	build_rhs_Vec(timestp, Cosb, Cosf_h, Cosmap, Cosu_h);
     
	phgFEMLoopTime3(crtime, FEMu_h, FEMu_p, CosC, CosG1, CosG2, Cosb);

/**********************************************************/
    for(i=0; i<N0; i++)
    {    for(j=0 ;j< N0; j++)
	    	PODmat[i*N0 + j] = PODG1[i*N0 + j] + eta_t(crtime)*PODG2[i*N0 + j];
	    	memcpy(VL->data, &U[i*maxnlocal], Cosmap->nlocal*sizeof(FLOAT));
            PODb[i] = phgVecDot(VL, 0, Cosb, 0, NULL);
	} 
	dgemv_("T", &N0, &N0, &done, PODC, &N0, PODu_p, &ONE, &dzero, y, &ONE);
    for(i=0; i< N0; i++)
		PODb[i] += y[i];
    
	dgesv_(&N0, &ONE, PODmat, &N0, IPIV, PODb, &N0, &info);
	phgPrintf("info :%d\n", info);
    /***************transform PODb to DOF****************/
	phgPODDofToFEMDof(Cosmap, N0, U, maxnlocal, PODb, Cosu_h);
	k++;
	Cosuhs[k] = phgDofCopy(Cosu_h, NULL, DOF_P1, "u_p1");
	Cosuhs[k]->userfunc = DofInterpolation;
	/******************computing error indicator******************/
    phgDofAXPBY(1.0, Cosu_h, 0, &error);
	phgDofAXPBY(-1.0, FEMu_h, 1.0, &error);
    phgPrintf("error:%f\n", phgDofNormL2(error));
	PODflag++;
	FEMflag++;
  }

	phgPrintf("k: %d\t , TotalN1: %d\n", k, TotalN1);
	phgMatDestroy(&CosG1);
	phgMatDestroy(&CosG2);
	phgMatDestroy(&CosC);
	for(i=0; i< TotalN1; i++)
		phgDofFree(Cosuhs + i);
	phgFree(Cosuhs);
	free(PODmat);  free(PODC); free(PODu_p);free(y);  free(IPIV);
	free(PODb);  free(PODG1); free(PODG2);
	phgDofFree(&FEMu_p);
	phgDofFree(&FEMu_h);
	phgDofFree(&Cosu_h);
	phgDofFree(&Cosu_p);
	phgDofFree(&Cosf_h);
	phgDofFree(&error);
	phgDofFree(&CosB1);
	phgDofFree(&CosB2);
	phgDofFree(&CosB3);
	phgDofFree(&CosC1);
	phgDofFree(&CosC2);
	phgDofFree(&CosC3);

	return 0;
}
static void
phgRefineFEMLoopTime(FLOAT timestp, GRID *g, INT TotalN,  DOF **uhs, INT rank, char *snD1)
{
    INT i, j, k, flag = 0;
	FLOAT err;
	MAT *AA, *GG1, *GG2, *CC;
	VEC *bb, *bbtmp, *bbtmp1;
	MAP *mmap;
	SOLVER *ssolver;
	ELEMENT *ee;
	DOF *uu_h, **uuhs, *ff_h;  /* numerical solution at time t_n */
	DOF *uu_p, *ggrad_u;  /* numerical solution at time t_{n-1} */
	DOF  *BB1, *BB2, *BB3, *CC1, *CC2, *CC3, *eerror;
    
	uuhs = phgAlloc(TotalN * sizeof(*uhs));
	uu_p = phgDofNew(g, DOF_P1, 1, "uu_p", DofInterpolation);
	ff_h = phgDofNew(g, DOF_P1, 1, "ff_h", DofInterpolation);
	phgDofSetDataByFunction(uu_p, func_u0);
	uu_h = phgDofNew(g, DOF_P1, 1, "uu_h", DofInterpolation);
	eerror = phgDofNew(g, DOF_P1, 1, "eerror", DofInterpolation);	
	uuhs[0] = phgDofCopy(uu_p, NULL, DOF_P1, "uu_p1");
	uuhs[0] ->userfunc = DofInterpolation;
	k = 0;
	BB1=phgDofNew(g, DOF_P1, 1, "BB1", func_B1);
	BB2=phgDofNew(g, DOF_P1, 1, "BB2", func_B2);
	BB3=phgDofNew(g, DOF_P1, 1, "BB3", func_B3);
	CC1=phgDofNew(g, DOF_P1, 1, "CC1", func_C1);
	CC2=phgDofNew(g, DOF_P1, 1, "CC2", func_C2);
	CC3=phgDofNew(g, DOF_P1, 1, "CC3", func_C3);

	mmap = phgMapCreate(uu_h, NULL);
	GG1 = phgMapCreateMat(mmap, mmap);
	GG2 = phgMapCreateMat(mmap, mmap);
	CC = phgMapCreateMat(mmap, mmap);

	bbtmp = phgMapCreateVec(mmap, 1);
	phgVecDisassemble(bbtmp);

	build_rhs_Mat(CC, mmap, uu_h);
	build_stiff_Mat1(timestp, D0, GG1, mmap, BB1, BB2, BB3, uu_h);
	build_stiff_Mat2(timestp, GG2, mmap, CC1, CC2, CC3, uu_h);
	crtime = 0;
	while(crtime < T - 1e-8)
	{
		/********************************************************************/	
		ptime = crtime;
		crtime += timestp;
		phgPrintf("\n/********* start new time layer *************/\n");
		phgPrintf("current time layer: [%lf, %lf]\n", (double)ptime, (double)crtime);
		flag++;
		if (flag > 1)
		{   /* update u_p */
			phgDofFree(&uu_p);
			uu_p = phgDofCopy(uu_h, NULL, DOF_P1, "uu_p");
			uu_p->userfunc = DofInterpolation;
		}
		phgMapDofArraysToVec(mmap, 1, FALSE, &bbtmp, &uu_p, NULL);
		bb = phgMapCreateVec(mmap, 1);
		bbtmp1 = phgMapCreateVec(mmap, 1);
		phgVecDisassemble(bb);
		phgVecDisassemble(bbtmp1);

		phgMatVec(0, 1.0, CC, bbtmp, 0.0, &bb);	
		phgDofSetDataByFunction(ff_h, func_f);

		build_rhs_Vec(timestp, bbtmp1, ff_h, mmap, uu_h);
		phgVecAXPBY(1.0, bbtmp1, 1.0, &bb);

		AA = phgMapCreateMat(mmap, mmap);
		phgMatAXPBY(1.0, GG1, 0.0, &AA);
		phgMatAXPBY(eta_t(crtime), GG2, 1.0, &AA);

		ssolver = phgSolverCreate(SOLVER_DEFAULT, uu_h, NULL);
		ssolver->mat = AA;
		ssolver->rhs = bb;
		phgSolverSolve(ssolver, TRUE, uu_h, NULL);
		phgSolverDestroy(&ssolver);
        k++;
		uuhs[k] = phgDofCopy(uu_h, NULL, DOF_P1, "uu_p1");
		uuhs[k] ->userfunc = DofInterpolation;
		phgPrintf("\ntime step: %lg\n\n", (double)timestp);
	}

	for(i=0; i< TotalN; i++)
	{    
		phgDofAXPBY(1.0, uhs[i], 0, &eerror);
		phgDofAXPBY(-1.0, uuhs[i], 1.0, &eerror);
		err =phgDofNormL2(eerror);
		if(rank==0)
	    	fprintf(snD1, "%f\t%e\n", i*timestp, err);
	}
	for(i=0 ;i< TotalN; i++)
		phgDofFree(uuhs+i);
	phgFree(uuhs);
	phgDofFree(&eerror);
	phgDofFree(&BB1);
	phgDofFree(&BB2);
	phgDofFree(&BB3);
	phgDofFree(&ff_h);

	phgDofFree(&CC1);
	phgDofFree(&CC2);
	phgDofFree(&CC3);
	phgDofFree(&uu_h);
	phgDofFree(&uu_p);
}
