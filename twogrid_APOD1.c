/*************************************************************************
	> File Name: APODheat.c
	> Author:kuangxiong 
	> Mail: kuangxiong@lsec.cc.ac.cn
	> Created Time: 2016年10月21日 星期五 21时26分17秒
 ************************************************************************/

#include "phg.h"
#include "fun.h"

#if (PHG_VERSION_MAJOR <= 0 && PHG_VERSION_MINOR < 9)
# undef ELEMENT
typedef SIMPLEX ELEMENT;
#endif

#include <math.h>
#include <stdlib.h>

int
main(int argc, char *argv[])
{
    int i, j, nT, k, kk, rank, numprocs, TotalN, TotalN1, numt1, tmpt;
	GRID *g, *Cosg;
	MAT *A, *G1, *G2, *CosG1, *CosG2, *C, *CosC;
	VEC *b, *Cosb, *V, *V1;
	MAP *map, *Cosmap;
    DOF *u_h, *uu_h, *Cosu_h, *f_h, *Cosf_h;  /* numerical solution at time t_n */
    DOF  *B1, *B2, *B3, *C1, *C2, *C3, *u_p, **uhs, **FEMuhs, **Cosuhs;
    DOF  *CosB1, *CosB2, *CosB3, *CosC1, *CosC2, *CosC3, *Cosu_p;
    FLOAT tt[3], tt1[3], *svdM;
	FLOAT *U, *S, value, *U1, *S1, tmp;
    int  context, flag = 0, cflag, *timeflag;
    long long   Numdof = 0; 
    long long Nelem = 0; /*  elements */
    char sn[40];

    FILE *snD;
	int maxnlocal, mapnlocal, M, N, info, ZERO, ONE , N0, N1;
   
	ZERO =0; ONE =1;
	snD = fopen("APODerr.text", "w");

    phgInit(&argc, &argv); 
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

	g = phgNewGrid(-1);
	Cosg = phgNewGrid(-1);

	phgSetPeriodicity(g, X_MASK | Y_MASK | Z_MASK);
	phgSetPeriodicity(Cosg, X_MASK | Y_MASK | Z_MASK);
    
	if (!phgImport(g, fn, TRUE))
        phgError(1, "can't read file \"%s\".\n", fn);
    if (!phgImport(Cosg, fn, TRUE))
        phgError(1, "can't read file \"%s\".\n", fn);

	phgRefineAllElements(g, RefineN);
	phgRefineAllElements(Cosg, CosRefineN);
/*********************fine grid information********************************/
	u_p = phgDofNew(g, DOF_P1, 1, "u_p", DofInterpolation);
	phgDofSetDataByFunction(u_p, func_u0);
    u_h = phgDofNew(g, DOF_P1, 1, "u_h", DofInterpolation);
    uu_h = phgDofNew(g, DOF_P1, 1, "uu_h", DofInterpolation);
    f_h = phgDofNew(g, DOF_P1, 1, "f_h", DofInterpolation);
	B1 = phgDofNew(g, DOF_P1, 1, "B1", func_B1);
	B2 = phgDofNew(g, DOF_P1, 1, "B2", func_B2);
	B3 = phgDofNew(g, DOF_P1, 1, "B3", func_B3);
	C1 = phgDofNew(g, DOF_P1, 1, "C1", func_C1);
	C2=phgDofNew(g, DOF_P1, 1, "C2", func_C2);
	C3=phgDofNew(g, DOF_P1, 1, "C3", func_C3);
	
	map = phgMapCreate(u_h, NULL);
	G1 = phgMapCreateMat(map, map);
	G2 = phgMapCreateMat(map, map);
	C = phgMapCreateMat(map, map);
	build_rhs_Mat(C, map, u_h);
    build_stiff_Mat1(stime, D0, G1, map, B1, B2, B3, u_h);
    build_stiff_Mat2(stime, G2, map, C1, C2, C3, u_h);
/***********************coarse grid information*****************************/
//	Cosu_p = phgDofNew(Cosg, DOF_P1, 1, "Cosu_p", DofInterpolation);
//	phgDofSetDataByFunction(Cosu_p, func_u0);
    Cosu_h = phgDofNew(Cosg, DOF_P1, 1, "Cosu_h", DofInterpolation);
//    Cosf_h = phgDofNew(Cosg, DOF_P1, 1, "Cosf_h", DofInterpolation);
//	CosB1 = phgDofNew(Cosg, DOF_P1, 1, "CosB1", func_B1);
//	CosB2 = phgDofNew(Cosg, DOF_P1, 1, "CosB2", func_B2);
//	CosB3 = phgDofNew(Cosg, DOF_P1, 1, "CosB3", func_B3);
//	CosC1 = phgDofNew(Cosg, DOF_P1, 1, "CosC1", func_C1);
//	CosC2 = phgDofNew(Cosg, DOF_P1, 1, "CosC2", func_C2);
//	CosC3 = phgDofNew(Cosg, DOF_P1, 1, "CosC3", func_C3);
//
	Cosmap = phgMapCreate(Cosu_h, NULL);
	CosG1 = phgMapCreateMat(Cosmap, Cosmap);
	CosG2 = phgMapCreateMat(Cosmap, Cosmap);
	CosC = phgMapCreateMat(Cosmap, Cosmap);
    build_stiff_Mat1(stime1, D0, CosG1, Cosmap, CosB1, CosB2, CosB3, Cosu_h);
    build_stiff_Mat2(stime1, CosG2, Cosmap, CosC1, CosC2, CosC3, Cosu_h);
//****************************************************************************/
    Nelem = g->nleaf_global;
    Numdof = DofGetDataCountGlobal(u_h);
     
    nT = PODT0 / stime + 1;
    nT1 = PODT0 /stime1 + 1;
	TotalN = T / stime + 1;
	TotalN1 = T / stime1 + 1;
    numt1 = CLT / stime1 + 1;
    uhs = phgAlloc(TotalN * sizeof(*uhs));
//    Cosuhs = phgAlloc(numt * sizeof(*Cosuhs));

//	for(i=0; i< numt1; i++)
//	    Cosuhs[i]->userfunc = DofInterpolation;
         
//	Cosuhs[0] = phgDofCopy(Cosu_p, NULL, DOF_P1, NULL);
//    Cosuhs[0]->userfunc = DofInterpolation;
/*******************************************************************************/
/***************************set timeflag****************************************/
	timeflag = (int *)malloc(TotalN1 * sizeof(int));
    for(i=0; i< TotalN1; i++)
		timeflag[i] = 0;
/******************** FEM **********************************************/	
//	crtime = 0;
//	for(i=0; i< numt1; i++)
//	{	
//		crtime = crtime + stime1;
//	    Cosb = phgMapCreateVec(Cosmap, 1);
//		phgVecDisassemble(Cosb);
//		phgDofSetDataByFunction(Cosf_h, func_f);
//		build_rhs_Vec(stime1, Cosb, Cosf_h, Cosmap, Cosu_h);
//		phgFEMLoopTime3(crtime, Cosu_h, Cosuhs[i], CosC, CosG1, CosG2, Cosb);
//       Cosuhs[i+1] = phgDofCopy(Cosu_h, NULL, DOF_P1, "u_p1");
//	}
/****************************end timeflag*************************************/

	uhs[0] = phgDofCopy(u_p, NULL, DOF_P1, "u_p1");
    uhs[0]->userfunc = DofInterpolation;
//  Cosuhs[0] = phgDofCopy(u_p, NULL, DOF_P1, "u_p1");
//  Cosuhs[0]->userfunc = DofInterpolation;
	kk = 0;  k = 0;
    phgGetTime(tt);
//   /**********************************************************************/
//  /* init parameters */
   flag =0;
   while(crtime < PODT0 - 1e-8)
   {
	/********************************************************************/	
		ptime = crtime;
        crtime += stime;
		phgPrintf("\n/********* start new time layer *************/\n");
        phgPrintf("current time layer: [%lf, %lf]\n", (FLOAT)ptime, (FLOAT)crtime);
        
        flag++;
        if (flag > 1)
		{   /* update u_p */
          phgDofFree(&u_p);
          u_p = phgDofCopy(u_h, NULL, DOF_P1, "u_p1");
          u_p->userfunc = DofInterpolation;
        }
        b = phgMapCreateVec(map, 1);
		phgVecDisassemble(b);
		phgDofSetDataByFunction(f_h, func_f);
		build_rhs_Vec(stime, b, f_h, map, u_h);
	    phgFEMLoopTime3(crtime, u_h, u_p, C, G1, G2, b);
	
        k++;
		uhs[k] = phgDofCopy(u_h, NULL, DOF_P1, "u_p1");
		uhs[k]->userfunc = DofInterpolation;
		phgPrintf("\ntime step: %lg\n\n", (FLOAT)stime);
    }
/*************************get Vector from Dofs**************************/
	V = phgMapCreateVec(map, nT);
	phgMapDofArraysToVec(map, nT, FALSE, &V, uhs, NULL);
/******************************SVD decomposition***********************/
	MPI_Reduce(&map->nlocal, &maxnlocal, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Bcast(&maxnlocal, 1, MPI_INT, 0, MPI_COMM_WORLD);
	M = numprocs * maxnlocal;
	N = nT;

    U = (FLOAT *)malloc((maxnlocal * M ) * sizeof(FLOAT));
    U1 = (FLOAT *)malloc((maxnlocal * M ) * sizeof(FLOAT));

    S = (FLOAT *)malloc(N*sizeof(FLOAT));
    phgMatSvd(V->data, map, nT, U, S, numprocs, rank);
//	for(i=0; i< N; i++)
//		phgPrintf("%f\n", S[i]);
///*************************************************************/
	N0 = phgGetPodNumbasis(gamma1, S, nT);
	free(S);
/*******************************POD********************************************/
/************************assemble POD basis matrix*******************/ 
    FLOAT *PODmat, *PODC, *PODu_p;
    VEC *PODu0, *VL, *VR, *Vtmp, *Vtf;
	FLOAT *PODG1, *PODG2, *PODb, *y;
    int  PODflag = 0, *IPIV;
	double  done, dzero;
    done = 1.0;  dzero = 0.0;

	Vtmp = phgMapCreateVec(map, 1);
	PODu0 = phgMapCreateVec(map, 1);
    VL = phgMapCreateVec(map, 1);
	VR = phgMapCreateVec(map, 1);
/*********************start  loop time*****************************/
   tmpt = 0 ; PODflag = 0;
   phgPrintf("crtime: %f\n", crtime);
   while(crtime < T - 1e-8)
   {
    if(tmpt == 0 || tmpt == numt1)
	{
		if(tmpt == numt1)
	    {
			free(PODmat);  free(PODC); free(PODu_p);free(y);  free(IPIV);
			free(PODb);  free(PODG1); free(PODG2);
		}
/*********************get Cosgrid dof************/
	    V1 = phgMapCreateVec(map, numt);

		for(i = 0; i< numt-1; i++)
		{	
			crtime = crtime + stime1;
		    Cosb = phgMapCreateVec(Cosmap, 1);
			phgVecDisassemble(Cosb);
			phgDofSetDataByFunction(Cosf_h, func_f);
			build_rhs_Vec(stime1, Cosb, Cosf_h, Cosmap, Cosu_h);
			phgFEMLoopTime3(crtime, Cosu_h, Cosuhs[i], CosC, CosG1, CosG2, Cosb);
		   	phgDofCopy(Cosuhs[i+1], &Cosu_h, DOF_P1, "u_p1");
		}
//*********************SVD1*************************/
		phgMapDofArraysToVec(map, numt, FALSE, &V1, &Cosuhs[0], NULL);
		S1 = (FLOAT *)malloc(numt * sizeof(FLOAT));
        phgMatSvd(V1->data, map, numt, U1, S1, numprocs, rank);
		N1 = phgGetPodNumbasis(gamma2, S1, numt);
		for(i=0; i< numt; i++)
			phgDofFree(Cosuhs+i);
		phgFree(Cosuhs);
		crtime = tmp;
/***************************SVD2*****************************/
		svdM =(FLOAT *) malloc(map->nlocal*(N0 + N1)*sizeof(FLOAT));
        for(i=0; i< N0 + N1; i++)
		{
			if(i< N0)
			memcpy(&svdM[i*map->nlocal], &U[i*maxnlocal], map->nlocal*sizeof(FLOAT));
	        if(i >= N0)
			memcpy(&svdM[i*map->nlocal], &U1[(i-N0)*maxnlocal], map->nlocal*sizeof(FLOAT));
		}
    	MPI_Barrier(MPI_COMM_WORLD);
		S = (FLOAT *)malloc((N0 + N1) * sizeof(FLOAT));
		phgMatSvd(svdM, map, N0+N1, U, S, numprocs, rank);
		N0 = phgGetPodNumbasis(gamma3, S, N0 + N1);
		free(S1);  free(S); free(svdM);
		phgPrintf("number of POD basis:%d\n", N0);
/*************************************************************/
	    PODmat = (FLOAT *)malloc(N0*N0*sizeof(FLOAT));
		PODC = (FLOAT *)malloc(N0*N0 * sizeof(FLOAT));
		PODu_p = (FLOAT *)malloc(N0 * sizeof(FLOAT));
    
		PODG1 = (FLOAT *)malloc(N0*N0*sizeof(FLOAT));
		PODG2 = (FLOAT *)malloc(N0*N0*sizeof(FLOAT));
		PODb = (FLOAT *)malloc(N0 * sizeof(FLOAT));
		y = (FLOAT *)malloc(N0 * sizeof(FLOAT));
		IPIV = (int *)malloc(N0 * sizeof(int));
/************build POD linear system**************************/	
    	phgMapDofArraysToVec(map, 1, FALSE, &PODu0, &uhs[k], NULL);
    	for(i=0; i< N0; i++)
    	{
	    	memcpy(VL->data, &U[i*maxnlocal], map->nlocal*sizeof(FLOAT));
	    	for(j=0; j< N0; j++)
	        {
		       memcpy(VR->data, &U[j*maxnlocal], map->nlocal*sizeof(FLOAT));
			   phgMatVec(0, 1.0, G1, VR, 0.0, &Vtmp);
	 	       PODG1[j*N0 + i] = phgVecDot(Vtmp, 0, VL, 0, NULL);
		       phgMatVec(0, 1.0, G2, VR, 0.0, &Vtmp);
		       PODG2[j*N0 + i] = phgVecDot(Vtmp, 0, VL, 0, NULL);
		       phgMatVec(0, 1.0, C, VR, 0.0, &Vtmp);
		       PODC[j*N0 + i]  = phgVecDot(Vtmp, 0, VL, 0, NULL);
		    }
           PODu_p[i] = phgVecDot(VL, 0, PODu0, 0, NULL);
  	   }
		tmpt = 1;
		PODflag = 0;
	}
	ptime = crtime;
    crtime += stime;
    phgPrintf("current time layer: [%lf, %lf]\n", (FLOAT)ptime, (FLOAT)crtime);
	if(PODflag > 0)
		memcpy(PODu_p, PODb, N0*sizeof(FLOAT));

    b = phgMapCreateVec(map, 1);
	phgVecDisassemble(b);
	phgDofSetDataByFunction(f_h, func_f);
	build_rhs_Vec(stime, b, f_h, map, u_h);
/**********************************************************/
    for(i=0; i < N0; i++)
    {    for(j=0 ;j< N0; j++)
	    	PODmat[i*N0 + j] = PODG1[i*N0 + j] + eta_t(crtime)*PODG2[i*N0 + j];
	    	memcpy(VL->data, &U[i*maxnlocal], map->nlocal*sizeof(FLOAT));
            PODb[i] = phgVecDot(VL, 0, b, 0, NULL);
	} 
	dgemv_("T", &N0, &N0, &done, PODC, &N0, PODu_p, &ONE, &dzero, y, &ONE);
    
	for(i=0; i< N0; i++)
		PODb[i] += y[i];
  
    dgesv_(&N0, &ONE, PODmat, &N0, IPIV, PODb, &N0, &info);

	phgPrintf("info :%d\n", info);
    /***************transform PODb to DOF****************/
	phgPODDofToFEMDof(map, N0, U, maxnlocal, PODb, u_h);
	k++;
    uhs[k] = phgDofCopy(u_h, NULL, DOF_P1, "u_p1");
	uhs[k]->userfunc = DofInterpolation;
    tmpt++;
	PODflag++;
 }
// /*****************computing D11*******************************/
    phgGetTime(tt1);
//*******************************computing err***********************************/
//    phgRefineAllElements(g,2);
	phgRefineFEMLoopTime(stime, g, TotalN, uhs, rank, snD);
///*************************************************************/	
	phgPrintf("number of POD basis:%d\n", N0);
	phgPrintf("\n Total time: %f \n ", (FLOAT)(tt1[2]-tt[2]));
	phgPrintf("\n DOF:%"dFMT", elements:%"dFMT"\n", 
  		DofGetDataCountGlobal(u_h), g->nleaf_global);

	free(PODmat);  free(PODC); free(PODu_p);free(y);  free(IPIV);
	free(PODb);  free(PODG1); free(PODG2);
	free(U); free(U1); 
	phgMatDestroy(&G1);
	phgMatDestroy(&G2);
	phgMatDestroy(&CosG1);
	phgMatDestroy(&CosG2);
	for(i=0; i< TotalN; i++)
	   phgDofFree(uhs+i);
    phgFree(uhs);
	phgDofFree(&B1);  phgDofFree(&CosB1);
    phgDofFree(&B2);  phgDofFree(&CosB2);
    phgDofFree(&B3);  phgDofFree(&CosB3);
    phgDofFree(&C1);  phgDofFree(&CosC1);
    phgDofFree(&C2);  phgDofFree(&CosC2);
    phgDofFree(&C3);  phgDofFree(&CosC3);
    phgDofFree(&u_h);  phgDofFree(&Cosu_h);
    phgDofFree(&uu_h); 
    phgDofFree(&f_h);  phgDofFree(&Cosf_h);
    phgDofFree(&u_p);  phgDofFree(&Cosu_p);
    phgFreeGrid(&g);
    phgFreeGrid(&Cosg);
	phgFinalize();
    return 0;
}
