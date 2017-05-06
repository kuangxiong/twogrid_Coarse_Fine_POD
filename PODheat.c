/* ************************************************************
 *   problem: *   u_{t} - \Delta{u} = func_f  (x, t) \in \Omega X (0, T)
 *   u = func_g (x, t) \in \partial\Omega X [0, T]
 *   u = func_u0       x \in \Omega  t = 0.0
 *
 *   The space-time adaptation algorithm is based on:
 *	Zhiming CHEN et. al,
 *	An adaptive finite element method with reliable and efficient error 
 *	control for linear parabolic problems, Math. Comp. 73 (2004), 1163-1197.
 *   
 * $Id: heat.c,v 1.56 2015/10/16 01:07:47 zlb Exp $
 **************************************************************/

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
    int i, j, nT, k, rank, numprocs, TotalN;
	GRID *g;
	MAT *A, *G1, *G2, *C;
	VEC *b, *g1, *g2, *btmp, *btmp1, *V;
	MAP *map, *map1;
	SOLVER *solver;
	ELEMENT *e;
    DOF *u_h, *err;  /* numerical solution at time t_n */
    DOF *u_p, *grad_u;  /* numerical solution at time t_{n-1} */
    DOF  *B1, *B2, *B3, *C1, *C2, *C3, **uhs , *f_h;
    double tt[3], tt1[3], tt2[3], t1[3], t2[3], FEMt, PODt;
	FLOAT *U, *S, value, etaspace;
    int  context, flag = 0, info;
    long long Numdof = 0; 
    long long Nelem = 0; /*  elements */
    char sn[40];/* output file name */

    FILE *snD , *snD1, *snD2;
	int maxnlocal, mapnlocal, M, N, ZERO =0, ONE =1, N0;
	static BOOLEAN export_mesh = FALSE;
	static BOOLEAN export_D11 = TRUE;
	FLOAT localnorm, globalnorm, D11;

 //   FEMt = 0;  PODt = 0;
//	snD = fopen("PODD11_0.5.text", "w");
	snD1 = fopen("PODerr.text", "w");
	snD2 = fopen("PODerrindicator.text", "w");

    phgInit(&argc, &argv); 
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	g = phgNewGrid(-1);

	phgSetPeriodicity(g, X_MASK | Y_MASK | Z_MASK);
    if (!phgImport(g, fn, TRUE))
        phgError(1, "can't read file \"%s\".\n", fn);
	phgRefineAllElements(g, RefineN);
    err = phgDofNew(g, DOF_P1, 1, "err", DofInterpolation);
	u_p = phgDofNew(g, DOF_P1, 1, "u_p", func_u0);
    u_h = phgDofNew(g, DOF_P1, 1, "u_h", DofInterpolation);
    f_h = phgDofNew(g, DOF_P1, 1, "f_h", DofInterpolation);
	B1=phgDofNew(g, DOF_P1, 1, "B1", func_B1);
	B2=phgDofNew(g, DOF_P1, 1, "B2", func_B2);
	B3=phgDofNew(g, DOF_P1, 1, "B3", func_B3);
	C1=phgDofNew(g, DOF_P1, 1, "C1", func_C1);
	C2=phgDofNew(g, DOF_P1, 1, "C2", func_C2);
	C3=phgDofNew(g, DOF_P1, 1, "C3", func_C3);
    
	map = phgMapCreate(u_h, NULL);
	G1 = phgMapCreateMat(map, map);
	G2 = phgMapCreateMat(map, map);
	C = phgMapCreateMat(map, map);

    btmp = phgMapCreateVec(map, 1);

    phgVecDisassemble(btmp);

	build_rhs_Mat(C, map, u_h);
    build_stiff_Mat1(stime, D0, G1, map, B1, B2, B3, u_h);
    build_stiff_Mat2(stime, G2, map, C1, C2, C3, u_h);


    Nelem = g->nleaf_global;
    Numdof = DofGetDataCountGlobal(u_h);
     
   nT = PODT0 / stime +1;
   TotalN = T / stime +1;
   phgPrintf("test:nT:%d\n", nT);
   uhs = phgAlloc(TotalN * sizeof(*uhs));

   uhs[0] = phgDofCopy(u_p, NULL, DOF_P1, "u_p1");
   uhs[0]->userfunc = DofInterpolation;
   k = 1;
	if (export_mesh) 
		{
		    *sn = '\0';
			 sprintf(sn, "../exportdata/T100t0.2Refine4_2pi/a4_%d.vtk", k);
			 phgPrintf("output mesh to %s:\n", sn);
			 phgExportVTK(g, sn, uhs[0], NULL);
		}

    phgGetTime(tt);
 while(crtime < PODT0 - 1e-8)
   {
	/********************************************************************/	
		ptime = crtime;
        crtime += stime;
		phgPrintf("\n/********* start new time layer *************/\n");
        phgPrintf("current time layer: [%lf, %lf]\n", (double)ptime, (double)crtime);
        
        if (crtime > T)
    	{
            crtime = T;
            stime = T - ptime;
            phgPrintf("current time layer: [%lf, %lf]\n", (double)ptime, (double)crtime);
        }
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
     
//		phgFEMLoopTime(crtime, u_h, u_p, C, G1, G2, g1, g2);
        phgFEMLoopTime3(crtime, u_h, u_p, C, G1, G2, b);

		uhs[k] = phgDofCopy(u_h, NULL, DOF_P1, "u_p1");
		uhs[k]->userfunc = DofInterpolation;
		phgPrintf("\ntime step: %lg\n\n", (double)stime);
		if (export_mesh) 
		 {
		    *sn = '\0';
			 sprintf(sn, "../exportdata/T100t0.2Refine4_2pi/a4_%d.vtk", k);
			 phgPrintf("output mesh to %s:\n", sn);
			 phgExportVTK(g, sn, u_h, NULL);
		}
		k++;
    }
/*************************get Vector from Dofs**************************/
    MPI_Barrier(MPI_COMM_WORLD);
    V = phgMapCreateVec(map, nT);
    phgMapDofArraysToVec(map, nT, FALSE, &V, uhs, NULL);
/******************************SVD decomposition***********************/ 
    MPI_Reduce(&map->nlocal, &maxnlocal, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Bcast(&maxnlocal, 1, MPI_INT, 0, MPI_COMM_WORLD);
    M = numprocs * maxnlocal;
    N = nT;
    U = (double *)malloc((maxnlocal * M) * sizeof(double));
    S = (double *)malloc(N*sizeof(double));
    phgMatSvd(V->data, map, nT, U, S, numprocs, rank);
//	for(i=0; i< N; i++)
//		phgPrintf("%f\n", S[i]);
    N0 = phgGetPodNumbasis(gamma1, S, nT);
    phgPrintf("number of sparse basis:\t%d\t\n", N0, gamma1);
///*******************************POD********************************************/
/**************assem-ble POD basis matrix******/ 
 //   N0 = nT;
   /*********************set N0************************/
    FLOAT PODmat[N0][N0], PODC[N0][N0], PODu_p[N0], PODu_h[N0], test1[N0], test2[N0];
    VEC *PODu0, *VL, *VR, *Vtmp, *Vtf;
    FLOAT PODG1[N0][N0], PODG2[N0][N0], PODg1[N0], PODg2[N0], PODb[N0], y[N0];
    int  PODflag=0, IPIV[N0];
    double  done, dzero;
    done = 1.0;   dzero = 0.0;
    Vtmp = phgMapCreateVec(map, 1);
    PODu0 = phgMapCreateVec(map, 1);
    VL = phgMapCreateVec(map, 1);
    VR = phgMapCreateVec(map, 1);

    phgVecAssemble(Vtmp);
    phgVecAssemble(VL);
    phgVecAssemble(VR);
/************build POD linear system***********************/	
    memcpy(PODu0->data,&V->data[(nT-1)*map->nlocal], map->nlocal*sizeof(*PODu0->data));
    for(i=0; i< N0; i++)
    {
    	memcpy(VL->data, &U[i*maxnlocal], map->nlocal*sizeof(FLOAT));
    	for(j=0; j< N0; j++)
    	{
    	  memcpy(VR->data, &U[j*maxnlocal], map->nlocal*sizeof(double));
    	  phgMatVec(0, 1.0, G1, VR, 0.0, &Vtmp);
     	  PODG1[j][i] = phgVecDot(Vtmp, 0, VL, 0, NULL);
    	  phgMatVec(0, 1.0, G2, VR, 0.0, &Vtmp);
    	  PODG2[j][i] = phgVecDot(Vtmp, 0, VL, 0, NULL);
    	  phgMatVec(0, 1.0, C, VR, 0.0, &Vtmp);
    	  PODC[j][i]  = phgVecDot(Vtmp, 0, VL, 0, NULL);
    	}
 //       PODg1[i] = phgVecDot(VL, 0, g1, 0, NULL);
  //      PODg2[i] = phgVecDot(VL, 0, g2, 0, NULL);
        PODu_p[i] = phgVecDot(VL, 0, PODu0, 0, NULL);
    }
//    for(i=0 ; i< N0; i++)
//    	phgPrintf("%f\n", PODu_p[i]);
/**********start  loop time********/
 while(crtime < T - 1e-8)
 {
    ptime = crtime;
    crtime += stime;
    phgPrintf("current time layer: [%lf, %lf]\n", (double)ptime, (double)crtime);
    
    if(PODflag > 0)
    	memcpy(PODu_p, PODb, N0*sizeof(double));

    b = phgMapCreateVec(map, 1);
    phgVecDisassemble(b);
    phgDofSetDataByFunction(f_h, func_f);
	build_rhs_Vec(b, f_h, map, u_h);
	for(i=0; i<N0; i++)
    {	for(j=0 ;j< N0; j++)
 		PODmat[i][j] = PODG1[i][j] + eta_t(crtime)*PODG2[i][j];
    	memcpy(VL->data, &U[i*maxnlocal], map->nlocal*sizeof(FLOAT));
        PODb[i] = phgVecDot(VL, 0, b, 0, NULL);
    } 
	/*************test**************************/
//    for(i=0; i< N0; i++)
//		test1[i] = PODb[i];
//	/******************************************/
   for(i=0; i< N0; i++)
 	   y[i] =0;
    dgemv_("T", &N0, &N0, &done, PODC, &N0, PODu_p, &ONE, &dzero, y, &ONE);
   for(i=0; i< N0; i++)
 	PODb[i] += y[i];
  
    dgesv_(&N0, &ONE, PODmat, &N0, IPIV, PODb, &N0, &info);
    phgPrintf("info :%d\n", info);
	/***************transform PODb to DOF****************/ 
    phgPODDofToFEMDof(map, N0, U, maxnlocal, PODb, u_h);
/******************test***********************/
//    for(i=0; i< N0; i++)
//		test2[i] = test1[i] - (PODb[i] - PODu_p[i])/stime;
//    phgPODDofToFEMDof(map, N0, U, maxnlocal, test2, err);
//    etaspace = phgDofNormL2(err);
//	if(rank==0 )
//	    fprintf(snD2, "%f\t%f\n",  crtime, etaspace);
///*********************************************/
	
    if (export_mesh) 
    {
    	*sn = '\0';
    	sprintf(sn, "../exportdata/T100t0.2Refine4_2pi/a4_%d.vtk", k);
    	phgPrintf("output mesh to %s:\n", sn);
    	phgExportVTK(g, sn, u_h, NULL);
    }
    uhs[k] = phgDofCopy(u_h, NULL, DOF_P1, "u_p1");
    uhs[k]->userfunc = DofInterpolation;
    PODflag++;
    k++;
 }
    phgGetTime(tt2);
    phgPrintf("\nTotal Time:%f\n",tt2[2] - tt[2]);
///************computing D11************/    
//   double time = 0;
//   for(i=0; i< TotalN; i++)
//   {
//// 	D11 = D0*(1 + GetD11(uhs[i]));
//    D11 = phgDofNormL2(uhs[i]);
//	if(rank==0 && export_D11)
//       	fprintf(snD,"%f\t%f\n", time, D11);
//	time = time + stime;
//  }
   /****************computing FEMerr****************************/
 //  phgRefineAllElements(g, 2);
    phgRefineFEMLoopTime(stime, g, TotalN, uhs, rank, snD1);
// /* cal total errors */
    phgPrintf("\n number of POD basis :%d\t%f\n", N0, gamma1);
    phgPrintf("\n DOF:%"dFMT", elements:%"dFMT"\n", 
    		DofGetDataCountGlobal(u_h), g->nleaf_global	);
    phgPrintf("\nTotal Time:%f\n",tt2[2] - tt[2]);

	for(i=0; i< TotalN; i++)
	   phgDofFree(uhs+i);
    phgFree(uhs);
    phgDofFree(&B1);
    phgDofFree(&B2);
    phgDofFree(&B3);
    phgDofFree(&C1);
    phgDofFree(&C2);
    phgDofFree(&C3);
    phgDofFree(&u_h);
    phgDofFree(&err);
    phgDofFree(&f_h);
    phgDofFree(&u_p);
    phgFreeGrid(&g);
    phgFinalize();
    return 0;
}
