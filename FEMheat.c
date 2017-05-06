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
	int i, k, nT, TotalN;
	GRID *g;
	MAT *A, *G1, *G2, *C;
	VEC *b, *g1, *g2, *btmp, *btmp1, *f_h;
	MAP *map;
	SOLVER *solver;
	ELEMENT *e;
	DOF *u_h, **uhs;  /* numerical solution at time t_n */
	DOF *u_p, *grad_u;  /* numerical solution at time t_{n-1} */
	DOF  *B1, *B2, *B3, *C1, *C2, *C3;

	double tt[3], tt1[3], tt2[3];
	int flag = 0;
	long long Nelem = 0; /*  elements */
	long long Ndof = 0;  /*  DOF */
	char sn[40];/* output file name */
	FILE *snD, *snD1;
	int rank;
	static BOOLEAN export_mesh = FALSE;
	static BOOLEAN export_D11 = TRUE;
	FLOAT localnorm, globalnorm, D11, etaspace;
	i = 0;
	snD = fopen("FEMD11_0.5.text", "w");
	snD1 = fopen("FEMerr.text", "w");
	phgInit(&argc, &argv); 
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	g = phgNewGrid(-1);

	phgSetPeriodicity(g, X_MASK | Y_MASK | Z_MASK);
	if (!phgImport(g, fn, TRUE))
		phgError(1, "can't read file \"%s\".\n", fn);
	phgRefineAllElements(g, RefineN);

	TotalN = T/ stime + 1;
	uhs = phgAlloc(TotalN * sizeof(*uhs));

	u_p = phgDofNew(g, DOF_P1, 1, "u_h", DofInterpolation);
	phgDofSetDataByFunction(u_p, func_u0);
	u_h = phgDofNew(g, DOF_P1, 1, "u_h", DofInterpolation);
	f_h = phgDofNew(g, DOF_P1, 1, "f_h", DofInterpolation);
	uhs[0] = phgDofCopy(u_p, NULL, DOF_P1, "u_p1");
	uhs[0] ->userfunc = DofInterpolation;
	k = 1;
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
	Ndof = DofGetDataCountGlobal(u_h);

	if (export_mesh) 
	{
		*sn = '\0';
		sprintf(sn, "../exportdata/T100t0.2Refine8_2pi/a8_%d.vtk", i);
		phgPrintf("output mesh to %s:\n", sn);
		phgExportVTK(g, sn, u_p, NULL);
	}
	/* init parameters */
	phgGetTime(tt);
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
			phgDofFree(&u_p);
			u_p = phgDofCopy(u_h, NULL, DOF_P1, "u_p");
			u_p->userfunc = DofInterpolation;

		}
		phgMapDofArraysToVec(map, 1, FALSE, &btmp, &u_p, NULL);
		b = phgMapCreateVec(map, 1);
	    btmp1 = phgMapCreateVec(map, 1);
		phgVecDisassemble(b);
		phgVecDisassemble(btmp1);

		phgMatVec(0, 1.0, C, btmp, 0.0, &b);
		phgDofSetDataByFunction(f_h, func_f);

		build_rhs_Vec(stime, btmp1, f_h, map, u_h);
		phgVecAXPBY(1.0, btmp1, 1.0, &b);

		A = phgMapCreateMat(map,map);
		phgMatAXPBY(1.0, G1, 0.0, &A);
		phgMatAXPBY(eta_t(crtime), G2, 1.0, &A);

		solver = phgSolverCreate(SOLVER_GMRES, u_h, NULL);
		solver->mat = A;
		solver->rhs = b;
		phgSolverSolve(solver, TRUE, u_h, NULL);
		phgSolverDestroy(&solver);

		uhs[k] = phgDofCopy(u_h, NULL, DOF_P1, "u_p1");
		uhs[k] ->userfunc = DofInterpolation;

		k++;
		phgPrintf("\ntime step: %lg\n\n", (double)stime);

		i++;
		if (export_mesh) 
		{
			*sn = '\0';
			sprintf(sn, "../exportdata/T100t0.2Refine8_2pi/a8_%d.vtk", i);
			phgPrintf("output mesh to %s:\n", sn);
			phgExportVTK(g, sn, u_h, NULL);
		}
	}
	phgGetTime(tt2);
 /********************************computing err******************************/
//    phgRefineAllElements(g, 2);
//  	phgRefineFEMLoopTime(stime, g, TotalN, uhs, rank, snD1);
	/****************************************************************************/
	phgPrintf("\nFEMTotal time:%lfs\n", (double)(tt2[2] - tt[2]));
	/* cal total errors */
	phgPrintf("\n time steps:%d\n", flag);
	phgPrintf("\n elements:%ld\n", Nelem);
	phgPrintf("\n DOF:%ld\n", Ndof);

	for(i=0 ;i< TotalN; i++)
		phgDofFree(uhs+i);
	phgFree(uhs);
	phgDofFree(&B1);
	phgDofFree(&B2);
	phgDofFree(&B3);
	phgDofFree(&C1);
	phgDofFree(&C2);
	phgDofFree(&C3);
	phgDofFree(&u_h);
	phgDofFree(&f_h);
	phgDofFree(&u_p);
	phgFreeGrid(&g);
	phgFinalize();
	return 0;
}
