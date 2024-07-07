/*! \file	IterPrepare.cpp
 *
 * \brief	Preparation for the calculation of stable interface structure;
 */

#include "Data.h"
#include "Head.h"
#include "DataOperators.h"
#include "Mytimer.h"
#include "functs.h"


/**
 * \brief	Initialization and preparation;
 */
void fn_iter_prepare		(	stu_system_param		*ssysparam,
									bool				isTest )
{
	if ( myrank == 0 )
	{
		printf(" <========== Preparation for the calculation of ");
		printf("stable interface structure ==========> \n\n");
	}
	mytimer_t timer;
	timer.reset();
	timer.start();

	/* calculate 'Gsquare'; */
	fn_obt_system_Gsquare ( ssysparam );

	/* calculate the interaction potential operator; */
	fn_obt_interact_grad ( ssysparam, isTest );

	/* calculate the interaction potential operator on the boundaries; */
	fn_obt_interact_bnd_grad ( ssysparam, isTest );

	/* memory allocation; */
	fn_system_fftw_memory_allocation ( ssysparam );

	MPI_Barrier ( MPI_COMM_WORLD );
	timer.pause();
	if ( myrank == 0 )
		printf("\n\t ***** time cost of preparation: %f seconds *****\n\n", 
				timer.get_current_time());
}


/**
 * \brief	Get Gsquare from 'ssysparam'; i.e. |R'PBk|^2;
 */
void fn_obt_system_Gsquare ( stu_system_param	*ssysparam )
{
	if ( myrank == 0 ) printf("Calculate Gsquare.\n");

	/* initialization; */
	int cplxReDofs = ssysparam->cplxReDofs;
	int dimRePhy   = ssysparam->dimRePhy;
	tmat<double> projPlaneTrans	= fn_tmat_init<double>	( dimRePhy, alloc_local_sys );
	ssysparam->sfftv->Gsquare = fn_tvec_init<double> ( cplxReDofs );
	fn_tvec_setZero<double> ( ssysparam->sfftv->Gsquare );

	/* transpose "projPlane"; */
	for ( int i = 0; i < dimRePhy; i++ )
		for ( int j = 0; j < alloc_local_sys; j++ )
			projPlaneTrans.val[i][j] = ssysparam->sfftv->projPlane.val[j][i];
	MPI_Barrier ( MPI_COMM_WORLD );

	/* the total amount of data may be greater than cplxReDofs; */
	int dataAmount = 0;
	for ( int i = 0; i < nprocs; i++ )
		dataAmount += recvCount_sys[i];
	
	/* collect all informations; */
	tmat<double> projPlaneAll	= fn_tmat_init<double>	( dimRePhy, dataAmount );
	int sendcount = alloc_local_sys;
	for ( int i = 0; i < dimRePhy; i++ )
		MPI_Allgatherv ( projPlaneTrans.val[i], sendcount, MPI_DOUBLE, 
					projPlaneAll.val[i], recvCount_sys, displs_sys, 
					MPI_DOUBLE, MPI_COMM_WORLD );

	/* calculate Gsquare; */
	for (int i = 0; i < cplxReDofs; i++)
	{
		for (int kk = 0; kk < dimRePhy; kk++)
			ssysparam->sfftv->Gsquare.val[i] += pow(projPlaneAll.val[kk][i], 2);
	}

	/* save data for comparison; */
//	char fwFile[FILELEN];
//	sprintf(fwFile, "%s/rank%d/gsquare.dat", rsltDir, myrank);
//	fn_tvec_save<double> ( ssysparam->sfftv->Gsquare, fwFile );

	/* release memory; */
	fn_tmat_free<double> ( projPlaneTrans );
	fn_tmat_free<double> ( projPlaneAll );
}


/**
 * \brief	Calculate the interaction potential operator;
 *			(\lap+1)^2 for LB;
 *			(\lap+1)^2 (\lap+q^2)^2 for LP;
 *
 */
void fn_obt_interact_grad ( stu_system_param	*ssysparam,
								bool			isTest )
{
	if ( myrank == 0 ) printf("Calculate the interaction potential operator.\n");

	int	cplxReDofs = ssysparam->cplxReDofs;
	int	nd	 = ssysparam->sGJPv->nd;	// innSMatdsJJ: [nd, nd];
	int	xlen = ssysparam->sGJPv->xlen;
	int glen = cplxReDofs * nd;

	tCCSmat<double> GinndsJJaddRslt;
	if ( strcmp( model_type, "LB" ) == 0 )
	{
		/* calculate (\lap+1)^2 with its expansion for LB model; */
		/* fourth order; */
		if ( myrank == 0 ) printf("---> fourth order.\n");
		tCCSmat<double> G0innd2JJ = fn_tensor_diag_dCCSmat ( 
				ssysparam->sfftv->Gsquare, ssysparam->sGJPv->innSMatd2JJ, 0 );
		tCCSmat<double> G1innd1JJ = fn_tensor_diag_dCCSmat ( 
				ssysparam->sfftv->Gsquare, ssysparam->sGJPv->innSMatd1JJ, 1 );
		tCCSmat<double> G2innd0JJ = fn_tensor_diag_dCCSmat ( 
				ssysparam->sfftv->Gsquare, ssysparam->sGJPv->innSMatd0JJ, 2 );
		/* second order; */
		if ( myrank == 0 ) printf("---> second order.\n");
		tCCSmat<double> G0innd1JJ = fn_tensor_diag_dCCSmat ( 
				ssysparam->sfftv->Gsquare, ssysparam->sGJPv->innSMatd1JJ, 0 );
		tCCSmat<double> G1innd0JJ = fn_tensor_diag_dCCSmat ( 
				ssysparam->sfftv->Gsquare, ssysparam->sGJPv->innSMatd0JJ, 1 );
		/* zero order; */
		if ( myrank == 0 ) printf("---> zero order.\n");
		tCCSmat<double> G0innd0JJ = fn_tensor_diag_dCCSmat ( 
				ssysparam->sfftv->Gsquare, ssysparam->sGJPv->innSMatd0JJ, 0 );

		/* addition result; */
		if ( myrank == 0 ) printf("---> add up.\n");
		tCCSmat<double> GinndsJJadd0 = fn_add_dCCSmat ( G0innd2JJ,	  G1innd1JJ, 1.0,  2.0 );
		tCCSmat<double> GinndsJJadd1 = fn_add_dCCSmat ( GinndsJJadd0, G2innd0JJ, 1.0,  1.0 );
		tCCSmat<double> GinndsJJadd2 = fn_add_dCCSmat ( GinndsJJadd1, G0innd1JJ, 1.0, -2.0 );
		tCCSmat<double> GinndsJJadd3 = fn_add_dCCSmat ( GinndsJJadd2, G1innd0JJ, 1.0, -2.0 );
		GinndsJJaddRslt				 = fn_add_dCCSmat ( GinndsJJadd3, G0innd0JJ, 1.0,  1.0 );

		/* release memory for temporary variables; */
		/* tensor product; */
		fn_tCCSmat_free<double> ( G0innd2JJ );
		fn_tCCSmat_free<double> ( G1innd1JJ );
		fn_tCCSmat_free<double> ( G2innd0JJ );
		fn_tCCSmat_free<double> ( G0innd1JJ );
		fn_tCCSmat_free<double> ( G1innd0JJ );
		fn_tCCSmat_free<double> ( G0innd0JJ );
		/* addition; */
		fn_tCCSmat_free<double> ( GinndsJJadd0 );
		fn_tCCSmat_free<double> ( GinndsJJadd1 );
		fn_tCCSmat_free<double> ( GinndsJJadd2 );
		fn_tCCSmat_free<double> ( GinndsJJadd3 );
	}
	else if ( strcmp( model_type, "LP" ) == 0 )
	{
		double q = ssysparam->scale_val[1];
		double q2 = pow(q, 2);
		double q4 = pow(q2, 2);

		/* calculate (\lap+1)^2 (\lap+q^2)^2 with its expansion for LP model; */
		/* eight order; */
		if ( myrank == 0 ) printf("---> eight order.\n");
		tCCSmat<double> G0innd4JJ = fn_tensor_diag_dCCSmat ( 
				ssysparam->sfftv->Gsquare, ssysparam->sGJPv->innSMatd4JJ, 0 );
		tCCSmat<double> G1innd3JJ = fn_tensor_diag_dCCSmat ( 
				ssysparam->sfftv->Gsquare, ssysparam->sGJPv->innSMatd3JJ, 1 );
		tCCSmat<double> G2innd2JJ = fn_tensor_diag_dCCSmat ( 
				ssysparam->sfftv->Gsquare, ssysparam->sGJPv->innSMatd2JJ, 2 );
		tCCSmat<double> G3innd1JJ = fn_tensor_diag_dCCSmat ( 
				ssysparam->sfftv->Gsquare, ssysparam->sGJPv->innSMatd1JJ, 3 );
		tCCSmat<double> G4innd0JJ = fn_tensor_diag_dCCSmat ( 
				ssysparam->sfftv->Gsquare, ssysparam->sGJPv->innSMatd0JJ, 4 );
		/* sixth order; */
		if ( myrank == 0 ) printf("---> sixth order.\n");
		tCCSmat<double> G0innd3JJ = fn_tensor_diag_dCCSmat ( 
				ssysparam->sfftv->Gsquare, ssysparam->sGJPv->innSMatd3JJ, 0 );
		tCCSmat<double> G1innd2JJ = fn_tensor_diag_dCCSmat ( 
				ssysparam->sfftv->Gsquare, ssysparam->sGJPv->innSMatd2JJ, 1 );
		tCCSmat<double> G2innd1JJ = fn_tensor_diag_dCCSmat ( 
				ssysparam->sfftv->Gsquare, ssysparam->sGJPv->innSMatd1JJ, 2 );
		tCCSmat<double> G3innd0JJ = fn_tensor_diag_dCCSmat ( 
				ssysparam->sfftv->Gsquare, ssysparam->sGJPv->innSMatd0JJ, 3 );
		/* fourth order; */
		if ( myrank == 0 ) printf("---> fourth order.\n");
		tCCSmat<double> G0innd2JJ = fn_tensor_diag_dCCSmat ( 
				ssysparam->sfftv->Gsquare, ssysparam->sGJPv->innSMatd2JJ, 0 );
		tCCSmat<double> G1innd1JJ = fn_tensor_diag_dCCSmat ( 
				ssysparam->sfftv->Gsquare, ssysparam->sGJPv->innSMatd1JJ, 1 );
		tCCSmat<double> G2innd0JJ = fn_tensor_diag_dCCSmat ( 
				ssysparam->sfftv->Gsquare, ssysparam->sGJPv->innSMatd0JJ, 2 );
		/* second order; */
		if ( myrank == 0 ) printf("---> second order.\n");
		tCCSmat<double> G0innd1JJ = fn_tensor_diag_dCCSmat ( 
				ssysparam->sfftv->Gsquare, ssysparam->sGJPv->innSMatd1JJ, 0 );
		tCCSmat<double> G1innd0JJ = fn_tensor_diag_dCCSmat ( 
				ssysparam->sfftv->Gsquare, ssysparam->sGJPv->innSMatd0JJ, 1 );
		/* zero order; */
		if ( myrank == 0 ) printf("---> zero order.\n");
		tCCSmat<double> G0innd0JJ = fn_tensor_diag_dCCSmat ( 
				ssysparam->sfftv->Gsquare, ssysparam->sGJPv->innSMatd0JJ, 0 );

		/* addition result; */
		if ( myrank == 0 ) printf("---> add up.\n");
		double coeffTmp0 = -2.0 * (q2 + 1.0);
		double coeffTmp1 =  3.0 * coeffTmp0;
		double coeffTmp2 = q4 + 4.0*q2 + 1.0;
		double coeffTmp3 =  2.0 * coeffTmp2;
		double coeffTmp4 = -2.0 * (q4 + q2);
		tCCSmat<double> GinndsJJadd0  = fn_add_dCCSmat ( G0innd4JJ,		G1innd3JJ, 1.0,  4.0 );
		tCCSmat<double> GinndsJJadd1  = fn_add_dCCSmat ( GinndsJJadd0,	G2innd2JJ, 1.0,  6.0 );
		tCCSmat<double> GinndsJJadd2  = fn_add_dCCSmat ( GinndsJJadd1,	G3innd1JJ, 1.0,  4.0 );
		tCCSmat<double> GinndsJJadd3  = fn_add_dCCSmat ( GinndsJJadd2,	G4innd0JJ, 1.0,  1.0 );
		tCCSmat<double> GinndsJJadd4  = fn_add_dCCSmat ( GinndsJJadd3,	G0innd3JJ, 1.0, coeffTmp0 );
		tCCSmat<double> GinndsJJadd5  = fn_add_dCCSmat ( GinndsJJadd4,	G1innd2JJ, 1.0, coeffTmp1 );
		tCCSmat<double> GinndsJJadd6  = fn_add_dCCSmat ( GinndsJJadd5,	G2innd1JJ, 1.0, coeffTmp1 );
		tCCSmat<double> GinndsJJadd7  = fn_add_dCCSmat ( GinndsJJadd6,	G3innd0JJ, 1.0, coeffTmp0 );
		tCCSmat<double> GinndsJJadd8  = fn_add_dCCSmat ( GinndsJJadd7,	G0innd2JJ, 1.0, coeffTmp2 );
		tCCSmat<double> GinndsJJadd9  = fn_add_dCCSmat ( GinndsJJadd8,  G1innd1JJ, 1.0, coeffTmp3 );
		tCCSmat<double> GinndsJJadd10 = fn_add_dCCSmat ( GinndsJJadd9,  G2innd0JJ, 1.0, coeffTmp2 );
		tCCSmat<double> GinndsJJadd11 = fn_add_dCCSmat ( GinndsJJadd10, G0innd1JJ, 1.0, coeffTmp4 );
		tCCSmat<double> GinndsJJadd12 = fn_add_dCCSmat ( GinndsJJadd11, G1innd0JJ, 1.0, coeffTmp4 );
		GinndsJJaddRslt				  = fn_add_dCCSmat ( GinndsJJadd12, G0innd0JJ, 1.0, q4 );

		/* release memory for temporary variables; */
		/* tensor product; */
		fn_tCCSmat_free<double> ( G0innd4JJ );
		fn_tCCSmat_free<double> ( G1innd3JJ );
		fn_tCCSmat_free<double> ( G2innd2JJ );
		fn_tCCSmat_free<double> ( G3innd1JJ );
		fn_tCCSmat_free<double> ( G4innd0JJ );
		fn_tCCSmat_free<double> ( G0innd3JJ );
		fn_tCCSmat_free<double> ( G1innd2JJ );
		fn_tCCSmat_free<double> ( G2innd1JJ );
		fn_tCCSmat_free<double> ( G3innd0JJ );
		fn_tCCSmat_free<double> ( G0innd2JJ );
		fn_tCCSmat_free<double> ( G1innd1JJ );
		fn_tCCSmat_free<double> ( G2innd0JJ );
		fn_tCCSmat_free<double> ( G0innd1JJ );
		fn_tCCSmat_free<double> ( G1innd0JJ );
		fn_tCCSmat_free<double> ( G0innd0JJ );
		/* addition; */
		fn_tCCSmat_free<double> ( GinndsJJadd0  );
		fn_tCCSmat_free<double> ( GinndsJJadd1  );
		fn_tCCSmat_free<double> ( GinndsJJadd2  );
		fn_tCCSmat_free<double> ( GinndsJJadd3  );
		fn_tCCSmat_free<double> ( GinndsJJadd4  );
		fn_tCCSmat_free<double> ( GinndsJJadd5  );
		fn_tCCSmat_free<double> ( GinndsJJadd6  );
		fn_tCCSmat_free<double> ( GinndsJJadd7  );
		fn_tCCSmat_free<double> ( GinndsJJadd8  );
		fn_tCCSmat_free<double> ( GinndsJJadd9  );
		fn_tCCSmat_free<double> ( GinndsJJadd10 );
		fn_tCCSmat_free<double> ( GinndsJJadd11 );
		fn_tCCSmat_free<double> ( GinndsJJadd12 );
	}
	else
	{
		GinndsJJaddRslt				  = fn_tCCSmat_init<double> ( glen, glen, 0 );
		if ( myrank == 0 )
		{
			printf("Error use 'fn_obt_interact_grad'\n");
			printf("'model_type' %s out of consideration.\n", model_type);
		}
	}

	/* copy the calculation result; */
	ssysparam->interact_grad = fn_tCCSmat_init<double> ( glen, glen, GinndsJJaddRslt.nnz );
	for ( int i = 0; i < glen+1; i++ )
		ssysparam->interact_grad.JA[i] = GinndsJJaddRslt.JA[i];
	for ( int i = 0; i < GinndsJJaddRslt.nnz; i++ )
		ssysparam->interact_grad.IA[i] = GinndsJJaddRslt.IA[i];
	for ( int i = 0; i < GinndsJJaddRslt.nnz; i++ )
		ssysparam->interact_grad.val[i] = GinndsJJaddRslt.val[i];

	/* release memory; */
	fn_tCCSmat_free<double> ( GinndsJJaddRslt );

	/* save data; */
	if ( isTest )
	{
		char fwFile[FILELEN];
		sprintf(fwFile, "%s/rank%d/interact_grad.dat", rsltDir, myrank);
		fn_tCCSmat_save<double> ( ssysparam->interact_grad, fwFile );
	}
}


/**
 * \brief	Calculate the interaction potential operator on the boundaries;
 *			(\lap+1)^2 for LB;
 *			(\lap+1)^2 (\lap+q^2)^2 for LP;
 *
 */
void fn_obt_interact_bnd_grad ( stu_system_param	*ssysparam,
									bool			isTest )
{
	if ( myrank == 0 ) printf("Calculate the interaction potential operator on the boundaries.\n");

	int xlen	 = ssysparam->scbndv->xlen;
	int cplxDofs = ssysparam->scbndv->cplxDofs;
	int bndLen	 = xlen * alloc_local_sys;
	ssysparam->interact_bnd_grad = fn_tvec_init<fftw_complex> ( bndLen );

	int ind = 0;
	if ( strcmp(model_type, "LB") == 0 )
	{
		/* calculate (\lap+1)^2 with its expansion; */
		double pk2, pk4;
		for ( int j = 0; j < alloc_local_sys; j++ )
		{
			for ( int i = 0; i < xlen; i++ )
			{
				ind = i * alloc_local_sys + j;
				int ind1 = displs_sys[myrank] + j;
				pk2 = ssysparam->sfftv->Gsquare.val[ind1];	// |R'PB k|^2;
				pk4 = pow(pk2, 2);
				ssysparam->interact_bnd_grad.val[ind][0] = 
					(2.0-2.0*pk2)	  * ssysparam->scbndv->d2bnd.val[ind][0] +
					(pk4-2.0*pk2+1.0) * ssysparam->scbndv->d0bnd.val[ind][0];
				ssysparam->interact_bnd_grad.val[ind][1] = 
					(2.0-2.0*pk2)	  * ssysparam->scbndv->d2bnd.val[ind][1] +
					(pk4-2.0*pk2+1.0) * ssysparam->scbndv->d0bnd.val[ind][1];
			}
		}
	}
	else if ( strcmp(model_type, "LP") == 0 )
	{
		double q = ssysparam->scale_val[1];
		double q2 = pow(q, 2);
		double q4 = pow(q2, 2);
		double coeffTmp0 = 2.0 * (q2 + 1.0);
		double coeffTmp1 = q4 + 4.0*q2 + 1.0;
		double coeffTmp2 = 2.0 * (q4 + q2);

		/* calculate (\lap+1)^2 (\lap+q^2)^2 with its expansion; */
		double pk2, pk4, pk6, pk8;
		double coeff0, coeff2, coeff4, coeff6;	// coefficients before p0, p2, p4, p6;
		for ( int j = 0; j < alloc_local_sys; j++ )
		{
			for ( int i = 0; i < xlen; i++ )
			{
				ind = i * alloc_local_sys + j;
				int ind1 = displs_sys[myrank] + j;
				pk2 = ssysparam->sfftv->Gsquare.val[ind1];	// |R'PB k|^2;
				pk4 = pow(pk2, 2);
				pk6 = pow(pk2, 3);
				pk8 = pow(pk2, 4);
				coeff6 = coeffTmp0 - 4.0*pk2;
				coeff4 = coeffTmp1 - 3.0*coeffTmp0*pk2 + 6.0*pk4;
				coeff2 = coeffTmp2 - 2.0*coeffTmp1*pk2 + 3.0*coeffTmp0*pk4 - 4.0*pk6;
				coeff0 = q4 - coeffTmp2*pk2 + coeffTmp1*pk4 - coeffTmp0*pk6 + pk8;
				/* (\lap+1)^2 (\lap+q^2)^2 on boundaries; */
				ssysparam->interact_bnd_grad.val[ind][0] = 
					coeff6 * ssysparam->scbndv->d6bnd.val[ind][0] +
					coeff4 * ssysparam->scbndv->d4bnd.val[ind][0] +
					coeff2 * ssysparam->scbndv->d2bnd.val[ind][0] +
					coeff0 * ssysparam->scbndv->d0bnd.val[ind][0];
				ssysparam->interact_bnd_grad.val[ind][1] = 
					coeff6 * ssysparam->scbndv->d6bnd.val[ind][1] +
					coeff4 * ssysparam->scbndv->d4bnd.val[ind][1] +
					coeff2 * ssysparam->scbndv->d2bnd.val[ind][1] +
					coeff0 * ssysparam->scbndv->d0bnd.val[ind][1];
			}
		}
	}
	else
	{
		if ( myrank == 0 )
		{
			printf("Error use 'fn_obt_interact_bnd_grad'\n");
			printf("'model_type' %s out of consideration.\n", model_type);
		}
	}

	/* save data; */
	if ( isTest )
	{
		char fwFile[FILELEN];
		sprintf(fwFile, "%s/rank%d/interact_bnd_grad.dat", rsltDir, myrank);
		FILE *fwPath = fopen(fwFile, "w");
		fprintf(fwPath, "%d\t%ld\t%d\n", xlen, alloc_local_sys, bndLen);
		for ( int j = 0; j < bndLen; j++ )
		{
			fprintf(fwPath, "%+.15E\t%+.15E\n", ssysparam->interact_bnd_grad.val[j][0], 
							ssysparam->interact_bnd_grad.val[j][1]);
		}
		fclose(fwPath);
	}
}


/**
 * \brief	Generate the iteration matrix for the SIS method;
 */
void fn_obt_iter_matrix		(	stu_system_param		*ssysparam,
									double				step_size,
									int					iterator,
									bool				isTest )
{
	double model_xi	 = ssysparam->model_xi; // has been squared;
	int cplxReDofs	 = ssysparam->cplxReDofs;

	/* model_xi * interact_grad + 1/dt * innSMatd0JJ; */
	tCCSmat<double> G0innd0JJ = fn_tensor_diag_dCCSmat ( 
			ssysparam->sfftv->Gsquare, ssysparam->sGJPv->innSMatd0JJ, 0 );
	ssysparam->iter_matrix	  = fn_add_dCCSmat ( 
			ssysparam->interact_grad, G0innd0JJ, model_xi, 1.0/step_size );

	/* release memory; */
	fn_tCCSmat_free<double> ( G0innd0JJ );

	if ( isTest )
	{
		char fwFile[FILELEN];
		sprintf(fwFile, "%s/rank%d/iter_matrix%d.dat", rsltDir, myrank, iterator);
		fn_tCCSmat_save ( ssysparam->iter_matrix, fwFile );
	}
}


/**
 * \brief	Project 'rhoJCplx' to Fourier space;
 *
 * \return	rhoCplx = rhoJCplx * d0JJ + d0bnd;
 */
tvec<fftw_complex> fn_back_Four_space	(	stu_system_param		*ssysparam,
											tvec<fftw_complex>		rhoJCplx )
{
	/* parameters; */
	int		nd		= ssysparam->scbndv->nd;
	int		xlen	= ssysparam->scbndv->xlen;
	int cplxReDofs  = ssysparam->scbndv->cplxDofs;
	
	/* memory allocation; */
	tvec<fftw_complex> rhoCplx = fn_tvec_init<fftw_complex> ( xlen * alloc_local_sys );

	/** 
	 * project to Fourier space;
	 *			 rhoCplx = rhoJCplx * d0JJ + d0bnd;
	 */
	int ind0, ind1;
	for ( int i = 0; i < xlen; i++ )
	{
		for ( int j0 = 0; j0 < alloc_local_sys; j0++ )
		{
			ind1 = i*alloc_local_sys + j0;
			fn_complex_setZero ( rhoCplx.val[ind1] );
			for ( int j1 = 0; j1 < nd; j1++ )
			{
				ind0 = j1*alloc_local_sys + j0;
				rhoCplx.val[ind1][0] += rhoJCplx.val[ind0][0] * ssysparam->sGJPv->d0JJ.val[i][j1];
				rhoCplx.val[ind1][1] += rhoJCplx.val[ind0][1] * ssysparam->sGJPv->d0JJ.val[i][j1];
			}
			rhoCplx.val[ind1][0] += ssysparam->scbndv->d0bnd.val[ind1][0];
			rhoCplx.val[ind1][1] += ssysparam->scbndv->d0bnd.val[ind1][1];
		}
	}

	return rhoCplx;
}


/**
 * \brief	Mass conservation by using 'fn_obt_mass' and 'fn_main_mass' many
 *			times to reach the desired precision;
 *
 * 'rhoJCplx':		given variable;
 * 'desir_mass':	the desired mass of the system;
 * 'TOL':			the maximal allowable tollerance;
 *
 * 'rhoJCplxMass':	return variable whose mass reaches 'desir_mass';
 * 'mass':			return the system mass;
 */
void fn_reach_mass		(	stu_system_param		*ssysparam,
							tvec<fftw_complex>		rhoJCplx,
							tvec<fftw_complex>		rhoJCplxMass,
							fftw_complex			mass,
							double					desir_mass,
							double					TOL	)
{
	/* copy rhoJCplx; */
	tvec<fftw_complex> rhoJCplxTmp  = fn_tvec_init<fftw_complex> ( rhoJCplx.len );
	memcpy ( rhoJCplxMass.val, rhoJCplx.val, sizeof(fftw_complex) * rhoJCplx.len );

	/* calculate the system mass; */
	fn_obt_mass ( ssysparam, rhoJCplxMass, mass );

	double err;
   	err = fn_complex_abs(mass) - desir_mass;
//	printf("mass: %.4e, err: %.4e\n", fn_complex_abs(mass), err);
	while ( fabs(err) > TOL )
	{
		mass[0] += desir_mass;
		fn_maintain_mass ( ssysparam, rhoJCplxMass, rhoJCplxTmp, mass );
		memcpy( rhoJCplxMass.val, rhoJCplxTmp.val, sizeof(fftw_complex) * rhoJCplx.len );
		fn_obt_mass ( ssysparam, rhoJCplxMass, mass );
		err = fn_complex_abs(mass) - desir_mass;
//		printf("mass: %.4e, err: %.4e\n", fn_complex_abs(mass), err);
	}

	/* release memory; */
	fn_tvec_free<fftw_complex> ( rhoJCplxTmp );
}


/**
 * \brief	Calculate the integral of 'rhoJCplx' to check whether the system is
 *			mass conservation;
 *
 * 'rhoJCplx':		given variable;
 *
 * 'mass':			the system mass;
 */
void fn_obt_mass			(	stu_system_param		*ssysparam,
								tvec<fftw_complex>		rhoJCplx,
								fftw_complex			mass )
{
	/* parameters; */
	int		nd		= ssysparam->scbndv->nd;
	int		xlen	= ssysparam->scbndv->xlen;
	int cplxReDofs  = ssysparam->scbndv->cplxDofs;
	
	/** 
	 * calculate the integral value;
	 */
	int ind0, ind1;
	fn_complex_setZero ( mass );
	for ( int i = 0; i < xlen; i++ )
	{
		/**
		 * rhoCplx = rhoJCplx * d0JJ + d0bnd;
		 *		only consider the first element along Fourier direction
		 *		since we are calculating the integral value;
		 */
		fftw_complex rhoCplx;
		fn_complex_setZero ( rhoCplx );
		for ( int j1 = 0; j1 < nd; j1++ )
		{
			ind0 = j1*alloc_local_sys;	// the first element along Fourier direction;
			rhoCplx[0] += rhoJCplx.val[ind0][0] * ssysparam->sGJPv->d0JJ.val[i][j1];
			rhoCplx[1] += rhoJCplx.val[ind0][1] * ssysparam->sGJPv->d0JJ.val[i][j1];
		}
		ind1 = i*alloc_local_sys;		// the first element along Fourier direction;
		rhoCplx[0] += ssysparam->scbndv->d0bnd.val[ind1][0];
		rhoCplx[1] += ssysparam->scbndv->d0bnd.val[ind1][1];
		/* multiply weight; */
		mass[0] += rhoCplx[0] * ssysparam->sGJPv->w.val[i];
		mass[1] += rhoCplx[1] * ssysparam->sGJPv->w.val[i];
	}
	mass[0] *= 0.5;
	mass[1] *= 0.5;

	/* obtain all information from different processes; */
	fftw_complex massRecv;
	massRecv[0] = mass[0];
	massRecv[1] = mass[1];
	MPI_Reduce ( &mass[0], &massRecv[0], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
	MPI_Reduce ( &mass[1], &massRecv[1], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
	MPI_Barrier ( MPI_COMM_WORLD );
	MPI_Bcast ( &mass[0], 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
	MPI_Bcast ( &mass[1], 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
}


/**
 * \brief	Maintain the system mass;
 *
 * 'rhoJCplx':		given variable;
 * 'mass':			the system mass;
 *
 * 'rhoJCplxMass':	return variable;
 */
void fn_maintain_mass		(	stu_system_param		*ssysparam,
								tvec<fftw_complex>		rhoJCplx,
								tvec<fftw_complex>		rhoJCplxMass,
								fftw_complex			mass	)
{
	/* parameters; */
	int		nd		= ssysparam->scbndv->nd;
	int		xlen	= ssysparam->scbndv->xlen;
	int cplxReDofs  = ssysparam->scbndv->cplxDofs;
	
	/* memory allocation; */
	tvec<fftw_complex> rhoCplxIgn	= fn_tvec_init<fftw_complex> ( xlen );
	tvec<fftw_complex> rhs			= fn_tvec_init<fftw_complex> ( nd );
	tvec<fftw_complex> rslt			= fn_tvec_init<fftw_complex> ( nd );

	/* copy rhoJCplx; */
	memcpy ( rhoJCplxMass.val, rhoJCplx.val, sizeof(fftw_complex) * rhoJCplx.len );

	/** 
	 * project to Fourier space but ignoring 'd0bnd';
	 *			 rhoCplxIgn(x,0) = d0JJ * rhoJCplx(x,0);
	 * meanwhile, subtract mass; rhoCplxIgn(x,0) - mass = rhoJCplxMass * d0JJ;
	 */
	int status;
	int ind0 = 0, ind1 = 0;
	/* rhoCplxIgn(x,0) - mass; */
	for ( int j = 0; j < xlen; j++ )
	{
		fn_complex_setZero ( rhoCplxIgn.val[j] );
		for ( int j0 = 0; j0 < nd; j0++ )
		{
			ind0 = j0 * alloc_local_sys;	// rhoJCplx(x,0);
			rhoCplxIgn.val[j][0] += ssysparam->sGJPv->d0JJ.val[j][j0] * rhoJCplx.val[ind0][0];
			rhoCplxIgn.val[j][1] += ssysparam->sGJPv->d0JJ.val[j][j0] * rhoJCplx.val[ind0][1];
		}
		rhoCplxIgn.val[j][0] -= mass[0];
		rhoCplxIgn.val[j][1] -= mass[1];
	}

	/* calculate the right term; */
	for ( int j = 0; j < nd; j++ )
	{
		fn_complex_setZero ( rhs.val[j] );
		for ( int j0 = 0; j0 < xlen; j0++ )
		{
			rhs.val[j][0] += ssysparam->sGJPv->d0JJ.val[j0][j] * 
							ssysparam->sGJPv->w.val[j0] * rhoCplxIgn.val[j0][0];
			rhs.val[j][1] += ssysparam->sGJPv->d0JJ.val[j0][j] * 
							ssysparam->sGJPv->w.val[j0] * rhoCplxIgn.val[j0][1];
		}
	}

	/* project onto the space of GJPs;
	 *	solve (u,J) = (J,J) \hat{u}; \hat{u} is coefficients;
	 */
	status = fn_umfpack_complex_solver ( ssysparam->sGJPv->innSMatd0JJ, rhs, rslt );
	if ( myrank == 0 )
		for (int j0 = 0; j0 < nd; j0++)
		{
			ind1 = j0 * alloc_local_sys; // the first element along Fourier direction;
			rhoJCplxMass.val[ind1][0] = rslt.val[j0][0];
			rhoJCplxMass.val[ind1][1] = rslt.val[j0][1];
		}

//
//	/* innSMatd0JJ * rhoJCplxMass; */
//	tCCSmat<double> innSMatd0JJTrans = fn_tCCSmat_init<double> (
//			ssysparam->sGJPv->innSMatd0JJ.col, 
//			ssysparam->sGJPv->innSMatd0JJ.row,
//			ssysparam->sGJPv->innSMatd0JJ.nnz );
//	tvec<fftw_complex> rhsCheck	 = fn_tvec_init<fftw_complex> ( nd );
//	innSMatd0JJTrans = fn_fast_trans_dCCSmat ( ssysparam->sGJPv->innSMatd0JJ );
//	fn_cvec_multiply_dCCSmat ( innSMatd0JJTrans, rslt, rhsCheck );
//	/* calculate error; */
//	double err = 100.0;
//	for ( int i = 0; i < nd; i++ )
//	{
//		fftw_complex rhsErrVec;
//		rhsErrVec[0] = rhsCheck.val[i][0] - rhs.val[i][0];
//		rhsErrVec[1] = rhsCheck.val[i][1] - rhs.val[i][1];
//		double rhsErr = fn_complex_abs ( rhsErrVec );
//		err = err < rhsErr ? err : rhsErr;
//	}
//	if ( myrank == 0 )
//		printf("err = %.10e\n", err);
//	fn_tCCSmat_free<double>		( innSMatd0JJTrans );
//	fn_tvec_free<fftw_complex>	( rhsCheck );
//

	/* release memory; */
	fn_tvec_free<fftw_complex> ( rhoCplxIgn );
	fn_tvec_free<fftw_complex> ( rhs  );
	fn_tvec_free<fftw_complex> ( rslt );
}


/**
 * \brief	memory allocation for FFTW in the calculation of the interface system;
 */
void fn_system_fftw_memory_allocation ( stu_system_param	*ssysparam )
{
	if ( myrank == 0 ) printf("Memory allocation for FFTW.\n");

	int	xlen	 = ssysparam->sGJPv->xlen;
	int cplxDofs = ssysparam->cplxReDofs;
	int dimCpt	 = ssysparam->dimReCpt;
	int	dimPhy	 = ssysparam->dimRePhy;
	int	Four_num = ssysparam->Four_num;
//	ssysparam->NCpt	= fn_tvec_init<ptrdiff_t> ( dimCpt );
//	for ( int i = 0; i < dimCpt; i++ )
//		ssysparam->NCpt.val[i] = Four_num;
	
	int	nd		= ssysparam->sGJPv->nd;
	int	rhsLen	= alloc_local_sys * nd;
	ssysparam->iter_rho_rhs			= fn_tvec_init<fftw_complex> ( rhsLen );
	ssysparam->iter_rho_rhs.row		= alloc_local_sys;
	ssysparam->iter_rho_rhs.col		= nd;
	ssysparam->iter_entropy_rhs		= fn_tvec_init<fftw_complex> ( rhsLen );
	ssysparam->iter_entropy_rhs.row = alloc_local_sys;
	ssysparam->iter_entropy_rhs.col = nd;
	ssysparam->grad_err			= fn_tvec_init<fftw_complex> ( rhsLen );
	ssysparam->grad_err.row		= alloc_local_sys;
	ssysparam->grad_err.col		= nd;
	fn_tvec_setZero_complex ( ssysparam->grad_err );
	ssysparam->hessian			= fn_tvec_init<fftw_complex> ( xlen*alloc_local_sys );
	ssysparam->hessian.row		= xlen;
	ssysparam->hessian.col		= alloc_local_sys;

//	ssysparam->sfftv->indKspace = fn_tmat_init<int>			 ( alloc_local_sys, dimCpt );
//	ssysparam->sfftv->projPlane = fn_tmat_init<double>		 ( alloc_local_sys, dimPhy );
//	ssysparam->sfftv->Gsquare   = fn_tvec_init<double>		 ( cplxDofs );
	ssysparam->sfftv->rhoCplx	= fn_tvec_init<fftw_complex> ( alloc_local_sys );
	ssysparam->sfftv->rhoReal	= fn_tvec_init<fftw_complex> ( alloc_local_sys );
	ssysparam->sfftv->fftw_Ctmp = fn_tvec_init<fftw_complex> ( alloc_local_sys );
	ssysparam->sfftv->fftw_Rtmp = fn_tvec_init<fftw_complex> ( alloc_local_sys );
	ssysparam->sfftv->cplxTmp   = fn_tvec_init<fftw_complex> ( alloc_local_sys );
	ssysparam->sfftv->quadTerm  = fn_tvec_init<fftw_complex> ( alloc_local_sys );
	ssysparam->sfftv->cubTerm   = fn_tvec_init<fftw_complex> ( alloc_local_sys );
	ssysparam->sfftv->quarTerm  = fn_tvec_init<fftw_complex> ( alloc_local_sys );
	ssysparam->sfftv->gradient  = fn_tvec_init<fftw_complex> ( xlen*alloc_local_sys );
	ssysparam->sfftv->gradient.row = xlen;
	ssysparam->sfftv->gradient.col = alloc_local_sys;

	/* prepare FFTW plan; */
	int offset = 0;
	char strNCpt[FILELEN-100];
	offset += sprintf(strNCpt, "%d_", dimCpt);
	for ( int i = 0; i < dimCpt; i++ )
	{
		offset += sprintf(strNCpt + offset, "%ld_", ssysparam->NCpt.val[i]);
	}
	offset += sprintf(strNCpt + offset, "%d.txt", flags);

	char fnr2c[FILELEN], fnc2r[FILELEN];
	sprintf(fnr2c, "plan/sys_forward_%s",  strNCpt);
	sprintf(fnc2r, "plan/sys_backward_%s", strNCpt);

	if ( access(fnr2c, F_OK) == 0 )
	{
		if ( myrank == 0 )
			printf("%s exist.\n", fnr2c);
		fftw_import_wisdom_from_filename( fnr2c );
	}
	if ( access(fnc2r, F_OK) == 0 )
	{
		if ( myrank == 0 )
			printf("%s exist.\n", fnc2r);
		fftw_import_wisdom_from_filename( fnc2r );
	}

	ssysparam->sfftv->Planc2cFord = fftw_mpi_plan_dft( dimCpt, ssysparam->NCpt.val, 
			ssysparam->sfftv->fftw_Rtmp.val, ssysparam->sfftv->fftw_Ctmp.val, 
			MPI_COMM_WORLD, FFTW_FORWARD,  flags );  // real to cplx
	ssysparam->sfftv->Planc2cBack = fftw_mpi_plan_dft( dimCpt, ssysparam->NCpt.val, 
			ssysparam->sfftv->fftw_Ctmp.val, ssysparam->sfftv->fftw_Rtmp.val, 
			MPI_COMM_WORLD, FFTW_BACKWARD, flags );  // cplx to real 

	if ( access(fnr2c, F_OK) != 0 )
	{
		fftw_export_wisdom_to_filename( fnr2c );
		if ( myrank == 0 )
			printf("%s save.\n", fnr2c);
	}
	if ( access(fnc2r, F_OK) != 0 )
	{
		fftw_export_wisdom_to_filename( fnc2r );
		if ( myrank == 0 )
			printf("%s save.\n", fnc2r);
	}
}


/**
 * \brief	Calculate the free energy and gradient term of the interface system;
 *
 * \param	energy:		three parts:
 *						interaction potential energy;
 *						entropy energy;
 *						hamilton energy;
 */
void fn_calc_system_energy		(	stu_system_param		*ssysparam,
									tvec<fftw_complex>		rhoJCplx,
									tvec<double>			energy,
										bool				withGrad,
										bool				isTest )
{
	/* parameters; */
	int		nd		= ssysparam->scbndv->nd;
	int		xlen	= ssysparam->scbndv->xlen;
	int cplxReDofs  = ssysparam->cplxReDofs;
	int	dimReCpt	= ssysparam->dimReCpt;
	
	/* memory allocation; */
	tvec<fftw_complex> rhoPotenCplx		= fn_tvec_init<fftw_complex> ( alloc_local_sys );
	tvec<fftw_complex> rhoPotenCplxQuad = fn_tvec_init<fftw_complex> ( alloc_local_sys );

	/* reset energy; */
	for ( int i = 0; i < 3; i++ )
		energy.val[i] = 0.0;

	/** 
	 * project to Fourier space;
	 *					rhoCplx = rhoJCplx * d0JJ + d0bnd;
	 * \partial_{x}^{2} rhoCplx = rhoJCplx * d2JJ + d2bnd;
	 * \partial_{x}^{4} rhoCplx = rhoJCplx * d4JJ + d4bnd;
	 * \partial_{x}^{6} rhoCplx = rhoJCplx * d6JJ + d6bnd;
	 *
	 * calculate free energies and entropy gradient;
	 */
	if ( strcmp(model_type, "LB") == 0 )
	{
		int ind0, ind1, ind2;
		double pk2;
		fftw_complex rhod0Cplx, rhod2Cplx;
		for ( int i = 0; i < xlen; i++ )
		{
			for ( int j0 = 0; j0 < alloc_local_sys; j0++ )
			{
				ind1 = i*alloc_local_sys + j0;
				/* calculate rhoCplx; */
				fn_complex_setZero ( rhod0Cplx );
				for ( int j1 = 0; j1 < nd; j1++ )
				{
					ind0 = j1*alloc_local_sys + j0;
					rhod0Cplx[0] += rhoJCplx.val[ind0][0] * ssysparam->sGJPv->d0JJ.val[i][j1];
					rhod0Cplx[1] += rhoJCplx.val[ind0][1] * ssysparam->sGJPv->d0JJ.val[i][j1];
				}
				rhod0Cplx[0] += ssysparam->scbndv->d0bnd.val[ind1][0];
				rhod0Cplx[1] += ssysparam->scbndv->d0bnd.val[ind1][1];
				/* calculate \partial_{x}^{2} rhoCplx; */
				fn_complex_setZero ( rhod2Cplx );
				for ( int j1 = 0; j1 < nd; j1++ )
				{
					ind0 = j1*alloc_local_sys + j0;
					rhod2Cplx[0] += rhoJCplx.val[ind0][0] * ssysparam->sGJPv->d2JJ.val[i][j1];
					rhod2Cplx[1] += rhoJCplx.val[ind0][1] * ssysparam->sGJPv->d2JJ.val[i][j1];
				}
				rhod2Cplx[0] += ssysparam->scbndv->d2bnd.val[ind1][0];
				rhod2Cplx[1] += ssysparam->scbndv->d2bnd.val[ind1][1];
				/* calculate (\lap+1) rho; */	
				ind2 = displs_sys[myrank] + j0;
				pk2 = ssysparam->sfftv->Gsquare.val[ind2];	// |R'PB k|^2;
				rhoPotenCplx.val[j0][0] = rhod2Cplx[0] + (1-pk2)*rhod0Cplx[0];
				rhoPotenCplx.val[j0][1] = rhod2Cplx[1] + (1-pk2)*rhod0Cplx[1];
				/* save rhoCplx for FFT; */
				ssysparam->sfftv->rhoCplx.val[j0][0] = rhod0Cplx[0];
				ssysparam->sfftv->rhoCplx.val[j0][1] = rhod0Cplx[1];
			}
			/* convolution calculation; */
			fn_convolution(rhoPotenCplx, rhoPotenCplxQuad, ssysparam->sfftv, displs_sys,
							alloc_local_sys, cplxReDofs, 2);
			fn_convolution(ssysparam->sfftv->rhoCplx, ssysparam->sfftv->quadTerm, 
							ssysparam->sfftv, displs_sys, alloc_local_sys, cplxReDofs, 2);
			fn_convolution(ssysparam->sfftv->rhoCplx, ssysparam->sfftv->cubTerm, 
							ssysparam->sfftv, displs_sys, alloc_local_sys, cplxReDofs, 3);
			fn_convolution(ssysparam->sfftv->rhoCplx, ssysparam->sfftv->quarTerm, 
							ssysparam->sfftv, displs_sys, alloc_local_sys, cplxReDofs, 4);

			if ( myrank == 0 )
			{
				/* energy calculation along Fourier directions; */
				/* interaction potential energy; */
				energy.val[0] += 0.5 * ssysparam->model_xi * rhoPotenCplxQuad.val[0][0] * 
									ssysparam->sGJPv->w.val[i];
				double entropyTmp = 
					-ssysparam->sfftv->quadTerm.val[0][0]	* ssysparam->model_tau/2.0 - 
					ssysparam->sfftv->cubTerm.val[0][0]		* ssysparam->model_gamma/3.0 + 
					ssysparam->sfftv->quarTerm.val[0][0]	* ssysparam->model_kappa/4.0; 
				energy.val[1] += entropyTmp * ssysparam->sGJPv->w.val[i];
			}
			/* entropy gradient; */
			if ( withGrad )
			{
				for ( int j0 = 0; j0 < alloc_local_sys; j0++ )
				{
					ind1 = i*alloc_local_sys + j0;
					ssysparam->sfftv->gradient.val[ind1][0] = 
						-ssysparam->model_tau	* ssysparam->sfftv->rhoCplx.val[j0][0] - 
						ssysparam->model_gamma	* ssysparam->sfftv->quadTerm.val[j0][0] +
						ssysparam->model_kappa	* ssysparam->sfftv->cubTerm.val[j0][0];
					ssysparam->sfftv->gradient.val[ind1][1] = 
						-ssysparam->model_tau	* ssysparam->sfftv->rhoCplx.val[j0][1] - 
						ssysparam->model_gamma	* ssysparam->sfftv->quadTerm.val[j0][1] +
						ssysparam->model_kappa	* ssysparam->sfftv->cubTerm.val[j0][1];
				}
				if ( isTest && i < 10 )
				{
					char fwFile[FILELEN];
					sprintf(fwFile, "%s/rank%d/rhoCplx%d.dat", rsltDir, myrank, i);
					fn_tvec_save_complex ( ssysparam->sfftv->rhoCplx, fwFile );
					sprintf(fwFile, "%s/rank%d/quadTerm%d.dat", rsltDir, myrank, i);
					fn_tvec_save_complex ( ssysparam->sfftv->quadTerm, fwFile );
				}
			}
		}
	}
	else if ( strcmp(model_type, "LP") == 0 )
	{
		double q = ssysparam->scale_val[1];
		double q2 = pow(q, 2);
		double coeffTmp = 1.0 + q2;
		int ind0, ind1, ind2;
		double pk2, pk4;
		fftw_complex rhod0Cplx, rhod2Cplx, rhod4Cplx;
		for ( int i = 0; i < xlen; i++ )
		{
			for ( int j0 = 0; j0 < alloc_local_sys; j0++ )
			{
				ind1 = i*alloc_local_sys + j0;
				/* calculate rhoCplx; */
				fn_complex_setZero ( rhod0Cplx );
				for ( int j1 = 0; j1 < nd; j1++ )
				{
					ind0 = j1*alloc_local_sys + j0;
					rhod0Cplx[0] += rhoJCplx.val[ind0][0] * ssysparam->sGJPv->d0JJ.val[i][j1];
					rhod0Cplx[1] += rhoJCplx.val[ind0][1] * ssysparam->sGJPv->d0JJ.val[i][j1];
				}
				rhod0Cplx[0] += ssysparam->scbndv->d0bnd.val[ind1][0];
				rhod0Cplx[1] += ssysparam->scbndv->d0bnd.val[ind1][1];
				/* calculate \partial_{x}^{2} rhoCplx; */
				fn_complex_setZero ( rhod2Cplx );
				for ( int j1 = 0; j1 < nd; j1++ )
				{
					ind0 = j1*alloc_local_sys + j0;
					rhod2Cplx[0] += rhoJCplx.val[ind0][0] * ssysparam->sGJPv->d2JJ.val[i][j1];
					rhod2Cplx[1] += rhoJCplx.val[ind0][1] * ssysparam->sGJPv->d2JJ.val[i][j1];
				}
				rhod2Cplx[0] += ssysparam->scbndv->d2bnd.val[ind1][0];
				rhod2Cplx[1] += ssysparam->scbndv->d2bnd.val[ind1][1];
				/* calculate \partial_{x}^{4} rhoCplx; */
				fn_complex_setZero ( rhod4Cplx );
				for ( int j1 = 0; j1 < nd; j1++ )
				{
					ind0 = j1*alloc_local_sys + j0;
					rhod4Cplx[0] += rhoJCplx.val[ind0][0] * ssysparam->sGJPv->d4JJ.val[i][j1];
					rhod4Cplx[1] += rhoJCplx.val[ind0][1] * ssysparam->sGJPv->d4JJ.val[i][j1];
				}
				rhod4Cplx[0] += ssysparam->scbndv->d4bnd.val[ind1][0];
				rhod4Cplx[1] += ssysparam->scbndv->d4bnd.val[ind1][1];
				/* calculate (\lap+1) (\lap+q^2) rho; */	
				ind2 = displs_sys[myrank] + j0;
				pk2 = ssysparam->sfftv->Gsquare.val[ind2];	// |R'PB k|^2;
				pk4 = pow(pk2, 2);
				rhoPotenCplx.val[j0][0] = rhod4Cplx[0] + (coeffTmp-2.0*pk2)*rhod2Cplx[0] +
					(q2-coeffTmp*pk2+pk4)*rhod0Cplx[0];
				rhoPotenCplx.val[j0][1] = rhod4Cplx[1] + (coeffTmp-2.0*pk2)*rhod2Cplx[1] +
					(q2-coeffTmp*pk2+pk4)*rhod0Cplx[1];
				/* save rhoCplx for FFT; */
				ssysparam->sfftv->rhoCplx.val[j0][0] = rhod0Cplx[0];
				ssysparam->sfftv->rhoCplx.val[j0][1] = rhod0Cplx[1];
			}
			/* convolution calculation; */
			fn_convolution(rhoPotenCplx, rhoPotenCplxQuad, ssysparam->sfftv, displs_sys, 
							alloc_local_sys, cplxReDofs, 2);
			fn_convolution(ssysparam->sfftv->rhoCplx, ssysparam->sfftv->quadTerm, 
							ssysparam->sfftv, displs_sys, alloc_local_sys, cplxReDofs, 2);
			fn_convolution(ssysparam->sfftv->rhoCplx, ssysparam->sfftv->cubTerm, 
							ssysparam->sfftv, displs_sys, alloc_local_sys, cplxReDofs, 3);
			fn_convolution(ssysparam->sfftv->rhoCplx, ssysparam->sfftv->quarTerm, 
							ssysparam->sfftv, displs_sys, alloc_local_sys, cplxReDofs, 4);

			if ( myrank == 0 )
			{
				/* energy calculation along Fourier directions; */
				/* interaction potential energy; */
				energy.val[0] += 0.5 * ssysparam->model_xi * rhoPotenCplxQuad.val[0][0] * 
									ssysparam->sGJPv->w.val[i];
				double entropyTmp = 
					-ssysparam->sfftv->quadTerm.val[0][0]	* ssysparam->model_tau/2.0 - 
					ssysparam->sfftv->cubTerm.val[0][0]		* ssysparam->model_gamma/3.0 + 
					ssysparam->sfftv->quarTerm.val[0][0]	* ssysparam->model_kappa/4.0; 
				energy.val[1] += entropyTmp * ssysparam->sGJPv->w.val[i];
			}
			/* entropy gradient; */
			if ( withGrad )
			{
				for ( int j0 = 0; j0 < alloc_local_sys; j0++ )
				{
					ind1 = i*alloc_local_sys + j0;
					ssysparam->sfftv->gradient.val[ind1][0] = 
						-ssysparam->model_tau	* ssysparam->sfftv->rhoCplx.val[j0][0] - 
						ssysparam->model_gamma	* ssysparam->sfftv->quadTerm.val[j0][0] +
						ssysparam->model_kappa	* ssysparam->sfftv->cubTerm.val[j0][0];
					ssysparam->sfftv->gradient.val[ind1][1] = 
						-ssysparam->model_tau	* ssysparam->sfftv->rhoCplx.val[j0][1] - 
						ssysparam->model_gamma	* ssysparam->sfftv->quadTerm.val[j0][1] +
						ssysparam->model_kappa	* ssysparam->sfftv->cubTerm.val[j0][1];
				}
			}
		}
	}
	else
	{
		if ( myrank == 0 )
		{
			printf("Error use 'fn_calc_system_energy'\n");
			printf("'model_type' %s out of consideration.\n", model_type);
		}
	}

	/* total energy; */
	if ( myrank == 0 )
	{
		energy.val[0] *= 0.5;
		energy.val[1] *= 0.5;
		energy.val[2] = energy.val[0] + energy.val[1];
	}
	MPI_Barrier ( MPI_COMM_WORLD );
	MPI_Bcast ( &energy.val[0], 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
	MPI_Bcast ( &energy.val[1], 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
	MPI_Bcast ( &energy.val[2], 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );

	/* release memory; */
	fn_tvec_free<fftw_complex> ( rhoPotenCplx );
	fn_tvec_free<fftw_complex> ( rhoPotenCplxQuad );
}


/**
 * \brief	Calculate the right term for iteration;
 *
 * \note!	'step_size' should be divided by 'iter_rho_rhs';
 */
void fn_obt_iter_rhs (	stu_system_param		*ssysparam,
						tvec<fftw_complex>		rhoJCplx,
						tvec<fftw_complex>		gradient,
							int					iterator,
							bool				isTest )
{
	int		nd		 = ssysparam->sGJPv->nd;
	int		xlen	 = ssysparam->sGJPv->xlen;
	int	cplxReDofs	 = ssysparam->cplxReDofs;
	int		rhsLen	 = nd * alloc_local_sys;
	double model_xi  = ssysparam->model_xi;

	/* memory allocation; */
	tvec<fftw_complex> rhoTmp1 = fn_tvec_init<fftw_complex> ( xlen*alloc_local_sys );
	tvec<fftw_complex> rhoTmp2 = fn_tvec_init<fftw_complex> ( xlen*alloc_local_sys );

	/* calculate the right term; */
	/**
	 * Jpoly.' * diag(w) * (rhoTmp1 - rhoTmp2);
	 *		rhoTmp1 = Jpoly * rhoJCplx / step_size
	 *		rhoTmp2 = gradient + model_xi*interact_bnd_grad;
	 *	Jpoly:				[xlen, nd];
	 *	rhoJCplx:			[nd, alloc_local_sys];	has been straightened;
	 *	gradient:			[xlen, alloc_local_sys];
	 *	interact_bnd_grad:	[xlen, alloc_local_sys];
	 */
	/* Jpoly*rhoJCplx/step_size - gradient - model_xi*interact_bnd_grad; */
	int ind0, ind1;
	for ( int j0 = 0; j0 < alloc_local_sys; j0++ )
	{
		for ( int i = 0; i < xlen; i++ )
		{
			ind1 = i*alloc_local_sys + j0;
			/* rhoTmp1 = Jpoly * rhoJCplx; */
			fn_complex_setZero ( rhoTmp1.val[ind1] ); // complex zero;
			for ( int j1 = 0; j1 < nd; j1++ )
			{
				ind0 = j1*alloc_local_sys + j0;
				rhoTmp1.val[ind1][0] += rhoJCplx.val[ind0][0] * 
											ssysparam->sGJPv->d0JJ.val[i][j1];
				rhoTmp1.val[ind1][1] += rhoJCplx.val[ind0][1] * 
											ssysparam->sGJPv->d0JJ.val[i][j1];
			}
			/* rhoTmp2 = gradient + model_xi*interact_bnd_grad; */
			rhoTmp2.val[ind1][0] = gradient.val[ind1][0] + 
					model_xi * ssysparam->interact_bnd_grad.val[ind1][0];
			rhoTmp2.val[ind1][1] = gradient.val[ind1][1] + 
					model_xi * ssysparam->interact_bnd_grad.val[ind1][1];
		}
	}
	/* Jpoly.' * diag(w) * rhoTmp1; */
	fn_tvec_setZero_complex ( ssysparam->iter_rho_rhs );
	fn_tvec_setZero_complex ( ssysparam->iter_entropy_rhs );
	for ( int j0 = 0; j0 < alloc_local_sys; j0++ )
	{
		for ( int i = 0; i < nd; i++ )
		{
			ind1 = j0*nd + i;
			for ( int j1 = 0; j1 < xlen; j1++ )
			{
				ind0 = j1*alloc_local_sys + j0;
				/* rhoTmp1; */
				ssysparam->iter_rho_rhs.val[ind1][0] += ssysparam->sGJPv->d0JJ.val[j1][i] *
						ssysparam->sGJPv->w.val[j1] * rhoTmp1.val[ind0][0];
				ssysparam->iter_rho_rhs.val[ind1][1] += ssysparam->sGJPv->d0JJ.val[j1][i] *
						ssysparam->sGJPv->w.val[j1] * rhoTmp1.val[ind0][1];
				/* rhoTmp2; */
				ssysparam->iter_entropy_rhs.val[ind1][0] += ssysparam->sGJPv->d0JJ.val[j1][i] *
						ssysparam->sGJPv->w.val[j1] * rhoTmp2.val[ind0][0];
				ssysparam->iter_entropy_rhs.val[ind1][1] += ssysparam->sGJPv->d0JJ.val[j1][i] *
						ssysparam->sGJPv->w.val[j1] * rhoTmp2.val[ind0][1];
			}
		}
	}

	if ( isTest )
	{
		char fwFile[FILELEN];
		sprintf(fwFile, "%s/rank%d/gradient%d.dat", rsltDir, myrank, iterator);
		fn_tvec_save_complex ( gradient, fwFile );
		sprintf(fwFile, "%s/rank%d/interact_bnd_grad%d.dat", rsltDir, myrank, iterator);
		fn_tvec_save_complex ( ssysparam->interact_bnd_grad, fwFile );
		/* save 'rhoTmp1'; */
		sprintf(fwFile, "%s/rank%d/iter_rhoTmp1%d.dat", rsltDir, myrank, iterator);
		fn_tvec_save_complex ( rhoTmp1, fwFile );
		/* save 'rhoTmp2'; */
		sprintf(fwFile, "%s/rank%d/iter_rhoTmp2%d.dat", rsltDir, myrank, iterator);
		fn_tvec_save_complex ( rhoTmp2, fwFile );
		/* save 'iter_rho_rhs'; */
		sprintf(fwFile, "%s/rank%d/iter_rho_rhs%d.dat", rsltDir, myrank, iterator);
		ssysparam->iter_rho_rhs.row = alloc_local_sys;
		ssysparam->iter_rho_rhs.col = nd;
		fn_tvec_save_complex ( ssysparam->iter_rho_rhs, fwFile );
		/* save 'iter_entropy_rhs'; */
		sprintf(fwFile, "%s/rank%d/iter_entropy_rhs%d.dat", rsltDir, myrank, iterator);
		ssysparam->iter_entropy_rhs.row = alloc_local_sys;
		ssysparam->iter_entropy_rhs.col = nd;
		fn_tvec_save_complex ( ssysparam->iter_entropy_rhs, fwFile );
	}

	/* release memory; */
	fn_tvec_free<fftw_complex> ( rhoTmp1 );
	fn_tvec_free<fftw_complex> ( rhoTmp2 );
}


/**
 * \brief	Obtain 'iter_rhs_all' containing all informations of 'iter_rhs';
 *
 */
void fn_obt_iter_rhs_all	(	stu_system_param		*ssysparam,
								tvec<fftw_complex>		iter_rhs,
								tvec<fftw_complex>		iter_rhs_all )
{
	/* load parameters; */
	int		nd		 = ssysparam->sGJPv->nd;
	int		xlen	 = ssysparam->sGJPv->xlen;
	int	cplxReDofs	 = ssysparam->cplxReDofs;

	/* offset and receive count; */
    int *displs    = (int *) malloc ( sizeof(int) * nprocs ); // offset of each process;
    int *recvCount = (int *) malloc ( sizeof(int) * nprocs ); // offset of each process;
	int sendcount;
	sendcount = 2 * nd * alloc_local_sys;
	MPI_Allgather ( &sendcount, 1, MPI_INT, recvCount, 1, MPI_INT, MPI_COMM_WORLD );
	displs[0] = 0;
	for ( int i = 1; i < nprocs; i++ )
		displs[i] = displs[i-1] + recvCount[i-1];
	
	/* collect all informations; */
	MPI_Allgatherv ( iter_rhs.val, sendcount, MPI_DOUBLE, 
				iter_rhs_all.val, recvCount, displs, 
				MPI_DOUBLE, MPI_COMM_WORLD );
	
	/* release memory; */
	free ( displs );
	free ( recvCount );
}


/**
 * \brief	Calculate the rhoJCplx weighted difference between two adjacent steps;
 *
 * \return	rhoJCplxDiff = (rhoJCplxNew - rhoJCplx) * d0JJ * w;
 *			rhoJCplxDiff: a 'tvec' vector with 'cplxReDofs' elements;
 */
void fn_diff_weight_rhoJCplx	(	stu_system_param		*ssysparam,
									tvec<fftw_complex>		rhoJCplxNew,
									tvec<fftw_complex>		rhoJCplx,
									tvec<fftw_complex>		rhoJCplxDiff )
{
	/* parameters; */
	int		nd		= ssysparam->scbndv->nd;
	int		xlen	= ssysparam->scbndv->xlen;
	int cplxReDofs  = ssysparam->scbndv->cplxDofs;

	/* temporary variable for w.' * d0JJ; */
	tvec<double> wd0JJ = fn_tvec_init<double> ( nd );
	fn_tvec_setZero<double> ( wd0JJ );
	for ( int i = 0; i < xlen; i++ )
		for ( int j1 = 0; j1 < nd; j1++ )
			wd0JJ.val[j1] += ssysparam->sGJPv->w.val[i] * 
							ssysparam->sGJPv->d0JJ.val[i][j1];

	/** 
	 * project to Fourier space;
	 *		rhoJCplxNew	 = rhoJCplxNew * d0JJ + d0bnd;
	 *		rhoJCplx	 = rhoJCplx	 * d0JJ + d0bnd;
	 *		rhoJCplxDiff = ( rhoJCplxNew - rhoJCplx ) * w
	 *					 = ( rhoJCplxNew - rhoJCplx ) * d0JJ * w;
	 */
	int ind0, ind1;
	fftw_complex tmp;
	fn_tvec_setZero_complex ( rhoJCplxDiff );
	for ( int j0 = 0; j0 < alloc_local_sys; j0++ )
	{
		for ( int j1 = 0; j1 < nd; j1++ )
		{
			ind0 = j1*alloc_local_sys + j0;
			tmp[0] = rhoJCplxNew.val[ind0][0] - rhoJCplx.val[ind0][0];
			tmp[1] = rhoJCplxNew.val[ind0][1] - rhoJCplx.val[ind0][1];
			rhoJCplxDiff.val[j0][0] += tmp[0] * wd0JJ.val[j1];
			rhoJCplxDiff.val[j0][1] += tmp[1] * wd0JJ.val[j1];
		}
	}
	fn_tvec_free<double> ( wd0JJ );
}


/**
 * \brief	Calculate the gradient weighted difference between two adjacent steps;
 *
 * \return	gradientW = (gradient - gradientOld) * w;
 *			gradientW: a 'tvec' vector with 'cplxReDofs' elements;
 */
void fn_diff_weight_gradient	(	stu_system_param		*ssysparam,
									tvec<fftw_complex>		gradientNew,
									tvec<fftw_complex>		gradient,
									tvec<fftw_complex>		gradientDiff )
{
	/* parameters; */
	int		nd		= ssysparam->scbndv->nd;
	int		xlen	= ssysparam->scbndv->xlen;
	int cplxReDofs  = ssysparam->scbndv->cplxDofs;

	/* gradientDiff = ( gradientNew - gradient ) w; */
	int ind0, ind1;
	fftw_complex tmp;
	fn_tvec_setZero_complex ( gradientDiff );
	for ( int j0 = 0; j0 < alloc_local_sys; j0++ )
	{
		for ( int j1 = 0; j1 < xlen; j1++ )
		{
			ind0 = j1*alloc_local_sys + j0;
			tmp[0] = gradientNew.val[ind0][0] - gradient.val[ind0][0];
			tmp[1] = gradientNew.val[ind0][1] - gradient.val[ind0][1];
			gradientDiff.val[j0][0] += tmp[0] * ssysparam->sGJPv->w.val[j1];
			gradientDiff.val[j0][1] += tmp[1] * ssysparam->sGJPv->w.val[j1];
		}
	}
}


/**
 * \brief	Calculate the Hessian matrix of the interface system;
 */
void fn_calc_system_hessian		(	stu_system_param		*ssysparam,
									tvec<fftw_complex>		rhoJCplx,
										int					iterator,
										bool				isTest )
{
	/* parameters; */
	int		nd		= ssysparam->scbndv->nd;
	int		xlen	= ssysparam->scbndv->xlen;
	int cplxReDofs  = ssysparam->cplxReDofs;
	int	dimReCpt	= ssysparam->dimReCpt;
	
	/** 
	 * project to Fourier space;
	 *		rhoCplx = rhoJCplx * d0JJ + d0bnd;
	 *
	 * calculate Hessian matrix;
	 */
	int ind0, ind1;
	fftw_complex rhod0Cplx;
	for ( int i = 0; i < xlen; i++ )
	{
		for ( int j0 = 0; j0 < alloc_local_sys; j0++ )
		{
			ind1 = i*alloc_local_sys + j0;
			/* calculate rhoCplx; */
			fn_complex_setZero ( rhod0Cplx );
			for ( int j1 = 0; j1 < nd; j1++ )
			{
				ind0 = j1*alloc_local_sys + j0;
				rhod0Cplx[0] += rhoJCplx.val[ind0][0] * ssysparam->sGJPv->d0JJ.val[i][j1];
				rhod0Cplx[1] += rhoJCplx.val[ind0][1] * ssysparam->sGJPv->d0JJ.val[i][j1];
			}
			rhod0Cplx[0] += ssysparam->scbndv->d0bnd.val[ind1][0];
			rhod0Cplx[1] += ssysparam->scbndv->d0bnd.val[ind1][1];
			/* save rhoCplx for FFT; */
			ssysparam->sfftv->rhoCplx.val[j0][0] = rhod0Cplx[0];
			ssysparam->sfftv->rhoCplx.val[j0][1] = rhod0Cplx[1];
		}
		/* convolution calculation; */
		fn_convolution(ssysparam->sfftv->rhoCplx, ssysparam->sfftv->quadTerm, 
						ssysparam->sfftv, displs_sys, alloc_local_sys, cplxReDofs, 2);
		/* Hessian matrix; */
		for ( int j0 = 0; j0 < alloc_local_sys; j0++ )
		{
			ind1 = i*alloc_local_sys + j0;
			ssysparam->hessian.val[ind1][0] = -ssysparam->model_tau - 
				2.0*ssysparam->model_gamma	* ssysparam->sfftv->rhoCplx.val[j0][0] +
				3.0*ssysparam->model_kappa	* ssysparam->sfftv->quadTerm.val[j0][0];
			ssysparam->hessian.val[ind1][1] = -ssysparam->model_tau	- 
				2.0*ssysparam->model_gamma	* ssysparam->sfftv->rhoCplx.val[j0][1] +
				3.0*ssysparam->model_kappa	* ssysparam->sfftv->quadTerm.val[j0][1];
		}
	}

	if ( isTest )
	{
		char fwFile[FILELEN];
		sprintf(fwFile, "%s/rank%d/hessian%d.dat", rsltDir, myrank, iterator);
		fn_tvec_save_complex ( ssysparam->hessian, fwFile );
	}
}


/**
 * \brief	Calculate the regularized Hessian matrix of the interface system;
 *
 * \note!	(\tilde{J}_{k} + \mu_{k}I) d_{k} = -\tilde{g}_{k};
 */
void fn_calc_system_reg_hessian	(	stu_system_param		*ssysparam,
									tvec<fftw_complex>		direction,
									tvec<fftw_complex>		regHessian,
									tCCSmat<double>			interact_grad_trans,
										double				mu )
{
	/* parameters; */
	int		nd		= ssysparam->scbndv->nd;
	int		xlen	= ssysparam->scbndv->xlen;
	int cplxReDofs  = ssysparam->cplxReDofs;
	int	dimReCpt	= ssysparam->dimReCpt;
	int	rhsLen		= alloc_local_sys * nd;
	double model_xi	= ssysparam->model_xi; // has been squared;

	/* temporary variables; */
	tvec<fftw_complex> tmp0		= fn_tvec_init<fftw_complex> ( alloc_local_sys );
	tvec<fftw_complex> tmp1		= fn_tvec_init<fftw_complex> ( alloc_local_sys );
	tvec<fftw_complex> convTmp	= fn_tvec_init<fftw_complex> ( alloc_local_sys );
	tvec<fftw_complex> interact_direction = fn_tvec_init<fftw_complex> ( rhsLen );

	/* calculate interact_grad * direction; */
	fn_cvec_multiply_dCCSmat ( interact_grad_trans, direction, interact_direction );
	
	/* calculate regularized Hessian matrix; */
	int ind0, ind1;
	for ( int i = 0; i < nd; i++ )
	{
		for ( int j0 = 0; j0 < alloc_local_sys; j0++ )
		{
			ind0 = i*alloc_local_sys + j0;
			ind1 = j0*nd + i;
			/* choose elements; */
			tmp0.val[j0][0] = ssysparam->hessian.val[ind0][0];
			tmp0.val[j0][1] = ssysparam->hessian.val[ind0][1];
			tmp1.val[j0][0] = direction.val[ind1][0];
			tmp1.val[j0][1] = direction.val[ind1][1];
		}
		/* convolution calculation of 'hessian' and 'direction'; */
		fn_convolution_general ( tmp0, tmp1, convTmp, ssysparam->sfftv, cplxReDofs );
		/* regularized Hessian matrix; */
		for ( int j0 = 0; j0 < alloc_local_sys; j0++ )
		{
			//ind1 = i*alloc_local_sys + j0;
			ind1 = j0*nd + i;
			regHessian.val[ind1][0] = model_xi * interact_direction.val[ind1][0] +
				mu * direction.val[ind1][0] + convTmp.val[j0][0];
			regHessian.val[ind1][1] = model_xi * interact_direction.val[ind1][1] +
				mu * direction.val[ind1][1] + convTmp.val[j0][1];
		}
	}

	/* release memory; */
	fn_tvec_free<fftw_complex> ( tmp0 );
	fn_tvec_free<fftw_complex> ( tmp1 );
	fn_tvec_free<fftw_complex> ( convTmp );
	fn_tvec_free<fftw_complex> ( interact_direction );
}


/**
 * \brief	pre-conditioner;
 *
 * \param	preconditioner:		has multiplied direction;
 */
void fn_calc_system_preconditioner ( stu_system_param		*ssysparam,
									tvec<fftw_complex>		direction,
									tvec<fftw_complex>		preconditioner,
										double				mu,
										double				delta,
										int					iterator,
										bool				isTest )
{
	/* parameters; */
	int		nd		= ssysparam->scbndv->nd;
	int		xlen	= ssysparam->scbndv->xlen;
	int cplxReDofs  = ssysparam->cplxReDofs;
	int	dimReCpt	= ssysparam->dimReCpt;
	int	rhsLen		= alloc_local_sys * nd;
	double model_xi	= ssysparam->model_xi; // has been squared;

	/* calculate preconditioner: Z (D(i,i) + delta I + mu I)^{-1} Z; */
	double coeff = mu + ssysparam->pcg_delta_coeff * delta;
	tCCSmat<double> precondMatrix;
	if ( ssysparam->pcg_type == 0 )
	{
		/* add 'coeff' to diagonal elements; */
		precondMatrix = fn_diag_add_dCCSmat ( 
						ssysparam->interact_grad, model_xi, coeff );
	}
	else if ( ssysparam->pcg_type == 1 )
	{
		/* only save diagonal elements and add 'coeff'; */
		precondMatrix = fn_obt_diag_add_dCCSmat ( 
						ssysparam->interact_grad, model_xi, coeff );
	}
	else if ( ssysparam->pcg_type == 2 )
	{
		/* add 'coeff' to all nonzero elements; */
		int row = ssysparam->interact_grad.row;
		int col = ssysparam->interact_grad.col;
		int nnz = ssysparam->interact_grad.nnz;
		precondMatrix = fn_tCCSmat_init<double> ( row, col, nnz );
		memcpy ( precondMatrix.JA,  ssysparam->interact_grad.JA,  sizeof(int) * (col+1) );
		memcpy ( precondMatrix.IA,  ssysparam->interact_grad.IA,  sizeof(int) * nnz );
		for ( int i = 0; i < nnz; i++ )
			precondMatrix.val[i] = model_xi * ssysparam->interact_grad.val[i] + coeff;
	}
	else
	{
		if ( myrank == 0 )
		{
			printf("Error use 'fn_calc_system_preconditioner'.\n");
			printf("Error 'pcg_type'.\n");
		}
	}

	/*  preconditioner \ direction; */
	int status = fn_umfpack_complex_solver ( precondMatrix, direction, preconditioner );

	if ( isTest )
	{
		char fwFile[FILELEN];
		/* precondMatrix; */
		sprintf(fwFile, "%s/rank%d/PCG_precondMatrix%d.dat", rsltDir, myrank, iterator);
		fn_tCCSmat_save<double> ( precondMatrix, fwFile );	
		/* preconditioner; */
		sprintf(fwFile, "%s/rank%d/PCG_preconditioner%d.dat", rsltDir, myrank, iterator);
		preconditioner.row = cplxReDofs;
		preconditioner.col = nd;
		fn_tvec_save_complex ( preconditioner, fwFile );
	}

	/* release memory; */
	fn_tCCSmat_free<double>		( precondMatrix );
}


/**
 * \brief	Calculate Newton direction directly;
 */
void fn_calc_system_newton_direction	(	stu_system_param		*ssysparam,
											tvec<fftw_complex>		direction,
											tCCSmat<double>			interact_grad_trans )
{
	/* parameters; */
	int		nd		= ssysparam->scbndv->nd;
	int		xlen	= ssysparam->scbndv->xlen;
	int cplxReDofs  = ssysparam->cplxReDofs;
	int	dimReCpt	= ssysparam->dimReCpt;
	int	rhsLen		= alloc_local_sys * nd;
	double model_xi	= ssysparam->model_xi; // has been squared;

	/* memory allocation; */
	tvec<fftw_complex> rhs  = fn_tvec_init<fftw_complex> ( rhsLen );

	/* obtain the diagonal elements; */
	tCCSmat<double>	precondMatrix = fn_obt_diag_add_dCCSmat ( 
								ssysparam->interact_grad, model_xi, 0.0 );
	for ( int i = 0; i < precondMatrix.nnz; i++ )
		precondMatrix.val[i] = 1.0 / precondMatrix.val[i];

	/* precondMatrix * interact_grad; */
	tCCSmat<double> preMgrad = fn_multiply_diag_dCCSmat ( 
								precondMatrix, ssysparam->interact_grad );

	/* -g; */
	fn_cvec_multiply_dCCSmat ( precondMatrix, ssysparam->grad_err, rhs );
	fn_tvec_constMultiply_complex (	rhs, -1.0 );

	/*  preconditioner \ direction; */
	int status = fn_umfpack_complex_solver ( preMgrad, rhs, direction );

	/* release memory; */
	fn_tvec_free<fftw_complex>	( rhs );
	fn_tCCSmat_free<double>		( preMgrad );
}


/**
 * \brief	Preconditioned conjugate gradient (PCG) method;
 *
 * \ref		Jiang, Si, Chen, and Bao, SIAM J. Sci. Comput., vol. 42, No. 6, 2020;
 *			Page 1364;
 */
int fn_calc_system_PCG			(	stu_system_param		*ssysparam,
									tvec<fftw_complex>		x,
									tCCSmat<double>			interact_grad_trans,
										double				mu,
										double				delta,
										int					iterator,
										bool				isTest )
{
	/* parameters; */
	int		nd		= ssysparam->scbndv->nd;
	int		xlen	= ssysparam->scbndv->xlen;
	int cplxReDofs  = ssysparam->cplxReDofs;
	int	dimReCpt	= ssysparam->dimReCpt;
	int	rhsLen		= alloc_local_sys * nd;
	double model_xi	= ssysparam->model_xi; // has been squared;

	/* memory allocation for temporary variables; */
	tvec<fftw_complex> p  = fn_tvec_init<fftw_complex> ( rhsLen ); // p^{i}:
	tvec<fftw_complex> Ap = fn_tvec_init<fftw_complex> ( rhsLen ); // Ap^{i}:
	tvec<fftw_complex> r  = fn_tvec_init<fftw_complex> ( rhsLen ); // r^{i};
	tvec<fftw_complex> Mr = fn_tvec_init<fftw_complex> ( rhsLen ); // Mr^{i};

	/* initialization of PCG; */
	fn_tvec_setZero_complex ( x ); // x^{0} = 0;
	/* r_{0] = A x^{0} - b = -b;	b = -\tilde{g}_{k}; */
//	fn_calc_system_reg_hessian ( ssysparam, x, Ap, interact_grad_trans, mu ); // A x^{0};
//	memcpy ( r.val, Ap.val, sizeof(fftw_complex) * rhsLen );
//	fn_tvec_add_complex ( r, ssysparam->grad_err, 1.0, 1.0 ); // r_{0} = A x^{0} - b;
	memcpy ( r.val, ssysparam->grad_err.val, sizeof(fftw_complex) * rhsLen ); // r_{0} = -b;

	/* calculate the pre-conditioner; M^{-1} r_{0}; */
	fn_calc_system_preconditioner ( ssysparam, r, Mr, mu, delta, iterator, isTest );

	/* initialization of x; p_{0} = -M^{-1} r_{0} */
	memcpy ( p.val, Mr.val, sizeof(fftw_complex) * rhsLen );
	fn_tvec_constMultiply_complex (	p, -1.0 );

	/* \|r_{i}\|_{M^{-1}}^{2}; */
	fftw_complex rMNormOld, rMNorm;
	fn_tvec_dotMultiplySum_complex ( r, Mr, r.len, rMNorm );
	memcpy ( rMNormOld, rMNorm, sizeof(fftw_complex) );

	double r2Norm = fn_tvec_norm_complex ( Mr ) / sqrt(cplxReDofs); // \|r_{i}\|_{M^{-1}};
	double tol	  = 1.0e-4 * r2Norm;	// maximal tolerance error of PCG; \eta;
	fftw_complex pANorm;	// \|p_{i}\|_{A}^{2};
	double alpha;			// \alpha_{i+1};
	double beta;			// \beta_{i+1};
	double r2NormOld = r2Norm;
	int	subIter	   = 0;
	int print_step = ssysparam->pcg_print_step;

	if ( isTest && myrank == 0 )
	{
		printf("\n");
		printf("rMNorm = %.10e, %.10e.\n", rMNorm[0], rMNorm[1]);
		printf("r2Norm = %.10e.\n", r2Norm);
		printf("tol = %.10e.\n", tol);
		printf("\n");
	}

	/* PCG iteration; */
	while ( r2Norm > tol && subIter < ssysparam->pcg_iter_max )
	{
		/* regularized Hessian matrix; A p_{i} */
		fn_calc_system_reg_hessian ( ssysparam, p, Ap, interact_grad_trans, mu );

		/* \|p_{i}\|_{A}^{2}; */
		fn_tvec_dotMultiplySum_complex ( p, Ap, p.len, pANorm );

		/* \alpha_{i+1} = \|r_{i}\|_{M^{-1}}^{2} / \|p_{i}\|_{A}^{2}; */
		fftw_complex alphaCplx;
		fn_complex_divide( rMNormOld, pANorm, alphaCplx );
		alpha = alphaCplx[0];

		/* x^{i+1} = x^{i} + \alpha_{i+1} p_{i}; */
		fn_tvec_add_complex ( x,  p, 1.0, alpha );

		/* r_{i+1} = r_{i} + \alpha_{i+1} Ap_{i}; */
		fn_tvec_add_complex ( r, Ap, 1.0, alpha );

		/* M^{-1} r_{i+1}; */
		fn_calc_system_preconditioner ( ssysparam, r, Mr, mu, delta, subIter, isTest );

		/* \|r_{i+1}\|; */
		r2Norm = fn_tvec_norm_complex ( r ) / sqrt(cplxReDofs);

		/* \|r_{i}\|_{M^{-1}}^{2}; */
		fn_tvec_dotMultiplySum_complex ( r, Mr, r.len, rMNorm );

		/* beta_{i+1} = \|r_{i+1}\|_{M^{-1}}^{2} / \|r_{i}\|_{M^{-1}}^{2}; */
		fftw_complex betaCplx;
		fn_complex_divide( rMNorm, rMNormOld, betaCplx );
		beta = betaCplx[0];

		/* p_{i+1} = -M^{-1} r_{i+1} + \beta_{i+1} p_{i}; */
		fn_tvec_add_complex ( p, Mr, beta, -1.0 );

		if ( subIter%print_step == 0 && myrank == 0 )
			printf("subIter = %d \t r2Norm = %.10e.\n\n", subIter, r2Norm);

		if ( isTest && myrank == 0 )
		{
			printf("subIter = %d \t r2Norm = %.10e.\n", subIter, r2Norm);
			printf("rMNorm = %.5e, %.5e \t pANorm = %.5e, %.5e.\n", 
					rMNorm[0], rMNorm[1], pANorm[0], pANorm[1]);
			printf("alpha = %.5e \t beta = %.5e.\n\n", alpha, beta);
		}

		/* save 'rMNorm'; */
		memcpy ( rMNormOld, rMNorm, sizeof(fftw_complex) );

		subIter ++ ;
	}


	if ( isTest )
	{
		char fwFile[FILELEN];
		/* x; */
		sprintf(fwFile, "%s/rank%d/PCG_x%d.dat", rsltDir, myrank, iterator);
		x.row = cplxReDofs;
		x.col = nd;
		fn_tvec_save_complex ( x, fwFile );
		/* p; */
		sprintf(fwFile, "%s/rank%d/PCG_p%d.dat", rsltDir, myrank, iterator);
		p.row = cplxReDofs;
		p.col = nd;
		fn_tvec_save_complex ( p, fwFile );
		/* Ap; */
		sprintf(fwFile, "%s/rank%d/PCG_Ap%d.dat", rsltDir, myrank, iterator);
		Ap.row = cplxReDofs;
		Ap.col = nd;
		fn_tvec_save_complex ( Ap, fwFile );
		/* r; */
		sprintf(fwFile, "%s/rank%d/PCG_r%d.dat", rsltDir, myrank, iterator);
		r.row = cplxReDofs;
		r.col = nd;
		fn_tvec_save_complex ( r, fwFile );
		/* Mr; */
		sprintf(fwFile, "%s/rank%d/PCG_Mr%d.dat", rsltDir, myrank, iterator);
		Mr.row = cplxReDofs;
		Mr.col = nd;
		fn_tvec_save_complex ( Mr, fwFile );
	}

	/* release memory; */
	fn_tvec_free<fftw_complex> ( p  );
	fn_tvec_free<fftw_complex> ( Ap );
	fn_tvec_free<fftw_complex> ( r  );
	fn_tvec_free<fftw_complex> ( Mr );

	return subIter;
}


/**
 * \brief	recover 'rhoJCplx' after screening spectra;
 */
int  fn_recover_rhoJCplx	(	stu_system_param		*ssysparam,
								tvec<fftw_complex>		rhoJCplx )
{
	printf("Recover 'rhoJCplx'\n");

	/* parameters; */
	double	TOL		= ssysparam->recoverTOL;
	int		nd		= ssysparam->scbndv->nd;
	int		xlen	= ssysparam->scbndv->xlen;
	int cplxReDofs  = ssysparam->scbndv->cplxDofs;
	int rhsLen		= nd * alloc_local_sys;
	
	/* memory allocation; */
	tvec<fftw_complex> rhoCplx		= fn_tvec_init<fftw_complex> ( xlen * alloc_local_sys );
	tvec<fftw_complex> rhoRecover	= fn_tvec_init<fftw_complex> ( rhsLen );
	tvec<fftw_complex> rhoJCplxTmp	= fn_tvec_init<fftw_complex> ( rhsLen );
	
	/** 
	 * project to Fourier space and screen spectra by 'TOL';
	 *			 rhoCplx = rhoJCplx * d0JJ + d0bnd;
	 */
	int ind0, ind1;
	int count = 0; // count the number of tmpCplx with absolute value greater than TOL;
	fftw_complex tmpCplx;
	for ( int i = 0; i < xlen; i++ )
	{
		for ( int j0 = 0; j0 < alloc_local_sys; j0++ )
		{
			ind1 = i*alloc_local_sys + j0;
			fn_complex_setZero ( rhoCplx.val[ind1] );
			tmpCplx[0] = 0.0, tmpCplx[1] = 0.0;
			for ( int j1 = 0; j1 < nd; j1++ )
			{
				ind0 = j1*alloc_local_sys + j0;
				tmpCplx[0] += rhoJCplx.val[ind0][0] * ssysparam->sGJPv->d0JJ.val[i][j1];
				tmpCplx[1] += rhoJCplx.val[ind0][1] * ssysparam->sGJPv->d0JJ.val[i][j1];
			}
			tmpCplx[0] += ssysparam->scbndv->d0bnd.val[ind1][0];
			tmpCplx[1] += ssysparam->scbndv->d0bnd.val[ind1][1];
			if ( fn_complex_abs(tmpCplx) > TOL )
			{
				count ++;
				rhoCplx.val[ind1][0] = tmpCplx[0] - ssysparam->scbndv->d0bnd.val[ind1][0];
				rhoCplx.val[ind1][1] = tmpCplx[1] - ssysparam->scbndv->d0bnd.val[ind1][1];
			}
		}
	}
	int countAll = xlen * alloc_local_sys;
	printf("The ratio of spectra with intensities greater than %.4e: %d/%d = %.4e\n",
			TOL, count, countAll, (double) count / countAll);

	/* Jpoly.' * diag(w) * rhoTmp1; */
	for ( int j0 = 0; j0 < alloc_local_sys; j0++ )
	{
		for ( int i = 0; i < nd; i++ )
		{
			ind1 = j0*nd + i;
			rhoRecover.val[ind1][0] = 0.0;
			rhoRecover.val[ind1][1] = 0.0;
			for ( int j1 = 0; j1 < xlen; j1++ )
			{
				ind0 = j1*alloc_local_sys + j0;
				rhoRecover.val[ind1][0] += ssysparam->sGJPv->d0JJ.val[j1][i] *
						ssysparam->sGJPv->w.val[j1] * rhoCplx.val[ind0][0];
				rhoRecover.val[ind1][1] += ssysparam->sGJPv->d0JJ.val[j1][i] *
						ssysparam->sGJPv->w.val[j1] * rhoCplx.val[ind0][1];
			}
		}
	}

	/* innSMatd0JJ; */
	tCCSmat<double> G0innd0JJ = fn_tensor_diag_dCCSmat ( 
			ssysparam->sfftv->Gsquare, ssysparam->sGJPv->innSMatd0JJ, 0 );

	/* calculate the new 'rhoJCplx'; */
	int status = fn_umfpack_complex_solver ( G0innd0JJ, rhoRecover, rhoJCplx );

	/* transpose 'rhoJCplx'; */
	memcpy ( rhoJCplxTmp.val, rhoJCplx.val, sizeof(fftw_complex) * rhsLen );
	for ( int j0 = 0; j0 < alloc_local_sys; j0++ )
	{
		for ( int j1 = 0; j1 < nd; j1++ )
		{
			ind0 = j1*alloc_local_sys + j0;
			ind1 = j0*nd + j1;
			rhoJCplx.val[ind0][0] = rhoJCplxTmp.val[ind1][0];
			rhoJCplx.val[ind0][1] = rhoJCplxTmp.val[ind1][1];
		}
	}

	/* release memory; */
	fn_tvec_free<fftw_complex>	( rhoCplx	  );
	fn_tvec_free<fftw_complex>	( rhoRecover  );
	fn_tvec_free<fftw_complex>	( rhoJCplxTmp );
	fn_tCCSmat_free<double>		( G0innd0JJ	  );

	return status;
}


/* calculate the error between density profiles at step0 and step1; */
void fn_recover_error	( stu_system_param	*ssysparam,
								int			step0,
								int			step1	)
{
	char fname0[FILELEN];
	char fname1[FILELEN];
	sprintf(fname0, "%s/sys_density%d.dat", rsltDir, step0);
	sprintf(fname1, "%s/sys_density%d.dat", rsltDir, step1);

	int		xlen	= ssysparam->scbndv->xlen;
	int		y_num	= ssysparam->y_num;
	int		z_num	= ssysparam->z_num;

	int status = 1;
	if ( ssysparam->dimRePhy == 1 )
	{
		/* read data from 'fname0'; */
		tvec<double> x0		= fn_tvec_init<double> ( xlen );
		tvec<double> y0		= fn_tvec_init<double> ( y_num );
		tvec<double> data0	= fn_tvec_init<double> ( xlen * y_num );
		status = fn_read_data_2d ( fname0, x0, y0, data0 );

 		/* read data from 'fname1'; */
		tvec<double> x1		= fn_tvec_init<double> ( xlen );
		tvec<double> y1		= fn_tvec_init<double> ( y_num );
		tvec<double> data1	= fn_tvec_init<double> ( xlen * y_num );
		status = fn_read_data_2d ( fname1, x1, y1, data1 );

		/* calculate error; */
		fn_tvec_add<double> ( x0,	 x1,	1.0, -1.0 ); // x0 - x1;
		fn_tvec_add<double> ( y0,	 y1,	1.0, -1.0 ); // y0 - y1;
		fn_tvec_add<double> ( data0, data1, 1.0, -1.0 ); // data0 - data1;
		double err_x	= fn_tvec_maxAbs<double> ( x0 );
		double err_y	= fn_tvec_maxAbs<double> ( y0 );
		double err_data = fn_tvec_maxAbs<double> ( data0 );
		double err_norm = fn_tvec_norm2 <double> ( data0 );
		if ( myrank == 0 )
		{
			printf("\t Error between density profiles at step %d and %d.\n", step0, step1);
			printf("\t ---> x error = %.6e\n", err_x);
			printf("\t ---> y error = %.6e\n", err_y);
			printf("\t ---> density error: l^inf = %.6e, l^2 = %.6e\n", err_data, err_norm);
		}

		/* release memory; */
		fn_tvec_free<double> ( x0	 );
		fn_tvec_free<double> ( x1	 );
		fn_tvec_free<double> ( y0	 );
		fn_tvec_free<double> ( y1	 );
		fn_tvec_free<double> ( data0 );
		fn_tvec_free<double> ( data1 );
	} 
	else if ( ssysparam->dimRePhy == 2 )
	{
		/* read data from 'fname0'; */
		tvec<double> x0		= fn_tvec_init<double> ( xlen );
		tvec<double> y0		= fn_tvec_init<double> ( y_num );
		tvec<double> z0		= fn_tvec_init<double> ( z_num );
		tvec<double> data0	= fn_tvec_init<double> ( xlen * y_num * z_num );
		fn_read_data_3d ( fname0, x0, y0, z0, data0 );

 		/* read data from 'fname1'; */
		tvec<double> x1		= fn_tvec_init<double> ( xlen );
		tvec<double> y1		= fn_tvec_init<double> ( y_num );
		tvec<double> z1		= fn_tvec_init<double> ( z_num );
		tvec<double> data1	= fn_tvec_init<double> ( xlen * y_num * z_num );
		fn_read_data_3d ( fname1, x1, y1, z1, data1 );

		/* calculate error; */
		fn_tvec_add<double> ( x0,	 x1,	1.0, -1.0 ); // x0 - x1;
		fn_tvec_add<double> ( y0,	 y1,	1.0, -1.0 ); // y0 - y1;
		fn_tvec_add<double> ( z0,	 z1,	1.0, -1.0 ); // z0 - z1;
		fn_tvec_add<double> ( data0, data1, 1.0, -1.0 ); // data0 - data1;
		double err_x	= fn_tvec_maxAbs<double> ( x0 );
		double err_y	= fn_tvec_maxAbs<double> ( y0 );
		double err_z	= fn_tvec_maxAbs<double> ( z0 );
		double err_data = fn_tvec_maxAbs<double> ( data0 );
		double err_norm = fn_tvec_norm2 <double> ( data0 );
		if ( myrank == 0 )
		{
			printf("\t Error between density profiles at step %d and %d.\n", step0, step1);
			printf("\t ---> x error = %.6e\n", err_x);
			printf("\t ---> y error = %.6e\n", err_y);
			printf("\t ---> z error = %.6e\n", err_z);
			printf("\t ---> density error: l^inf = %.6e, l^2 = %.6e\n", err_data, err_norm);
		}

		/* release memory; */
		fn_tvec_free<double> ( x0	 );
		fn_tvec_free<double> ( x1	 );
		fn_tvec_free<double> ( y0	 );
		fn_tvec_free<double> ( y1	 );
		fn_tvec_free<double> ( z0	 );
		fn_tvec_free<double> ( z1	 );
		fn_tvec_free<double> ( data0 );
		fn_tvec_free<double> ( data1 );
	}
}


/**
 * \brief	The integration of functions about saving data;
 *			save Fourier k, densities, R'PBk;
 */
void fn_save_system_phase	(	stu_system_param		*ssysparam,
									int					iterator )
{	
	/* save rhoJCplx; */
	if ( ssysparam->save_type[0] == 'y' || ssysparam->save_type[0] == 'Y' )
	{
		char fwFile[FILELEN];
		sprintf(fwFile, "%s/rank%d/sys_rhoJCplx%d.dat", rsltDir, myrank, iterator);
		ssysparam->scbndv->rhoJCplx.row = ssysparam->sGJPv->nd;
		ssysparam->scbndv->rhoJCplx.col = ssysparam->cplxReDofs;
		fn_tvec_save_complex ( ssysparam->scbndv->rhoJCplx, fwFile );
	}
	/* save density calculating by rhoJCplx (in scbndv); */
	if ( ssysparam->save_type[1] == 'y' || ssysparam->save_type[1] == 'Y' )
	{
		fn_disp_system_density ( ssysparam->scbndv->rhoJCplx, ssysparam, iterator );
	}
	if ( ssysparam->save_type[2] == 'y' || ssysparam->save_type[2] == 'Y' )
	{
	}
}
