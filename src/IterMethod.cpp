/*! \file	IterMethod.cpp
 *
 * \brief	Iteration methods for the calculation of stable interface structure;
 */

#include "Data.h"
#include "Head.h"
#include "DataOperators.h"
#include "Mytimer.h"
#include "functs.h"

/**
 * \brief	Choose iteration method;
 */
double fn_choose_method		(	stu_system_param		*ssysparam,
									bool				isTest,
									int					iter_load	)
{
	if ( myrank == 0 )
	{
		printf(" <========== Iteration for the calculation of ");
		printf("stable interface structure ==========> \n\n");
	}
	mytimer_t timer;
	timer.reset();
	timer.start();
	tvec<double> energy = fn_tvec_init<double> ( 3 );
	fn_tvec_setZero<double> ( energy );
	int iterator = 0;

	/* parameters; */
	int		nd	   = ssysparam->sGJPv->nd;
	int cplxReDofs = ssysparam->cplxReDofs;	
	int	rhsLen	   = nd * alloc_local_sys;

	/* memory allocation; */
	fftw_complex mass;
	ssysparam->massVec = fn_tvec_init<fftw_complex> ( nd );

	/* preparation for mass conservation; */
	char massFile[FILELEN];
	sprintf(massFile, "%s/massVec.dat", rsltDir);
	if ( iter_load <= 0 ) // generate new mass vector;
	{
		/* mass constraint; */
		fn_obt_mass ( ssysparam, ssysparam->scbndv->rhoJCplx, mass );
		if ( myrank == 0 )
			printf("---> Mass: %.10e\n", fn_complex_abs(mass));

		/* force mass constraint; */
		/*
		tvec<fftw_complex> rhoJCplxMass	 = fn_tvec_init<fftw_complex> ( rhsLen );
		fn_obt_mass ( ssysparam, ssysparam->scbndv->rhoJCplx, mass );
		if ( myrank == 0 )
			printf("---> Mass before maintaining process: %.10e\n", fn_complex_abs(mass));
		fn_reach_mass ( ssysparam, ssysparam->scbndv->rhoJCplx, rhoJCplxMass, mass, 0.0, 1e-8 );
		memcpy( ssysparam->scbndv->rhoJCplx.val, rhoJCplxMass.val, sizeof(fftw_complex) * rhsLen );
		if ( myrank == 0 )
			printf("---> Mass after  maintaining process: %.10e\n", fn_complex_abs(mass));
		fn_tvec_free<fftw_complex>	( rhoJCplxMass	);
		*/

		/* save mass vector; */
		if ( myrank == 0 )
			for ( int j1 = 0; j1 < nd; j1++ )
			{
				int ind0 = j1*alloc_local_sys;	// the first element along Fourier direction;
				ssysparam->massVec.val[j1][0] = ssysparam->scbndv->rhoJCplx.val[ind0][0];
				ssysparam->massVec.val[j1][1] = ssysparam->scbndv->rhoJCplx.val[ind0][1];
			}

		/* save massVec; */
		if ( myrank == 0 )
			fn_tvec_save_complex ( ssysparam->massVec, massFile );
		MPI_Barrier ( MPI_COMM_WORLD );
	}

	/* load mass vector to avoid error; */
	FILE *fp = fopen(massFile, "r");
	if ( fp == NULL )	// invalid file;
	{
		if ( myrank == 0 )
		{
			printf("Error using 'fn_choose_method'\n");
			printf("Unable to read file '%s'. No such file or directory.\n", massFile);
		}
		return 0.0;
	}
	else
	{
		/* load 'massVec'; */
		int val, len;
		val = fscanf(fp, "%d", &(len));
		for ( int i = 0; i < nd; i++ )
		{
			val = fscanf(fp, "%lf", &(ssysparam->massVec.val[i][0]));
			val = fscanf(fp, "%lf", &(ssysparam->massVec.val[i][1]));
		}

		/* mass constraint; */
		if ( myrank == 0 )
			for ( int j1 = 0; j1 < nd; j1++ )
			{
				int ind0 = j1*alloc_local_sys;	// the first element along Fourier direction;
				ssysparam->scbndv->rhoJCplx.val[ind0][0] = ssysparam->massVec.val[j1][0]; 
				ssysparam->scbndv->rhoJCplx.val[ind0][1] = ssysparam->massVec.val[j1][1];
			}

		/* mass; */
		fn_obt_mass ( ssysparam, ssysparam->scbndv->rhoJCplx, mass );
		if ( myrank == 0 )
			printf("---> Mass with loading mass vector: %.10e\n", fn_complex_abs(mass));
	}
	fclose(fp);	
	
	/* choose iteration method; */
	if ( myrank == 0 )
		printf("\t Iteration method is %s.\n", ssysparam->iter_method);
	if ( strcmp(ssysparam->iter_method, "SIS") == 0 )
		iterator = fn_sis_method ( ssysparam, energy, isTest, iter_load );
	else if ( strcmp(ssysparam->iter_method, "rSIS") == 0 )
		iterator = fn_sis_rand_method ( ssysparam, energy, isTest, iter_load );
	else if ( strcmp(ssysparam->iter_method, "APG") == 0 )
		iterator = fn_apg_method ( ssysparam, energy, ssysparam->tol, isTest, iter_load );
	else if ( strcmp(ssysparam->iter_method, "nAPG") == 0 )
	{
		if ( myrank == 0 )
			printf("\n===> APG step: \n");
		iterator = fn_apg_method	( ssysparam, energy, ssysparam->newton_tol, isTest, iter_load );
		if ( myrank == 0 )
			printf("\n===> Newton step: \n");
		iterator = fn_newton_method ( ssysparam, energy, iterator, isTest );
	}
	else
	{
		if ( myrank == 0 )
		{
			printf("Error use 'fn_choose_method'\n");
			printf("Iteration method %s is out of consideration.\n", ssysparam->iter_method);
		}
		return 0.0;
	}

	/* save stationary value; */
	if ( ssysparam->iter_max > 0 )
		fn_save_system_phase (ssysparam, -1);

	double hamilton = energy.val[2];
	fn_tvec_free<double>		( energy );
	fn_tvec_free<fftw_complex>	( ssysparam->massVec );

	MPI_Barrier ( MPI_COMM_WORLD );
	timer.pause();
	if ( myrank == 0 )
		printf("\n\t ***** time cost of iteration: %f seconds *****\n\n", timer.get_current_time());

	return hamilton;
}


/**
 * \brief	Semi-implicit method;
 */
int fn_sis_method		(	stu_system_param		*ssysparam,
							tvec<double>			energy,
								bool				isTest,
								int					iter_load	)
{
	double	tol			= ssysparam->tol;
	double	tolham		= ssysparam->tolham;
	double	step_size	= ssysparam->step_size;
	int		iter_max	= ssysparam->iter_max;
	int		print_step	= ssysparam->print_step;
	int		save_step	= ssysparam->save_step;

	double	model_xi	= ssysparam->model_xi;
	double	model_tau   = ssysparam->model_tau;
	double	model_gamma = ssysparam->model_gamma;
	double	model_kappa = ssysparam->model_kappa;
	fftw_complex mass;

	/* transpose 'interact_grad'; */
	mytimer_t timer;
	timer.reset();
	timer.start();
	tCCSmat<double> interact_grad_trans = fn_fast_trans_dCCSmat ( ssysparam->interact_grad );

	MPI_Barrier ( MPI_COMM_WORLD );
	timer.pause();
	if ( myrank == 0 )
		printf("\n\t ***** time cost of transpose dCCSmat: %f seconds *****\n\n", 
				timer.get_current_time());

	/* parameters; */
	int		nd	   = ssysparam->sGJPv->nd;
	int cplxReDofs = ssysparam->cplxReDofs;	
	int	rhsLen	   = nd * alloc_local_sys;
	int dataLen;
	if ( myrank == nprocs-1 )
		dataLen = cplxReDofs - displs_sys[myrank];
	else
		dataLen = alloc_local_sys;

	/* offset and receive count; */
    int *displs    = (int *) malloc ( sizeof(int) * nprocs ); // offset of each process;
    int *recvCount = (int *) malloc ( sizeof(int) * nprocs ); // offset of each process;
	int sendcount;
	sendcount = nd * alloc_local_sys;
	MPI_Allgather ( &sendcount, 1, MPI_INT, recvCount, 1, MPI_INT, MPI_COMM_WORLD );
	displs[0] = 0;
	for ( int i = 1; i < nprocs; i++ )
		displs[i] = displs[i-1] + recvCount[i-1];
	int dataAmount = 0;
	for ( int i = 0; i < nprocs; i++ )
		dataAmount += recvCount[i];

	/* memory allocation; */
	tvec<fftw_complex> rhoJCplxNew	 = fn_tvec_init<fftw_complex> ( rhsLen );
	tvec<fftw_complex> rhoJCplxTrans = fn_tvec_init<fftw_complex> ( rhsLen );
	tvec<fftw_complex> iter_rhs		 = fn_tvec_init<fftw_complex> ( rhsLen );
	tvec<fftw_complex> rhoJCplxAll	 = fn_tvec_init<fftw_complex> ( dataAmount );
	tvec<fftw_complex> iter_rhs_all	 = fn_tvec_init<fftw_complex> ( dataAmount );
	tvec<fftw_complex> grad_err_all	 = fn_tvec_init<fftw_complex> ( dataAmount );
	rhoJCplxNew.row   = nd;		rhoJCplxNew.col   = alloc_local_sys;
	rhoJCplxTrans.col = nd;		rhoJCplxTrans.row = alloc_local_sys;
	rhoJCplxAll.row = nd;		rhoJCplxAll.col = cplxReDofs;

	/* obtain the iteration matrix; */
	fn_obt_iter_matrix ( ssysparam, step_size, -1, isTest );

	/* initialization; */
	int iterator = iter_load;
	int	status = 1;
	double res = 1.0;
	double hamilton, hamiltonOld, diffham;
	char dataName[FILELEN];
	sprintf(dataName, "%s/sys_energy_error.dat", rsltDir);

	int isDecay = 1;
	hamiltonOld = 100.0;
	diffham		= 1.0;

	/* data for reloading; */
	FILE *fdata, *tmp_fdata; 
	if ( iter_load > 0 )
	{
		/* load step size; */
		if ( access(dataName, F_OK) == 0 ) // exist;
		{
			/* memory allocation for existing data; */
			tvec<double> tmp_step	 = fn_tvec_init<double> ( iter_load );
			tvec<double> tmp_energy0 = fn_tvec_init<double> ( iter_load );
			tvec<double> tmp_energy1 = fn_tvec_init<double> ( iter_load );
			tvec<double> tmp_energy2 = fn_tvec_init<double> ( iter_load );
			tvec<double> tmp_res	 = fn_tvec_init<double> ( iter_load );
			tvec<int>	 tmp_isDacy	 = fn_tvec_init<int>	( iter_load );

			/* load data before 'iter_load'; */
			int val;
			tmp_fdata = fopen(dataName, "r");
			for ( int i = 0; i < iter_load; i++ )
			{
				val = fscanf(tmp_fdata, "%lf", &(tmp_step.val[i]));
				val = fscanf(tmp_fdata, "%lf", &(tmp_energy0.val[i]));
				val = fscanf(tmp_fdata, "%lf", &(tmp_energy1.val[i]));
				val = fscanf(tmp_fdata, "%lf", &(tmp_energy2.val[i]));
				val = fscanf(tmp_fdata, "%lf", &(tmp_res.val[i]));
				val = fscanf(tmp_fdata, "%d",  &(tmp_isDacy.val[i]));
			}
			fclose(tmp_fdata);
			if ( iter_load == 1 )
				diffham = -10.0;
			else
				diffham = tmp_energy2.val[iter_load-1] - tmp_energy2.val[iter_load-2];
			hamiltonOld = tmp_energy2.val[iter_load-1];

			/* rewrite the existing data; */
			fdata  = fopen(dataName, "w");
			if ( myrank == 0 )
			{
				for ( int i = 0; i < iter_load; i++ )
				{
					fprintf(fdata, "%.6E\t%+.15E\t%+.15E\t%+.20E\t%+.15E\t%d\n", 
							tmp_step.val[i], tmp_energy0.val[i], tmp_energy1.val[i], 
							tmp_energy2.val[i], tmp_res.val[i], tmp_isDacy.val[i]);
				}
			}

			/* release memory of temporay variables; */
			fn_tvec_free<double> ( tmp_step );
			fn_tvec_free<double> ( tmp_energy0 );
			fn_tvec_free<double> ( tmp_energy1 );
			fn_tvec_free<double> ( tmp_energy2 );
			fn_tvec_free<double> ( tmp_res );
			fn_tvec_free<int>	 ( tmp_isDacy );
		}
		else
		{
			fdata = fopen(dataName, "w");
		}
	}
	else
	{
		fdata = fopen(dataName, "w");
	}

	if ( myrank == 0 )
		printf("---> start iteration.\n");
	while ( iterator < iter_max )
	{
		/* calculate free energy and entropy gradient; */
		fn_calc_system_energy ( ssysparam, ssysparam->scbndv->rhoJCplx, energy, true, isTest );
		hamilton = energy.val[2];

		/* save data; */
		if ( myrank == 0 )
			fprintf(fdata, "%.6E\t%+.15E\t%+.15E\t%+.20E\t%+.15E\t%d\n", step_size,
						energy.val[0], energy.val[1], energy.val[2], res, isDecay);
		if ( iterator%print_step == 0 )
		{
			fn_obt_mass ( ssysparam, ssysparam->scbndv->rhoJCplx, mass );
			if ( myrank == 0 )
			{
				timer.pause();
				printf("\nIterator %d: step_size = %.6e \t res = %.10e", iterator, step_size, res);
				printf(" \t diffham = %.10e \t isDecay = %d\n", diffham, isDecay);
				printf("\t Laplace term: % .10e \t Nonlinear term: % .10e ", 
						energy.val[0], energy.val[1]);
				printf("\t hamilton = % .10e\n", energy.val[2]);
				printf("\t Mass: %.10e \t time cost: %f seconds\n", 
						fn_complex_abs(mass), timer.get_current_time());
			}
		}

		if ( (iterator > iter_load ) && (iterator % save_step) == 0 )
			fn_save_system_phase (ssysparam, iterator);		// save data;

		diffham = hamilton - hamiltonOld;
		hamiltonOld = hamilton;

		/* set flag to check if the hamilton energy is decay; */
		if ( diffham < 0 )
			isDecay = 1;
		else
			isDecay = 0;
	
		/* generate the right term for iteration; */
		fn_obt_iter_rhs ( ssysparam, ssysparam->scbndv->rhoJCplx,
						ssysparam->sfftv->gradient, iterator, isTest );
		memcpy ( iter_rhs.val, ssysparam->iter_rho_rhs.val, sizeof(fftw_complex) * rhsLen );
		fn_tvec_add_complex ( iter_rhs, ssysparam->iter_entropy_rhs, 1.0/step_size, -1.0 );
		fn_obt_iter_rhs_all ( ssysparam, iter_rhs, iter_rhs_all );

		/* calculate the new 'rhoJCplx'; */
		status = fn_umfpack_complex_solver ( ssysparam->iter_matrix, iter_rhs_all, rhoJCplxAll );
		/* mass conservation; */
		if ( ssysparam->massFlag == 1 )
			for ( int j1 = 0; j1 < nd; j1++ )
			{
				int ind0 = j1; // the first element along Fourier direction;
				rhoJCplxAll.val[ind0][0] = ssysparam->massVec.val[j1][0];
				rhoJCplxAll.val[ind0][1] = ssysparam->massVec.val[j1][1];
			}

		/* calculate gradient error; */
		fn_cvec_multiply_dCCSmat ( interact_grad_trans, rhoJCplxAll, grad_err_all );
		for ( int j0 = 0; j0 < dataLen; j0++ )
		{
			for ( int j1 = 0; j1 < nd; j1++ )
			{
				int ind0 = j0 * nd + j1;
				int ind1 = displs[myrank] + j0*nd + j1;
				ssysparam->grad_err.val[ind0][0] = grad_err_all.val[ind1][0];
				ssysparam->grad_err.val[ind0][1] = grad_err_all.val[ind1][1];
			}
		}
		fn_tvec_add_complex ( ssysparam->grad_err, ssysparam->iter_entropy_rhs, model_xi, 1.0 );
		if ( myrank == 0 )
			for ( int i = 0; i < nd; i++ )
				fn_complex_setZero ( ssysparam->grad_err.val[i] );
		res = fn_tvec_maxAbs_complex ( ssysparam->grad_err );

		/* update 'rhoJCplx'; transpose for consistency; */
		for ( int j0 = 0; j0 < alloc_local_sys; j0++ )
		{
			for ( int j1 = 0; j1 < nd; j1++ )
			{
				//int ind0 = j0 * nd + j1; // no transpose;
				int ind0 = j1 * alloc_local_sys + j0;
				int ind1 = displs[myrank] + j0*nd + j1;
				ssysparam->scbndv->rhoJCplx.val[ind0][0] = rhoJCplxAll.val[ind1][0];
				ssysparam->scbndv->rhoJCplx.val[ind0][1] = rhoJCplxAll.val[ind1][1];
			}
		}
		if ( isTest )
		{
			char fwFile[FILELEN];
			sprintf(fwFile, "%s/rank%d/iterMat%d.dat", rsltDir, myrank, iterator);
			fn_tCCSmat_save ( ssysparam->iter_matrix, fwFile );
			sprintf(fwFile, "%s/rank%d/iterRhsAll%d.dat", rsltDir, myrank, iterator);
			fn_tvec_save_complex ( iter_rhs_all, fwFile );
			sprintf(fwFile, "%s/rank%d/rhoJCplxAll%d.dat", rsltDir, myrank, iterator);
			fn_tvec_save_complex ( rhoJCplxAll, fwFile );
			sprintf(fwFile, "%s/rank%d/rhoJCplx%d.dat", rsltDir, myrank, iterator);
			fn_tvec_save_complex ( ssysparam->scbndv->rhoJCplx, fwFile );
			sprintf(fwFile, "%s/rank%d/gradErr%d.dat", rsltDir, myrank, iterator);
			fn_tvec_save_complex ( ssysparam->grad_err, fwFile );
		}

		iterator ++;
//		if ( (res < tol && iterator > iter_load+2) || fabs(diffham) < 1.0e-15 )
//		if ( res < tol && iterator > iter_load+2 )
//		if ( res < tol )
		if ( res < tol || fabs(diffham) < tolham )
			break;
	}

	/* output data about the stationary state; */
	fn_calc_system_energy ( ssysparam, ssysparam->scbndv->rhoJCplx, energy, false, isTest );
	hamilton = energy.val[2];
	if ( myrank == 0 )
		fprintf(fdata, "%.6E\t%+.15E\t%+.15E\t%+.20E\t%+.15E\t%d\n", step_size,
				energy.val[0], energy.val[1], energy.val[2], res, isDecay);
	if ( myrank == 0 )
	{
		timer.pause();
		printf("\n\t ***** the interface system\n");
		printf("\t ***** iterator = %d \t step size = %.6e \t time cost: %f seconds\n", 
				iterator, step_size, timer.get_current_time());
		printf("\t ***** res = % .10e \t diffham = % .10e\n", res, diffham);
		printf("\t ***** Laplace term: % .15e,\t Nonlinear term: % .15e\n", 
				energy.val[0], energy.val[1]);
		printf("\t ***** hamilton = % .20e\n\n", energy.val[2]);
	}
	fclose(fdata);

	/* Releases memory; */
	fn_tCCSmat_free<double>		( interact_grad_trans );
	fn_tvec_free<fftw_complex>	( rhoJCplxNew	);
	fn_tvec_free<fftw_complex>	( rhoJCplxTrans );
	fn_tvec_free<fftw_complex>	( iter_rhs		);
	fn_tvec_free<fftw_complex>	( rhoJCplxAll	);
	fn_tvec_free<fftw_complex>	( iter_rhs_all	);
	fn_tvec_free<fftw_complex>	( grad_err_all	);
	free ( displs );		free ( recvCount );

	return iterator;
}


/**
 * \brief	Semi-implicit method with random disturbance;
 */
int fn_sis_rand_method	(	stu_system_param		*ssysparam,
							tvec<double>			energy,
								bool				isTest,
								int					iter_load	)
{
	double	tol			= ssysparam->tol;
	double	tolham		= ssysparam->tolham;
	double	step_size	= ssysparam->step_size;
	int		iter_max	= ssysparam->iter_max;
	int		print_step	= ssysparam->print_step;
	int		save_step	= ssysparam->save_step;

	double	model_xi	= ssysparam->model_xi;
	double	model_tau   = ssysparam->model_tau;
	double	model_gamma = ssysparam->model_gamma;
	double	model_kappa = ssysparam->model_kappa;
	fftw_complex mass;

	double	rand_coeff	= 1e-4;

	/* transpose 'interact_grad'; */
	mytimer_t timer;
	timer.reset();
	timer.start();
	tCCSmat<double> interact_grad_trans = fn_fast_trans_dCCSmat ( ssysparam->interact_grad );

	MPI_Barrier ( MPI_COMM_WORLD );
	timer.pause();
	if ( myrank == 0 )
		printf("\n\t ***** time cost of transpose dCCSmat: %f seconds *****\n\n", 
				timer.get_current_time());

	/* parameters; */
	int		nd	   = ssysparam->sGJPv->nd;
	int cplxReDofs = ssysparam->cplxReDofs;	
	int	rhsLen	   = nd * alloc_local_sys;
	int dataLen;
	if ( myrank == nprocs-1 )
		dataLen = cplxReDofs - displs_sys[myrank];
	else
		dataLen = alloc_local_sys;

	/* offset and receive count; */
    int *displs    = (int *) malloc ( sizeof(int) * nprocs ); // offset of each process;
    int *recvCount = (int *) malloc ( sizeof(int) * nprocs ); // offset of each process;
	int sendcount;
	sendcount = nd * alloc_local_sys;
	MPI_Allgather ( &sendcount, 1, MPI_INT, recvCount, 1, MPI_INT, MPI_COMM_WORLD );
	displs[0] = 0;
	for ( int i = 1; i < nprocs; i++ )
		displs[i] = displs[i-1] + recvCount[i-1];
	int dataAmount = 0;
	for ( int i = 0; i < nprocs; i++ )
		dataAmount += recvCount[i];

	/* memory allocation; */
	tvec<fftw_complex> rhoJCplxNew	 = fn_tvec_init<fftw_complex> ( rhsLen );
	tvec<fftw_complex> rhoJCplxTrans = fn_tvec_init<fftw_complex> ( rhsLen );
	tvec<fftw_complex> iter_rhs		 = fn_tvec_init<fftw_complex> ( rhsLen );
	tvec<fftw_complex> rhoJCplxAll	 = fn_tvec_init<fftw_complex> ( dataAmount );
	tvec<fftw_complex> iter_rhs_all	 = fn_tvec_init<fftw_complex> ( dataAmount );
	tvec<fftw_complex> grad_err_all	 = fn_tvec_init<fftw_complex> ( dataAmount );
	rhoJCplxNew.row   = nd;		rhoJCplxNew.col   = alloc_local_sys;
	rhoJCplxTrans.col = nd;		rhoJCplxTrans.row = alloc_local_sys;
	rhoJCplxAll.row = nd;		rhoJCplxAll.col = cplxReDofs;

	/* read 'rhoJCplx'; */
	char fname[FILELEN];
	sprintf(fname, "%s/rank%d/sys_rhoJCplx%d.dat", "result/LP_HEX_HEX_45s71a", myrank, -1);
	if ( myrank == 0 )
		printf("Read file '%s'.\n", fname);
	ssysparam->scbndv->rhoJCplx = fn_tvec_read_complex ( fname );

	/* mass constraint; */
   	fn_obt_mass ( ssysparam, ssysparam->scbndv->rhoJCplx, mass );
	if ( myrank == 0 )
		printf("---> Mass: %.10e\n", fn_complex_abs(mass));

	/* save mass vector; */
	ssysparam->massVec = fn_tvec_init<fftw_complex> ( nd );
	if ( myrank == 0 )
		for ( int j1 = 0; j1 < nd; j1++ )
		{
			int ind0 = j1 * alloc_local_sys;	// the first element along Fourier direction;
			ssysparam->massVec.val[j1][0] = ssysparam->scbndv->rhoJCplx.val[ind0][0];
			ssysparam->massVec.val[j1][1] = ssysparam->scbndv->rhoJCplx.val[ind0][1];
		}

	/* add random disturbance; */
	tvec<fftw_complex> randomReal	 = fn_tvec_init<fftw_complex> ( alloc_local_sys );
	tvec<fftw_complex> randomCplx	 = fn_tvec_init<fftw_complex> ( alloc_local_sys );
	for ( int j1 = 0; j1 < nd; j1++ )
	{
		srand(time(NULL)); // random seed;
		for ( int j0 = 0; j0 < alloc_local_sys; j0++ )
		{
			randomReal.val[j0][0] = rand()/(RAND_MAX+1.0) * rand_coeff;
			randomReal.val[j0][1] = 0.0;
		}
		memcpy(ssysparam->sfftv->fftw_Rtmp.val, randomReal.val, sizeof(fftw_complex)*alloc_local_sys);
		fftw_execute(ssysparam->sfftv->Planc2cFord);  // real to cplx;
		fn_tvec_constMultiply_complex ( ssysparam->sfftv->fftw_Ctmp, 1.0/cplxReDofs );
		memcpy(randomCplx.val, ssysparam->sfftv->fftw_Ctmp.val, sizeof(fftw_complex)*alloc_local_sys);

		for ( int j0 = 0; j0 < alloc_local_sys; j0++ )
		{
			int ind1 = j0*nd + j1;
			ssysparam->scbndv->rhoJCplx.val[ind1][0] += randomCplx.val[j0][0];
			ssysparam->scbndv->rhoJCplx.val[ind1][1] += randomCplx.val[j0][1];
		}
	}
	fn_tvec_free<fftw_complex>	( randomReal );
	fn_tvec_free<fftw_complex>	( randomCplx );


	/* obtain the iteration matrix; */
	fn_obt_iter_matrix ( ssysparam, step_size, -1, isTest );

	/* initialization; */
	int iterator = iter_load;
	int	status = 1;
	double res = 1.0;
	double hamilton, hamiltonOld, diffham;
	char dataName[FILELEN];
	sprintf(dataName, "%s/sys_energy_error.dat", rsltDir);
	FILE *fdata = fopen(dataName, "w");

	int isDecay = 1;
	hamiltonOld = 100.0;
	diffham		= 1.0;

	if ( myrank == 0 )
		printf("---> start iteration.\n");
	while ( iterator < iter_max )
	{
		/* calculate free energy and entropy gradient; */
		fn_calc_system_energy ( ssysparam, ssysparam->scbndv->rhoJCplx, energy, true, isTest );
		hamilton = energy.val[2];

		/* save data; */
		if ( myrank == 0 )
			fprintf(fdata, "%.6E\t%+.15E\t%+.15E\t%+.20E\t%+.15E\t%d\n", step_size,
						energy.val[0], energy.val[1], energy.val[2], res, isDecay);
		if ( iterator%print_step == 0 )
		{
			fn_obt_mass ( ssysparam, ssysparam->scbndv->rhoJCplx, mass );
			if ( myrank == 0 )
			{
				timer.pause();
				printf("\nIterator %d: step_size = %.6e \t res = %.10e", iterator, step_size, res);
				printf(" \t diffham = %.10e \t isDecay = %d\n", diffham, isDecay);
				printf("\t Laplace term: % .10e \t Nonlinear term: % .10e ", 
						energy.val[0], energy.val[1]);
				printf("\t hamilton = % .10e\n", energy.val[2]);
				printf("\t Mass: %.10e \t time cost: %f seconds\n", 
						fn_complex_abs(mass), timer.get_current_time());
			}
		}

		if ( (iterator > iter_load ) && (iterator % save_step) == 0 )
			fn_save_system_phase (ssysparam, iterator);		// save data;

		diffham = hamilton - hamiltonOld;
		hamiltonOld = hamilton;

		/* set flag to check if the hamilton energy is decay; */
		if ( diffham < 0 )
			isDecay = 1;
		else
			isDecay = 0;

		/* generate the right term for iteration; */
		fn_obt_iter_rhs ( ssysparam, ssysparam->scbndv->rhoJCplx,
						ssysparam->sfftv->gradient, iterator, isTest );
		memcpy ( iter_rhs.val, ssysparam->iter_rho_rhs.val, sizeof(fftw_complex) * rhsLen );
		fn_tvec_add_complex ( iter_rhs, ssysparam->iter_entropy_rhs, 1.0/step_size, -1.0 );
		fn_obt_iter_rhs_all ( ssysparam, iter_rhs, iter_rhs_all );

		/* calculate the new 'rhoJCplx'; */
		status = fn_umfpack_complex_solver ( ssysparam->iter_matrix, iter_rhs, rhoJCplxAll );

		/* mass conservation; */
		if ( ssysparam->massFlag == 1 )
			for ( int j1 = 0; j1 < nd; j1++ )
			{
				int ind0 = j1; // the first element along Fourier direction;
				rhoJCplxAll.val[ind0][0] = ssysparam->massVec.val[j1][0];
				rhoJCplxAll.val[ind0][1] = ssysparam->massVec.val[j1][1];
			}

		/* calculate gradient error; */
		fn_cvec_multiply_dCCSmat ( interact_grad_trans, rhoJCplxAll, grad_err_all );
		for ( int j0 = 0; j0 < dataLen; j0++ )
		{
			for ( int j1 = 0; j1 < nd; j1++ )
			{
				int ind0 = j0 * nd + j1;
				int ind1 = displs[myrank] + j0*nd + j1;
				ssysparam->grad_err.val[ind0][0] = grad_err_all.val[ind1][0];
				ssysparam->grad_err.val[ind0][1] = grad_err_all.val[ind1][1];
			}
		}
		fn_tvec_add_complex ( ssysparam->grad_err, ssysparam->iter_entropy_rhs, model_xi, 1.0 );
		if ( myrank == 0 )
			for ( int i = 0; i < nd; i++ )
				fn_complex_setZero ( ssysparam->grad_err.val[i] );
		res = fn_tvec_maxAbs_complex ( ssysparam->grad_err );

		/* update 'rhoJCplx'; * transpose for consistency; */
		for ( int j0 = 0; j0 < alloc_local_sys; j0++ )
		{
			for ( int j1 = 0; j1 < nd; j1++ )
			{
				//int ind0 = j0 * nd + j1; // no transpose;
				int ind0 = j1 * alloc_local_sys + j0;
				int ind1 = displs[myrank] + j0*nd + j1;
				ssysparam->scbndv->rhoJCplx.val[ind0][0] = rhoJCplxAll.val[ind1][0];
				ssysparam->scbndv->rhoJCplx.val[ind0][1] = rhoJCplxAll.val[ind1][1];
			}
		}

		iterator ++;
//		if ( (res < tol && iterator > iter_load+2) || fabs(diffham) < 1.0e-15 )
//		if ( res < tol && iterator > iter_load+2 )
//		if ( res < tol )
		if ( res < tol || fabs(diffham) < tolham )
			break;
	}

	/* output data about the stationary state; */
	fn_calc_system_energy ( ssysparam, ssysparam->scbndv->rhoJCplx, energy, false, isTest );
	hamilton = energy.val[2];
	if ( myrank == 0 )
		fprintf(fdata, "%.6E\t%+.15E\t%+.15E\t%+.20E\t%+.15E\t%d\n", step_size,
				energy.val[0], energy.val[1], energy.val[2], res, isDecay);
	if ( myrank == 0 )
	{
		timer.pause();
		printf("\n\t ***** the interface system\n");
		printf("\t ***** iterator = %d \t step size = %.6e \t time cost: %f seconds\n", 
				iterator, step_size, timer.get_current_time());
		printf("\t ***** res = % .10e \t diffham = % .10e\n", res, diffham);
		printf("\t ***** Laplace term: % .15e,\t Nonlinear term: % .15e\n", 
				energy.val[0], energy.val[1]);
		printf("\t ***** hamilton = % .20e\n\n", energy.val[2]);
	}
	fclose(fdata);

	/* Releases memory; */
	fn_tCCSmat_free<double>		( interact_grad_trans );
	fn_tvec_free<fftw_complex>	( rhoJCplxNew	);
	fn_tvec_free<fftw_complex>	( rhoJCplxTrans );
	fn_tvec_free<fftw_complex>	( iter_rhs	);
	fn_tvec_free<fftw_complex>	( rhoJCplxAll	);
	fn_tvec_free<fftw_complex>	( iter_rhs_all	);
	fn_tvec_free<fftw_complex>	( grad_err_all	);
	free ( displs );		free ( recvCount );
	
	return iterator;
}


/**
 * \brief	Adaptive accelerated Bregman proximal gradient method (P2 case);
 *			Unpublished name is accelerated proximal gradient method (APG);
 */
int fn_apg_method		(	stu_system_param		*ssysparam,
							tvec<double>			energy,
								double				tol,
								bool				isTest,
								int					iter_load	)
{
	double	tolham		= ssysparam->tolham;
	int		bbType	    = 1;					// the type of BB stepping size;
	double	step0	    = ssysparam->step_size;	// initial step size;
	double	step_min    = ssysparam->step_min;	// the lower bound of step size;
	double	step_max    = ssysparam->step_max;	// the upper bound of step size;
	int		iter_max    = ssysparam->iter_max;
	int		print_step  = ssysparam->print_step;
	int		save_step   = ssysparam->save_step;

	double	model_xi	= ssysparam->model_xi;
	double	model_tau   = ssysparam->model_tau;
	double	model_gamma = ssysparam->model_gamma;
	double	model_kappa = ssysparam->model_kappa;
	fftw_complex mass;

	/* transpose 'interact_grad'; */
	mytimer_t timer;
	timer.reset();
	timer.start();
	tCCSmat<double> interact_grad_trans = fn_fast_trans_dCCSmat ( ssysparam->interact_grad );

	MPI_Barrier ( MPI_COMM_WORLD );
	timer.pause();
	if ( myrank == 0 )
		printf("\n\t ***** time cost of transpose dCCSmat: %f seconds *****\n\n", 
				timer.get_current_time());

	/* parameters; */
	int		nd	   = ssysparam->sGJPv->nd;
	int		xlen   = ssysparam->sGJPv->xlen;
	int cplxReDofs = ssysparam->cplxReDofs;	
	int	rhsLen	   = nd * alloc_local_sys;
	int dataLen;
	if ( myrank == nprocs-1 )
		dataLen = cplxReDofs - displs_sys[myrank];
	else
		dataLen = alloc_local_sys;

	/* offset and receive count; */
    int *displs    = (int *) malloc ( sizeof(int) * nprocs ); // offset of each process;
    int *recvCount = (int *) malloc ( sizeof(int) * nprocs ); // offset of each process;
	int sendcount;
	sendcount = nd * alloc_local_sys;
	MPI_Allgather ( &sendcount, 1, MPI_INT, recvCount, 1, MPI_INT, MPI_COMM_WORLD );
	displs[0] = 0;
	for ( int i = 1; i < nprocs; i++ )
		displs[i] = displs[i-1] + recvCount[i-1];
	int dataAmount = 0;
	for ( int i = 0; i < nprocs; i++ )
		dataAmount += recvCount[i];

	/* memory allocation; */
	tvec<fftw_complex> rhoJCplxNew   = fn_tvec_init<fftw_complex> ( rhsLen );
	tvec<fftw_complex> rhoJCplxTrans = fn_tvec_init<fftw_complex> ( rhsLen );
	tvec<fftw_complex> iter_rhs		 = fn_tvec_init<fftw_complex> ( rhsLen );
	tvec<fftw_complex> gradientOld	 = fn_tvec_init<fftw_complex> ( xlen * alloc_local_sys );
	tvec<fftw_complex> rhoJCplxAll	 = fn_tvec_init<fftw_complex> ( dataAmount );
	tvec<fftw_complex> iter_rhs_all	 = fn_tvec_init<fftw_complex> ( dataAmount );
	tvec<fftw_complex> grad_err_all	 = fn_tvec_init<fftw_complex> ( dataAmount );
	rhoJCplxNew.row   = nd;		rhoJCplxNew.col   = alloc_local_sys;
	rhoJCplxTrans.col = nd;		rhoJCplxTrans.row = alloc_local_sys;
	gradientOld.row	  = xlen;	gradientOld.col	  = alloc_local_sys;
	rhoJCplxAll.row = nd;		rhoJCplxAll.col = cplxReDofs;

	/* temporary variables for APG method; */
	tvec<fftw_complex> rhoJCplxTmp0  = fn_tvec_init<fftw_complex> ( rhsLen );
	tvec<fftw_complex> rhoJCplxTmp1  = fn_tvec_init<fftw_complex> ( rhsLen );
	tvec<fftw_complex> rhoJCplxInp	 = fn_tvec_init<fftw_complex> ( rhsLen );
	tvec<fftw_complex> rhoJCplxDiff  = fn_tvec_init<fftw_complex> ( alloc_local_sys );
	tvec<fftw_complex> gradientDiff  = fn_tvec_init<fftw_complex> ( alloc_local_sys );
	memcpy( rhoJCplxTmp0.val, ssysparam->scbndv->rhoJCplx.val, sizeof(fftw_complex) * rhsLen );
	memcpy( rhoJCplxTmp1.val, ssysparam->scbndv->rhoJCplx.val, sizeof(fftw_complex) * rhsLen );
	memcpy( rhoJCplxInp.val,  ssysparam->scbndv->rhoJCplx.val, sizeof(fftw_complex) * rhsLen );
	rhoJCplxTmp0.row   = nd;		rhoJCplxTmp0.col   = alloc_local_sys;
	rhoJCplxTmp1.row   = nd;		rhoJCplxTmp1.col   = alloc_local_sys;
	rhoJCplxInp.row    = nd;		rhoJCplxInp.col    = alloc_local_sys;

	/* initialization; */
	int iterator = iter_load;
	double step_size = step0;
	int	status = 1;
	double res = 1.0;
	double hamilton, hamiltonOld, diffham;
	char dataName[FILELEN];
	sprintf(dataName, "%s/sys_energy_error.dat", rsltDir);

	/* energy decay; */
	int reCount = 0;
	int isDecay = 1;
	hamiltonOld = 100.0;
	diffham		= 1.0;
	tvec<double> energyNew	= fn_tvec_init<double> ( 3 );
	fn_tvec_setZero<double> ( energyNew );

	/* innSMatd0JJ; */
	tCCSmat<double> G0innd0JJ = fn_tensor_diag_dCCSmat ( 
			ssysparam->sfftv->Gsquare, ssysparam->sGJPv->innSMatd0JJ, 0 );
	ssysparam->iter_matrix	  = fn_tCCSmat_init<double> (
			G0innd0JJ.row, G0innd0JJ.col, ssysparam->interact_grad.nnz + G0innd0JJ.nnz );

	/* parameters for Nesterov accelerated method; */
	double	theta = 1.0;
	double	q	  = 0.0;

	/* parameters for inexact linear search; */
	bool	isBreak = true;
	double	rho		= 0.5 * (sqrt(5) - 1.0);
	double	delta	= 1.0e-14;

	/* data for reloading; */
	FILE *fdata, *tmp_fdata; 
	if ( iter_load > 0 )
	{
		/* load step size; */
		if ( access(dataName, F_OK) == 0 ) // exist;
		{
			/* memory allocation for existing data; */
			tvec<double> tmp_step	 = fn_tvec_init<double> ( iter_load );
			tvec<double> tmp_theta	 = fn_tvec_init<double> ( iter_load );
			tvec<double> tmp_energy0 = fn_tvec_init<double> ( iter_load );
			tvec<double> tmp_energy1 = fn_tvec_init<double> ( iter_load );
			tvec<double> tmp_energy2 = fn_tvec_init<double> ( iter_load );
			tvec<double> tmp_res	 = fn_tvec_init<double> ( iter_load );
			tvec<int>	 tmp_isDacy	 = fn_tvec_init<int>	( iter_load );

			/* load data before 'iter_load'; */
			int val;
			tmp_fdata = fopen(dataName, "r");
			for ( int i = 0; i < iter_load; i++ )
			{
				val = fscanf(tmp_fdata, "%lf", &(tmp_step.val[i]));
				val = fscanf(tmp_fdata, "%lf", &(tmp_theta.val[i]));
				val = fscanf(tmp_fdata, "%lf", &(tmp_energy0.val[i]));
				val = fscanf(tmp_fdata, "%lf", &(tmp_energy1.val[i]));
				val = fscanf(tmp_fdata, "%lf", &(tmp_energy2.val[i]));
				val = fscanf(tmp_fdata, "%lf", &(tmp_res.val[i]));
				val = fscanf(tmp_fdata, "%d",  &(tmp_isDacy.val[i]));
			}
			fclose(tmp_fdata);
			step_size = tmp_step.val[iter_load-1]; // obtain the step size at 'iter_load';
			theta = tmp_theta.val[iter_load-1]; // obtain the parameter 'theta' at 'iter_load';
			if ( iter_load == 1 )
				diffham = -10.0;
			else
				diffham = tmp_energy2.val[iter_load-1] - tmp_energy2.val[iter_load-2];

			/* count of restarting operators; */
			for ( int j1 = 0; j1 < iter_load; j1++ )
				if ( isDecay == 0 )
					reCount ++ ;

			/* rewrite the existing data; */
			fdata  = fopen(dataName, "w");
			if ( myrank == 0 )
				for ( int i = 0; i < iter_load; i++ )
				{
					fprintf(fdata, "%.6E\t%.15E\t%+.15E\t%+.15E\t%+.20E\t%+.15E\t%d\n", 
							tmp_step.val[i], tmp_theta.val[i], tmp_energy0.val[i], 
							tmp_energy1.val[i], tmp_energy2.val[i], tmp_res.val[i], 
							tmp_isDacy.val[i]);
					printf("%.6E\t%.15E\t%+.15E\t%+.15E\t%+.20E\t%+.15E\t%d\n", 
							tmp_step.val[i], tmp_theta.val[i], tmp_energy0.val[i], 
							tmp_energy1.val[i], tmp_energy2.val[i], tmp_res.val[i], 
							tmp_isDacy.val[i]);
				}

			/* release memory of temporay variables; */
			fn_tvec_free<double> ( tmp_step );
			fn_tvec_free<double> ( tmp_theta );
			fn_tvec_free<double> ( tmp_energy0 );
			fn_tvec_free<double> ( tmp_energy1 );
			fn_tvec_free<double> ( tmp_energy2 );
			fn_tvec_free<double> ( tmp_res );
			fn_tvec_free<int>	 ( tmp_isDacy );
		}
		else
		{
			fdata = fopen(dataName, "w");
		}

		/* load previous 'rhoJCplx'; */
		char frFile[FILELEN];
		sprintf(frFile, "%s/rank%d/sys_rhoJCplx%d.dat", rsltDir, myrank, iter_load-1);
		FILE *fp = fopen(frFile, "r");
		if ( fp == NULL )	// invalid file;
		{
			if ( myrank == 0 )
				printf("Unable to read file '%s'. No such file or directory.\n", frFile);
			return 1;
		}
		else
		{
			int val;
			val = fscanf(fp, "%d", &(rhoJCplxTmp0.row));
			val = fscanf(fp, "%d", &(rhoJCplxTmp0.col));
			val = fscanf(fp, "%d", &(rhoJCplxTmp0.len));
			for ( int i = 0; i < rhsLen; i++ )
			{
				val = fscanf(fp, "%lf", &(rhoJCplxTmp0.val[i][0]));
				val = fscanf(fp, "%lf", &(rhoJCplxTmp0.val[i][1]));
			}
		}
		fclose(fp);	

		/* calculate 'gradientOld'; */
		fn_calc_system_energy ( ssysparam, rhoJCplxTmp0, energy, true, isTest );
		memcpy ( gradientOld.val, ssysparam->sfftv->gradient.val, 
					sizeof(fftw_complex) * gradientOld.len );
		hamiltonOld = energy.val[2];

		/* restart if energy is not decay; */
		if ( diffham < 0 )
		{
			/* set flag to check if the hamilton energy is decay; */
			isDecay = 1;
            /* the parameter for Lagrange extrapolation; */
			double theta2 = pow(theta, 2);
            double thetaTmp = -0.5*(theta2-q) + sqrt(0.25*pow(theta2-q,2) + theta2);
            double beta = theta*(1.0-theta) / (theta2+thetaTmp);
            theta = thetaTmp;
		}
		else
		{
			/* set flag to check if the hamilton energy is decay; */
			isDecay = 0;
			/* restart; */
            theta = 1.0;
			if ( myrank == 0 )
				printf("--> restart\n");
			reCount ++ ;
		}
	}
	else
	{
		fdata = fopen(dataName, "w");
	}

	if ( myrank == 0 )
		printf("---> start iteration.\n");
	while ( iterator < iter_max )
	{
		/* calculate free energy and entropy gradient; */
		fn_calc_system_energy ( ssysparam, ssysparam->scbndv->rhoJCplx, energy, true, isTest );
		hamilton = energy.val[2];

		/* estimate the time stepping size by BB method; */
		if ( iterator > 0 )
		{
			/* difference of order parameters; */
			fn_diff_weight_rhoJCplx ( ssysparam, rhoJCplxTmp1, 
										rhoJCplxTmp0, rhoJCplxDiff );
			fn_diff_weight_gradient ( ssysparam, ssysparam->sfftv->gradient, 
										gradientOld, gradientDiff );
			if ( isTest && iterator == iter_load )
			//if ( isTest && iterator == 1 )
			{
				char fwFile[FILELEN];
				sprintf(fwFile, "%s/rank%d/rhoJCplxTmp0.dat", rsltDir, myrank);
				fn_tvec_save_complex ( rhoJCplxTmp0, fwFile );
				sprintf(fwFile, "%s/rank%d/rhoJCplxTmp1.dat", rsltDir, myrank);
				fn_tvec_save_complex ( rhoJCplxTmp1, fwFile );
				sprintf(fwFile, "%s/rank%d/rhoJCplxDiff.dat", rsltDir, myrank);
				fn_tvec_save_complex ( rhoJCplxDiff, fwFile );

				sprintf(fwFile, "%s/rank%d/gradientOld.dat", rsltDir, myrank);
				fn_tvec_save_complex ( gradientOld, fwFile );
				sprintf(fwFile, "%s/rank%d/gradient.dat", rsltDir, myrank);
				fn_tvec_save_complex ( ssysparam->sfftv->gradient, fwFile );
				sprintf(fwFile, "%s/rank%d/gradientDiff.dat", rsltDir, myrank);
				fn_tvec_save_complex ( gradientDiff, fwFile );
			}
	
			/* calculate BB stepping size; */
			fftw_complex bbTmp0, bbTmp1, bbStep;
			if ( bbType == 1 )
			{
				/* BB stepping size (type 1); */
				/* first part; */
				fn_tvec_dotMultiplySum_complex ( rhoJCplxDiff, rhoJCplxDiff, dataLen, bbTmp0 );
				/* second part; */
				fn_tvec_dotMultiplySum_complex ( rhoJCplxDiff, gradientDiff, dataLen, bbTmp1 );
			}
			else
			{
				/* BB stepping size (type 2); */
				/* first part; */
				fn_tvec_dotMultiplySum_complex ( rhoJCplxDiff, gradientDiff, dataLen, bbTmp0 );
				/* second part; */
				fn_tvec_dotMultiplySum_complex ( gradientDiff, gradientDiff, dataLen, bbTmp1 );
			}
			if ( myrank == 0 )
			{
				fn_complex_divide ( bbTmp0, bbTmp1, bbStep ); // bbStep = bbTmp0 / bbTmp1;
				step_size  = fabs ( bbStep[0] );
				/* 'step_size' must belong to the region [step_min, step_max]; */
				step_size  = ( step_size > step_min ? step_size : step_min );
				step_size  = ( step_size < step_max ? step_size : step_max );
			}
			MPI_Barrier ( MPI_COMM_WORLD );
			MPI_Bcast ( &step_size, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
			isBreak = false;
		}

		/* generate the right term for iteration; */
		fn_obt_iter_rhs ( ssysparam, ssysparam->scbndv->rhoJCplx,
						ssysparam->sfftv->gradient, iterator, isTest );

		/* search time stepping size for energy decay; */
		while ( true )
		{
			/* adopt the minimal stepping size and then break recurrence; */
			if ( step_size <= step_min )
			{
				step_size = step_min;
				isBreak = true;
			}

			/* obtain the iteration matrix; */
			//fn_obt_iter_matrix ( ssysparam, step_size, iterator, isTest );
			fn_tCCSmat_free<double> ( ssysparam->iter_matrix );
			ssysparam->iter_matrix	  = fn_add_dCCSmat ( 
					ssysparam->interact_grad, G0innd0JJ, model_xi, 1.0/step_size );

			/* generate the right term for iteration; */
			memcpy ( iter_rhs.val, ssysparam->iter_rho_rhs.val, sizeof(fftw_complex) * rhsLen );
			fn_tvec_add_complex ( iter_rhs, ssysparam->iter_entropy_rhs, 1.0/step_size, -1.0 );
			fn_obt_iter_rhs_all ( ssysparam, iter_rhs, iter_rhs_all );

			/* calculate the new 'rhoJCplx'; */
			status = fn_umfpack_complex_solver ( 
						ssysparam->iter_matrix, iter_rhs_all, rhoJCplxAll );

			/* mass conservation; */
			if ( ssysparam->massFlag == 1 )
				for ( int j1 = 0; j1 < nd; j1++ )
				{
					int ind0 = j1; // the first element along Fourier direction;
					rhoJCplxAll.val[ind0][0] = ssysparam->massVec.val[j1][0];
					rhoJCplxAll.val[ind0][1] = ssysparam->massVec.val[j1][1];
				}

			/* obtain 'rhoJCplx' on each process; */
			for ( int j0 = 0; j0 < alloc_local_sys; j0++ )
			{
				for ( int j1 = 0; j1 < nd; j1++ )
				{
					//int ind0 = j0 * nd + j1; // no transpose;
					int ind0 = j1 * alloc_local_sys + j0;
					int ind1 = displs[myrank] + j0*nd + j1;
					rhoJCplxTrans.val[ind0][0] = rhoJCplxAll.val[ind1][0];
					rhoJCplxTrans.val[ind0][1] = rhoJCplxAll.val[ind1][1];
				}
			}

			/* compute energy about the temporary order parameter; */
			fn_calc_system_energy ( ssysparam, rhoJCplxTrans, energyNew, false, isTest );
			
			/* ensure enough energy decay; */
			fn_diff_weight_rhoJCplx ( ssysparam, rhoJCplxTrans, 
								ssysparam->scbndv->rhoJCplx, rhoJCplxDiff );
			double energy_decay = fn_tvec_norm_complex ( rhoJCplxDiff ) / sqrt(cplxReDofs);
			energy_decay = energy.val[2] - energyNew.val[2] - delta * energy_decay;
//			printf("rank%d: energy = %.4e, energyNew = %.4e, energy_decay = %.4e\n",
//					myrank, energy.val[2], energyNew.val[2], energy_decay);
//			printf("step size = %.4e, isBreak = %d\n", step_size, isBreak);

			if ( isTest )
			{
				char fwFile[FILELEN];
				sprintf(fwFile, "%s/rank%d/rhoJCplxTrans%d.dat", rsltDir, myrank, iterator);
				fn_tvec_save_complex ( rhoJCplxTrans, fwFile );
				if ( myrank == 0 )
				{
					printf("energy: %.6e, %.6e\n", energy.val[0], energy.val[1]);
					printf("energyNew: %.6e, %.6e\n", energyNew.val[0], energyNew.val[1]);
					printf("step size = %.6e.\n", step_size);
					printf("hamilton = %.6e \t hamilton temp = %.6e\n", 
							energy.val[2], energyNew.val[2]);
					printf("energy decay = %.6e.\n\n", energy_decay);
				}
			}
			if ( energy_decay <= 0 && !isBreak )
				step_size *= rho;
			else
				break;
		}

		/* save data; */
		if ( iterator > iter_load )
		{
			/* save data of current step; */
			if ( (iterator % save_step) == 0 )
				fn_save_system_phase (ssysparam, iterator);

			/* only save 'rhoJCplx' of previous step; */
			if ( (save_step > 1) && ((iterator+1) % save_step) == 0 )
			{
				char fwFile[FILELEN];
				sprintf(fwFile, "%s/rank%d/sys_rhoJCplx%d.dat", rsltDir, myrank, iterator);
				fn_tvec_save_complex ( ssysparam->scbndv->rhoJCplx, fwFile );
			}
		}

		/* copy gradient values; must be before gradient error calculation; */
		memcpy ( gradientOld.val, ssysparam->sfftv->gradient.val, 
					sizeof(fftw_complex) * gradientOld.len );

		/* calculate gradient error; */
		fn_cvec_multiply_dCCSmat ( interact_grad_trans, rhoJCplxAll, grad_err_all );
		for ( int j0 = 0; j0 < dataLen; j0++ )
		{
			for ( int j1 = 0; j1 < nd; j1++ )
			{
				int ind0 = j0 * nd + j1;
				int ind1 = displs[myrank] + j0*nd + j1;
				ssysparam->grad_err.val[ind0][0] = grad_err_all.val[ind1][0];
				ssysparam->grad_err.val[ind0][1] = grad_err_all.val[ind1][1];
			}
		}
		fn_tvec_add_complex ( ssysparam->grad_err, ssysparam->iter_entropy_rhs, model_xi, 1.0 );
		if ( myrank == 0 )
			for ( int i = 0; i < nd; i++ )
				fn_complex_setZero ( ssysparam->grad_err.val[i] );
		res = fn_tvec_maxAbs_complex ( ssysparam->grad_err );

		if ( isTest )
		{
			char fwFile[FILELEN];
			sprintf(fwFile, "%s/rank%d/grad_err%d.dat", rsltDir, myrank, iterator);
			fn_tvec_save_complex ( ssysparam->grad_err, fwFile );
		}
		
		/* compute the difference between the hamilton energies of two adjoin steps; */
		diffham = hamilton - hamiltonOld;
		hamiltonOld = hamilton;

		/* save data; */
		if ( myrank == 0 )
			fprintf(fdata, "%.6E\t%.15E\t%+.15E\t%+.15E\t%+.20E\t%+.15E\t%d\n", 
						step_size, theta, energy.val[0], energy.val[1], 
						energy.val[2], res, isDecay);
		if ( iterator%print_step == 0 )
		{
			//fn_obt_mass ( ssysparam, ssysparam->scbndv->rhoJCplx, mass );
			fn_obt_mass ( ssysparam, rhoJCplxTrans, mass );
			if ( myrank == 0 )
			{
				timer.pause();
				printf("\nIterator %d: step_size = %.6e \t res = %.10e", iterator, step_size, res);
				printf(" \t diffham = %.10e \t isDecay = %d\n", diffham, isDecay);
				printf("\t Laplace term: % .10e \t Nonlinear term: % .10e ", 
						energy.val[0], energy.val[1]);
				printf("\t hamilton = % .10e\n", energy.val[2]);
				printf("\t Mass: %.10e \t theta: %.10e \t time cost: %f seconds\n", 
						fn_complex_abs(mass), theta, timer.get_current_time());
			}
		}

		/* restart if energy is not decay; */
		if ( diffham < 0 )
		{
			/* set flag to check if the hamilton energy is decay; */
			isDecay = 1;
            /* the parameter for Lagrange extrapolation; */
			double theta2 = pow(theta, 2);
            double thetaTmp = -0.5*(theta2-q) + sqrt(0.25*pow(theta2-q,2) + theta2);
            double beta = theta*(1.0-theta) / (theta2+thetaTmp);
            theta = thetaTmp;
            /* Interpolation; */
			memcpy ( rhoJCplxInp.val, rhoJCplxTrans.val, sizeof(fftw_complex) * rhoJCplxTrans.len );
			fn_tvec_add_complex ( rhoJCplxInp, ssysparam->scbndv->rhoJCplx, 1.0+beta, -beta ); 
			/* test data; */
			if ( isTest )
			{
				char fwFile[FILELEN];
				sprintf(fwFile, "%s/rank%d/rhoJCplxInp%d.dat", rsltDir, myrank, iterator);
				fn_tvec_save_complex ( rhoJCplxInp, fwFile );
			}
		}
		else
		{
			/* set flag to check if the hamilton energy is decay; */
			isDecay = 0;
			/* restart; */
            theta = 1.0;
			memcpy ( rhoJCplxInp.val, rhoJCplxTrans.val, sizeof(fftw_complex) * rhoJCplxTrans.len );
			if ( myrank == 0 )
				printf("--> restart\n");
			reCount ++ ;
		}

		/* update 'rhoJCplx'; * transpose for consistency; */
		memcpy ( ssysparam->scbndv->rhoJCplx.val, rhoJCplxTrans.val,
					sizeof(fftw_complex) * rhoJCplxTrans.len );

		/* update 'rhoJCplxTmp0', 'rhoJCplxTmp1'; */
		if ( iterator > iter_load )
			memcpy ( rhoJCplxTmp0.val, rhoJCplxTmp1.val, sizeof(fftw_complex) * rhoJCplxTmp0.len );
		memcpy ( rhoJCplxTmp1.val, rhoJCplxInp.val, sizeof(fftw_complex) * rhoJCplxTmp1.len );

		iterator ++;
//		if ( (res < tol && iterator > iter_load+2) || fabs(diffham) < 1.0e-15 )
//		if ( res < tol && iterator > iter_load+2 )
		if ( res < tol || fabs(diffham) < tolham )
			break;
		if ( reCount > 5 )
			break;
	}

	/* output data about the stationary state; */
	fn_calc_system_energy ( ssysparam, ssysparam->scbndv->rhoJCplx, energy, false, isTest );
	hamilton = energy.val[2];
	if ( myrank == 0 )
		fprintf(fdata, "%.6E\t%.15E\t%+.15E\t%+.15E\t%+.20E\t%+.15E\t%d\n", 
					step_size, theta, energy.val[0], energy.val[1], 
					energy.val[2], res, isDecay);
	if ( myrank == 0 )
	{
		timer.pause();
		printf("\n\t ***** the interface system\n");
		printf("\t ***** iterator = %d \t time cost: %f seconds\n", 
				iterator, timer.get_current_time());
		printf("\t ***** res = % .10e \t diffham = % .10e\n", res, diffham);
		printf("\t ***** Laplace term: % .15e,\t Nonlinear term: % .15e\n", 
				energy.val[0], energy.val[1]);
		printf("\t ***** hamilton = % .20e\n\n", energy.val[2]);
	}
	fclose(fdata);

	/* Releases memory; */
	fn_tvec_free<double>		( energyNew );
	fn_tCCSmat_free<double>		( interact_grad_trans );
	fn_tvec_free<fftw_complex>	( rhoJCplxNew	);
	fn_tvec_free<fftw_complex>	( rhoJCplxTrans );
	fn_tvec_free<fftw_complex>	( iter_rhs		);
	fn_tvec_free<fftw_complex>	( gradientOld	);
	fn_tvec_free<fftw_complex>	( rhoJCplxAll	);
	fn_tvec_free<fftw_complex>	( iter_rhs_all	);
	fn_tvec_free<fftw_complex>	( grad_err_all	);
	fn_tCCSmat_free<double>		( G0innd0JJ		);
	/* temporary variables for APG method; */
	fn_tvec_free<fftw_complex>	( rhoJCplxTmp0 );
	fn_tvec_free<fftw_complex>	( rhoJCplxTmp1 );
	fn_tvec_free<fftw_complex>	( rhoJCplxInp  );
	fn_tvec_free<fftw_complex>	( rhoJCplxDiff );
	fn_tvec_free<fftw_complex>	( gradientDiff );
	free ( displs );		free ( recvCount );

	return iterator;
}


/**
 * \brief	Newton PCG method;
 */
int fn_newton_method	(	stu_system_param		*ssysparam,
							tvec<double>			energy,
								int					iterator,
								bool				isTest )
{
	double	tol			= ssysparam->tol;
	double	tolham		= ssysparam->tolham;
	double	step0		= ssysparam->newton_step_size;
	int		iter_max	= ssysparam->iter_max;
	int		print_step	= ssysparam->print_step;
	int		save_step	= ssysparam->save_step;
	double	mu_para		= ssysparam->pcg_mu_para;

	double	model_xi	= ssysparam->model_xi;
	double	model_tau   = ssysparam->model_tau;
	double	model_gamma = ssysparam->model_gamma;
	double	model_kappa = ssysparam->model_kappa;
	fftw_complex mass;

	/* transpose 'interact_grad'; */
	mytimer_t timer;
	timer.reset();
	timer.start();
	tCCSmat<double> interact_grad_trans = fn_fast_trans_dCCSmat ( ssysparam->interact_grad );

	MPI_Barrier ( MPI_COMM_WORLD );
	timer.pause();
	if ( myrank == 0 )
		printf("\n\t ***** time cost of transpose dCCSmat: %f seconds *****\n\n", 
				timer.get_current_time());

	/* parameters; */
	int		nd	   = ssysparam->sGJPv->nd;
	int cplxReDofs = ssysparam->cplxReDofs;	
	int	rhsLen	   = nd * alloc_local_sys;
	int dataLen;
	if ( myrank == nprocs-1 )
		dataLen = cplxReDofs - displs_sys[myrank];
	else
		dataLen = alloc_local_sys;

	/* memory allocation; */
	tvec<fftw_complex> rhoJCplxNew		  = fn_tvec_init<fftw_complex> ( rhsLen );
	tvec<fftw_complex> PCG_direction	  = fn_tvec_init<fftw_complex> ( rhsLen ); // x^{i} of PCG;
	tvec<fftw_complex> PCG_directionTrans = fn_tvec_init<fftw_complex> ( rhsLen );
	rhoJCplxNew.col			= nd;		rhoJCplxNew.row			= alloc_local_sys;
	PCG_direction.col		= nd;		PCG_direction.row		= alloc_local_sys;
	PCG_directionTrans.row	= nd;		PCG_directionTrans.col	= alloc_local_sys;

	tvec<fftw_complex> rhoJCplxTrans	  = fn_tvec_init<fftw_complex> ( rhsLen );
	rhoJCplxTrans.row		= nd;		rhoJCplxTrans.col		= alloc_local_sys;

	/* initialization; */
	double step_size = step0;
	int subIter = 0;
	int	status = 1;
	double mu, delta;
	double hamilton, diffham;
	double hamiltonOld = energy.val[2];
	double res = fn_tvec_maxAbs_complex ( ssysparam->grad_err );
	char dataName[FILELEN];
	sprintf(dataName, "%s/sys_energy_error.dat", rsltDir);
	FILE *fdata = fopen(dataName, "a");

	int isDecay = 3;
	diffham		= 1.0;

	if ( isTest )
	{
		if ( myrank == 0 )
			printf("mu_para = %.10e \t res = %.10e.\n", mu_para, res);
		char fwFile[FILELEN];
		sprintf(fwFile, "%s/rank%d/newton_rhoJCplx%d.dat", rsltDir, myrank, 0);
		fn_tvec_save_complex ( ssysparam->scbndv->rhoJCplx, fwFile );
		sprintf(fwFile, "%s/rank%d/newton_grad_err%d.dat", rsltDir, myrank, 0);
		fn_tvec_save_complex ( ssysparam->grad_err, fwFile );
	}

	/* obtain Hessian matrix; */
	fn_calc_system_hessian ( ssysparam, ssysparam->scbndv->rhoJCplx, iterator, isTest );

	if ( myrank == 0 )
		printf("---> start iteration.\n");
	while ( iterator < iter_max )
	{
		/* precondition; */
		double tmp = mu_para * res;	
		mu = tmp < 100 ? tmp : 100;
		mu = mu > 1.0e-20 ? mu : 1.0e-20;
		delta = fn_tvec_maxAbs_complex ( ssysparam->hessian );

		/* Pre-conditional Projected Conjugate Gradient method; */
		subIter = fn_calc_system_PCG ( ssysparam, PCG_direction, 
					interact_grad_trans, mu, delta, iterator, false );

		/* transpose; */
		fn_tvec_trans_complex ( PCG_directionTrans, PCG_direction, nd, alloc_local_sys );

		/* rho = - \<x^{i}, b\> / \|x^{i}\|^{2};		b = grad_err; */
		fftw_complex rhoTmp0;
		fn_tvec_dotMultiplySum_complex ( ssysparam->grad_err, PCG_direction, dataLen, rhoTmp0 );
		double rhoTmp1 = fn_tvec_norm_complex ( PCG_direction ) / sqrt(cplxReDofs);
		double rho = -rhoTmp0[0] / pow(rhoTmp1,2);

		/* update 'mu_para'; */
		if ( fabs(rho) > 0.1 )
			mu_para = 0.1 * mu_para;
		else if ( fabs(rho) < 0.01 )
			mu_para = 10.0 * mu_para;
		else
			mu_para = 0.5 * mu_para;

		if ( isTest )
		{
			if ( myrank == 0 )
			{
				printf("mu = %.10e \t delta = %.10e \t rho = %.10e.\n", mu, delta, rho);
				printf("rhoTmp0 = %.5e, %.5e \t rhoTmp1 = %.5e.\n",
						rhoTmp0[0], rhoTmp0[1], rhoTmp1);
				printf("hamilton old = %.10e.\n\n", hamiltonOld);
			}
			char fwFile[FILELEN];
			sprintf(fwFile, "%s/rank%d/PCG_direction%d.dat", rsltDir, myrank, iterator);
			fn_tvec_save_complex ( PCG_direction, fwFile );
		}

		/* update 'rhoJCplx' and free energies; */
		step_size = step0;
		while ( step_size > 1.0e-6 )
		{
			/* update 'rhoJCplx' : rhoJCplx + step_size * PCG_directionTrans; */
			memcpy ( rhoJCplxNew.val, ssysparam->scbndv->rhoJCplx.val, 
						sizeof(fftw_complex) * rhsLen );
			fn_tvec_add_complex ( rhoJCplxNew, PCG_directionTrans, 1.0, step_size );

			/* transpose; */
			fn_tvec_trans_complex ( rhoJCplxTrans, rhoJCplxNew, nd, alloc_local_sys );

			/* mass conservation; */
			//fn_obt_mass ( ssysparam, rhoJCplxTrans, mass );
			//if ( myrank == 0 )
			//	printf("Mass: %.10e\n", fn_complex_abs(mass));
			if ( myrank == 0 && ssysparam->massFlag == 1 )
				for ( int j1 = 0; j1 < nd; j1++ )
				{
					int ind0 = j1 * alloc_local_sys; // the first element along Fourier direction;
					rhoJCplxTrans.val[ind0][0] = ssysparam->massVec.val[j1][0];
					rhoJCplxTrans.val[ind0][1] = ssysparam->massVec.val[j1][1];
				}

			/* transpose; */
			fn_tvec_trans_complex ( rhoJCplxNew, rhoJCplxTrans, alloc_local_sys, nd );

			/* transpose; */
			fn_tvec_trans_complex ( rhoJCplxNew, rhoJCplxTrans, alloc_local_sys, nd );

			/* calculate free energy and entropy gradient; */
			fn_calc_system_energy ( ssysparam, rhoJCplxNew, energy, true, isTest );
			hamilton = energy.val[2];

			/* if decay; */
			diffham = hamilton - hamiltonOld;
			double diffhamCheck = 1.0e-6 * step_size * rhoTmp0[0];
			if ( diffham < diffhamCheck )
				break;
			else
				step_size = 0.5 * step_size;
//			if ( myrank == 0 )
//			{
//				printf("energy: %.10e \t %.10e \t %.10e.\n", 
//						energy.val[0], energy.val[1], energy.val[2]);
//				printf("diffhamCheck = %.5e \t step_size = %.5e \t diffham = %.5e.\n\n", 
//						diffhamCheck, step_size, diffham);
//			}
		}
		hamilton = energy.val[2];
		diffham = hamilton - hamiltonOld;
		hamiltonOld = hamilton;
		if ( diffham < 0 )
			isDecay = 3;
		else
			isDecay = 2;

		/* calculate gradient error; */
		/*
		memcpy ( ssysparam->grad_err.val, rhoJCplxNew.val, 
					sizeof(fftw_complex) * rhoJCplxNew.len );
		fn_tvec_add_complex ( ssysparam->grad_err, ssysparam->scbndv->rhoJCplx, 1.0, -1.0 );
		res = fn_tvec_maxAbs_complex ( ssysparam->grad_err );
		res /= step_size;
		*/

		/* calculate gradient error; */
		//fn_calc_system_energy ( ssysparam, rhoJCplxNew, energy, true, isTest );
		fn_obt_iter_rhs	( ssysparam, rhoJCplxNew, ssysparam->sfftv->gradient, iterator, isTest );
		fn_tvec_trans_complex ( rhoJCplxTrans, rhoJCplxNew, alloc_local_sys, nd );
		fn_cvec_multiply_dCCSmat ( interact_grad_trans, rhoJCplxTrans, ssysparam->grad_err );
		fn_tvec_add_complex ( ssysparam->grad_err, ssysparam->iter_entropy_rhs, model_xi, 1.0 );

		for ( int i = 0; i < nd; i++ )
			fn_complex_setZero ( ssysparam->grad_err.val[i] );
		res = fn_tvec_maxAbs_complex ( ssysparam->grad_err );

		/* update 'rhoJCplx';	transpose for consistency; */
		memcpy ( ssysparam->scbndv->rhoJCplx.val, rhoJCplxNew.val,
					sizeof(fftw_complex) * rhoJCplxNew.len );

		/* update Hessian matrix;	interact_grad_trans; */
		fn_calc_system_hessian ( ssysparam, rhoJCplxNew, iterator, isTest );

		if ( isTest )
		{
			char fwFile[FILELEN];
			sprintf(fwFile, "%s/rank%d/Newton_rhoJCplxNew%d.dat", rsltDir, myrank, iterator);
			fn_tvec_save_complex ( rhoJCplxNew, fwFile );
			sprintf(fwFile, "%s/rank%d/Newton_iter_entropy_rhs%d.dat", rsltDir, myrank, iterator);
			fn_tvec_save_complex ( ssysparam->iter_entropy_rhs, fwFile );
			sprintf(fwFile, "%s/rank%d/Newton_grad_err%d.dat", rsltDir, myrank, iterator);
			fn_tvec_save_complex ( ssysparam->grad_err, fwFile );
		}

		/* save data; */
		if ( myrank == 0 )
			fprintf(fdata, "%.6E\t%+.15E\t%+.15E\t%+.20E\t%+.15E\t%d\n", step_size,
						energy.val[0], energy.val[1], energy.val[2], res, isDecay);
		if ( iterator%print_step == 0 )
		{
			fn_obt_mass ( ssysparam, ssysparam->scbndv->rhoJCplx, mass );
			if ( myrank == 0 )
			{
				timer.pause();
				printf("\nIterator %d: step_size = %.6e \t res = %.10e", 
						iterator, step_size, res);
				printf(" \t diffham = %.10e \t isDecay = %d\n", diffham, isDecay);
				printf("\t Laplace term: % .10e \t Nonlinear term: % .10e ", 
						energy.val[0], energy.val[1]);
				printf("\t hamilton = % .10e\n", energy.val[2]);
				printf("\t Mass: %.10e \t time cost: %f seconds\n", 
						fn_complex_abs(mass), timer.get_current_time());
			}
		}

		if ( (iterator > 0 ) && (iterator % save_step) == 0 )
			fn_save_system_phase (ssysparam, iterator);		// save data;

		iterator ++;
//		if ( (res < tol && iterator > 2) || fabs(diffham) < 1.0e-15 )
//		if ( res < tol && iterator > 2 )
//		if ( res < tol )
		if ( res < tol || fabs(diffham) < tolham )
			break;
	}

	/* output data about the stationary state; */
	fn_calc_system_energy ( ssysparam, ssysparam->scbndv->rhoJCplx, energy, false, isTest );
	hamilton = energy.val[2];
	if ( myrank == 0 )
		fprintf(fdata, "%.6E\t%+.15E\t%+.15E\t%+.20E\t%+.15E\t%d\n", step_size,
				energy.val[0], energy.val[1], energy.val[2], res, isDecay);
	if ( myrank == 0 )
	{
		timer.pause();
		printf("\n\t ***** the interface system\n");
		printf("\t ***** iterator = %d \t step size = %.6e \t time cost: %f seconds\n", 
				iterator, step_size, timer.get_current_time());
		printf("\t ***** res = % .10e \t diffham = % .10e\n", res, diffham);
		printf("\t ***** Laplace term: % .15e,\t Nonlinear term: % .15e\n", 
				energy.val[0], energy.val[1]);
		printf("\t ***** hamilton = % .20e\n\n", energy.val[2]);
	}
	fclose(fdata);

	/* Releases memory; */
	fn_tCCSmat_free<double>		( interact_grad_trans );
	fn_tvec_free<fftw_complex>	( rhoJCplxNew	);
	fn_tvec_free<fftw_complex>	( rhoJCplxTrans );
	fn_tvec_free<fftw_complex>	( PCG_direction );
	fn_tvec_free<fftw_complex>	( PCG_directionTrans );

	return iterator;
}
