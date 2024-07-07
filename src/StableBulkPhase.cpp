/*! \file	StableBulkPhase.cpp
 *
 * \brief	obtain the stable state of the bulk phase;
 *
 */

#include "Head.h"
#include "DataOperators.h"
#include "Mytimer.h"
#include "functs.h"


/**
 * \brief	Initialization and preparation;
 *			update bulk phase to obtain the stable state;
 */
double fn_obt_stable_bulk_phase (	stu_bulk_param	*sbulkparam )
{
	/* directory for the truncation data; */
	char bulkDir[FILELEN0];
	sprintf(bulkDir, "./bulk/%s_%s", model_type, sbulkparam->phase);
	mkdir("./bulk/", 0755);
	mkdir(bulkDir, 0755); // create fold;

	/* file for saving; */
	if ( strcmp( main_type, "stable_bulk1" ) != 0 && 
			strcmp( main_type, "stable_bulk2" ) != 0 ) // not stable bulk phases;
	{
		char bulkFile[FILELEN];
		int paraOrder = fn_paraOrder_bulk_phase	( sbulkparam );
		sprintf(bulkFile, "./bulk/%s_%s/para%d.dat", model_type, sbulkparam->phase, paraOrder);
		if ( access(bulkFile, F_OK) == 0 )
		{
			if ( myrank == 0 )
			{
				printf(" <========== load %s of the truncation of stable ", bulkFile);
				if ( sbulkparam->sflag == 1 )		printf("left bulk phase.\n");
				else if ( sbulkparam->sflag == 2 )	printf("right bulk phase.\n");
				printf("\n\n");
			}
			return fn_load_trunc_bulk_phase	( sbulkparam );
		}
		else
		{
			if ( myrank == 0 )
				printf("%s not exist.\n", bulkFile);
		}
	}


	bool isTest = false;

	if ( myrank == 0 )
	{
		printf(" <========== Fixed box, generate equilibrium state of ");
		if ( sbulkparam->sflag == 1 )
			printf("left bulk phase ==========> \n\n");
		else if ( sbulkparam->sflag == 2 )
			printf("right bulk phase ==========> \n\n");
	}
	mytimer_t timer;
	timer.reset();
	timer.start();

	/* initialization and preparation; */
	fn_bulk_memory_allocation	( sbulkparam );			// memory allocation;
	fn_serial_obt_kIndex ( sbulkparam->sfftv->indKspace, sbulkparam->NCpt ); // obtain Fourier k;
	fn_get_projection_plane	( sbulkparam );				// project plane; R'PBk;
	fn_get_Gsquare			( sbulkparam );				// Gsquare; |R'PBk|^2;
	fn_get_init_bulk_phase	( sbulkparam );				// obtain the initial rho in Fourier space;
	fn_save_bulk_phase		( sbulkparam, 0, 0 );		// save initial value;

	/* file for saving free energies; */
	char dataName[FILELEN];
	sprintf(dataName, "%s/bulk%d_energy_error.dat", rsltDir, sbulkparam->sflag);
	FILE *fdata = fopen(dataName, "w");

	/* optimize computational box; */
	int opt_iter_max = sbulkparam->opt_iter_max;
	double opt_tol	 = sbulkparam->opt_tol;
	int opt_iter = 0;
	double opt_err = 1.0;
	double hamilton, hamiltonOld;

	/* update bulk phase to obtain the stable state; */
	hamilton = fn_update_bulk_phase	( sbulkparam, fdata, opt_iter );
	hamiltonOld = hamilton;

	while ( opt_err > opt_tol && opt_iter < opt_iter_max )
	{
		fn_opt_rcpBox ( sbulkparam );	// optimize computational box;
		hamilton = fn_update_bulk_phase	( sbulkparam, fdata, opt_iter ); // calculate stable state;

		/* calculate error; */
		opt_err = fabs( hamilton - hamiltonOld );
		if ( myrank == 0 )
			printf("\t ***** Iter %d: hamilton = % .20e\terr = % .10e.\n\n", 
					opt_iter, hamilton, opt_err);
		
		/* update energy; */
		hamiltonOld = hamilton;
		opt_iter ++ ;
	}
	fclose(fdata);

	/* update 'dirBox'; */
	matrix_inverse ( sbulkparam->rcpBox.val, sbulkparam->dirBox.val, sbulkparam->dimCpt	);
	for ( int i = 0; i < sbulkparam->dimCpt; i++ )
		for ( int j = 0; j < sbulkparam->dimCpt; j++ )
			sbulkparam->dirBox.val[i][j] *= 2*PI;
	/* update 'translProjBoxVec'; */
	fn_obt_translProjBox_matrix ( sbulkparam );
	if ( myrank == 0 )
	{
		printf("rcpBox: \n");
		fn_tmat_print<double> ( sbulkparam->rcpBox );
		printf("\n");
	}

	/* take translation in the bulk phase; */
	fn_save_bulk_phase ( sbulkparam, -1, -1 ); // save stationary value before translation;
	fn_translate_bulk_phase	( sbulkparam );
	fn_save_bulk_phase ( sbulkparam, -1, -2 ); // save stationary value after translation;

	/* only translation without rotation; */
	if ( isTest )
	{
		int dimPhy = sbulkparam->dimPhy;
		int dimCpt = sbulkparam->dimCpt;

		/* save rotation matrix; */
		tmat<double> rotateTmp = fn_tmat_init<double> ( dimPhy, dimPhy );
		for ( int j1 = 0; j1 < dimPhy; j1++ )
		{
			for ( int j2 = 0; j2 < dimPhy; j2++ )
			{
				rotateTmp.val[j1][j2] = sbulkparam->rotateMat.val[j1][j2];
				if ( j1 == j2 )
					sbulkparam->rotateMat.val[j1][j2] = 1.0;
				else
					sbulkparam->rotateMat.val[j1][j2] = 0.0;
			}
		}
		fn_obt_rotateProjBox_matrix ( sbulkparam );
		if ( myrank == 0 )
			fn_matrix_print ( sbulkparam->rotateProjBoxMat );
		fn_get_projection_plane	( sbulkparam );				// project plane; R'PBk;
		fn_save_bulk_phase ( sbulkparam, -1, -3 ); // save stationary value after translation;

		/* recover rotation matrix; */
		for ( int j1 = 0; j1 < dimPhy; j1++ )
			for ( int j2 = 0; j2 < dimPhy; j2++ )
				sbulkparam->rotateMat.val[j1][j2] = rotateTmp.val[j1][j2];
		fn_obt_rotateProjBox_matrix ( sbulkparam );
		fn_get_projection_plane	( sbulkparam );				// project plane; R'PBk;
		fn_tmat_free<double> ( rotateTmp );
	}

	/* save optimized reciprocal box; */
	if ( myrank == 0 )
	{
		char fname[FILELEN];
		sprintf(fname, "%s/parameter_init.dat", rsltDir);
		FILE *fp = fopen(fname, "a+");
		fprintf(fp, "\n\n");
		fprintf(fp, "bulk%d_opt_rcpBox\t\t= \n", sbulkparam->sflag);
		fn_matrix_save ( sbulkparam->rcpBox, fp);
		fclose(fp);
	}

	/* costing time; */
	timer.pause();
	if ( myrank == 0 )
	{
		printf("\n\t ***** time cost of the calculation of stable ");
		if ( sbulkparam->sflag == 1 )
			printf("left  bulk phase: %f seconds *****\n\n", timer.get_current_time());
		else if ( sbulkparam->sflag == 2 )
			printf("right bulk phase: %f seconds *****\n\n", timer.get_current_time());
	}

	/* release memory; only save rhoCplx; */
	fn_bulk_memory_free	(sbulkparam, "stable bulk");

	/* truncation to reduce computational cost; */
	fn_truncation_bulk_phase ( sbulkparam, hamilton );

	return hamilton;
}


/**
 * \brief	Allocation memory for calculating stable bulk phase;
 */
void fn_bulk_memory_allocation ( stu_bulk_param		*sbulkparam )
{
	/* load parameters; */
	int cplxDofs = sbulkparam->cplxDofs;
	int dimPhy	 = sbulkparam->dimPhy;
	int dimCpt	 = sbulkparam->dimCpt;

	/* memory allocation; */
	sbulkparam->sfftv->indKspace = fn_tmat_init<int>			( cplxDofs, dimCpt );
	sbulkparam->sfftv->projPlane = fn_tmat_init<double>			( cplxDofs, dimPhy );
	sbulkparam->sfftv->Gsquare   = fn_tvec_init<double>			( cplxDofs );

	sbulkparam->sfftv->rhoCplx   = fn_tvec_init<fftw_complex>	( cplxDofs );
	sbulkparam->sfftv->rhoReal	 = fn_tvec_init<fftw_complex>	( cplxDofs );

	sbulkparam->sfftv->fftw_Ctmp = fn_tvec_init<fftw_complex>	( cplxDofs );
	sbulkparam->sfftv->fftw_Rtmp = fn_tvec_init<fftw_complex>	( cplxDofs );
	sbulkparam->sfftv->gradient  = fn_tvec_init<fftw_complex>	( cplxDofs );
	sbulkparam->sfftv->cplxTmp   = fn_tvec_init<fftw_complex>	( cplxDofs );
	sbulkparam->sfftv->quadTerm  = fn_tvec_init<fftw_complex>	( cplxDofs );
	sbulkparam->sfftv->cubTerm   = fn_tvec_init<fftw_complex>	( cplxDofs );
	sbulkparam->sfftv->quarTerm  = fn_tvec_init<fftw_complex>	( cplxDofs );

	/* initialization; */
	fn_tmat_setZero<int>	( sbulkparam->sfftv->indKspace );
	fn_tmat_setZero<double> ( sbulkparam->sfftv->projPlane );
	fn_tvec_setZero<double> ( sbulkparam->sfftv->Gsquare  );
	fn_tvec_setZero_complex ( sbulkparam->sfftv->rhoCplx  );
	fn_tvec_setZero_complex ( sbulkparam->sfftv->rhoReal  );
	fn_tvec_setZero_complex ( sbulkparam->sfftv->gradient );
	fn_tvec_setZero_complex ( sbulkparam->sfftv->cplxTmp  );

	/* prepare FFTW plan; */
	int offset = 0;
	char strNCpt[FILELEN-100];
	offset += sprintf(strNCpt, "%d_", dimCpt);
	for ( int i = 0; i < dimCpt; i++ )
	{
		offset += sprintf(strNCpt + offset, "%d_", sbulkparam->NCpt.val[i]);
	}
	offset += sprintf(strNCpt + offset, "%d.txt", flags);

	char fnr2c[FILELEN], fnc2r[FILELEN];
	sprintf(fnr2c, "plan/bulk%d_forward_%s",  sbulkparam->sflag, strNCpt);
	sprintf(fnc2r, "plan/bulk%d_backward_%s", sbulkparam->sflag, strNCpt);

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

	sbulkparam->sfftv->Planc2cFord = fftw_plan_dft ( dimCpt, 
			sbulkparam->NCpt.val, sbulkparam->sfftv->fftw_Rtmp.val, 
			sbulkparam->sfftv->fftw_Ctmp.val, FFTW_FORWARD,  flags); // real to cplx
	sbulkparam->sfftv->Planc2cBack = fftw_plan_dft ( dimCpt, 
			sbulkparam->NCpt.val, sbulkparam->sfftv->fftw_Ctmp.val, 
			sbulkparam->sfftv->fftw_Rtmp.val, FFTW_BACKWARD, flags); // cplx to real 

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
 * \brief	Get projection plane; i.e. R'PBk;
 *			R is rotateMat;
 *			P is projMat;
 *			B is rcpBox;
 *			k is Fourier index;
 *			R'PB: [dimPhy, dimCpt];		k: [dimCpt, 1];
 */
void fn_get_projection_plane	(	stu_bulk_param		*sbulkparam	)
{
	double RPBKtmp;
	for (int i = 0; i < sbulkparam->cplxDofs; i++)
	{
		for (int kk = 0; kk < sbulkparam->dimPhy; kk++)
		{
			RPBKtmp = 0.0;
			for (int jj = 0; jj < sbulkparam->dimCpt; jj++)
				RPBKtmp += sbulkparam->rotateProjBoxMat.val[kk][jj] * 
							sbulkparam->sfftv->indKspace.val[i][jj]; // R'PB * k;
			sbulkparam->sfftv->projPlane.val[i][kk] = RPBKtmp;
		}
	}
}


/**
 * \brief	Get Gsquare; i.e. |R'PBk|^2;
 *			based on 'fn_get_projection_plane';
 */
void fn_get_Gsquare			(	stu_bulk_param		*sbulkparam )
{
	for (int i = 0; i < sbulkparam->cplxDofs; i++)
	{
		sbulkparam->sfftv->Gsquare.val[i] = 0.0;
		for (int kk = 0; kk < sbulkparam->dimPhy; kk++)
		{
			sbulkparam->sfftv->Gsquare.val[i] += 
				pow(sbulkparam->sfftv->projPlane.val[i][kk], 2);
		}
	}
}


/**
 * \brief	Get the initial state of bulk phase;
 */
void fn_get_init_bulk_phase	(	stu_bulk_param	*sbulkparam )
{
	/* parameters; */
	int dimCpt = sbulkparam->dimCpt;
	int kk;
	double translVal;
	double tmpReal, tmpCplx;
	int cplxDofs = sbulkparam->cplxDofs;

	/* memory allocation; */
	tvec<int> kIndex = fn_tvec_init<int> ( dimCpt );

	/* obtain all information of rhoCplx; */
	for (int i = 0; i < sbulkparam->initNum; i++)
	{
		/* k to global index; */
		kk = fn_serial_kIndex_to_gIndex ( sbulkparam->initIndex.val[i], sbulkparam->NCpt );
		/* assign to initial value; */
		sbulkparam->sfftv->rhoCplx.val[kk][0] = sbulkparam->initCoeff.val[i][0];
		sbulkparam->sfftv->rhoCplx.val[kk][1] = sbulkparam->initCoeff.val[i][1];

		/* symmetry; */
		if ( sbulkparam->initIndex.val[i][dimCpt-1] > 0 )
		{
			for ( int j = 0; j < dimCpt; j++ )
				kIndex.val[j] = -1 * sbulkparam->initIndex.val[i][j];
			/* k to global index; */
			kk = fn_serial_kIndex_to_gIndex ( kIndex.val, sbulkparam->NCpt );
			/* assign to initial value; */
			sbulkparam->sfftv->rhoCplx.val[kk][0] = sbulkparam->initCoeff.val[i][0];
			sbulkparam->sfftv->rhoCplx.val[kk][1] = sbulkparam->initCoeff.val[i][1];
		}
	}
	fn_complex_setZero ( sbulkparam->sfftv->rhoCplx.val[0] );	// mass conservation;
	fn_tvec_free<int>  ( kIndex );
}


/**
 * \brief	Update bulk phase and obtain the stable bulk phase;
 */
double fn_update_bulk_phase (	stu_bulk_param		*sbulkparam,
									FILE			*fdata,
									int				opt_iter )
{
	/* load iteration parameters; */
	double tol		= sbulkparam->tol;
	double stepSize = sbulkparam->step_size;
	int itMax		= sbulkparam->iter_max;
	int optSave		= sbulkparam->opt_save_step;
	int stepSave	= sbulkparam->save_step;

	/* initialization; */
	int iterator = 0;
	double res = 1.0;
	double tmp, temp;
	clock_t start, finish;
	double duration, hamilton, oldHamilton, diffham;

	oldHamilton = 100.0;
	diffham		= 0.0;
	tvec<double> energy = fn_tvec_init<double> ( 3 );
	fn_tvec_setZero<double> ( energy );

	while (res > tol && iterator < itMax)
	{
		/* calculate hamilton energy; */
		fn_calc_bulk_energy ( sbulkparam, energy );
		hamilton = energy.val[2];

		/* save data; */
		if ( myrank == 0 )
		{
			fprintf(fdata, "%d\t%d\t%+.15E\t%+.15E\t%+.20E\t%+.15E\n", 
					opt_iter, iterator, energy.val[0], energy.val[1], energy.val[2], res);
			if ( iterator%sbulkparam->print_step == 0 )
			{
				printf("\nIter %d: stepSize = %.5e \t res = %.10e \t diffham = %.10e\n", 
						iterator, stepSize, res, diffham);
				printf("\t Laplace term: % .10e \t Nonlinear term: % .10e,", 
						energy.val[0], energy.val[1]);
				printf("\t hamilton = % .10e\n", energy.val[2]);
			}
		}

		if ( (opt_iter % optSave) == 0 && (iterator > 0) && (iterator % stepSave) == 0 )
			fn_save_bulk_phase (sbulkparam, opt_iter, iterator);	// save data;

		diffham = hamilton - oldHamilton;
		oldHamilton = hamilton;

		for (int i = 0; i < sbulkparam->cplxDofs; i++)
		{
			/* [(q0^2-|G|^2)*(q1^2-|G|^2)*...]^2; */
			temp = 1.0;
			for (int j = 0; j < sbulkparam->scale_num; j++)
			{
				tmp = pow(sbulkparam->scale_val[j], 2) - sbulkparam->sfftv->Gsquare.val[i];
				temp *= pow(tmp, 2);
			}

			/* calculate gradient part; */
			sbulkparam->sfftv->gradient.val[i][0] = 
				-sbulkparam->model_tau	* sbulkparam->sfftv->rhoCplx.val[i][0] - 
				sbulkparam->model_gamma	* sbulkparam->sfftv->quadTerm.val[i][0] +
				sbulkparam->model_kappa * sbulkparam->sfftv->cubTerm.val[i][0];
			sbulkparam->sfftv->gradient.val[i][1] = 
				-sbulkparam->model_tau	* sbulkparam->sfftv->rhoCplx.val[i][1] - 
				sbulkparam->model_gamma	* sbulkparam->sfftv->quadTerm.val[i][1] + 
				sbulkparam->model_kappa	* sbulkparam->sfftv->cubTerm.val[i][1];

			/* semi-implicit scheme to update rhoCplx; */
			sbulkparam->sfftv->rhoCplx.val[i][0] -= stepSize * sbulkparam->sfftv->gradient.val[i][0];
			sbulkparam->sfftv->rhoCplx.val[i][1] -= stepSize * sbulkparam->sfftv->gradient.val[i][1];
			sbulkparam->sfftv->rhoCplx.val[i][0] /= (1.0 + stepSize*temp*sbulkparam->model_xi);
			sbulkparam->sfftv->rhoCplx.val[i][1] /= (1.0 + stepSize*temp*sbulkparam->model_xi);

//			if ( myrank == 0 )
//				printf("%d : temp = %.4e \t gradient = %.4e, %.4e \t rhoCplx = %.4e, %.4e.\n", i, temp,
//					sbulkparam->sfftv->gradient.val[i][0], sbulkparam->sfftv->gradient.val[i][1],
//					sbulkparam->sfftv->rhoCplx.val[i][0],  sbulkparam->sfftv->rhoCplx.val[i][1]);

			/* update gradient error; */
			sbulkparam->sfftv->gradient.val[i][0] += 
				sbulkparam->model_xi * temp * sbulkparam->sfftv->rhoCplx.val[i][0];
			sbulkparam->sfftv->gradient.val[i][1] += 
				sbulkparam->model_xi * temp * sbulkparam->sfftv->rhoCplx.val[i][1];
		}
		fn_complex_setZero ( sbulkparam->sfftv->rhoCplx.val[0]  );
		fn_complex_setZero ( sbulkparam->sfftv->gradient.val[0] );
		res = fn_tvec_maxAbs_complex ( sbulkparam->sfftv->gradient );

		iterator ++;
	}

	if ( (opt_iter % optSave) == 0 && itMax > 0 )
		fn_save_bulk_phase (sbulkparam, opt_iter, -1);		// save stationary value;

	// output data about the stationary state;
	fn_calc_bulk_energy ( sbulkparam, energy );
	hamilton = energy.val[2];
	if ( myrank == 0 )
	{
		fprintf(fdata, "%d\t%d\t%+.15E\t%+.15E\t%+.20E\t%+.15E\n", 
				opt_iter, iterator, energy.val[0], energy.val[1], energy.val[2], res);
		printf("\n\t ***** the bulk phase is %s\n", sbulkparam->phase);
		printf("\t ***** iterator = %d \t step size = %.3e\n", iterator, stepSize);
		printf("\t ***** res = % .10e \t diffham = % .10e\n", res, diffham);
		printf("\t ***** Laplace term: % .15e,\t Nonlinear term: % .15e\n", 
				energy.val[0], energy.val[1]);
		printf("\t ***** hamilton = % .20e\n\n", energy.val[2]);
	}

	fn_tvec_free<double> ( energy );
	return hamilton;
}


/**
 * \brief	Calculate the free energy of bulk phase;
 */
void fn_calc_bulk_energy (	stu_bulk_param		*sbulkparam, 
							tvec<double>		energy )
{
	int cplxDofs = sbulkparam->cplxDofs;

	/* initialization; */
	fn_tvec_setZero<double> ( energy );
	tvec<fftw_complex> DiffTerm		= fn_tvec_init<fftw_complex> ( cplxDofs );
	tvec<fftw_complex> DiffTermTmp	= fn_tvec_init<fftw_complex> ( cplxDofs );

	/* (q0^2-|G|^2) * (q1^2-|G|^2) * ... * rhoCplx; */
	double Difftmp;
	for (int i = 0; i < cplxDofs; i++)
	{
		Difftmp = 1.0;
		for (int j=0; j < sbulkparam->scale_num; j++)
		{
			Difftmp *= pow(sbulkparam->scale_val[j], 2) - sbulkparam->sfftv->Gsquare.val[i];  
		}
		DiffTermTmp.val[i][0] = Difftmp*sbulkparam->sfftv->rhoCplx.val[i][0];
		DiffTermTmp.val[i][1] = Difftmp*sbulkparam->sfftv->rhoCplx.val[i][1];
	}

	/* convolution calculation; */
	fn_serial_convolution(DiffTermTmp, DiffTerm, sbulkparam->sfftv, 2); // [(nabla^2+q0^2)*(nabla^2+q1^2)*...*rho]^2;
	fn_serial_convolution(sbulkparam->sfftv->rhoCplx, sbulkparam->sfftv->quadTerm, 
					sbulkparam->sfftv, 2); // rho^2;
	fn_serial_convolution(sbulkparam->sfftv->rhoCplx, sbulkparam->sfftv->cubTerm, 
					sbulkparam->sfftv, 3); // rho^3;
	fn_serial_convolution(sbulkparam->sfftv->rhoCplx, sbulkparam->sfftv->quarTerm, 
					sbulkparam->sfftv, 4); // rho^4;

	/* interaction energy; operator part; */
	energy.val[0] = DiffTerm.val[0][0]*sbulkparam->model_xi/2.0;
	/* entropy energy; nonlinear part; */
	energy.val[1] = 
		-sbulkparam->sfftv->quadTerm.val[0][0]	* sbulkparam->model_tau/2.0 - 
		sbulkparam->sfftv->cubTerm.val[0][0]	* sbulkparam->model_gamma/3.0 + 
		sbulkparam->sfftv->quarTerm.val[0][0]	* sbulkparam->model_kappa/4.0; 
	/* hamilton energy; */
	energy.val[2] = energy.val[0] + energy.val[1];

	fn_tvec_free<fftw_complex> ( DiffTerm );
	fn_tvec_free<fftw_complex> ( DiffTermTmp );
}


/**
 * \brief	The integration of functions about saving data;
 *			save Fourier coefficients, densities, plane wave;
 */
void fn_save_bulk_phase (	stu_bulk_param		*sbulkparam,
								int				opt_iter,
								int				iterator )
{	
	if ( myrank == 0 )
	{
		if ( sbulkparam->save_type[0] == 'y' || sbulkparam->save_type[0] == 'Y' )
			fn_disp_Fourier_coeff  ( sbulkparam->sfftv->rhoCplx, sbulkparam, opt_iter, iterator );
		if ( sbulkparam->save_type[1] == 'y' || sbulkparam->save_type[1] == 'Y' )
			fn_disp_bulk_density   ( sbulkparam->sfftv->rhoCplx, sbulkparam, opt_iter, iterator );
		if ( sbulkparam->save_type[2] == 'y' || sbulkparam->save_type[2] == 'Y' )
			fn_disp_bulk_plane_wave( sbulkparam->sfftv->rhoCplx, sbulkparam, opt_iter, iterator );
	}
}


/**
 * \brief	Optimize the computational box;
 */
void fn_opt_rcpBox		( stu_bulk_param		*sbulkparam )
{
	int bbType	 = sbulkparam->box_bbType;
	int	iter_max = sbulkparam->box_iter_max;
	double step0 = sbulkparam->box_step_size;
	double  tol	 = sbulkparam->box_tol;

	int optType	 = 0; // just for cuboid;

	/* memory allocation; */
	int row = sbulkparam->rcpBox.row;
	int col = sbulkparam->rcpBox.col;
	tmat<double> box0		= fn_tmat_init<double> ( row, col );
	tmat<double> box1		= fn_tmat_init<double> ( row, col );
	tmat<double> box		= fn_tmat_init<double> ( row, col );
	tmat<double> boxDiff	= fn_tmat_init<double> ( row, col );
	tmat<double> gradBOld	= fn_tmat_init<double> ( row, col );
	tmat<double> gradB		= fn_tmat_init<double> ( row, col );
	tmat<double> gradBDiff	= fn_tmat_init<double> ( row, col );
	fn_tmat_copy<double> ( box0,	sbulkparam->rcpBox );
	fn_tmat_copy<double> ( box,		sbulkparam->rcpBox );

	/* the step length for computing difference quotient; */
	double dh = 1.0e-4;

	/* initialization; */
	int		iter = 0;
	double  err	 = 1.0;
	double	step_size = step0;
	double  lham, rham;
	tvec<double> energy = fn_tvec_init<double> ( 3 );
	fn_tvec_setZero<double> ( energy );

	/* recurrence for updating the computational box; */
	while ( err > tol && iter < iter_max )
	{
		/* compute the gradient box; */
		fn_tmat_setZero<double> ( gradB );
		lham = 0.0, rham = 0.0;
		if ( strcmp ( sbulkparam->boxType, "cube" ) == 0 )
		{
			/* calculate the energy of boxOld + dh; */
			fn_tmat_copy<double> ( sbulkparam->rcpBox, box );
			for ( int j = 0; j < col; j++ )
				sbulkparam->rcpBox.val[j][j] += dh;
			fn_obt_rotateProjBox_matrix ( sbulkparam );
			fn_get_projection_plane	( sbulkparam );	// project plane; R'PBk;
			fn_get_Gsquare			( sbulkparam );	// Gsquare; |R'PBk|^2;
			fn_calc_bulk_energy ( sbulkparam, energy );
			rham = energy.val[2];

			/* calculate the energy of boxOld - dh; */
			for ( int j = 0; j < col; j++ )
				sbulkparam->rcpBox.val[j][j] -= 2.0 * dh;
			fn_obt_rotateProjBox_matrix ( sbulkparam );
			fn_get_projection_plane	( sbulkparam );	// project plane; R'PBk;
			fn_get_Gsquare			( sbulkparam );	// Gsquare; |R'PBk|^2;
			fn_calc_bulk_energy ( sbulkparam, energy );
			lham = energy.val[2];

			/* gradient box; */
			for ( int j = 0; j < col; j++ )
				gradB.val[j][j] = ( rham - lham ) / ( 2.0*dh );
		}
		else if ( strcmp ( sbulkparam->boxType, "cuboid" ) == 0 )
		{
			if ( optType == 0 )
			{
				/* calculate the energy of boxOld + dh; */
				fn_tmat_copy<double> ( sbulkparam->rcpBox, box );
				for ( int j = 0; j < col; j++ )
					sbulkparam->rcpBox.val[j][j] += dh * box.val[j][j];
				fn_obt_rotateProjBox_matrix ( sbulkparam );
				fn_get_projection_plane	( sbulkparam );	// project plane; R'PBk;
				fn_get_Gsquare			( sbulkparam );	// Gsquare; |R'PBk|^2;
				fn_calc_bulk_energy ( sbulkparam, energy );
				rham = energy.val[2];

				/* calculate the energy of boxOld - dh; */
				for ( int j = 0; j < col; j++ )
					sbulkparam->rcpBox.val[j][j] -= 2.0 * dh * box.val[j][j];
				fn_obt_rotateProjBox_matrix ( sbulkparam );
				fn_get_projection_plane	( sbulkparam );	// project plane; R'PBk;
				fn_get_Gsquare			( sbulkparam );	// Gsquare; |R'PBk|^2;
				fn_calc_bulk_energy ( sbulkparam, energy );
				lham = energy.val[2];

				/* gradient box; */
				for ( int j = 0; j < col; j++ )
					gradB.val[j][j] = ( rham - lham ) / ( 2.0*dh ) * box.val[j][j];
			}
			else
			{
				for ( int j = 0; j < col; j++ )
				{
					/* calculate the energy of boxOld + dh; */
					fn_tmat_copy<double> ( sbulkparam->rcpBox, box );
					double dhBox = sbulkparam->rcpBox.val[j][j] * dh;
					sbulkparam->rcpBox.val[j][j] += dhBox;
					fn_obt_rotateProjBox_matrix ( sbulkparam );
					fn_get_projection_plane	( sbulkparam );	// project plane; R'PBk;
					fn_get_Gsquare			( sbulkparam );	// Gsquare; |R'PBk|^2;
					fn_calc_bulk_energy ( sbulkparam, energy );
					rham = energy.val[2];

					/* calculate the energy of boxOld - dh; */
					sbulkparam->rcpBox.val[j][j] -= 2.0 * dhBox;
					fn_obt_rotateProjBox_matrix ( sbulkparam );
					fn_get_projection_plane	( sbulkparam );	// project plane; R'PBk;
					fn_get_Gsquare			( sbulkparam );	// Gsquare; |R'PBk|^2;
					fn_calc_bulk_energy ( sbulkparam, energy );
					lham = energy.val[2];

					/* gradient box; */
					gradB.val[j][j] = ( rham - lham ) / ( 2.0*dh );
				}
			}
		}
		else if ( strcmp ( sbulkparam->boxType, "hex" ) == 0 )
		{
			for ( int j = 0; j < col; j++ )
			{
				for ( int i = 0; i < j; i++ )
				{
					/* calculate the energy of boxOld + dh; */
					fn_tmat_copy<double> ( sbulkparam->rcpBox, box );
					sbulkparam->rcpBox.val[i][j] += dh;
					fn_obt_rotateProjBox_matrix ( sbulkparam );
					fn_get_projection_plane	( sbulkparam );	// project plane; R'PBk;
					fn_get_Gsquare			( sbulkparam );	// Gsquare; |R'PBk|^2;
					fn_calc_bulk_energy ( sbulkparam, energy );
					rham = energy.val[2];

					/* calculate the energy of boxOld + dh; */
					sbulkparam->rcpBox.val[i][j] -= 2.0 * dh;
					fn_obt_rotateProjBox_matrix ( sbulkparam );
					fn_get_projection_plane	( sbulkparam );	// project plane; R'PBk;
					fn_get_Gsquare			( sbulkparam );	// Gsquare; |R'PBk|^2;
					fn_calc_bulk_energy ( sbulkparam, energy );
					lham = energy.val[2];

					/* gradient box; */
					gradB.val[i][j] = ( rham - lham ) / ( 2.0*dh );
				}
			}
		}

		/* stepping size; */
		if ( iter == 0 )
		{
			step_size = step0;	// the initial stepping size;
		}
		else
		{
			/* 'box' - 'boxOld'; */
			fn_tmat_copy<double> ( boxDiff, box1 );
			fn_tmat_add<double>	 ( boxDiff, box0, 1.0, -1.0 );

			/* 'gradB' - 'gradBOld'; */
			fn_tmat_copy<double> ( gradBDiff, gradB );
			fn_tmat_add<double>  ( gradBDiff, gradBOld, 1.0, -1.0 );

			/* BB stepping size; */
			double bbTmp0, bbTmp1, bbStep;
			if ( bbType == 1 )
			{
				/* BB stepping size (type 1); */
				bbTmp0 = fn_tmat_inner ( boxDiff,   boxDiff );
				bbTmp1 = fn_tmat_inner ( boxDiff, gradBDiff );
			}
			else
			{
				/* BB stepping size (type 2); */
				bbTmp0 = fn_tmat_inner (   boxDiff, gradBDiff );
				bbTmp1 = fn_tmat_inner ( gradBDiff, gradBDiff );
			}
			step_size  = bbTmp0 / bbTmp1;
		}

		/* update computational box; */
		fn_tmat_add<double>  ( box, gradB, 1.0, -step_size );

		/* update; */
		if ( iter > 0 )
			fn_tmat_copy<double> ( box0, box1 );
		fn_tmat_copy<double> ( box1, box );
		fn_tmat_copy<double> ( gradBOld, gradB );

		/* calculate gradient error; */
		err = fn_tmat_maxAbs ( gradB );
		if ( myrank == 0 )
			printf("\t\t === > iter %d: step_size = %.5e \t res = % .15e\n", iter, step_size, err );
		iter ++ ;
	}

	if ( myrank == 0 )
	{
		printf("rcpBox: \n");
		fn_tmat_print<double> ( box );
		printf("\n");
	}

	/* update the optimized result; */
	fn_tmat_copy<double> ( sbulkparam->rcpBox, box );
	fn_obt_rotateProjBox_matrix ( sbulkparam );
	fn_get_projection_plane	( sbulkparam );	// project plane; R'PBk;
	fn_get_Gsquare			( sbulkparam );	// Gsquare; |R'PBk|^2;

	/* release memory; */
	fn_tmat_free<double> ( box0 );
	fn_tmat_free<double> ( box1 );
	fn_tmat_free<double> ( box );
	fn_tmat_free<double> ( boxDiff );
	fn_tmat_free<double> ( gradBOld );
	fn_tmat_free<double> ( gradB );
	fn_tmat_free<double> ( gradBDiff );
	fn_tvec_free<double> ( energy );
}


/**
 * \brief	Take translation in the bulk phase;
 */
void fn_translate_bulk_phase	(	stu_bulk_param		*sbulkparam )
{
	/* parameters; */
	int dimCpt	 = sbulkparam->dimCpt;
	int cplxDofs = sbulkparam->cplxDofs;

	/* take translation in the bulk phase; */
	double translVal;
	double tmpReal, tmpCplx;
	for (int i = 0; i < cplxDofs; i++)
	{
		/* consider translation; calculate (PBk)'t = k'B'P't; superscript ' means transpose; */
		/* t'PB * k; */
		translVal = 0.0;
		for (int j = 0; j < dimCpt; j++)
		{
			translVal += sbulkparam->translProjBoxVec.val[j] * 
				sbulkparam->sfftv->indKspace.val[i][j]; // t'PB * k;
		}

		/* rhoCplx * exp(i(PBk)'t); t is translation matrix; */
		tmpReal = sbulkparam->sfftv->rhoCplx.val[i][0];
		tmpCplx = sbulkparam->sfftv->rhoCplx.val[i][1];
		sbulkparam->sfftv->rhoCplx.val[i][0] = 
			tmpReal * cos(translVal) - tmpCplx * sin(translVal);
		sbulkparam->sfftv->rhoCplx.val[i][1] =
			tmpReal * sin(translVal) + tmpCplx * cos(translVal);
	}

	/*
	if ( myrank == 0 )
	{
		printf("translProjBoxVec:\n");
		for ( int j = 0; j < dimCpt; j++ )
			printf("%.10e\n", sbulkparam->translProjBoxVec.val[j]);
	}
	*/
}


/**
 * \brief	Truncation to reduce the computational cost;
 */
void fn_truncation_bulk_phase	(	stu_bulk_param		*sbulkparam,
										double			hamilton	)
{

	/* truncation; */
	double trunc_tol = sbulkparam->trunc_tol; // the tolerance for truncation;
	int cplxDofs = sbulkparam->cplxDofs;
	int dimPhy	 = sbulkparam->dimPhy;
	int dimCpt	 = sbulkparam->dimCpt;

	/* memory allocation for temporary variables; */
	tvec<fftw_complex>  rhoCplx   = fn_tvec_init<fftw_complex>	( cplxDofs );
	tmat<int>			indKspace = fn_tmat_init<int>			( cplxDofs, dimCpt );
	tmat<double>		projPlane = fn_tmat_init<double>		( cplxDofs, dimPhy );

	int ind = 0;
	double trunc_err = 0;
	for ( int j1 = 0; j1 < cplxDofs; j1++ )
	{
		double val = fn_complex_abs ( sbulkparam->sfftv->rhoCplx.val[j1] );
		if ( val > trunc_tol )
		{
			rhoCplx.val[ind][0] = sbulkparam->sfftv->rhoCplx.val[j1][0];
			rhoCplx.val[ind][1] = sbulkparam->sfftv->rhoCplx.val[j1][1];
			for ( int j2 = 0; j2 < dimCpt; j2++ )
				indKspace.val[ind][j2] = sbulkparam->sfftv->indKspace.val[j1][j2];
			for ( int j2 = 0; j2 < dimPhy; j2++ )
				projPlane.val[ind][j2] = sbulkparam->sfftv->projPlane.val[j1][j2];
			ind ++ ;
		}
		else
			trunc_err += pow(val,2);
	}
	sbulkparam->cplxDofs = ind + 1; // update cplxDofs but save rhoCplx.val[0][0] for mass;
	trunc_err = sqrt(trunc_err);
	if ( myrank == 0 )
		printf("truncation length is %d and error is %.2e.\n", sbulkparam->cplxDofs, trunc_err);

	/* obtain index to sort descend; */
	double tmp;
	vector <stu_sort> sort_array ( ind );
	for ( int i = 0; i < ind; i++ )
	{
		sort_array[i].ind = i;
		sort_array[i].val = fn_complex_abs ( rhoCplx.val[i] );
	}
	sort(sort_array.begin(), sort_array.end(), fn_compare_descend);

	/* release memory and memory allocation again for sbulkparam; */
	fn_tvec_free<fftw_complex>  ( sbulkparam->sfftv->rhoCplx );	
	fn_tmat_free<int>			( sbulkparam->sfftv->indKspace );
	fn_tmat_free<double>		( sbulkparam->sfftv->projPlane );
	sbulkparam->sfftv->rhoCplx	 = fn_tvec_init<fftw_complex>	( sbulkparam->cplxDofs );
	sbulkparam->sfftv->indKspace = fn_tmat_init<int>			( sbulkparam->cplxDofs, dimCpt );
	sbulkparam->sfftv->projPlane = fn_tmat_init<double>			( sbulkparam->cplxDofs, dimPhy );
	sbulkparam->sfftv->ind		 = fn_tvec_init<int>			( sbulkparam->cplxDofs );

	/* assign values for [0]; */
	sbulkparam->sfftv->ind.val[0] = 0;
	sbulkparam->sfftv->rhoCplx.val[0][0] = 0.0;
	sbulkparam->sfftv->rhoCplx.val[0][1] = 0.0;
	for ( int j2 = 0; j2 < dimCpt; j2++ )
		sbulkparam->sfftv->indKspace.val[0][j2] = 0;
	for ( int j2 = 0; j2 < dimPhy; j2++ )
		sbulkparam->sfftv->projPlane.val[0][j2] = 0.0;

	/* assign values for others; */
	for ( int j1 = 1; j1 < sbulkparam->cplxDofs; j1++ )
	{
		int ind0 = sort_array[j1-1].ind;
		sbulkparam->sfftv->ind.val[j1] = ind0;
		sbulkparam->sfftv->rhoCplx.val[j1][0] = rhoCplx.val[ind0][0];
		sbulkparam->sfftv->rhoCplx.val[j1][1] = rhoCplx.val[ind0][1];
		for ( int j2 = 0; j2 < dimCpt; j2++ )
			sbulkparam->sfftv->indKspace.val[j1][j2] = indKspace.val[ind0][j2];
		for ( int j2 = 0; j2 < dimPhy; j2++ )
			sbulkparam->sfftv->projPlane.val[j1][j2] = projPlane.val[ind0][j2];
	}

	/* release memory; */
	fn_tvec_free<fftw_complex>  ( rhoCplx );	
	fn_tmat_free<int>			( indKspace );
	fn_tmat_free<double>		( projPlane );

	/* save the truncation data for bulk phase; */
	if ( myrank == 0 )
	{
		if ( strcmp( main_type, "stable_bulk1" ) == 0 ||
			 strcmp( main_type, "stable_bulk2" ) == 0 ) // stable bulk phases;
			fn_save_trunc_bulk_phase ( sbulkparam, hamilton, false );
		else
			fn_save_trunc_bulk_phase ( sbulkparam, hamilton, true );
	}
}


/**
 * \brief	obtain the serial number of model parameters for bulk phase 
 *			or save the model parameters for bulk phase;
 */
int fn_paraOrder_bulk_phase	(	stu_bulk_param		*sbulkparam )
{
	char paraBulkFile[FILELEN];
	sprintf(paraBulkFile, "./bulk/%s_%s/paralist.dat", model_type, sbulkparam->phase);

	/* search the serial number of model parameters; */
	int paraOrder = -1;
	int nline = 0;
	int Four_num;
	double tmp_xi, tmp_tau, tmp_gamma, tmp_kappa, tmp_trunc_tol;
	if ( access(paraBulkFile, F_OK) == 0 )
	{
		char    buffer[STRLEN]; // Note: max number of char for each line!
		int     val;
		int		status = 1;		// 1: success; 0: failure;
		FILE *fp = fopen(paraBulkFile, "r");

		while ( status == 1 )
		{
			val = fscanf(fp, "%s", buffer);
			if (val==EOF) break;
			if (val!=1) { status = 0; break; }
			if (buffer[0]=='#' || buffer[0]=='%' || buffer[0]=='|' || buffer[0]=='/')
			{
				if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
				continue;
			}

			val = fscanf(fp, "%lf", &(tmp_xi));
			val = fscanf(fp, "%lf", &(tmp_tau));
			val = fscanf(fp, "%lf", &(tmp_gamma));
			val = fscanf(fp, "%lf", &(tmp_kappa));
			val = fscanf(fp, "%d", &(Four_num));
			val = fscanf(fp, "%lf", &(tmp_trunc_tol));
			if ( abs(sbulkparam->model_xi - tmp_xi) < ZEROTOL && 
				 abs(sbulkparam->model_tau - tmp_tau) < ZEROTOL &&
				 abs(sbulkparam->model_gamma - tmp_gamma) < ZEROTOL && 
				 abs(sbulkparam->model_kappa - tmp_kappa) < ZEROTOL &&
				 abs(sbulkparam->Four_num - Four_num) < ZEROTOL &&
				 abs(sbulkparam->trunc_tol - tmp_trunc_tol) < ZEROTOL )
			{
				paraOrder = nline;
				status = 0;		break;
			}
			nline ++ ;
		}
		fclose(fp);	
	}
	else
	{
		paraOrder = 0;
		FILE *fp = fopen(paraBulkFile, "w");
		fprintf(fp, "para\t%+.10E\t%+.10E\t%+.10E\t%+.10E\t%d\t%+.10E\n", sbulkparam->model_xi, 
				sbulkparam->model_tau, sbulkparam->model_gamma, sbulkparam->model_kappa, 
				sbulkparam->Four_num, sbulkparam->trunc_tol);
		fclose(fp);
	}

	if ( paraOrder == -1 ) // paraBulkFile exist but not save the current parameters;
	{
		paraOrder = nline;
		FILE *fp = fopen(paraBulkFile, "a+");
		fprintf(fp, "para\t%+.10E\t%+.10E\t%+.10E\t%+.10E\t%d\t%+.10E\n", sbulkparam->model_xi, 
				sbulkparam->model_tau, sbulkparam->model_gamma, sbulkparam->model_kappa, 
				sbulkparam->Four_num, sbulkparam->trunc_tol);
		fclose(fp);
	}
	return paraOrder;
}


/**
 * \brief	save the truncation data for bulk phase;
 */
void fn_save_trunc_bulk_phase	(	stu_bulk_param		*sbulkparam,
										double			hamilton,
										bool			ischeck )
{
	/* obtain the serial number of model parameters for bulk phase; */
	int paraOrder = fn_paraOrder_bulk_phase	( sbulkparam );
	if ( myrank == 0 )
		printf("the serial number of model parameters is %d.\n", paraOrder);

	/* file for saving; */
	char bulkFile[FILELEN];
	sprintf(bulkFile, "./bulk/%s_%s/para%d.dat", model_type, sbulkparam->phase, paraOrder);

	/* save data; */
	if ( ischeck && access(bulkFile, F_OK) == 0 )
	{
		printf("%s exist.\n", bulkFile);
	}
	else
	{
		FILE *fdata = fopen(bulkFile, "w");
		fprintf(fdata, "%d\t%d\t%d\t%+.15E\n", sbulkparam->cplxDofs, sbulkparam->dimCpt, sbulkparam->dimPhy, hamilton);
		for ( int j1 = 0; j1 < sbulkparam->cplxDofs; j1++ )
		{
			fprintf(fdata, "%+.15E\t%+.15E", sbulkparam->sfftv->rhoCplx.val[j1][0], sbulkparam->sfftv->rhoCplx.val[j1][1]);
			for ( int j2 = 0; j2 < sbulkparam->dimCpt; j2++ )
				fprintf(fdata, "\t%d", sbulkparam->sfftv->indKspace.val[j1][j2]);
			fprintf(fdata, "\t%d", sbulkparam->sfftv->ind.val[j1]);
			fprintf(fdata, "\n");
		}
		fclose(fdata);
	}
}


/**
 * \brief	load the truncation data for bulk phase;
 */
double fn_load_trunc_bulk_phase	(	stu_bulk_param		*sbulkparam	)
{
	/* obtain the serial number of model parameters for bulk phase; */
	int paraOrder = fn_paraOrder_bulk_phase	( sbulkparam );

	/* file for saving; */
	char bulkFile[FILELEN];
	sprintf(bulkFile, "./bulk/%s_%s/para%d.dat", model_type, sbulkparam->phase, paraOrder);

	/* load data; */
	double hamilton = 0.0;
	if ( access(bulkFile, F_OK) == 0 )
	{
		int		val;
		FILE *fp = fopen(bulkFile, "r");
		val = fscanf(fp, "%d", &(sbulkparam->cplxDofs));
		val = fscanf(fp, "%d", &(sbulkparam->dimCpt));
		val = fscanf(fp, "%d", &(sbulkparam->dimPhy));
		val = fscanf(fp, "%lf", &(hamilton));

		/* memory allocation for truncation data; */
		sbulkparam->sfftv->rhoCplx	 = fn_tvec_init<fftw_complex>	( sbulkparam->cplxDofs );
		sbulkparam->sfftv->indKspace = fn_tmat_init<int>			( sbulkparam->cplxDofs, sbulkparam->dimCpt );
		sbulkparam->sfftv->projPlane = fn_tmat_init<double>			( sbulkparam->cplxDofs, sbulkparam->dimPhy );
		sbulkparam->sfftv->ind		 = fn_tvec_init<int>			( sbulkparam->cplxDofs );

		for ( int j1 = 0; j1 < sbulkparam->cplxDofs; j1++ )
		{
			val = fscanf(fp, "%lf", &(sbulkparam->sfftv->rhoCplx.val[j1][0]));
			val = fscanf(fp, "%lf", &(sbulkparam->sfftv->rhoCplx.val[j1][1]));
			for ( int j2 = 0; j2 < sbulkparam->dimCpt; j2++ )
				val = fscanf(fp, "%d", &(sbulkparam->sfftv->indKspace.val[j1][j2]));
			val = fscanf(fp, "%d", &(sbulkparam->sfftv->ind.val[j1]));
			for ( int j2 = 0; j2 < sbulkparam->dimPhy; j2++ )
			{
				double tmp = 0.0;
				for ( int j3 = 0; j3 < sbulkparam->dimCpt; j3++ )
					tmp += sbulkparam->rotateProjBoxMat.val[j2][j3] * sbulkparam->sfftv->indKspace.val[j1][j3];
				sbulkparam->sfftv->projPlane.val[j1][j2] = tmp;
			}
		}
		fclose(fp);
		if ( myrank == 0 )
		{
			printf("truncation length is %d.\n", sbulkparam->cplxDofs);
			printf("convergence hamilton is %+.15E.\n", hamilton);
		}
	}
	else
	{
		if ( myrank == 0 )
			printf("%s not exist.\n", bulkFile);
	}
	return hamilton;
}
