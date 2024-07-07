/*! \file	CommonFourGJP.cpp
 *
 * \brief	Transform rhoJCplx (in the space of Fourier and GJPs) to 
 *			the common space (Fourier and GJPs);
 */

#include "Data.h"
#include "DataOperators.h"
#include "Head.h"
#include "Mytimer.h"
#include "functs.h"


/**
 * \brief	Initialization and preparation;
 *			Transform rhoJCplx of the two bulk phases in a common space;	
 *
 * \param	sbulkparam:		the structure body for bulk phases;
 * \param	ssysparam:		the structure body for the interface system;
 *
 */
void fn_common_Fourier_GJP (	stu_bulk_param		*sbulkparam1,
								stu_bulk_param		*sbulkparam2,
								stu_system_param	*ssysparam )
{
	if ( myrank == 0 )
		printf(" <========== Transform to the common space ==========> \n\n");
	mytimer_t timer;
	timer.reset();
	timer.start();

	/* obtain the projection plane about the common 'rotateProjBoxMat'; */
	fn_get_system_projection_plane ( ssysparam );

//	sbulkparam1->srebndv->isSave  = true;
//	sbulkparam2->srebndv->isSave  = true;
//	sbulkparam1->srebndv->isTest  = true;
//	sbulkparam2->srebndv->isTest  = true;

	/**
	 * re-represent by the common 'rotateProjBoxMat';
	 *	including 'rhoJCplx', 'd0bnd', ...;
	 */
	fn_total_re_represent_common ( sbulkparam1, ssysparam );
	fn_total_re_represent_common ( sbulkparam2, ssysparam );

	/* release memory about 'sbndv'; */
	fn_bnd_memory_free ( sbulkparam1->sbndv, "total" );
	fn_bnd_memory_free ( sbulkparam2->sbndv, "total" );

	MPI_Barrier ( MPI_COMM_WORLD );
	timer.pause();
	if ( myrank == 0 )
	{
		printf("\n\t ***** time cost of transformation to the common space: ");
		printf("%f seconds *****\n\n", timer.get_current_time());
	}
}


/**
 * \brief	Get projection plane in the interface system; i.e. R'PBk;
 *			R is rotateMat;
 *			P is projMat;
 *			B is rcpBox;
 *			k is Fourier index;
 *			R'PB: [dimRePhy, dimReCpt];		k: [dimReCpt, 1];
 */
void fn_get_system_projection_plane ( stu_system_param		*ssysparam )
{
	/* parameters; */
	int dimRePhy	  = ssysparam->dimRePhy;
	int dimReCpt	  = ssysparam->dimReCpt;
	int Four_num	  = ssysparam->Four_num;
	int cplxReDofs    = 1;
	ssysparam->NCpt = fn_tvec_init<ptrdiff_t> ( dimReCpt );
	for ( int i = 0; i < dimReCpt; i++ )
	{
		ssysparam->NCpt.val[i] = Four_num;
		cplxReDofs *= Four_num;
	}
	ssysparam->cplxReDofs = cplxReDofs;

	/* local memory of mpi processes; */
    alloc_local_sys = fftw_mpi_local_size( dimReCpt, ssysparam->NCpt.val, 
				MPI_COMM_WORLD, &local_n0_sys, &local_0_start_sys);
	ssysparam->sfftv->indKspace = fn_tmat_init<int>		( alloc_local_sys, dimReCpt );
	ssysparam->sfftv->projPlane	= fn_tmat_init<double>	( alloc_local_sys, dimRePhy );

	/* offset and receive count; */
    displs_sys    = (int *) malloc ( sizeof(int) * nprocs ); // offset of each process;
    recvCount_sys = (int *) malloc ( sizeof(int) * nprocs ); // offset of each process;
	int sendcount;
	sendcount = alloc_local_sys;
	MPI_Allgather ( &sendcount, 1, MPI_INT, recvCount_sys, 1, MPI_INT, MPI_COMM_WORLD );
	displs_sys[0] = 0;
	for ( int i = 1; i < nprocs; i++ )
		displs_sys[i] = displs_sys[i-1] + recvCount_sys[i-1];
	if ( myrank == 0 )
	{
		printf("recvCount_sys: ");
		for ( int i = 0; i < nprocs; i++ )
			printf("%d ", recvCount_sys[i]);
		printf("\n");
	}

	/* obtain k index; */
	fn_obt_kIndex ( ssysparam->sfftv->indKspace, ssysparam->NCpt, displs_sys, alloc_local_sys );
//	char fn[FILELEN];
//	sprintf( fn, "k%d.txt", myrank );
//	fn_tmat_save_int  ( ssysparam->sfftv->indKspace, fn );

	/* obtain projection wave; */
	double RPBKtmp;
	for (int i = 0; i < alloc_local_sys; i++)
	{
		for (int kk = 0; kk < dimRePhy; kk++)
		{
			RPBKtmp = 0.0;
			for (int jj = 0; jj < dimReCpt; jj++) // R'PB * k;
				RPBKtmp += ssysparam->rotateProjBoxMat.val[kk][jj] * 
							ssysparam->sfftv->indKspace.val[i][jj];
			ssysparam->sfftv->projPlane.val[i][kk] = RPBKtmp;
		}
	}
}


/**
 * \brief	summarization of 'fn_re_represent_common' or 'fn_new_re_represent_common';
 */
void fn_total_re_represent_common (		stu_bulk_param		*sbulkparam,
									   stu_system_param		*ssysparam )
{
	int type = 0;
	if ( myrank == 0 )
	{
		if ( sbulkparam->sflag == 1 )
			printf("re-represent data of left bulk phase.\n");
		else if ( sbulkparam->sflag == 2 )
			printf("re-represent data of right bulk phase.\n");
	}
	/* memory allocation; */
	fn_common_memory_alloc (sbulkparam->sbndv, sbulkparam->srebndv, 
							ssysparam->cplxReDofs, "total");

	if ( type == 0 )
	{
		/* re-represent 'rhoJCplx' by the common 'rotateProjBoxMat'; */
		fn_re_represent_common (sbulkparam, ssysparam, sbulkparam->sbndv->rhoJCplx, 
								sbulkparam->srebndv->rhoJCplx, "rhoReJCplx");

		/* re-represent 'd0bnd', ..., by the common 'rotateProjBoxMat'; */
		fn_re_represent_common (sbulkparam, ssysparam, sbulkparam->sbndv->d0bnd, 
								sbulkparam->srebndv->d0bnd, "d0rebnd");
		fn_re_represent_common (sbulkparam, ssysparam, sbulkparam->sbndv->d1bnd, 
								sbulkparam->srebndv->d1bnd, "d1rebnd");
		fn_re_represent_common (sbulkparam, ssysparam, sbulkparam->sbndv->d2bnd, 
								sbulkparam->srebndv->d2bnd, "d2rebnd");
		if ( strcmp(model_type, "LP") == 0 )
		{
			fn_re_represent_common (sbulkparam, ssysparam, sbulkparam->sbndv->d3bnd, 
									sbulkparam->srebndv->d3bnd, "d3rebnd");
			fn_re_represent_common (sbulkparam, ssysparam, sbulkparam->sbndv->d4bnd, 
									sbulkparam->srebndv->d4bnd, "d4rebnd");
			fn_re_represent_common (sbulkparam, ssysparam, sbulkparam->sbndv->d5bnd, 
									sbulkparam->srebndv->d5bnd, "d5rebnd");
			fn_re_represent_common (sbulkparam, ssysparam, sbulkparam->sbndv->d6bnd, 
									sbulkparam->srebndv->d6bnd, "d6rebnd");
		}
	}
	else
	{
		/* re-represent 'rhoJCplx' by the common 'rotateProjBoxMat'; */
		fn_new_re_represent_common (sbulkparam, ssysparam, sbulkparam->sbndv->rhoJCplx, 
								sbulkparam->srebndv->rhoJCplx, "rhoReJCplx");

		/* re-represent 'd0bnd', ..., by the common 'rotateProjBoxMat'; */
		fn_new_re_represent_common (sbulkparam, ssysparam, sbulkparam->sbndv->d0bnd, 
								sbulkparam->srebndv->d0bnd, "d0rebnd");
		fn_new_re_represent_common (sbulkparam, ssysparam, sbulkparam->sbndv->d1bnd, 
								sbulkparam->srebndv->d1bnd, "d1rebnd");
		fn_new_re_represent_common (sbulkparam, ssysparam, sbulkparam->sbndv->d2bnd, 
								sbulkparam->srebndv->d2bnd, "d2rebnd");
		if ( strcmp(model_type, "LP") == 0 )
		{
			fn_new_re_represent_common (sbulkparam, ssysparam, sbulkparam->sbndv->d3bnd, 
									sbulkparam->srebndv->d3bnd, "d3rebnd");
			fn_new_re_represent_common (sbulkparam, ssysparam, sbulkparam->sbndv->d4bnd, 
									sbulkparam->srebndv->d4bnd, "d4rebnd");
			fn_new_re_represent_common (sbulkparam, ssysparam, sbulkparam->sbndv->d5bnd, 
									sbulkparam->srebndv->d5bnd, "d5rebnd");
			fn_new_re_represent_common (sbulkparam, ssysparam, sbulkparam->sbndv->d6bnd, 
									sbulkparam->srebndv->d6bnd, "d6rebnd");
		}
	}

	/* save density calculating by re-represented rhoJCplx; */
	if ( sbulkparam->srebndv->isSave )
	{
		fn_disp_bulk_common_density (sbulkparam->srebndv->rhoJCplx, sbulkparam, ssysparam);
		/* calculate error; */
		if ( myrank == 0 )
			fn_proj_common_error ( sbulkparam, ssysparam );
	}

	/* compute the mass; */
	fftw_complex mass;
	fn_obt_bulk_mass_common ( sbulkparam, ssysparam, mass );
	if ( myrank == 0 )
	{
		if ( sbulkparam->sflag == 1 )
			printf("\n\t ---> Mass of re-representing left  bulk phase: %.10e\n", fn_complex_abs(mass));
		else if ( sbulkparam->sflag == 2 )
			printf("\n\t ---> Mass of re-representing right bulk phase: %.10e\n", fn_complex_abs(mass));
	}

	/* valid data; */
	int reComLen;
	if ( myrank == nprocs-1 )
		reComLen = ssysparam->cplxReDofs - displs_sys[myrank];
	else
		reComLen = alloc_local_sys;

	/* skip [0] to maintain mass conservation; error examination; */
	/*
	if ( myrank == 0 )
		printf("\t ---> Set mass %.10e equal to zero in CommonFourGJP.cpp\n", fn_complex_abs(mass));
	for ( int j1 = 0; j1 < sbulkparam->srebndv->nd; j1++ )
	{
		int ind0 = j1*reComLen;	// the first element along Fourier direction;
		fn_complex_setZero ( sbulkparam->srebndv->rhoJCplx.val[ind0] );
	}
	for ( int j1 = 0; j1 < sbulkparam->srebndv->xlen; j1++ )
	{
		int ind0 = j1*reComLen;	// the first element along Fourier direction;
		fn_complex_setZero ( sbulkparam->srebndv->d0bnd.val[ind0] );
		fn_complex_setZero ( sbulkparam->srebndv->d1bnd.val[ind0] );
		fn_complex_setZero ( sbulkparam->srebndv->d2bnd.val[ind0] );
		if ( strcmp(model_type, "LP") == 0 )
		{
			fn_complex_setZero ( sbulkparam->srebndv->d3bnd.val[ind0] );
			fn_complex_setZero ( sbulkparam->srebndv->d4bnd.val[ind0] );
			fn_complex_setZero ( sbulkparam->srebndv->d5bnd.val[ind0] );
			fn_complex_setZero ( sbulkparam->srebndv->d6bnd.val[ind0] );
		}
	}
	*/
}


/**
 * \brief	memory allocation;
 *
 * \param	sbndv:		provide parameters;
 * \param	sobjbndv:	the objective structure body;
 * \param	cplxReDofs:	a parameter provided by 'ssysparam';
 */
void fn_common_memory_alloc (	stu_bnd_var			*sbndv,
								stu_bnd_var			*sobjbndv,
									int				cplxReDofs,
								const char			*memType)
{
	sobjbndv->polyDegree = sbndv->polyDegree;
	sobjbndv->nd		 = sbndv->nd;
	sobjbndv->xlen		 = sbndv->xlen;
	sobjbndv->cplxDofs	 = cplxReDofs;
	sobjbndv->x_range	 = sbndv->x_range;

	int rhoLen = sobjbndv->nd   * alloc_local_sys;
	sobjbndv->rhoJCplx	= fn_tvec_init<fftw_complex> ( rhoLen );

	int bndLen = sobjbndv->xlen * alloc_local_sys;
	sobjbndv->d0bnd		= fn_tvec_init<fftw_complex> ( bndLen );
	sobjbndv->d2bnd		= fn_tvec_init<fftw_complex> ( bndLen );
	if ( strcmp(model_type, "LP") == 0 )
	{
		sobjbndv->d4bnd		= fn_tvec_init<fftw_complex> ( bndLen );
		sobjbndv->d6bnd		= fn_tvec_init<fftw_complex> ( bndLen );
	}

	if ( strcmp(memType, "total") == 0 )
	{
		sobjbndv->d1bnd		= fn_tvec_init<fftw_complex> ( bndLen );
		if ( strcmp(model_type, "LP") == 0 )
		{
			sobjbndv->d3bnd		= fn_tvec_init<fftw_complex> ( bndLen );
			sobjbndv->d5bnd		= fn_tvec_init<fftw_complex> ( bndLen );
		}
	}
}


/**
 * \brief	Re-represent 'rhoJCplx' in the common 'rotateProjBoxMat';
 *			Re-represent 'd0bnd', ... in the common 'rotateProjBoxMat';
 *
 * \param	sbulkparam:		provide projection plane about bulk phases;
 * \param	ssysparam:		provide projection plane about the interface system;
 * \param	origBulk:		the original bulk phase (projected in the system);
 * \param	reComBulk:		the bulk phase transformed by the common 'rotateProjBoxMat';
 * \param	saveStr:		a string for saving data;
 *
 * \note!	The essential spectra are consistent;
 */
void fn_re_represent_common	(	stu_bulk_param		*sbulkparam,
								stu_system_param	*ssysparam,
								tvec<fftw_complex>	origBulk,
								tvec<fftw_complex>	reComBulk,
									const char		*saveStr )
{
	if ( myrank == 0 )
		printf("re-represent %s.\n", saveStr);
	/* load parameters; */
	double tol = 1e-10;
	int nd		   = sbulkparam->sbndv->nd; 
	int xlen	   = sbulkparam->sbndv->xlen;
	int cplxDofs   = sbulkparam->sbndv->cplxDofs; // the number of spectra in bulk phase;
	int dimRePhy   = ssysparam->dimRePhy;
	int dimReCpt   = ssysparam->dimReCpt;
	int cplxReDofs = ssysparam->cplxReDofs;		 // the number of spectra for re-representation;
	int ind0, ind1;
	int sendcount;

	/* considering two cases: 'rhoJCplx'; 'd0bnd',...; */
	int nj;
	if ( strcmp(saveStr, "rhoReJCplx") == 0 ) nj = nd;
	else nj = xlen;

	/* initialization of results; */
	int lenRe = nj * alloc_local_sys;
	for ( int i = 0; i < lenRe; i++ )
		fn_complex_setZero ( reComBulk.val[i] );

	/* valid data; */
	int reComLen;
	if ( myrank == nprocs-1 )
		reComLen = cplxReDofs - displs_sys[myrank];
	else
		reComLen = alloc_local_sys;

	/* re-representation of 'rhoJCplx'; */
	if ( myrank == 0 )
		printf("progress: ");
	int pmod = reComLen / 10;
	double res, tmp0, tmp1;
	for ( int j = 0; j < reComLen; j++ )
	{
		for ( int i = 0; i < cplxDofs; i++ )
		{
			res = 0.0;
			for ( int j0 = 0; j0 < dimRePhy; j0++ )	// calculate residual; norm 2;
			{
				tmp0 = sbulkparam->sfftv->projPlane.val[i][j0+1];
				tmp1 = ssysparam->sfftv->projPlane.val[j][j0];
				res += pow(tmp0-tmp1, 2);
			}
			res = sqrt(res);
			if ( res < tol )
			{
				for ( int j1 = 0; j1 < nj; j1++ )	// same spectra;
				{
					ind0 = j1 * reComLen + j;
					ind1 = j1 * cplxDofs + i;
					reComBulk.val[ind0][0] += origBulk.val[ind1][0];
					reComBulk.val[ind0][1] += origBulk.val[ind1][1];
				}
			}
		}
		if ( myrank == 0 )
			if ( j % pmod == 0 )
				printf("%.2lf ", (double) j/reComLen);
	}
	if ( myrank == 0 )
		printf("finish\n");

	/* save results; */
	if ( sbulkparam->srebndv->isTest )
	{
	}
}


/**
 * \brief	Re-represent 'rhoJCplx' in the common 'rotateProjBoxMat';
 *			Re-represent 'd0bnd', ... in the common 'rotateProjBoxMat';
 *
 * \param	sbulkparam:		provide projection plane about bulk phases;
 * \param	ssysparam:		provide projection plane about the interface system;
 * \param	origBulk:		the original bulk phase (projected in the system);
 * \param	reComBulk:		the bulk phase transformed by the common 'rotateProjBoxMat';
 * \param	saveStr:		a string for saving data;
 *
 * \note!	The essential spectra are consistent;
 */
void fn_new_re_represent_common	(	stu_bulk_param		*sbulkparam,
									stu_system_param	*ssysparam,
									tvec<fftw_complex>	origBulk,
									tvec<fftw_complex>	reComBulk,
										const char		*saveStr )
{
	if ( myrank == 0 )
		printf("re-represent %s.\n", saveStr);
	/* load parameters; */
	double tol = 1e-10;
	int nd		   = sbulkparam->sbndv->nd; 
	int xlen	   = sbulkparam->sbndv->xlen;
	int cplxDofs   = sbulkparam->sbndv->cplxDofs; // the number of spectra in bulk phase;
	int dimCpt	   = sbulkparam->dimCpt;
	int dimRePhy   = ssysparam->dimRePhy;
	int dimReCpt   = ssysparam->dimReCpt;
	int cplxReDofs = ssysparam->cplxReDofs;		 // the number of spectra for re-representation;
	int ind0, ind1;
	int sendcount;

	/* load coefficient matrix; */
	tmat<int> coeffmat	= fn_tmat_init<int> ( dimCpt, ssysparam->coeffmat.col );
	tvec<int> coeffrslt = fn_tvec_init<int>	( coeffmat.col );
	int bulk1len = ssysparam->coeffmat.row - dimCpt; 
	for ( int i = 0; i < coeffmat.row; i++ )
	{
		for ( int j = 0; j < coeffmat.col; j++ )
		{
			if ( sbulkparam->sflag == 1 )
				coeffmat.val[i][j] = ssysparam->coeffmat.val[i][j];
			else if ( sbulkparam->sflag == 2 )
				coeffmat.val[i][j] = ssysparam->coeffmat.val[bulk1len+i][j];
		}
	}

	/* considering two cases: 'rhoJCplx'; 'd0bnd',...; */
	int nj;
	if ( strcmp(saveStr, "rhoReJCplx") == 0 ) nj = nd;
	else nj = xlen;

	/* memory allocation for temporal variables; */
	int lenReAll = nj * cplxReDofs;

	/* valid data; */
	int reComLen;
	if ( myrank == nprocs-1 )
		reComLen = cplxReDofs - displs_sys[myrank];
	else
		reComLen = alloc_local_sys;

	/* initialization of results; */
	int lenRe = nj * reComLen;
	for ( int i = 0; i < lenRe; i++ )
		fn_complex_setZero ( reComBulk.val[i] );

	/* re-representation of 'rhoJCplx'; */
	if ( myrank == 0 )
		printf("progress: ");
	int pmod = cplxDofs / 10;
	for ( int i = 0; i < cplxDofs; i++ )
	{
		/* coeffmat * indKspace; */
		for ( int j1 = 0; j1 < coeffmat.col; j1++ )
		{
			int tmpVal = 0;
			for ( int j2 = 0; j2 < dimCpt; j2++ )
				tmpVal += coeffmat.val[j2][j1] * sbulkparam->sfftv->indKspace.val[i][j2];
			coeffrslt.val[j1] = tmpVal;
		}

		/* local index transform to global index; */
		int kk = fn_kIndex_to_gIndex ( coeffrslt.val, ssysparam->NCpt );
		int kkdispls = kk - displs_sys[myrank];

		/* assign to represent value; */
		if ( kkdispls < reComLen && kkdispls >= 0 )
		{
			for ( int j1 = 0; j1 < nj; j1++ )	// same spectra;
			{
				ind0 = j1 * reComLen + kkdispls;
				ind1 = j1 * cplxDofs + i;
				reComBulk.val[ind0][0] += origBulk.val[ind1][0];
				reComBulk.val[ind0][1] += origBulk.val[ind1][1];
			}
		}

		if ( myrank == 0 )
			if ( i % pmod == 0 )
				printf("%.2lf ", (double) i/cplxDofs);
	}

	if ( myrank == 0 )
		printf("finish\n");

	/* save results; */
	if ( sbulkparam->srebndv->isTest )
	{
	}

	/* release memory; */
	fn_tmat_free<int>	 ( coeffmat );
	fn_tvec_free<int>	 ( coeffrslt );
}


/**
 * \brief	Calculate the error between 'bulk1_FGJP_density' and 'bulk1_FGJPre_density';
 *			also for 'bulk2_FGJP_density' and 'bulk2_FGJPre_density';
 */
void fn_proj_common_error ( stu_bulk_param		*sbulkparam,
							stu_system_param	*ssysparam )
{
	char fname0[FILELEN];
	char fname1[FILELEN];
	sprintf(fname0, "%s/bulk%d_FGJP_density.dat", rsltDir, sbulkparam->sflag);
	sprintf(fname1, "%s/bulk%d_FGJPre_density.dat", rsltDir, sbulkparam->sflag);

	int		xlen0	= sbulkparam->sbndv->xlen;
	int		xlen1	= sbulkparam->srebndv->xlen;
	int		y_num	= ssysparam->err_y_num;
	int		z_num	= ssysparam->err_z_num;

	int status = 1;
	if ( sbulkparam->dimPhy == 2 )
	{
		/* read data from 'fname0'; */
		tvec<double> x0		= fn_tvec_init<double> ( xlen0 );
		tvec<double> y0		= fn_tvec_init<double> ( y_num );
		tvec<double> data0	= fn_tvec_init<double> ( xlen0 * y_num );
		status = fn_read_data_2d ( fname0, x0, y0, data0 );

 		/* read data from 'fname1'; */
		tvec<double> x1		= fn_tvec_init<double> ( xlen1 );
		tvec<double> y1		= fn_tvec_init<double> ( y_num );
		tvec<double> data1	= fn_tvec_init<double> ( xlen1 * y_num );
		status = fn_read_data_2d ( fname1, x1, y1, data1 );

		/* calculate error; */
		fn_tvec_add<double> ( x0,	 x1,	1.0, -1.0 ); // x0 - x1;
		fn_tvec_add<double> ( y0,	 y1,	1.0, -1.0 ); // y0 - y1;
		fn_tvec_add<double> ( data0, data1, 1.0, -1.0 ); // data0 - data1;
		double err_x	= fn_tvec_maxAbs<double> ( x0 );
		double err_y	= fn_tvec_maxAbs<double> ( y0 );
		double err_data = fn_tvec_maxAbs<double> ( data0 );
		if ( myrank == 0 )
		{
			printf("\t Error between the 'phi' before and after changing space.\n");
			printf("\t ---> x error = %.6e\n", err_x);
			printf("\t ---> y error = %.6e\n", err_y);
			printf("\t ---> density error = %.6e\n", err_data);
		}

		/* release memory; */
		fn_tvec_free<double> ( x0	 );
		fn_tvec_free<double> ( x1	 );
		fn_tvec_free<double> ( y0	 );
		fn_tvec_free<double> ( y1	 );
		fn_tvec_free<double> ( data0 );
		fn_tvec_free<double> ( data1 );
	} 
	else if ( sbulkparam->dimPhy == 3 )
	{
		/* read data from 'fname0'; */
		tvec<double> x0		= fn_tvec_init<double> ( xlen0 );
		tvec<double> y0		= fn_tvec_init<double> ( y_num );
		tvec<double> z0		= fn_tvec_init<double> ( z_num );
		tvec<double> data0	= fn_tvec_init<double> ( xlen0 * y_num * z_num );
		fn_read_data_3d ( fname0, x0, y0, z0, data0 );

 		/* read data from 'fname1'; */
		tvec<double> x1		= fn_tvec_init<double> ( xlen1 );
		tvec<double> y1		= fn_tvec_init<double> ( y_num );
		tvec<double> z1		= fn_tvec_init<double> ( z_num );
		tvec<double> data1	= fn_tvec_init<double> ( xlen1 * y_num * z_num );
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
		if ( myrank == 0 )
		{
			printf("\t Error between the 'phi' before and after changing space.\n");
			printf("\t ---> x error = %.6e\n", err_x);
			printf("\t ---> y error = %.6e\n", err_y);
			printf("\t ---> z error = %.6e\n", err_z);
			printf("\t ---> density error = %.6e\n", err_data);
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
 * \brief	Read file for 2D data;
 */
int fn_read_data_2d		(		char		fname[FILELEN],		
							tvec<double>	x,
							tvec<double>	y,
							tvec<double>	data )
{
	FILE	*fp = fopen(fname, "r");
	int		status;
    if ( fp == NULL )	// invalid file;
	{
		if ( myrank == 0 )
			printf("Unable to read file '%s'. No such file or directory.\n", fname);
		return 1;
	}
    
	int nx, ny, len;
    status = fscanf(fp, "%d", &(nx));
	status = fscanf(fp, "%d", &(ny));
	status = fscanf(fp, "%d", &(len));

	for ( int i = 0; i < nx; i++ )
		status = fscanf(fp, "%lf", &(x.val[i]));
	for ( int i = 0; i < ny; i++ )
		status = fscanf(fp, "%lf", &(y.val[i]));
	for ( int i = 0; i < len; i++ )
		status = fscanf(fp, "%lf", &(data.val[i]));
	fclose(fp);
	return 0;
}


/**
 * \brief	Read file for 3D data;
 */
int fn_read_data_3d		(		char		fname[FILELEN],		
							tvec<double>	x,
							tvec<double>	y,
							tvec<double>	z,
							tvec<double>	data )
{
	FILE	*fp = fopen(fname, "r");
	int		status;
    if ( fp == NULL )	// invalid file;
	{
		if ( myrank == 0 )
			printf("Unable to read file '%s'. No such file or directory.\n", fname);
		return 1;
	}
    
	int nx, ny, nz, len;
    status = fscanf(fp, "%d", &(nx));
	status = fscanf(fp, "%d", &(ny));
	status = fscanf(fp, "%d", &(nz));
	status = fscanf(fp, "%d", &(len));

	for ( int i = 0; i < nx; i++ )
		status = fscanf(fp, "%lf", &(x.val[i]));
	for ( int i = 0; i < ny; i++ )
		status = fscanf(fp, "%lf", &(y.val[i]));
	for ( int i = 0; i < nz; i++ )
		status = fscanf(fp, "%lf", &(z.val[i]));
	for ( int i = 0; i < len; i++ )
		status = fscanf(fp, "%lf", &(data.val[i]));
	fclose(fp);
	return 0;
}


/**
 * \brief	compute the mass of left and right bulk phases after re-representation;
 */
void fn_obt_bulk_mass_common	(	stu_bulk_param			*sbulkparam,
									stu_system_param		*ssysparam,
									fftw_complex				mass	)
{
	/* parameters; */
	int		nd		= sbulkparam->srebndv->nd;
	int		xlen	= sbulkparam->srebndv->xlen;

	/* valid data; */
	int reComLen;
	if ( myrank == nprocs-1 )
		reComLen = ssysparam->cplxReDofs - displs_sys[myrank];
	else
		reComLen = alloc_local_sys;

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
			ind0 = j1*reComLen;	// the first element along Fourier direction;
			rhoCplx[0] += sbulkparam->srebndv->rhoJCplx.val[ind0][0] * ssysparam->sGJPv->d0JJ.val[i][j1];
			rhoCplx[1] += sbulkparam->srebndv->rhoJCplx.val[ind0][1] * ssysparam->sGJPv->d0JJ.val[i][j1];
		}
		ind1 = i*reComLen;		// the first element along Fourier direction;
		rhoCplx[0] += sbulkparam->srebndv->d0bnd.val[ind1][0];
		rhoCplx[1] += sbulkparam->srebndv->d0bnd.val[ind1][1];
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
