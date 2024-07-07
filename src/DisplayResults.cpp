/*! \file	DisplayResults.cpp
 *
 *  \brief	Display results;
 *			calculate the corresponding drawing data for displaying results;
 *
 */

#include "Data.h"
#include "Head.h"
#include "DataOperators.h"
#include "umfpack.h"
#include "Mytimer.h"
#include "functs.h"


/**
 * \brief	display densities of bulk phase;
 */
void fn_disp_bulk_density	(	tvec<fftw_complex>		src, 
								stu_bulk_param			*sbulkparam, 
									int					opt_iter,
									int					step )
{
	if ( myrank == 0 )
	{
		if ( sbulkparam->sflag == 1 )
			printf("\n===> Output plotted data (density of left bulk phase): ");
		else if ( sbulkparam->sflag == 2 )
			printf("\n===> Output plotted data (density of right bulk phase): ");
	}
	int	cplxDofs = sbulkparam->cplxDofs;

	mytimer_t timer;
	timer.reset();
	timer.start();

	/* calculate 'dirBox'; */
	matrix_inverse ( sbulkparam->rcpBox.val, sbulkparam->dirBox.val, sbulkparam->dimCpt	);
	for ( int i = 0; i < sbulkparam->dimCpt; i++ )
		for ( int j = 0; j < sbulkparam->dimCpt; j++ )
			sbulkparam->dirBox.val[i][j] *= 2*PI;

	FILE *densityfile;
	char densityName[FILELEN];
	sprintf(densityName, "%s/bulk%d_opt%d_density%d.dat", 
			rsltDir, sbulkparam->sflag, opt_iter, step);
	densityfile = fopen(densityName, "w");

	if ( sbulkparam->dimPhy == 2 )
	{
		double enlarge = sbulkparam->enlarge;
		/* the number of discrete points in real space; */
		int ny = (int) sbulkparam->plot_num, nz = ny; 
		double sizey = enlarge*sbulkparam->dirBox.val[0][0];
		double sizez = enlarge*sbulkparam->dirBox.val[1][1];
		double dy = sizey / (double)(ny-1);
		double dz = sizez / (double)(nz-1);

		/* save the number of discrete points; */
		fprintf(densityfile, "%d\t%d\t%d\n", ny, nz, ny*nz);

		// save discrete grids along each direction;
		for (int ky = 0; ky < ny; ky++) 
			fprintf(densityfile, "%+.15E\n", dy*ky - sizey/2.0);
		for (int kz = 0; kz < nz; kz++)
			fprintf(densityfile, "%+.15E\t", dz*kz - sizez/2.0);

		// save density; projection;
		for (int ky = 0; ky < ny; ky++)
		{
			double yVal = dy*ky - sizey/2.0;
			for (int kz = 0; kz < nz; kz++)
			{
				double zVal = dz*kz - sizez/2.0;
				double rho = 0.0;
				double tmpphase;

				// rhoCplx * exp(i (k1*y+k2*z));
				for (int i = 0; i < cplxDofs; i++)
				{
					double elm = fn_complex_abs ( src.val[i] );
					if ( elm > SPECTOL )
					{
						tmpphase = sbulkparam->sfftv->projPlane.val[i][0]*yVal + 
								sbulkparam->sfftv->projPlane.val[i][1]*zVal;
						rho += (src.val[i][0]*cos(tmpphase) - src.val[i][1]*sin(tmpphase));
					}
				}
				if ( myrank == 0 )
					fprintf(densityfile, "%+.15E\n", rho);
			}
		}
	}
	else if ( sbulkparam->dimPhy == 3 )
	{
		double enlarge = sbulkparam->enlarge;
		/* the number of discrete points in real space; */
		int nx = (int) sbulkparam->plot_num, ny = nx, nz = nx; 
		double sizex = enlarge*sbulkparam->dirBox.val[0][0];
		double sizey = enlarge*sbulkparam->dirBox.val[1][1];
		double sizez = enlarge*sbulkparam->dirBox.val[2][2];
		double dx = sizex / (double)(nx-1);
		double dy = sizey / (double)(ny-1);
		double dz = sizez / (double)(nz-1);

		/* save the number of discrete points; */
		fprintf(densityfile, "%d\t%d\t%d\t%d\n", nx, ny, nz, nx*ny*nz);

		/* save discrete grids along each direction; */
		for (int kx = 0; kx < nx; kx++)
			fprintf(densityfile, "%+.15E\n", dx*kx - sizex/2.0);
		for (int ky = 0; ky < ny; ky++)
			fprintf(densityfile, "%+.15E\n", dy*ky - sizey/2.0);
		for (int kz = 0; kz < nz; kz++)
			fprintf(densityfile, "%+.15E\n", dz*kz - sizez/2.0);

		/* save density; projection; */
		for (int kx = 0; kx < nx; kx++)
		{
			double xVal = dx*kx - sizex/2.0;
			for (int ky = 0; ky < ny; ky++)
			{
				double yVal = dy*ky - sizey/2.0;
				for (int kz = 0; kz < nz; kz++)
				{
					double zVal = dz*kz - sizez/2.0;
					double rho = 0.0;
					double tmpphase;

					/* rhoCplx * exp(i (k1*y+k2*z)); */
					for (int i = 0; i < cplxDofs; i++)
					{
						double elm = fn_complex_abs ( src.val[i] );
						if ( elm > SPECTOL )
						{
							tmpphase = sbulkparam->sfftv->projPlane.val[i][0]*xVal + 
								sbulkparam->sfftv->projPlane.val[i][1]*yVal + 
								sbulkparam->sfftv->projPlane.val[i][2]*zVal;
							rho += (src.val[i][0]*cos(tmpphase) - src.val[i][1]*sin(tmpphase));
						}
					}
					if ( myrank == 0 )
						fprintf(densityfile, "%+.15E\n", rho);
				}
			}
		}
	}
	timer.pause();
	if ( myrank == 0 )
		printf("Step %d, %f seconds\n\n", step, timer.get_current_time());
	fclose(densityfile);
}


bool myComp(mySortVec a, mySortVec b)
{
	return (a.Data[0]*a.Data[0]+a.Data[1]*a.Data[1] > b.Data[0]*b.Data[0]+b.Data[1]*b.Data[1]);
}


bool scaleComp(mySortWave a, mySortWave b)
{
	return (a.Data > b.Data);
}


/**
 * \brief	display Fourier coefficients of bulk phase;
 */
void fn_disp_Fourier_coeff	(	tvec<fftw_complex>		src, 
								stu_bulk_param			*sbulkparam, 
									int					opt_iter,
									int					step )
{
	int	cplxDofs = sbulkparam->cplxDofs;

	FILE *fFourCoeff;
	char fname[FILELEN];
	sprintf(fname, "%s/rank%d/bulk%d_opt%d_field%d.dat", 
			rsltDir, myrank, sbulkparam->sflag, opt_iter, step);
	/* clear existing contents; */
	fFourCoeff = fopen(fname, "w");

	/* save the Fourier coefficients whose values are greater than 1.0e-16; */
	for (int k = 0; k < cplxDofs; k++)
	{
		double tmp = 0.0;
		for ( int j = 0; j < 2; j++ )
		{
			double tmpVal = src.val[k][j];
			tmp += tmpVal * tmpVal;
		}
		tmp = sqrt(tmp);
		if (tmp > 1.0e-16)
		{
			for (int i = 0; i < sbulkparam->dimCpt; i++)
			{
				fprintf(fFourCoeff, "% d\t", sbulkparam->sfftv->indKspace.val[k][i]);
			}
			fprintf(fFourCoeff, "% e\t% e\n", src.val[k][0], src.val[k][1]);
		}
	}
	fclose(fFourCoeff);
}


/**
 * \brief	display plane waves of bulk phase;
 */
void fn_disp_bulk_plane_wave	(	tvec<fftw_complex>		src, 
									stu_bulk_param			*sbulkparam, 
										int					opt_iter,
										int					step)
{
	int	cplxDofs = sbulkparam->cplxDofs;

	FILE *fprojPlane;
	char projPlaneName[FILELEN];
	sprintf(projPlaneName, "%s/rank%d/bulk%d_opt%d_planeWave%d.dat", 
			rsltDir, myrank, sbulkparam->sflag, opt_iter, step);
	fprojPlane = fopen(projPlaneName, "w");

	mySortWave *myPlaneWave = (mySortWave*) malloc(sizeof(mySortWave) * cplxDofs);
	for (int k = 0; k < cplxDofs; k++)
	{
		myPlaneWave[k].Wave = (double *)malloc(sizeof(double) * sbulkparam->dimPhy);
		for (int i = 0; i < sbulkparam->dimPhy; i++)
			myPlaneWave[k].Wave[i] = sbulkparam->sfftv->projPlane.val[k][i];
		myPlaneWave[k].Data = sqrt(src.val[k][0]*src.val[k][0] + src.val[k][1]*src.val[k][1]);
	}

	std::sort(myPlaneWave, myPlaneWave+cplxDofs, scaleComp);

	// save the Fourier coefficients whose values are greater than 1.0e-16;
	for (int i = 0; i < cplxDofs; i ++)
	{
		if (myPlaneWave[i].Data > 1.0e-16)
		{
			for (int j = 0; j < sbulkparam->dimPhy; j++)
				fprintf(fprojPlane, "% f\t ", myPlaneWave[i].Wave[j]);
			fprintf(fprojPlane, "% e\n", myPlaneWave[i].Data);
		}
	}
	
	for (int i = 0; i < cplxDofs; i++)
	{
		free(myPlaneWave[i].Wave);
	}
	free(myPlaneWave);
	fclose(fprojPlane);
}



/**
 * \brief	display densities of bulk phase by a parallel code;
 */
void fn_disp_bulk_density_parallel	(	tvec<fftw_complex>		src, 
										stu_bulk_param			*sbulkparam, 
											int					opt_iter,
											int					step )
{
	if ( myrank == 0 )
	{
		if ( sbulkparam->sflag == 1 )
			printf("\n===> Output plotted data (density of left bulk phase): ");
		else if ( sbulkparam->sflag == 2 )
			printf("\n===> Output plotted data (density of right bulk phase): ");
	}
	
	mytimer_t timer;
	timer.reset();
	timer.start();

	/* calculate 'dirBox'; */
	matrix_inverse ( sbulkparam->rcpBox.val, sbulkparam->dirBox.val, sbulkparam->dimCpt	);
	for ( int i = 0; i < sbulkparam->dimCpt; i++ )
		for ( int j = 0; j < sbulkparam->dimCpt; j++ )
			sbulkparam->dirBox.val[i][j] *= 2*PI;

	FILE *densityfile;
	char densityName[FILELEN];
	sprintf(densityName, "%s/bulk%d_opt%d_density%d.dat", 
			rsltDir, sbulkparam->sflag, opt_iter, step);
	densityfile = fopen(densityName, "w");

	if ( sbulkparam->dimPhy == 2 )
	{
		double enlarge = sbulkparam->enlarge;
		/* the number of discrete points in real space; */
		int ny = (int) sbulkparam->plot_num, nz = ny; 
		double sizey = enlarge*sbulkparam->dirBox.val[0][0];
		double sizez = enlarge*sbulkparam->dirBox.val[1][1];
		double dy = sizey / (double)(ny-1);
		double dz = sizez / (double)(nz-1);

		/* save the number of discrete points; */
		fprintf(densityfile, "%d\t%d\t%d\n", ny, nz, ny*nz);

		// save discrete grids along each direction;
		for (int ky = 0; ky < ny; ky++) 
			fprintf(densityfile, "%+.15E\n", dy*ky - sizey/2.0);
		for (int kz = 0; kz < nz; kz++)
			fprintf(densityfile, "%+.15E\t", dz*kz - sizez/2.0);

		// save density; projection;
		for (int ky = 0; ky < ny; ky++)
		{
			double yVal = dy*ky - sizey/2.0;
			for (int kz = 0; kz < nz; kz++)
			{
				double zVal = dz*kz - sizez/2.0;
				double rho = 0.0;
				double tmpphase;

				// rhoCplx * exp(i (k1*y+k2*z));
				for (int i = 0; i < alloc_local; i++)
				{
					double elm = fn_complex_abs ( src.val[i] );
					if ( elm > SPECTOL )
					{
						tmpphase = sbulkparam->sfftv->projPlane.val[i][0]*yVal + 
								sbulkparam->sfftv->projPlane.val[i][1]*zVal;
						rho += (src.val[i][0]*cos(tmpphase) - src.val[i][1]*sin(tmpphase));
					}
				}
				double rhoRecv = rho;
				MPI_Reduce ( &rho, &rhoRecv, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
				if ( myrank == 0 )
					fprintf(densityfile, "%+.15E\n", rhoRecv);
			}
		}
	}
	else if ( sbulkparam->dimPhy == 3 )
	{
		double enlarge = sbulkparam->enlarge;
		/* the number of discrete points in real space; */
		int nx = (int) sbulkparam->plot_num, ny = nx, nz = nx; 
		double sizex = enlarge*sbulkparam->dirBox.val[0][0];
		double sizey = enlarge*sbulkparam->dirBox.val[1][1];
		double sizez = enlarge*sbulkparam->dirBox.val[2][2];
		double dx = sizex / (double)(nx-1);
		double dy = sizey / (double)(ny-1);
		double dz = sizez / (double)(nz-1);

		/* save the number of discrete points; */
		fprintf(densityfile, "%d\t%d\t%d\t%d\n", nx, ny, nz, nx*ny*nz);

		/* save discrete grids along each direction; */
		for (int kx = 0; kx < nx; kx++)
			fprintf(densityfile, "%+.15E\n", dx*kx - sizex/2.0);
		for (int ky = 0; ky < ny; ky++)
			fprintf(densityfile, "%+.15E\n", dy*ky - sizey/2.0);
		for (int kz = 0; kz < nz; kz++)
			fprintf(densityfile, "%+.15E\n", dz*kz - sizez/2.0);

		/* save density; projection; */
		for (int kx = 0; kx < nx; kx++)
		{
			double xVal = dx*kx - sizex/2.0;
			for (int ky = 0; ky < ny; ky++)
			{
				double yVal = dy*ky - sizey/2.0;
				for (int kz = 0; kz < nz; kz++)
				{
					double zVal = dz*kz - sizez/2.0;
					double rho = 0.0;
					double tmpphase;

					/* rhoCplx * exp(i (k1*y+k2*z)); */
					for (int i = 0; i < alloc_local; i++)
					{
						double elm = fn_complex_abs ( src.val[i] );
						if ( elm > SPECTOL )
						{
							tmpphase = sbulkparam->sfftv->projPlane.val[i][0]*xVal + 
								sbulkparam->sfftv->projPlane.val[i][1]*yVal + 
								sbulkparam->sfftv->projPlane.val[i][2]*zVal;
							rho += (src.val[i][0]*cos(tmpphase) - src.val[i][1]*sin(tmpphase));
						}
					}
					double rhoRecv = rho;
					MPI_Reduce ( &rho, &rhoRecv, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
					if ( myrank == 0 )
						fprintf(densityfile, "%+.15E\n", rhoRecv);
				}
			}
		}
	}
	MPI_Barrier ( MPI_COMM_WORLD );
	timer.pause();
	if ( myrank == 0 )
		printf("Step %d, %f seconds\n\n", step, timer.get_current_time());
	fclose(densityfile);
}


/**
 * \brief	display Fourier coefficients of bulk phase by a parallel code;
 */
void fn_disp_Fourier_coeff_parallel	(	tvec<fftw_complex>		src, 
										stu_bulk_param			*sbulkparam, 
											int					opt_iter,
											int					step )
{
	FILE *fFourCoeff;
	char fname[FILELEN];
	sprintf(fname, "%s/rank%d/bulk%d_opt%d_field%d.dat", 
			rsltDir, myrank, sbulkparam->sflag, opt_iter, step);
	/* clear existing contents; */
	fFourCoeff = fopen(fname, "w");

	/* save the Fourier coefficients whose values are greater than 1.0e-16; */
	for (int k = 0; k < alloc_local; k++)
	{
		double tmp = 0.0;
		for ( int j = 0; j < 2; j++ )
		{
			double tmpVal = src.val[k][j];
			tmp += tmpVal * tmpVal;
		}
		tmp = sqrt(tmp);
		if (tmp > 1.0e-16)
		{
			for (int i = 0; i < sbulkparam->dimCpt; i++)
			{
				fprintf(fFourCoeff, "% d\t", sbulkparam->sfftv->indKspace.val[k][i]);
			}
			fprintf(fFourCoeff, "% e\t% e\n", src.val[k][0], src.val[k][1]);
		}
	}
	fclose(fFourCoeff);
}


/**
 * \brief	display plane waves of bulk phase by a parallel code;
 */
void fn_disp_bulk_plane_wave_parallel	(	tvec<fftw_complex>		src, 
											stu_bulk_param			*sbulkparam, 
												int					opt_iter,
												int					step)
{
	FILE *fprojPlane;
	char projPlaneName[FILELEN];
	sprintf(projPlaneName, "%s/rank%d/bulk%d_opt%d_planeWave%d.dat", 
			rsltDir, myrank, sbulkparam->sflag, opt_iter, step);
	fprojPlane = fopen(projPlaneName, "w");

	mySortWave *myPlaneWave = (mySortWave*) malloc(sizeof(mySortWave) * alloc_local);
	for (int k = 0; k < alloc_local; k++)
	{
		myPlaneWave[k].Wave = (double *)malloc(sizeof(double) * sbulkparam->dimPhy);
		for (int i = 0; i < sbulkparam->dimPhy; i++)
			myPlaneWave[k].Wave[i] = sbulkparam->sfftv->projPlane.val[k][i];
		myPlaneWave[k].Data = sqrt(src.val[k][0]*src.val[k][0] + src.val[k][1]*src.val[k][1]);
	}

	std::sort(myPlaneWave, myPlaneWave+alloc_local, scaleComp);

	// save the Fourier coefficients whose values are greater than 1.0e-16;
	for (int i = 0; i < alloc_local; i ++)
	{
		if (myPlaneWave[i].Data > 1.0e-16)
		{
			for (int j = 0; j < sbulkparam->dimPhy; j++)
				fprintf(fprojPlane, "% f\t ", myPlaneWave[i].Wave[j]);
			fprintf(fprojPlane, "% e\n", myPlaneWave[i].Data);
		}
	}
	
	for (int i = 0; i < alloc_local; i++)
	{
		free(myPlaneWave[i].Wave);
	}
	free(myPlaneWave);
	fclose(fprojPlane);
}



/*
 * \brief	display density of bulk phases after projecting to the space of GJPs;
 *
 * \param	rhoJCplx:		coefficients in Fourier and GJP space at the same time;
 * \param	sbulkparam:		structure body for bulk phase;
 * \param	ssysparam:		structure body for the interface system;
 *							including some plotting parameters and sGJPv;
 */
void fn_disp_bulk_proj_density (	tvec<fftw_complex>	rhoJCplx, 
									stu_bulk_param		*sbulkparam, 
									stu_system_param	*ssysparam )
{
	if ( myrank == 0 )
	{
		if ( sbulkparam->sflag == 1 )
			printf("\n===> Output plotted data (density of left bulk phase after projecting): \n");
		else if ( sbulkparam->sflag == 2 )
			printf("\n===> Output plotted data (density of right bulk phase after projecting): \n");
	}
	
	/* parameters; */
	int		nd		= sbulkparam->sbndv->nd;
	int		xlen	= sbulkparam->sbndv->xlen;
	int	cplxDofs	= sbulkparam->sbndv->cplxDofs;
	double	y_start	= ssysparam->err_y_start;
	double	z_start	= ssysparam->err_z_start;
	double	y_range = ssysparam->err_y_range;
	double	z_range = ssysparam->err_z_range;
	int		y_num	= ssysparam->err_y_num;
	int		z_num	= ssysparam->err_z_num;

	double dy = y_range / y_num;		// step along y-direction;
	double dz = z_range / z_num;		// step along z-direction;

	mytimer_t timer;
	timer.reset();
	timer.start();

	/* memory allocation for temporary variables; */
	int ind0, ind1;
	tvec<fftw_complex> tmp = fn_tvec_init<fftw_complex> ( xlen * cplxDofs );
	tmp.row = xlen;
	tmp.col = cplxDofs;
	fn_tvec_setZero_complex ( tmp );

	/** 
	 * project to Fourier space;
	 *		rhoJCplx * d0JJ + d0bndlr;
	 */
	for ( int i = 0; i < xlen; i++ )
	{
		for ( int j0 = 0; j0 < cplxDofs; j0++ )
		{
			ind1 = i*cplxDofs + j0;
			for ( int j1 = 0; j1 < nd; j1++ )
			{
				ind0 = j1*cplxDofs + j0;
				tmp.val[ind1][0] += rhoJCplx.val[ind0][0] * ssysparam->sGJPv->d0JJ.val[i][j1];
				tmp.val[ind1][1] += rhoJCplx.val[ind0][1] * ssysparam->sGJPv->d0JJ.val[i][j1];
			}
			tmp.val[ind1][0] += sbulkparam->sbndv->d0bnd.val[ind1][0];
			tmp.val[ind1][1] += sbulkparam->sbndv->d0bnd.val[ind1][1];
		}
	}

	/* obtain global indexes of special positions (N/2); */
	tvec<int> gIndex = fn_serial_obt_half_gIndex ( sbulkparam->NCpt, cplxDofs );

	/* the start index of each process; */
	int gIndexRank = 0;
	for ( int j = 0; j < gIndex.len; j++ )
	{
		if ( gIndex.val[j] >= 0 )
		{
			gIndexRank = j;
			break;
		}
	}

	/* open file; */
	FILE *densityfile;
	char densityName[FILELEN];
	sprintf(densityName, "%s/bulk%d_FGJP_density.dat", rsltDir, sbulkparam->sflag);
	densityfile = fopen(densityName, "w");

	/**
	 * calculate density by 
	 *		sum_{k} sum_{j} rhoJCplx d0JJ exp(i(tilde{R}'PBk)'tilde{r});
	 */
	if ( myrank == 0 )
		printf("progress: ");
	int pmod = xlen / 10;
	if ( sbulkparam->dimPhy == 2 )
	{
		/* save the number of discrete points; */
		fprintf(densityfile, "%d\t%d\t%d\n", xlen, y_num, xlen*y_num);

		/* save discrete grids along each direction; */
		for ( int i = 0; i < xlen; i++ )
			fprintf(densityfile, "%+.15E\n", ssysparam->sGJPv->x.val[i]);
		for ( int i = 0; i < y_num; i++ ) 
			fprintf(densityfile, "%+.15E\n", y_start + dy*i);

		/**
		 * save density rho;
		 *		rho:	size is xlen*y_num;
		 */
		double tmpPhase, rho;
		for ( int i = 0; i < xlen; i++ )
		{
			for ( int j = 0; j < y_num; j++ )
			{
				rho = 0.0;
				int ind = gIndexRank;
				for ( int j0 = 0; j0 < cplxDofs; j0++ )
				{
					if ( sbulkparam->skip_half == 1 && j0 == gIndex.val[ind] ) // skip the special index;
						ind ++ ;
					else
					{
						ind1 = i*cplxDofs + j0;
						/* save elements whose modulus are greater than SPECTOL; */
						double elm = fn_complex_abs ( tmp.val[ind1] );
						if ( elm > SPECTOL )
						{
							tmpPhase = sbulkparam->sfftv->projPlane.val[j0][1]*(y_start+dy*j);
							rho += tmp.val[ind1][0]*cos(tmpPhase) - tmp.val[ind1][1]*sin(tmpPhase);
						}
					}
				}
				if ( myrank == 0 )
					fprintf(densityfile, "%+.15E\n", rho);
			}
			if ( myrank == 0 )
				if ( i % pmod == 0 )
					printf("%.2lf ", (double) i/xlen);
		}
	}
	else if ( sbulkparam->dimPhy == 3 )
	{
		/* save the number of discrete points; */
		fprintf(densityfile, "%d\t%d\t%d\t%d\n", xlen, y_num, z_num, xlen*y_num*z_num);

		/* save discrete grids along each direction; */
		for ( int i = 0; i < xlen; i++ )
			fprintf(densityfile, "%+.15E\n", ssysparam->sGJPv->x.val[i]);
		for ( int i = 0; i < y_num; i++ ) 
			fprintf(densityfile, "%+.15E\n", y_start + dy*i);
		for ( int i = 0; i < z_num; i++ ) 
			fprintf(densityfile, "%+.15E\n", z_start + dz*i);

		/**
		 * save density rho;
		 *		rho:	size is xlen*y_num;
		 */
		double tmpPhase, rho;
		for ( int i = 0; i < xlen; i++ )
		{
			for ( int j = 0; j < y_num; j++ )
			{
				for ( int k = 0; k < z_num; k++ )
				{
					rho = 0.0;
					int ind = gIndexRank;
					for ( int j0 = 0; j0 < cplxDofs; j0++ )
					{
						if ( sbulkparam->skip_half == 1 && j0 == gIndex.val[ind] ) // skip the special index;
							ind ++ ;
						else
						{
							ind1 = i*cplxDofs + j0;
							/* save elements whose modulus are greater than SPECTOL; */
							double elm = fn_complex_abs ( tmp.val[ind1] );
							if ( elm > SPECTOL )
							{
								tmpPhase = sbulkparam->sfftv->projPlane.val[j0][1]*(y_start+dy*j) +
										   sbulkparam->sfftv->projPlane.val[j0][2]*(z_start+dz*k);
								rho += tmp.val[ind1][0]*cos(tmpPhase) - tmp.val[ind1][1]*sin(tmpPhase);
							}
						}
					}
					if ( myrank == 0 )
						fprintf(densityfile, "%+.15E\n", rho);
				}
			}
			if ( myrank == 0 )
				if ( i % pmod == 0 )
					printf("%.2lf ", (double) i/xlen);
		}
	}
	if ( myrank == 0 )
		printf("finish\n");
	timer.pause();
	if ( myrank == 0 )
	{
		if ( sbulkparam->sflag == 1 )
			printf("left bulk density after projecting: ");
		else if ( sbulkparam->sflag == 2 )
			printf("right bulk density after projecting: ");
		printf("%f seconds\n\n", timer.get_current_time());
	}
	fclose(densityfile);

	/* release memory; */
	fn_tvec_free<fftw_complex> ( tmp );
	fn_tvec_free<int> ( gIndex );
}


/*
 * \brief	display density of bulk phases in the common space;
 *
 * \param	rhoReJCplx:		coefficients in the common space;
 * \param	sbulkparam:		structure body for bulk phase;
 * \param	ssysparam:		structure body for the interface system;
 *							including some plotting parameters and sGJPv;
 */
void fn_disp_bulk_common_density (	tvec<fftw_complex>	rhoReJCplx, 
									stu_bulk_param		*sbulkparam, 
									stu_system_param	*ssysparam )
{
	if ( myrank == 0 )
	{
		if ( sbulkparam->sflag == 1 )
			printf("\n===> Output plotted data (density of left bulk phase in the common space): \n");
		else if ( sbulkparam->sflag == 2 )
			printf("\n===> Output plotted data (density of right bulk phase in the common space): \n");
	}
	
	/* parameters; */
	int		nd		= sbulkparam->srebndv->nd;
	int		xlen	= sbulkparam->srebndv->xlen;
	int cplxReDofs  = sbulkparam->srebndv->cplxDofs;
	double	y_start	= ssysparam->err_y_start;
	double	z_start	= ssysparam->err_z_start;
	double	y_range = ssysparam->err_y_range;
	double	z_range = ssysparam->err_z_range;
	int		y_num	= ssysparam->err_y_num;
	int		z_num	= ssysparam->err_z_num;

	double dy = y_range / y_num;		// step along y-direction;
	double dz = z_range / z_num;		// step along z-direction;

	mytimer_t timer;
	timer.reset();
	timer.start();

	/* memory allocation for temporary variables; */
	int ind0, ind1;
	tvec<fftw_complex> tmp = fn_tvec_init<fftw_complex> ( xlen * alloc_local_sys );
	tmp.row = xlen;
	tmp.col = alloc_local_sys;
	fn_tvec_setZero_complex ( tmp );

	/** 
	 * project to Fourier space;
	 *		rhoReJCplx * d0JJ + d0rebndlr;
	 */
	for ( int i = 0; i < xlen; i++ )
	{
		for ( int j0 = 0; j0 < alloc_local_sys; j0++ )
		{
			ind1 = i*alloc_local_sys + j0;
			for ( int j1 = 0; j1 < nd; j1++ )
			{
				ind0 = j1*alloc_local_sys + j0;
				tmp.val[ind1][0] += rhoReJCplx.val[ind0][0] * ssysparam->sGJPv->d0JJ.val[i][j1];
				tmp.val[ind1][1] += rhoReJCplx.val[ind0][1] * ssysparam->sGJPv->d0JJ.val[i][j1];
			}
			tmp.val[ind1][0] += sbulkparam->srebndv->d0bnd.val[ind1][0];
			tmp.val[ind1][1] += sbulkparam->srebndv->d0bnd.val[ind1][1];
		}
	}

	/* obtain global indexes of special positions (N/2); */
	tvec<int> gIndex = fn_obt_half_gIndex ( ssysparam->NCpt, cplxReDofs );
	for ( int j = 0; j < gIndex.len; j++ ) gIndex.val[j] -= displs_sys[myrank];

	/* the start index of each process; */
	int gIndexRank = 0;
	for ( int j = 0; j < gIndex.len; j++ )
	{
		if ( gIndex.val[j] >= 0 )
		{
			gIndexRank = j;
			break;
		}
	}

	/* open file; */
	FILE *densityfile;
	char densityName[FILELEN];
	sprintf(densityName, "%s/bulk%d_FGJPre_density.dat", rsltDir, sbulkparam->sflag);
	densityfile = fopen(densityName, "w");

	/**
	 * calculate density by 
	 *		sum_{k} sum_{j} rhoReJCplx d0JJ exp(i(common rotateProjBoxMat)'tilde{r});
	 */
	if ( myrank == 0 )
		printf("progress: ");
	int pmod = xlen / 10;
	if ( sbulkparam->dimPhy == 2 )
	{
		/* save the number of discrete points; */
		fprintf(densityfile, "%d\t%d\t%d\n", xlen, y_num, xlen*y_num);

		/* save discrete grids along each direction; */
		for ( int i = 0; i < xlen; i++ )
			fprintf(densityfile, "%+.15E\n", ssysparam->sGJPv->x.val[i]);
		for ( int i = 0; i < y_num; i++ ) 
			fprintf(densityfile, "%+.15E\n", y_start + dy*i);

		/**
		 * save density rho;
		 *		rho:	size is xlen*y_num;
		 */
		double tmpPhase, rho;
		for ( int i = 0; i < xlen; i++ )
		{
			for ( int j = 0; j < y_num; j++ )
			{
				rho = 0.0;
				int ind = gIndexRank;
				for ( int j0 = 0; j0 < alloc_local_sys; j0++ )
				{
					if ( ssysparam->skip_half == 1 && j0 == gIndex.val[ind] ) // skip the special index;
						ind ++ ;
					else
					{
						ind1 = i*alloc_local_sys + j0;
						/* save elements whose modulus are greater than SPECTOL; */
						double elm = fn_complex_abs ( tmp.val[ind1] );
						if ( elm > SPECTOL )
						{
							tmpPhase = ssysparam->sfftv->projPlane.val[j0][0]*(y_start+dy*j);
							rho += tmp.val[ind1][0]*cos(tmpPhase) - tmp.val[ind1][1]*sin(tmpPhase);
						}
					}
				}
				double rhoRecv = rho;
				MPI_Reduce ( &rho, &rhoRecv, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
				if ( myrank == 0 )
					fprintf(densityfile, "%+.15E\n", rhoRecv);
			}
			if ( myrank == 0 )
				if ( i % pmod == 0 )
					printf("%.2lf ", (double) i/xlen);
		}
	}
	else if ( sbulkparam->dimPhy == 3 )
	{
		/* save the number of discrete points; */
		fprintf(densityfile, "%d\t%d\t%d\t%d\n", xlen, y_num, z_num, xlen*y_num*z_num);

		/* save discrete grids along each direction; */
		for ( int i = 0; i < xlen; i++ )
			fprintf(densityfile, "%+.15E\n", ssysparam->sGJPv->x.val[i]);
		for ( int i = 0; i < y_num; i++ ) 
			fprintf(densityfile, "%+.15E\n", y_start + dy*i);
		for ( int i = 0; i < z_num; i++ ) 
			fprintf(densityfile, "%+.15E\n", z_start + dz*i);
	
		/**
		 * save density rho;
		 *		rho:	size is xlen*y_num;
		 */
		double tmpPhase, rho;
		for ( int i = 0; i < xlen; i++ )
		{
			for ( int j = 0; j < y_num; j++ )
			{
				for ( int k = 0; k < z_num; k++ )
				{
					rho = 0.0;
					int ind = gIndexRank;
					for ( int j0 = 0; j0 < alloc_local_sys; j0++ )
					{
						if ( ssysparam->skip_half == 1 && j0 == gIndex.val[ind] ) // skip the special index;
							ind ++ ;
						else
						{
							ind1 = i*alloc_local_sys + j0;
							/* save elements whose modulus are greater than SPECTOL; */
							double elm = fn_complex_abs ( tmp.val[ind1] );
							if ( elm > SPECTOL )
							{
								tmpPhase = ssysparam->sfftv->projPlane.val[j0][0]*(y_start+dy*j) +
										   ssysparam->sfftv->projPlane.val[j0][1]*(z_start+dz*k);
								rho += tmp.val[ind1][0]*cos(tmpPhase) - tmp.val[ind1][1]*sin(tmpPhase);
							}
						}
					}
					double rhoRecv = rho;
					MPI_Reduce ( &rho, &rhoRecv, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
					if ( myrank == 0 )
						fprintf(densityfile, "%+.15E\n", rhoRecv);
				}
			}
			if ( myrank == 0 )
				if ( i % pmod == 0 )
					printf("%.2lf ", (double) i/xlen);
		}
	}
	if ( myrank == 0 )
		printf("finish\n");
	MPI_Barrier ( MPI_COMM_WORLD );
	timer.pause();
	if ( myrank == 0 )
	{
		if ( sbulkparam->sflag == 1 )
			printf("left bulk density in the common space: ");
		else if ( sbulkparam->sflag == 2 )
			printf("right bulk density in the common space: ");
		printf("%f seconds\n\n", timer.get_current_time());
	}
	fclose(densityfile);

	/* release memory; */
	fn_tvec_free<fftw_complex> ( tmp );
	fn_tvec_free<int> ( gIndex );
}


/*
 * \brief	display density of the interfce system;
 *
 * \param	rhoJCplx:		coefficients of the interface system;
 * \param	ssysparam:		structure body for the interface system;
 *							including some plotting parameters and sGJPv;
 */
void fn_disp_system_density		(	tvec<fftw_complex>		rhoJCplx, 
									stu_system_param		*ssysparam,
										int					step )
{
	if ( myrank == 0 )
		printf("\n===> Output plotted data (density of interface system): \n");

	/* parameters; */
	int		nd		= ssysparam->scbndv->nd;
	int		xlen	= ssysparam->scbndv->xlen;
	int cplxReDofs  = ssysparam->scbndv->cplxDofs;
	double	y_start	= ssysparam->y_start;
	double	z_start	= ssysparam->z_start;
	double	y_range = ssysparam->y_range;
	double	z_range = ssysparam->z_range;
	int		y_num	= ssysparam->y_num;
	int		z_num	= ssysparam->z_num;

	double dy = y_range / y_num;		// step along y-direction;
	double dz = z_range / z_num;		// step along z-direction;

	mytimer_t timer;
	timer.reset();
	timer.start();

	/* memory allocation for temporary variables; */
	int ind0, ind1;
	tvec<fftw_complex> tmp = fn_tvec_init<fftw_complex> ( xlen * alloc_local_sys );
	tmp.row = xlen;
	tmp.col = alloc_local_sys;
	fn_tvec_setZero_complex ( tmp );

	/** 
	 * project to Fourier space;
	 *		rhoJCplx * d0JJ + d0rebndlr;
	 */
	for ( int i = 0; i < xlen; i++ )
	{
		for ( int j0 = 0; j0 < alloc_local_sys; j0++ )
		{
			ind1 = i*alloc_local_sys + j0;
			for ( int j1 = 0; j1 < nd; j1++ )
			{
				ind0 = j1*alloc_local_sys + j0;
				tmp.val[ind1][0] += rhoJCplx.val[ind0][0] * ssysparam->sGJPv->d0JJ.val[i][j1];
				tmp.val[ind1][1] += rhoJCplx.val[ind0][1] * ssysparam->sGJPv->d0JJ.val[i][j1];
			}
			tmp.val[ind1][0] += ssysparam->scbndv->d0bnd.val[ind1][0];
			tmp.val[ind1][1] += ssysparam->scbndv->d0bnd.val[ind1][1];
		}
	}

	/* obtain global indexes of special positions (N/2); */
	tvec<int> gIndex = fn_obt_half_gIndex ( ssysparam->NCpt, cplxReDofs );
	for ( int j = 0; j < gIndex.len; j++ ) gIndex.val[j] -= displs_sys[myrank];

	/* the start index of each process; */
	int gIndexRank = 0;
	for ( int j = 0; j < gIndex.len; j++ )
	{
		if ( gIndex.val[j] >= 0 )
		{
			gIndexRank = j;
			break;
		}
	}

	/* open file; */
	FILE *densityfile;
	char densityName[FILELEN];
	sprintf(densityName, "%s/sys_density%d.dat", rsltDir, step);
	densityfile = fopen(densityName, "w");

	/**
	 * calculate density by 
	 *		sum_{k} sum_{j} rhoJCplx d0JJ exp(i(common rotateProjBoxMat)'tilde{r});
	 */
	if ( myrank == 0 )
		printf("progress: ");
	int pmod = xlen / 10;
	if ( ssysparam->dimRePhy == 1 )
	{
		/* save the number of discrete points; */
		fprintf(densityfile, "%d\t%d\t%d\n", xlen, y_num, xlen*y_num);

		/* save discrete grids along each direction; */
		for ( int i = 0; i < xlen; i++ )
			fprintf(densityfile, "%+.15E\n", ssysparam->sGJPv->x.val[i]);
		for ( int i = 0; i < y_num; i++ ) 
			fprintf(densityfile, "%+.15E\n", y_start + dy*i);

		/**
		 * save density rho;
		 *		rho:	size is xlen*y_num;
		 */
		double tmpPhase, rho;
		for ( int i = 0; i < xlen; i++ )
		{
			for ( int j = 0; j < y_num; j++ )
			{
				rho = 0.0;
				int ind = gIndexRank;
				for ( int j0 = 0; j0 < alloc_local_sys; j0++ )
				{
					if ( ssysparam->skip_half == 1 && j0 == gIndex.val[ind] ) // skip the special index;
						ind ++ ;
					else
					{
						ind1 = i*alloc_local_sys + j0;
						/* save elements whose modulus are greater than SPECTOL; */
						double elm = fn_complex_abs ( tmp.val[ind1] );
						if ( elm > SPECTOL )
						{
							tmpPhase = ssysparam->sfftv->projPlane.val[j0][0]*(y_start+dy*j);
							rho += tmp.val[ind1][0]*cos(tmpPhase) - tmp.val[ind1][1]*sin(tmpPhase);
						}
					}
				}
				double rhoRecv = rho;
				MPI_Reduce ( &rho, &rhoRecv, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
				if ( myrank == 0 )
					fprintf(densityfile, "%+.15E\n", rhoRecv);
			}
			if ( myrank == 0 )
				if ( i % pmod == 0 )
					printf("%.2lf ", (double) i/xlen);
		}
	}
	else if ( ssysparam->dimRePhy == 2 )
	{
		/* save the number of discrete points; */
		fprintf(densityfile, "%d\t%d\t%d\t%d\n", xlen, y_num, z_num, xlen*y_num*z_num);

		/* save discrete grids along each direction; */
		for ( int i = 0; i < xlen; i++ )
			fprintf(densityfile, "%+.15E\n", ssysparam->sGJPv->x.val[i]);
		for ( int i = 0; i < y_num; i++ ) 
			fprintf(densityfile, "%+.15E\n", y_start + dy*i);
		for ( int i = 0; i < z_num; i++ ) 
			fprintf(densityfile, "%+.15E\n", z_start + dz*i);

		/**
		 * save density rho;
		 *		rho:	size is xlen*y_num;
		 */
		double tmpPhase, rho;
		for ( int i = 0; i < xlen; i++ )
		{
			for ( int j = 0; j < y_num; j++ )
			{
				for ( int k = 0; k < z_num; k++ )
				{
					rho = 0.0;
					int ind = gIndexRank;
					for ( int j0 = 0; j0 < alloc_local_sys; j0++ )
					{
						if ( ssysparam->skip_half == 1 && j0 == gIndex.val[ind] ) // skip the special index;
							ind ++ ;
						else
						{
							ind1 = i*alloc_local_sys + j0;
							/* save elements whose modulus are greater than SPECTOL; */
							double elm = fn_complex_abs ( tmp.val[ind1] );
							if ( elm > SPECTOL )
							{
								tmpPhase = ssysparam->sfftv->projPlane.val[j0][0]*(y_start+dy*j) +
										   ssysparam->sfftv->projPlane.val[j0][1]*(z_start+dz*k);
								rho += tmp.val[ind1][0]*cos(tmpPhase) - tmp.val[ind1][1]*sin(tmpPhase);
							}
						}
					}
					double rhoRecv = rho;
					MPI_Reduce ( &rho, &rhoRecv, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
					if ( myrank == 0 )
						fprintf(densityfile, "%+.15E\n", rhoRecv);
				}
			}
			if ( myrank == 0 )
				if ( i % pmod == 0 )
					printf("%.2lf ", (double) i/xlen);
		}
	}
	if ( myrank == 0 )
		printf("finish\n");
	MPI_Barrier ( MPI_COMM_WORLD );
	timer.pause();
	if ( myrank == 0 )
		printf("interface density: %f seconds\n\n", timer.get_current_time());
	fclose(densityfile);

	/* release memory; */
	fn_tvec_free<fftw_complex> ( tmp );
	fn_tvec_free<int> ( gIndex );
}
