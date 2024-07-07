/*! \file	ProjFourGJP.cpp
 *
 * \brief	project rhoCplx to the space of general Jacobi polynomials (GJPs);
 *			introduce polynomials to make boundaries homogeneous;
 *
 */

#include "Data.h"
#include "Head.h"
#include "DataOperators.h"
#include "Mytimer.h"
#include "functs.h"


/**
 * \brief	Initialization and preparation;
 *			introduce polynomials to make boundaries homogeous;
 *			project to GJPs;
 *
 * \param	sbulkparam:		the structure body for bulk phases;
 * \param	ssysparam:		the structure body for the interface system;
 *
 */
void fn_project_Fourier_GJP (	stu_bulk_param		*sbulkparam,
								stu_system_param	*ssysparam )
{
	if ( myrank == 0 )
	{
		if ( sbulkparam->sflag == 1 )
			printf(" <========== Project to GJP space with left bulk phase ==========> \n\n");
		else if ( sbulkparam->sflag == 2 )
			printf(" <========== Project to GJP space with right bulk phase ==========> \n\n");
	}
	mytimer_t timer;
	timer.reset();
	timer.start();

	/* memory allocation; */
	fn_sbndv_memory_allocation(sbulkparam, ssysparam->sGJPv);
//	sbulkparam->sbndv->isSave = true;		// true: save; false: not save;
//	sbulkparam->sbndv->isTest = true;		// true: save; false: not save;

	/* obtain polynomials to deal with inhomogeneous boundary conditions; */
	if ( myrank == 0 )
		printf("deal with the inhomogeneous boundary conditions of order = ");
	fn_obt_poly_homo_bnd (sbulkparam, ssysparam->sGJPv, sbulkparam->sbndv->d0bnd, 0);
	fn_obt_poly_homo_bnd (sbulkparam, ssysparam->sGJPv, sbulkparam->sbndv->d1bnd, 1);
	fn_obt_poly_homo_bnd (sbulkparam, ssysparam->sGJPv, sbulkparam->sbndv->d2bnd, 2);
	if ( strcmp(model_type, "LP") == 0 )
	{
		fn_obt_poly_homo_bnd (sbulkparam, ssysparam->sGJPv, sbulkparam->sbndv->d3bnd, 3);
		fn_obt_poly_homo_bnd (sbulkparam, ssysparam->sGJPv, sbulkparam->sbndv->d4bnd, 4);
		fn_obt_poly_homo_bnd (sbulkparam, ssysparam->sGJPv, sbulkparam->sbndv->d5bnd, 5);
		fn_obt_poly_homo_bnd (sbulkparam, ssysparam->sGJPv, sbulkparam->sbndv->d6bnd, 6);
	}
	if ( myrank == 0 ) printf("\n");

	/* projection original space to GJP; */
	fn_project_GJP ( sbulkparam, ssysparam->sGJPv );

	/* save density calculating by rhoJCplx; */
	if ( myrank == 0 && sbulkparam->sbndv->isSave )
		fn_disp_bulk_proj_density ( sbulkparam->sbndv->rhoJCplx, sbulkparam, ssysparam );

	/* compute the mass; */
	fftw_complex mass;
	fn_obt_bulk_mass_project ( sbulkparam, ssysparam, mass );
	if ( myrank == 0 )
	{
		if ( sbulkparam->sflag == 1 )
			printf("\n\t ---> Mass of projecting left  bulk phase: %.10e\n", fn_complex_abs(mass));
		else if ( sbulkparam->sflag == 2 )
			printf("\n\t ---> Mass of projecting right bulk phase: %.10e\n", fn_complex_abs(mass));
	}

	timer.pause();
	if ( myrank == 0 )
		printf("\n\t ***** time cost of projection to GJPs: %f seconds *****\n\n", 
				timer.get_current_time());
}


/**
 * \brief	Memory allocation for boundary;
 */
void fn_sbndv_memory_allocation	(	stu_bulk_param		*sbulkparam, 
									stu_GJP_var			*sGJPv )
{	
	sbulkparam->sbndv->nd			= sGJPv->nd;
	sbulkparam->sbndv->xlen			= sGJPv->xlen;
	sbulkparam->sbndv->cplxDofs		= sbulkparam->cplxDofs;
	sbulkparam->sbndv->polyDegree	= sGJPv->polyDegree;
	sbulkparam->sbndv->x_range		= sGJPv->x_range;

	/* memory allocation; */
	int len	= sGJPv->xlen * sbulkparam->cplxDofs;
	sbulkparam->sbndv->d0bnd	= fn_tvec_init<fftw_complex> ( len );
	sbulkparam->sbndv->d1bnd	= fn_tvec_init<fftw_complex> ( len );
	sbulkparam->sbndv->d2bnd	= fn_tvec_init<fftw_complex> ( len );
	if ( strcmp(model_type, "LP") ==  0 )
	{
		sbulkparam->sbndv->d3bnd	= fn_tvec_init<fftw_complex> ( len );
		sbulkparam->sbndv->d4bnd	= fn_tvec_init<fftw_complex> ( len );
		sbulkparam->sbndv->d5bnd	= fn_tvec_init<fftw_complex> ( len );
		sbulkparam->sbndv->d6bnd	= fn_tvec_init<fftw_complex> ( len );
	}
}


/**
 * \brief	deal with the inhomogeneous boundary conditions;
 *			introduce polynomial to make boundary conditions homogeneous;
 *
 * \param	sbulkparam:		structure body for bulk phases;
 * \param	sGJPv:			structure body for GJP;
 * \param	bndlr;			cvector (also structure body variable);
 * \param	order;			the order of derivative;
 */
void fn_obt_poly_homo_bnd	(	stu_bulk_param			*sbulkparam,
								stu_GJP_var				*sGJPv, 
								tvec<fftw_complex>		bndlr,
									int					order )
{
	int nd		 = sbulkparam->sbndv->nd; 
	int xlen	 = sbulkparam->sbndv->xlen;
	int	cplxDofs = sbulkparam->sbndv->cplxDofs;
	double L	 = sbulkparam->sbndv->x_range;

	/* memory allocation for calculating the results of interpolating; */
	fftw_complex RPBkCplx, RPBkxCplx, rhoRPBkxCplx0, rhoRPBkxCplx1;
	tvec<fftw_complex> BC	= fn_tvec_init<fftw_complex> ( sGJPv->polyDegree ); 
	tvec<fftw_complex> pBC	= fn_tvec_init<fftw_complex> ( xlen );

	/* make boundaries homogeneous; */
	int ind = 0;
	for (int i = 0; i < cplxDofs; i++)
	{
		/* exp(iR_x'PBk); */
		RPBkCplx[0] = 0.0;
		RPBkCplx[1] = sbulkparam->sfftv->projPlane.val[i][0];
		for (int j = 0; j < sGJPv->polyDegree; j++)
		{
			BC.val[j][0] = 0.0;
			BC.val[j][1] = 0.0;
		}
		for (int s = 1; s <= 2; s++)			// s=1 for left; s=2 for right;
		{
			/* exp(i(R_x'PBk)'x); */
			RPBkxCplx[0] = cos( sbulkparam->sfftv->projPlane.val[i][0] * pow(-1,s) * L );
			RPBkxCplx[1] = sin( sbulkparam->sfftv->projPlane.val[i][0] * pow(-1,s) * L );
			/* rhoCplx * exp(i(R_x'PBk)'x); */
			fn_complex_multiply ( sbulkparam->sfftv->rhoCplx.val[i], RPBkxCplx, rhoRPBkxCplx0 );
			/**
			 * BC: the boundary values in the following order;
			 * [p(-L), p(L), p'(-L), p'(L), p''(-L), p''(L), p'''(-L), p'''(L)] for LP model;
			 * [p(-L), p(L), p'(-L), p'(L)] for LB model;
			 * (iR_x'PBk)^j * rhoCplx * exp(i(R_x'PBk)'x); j = 0,1,2,...
			 */
			for (int j = 0; j < sGJPv->polyDegree/2; j++)
			{
				BC.val[2*j+s-1][0] = rhoRPBkxCplx0[0];
				BC.val[2*j+s-1][1] = rhoRPBkxCplx0[1];
				fn_complex_multiply ( RPBkCplx, rhoRPBkxCplx0, rhoRPBkxCplx1 );
				rhoRPBkxCplx0[0] = rhoRPBkxCplx1[0];
				rhoRPBkxCplx0[1] = rhoRPBkxCplx1[1];
			}
		}
		fn_interpl_complex_bnd (sGJPv->x, pBC, BC, L, order);

		for ( int j = 0; j < xlen; j++ )
		{
			ind = j * cplxDofs + i;
			bndlr.val[ind][0] = pBC.val[j][0];
			bndlr.val[ind][1] = pBC.val[j][1];
		}
	}
	if ( myrank == 0 ) printf("%d ", order);

	/* release memory; */
	fn_tvec_free<fftw_complex> ( BC );
	fn_tvec_free<fftw_complex> ( pBC );

	/* save data; */
	if ( sbulkparam->sbndv->isTest )
	{
		char fwFile[FILELEN];
		sprintf(fwFile, "%s/bulk%d_d%dbnd.dat", rsltDir, sbulkparam->sflag, order);
		bndlr.row = xlen;
		bndlr.col = cplxDofs;
		fn_tvec_save_complex ( bndlr, fwFile );
	}
}
 

/**
 * \brief	construct a polynomial satisfying boundary conditions up to 3rd order derivative;
 *			complex type;
 *
 * \param	x:		the LGL nodes;
 * \param	fun:	the order-th order derivative of the constructed polynomial
 * \param	BC:		the boundary values in the following order;
 *				[p(-L), p(L), p'(-L), p'(L), p''(-L), p''(L), p'''(-L), p'''(L)] for LP model;
 *				[p(-L), p(L), p'(-L), p'(L)] for LB model;
 * \param	L:		the distance between two anchoring planes; x_range;
 * \param	order:	specifies the derivative order of returned function;
 *
 */
void fn_interpl_complex_bnd		(	tvec<double>		x,
									tvec<fftw_complex>	fun,
									tvec<fftw_complex>	BC, 
										double			L,
										int				order )
{
	/* memory allocation; */
	tvec<double> funReal = fn_tvec_init<double> ( fun.len );
	tvec<double> funCplx = fn_tvec_init<double> ( fun.len );
	tvec<double> BCReal	 = fn_tvec_init<double> ( BC.len  );
	tvec<double> BCCplx	 = fn_tvec_init<double> ( BC.len  );
	for ( int i = 0; i < BC.len; i++ )
	{
		BCReal.val[i] = BC.val[i][0];
		BCCplx.val[i] = BC.val[i][1];
	}
	
	/* construct a complex polynomial from real and complex space respectively; */
	fn_interpl_bnd ( x, funReal, BCReal, L, order );
	fn_interpl_bnd ( x, funCplx, BCCplx, L, order );
	for ( int i = 0; i < fun.len; i++ )
	{
		fun.val[i][0] = funReal.val[i];
		fun.val[i][1] = funCplx.val[i];
	}

	/* release memory; */
	fn_tvec_free<double> ( funReal );
	fn_tvec_free<double> ( funCplx );
	fn_tvec_free<double> ( BCReal  );
	fn_tvec_free<double> ( BCCplx  );
}


/**
 * \brief	construct a polynomial satisfying boundary conditions up to 3rd order derivative;
 *			double type;
 *
 * \param	x:		the LGL nodes;
 * \param	fun:	the order-th order derivative of the constructed polynomial
 * \param	BC:		the boundary values in the following order;
 *				[p(-L), p(L), p'(-L), p'(L), p''(-L), p''(L), p'''(-L), p'''(L)] for LP model;
 *				[p(-L), p(L), p'(-L), p'(L)] for LB model;
 * \param	L:		the distance between two anchoring planes; x_range;
 * \param	order:	specifies the derivative order of returned function;
 *
 */
void fn_interpl_bnd		(	tvec<double>	x,
							tvec<double>	fun,
							tvec<double>	BC, 
								double		L,
								int			order )
{
	int xlen = x.len;

	// compute the coefficients of polynomials and return the derivative function;
	if ( strcmp(model_type, "LB") == 0 )
	{
		double f0, f1, f2, f3;
		f0 = BC.val[0];	f1 = BC.val[1];	f2 = BC.val[2];	f3 = BC.val[3];
		double C0 = (f0-f1)/(4*pow(L,3)) + (f2+f3)/(4*pow(L,2));
		double C1 = (-f2+f3)/(4*L);
		double C2 = (-3*f0+3*f1)/(4*L) + (-f2-f3)/4;
		double C3 = (f0+f1)/2 + (f2-f3)*L/4;
		if ( order == 0 )
			for (int i = 0; i < xlen; i++)
				fun.val[i] = C0*pow(x.val[i],3) + C1*pow(x.val[i],2) + C2*x.val[i] + C3;
		else if ( order == 1 )
			for (int i = 0; i < xlen; i++)
				fun.val[i] = 3*C0*pow(x.val[i],2) + 2*C1*x.val[i] + C2;
		else if ( order == 2 )
			for (int i = 0; i < xlen; i++)
				fun.val[i] = 6*C0*x.val[i] + 2*C1;
		else if ( order == 3 )
			for (int i = 0; i < xlen; i++)
				fun.val[i] = 6*C0;
		else
			for (int i = 0; i < xlen; i++)
				fun.val[i] = 0.0;
	}
	else if ( strcmp(model_type, "LP") == 0 )
	{
		double f0, f1, f2, f3, f4, f5, f6, f7;
		f0 = BC.val[0];	f1 = BC.val[1];	f2 = BC.val[2];	f3 = BC.val[3];
		f4 = BC.val[4];	f5 = BC.val[5];	f6 = BC.val[6];	f7 = BC.val[7];
		double C0 = (5*f0-5*f1)/(32*pow(L,7)) + (5*f2+5*f3)/(32*pow(L,6)) + 
			(f4-f5)/(16*pow(L,5)) + (f6+f7)/(96*pow(L,4));
		double C1 = (f3-f2)/(32*pow(L,5)) + (-f4-f5)/(32*pow(L,4)) + 
			(-f6+f7)/(96*pow(L,3));
		double C2 = (-21*f0+21*f1)/(32*pow(L,5)) + (-21*f2-21*f3)/(32*pow(L,4)) + 
			(-f4+f5)/(4*pow(L,3)) + (-f6-f7)/(32*pow(L,2));
		double C3 = (5*f2-5*f3)/(32*pow(L,3)) + (5*f4+5*f5)/(32*pow(L,2)) + 
			(f6 - f7)/(32*L);
		double C4 = (35*f0-35*f1)/(32*pow(L,3)) + (35*f2+35*f3)/(32*pow(L,2)) + 
			(5*f4-5*f5)/(16*L) + (f6+f7)/32;
		double C5 = (-7*f4-7*f5)/32 + (-15*f2+15*f3)/(32*L) + (-f6+f7)*L/32;
		double C6 = (-19*f2-19*f3)/32 + (-35*f0+35*f1)/(32*L) + (-f4+f5)*L/8 + 
			(-f6-f7)*pow(L,2)/96;
		double C7 = (f0+f1)/2 + (11*f2-11*f3)*L/32 + (3*f4+3*f5)*pow(L,2)/32 + 
			(f6-f7)*pow(L,3)/96;
		if ( order == 0 )
			for (int i = 0; i < xlen; i++)
				fun.val[i] = C0*pow(x.val[i],7) + C1*pow(x.val[i],6) + 
					C2*pow(x.val[i],5) + C3*pow(x.val[i],4) + 
					C4*pow(x.val[i],3) + C5*pow(x.val[i],2) + 
					C6*x.val[i] + C7;
		else if ( order == 1 )
			for (int i = 0; i < xlen; i++)
				fun.val[i] = 7*C0*pow(x.val[i],6) + 6*C1*pow(x.val[i],5) + 
					5*C2*pow(x.val[i],4) + 4*C3*pow(x.val[i],3) + 
					3*C4*pow(x.val[i],2) + 2*C5*x.val[i] + C6;
		else if ( order == 2 )
			for (int i = 0; i < xlen; i++)
				fun.val[i] = 42*C0*pow(x.val[i],5) + 30*C1*pow(x.val[i],4) + 
					20*C2*pow(x.val[i],3) + 12*C3*pow(x.val[i],2) + 
					6*C4*x.val[i] + 2*C5;
		else if ( order == 3 )
			for (int i = 0; i < xlen; i++)
				fun.val[i] = 210*C0*pow(x.val[i],4) + 120*C1*pow(x.val[i],3) + 
					60*C2*pow(x.val[i],2) + 24*C3*x.val[i] + 6*C4;
		else if ( order == 4 )
			for (int i = 0; i < xlen; i++)
				fun.val[i] = 840*C0*pow(x.val[i],3) + 360*C1*pow(x.val[i],2) + 
					120*C2*x.val[i] + 24*C3;
		else if ( order == 5 )
			for (int i = 0; i < xlen; i++)
				fun.val[i] = 2520*C0*pow(x.val[i],2) + 720*C1*x.val[i] + 120*C2;
		else if ( order == 6 )
			for (int i = 0; i < xlen; i++)
				fun.val[i] = 5040*C0*x.val[i] + 720*C1;
		else if ( order == 7 )
			for (int i = 0; i < xlen; i++)
				fun.val[i] = 5040*C0;
		else
			for (int i = 0; i < xlen; i++)
				fun.val[i] = 0.0;
	}
}


/**
 * \brief	project rhoCplx onto the GJP space;
 *			rhoCplx: has contained the translation operator;
 *
 * \param	sbulkparam:		structure body for bulk phases;
 * \param	sGJPv:			structure body for GJP;
 */
void fn_project_GJP			(	stu_bulk_param		*sbulkparam,
									stu_GJP_var		*sGJPv	)
{
	if ( myrank == 0 )
		printf("project rhoCplx onto the GJP space.\n");

	int nd		 = sbulkparam->sbndv->nd; 
	int xlen	 = sbulkparam->sbndv->xlen;
	int	cplxDofs = sbulkparam->sbndv->cplxDofs;

	/* memory allocation; */
	fftw_complex RPBkxCplx, rhoRPBkxCplx;
	tvec<fftw_complex> xInterp		= fn_tvec_init<fftw_complex> ( xlen );
	tvec<fftw_complex> rhs			= fn_tvec_init<fftw_complex> ( nd );
	tvec<fftw_complex> rslt			= fn_tvec_init<fftw_complex> ( nd );

	/* rho on the GJP space; */
	int len = nd * cplxDofs;
	sbulkparam->sbndv->rhoJCplx = fn_tvec_init<fftw_complex> ( len );

	/* project; */
	if ( myrank == 0 )
		printf("progress: ");
	int pmod = cplxDofs / 20;
	int ind0 = 0, ind1 = 0;
	int status;
	int iter = 0;
	for (int i = 0; i < cplxDofs; i++)
	{
		for (int j = 0; j < xlen; j++)
		{
			/* exp(i(R_x'PBk)'x); */
			RPBkxCplx[0] = cos( sbulkparam->sfftv->projPlane.val[i][0] * sGJPv->x.val[j] );
			RPBkxCplx[1] = sin( sbulkparam->sfftv->projPlane.val[i][0] * sGJPv->x.val[j] );
			/* rhoCplx * exp(i(R_x'PBk)'x); */
			fn_complex_multiply ( sbulkparam->sfftv->rhoCplx.val[i], RPBkxCplx, rhoRPBkxCplx );
			/**
			 * project onto the general Jacobi polynomial space;
			 * u = phi - p; phi = rho;
			 */
			ind0 = j * cplxDofs + i;
			xInterp.val[j][0] = rhoRPBkxCplx[0] - 
									sbulkparam->sbndv->d0bnd.val[ind0][0];
			xInterp.val[j][1] = rhoRPBkxCplx[1] - 
									sbulkparam->sbndv->d0bnd.val[ind0][1];
		}

		/* calculate the inner product: (u,J); w is weight; */
		for (int j0 = 0; j0 < nd; j0++)
		{
			rhs.val[j0][0] = 0.0;		rhs.val[j0][1] = 0.0;
			for (int j1 = 0; j1 < xlen; j1++)
			{
				rhs.val[j0][0] += sGJPv->d0JJ.val[j1][j0] * 
									sGJPv->w.val[j1] * xInterp.val[j1][0];
				rhs.val[j0][1] += sGJPv->d0JJ.val[j1][j0] * 
									sGJPv->w.val[j1] * xInterp.val[j1][1];
			}
		}

		/* temporary data; */
		if ( sbulkparam->sbndv->isTest )
		{
			char dirname[FILELEN], fname[FILELEN+64];
			FILE *fwPath;
			/* xInterp; */
			sprintf(dirname, "%s/rank%d/xInterp_bulk%d", rsltDir, myrank, sbulkparam->sflag);
			mkdir(dirname, 0755);
			sprintf(fname, "%s/%d.dat", dirname, i);
			fwPath = fopen(fname, "w");
			for ( int j0 = 0; j0 < xlen; j0++ )
				fprintf(fwPath, "%+.15E\t%+.15E\n", 
						xInterp.val[j0][0], xInterp.val[j0][1]);
			fclose(fwPath);
			/* rhs; */
			sprintf(dirname, "%s/rank%d/rhs_bulk%d", rsltDir, myrank, sbulkparam->sflag);
			mkdir(dirname, 0755);
			sprintf(fname, "%s/%d.dat", dirname, i);
			fwPath = fopen(fname, "w");
			for ( int j0 = 0; j0 < nd; j0++ )
				fprintf(fwPath, "%+.15E\t%+.15E\n", rhs.val[j0][0], rhs.val[j0][1]);
			fclose(fwPath);
		}

		/* solve (u,J) = (J,J) \hat{u}; \hat{u} is coefficients; */
		status = fn_umfpack_complex_solver ( sGJPv->innSMatd0JJ, rhs, rslt );
		for (int j0 = 0; j0 < nd; j0++)
		{
			ind1 = j0 * cplxDofs + i;
			sbulkparam->sbndv->rhoJCplx.val[ind1][0] = rslt.val[j0][0];
			sbulkparam->sbndv->rhoJCplx.val[ind1][1] = rslt.val[j0][1];
		}
		if ( myrank == 0 )
			if ( i % pmod == 0 )
				printf("%.2lf ", (double) i/cplxDofs);
	}
	if ( myrank == 0 )
		printf("finish\n");

	/* save rhoJCplx; */
	if ( sbulkparam->sbndv->isSave )
	{
		char fwFile[FILELEN];
		sprintf(fwFile, "%s/rank%d/bulk%d_rhoJCplx.dat", rsltDir, myrank, sbulkparam->sflag);
		FILE *fwPath = fopen(fwFile, "w");
		fprintf(fwPath, "%d\t%d\t%d\n", nd, cplxDofs, len);
		for ( int j = 0; j < len; j++ )
		{
			fprintf(fwPath, "%+.15E\t%+.15E\n", sbulkparam->sbndv->rhoJCplx.val[j][0], 
												sbulkparam->sbndv->rhoJCplx.val[j][1]);
		}
		fclose(fwPath);
	}

	/* free memory; */
	fn_tvec_free<fftw_complex> ( xInterp );
	fn_tvec_free<fftw_complex> ( rhs );
	fn_tvec_free<fftw_complex> ( rslt );
}


/**
 * \brief	compute the mass of left and right bulk phases after projecting;
 */
void fn_obt_bulk_mass_project	(	stu_bulk_param			*sbulkparam,
									stu_system_param		*ssysparam,
									fftw_complex				mass	)
{
	/* parameters; */
	int		nd		= sbulkparam->sbndv->nd;
	int		xlen	= sbulkparam->sbndv->xlen;
	int	  cplxDofs  = sbulkparam->sbndv->cplxDofs;

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
			ind0 = j1*cplxDofs;	// the first element along Fourier direction;
			rhoCplx[0] += sbulkparam->sbndv->rhoJCplx.val[ind0][0] * ssysparam->sGJPv->d0JJ.val[i][j1];
			rhoCplx[1] += sbulkparam->sbndv->rhoJCplx.val[ind0][1] * ssysparam->sGJPv->d0JJ.val[i][j1];
		}
		ind1 = i*cplxDofs;		// the first element along Fourier direction;
		rhoCplx[0] += sbulkparam->sbndv->d0bnd.val[ind1][0];
		rhoCplx[1] += sbulkparam->sbndv->d0bnd.val[ind1][1];
		/* multiply weight; */
		mass[0] += rhoCplx[0] * ssysparam->sGJPv->w.val[i];
		mass[1] += rhoCplx[1] * ssysparam->sGJPv->w.val[i];
	}
	mass[0] *= 0.5;
	mass[1] *= 0.5;
}
