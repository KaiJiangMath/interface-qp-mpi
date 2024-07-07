/*! \file	SysInitVal.cpp
 *
 * \brief	Construct the initial value of interface system
 *			by connecting two bulk phases;
 *
 */

#include "Data.h"
#include "Head.h"
#include "DataOperators.h"
#include "Mytimer.h"
#include "functs.h"


/**
 * \brief	Initialization and preparation;
 *			Connect two bulk phases to construct the system initial value;	
 *
 * \param	sbulkparam:		the structure body for bulk phases;
 * \param	ssysparam:		the structure body for the interface system;
 *
 */
void fn_system_initial_value (	stu_bulk_param			*sbulkparam1,
								stu_bulk_param			*sbulkparam2,
								stu_system_param		*ssysparam )
{
	if ( myrank == 0 )
		printf(" <========== Generate the system initial value ==========> \n\n");
	mytimer_t timer;
	timer.reset();
	timer.start();

	/* memory allocation for 'scbndv'; */
//	ssysparam->scbndv->isTest = true;		// true: save; false: not save;
	fn_common_memory_alloc (sbulkparam1->srebndv, ssysparam->scbndv, 
							ssysparam->cplxReDofs, "part");

	/* construct initial value of the interface system; */
	fn_connect_init_value ( sbulkparam1, sbulkparam2, ssysparam );

	/* connect boundary polynomials of the left/right bulk phases; */
	fn_connect_rebnd ( sbulkparam1->srebndv, sbulkparam2->srebndv, ssysparam );

	/* save data and plot; */
	fn_save_system_phase  ( ssysparam, 0 );

	/* release memory about 'srebndv'; */
	fn_total_bulk_memory_free ( sbulkparam1 );
	fn_total_bulk_memory_free ( sbulkparam2 );

	MPI_Barrier ( MPI_COMM_WORLD );
	timer.pause();
	if ( myrank == 0 )
	{
		printf("\t ***** time cost of generation of system initial value: ");
		printf("%f seconds *****\n\n", timer.get_current_time());
	}
}


/**
 * \brief	Construct initial value;
 *			connect the two bulk phases (in the common space);
 */
void fn_connect_init_value (	stu_bulk_param		*sbulkparam1,
								stu_bulk_param		*sbulkparam2,
								stu_system_param	*ssysparam )
{
	/* load parameters; */
	int nd			= ssysparam->scbndv->nd; 
	int xlen		= ssysparam->scbndv->xlen;
	int cplxReDofs	= ssysparam->cplxReDofs;
	int dimRePhy	= ssysparam->dimRePhy;
	int dimReCpt	= ssysparam->dimReCpt;
	double smooth	= ssysparam->smooth;
	double initDist1 = ssysparam->initDist1;	// the end position of left bulk phase;
	double initDist2 = ssysparam->initDist2;	// the end position of right bulk phase;

	/* structure body; */
	tvec<fftw_complex> xInterp	= fn_tvec_init<fftw_complex> ( xlen );
	tvec<fftw_complex> rhs		= fn_tvec_init<fftw_complex> ( nd );
	tvec<fftw_complex> rslt		= fn_tvec_init<fftw_complex> ( nd );

	/* connecting; */
	int status;
	int ind0 = 0, ind1 = 0;
	double maskL, maskR;	// temporary variables for the smooth function;
	double tmpReal1, tmpImag1, tmpReal2, tmpImag2;
	for ( int i = 0; i < alloc_local_sys; i++ )
	{
		for ( int j = 0; j < xlen; j++ )
		{
			/* calculate connection functions; */
			maskL = 0.5 * ( 1.0 - tanh( smooth*(ssysparam->sGJPv->x.val[j] - initDist1) ) );
			maskR = 0.5 * ( 1.0 + tanh( smooth*(ssysparam->sGJPv->x.val[j] - initDist2) ) );
			
			/* connect initial values; */
			/* (1-b(x)) * u1(x); */
			tmpReal1 = 0.0;		tmpImag1 = 0.0;
			for ( int j0 = 0; j0 < nd; j0++ )
			{
				ind0 = j0 * alloc_local_sys + i;
				tmpReal1 += ssysparam->sGJPv->d0JJ.val[j][j0] * 
							sbulkparam1->srebndv->rhoJCplx.val[ind0][0];
				tmpImag1 += ssysparam->sGJPv->d0JJ.val[j][j0] * 
							sbulkparam1->srebndv->rhoJCplx.val[ind0][1];
			}
			tmpReal1 *= maskL;
			tmpImag1 *= maskL;
			/* b(x) * u2(x); */
			tmpReal2 = 0.0;		tmpImag2 = 0.0;
			for ( int j0 = 0; j0 < nd; j0++ )
			{
				ind0 = j0 * alloc_local_sys + i;
				tmpReal2 += ssysparam->sGJPv->d0JJ.val[j][j0] * 
							sbulkparam2->srebndv->rhoJCplx.val[ind0][0];
				tmpImag2 += ssysparam->sGJPv->d0JJ.val[j][j0] * 
							sbulkparam2->srebndv->rhoJCplx.val[ind0][1];
			}
			tmpReal2 *= maskR;
			tmpImag2 *= maskR;
			/* (1-b(x)) * u1(x) + b(x) * u2(x); */
			xInterp.val[j][0] = tmpReal1 + tmpReal2;
			xInterp.val[j][1] = tmpImag1 + tmpImag2;
		}

		/* calculate the right term; */
		for ( int j = 0; j < nd; j++ )
		{
			fn_complex_setZero ( rhs.val[j] );
			for ( int j0 = 0; j0 < xlen; j0++ )
			{
				rhs.val[j][0] += ssysparam->sGJPv->d0JJ.val[j0][j] * 
								ssysparam->sGJPv->w.val[j0] * xInterp.val[j0][0];
				rhs.val[j][1] += ssysparam->sGJPv->d0JJ.val[j0][j] * 
								ssysparam->sGJPv->w.val[j0] * xInterp.val[j0][1];
			}
		}

		/* project onto the space of GJPs;
		 *	solve (u,J) = (J,J) \hat{u}; \hat{u} is coefficients;
		 */
		status = fn_umfpack_complex_solver ( ssysparam->sGJPv->innSMatd0JJ, rhs, rslt );
		for (int j0 = 0; j0 < nd; j0++)
		{
			ind1 = j0 * alloc_local_sys + i;
			ssysparam->scbndv->rhoJCplx.val[ind1][0] = rslt.val[j0][0];
			ssysparam->scbndv->rhoJCplx.val[ind1][1] = rslt.val[j0][1];
		}
	}

	/* Release memory; */
	fn_tvec_free<fftw_complex> ( xInterp );
	fn_tvec_free<fftw_complex> ( rhs );
	fn_tvec_free<fftw_complex> ( rslt );
}


/**
 * \brief	connect the polynomials (left/right) by the smooth function;
 */
void fn_connect_rebnd  (	stu_bnd_var			*srebndv1,
							stu_bnd_var			*srebndv2,
							stu_system_param	*ssysparam )
{
	/* load parameters; */
	int nd			= ssysparam->scbndv->nd; 
	int xlen		= ssysparam->scbndv->xlen;
	int cplxReDofs	= ssysparam->cplxReDofs;
	int dimRePhy	= ssysparam->dimRePhy;
	int dimReCpt	= ssysparam->dimReCpt;
	double smooth	= ssysparam->smooth;
	double initDist1 = ssysparam->initDist1;	// the end position of left bulk phase;
	double initDist2 = ssysparam->initDist2;	// the end position of right bulk phase;

	/* temporary variables for the smooth function; */
	int maskLen = 0;
	if ( strcmp(model_type, "LB") == 0 )
		maskLen = 3;
	else if ( strcmp(model_type, "LP") == 0 )
		maskLen = 7;
	else
	{
		if ( myrank == 0 )
		{
			printf("Error using 'fn_connect_rebnd'");
			printf("maskLen error.\n");
		}
	}
	tvec<double> maskL = fn_tvec_init<double> ( maskLen );
	tvec<double> maskR = fn_tvec_init<double> ( maskLen );

	/* connecting; */
	int status;
	int ind0 = 0, ind1 = 0;
	double tmpReal1, tmpImag1, tmpReal2, tmpImag2;
	for ( int i = 0; i < alloc_local_sys; i++ )
	{
		for ( int j = 0; j < xlen; j++ )
		{
			/* calculate connection functions; */
			fn_obt_connect_fun ( ssysparam->sGJPv->x.val[j], smooth, initDist1, 1, maskL );
			fn_obt_connect_fun ( ssysparam->sGJPv->x.val[j], smooth, initDist2, 2, maskR );
			
			/**
			 * calculate 'd0bnd', 'd2bnd', 'd4bnd', 'd6bnd' in 'ssysparam';
			 *		'd1bnd', 'd3bnd',... will not be used in 'LB'/'LP';
			 * connect 'd0bnd',... of 'srebnd1' and 'srebnd2';
			 */
			ind0 = j * alloc_local_sys + i;
			/* d0bnd; */
			ssysparam->scbndv->d0bnd.val[ind0][0] = 
					maskL.val[0]*srebndv1->d0bnd.val[ind0][0] + 
					maskR.val[0]*srebndv2->d0bnd.val[ind0][0];
			ssysparam->scbndv->d0bnd.val[ind0][1] = 
					maskL.val[0]*srebndv1->d0bnd.val[ind0][1] + 
					maskR.val[0]*srebndv2->d0bnd.val[ind0][1];
			/* d2bnd; */
			ssysparam->scbndv->d2bnd.val[ind0][0] = 
					maskL.val[2]*srebndv1->d0bnd.val[ind0][0] + 
				2.0*maskL.val[1]*srebndv1->d1bnd.val[ind0][0] + 
					maskL.val[0]*srebndv1->d2bnd.val[ind0][0] +
					maskR.val[2]*srebndv2->d0bnd.val[ind0][0] + 
				2.0*maskR.val[1]*srebndv2->d1bnd.val[ind0][0] + 
					maskR.val[0]*srebndv2->d2bnd.val[ind0][0];
			ssysparam->scbndv->d2bnd.val[ind0][1] = 
					maskL.val[2]*srebndv1->d0bnd.val[ind0][1] + 
				2.0*maskL.val[1]*srebndv1->d1bnd.val[ind0][1] + 
					maskL.val[0]*srebndv1->d2bnd.val[ind0][1] +
					maskR.val[2]*srebndv2->d0bnd.val[ind0][1] + 
				2.0*maskR.val[1]*srebndv2->d1bnd.val[ind0][1] + 
					maskR.val[0]*srebndv2->d2bnd.val[ind0][1];
			if ( strcmp(model_type, "LP") == 0 )
			{
				/* d4bnd; */
				ssysparam->scbndv->d4bnd.val[ind0][0] =
						maskL.val[4]*srebndv1->d0bnd.val[ind0][0] + 
					4.0*maskL.val[3]*srebndv1->d1bnd.val[ind0][0] +
					6.0*maskL.val[2]*srebndv1->d2bnd.val[ind0][0] + 
					4.0*maskL.val[1]*srebndv1->d3bnd.val[ind0][0] +
						maskL.val[0]*srebndv1->d4bnd.val[ind0][0] +
						maskR.val[4]*srebndv2->d0bnd.val[ind0][0] + 
					4.0*maskR.val[3]*srebndv2->d1bnd.val[ind0][0] +
					6.0*maskR.val[2]*srebndv2->d2bnd.val[ind0][0] + 
					4.0*maskR.val[1]*srebndv2->d3bnd.val[ind0][0] +
						maskR.val[0]*srebndv2->d4bnd.val[ind0][0];
				ssysparam->scbndv->d4bnd.val[ind0][1] =
						maskL.val[4]*srebndv1->d0bnd.val[ind0][1] + 
					4.0*maskL.val[3]*srebndv1->d1bnd.val[ind0][1] +
					6.0*maskL.val[2]*srebndv1->d2bnd.val[ind0][1] + 
					4.0*maskL.val[1]*srebndv1->d3bnd.val[ind0][1] +
						maskL.val[0]*srebndv1->d4bnd.val[ind0][1] +
						maskR.val[4]*srebndv2->d0bnd.val[ind0][1] + 
					4.0*maskR.val[3]*srebndv2->d1bnd.val[ind0][1] +
					6.0*maskR.val[2]*srebndv2->d2bnd.val[ind0][1] + 
					4.0*maskR.val[1]*srebndv2->d3bnd.val[ind0][1] +
						maskR.val[0]*srebndv2->d4bnd.val[ind0][1];
				/* d6bnd; */
				ssysparam->scbndv->d6bnd.val[ind0][0] =
						 maskL.val[6]*srebndv1->d0bnd.val[ind0][0] +  
					 6.0*maskL.val[5]*srebndv1->d1bnd.val[ind0][0] +
					15.0*maskL.val[4]*srebndv1->d2bnd.val[ind0][0] + 
					20.0*maskL.val[3]*srebndv1->d3bnd.val[ind0][0] +
					15.0*maskL.val[2]*srebndv1->d4bnd.val[ind0][0] +
					 6.0*maskL.val[1]*srebndv1->d5bnd.val[ind0][0] +
						 maskL.val[0]*srebndv1->d6bnd.val[ind0][0] +
						 maskR.val[6]*srebndv2->d0bnd.val[ind0][0] +
					 6.0*maskR.val[5]*srebndv2->d1bnd.val[ind0][0] +
					15.0*maskR.val[4]*srebndv2->d2bnd.val[ind0][0] + 
					20.0*maskR.val[3]*srebndv2->d3bnd.val[ind0][0] +
					15.0*maskR.val[2]*srebndv2->d4bnd.val[ind0][0] +
					 6.0*maskR.val[1]*srebndv2->d5bnd.val[ind0][0] +
						 maskR.val[0]*srebndv2->d6bnd.val[ind0][0];
				ssysparam->scbndv->d6bnd.val[ind0][1] =
						 maskL.val[6]*srebndv1->d0bnd.val[ind0][1] +
					 6.0*maskL.val[5]*srebndv1->d1bnd.val[ind0][1] +
					15.0*maskL.val[4]*srebndv1->d2bnd.val[ind0][1] +
					20.0*maskL.val[3]*srebndv1->d3bnd.val[ind0][1] +
					15.0*maskL.val[2]*srebndv1->d4bnd.val[ind0][1] +
					 6.0*maskL.val[1]*srebndv1->d5bnd.val[ind0][1] +
						 maskL.val[0]*srebndv1->d6bnd.val[ind0][1] +
						 maskR.val[6]*srebndv2->d0bnd.val[ind0][1] +
					 6.0*maskR.val[5]*srebndv2->d1bnd.val[ind0][1] +
					15.0*maskR.val[4]*srebndv2->d2bnd.val[ind0][1] +
					20.0*maskR.val[3]*srebndv2->d3bnd.val[ind0][1] +
					15.0*maskR.val[2]*srebndv2->d4bnd.val[ind0][1] +
					 6.0*maskR.val[1]*srebndv2->d5bnd.val[ind0][1] +
						 maskR.val[0]*srebndv2->d6bnd.val[ind0][1];
			}
		}
	}

	/* save results; */
	fn_save_cbnd ( ssysparam->scbndv->d0bnd, 0, xlen, alloc_local_sys );
	if ( ssysparam->scbndv->isTest )
	{
		fn_save_cbnd ( ssysparam->scbndv->d2bnd, 2, xlen, alloc_local_sys );
		if ( strcmp(model_type, "LP") == 0 )
		{
			fn_save_cbnd ( ssysparam->scbndv->d4bnd, 4, xlen, alloc_local_sys );
			fn_save_cbnd ( ssysparam->scbndv->d6bnd, 6, xlen, alloc_local_sys );
		}
	}

	/* Release memory; */
	fn_tvec_free<double> ( maskL );
	fn_tvec_free<double> ( maskR );
}


/**
 * \brief	Generate the connection function;
 *
 * \param	xVal:		the value of x;
 * \param	smooth, initDist:		parameters controlling the smooth function;
 * \param	sflag:		1: left; 2: right;
 * \param	masklr:		return value;
 */
void fn_obt_connect_fun (	double			xVal,
							double			smooth,
							double			initDist, 
							int				sflag,
							tvec<double>	masklr )
{
	int coeff	= pow(-1, sflag);	// minus or plus;
	double xlr	= smooth * (xVal - initDist);
	double mxlr = tanh(xlr);
	masklr.val[0] = 0.5 * ( 1.0 + coeff*mxlr );
	masklr.val[1] = coeff * 0.5*smooth * ( 1.0 - pow(tanh(xlr),2) );
	masklr.val[2] = coeff * pow(smooth,2) * ( pow(mxlr,3) - mxlr );
	if ( strcmp(model_type, "LP") == 0 )
	{
		masklr.val[3] = -coeff * pow(smooth,3) * 
			( 3.0*pow(mxlr,4) - 4.0*pow(mxlr,2) + 1.0 );
		masklr.val[4] =  coeff * 4.0*pow(smooth,4) * mxlr * 
			( 3.0*pow(mxlr,4) - 5.0*pow(mxlr,2) + 2.0 );
		masklr.val[5] = -coeff * pow(smooth,5) * ( 60.0*pow(mxlr,6) - 
				120.0*pow(mxlr,4) + 68.0*pow(mxlr,2) - 8.0 );
		masklr.val[6] =  coeff * 8.0*pow(smooth,6) * mxlr * 
			( 45.0*pow(mxlr,6) - 105.0*pow(mxlr,4) + 77.0*pow(mxlr,2) - 17.0 );
	}
}


/**
 * \brief	save 'd0bnd', 'd2bnd', ... of the structure body 'scbndv';
 *
 * \param	order:		the 'order'-th derivative;
 * \param	xlen:		the number of rows;
 * \param	cplxReDofs:	the number of columns;
 */
void fn_save_cbnd (	tvec<fftw_complex>	dsbnd,	int order,	int xlen, int cplxReDofs)
{
	int len = dsbnd.len;

	char fwFile[FILELEN];
	sprintf(fwFile, "%s/rank%d/sys_d%dbnd.dat", rsltDir, myrank, order);

	FILE *fwPath = fopen(fwFile, "w");
	fprintf(fwPath, "%d\t%d\t%d\n", xlen, cplxReDofs, len);
	for ( int j = 0; j < len; j++ )
		fprintf(fwPath, "%+.15E\t%+.15E\n", dsbnd.val[j][0], dsbnd.val[j][1]);
	fclose(fwPath);
}


/**
 * \brief	Release memory about bulk phases;
 */
void fn_total_bulk_memory_free	(stu_bulk_param		*sbulkparam)
{
	fn_bnd_memory_free  ( sbulkparam->srebndv, "total" );
	fn_bulk_memory_free ( sbulkparam, "project bulk" );
	fn_bulk_memory_free ( sbulkparam, "bulk param" );
}
