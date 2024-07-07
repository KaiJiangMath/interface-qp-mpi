/*! \file	MemFree.cpp
 *
 * \brief	Release memory;
 */

#include "Head.h"
#include "Data.h"
#include "DataOperators.h"
#include "functs.h"


/**
 * \brief	Release memory about bulk phases;
 *
 * \param	sbulkparam:		the structure body for bulk phases;
 * \param	memType:		the type of memory;
 */
void fn_bulk_memory_free	( stu_bulk_param	*sbulkparam,
								const char		*memType )
{
	if ( strcmp(memType, "bulk param") == 0 )
	{
		fn_tvec_free<int>		 ( sbulkparam->NCpt );
		fn_tvec_free<ptrdiff_t>	 ( sbulkparam->pNCpt );
		fn_tvec_free<double> ( sbulkparam->rotate_angle );
		fn_tvec_free<double> ( sbulkparam->transl_var );

		fn_tmat_free<int>	 ( sbulkparam->initIndex );
		fn_tvec_free<fftw_complex> ( sbulkparam->initCoeff );

		fn_tmat_free<double> ( sbulkparam->rcpBox	 );
		fn_tmat_free<double> ( sbulkparam->dirBox	 );
		fn_tmat_free<double> ( sbulkparam->projMat	 );
		fn_tmat_free<double> ( sbulkparam->rotateMat );
		fn_tmat_free<double> ( sbulkparam->rotateProjBoxMat );

		fn_tvec_free<double> ( sbulkparam->translVec );
		fn_tvec_free<double> ( sbulkparam->translProjBoxVec );
	}
	else if ( strcmp(memType, "stable bulk") == 0 )
	{
		fn_tvec_free<fftw_complex> ( sbulkparam->sfftv->rhoReal );
		fn_tvec_free<fftw_complex> ( sbulkparam->sfftv->fftw_Ctmp );
		fn_tvec_free<fftw_complex> ( sbulkparam->sfftv->fftw_Rtmp );
		fn_tvec_free<fftw_complex> ( sbulkparam->sfftv->gradient );
		fn_tvec_free<fftw_complex> ( sbulkparam->sfftv->cplxTmp );
		fn_tvec_free<fftw_complex> ( sbulkparam->sfftv->quadTerm );
		fn_tvec_free<fftw_complex> ( sbulkparam->sfftv->cubTerm );
		fn_tvec_free<fftw_complex> ( sbulkparam->sfftv->quarTerm );

		fn_tvec_free<double> ( sbulkparam->sfftv->Gsquare );

		fftw_destroy_plan	 ( sbulkparam->sfftv->Planc2cFord );
		fftw_destroy_plan	 ( sbulkparam->sfftv->Planc2cBack );
	}
	else if ( strcmp(memType, "project bulk") == 0 )
	{
		fn_tvec_free<fftw_complex>	( sbulkparam->sfftv->rhoCplx );
		fn_tmat_free<int>			( sbulkparam->sfftv->indKspace );
		fn_tmat_free<double>		( sbulkparam->sfftv->projPlane );
		fn_tvec_free<int>			( sbulkparam->sfftv->ind );
	}
}


/**
 * \brief	Release memory after 'CommonFourGJP';
 */
void fn_bnd_memory_free		   (	stu_bnd_var		*sbndv,
									const char		*memType )
{
	fn_tvec_free<fftw_complex> ( sbndv->rhoJCplx );

	fn_tvec_free<fftw_complex> ( sbndv->d0bnd );
	fn_tvec_free<fftw_complex> ( sbndv->d2bnd );
	if ( strcmp(model_type, "LP") == 0 )
	{
		fn_tvec_free<fftw_complex> ( sbndv->d4bnd );
		fn_tvec_free<fftw_complex> ( sbndv->d6bnd );
	}

	if ( strcmp(memType, "total") == 0 )
	{
		fn_tvec_free<fftw_complex> ( sbndv->d1bnd );
		if ( strcmp(model_type, "LP") == 0 )
		{
			fn_tvec_free<fftw_complex> ( sbndv->d3bnd );
			fn_tvec_free<fftw_complex> ( sbndv->d5bnd );
		}
	}
}


/**
 * \brief	Release memory about the interface system;
 */
void fn_system_memory_free ( stu_system_param		*ssysparam )
{
	fn_tvec_free<ptrdiff_t>			( ssysparam->NCpt );

	fn_tmat_free<double>		( ssysparam->rotateProjBoxMat );
	fn_tmat_free<int>			( ssysparam->coeffmat );

	fn_tmat_free<int>			( ssysparam->sfftv->indKspace );
	fn_tmat_free<double>		( ssysparam->sfftv->projPlane );
	fn_tvec_free<double>		( ssysparam->sfftv->Gsquare );

	fn_tvec_free<fftw_complex>	( ssysparam->sfftv->rhoCplx );
	fn_tvec_free<fftw_complex>	( ssysparam->sfftv->rhoReal );
	fn_tvec_free<fftw_complex>	( ssysparam->sfftv->fftw_Ctmp );
	fn_tvec_free<fftw_complex>	( ssysparam->sfftv->fftw_Rtmp );
	fn_tvec_free<fftw_complex>	( ssysparam->sfftv->gradient );
	fn_tvec_free<fftw_complex>	( ssysparam->sfftv->cplxTmp );
	fn_tvec_free<fftw_complex>	( ssysparam->sfftv->quadTerm );
	fn_tvec_free<fftw_complex>	( ssysparam->sfftv->cubTerm );
	fn_tvec_free<fftw_complex>	( ssysparam->sfftv->quarTerm );

	fftw_destroy_plan(ssysparam->sfftv->Planc2cFord);
	fftw_destroy_plan(ssysparam->sfftv->Planc2cBack);

	fn_tvec_free<fftw_complex> 	( ssysparam->iter_rho_rhs );
	fn_tvec_free<fftw_complex> 	( ssysparam->iter_entropy_rhs );
	fn_tvec_free<fftw_complex> 	( ssysparam->interact_bnd_grad );
	fn_tvec_free<fftw_complex>	( ssysparam->grad_err );
	fn_tvec_free<fftw_complex>	( ssysparam->hessian );

	fn_tCCSmat_free<double>		( ssysparam->iter_matrix );
	fn_tCCSmat_free<double>		( ssysparam->interact_grad );
}


/**
 * \brief	Release memory about GJPs;
 */
void fn_GJP_memory_free ( stu_GJP_var	*sGJPv )
{
	fn_tvec_free<double> ( sGJPv->x );
	fn_tvec_free<double> ( sGJPv->w );

	fn_tmat_free<double> ( sGJPv->d0JJ );
	fn_tmat_free<double> ( sGJPv->d1JJ );
	fn_tmat_free<double> ( sGJPv->d2JJ );
	fn_tmat_free<double> ( sGJPv->d3JJ );
	fn_tmat_free<double> ( sGJPv->d4JJ );

	fn_tCCSmat_free<double> ( sGJPv->innSMatd0JJ );
	fn_tCCSmat_free<double> ( sGJPv->innSMatd1JJ );
	fn_tCCSmat_free<double> ( sGJPv->innSMatd2JJ );
	if ( strcmp(model_type, "LP") == 0 )
	{
		fn_tCCSmat_free<double> ( sGJPv->innSMatd3JJ );
		fn_tCCSmat_free<double> ( sGJPv->innSMatd4JJ );
	}
}
