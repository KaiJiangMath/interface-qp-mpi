/*! \file	MainCollect.cpp
 *
 *  \brief	Collection of main code to realize different functions;
 *
 */

#include "Head.h"
#include "Data.h"
#include "DataOperators.h"
#include "Mytimer.h"
#include "functs.h"
#include "umfpack.h"

/*---------------------------------*/
/*--       Main Functions        --*/
/*---------------------------------*/

/**
 * \brief	Calculation of the stable interface structure;
 */
double fn_main_stable_interface	( )
{
	/* introduce structure bodies; */
	stu_bulk_param sbulkparam1Type;
	stu_bulk_param *sbulkparam1 = &sbulkparam1Type;
	stu_bulk_param sbulkparam2Type;
	stu_bulk_param *sbulkparam2 = &sbulkparam2Type;
	stu_system_param ssysparamType;
	stu_system_param *ssysparam = &ssysparamType;

	/* true: save;	false: not save; */
	/* whether save densities; */
	sbulkparam1->sbndv->isSave	 = true; // 'fn_project_Fourier_GJP';
	sbulkparam2->sbndv->isSave	 = true; // 'fn_project_Fourier_GJP';
	sbulkparam1->srebndv->isSave = true; // 'fn_common_Fourier_GJP';
	sbulkparam2->srebndv->isSave = true; // 'fn_common_Fourier_GJP';
	/* variables for debug; */
	ssysparam->sGJPv->isTest	 = false; // 'fn_obt_system_gen_Jac_poly';
	sbulkparam1->sbndv->isTest	 = false; // 'fn_project_Fourier_GJP';
	sbulkparam2->sbndv->isTest	 = false; // 'fn_project_Fourier_GJP';
	sbulkparam1->srebndv->isTest = false; // 'fn_common_Fourier_GJP';
	sbulkparam2->srebndv->isTest = false; // 'fn_common_Fourier_GJP';
	ssysparam->scbndv->isTest	 = false; // 'fn_system_initial_value';
	bool iter_prepare_isTest	 = false; // 'fn_iter_prepare';
	bool iter_method_isTest		 = false; // 'fn_choose_method';

	/* input and save parameters; */
	fn_param_input ( sbulkparam1, sbulkparam2, ssysparam );

	/* obtain the stable bulk phases; */
	double hamilton1 = fn_obt_stable_bulk_phase ( sbulkparam1 );
	double hamilton2 = fn_obt_stable_bulk_phase ( sbulkparam2 );

	/* update 'x_range'; */
	if ( strcmp ( ssysparam->x_range_type, "direct" ) == 0 )
	{
		if ( myrank == 0 )
			printf("x_range_type = %s.\n", ssysparam->x_range_type);
	}
	else if ( strcmp ( ssysparam->x_range_type, "cubePlane" ) == 0 )
	{
		double angle = fabs ( sbulkparam1->rotate_angle.val[0] ) / 180.0 * PI;
		double tmp1 = sin(angle) + cos(angle);
		tmp1 *= sbulkparam1->dirBox.val[0][0];
		ssysparam->x_range *= tmp1;
	}
	else
	{
		if ( myrank == 0 )
			printf("Error x_range_type %s.\n", ssysparam->x_range_type);
	}
	ssysparam->sGJPv->x_range = ssysparam->x_range;

	/* prepare general Jacobi polynomials; */
	fn_obt_system_gen_Jac_poly ( ssysparam->sGJPv );

	/* project rhoCplx to the space of GJPs; */
	fn_project_Fourier_GJP ( sbulkparam1, ssysparam );
	fn_project_Fourier_GJP ( sbulkparam2, ssysparam );

	/* obtain the common 'rotateProjBoxMat'; */
	int status = fn_obt_commom_rotateProjBoxMat ( sbulkparam1, sbulkparam2, ssysparam );
	if ( status != 0 )
		return 0.0;

	/* save parameters (optimal computational box); */
	fn_param_save ( sbulkparam1, sbulkparam2, ssysparam, "opt" );

	/* to a common space; */
	fn_common_Fourier_GJP   ( sbulkparam1, sbulkparam2, ssysparam	);

	/* connect two bulk phases to construct initial value; */
	fn_system_initial_value ( sbulkparam1, sbulkparam2, ssysparam	);

	/* iteration preparation; */
	fn_iter_prepare	( ssysparam, iter_prepare_isTest );

	/* iteration by numerical methods; */
	double hamilton = fn_choose_method ( ssysparam, iter_method_isTest, 0 );

	/* memory free about the interface system; */
	fn_GJP_memory_free	  ( ssysparam->sGJPv );
	fn_bnd_memory_free	  ( ssysparam->scbndv, "part" );
	fn_system_memory_free ( ssysparam );
	free ( displs_sys );		free ( recvCount_sys );

	return hamilton;
}


/**
 * \brief	Calculation of the stable interface structure 
 *			by loading existing files;
 */
double fn_main_load_stable_interface	(  )
{
	/* introduce structure bodies; */
	stu_bulk_param sbulkparam1Type;
	stu_bulk_param *sbulkparam1 = &sbulkparam1Type;
	stu_bulk_param sbulkparam2Type;
	stu_bulk_param *sbulkparam2 = &sbulkparam2Type;
	stu_system_param ssysparamType;
	stu_system_param *ssysparam = &ssysparamType;

	/* true: save;	false: not save; */
	/* whether save densities; */
	sbulkparam1->sbndv->isSave	 = true; // 'fn_project_Fourier_GJP';
	sbulkparam2->sbndv->isSave	 = true; // 'fn_project_Fourier_GJP';
	sbulkparam1->srebndv->isSave = true; // 'fn_common_Fourier_GJP';
	sbulkparam2->srebndv->isSave = true; // 'fn_common_Fourier_GJP';
	/* variables for debug; */
	ssysparam->sGJPv->isTest	 = false; // 'fn_obt_system_gen_Jac_poly';
	sbulkparam1->sbndv->isTest	 = false; // 'fn_project_Fourier_GJP';
	sbulkparam2->sbndv->isTest	 = false; // 'fn_project_Fourier_GJP';
	sbulkparam1->srebndv->isTest = false; // 'fn_common_Fourier_GJP';
	sbulkparam2->srebndv->isTest = false; // 'fn_common_Fourier_GJP';
	ssysparam->scbndv->isTest	 = false; // 'fn_system_initial_value';
	bool iter_prepare_isTest	 = false; // 'fn_iter_prepare';
	bool iter_method_isTest		 = false; // 'fn_choose_method';


	/* input and save parameters; */
	fn_param_input ( sbulkparam1, sbulkparam2, ssysparam );

	/* obtain the stable bulk phases; */
	double hamilton1 = fn_obt_stable_bulk_phase ( sbulkparam1 );
	double hamilton2 = fn_obt_stable_bulk_phase ( sbulkparam2 );

	/* update 'x_range'; */
	if ( strcmp ( ssysparam->x_range_type, "direct" ) == 0 )
	{
		if ( myrank == 0 )
			printf("x_range_type = %s.\n", ssysparam->x_range_type);
	}
	else if ( strcmp ( ssysparam->x_range_type, "cubePlane" ) == 0 )
	{
		double angle = fabs ( sbulkparam1->rotate_angle.val[0] ) / 180.0 * PI;
		double tmp1 = sin(angle) + cos(angle);
		tmp1 *= sbulkparam1->dirBox.val[0][0];
		ssysparam->x_range *= tmp1;
	}
	else
	{
		if ( myrank == 0 )
			printf("Error x_range_type %s.\n", ssysparam->x_range_type);
	}
	ssysparam->sGJPv->x_range = ssysparam->x_range;

	/* prepare general Jacobi polynomials; */
	fn_obt_system_gen_Jac_poly ( ssysparam->sGJPv );

	/* project rhoCplx to the space of GJPs; */
	fn_project_Fourier_GJP ( sbulkparam1, ssysparam );
	fn_project_Fourier_GJP ( sbulkparam2, ssysparam );

	/* obtain the common 'rotateProjBoxMat'; */
	int status = fn_obt_commom_rotateProjBoxMat ( sbulkparam1, sbulkparam2, ssysparam );
	if ( status != 0 )
		return 0.0;

	/* save parameters (optimal computational box); */
	fn_param_save ( sbulkparam1, sbulkparam2, ssysparam, "opt" );

	/* to a common space; */
	fn_common_Fourier_GJP   ( sbulkparam1, sbulkparam2, ssysparam	);

	/* connect two bulk phases to construct initial value; */
	fn_system_initial_value ( sbulkparam1, sbulkparam2, ssysparam	);

	/* iteration preparation; */
	fn_iter_prepare	( ssysparam, iter_prepare_isTest );

	/* load 'rhoJCplx'; */
	char frFile[FILELEN];
	sprintf(frFile, "%s/rank%d/sys_rhoJCplx%d.dat", rsltDir, myrank, ssysparam->iter_load);
    FILE *fp = fopen(frFile,"r");
    if ( fp == NULL )	// invalid file;
	{
		if ( myrank == 0 )
		{
			printf("Error using 'fn_main_load_stable_interface'\n");
			printf("Unable to read file '%s'. No such file or directory.\n", frFile);
		}
		return 0.0;
	}
	else
	{
		int val;
		val = fscanf(fp, "%d", &(ssysparam->scbndv->rhoJCplx.row));
		val = fscanf(fp, "%d", &(ssysparam->scbndv->rhoJCplx.col));
		val = fscanf(fp, "%d", &(ssysparam->scbndv->rhoJCplx.len));
		for ( int i = 0; i < ssysparam->scbndv->rhoJCplx.len; i++ )
		{
			val = fscanf(fp, "%lf", &(ssysparam->scbndv->rhoJCplx.val[i][0]));
			val = fscanf(fp, "%lf", &(ssysparam->scbndv->rhoJCplx.val[i][1]));
		}
	}
	fclose(fp);	

	/* iteration by numerical methods; */
	int	iter_load = ssysparam->iter_load;
	double hamilton = fn_choose_method ( ssysparam, iter_method_isTest, iter_load );

	/* memory free about the interface system; */
	fn_GJP_memory_free	  ( ssysparam->sGJPv );
	fn_bnd_memory_free	  ( ssysparam->scbndv, "part" );
	fn_system_memory_free ( ssysparam );
	free ( displs_sys );		free ( recvCount_sys );

	return hamilton;
}


/**
 * \brief	Calculation of the stable bulk phase;
 */
double fn_main_stable_bulk ( )
{
	/* introduce structure bodies; */
	stu_bulk_param sbulkparam1Type;
	stu_bulk_param *sbulkparam1 = &sbulkparam1Type;
	stu_bulk_param sbulkparam2Type;
	stu_bulk_param *sbulkparam2 = &sbulkparam2Type;
	stu_system_param ssysparamType;
	stu_system_param *ssysparam = &ssysparamType;

	/* input and save parameters; */
	fn_param_input ( sbulkparam1, sbulkparam2, ssysparam );

	/* obtain the stable bulk phases; */
	double hamilton = 0.0;
	if ( strcmp( main_type, "stable_bulk1" ) == 0 )
	{
		hamilton = fn_obt_stable_bulk_phase ( sbulkparam1 );
		fn_bulk_memory_free ( sbulkparam1, "project bulk" );
	}
	else if ( strcmp( main_type, "stable_bulk2" ) == 0 )
	{
		hamilton = fn_obt_stable_bulk_phase ( sbulkparam2 );
		fn_bulk_memory_free ( sbulkparam2, "project bulk" );
	}
	else
	{
		if ( myrank == 0 )
		{
			printf("Error use 'fn_main_stable_bulk'.\n");
			printf("'main_type' should be 'stable_bulk1' or 'stable_bulk2'.\n");
		}
	}

	/* release memory; */
	fn_bulk_memory_free	( sbulkparam1, "bulk param" );
	fn_bulk_memory_free	( sbulkparam2, "bulk param" );

	return hamilton;
}


/**
 * \brief	Calculation of the stable bulk phase by a parallel code;
 */
double fn_main_stable_bulk_parallel ( )
{
	/* introduce structure bodies; */
	stu_bulk_param sbulkparam1Type;
	stu_bulk_param *sbulkparam1 = &sbulkparam1Type;
	stu_bulk_param sbulkparam2Type;
	stu_bulk_param *sbulkparam2 = &sbulkparam2Type;
	stu_system_param ssysparamType;
	stu_system_param *ssysparam = &ssysparamType;

	/* input and save parameters; */
	fn_param_input ( sbulkparam1, sbulkparam2, ssysparam );

	/* obtain the stable bulk phases; */
    displs_bulk    = (int *) malloc ( sizeof(int) * nprocs ); // offset of each process;
    recvCount_bulk = (int *) malloc ( sizeof(int) * nprocs ); // offset of each process;
	double hamilton = 0.0;
	if ( strcmp( main_type, "stable_bulk1_parallel" ) == 0 )
	{
		hamilton = fn_obt_stable_bulk_parallel ( sbulkparam1 );
		if ( myrank == 0 )
			fn_bulk_memory_free ( sbulkparam1, "project bulk" );
	}
	else if ( strcmp( main_type, "stable_bulk2_parallel" ) == 0 )
	{
		hamilton = fn_obt_stable_bulk_parallel ( sbulkparam2 );
		if ( myrank == 0 )
			fn_bulk_memory_free ( sbulkparam2, "project bulk" );
	}
	else
	{
		if ( myrank == 0 )
		{
			printf("Error use 'fn_main_stable_bulk_parallel'.\n");
			printf("'main_type' should be 'stable_bulk1_parallel' or 'stable_bulk2_parallel'.\n");
		}
	}

	/* release memory; */
	fn_bulk_memory_free	( sbulkparam1, "bulk param" );
	fn_bulk_memory_free	( sbulkparam2, "bulk param" );

	return hamilton;
}


/**
 * \brief	Projection of the stable bulk phase;
 */
void fn_main_project_bulk ( )
{
	/* introduce structure bodies; */
	stu_bulk_param sbulkparam1Type;
	stu_bulk_param *sbulkparam1 = &sbulkparam1Type;
	stu_bulk_param sbulkparam2Type;
	stu_bulk_param *sbulkparam2 = &sbulkparam2Type;
	stu_system_param ssysparamType;
	stu_system_param *ssysparam = &ssysparamType;

	/* true: save;	false: not save; */
	/* whether save densities; */
	sbulkparam1->sbndv->isSave	 = true; // 'fn_project_Fourier_GJP';
	sbulkparam2->sbndv->isSave	 = true; // 'fn_project_Fourier_GJP';
	/* variables for debug; */
	ssysparam->sGJPv->isTest	 = false; // 'fn_obt_system_gen_Jac_poly';
	sbulkparam1->sbndv->isTest	 = false; // 'fn_project_Fourier_GJP';
	sbulkparam2->sbndv->isTest	 = false; // 'fn_project_Fourier_GJP';


	/* input and save parameters; */
	fn_param_input ( sbulkparam1, sbulkparam2, ssysparam );

	double hamilton = 0.0;
	if ( strcmp( main_type, "proj_bulk1" ) == 0 )
	{
		/* obtain the stable bulk phases; */
		hamilton = fn_obt_stable_bulk_phase ( sbulkparam1 );

		/* prepare general Jacobi polynomials; */
		fn_obt_system_gen_Jac_poly ( ssysparam->sGJPv );

		/* project rhoCplx to the space of GJPs; */
		fn_project_Fourier_GJP ( sbulkparam1, ssysparam );

		/* release memory; */
		fn_bulk_memory_free ( sbulkparam1, "project bulk" );
		fn_bnd_memory_free  ( sbulkparam1->sbndv, "total" );
	}
	else if ( strcmp( main_type, "proj_bulk2" ) == 0 )
	{
		/* obtain the stable bulk phases; */
		hamilton = fn_obt_stable_bulk_phase ( sbulkparam2 );

		/* prepare general Jacobi polynomials; */
		fn_obt_system_gen_Jac_poly ( ssysparam->sGJPv );

		/* project rhoCplx to the space of GJPs; */
		fn_project_Fourier_GJP ( sbulkparam2, ssysparam );

		/* release memory; */
		fn_bulk_memory_free ( sbulkparam2, "project bulk" );
		fn_bnd_memory_free  ( sbulkparam2->sbndv, "total" );
	}
	else
	{
		if ( myrank == 0 )
		{
			printf("Error use 'fn_main_proj_bulk'.\n");
			printf("'main_type' should be 'proj_bulk1' or 'proj_bulk2'.\n");
		}
	}

	/* release memory; */
	fn_bulk_memory_free	  ( sbulkparam1, "bulk param" );
	fn_bulk_memory_free	  ( sbulkparam2, "bulk param" );
	fn_GJP_memory_free	  ( ssysparam->sGJPv );
}


/**
 * \brief	Transform rhoJCplx (in the space of Fourier and GJPs) to 
 *			the common space (Fourier and GJPs);
 */
void fn_main_common_bulk ( )
{
	/* introduce structure bodies; */
	stu_bulk_param sbulkparam1Type;
	stu_bulk_param *sbulkparam1 = &sbulkparam1Type;
	stu_bulk_param sbulkparam2Type;
	stu_bulk_param *sbulkparam2 = &sbulkparam2Type;
	stu_system_param ssysparamType;
	stu_system_param *ssysparam = &ssysparamType;

	/* true: save;	false: not save; */
	/* whether save densities; */
	sbulkparam1->sbndv->isSave	 = true; // 'fn_project_Fourier_GJP';
	sbulkparam2->sbndv->isSave	 = true; // 'fn_project_Fourier_GJP';
	sbulkparam1->srebndv->isSave = true; // 'fn_common_Fourier_GJP';
	sbulkparam2->srebndv->isSave = true; // 'fn_common_Fourier_GJP';
	/* variables for debug; */
	ssysparam->sGJPv->isTest	 = false; // 'fn_obt_system_gen_Jac_poly';
	sbulkparam1->sbndv->isTest	 = false; // 'fn_project_Fourier_GJP';
	sbulkparam2->sbndv->isTest	 = false; // 'fn_project_Fourier_GJP';
	sbulkparam1->srebndv->isTest = false; // 'fn_common_Fourier_GJP';
	sbulkparam2->srebndv->isTest = false; // 'fn_common_Fourier_GJP';

	
	/* input and save parameters; */
	fn_param_input ( sbulkparam1, sbulkparam2, ssysparam );

	/* obtain the stable bulk phases; */
	double hamilton1 = fn_obt_stable_bulk_phase ( sbulkparam1 );
	double hamilton2 = fn_obt_stable_bulk_phase ( sbulkparam2 );

	/* obtain the common 'rotateProjBoxMat'; */
	int status = fn_obt_commom_rotateProjBoxMat ( sbulkparam1, sbulkparam2, ssysparam );

	double hamilton = 0.0;
	if ( strcmp( main_type, "com_bulk1" ) == 0 )
	{
		/* prepare general Jacobi polynomials; */
		fn_obt_system_gen_Jac_poly ( ssysparam->sGJPv );

		/* project rhoCplx to the space of GJPs; */
		fn_project_Fourier_GJP ( sbulkparam1, ssysparam );

		if ( myrank == 0 )
			printf(" <========== Transform to the common space ==========> \n\n");
		mytimer_t timer;
		timer.reset();
		timer.start();

		/* obtain the projection plane about the common 'rotateProjBoxMat'; */
		fn_get_system_projection_plane ( ssysparam );

		/**
		 * re-represent by the common 'rotateProjBoxMat';
		 *	including 'rhoJCplx', 'd0bnd', ...;
		 */
		fn_total_re_represent_common ( sbulkparam1, ssysparam );

		MPI_Barrier ( MPI_COMM_WORLD );
		timer.pause();
		if ( myrank == 0 )
		{
			printf("\n\t ***** time cost of transformation to the common space: ");
			printf("%f seconds *****\n\n", timer.get_current_time());
		}

		/* release memory; */
		fn_bnd_memory_free  ( sbulkparam1->sbndv, "total" );
		fn_bnd_memory_free  ( sbulkparam1->srebndv, "total" );
	}
	else if ( strcmp( main_type, "com_bulk2" ) == 0 )
	{
		/* prepare general Jacobi polynomials; */
		fn_obt_system_gen_Jac_poly ( ssysparam->sGJPv );

		/* project rhoCplx to the space of GJPs; */
		fn_project_Fourier_GJP ( sbulkparam2, ssysparam );

		if ( myrank == 0 )
			printf(" <========== Transform to the common space ==========> \n\n");
		mytimer_t timer;
		timer.reset();
		timer.start();

		/* obtain the projection plane about the common 'rotateProjBoxMat'; */
		fn_get_system_projection_plane ( ssysparam );

		/**
		 * re-represent by the common 'rotateProjBoxMat';
		 *	including 'rhoJCplx', 'd0bnd', ...;
		 */
		fn_total_re_represent_common ( sbulkparam2, ssysparam );

		MPI_Barrier ( MPI_COMM_WORLD );
		timer.pause();
		if ( myrank == 0 )
		{
			printf("\n\t ***** time cost of transformation to the common space: ");
			printf("%f seconds *****\n\n", timer.get_current_time());
		}

		/* release memory; */
		fn_bnd_memory_free  ( sbulkparam2->sbndv, "total" );
		fn_bnd_memory_free  ( sbulkparam2->srebndv, "total" );
	}
	else
	{
		if ( myrank == 0 )
		{
			printf("Error use 'fn_main_proj_bulk'.\n");
			printf("'main_type' should be 'proj_bulk1' or 'proj_bulk2'.\n");
		}
	}

	/* release memory; */
	fn_bulk_memory_free		( sbulkparam1, "project bulk" );
	fn_bulk_memory_free		( sbulkparam2, "project bulk" );
	fn_bulk_memory_free		( sbulkparam1, "bulk param" );
	fn_bulk_memory_free		( sbulkparam2, "bulk param" );
	fn_GJP_memory_free		( ssysparam->sGJPv );
	fn_tvec_free<ptrdiff_t>	( ssysparam->NCpt );
	fn_tmat_free<double>	( ssysparam->rotateProjBoxMat );
	fn_tmat_free<int>		( ssysparam->coeffmat );
	fn_tmat_free<int>		( ssysparam->sfftv->indKspace );
	fn_tmat_free<double>	( ssysparam->sfftv->projPlane );
	free ( displs_sys );		free ( recvCount_sys );
}


/**
 * \brief	Calculation of the common projection matrix;
 */
void fn_main_common_rotateProjBoxMat ( )
{
	/* introduce structure bodies; */
	stu_bulk_param sbulkparam1Type;
	stu_bulk_param *sbulkparam1 = &sbulkparam1Type;
	stu_bulk_param sbulkparam2Type;
	stu_bulk_param *sbulkparam2 = &sbulkparam2Type;
	stu_system_param ssysparamType;
	stu_system_param *ssysparam = &ssysparamType;
	
	/* input and save parameters; */
	fn_param_input ( sbulkparam1, sbulkparam2, ssysparam );

	/* obtain the common 'rotateProjBoxMat'; */
	int status = fn_obt_commom_rotateProjBoxMat ( sbulkparam1, sbulkparam2, ssysparam );

	/* release memory; */
	fn_bulk_memory_free		( sbulkparam1, "bulk param" );
	fn_bulk_memory_free		( sbulkparam2, "bulk param" );
	fn_tmat_free<double>	( ssysparam->rotateProjBoxMat );
	fn_tmat_free<int>		( ssysparam->coeffmat );
}


/**
 * \brief	Display density according to the file about 'rhoJCplx';
 */
int fn_main_disp_system_density ( )
{
	/* introduce structure bodies; */
	stu_bulk_param sbulkparam1Type;
	stu_bulk_param *sbulkparam1 = &sbulkparam1Type;
	stu_bulk_param sbulkparam2Type;
	stu_bulk_param *sbulkparam2 = &sbulkparam2Type;
	stu_system_param ssysparamType;
	stu_system_param *ssysparam = &ssysparamType;

	/* input and save parameters; */
	fn_param_input ( sbulkparam1, sbulkparam2, ssysparam );

	FILE *fp;

	/* read 'd0JJ'; */
	char fname[FILELEN];
	sprintf(fname, "%s/d0GJP.dat", rsltDir); // file for saving parameters;
	fp = fopen(fname, "r");
	if ( fp == NULL )	// invalid file;
	{
		if ( myrank == 0 )
			printf("Unable to read file '%s'. No such file or directory.\n", fname);
		return 0;
	}
	else
	{
		ssysparam->sGJPv->d0JJ = fn_tmat_read_double ( fname );
	}
	fclose(fp);

	/* read 'd0bnd'; */
	sprintf(fname, "%s/rank%d/sys_d0bnd.dat", rsltDir, myrank); // file for saving parameters;
	fp = fopen(fname, "r");
	if ( fp == NULL )	// invalid file;
	{
		if ( myrank == 0 )
			printf("Unable to read file '%s'. No such file or directory.\n", fname);
		return 0;
	}
	else	
	{
		ssysparam->scbndv->d0bnd = fn_tvec_read_complex ( fname );
	}
	fclose(fp);

	/* read 'x'; */
	sprintf(fname, "%s/x.dat", rsltDir); // file for saving parameters;
	fp = fopen(fname, "r");
	if ( fp == NULL )	// invalid file;
	{
		if ( myrank == 0 )
			printf("Unable to read file '%s'. No such file or directory.\n", fname);
		return 0;
	}
	else
	{
		ssysparam->sGJPv->x = fn_tvec_read_double ( fname );
	}
	fclose(fp);

	/* read the common 'rotateProjBoxMat'; */
	char paraFile[FILELEN];
	sprintf(paraFile, "%s/parameter_opt.dat", rsltDir);
	if ( myrank == 0 )
		printf("parameter file: %s.\n\n", paraFile);
	fn_obt_com_projmat ( ssysparam, paraFile );	// read file;
	if ( myrank == 0 )
	{
		printf("\t ---> dimPhy = %d \t dimCpt = %d.\n", 
				ssysparam->dimRePhy, ssysparam->dimReCpt);
		printf("\t ---> The common projection matrix (read from %s) is \n", paraFile);
		fn_matrix_print ( ssysparam->rotateProjBoxMat );
		printf("\n");
	}

	/* parameters; */
	ssysparam->scbndv->xlen		= ssysparam->sGJPv->d0JJ.row;
	ssysparam->scbndv->nd		= ssysparam->sGJPv->d0JJ.col;

	/* generate 'projPlane'; */
	fn_get_system_projection_plane ( ssysparam );
	ssysparam->scbndv->cplxDofs = ssysparam->cplxReDofs;

	/* calculate the maximal iterator; */
	int line = 0;
	sprintf(fname, "%s/sys_energy_error.dat", rsltDir);
	if ( access(fname, F_OK) == 0 ) // exist;
	{
		int c, lc=0;
		fp = fopen(fname, "r");
		while ( (c = fgetc(fp)) != EOF )
		{
			if ( c == '\n' ) line ++;
			lc = c;
		}
		fclose(fp);
		if ( lc != '\n' ) line ++;
	}
	int iter_max = line - 1;
//	int iter_max = 0; // 0 for only plotting "sys_rhoJCplx-1.dat";

	/* read 'rhoJCplx'; */
	int	save_step	= ssysparam->save_step;
	int	iterator	= -1;
	tvec<fftw_complex> rhoJCplx;
	while ( iterator < iter_max )
	{
//		if ( (iterator % save_step) == 0 )
		if ( (iterator+1) % save_step == 0 || iterator == -1 || iterator == 0 )
		{
			sprintf(fname, "%s/rank%d/sys_rhoJCplx%d.dat", rsltDir, myrank, iterator);
			fp = fopen(fname, "r");
			if ( fp == NULL )	// invalid file;
			{
				if ( myrank == 0 )
					printf("Unable to read file '%s'. No such file or directory.\n", fname);
			}
			else
			{
				if ( myrank == 0 )
					printf("Read file '%s'.\n", fname);
				rhoJCplx = fn_tvec_read_complex ( fname );
				/* deal with NAN; */
				for ( int i = 0; i < rhoJCplx.len; i++ )
				{
					if ( isnan(rhoJCplx.val[i][0]) )
						rhoJCplx.val[i][0] = 0.0;
					if ( isnan(rhoJCplx.val[i][1]) )
						rhoJCplx.val[i][1] = 0.0;
				}
				fn_disp_system_density	( rhoJCplx, ssysparam, iterator );
			}
			fclose(fp);
		}
		iterator ++ ;
	}

	/* memory free about the interface system; */
	fn_bulk_memory_free			( sbulkparam1, "bulk param" );
	fn_bulk_memory_free			( sbulkparam2, "bulk param" );
	fn_tmat_free<double>		( ssysparam->sGJPv->d0JJ );
	fn_tvec_free<fftw_complex>	( ssysparam->scbndv->d0bnd );
	fn_tvec_free<fftw_complex>	( rhoJCplx );
	fn_tvec_free<ptrdiff_t>		( ssysparam->NCpt );
	fn_tvec_free<double>		( ssysparam->sGJPv->x );
	fn_tmat_free<double>		( ssysparam->rotateProjBoxMat );
	fn_tmat_free<double>		( ssysparam->sfftv->projPlane );
	fn_tmat_free<int>			( ssysparam->sfftv->indKspace );
	free ( displs_sys );		free ( recvCount_sys );
	return 1;
}


/**
 * \brief	Collect all data for interface system;
 */
int fn_main_collect_system_data ( )
{
	if ( nprocs > 1 )
	{
		printf("only one process, please.\n\n");
		return 0;
	}

	/* introduce structure bodies; */
	stu_bulk_param sbulkparam1Type;
	stu_bulk_param *sbulkparam1 = &sbulkparam1Type;
	stu_bulk_param sbulkparam2Type;
	stu_bulk_param *sbulkparam2 = &sbulkparam2Type;
	stu_system_param ssysparamType;
	stu_system_param *ssysparam = &ssysparamType;

	/* input and save parameters; */
	fn_param_input ( sbulkparam1, sbulkparam2, ssysparam );

	/* parameters; */
	int		nd	   = ssysparam->sGJPv->nd;
	int		xlen   = ssysparam->sGJPv->xlen;
	FILE *fpr;
	FILE *fpw;
	int	nprocsData = nprocs;

	/* the directory for total data; */
	char rsltDirAll[FILELEN0+10];
	sprintf(rsltDirAll, "%s/all", rsltDir);
	mkdir(rsltDirAll, 0755);

	/* read the number of process 'nprocsData' when users obtained the data; */
    char    buffer[STRLEN]; // Note: max number of char for each line!
    int     val;
	int		status = 1;		// 1: success; 0: failure;
	char paraFile[FILELEN];
	sprintf(paraFile, "%s/parameter_opt.dat", rsltDir);

    fpr = fopen(paraFile, "r");
    if ( fpr == NULL )	// invalid file;
	{
		if ( myrank == 0 )
			printf("Unable to read file '%s'. No such file or directory.\n", paraFile);
		return 0;
	}

    while ( status == 1 )
	{
		int     ibuff;
    
        val = fscanf(fpr, "%s", buffer);
        if (val==EOF) break;
        if (val!=1) { status = 0; break; }
        if (buffer[0]=='#' || buffer[0]=='%' || buffer[0]=='|' || buffer[0]=='/')
		{
            if (fscanf(fpr, "%*[^\n]")) {/* skip rest of line and do nothing */ };
            continue;
        }
    
        /* match keyword and scan for value; */
		if ( strcmp ( buffer, "nprocs" ) == 0 )
		{
			val = fscanf(fpr, "%s", buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			val = fscanf(fpr, "%d", &(ibuff));
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			nprocsData = ibuff;
			if (fscanf(fpr, "%*[^\n]")) {/* skip rest of line and do nothing */ };
		}
		else if ( strcmp ( buffer, "com_projmat_size" ) == 0 )
		{
			val = fscanf(fpr, "%s", buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			val = fscanf(fpr, "%d", &(ibuff));
			val = fscanf(fpr, "%d", &(ibuff));
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			ssysparam->dimReCpt = ibuff;
			if (fscanf(fpr, "%*[^\n]")) {/* skip rest of line and do nothing */ };
		}
	}
    fclose(fpr);
	if ( myrank == 0 )
		printf("The number of process when users obtained the data is %d.\n", nprocsData);
	int cplxReDofs = 1;	
	for ( int i = 0; i < ssysparam->dimReCpt; i++ )
		cplxReDofs *= ssysparam->Four_num;

	/* write 'd0bnd' which contains complete information; */
	tvec<fftw_complex> d0bnd = fn_tvec_init<fftw_complex> ( xlen*cplxReDofs );
	d0bnd.row = xlen;		d0bnd.col = cplxReDofs;
	int size, ind;
    int *displs    = (int *) malloc ( sizeof(int) * nprocsData ); // offset of each process;
    int *recvCount = (int *) malloc ( sizeof(int) * nprocsData ); // offset of each process;
	displs[0] = 0;
	double dataReal, dataCplx;
	char fname[FILELEN];

	/* read 'd0bnd'; */
	for ( int n = 0; n < nprocsData; n++ )
	{
		sprintf(fname, "%s/rank%d/sys_d0bnd.dat", rsltDir, n);
		fpr = fopen(fname, "r");
		val = fscanf(fpr, "%d", &(size));		// read 'row';
		val = fscanf(fpr, "%d", &(size));		// read 'col';
		recvCount[n] = size;
		if ( n > 0 )
			displs[n] = displs[n-1] + recvCount[n-1];
		val = fscanf(fpr, "%d", &(size));		// read 'len';
		for ( int s = 0; s < xlen; s++ )
		{
			for ( int i = 0; i < recvCount[n]; i++ )
			{
				val = fscanf(fpr, "%lf", &(dataReal));
				val = fscanf(fpr, "%lf", &(dataCplx));
				ind = s*cplxReDofs + i + displs[n];
				d0bnd.val[ind][0] = dataReal;
				d0bnd.val[ind][1] = dataCplx;
			}
		}
		fclose(fpr);
		if ( myrank == 0 )
			printf("Read sys_d0bnd.dat --> progress: %d/%d.\n", n, nprocsData);
	}

	/* write 'd0bnd'; */
	sprintf(fname, "%s/all/sys_d0bnd.dat", rsltDir); // file for saving parameters;
	fpw = fopen(fname, "w");
	fprintf(fpw, "%d\t%d\t%d\n", d0bnd.row, d0bnd.col, d0bnd.len);
	for ( int i = 0; i < d0bnd.len; i++ )
		fprintf(fpw, "%+.15E\t%+.15E\n", d0bnd.val[i][0], d0bnd.val[i][1]);
	fclose(fpw);
	if ( myrank == 0 )
		printf("Write sys_d0bnd.dat.\n\n");
	fn_tvec_free<fftw_complex>	( d0bnd );

	/* calculate the maximal iterator; */
	int line = 0;
	sprintf(fname, "%s/sys_energy_error.dat", rsltDir);
	if ( access(fname, F_OK) == 0 ) // exist;
	{
		int c, lc=0;
		fpr = fopen(fname, "r");
		while ( (c = fgetc(fpr)) != EOF )
		{
			if ( c == '\n' ) line ++;
			lc = c;
		}
		fclose(fpr);
		if ( lc != '\n' ) line ++;
	}
	int iter_max = line - 1;
//	int iter_max = 0; // 0 for only plotting "sys_rhoJCplx-1.dat";

	/* write 'rhoJCplx' which contains complete information; */
	tvec<fftw_complex> rhoJCplx = fn_tvec_init<fftw_complex> ( nd*cplxReDofs );
	rhoJCplx.row = nd;		rhoJCplx.col = cplxReDofs;
	int	save_step	= ssysparam->save_step;
	int	iterator	= -1;
	if ( myrank == 0 )
		printf("save_step = %d, iter_max = %d\n\n", save_step, iter_max);
	while ( iterator < iter_max )
	{
		if ( (iterator+1) % save_step == 0 || iterator == -1 || iterator == 0 )
		{
			/* read 'rhoJCplx'; */
			for ( int n = 0; n < nprocsData; n++ )
			{
				sprintf(fname, "%s/rank%d/sys_rhoJCplx%d.dat", rsltDir, n, iterator);
				fpr = fopen(fname, "r");
				/* old code only saves 'len', please note that; */
				val = fscanf(fpr, "%d", &(size));	// read 'row';
				val = fscanf(fpr, "%d", &(size));	// read 'col';
				val = fscanf(fpr, "%d", &(size));	// read 'len';
				for ( int s = 0; s < nd; s++ )
				{
					for ( int i = 0; i < recvCount[n]; i++ )
					{
						val = fscanf(fpr, "%lf", &(dataReal));
						val = fscanf(fpr, "%lf", &(dataCplx));
						ind = s*cplxReDofs + i + displs[n];
						rhoJCplx.val[ind][0] = dataReal;
						rhoJCplx.val[ind][1] = dataCplx;
					}
				}
				fclose(fpr);
				if ( myrank == 0 )
					printf("read sys_rhoJCplx%d.dat --> progress: %d/%d.\n", 
						iterator, n, nprocsData);
			}

			/* write 'rhoJCplx'; */
			sprintf(fname, "%s/all/sys_rhoJCplx%d.dat", rsltDir, iterator);
			fpw = fopen(fname, "w");
			fprintf(fpw, "%d\t%d\t%d\n", rhoJCplx.row, rhoJCplx.col, rhoJCplx.len);
			for ( int i = 0; i < rhoJCplx.len; i++ )
				fprintf(fpw, "%+.15E\t%+.15E\n", rhoJCplx.val[i][0], rhoJCplx.val[i][1]);
			fclose(fpw);
			if ( myrank == 0 )
				printf("Write sys_rhoJCplx%d.dat.\n\n", iterator);
		}
		iterator ++ ;
	}
	fn_tvec_free<fftw_complex>	( rhoJCplx );

	/* memory free about the interface system; */
	fn_bulk_memory_free			( sbulkparam1, "bulk param" );
	fn_bulk_memory_free			( sbulkparam2, "bulk param" );
	free ( displs );		free ( recvCount );
	return 1;
}


/**
 * \brief	recover density profile by the screened 'rhoJCplx';
 */
int fn_main_recover_system_density ( )
{
	/* introduce structure bodies; */
	stu_bulk_param sbulkparam1Type;
	stu_bulk_param *sbulkparam1 = &sbulkparam1Type;
	stu_bulk_param sbulkparam2Type;
	stu_bulk_param *sbulkparam2 = &sbulkparam2Type;
	stu_system_param ssysparamType;
	stu_system_param *ssysparam = &ssysparamType;

	/* input and save parameters; */
	fn_param_input ( sbulkparam1, sbulkparam2, ssysparam );

	FILE *fp;

	/* read 'd0JJ'; */
	char fname[FILELEN];
	sprintf(fname, "%s/d0GJP.dat", rsltDir); // file for saving parameters;
	fp = fopen(fname, "r");
	if ( fp == NULL )	// invalid file;
	{
		if ( myrank == 0 )
			printf("Unable to read file '%s'. No such file or directory.\n", fname);
		return 0;
	}
	else
	{
		ssysparam->sGJPv->d0JJ = fn_tmat_read_double ( fname );
	}
	fclose(fp);

	/* read 'd0bnd'; */
	sprintf(fname, "%s/rank%d/sys_d0bnd.dat", rsltDir, myrank); // file for saving parameters;
	fp = fopen(fname, "r");
	if ( fp == NULL )	// invalid file;
	{
		if ( myrank == 0 )
			printf("Unable to read file '%s'. No such file or directory.\n", fname);
		return 0;
	}
	else	
	{
		ssysparam->scbndv->d0bnd = fn_tvec_read_complex ( fname );
	}
	fclose(fp);

	/* read 'x'; */
	sprintf(fname, "%s/x.dat", rsltDir); // file for saving parameters;
	fp = fopen(fname, "r");
	if ( fp == NULL )	// invalid file;
	{
		if ( myrank == 0 )
			printf("Unable to read file '%s'. No such file or directory.\n", fname);
		return 0;
	}
	else
	{
		ssysparam->sGJPv->x = fn_tvec_read_double ( fname );
	}
	fclose(fp);

	/* read 'w'; */
	sprintf(fname, "%s/w.dat", rsltDir); // file for saving parameters;
	fp = fopen(fname, "r");
	if ( fp == NULL )	// invalid file;
	{
		if ( myrank == 0 )
			printf("Unable to read file '%s'. No such file or directory.\n", fname);
		return 0;
	}
	else
	{
		ssysparam->sGJPv->w = fn_tvec_read_double ( fname );
	}
	fclose(fp);

	/* read the common 'rotateProjBoxMat'; */
	char paraFile[FILELEN];
	sprintf(paraFile, "%s/parameter_opt.dat", rsltDir);
	if ( myrank == 0 )
		printf("parameter file: %s.\n\n", paraFile);
	fn_obt_com_projmat ( ssysparam, paraFile );	// read file;
	if ( myrank == 0 )
	{
		printf("\t ---> dimPhy = %d \t dimCpt = %d.\n", 
				ssysparam->dimRePhy, ssysparam->dimReCpt);
		printf("\t ---> The common projection matrix (read from %s) is \n", paraFile);
		fn_matrix_print ( ssysparam->rotateProjBoxMat );
		printf("\n");
	}

	/* parameters; */
	ssysparam->scbndv->xlen		= ssysparam->sGJPv->d0JJ.row;
	ssysparam->scbndv->nd		= ssysparam->sGJPv->d0JJ.col;

	/* generate 'projPlane'; */
	fn_get_system_projection_plane ( ssysparam );
	ssysparam->scbndv->cplxDofs = ssysparam->cplxReDofs;

	/* calculate 'Gsquare'; */
	fn_obt_system_Gsquare ( ssysparam );

	/* obtain the inner product matrices; */
	ssysparam->sGJPv->innSMatd0JJ = fn_innSMat_gen_Jac_poly(
			ssysparam->sGJPv, ssysparam->sGJPv->d0JJ, 0, 1e-8);

	/* read 'rhoJCplx'; */
	int	save_step	= ssysparam->save_step;
	int	iterator	= -1;
	tvec<fftw_complex> rhoJCplx;
	sprintf(fname, "%s/rank%d/sys_rhoJCplx%d.dat", rsltDir, myrank, iterator);
	fp = fopen(fname, "r");
	if ( fp == NULL )	// invalid file;
	{
		if ( myrank == 0 )
			printf("Unable to read file '%s'. No such file or directory.\n", fname);
	}
	else
	{
		if ( myrank == 0 )
			printf("Read file '%s'.\n", fname);
		rhoJCplx = fn_tvec_read_complex ( fname );
		/* deal with NAN; */
		for ( int i = 0; i < rhoJCplx.len; i++ )
		{
			if ( isnan(rhoJCplx.val[i][0]) )
				rhoJCplx.val[i][0] = 0.0;
			if ( isnan(rhoJCplx.val[i][1]) )
				rhoJCplx.val[i][1] = 0.0;
		}
		/* recover the density profile by the screened spectra; */
		int status = fn_recover_rhoJCplx ( ssysparam, rhoJCplx );
		if ( status == 0 )
		{
			/* save 'rhoJCplx'; */
			sprintf(fname, "%s/rank%d/sys_rhoJCplx%d.dat", rsltDir, myrank, -2);
			fn_tvec_save_complex ( rhoJCplx, fname );
			/* save density calculated by 'rhoJCplx'; */
			fn_disp_system_density	( rhoJCplx, ssysparam, -2 );
		}
		/* calculate the error between density profiles 
		 * before and after the recover operation; */
		fn_recover_error ( ssysparam, -1, -2 );
	}
	fclose(fp);

	/* memory free about the interface system; */
	fn_bulk_memory_free			( sbulkparam1, "bulk param" );
	fn_bulk_memory_free			( sbulkparam2, "bulk param" );
	fn_tvec_free<fftw_complex>	( rhoJCplx );
	fn_tvec_free<fftw_complex>	( ssysparam->scbndv->d0bnd );
	fn_tvec_free<ptrdiff_t>		( ssysparam->NCpt );
	fn_tmat_free<double>		( ssysparam->sGJPv->d0JJ );
	fn_tvec_free<double>		( ssysparam->sGJPv->x );
	fn_tvec_free<double>		( ssysparam->sGJPv->w );
	fn_tCCSmat_free<double>		( ssysparam->sGJPv->innSMatd0JJ );
	fn_tmat_free<double>		( ssysparam->rotateProjBoxMat );
	fn_tmat_free<int>			( ssysparam->coeffmat );
	fn_tmat_free<double>		( ssysparam->sfftv->projPlane );
	fn_tmat_free<int>			( ssysparam->sfftv->indKspace );
	fn_tvec_free<double>		( ssysparam->sfftv->Gsquare );
	free ( displs_sys );		free ( recvCount_sys );
	return 1;
}


/**
 * \brief	Generate parameter files in batches;
 */
void fn_param_batch ( )
{
	/* the path of parameter file; */
	if ( myrank == 0 )
		printf(" <========== Generate parameter files in batches ==========> \n\n");
	char inFile[FILELEN], outFile[FILELEN];
	sprintf(inFile, "./para/%s/input%s.dat", paraDir, para_flag);
	if ( myrank == 0 )
		printf("inFile: %s.\n", inFile);

	for ( int theta = 0; theta <= 2; theta++ )
	{
		sprintf(outFile, "./para/%s/input%d%s.dat", paraDir, theta, "copy");
		if ( myrank == 0 )
			printf("outFile: %s.\n", outFile);

		/* what needs to be modified; */
		int mLine0 = 54;
		char mContent0[1024];
		sprintf(mContent0, "bulk1_rotate_angle  = 0.0  -45.0  -%d", theta);
		int mLine1 = 91;
		char mContent1[1024];
		sprintf(mContent1, "bulk2_rotate_angle  = 0.0  -45.0  %d", theta);
		int mLine2 = 144;
		char mContent2[1024];
		sprintf(mContent2, "x_range_type        = cubePlane%d", 1);

		/* read data from "inFile" and write data into "outFile"; */
		ifstream in;
		in.open(inFile);
		char line[1024] = {'\0'};
		int i = 0;
		string tempStr;
		while ( in.getline(line, sizeof(line)) )
		{
			i++;
			if ( i == mLine0 )
			{
				tempStr += charToStr(mContent0);
			}
			else if ( i == mLine1 )
			{
				tempStr += charToStr(mContent1);
			}
			else if ( i == mLine2 )
			{
				tempStr += charToStr(mContent2);
			}
			else
			{
				tempStr += charToStr(line);
			}
			tempStr += '\n';
		}
		in.close();

		/* write data into "outFile"; */
		ofstream out;
		out.open(outFile);
		out.flush();
		out << tempStr;
		out.close();
	}
}


/**
 * \brief	char to str;
 */
string charToStr ( char *contentChar )
{
	string tempStr;
	for ( int i=0; contentChar[i]!='\0'; i++ )
	{
		tempStr += contentChar[i];
	}
	return tempStr;
}
