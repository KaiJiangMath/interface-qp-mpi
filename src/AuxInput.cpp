/*! \file	AuxInput.cpp
 *
 *  \brief	Read input parameters
 *
 */

#include "Head.h"
#include "Data.h"
#include "Mytimer.h"
#include "DataOperators.h"
#include "functs.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/


/**
 * \brief	The main code of input parameters;
 *
 * \param	sbulkparam1		Input parameters about the left bulk phase;
 * \param	sbulkparam2		Input parameters about the right bulk phase;
 * \param	ssysparam		Input parameters about the interface system;
 *
 */
void fn_param_input (	stu_bulk_param			*sbulkparam1,
						stu_bulk_param			*sbulkparam2,
						stu_system_param		*ssysparam )
{
	if ( myrank == 0 )
		printf(" <========== Input parameters ==========> \n\n");
	mytimer_t timer;
	timer.reset();
	timer.start();

	char paraFile[FILELEN];
	sprintf(paraFile, "./para/%s/input%s.dat", paraDir, para_flag);

	/* parameters about the left/right bulk phase; */
	sbulkparam1->sflag = 1;
	sbulkparam2->sflag = 2;
	fn_stu_bulk_param_input( paraFile, sbulkparam1 );
	fn_stu_bulk_param_input( paraFile, sbulkparam2 );

	/* parameters about the interface system; */
	fn_stu_system_param_input( paraFile, ssysparam );

	/**
	 * the parameters about GJPs;
	 * correct for LB, LP,...
	 */
	ssysparam->sGJPv->polyDegree = 4 * ssysparam->scale_num;
	ssysparam->sGJPv->alpha		 = 2.0 * ssysparam->scale_num;
	ssysparam->sGJPv->beta		 = 2.0 * ssysparam->scale_num;
	ssysparam->sGJPv->nd		 = ssysparam->GJP_degree + 1;
	ssysparam->sGJPv->xlen		 = ssysparam->LGL_num	+ 1;
	ssysparam->sGJPv->x_range	 = ssysparam->x_range;

	/* directory for saving results; */
	sprintf(rsltDir, "./result/%s_%s_%s_%s", model_type,
			sbulkparam1->phase, sbulkparam2->phase, para_flag);

	/* copy parameter from stu_system_param to stu_bulk_param; */
	fn_param_copy_system_to_bulk( ssysparam, sbulkparam1 );
	fn_param_copy_system_to_bulk( ssysparam, sbulkparam2 );

	/* save parameters;
	 *	if the computational box will be optimized, so we save it again;
	 */
	mkdir("./result/", 0755);
	mkdir(rsltDir, 0755);		// create fold;
	/* the directory for each process; */
	char rsltDirRank[FILELEN0+10];
	sprintf(rsltDirRank, "%s/rank%d", rsltDir, myrank);
	mkdir(rsltDirRank, 0755);
	if ( myrank == 0 )
		fn_param_save ( sbulkparam1, sbulkparam2, ssysparam, "init" );

	/* print parameters; */
	if ( myrank == 0 )
	{
		fn_stu_bulk_param_print( sbulkparam1 );
		fn_stu_bulk_param_print( sbulkparam2 );
		if ( strcmp( main_type, "stable_system" ) == 0 )
			fn_stu_system_param_print( ssysparam );
	}

	/* adjust model parameters; */
	fn_adjust_bulk_model_param	 ( sbulkparam1 );
	fn_adjust_bulk_model_param	 ( sbulkparam2 );
	fn_adjust_system_model_param ( ssysparam );

	timer.pause();
	if ( myrank == 0 )
		printf("\n\t ***** time cost of loading parameters: %f seconds *****\n\n", 
				timer.get_current_time());
}


/**
 * \brief	Default input parameters about bulk phases;
 *
 * \param	sbulkparam		Input parameters about bulk phases;
 *
 */
void fn_stu_bulk_param_input_init (stu_bulk_param	*sbulkparam)
{
	strcpy(sbulkparam->phase, "HEX");
	sbulkparam->print_level			= 0;		// 0 is not displayed;

	/* the rotation order and rotation angles around x,y,z-axis; */
	strcpy(sbulkparam->motion_order, "rotate_transl");
	strcpy(sbulkparam->rotate_order, "iiz");
	fn_tvec_setZero<double> ( sbulkparam->rotate_angle );

	/* the translation vector; */
	fn_tvec_setZero<double> ( sbulkparam->transl_var );

	/* the calculation and plotting parameters; */
	sbulkparam->Four_num				= 4;
	sbulkparam->enlarge					= 2;
	sbulkparam->plot_num				= 16;
	sbulkparam->skip_half				= 0;

	/* the iteration parameters; */
	sbulkparam->tol						= 1e-6;
	sbulkparam->step_size				= 0.1;
	sbulkparam->iter_max				= 10000;
	sbulkparam->print_step				= 100;
	sbulkparam->save_step				= 100;
	strcpy(sbulkparam->save_type, "yyn");

	/* the parameters in the optimization of computational box; */
	sbulkparam->opt_tol					= 1e-6;
	sbulkparam->opt_iter_max			= 10;
	sbulkparam->opt_save_step			= 10;
	sbulkparam->box_bbType				= 1;
	sbulkparam->box_iter_max			= 50;
	sbulkparam->box_step_size			= 0.1;
	sbulkparam->box_tol					= 1e-6;

	/* the parameters for truncation; */
	sbulkparam->trunc_tol				= 1e-10;
}


/**
 * \brief	Default input parameters about bulk phases;
 *			(natural parameters of bulk phases);
 *
 * \param	sbulkparam		Input natural parameters about bulk phases;
 *
 */
void fn_bulk_nature_param_input_init (stu_bulk_param	*sbulkparam)
{
	/* parameters depending on the bulk phase; */
	sbulkparam->dimPhy				= 2;
	sbulkparam->dimCpt				= 2;
	sbulkparam->nfold				= 6;		// only work for quasicrystals;
}


/**
 * \brief	Default input parameters about the interface;
 *
 * \param	ssysparam		Input parameters about the interface;
 *
 */
void fn_stu_system_param_input_init (stu_system_param	*ssysparam)
{
	strcpy ( model_type, "LB" );
	strcpy ( ssysparam->com_projmat_way, "calculate" );
	strcpy ( ssysparam->iter_method, "SIS" );
	ssysparam->print_level			= 0;	// 0 is not displayed;

	/* model parameters; */
	ssysparam->scale_num			= 1;
	ssysparam->scale_val[0]			= 1.0;
	ssysparam->model_xi				= 1.0;
	ssysparam->model_tau			= 0.0;
	ssysparam->model_gamma			= 1.0;
	ssysparam->model_kappa			= 1.0;

	/* spatial discrete parameters; */
	ssysparam->GJP_degree			= 256;
	ssysparam->LGL_num				= 512;
	ssysparam->Four_num				= 20;
	strcpy ( ssysparam->x_range_type, "direct" );
	ssysparam->x_range				= 1.0;
	ssysparam->smooth				= 0.1;
	ssysparam->initDist1			= 0.0;
	ssysparam->initDist2			= 0.0;
	ssysparam->searchReg			= 40;
	ssysparam->adjustReg			= 40;

	/* iteration parameters; */
	ssysparam->massFlag				= 1;
	ssysparam->tol					= 1e-6;
	ssysparam->tolham				= 1e-16;
	ssysparam->step_size			= 0.1;
	ssysparam->step_min				= 0.01;
	ssysparam->step_max				= 1.0;
	ssysparam->iter_load			= 0;
	ssysparam->iter_max				= 10000;
	ssysparam->print_step			= 100;
	ssysparam->save_step			= 1000;
	strcpy(ssysparam->save_type, "yyn");

	/* some iteration parameters for Newton-PCG method; */
	ssysparam->newton_tol			= 1e-4;
	ssysparam->newton_step_size		= 0.5;
	ssysparam->pcg_type				= 0;
	ssysparam->pcg_iter_max			= 1000;
	ssysparam->pcg_print_step		= 100;
	ssysparam->pcg_delta_coeff		= 0.1;
	ssysparam->pcg_mu_para			= 2.0;

	/* plotting parameters; */
	ssysparam->y_start				= 0.0;
	ssysparam->z_start				= 0.0;
	ssysparam->y_range				= 10.0;
	ssysparam->z_range				= 0.0;
	ssysparam->y_num				= 128;
	ssysparam->z_num				= 0;
	ssysparam->skip_half			= 0;

	/* parameters for representing error; */
	ssysparam->err_y_start			= 0.0;
	ssysparam->err_z_start			= 0.0;
	ssysparam->err_y_range			= 10.0;
	ssysparam->err_z_range			= 10.0;
	ssysparam->err_y_num			= 8;
	ssysparam->err_z_num			= 8;

	/* parameter to recover; */
	ssysparam->recoverTOL			= 1e-4;
}


/**
 * \brief	Read input parameters from disk file;
 *
 * \param	fname			File name for input file;
 * \param	sbulkparam		Input parameters about bulk phases;
 *
 */
void fn_stu_bulk_param_input (	const char		*fname,
								stu_bulk_param	*sbulkparam )
{
	int		sflag = sbulkparam->sflag;
    char    buffer[STRLEN]; // Note: max number of char for each line!
    int     val;
	int		status = 1;		// 1: success; 0: failure;
    FILE    *fp;

	/* memory allocation of 'rotate_angle' and 'transl_var'; */
	sbulkparam->rotate_angle = fn_tvec_init<double> ( 3 ); // x,y,z;
	sbulkparam->transl_var	 = fn_tvec_init<double> ( 3 ); // x,y,z;

    fn_stu_bulk_param_input_init(sbulkparam);	// set default input parameters
   
    fp = fopen(fname,"r");
    if ( fp == NULL )	// invalid file;
	{
		if ( myrank == 0 )
		{
			printf("Error using 'fn_stu_bulk_param_input'\n");
			printf("Unable to read file '%s'. No such file or directory.\n", fname);
		}
		return;			// return the default values;
	}
    
    while ( status == 1 )
	{
        int     ibuff;
        double  dbuff;
        char    sbuff[STRLEN];
		char	keyword[WORDLEN];
    
        val = fscanf(fp, "%s", buffer);
        if (val==EOF) break;
        if (val!=1) { status = 0; break; }
        if (buffer[0]=='#' || buffer[0]=='%' || buffer[0]=='|' || buffer[0]=='/')
		{
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
            continue;
        }
    
        /* match keyword and scan for value; */
		sprintf(keyword, "bulk%d_phase", sflag);
		if ( strcmp(buffer, keyword) == 0 )
		{
			val = fscanf(fp, "%s", buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			val = fscanf(fp, "%s", sbuff);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			strcpy(sbulkparam->phase, sbuff);
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
		}
		sprintf(keyword, "bulk%d_print_level", sflag);
		if ( strcmp(buffer, keyword) == 0 )
		{
			val = fscanf(fp, "%s", buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			val = fscanf(fp, "%d", &(ibuff));
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			sbulkparam->print_level = ibuff;
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
		}
		/* motion order: rotate, then translate; translate, then rotate; */
		sprintf(keyword, "bulk%d_motion_order", sflag);
		if ( strcmp(buffer, keyword) == 0 )
		{
			val = fscanf(fp, "%s", buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			val = fscanf(fp, "%s", sbuff);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			strcpy(sbulkparam->motion_order, sbuff);
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
		}
		/* the rotation order and rotation angles around x,y,z-axis; */
		sprintf(keyword, "bulk%d_rotate_order", sflag);
		if ( strcmp(buffer, keyword) == 0 )
		{
			val = fscanf(fp, "%s", buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			val = fscanf(fp, "%s", sbuff);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			strcpy(sbulkparam->rotate_order, sbuff);
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
		}
		sprintf(keyword, "bulk%d_rotate_angle", sflag);
		if ( strcmp(buffer, keyword) == 0 )
		{
			val = fscanf(fp, "%s", buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			for ( int i = 0; i < sbulkparam->rotate_angle.len; i++ )
			{
				val = fscanf(fp, "%lf", &(dbuff));
				if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
				sbulkparam->rotate_angle.val[i] = dbuff;
			}
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
		}
		/* the translation vector; */
		sprintf(keyword, "bulk%d_transl_var", sflag);
		if ( strcmp(buffer, keyword) == 0 )
		{
			val = fscanf(fp, "%s", buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			for ( int i = 0; i < sbulkparam->transl_var.len; i++ )
			{
				val = fscanf(fp, "%lf", &(dbuff));
				if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
				sbulkparam->transl_var.val[i] = dbuff;
			}
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
		}
		/* the calculation and plotting parameters; */
		sprintf(keyword, "bulk%d_Four_num", sflag);
		if ( strcmp(buffer, keyword) == 0 )
		{
			val = fscanf(fp, "%s", buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			val = fscanf(fp, "%d", &(ibuff));
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			sbulkparam->Four_num = ibuff;
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
		}
		sprintf(keyword, "bulk%d_enlarge", sflag);
		if ( strcmp(buffer, keyword) == 0 )
		{
			val = fscanf(fp, "%s", buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			val = fscanf(fp, "%lf", &(dbuff));
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			sbulkparam->enlarge = dbuff;
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
		}
		sprintf(keyword, "bulk%d_plot_num", sflag);
		if ( strcmp(buffer, keyword) == 0 )
		{
			val = fscanf(fp, "%s", buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			val = fscanf(fp, "%d", &(ibuff));
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			sbulkparam->plot_num = ibuff;
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
		}
		sprintf(keyword, "bulk%d_is_half", sflag);
		if ( strcmp(buffer, keyword) == 0 )
		{
			val = fscanf(fp, "%s", buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			val = fscanf(fp, "%d", &(ibuff));
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			sbulkparam->skip_half = ibuff;
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
		}
		sprintf(keyword, "bulk%d_tol", sflag);
		if ( strcmp(buffer, keyword) == 0 )
		{
			val = fscanf(fp, "%s", buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			val = fscanf(fp, "%lf", &(dbuff));
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			sbulkparam->tol = dbuff;				// load coefficient;
			val = fscanf(fp, "%lf", &(dbuff));
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			sbulkparam->tol *= pow(10, dbuff);		// load index;
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
		}
		sprintf(keyword, "bulk%d_step_size", sflag);
		if ( strcmp(buffer, keyword) == 0 )
		{
			val = fscanf(fp, "%s", buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			val = fscanf(fp, "%lf", &(dbuff));
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			sbulkparam->step_size = dbuff;
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
		}
		sprintf(keyword, "bulk%d_iter_max", sflag);
		if ( strcmp(buffer, keyword) == 0 )
		{
			val = fscanf(fp, "%s", buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			val = fscanf(fp, "%d", &(ibuff));
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			sbulkparam->iter_max = ibuff;
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
		}
		sprintf(keyword, "bulk%d_print_step", sflag);
		if ( strcmp(buffer, keyword) == 0 )
		{
			val = fscanf(fp, "%s", buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			val = fscanf(fp, "%d", &(ibuff));
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			sbulkparam->print_step = ibuff;
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
		}
		sprintf(keyword, "bulk%d_save_step", sflag);
		if ( strcmp(buffer, keyword) == 0 )
		{
			val = fscanf(fp, "%s", buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			val = fscanf(fp, "%d", &(ibuff));
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			sbulkparam->save_step = ibuff;
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
		}
		sprintf(keyword, "bulk%d_save_type", sflag);
		if ( strcmp(buffer, keyword) == 0 )
		{
			val = fscanf(fp, "%s", buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			val = fscanf(fp, "%s", sbuff);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			strcpy(sbulkparam->save_type, sbuff);
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
		}
		/* parameters to optimize computational box; */
		sprintf(keyword, "bulk%d_opt_tol", sflag);
		if ( strcmp(buffer, keyword) == 0 )
		{
			val = fscanf(fp, "%s", buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			val = fscanf(fp, "%lf", &(dbuff));
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			sbulkparam->opt_tol = dbuff;				// load coefficient;
			val = fscanf(fp, "%lf", &(dbuff));
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			sbulkparam->opt_tol *= pow(10, dbuff);		// load index;
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
		}
		sprintf(keyword, "bulk%d_opt_iter_max", sflag);
		if ( strcmp(buffer, keyword) == 0 )
		{
			val = fscanf(fp, "%s", buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			val = fscanf(fp, "%d", &(ibuff));
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			sbulkparam->opt_iter_max = ibuff;
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
		}
		sprintf(keyword, "bulk%d_opt_save_step", sflag);
		if ( strcmp(buffer, keyword) == 0 )
		{
			val = fscanf(fp, "%s", buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			val = fscanf(fp, "%d", &(ibuff));
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			sbulkparam->opt_save_step = ibuff;
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
		}
		sprintf(keyword, "bulk%d_box_bbType", sflag);
		if ( strcmp(buffer, keyword) == 0 )
		{
			val = fscanf(fp, "%s", buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			val = fscanf(fp, "%d", &(ibuff));
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			sbulkparam->box_bbType = ibuff;
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
		}
		sprintf(keyword, "bulk%d_box_iter_max", sflag);
		if ( strcmp(buffer, keyword) == 0 )
		{
			val = fscanf(fp, "%s", buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			val = fscanf(fp, "%d", &(ibuff));
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			sbulkparam->box_iter_max = ibuff;
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
		}
		sprintf(keyword, "bulk%d_box_step_size", sflag);
		if ( strcmp(buffer, keyword) == 0 )
		{
			val = fscanf(fp, "%s", buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			val = fscanf(fp, "%lf", &(dbuff));
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			sbulkparam->box_step_size = dbuff;
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
		}
		sprintf(keyword, "bulk%d_box_tol", sflag);
		if ( strcmp(buffer, keyword) == 0 )
		{
			val = fscanf(fp, "%s", buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			val = fscanf(fp, "%lf", &(dbuff));
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			sbulkparam->box_tol = dbuff;				// load coefficient;
			val = fscanf(fp, "%lf", &(dbuff));
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			sbulkparam->box_tol *= pow(10, dbuff);		// load index;
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
		}
		sprintf(keyword, "bulk%d_trunc_tol", sflag);
		if ( strcmp(buffer, keyword) == 0 )
		{
			val = fscanf(fp, "%s", buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			val = fscanf(fp, "%lf", &(dbuff));
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			sbulkparam->trunc_tol = dbuff;				// load coefficient;
			val = fscanf(fp, "%lf", &(dbuff));
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			sbulkparam->trunc_tol *= pow(10, dbuff);		// load index;
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
		}
    }
    fclose(fp);

	/* read the natural parameters of bulk phases; */
	fn_bulk_nature_param_input (sbulkparam);

	/* discrete grids; */
	sbulkparam->NCpt  = fn_tvec_init<int> ( sbulkparam->dimCpt );
	sbulkparam->pNCpt = fn_tvec_init<ptrdiff_t> ( sbulkparam->dimCpt );
	sbulkparam->realDofs = 1;
	sbulkparam->cplxDofs = 1;
	for ( int i = 0; i < sbulkparam->dimCpt; i++ )
	{
		sbulkparam->NCpt.val[i]  = sbulkparam->Four_num;
		sbulkparam->pNCpt.val[i] = sbulkparam->Four_num;
		sbulkparam->realDofs *= sbulkparam->Four_num;
		sbulkparam->cplxDofs *= sbulkparam->Four_num;
	}

	/* direction box in real space; */
	sbulkparam->dirBox = fn_tmat_init<double> ( sbulkparam->dimCpt, sbulkparam->dimCpt );
	matrix_inverse ( sbulkparam->rcpBox.val, sbulkparam->dirBox.val, sbulkparam->dimCpt	);
	for ( int i = 0; i < sbulkparam->dimCpt; i++ )
		for ( int j = 0; j < sbulkparam->dimCpt; j++ )
			sbulkparam->dirBox.val[i][j] *= 2*PI;

	/* projection matrix; */
	fn_obt_projection_matrix ( sbulkparam );

	/* rotation matrix; */
	fn_obt_rotation_matrix ( fname, sbulkparam );

	/* translation matrix; */
	fn_obt_translation_vector ( sbulkparam );

	/* calculate rotateMat' * projMat * rcpBox; */
	sbulkparam->rotateProjBoxMat = fn_tmat_init<double> ( sbulkparam->dimPhy, sbulkparam->dimCpt );
	fn_obt_rotateProjBox_matrix ( sbulkparam );

	/* calculate translVec' * projMat * rcpBox; */
	sbulkparam->translProjBoxVec = fn_tvec_init<double> ( sbulkparam->dimCpt );
	fn_obt_translProjBox_matrix ( sbulkparam );

#if DEBUG_MODE > 1
	if ( myrank == 0 )
		printf("### DEBUG: Reading input (bulk %d) status = %d\n", sflag, status);
#endif
}


/**
 * \brief	Read input parameters from disk file; 
 *			(natural parameters of bulk phases);
 *
 * \param	sbulkparam		Input natural parameters about bulk phases;
 *
 */
void fn_bulk_nature_param_input (stu_bulk_param	*sbulkparam)
{
	char    buffer[STRLEN]; // Note: max number of char for each line!
    int     val;
	int		status = 1;		// 1: success; 0: failure;
    FILE    *fp;

    fn_bulk_nature_param_input_init(sbulkparam);	// set default input parameters
   
	char fname[FILELEN];
	sprintf(fname, "./initData/%s", sbulkparam->phase); // file path;
    fp = fopen(fname,"r");
    if ( fp == NULL )	// invalid file;
	{
		if ( myrank == 0 )
		{
			printf("Error using 'bulk_nature_param_input'\n");
			printf("Unable to read file '%s'. No such file or directory.\n", fname);
		}
		return;			// return the default values;
	}
    
    while ( status == 1 )
	{
        int     ibuff;
        double  dbuff;
        char    sbuff[STRLEN];
		char	keyword[WORDLEN];
    
        val = fscanf(fp, "%s", buffer);
        if (val==EOF) break;
        if (val!=1) { status = 0; break; }
        if (buffer[0]=='#' || buffer[0]=='%' || buffer[0]=='|' || buffer[0]=='/')
		{
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
            continue;
        }
    
        /* match keyword and scan for value; */
		/* parameters depending on the bulk phase; */
		if ( strcmp(buffer, "dimPhy") == 0 )
		{
			val = fscanf(fp, "%s", buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			val = fscanf(fp, "%d", &(ibuff));
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			sbulkparam->dimPhy = ibuff;
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
		}
		else if ( strcmp(buffer, "dimCpt") == 0 )
		{
			val = fscanf(fp, "%s", buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			val = fscanf(fp, "%d", &(ibuff));
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			sbulkparam->dimCpt = ibuff;
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
		}
		else if ( strcmp(buffer, "nfold") == 0 )
		{
			val = fscanf(fp, "%s", buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			val = fscanf(fp, "%d", &(ibuff));
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			sbulkparam->nfold = ibuff;
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
		}
		else if ( strcmp(buffer, "boxType") == 0 )
		{
			val = fscanf(fp, "%s", buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			val = fscanf(fp, "%s", sbuff);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			strcpy(sbulkparam->boxType, sbuff);
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
		}
		else if ( strcmp(buffer, "boxVal") == 0 ) // reciprocal box (computational box);
		{
			val = fscanf(fp, "%s", buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			/* initialization; */
			sbulkparam->rcpBox = fn_tmat_init<double> ( sbulkparam->dimCpt, sbulkparam->dimCpt );
			fn_tmat_setZero<double> ( sbulkparam->rcpBox );
			/* cubic box or hexagonal box; */
			if ( strcmp(sbulkparam->boxType, "cube") == 0 )			// same elements on the diagonal;
			{
				val = fscanf(fp, "%lf", &(dbuff));
				if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
				for ( int i = 0; i < sbulkparam->rcpBox.row; i++ )
				{
					sbulkparam->rcpBox.val[i][i] = dbuff;
				}
			}
			else if ( strcmp(sbulkparam->boxType, "cuboid") == 0 )	// only diagonal elements;
			{
				for ( int i = 0; i < sbulkparam->rcpBox.row; i++ )
				{
					val = fscanf(fp, "%lf", &(dbuff));
					if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
					sbulkparam->rcpBox.val[i][i] = dbuff;
				}
			}
			else if ( strcmp(sbulkparam->boxType, "hex") == 0 )		// all elements;
			{
				for ( int i = 0; i < sbulkparam->rcpBox.row; i++ )
				{
					for ( int j = 0; j < sbulkparam->rcpBox.col; j++ )
					{
						val = fscanf(fp, "%lf", &(dbuff));
						if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
						sbulkparam->rcpBox.val[i][j] = dbuff;
					}
				}
			}
			else
			{
				if ( myrank == 0 )
					printf("Error setting 'boxType' or 'boxVal' in %s.\n", fname);
				return;			// return the default values;
			}
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
	    }
		else if ( strcmp(buffer, "initNum") == 0 )
		{
			val = fscanf(fp, "%s", buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			val = fscanf(fp, "%d", &(ibuff));
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			sbulkparam->initNum = ibuff;
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
		}
		else if ( strcmp(buffer, "initVal") == 0 )
		{
			val = fscanf(fp, "%s", buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			/* initialization; */
			sbulkparam->initIndex = fn_tmat_init<int> ( sbulkparam->initNum, sbulkparam->dimCpt );
			sbulkparam->initCoeff = fn_tvec_init<fftw_complex> ( sbulkparam->initNum );
			fn_tvec_setZero_complex ( sbulkparam->initCoeff );
			/* read file to obtain primary spectra; */
			for ( int i = 0; i < sbulkparam->initNum; i++ )
			{
				for ( int j = 0; j < sbulkparam->dimCpt; j++ )
				{
					val = fscanf(fp, "%d", &(ibuff));
					if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
					sbulkparam->initIndex.val[i][j] = ibuff;
				}
				for ( int j = 0; j < 2; j++ )
				{
					val = fscanf(fp, "%lf", &(dbuff));
					if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
					sbulkparam->initCoeff.val[i][j] = dbuff;
				}
			}
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
		}
	}
    fclose(fp);

#if DEBUG_MODE > 1
	if ( myrank == 0 )
		printf("### DEBUG: Reading input (bulk %d) status = %d\n", sbulkparam->sflag, status);
#endif
}


/**
 * \brief	Set the projection matrix;
 *
 * \param	sbulkparam		Return projection matrix;
 *
 */
void fn_obt_projection_matrix (	stu_bulk_param		*sbulkparam )
{
	// memory allocation of projection matrix;
	sbulkparam->projMat = fn_tmat_init<double> ( sbulkparam->dimPhy, sbulkparam->dimCpt );
	fn_tmat_setZero<double> ( sbulkparam->projMat );

	// projection matrix for crystals or quasicrystals;
	if ( sbulkparam->projMat.row == sbulkparam->projMat.col )
	{
		for ( int i = 0; i < sbulkparam->projMat.row; i ++ )
			sbulkparam->projMat.val[i][i] = 1.0;
	}
	else if ( sbulkparam->projMat.row < sbulkparam->projMat.col )
	{
		if ( sbulkparam->projMat.row == 2 )
		{
			for ( int j = 0; j < sbulkparam->projMat.col; j++ )
			{
				sbulkparam->projMat.val[0][j] = cos(2*PI*j/sbulkparam->nfold);
				sbulkparam->projMat.val[1][j] = sin(2*PI*j/sbulkparam->nfold);
			}
		}
		else if ( sbulkparam->projMat.row == 3 && sbulkparam->projMat.col == 6 )
		{
			double qq = 2*cos(PI/5);
			sbulkparam->projMat.val[0][0] = 1.0;
			sbulkparam->projMat.val[0][1] = 0.5*qq;
			sbulkparam->projMat.val[0][2] = 0.5*qq;
			sbulkparam->projMat.val[0][3] = 0.5;
			sbulkparam->projMat.val[0][4] = 0.0;
			sbulkparam->projMat.val[0][5] = 0.0;
			sbulkparam->projMat.val[1][0] = 0.0;
			sbulkparam->projMat.val[1][1] = 0.5;
			sbulkparam->projMat.val[1][2] = -0.5;
			sbulkparam->projMat.val[1][3] = 0.5*(qq-1.0);
			sbulkparam->projMat.val[1][4] = 1.0;
			sbulkparam->projMat.val[1][5] = 0.0;
			sbulkparam->projMat.val[2][0] = 0.0;
			sbulkparam->projMat.val[2][1] = 0.5*(qq-1.0);
			sbulkparam->projMat.val[2][2] = 0.5*(1.0-qq);
			sbulkparam->projMat.val[2][3] = -0.5*qq;
			sbulkparam->projMat.val[2][4] = 0.0;
			sbulkparam->projMat.val[2][5] = 1.0;
		}
		else
		{
			if ( myrank == 0 )
			{
				printf("Error using 'fn_obt_projection_matrix'\n");
				printf("Only consider 'dimPhy = 2' for 'dimPhy < dimCpt'.\n");
			}
		}
	}
	else
	{
		if ( myrank == 0 )
		{
			printf("Error using 'fn_obt_projection_matrix'\n");
			printf("'dimPhy' cannot be smaller than 'dimCpt'.\n");
		}
	}
}


/**
 * \brief	Generate the rotation matrix;
 *
 * \param	sbulkparam		Input the rotation angles;
 *
 */
void fn_obt_rotation_matrix (	const char		*fname,
								stu_bulk_param	*sbulkparam )
{
	/* initialization; */
	double theta, cosVal, sinVal;
	sbulkparam->rotateMat = fn_tmat_init<double> ( sbulkparam->dimPhy, sbulkparam->dimPhy );
	
	/* generate rotation matrix; */
	if ( strcmp(sbulkparam->rotate_order, "direct") == 0 )
	{
		int		sflag = sbulkparam->sflag;
		char    buffer[STRLEN]; // Note: max number of char for each line!
		int     val;
		int		status = 1;		// 1: success; 0: failure;
		FILE    *fp;

		fp = fopen(fname, "r");

		if ( fp == NULL )	// invalid file;
		{
			if ( myrank == 0 )
			{
				printf("Error using 'fn_obt_rotation_matrix'\n");
				printf("Unable to read file '%s'. No such file or directory.\n", fname);
			}
			return;			// return the default values;
		}
	 
		while ( status == 1 )
		{
			int     ibuff;
			double  dbuff;
			char    sbuff[STRLEN];
			char	keyword[WORDLEN];
		
			val = fscanf(fp, "%s", buffer);
			if (val==EOF) break;
			if (val!=1) { status = 0; break; }
			if (buffer[0]=='#' || buffer[0]=='%' || buffer[0]=='|' || buffer[0]=='/')
			{
				if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
				continue;
			}
		
			/* match keyword and scan for value; */
			sprintf(keyword, "bulk%d_rotateMat", sflag);
			if ( strcmp(buffer, keyword) == 0 )
			{
				val = fscanf(fp, "%s", buffer);
				if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
				/* directly read rotation matrix; */	
				for ( int i = 0; i < sbulkparam->rotateMat.row; i++ )
				{
					for ( int j = 0; j < sbulkparam->rotateMat.col; j++ )
					{
						val = fscanf(fp, "%lf", &(dbuff));
						if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
						sbulkparam->rotateMat.val[i][j] = dbuff;
					}
				}
				if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
			}
		}
		fclose(fp);

#if DEBUG_MODE > 1
		if ( myrank == 0 )
			printf("### DEBUG: Reading input (bulk %d) status = %d\n", 
					sbulkparam->sflag, status);
#endif
	}
	else
	{
		if ( sbulkparam->dimPhy == 1 )
			sbulkparam->rotateMat.val[0][0] = 1.0;
		else if ( sbulkparam->dimPhy == 2 )
		{
			theta = sbulkparam->rotate_angle.val[0] / 180.0 * PI;	// angle to radians;
			cosVal = cos(theta), sinVal = sin(theta);
			sbulkparam->rotateMat.val[0][0] = cosVal;
			sbulkparam->rotateMat.val[0][1] = -sinVal;
			sbulkparam->rotateMat.val[1][0] = sinVal;
			sbulkparam->rotateMat.val[1][1] = cosVal;
		}
		else if ( sbulkparam->dimPhy == 3 )
		{
			/* set temporary variables and initialization; */
			tmat<double> rotateMatXYZ  = fn_tmat_init<double> ( 3, 3 ); // temporary matrix;
			tmat<double> rotateMatRslt = fn_tmat_init<double> ( 3, 3 ); // product of two matrices;
			fn_tmat_setIdentity<double> ( sbulkparam->rotateMat );

			/* calculate rotation matrix by the rule 'rotate_order'; */
			for ( int i = 0; i < 3; i++ ) // Rx, Ry, Rz;
			{
				fn_tmat_setIdentity<double> ( rotateMatXYZ );
				theta = sbulkparam->rotate_angle.val[i] / 180.0 * PI;	// angle to radian;
				cosVal = cos(theta), sinVal = sin(theta);
				if ( sbulkparam->rotate_order[i] == 'x' || sbulkparam->rotate_order[i] == 'X' )
				{
					rotateMatXYZ.val[1][1] = cosVal;
					rotateMatXYZ.val[1][2] = -sinVal;
					rotateMatXYZ.val[2][1] = sinVal;
					rotateMatXYZ.val[2][2] = cosVal;
				}
				else if ( sbulkparam->rotate_order[i] == 'y' || sbulkparam->rotate_order[i] == 'Y' )
				{
					rotateMatXYZ.val[0][0] = cosVal;
					rotateMatXYZ.val[0][2] = sinVal;
					rotateMatXYZ.val[2][0] = -sinVal;
					rotateMatXYZ.val[2][2] = cosVal;
				}
				else if ( sbulkparam->rotate_order[i] == 'z' || sbulkparam->rotate_order[i] == 'Z' )
				{
					rotateMatXYZ.val[0][0] = cosVal;
					rotateMatXYZ.val[0][1] = -sinVal;
					rotateMatXYZ.val[1][0] = sinVal;
					rotateMatXYZ.val[1][1] = cosVal;
				}
				/**
				 * rotateMatRslt = rotateMat * rotateMatXYZ;
				 * rotateMat = rotateMatRslt to update;
				 */
				fn_tmat_multiply<double> ( sbulkparam->rotateMat, rotateMatXYZ, rotateMatRslt );
				fn_tmat_copy<double>	 ( sbulkparam->rotateMat, rotateMatRslt );
			}

			/* memory free; */
			fn_tmat_free<double> ( rotateMatXYZ );
			fn_tmat_free<double> ( rotateMatRslt );
		}
	}
}


/**
 * \brief	Set the translation vector;
 *
 * \param	stu_bulk_param		Input translation in each direction;
 *
 */
void fn_obt_translation_vector (stu_bulk_param	*sbulkparam)
{
 	sbulkparam->translVec = fn_tvec_init<double> ( sbulkparam->dimPhy );
	for ( int i = 0; i < sbulkparam->dimPhy; i++ ) // dimPhy <= 3;
	{
		sbulkparam->translVec.val[i] = sbulkparam->transl_var.val[i];
	}
}


/**
 * \brief	Calculate the products of rotateMat, projMat, and rcpBox;
 *
 * \param	sbulkparam		Input rotateMat, projMat, rcpBox;
 *
 */
void fn_obt_rotateProjBox_matrix ( stu_bulk_param	*sbulkparam )
{
	/* load parameters; */
	int dimPhy = sbulkparam->dimPhy;
	int dimCpt = sbulkparam->dimCpt;

	/* initialization; */
	tmat<double> rpTmp = fn_tmat_init<double> ( dimPhy, dimCpt );

	/**
	 * calculate rotateMat'*projMat;
	 *		rotateMat:	size is [dimPhy, dimPhy];
	 *		projMat:	size is [dimPhy, dimCpt];
	 * superscript ' means transpose;
	 * the product is denoted as rpTmp;
	 */
	for (int i = 0; i < dimPhy; i++)
	{
		for (int j = 0; j < dimCpt; j++)
		{
			rpTmp.val[i][j] = 0.0;
			for (int jj = 0; jj < dimPhy; jj++)
				rpTmp.val[i][j] += sbulkparam->rotateMat.val[jj][i] * 
								sbulkparam->projMat.val[jj][j];
		}
	}

	/**
	 * calculate rotateMat'*projMat*rcpBox;
	 *		i.e. rpTmp * rcpBox;
	 *			rpTmp:	size is [dimPhy, dimCpt];
	 *			rcpBox:	size is [dimCpt, dimCpt];
	 */
	for (int i = 0; i < dimPhy; i++)
	{
		for (int j = 0; j < dimCpt; j++)
		{
			sbulkparam->rotateProjBoxMat.val[i][j] = 0.0;
			for (int jj = 0; jj < dimCpt; jj++)
				sbulkparam->rotateProjBoxMat.val[i][j] += 
					rpTmp.val[i][jj] * sbulkparam->rcpBox.val[jj][j];
		}
	}

	/* memory free; */
	fn_tmat_free<double> ( rpTmp );
}


/**
 * \brief	Calculate the product of translVec, projMat, rcpBox;
 *
 * \param	sbulkparam		Input translVec, projMat, rcpBox;
 *
 */
void fn_obt_translProjBox_matrix (stu_bulk_param	*sbulkparam)
{
	if ( strcmp(sbulkparam->motion_order, "transl_rotate") == 0 )
	{
		/* initialization; */
		tvec<double> tpTmp = fn_tvec_init<double> ( sbulkparam->dimCpt );

		/**
		 * calculate -translVec' * projMat;
		 *			translVec:	size is [dimPhy, 1];
		 *			projMat:	size is [dimPhy, dimCpt];
		 * ' means transpose;
		 * the product result is denoted as tpTmp;
		 */
		for (int i = 0; i < sbulkparam->dimCpt; i++)
		{
			tpTmp.val[i] = 0.0;
			for (int jj = 0; jj < sbulkparam->dimPhy; jj++)
			{
				tpTmp.val[i] -= sbulkparam->translVec.val[jj] * 
					sbulkparam->projMat.val[jj][i];
			}
		}

		/**
		 * calculate -translVec' * projMat * rcpBox;
		 *			i.e. tpTmp * rcpBox;
		 *				tpTmp:	size is [1, dimCpt];
		 *				rcpBox:	size is [dimCpt, dimCpt];
		 */
		for (int i = 0; i < sbulkparam->dimCpt; i++)
		{
			sbulkparam->translProjBoxVec.val[i] = 0.0;
			for (int jj = 0; jj < sbulkparam->dimCpt; jj++)
			{
				sbulkparam->translProjBoxVec.val[i] += 
					tpTmp.val[jj] * sbulkparam->rcpBox.val[jj][i];
			}
		}
		fn_tvec_free<double> ( tpTmp );
	}
	else
	{
		/**
		 * calculate -translVec' * rotateProjBoxMat;
		 *			rotateProjBoxMat = rotateMat'*projMat*rcpBox;
		 *			translVec:	size is [dimPhy, 1];
		 *			rotateProjBoxMat:	size is [dimPhy, dimCpt];
		 * ' means transpose;
		 */
		for (int i = 0; i < sbulkparam->dimCpt; i++)
		{
			sbulkparam->translProjBoxVec.val[i] = 0.0;
			for (int jj = 0; jj < sbulkparam->dimPhy; jj++)
			{
				sbulkparam->translProjBoxVec.val[i] -= 
					sbulkparam->translVec.val[jj] * 
					sbulkparam->rotateProjBoxMat.val[jj][i];
			}
		}
	}
}


/**
 * \brief	Read input parameters from disk file;
 *
 * \param	fname			File name for input file;
 * \param	ssysparam		Input parameters about the interface;
 *
 */
void fn_stu_system_param_input (	const char		*fname,
								stu_system_param	*ssysparam)
{
    char    buffer[STRLEN]; // Note: max number of char for each line!
    int     val;
	int		status = 1;		// 1: success; 0: failure;
    FILE    *fp;

    fn_stu_system_param_input_init(ssysparam);	// set default input parameters
   
    fp = fopen(fname,"r");
    if ( fp == NULL )	// invalid file;
	{
		if ( myrank == 0 )
		{
			printf("Error using 'fn_stu_system_param_input'\n");
			printf("Unable to read file '%s'. No such file or directory.\n", fname);
		}
		return;			// return the default values;
	}
    
    while ( status == 1 )
	{
        int     ibuff;
        double  dbuff;
        char    sbuff[STRLEN];
		char	keyword[WORDLEN];
    
        val = fscanf(fp, "%s", buffer);
        if (val==EOF) break;
        if (val!=1) { status = 0; break; }
        if (buffer[0]=='#' || buffer[0]=='%' || buffer[0]=='|' || buffer[0]=='/')
		{
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
            continue;
        }
    
        /* match keyword and scan for value; */
		if ( strcmp(buffer, "model_type") == 0 )
		{
			val = fscanf(fp, "%s", buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			val = fscanf(fp, "%s", sbuff);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			strcpy(model_type, sbuff);
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
		}
		else if ( strcmp(buffer, "com_projmat_way") == 0 )
		{
			val = fscanf(fp, "%s", buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			val = fscanf(fp, "%s", sbuff);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			strcpy(ssysparam->com_projmat_way, sbuff);
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
		}
		else if ( strcmp(buffer, "iter_method") == 0 )
		{
			val = fscanf(fp, "%s", buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			val = fscanf(fp, "%s", sbuff);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			strcpy(ssysparam->iter_method, sbuff);
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
		}
		else if ( strcmp(buffer, "print_level") == 0 )
		{
			val = fscanf(fp, "%s", buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			val = fscanf(fp, "%d", &(ibuff));
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			ssysparam->print_level = ibuff;
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
		}
		/* model parameters; */
		else if ( strcmp(buffer, "scale_num") == 0 )
		{
			val = fscanf(fp, "%s", buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			val = fscanf(fp, "%d", &(ibuff));
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			ssysparam->scale_num = ibuff;
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
		}
		else if ( strcmp(buffer, "scale_val") == 0 )
		{
			val = fscanf(fp, "%s", buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			for ( int i = 0; i < ssysparam->scale_num; i++ )
			{
				val = fscanf(fp, "%lf", &(dbuff));
				if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
				ssysparam->scale_val[i] = dbuff;
			}
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
		}
		else if ( strcmp(buffer, "model_xi") == 0 )
		{
			val = fscanf(fp, "%s", buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			val = fscanf(fp, "%lf", &(dbuff));
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			ssysparam->model_xi = dbuff;
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
		}
		else if ( strcmp(buffer, "model_tau") == 0 )
		{
			val = fscanf(fp, "%s", buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			val = fscanf(fp, "%lf", &(dbuff));
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			ssysparam->model_tau = dbuff;
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
		}
		else if ( strcmp(buffer, "model_gamma") == 0 )
		{
			val = fscanf(fp, "%s", buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			val = fscanf(fp, "%lf", &(dbuff));
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			ssysparam->model_gamma = dbuff;
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
		}
		else if ( strcmp(buffer, "model_kappa") == 0 )
		{
			val = fscanf(fp, "%s", buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			val = fscanf(fp, "%lf", &(dbuff));
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			ssysparam->model_kappa = dbuff;
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
		}
		/* spatial discrete parameters; */
		else if ( strcmp(buffer, "GJP_degree") == 0 )
		{
			val = fscanf(fp, "%s", buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			val = fscanf(fp, "%d", &(ibuff));
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			ssysparam->GJP_degree = ibuff;
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
		}
		else if ( strcmp(buffer, "LGL_num") == 0 )
		{
			val = fscanf(fp, "%s", buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			val = fscanf(fp, "%d", &(ibuff));
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			ssysparam->LGL_num = ibuff;
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
		}
		else if ( strcmp(buffer, "Four_num") == 0 )
		{
			val = fscanf(fp, "%s", buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			val = fscanf(fp, "%d", &(ibuff));
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			ssysparam->Four_num = ibuff;
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
		}
		else if ( strcmp(buffer, "x_range_type") == 0 )
		{
			val = fscanf(fp, "%s", buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			val = fscanf(fp, "%s", sbuff);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			strcpy(ssysparam->x_range_type, sbuff);
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
		}
		else if ( strcmp(buffer, "x_range") == 0 )
		{
			val = fscanf(fp, "%s", buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			val = fscanf(fp, "%lf", &(dbuff));
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			ssysparam->x_range = dbuff;
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
		}
		else if ( strcmp(buffer, "smooth") == 0 )
		{
			val = fscanf(fp, "%s", buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			val = fscanf(fp, "%lf", &(dbuff));
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			ssysparam->smooth = dbuff;
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
		}
		else if ( strcmp(buffer, "initDist1") == 0 )
		{
			val = fscanf(fp, "%s", buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			val = fscanf(fp, "%lf", &(dbuff));
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			ssysparam->initDist1 = dbuff;
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
		}
		else if ( strcmp(buffer, "initDist2") == 0 )
		{
			val = fscanf(fp, "%s", buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			val = fscanf(fp, "%lf", &(dbuff));
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			ssysparam->initDist2 = dbuff;
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
		}
		else if ( strcmp(buffer, "searchReg") == 0 )
		{
			val = fscanf(fp, "%s", buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			val = fscanf(fp, "%d", &(ibuff));
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			ssysparam->searchReg = ibuff;
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
		}
		else if ( strcmp(buffer, "adjustReg") == 0 )
		{
			val = fscanf(fp, "%s", buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			val = fscanf(fp, "%d", &(ibuff));
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			ssysparam->adjustReg = ibuff;
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
		}
		/* iteration parameters; */
		else if ( strcmp(buffer, "massFlag") == 0 )
		{
			val = fscanf(fp, "%s", buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			val = fscanf(fp, "%d", &(ibuff));
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			ssysparam->massFlag = ibuff;
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
		}
		else if ( strcmp(buffer, "tol") == 0 )
		{
			val = fscanf(fp, "%s", buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			val = fscanf(fp, "%lf", &(dbuff));
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			ssysparam->tol = dbuff;				// load coefficient;
			val = fscanf(fp, "%lf", &(dbuff));
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			ssysparam->tol *= pow(10, dbuff);		// load index;
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
		}
		else if ( strcmp(buffer, "tolham") == 0 )
		{
			val = fscanf(fp, "%s", buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			val = fscanf(fp, "%lf", &(dbuff));
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			ssysparam->tolham = dbuff;				// load coefficient;
			val = fscanf(fp, "%lf", &(dbuff));
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			ssysparam->tolham *= pow(10, dbuff);		// load index;
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
		}
		else if ( strcmp(buffer, "step_size") == 0 )
		{
			val = fscanf(fp, "%s", buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			val = fscanf(fp, "%lf", &(dbuff));
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			ssysparam->step_size = dbuff;
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
		}
		else if ( strcmp(buffer, "step_min") == 0 )
		{
			val = fscanf(fp, "%s", buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			val = fscanf(fp, "%lf", &(dbuff));
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			ssysparam->step_min = dbuff;
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
		}
		else if ( strcmp(buffer, "step_max") == 0 )
		{
			val = fscanf(fp, "%s", buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			val = fscanf(fp, "%lf", &(dbuff));
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			ssysparam->step_max = dbuff;
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
		}
		else if ( strcmp(buffer, "iter_load") == 0 )
		{
			val = fscanf(fp, "%s", buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			val = fscanf(fp, "%d", &(ibuff));
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			ssysparam->iter_load = ibuff;
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
		}
		else if ( strcmp(buffer, "iter_max") == 0 )
		{
			val = fscanf(fp, "%s", buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			val = fscanf(fp, "%d", &(ibuff));
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			ssysparam->iter_max = ibuff;
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
		}
		else if ( strcmp(buffer, "print_step") == 0 )
		{
			val = fscanf(fp, "%s", buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			val = fscanf(fp, "%d", &(ibuff));
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			ssysparam->print_step = ibuff;
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
		}
		else if ( strcmp(buffer, "save_step") == 0 )
		{
			val = fscanf(fp, "%s", buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			val = fscanf(fp, "%d", &(ibuff));
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			ssysparam->save_step = ibuff;
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
		}
		else if ( strcmp(buffer, "save_type") == 0 )
		{
			val = fscanf(fp, "%s", buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			val = fscanf(fp, "%s", sbuff);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			strcpy(ssysparam->save_type, sbuff);
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
		}
		/* some iteration parameters for Newton-PCG method; */
		else if ( strcmp(buffer, "newton_tol") == 0 )
		{
			val = fscanf(fp, "%s", buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			val = fscanf(fp, "%lf", &(dbuff));
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			ssysparam->newton_tol = dbuff;				// load coefficient;
			val = fscanf(fp, "%lf", &(dbuff));
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			ssysparam->newton_tol *= pow(10, dbuff);	// load index;
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
		}
		else if ( strcmp(buffer, "newton_step_size") == 0 )
		{
			val = fscanf(fp, "%s", buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			val = fscanf(fp, "%lf", &(dbuff));
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			ssysparam->newton_step_size = dbuff;
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
		}
		else if ( strcmp(buffer, "pcg_type") == 0 )
		{
			val = fscanf(fp, "%s", buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			val = fscanf(fp, "%d", &(ibuff));
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			ssysparam->pcg_type = ibuff;
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
		}
		else if ( strcmp(buffer, "pcg_iter_max") == 0 )
		{
			val = fscanf(fp, "%s", buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			val = fscanf(fp, "%d", &(ibuff));
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			ssysparam->pcg_iter_max = ibuff;
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
		}
		else if ( strcmp(buffer, "pcg_print_step") == 0 )
		{
			val = fscanf(fp, "%s", buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			val = fscanf(fp, "%d", &(ibuff));
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			ssysparam->pcg_print_step = ibuff;
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
		}
		else if ( strcmp(buffer, "pcg_delta_coeff") == 0 )
		{
			val = fscanf(fp, "%s", buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			val = fscanf(fp, "%lf", &(dbuff));
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			ssysparam->pcg_delta_coeff = dbuff;
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
		}
		else if ( strcmp(buffer, "pcg_mu_para") == 0 )
		{
			val = fscanf(fp, "%s", buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			val = fscanf(fp, "%lf", &(dbuff));
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			ssysparam->pcg_mu_para = dbuff;
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
		}
		/* plotting parameters; */
		else if ( strcmp(buffer, "y_start") == 0 )
		{
			val = fscanf(fp, "%s", buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			val = fscanf(fp, "%lf", &(dbuff));
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			ssysparam->y_start = dbuff;
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
		}
		else if ( strcmp(buffer, "z_start") == 0 )
		{
			val = fscanf(fp, "%s", buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			val = fscanf(fp, "%lf", &(dbuff));
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			ssysparam->z_start = dbuff;
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
		}
		else if ( strcmp(buffer, "y_range") == 0 )
		{
			val = fscanf(fp, "%s", buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			val = fscanf(fp, "%lf", &(dbuff));
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			ssysparam->y_range = dbuff;
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
		}
		else if ( strcmp(buffer, "z_range") == 0 )
		{
			val = fscanf(fp, "%s", buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			val = fscanf(fp, "%lf", &(dbuff));
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			ssysparam->z_range = dbuff;
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
		}
		else if ( strcmp(buffer, "y_num") == 0 )
		{
			val = fscanf(fp, "%s", buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			val = fscanf(fp, "%d", &(ibuff));
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			ssysparam->y_num = ibuff;
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
		}
		else if ( strcmp(buffer, "z_num") == 0 )
		{
			val = fscanf(fp, "%s", buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			val = fscanf(fp, "%d", &(ibuff));
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			ssysparam->z_num = ibuff;
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
		}
		else if ( strcmp(buffer, "skip_half") == 0 )
		{
			val = fscanf(fp, "%s", buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			val = fscanf(fp, "%d", &(ibuff));
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			ssysparam->skip_half = ibuff;
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
		}
		/* parameters for representing error; */
		else if ( strcmp(buffer, "err_y_start") == 0 )
		{
			val = fscanf(fp, "%s", buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			val = fscanf(fp, "%lf", &(dbuff));
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			ssysparam->err_y_start = dbuff;
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
		}
		else if ( strcmp(buffer, "err_z_start") == 0 )
		{
			val = fscanf(fp, "%s", buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			val = fscanf(fp, "%lf", &(dbuff));
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			ssysparam->err_z_start = dbuff;
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
		}
		else if ( strcmp(buffer, "err_y_range") == 0 )
		{
			val = fscanf(fp, "%s", buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			val = fscanf(fp, "%lf", &(dbuff));
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			ssysparam->err_y_range = dbuff;
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
		}
		else if ( strcmp(buffer, "err_z_range") == 0 )
		{
			val = fscanf(fp, "%s", buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			val = fscanf(fp, "%lf", &(dbuff));
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			ssysparam->err_z_range = dbuff;
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
		}
		else if ( strcmp(buffer, "err_y_num") == 0 )
		{
			val = fscanf(fp, "%s", buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			val = fscanf(fp, "%d", &(ibuff));
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			ssysparam->err_y_num = ibuff;
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
		}
		else if ( strcmp(buffer, "err_z_num") == 0 )
		{
			val = fscanf(fp, "%s", buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			val = fscanf(fp, "%d", &(ibuff));
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			ssysparam->err_z_num = ibuff;
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
		}
		else if ( strcmp(buffer, "recoverTOL") == 0 )
		{
			val = fscanf(fp, "%s", buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			val = fscanf(fp, "%lf", &(dbuff));
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			ssysparam->recoverTOL = dbuff;				// load coefficient;
			val = fscanf(fp, "%lf", &(dbuff));
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			ssysparam->recoverTOL *= pow(10, dbuff);		// load index;
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
		}
    }
    
    fclose(fp);

#if DEBUG_MODE > 1
	if ( myrank == 0 )
		printf("### DEBUG: Reading input (interface) status = %d\n", status);
#endif
}


/**
 * \brief	Copy data from stu_system_param to stu_bulk_param;
 *
 * \param	ssysparam			Input data;
 * \param	sbulkparam			Save data;
 */
void fn_param_copy_system_to_bulk (stu_system_param	*ssysparam,
								stu_bulk_param		*sbulkparam)
{
	/* copy model parameters; */
	sbulkparam->scale_num	= ssysparam->scale_num;
	for ( int i = 0; i < ssysparam->scale_num; i++ )
		sbulkparam->scale_val[i]		= ssysparam->scale_val[i];
	sbulkparam->model_xi		= ssysparam->model_xi;
	sbulkparam->model_tau	= ssysparam->model_tau;
	sbulkparam->model_gamma	= ssysparam->model_gamma;
	sbulkparam->model_kappa	= ssysparam->model_kappa;
}


/**
 * \brief	Save input parameters;
 *
 * \param	sbulkparam1		Input parameters about the left bulk phase;
 * \param	sbulkparam2		Input parameters about the right bulk phase;
 * \param	ssysparam		Input parameters about the interface system;
 *
 */
void fn_param_save (	stu_bulk_param		*sbulkparam1,
						stu_bulk_param		*sbulkparam2,
						stu_system_param	*ssysparam,
						const	char		*key )
{
	/* file for saving parameters; */
	char fname[FILELEN];
	sprintf(fname, "%s/parameter_%s.dat", rsltDir, key);

	/* save parameters; */
	FILE *fp = fopen(fname, "w");
	fprintf(fp, "#+++++++++++++++++++++++++++++++++++++++++++++++++++++++#\n");
	fprintf(fp, "#		input parameters: interface framework			#\n");
	fprintf(fp, "#		lines starting with # are comments				#\n");
	fprintf(fp, "#		must have spaces around the equal sign \"=\"		#\n");
	fprintf(fp, "#+++++++++++++++++++++++++++++++++++++++++++++++++++++++#\n\n");
	fprintf(fp, "## the type of main code, including\n");
	fprintf(fp, "##      stable_system:  stable interface system;\n");
	fprintf(fp, "##      stable_bulk1:   stable left  bulk phase;\n");
	fprintf(fp, "##      stable_bulk2:   stable right bulk phase;\n");
	fprintf(fp, "##      proj_bulk1:     project left bulk phase;\n");
	fprintf(fp, "##      proj_bulk2:     project right bulk phase;\n");
	fprintf(fp, "##      com_bulk1:      left  bulk phase in the common space;\n");
	fprintf(fp, "##      com_bulk2:      right bulk phase in the common space;\n");
	fprintf(fp, "##      com_projmat:    the common projection matrix;\n");
	fprintf(fp, "##      disp_density:   display densities according the files;\n\n");
	fprintf(fp, "nprocs\t\t\t\t= %d\n", nprocs);
	fprintf(fp, "main_type\t\t\t= %s\n\n\n", main_type);
	fprintf(fp, "#+++++++++++++++++++++++++++++++++++++++++++++++++++++++#\n");
	fprintf(fp, "#		model parameters including length scales	    #\n");
	fprintf(fp, "#+++++++++++++++++++++++++++++++++++++++++++++++++++++++#\n\n");
	fprintf(fp, "model_type\t\t\t= %s", model_type);
	fprintf(fp, "\t\t\t\t\t# the type of Landau model;\n");
	fprintf(fp, "print_level\t\t\t= %d", ssysparam->print_level);
	fprintf(fp, "\t\t\t\t\t\t# 0 is not displayed;\n\n");
	fprintf(fp, "scale_num\t\t\t= %d", ssysparam->scale_num);
	fprintf(fp, "\t\t\t\t\t\t# the number of length scales;\n");
	fprintf(fp, "## the values of length scale;\n");
	fprintf(fp, "scale_val\t\t\t= ");
	for ( int i = 0; i < ssysparam->scale_num; i++ )
		fprintf(fp, "% .15f\t", ssysparam->scale_val[i]);
	fprintf(fp, "\n\n");
	fprintf(fp, "model_xi\t\t\t= % .15f", ssysparam->model_xi);
	fprintf(fp, "\t# the penalty factor;\n");
	fprintf(fp, "model_tau\t\t\t= % .15f", ssysparam->model_tau);
	fprintf(fp, "\t# the coefficient before quadratic term;\n");
	fprintf(fp, "model_gamma\t\t\t= % .15f", ssysparam->model_gamma);
	fprintf(fp, "\t# the coefficient before cubic term;\n");
	fprintf(fp, "model_kappa\t\t\t= % .15f", ssysparam->model_kappa);
	fprintf(fp, "\t# the coefficient before quartic term;\n\n\n");
	fn_stu_bulk_param_save (sbulkparam1, fp);		// the left bulk phase;
	fn_stu_bulk_param_save (sbulkparam2, fp);		// the right bulk phase;
	fprintf(fp, "#+++++++++++++++++++++++++++++++++++++++++++++++++++++++#\n");
	fprintf(fp, "#			spatial discrete parameters					#\n");
	fprintf(fp, "#		parameters for interface framework				#\n");
	fprintf(fp, "#+++++++++++++++++++++++++++++++++++++++++++++++++++++++#\n\n");
	fprintf(fp, "GJP_degree\t\t= %d", ssysparam->GJP_degree);
	fprintf(fp, "\t\t\t\t\t# the degree of General Jacobi Polynomial up to GJP_degree;\n");
	fprintf(fp, "LGL_num\t\t\t= %d", ssysparam->LGL_num);
	fprintf(fp, "\t\t\t\t\t# the number of Legendre Gauss+Lobatto points; i.e. discrete x;\n");
	fprintf(fp, "Four_num\t\t= %d", ssysparam->Four_num);
	fprintf(fp, "\t\t\t\t\t# the number of Fourier discrete points along each direction;\n\n");
	fprintf(fp, "x_range\t\t\t= %.15e", ssysparam->x_range);
	fprintf(fp, "\t# the distance between the two anchoring planes; (* 2);\n");
	fprintf(fp, "smooth\t\t\t= %.15e", ssysparam->smooth);
	fprintf(fp, "\t# the smooth constant to connect the two bulk phases;\n");
	fprintf(fp, "initDist1\t\t= % .15e", ssysparam->initDist1);
	fprintf(fp, "\t# the end position of the left bulk phase (initial value);\n");
	fprintf(fp, "initDist2\t\t= % .15e", ssysparam->initDist2);
	fprintf(fp, "\t# the end position of the right bulk phase (initial value);\n");
	fprintf(fp, "searchReg\t\t= %d", ssysparam->searchReg);
	fprintf(fp, "\t# the searching range in the computation of common 'rotateProjBoxMat';\n");
	fprintf(fp, "adjustReg\t\t= %d", ssysparam->adjustReg);
	fprintf(fp, "\t# the adjust range in the test of common 'rotateProjBoxMat';\n\n\n");
	fprintf(fp, "#+++++++++++++++++++++++++++++++++++++++++++++++++++++++#\n");
	fprintf(fp, "#		iteration parameters for interface framework	#\n");
	fprintf(fp, "#+++++++++++++++++++++++++++++++++++++++++++++++++++++++#\n\n");
	fprintf(fp, "massFlag\t\t= %d", ssysparam->massFlag);
	fprintf(fp, "\t# the flag of mass conservation; default 1 for forced mass conservation;\n");
	fprintf(fp, "tol\t\t\t\t= %.15e", ssysparam->tol);
	fprintf(fp, "\t# the maximal tolerance error;\n");
	fprintf(fp, "tolham\t\t\t\t= %.15e", ssysparam->tolham);
	fprintf(fp, "\t# the maximal tolerance error of Hamilton;\n");
	fprintf(fp, "step_size\t\t= %.15e", ssysparam->step_size);
	fprintf(fp, "\t# the time steeping size;\n");
	fprintf(fp, "step_min\t\t= %.15e", ssysparam->step_min);
	fprintf(fp, "\t# the lower bound of time steeping size;\n");
	fprintf(fp, "step_max\t\t= %.15e", ssysparam->step_max);
	fprintf(fp, "\t# the upper bound of time steeping size;\n");
	fprintf(fp, "iter_load\t\t= %d", ssysparam->iter_load);
	fprintf(fp, "\t\t\t\t\t# the loading iterator;\n");
	fprintf(fp, "iter_max\t\t= %d", ssysparam->iter_max);
	fprintf(fp, "\t\t\t\t\t# the maximal iterator;\n");
	fprintf(fp, "print_step\t\t= %d", ssysparam->print_step);
	fprintf(fp, "\t\t\t\t\t# the stepping size for printing data;\n");
	fprintf(fp, "save_step\t\t= %d", ssysparam->save_step);
	fprintf(fp, "\t\t\t\t\t# the stepping size for saving data;\n");
	fprintf(fp, "save_type\t\t= %s", ssysparam->save_type);
	fprintf(fp, "\t\t\t\t\t# the type of saving data; \n");
	fprintf(fp, "\t\t\t\t\t\t\t\t\t\t# followed by Fourier coefficient, density, plane wave;\n");
	fprintf(fp, "\t\t\t\t\t\t\t\t\t\t# y: yes; n: no;\n\n\n");
	fprintf(fp, "#+++++++++++++++++++++++++++++++++++++++++++++++++++++++#\n");
	fprintf(fp, "#		iteration parameters for Newton-PCG method		#\n");
	fprintf(fp, "#+++++++++++++++++++++++++++++++++++++++++++++++++++++++#\n\n");
	fprintf(fp, "newton_tol\t\t\t= %.15e", ssysparam->newton_tol);
	fprintf(fp, "\t# the maximal tolerance error of end the first-order method;\n");
	fprintf(fp, "newton_step_size\t= %.15e", ssysparam->newton_step_size);
	fprintf(fp, "\t# the time steeping size of Newton method;\n");
	fprintf(fp, "pcg_type\t\t\t= %d", ssysparam->pcg_type);
	fprintf(fp, "\t# the type of pre-conditioner;\n");
	fprintf(fp, "pcg_iter_max\t\t= %d", ssysparam->pcg_iter_max);
	fprintf(fp, "\t# the maximal iteration of PCG;\n");
	fprintf(fp, "pcg_print_step\t\t= %d", ssysparam->pcg_print_step);
	fprintf(fp, "\t# the stepping size for printing data in the calculation of PCG;\n");
	fprintf(fp, "pcg_delta_coeff\t\t= %.15e", ssysparam->pcg_delta_coeff);
	fprintf(fp, "\t# the coefficient of delta in the calculation of pre-conditioner;\n");
	fprintf(fp, "pcg_mu_para\t\t\t= %.15e", ssysparam->pcg_mu_para);
	fprintf(fp, "\t# the initial 'mu_para' for pre-conditioner;\n\n\n");
	fprintf(fp, "#+++++++++++++++++++++++++++++++++++++++++++++++++++++++#\n");
	fprintf(fp, "#		plotting parameters for interface framework		#\n");
	fprintf(fp, "#+++++++++++++++++++++++++++++++++++++++++++++++++++++++#\n\n");
	fprintf(fp, "y_start\t\t\t= %.15e", ssysparam->y_start);
	fprintf(fp, "\t# the start of plotting interface along y-direction;\n");
	fprintf(fp, "z_start\t\t\t= %.15e", ssysparam->z_start);
	fprintf(fp, "\t# the start of plotting interface along z-direction;\n");
	fprintf(fp, "y_range\t\t\t= %.15e", ssysparam->y_range);
	fprintf(fp, "\t# the plot range along y-direction;\n");
	fprintf(fp, "z_range\t\t\t= %.15e", ssysparam->z_range);
	fprintf(fp, "\t# the plot range along z-direction for 3D phases;\n");
	fprintf(fp, "y_num\t\t\t= %d", ssysparam->y_num);
	fprintf(fp, "\t\t\t\t\t# the number of discrete points along y-direction;\n");
	fprintf(fp, "z_num\t\t\t= %d", ssysparam->z_num);
	fprintf(fp, "\t\t\t\t\t# the number of discrete points along z-direction;\n");
	fprintf(fp, "skip_half\t\t= %d", ssysparam->skip_half);
	fprintf(fp, "\t\t\t\t\t# 1 to skip the half of Four_num in plotting, default 0;\n\n\n");
	fprintf(fp, "#+++++++++++++++++++++++++++++++++++++++++++++++++++++++#\n");
	fprintf(fp, "#			parameters for representing error			#\n");
	fprintf(fp, "#+++++++++++++++++++++++++++++++++++++++++++++++++++++++#\n\n");
	fprintf(fp, "err_y_start\t\t\t= %.15e", ssysparam->err_y_start);
	fprintf(fp, "\t# the start of plotting interface along y-direction;\n");
	fprintf(fp, "err_z_start\t\t\t= %.15e", ssysparam->err_z_start);
	fprintf(fp, "\t# the start of plotting interface along z-direction;\n");
	fprintf(fp, "err_y_range\t\t\t= %.15e", ssysparam->err_y_range);
	fprintf(fp, "\t# the plot range along y-direction;\n");
	fprintf(fp, "err_z_range\t\t\t= %.15e", ssysparam->err_z_range);
	fprintf(fp, "\t# the plot range along z-direction for 3D phases;\n");
	fprintf(fp, "err_y_num\t\t\t= %d", ssysparam->err_y_num);
	fprintf(fp, "\t\t\t\t\t# the number of discrete points along y-direction;\n");
	fprintf(fp, "err_z_num\t\t\t= %d", ssysparam->err_z_num);
	fprintf(fp, "\t\t\t\t\t# the number of discrete points along z-direction;\n\n\n");
	fprintf(fp, "#+++++++++++++++++++++++++++++++++++++++++++++++++++++++#\n");
	fprintf(fp, "#					parameter to recover				#\n");
	fprintf(fp, "#+++++++++++++++++++++++++++++++++++++++++++++++++++++++#\n\n");
	fprintf(fp, "recoverTOL\t\t\t\t= %.15e", ssysparam->recoverTOL);
	fprintf(fp, "\t# the maximal spectra mode;\n");

	if ( strcmp( key, "opt" ) == 0 )
	{
		fprintf(fp, "\n\n#+++++++++++++++++++++++++++++++++++++++++++++++++++++++#\n");
		fprintf(fp, "#		parameters for the common projection matrix		#\n");
		fprintf(fp, "#+++++++++++++++++++++++++++++++++++++++++++++++++++++++#\n\n");
		fprintf(fp, "# com_projmat_way = calculate: calculate the common projection matrix;\n");
		fprintf(fp, "# com_projmat_way = direct:    directly input the common projection matrix;\n\n");
		fprintf(fp, "com_projmat_way\t\t= %s\n", ssysparam->com_projmat_way);
		fprintf(fp, "com_projmat_size\t= %d\t%d\n", ssysparam->rotateProjBoxMat.row, 
				ssysparam->rotateProjBoxMat.col);
		fprintf(fp, "com_projmat_mat\t\t= \n");
		fn_matrix_save ( ssysparam->rotateProjBoxMat, fp);
		fprintf(fp, "com_projmat_coeff\t= \n");
		for ( int i = 0; i < ssysparam->coeffmat.row; i++ )
		{
			fprintf(fp, "\t");
			for ( int j = 0; j < ssysparam->coeffmat.col; j++ )
				fprintf(fp, "\t% d", ssysparam->coeffmat.val[i][j]);
			fprintf(fp, "\n");
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
}


/**
 * \brief	Save input parameters about bulk phases;
 *
 * \param	sbulkparam		Input parameters about the bulk phase;
 * \param	fp				The pointer of saving file;
 *
 */
void fn_stu_bulk_param_save (	stu_bulk_param		*sbulkparam,
									FILE			*fp )
{
	int sflag = sbulkparam->sflag;
	fprintf(fp, "#+++++++++++++++++++++++++++++++++++++++++++++++++++++++#\n");
	if ( sflag == 1 )
		fprintf(fp, "#			parameters for left bulk phase				#\n");
	else if ( sflag == 2 )
		fprintf(fp, "#			parameters for right bulk phase				#\n");
	fprintf(fp, "#+++++++++++++++++++++++++++++++++++++++++++++++++++++++#\n\n");
	fprintf(fp, "bulk%d_phase\t\t\t= %s", sflag, sbulkparam->phase);
	fprintf(fp, "\t\t\t\t\t# the type of the left bulk phase;\n");
	fprintf(fp, "bulk%d_print_level\t= %d", sflag, sbulkparam->print_level);
	fprintf(fp, "\t\t\t\t\t\t# 0 is not displayed;\n\n");
	fprintf(fp, "bulk%d_motion_order\t= %s", sflag, sbulkparam->motion_order);
	fprintf(fp, "\t\t\t\t\t# the motion order including rotate_transl, transl_rotate;\n");
	fprintf(fp, "bulk%d_rotate_order\t= %s", sflag, sbulkparam->rotate_order);
	fprintf(fp, "\t\t\t\t\t# the rotation order; such as xyz;\n");
	fprintf(fp, "bulk%d_rotate_angle\t= ", sflag);
	for ( int i = 0; i < sbulkparam->rotate_angle.len; i++ )
		fprintf(fp, "% .15f\t", sbulkparam->rotate_angle.val[i]);
	fprintf(fp, "# rotate bulk phase by 'rotate_order';\n");
	fprintf(fp, "bulk%d_transl_var\t= ", sflag);
	for ( int i = 0; i < sbulkparam->transl_var.len; i++ )
		fprintf(fp, "% .15f\t", sbulkparam->transl_var.val[i]);
	fprintf(fp, "# translation along x,y,z-direction;\n");
	fprintf(fp, "bulk%d_Four_num\t\t= %d", sflag, sbulkparam->Four_num);
	fprintf(fp, "\t\t\t\t\t\t# the number of Fourier discrete points along each direction;\n");
	fprintf(fp, "bulk%d_enlarge\t\t= %.15f", sflag, sbulkparam->enlarge);
	fprintf(fp, "\t# enlarge the plotting range for bulk phase;\n");
	fprintf(fp, "bulk%d_plot_num\t\t= %d", sflag, sbulkparam->plot_num);
	fprintf(fp, "\t\t\t\t\t# the number of plotting points along each direction;\n");
	fprintf(fp, "bulk%d_skip_half\t\t= %d", sflag, sbulkparam->skip_half);
	fprintf(fp, "\t\t\t\t\t# 1 to skip the half of Four_num in plotting, default 0;\n\n");
	fprintf(fp, "bulk%d_tol\t\t\t= %.15e", sflag, sbulkparam->tol);
	fprintf(fp, "\t# the maximal tolerance error;\n");
	fprintf(fp, "bulk%d_step_size\t\t= %.15e", sflag, sbulkparam->step_size);
	fprintf(fp, "\t# the time stepping size;\n");
	fprintf(fp, "bulk%d_iter_max\t\t= %d", sflag, sbulkparam->iter_max);
	fprintf(fp, "\t\t\t\t\t# the maximal iterator;\n");
	fprintf(fp, "bulk%d_print_step\t= %d", sflag, sbulkparam->print_step);
	fprintf(fp, "\t\t\t\t\t# the stepping size for printing data;\n");
	fprintf(fp, "bulk%d_save_step\t\t= %d", sflag, sbulkparam->save_step);
	fprintf(fp, "\t\t\t\t\t# the stepping size for saving data;\n");
	fprintf(fp, "bulk%d_save_type\t\t= %s", sflag, sbulkparam->save_type);
	fprintf(fp, "\t\t\t\t\t# the type of saving data;\n");
	fprintf(fp, "\t\t\t\t\t\t\t\t\t\t# followed by Fourier coefficient, density, plane wave;\n");
	fprintf(fp, "\t\t\t\t\t\t\t\t\t\t# y: yes; n: no;\n\n");
	fprintf(fp, "bulk%d_opt_tol\t\t= %.15e", sflag, sbulkparam->opt_tol);
	fprintf(fp, "\t\t# the maximal iterator of optimization;\n");
	fprintf(fp, "bulk%d_opt_iter_max\t= %d", sflag, sbulkparam->opt_iter_max);
	fprintf(fp, "\t\t# the maximal tolerance energy error of optimization;\n");
	fprintf(fp, "bulk%d_opt_save_step\t= %d", sflag, sbulkparam->opt_save_step);
	fprintf(fp, "\t\t# the stepping size for saving data in optimization;\n");
	fprintf(fp, "bulk%d_box_bbType\t= %d", sflag, sbulkparam->box_bbType);
	fprintf(fp, "\t\t# the type of BB stepping size in box optimization;\n");
	fprintf(fp, "bulk%d_box_iter_max\t= %d", sflag, sbulkparam->box_iter_max);
	fprintf(fp, "\t\t# the maximal iterator in box optimization;\n");
	fprintf(fp, "bulk%d_box_step_size\t= %.15e", sflag, sbulkparam->box_step_size);
	fprintf(fp, "\t\t# the initial stepping size in box optimization;\n");
	fprintf(fp, "bulk%d_box_tol\t\t= %.15e", sflag, sbulkparam->box_tol);
	fprintf(fp, "\t\t# the maximal tolerance error in box optimization;\n\n");
	fprintf(fp, "bulk%d_dimPhy\t\t= %d", sflag, sbulkparam->dimPhy);	
	fprintf(fp, "\t\t\t\t\t\t# the dimensionality of physical space;\n");
	fprintf(fp, "bulk%d_dimCpt\t\t= %d", sflag, sbulkparam->dimCpt);	
	fprintf(fp, "\t\t\t\t\t\t# the dimensionality of computational space;\n\n");
	fprintf(fp, "bulk%d_boxType\t\t= %s", sflag, sbulkparam->boxType);
	fprintf(fp, "\t\t\t\t\t# the type of computational box including 'cube', 'cuboid', 'hex';\n");
	fprintf(fp, "bulk%d_rcpBox\t\t= \n", sflag);
	fn_matrix_save ( sbulkparam->rcpBox, fp);
	fprintf(fp, "bulk%d_dirBox\t\t= \n", sflag);
	fn_matrix_save(sbulkparam->dirBox, fp);
	fprintf(fp, "# the projection matrix \n");
	fprintf(fp, "bulk%d_projMat\t\t= \n", sflag);
	fn_matrix_save(sbulkparam->projMat, fp);
	fprintf(fp, "# the rotation matrix \n");
	fprintf(fp, "bulk%d_rotateMat\t\t= \n", sflag);
	fn_matrix_save(sbulkparam->rotateMat, fp);
	fprintf(fp, "# the translation vector \n");
	fprintf(fp, "bulk%d_translVec\t\t= \n", sflag);
	for ( int i = 0; i < sbulkparam->dimPhy; i++ )
		fprintf(fp, "\t\t% .15f\n", sbulkparam->translVec.val[i]);
	fprintf(fp, "\n\n# the result of rotateMat' * projMat * rcpBox \n");
	fprintf(fp, "bulk%d_rotateProjBoxMat\t\t= \n", sflag);
	fn_matrix_save(sbulkparam->rotateProjBoxMat, fp);
	fprintf(fp, "# the result of translVec' * projMat * rcpBox \n");
	fprintf(fp, "bulk%d_translProjBoxMat\t\t= \n", sflag);
	for ( int i = 0; i < sbulkparam->dimPhy; i++ )
		fprintf(fp, "\t\t% .15f\n", sbulkparam->translProjBoxVec.val[i]);
	fprintf(fp, "\n");
	fprintf(fp, "bulk%d_trunc_tol\t\t= %.15e", sflag, sbulkparam->trunc_tol);
	fprintf(fp, "\t# the maximal tolerance error for truncation;\n");
	fprintf(fp, "\n\n");
}


/**
 * \brief	print input parameters about bulk phases;
 *
 * \param	sbulkparam		Input parameters about bulk phases;
 *
 */
void fn_stu_bulk_param_print (stu_bulk_param	*sbulkparam)
{
	int sflag = sbulkparam->sflag;
	if ( sbulkparam->print_level > 0 )
	{
		if ( sflag == 1 )
			printf("\t The left bulk phase:  %s.\n", sbulkparam->phase);
		else if ( sflag == 2 )
			printf("\t The right bulk phase: %s.\n", sbulkparam->phase);
		else
		{
			printf("Error using 'fn_stu_bulk_param_print'\n");
			printf("sflag must be 1 or 2.\n");
		}
		printf("\t The print level is %d.\n", sbulkparam->print_level);
		printf("\t Dimensionality and symmetry.\n");
		printf("\t ---> The dimensionality of physical space is %d.\n", sbulkparam->dimPhy);
		printf("\t ---> The dimensionality of computational space is %d.\n", sbulkparam->dimCpt);
		printf("\t ---> The symmetry is %d-fold.\n", sbulkparam->nfold);
		printf("\t Rotation and translation.\n");
		printf("\t ---> Motion order is '%s'.\n", sbulkparam->motion_order);
		printf("\t ---> Rotation order is '%s'.\n", sbulkparam->rotate_order);
		printf("\t ---> Counterclockwise rotation angles are \n\t\t ");
		for ( int i = 0; i < sbulkparam->rotate_angle.len; i++ )
			printf("% .15f\t", sbulkparam->rotate_angle.val[i]);
		printf("\n");
		printf("\t ---> Translation variables are \n\t\t ");
		for ( int i = 0; i < sbulkparam->transl_var.len; i++ )
			printf("% .15f\t", sbulkparam->transl_var.val[i]);
		printf("\n");
		printf("\t Spatial discretization.\n");
		printf("\t ---> The number of Fourier discrete points along each direction is %d.\n", sbulkparam->Four_num);
		printf("\t ---> The plotting range for bulk phase is [0, %.4f * 2PI].\n", sbulkparam->enlarge);
		printf("\t ---> The number of plotting points along each direction is %d.\n", sbulkparam->plot_num);
		printf("\t ---> The flag is %d to determine whether skipping the half of Four_num in plotting.\n", sbulkparam->skip_half);
		printf("\t Iteration parameters.\n");
		printf("\t ---> The maximal tolerance error is %.15e.\n", sbulkparam->tol);
		printf("\t ---> The time stepping size is %.15e.\n", sbulkparam->step_size);
		printf("\t ---> The maximal iterator is %d.\n", sbulkparam->iter_max);
		printf("\t ---> The stepping size for printing data is %d.\n", sbulkparam->print_step);
		printf("\t ---> The stepping size for saving data is %d.\n", sbulkparam->save_step);
		printf("\t ---> The type of saving data is %s. (y: yes; n: no)\n", sbulkparam->save_type);
		printf("\t Parameters in the box optimization.\n");
		printf("\t ---> The maximal tolerance energy error of optimization is %.15e.\n", sbulkparam->opt_tol);
		printf("\t ---> The maximal iterator of optimization is %d.\n", sbulkparam->opt_iter_max);
		printf("\t ---> The stepping size for saving data in optimization is %d.\n", sbulkparam->opt_save_step);
		printf("\t ---> The type of BB stepping size in box optimization is %d.\n", sbulkparam->box_bbType);
		printf("\t ---> The maximal iterator in box optimization is %d.\n", sbulkparam->box_iter_max);
		printf("\t ---> The initial stepping size in box optimization is %.15e.\n", sbulkparam->box_step_size);
		printf("\t ---> The maximal tolerance error in box optimization is %.15e.\n", sbulkparam->box_tol);
		printf("\t Parameters for truncation.\n");
		printf("\t ---> The maximal tolerance error for truncation is %.15e.\n", sbulkparam->trunc_tol);
	}
	if ( sbulkparam->print_level > 1 )
	{
		printf("\t Rotation, translation, projection, computational box.\n");
		printf("\t ---> The type of reciprocal box is %s.\n", sbulkparam->boxType);
		printf("\t ---> The reciprocal box is \n");
		fn_matrix_print ( sbulkparam->rcpBox );
		printf("\t ---> The direct box is \n");
		fn_matrix_print ( sbulkparam->dirBox );
		printf("\t ---> The projection matrix is \n");
		fn_matrix_print ( sbulkparam->projMat );
		printf("\t ---> The rotation matrix is \n");
		fn_matrix_print ( sbulkparam->rotateMat );
		printf("\t ---> The translation vector is \n");
		for ( int i = 0; i < sbulkparam->dimPhy; i++ )
			printf("\t\t% .15f\n", sbulkparam->translVec.val[i]);
	}
	if ( sbulkparam->print_level > 2 )
	{
		printf("\t ---> The result of rotateMat' * projMat * rcpBox is \n");
		fn_matrix_print ( sbulkparam->rotateProjBoxMat );
		printf("\t ---> The result of translVec' * projMat * rcpBox is \n");
		printf("\t");
		for ( int i = 0; i < sbulkparam->dimPhy; i++ )
			printf("\t% .15f", sbulkparam->translProjBoxVec.val[i]);
		printf("\n");
		printf("\t The initial values is \n");
		for ( int i = 0; i < sbulkparam->initNum; i++ )
		{
			printf("\t");
			for ( int j = 0; j < sbulkparam->dimCpt; j++ )
				printf("\t% d", sbulkparam->initIndex.val[i][j]);
			for ( int j = 0; j < 2; j++ )
				printf("\t% .5f", sbulkparam->initCoeff.val[i][j]);
			printf("\n");
		}
	}
	if ( sbulkparam->print_level > 3 )
	{
		printf("\t ---> The model type is %s.\n", model_type);
		printf("\t ---> The number of length scales is %d.\n", sbulkparam->scale_num);
		printf("\t ---> The values of length scales: ");
		for ( int i = 0; i < sbulkparam->scale_num; i++ )
			printf("\t% .15f", sbulkparam->scale_val[i]);
		printf("\n");
		printf("\t ---> model_xi    = % .15f.\n", sbulkparam->model_xi);
		printf("\t ---> model_tau   = % .15f.\n", sbulkparam->model_tau);
		printf("\t ---> model_gamma = % .15f.\n", sbulkparam->model_gamma);
		printf("\t ---> model_kappa = % .15f.\n", sbulkparam->model_kappa);
	}
	printf("\n\n");
}


/**
 * \brief	print input parameters about interface system;
 *
 * \param	ssysparam		Input parameters about the interface;
 *
 */
void fn_stu_system_param_print (stu_system_param	*ssysparam)
{
	if ( ssysparam->print_level > 0 )
	{
		printf("\t The interface system.\n");
		printf("\t The print level is %d.\n", ssysparam->print_level);
		printf("\t The model parameters:\n");
		printf("\t ---> The model type is %s.\n", model_type);
		printf("\t ---> The number of length scales is %d.\n", ssysparam->scale_num);
		printf("\t ---> The values of length scales: ");
		for ( int i = 0; i < ssysparam->scale_num; i++ )
			printf("\t% .15f", ssysparam->scale_val[i]);
		printf("\n");
		printf("\t ---> model_xi    = % .15f.\n", ssysparam->model_xi);
		printf("\t ---> model_tau   = % .15f.\n", ssysparam->model_tau);
		printf("\t ---> model_gamma = % .15f.\n", ssysparam->model_gamma);
		printf("\t ---> model_kappa = % .15f.\n", ssysparam->model_kappa);
		printf("\t Spatial discrete parameters:\n");
		printf("\t ---> The degree of General Jacobi Polynomial up to %d.\n", ssysparam->GJP_degree);
		printf("\t ---> The number of Legendre Gauss-Lobatto points is %d.\n", ssysparam->LGL_num);
		printf("\t ---> The number of Fourier discrete points along each direction is %d.\n", ssysparam->Four_num);
		printf("\t ---> The distance between the two anchoring planes is % .5e * 2.\n", ssysparam->x_range);
		printf("\t ---> The smooth constant to connect the two bulk phases is % .5e.\n", ssysparam->smooth);
		printf("\t ---> The end position of left bulk phase (initial value) is % .5e.\n", ssysparam->initDist1);
		printf("\t ---> The end position of right bulk phase (initial value) is % .5e.\n", ssysparam->initDist2);
		printf("\t ---> The searching range in the computation of common matrix is %d.\n", ssysparam->searchReg);
		printf("\t ---> The adjust range in the test of common matrix is %d.\n", ssysparam->adjustReg);
		printf("\t Iteration parameters:\n");
		printf("\t ---> The flag of mass conservation is %d.\n", ssysparam->massFlag);
		printf("\t ---> The maximal tolerance error is %.6e.\n", ssysparam->tol);
		printf("\t ---> The maximal tolerance error of Hamilton is %.6e.\n", ssysparam->tolham);
		printf("\t ---> The time stepping size is %.6e.\n", ssysparam->step_size);
		printf("\t ---> The lower bound of time stepping size is %.6e.\n", ssysparam->step_min);
		printf("\t ---> The upper bound of time stepping size is %.6e.\n", ssysparam->step_max);
		printf("\t ---> The loading iterator is %d.\n", ssysparam->iter_load);
		printf("\t ---> The maximal iterator is %d.\n", ssysparam->iter_max);
		printf("\t ---> The stepping size for printing data is %d.\n", ssysparam->print_step);
		printf("\t ---> The stepping size for saving data is %d.\n", ssysparam->save_step);
		printf("\t ---> The type of saving data is %s. (y: yes; n: no)\n", ssysparam->save_type);
		printf("\t Iteration parameters of Newton-PCG method.\n");
		printf("\t ---> The maximal tolerance error to end the first-order method is %.15e.\n", 
				ssysparam->newton_tol);
		printf("\t ---> The time stepping size of Newton method is %.15e.\n", ssysparam->newton_step_size);
		printf("\t ---> The type of PCG is %d.\n", ssysparam->pcg_type);
		printf("\t ---> The maximal iteration of PCG is %d.\n", ssysparam->pcg_iter_max);
		printf("\t ---> The stepping size for printing data in the calculation of PCG is %d.\n", 
				ssysparam->pcg_print_step);
		printf("\t ---> The coefficient of delta in the calculation of pre-conditioner is %.15e.\n", 
				ssysparam->pcg_delta_coeff);
		printf("\t ---> The initial 'mu_para' for pre-conditioner is %.15e.\n", ssysparam->pcg_mu_para);
		printf("\t Plotting parameters:\n");
		printf("\t ---> The start of plotting interface along y-direction is % .5e.\n", ssysparam->y_start);
		printf("\t ---> The start of plotting interface along z-direction is % .5e.\n", ssysparam->z_start);
		printf("\t ---> The plot range along y-direction is % .5e.\n", ssysparam->y_range);
		printf("\t ---> The plot range along z-direction is % .5e.\n", ssysparam->z_range);
		printf("\t ---> The number of discrete points along y-direction is %d.\n", ssysparam->y_num);
		printf("\t ---> The number of discrete points along z-direction is %d.\n", ssysparam->z_num);
		printf("\t ---> The flag is %d to determine whether skipping the half of Four_num in plotting.\n", ssysparam->skip_half);
		printf("\t Parameters for representing error:\n");
		printf("\t ---> The start of plotting interface along y-direction is % .5e.\n", ssysparam->err_y_start);
		printf("\t ---> The start of plotting interface along z-direction is % .5e.\n", ssysparam->err_z_start);
		printf("\t ---> The plot range along y-direction is % .5e.\n", ssysparam->err_y_range);
		printf("\t ---> The plot range along z-direction is % .5e.\n", ssysparam->err_z_range);
		printf("\t ---> The number of discrete points along y-direction is %d.\n", ssysparam->err_y_num);
		printf("\t ---> The number of discrete points along z-direction is %d.\n", ssysparam->err_z_num);
		printf("\t Parameter to recover:\n");
		printf("\t ---> The maximal spectra mode is %.6e.\n", ssysparam->recoverTOL);
	}
	if ( ssysparam->print_level > 1 )
	{
		printf("\t Parameters about GJPs:\n");
		printf("\t ---> The degree of general Jacobi polynomials is %d.\n", ssysparam->sGJPv->polyDegree);
		printf("\t ---> In GJPs, alpha = % .4f, beta = % .4f.\n", ssysparam->sGJPv->alpha, ssysparam->sGJPv->beta);
		printf("\t ---> Some parameters for GJPs: nd = %d, xlen = %d, x_range = % .5e.\n",
				ssysparam->sGJPv->nd, ssysparam->sGJPv->xlen, ssysparam->sGJPv->x_range);
	}
	printf("\n\n");
}


/**
 * \brief	save matrix;
 *
 * \param	mat			Printing matrix;
 * \param	fp			The pointer of saving file;
 *
 */
void fn_matrix_save  (	tmat<double>	src,	FILE	*fp )
{
	for ( int i = 0; i < src.row; i++ )
	{
		fprintf(fp, "\t");
		for ( int j = 0; j < src.col; j++ )
			fprintf(fp, "\t% .15f", src.val[i][j]);
		fprintf(fp, "\n");
	}
	fprintf(fp, "\n");
}


/**
 * \brief	print matrix;
 *
 * \param	src:		Printing matrix;
 *
 */
void fn_matrix_print ( tmat<double> src )
{
	for ( int i = 0; i < src.row; i++ )
	{
		printf("\t");
		for ( int j = 0; j < src.col; j++ )
			printf("\t% .15f", src.val[i][j]);
		printf("\n");
	}
}


/**
 * \brief	Adjust model parameters to match different models;
 *
 * \param	sbulkparam		Input the model parameters and save the modified values;
 *
 */
void fn_adjust_bulk_model_param (stu_bulk_param			*sbulkparam)
{
	sbulkparam->model_xi = pow(sbulkparam->model_xi, 2);	// xi^2;
	if ( strcmp(model_type, "LB") == 0 )
	{
		sbulkparam->model_gamma /= 2.0;
		sbulkparam->model_kappa /= 6.0;
	}
}


/**
 * \brief	Adjust model parameters to match different models;
 *
 * \param	ssysparam		Input the model parameters and save the modified values;
 *
 */
void fn_adjust_system_model_param	(	stu_system_param		*ssysparam	)
{
	ssysparam->model_xi = pow(ssysparam->model_xi, 2);	// xi^2;
	if ( strcmp(model_type, "LB") == 0 )
	{
		ssysparam->model_gamma /= 2.0;
		ssysparam->model_kappa /= 6.0;
	}
//	if ( strcmp ( ssysparam->x_range_type, "direct" ) == 0 )
//	{
//		ssysparam->x_range *= 2*PI;
//	}
	ssysparam->sGJPv->x_range = ssysparam->x_range;
}


/**
 * \brief	Read 'main_type' from the parameter file;
 */
void fn_obt_main_type (  )
{
	/* parameter file; */
	char paraFile[FILELEN];
	sprintf(paraFile, "./para/%s/input%s.dat", paraDir, para_flag);
	if ( myrank == 0 )
		printf("parameter file: %s.\n\n", paraFile);

    char    buffer[STRLEN]; // Note: max number of char for each line!
    int     val;
	int		status = 1;		// 1: success; 0: failure;
    FILE    *fp;

    fp = fopen(paraFile, "r");
    if ( fp == NULL )	// invalid file;
	{
		if ( myrank == 0 )
		{
			printf("Error using 'fn_obt_main_type'.\n");
			printf("Unable to read file '%s'. No such file or directory.\n", paraFile);
		}
	}

    while ( status == 1 )
	{
        char    sbuff[STRLEN];
    
        val = fscanf(fp, "%s", buffer);
        if (val==EOF) break;
        if (val!=1) { status = 0; break; }
        if (buffer[0]=='#' || buffer[0]=='%' || buffer[0]=='|' || buffer[0]=='/')
		{
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
            continue;
        }
    
        /* match keyword and scan for value; */
		if ( strcmp(buffer, "main_type") == 0 )
		{
			val = fscanf(fp, "%s", buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			val = fscanf(fp, "%s", sbuff);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			strcpy(main_type, sbuff);
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
		}
	}
    fclose(fp);
	if ( myrank == 0 )
		printf("\t The type of main code is %s.\n\n", main_type);

#if DEBUG_MODE > 1
	if ( myrank == 0 )
		printf("### DEBUG: Reading input (interface) status = %d\n", status);
#endif
}


/**
 * \brief	Obtain 'rsltDir' from the parameter file;
 */
void fn_obt_rsltDir (  )
{
	char	bulk1_phase[STRLEN];	// the type of the left bulk phase;
	char	bulk2_phase[STRLEN];	// the type of the right bulk phase;

	/* parameter file; */
	char paraFile[FILELEN];
	sprintf(paraFile, "./para/%s/input%s.dat", paraDir, para_flag);
	if ( myrank == 0 )
		printf("parameter file: %s.\n\n", paraFile);

    char    buffer[STRLEN]; // Note: max number of char for each line!
    int     val;
	int		status = 1;		// 1: success; 0: failure;
    FILE    *fp;

    fp = fopen(paraFile, "r");
    if ( fp == NULL )	// invalid file;
	{
		if ( myrank == 0 )
		{
			printf("Error using 'fn_obt_rsltDir'.\n");
			printf("Unable to read file '%s'. No such file or directory.\n", paraFile);
		}
	}

    while ( status == 1 )
	{
        char    sbuff[STRLEN];
    
        val = fscanf(fp, "%s", buffer);
        if (val==EOF) break;
        if (val!=1) { status = 0; break; }
        if (buffer[0]=='#' || buffer[0]=='%' || buffer[0]=='|' || buffer[0]=='/')
		{
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
            continue;
        }
    
        /* match keyword and scan for value; */
		if ( strcmp(buffer, "model_type") == 0 )
		{
			val = fscanf(fp, "%s", buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			val = fscanf(fp, "%s", sbuff);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			strcpy(model_type, sbuff);
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
		}
		else if ( strcmp(buffer, "bulk1_phase") == 0 )
		{
			val = fscanf(fp, "%s", buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			val = fscanf(fp, "%s", sbuff);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			strcpy(bulk1_phase, sbuff);
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
		}
		else if ( strcmp(buffer, "bulk2_phase") == 0 )
		{
			val = fscanf(fp, "%s", buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			val = fscanf(fp, "%s", sbuff);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			strcpy(bulk2_phase, sbuff);
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
		}
	}
    fclose(fp);

	sprintf(rsltDir, "./result/%s_%s_%s_%s", model_type, bulk1_phase, bulk2_phase, para_flag);

#if DEBUG_MODE > 1
	if ( myrank == 0 )
		printf("### DEBUG: Reading input (interface) status = %d\n", status);
#endif
}


/**
 * \brief	Read 'com_projmat' from the parameter file;
 */
void fn_obt_com_projmat ( stu_system_param		*ssysparam,
								char			paraFile[FILELEN] )
{
    char    buffer[STRLEN]; // Note: max number of char for each line!
    int     val;
	int		status = 1;		// 1: success; 0: failure;
    FILE    *fp;

    fp = fopen(paraFile, "r");
    if ( fp == NULL )	// invalid file;
	{
		if ( myrank == 0 )
		{
			printf("Error using 'fn_obt_com_projmat'.\n");
			printf("Unable to read file '%s'. No such file or directory.\n", paraFile);
		}
	}

    while ( status == 1 )
	{
		int     ibuff;
        double  dbuff;
        char    sbuff[STRLEN];
    
        val = fscanf(fp, "%s", buffer);
        if (val==EOF) break;
        if (val!=1) { status = 0; break; }
        if (buffer[0]=='#' || buffer[0]=='%' || buffer[0]=='|' || buffer[0]=='/')
		{
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
            continue;
        }
    
        /* match keyword and scan for value; */
		if ( strcmp ( buffer, "com_projmat_size" ) == 0 )
		{
			val = fscanf(fp, "%s", buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			val = fscanf(fp, "%d", &(ibuff));
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			ssysparam->dimRePhy = ibuff;
			val = fscanf(fp, "%d", &(ibuff));
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			ssysparam->dimReCpt = ibuff;
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
		}
		else if ( strcmp ( buffer, "com_projmat_mat" ) == 0 )
		{
			ssysparam->rotateProjBoxMat = fn_tmat_init<double> ( 
					ssysparam->dimRePhy, ssysparam->dimReCpt );
			val = fscanf(fp, "%s", buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
			for ( int i = 0; i < ssysparam->dimRePhy; i++ )
			{
				for ( int j = 0; j < ssysparam->dimReCpt; j++ )
				{
					val = fscanf(fp, "%lf", &(dbuff));
					if (val!=1 || strcmp(buffer,"=")!=0) { status = 0;	break; }
					ssysparam->rotateProjBoxMat.val[i][j] = dbuff;
				}
			}
            if (fscanf(fp, "%*[^\n]")) {/* skip rest of line and do nothing */ };
		}
	}
    fclose(fp);

#if DEBUG_MODE > 1
	if ( myrank == 0 )
		printf("### DEBUG: Reading input (com_projmat) status = %d\n", status);
#endif
}


/**
 * \brief	Read double tmat matrix;
 *
 */
tmat<double> fn_tmat_read_double ( char filename[] )
{
    int row;	// number of rows;
	int col;	// number of columns;
	int ele;	// number of element; equal to row * col;
	tmat<double> rslt;

    int     val;
    FILE    *fp;

    fp = fopen(filename, "r");
    if ( fp == NULL )	// invalid file;
	{
		if ( myrank == 0 )
		{
			printf("Error using 'fn_tmat_read'.\n");
			printf("Unable to read file '%s'. No such file or directory.\n", filename);
		}
	}

	val = fscanf(fp, "%d", &(row));		// read 'row';
	val = fscanf(fp, "%d", &(col));		// read 'col';
	val = fscanf(fp, "%d", &(ele));		// read 'ele';

	/* read 'rslt'; */
	rslt = fn_tmat_init<double> ( row, col );
	for ( int i = 0; i < row; i++ )
	{
		for ( int j = 0; j < col; j++ )
		{
			val = fscanf(fp, "%lf", &(rslt.val[i][j]));
		}
	}
    fclose(fp);

	return rslt;
}


/**
 * \brief	Read double tvec matrix;
 *
 */
tvec<double> fn_tvec_read_double ( char filename[] )
{
	int row;	// number of rows; default value is 1;
	int col;	// number of cols; default value is 1;
    int len;	// number of elements;
	tvec<double> rslt;

    int     val;
    FILE    *fp;

    fp = fopen(filename, "r");
    if ( fp == NULL )	// invalid file;
	{
		if ( myrank == 0 )
		{
			printf("Error using 'fn_tmat_read'.\n");
			printf("Unable to read file '%s'. No such file or directory.\n", filename);
		}
	}

	val = fscanf(fp, "%d", &(row));		// read 'row';
	val = fscanf(fp, "%d", &(col));		// read 'col';
	val = fscanf(fp, "%d", &(len));		// read 'len';

	/* read 'rslt'; */
	rslt = fn_tvec_init<double> ( len );
	for ( int i = 0; i < len; i++ )
	{
		val = fscanf(fp, "%lf", &(rslt.val[i]));
	}
    fclose(fp);
	rslt.row = row;
	rslt.col = col;

	return rslt;
}


/**
 * \brief	Read fftw_complex tvec matrix;
 *
 */
tvec<fftw_complex> fn_tvec_read_complex ( char filename[] )
{
	int row;	// number of rows; default value is 1;
	int col;	// number of cols; default value is 1;
    int len;	// number of elements;
	tvec<fftw_complex> rslt;

    int     val;
    FILE    *fp;

    fp = fopen(filename, "r");
    if ( fp == NULL )	// invalid file;
	{
		if ( myrank == 0 )
		{
			printf("Error using 'fn_tmat_read'.\n");
			printf("Unable to read file '%s'. No such file or directory.\n", filename);
		}
	}

	val = fscanf(fp, "%d", &(row));		// read 'row';
	val = fscanf(fp, "%d", &(col));		// read 'col';
	val = fscanf(fp, "%d", &(len));		// read 'len';

	/* read 'rslt'; */
	rslt = fn_tvec_init<fftw_complex> ( len );
	for ( int i = 0; i < len; i++ )
	{
		val = fscanf(fp, "%lf", &(rslt.val[i][0]));
		val = fscanf(fp, "%lf", &(rslt.val[i][1]));
	}
    fclose(fp);
	rslt.row = row;
	rslt.col = col;

	return rslt;
}


/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
