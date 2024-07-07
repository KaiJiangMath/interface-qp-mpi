/*! \file	CommonRotateProjBoxMat.cpp
 *
 * \brief	Generate the common 'rotateProjBoxMat';
 */

#include "Data.h"
#include "Head.h"
#include "DataOperators.h"
#include "Mytimer.h"
#include "functs.h"


/**
 * \brief	Calculate the common 'rotateProjBoxMat';
 *
 * \param	sbulkparam1:		structure body for the left   'rotateProjBoxMat';
 * \param	sbulkparam2:		structure body for the right  'rotateProjBoxMat';
 * \param	ssysparam:			structure body for the common 'rotateProjBoxMat';
 *
 * \note!	here all 'rotateProjBoxMat' not consider the first row;
 *			the first row is corresponding to x (GJP);
 */
int fn_obt_commom_rotateProjBoxMat  (	stu_bulk_param		*sbulkparam1,
										stu_bulk_param		*sbulkparam2,
										stu_system_param	*ssysparam )
{
	if ( myrank == 0 )
		printf(" <========== Calculation of the common projection matrix ==========> \n\n");
	mytimer_t timer;
	timer.reset();
	timer.start();

	double tol		= 1e-6;
	int int_reg		= ssysparam->searchReg;
	int adjustReg	= ssysparam->adjustReg;

	if ( sbulkparam1->dimPhy != sbulkparam2->dimPhy )
	{
		if ( myrank == 0 )
		{
			printf("Error using 'fn_obt_commom_rotateProjBoxMat'\n");
			printf("'dimPhy' of the left/right bulk phases must be equal.\n");
		}
		return 1;
	}
	int nr = sbulkparam1->dimPhy - 1;

	/* matrix of the left bulk phase; */
	tmat<double> xmat;
	xmat.row = nr;
	xmat.col = sbulkparam1->dimCpt;
	xmat.val = (double **) malloc( sizeof(double *) * xmat.col );
	for ( int i = 0; i < xmat.col; i++ )
		xmat.val[i] = (double *) malloc( sizeof(double) * xmat.row );
	for ( int i = 0; i < xmat.col; i++ )
		for ( int j = 0; j < xmat.row; j++ )
			xmat.val[i][j] = sbulkparam1->rotateProjBoxMat.val[j+1][i];

	/* matrix of the right bulk phase; */
	tmat<double> ymat;
	ymat.row = nr;
	ymat.col = sbulkparam2->dimCpt;
	ymat.val = (double **) malloc( sizeof(double *) * ymat.col );
	for ( int i = 0; i < ymat.col; i++ )
		ymat.val[i] = (double *) malloc( sizeof(double) * ymat.row );
	for ( int i = 0; i < ymat.col; i++ )
		for ( int j = 0; j < ymat.row; j++ )
			ymat.val[i][j] = sbulkparam2->rotateProjBoxMat.val[j+1][i];

	/* connect xmat and ymat; */
	tmat<double> xymat = fn_matrix_connect ( xmat, ymat );
	if ( myrank == 0 )
	{
		printf("xmat:\n");
		for ( int i = 0; i < xmat.row; i++ )
		{
			for ( int j = 0; j < xmat.col; j++ )
				printf("% .15f\t", xmat.val[j][i]);
			printf("\n");
		}
		printf("\n");
		printf("ymat:\n");
		for ( int i = 0; i < ymat.row; i++ )
		{
			for ( int j = 0; j < ymat.col; j++ )
				printf("% .15f\t", ymat.val[j][i]);
			printf("\n");
		}
		printf("\n");
	}


	if ( strcmp ( ssysparam->com_projmat_way, "direct" ) == 0 )
	{
		/* parameter file; */
		char paraFile[FILELEN];
		sprintf(paraFile, "./para/%s/input%s.dat", paraDir, para_flag);
		if ( myrank == 0 )
		{
			printf("parameter file: %s.\n\n", paraFile);
			printf("\t\tdirect input:\n");
		}
		fn_obt_com_projmat ( ssysparam, paraFile );

		int dimRePhy = ssysparam->dimRePhy;
		int dimReCpt = ssysparam->dimReCpt;
		int cplxDofs = pow(int_reg, dimReCpt);
		int **intSpace = (int **) malloc( sizeof(int *) * cplxDofs );
		for ( int i = 0; i < cplxDofs; i++ )
			intSpace[i] = (int *) malloc( sizeof(int) * dimReCpt );
		ssysparam->coeffmat = fn_tmat_init<int> ( xymat.col, dimReCpt );
		fn_tmat_setZero<int> ( ssysparam->coeffmat );

		/* generate all cases of integer coefficients; */
		fn_obt_int_coeff (intSpace, dimReCpt, int_reg, cplxDofs);

		/* calculate the representation error; */
		double errMax = 0.0;
		double err = 0.0;
		for ( int xyInd = 0; xyInd < xymat.col; xyInd++ )
		{
			for ( int j0 = 0; j0 < cplxDofs; j0++ )
			{
				err = 0.0;
				for ( int j = 0; j < dimRePhy; j++ )
				{
					double errTmp = 0.0;
					for ( int i = 0; i < dimReCpt; i++ )
						errTmp += intSpace[j0][i] * ssysparam->rotateProjBoxMat.val[j][i];
					errTmp = fabs( errTmp - xymat.val[xyInd][j] );
					err = ( err > errTmp ? err : errTmp );
				}
				if ( fabs(err) < tol )
				{
					errMax = ( err > errMax ? err : errMax );
					for ( int i = 0; i < dimReCpt; i++ )
						ssysparam->coeffmat.val[xyInd][i] = intSpace[j0][i];
					break;
				}
			}
		}

		if ( myrank == 0 )
		{
			printf("\t ---> the coefficient matrix is \n");
			for ( int i = 0; i < ssysparam->coeffmat.row; i++ )
			{
				printf("\t    ");
				for ( int j0 = 0; j0 < ssysparam->coeffmat.col; j0++ )
					printf("%d\t", ssysparam->coeffmat.val[i][j0]);
				printf("\n");
			}
			printf("\t ---> the error between the original matrix and ");
			printf("the new representation is % .5e\n", errMax);
		}

		/* release memory; */
		for ( int i = 0; i < cplxDofs; i++ )
			free(intSpace[i]);
		free(intSpace);
	}
	else
	{
		/* calculate minimal value of 'xymat'; */
		double minVal = 1.0e9;
		for ( int i = 0; i < xymat.row; i++ )
		{
			for ( int j = 0; j < xymat.col; j++ )
			{
				double tmp = fabs(xymat.val[j][i]);
				minVal = ( minVal < tmp ? minVal : tmp );
			}
		}

		/* check if the common matrix is a identical matrix 
		 * multiplying a constant; */
		int errFlag = 0;
		double err = 0.0;
		for ( int i = 0; i < xymat.row; i++ )
		{
			for ( int j = 0; j < xymat.col; j++ )
			{
				int intCoeff = round( xymat.val[j][i] / minVal );
				err = xymat.val[j][i] - intCoeff * minVal;
//				printf("err = %.4e\n", err);
				if ( fabs(err) > tol )
				{
					errFlag ++;
					break;
				}
			}
		}

//		if ( myrank == 0 ) printf("minVal = %.6f, errFlag = %d\n", minVal, errFlag);

		/* memory allocation for the integer coefficients; */
		tmat<int> coeffmat;
		coeffmat.row = xymat.col;
		coeffmat.col = xymat.col;
		coeffmat.val = (int **) malloc( sizeof(int *) * coeffmat.col );
		for ( int i = 0; i < coeffmat.col; i++ )
			coeffmat.val[i] = (int *) malloc( sizeof(int) * coeffmat.row );
		for ( int i = 0; i < coeffmat.col; i++ )
			for ( int j = 0; j < coeffmat.row; j++ )
				coeffmat.val[i][j] = 0;
	
		/* calculate the common 'rotateProjBoxMat'; */
		tmat<double> rsltmat;
		if ( errFlag > 0 )
		{
			rsltmat = fn_matrix_rational_dependent ( xymat, coeffmat, tol, int_reg, adjustReg );
			coeffmat.col = rsltmat.col;	
		}
		else
		{
			rsltmat.row = nr;
			rsltmat.col = nr;
			rsltmat.val = (double **) malloc( sizeof(double *) * rsltmat.col );
			for ( int i = 0; i < rsltmat.col; i++ )
				rsltmat.val[i] = (double *) malloc( sizeof(double) * nr );
			for ( int i = 0; i < rsltmat.col; i++ )
			{
				for ( int j = 0; j < nr; j++ )
				{
					if ( i == j )
						rsltmat.val[i][j] = minVal;
					else
						rsltmat.val[i][j] = 0.0;
				}
			}

			/* release the extra memory; */
			for ( int i = nr; i < coeffmat.col; i++ )
				free(coeffmat.val[i]);
			coeffmat.col = nr;
			/* identical matrix; */
			for ( int i = 0; i < coeffmat.row; i++ )
			{
				for ( int j = 0; j < nr; j++ )
				{
					if ( i == j || i == j + nr )
						coeffmat.val[i][j] = 1;
				}
			}
		}

		/* save results; */
		ssysparam->dimRePhy	= rsltmat.row;
		ssysparam->dimReCpt	= rsltmat.col;
		/* save roateProjBoxMat; */
		ssysparam->rotateProjBoxMat = fn_tmat_init<double> ( 
									ssysparam->dimRePhy, ssysparam->dimReCpt );
		for ( int i = 0; i < ssysparam->dimRePhy; i++ )	// transpose;
			for ( int j = 0; j < ssysparam->dimReCpt; j++ )
				ssysparam->rotateProjBoxMat.val[i][j] = rsltmat.val[j][i];
		/* save coeffmat; */
		int coeffmatlen = sbulkparam1->dimCpt + sbulkparam2->dimCpt;
		ssysparam->coeffmat = fn_tmat_init<int> ( coeffmatlen, coeffmat.col );
		for ( int i = 0; i < coeffmatlen; i++ ) // transpose;
			for ( int j = 0; j < coeffmat.col; j++ )
				ssysparam->coeffmat.val[i][j] = 0;
		int ind0 = 0;
		for ( int i = 0; i < xmat.col; i++ )	// xmat part;
		{
			if ( normRealInfty(xmat.val[i], nr) > ZEROTOL )	// start position;
			{
				for ( int j = 0; j < coeffmat.col; j++ )
					ssysparam->coeffmat.val[i][j] = coeffmat.val[j][ind0];
				ind0 ++ ;
			}
		}
		for ( int i = 0; i < ymat.col; i++ )	// ymat part;
		{
			if ( normRealInfty(ymat.val[i], nr) > ZEROTOL )	// start position;
			{
				for ( int j = 0; j < coeffmat.col; j++ )
					ssysparam->coeffmat.val[xmat.col+i][j] = coeffmat.val[j][ind0];
				ind0 ++ ;
			}
		}

		if ( myrank == 0 )
		{
			printf("\t ---> the coefficient matrix is \n");
			for ( int i = 0; i < ssysparam->coeffmat.row; i++ )
			{
				printf("\t    ");
				for ( int j0 = 0; j0 < ssysparam->coeffmat.col; j0++ )
					printf("%d\t", ssysparam->coeffmat.val[i][j0]);
				printf("\n");
			}
		}

		/* release memory; */
		for ( int i = 0; i < rsltmat.col; i++ )
			free(rsltmat.val[i]);
		free(rsltmat.val);
		for ( int i = 0; i < coeffmat.col; i++ )
			free(coeffmat.val[i]);
		free(coeffmat.val);

		MPI_Barrier ( MPI_COMM_WORLD );
		timer.pause();
		if ( myrank == 0 )
		{
			printf("\t ***** time cost of calculation of the common projection matrix: ");
			printf("%f seconds *****\n\n", timer.get_current_time());
		}
	}

	if ( myrank == 0 )
	{
		printf("\t ---> The common projection matrix (%s) is \n", ssysparam->com_projmat_way);
		fn_matrix_print ( ssysparam->rotateProjBoxMat );
		printf("\n");
	}

	/* release memory; */
	for ( int i = 0; i < xmat.col; i++ )
		free(xmat.val[i]);
	free(xmat.val);
	for ( int i = 0; i < ymat.col; i++ )
		free(ymat.val[i]);
	free(ymat.val);
	for ( int i = 0; i < xymat.col; i++ )
		free(xymat.val[i]);
	free(xymat.val);

	/* create fold; */
	/*
	mkdir("./result/", 0755);
	mkdir(rsltDir, 0755);
	char fname[FILELEN];
	sprintf(fname, "%s/parameter_opt.dat", rsltDir);	// file for saving parameters;
	*/

	/* save parameters; */
	/*
	if ( myrank ==  0 )
	{
		FILE *fp = fopen(fname, "a");
		fprintf(fp, "\n\n#+++++++++++++++++++++++++++++++++++++++++++++++++++++++#\n");
		fprintf(fp, "#		the way to get the common projection matrix		#\n");
		fprintf(fp, "#+++++++++++++++++++++++++++++++++++++++++++++++++++++++#\n\n");
		fprintf(fp, "# com_projmat_way = calculate: calculate the common projection matrix;\n");
		fprintf(fp, "# com_projmat_way = direct:    directly input the common projection matrix;\n\n");
		fprintf(fp, "com_projmat_way\t\t= %s\n", ssysparam->com_projmat_way);
		fprintf(fp, "com_projmat_size\t= %d\t%d\n", ssysparam->dimRePhy, ssysparam->dimReCpt);
		fprintf(fp, "com_projmat_mat\t\t= \n");
		fn_matrix_save ( ssysparam->rotateProjBoxMat, fp );
		fclose(fp);
	}
	*/

	/* file for saving coeffient matrix; */
	if ( myrank == 0 )
	{
		char fname[FILELEN];
		sprintf(fname, "%s/parameter_coeff.dat", rsltDir);
		FILE *fp = fopen(fname, "w");
		for ( int i = 0; i < ssysparam->coeffmat.row; i++ )
		{
			for ( int j = 0; j < ssysparam->coeffmat.col; j++ )
				fprintf(fp, "% d\t", ssysparam->coeffmat.val[i][j]);
			fprintf(fp, "\n");
		}
//		fprintf(fp, "\n");
		fclose(fp);
	}

	return 0;
}


/**
 * \brief	Connect 'rotateProjBoxMat' of two bulk phases;
 *			delete the zero vector;
 */
tmat<double> fn_matrix_connect	(	tmat<double>	xmat, 
									tmat<double>	ymat )
{
	int nr = xmat.row; // or ymat.row;

	/* initialize result; */
	tmat<double> xymat;
	xymat.row = nr;
	xymat.col = xmat.col + ymat.col;
	xymat.val = (double **) malloc( sizeof(double *) * xymat.col );
	for ( int i = 0; i < xymat.col; i++ )
		xymat.val[i] = (double *) malloc( sizeof(double) * nr );

	/* connect; */
	int ind0 = 0;
	for ( int i = 0; i < xmat.col; i++ )	// xmat part;
	{
		if ( normRealInfty(xmat.val[i], nr) > ZEROTOL )	// start position;
		{
			for ( int j = 0; j < nr; j++ )
				xymat.val[ind0][j] = xmat.val[i][j];
			ind0 ++ ;
		}
	}
	for ( int i = 0; i < ymat.col; i++ )	// ymat part;
	{
		if ( normRealInfty(ymat.val[i], nr) > ZEROTOL )	// start position;
		{
			for ( int j = 0; j < nr; j++ )
				xymat.val[ind0][j] = ymat.val[i][j];
			ind0 ++ ;
		}
	}
	for ( int i = ind0; i < xymat.col; i++ ) // release the extra memory;
		free(xymat.val[i]);
	xymat.col = ind0;						 // update the number of columns;

	return xymat;
}


/**
 * \brief	Numerically check whether two vectors are linear dependent
 *				over the rational number domain;
 *
 * \param	xymat:			real matrix;
 *							have been transposed;
 *			coeffmat:		integer matrix for storing coefficients;
 *
 * \return	rsltmat:		return a real matrix;
 */
tmat<double> fn_matrix_rational_dependent (	tmat<double>	xymat, 
											tmat<int>		coeffmat,
												double		tol, 
												int			int_reg,
												int			adjustReg )
{
	int nr = xymat.row;	// invariant; or ymat.row;
//	bool isPrint = true;
	bool isPrint = false;

	/* memory allocation for result; */
	tmat<double> rsltmat;
	rsltmat.row = nr;
	rsltmat.col = xymat.col;
	rsltmat.val = (double **) malloc( sizeof(double *) * rsltmat.col );
	for ( int i = 0; i < rsltmat.col; i++ )
		rsltmat.val[i] = (double *) malloc( sizeof(double) * nr );
	for ( int i = 0; i < rsltmat.col; i++ )	// initialization;
		for ( int j = 0; j < nr; j++ )
			rsltmat.val[i][j] = 0.0;


	/* -------------------------------------------------------------------- */
	/**
	 * find the maximal linearly independent group over the rational number domain;
	 */
	int xyIntCoeff, dpdflag;
	double xyTmp, err, errTmp;
	int intSpaceFlag = 0;					// 0: create 'intSpace';
	int rsltLen = 1;
	for ( int j = 0; j < nr; j++ )			// select the first column vector;
		rsltmat.val[0][j] = xymat.val[0][j];
	coeffmat.val[0][0] = 1;

	/**
	 * check each column vector of 'xymat';
	 *	here is row vector since 'xymat' has been transposed;
	 *	transposing matrix is more convenient;
	 */
	int cplxDofs, comDiv;
	int **intSpace;
	for ( int xyInd = 1; xyInd < xymat.col; xyInd++ )	// the first one has been used;
	{
		if ( intSpaceFlag == 0 )
		{
			cplxDofs = pow(int_reg, rsltLen);
			intSpace = (int **) malloc( sizeof(int *) * cplxDofs );
			for ( int i = 0; i < cplxDofs; i++ )
				intSpace[i] = (int *) malloc( sizeof(int) * rsltLen );
			/* generate all cases of integer coefficients; */
			fn_obt_int_coeff (intSpace, rsltLen, int_reg, cplxDofs);
		}

		/**
		 * find a nonzero element in the 'xyInd'-th column vector of 'xymat';
		 * it should exist since we have deleted all zero column vectors;
		 */
		int nonzeroInd = 0;
		for ( int j = 0; j < nr; j++ )
		{
			if ( fabs(xymat.val[xyInd][j]) > tol )
			{
				nonzeroInd = j;
				break;
			}
		}

		/* search representation; */
		dpdflag = 0;	// the flag to check whether they are linearly independent;
		for ( int j0 = 0; j0 < cplxDofs; j0++ )
		{
			/* calculate the coefficient by one element; */
			xyTmp = 0.0;
			for ( int i = 0; i < rsltLen; i++ )	// 'rsltLen' is the dimensionality of 'intSpace';
				xyTmp += intSpace[j0][i] * rsltmat.val[i][nonzeroInd];
			xyIntCoeff = round( xyTmp / xymat.val[xyInd][nonzeroInd] );
			/* xyIntCoeff cannot be 0 and cannot be greater than int_reg; */
			if ( xyIntCoeff == 0 || fabs(xyIntCoeff) > int_reg ) continue;

			/* calculate the representation error; */
			err = 0.0;
			for ( int j = 0; j < nr; j++ )
			{
				errTmp = 0.0;
				for ( int i = 0; i < rsltLen; i++ )
					errTmp += intSpace[j0][i] * rsltmat.val[i][j];
				errTmp = fabs( errTmp - xyIntCoeff*xymat.val[xyInd][j] );
				err = ( err > errTmp ? err : errTmp );
			}

			if ( err < tol )
			{
				/* linearly dependent; */
				for ( int i = 0; i < rsltLen; i++ )
				{
					coeffmat.val[i][xyInd] = intSpace[j0][i];		 // save coefficients;
					if ( intSpace[j0][i] != 0 )
					{
						comDiv = __gcd(xyIntCoeff, intSpace[j0][i]); // common divisor;
						for ( int j = 0; j < rsltLen; j++ )
							coeffmat.val[i][j] *= (xyIntCoeff/comDiv); // update the calculated 'rsltmat';
						coeffmat.val[i][xyInd] /= comDiv;			 // update coefficients;
						for ( int j = 0; j < nr; j++ )
						{
							rsltmat.val[i][j] = rsltmat.val[i][j] / (xyIntCoeff/comDiv);
						}
					}
				}
				dpdflag = 1;
				break;
			}
		}

		if ( isPrint && myrank == 0 )
		{
			printf("xyInd = %d \t coeffmat: \n", xyInd);
			for ( int i = 0; i < coeffmat.col; i++ )
			{
				for ( int j = 0; j < coeffmat.row; j++ )
					printf("% d\t", coeffmat.val[i][j]);
				printf("\n");
			}
			printf("\n");
		}
		
		/* add the linearly independent column vector; */
		if ( dpdflag == 0 )
		{
			coeffmat.val[rsltLen][xyInd] = 1;
			for ( int j = 0; j < nr; j++ )
				rsltmat.val[rsltLen][j] = xymat.val[xyInd][j];
			rsltLen ++ ;
		}

		/* release memory of 'intSpace'; */
		if ( dpdflag == 0 || xyInd == xymat.col-1 )
		{
			for ( int i = 0; i < cplxDofs; i++ )
				free(intSpace[i]);
			free(intSpace);
			intSpaceFlag = 0;
		}
		else
			intSpaceFlag = 1;
	}

	/* release the extra memory; */
	for ( int i = rsltLen; i < rsltmat.col; i++ )
		free(rsltmat.val[i]);
	rsltmat.col = rsltLen;
	for ( int i = rsltLen; i < coeffmat.col; i++ )
		free(coeffmat.val[i]);
	coeffmat.col = rsltLen;

	/* the maximal absolute value as the check criterion; */
	int check = fn_max_abs ( coeffmat.val, coeffmat.col, coeffmat.row );

	/* calculate error; */
	err = 0.0;
	for ( int i = 0; i < coeffmat.row; i++ )
	{
		for ( int j = 0; j < rsltmat.row; j++ )
		{
			double chetmp = 0.0;
			for ( int j0 = 0; j0 < coeffmat.col; j0++ )
			{
				chetmp += coeffmat.val[j0][i] * rsltmat.val[j0][j];
			}
			chetmp -= xymat.val[i][j];
			chetmp  = fabs(chetmp);
			err = ( err > chetmp ? err : chetmp );
		}
	}

	if ( ( check > 2 && rsltLen > 1 ) || ( err > 1e-2 ) )
	{
		/* calculate coeffmat; */
		int searchLen = pow( int_reg, coeffmat.col ); // positive, negative, and zero;
		int **intSpace = (int **) malloc( sizeof(int *) * searchLen );
		for ( int i = 0; i < searchLen; i++ )
			intSpace[i] = (int *) malloc( sizeof(int) * coeffmat.col );
		/* generate all cases of integer coefficients; */
		fn_obt_int_coeff (intSpace, coeffmat.col, int_reg, searchLen);

		/* calculate coeffmat and its maximal absolute value; */
		int adjustMax = 0, adjustInd = 0, adjustTmp;
		double errNew = 0.0;
		double errTmp0, errTmp1;
		for ( int j0 = 0; j0 < coeffmat.row; j0++ )
		{
			errTmp0 = 0.0;
			for ( int i = 0; i < searchLen; i++ )
			{
				errTmp0 = 0.0;
				for ( int j1 = 0; j1 < rsltmat.row; j1++ )
				{
					errTmp1 = 0.0;
					for ( int j2 = 0; j2 < rsltmat.col; j2++ )
					{
						errTmp1 += intSpace[i][j2] * rsltmat.val[j2][j1];
					}
					errTmp1 = fabs( errTmp1 - xymat.val[j0][j1] );
					errTmp0 = ( errTmp0 > errTmp1 ? errTmp0 : errTmp1 );
				}
				if ( errTmp0 < err )
				{
					for ( int j2 = 0; j2 < rsltmat.col; j2++ )
					{
						coeffmat.val[j2][j0] = intSpace[i][j2];
						adjustTmp = fabs(coeffmat.val[j2][j0]);
						if ( adjustTmp > adjustMax )
						{
							adjustMax = adjustTmp;
							adjustInd = j2;
						}
					}
					if ( errTmp0 < tol ) break;
				}
			}
			errNew += errTmp0;
		}
		if ( errNew < err ) err = errNew; // update error;

		/* print coeffmat before the adjusting operator; */
		if ( isPrint && myrank == 0 )
		{
			printf("\t ---> the coefficient matrix before adjusting is \n");
			for ( int i = 0; i < coeffmat.row; i++ )
			{
				printf("\t    ");
				for ( int j0 = 0; j0 < coeffmat.col; j0++ )
					printf("%d\t", coeffmat.val[j0][i]);
				printf("\n");
			}
			printf("\t ---> the error between the original matrix and ");
			printf("the new representation is % .5e\n", err);
			printf("\t ---> adjust: max value = %d, index = %d\n", adjustMax, adjustInd);
		}

		/* --------------------------- adjust rsltmat --------------------------- */

		/* allocate memory for optimal coeffmat and rsltmat; */
		tmat<int> coeffmatOpt;
		coeffmatOpt.row = coeffmat.row;
		coeffmatOpt.col = coeffmat.col;
		coeffmatOpt.val = (int **) malloc( sizeof(int *) * coeffmatOpt.col );
		for ( int i = 0; i < coeffmatOpt.col; i++ )
			coeffmatOpt.val[i] = (int *) malloc( sizeof(int) * coeffmatOpt.row );
		tmat<double> rsltmatOpt;
		rsltmatOpt.row = rsltmat.row;
		rsltmatOpt.col = rsltmat.col;
		rsltmatOpt.val = (double **) malloc( sizeof(double *) * rsltmatOpt.col );
		for ( int i = 0; i < rsltmatOpt.col; i++ )
			rsltmatOpt.val[i] = (double *) malloc( sizeof(double) * rsltmatOpt.row );

		/* adjust coeffmat and rsltmat; */
		int adjustLen = pow( adjustReg, coeffmat.col ); // positive, negative, and zero;
		int **adjustSpace = (int **) malloc( sizeof(int *) * adjustLen );
		for ( int i = 0; i < adjustLen; i++ )
			adjustSpace[i] = (int *) malloc( sizeof(int) * coeffmat.col );
		/* generate all cases of integer coefficients; */
		fn_opt_int_coeff (adjustSpace, coeffmat.col, adjustReg, adjustLen);

		/* integer matrix to obtain square coefficient matrix; */
		int squareLen = pow( coeffmat.row, coeffmat.col );
		int **squareSpace = (int **) malloc( sizeof(int *) * squareLen );
		for ( int i = 0; i < squareLen; i++ )
			squareSpace[i] = (int *) malloc( sizeof(int) * coeffmat.col );
		/* generate elements of squareSpace; */
		int *k = (int *) malloc( sizeof(int) * coeffmat.col );
		for ( int i = 0; i < coeffmat.col; i++ ) k[i] = 0;
		int squareInd = 0;
		for ( int i = 0; i < squareLen; i++ )
		{
			bool isSame	 = false;
			bool isBreak = true;
			for ( int j0 = 0; j0 < coeffmat.col && isBreak; j0++ )
				for ( int j1 = j0+1; j1 < coeffmat.col && isBreak; j1++ )
					if ( k[j0] == k[j1] )
					{
						isSame  = true;
						isBreak = false;
					}
			if ( !isSame )
			{
				for ( int j = 0; j < coeffmat.col; j++ )
					squareSpace[squareInd][j] = k[j];
				squareInd ++;
			}
			k[coeffmat.col-1] ++;
			for ( int jj = coeffmat.col-1; jj > 0; jj-- )
			{
				if (k[jj] >= coeffmat.col)
				{
					k[jj] = 0;
					k[jj-1] ++;
				}
			}
		}
		free(k);
		if ( isPrint && myrank == 0 )
			printf("squareInd = %d\n", squareInd);
		
		/* adjust coeffmat; */
		int adjustMaxOld = adjustMax;
		int adjustIndOld = adjustInd;
		bool adjustFlag	 = false;
		int *adjustMaxVec = (int *) malloc( sizeof(int) * coeffmat.row );
		for ( int j0 = 0; j0 < coeffmat.row; j0++ )
			adjustMaxVec[j0] = 0;
		/* calculate the maximal absolute coefficient of coeffmat with different cases; */
		for ( int i = 0; i < adjustLen; i++ )
		{
			for ( int j0 = 0; j0 < coeffmat.row; j0++ )
			{
				adjustTmp = 0;
				for ( int j1 = 0; j1 < coeffmat.col; j1++ )
					adjustTmp += adjustSpace[i][j1] * coeffmat.val[j1][j0];
				adjustMaxVec[j0] = adjustTmp;
			}
			/* common divisor; */
			comDiv = __gcd( adjustMaxVec[0], adjustMaxVec[1] );
			for ( int j0 = 1; j0 < coeffmat.row; j0++ )
				comDiv = __gcd( comDiv, adjustMaxVec[j0] );
			/* divide the common divisor; */
			for ( int j0 = 0; j0 < coeffmat.row; j0++ )
				adjustMaxVec[j0] /= comDiv;

			/* calculate the maximal absolute of 'adjustMaxVec'; */
			adjustMax = 0;
			for ( int j0 = 0; j0 < coeffmat.row; j0++ )
			{
				adjustTmp = fabs(adjustMaxVec[j0]);
				if ( adjustTmp > adjustMax )
					adjustMax = adjustTmp; 
			}

			/* update adjustMax and adjustInd; */
			if ( adjustMax < adjustMaxOld )
			{
				if ( isPrint && myrank == 0 )
					printf("\t ---> adjust: max value = %d, index = %d\n", adjustMax, adjustInd);
				adjustFlag	 = true;
				adjustMaxOld = adjustMax;
				adjustIndOld = adjustInd;

				/* copy coeffmat and rsltmat; */
				for ( int j0 = 0; j0 < coeffmat.col; j0++ )
					for ( int j1 = 0; j1 < coeffmat.row; j1++ )
						coeffmatOpt.val[j0][j1] = coeffmat.val[j0][j1];
				for ( int j0 = 0; j0 < rsltmat.col; j0++ )
					for ( int j1 = 0; j1 < rsltmat.row; j1++ )
						rsltmatOpt.val[j0][j1] = rsltmat.val[j0][j1];

				/* adjust coeffmat; */
				for ( int j0 = 0; j0 < coeffmat.row; j0++ )
				{
					adjustTmp = 0;
					for ( int j1 = 0; j1 < coeffmat.col; j1++ )
						adjustTmp += adjustSpace[i][j1] * coeffmat.val[j1][j0] / comDiv;
					coeffmatOpt.val[adjustInd][j0] = adjustTmp;
				}
				if ( isPrint && myrank == 0 )
				{
					printf("adjustSpace[%d]: ", i);
					for ( int j1 = 0; j1 < coeffmat.col; j1++ )
						printf("%d ", adjustSpace[i][j1]);
					printf("\n");	
					printf("coeffmatOpt: \n");
					for ( int j0 = 0; j0 < coeffmat.row; j0++ )
					{
						for ( int j1 = 0; j1 < coeffmat.col; j1++ )
							printf("%d ", coeffmatOpt.val[j1][j0]);
						printf("\n");
					}
				}

				/* calculate rsltmat; */
				/* find inverse coefficient matrix; */
				double **coeffSquare = (double **) malloc( sizeof(double *) * coeffmat.col );
				for ( int j0 = 0; j0 < coeffmat.col; j0++ )
					coeffSquare[j0] = (double *) malloc( sizeof(double) * coeffmat.col );
				double **coeffSquareInv = (double **) malloc( sizeof(double *) * coeffmat.col );
				for ( int j0 = 0; j0 < coeffmat.col; j0++ )
					coeffSquareInv[j0] = (double *) malloc( sizeof(double) * coeffmat.col );
				/* calculate det to check if the matrix can inverse; */
				int obtInd = 0;
				double squareDet = 0.0;
				for ( int j0 = 0; j0 < squareInd; j0++ )
				{
					for ( int j1 = 0; j1 < coeffmat.col; j1++ ) // row;
					{
						int rowNum = squareSpace[j0][j1];
						for ( int j2 = 0; j2 < coeffmat.col; j2++ ) // col;
							coeffSquare[j1][j2] = (double) coeffmatOpt.val[j2][rowNum];
					}
					squareDet = det ( coeffSquare, coeffmat.col );
					if ( fabs(squareDet) > 0 )
					{
						obtInd = j0;
						break;
					}
				}
				if ( fabs(squareDet) > 0 )
				{
					/* calculate inverse matrix of coeffSquare; */
					matrix_inverse ( coeffSquare, coeffSquareInv, coeffmat.col );
					/* calculate rsltmat; */
					for ( int j1 = 0; j1 < coeffmat.col; j1++ ) // row;
					{
						for ( int j3 = 0; j3 < rsltmat.row; j3++ )
						{
							double rsltTmp = 0.0;
							for ( int j2 = 0; j2 < coeffmat.col; j2++ ) // col;
							{
								int rowNum = squareSpace[obtInd][j3];
								rsltTmp += coeffSquareInv[j1][j2] * xymat.val[j2][rowNum];
							}
							rsltmatOpt.val[j1][j3] = rsltTmp;
						}
					}
					if ( isPrint && myrank == 0 )
					{
						printf("coeffSquareInv:\n");
						for ( int j1 = 0; j1 < coeffmat.col; j1++ )
						{
							for ( int j2 = 0; j2 < coeffmat.col; j2++ )
								printf("%.4f ", coeffSquareInv[j1][j2]);
							printf("\n");
						}
					}
				}

				if ( isPrint && myrank == 0 )
				{
					printf("xymat:\n");
					for ( int j1 = 0; j1 < rsltmat.row; j1++ )
					{
						int rowNum = squareSpace[obtInd][j1];
						for ( int j2 = 0; j2 < coeffmat.col; j2++ )
							printf("%.4f ", xymat.val[j2][rowNum]);
						printf("\n");
					}
					printf("rsltmat:\n");
					for ( int j1 = 0; j1 < coeffmat.col; j1++ )
					{
						for ( int j2 = 0; j2 < rsltmat.row; j2++ )
							printf("%.4f ", rsltmatOpt.val[j1][j2]);
						printf("\n");
					}
				}

				/* release memory; */
				for ( int j0 = 0; j0 < coeffmat.col; j0++ ) free ( coeffSquare[j0] );
				free ( coeffSquare );
				for ( int j0 = 0; j0 < coeffmat.col; j0++ ) free ( coeffSquareInv[j0] );
				free ( coeffSquareInv );

				if ( isPrint && myrank == 0 )
				{
					printf("\t ---> adjustSpace[%d]: ", i);
					for ( int j1 = 0; j1 < rsltmat.col; j1++ )
						printf("%d ", adjustSpace[i][j1]);
					printf("\n");
				}
			}
		}

		if ( adjustFlag )
		{
			/* copy coeffmatOpt to coeffmat; */
			for ( int j0 = 0; j0 < coeffmat.col; j0++ )
				for ( int j1 = 0; j1 < coeffmat.row; j1++ )
					coeffmat.val[j0][j1] = coeffmatOpt.val[j0][j1];
			/* copy rsltmatOpt to rsltmat; */
			for ( int j0 = 0; j0 < rsltmat.col; j0++ )
				for ( int j1 = 0; j1 < rsltmat.row; j1++ )
					rsltmat.val[j0][j1] = rsltmatOpt.val[j0][j1];
		}
		/* the maximal absolute value as the check criterion; */
		check = fn_max_abs ( coeffmat.val, coeffmat.col, coeffmat.row );

		/* calculate error; */
		err = 0.0;
		for ( int i = 0; i < coeffmat.row; i++ )
		{
			for ( int j = 0; j < rsltmat.row; j++ )
			{
				double chetmp = 0.0;
				for ( int j0 = 0; j0 < coeffmat.col; j0++ )
				{
					chetmp += coeffmat.val[j0][i] * rsltmat.val[j0][j];
				}
				chetmp -= xymat.val[i][j];
				chetmp  = fabs(chetmp);
				err = ( err > chetmp ? err : chetmp );
			}
		}

		/* release memory; */
		for ( int i = 0; i < searchLen; i++ )
			free( intSpace[i] );
		free( intSpace );
		for ( int i = 0; i < adjustLen; i++ )
			free( adjustSpace[i] );
		free( adjustSpace );
		for ( int i = 0; i < squareLen; i++ )
			free( squareSpace[i] );
		free( squareSpace );
		for ( int i = 0; i < coeffmatOpt.col; i++ )
			free( coeffmatOpt.val[i] );
		free( coeffmatOpt.val );
		for ( int i = 0; i < rsltmatOpt.col; i++ )
			free( rsltmatOpt.val[i] );
		free( rsltmatOpt.val );
		free( adjustMaxVec );
	}

	/* print result; */
	if ( myrank == 0 )
	{
		printf("\t Common projection matrix: \n");
		printf("\t ---> the maximal absolute coefficient is %d\n", check);
		printf("\t ---> the error between the original matrix and ");
		printf("the new representation is % .5e\n", err);
	}

	return rsltmat;
}


/**
 * \brief	Generate all possible cases of the integer coefficients;
 *
 * \param	intSpace:		return value; all possible cases;
 * \param	dim:		dimensionality;
 * \param	N:			the value range; 1,2,...,N;
 */
void fn_obt_int_coeff (int **intSpace, int dim, int N, int cplxDofs)
{
	int *k = (int *)malloc( sizeof(int) * dim);
	for ( int i = 0; i < dim; i++ ) k[i] = -N/2+1;

	for ( int i = 0; i < cplxDofs; i++ )
	{
		for ( int j = 0; j < dim; j++ )
			intSpace[i][j] = k[j];
		k[dim-1] ++;
		for ( int jj = dim-1; jj > 0; jj-- )
		{
			if (k[jj] > N/2)
			{
				k[jj] = -N/2+1;
				k[jj-1] ++;
			}
		}
	}
	free(k);

	/* sort 'intSpace' by modules; */
	double tmp;
	vector <stu_sort> sort_array(cplxDofs);
	for ( int i = 0; i < cplxDofs; ++i )
	{
		tmp = 0.0;
		for ( int j = 0; j < dim; j++ )
			tmp += pow(intSpace[i][j], 2);
		sort_array[i].ind = i;
		sort_array[i].val = tmp;
	}
	sort(sort_array.begin(), sort_array.end(), fn_compare);

	/* obtain the sorted 'intSpace'; */
	int **rslt = (int **) malloc( sizeof(int *) * cplxDofs );
	for ( int i = 0; i < cplxDofs; i++ )
		rslt[i] = (int *) malloc( sizeof(int) * dim );
	for ( int i = 0; i < cplxDofs; i++ )
		for ( int j = 0; j < dim; j++ )
			rslt[i][j] = intSpace[sort_array[i].ind][j];
	for ( int i = 0; i < cplxDofs; i++ ) // copy;
		for ( int j = 0; j < dim; j++ )
			intSpace[i][j] = rslt[i][j];

	/* release memory; */
	for ( int i = 0; i < cplxDofs; i++ )
		free(rslt[i]);
	free(rslt);
}


/**
 * \brief	Generate all possible cases to optimize the integer coefficients;
 *
 * \param	intSpace:		return value; all possible cases;
 * \param	dim:		dimensionality;
 * \param	N:			the value range; 1,2,...,N;
 */
void fn_opt_int_coeff (int **intSpace, int dim, int N, int cplxDofs)
{
	int *k = (int *)malloc( sizeof(int) * dim);
	for ( int i = 0; i < dim; i++ ) k[i] = -N/2;

	for ( int i = 0; i < cplxDofs; i++ )
	{
		for ( int j = 0; j < dim; j++ )
		{
			if ( fabs(k[j]) == 0 )
				k[j] ++;			// skip 0;
			intSpace[i][j] = k[j];
		}
		k[dim-1] ++;
		for ( int jj = dim-1; jj > 0; jj-- )
		{
			if (k[jj] > N/2)
			{
				k[jj] = -N/2;
				k[jj-1] ++;
			}
		}
	}
	free(k);

	/* sort 'intSpace' by modules; */
	double tmp;
	vector <stu_sort> sort_array(cplxDofs);
	for ( int i = 0; i < cplxDofs; ++i )
	{
		tmp = 0.0;
		for ( int j = 0; j < dim; j++ )
			tmp += pow(intSpace[i][j], 2);
		sort_array[i].ind = i;
		sort_array[i].val = tmp;
	}
	sort(sort_array.begin(), sort_array.end(), fn_compare);

	/* obtain the sorted 'intSpace'; */
	int **rslt = (int **) malloc( sizeof(int *) * cplxDofs );
	for ( int i = 0; i < cplxDofs; i++ )
		rslt[i] = (int *) malloc( sizeof(int) * dim );
	for ( int i = 0; i < cplxDofs; i++ )
		for ( int j = 0; j < dim; j++ )
			rslt[i][j] = intSpace[sort_array[i].ind][j];
	for ( int i = 0; i < cplxDofs; i++ ) // copy;
		for ( int j = 0; j < dim; j++ )
			intSpace[i][j] = rslt[i][j];

	/* release memory; */
	for ( int i = 0; i < cplxDofs; i++ )
		free(rslt[i]);
	free(rslt);
}


/**
 * \brief	Numerically check whether two vectors are linear dependent
 *				over the rational number domain;
 * \param	x, y:			two real vectors (nonzero);
 * \param	z:				return a real vector;
 * \param	n:				the length of x,y;
 *
 * \return	an integer:		0: x,y are nonzero and independent;
 *							1: x,y are nonzero and dependent, select x;
 *							2: x,y are nonzero and dependent, select y;
 *							-1: |y| = 0, x is nonzero;
 *							-2: |x| = 0, y is nonzero;
 *							-3: |x| = |y| = 0;
 */
int fn_check_vector_rational_dependent (double *x,	double *y,	double *z,	
										int		n,	double tol,	int	   int_reg)
{
	if ( normRealInfty(x, n) < tol )
	{
		if ( normRealInfty(y, n) < tol )	// |x| = |y| = 0;
			return -3;
		else
			return -2;						// |x| = 0, y is nonzero;
	}
	else if ( normRealInfty(y, n) < tol )
		return -1;							// |y| = 0, x is nonzero;

	/* check zero; */
	for ( int i = 0; i < n; i++ )
	{
		if ( ( fabs(x[i]) < tol && fabs(y[i]) > tol ) ||
		   ( fabs(x[i]) > tol && fabs(y[i]) < tol ) )
			   return 0;					// x[i], y[i]: only one zero;
	}

	/* check that x[i] and y[i] both are nonzero; */
	int *ind = (int *) malloc( sizeof(int) * n );
	int len = 0;
	for ( int i = 0; i < n; i++ )
	{
		if ( fabs(x[i]) > tol && fabs(y[i]) > tol ) // x[i], y[i]: both nonzero;
		{
			ind[len] = i;
			len++;
		}
	}

	/** 
	 * check whether x,y are linearly dependent over the rational number domain;
	 *		x, y are nonzero;
	 *	 dependent: ax=by, a,b are integers belonging to [1, int_reg);
	 * independent: cannot find two integers a,b belonging to [1, int_reg) 
	 *				to satisfy ax=by;
	 */
	double b, err, tmp;
	int	*intb = (int *) malloc( sizeof(int) * len );
	int  intb0;
	for ( int i = 0; i < n; i++ )
		z[i] = 0.0;
	for ( int a = 1; a < int_reg; a++ )
	{
		err = 0.0;
		for ( int i = 0; i < len; i++ )	// search b;
		{
			b		= a * x[ind[i]] / y[ind[i]];
			intb[i] = round(b);
			tmp		= fabs( b - intb[i] );
			err		= (err > tmp ? err : tmp);
		}
		intb0 = intb[0];
		for ( int i = 1; i < len; i++ )	// all elements of b should be equal;
		{
			if ( fabs(intb0) < tol || fabs(intb[i] - intb0) > tol )
			{
				err = 1.0;
				break;
			}
		}
		if ( err < tol )
		{
			if ( fabs(a) > fabs(intb0) )
			{
				for ( int i = 0; i < len; i++ ) // return x;
					z[ind[i]] = x[ind[i]] / intb0;
				free(ind);	free(intb);
				return 1;		// x,y are nonzero and dependent, select x;
			}
			else
			{
				for ( int i = 0; i < len; i++ ) // return y;
					z[ind[i]] = y[ind[i]] / a;
				free(ind);	free(intb);
				return 2;		// x,y are nonzero and dependent, select y;
			}
		}
	}
	free(ind);	free(intb);
	return 0;					// x,y are nonzero and independent;
}


/**
 * \brief	Numerically check whether two real numbers are linearly dependent 
 *				over the rational number domain;
 *
 * \param	x, y:			two real numbers;
 *
 * \return	an integer:		0: x,y are nonzero and independent;
 *							1: x,y are nonzero and dependent, select x;
 *							2: x,y are nonzero and dependent, select y;
 *							-1: y = 0, x is nonzero;
 *							-2: x = 0, y is nonzero;
 *							-3: x = y = 0;
 *							
 */
int fn_check_number_rational_dependent (double x,	double y,
									    double tol,	int	   int_reg)
{
	if ( fabs(x) < tol )
	{
		if ( fabs(y) < tol )	// x = y = 0;
			return -3;
		else
			return -2;			// x = 0, y is nonzero;
	}
	else if ( fabs(y) < tol )
		return -1;				// y = 0, x is nonzero;

	/** 
	 * check whether x,y are linearly dependent over the rational number domain;
	 *		x, y are nonzero;
	 *	 dependent: ax=by, a,b are integers belonging to [1, int_reg);
	 * independent: cannot find two integers a,b belonging to [1, int_reg) 
	 *				to satisfy ax=by;
	 */
	double b, err;
	int intb;
	for ( int a = 1; a < int_reg; a++ )
	{
		b	 = a * x / y;
		intb = round(b);
		err  = fabs( b - round(b) );
		if ( err < tol )
		{
			if ( fabs(a) > fabs(intb) )
				return 1;		// x,y are nonzero and dependent, select x;
			else
				return 2;		// x,y are nonzero and dependent, select y;
		}
	}
	return 0;					// x,y are nonzero and independent;
}
