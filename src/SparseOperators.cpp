/*! \file	SparseOperators.cpp
 *
 * \brief	Basic operators of sparse matrix (CCS);
 */

#include "Data.h"
#include "Head.h"
#include "DataOperators.h"
#include "functs.h"
#include "umfpack.h"


/**
 * \brief	Save sparse matrix of int type in CSR format;
 *
 */
void fn_tCSRmat_save_int ( tCSRmat<int> src, char filename[] )
{
	/* open file; */
	FILE *fwPath = fopen(filename, "w");
	if ( fwPath == NULL && myrank == 0 )
		printf("%s error.\n", filename);

	/* save IA and JA */
	fprintf(fwPath, "%d\t%d\t%d\n", src.row, src.col, src.nnz);
	for ( int i = 0; i <= src.row; i++ )
		fprintf(fwPath, "%d\n", src.IA[i]);
	for ( int i = 0; i < src.nnz; i++ )
		fprintf(fwPath, "%d\n", src.JA[i]);

	/* save data; */
	for ( int i = 0; i < src.nnz; i++ )
		fprintf(fwPath, "%+d\n", src.val[i]);
	fclose(fwPath);
}


/**
 * \brief	Save sparse matrix of int type in CCS format;
 *
 */
void fn_tCCSmat_save_int ( tCCSmat<int> src, char filename[] )
{
	/* open file; */
	FILE *fwPath = fopen(filename, "w");
	if ( fwPath == NULL && myrank == 0 )
		printf("%s error.\n", filename);

	/* save IA and JA */
	fprintf(fwPath, "%d\t%d\t%d\n", src.row, src.col, src.nnz);
	for ( int i = 0; i <= src.col; i++ )
		fprintf(fwPath, "%d\n", src.JA[i]);
	for ( int i = 0; i < src.nnz; i++ )
		fprintf(fwPath, "%d\n", src.IA[i]);

	/* save data; */
	for ( int i = 0; i < src.nnz; i++ )
		fprintf(fwPath, "%+d\n", src.val[i]);
	fclose(fwPath);
}


/**
 * \brief	Generate block diagonal sparse matrix;
 *
 * \param	ele:		the element on the block diagonal line;
 * \param	d:			the distance between two adjacent block diagonal lines;
 *
 * \return	rslt:		the block diagonal sparse matrix;
 *
 * \note!	Only consider a square matrix;
 */
tCCSmat<double> fn_block_diag_dCCSmat	( int row, int col,	int	d )
{
	/* initialization; */
	int blockNum = row / d;
	tCCSmat<double> rslt = fn_tCCSmat_init<double> ( row, col, blockNum*col );

	/* generate CCS matrix; */
	int nnz = 0;
	for ( int i = 0; i < col; ++i )
	{
		int j = i%d;
		while ( j < row )
		{
			rslt.IA[nnz]  = j;
			rslt.val[nnz] = 1.0;
			j += d;
			nnz ++ ;
		}
		rslt.JA[i+1] = nnz;
	}
	rslt.nnz = nnz;

	return rslt;
}


/**
 * \brief	Calculate the tensor product of diag(diagSrc) and dccsSrc;
 *
 * \param	diagSrc:	a dense vector;
 * \param	dccsSrc:	a structure body with double CCS type;
 * \param	order:		a parameter for 'diagSrc';
 *
 * \return	a structure body with doulbe CCS type;
 *			the return matrix is
 *			diag( diagSrc[0]*dccsSrc, ..., diagSrc[-1]*dccsSrc );
 */
tCCSmat<double> fn_tensor_diag_dCCSmat	(	tvec<double>		diagSrc,
											tCCSmat<double>		dccsSrc,
												int				order )
{
	/* Initialization of result; */
	int row = diagSrc.len * dccsSrc.row;
	int col = diagSrc.len * dccsSrc.col;
	int nnz = diagSrc.len * dccsSrc.nnz;	// big enough;
	tCCSmat<double> tensorRslt = fn_tCCSmat_init<double> ( row, col, nnz );

	double eleDiagTmp, eleTensorTmp;
	int indJA = 0, indnnz = 0;
	int startJA, endJA;

	/* calculate the tensor product; */
	for ( int i = 0; i < diagSrc.len; i++ )
	{
		eleDiagTmp = pow(diagSrc.val[i], order);
		/* each block is 'dccsSrc'; */
		for ( int j0 = 0; j0 < dccsSrc.col; j0++ )
		{
			startJA = dccsSrc.JA[j0];
			endJA	= dccsSrc.JA[j0+1];
			for ( int j1 = startJA; j1 < endJA; j1++ )
			{
				eleTensorTmp = eleDiagTmp * dccsSrc.val[j1];
				if ( fabs(eleTensorTmp) > ZEROTOL )
				{
					/* matrix block along the diagonal line, so we add i*dccsSrc.row; */
					tensorRslt.IA[indnnz] = i*dccsSrc.row + dccsSrc.IA[j1];
					tensorRslt.val[indnnz] = eleTensorTmp;
					indnnz ++ ;
				}
			}
			indJA ++ ;
			tensorRslt.JA[indJA] = indnnz;
		}
	}
	tensorRslt.nnz = indnnz;

	return tensorRslt;
}


/**
 * \brief	Calculate the product of dccsSrc * diagSrc;
 *
 * \param	diagSrc:	a diagonal CCS matrix;
 * \param	dccsSrc:	a structure body with double CCS type;
 *
 * \return	a structure body with doulbe CCS type;
 */
tCCSmat<double> fn_multiply_diag_dCCSmat(	tCCSmat<double>		diagSrc,
											tCCSmat<double>		dccsSrc )
{
	/* Initialization of result; */
	int row = dccsSrc.row;
	int col = dccsSrc.col;
	int nnz = row * col;
	tCCSmat<double> rsltTmp = fn_tCCSmat_init<double> ( row, col, nnz );
	rsltTmp.JA[0] = 0;

	/* calculate result; */
	int indnnz = 0;
	for ( int i = 0; i < col; i ++ )
	{
		/* there is a nonzero element in 'diagSrc'; */
		if ( diagSrc.JA[i] < diagSrc.JA[i+1] )
		{
			double coeff = diagSrc.val[ diagSrc.JA[i] ];
			for ( int j = dccsSrc.JA[i]; j < dccsSrc.JA[i+1]; j++ )
			{
				rsltTmp.IA[indnnz]  = dccsSrc.IA[j];
				rsltTmp.val[indnnz] = coeff * dccsSrc.val[j];
				indnnz ++ ;
			}
		}
		rsltTmp.JA[i+1] = indnnz;
	}

	/* return result; */
	tCCSmat<double> rslt = fn_tCCSmat_init<double> ( row, col, indnnz );
	memcpy ( rslt.JA, rsltTmp.JA, sizeof(int) * (col+1) );
	for ( int i = 0; i < indnnz; i ++ )
	{
		rslt.IA[i]  = rsltTmp.IA[i];
		rslt.val[i] = rsltTmp.val[i];
	}
	
	/* release memory; */
	fn_tCCSmat_free<double> ( rsltTmp );

	return rslt;
}



/**
 * \brief	Calculate the result of coeff * dccsSrc;
 *
 * \param	dccsSrc:	a structure body with double CCS type;
 * \param	coeff:		a constant real value;
 */
tCCSmat<double> fn_const_multiply_dCCSmat	(	tCCSmat<double>		dccsSrc,
													double			coeff )
{
	/* Initialization of result; */
	int row = dccsSrc.row;
	int col = dccsSrc.col;
	int nnz = dccsSrc.nnz;
	tCCSmat<double> rslt = fn_tCCSmat_init<double> ( row, col, nnz );

	/* coeff * dccsSrc; */
	memcpy ( rslt.JA, dccsSrc.JA, sizeof(int) * (col+1) );
	if ( fabs(coeff) > ZEROTOL )
	{
		memcpy ( rslt.IA, dccsSrc.IA, sizeof(int) * nnz );
		for ( int i = 0; i < nnz; i++ )
			rslt.val[i] = coeff * dccsSrc.val[i];
	}

	/*
	double tmp;
	int indnnz = 0;
	for ( int i = 0; i < dccsSrc.col; ++i )
	{
		for ( int j = dccsSrc.JA[i]; j < dccsSrc.JA[i+1]; j++ )
		{
			tmp = coeff * dccsSrc.val[j];
			if ( fabs(tmp) > ZEROTOL )
			{
				rslt.val[j] = tmp;
				rslt.IA[j]  = dccsSrc.IA[j];
				indnnz ++ ;
			}
		}
		rslt.JA[i+1] = indnnz;
	}
	rslt.nnz = indnnz;
	*/

	return rslt;
}


/**
 * \brief	Obtain the diagonal elements of dccsSrc and add coeff*I on them;
 *
 * \param	dccsSrc:	a structure body with double CCS type;
 * \param	coeff1:		a constant real value on 'dccsSrc';
 * \param	coeff2:		a constant real value on the identity matrix;
 *
 * \return	rslt:		a diagonal sparse matrix;
 */
tCCSmat<double> fn_obt_diag_add_dCCSmat	(	tCCSmat<double>		dccsSrc,
												double			coeff1,
												double			coeff2 )
{
	/* Initialization of result; */
	int len = dccsSrc.row < dccsSrc.col ? dccsSrc.row : dccsSrc.col;
	tCCSmat<double> rslt = fn_tCCSmat_init<double> ( 
			dccsSrc.row, dccsSrc.col, dccsSrc.nnz + len );

	double tmp;
	int	 indnnz = 0;
	bool isDiag = false;
	for ( int i = 0; i < dccsSrc.col; ++i )
	{
		isDiag = false;
		for ( int j = dccsSrc.JA[i]; j < dccsSrc.JA[i+1]; j++ )
		{
			if ( dccsSrc.IA[j] == i )	// diagonal part;
			{
				isDiag = true;
				tmp = coeff1 * dccsSrc.val[j] + coeff2;
				if ( fabs(tmp) > ZEROTOL )
				{
					rslt.val[indnnz] = tmp;
					rslt.IA[indnnz]  = dccsSrc.IA[j];
					indnnz ++ ;
				}
			}
		}
		if ( !isDiag && i < len )
		{
			if ( fabs(coeff2) > ZEROTOL )
			{
				rslt.val[indnnz] = coeff2;
				rslt.IA[indnnz]  = i;
				indnnz ++ ;
			}
		}
		rslt.JA[i+1] = indnnz;
	}
	rslt.nnz = indnnz;

	return rslt;
}


/**
 * \brief	Calculate the result of dccsSrc + coeff*I;
 *
 * \param	dccsSrc:	a structure body with double CCS type;
 * \param	coeff1:		a constant real value on 'dccsSrc';
 * \param	coeff2:		a constant real value on the identity matrix;
 *
 * \return	rslt:		only the diagonal elements are changed;
 */
tCCSmat<double> fn_diag_add_dCCSmat		(	tCCSmat<double>		dccsSrc,
												double			coeff1,
												double			coeff2 )
{
	/* Initialization of result; */
	int len = dccsSrc.row < dccsSrc.col ? dccsSrc.row : dccsSrc.col;
	tCCSmat<double> rslt = fn_tCCSmat_init<double> ( 
			dccsSrc.row, dccsSrc.col, dccsSrc.nnz + len );

	double tmp;
	int	 indnnz = 0;
	bool isDiag = false;
	for ( int i = 0; i < dccsSrc.col; ++i )
	{
		isDiag = false;
		for ( int j = dccsSrc.JA[i]; j < dccsSrc.JA[i+1]; j++ )
		{
			if ( dccsSrc.IA[j] == i )	// diagonal part;
			{
				isDiag = true;
				tmp = coeff1 * dccsSrc.val[j] + coeff2;
			}
			else	// not diagonal part;
			{
				tmp = coeff1 * dccsSrc.val[j];
			}
			if ( fabs(tmp) > ZEROTOL )
			{
				rslt.val[indnnz] = tmp;
				rslt.IA[indnnz]  = dccsSrc.IA[j];
				indnnz ++ ;
			}
		}
		if ( !isDiag && i < len )
		{
			if ( fabs(coeff2) > ZEROTOL )
			{
				rslt.val[indnnz] = coeff2;
				rslt.IA[indnnz]  = i;
				indnnz ++ ;
			}
		}
		rslt.JA[i+1] = indnnz;
	}
	rslt.nnz = indnnz;

	return rslt;
}


/**
 * \brief	Calculate the result of dccsSrc + coeff*blockI;
 *
 * \param	dccsSrc:	a structure body with double CCS type;
 * \param	coeff1:		a constant real value on 'dccsSrc';
 * \param	coeff2:		a constant real value on the identity matrix;
 * \param	d:			the distance between two adjacent block diagonal lines;
 */
tCCSmat<double> fn_block_diag_add_dCCSmat (	tCCSmat<double>		dccsSrc,
												double			coeff1,
												double			coeff2,
												int				d )
{
	/* Initialization of result; */
	int len = dccsSrc.row < dccsSrc.col ? dccsSrc.row : dccsSrc.col;
	tCCSmat<double> rslt = fn_tCCSmat_init<double> ( 
			dccsSrc.row, dccsSrc.col, dccsSrc.nnz + len );

	double tmp;
	int	 indnnz = 0;
	bool isDiag = false;
	for ( int i = 0; i < dccsSrc.col; ++i )
	{
		isDiag = false;
		for ( int j = dccsSrc.JA[i]; j < dccsSrc.JA[i+1]; j++ )
		{
			int check = fabs( dccsSrc.IA[j] - i );
			if ( check%d == 0 )	// diagonal part;
			{
				isDiag = true;
				tmp = coeff1 * dccsSrc.val[j] + coeff2;
			}
			else	// not diagonal part;
			{
				tmp = coeff1 * dccsSrc.val[j];
			}
			if ( fabs(tmp) > ZEROTOL )
			{
				rslt.val[indnnz] = tmp;
				rslt.IA[indnnz]  = dccsSrc.IA[j];
				indnnz ++ ;
			}
		}
		if ( !isDiag && i < len )
		{
			if ( fabs(coeff2) > ZEROTOL )
			{
				rslt.val[indnnz] = coeff2;
				rslt.IA[indnnz]  = i;
				indnnz ++ ;
			}
		}
		rslt.JA[i+1] = indnnz;
	}
	rslt.nnz = indnnz;

	return rslt;
}


/**
 * \brief	Transpose the CCS 'dccsSrc';
 *
 * \param	dccsSrc:	a structure body with double CCS type;
 */
tCCSmat<double> fn_trans_dCCSmat ( tCCSmat<double> dccsSrc )
{
	/* memory allocation of results; */
	tCCSmat<double> transRslt = fn_tCCSmat_init<double> ( 
			dccsSrc.col, dccsSrc.row, dccsSrc.nnz );

	/* calculate 'transRslt'; */
	int nnz = 0;
	for ( int i = 0; i < transRslt.col; i++ )
	{
		for ( int j = 0; j < dccsSrc.nnz; j++ )
		{
			if ( dccsSrc.IA[j] == i )
			{
				transRslt.val[nnz] = dccsSrc.val[j]; // corresponding values;
				/* search 'transRslt.IA'; */
				for ( int j0 = 1; j0 <= dccsSrc.col; j0++ )
				{
					if ( dccsSrc.JA[j0] > j )
					{
						transRslt.IA[nnz] = j0 - 1;
						break;
					}
				}
				nnz ++ ;
			}
		}
		transRslt.JA[i+1] = nnz;	// nonzero elements on each row;
	}

	return transRslt;
}


/**
 * \brief	Transpose the CCS 'dccsSrc' by a faster way;
 *
 * \param	dccsSrc:	a structure body with double CCS type;
 */
tCCSmat<double> fn_fast_trans_dCCSmat ( tCCSmat<double> dccsSrc )
{
	/* memory allocation of results; */
	tCCSmat<double> transRslt = fn_tCCSmat_init<double> ( 
			dccsSrc.col, dccsSrc.row, dccsSrc.nnz );

	/* obtain the row index, col index, and value of all nonzero elements; */
	tvec<int> rowInd	= fn_tvec_init<int>		( dccsSrc.nnz );
	tvec<int> colInd	= fn_tvec_init<int>		( dccsSrc.nnz );
	tvec<double> values = fn_tvec_init<double>	( dccsSrc.nnz );
	int nnz = 0;
	for ( int i = 0; i < dccsSrc.col; ++i )
	{
		for ( int j = dccsSrc.JA[i]; j < dccsSrc.JA[i+1]; j++ )
		{
			rowInd.val[nnz] = dccsSrc.IA[j];
			colInd.val[nnz] = i;
			values.val[nnz] = dccsSrc.val[j];
			nnz ++ ;
		}
	}

	/* sort 'rowInd'; */
	double tmp;
	vector <stu_sort> sort_array ( dccsSrc.nnz );
	for ( int i = 0; i < dccsSrc.nnz; ++i )
	{
		sort_array[i].ind = i;
		sort_array[i].val = rowInd.val[i];
	}
	sort(sort_array.begin(), sort_array.end(), fn_compare);

	/* transpose 'dccsSrc'; */
	int sortInd, indOld, indNew;
	indOld = rowInd.val[ sort_array[0].ind ];
	nnz = 0;
	for ( int i = 0; i < transRslt.col+1; i++ )
		transRslt.JA[i] = 0;
	for ( int i = 0; i < dccsSrc.nnz; i++ )
	{
		sortInd = sort_array[i].ind;
		indNew = rowInd.val[ sortInd ];
		if ( i > 0 && indNew != indOld )
		{
			transRslt.JA[indNew] = nnz;
		}
		indOld = indNew;
		transRslt.IA[i]  = colInd.val[ sortInd ];
		transRslt.val[i] = values.val[ sortInd ];
		nnz ++ ;
	}
	transRslt.JA[ transRslt.col ] = nnz;

	/* release memory; */
	fn_tvec_free<int>	 ( rowInd );
	fn_tvec_free<int>	 ( colInd );
	fn_tvec_free<double> ( values );

	return transRslt;
}


/**
 * \brief	Calculate the result of dccsSrc * rhoJCplx;
 *			considering 'dccsSrc' is CCS, we first transpose
 *			'dccsSrc' and use its each compressed column to 
 *			multiply 'rhoJCplx';
 *			i.e., dccsSrcTrans' * rhoJCplx;
 *
 * \param	dccsSrc:	a structure body with double CCS type;
 * \param	rhoJCplx:	a vector with complex type;
 * 
 * \note!	'dccsSrc.row' == 'rhoJCplx.row';
 */
void fn_cvec_multiply_dCCSmat	(	tCCSmat<double>			dccsSrcTrans,
									tvec<fftw_complex>		rhoJCplx,
									tvec<fftw_complex>		rslt )
{
	/* calculate dccsSrc * rhoJCplx by transposing 'dccsSrc'; */
	for ( int i = 0; i < dccsSrcTrans.col; ++i )
	{
		fn_complex_setZero ( rslt.val[i] );
		for ( int j = dccsSrcTrans.JA[i]; j < dccsSrcTrans.JA[i+1]; j++ )
		{
			rslt.val[i][0] += dccsSrcTrans.val[j] * 
								rhoJCplx.val[ dccsSrcTrans.IA[j] ][0];
			rslt.val[i][1] += dccsSrcTrans.val[j] * 
								rhoJCplx.val[ dccsSrcTrans.IA[j] ][1];
		}
	}
}


/**
 * \brief	Calculate the result of dccsSrc1 add dccsSrc2;
 *			coeff1 * dccsSrc1 + coeff2 * dccsSrc2;
 *
 * \param	dccsSrc1:	a structure body with double CCS type;
 * \param	dccsSrc2:	a structure body with double CCS type;
 * \param	coeff1:		the coefficient before 'dccsSrc1';
 * \param	coeff2:		the coefficient before 'dccsSrc2';
 *
 * \note!	'dccsSrc1' and 'dccsSrc2' must have the same number of 
 *			rows and columns;
 *
 * \return	a structure body with double CCS type;
 */
tCCSmat<double> fn_add_dCCSmat		(	tCCSmat<double>		dccsSrc1,
										tCCSmat<double>		dccsSrc2,
											double			coeff1,
											double			coeff2 )
{
	/**
	 * dccsSrc1.row must be equal to dccsSrc2.row;
	 * dccsSrc1.col must be equal to dccsSrc2.col;
	 */
	tCCSmat<double> sumRslt = fn_tCCSmat_init<double> ( 
			dccsSrc1.row, dccsSrc1.col, dccsSrc1.nnz + dccsSrc2.nnz );

	double sumTmp;
	int indJA = 0, indnnz = 0;
	int indIA1, endJA1;
	int indIA2, endJA2;

	/* calculate the result of 'dccsSrc1' add 'dccsSrc2'; */
	for ( int j0 = 0; j0 < sumRslt.col; j0++ )
	{
		int j1 = dccsSrc1.JA[j0];
		int j2 = dccsSrc2.JA[j0];
		endJA1	= dccsSrc1.JA[j0+1];
		endJA2	= dccsSrc2.JA[j0+1];
		while ( j1 < endJA1 || j2 < endJA2 )
		{
			/* index for left matrix; */
			if ( j1 == endJA1 )
				indIA1 = dccsSrc1.row;
			else
				indIA1 = dccsSrc1.IA[j1];

			/* index for right matrix; */
			if ( j2 == endJA2 )
				indIA2 = dccsSrc2.row;
			else
				indIA2 = dccsSrc2.IA[j2];

			if ( indIA1 < indIA2 )			// zero in 'dccsSrc2';
			{
				sumRslt.IA[indnnz] = dccsSrc1.IA[j1];
				sumRslt.val[indnnz] = coeff1 * dccsSrc1.val[j1];
				indnnz ++ ;
				j1 ++ ;
			}
			else if ( indIA1 == indIA2 )	// two elements both are nonzero;
			{
				sumTmp = coeff1 * dccsSrc1.val[j1] + coeff2 * dccsSrc2.val[j2];
				if ( fabs(sumTmp) > ZEROTOL )
				{
					sumRslt.IA[indnnz] = dccsSrc1.IA[j1];	// or dccsSrc2.IA[j2];
					sumRslt.val[indnnz] = sumTmp;
					indnnz ++ ;
				}
				j1 ++ ;		j2 ++ ;
			}
			else if ( indIA1 > indIA2 )		// zero in 'dccsSrc1';
			{
				sumRslt.IA[indnnz] = dccsSrc2.IA[j2];
				sumRslt.val[indnnz] = coeff2 * dccsSrc2.val[j2];
				indnnz ++ ;
				j2 ++ ;
			}
		}
		indJA ++ ;
		sumRslt.JA[indJA] = indnnz;
	}
	sumRslt.nnz = indnnz;

	return sumRslt;
}


/**
 * \brief	Solve Ax = b by umfpack.h; (double type);
 *
 * \param	A:		a sparse matrix whose format is double CCS type;
 * \param	b:		a dense matrix with double type;
 *
 */
int fn_umfpack_solver	(	tCCSmat<double>		A, 
							tvec<double>		b, 
							tvec<double>		x )
{
	if ( A.row != A.col )
	{
		if ( myrank == 0 )
		{
			printf("Error using 'fn_umfpack_solver'\n");
			printf("Coefficient matrix must be square.\n");
		}
		return 1;
	}

	/* initialization; */
	int status;
	double Control [UMFPACK_CONTROL], Info [UMFPACK_INFO];
	umfpack_di_defaults (Control) ;
//	Control [UMFPACK_IRSTEP] = 5;
//	Control [UMFPACK_PIVOT_TOLERANCE] = 1e-2;

	/* symbols; */
	void *Symbolic, *Numeric ;
	status = umfpack_di_symbolic (A.row, A.col, A.JA, A.IA, A.val, 
									&Symbolic, Control, Info) ;

	/* numerical; */
	status = umfpack_di_numeric (A.JA, A.IA, A.val, Symbolic, 
									&Numeric, Control, Info) ;
	umfpack_di_free_symbolic (&Symbolic) ;

	/* solve Ax = b; */
	status = umfpack_di_solve (UMFPACK_A, A.JA, A.IA, A.val, 
									x.val, b.val, Numeric, Control, Info) ;
	umfpack_di_free_numeric (&Numeric) ;
	return 0;
}


/**
 * \brief	Solve Ax = b by umfpack.h; ( complex type );
 *
 * \param	A:		a sparse matrix whose format is double CCS type;
 * \param	b:		a dense matrix with fftw_complex type;
 *
 */
int fn_umfpack_complex_solver	(	tCCSmat<double>			A, 
									tvec<fftw_complex>		b, 
									tvec<fftw_complex>		x )
{
	if ( A.row != A.col )
	{
		if ( myrank == 0 )
		{
			printf("Error using 'fn_umfpack_solver'\n");
			printf("Coefficient matrix must be square.\n");
		}
		return 1;
	}

	/* initialization; */
	int status;
	double Control [UMFPACK_CONTROL], Info [UMFPACK_INFO];
	umfpack_di_defaults (Control) ;
//	Control [UMFPACK_IRSTEP] = 5;
//	Control [UMFPACK_PIVOT_TOLERANCE] = 1e-2;

	/* symbols; */
	void *Symbolic, *Numeric ;
	status = umfpack_di_symbolic (A.row, A.col, A.JA, A.IA, A.val, 
									&Symbolic, Control, Info) ;

	/* numerical; */
	status = umfpack_di_numeric (A.JA, A.IA, A.val, Symbolic, 
									&Numeric, Control, Info) ;
	umfpack_di_free_symbolic (&Symbolic) ;

	/* memory allocation; */
	tvec<double> xReal = fn_tvec_init<double> ( x.len );
	tvec<double> xCplx = fn_tvec_init<double> ( x.len );
	tvec<double> bReal = fn_tvec_init<double> ( b.len );
	tvec<double> bCplx = fn_tvec_init<double> ( b.len );
	for ( int i = 0; i < b.len; i++ )
	{
		bReal.val[i] = b.val[i][0];
		bCplx.val[i] = b.val[i][1];
	}

	/* solve Ax = b; */
	status = umfpack_di_solve (UMFPACK_A, A.JA, A.IA, A.val, 
									xReal.val, bReal.val, Numeric, Control, Info) ;
	status = umfpack_di_solve (UMFPACK_A, A.JA, A.IA, A.val, 
									xCplx.val, bCplx.val, Numeric, Control, Info) ;

	/* save and release memory; */
	for ( int i = 0; i < x.len; i++ )
	{
		x.val[i][0] = xReal.val[i];
		x.val[i][1] = xCplx.val[i];
	}
	fn_tvec_free<double> ( xReal );
	fn_tvec_free<double> ( xCplx );
	fn_tvec_free<double> ( bReal );
	fn_tvec_free<double> ( bCplx );

	umfpack_di_free_numeric (&Numeric) ;
	return 0;
}
