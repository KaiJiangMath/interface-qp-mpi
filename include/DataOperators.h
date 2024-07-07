#ifndef __DataOperators_h
#define __DataOperators_h

#include "Head.h"
#include "Data.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/


/* --------------------		Sparse matrix (CSR)		--------------------*/

/**
 * \brief	Initialization of sparse matrix of <typename T> type in CSR format;
 *
 * <typename T>:	int, float, double, fftw_complex;
 *
 * \param	src:	sparse matrix of <typename T> type in CSR format;
 * \param	row:	row number of matrix 'src';
 * \param	col:	column of matrix 'src';
 * \param	nnz:	number of nonzero entries;
 */
template <typename T>
tCSRmat<T> fn_tCSRmat_init ( int row, int col, int nnz )
{
	tCSRmat<T> src;
	src.row = row;
	src.col = col;
	src.nnz = nnz;
	src.IA	= (int *) malloc( sizeof(T) * (row+1) );
	src.JA	= (int *) malloc( sizeof(T) * nnz );
	src.val = ( T *)  malloc( sizeof(T) * nnz );
	src.IA[0] = 0;
	return src;
}


/**
 * \brief	Save sparse matrix of <typename T> type in CSR format;
 *
 * <typename T>:	float, double;
 *
 */
template <typename T>
void fn_tCSRmat_save ( tCSRmat<T> src, char filename[] )
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
		fprintf(fwPath, "%+.15E\n", src.val[i]);
	fclose(fwPath);
}


/**
 * \brief	Save sparse matrix of <typename T> type in CSR format;
 *
 * <typename T>:	int, float, double;
 *
 */
template <typename T>
void fn_tCSRmat_print ( tCSRmat<T> src )
{
	/* print data; */
	printf("nrow = %d, ncol = %d, nnz = %d\n", src.row, src.col, src.nnz);
	for ( int i = 0; i < src.row; ++i )
	{
		for ( int j = src.IA[i]; j < src.IA[i+1]; j++ )
			printf("(%d,%d) = %+.10E\n", i, src.JA[j], src.val[j]);
	}
}


/**
 * \brief	Release memory of sparse matrix in CSR format;
 *
 * <typename T>:	int, float, double, fftw_complex;
 *
 */
template <typename T>
void fn_tCSRmat_free ( tCSRmat<T> src )
{
	free ( src.IA  );
	free ( src.JA  );
	free ( src.val );
}

/* --------------------		Sparse matrix (CSR)		--------------------*/


/* --------------------		Sparse matrix (CCS)		--------------------*/

/**
 * \brief	Initialization of sparse matrix of <typename T> type in CCS format;
 *
 * <typename T>:	int, float, double, fftw_complex;
 *
 * \param	src:	sparse matrix of <typename T> type in CCS format;
 * \param	row:	row number of matrix 'src';
 * \param	col:	column of matrix 'src';
 * \param	nnz:	number of nonzero entries;
 */
template <typename T>
tCCSmat<T> fn_tCCSmat_init ( int row, int col, int nnz )
{
	tCCSmat<T> src;
	src.row = row;
	src.col = col;
	src.nnz = nnz;
	src.JA	= (int *) malloc( sizeof(T) * (col+1) );
	src.IA	= (int *) malloc( sizeof(T) * nnz );
	src.val = ( T *)  malloc( sizeof(T) * nnz );
	src.JA[0] = 0;
	return src;
}


/**
 * \brief	Save sparse matrix of <typename T> type in CCS format;
 *
 * <typename T>:	float, double;
 *
 */
template <typename T>
void fn_tCCSmat_save ( tCCSmat<T> src, char filename[] )
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
		fprintf(fwPath, "%+.15E\n", src.val[i]);
	fclose(fwPath);
}


/**
 * \brief	Save sparse matrix of <typename T> type in CCS format;
 *
 * <typename T>:	int, float, double;
 *
 */
template <typename T>
void fn_tCCSmat_print ( tCCSmat<T> src )
{
	/* print data; */
	printf("nrow = %d, ncol = %d, nnz = %d\n", src.row, src.col, src.nnz);
	for ( int i = 0; i < src.col; ++i )
	{
		for ( int j = src.JA[i]; j < src.JA[i+1]; j++ )
			printf("(%d,%d) = %+.10E\n", src.IA[j], i, src.val[j]);
	}
}


/**
 * \brief	Release memory of sparse matrix in CCS format;
 *
 * <typename T>:	int, float, double, fftw_complex;
 *
 */
template <typename T>
void fn_tCCSmat_free ( tCCSmat<T> src )
{
	free ( src.JA  );
	free ( src.IA  );
	free ( src.val );
}

/* --------------------		Sparse matrix (CCS)	--------------------*/


/* --------------------		Dense matrix (vector)	--------------------*/

/**
 * \brief	Initialization of 'tvec' vector;
 *
 * <typename T>:	int, float, double, fftw_complex;
 *
 */
template <typename T>
tvec<T> fn_tvec_init ( int len )
{
	tvec<T> src;
	src.row = 1;
	src.col = 1;
	src.len = len;
	src.val = (T *) malloc( sizeof(T) * len );
	return src;
}


/**
 * \brief	Set zero of a 'tvec' vector;
 *
 * <typename T>:	int, float, double;
 *
 */
template <typename T>
void fn_tvec_setZero ( tvec<T> src )
{
	for ( int i = 0; i < src.len; i++ )
		src.val[i] = 0;
}


/**
 * \brief	Transpose a matrix whose is straightened as a 'tvec' vector;
 *
 * <typename T>:	int, float, double;
 *
 * \param	row:	the number of rows of 'src';
 *			col:	the number of columns of 'src';
 *
 * \note!	rslt.len == src.len;
 *
 */
template <typename T>
void fn_tvec_trans ( tvec<T> rslt, tvec<T> src, int row, int col )
{
	for ( int i = 0; i < row; i++ )
	{
		for ( int j = 0; j < col; j++ )
		{
			int srcInd  = i*col + j;
			int rsltInd = j*row + i;
			rslt.val[rsltInd] = src.val[srcInd];
		}
	}
}


/**
 * \brief	Result of adding two 'tvec' vectors;
 *
 * <typename T>:	int, float, double;
 *
 */
template <typename T>
void fn_tvec_add ( tvec<T> dst, tvec<T> src, T dstCoeff, T srcCoeff )
{
	for ( int i = 0; i < dst.len; i++ )
		dst.val[i] = dstCoeff*dst.val[i] + srcCoeff*src.val[i];
}


/**
 * \brief	The inner product of a 'tvec' vector;
 *
 * <typename T>:	int, float, double;
 *
 */
template <typename T>
T fn_tvec_inner ( tvec<T> lft, tvec<T> rht )
{
	T rslt = 0;
	for ( int i = 0; i < lft.len; i++ )
		rslt += lft.val[i] * rht.val[i];
	return rslt;
}


/**
 * \brief	Maximal absolute value of a 'tvec' vector;
 *
 * <typename T>:	int, float, double;
 *
 */
template <typename T>
T fn_tvec_maxAbs ( tvec<T> src )
{
	T tmp, rslt;
	rslt = 0;
	for ( int i = 0; i < src.len; i++ )
	{
		tmp = fabs(src.val[i]);
		rslt = (rslt > tmp ? rslt : tmp);
	}
	return rslt;
}


/**
 * \brief	norm-2 value of a 'tvec' vector;
 *
 * <typename T>:	int, float, double;
 *
 */
template <typename T>
T fn_tvec_norm2 ( tvec<T> src )
{
	T rslt;
	rslt = 0;
	for ( int i = 0; i < src.len; i++ )
	{
		rslt += pow(src.val[i], 2);
	}
	rslt /= src.len;
	return sqrt(rslt);
}


/**
 * \brief	Save 'tvec' vector;
 *
 * <typename T>:	float, double;
 *
 */
template <typename T>
void fn_tvec_save ( tvec<T> src, char filename[] )
{
	/* open file; */
	FILE *fwPath = fopen(filename, "w");
	if ( fwPath == NULL && myrank == 0 )
		printf("%s error.\n", filename);

	/* save data; */
//	if ( src.row == 1 && src.col == 1 )
//		fprintf(fwPath, "%d\n", src.len);
//	else
		fprintf(fwPath, "%d\t%d\t%d\n", src.row, src.col, src.len);
	for ( int i = 0; i < src.len; i++ )
		fprintf(fwPath, "%+.15E\n", src.val[i]);
	fclose(fwPath);
}


/**
 * \brief	Print 'tvec' vector;
 *
 * <typename T>:	int, float, double;
 *
 */
template <typename T>
void fn_tvec_print ( tvec<T> src )
{
	/* print data; */
	if ( src.row == 1 && src.col == 1 )
		printf("len = %d\n", src.len);
	else
		printf("nrow = %d, ncol = %d, nnz = %d\n", src.row, src.col, src.len);
	for ( int i = 0; i < src.len; i++ )
		printf("%+.10E\n", src.val[i]);
}


/**
 * \brief	Release memory of 'tvec' vector;
 *
 * <typename T>:	int, float, double, fftw_complex;
 *
 */
template <typename T>
void fn_tvec_free ( tvec<T> src )
{
	free ( src.val );
}

/* --------------------		Dense matrix (vector)	--------------------*/

/* --------------------		Dense matrix (matrix)	--------------------*/

/**
 * \brief	Initialization of 'tmat' matrix;
 *
 * <typename T>:	int, float, double, fftw_complex;
 *
 */
template <typename T>
tmat<T> fn_tmat_init ( int row, int col )
{
	tmat<T> src;
	src.row = row;
	src.col = col;
	src.ele = row * col;
	src.val = (T **) malloc( sizeof(T *) * row );
	for ( int i = 0; i < row; i++ )
		src.val[i] = (T *) malloc( sizeof(T) * col );
	return src;
}


/**
 * \brief	Set zero of 'tmat' matrix;
 *
 * <typename T>:	int, float, double;
 *
 */
template <typename T>
void fn_tmat_setZero ( tmat<T> src )
{
	for ( int i = 0; i < src.row; i++ )
		for ( int j = 0; j < src.col; j++ )
			src.val[i][j] = 0;
}


/**
 * \brief	Set identity matrix of 'tmat' matrix;
 *
 * <typename T>:	int, float, double;
 *
 */
template <typename T>
void fn_tmat_setIdentity ( tmat<T> src )
{
	if ( src.row == src.col )
	{
		for ( int i = 0; i < src.row; i++ )
		{
			for ( int j = 0; j < src.col; j++ )
			{
				if ( i == j )
					src.val[i][j] = 1;
				else
					src.val[i][j] = 0;
			}
		}
	}
	else
		printf("Error using 'fn_tmat_setIdentity'.\n");
}


/**
 * \brief	Copy 'tmat' matrix;
 *
 * <typename T>:	int, float, double;
 *
 */
template <typename T>
void fn_tmat_copy ( tmat<T> dst, tmat<T> src )
{
	for ( int i = 0; i < src.row; i++ )
		for ( int j = 0; j < src.col; j++ )
			dst.val[i][j] = src.val[i][j];
}


/**
 * \brief	Transpose 'tmat' matrix;
 *
 * <typename T>:	int, float, double;
 *
 */
template <typename T>
tmat<T> fn_tmat_trans ( tmat<T> src )
{
	tmat<T> dst = fn_tmat_init<T> ( src.col, src.row );
	for ( int i = 0; i < src.row; i++ )
		for ( int j = 0; j < src.col; j++ )
			dst.val[j][i] = src.val[i][j];
	return dst;
}


/**
 * \brief	Result of adding two 'tmat' matrices;
 *
 * <typename T>:	int, float, double;
 *
 */
template <typename T>
void fn_tmat_add ( tmat<T> dst, tmat<T> src, T dstCoeff, T srcCoeff )
{
	for ( int i = 0; i < dst.row; i++ )
		for ( int j = 0; j < dst.col; j++ )
			dst.val[i][j] = dstCoeff*dst.val[i][j] + srcCoeff*src.val[i][j];
}


/**
 * \brief	matrix 'lhs' multiplies matrix 'rhs'; result is lrhs;
 *
 * <typename T>:	int, float, double;
 *
 * \note!	lhs.col must be equal to rhs.row;
 *
 */
template <typename T>
void fn_tmat_multiply ( tmat<T> lhs, tmat<T> rhs, tmat<T> lrhs )
{
	for ( int i = 0; i < lhs.row; i++ )
	{
		for ( int j = 0; j < rhs.col; j++ )
		{
			lrhs.val[i][j] = 0;
			for ( int k = 0; k < lhs.col; k++ )
				lrhs.val[i][j] += lhs.val[i][k] * rhs.val[k][j];
		}
	}
}


/**
 * \brief	The inner product of two 'tmat' matrices;
 *
 * <typename T>:	int, float, double;
 *
 */
template <typename T>
T fn_tmat_inner ( tmat<T> lhs, tmat<T> rhs )
{
	T rslt = 0;
	for ( int i = 0; i < lhs.row; i++ )
		for ( int j = 0; j < rhs.col; j++ )
			rslt += lhs.val[i][j] * rhs.val[i][j];
	return rslt;
}


/**
 * \brief	Maximal absolute value of a 'tmat' matrix;
 *
 * <typename T>:	int, float, double;
 *
 */
template <typename T>
T fn_tmat_maxAbs ( tmat<T> src )
{
	T tmp, rslt;
	rslt = 0;
	for ( int i = 0; i < src.row; i++ )
	{
		for ( int j = 0; j < src.col; j++ )
		{
			tmp = fabs(src.val[i][j]);
			rslt = ( rslt > tmp ? rslt : tmp );
		}
	}
	return rslt;
}


/**
 * \brief	norm-2 value of a 'tmat' matrix;
 *
 * <typename T>:	int, float, double;
 *
 */
template <typename T>
T fn_tmat_norm2 ( tmat<T> src )
{
	T rslt;
	rslt = 0;
	for ( int i = 0; i < src.row; i++ )
	{
		for ( int j = 0; j < src.col; j++ )
		{
			rslt += pow(src.val[i][j], 2);
		}
	}
	rslt /= (src.row * src.col);
	return sqrt(rslt);
}


/**
 * \brief	Save 'tmat' matrix;
 *
 * <typename T>:	float, double;
 *
 */
template <typename T>
void fn_tmat_save ( tmat<T> src, char filename[] )
{
	/* open file; */
	FILE *fwPath = fopen(filename, "w");
	if ( fwPath == NULL && myrank == 0 )
		printf("%s error.\n", filename);

	/* save data; */
	fprintf(fwPath, "%d\t%d\t%d\n", src.row, src.col, src.ele);
	for ( int i = 0; i < src.row; i++ )
		for ( int j = 0; j < src.col; j++ )
			fprintf(fwPath, "%+.15E\n", src.val[i][j]);
	fclose(fwPath);
}


/**
 * \brief	Print 'tmat' matrix;
 *
 * <typename T>:	float, double;
 *
 */
template <typename T>
void fn_tmat_print ( tmat<T> src )
{
	/* print data; */
	printf("nrow = %d, ncol = %d, nnz = %d\n", src.row, src.col, src.ele);
	for ( int i = 0; i < src.row; i++ )
	{
		for ( int j = 0; j < src.col; j++ )
			printf("%+.10E\t", src.val[i][j]);
		printf("\n");
	}
}


/**
 * \brief	Release memory of 'tmat' matrix;
 *
 * <typename T>:	int, float, double, fftw_complex;
 *
 */
template <typename T>
void fn_tmat_free ( tmat<T> src )
{
	for ( int i = 0; i < src.row; i++ )
		free ( src.val[i] );
	free ( src.val );
}


/* --------------------		Dense matrix (matrix)	--------------------*/

#endif
