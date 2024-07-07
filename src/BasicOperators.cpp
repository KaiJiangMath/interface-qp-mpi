/*! \file	BasicOperators.cpp
 *
 *  \brief	Basic operators;
 *
 */


#include "Head.h"
#include "Data.h"
#include "DataOperators.h"
#include "functs.h"


/**
 * \brief	Maximal absolute value;
 */
double normRealInfty ( double *src, int n )
{
    double tmp;
    double rslt = 0.0;
    for (int i = 0; i < n; i++)
    {
        tmp = fabs(src[i]);
        rslt = (rslt > tmp ? rslt : tmp);
    }
    return rslt;
}


/**
 * \brief	Maximal absolute value of an integer matrix;
 */
int fn_max_abs ( int **src, int row, int col )
{
	int tmp, rslt = 0;
	for ( int i = 0; i < row; i++ )
	{
		for ( int j = 0; j < col; j++ )
		{
			tmp = fabs(src[i][j]);
			rslt = ( rslt > tmp ? rslt : tmp );
		}
	}
	return rslt;
}


/**
 * \brief	Define the sort rule;
 */
bool fn_compare (stu_sort a, stu_sort b)
{
	return a.val < b.val;	// increase;
}

/**
 * \brief	Define the descend sort rule;
 */
bool fn_compare_descend (stu_sort a, stu_sort b)
{
	return a.val > b.val;	// decrease;
}

/* --------------------		Value (complex)	--------------------*/

/**
 * \brief	Set zero for a complex value;
 */
void fn_complex_setZero ( fftw_complex rslt )
{
	rslt[0] = 0.0; rslt[1] = 0.0;
}


/**
 * \brief	Absolute value for a complex value;
 */
double fn_complex_abs ( fftw_complex src )
{
	double rslt;
	rslt = src[0]*src[0] + src[1]*src[1];
	return sqrt( rslt );
}


/**
 * \brief	complex F1 multiplies complex F2;
 */
void fn_complex_multiply	(	fftw_complex	F1, 
								fftw_complex	F2,
								fftw_complex	rslt )
{
	// (a+ib)(c+id) = (ac-bd) + i(ad+bc);
	rslt[0] = F1[0]*F2[0] - F1[1]*F2[1];
	rslt[1] = F1[0]*F2[1] + F1[1]*F2[0];
}


/**
 * \brief	complex F1 divides complex F2;
 */
void fn_complex_divide		(	fftw_complex	F1, 
								fftw_complex	F2,
								fftw_complex	rslt )
{
	// (a+ib)/(c+id) = (ac+bd)/(c^2+d^2) + i(bc-ad)/(c^2+d^2);
	double F2norm = F2[0]*F2[0] + F2[1]*F2[1];
	rslt[0] = ( F1[0]*F2[0] + F1[1]*F2[1] ) / F2norm;
	rslt[1] = ( F1[1]*F2[0] - F1[0]*F2[1] ) / F2norm;
}

/* --------------------		Value (complex)	--------------------*/

/* --------------------		tvec (int)	--------------------*/

/**
 * \brief	Save 'tvec' vector;
 */
void fn_tvec_save_int ( tvec<int> src, char filename[] )
{
	/* open file; */
	FILE *fwPath = fopen(filename, "w");
	if ( fwPath == NULL && myrank == 0 )
		printf("%s error.\n", filename);

	/* save data; */
	if ( src.row == 1 && src.col == 1 )
		fprintf(fwPath, "%d\n", src.len);
	else
		fprintf(fwPath, "%d\t%d\t%d\n", src.row, src.col, src.len);
	for ( int i = 0; i < src.len; i++ )
		fprintf(fwPath, "%+d\n", src.val[i]);
	fclose(fwPath);
}


/**
 * \brief	Print 'tvec' vector;
 */
void fn_tvec_print_int ( tvec<int> src )
{
	/* print data; */
	if ( myrank == 0 )
	{
		if ( src.row == 1 && src.col == 1 )
			printf("len = %d\n", src.len);
		else
			printf("nrow = %d, ncol = %d, nnz = %d\n", src.row, src.col, src.len);
		for ( int i = 0; i < src.len; i++ )
			printf("%+d\n", src.val[i]);
	}
}


/* --------------------		tvec (int)	--------------------*/

/* --------------------		tvec (complex)	--------------------*/

/**
 * \brief	Set zero of a 'tvec' vector;
 */
void fn_tvec_setZero_complex ( tvec<fftw_complex> src )
{
	for ( int i = 0; i < src.len; i++ )
		src.val[i][0] = 0, src.val[i][1] = 0;
}


/**
 * \brief	Transpose a matrix whose is straightened as a 'tvec' vector;
 *
 * \param	row:	the number of rows of 'src';
 *			col:	the number of columns of 'src';
 
 * \note!	rslt.len == src.len;
 */
void fn_tvec_trans_complex (	tvec<fftw_complex>		rslt,
								tvec<fftw_complex>		src,
										int				row,
										int				col )
{
	for ( int i = 0; i < row; i++ )
	{
		for ( int j = 0; j < col; j++ )
		{
			int rsltInd = i*col + j;
			int srcInd	= j*row + i;
			rslt.val[rsltInd][0] = src.val[srcInd][0];
			rslt.val[rsltInd][1] = src.val[srcInd][1];
		}
	}
}


/**
 * \brief	Maximal absolute value of a complex 'tvec' vector;
 */
double fn_tvec_maxAbs_complex			( tvec<fftw_complex> src )
{
	double tmp, rsltAll;
	double rslt = 0.0;
	for ( int i = 0; i < src.len; i++ )
	{
		tmp = fn_complex_abs ( src.val[i] );
		rslt = (rslt > tmp ? rslt : tmp);
	}
	MPI_Barrier ( MPI_COMM_WORLD );
	MPI_Allreduce ( &rslt, &rsltAll, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );
	return rsltAll;
}


/**
 * \brief	Norm 2 of a complex 'tvec' vector;
 */
double fn_tvec_norm_complex				( tvec<fftw_complex>	src )
{
	double rslt = 0.0;
	double rsltAll;
	for ( int i = 0; i < src.len; i++ )
	{
		rslt += pow(src.val[i][0],2) + pow(src.val[i][1],2);
	}
	MPI_Barrier ( MPI_COMM_WORLD );
	MPI_Allreduce ( &rslt, &rsltAll, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
	return sqrt(rsltAll);
}


/**
 * \brief	The adding result of a complex 'tvec' vector and a double value;
 */
void fn_tvec_constAdd_complex			( tvec<fftw_complex>		rslt,
												double				a )
{
	for ( int i = 0; i < rslt.len; i++ )
	{
		rslt.val[i][0] += a;
		rslt.val[i][1] += a;
	}
}


/**
 * \brief	The multiplying result of a complex 'tvec' vector and a double value;
 */
void fn_tvec_constMultiply_complex		(	tvec<fftw_complex>		src,
												double				a )
{
	for ( int i = 0; i < src.len; i++ )
	{
		src.val[i][0] *= a;
		src.val[i][1] *= a;
	}
}


/**
 * \brief	Result of adding two complex 'tvec' vectors;
 */
void fn_tvec_add_complex (	tvec<fftw_complex>		dst, 
							tvec<fftw_complex>		src, 
								double				dstCoeff, 
								double				srcCoeff )
{
	for ( int i = 0; i < dst.len; i++ )
	{
		dst.val[i][0] = dstCoeff*dst.val[i][0] + srcCoeff*src.val[i][0];
		dst.val[i][1] = dstCoeff*dst.val[i][1] + srcCoeff*src.val[i][1];
	}
}


/**
 * \brief	The dot multiplying result of two complex 'tvec' vectors;
 *
 * \note!	F1.len == F2.len == rslt.len;
 */
void fn_tvec_dotMultiply_complex (	tvec<fftw_complex>		F1,
									tvec<fftw_complex>		F2,
									tvec<fftw_complex>		rslt )
{
	for (int i = 0; i < rslt.len; i++)
	{
		// (a+ib)(c+id) = (ac-bd) + i(ad+bc);
		rslt.val[i][0] = F1.val[i][0]*F2.val[i][0] - F1.val[i][1]*F2.val[i][1];
		rslt.val[i][1] = F1.val[i][0]*F2.val[i][1] + F1.val[i][1]*F2.val[i][0];
	}
}


/**
 * \brief	The summation of the dot multiplying result of two complex 'tvec' vectors;
 * 
 * \note!	F1.len == F2.len >= len;
 */
void fn_tvec_dotMultiplySum_complex	(	tvec<fftw_complex>		F1,
										tvec<fftw_complex>		F2,
												int				len,
											fftw_complex		rslt )
{
	fn_complex_setZero ( rslt );
	for ( int i = 0; i < len; i++ )
	{
		// (a+ib)(c+id) = (ac-bd) + i(ad+bc);
		rslt[0] += F1.val[i][0]*F2.val[i][0] - F1.val[i][1]*F2.val[i][1];
		rslt[1] += F1.val[i][0]*F2.val[i][1] + F1.val[i][1]*F2.val[i][0];
	}
	fftw_complex rsltSend;
	rsltSend[0] = rslt[0];		rsltSend[1] = rslt[1];
	MPI_Reduce ( &rsltSend[0], &rslt[0], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
	MPI_Reduce ( &rsltSend[1], &rslt[1], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
}


/**
 * \brief	Save 'tvec' vector;
 */
void fn_tvec_save_complex ( tvec<fftw_complex> src, char filename[] )
{
	/* open file; */
	FILE *fwPath = fopen(filename, "w");
	if ( fwPath == NULL && myrank == 0 )
		printf("%s error.\n", filename);

	/* save data; */
	if ( src.row == 1 && src.col == 1 )
		fprintf(fwPath, "%d\n", src.len);
	else
		fprintf(fwPath, "%d\t%d\t%d\n", src.row, src.col, src.len);
	for ( int i = 0; i < src.len; i++ )
		fprintf(fwPath, "%+.15E\t%+.15E\n", src.val[i][0], src.val[i][1]);
	fclose(fwPath);
}


/**
 * \brief	Print 'tvec' vector;
 */
void fn_tvec_print_complex ( tvec<fftw_complex> src )
{
	/* print data; */
	if ( myrank == 0 )
	{
		if ( src.row == 1 && src.col == 1 )
			printf("len = %d\n", src.len);
		else
			printf("nrow = %d, ncol = %d, nnz = %d\n", src.row, src.col, src.len);
		for ( int i = 0; i < src.len; i++ )
			printf("%+.10E\t%+.10E\n", src.val[i][0], src.val[i][1]);
	}
}

/* --------------------		tvec (complex)	--------------------*/

/* --------------------		tmat (int)	--------------------*/

/**
 * \brief	Save 'tmat' matrix;
 */
void fn_tmat_save_int ( tmat<int> src, char filename[] )
{
	/* open file; */
	FILE *fwPath = fopen(filename, "w");
	if ( fwPath == NULL && myrank == 0 )
		printf("%s error.\n", filename);

	/* save data; */
	fprintf(fwPath, "%d\t%d\t%d\n", src.row, src.col, src.ele);
	for ( int i = 0; i < src.row; i++ )
		for ( int j = 0; j < src.col; j++ )
			fprintf(fwPath, "%+d\n", src.val[i][j]);
	fclose(fwPath);
}


/**
 * \brief	Print 'tmat' matrix;
 */
void fn_tmat_print_int ( tmat<int> src )
{
	/* print data; */
	if ( myrank == 0 )
	{
		printf("nrow = %d, ncol = %d, nnz = %d\n", src.row, src.col, src.ele);
		for ( int i = 0; i < src.row; i++ )
			for ( int j = 0; j < src.col; j++ )
				printf("%+d\n", src.val[i][j]);
	}
}

/* --------------------		tmat (int)	--------------------*/

/* --------------------		tmat (complex)	--------------------*/

/**
 * \brief	Set zero of 'tmat' matrix;
 */
void fn_tmat_setZero_complex ( tmat<fftw_complex> src )
{
	for ( int i = 0; i < src.row; i++ )
		for ( int j = 0; j < src.col; j++ )
			src.val[i][j][0] = 0, src.val[i][j][1] = 0;
}


/**
 * \brief	Save 'tmat' matrix;
 */
void fn_tmat_save_complex ( tmat<fftw_complex> src, char filename[] )
{
	/* open file; */
	FILE *fwPath = fopen(filename, "w");
	if ( fwPath == NULL && myrank == 0 )
		printf("%s error.\n", filename);

	/* save data; */
	fprintf(fwPath, "%d\t%d\t%d\n", src.row, src.col, src.ele);
	for ( int i = 0; i < src.row; i++ )
		for ( int j = 0; j < src.col; j++ )
			fprintf(fwPath, "%+.15E\t%+.15E\n", 
					src.val[i][j][0], src.val[i][j][1]);
	fclose(fwPath);
}


/**
 * \brief	Print 'tmat' matrix;
 */
void fn_tmat_print_complex ( tmat<fftw_complex> src )
{
	/* print data; */
	if ( myrank == 0 )
	{
		printf("nrow = %d, ncol = %d, nnz = %d\n", src.row, src.col, src.ele);
		for ( int i = 0; i < src.row; i++ )
			for ( int j = 0; j < src.col; j++ )
				printf("%+.10E\t%+.10E\n", src.val[i][j][0], src.val[i][j][1]);
	}
}

/* --------------------		tmat (complex)	--------------------*/


/* --------------------		matrix det, inverse --------------------*/


/**
 * \brief	calculate det of matrix;
 *
 * \param	D:		the input matrix;
 * \param	n:		the order of matrix;
 *
 * \return:	the det of matrix;
 */
double det	(	double **D,		int n	)
{
	double d=0;
	
	if(n==1)d=D[0][0];
	if(n==2)d=D[0][0]*D[1][1]-D[0][1]*D[1][0];
	else{
	for(int k=0;k<n;k++){
		/* allocate memory for algebraic complement; */
		double **M;
		M=(double**)malloc((n-1)*sizeof(double*));
		for(int i=0;i<n-1;i++)
			M[i]=(double*)malloc((n-1)*sizeof(double));
			
		/* assign to algebraic complement; */
		for(int i=0;i<n-1;i++)
			for(int j=0;j<n-1;j++)
				M[i][j]=D[i+1][j<k?j:j+1];
				
		/* expand by the first row;
		 * calculate algebraic complement by recurrence;
		 * one can ignore the zero element to accelerate computation; */
		if(D[0][k])
			d+=D[0][k]*det(M,n-1)*(((2+k)%2)?-1:1);
		
		/* release memory; */
		for(int i=0;i<n-1;i++)free(M[i]);
		free(M);
		}
	}
	return d;        
}


/**
 * \brief	calculate the inverse matrix;
 *
 * \param	a:		the input matrix;
 * \param	b:		the inverse matrix of a; the return value;
 * \param	N:		the number of row/column of matrix;
 */
void matrix_inverse	(	double **a,		double **b,		int N	)
{
	using namespace std;
	int i, j, k;
	double max, temp;
	/* define a temporary matrix t; */
	double t[N][N];
	/* save matrix a into matrix t; */
	for (i = 0; i < N; i++)
	{
		for (j = 0; j < N; j++)
		{
			t[i][j] = a[i][j];
		}
	}
	/* initialize matrix b as an identity matrix; */
	for (i = 0; i < N; i++)
	{
		for (j = 0; j < N; j++)
		{
			b[i][j] = (i == j) ? (double)1 : 0;
		}
	}
	/* search the main element on each column; */
	for (i = 0; i < N; i++)
	{
		max = t[i][i];
		/* remember the index of main element on each column; */
		k = i;
		/* search the main element on each column; */
		for (j = i + 1; j < N; j++)
		{
			if (fabs(t[j][i]) > fabs(max))
			{
				max = t[j][i];
				k = j;
			}
		}
		//cout<<"the max number is "<<max<<endl;
		/* exchange if the main element is not on i-th column; */
		if (k != i)
		{
			/* exchange the elements on two columns; */
			for (j = 0; j < N; j++)
			{
				temp = t[i][j];
				t[i][j] = t[k][j];
				t[k][j] = temp;
				/* exchange the elements on two columns for adjoint matrix b; */
				temp = b[i][j];
				b[i][j] = b[k][j];
				b[k][j] = temp;
			}
		}
		if (t[i][i] == 0)
		{
			cout << "\nthe matrix does not exist inverse matrix\n";
			break;
		}
		/* obtain the main elements on columns; */
		temp = t[i][i];
		/* identify the column where the main element is located; */
		//cout<<"\nthe temp is "<<temp<<endl;
		for (j = 0; j < N; j++)
		{
			t[i][j] = t[i][j] / temp;
			b[i][j] = b[i][j] / temp;
		}
		for (j = 0; j < N; j++)
		{
			if (j != i)
			{
				temp = t[j][i];
				/* delete the other elements on this column; */
				for (k = 0; k < N; k++)
				{
					t[j][k] = t[j][k] - temp * t[i][k];
					b[j][k] = b[j][k] - temp * b[i][k];
				}
			}
 
		}
	}
}

/* --------------------		matrix det, inverse --------------------*/
