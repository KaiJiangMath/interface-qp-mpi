/*! \file	GenJacPoly.cpp
 *
 * \brief	See the book: J. Shen, T. Tang and L. Wang, Spectral Methods:
 *				Algorithms, Analysis and Applications, Springer Series in
 *				Compuational Mathematics, 41, Springer, 2011.
 */

#include "Head.h"
#include "Data.h"
#include "DataOperators.h"
#include "Mytimer.h"
#include "functs.h"

/**
 * \brief	obtain general Jacobi polynomials and inner product matrices;
 */
void fn_obt_system_gen_Jac_poly ( stu_GJP_var	*sGJPv )
{
	if ( myrank == 0 )
		printf(" <========== Generate GJPs ==========> \n\n");
	mytimer_t timer;
	timer.reset();
	timer.start();

	/* memory allocation; */
	fn_GJP_memory_allocation (sGJPv);
//	sGJPv->isTest = true;		// true: save; false: not save;

	/* obtain LGL nodes and weights; */
	fn_obt_Leg_Gau_Lob (sGJPv->x, sGJPv->w, sGJPv->xlen);

	// obtain general Jacobi polynomial;
	fn_obt_gen_Jac_poly (sGJPv, sGJPv->d0JJ, 0);
	fn_obt_gen_Jac_poly (sGJPv, sGJPv->d1JJ, 1);
	fn_obt_gen_Jac_poly (sGJPv, sGJPv->d2JJ, 2);
	if ( strcmp(model_type, "LP") == 0 )
	{
		fn_obt_gen_Jac_poly (sGJPv, sGJPv->d3JJ, 3);
		fn_obt_gen_Jac_poly (sGJPv, sGJPv->d4JJ, 4);
	}

	/**
	 * obtain the inner product matrices; 
	 *		CCS sparse matrix;
	 * meanwhile save data;
	 */
	if ( strcmp(model_type, "LB") == 0 )
	{
		sGJPv->innSMatd0JJ = fn_innSMat_gen_Jac_poly(sGJPv, sGJPv->d0JJ, 0, 1e-8);
		sGJPv->innSMatd1JJ = fn_innSMat_gen_Jac_poly(sGJPv, sGJPv->d1JJ, 1, 1e-8);
		sGJPv->innSMatd2JJ = fn_diag_innSMat_gen_Jac_poly(sGJPv, sGJPv->d2JJ, 2, 1e-16);
	}
	else if ( strcmp(model_type, "LP") == 0 )
	{
		sGJPv->innSMatd0JJ = fn_innSMat_gen_Jac_poly(sGJPv, sGJPv->d0JJ, 0, 1e-8);
		sGJPv->innSMatd1JJ = fn_innSMat_gen_Jac_poly(sGJPv, sGJPv->d1JJ, 1, 1e-8);
		sGJPv->innSMatd2JJ = fn_innSMat_gen_Jac_poly(sGJPv, sGJPv->d2JJ, 2, 1e-8);
		sGJPv->innSMatd3JJ = fn_innSMat_gen_Jac_poly(sGJPv, sGJPv->d3JJ, 3, 1e-7);
		sGJPv->innSMatd4JJ = fn_diag_innSMat_gen_Jac_poly(sGJPv, sGJPv->d4JJ, 4, 1e-16);
	}

	/* zoom x; */
	for (int i = 0; i < sGJPv->xlen; i++)
		sGJPv->x.val[i] *= sGJPv->x_range;

	/* save x and w; */
	if ( myrank == 0 )
	{
		char fname[FILELEN];
		/* save x; */
		sprintf(fname, "%s/x.dat", rsltDir);
		fn_tvec_save<double> ( sGJPv->x, fname );
		/* save w; */
		sprintf(fname, "%s/w.dat", rsltDir);
		fn_tvec_save<double> ( sGJPv->w, fname );
	}

	timer.pause();
	if ( myrank == 0 )
		printf("\n\t ***** time cost of GJPs: %f seconds *****\n\n", 
				timer.get_current_time());
}


/**
 * \brief	Allocation memory for calculating general Jacobi polynomials;
 */
void fn_GJP_memory_allocation (stu_GJP_var *sGJPv)
{
	/* load parameters; */
	int	xlen = sGJPv->xlen;
	int nd	 = sGJPv->nd;

	/* Initialization; */
	sGJPv->x	= fn_tvec_init<double> ( xlen );
	sGJPv->w	= fn_tvec_init<double> ( xlen );
	sGJPv->d0JJ = fn_tmat_init<double> ( xlen, nd );
	sGJPv->d1JJ = fn_tmat_init<double> ( xlen, nd );
	sGJPv->d2JJ = fn_tmat_init<double> ( xlen, nd );
	sGJPv->d3JJ = fn_tmat_init<double> ( xlen, nd );
	sGJPv->d4JJ = fn_tmat_init<double> ( xlen, nd );
}

	
/**
 * \brief	Construct General Jacobi polynomial, see Shen, Wang and Tang book Chp 6;
 *
 * \param	sGJPv	the structure body about GJP;
 * \param	dsJJ	the s-th derivative of GJP;
 * \param	order	the derivative order; 0, 1, 2,...
 */
void fn_obt_gen_Jac_poly	(	stu_GJP_var		*sGJPv,
								tmat<double>	dsJJ, 
									int			order )
{
	double	jAlpha	= sGJPv->alpha;
	double	jBeta	= sGJPv->beta;
	int		nd		= sGJPv->nd;
	int		xlen	= sGJPv->xlen;

	/* obtain Jacobi polynomial; */
	tmat<double> dsJpoly = fn_tmat_init<double> ( xlen, nd+order );
	fn_obt_Jac_poly ( sGJPv->x, dsJpoly, jAlpha-order, jBeta-order );

	/* obtain general Jacobi polynomial; */
	double xTmp, coefTmp;
	for (int i = 0; i < xlen; i++)
	{
		/**
		 * (1-x)^{k-m} * (1+x)^{l-m}; 
		 *		see README_interface;
		 */
		xTmp =	pow(1.0 - sGJPv->x.val[i], jAlpha-order) * 
				pow(1.0 + sGJPv->x.val[i], jBeta-order);
		for (int j = 0; j < nd; j++)
		{
			/* calculate coefficient; */
			coefTmp = pow(-2.0, order);
			for (int k = j+1; k <= j+order; k++)
				coefTmp *= k;
			coefTmp /= pow(sGJPv->x_range, order);
			/* calculate general Jacobi polynomial; */
			dsJJ.val[i][j] = coefTmp * xTmp * dsJpoly.val[i][j+order];
		}
	}

	/* save GJP and its derivatives; */
	if ( myrank == 0 )
	{
		if ( order == 0 )
		{
			char dsJJfwFile[FILELEN];
			sprintf(dsJJfwFile, "%s/d%dGJP.dat", rsltDir, order);
			fn_tmat_save<double> ( dsJJ, dsJJfwFile );
		}
		else
		{
			if ( sGJPv->isTest )
			{
				char dsJJfwFile[FILELEN];
				sprintf(dsJJfwFile, "%s/d%dGJP.dat", rsltDir, order);
				fn_tmat_save<double> ( dsJJ, dsJJfwFile );
			}
		}
	}

	/* release memory; */
	fn_tmat_free<double> ( dsJpoly );
}


/**
 * \brief	Construct the inner product matrices;
 *
 * \param	sGJPv			the structure body for GJP;
 * \param	dsJJ			the general Jacobi polynomial; can be replaced by d1JJ, d2JJ, d3JJ, d4JJ;
 * \param	TOL				the maximal allowable tolerance for sparse CCS matrix;
 *
 * \return	innSMatdsJJ		(dsJJ, dsJJ)_{w} = dsJJ.' * diag(w) * dsJJ; a sparse matrix;
 *							dsJJ.': (nd)x(xlen); 
 *							diag(w): (xlen)x(xlen); 
 *							dsJJ: (xlen)x(nd);
 */
tCCSmat<double> fn_innSMat_gen_Jac_poly (	stu_GJP_var		*sGJPv,
											tmat<double>	dsJJ, 
												int			order,
												double		TOL )
{
	int		nd		= sGJPv->nd;
	int		xlen	= sGJPv->xlen;

	/* construct temporary vectors; */
	int *JAvector = (int *)malloc(sizeof(int) * (nd+1));
	int *IAvector = (int *)malloc(sizeof(int) * (nd*nd));
	double *valvector = (double *)malloc(sizeof(double) * (nd*nd));

	/**
	 * calculate (dsJJ, dsJJ)_{w} = dsJJ.' * diag(w) * dsJJ;
	 * firstly calculate dsJJ.' * diag(w);
	 *		corresponding to column transformation;
	 */
	JAvector[0] = 0;
	double valTmp;
	int nnzTmp = 0;
	for ( int j = 0; j < nd; j++ )
	{
		for ( int i = 0; i < nd; i++ )
		{
			valTmp = 0.0;
			for ( int k = 0; k < xlen; k++ )	// (dsJJ, dsJJ)_{w}; (i,j) position;
			{
				valTmp += dsJJ.val[k][i] * sGJPv->w.val[k] * dsJJ.val[k][j];
			}
			if ( fabs(valTmp) > TOL )			// nonzero element;
			{
				IAvector[nnzTmp] = i;
				valvector[nnzTmp] = valTmp;
				nnzTmp++;
			}
		}
		JAvector[j+1] = nnzTmp;
	}

	/* construct the sparse CCS matrix; */
	tCCSmat<double> innSMatdsJJ = fn_tCCSmat_init<double> ( nd, nd, nnzTmp );
	for ( int i = 0; i < nd+1; i++ )
		innSMatdsJJ.JA[i] = JAvector[i];
	for ( int i = 0; i < nnzTmp; i++ )
	{
		innSMatdsJJ.IA[i]	= IAvector[i];
		innSMatdsJJ.val[i]	= valvector[i];
	}

	/* save results; */
	if ( sGJPv->isTest && myrank == 0 )
	{
		char fwFile[FILELEN];
		sprintf(fwFile, "%s/d%dGJPsmat.dat", rsltDir, order);
		fn_tCCSmat_save<double> ( innSMatdsJJ, fwFile );
		printf("\t derivative order: %d \t element number: %d\n", order, innSMatdsJJ.nnz);
	}

	/* free memory; */
	free(JAvector);
	free(IAvector);
	free(valvector);
	return innSMatdsJJ;
}


/**
 * \brief	Construct the inner product matrices; only keep the diagonal elements;
 *
 * \param	sGJPv			the structure body for GJP;
 * \param	dsJJ			the general Jacobi polynomial; can be replaced by d1JJ, d2JJ, d3JJ, d4JJ;
 * \param	TOL				the maximal allowable tolerance for sparse CCS matrix;
 *
 * \return	innSMatdsJJ		(dsJJ, dsJJ)_{w} = dsJJ.' * diag(w) * dsJJ; a sparse matrix;
 *							dsJJ.': (nd)x(xlen); 
 *							diag(w): (xlen)x(xlen); 
 *							dsJJ: (xlen)x(nd);
 */
tCCSmat<double> fn_diag_innSMat_gen_Jac_poly	(	stu_GJP_var		*sGJPv,
													tmat<double>	dsJJ, 
														int			order,
														double		TOL)
{
	int		nd		= sGJPv->nd;
	int		xlen	= sGJPv->xlen;

	/* construct temporary vectors; */
	int *JAvector = (int *)malloc(sizeof(int) * (nd+1));
	int *IAvector = (int *)malloc(sizeof(int) * (nd*nd));
	double *valvector = (double *)malloc(sizeof(double) * (nd*nd));

	/**
	 * calculate (dsJJ, dsJJ)_{w} = dsJJ.' * diag(w) * dsJJ;
	 * firstly calculate dsJJ.' * diag(w);
	 *		corresponding to column transformation;
	 */
	JAvector[0] = 0;
	double valTmp;
	int nnzTmp = 0;
	for ( int i = 0; i < nd; i++ )
	{
		valTmp = 0.0;
		for ( int k = 0; k < xlen; k++ )	// (dsJJ, dsJJ)_{w}; (i,i) position;
		{
			valTmp += dsJJ.val[k][i] * sGJPv->w.val[k] * dsJJ.val[k][i];
		}
		if ( fabs(valTmp) > TOL )			// nonzero element;
		{
			IAvector[nnzTmp] = i;
			valvector[nnzTmp] = valTmp;
			nnzTmp++;
		}
		JAvector[i+1] = nnzTmp;
	}

	/* construct the sparse CCS matrix; */
	tCCSmat<double> innSMatdsJJ = fn_tCCSmat_init<double> ( nd, nd, nnzTmp );
	for ( int i = 0; i < nd+1; i++ )
		innSMatdsJJ.JA[i] = JAvector[i];
	for ( int i = 0; i < nnzTmp; i++ )
	{
		innSMatdsJJ.IA[i] = IAvector[i];
		innSMatdsJJ.val[i] = valvector[i];
	}

	/* save results; */
	if ( sGJPv->isTest && myrank == 0 )
	{
		char fwFile[FILELEN];
		sprintf(fwFile, "%s/d%dGJPsmat.dat", rsltDir, order);
		fn_tCCSmat_save<double> ( innSMatdsJJ, fwFile );
		printf("\t derivative order: %d \t element number: %d\n", order, innSMatdsJJ.nnz);
	}

	// free memory;
	free(JAvector);
	free(IAvector);
	free(valvector);
	return innSMatdsJJ;
}


/**
 * \brief	obtain Legendre polynomial; 
 *			see page 94;
 *
 * \param	y		the Legendre polynomial about x;
 * \param	n		the degree of the Legendre polynomial;
 */
void fn_obt_Leg_poly	(	tvec<double>	x,
							tvec<double>	y, 
								int			n )
{
	int xlen = x.len;	// or y.len;

	if ( n == 0 )
		for (int i = 0; i < xlen; i++)	y.val[i] = 1.0;			// L_{0}(x) = 1;
	else if ( n == 1 )
		for (int i = 0; i < xlen; i++)	y.val[i] = x.val[i];	// L_{1}(x) = x;
	else
	{
		double L0, L1, L2;
		for (int i = 0; i < xlen; i++)
		{
			L0 = 1.0, L1 = x.val[i], L2 = 0.0;
			for (int k = 1; k < n; k++)
			{
				/**
				 * three-term recurrence relation; 
				 * (k+1)L_{k+1}(x) = (2k+1)x L_{k}(x) - k L_{k-1}(x), k = 1,2,...
				 */
				L2 = ( (2*k+1)*x.val[i]*L1 - k*L0 ) / (k+1);
				L0 = L1, L1 = L2;
			}
			y.val[i] = L2;
		}
	}
}


/**
 * \brief	obtain the Legendre polynomial and its first-order derivative; 
 *			see page 94, 95;
 *
 *	\param	y		the Legendre polynomial about x;
 *	\param	dy		the first-order derivative of the Legendre polynomial about x;
 *	\param	n		the degree of the Legendre polynomial;
 */
void fn_obt_first_deri_Leg_poly		(	tvec<double>	x,
										tvec<double>	y, 
										tvec<double>	dy, 
											int			n )
{
	int xlen = x.len;	// or y.len; or dy.len;

	if ( n == 0 )
		for (int i = 0; i < xlen; i++)
		{
			y.val[i] = 1.0, dy.val[i] = 0.0;		// L_{0}(x) = 1, L_{0}'(x) = 0;
		}
	else if ( n == 1 )
		for (int i = 0; i < xlen; i++)
		{
			y.val[i] = x.val[i], dy.val[i] = 1.0;	// L_{1}(x) = x, L_{1}'(x) = 1;
		}
	else
	{
		double L0, L1, L2;
		double dL0, dL1, dL2;
		for (int i = 0; i < xlen; i++)
		{
			L0 = 1.0, L1 = x.val[i], L2 = 0.0;		// L_{0}(x) = 1, L_{1}(x) = x;
			dL0 = 0.0, dL1 = 1.0, dL2 = 0.0;		// L_{0}'(x) = 0, L_{1}'(x) = 1;
			for (int k = 1; k < n; k++)
			{
				/**
				 * three-term recurrence relation;
				 * (k+1)L_{k+1}(x) = (2k+1)x L_{k}(x) - k L_{k-1}(x), k = 1,2,...
				 */
				L2 = ( (2*k+1)*x.val[i]*L1 - k*L0 ) / (k+1);
				/* L_{n+1}'(x) = L_{n-1}'(x) + (2n+1)L_{n}(x), n = 1,2,... */
				dL2 = dL0 + (2*k+1)*L1;
				L0 = L1, L1 = L2, dL0 = dL1, dL1 = dL2;
			}
			y.val[i] = L2, dy.val[i] = dL2;
		}
	}
}


/**
 * \brief	obtain Legendre-Gauss-Lobatto (LGL) nodes and weights by Newton method;
 *			see page 99;
 * \param	x:		the LGL nodes;
 * \param	w:		the corresponding weights;
 */
void fn_obt_Leg_Gau_Lob		(	tvec<double>	x,
								tvec<double>	w,
									int			n )
{
	/* Compute the initial guess of the interior LGL points; */
	double theta;
	int N = n-1, NN = N-1;
	tvec<double> sigma = fn_tvec_init<double> ( N );
	for (int k = 1; k <= N; k++)
	{
		/* theta_{k} = (4k-1)/(4N+2) * pi; */
		theta = (4.0*k-1.0) / (4.0*N+2.0) * PI;
		/* sigma_{k} = [ 1 - (N-1)/(8N^3) - (39-28/sin^2\theta_{k})/(384N^4) ] * cos\theta_{k}; */
		sigma.val[k-1] = 1.0 - (N-1.0)/8.0/pow(N,3);
		sigma.val[k-1] -= (39.0-28.0/pow(sin(theta),2)) /384.0/pow(N,4);
	   	sigma.val[k-1] *= -cos(theta);
	}

	/* memory allocation; */
	tvec<double> x0 = fn_tvec_init<double> ( NN );
	tvec<double> y  = fn_tvec_init<double> ( NN );
	tvec<double> dy = fn_tvec_init<double> ( NN );
	for (int i = 0; i < NN; i++) x0.val[i] = 0.0;
	for (int i = 0; i < NN; i++) y.val[i] = 0.0;
	for (int i = 0; i < NN; i++) dy.val[i] = 0.0;

	/* boundary values; */
	x.val[0] = -1.0;
	x.val[N] = 1.0;
	w.val[0] = 2.0/N/(N+1.0);
	w.val[N] = w.val[0];

	/** 
	 * the initial guess: 
	 *		x_{j}^{0} = (sigma_{j} + sigma_{j+1}) / 2, 1\leq j\leq (N-1);
	 */
	for (int j = 0; j < NN; j++)
		x0.val[j] = 0.5 * (sigma.val[j] + sigma.val[j+1]);

	/* Newton method; */
	double x0j, x1j, yj, dyj, res = 1.0, tmp;
	double TOL = 10*EPS;		// error tolerance for stopping iteration;
	while ( res >= TOL )
	{
		/* obtain the Legendre polynomial and its first-order derivative; */
		fn_obt_first_deri_Leg_poly(x0, y, dy, N);
		/* update; */
		res = 0.0;
		for (int j = 0; j < NN; j++)
		{
			x0j = x0.val[j], yj = y.val[j], dyj = dy.val[j];
			x1j = x0j - (1.0-x0j*x0j)*dyj / ( 2.0*x0j*dyj - N*(N+1.0)*yj );
			tmp = fabs(x1j-x0j);
			res = tmp > res ? tmp : res;
			x0.val[j] = x1j;
		}
	}
	for (int j = 1; j < N; j++)
	{
		x.val[j] = x0.val[j-1];
		/* Use the weight expression (3.188) to compute the weights; */
		w.val[j] = w.val[0] / (y.val[j-1]*y.val[j-1]);
	}

	/* release memory; */
	fn_tvec_free<double> ( sigma );
	fn_tvec_free<double> ( x0 );
	fn_tvec_free<double> ( y  );
	fn_tvec_free<double> ( dy );
}


/**
 * \brief	obtain Jacobi polynomial;
 *			see page 74;
 *
 * \param	y		the Jacobi polynomial about x;
 *
 * \note!	x.len = y.row;
 */
void fn_obt_Jac_poly	(	tvec<double>	x,
							tmat<double>	y, 
								double		jAlpha,
								double		jBeta )
{
	int xlen = y.row;	// or x.len;
	int n	 = y.col-1;

	if ( n == 0 )
		// J_{0}(x) = 1;
		for (int i = 0; i < xlen; i++)	y.val[i][0] = 1.0;
	else if ( n == 1 )
		// J_{1}(x) = 0.5(alpha+beta+2)x + 0.5(alpha-beta);
		for (int i = 0; i < xlen; i++)
		{
			y.val[i][0] = 1.0;
			y.val[i][1] = 0.5*(jAlpha+jBeta+2.0)*x.val[i] + 0.5*(jAlpha-jBeta);
		}
	else
	{
		double J0, J1, J2;
		double abc, a, b, c;
		double jAB = jAlpha + jBeta;
		for (int i = 0; i < xlen; i++)
		{
			J0 = 1.0, J2 = 0.0;
			J1 = 0.5*(jAlpha+jBeta+2.0)*x.val[i] + 0.5*(jAlpha-jBeta);
			y.val[i][0] = J0, y.val[i][1] = J1;
			for (int k = 1; k < n; k++)
			{
				/**
				 * three-term recurrence relation; (ignoring alpha, beta);
				 * J_{k+1}(x) = (a_{k}x - b_{k}) J_{k}(x) - c_{k} J_{k-1}(x), k = 1,2,...  
				 * abc = 2(k+1)(k+alpha+beta+1)(2k+alpha+beta); 
				 * a_{k} = (2k+alpha+beta)(2k+alpha+beta+1)(2k+alpha+beta+2)/abc; 
				 * b_{k} = (beta^2-alpha^2)(2k+alpha+beta+1)/abc; 
				 * c_{k} = 2(k+alpha)(k+beta)(2k+alpha+beta+2)/abc;
				 */
				/* coefficients; */
				abc = 2.0*(k+1.0)*(k+jAB+1.0)*(2.0*k+jAB);
				a = (2.0*k+jAB)*(2.0*k+jAB+1.0)*(2.0*k+jAB+2.0);
				b = jAB*(jBeta-jAlpha)*(2.0*k+jAB+1.0);
				c = 2.0*(k+jAlpha)*(k+jBeta)*(2.0*k+jAB+2.0);
				/* update Jacobi polynomial; */
				J2 = (a*x.val[i]-b) * J1 - c*J0;
				J2 /= abc;
				J0 = J1, J1 = J2;
				y.val[i][k+1] = J2;
			}
		}
	}
}


/**
 * \brief	obtain the Jacobi polynomial and its first-order derivative;
 *			its first-order derivative can be obtained directly;
 *			see page 74;
 *
 * \param	y		the Legendre polynomial about x;
 * \param	dy		the first-order derivative of the Legendre polynomial about x;
 *
 * \note!	x.len = y.row = dy.row;
 */
void fn_obt_first_deri_Jac_poly		(	tvec<double>	x,
										tmat<double>	y,
										tmat<double>	dy, 
											double		jAlpha,
											double		jBeta )
{
	int xlen = y.row;	// or x.len; or dy.row;
	int n	 = y.col-1;

	if ( n == 0 )
		for (int i = 0; i < xlen; i++)
		{
			// J_{0}(x) = 1, J_{0}'(x) = 0;
			y.val[i][0] = 1.0, dy.val[i][0] = 0.0;
		}
	else if ( n == 1 )
		for (int i = 0; i < xlen; i++)
		{
			// J_{0}(x) = 1, J_{0}'(x) = 0;
			y.val[i][0] = 1.0, dy.val[i][0] = 0.0;
			// J_{1}(x) = 0.5(alpha+beta+2)x + 0.5(alpha-beta);
			// J_{1}'(x) = 0.5(alpha+beta+2);
			y.val[i][1] = 0.5*(jAlpha+jBeta+2.0)*x.val[i] + 0.5*(jAlpha-jBeta);
			dy.val[i][1] = 0.5*(jAlpha+jBeta+2.0);
		}
	else
	{
		double J0, J1, J2;
		double dJ0, dJ1, dJ2;
		double abc, a, b, c;
		double jAB = jAlpha + jBeta;
		for (int i = 0; i < xlen; i++)
		{
			J0 = 1.0, dJ0 = 0.0, J2 = 0.0, dJ2 = 0.0;
			J1 = 0.5*(jAlpha+jBeta+2.0)*x.val[i] + 0.5*(jAlpha-jBeta);
			dJ1 = 0.5*(jAlpha+jBeta+2.0);
			y.val[i][0] = J0, dy.val[i][0] = dJ0;
			y.val[i][1] = J1, dy.val[i][1] = dJ1;
			for (int k = 1; k < n; k++)
			{
				/**
				 * coefficients; 
				 * abc = 2(k+1)(k+alpha+beta+1)(2k+alpha+beta); 
				 * a_{k} = (2k+alpha+beta)(2k+alpha+beta+1)(2k+alpha+beta+2)/abc; 
				 * b_{k} = (beta^2-alpha^2)(2k+alpha+beta+1)/abc; 
				 * c_{k} = 2(k+alpha)(k+beta)(2k+alpha+beta+2)/abc;
				 */
				abc = 2.0*(k+1.0)*(k+jAB+1.0)*(2.0*k+jAB);
				a = (2.0*k+jAB)*(2.0*k+jAB+1.0)*(2.0*k+jAB+2.0);
				b = jAB*(jBeta-jAlpha)*(2.0*k+jAB+1.0);
				c = 2.0*(k+jAlpha)*(k+jBeta)*(2.0*k+jAB+2.0);
				/**
				 * update Jacobi polynomial; 
				 * three-term recurrence relation; (ignoring alpha, beta);
				 * J_{k+1}(x) = (a_{k}x - b_{k}) J_{k}(x) - c_{k} J_{k-1}(x), k = 1,2,...
				 */
				J2 = (a*x.val[i]-b) * J1 - c*J0;
				J2 /= abc;
				/**
				 * update its first-order derivative;
				 * three-term recurrence relation; (ignoring alpha, beta);
				 * J_{k+1}'(x) = (a_{k}x - b_{k}) J_{k}'(x) - c_{k} J_{k-1}'(x) 
				 *					+ a_{k} J_{k}(x), k = 1,2,...
				 */
				dJ2 = (a*x.val[i]-b) * dJ1 - c*dJ0 + a*J1;
				dJ2 /= abc;
				J0 = J1, J1 = J2, dJ0 = dJ1, dJ1 = dJ2;
				y.val[i][k+1] = J2, dy.val[i][k+1] = dJ2;
			}
		}
	}
}
