#ifndef __Data_h
#define __Data_h

#include "Head.h"

#define DATALEN 32
#define STRLEN	64
#define WORDLEN 128
#define FILELEN0 256
#define FILELEN 512
#define ZEROTOL 1e-10
#define SPECTOL 1e-10	// the lowest spectra intensity of allowed density modes;

extern char main_type[STRLEN];		// the type of the main code;
extern char model_type[STRLEN];		// the type of the Landau model; LB/LP;
extern const char *para_flag;		// the flag of parameters;
extern char paraDir[FILELEN0];		// the directory about parameters;
extern char rsltDir[FILELEN0];		// the directory to save results;
extern unsigned flags;				// the flag of FFTW;

extern int nprocs, myrank;			// the number of process and the rank of this process;
extern int *displs_bulk,	*recvCount_bulk; // offset and receive count of bulk phases;
extern int *displs_sys,		*recvCount_sys;  // offset and receive count of interfaces;
extern ptrdiff_t alloc_local,	   local_n0,	 local_0_start;		// local memory and index of bulk phases;
extern ptrdiff_t alloc_local_sys, local_n0_sys, local_0_start_sys; // local memory and index of interfaces;

/* ------------------------------------------------------------------*/

typedef struct myVec
{
	double* Data;
	int* Index;
} mySortVec;


typedef struct myWave
{
	double* Wave;
	double Data;
} mySortWave;


/**
 * \struct	stu_sort;
 * \brief	define the structure body to return indexes after sorting;
 */
typedef struct stu_sort
{
	int		ind;
	double	val;
} stu_sort;

/* ------------------------------------------------------------------*/

/* --------------------		Sparse matrix	--------------------*/

/**
 * \brief	Sparse matrix of <typename T> type in CSR format;
 *
 * CSR Format (IA,JA,A) in <typename T>;
 *
 * <typename T>:	int, float, double;
 *
 * \note	The starting index of A is 0;
 */
template <typename T>
struct tCSRmat
{
    int row;	// row number of matrix A;
    int col;	// column of matrix A;
    int nnz;	// number of nonzero entries;
    int *IA;	// integer array of row pointers, the size is row+1;
    int *JA;	// integer array of column indexes, the size is nnz;
     T	*val;	// nonzero entries of A;
};


/**
 * \brief	Sparse matrix of <typename T> type in CCS format;
 *
 * CCS Format (JA,IA,A) in <typename T>;
 *
 * <typename T>:	int, float, double;
 *
 * \note	The starting index of A is 0.
 */
template <typename T>
struct tCCSmat
{
    int row;	// row number of matrix A;
    int col;	// column of matrix A;
    int nnz;	// number of nonzero entries;
    int *JA;	// integer array of column pointers, the size is col+1;
    int *IA;	// integer array of row indexes, the size is nnz;
     T	*val;	// nonzero entries of A;
};

/* --------------------		Sparse matrix	--------------------*/

/* --------------------		Dense matrix	--------------------*/

/**
 * \brief	Vector with 'len' entries of <typename T> type;
 *
 * <typename T>:	int, float, double;
 *
 * \note!	row = 1, col = 1;	for a vector;
 *
 */
template <typename T>
struct tvec
{
	int row;	// number of rows; default value is 1;
	int col;	// number of cols; default value is 1;
    int len;	// number of elements;
     T	*val;	// actual vector entries;
};


/**
 * \brief	Matrix with 'row x col' entries of <typename T> type;
 *
 * <typename T>:	int, float, double;
 *
 */
template <typename T>
struct tmat
{
    int row;	// number of rows;
	int col;	// number of columns;
	int ele;	// number of element; equal to row * col;
	 T	**val;	// actual matrix entries;
};

/* --------------------		Dense matrix	--------------------*/

/* --------------------		FFTW, GJP, bnd	--------------------*/

/**
 * \brief	Some variables for convenient FFTW;
 *
 */
typedef struct stu_fftw_var
{
	tmat<int>	 indKspace;			// k in Fourier space;	row = alloc_local, col = dimCpt;
	tmat<double> projPlane;			// R'PBk;				row = alloc_local, col = dimPhy;
	tvec<double> Gsquare;			// |R'PBk|^2;			len = alloc_local;
	tvec<int>	 ind;				// the global index for truncation;

	tmat<double> projPlaneAll;		// projPlane with all informations;	row = cplxDofs, col = dimPhy;

	tvec<fftw_complex> rhoReal;		// rho in real space;						len = alloc_local;
	tvec<fftw_complex> rhoCplx;		// rho in reciprocal space;					len = alloc_local;
	tvec<fftw_complex> gradient;	// gradient part of the considering model;	len = alloc_local;
	tvec<fftw_complex> fftw_Ctmp;	// temporal variables for FFTW;				len = alloc_local;
	tvec<fftw_complex> fftw_Rtmp;	// temporal variables for FFTW;				len = alloc_local;
	tvec<fftw_complex> cplxTmp;		// temporal variables for FFTW;				len = alloc_local;
	tvec<fftw_complex> quadTerm;	// quadratic term in reciprocal space;		len = alloc_local;
	tvec<fftw_complex> cubTerm;		// cubic term in reciprocal space;			len = alloc_local;
	tvec<fftw_complex> quarTerm;	// quartic term in reciprocal space;		len = alloc_local;

	fftw_plan Planc2cFord;		// FFTW plan; forward transformation;
	fftw_plan Planc2cBack;		// FFTW plan; backward transformation;
} stu_fftw_var;


/**
 * \struct stu_GJP_var;
 * \brief  Some variables for convenient General Jacobi Polynomials;
 *
 */
typedef struct stu_GJP_var
{
	int				polyDegree;		// the degree of general Jacobi polynomials;
	bool			isTest;			// true: save; false: not save;

	double			alpha, beta;	// J^{alpha, beta}; see Jacobi polynomials;
	int				nd;				// the degree of the Legendre polynomial;
	int				xlen;			// the length of x, w;
	double			x_range;		// the distance between the two anchoring planes; (* 2);
	tvec<double>	x;				// the Legendre-Gauss-Lobatto (LGL) nodes;	len = xlen;
	tvec<double>	w;				// the corresponding weights;				len = xlen;

	tmat<double>	d0JJ;			// general Jacobi polynomial (GJP);		row = xlen; col = nd;
	tmat<double>	d1JJ;			// the first-order derivative of GJP;	row = xlen; col = nd;
	tmat<double>	d2JJ;			// the second-order derivative of GJP;	row = xlen; col = nd;	
	tmat<double>	d3JJ;			// the third-order derivative of GJP;	row = xlen; col = nd;
	tmat<double>	d4JJ;			// the fourth-order derivative of GJP;	row = xlen; col = nd;	

	tCCSmat<double> innSMatd0JJ;	// (d0JJ, d0JJ)_{w};		row = col = nd;
	tCCSmat<double> innSMatd1JJ;	// (d1JJ, d1JJ)_{w};		row = col = nd;
	tCCSmat<double> innSMatd2JJ;	// (d2JJ, d2JJ)_{w};		row = col = nd;
	tCCSmat<double> innSMatd3JJ;	// (d3JJ, d3JJ)_{w};		row = col = nd;
	tCCSmat<double> innSMatd4JJ;	// (d4JJ, d4JJ)_{w};		row = col = nd;

} stu_GJP_var;


/**
 * \brief	The introducted polynomials to make boundaries homogeous;
 */
typedef struct stu_bnd_var
{
	int		polyDegree;		// the degree of general Jacobi polynomials;
	bool	isSave;			// true: save; false: not save;
	bool	isTest;			// true: save; false: not save;

	int		nd;				// the degree of the Legendre polynomial;
	int		xlen;			// the length of x, w;
	int		cplxDofs;		// the number of discrete points in computational space;
	double	x_range;		// the distance between the two anchoring planes; (* 2);

	/* project rhoCplx on the GJP space; */
	tvec<fftw_complex> rhoJCplx;	// rho on each process;

	/**
	 * the polynomials for deal with inhomogeneous boundary conditions;
	 *	LB: d0bnd, d1bnd, d2bnd;		total 3;
	 *	LP: d0bnd, d1bnd, ..., d6bnd;	total 7;
	 */
	tvec<fftw_complex> d0bnd;
	tvec<fftw_complex> d1bnd;
	tvec<fftw_complex> d2bnd;
	tvec<fftw_complex> d3bnd;
	tvec<fftw_complex> d4bnd;
	tvec<fftw_complex> d5bnd;
	tvec<fftw_complex> d6bnd;

} stu_bnd_var;

/* --------------------		FFTW, GJP, bnd		--------------------*/


/* --------------------		bulk and system		--------------------*/

/**
 * \struct stu_bulk_param;
 * \brief  Input the parameters of the bulk phase;
 *
 * Input the parameters of the bulk phase, reading from disk file;
 */
typedef struct stu_bulk_param
{
	char	phase[STRLEN];		// the type of the bulk phase;
	int		sflag;				// 1: left; 2: right;
	int		print_level;		// the printing level; 0 is not displayed;

	/* the rotation order and rotation angles around x,y,z-axis; */
	char motion_order[STRLEN];	// the motion order including rotate_transl, transl_rotate;
	char rotate_order[STRLEN];	// the rotation order; work for 3D; such as xyz;
	tvec<double> rotate_angle;	// counterclockwise rotation angle;
								// rotation order and axis are determined by rotate_order;

	/* the translation vector; */
	tvec<double> transl_var;	// translation along x,y,z-direction;

	/* the calculation and plotting parameters; */
	int		Four_num;			// the number of Fourier discrete points along each direction;
	double	enlarge;			// enlarge the plotting range;
	int		plot_num;			// the number of plotting points along each direction;
	int		skip_half;			// 1 to skip the half of Four_num in plotting, default 0;

	/* parameters depending on the bulk phase; */
	int		dimPhy;				// dimensionality of physical space;
	int		dimCpt;				// dimensionality of computational space;
	int		nfold;				// only for the projection matrix of quasicrystals;
	/* initial reciprocal box; */
	char	boxType[STRLEN];	// the type of computational box including 'cube', 'cuboid', 'hex';
	tmat<double> rcpBox;		// computational box in reciprocal space;
	/* initial values; */
	int		initNum;			// the number of initial spectral points;
	tmat<int> initIndex;		// the indexes of initial values;
	tvec<fftw_complex> initCoeff; // the Fourier coefficients of initial values;

	/* spatial grids; */
	tvec<int> NCpt;				// space grid;
	tvec<ptrdiff_t> pNCpt;		// space grid for parallel code;
	int		realDofs;			// the number of discrete points in physical space;
	int		cplxDofs;			// the number of discrete points in computational space;

	/* computational box, projection matrix, rotational matrix, and translation vector; */
	tmat<double> dirBox;			// computational box in direct space;
	tmat<double> projMat;			// projection matrix; denote P;
	tmat<double> rotateMat;			// rotational matrix; denote R;
	tmat<double> rotateProjBoxMat;	// R'PB;
	tvec<double> translVec;			// translation vector; denote t;
	tvec<double> translProjBoxVec;	// t'PB;

	/* model parameters are given by stu_system_param; */
	int		scale_num;			// the number of length scales;
	double	scale_val[DATALEN];	// the values of length scale;
	double	model_xi;			// the penalty factor;
	double	model_tau;			// the coefficient before quadratic term;
	double	model_gamma;		// the coefficient before cubic term;
	double	model_kappa;		// the coefficient before quartic term;

	/* iteration parameters; */
	double	trunc_tol;			// the maximal tolerance error for truncation;
	double	tol;				// the maximal tolerance error;
	double	step_size;			// the time steeping size;
	int		iter_max;			// the maximal iterator;
	int		print_step;			// the stepping size for printing data;
	int		save_step;			// the stepping size for saving data;
	char save_type[STRLEN];		// the type of saving data; 
								// followed by Fourier coefficient, density, plane wave; 
								// y: yes; n: no;
	
	/* parameters to optimize computational box; */
	double	opt_tol;			// the maximal tolerance energy error of optimization;
	int		opt_iter_max;		// the maximal iterator of optimization;
	int		opt_save_step;		// the stepping size for saving data in optimization;
	int		box_bbType;			// the type of BB stepping size in box optimization;
	int		box_iter_max;		// the maximal iterator in box optimization;
	double	box_step_size;		// the initial stepping size in box optimization;
	double	box_tol;			// the maximal tolerance error in box optimization;


	/* structure body 'stu_fftw_var' for calculating stable bulk phase; */
	stu_fftw_var sfftvType;
	stu_fftw_var *sfftv = &sfftvType;

	/* structure body 'stu_bnd_var' to make boundaries homogeneous; */
	stu_bnd_var sbndvType;		// in the original space;
	stu_bnd_var *sbndv = &sbndvType;

	/* structure body 'stu_bnd_var' under the common 'rotateProjBoxMat'; */
	stu_bnd_var srebndvType;	// re-representation;
	stu_bnd_var *srebndv = &srebndvType;

} stu_bulk_param;


/**
 * \struct stu_system_param;
 * \brief  Input the interface parameters;
 *
 * Input interface parameters, reading from disk file;
 */
typedef struct stu_system_param
{
	char iter_method[STRLEN];	// iteration method;
	int		print_level;		// the printing level; 0 is not displayed;

	/* model parameters; */
	int		scale_num;			// the number of length scales;
	double	scale_val[DATALEN];	// the values of length scale;
	double	model_xi;			// the penalty factor;
	double	model_tau;			// the coefficient before quadratic term;
	double	model_gamma;		// the coefficient before cubic term;
	double	model_kappa;		// the coefficient before quartic term;

	/* spatial discrete parameters; */
	int		GJP_degree;			// the degree of General Jacobi Polynomial up to GJP_degree; (Nx)
	int		LGL_num;			// the number of Legendre Gauss-Lobatto points; i.e. discrete x; (Mx)
	int		Four_num;			// the number of Fourier discrete points along each direction;

	char x_range_type[STRLEN];	// the calculation type of x_range;
	double	x_range;			// the distance between the two anchoring planes; (* 2);

	double	smooth;				// the smooth constant to connect the two bulk phases;
	double	initDist1;			// the end position of the left bulk phase when we construct initial value;
	double	initDist2;			// the end position of the rigth bulk phase when we construct initial value;

	/* iteration parameters; */
	double	tol;				// the maximal tolerance error;
	double	tolham;				// the maximal tolerance error of Hamilton;
	double	step_size;			// the time stepping size;
	double	step_min;			// the lower bound of time steeping size;
	double	step_max;			// the upper bound of time steeping size;
	int		iter_load;			// the loading iterator;
	int		iter_max;			// the maximal iterator;
	int		print_step;			// the stepping size for printing data;
	int		save_step;			// the stepping size for saving data;
	char save_type[STRLEN];		// the type of saving data; 
								// followed by Fourier coefficient, density, plane wave; 
								// y: yes; n: no;

	/* some iteration parameters for Newton-PCG method; */
	double	newton_tol;			// the maximal tolerance error to end the first-order method;
	double  newton_step_size;	// the time stepping size of Newton method;
	int		pcg_type;			// the type of pre-conditioner;
	int		pcg_iter_max;		// the maximal iteration of PCG;
	int		pcg_print_step;		// the stepping size for printing data in the calculation of PCG;
	double	pcg_delta_coeff;	// the coefficient of delta in the calculation of pre-conditioner;
	double	pcg_mu_para;		// the initial 'mu_para' for pre-conditioner;

	/* plotting parameters; */
	double	y_start;			// the start of plotting interface along y-direction;
	double	z_start;			// the start of plotting interface along z-direction;
	double	y_range;			// the plot range along y-direction;
	double	z_range;			// the plot range along z-direction for 3D phases;
	int		y_num;				// the number of discrete points along y-direction;
	int		z_num;				// the number of discrete points along z-direction;
	int		skip_half;			// 1 to skip the half of Four_num in plotting, default 0;

	/* parameters for representing error; */
	double	err_y_start;		// the start of plotting interface along y-direction;
	double	err_z_start;		// the start of plotting interface along z-direction;
	double	err_y_range;		// the plot range along y-direction;
	double	err_z_range;		// the plot range along z-direction for 3D phases;
	int		err_y_num;			// the number of discrete points along y-direction;
	int		err_z_num;			// the number of discrete points along z-direction;

	/* parameter to recover; */
	double recoverTOL;			// the maximal spectra mode;

	/* common 'rotateProjBoxMat'; */
	int		searchReg;			// the searching range in the computation of common 'rotateProjBoxMat';
	int		adjustReg;			// the adjust range in the test of common 'rotateProjBoxMat';
	int		dimRePhy;			// the number of row of 'rotateProjBoxMat';
	int		dimReCpt;			// the number of column of 'rotateProjBoxMat';
	tvec<ptrdiff_t> NCpt;		// space grid of the Fourier direction;
	int		cplxReDofs;			// i.e. pow(Four_num, dimReCplx);
	char com_projmat_way[STRLEN];	// the way to get the common 'rotateProjBoxMat';
	tmat<double> rotateProjBoxMat;	// the common 'rotateProjBoxMat';
	tmat<int>	coeffmat;			// the coefficient matrix containing two bulk phases;

	/* iteration variables; */
	int		massFlag;			// the flag of mass conservation;
								// 0: not forced, only show mass; might be good for some case, such as 12fold symmetric tilt grain boundary;
								// 1: forced mass conservation, default value;
	tCCSmat<double>	   interact_grad;		// the interaction potential operator;
	tvec<fftw_complex> interact_bnd_grad;	// the interaction potential operator on boundaries;
	tCCSmat<double>	   iter_matrix;			// the sparse matrix for iteration;
	tvec<fftw_complex> iter_rho_rhs;		// the right term for iteration; u^{n}/step_size;
	tvec<fftw_complex> iter_entropy_rhs;	// the right term for iteration; f(phi) + xi^2(\Delta+1)^2 p;
	tvec<fftw_complex> grad_err;			// the gradient error;
	tvec<fftw_complex> hessian;				// the Hessian matrix;

	/* structure body 'stu_fftw_var' for calculating the interface system; */
	stu_fftw_var sfftvType;
	stu_fftw_var *sfftv = &sfftvType;

	/* structure body 'stu_GJP_var' for projecting bulk phase; */
	stu_GJP_var sGJPvType;
	stu_GJP_var *sGJPv = &sGJPvType;

	/* structure body 'stu_bnd_var' for calculating the interface system; */
	stu_bnd_var scbndvType;		// under the common 'rotateProjBoxMat';
	stu_bnd_var *scbndv = &scbndvType;

	/* massVec; */
	tvec<fftw_complex> massVec;				// for mass conservation;

} stu_system_param;

/* --------------------		bulk and system		--------------------*/

#endif
