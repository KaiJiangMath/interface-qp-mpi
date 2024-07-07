#ifndef __functs_h
#define __functs_h

#include "Data.h"
#include "Head.h"

/*----------------- In file: FunctCollect.cpp -----------------*/

/* Calculation of the stable interface structure; */
extern double fn_main_stable_interface ( );

/* Calculation of the stable interface structure
 * by loading existing files; */
extern double fn_main_load_stable_interface ( );

/* Calculation of the stable bulk phase; */
extern double fn_main_stable_bulk ( );

/* Calculation of the stable bulk phase by a parallel code; */
extern double fn_main_stable_bulk_parallel ( );

/* Projection of the stable bulk phase; */
extern void fn_main_project_bulk  ( );

/* Transform rhoJCplx (in the space of Fourier and GJPs) to 
 * the common space (Fourier and GJPs); */
extern void fn_main_common_bulk	  ( );

/* Calculation of the common projection matrix; */
extern void fn_main_common_rotateProjBoxMat ( );

/* Display density according to the file about 'rhoJCplx'; */
extern int fn_main_disp_system_density ( );

/* Collect all data for interface system; */
extern int fn_main_collect_system_data ( );

/* Recover density profile by the screened 'rhoJCplx'; */
extern int fn_main_recover_system_density ( );

/* Generate parameter files in batches; */
extern void fn_param_batch ( );

/* char to str; */
string charToStr (char *contentChar);

/*----------------- In file: FunctCollect.cpp -----------------*/

/*----------------- In file: BasicOperators.cpp -----------------*/

/* Maximal absolute value; */
extern double normRealInfty ( double *src, int n );

/* Maximal absolute value of an integer matrix; */
extern int fn_max_abs ( int **src, int row, int col );

/* Define the sort rule; */
extern bool fn_compare (stu_sort a, stu_sort b);

/* Define the descend sort rule; */
extern bool fn_compare_descend (stu_sort a, stu_sort b);

/* Set zero for a complex value; */
extern void fn_complex_setZero					(	fftw_complex			rslt );

/* The absolute value for a complex value; */
extern double fn_complex_abs					(	fftw_complex			src );

/* Complex F1 multiply complex F2; */
extern void fn_complex_multiply					(	fftw_complex			F1, 
													fftw_complex			F2,
													fftw_complex			rslt );
/* complex F1 divides complex F2; */
extern void fn_complex_divide					(	fftw_complex			F1, 
													fftw_complex			F2,
													fftw_complex			rslt );
/* save 'tvec' vector; */
extern void fn_tvec_save_int					(	tvec<int>				src, 
														char				filename[] );

/* print 'tvec' vector; */
extern void fn_tvec_print_int					(	tvec<int>				src );

/* Set zero of a 'tvec' vector; */
extern void fn_tvec_setZero_complex				(	tvec<fftw_complex>		src );

/* Transpose a matrix whose is straightened as a 'tvec' vector; */
extern void fn_tvec_trans_complex				(	tvec<fftw_complex>		src,
													tvec<fftw_complex>		rslt,
															int				row,
															int				col );

/* Maximal absolute value of a complex 'tvec' vector; */
extern double fn_tvec_maxAbs_complex			(	tvec<fftw_complex>		src );

/* Norm 2 of a complex 'tvec' vector; */
extern double fn_tvec_norm_complex				( tvec<fftw_complex>		src );

/* The adding result of a complex 'tvec' vector and a double value; */
extern void fn_tvec_constAdd_complex			(	tvec<fftw_complex>		rslt,
														double				a );

/* The multiplying result of a complex 'tvec' vector and a double value; */
extern void fn_tvec_constMultiply_complex		(	tvec<fftw_complex>		src,
														double				a );

/* Result of adding two complex 'tvec' vectors; */
extern void fn_tvec_add_complex					(	tvec<fftw_complex>		dst, 
													tvec<fftw_complex>		src, 
														double				dstCoeff, 
														double				srcCoeff );

/* The dot multiplying result of two complex 'tvec' vectors; */
extern void fn_tvec_dotMultiply_complex			(	tvec<fftw_complex>		F1,
													tvec<fftw_complex>		F2,
													tvec<fftw_complex>		rslt );

/* The summation of the dot multiplying result of two complex 'tvec' vectors; */
extern void fn_tvec_dotMultiplySum_complex		(	tvec<fftw_complex>		F1,
													tvec<fftw_complex>		F2,
															int				len,
														fftw_complex		rslt );

/* Save 'tvec' vector; */
extern void fn_tvec_save_complex				(	tvec<fftw_complex>		src, 
															char			filename[] );

/* Print 'tvec' vector; */
extern void fn_tvec_print_complex				(	tvec<fftw_complex>		src );

/* Save 'tmat' matrix; */
extern void fn_tmat_save_int					(	tmat<int>				src, 
														char				filename[] );
/* Print 'tmat' matrix; */
extern void fn_tmat_print_int					(	tmat<int>				src );

/* Set zero of 'tmat' matrix; */
extern void fn_tmat_setZero_complex				(	tmat<fftw_complex>		src );

/* Save 'tmat' matrix; */
extern void fn_tmat_save_complex				(	tmat<fftw_complex>		src, 
															char			filename[] );

/* Print 'tmat' matrix; */
extern void fn_tmat_print_complex				(	tmat<fftw_complex>		src );

/* calculate the det of matrix; */
extern double det								(		double				**D, 
														int					n );
/* calculate the inverse matrix; */
extern void matrix_inverse						(		double				**a,
														double				**b,
														int					N	);

/*----------------- In file: BasicOperators.cpp -----------------*/

/*----------------- In file: SparseOperators.cpp -----------------*/

/* Save sparse matrix of int type in CSR format; */
extern void fn_tCSRmat_save_int					(	tCSRmat<int>			src, 
														char				filename[] );

/* Save sparse matrix of int type in CCS format; */
extern void fn_tCCSmat_save_int					(	tCCSmat<int>			src, 
														char				filename[] );

/* Generate block diagonal sparse matrix; */
extern tCCSmat<double> fn_block_diag_dCCSmat	( int row, int col,	int	d );

/* calculate the tensor product of diag(diagSrc) and dccsSrc; */
extern tCCSmat<double> fn_tensor_diag_dCCSmat	(	tvec<double>		diagSrc,
													tCCSmat<double>		dccsSrc,
														int				order );

/* Calculate the product of dccsSrc *  diagSrc;
 * 'diagSrc':	a diagonal CCS matrix; */
extern tCCSmat<double> fn_multiply_diag_dCCSmat	(	tCCSmat<double>		diagSrc,
													tCCSmat<double>		dccsSrc );

/* calculate the result of coeff * dccsSrc; */
extern tCCSmat<double> fn_const_multiply_dCCSmat (	tCCSmat<double>		dccsSrc,
														double			coeff );

/* Obtain the diagonal elements of dccsSrc and add coeff*I on them; */
extern tCCSmat<double> fn_obt_diag_add_dCCSmat	(	tCCSmat<double>		dccsSrc,
														double			coeff1,
														double			coeff2 );

/* Calculate the result of dccsSrc + coeff*I; */
extern tCCSmat<double> fn_diag_add_dCCSmat		(	tCCSmat<double>		dccsSrc,
														double			coeff1,
														double			coeff2 );

/* Calculate the result of dccsSrc + coeff*blockI; */
extern tCCSmat<double> fn_block_diag_add_dCCSmat (	tCCSmat<double>		dccsSrc,
														double			coeff1,
														double			coeff2,
														int				d );

/* Transpose the CCS 'dccsSrc'; */
extern tCCSmat<double> fn_trans_dCCSmat			 (	tCCSmat<double>		dccsSrc );

/* Transpose the CCS 'dccsSrc' by a faster way; */
extern tCCSmat<double> fn_fast_trans_dCCSmat	 (	tCCSmat<double>		dccsSrc );

/* calculate the result of dccsSrc * rhoJCplx; */
extern void fn_cvec_multiply_dCCSmat			 (	tCCSmat<double>		dccsSrcTrans,
													tvec<fftw_complex>	rhoJCplx,
													tvec<fftw_complex>	rslt );

/* calculate the result of dccsSrc1 add dccsSrc2; */
extern tCCSmat<double> fn_add_dCCSmat			 (	tCCSmat<double>		dccsSrc1,
													tCCSmat<double>		dccsSrc2,
														double			coeff1,
														double			coeff2 );

/* Solve Ax = b by umfpack.h; ( double type ); */
extern int fn_umfpack_solver					 (	tCCSmat<double>		A, 
													tvec<double>		b, 
													tvec<double>		x );

/* solve Ax = b by umfpack.h; ( complex type ); */
extern int fn_umfpack_complex_solver			 (	tCCSmat<double>		A, 
													tvec<fftw_complex>	b, 
													tvec<fftw_complex>	x );

/*----------------- In file: SparseOperators.cpp -----------------*/

/*----------------- In file: AuxInput.cpp -----------------*/

/* The main code of input parameters; */
extern void fn_param_input					(	stu_bulk_param		*sbulkparam1,
												stu_bulk_param		*sbulkparam2,
												stu_system_param	*ssysparam );

/* Default input parameters about bulk phases; */
extern void fn_bulk_nature_param_input_init (	stu_bulk_param		*sbulkparam );

/* Default input parameters about bulk phases; */
extern void fn_stu_bulk_param_input_init	(	stu_bulk_param		*sbulkparam );

/* Default input parameters about the interface; */
extern void fn_stu_system_param_input_init  (	stu_system_param	*ssysparam );

/* Read input parameters from disk file; */
extern void fn_bulk_nature_param_input		(	stu_bulk_param		*sbulkparam );

/* Read input parameters from disk file; */
extern void fn_stu_bulk_param_input			(	const char			*fname,
												stu_bulk_param		*sbulkparam );

/* Set the projection matrix; */
extern void fn_obt_projection_matrix		(	stu_bulk_param		*sbulkparam );

/* Generate the rotation matrix; */
extern void fn_obt_rotation_matrix			(	const char			*fname,
												stu_bulk_param		*sbulkparam );

/* Set the translation vector; */
extern void fn_obt_translation_vector		(	stu_bulk_param		*sbulkparam );

/* Calculate the products of rotateMat, projMat, and rcpBox; */
extern void fn_obt_rotateProjBox_matrix		(	stu_bulk_param		*sbulkparam );

/* Calculate the product of translVec, projMat, rcpBox; */
extern void fn_obt_translProjBox_matrix		(	stu_bulk_param		*sbulkparam );

/* Read input parameters from disk file; */
extern void fn_stu_system_param_input		(	const char			*fname,
												stu_system_param	*ssysparam );

/* Copy data from stu_system_param to stu_bulk_param; */
extern void fn_param_copy_system_to_bulk	(	stu_system_param	*ssysparam,
												stu_bulk_param		*sbulkparam );

/* Save input parameters; */
extern void fn_param_save					(	stu_bulk_param		*sbulkparam1,
												stu_bulk_param		*sbulkparam2,
												stu_system_param	*ssysparam,
												const char			*key );

/* Save input parameters about bulk phases; */
extern void fn_stu_bulk_param_save			(	stu_bulk_param		*sbulkparam,
													FILE			*fp );

/* print input parameters about bulk phases; */
extern void fn_stu_bulk_param_print			(	stu_bulk_param		*sbulkparam );

/* print input parameters about interface system; */
extern void fn_stu_system_param_print		(	stu_system_param	*ssysparam );

/* save matrix for parameters; */
extern void fn_matrix_save					(	tmat<double>		src,
													FILE			*fp );

/* print matrix for parameters; */
extern void fn_matrix_print					(	tmat<double>		src );

/* Adjust model parameters to match different models; */
extern void fn_adjust_bulk_model_param		(	stu_bulk_param		*sbulkparam );

/* Adjust model parameters to match different models; */
extern void fn_adjust_system_model_param	(	stu_system_param	*ssysparam	);

/* Read 'main_type' from the parameter file; */
extern void fn_obt_main_type (  );

/* Obtain 'rsltDir' from the parameter file; */
extern void fn_obt_rsltDir (  );

/* Read 'com_projmat' from the parameter file; */
extern void fn_obt_com_projmat				( stu_system_param		*ssysparam,
												char			paraFile[FILELEN] );

/* Read double tmat matrix; */
extern tmat<double>		  fn_tmat_read_double	( char filename[] );

/* Read double tvec matrix; */
extern tvec<double>		  fn_tvec_read_double	( char filename[] );

/* Read fftw_complex tvec matrix; */
extern tvec<fftw_complex> fn_tvec_read_complex	( char filename[] );

/*----------------- In file: AuxInput.cpp -----------------*/

/*----------------- In file: FftwToolkit.cpp -----------------*/

/* Obtain Fourier k for a serial code; */
extern void fn_serial_obt_kIndex			(	tmat<int>			kspace,
												tvec<int>			NCpt );

/* Obtain global indexes of special positions (N/2); */
extern tvec<int> fn_serial_obt_half_gIndex	(	tvec<int>			NCpt,
													int				cplxDofs );

/* Transform Fourier k into the global index; */
extern int fn_serial_kIndex_to_gIndex		(		int				*k, 
												tvec<int>			Ndof );

/* Calculation of convolution of pow('src', order) for a serial code; */
extern void fn_serial_convolution			(	tvec<fftw_complex>	src,
												tvec<fftw_complex>	orig, 
												stu_fftw_var		*sfftv,
													int				order );


/* Obtain Fourier k for a parallel code; */
extern void fn_obt_kIndex					(	tmat<int>			kspace,
												tvec<ptrdiff_t>		NCpt,
													int				*displs,
													int				local_len );

/* Obtain global indexes of special positions (N/2); */
extern tvec<int> fn_obt_half_gIndex			(	tvec<ptrdiff_t>		NCpt,
													int				cplxDofs );

/* Transform Fourier k into the global index; */
extern int fn_kIndex_to_gIndex				(		int				*k, 
												tvec<ptrdiff_t>		Ndof );

/* Calculation of convolution of pow('src', order) for a parallel code; */
extern void fn_convolution					(	tvec<fftw_complex>	src,
												tvec<fftw_complex>	orig, 
												stu_fftw_var		*sfftv,
													int				*displs,
													int				vecLen,
													int				cplxDofs,
													int				order );

/* Calculation of convolution of 'src1' and 'src2'; */
extern void fn_convolution_general			(	tvec<fftw_complex>	src1,
												tvec<fftw_complex>	src2,
												tvec<fftw_complex>	orig, 
												stu_fftw_var		*sfftv,
													int				cplxDofs );

/*----------------- In file: FftwToolkit.cpp -----------------*/

/*----------------- In file: StableBulkPhase.cpp -----------------*/

/* Update bulk phase to obtain the stable state; */
extern double fn_obt_stable_bulk_phase	(	stu_bulk_param	*sbulkparam );

/* Allocation memory for calculating stable bulk phase; */
extern void fn_bulk_memory_allocation	(	stu_bulk_param	*sbulkparam );

/* Get projection plane; i.e. R'PBk; */
extern void fn_get_projection_plane		(	stu_bulk_param	*sbulkparam );

/* Get Gsquare; i.e. |R'PBk|^2; based on 'fn_get_projection_plane'; */
extern void fn_get_Gsquare				(	stu_bulk_param	*sbulkparam );

/* Get the initial state of bulk phase; */
extern void	fn_get_init_bulk_phase		(	stu_bulk_param	*sbulkparam );

/* Update bulk phase and obtain the stable bulk phase; */
extern double fn_update_bulk_phase		(	stu_bulk_param	*sbulkparam,
												FILE		*fdata,
												int			opt_iter );

/* Calculate the free energy of bulk phase; */
extern void fn_calc_bulk_energy			(	stu_bulk_param	*sbulkparam,
											tvec<double>	energy );

/* The integration of functions about saving data;
 *	save Fourier coefficients, densities, plane wave;
 */
extern void fn_save_bulk_phase			(	stu_bulk_param	*sbulkparam,
												int			opt_iter,
												int			iterator );

/* Optimize the computational box; */
extern void fn_opt_rcpBox				(	stu_bulk_param	*sbulkparam );

/* Take translation in the bulk phase; */
extern void fn_translate_bulk_phase		(	stu_bulk_param	*sbulkparam );

/* Truncation to reduce the computational cost; */
extern void fn_truncation_bulk_phase	(	stu_bulk_param	*sbulkparam,
												double		hamilton	);

/* Obtain the serial number of model parameters for bulk phase 
 * or save the model parameters for bulk phase; */
extern int fn_paraOrder_bulk_phase		(	stu_bulk_param	*sbulkparam );
	
/* Save the truncation data for bulk phase; */
extern void fn_save_trunc_bulk_phase	(	stu_bulk_param	*sbulkparam,
												double		hamilton,
												bool		ischeck );

/* Load the truncation data for bulk phase; */
extern double fn_load_trunc_bulk_phase	(	stu_bulk_param	*sbulkparam	);

/*----------------- In file: StableBulkPhase.cpp -----------------*/

/*----------------- In file: StableBulkParallel.cpp -----------------*/

/* Update bulk phase to obtain the stable state; */
extern double fn_obt_stable_bulk_parallel		(	stu_bulk_param	*sbulkparam );

/* Allocation memory for calculating stable bulk phase; */
extern void fn_bulk_memory_allocation_parallel	(	stu_bulk_param	*sbulkparam );

/* Get projection plane; i.e. R'PBk; */
extern void fn_get_projection_plane_parallel	(	stu_bulk_param	*sbulkparam );

/* Get Gsquare; i.e. |R'PBk|^2; based on 'fn_get_projection_plane'; */
extern void fn_get_Gsquare_parallel				(	stu_bulk_param	*sbulkparam );

/* Get the initial state of bulk phase; */
extern void	fn_get_init_bulk_parallel			(	stu_bulk_param	*sbulkparam );

/* Update bulk phase and obtain the stable bulk phase; */
extern double fn_update_bulk_parallel			(	stu_bulk_param	*sbulkparam,
														FILE		*fdata,
														int			opt_iter );

/* Calculate the free energy of bulk phase; */
extern void fn_calc_bulk_energy_parallel		(	stu_bulk_param	*sbulkparam,
													tvec<double>	energy );

/* The integration of functions about saving data;
 *	save Fourier coefficients, densities, plane wave;
 */
extern void fn_save_bulk_parallel				(	stu_bulk_param	*sbulkparam,
														int			opt_iter,
														int			iterator );

/* Optimize the computational box; */
extern void fn_opt_rcpBox_parallel				(	stu_bulk_param	*sbulkparam );

/* Take translation in the bulk phase; */
extern void fn_translate_bulk_parallel			(	stu_bulk_param	*sbulkparam );

/* Truncation to reduce the computational cost; */
extern void fn_truncation_bulk_parallel			(	stu_bulk_param	*sbulkparam,
														double		hamilton	);

/* Obtain the serial number of model parameters for bulk phase 
 * or save the model parameters for bulk phase; */
extern int fn_paraOrder_bulk_parallel			(	stu_bulk_param	*sbulkparam );
	
/* Save the truncation data for bulk phase; */
extern void fn_save_trunc_bulk_parallel			(	stu_bulk_param	*sbulkparam,
														double		hamilton,
														bool		ischeck );

/* Load the truncation data for bulk phase; */
extern double fn_load_trunc_bulk_parallel		(	stu_bulk_param	*sbulkparam	);

/*----------------- In file: StableBulkParallel.cpp -----------------*/

/*----------------- In file: GenJacPoly.cpp -----------------*/

extern void fn_obt_system_gen_Jac_poly		(	stu_GJP_var		*sGJPv );

extern void fn_GJP_memory_allocation		(	stu_GJP_var		*sGJPv );

extern void fn_obt_gen_Jac_poly				(	stu_GJP_var		*sGJPv,
												tmat<double>	dsJJ, 
													int			order );

extern tCCSmat<double> fn_innSMat_gen_Jac_poly		(	stu_GJP_var		*sGJPv,
														tmat<double>	dsJJ, 
															int			order,
															double		TOL );

extern tCCSmat<double> fn_diag_innSMat_gen_Jac_poly	(	stu_GJP_var		*sGJPv,
														tmat<double>	dsJJ, 
															int			order,
															double		TOL);

extern void fn_obt_Leg_poly					(	tvec<double>	x,
												tvec<double>	y, 
													int			n );

extern void fn_obt_first_deri_Leg_poly		(	tvec<double>	x,
												tvec<double>	y, 
												tvec<double>	dy, 
													int			n );

extern void fn_obt_Leg_Gau_Lob				(	tvec<double>	x,
												tvec<double>	w,
													int			n );

extern void fn_obt_Jac_poly					(	tvec<double>	x,
												tmat<double>	y, 
													double		jAlpha,
													double		jBeta );

extern void fn_obt_first_deri_Jac_poly		(	tvec<double>	x,
												tmat<double>	y,
												tmat<double>	dy, 
													double		jAlpha,
													double		jBeta );

/*----------------- In file: GenJacPoly.cpp -----------------*/

/*----------------- In file: ProjFourGJP.cpp -----------------*/

extern void fn_project_Fourier_GJP		(stu_bulk_param		*sbulkparam,
										stu_system_param	*ssysparam);

extern void fn_sbndv_memory_allocation	(stu_bulk_param		*sbulkparam, 
											stu_GJP_var		*sGJPv);

/* introduce polynomial to make boundary conditions homogeneous; */
extern void fn_obt_poly_homo_bnd		(	stu_bulk_param		*sbulkparam,
											stu_GJP_var			*sGJPv, 
											tvec<fftw_complex>	bndlr,
												int				order );

/* construct a complex polynomial satisfying boundary conditions; */
extern void fn_interpl_complex_bnd		(	tvec<double>		x,
											tvec<fftw_complex>	fun,
											tvec<fftw_complex>	BC, 
												double			L,
												int				order );

/* construct a real polynomial satisfying boundary conditions; */
extern void fn_interpl_bnd				(	tvec<double>		x,
											tvec<double>		fun,
											tvec<double>		BC, 
												double			L,
												int				order );

/* project rhoCplx onto the GJP space; */
extern void fn_project_GJP				(	stu_bulk_param		*sbulkparam,
											stu_GJP_var			*sGJPv );

/* compute the mass of left and right bulk phases after projecting; */
extern void fn_obt_bulk_mass_project	(	stu_bulk_param		*sbulkparam,
											stu_system_param	*ssysparam,
											fftw_complex			mass	);

/*----------------- In file: ProjFourGJP.cpp -----------------*/

/*----------------- In file: CommonRotateProjBoxMat.cpp -----------------*/

extern int fn_obt_commom_rotateProjBoxMat (	 stu_bulk_param		*sbulkparam1,
											 stu_bulk_param		*sbulkparam2,
											 stu_system_param	*ssysparam);

extern tmat<double> fn_matrix_connect			(	tmat<double>	xmat, 
													tmat<double>	ymat );

extern tmat<double> fn_matrix_rational_dependent (	tmat<double>	xymat,
													tmat<int>		coeffmat,
													double			tol,
													int				int_reg,
													int				adjustReg );

extern void fn_obt_int_coeff (int **intSpace, int dim, int N, int cplxDofs);

extern void fn_opt_int_coeff (int **intSpace, int dim, int N, int cplxDofs);

extern int  fn_check_vector_rational_dependent (double *x,	double *y,	double *z,
												int		n,	double tol,	int	   int_reg);

extern int  fn_check_number_rational_dependent (double x,		double y,
												double tol,		int	   int_reg);

/*----------------- In file: CommonRotateProjBoxMat.cpp -----------------*/

/*----------------- In file: CommonFourGJP.cpp -----------------*/

/* Transform rhoJCplx of the two bulk phases in a common space;	 */
extern void fn_common_Fourier_GJP			( stu_bulk_param		*sbulkparam1,
											  stu_bulk_param		*sbulkparam2,
											  stu_system_param		*ssysparam );

/* obtain projPlane on each process; */
extern void fn_get_projPlaneAll				( stu_bulk_param		*sbulkparam );

/* obtain projPlane on each process; */
extern void fn_get_system_projection_plane  ( stu_system_param		*ssysparam);

/* Get projection plane in the interface system; i.e. R'PBk; */
extern void fn_total_re_represent_common	( stu_bulk_param		*sbulkparam,
											  stu_system_param		*ssysparam);

extern void fn_common_memory_alloc			( stu_bnd_var			*sbndv,
											  stu_bnd_var			*sobjbndv,
												int					cplxReDofs,
												const char			*memType);

extern void fn_re_represent_common		(	stu_bulk_param		*sbulkparam,
											stu_system_param	*ssysparam,
											tvec<fftw_complex>	origBulk,
											tvec<fftw_complex>	reComBulk,
											const char			*saveStr );

extern void fn_new_re_represent_common	(	stu_bulk_param		*sbulkparam,
											stu_system_param	*ssysparam,
											tvec<fftw_complex>	origBulk,
											tvec<fftw_complex>	reComBulk,
											const char			*saveStr );

/* calculate the error between 'bulk1_FGJP_density' and 'bulk1_FGJPre_density'; */
extern void fn_proj_common_error		(	stu_bulk_param		*sbulkparam,
											stu_system_param	*ssysparam );

/* read file for 2D data; */
extern int fn_read_data_2d				(		char			fname[FILELEN],		
											tvec<double>		x,
											tvec<double>		y,
											tvec<double>		data );

/* read file for 3D data; */
extern int fn_read_data_3d				(		char			fname[FILELEN],		
											tvec<double>		x,
											tvec<double>		y,
											tvec<double>		z,
											tvec<double>		data );

/* compute the mass of left and right bulk phases; */
extern void fn_obt_bulk_mass_common		(	stu_bulk_param		*sbulkparam,
											stu_system_param	*ssysparam,
											fftw_complex			mass	);

/*----------------- In file: CommonFourGJP.cpp -----------------*/

/*----------------- In file: SysInitVal.cpp -----------------*/

extern void fn_system_initial_value			(stu_bulk_param		*sbulkparam1,
										    stu_bulk_param		*sbulkparam2,
										    stu_system_param	*ssysparam);

extern void fn_connect_init_value			(stu_bulk_param		*sbulkparam1,
											stu_bulk_param		*sbulkparam2,
											stu_system_param	*ssysparam);

extern void fn_connect_rebnd			(	stu_bnd_var			*srebndv1,
											stu_bnd_var			*srebndv2,
											stu_system_param	*ssysparam);

/* generate the connection function; */
extern void fn_obt_connect_fun			(		double			xVal,
												double			smooth,
												double			initDist, 
												int				sflag,
											tvec<double>		masklr );

/* save 'd0bnd', 'd2bnd', ... of the structure body 'scbndv'; */
extern void fn_save_cbnd				(	tvec<fftw_complex>	dsbnd,
												int				order,
												int				xlen,
												int				cplxReDofs );

extern void fn_total_bulk_memory_free	(	stu_bulk_param		*sbulkparam);

/*----------------- In file: SysInitVal.cpp -----------------*/

/*----------------- In file: IterPrepare.cpp -----------------*/

extern void fn_iter_prepare			(	stu_system_param	*ssysparam,
											bool			isTest );

extern void fn_obt_system_Gsquare	(	stu_system_param	*ssysparam );

extern void fn_obt_interact_grad	(	stu_system_param	*ssysparam,
											bool			isTest );

extern void fn_obt_interact_bnd_grad(	stu_system_param	*ssysparam,
											bool			isTest );

extern void fn_obt_iter_matrix		(	stu_system_param	*ssysparam,
											double			step_size,
											int				iterator,
											bool			isTest );

/* Project 'rhoJCplx' to Fourier space; */
extern tvec<fftw_complex> fn_back_Four_space (stu_system_param	 *ssysparam,
											  tvec<fftw_complex> rhoJCplx );

/* Mass conservation by using 'fn_obt_mass' and 'fn_main_mass' many
 *			times to reach the desired precision;
 */
extern void fn_reach_mass				(	stu_system_param		*ssysparam,
											tvec<fftw_complex>		rhoJCplx,
											tvec<fftw_complex>		rhoJCplxMass,
											fftw_complex			mass,
											double					desir_mass,
											double					TOL	);

/* Calculate the integral of 'rhoJCplx' (check mass conservation); */
extern void fn_obt_mass					(	stu_system_param		*ssysparam,
											tvec<fftw_complex>		rhoJCplx,
											fftw_complex			mass );

/* Maintain the system mass; */
extern void fn_maintain_mass			(	stu_system_param		*ssysparam,
											tvec<fftw_complex>		rhoJCplx,
											tvec<fftw_complex>		rhoJCplxMass,
											fftw_complex			mass	);

extern void fn_system_fftw_memory_allocation (stu_system_param	*ssysparam);

/* Calculate the free energy and gradient term of the interface system; */
extern void fn_calc_system_energy	(	stu_system_param	*ssysparam,
										tvec<fftw_complex>	rhoJCplx,
										tvec<double>		energy, 
											bool			withGrad,
											bool			isTest );

/* Calculate the right term for iteration; */
extern void fn_obt_iter_rhs			(	stu_system_param	*ssysparam,
										tvec<fftw_complex>	rhoJCplx,
										tvec<fftw_complex>	gradient,
											int				iterator,
											bool			isTest );

/* Obtain 'iter_rhs_all' containing all informations of 'iter_rhs'; */
extern void fn_obt_iter_rhs_all		(	stu_system_param	*ssysparam,
										tvec<fftw_complex>	iter_rhs,
										tvec<fftw_complex>	iter_rhs_all );

/* Calculate the rhoJCplx weighted difference between two adjacent steps; */
extern void fn_diff_weight_rhoJCplx	(	stu_system_param	*ssysparam,
										tvec<fftw_complex>	rhoJCplxNew,
										tvec<fftw_complex>	rhoJCplx,
										tvec<fftw_complex>	rhoJCplxDiff );

/* Calculate the gradient weighted difference between two adjacent steps; */
extern void fn_diff_weight_gradient	(	stu_system_param	*ssysparam,
										tvec<fftw_complex>	gradientNew,
										tvec<fftw_complex>	gradient,
										tvec<fftw_complex>	gradientDiff );

/* Calculate the Hessian matrix of the interface system; */
extern void fn_calc_system_hessian	(	stu_system_param	*ssysparam,
										tvec<fftw_complex>	rhoJCplx,
											int				iterator,
											bool			isTest );

/* Calculate the regularized Hessian matrix of the interface system; */
extern void fn_calc_system_reg_hessian(	stu_system_param	*ssysparam,
										tvec<fftw_complex>	direction,
										tvec<fftw_complex>	regHessian,
										tCCSmat<double>		interact_grad_trans,
											double			mu );

/* pre-conditioner; */
extern void fn_calc_system_preconditioner ( stu_system_param	*ssysparam,
											tvec<fftw_complex>	direction,
											tvec<fftw_complex>	preconditioner,
												double			mu,
												double			delta,
												int				iterator,
												bool			isTest );

/* Calculate Newton direction directly; */
extern void fn_calc_system_newton_direction(stu_system_param	*ssysparam,
											tvec<fftw_complex>	direction,
											tCCSmat<double>		interact_grad_trans );

/* Preconditioned conjugate gradient (PCG) method; */
extern int fn_calc_system_PCG		(	stu_system_param		*ssysparam,
										tvec<fftw_complex>		x,
										tCCSmat<double>			interact_grad_trans,
											double				mu,
											double				delta,
											int					iterator,
											bool				isTest );

/* Recover 'rhoJCplx' after screening spectra; */
extern int fn_recover_rhoJCplx		(	stu_system_param		*ssysparam,
										tvec<fftw_complex>		rhoJCplx );

/* calculate the error between density profiles at step0 and step1; */
extern void fn_recover_error		(	stu_system_param	*ssysparam,
												int			step0,
												int			step1	);

/* The integration of functions about saving data; */
extern void fn_save_system_phase	(	stu_system_param	*ssysparam,
											int				iterator );

/*----------------- In file: IterPrepare.cpp -----------------*/

/*----------------- In file: IterMethod.cpp -----------------*/

/* Choose iteration method; */
extern double fn_choose_method		( stu_system_param		*ssysparam,
											bool			isTest,
											int				iter_load	);

/* Semi-implicit scheme (SIS); */
extern int fn_sis_method			(	stu_system_param	*ssysparam,
										tvec<double>		energy,
											bool			isTest,
											int				iter_load	);

/* Semi-implicit scheme (SIS) with random disturbance; */
extern int fn_sis_rand_method		(	stu_system_param	*ssysparam,
										tvec<double>		energy,
											bool			isTest,
											int				iter_load	);

/* Accelerated proximal gradient method (APG); */
extern int fn_apg_method			(	stu_system_param	*ssysparam,
										tvec<double>		energy,
											double			tol,
											bool			isTest,
											int				iter_load	);
/* Newton PCG method; */
extern int fn_newton_method			(	stu_system_param	*ssysparam,
										tvec<double>		energy,
											int				iterator,
											bool			isTest	);

/*----------------- In file: IterMethod.cpp -----------------*/

/*----------------- In file: DisplayResults.cpp -----------------*/

/* display densities of bulk phase; */
extern void fn_disp_bulk_density	(	tvec<fftw_complex>	rhoCplx, 
										 stu_bulk_param		*sbulkparam, 
											int				opt_iter,
											int				step );

extern bool myComp		(mySortVec a,	mySortVec b);

extern bool scaleComp	(mySortWave a,	mySortWave b);

/* display Fourier coefficients of bulk phase; */
extern void fn_disp_Fourier_coeff	(	tvec<fftw_complex>	src, 
										stu_bulk_param		*sbulkparam, 
											int				opt_iter,
											int				step );

/* display plane waves of bulk phase; */
extern void fn_disp_bulk_plane_wave	(	tvec<fftw_complex>	src, 
										stu_bulk_param		*sbulkparam, 
											int				opt_iter,
											int				step );

/* display densities of bulk phase by a parallel code; */
extern void fn_disp_bulk_density_parallel	(	tvec<fftw_complex>	rhoCplx, 
												 stu_bulk_param		*sbulkparam, 
													int				opt_iter,
													int				step );

/* display Fourier coefficients of bulk phase by a parallel code; */
extern void fn_disp_Fourier_coeff_parallel	(	tvec<fftw_complex>	src, 
												stu_bulk_param		*sbulkparam, 
													int				opt_iter,
													int				step );

/* display plane waves of bulk phase by a parallel code; */
extern void fn_disp_bulk_plane_wave_parallel (	tvec<fftw_complex>	src, 
												stu_bulk_param		*sbulkparam, 
													int				opt_iter,
													int				step );

/* display density of bulk phases after projecting to the space of GJPs; */
extern void fn_disp_bulk_proj_density		(	tvec<fftw_complex>	rhoJCplx, 
												stu_bulk_param		*sbulkparam, 
												stu_system_param	*ssysparam );

/* display density of bulk phases in the common space; */
extern void fn_disp_bulk_common_density		(	tvec<fftw_complex>	rhoReJCplx, 
												stu_bulk_param		*sbulkparam, 
												stu_system_param	*ssysparam );


/* display density of the interfce system; */
extern void fn_disp_system_density			(	tvec<fftw_complex>	rhoJCplx, 
												stu_system_param	*ssysparam,
													int				step );

/*----------------- In file: DisplayResults.cpp -----------------*/

/*----------------- In file: MemFree.cpp -----------------*/

extern void fn_bulk_memory_free			(stu_bulk_param		*sbulkparam,
										const char			*memType);

extern void fn_bnd_memory_free			(stu_bnd_var		*sbndv,
										const char			*memType);

extern void fn_system_memory_free		(stu_system_param	*ssysparam);

extern void fn_GJP_memory_free			(stu_GJP_var		*sGJPv);

/*----------------- In file: MemFree.cpp -----------------*/
#endif
