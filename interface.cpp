#include "Head.h"
#include "Data.h"
#include "DataOperators.h"
#include "Mytimer.h"
#include "functs.h"
#include "umfpack.h"

int main(int argc, char* argv[])
{
	/* get the number of process and the rank of this process; */
	//MPI_Init(NULL, NULL);
	MPI_Init ( &argc, &argv );
	MPI_Comm_size( MPI_COMM_WORLD, &nprocs );
	MPI_Comm_rank( MPI_COMM_WORLD, &myrank );

	if ( myrank == 0 )
		printf("\n\n ====================================   START PROGRAM  ======================================\n");

	pid_t pid = getpid();	// getpid: self; getppid: parent; getpgid: group;
	printf("\t\t myrank : %d  \t  program PID : %d\n\n", myrank, pid);
	MPI_Barrier ( MPI_COMM_WORLD );	
	
	mytimer_t timer;
	timer.reset();
	timer.start();
	if ( myrank == 0  )
		mkdir("./plan/", 0755);

//	flags = FFTW_MEASURE;
//	flags = FFTW_ESTIMATE;
	flags = FFTW_PATIENT;
//	flags = FFTW_EXHAUSTIVE;

	/* the flag of parameters; */
	if ( argc == 1 )
	{
		strcpy(paraDir, "test");	// default folder for parameters;
//		para_flag = "12fold0";		// default parameters;
		para_flag = "fcc40s91";		// default parameters;
	}
	else if ( argc == 2 )
	{
		strcpy(paraDir, "test");	// default folder for parameters;
		para_flag = argv[1];
	}
	else
	{
		strcpy(paraDir, argv[1]);
		para_flag = argv[2];
	}
//	bool paramGenerate = true;		// true for generating parameter files in batches;
	bool paramGenerate = false;		// flase for turning off the function;


	/* generate parameter files in batches; */
	if ( myrank == 0 && paramGenerate )
	{
		fn_param_batch ( );
		return 1;
	}


	/* obtain 'main_type'; */
	fn_obt_main_type ( );

    fftw_mpi_init();


	/* different functions; */
	double hamilton;
	if		( strcmp( main_type, "stable_system" ) == 0 ) // stable interface system;
	{
		hamilton = fn_main_stable_interface ( );
	}
	else if	( strcmp( main_type, "load_stable_system" ) == 0 ) // stable interface system by loading existing files;
	{
		hamilton = fn_main_load_stable_interface ( );
	}
	else if ( strcmp( main_type, "stable_bulk1" ) == 0 ||
			  strcmp( main_type, "stable_bulk2" ) == 0 ) // stable bulk phases;
	{
		hamilton = fn_main_stable_bulk ( );
	}
	else if ( strcmp( main_type, "stable_bulk1_parallel" ) == 0 ||
			  strcmp( main_type, "stable_bulk2_parallel" ) == 0 ) // stable bulk phases;
	{
		hamilton = fn_main_stable_bulk_parallel ( );
	}
	else if ( strcmp( main_type, "proj_bulk1" ) == 0 ||
			  strcmp( main_type, "proj_bulk2" ) == 0 )	// project bulk phases;
	{
		fn_main_project_bulk ( );
	}
	else if ( strcmp( main_type, "com_bulk1" ) == 0 ||
			  strcmp( main_type, "com_bulk2" ) == 0 )  // bulk phases in the common space;
	{
		fn_main_common_bulk	 ( );
	}
	else if ( strcmp( main_type, "com_projmat" ) == 0 ) // the common projection matrix;
	{
		fn_main_common_rotateProjBoxMat ( );
	}
	else if ( strcmp( main_type, "disp_density" ) == 0 ) // display densities according to files;
	{
		int status = fn_main_disp_system_density ( );
	}
	else if ( strcmp( main_type, "collect" ) == 0 )		// collect all data for interface system;
	{
		fn_main_collect_system_data ( );
	}
	else if ( strcmp( main_type, "recover" ) == 0 )		// recover density profile by the screened 'rhoJCplx'
	{
		fn_main_recover_system_density ( );
	}
	else
	{
		printf("Error 'main_type': %s.\n", main_type);
	}


	timer.pause();
	if ( myrank == 0  )
	{
		printf("\n\n\n\t\t time cost of program : %f seconds\n", timer.get_current_time());
		printf("\n\n ======================================   END PROGRAM  ======================================\n\n");
		char finishFile[FILELEN];
		sprintf(finishFile, "%s/finish.txt", rsltDir);
		FILE *fp = fopen(finishFile, "w");
		fprintf(fp, "done\n");
		fclose(fp);
	}

	fftw_mpi_cleanup(); // deallocate all persistent data and reset FFTW
    MPI_Finalize();

	return 0;
}
