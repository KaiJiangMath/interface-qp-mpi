#include "Data.h"

char main_type[STRLEN];		// the type of the main code;
char model_type[STRLEN];	// the type of the Landau model; LB/LP;
const char *para_flag;		// the serial number of parameters;
char paraDir[FILELEN0];		// the directory about parameters;
char rsltDir[FILELEN0];		// the directory for saving results;
unsigned flags;				// the flag of FFTW;

int nprocs, myrank;			// the number of process and the rank of this process;
int *displs_bulk,	*recvCount_bulk; // offset and receive count of bulk phases;
int *displs_sys,	*recvCount_sys;  // offset and receive count of interfaces; 
ptrdiff_t alloc_local,	   local_n0,	 local_0_start;		// local memory and index of bulk phases;
ptrdiff_t alloc_local_sys, local_n0_sys, local_0_start_sys; // local memory and index of interfaces;
