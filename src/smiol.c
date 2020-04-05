#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include "smiol.h"
#include "smiol_utils.h"

#ifdef SMIOL_PNETCDF
#include "pnetcdf.h"
#define PNETCDF_DEFINE_MODE 0
#define PNETCDF_DATA_MODE 1
#endif


/********************************************************************************
 *
 * SMIOL_fortran_init
 *
 * Initialize a SMIOL context from Fortran.
 *
 * This function is a simply a wrapper for the SMOIL_init routine that is intended
 * to be called from Fortran. Accordingly, the first argument is of type MPI_Fint
 * (a Fortran integer) rather than MPI_Comm.
 *
 ********************************************************************************/
int SMIOL_fortran_init(MPI_Fint comm, struct SMIOL_context **context)
{
	return SMIOL_init(MPI_Comm_f2c(comm), context);
}


/********************************************************************************
 *
 * SMIOL_init
 *
 * Initialize a SMIOL context.
 *
 * Initializes a SMIOL context, within which decompositions may be defined and
 * files may be read and written. At present, the only input argument is an MPI
 * communicator.
 *
 * Upon successful return the context argument points to a valid SMIOL context;
 * otherwise, it is NULL and an error code other than MPI_SUCCESS is returned.
 *
 * Note: It is assumed that MPI_Init has been called prior to this routine, so
 *       that any use of the provided MPI communicator will be valid.
 *
 ********************************************************************************/
int SMIOL_init(MPI_Comm comm, struct SMIOL_context **context)
{
	MPI_Comm smiol_comm;

	/*
	 * Before dereferencing context below, ensure that the pointer
	 * the context pointer is not NULL
	 */
	if (context == NULL) {
		return SMIOL_INVALID_ARGUMENT;
	}

	/*
	 * We cannot check for every possible invalid argument for comm, but
	 * at least we can verify that the communicator is not MPI_COMM_NULL
	 */
	if (comm == MPI_COMM_NULL) {
		/* Nullifying (*context) here may result in a memory leak, but this
		 * seems better than disobeying the stated behavior of returning
		 * a NULL context upon failure
		 */
		(*context) = NULL;

		return SMIOL_INVALID_ARGUMENT;
	}

	*context = (struct SMIOL_context *)malloc(sizeof(struct SMIOL_context));
	if ((*context) == NULL) {
		return SMIOL_MALLOC_FAILURE;
	}

	/*
	 * Initialize context
	 */
	(*context)->lib_ierr = 0;
	(*context)->lib_type = SMIOL_LIBRARY_UNKNOWN;

	/*
	 * Make a duplicate of the MPI communicator for use by SMIOL
	 */
	if (MPI_Comm_dup(comm, &smiol_comm) != MPI_SUCCESS) {
		free((*context));
		(*context) = NULL;
		return SMIOL_MPI_ERROR;
	}
	(*context)->fcomm = MPI_Comm_c2f(smiol_comm);

	if (MPI_Comm_size(smiol_comm, &((*context)->comm_size)) != MPI_SUCCESS) {
		free((*context));
		(*context) = NULL;
		return SMIOL_MPI_ERROR;
	}

	if (MPI_Comm_rank(smiol_comm, &((*context)->comm_rank)) != MPI_SUCCESS) {
		free((*context));
		(*context) = NULL;
		return SMIOL_MPI_ERROR;
	}

	return SMIOL_SUCCESS;
}


/********************************************************************************
 *
 * SMIOL_finalize
 *
 * Finalize a SMIOL context.
 *
 * Finalizes a SMIOL context and frees all memory in the SMIOL_context instance.
 * After this routine is called, no other SMIOL routines that make reference to
 * the finalized context should be called.
 *
 ********************************************************************************/
int SMIOL_finalize(struct SMIOL_context **context)
{
	MPI_Comm smiol_comm;

	/*
	 * If the pointer to the context pointer is NULL, assume we have nothing
	 * to do and declare success
	 */
	if (context == NULL) {
		return SMIOL_SUCCESS;
	}

	if ((*context) == NULL) {
		return SMIOL_SUCCESS;
	}

	smiol_comm = MPI_Comm_f2c((*context)->fcomm);
	if (MPI_Comm_free(&smiol_comm) != MPI_SUCCESS) {
		free((*context));
		(*context) = NULL;
		return SMIOL_MPI_ERROR;
	}

	free((*context));
	(*context) = NULL;

	return SMIOL_SUCCESS;
}


/********************************************************************************
 *
 * SMIOL_inquire
 *
 * Inquire about a SMIOL context.
 *
 * Detailed description.
 *
 ********************************************************************************/
int SMIOL_inquire(void)
{
	return SMIOL_SUCCESS;
}


/********************************************************************************
 *
 * SMIOL_open_file
 *
 * Opens a file within a SMIOL context.
 *
 * Depending on the specified file mode, creates or opens the file specified
 * by filename within the provided SMIOL context.
 *
 * Upon successful completion, SMIOL_SUCCESS is returned, and the file handle
 * argument will point to a valid file handle and the current frame for the
 * file will be set to zero. Otherwise, the file handle is NULL and an error
 * code other than SMIOL_SUCCESS is returned.
 *
 ********************************************************************************/
int SMIOL_open_file(struct SMIOL_context *context, const char *filename, int mode, struct SMIOL_file **file)
{
#ifdef SMIOL_PNETCDF
	int ierr;
#endif

	/*
	 * Before dereferencing file below, ensure that the pointer
	 * the file pointer is not NULL
	 */
	if (file == NULL) {
		return SMIOL_INVALID_ARGUMENT;
	}

	/*
	 * Check that context is valid
	 */
	if (context == NULL) {
		return SMIOL_INVALID_ARGUMENT;
	}

	*file = (struct SMIOL_file *)malloc(sizeof(struct SMIOL_file));
	if ((*file) == NULL) {
		return SMIOL_MALLOC_FAILURE;
	}

	/*
	 * Save pointer to context for this file
	 */
	(*file)->context = context;
	(*file)->frame = (SMIOL_Offset) 0;

	if (mode & SMIOL_FILE_CREATE) {
#ifdef SMIOL_PNETCDF
		if ((ierr = ncmpi_create(MPI_Comm_f2c(context->fcomm), filename,
					(NC_64BIT_DATA | NC_CLOBBER), MPI_INFO_NULL,
					&((*file)->ncidp))) != NC_NOERR) {
			free((*file));
			(*file) = NULL;
			context->lib_type = SMIOL_LIBRARY_PNETCDF;
			context->lib_ierr = ierr;
			return SMIOL_LIBRARY_ERROR;
		} else {
			(*file)->state = PNETCDF_DEFINE_MODE;
		}
#endif
	}
	else if (mode & SMIOL_FILE_WRITE) {
#ifdef SMIOL_PNETCDF
		if ((ierr = ncmpi_open(MPI_Comm_f2c(context->fcomm), filename,
				NC_WRITE, MPI_INFO_NULL, &((*file)->ncidp))) != NC_NOERR) {
			free((*file));
			(*file) = NULL;
			context->lib_type = SMIOL_LIBRARY_PNETCDF;
			context->lib_ierr = ierr;
			return SMIOL_LIBRARY_ERROR;
		} else {
			(*file)->state = PNETCDF_DATA_MODE;
		}
#endif
	}
	else if (mode & SMIOL_FILE_READ) {
#ifdef SMIOL_PNETCDF
		if ((ierr = ncmpi_open(MPI_Comm_f2c(context->fcomm), filename,
				NC_NOWRITE, MPI_INFO_NULL, &((*file)->ncidp))) != NC_NOERR) {
			free((*file));
			(*file) = NULL;
			context->lib_type = SMIOL_LIBRARY_PNETCDF;
			context->lib_ierr = ierr;
			return SMIOL_LIBRARY_ERROR;
		} else {
			(*file)->state = PNETCDF_DATA_MODE;
		}
#endif
	}
	else {
		free((*file));
		(*file) = NULL;
		return SMIOL_INVALID_ARGUMENT;
	}

	return SMIOL_SUCCESS;
}


/********************************************************************************
 *
 * SMIOL_close_file
 *
 * Closes a file within a SMIOL context.
 *
 * Closes the file associated with the provided file handle. Upon successful
 * completion, SMIOL_SUCCESS is returned, the file will be closed, and all memory
 * that is uniquely associated with the file handle will be deallocated.
 * Otherwise, an error code other than SMIOL_SUCCESS will be returned.
 *
 ********************************************************************************/
int SMIOL_close_file(struct SMIOL_file **file)
{
#ifdef SMIOL_PNETCDF
	int ierr;
#endif

	/*
	 * If the pointer to the file pointer is NULL, assume we have nothing
	 * to do and declare success
	 */
	if (file == NULL) {
		return SMIOL_SUCCESS;
	}

#ifdef SMIOL_PNETCDF
	if ((ierr = ncmpi_close((*file)->ncidp)) != NC_NOERR) {
		((*file)->context)->lib_type = SMIOL_LIBRARY_PNETCDF;
		((*file)->context)->lib_ierr = ierr;
		free((*file));
		(*file) = NULL;
		return SMIOL_LIBRARY_ERROR;
	}
#endif

	free((*file));
	(*file) = NULL;

	return SMIOL_SUCCESS;
}


/********************************************************************************
 *
 * SMIOL_define_dim
 *
 * Defines a new dimension in a file.
 *
 * Defines a dimension with the specified name and size in the file associated
 * with the file handle. If a negative value is provided for the size argument,
 * the dimension will be defined as an unlimited or record dimension.
 *
 * Upon successful completion, SMIOL_SUCCESS is returned; otherwise, an error
 * code is returned.
 *
 ********************************************************************************/
int SMIOL_define_dim(struct SMIOL_file *file, const char *dimname, SMIOL_Offset dimsize)
{
#ifdef SMIOL_PNETCDF
	int dimidp;
	int ierr;
	MPI_Offset len;
#endif

	/*
	 * Check that file handle is valid
	 */
	if (file == NULL) {
		return SMIOL_INVALID_ARGUMENT;
	}

	/*
	 * Check that dimension name is valid
	 */
	if (dimname == NULL) {
		return SMIOL_INVALID_ARGUMENT;
	}

#ifdef SMIOL_PNETCDF
	/*
	 * The parallel-netCDF library does not permit zero-length dimensions
	 */
	if (dimsize == (SMIOL_Offset)0) {
		return SMIOL_INVALID_ARGUMENT;
	}

	/*
	 * Handle unlimited / record dimension specifications
	 */
	if (dimsize < (SMIOL_Offset)0) {
		len = NC_UNLIMITED;
	}
	else {
		len = (MPI_Offset)dimsize;
	}

	/*
	 * If the file is in data mode, then switch it to define mode
	 */
	if (file->state == PNETCDF_DATA_MODE) {
		if ((ierr = ncmpi_redef(file->ncidp)) != NC_NOERR) {
			file->context->lib_type = SMIOL_LIBRARY_PNETCDF;
			file->context->lib_ierr = ierr;
			return SMIOL_LIBRARY_ERROR;
		}
		file->state = PNETCDF_DEFINE_MODE;
	}

	if ((ierr = ncmpi_def_dim(file->ncidp, dimname, len, &dimidp)) != NC_NOERR) {
		file->context->lib_type = SMIOL_LIBRARY_PNETCDF;
		file->context->lib_ierr = ierr;
		return SMIOL_LIBRARY_ERROR;
	}
#endif

	return SMIOL_SUCCESS;
}


/********************************************************************************
 *
 * SMIOL_inquire_dim
 *
 * Inquires about an existing dimension in a file.
 *
 * Inquire about the size of an existing dimension and optionally inquire if the
 * given dimension is the unlimited dimension or not. If dimsize is a non-NULL
 * pointer then the dimension size will be returned in dimsize. For unlimited
 * dimensions, the current size of the dimension is returned; future writes of
 * additional records to a file can lead to different return sizes for
 * unlimited dimensions.
 *
 * If is_unlimited is a non-NULL pointer and if the inquired dimension is the
 * unlimited dimension, is_unlimited will be set to 1; if the inquired
 * dimension is not the unlimited dimension then is_unlimited will be set to 0.
 *
 * Upon successful completion, SMIOL_SUCCESS is returned; otherwise, an error
 * code is returned.
 *
 ********************************************************************************/
int SMIOL_inquire_dim(struct SMIOL_file *file, const char *dimname,
                      SMIOL_Offset *dimsize, int *is_unlimited)
{
#ifdef SMIOL_PNETCDF
	int dimidp;
	int ierr;
	MPI_Offset len;
#endif
	/*
	 * Check that file handle is valid
	 */
	if (file == NULL) {
		return SMIOL_INVALID_ARGUMENT;
	}

	/*
	 * Check that dimension name is valid
	 */
	if (dimname == NULL) {
		return SMIOL_INVALID_ARGUMENT;
	}

	/*
	 * Check that dimension size is not NULL
	 */
	if (dimsize == NULL && is_unlimited == NULL) {
		return SMIOL_INVALID_ARGUMENT;
	}

	if (dimsize != NULL) {
		(*dimsize) = (SMIOL_Offset)0;   /* Default dimension size if no library provides a value */
	}

	if (is_unlimited != NULL) {
		(*is_unlimited) = 0; /* Return 0 if no library provides a value */
	}

#ifdef SMIOL_PNETCDF
	if ((ierr = ncmpi_inq_dimid(file->ncidp, dimname, &dimidp)) != NC_NOERR) {
		(*dimsize) = (SMIOL_Offset)(-1);  /* TODO: should there be a well-defined invalid size? */
		file->context->lib_type = SMIOL_LIBRARY_PNETCDF;
		file->context->lib_ierr = ierr;
		return SMIOL_LIBRARY_ERROR;
	}

	/*
	 * Inquire about dimsize
	 */
	if (dimsize != NULL) {
		ierr = ncmpi_inq_dimlen(file->ncidp, dimidp, &len);
		if (ierr != NC_NOERR) {
			(*dimsize) = (SMIOL_Offset)(-1);  /* TODO: should there be a well-defined invalid size? */
			file->context->lib_type = SMIOL_LIBRARY_PNETCDF;
			file->context->lib_ierr = ierr;
			return SMIOL_LIBRARY_ERROR;
		}

		(*dimsize) = (SMIOL_Offset)len;
	}


	/*
	 * Inquire if this dimension is the unlimited dimension
	 */
	if (is_unlimited != NULL) {
		int unlimdimidp;
		ierr = ncmpi_inq_unlimdim(file->ncidp, &unlimdimidp);
		if (ierr != NC_NOERR) {
			file->context->lib_type = SMIOL_LIBRARY_PNETCDF;
			file->context->lib_ierr = ierr;
			return SMIOL_LIBRARY_ERROR;
		}
		if (unlimdimidp == dimidp) {
			(*is_unlimited) = 1;
		} else {
			(*is_unlimited) = 0; // Not the unlimited dim
		}
	}
#endif

	return SMIOL_SUCCESS;
}


/********************************************************************************
 *
 * SMIOL_define_var
 *
 * Defines a new variable in a file.
 *
 * Defines a variable with the specified name, type, and dimensions in an open
 * file pointed to by the file argument. The varname and dimnames arguments
 * are expected to be null-terminated strings, except if the variable has zero
 * dimensions, in which case the dimnames argument may be a NULL pointer.
 *
 * Upon successful completion, SMIOL_SUCCESS is returned; otherwise, an error
 * code is returned.
 *
 ********************************************************************************/
int SMIOL_define_var(struct SMIOL_file *file, const char *varname, int vartype, int ndims, const char **dimnames)
{
#ifdef SMIOL_PNETCDF
	int *dimids;
	int ierr;
	int i;
	nc_type xtype;
	int varidp;
#endif

	/*
	 * Check that file handle is valid
	 */
	if (file == NULL) {
		return SMIOL_INVALID_ARGUMENT;
	}

	/*
	 * Check that variable name is valid
	 */
	if (varname == NULL) {
		return SMIOL_INVALID_ARGUMENT;
	}

	/*
	 * Check that the variable type is valid - handled below in a library-specific way...
	 */

	/*
	 * Check that variable dimension names are valid
	 */
	if (dimnames == NULL && ndims > 0) {
		return SMIOL_INVALID_ARGUMENT;
	}

#ifdef SMIOL_PNETCDF
	dimids = (int *)malloc(sizeof(int) * (size_t)ndims);
	if (dimids == NULL) {
		return SMIOL_MALLOC_FAILURE;
	}

	/*
	 * Build a list of dimension IDs
	 */
	for (i=0; i<ndims; i++) {
		if ((ierr = ncmpi_inq_dimid(file->ncidp, dimnames[i], &dimids[i])) != NC_NOERR) {
			free(dimids);
			file->context->lib_type = SMIOL_LIBRARY_PNETCDF;
			file->context->lib_ierr = ierr;
			return SMIOL_LIBRARY_ERROR;
		}
	}

	/*
	 * Translate SMIOL variable type to parallel-netcdf type
	 */
	switch (vartype) {
		case SMIOL_REAL32:
			xtype = NC_FLOAT;
			break;
		case SMIOL_REAL64:
			xtype = NC_DOUBLE;
			break;
		case SMIOL_INT32:
			xtype = NC_INT;
			break;
		case SMIOL_CHAR:
			xtype = NC_CHAR;
			break;
		default:
			free(dimids);
			return SMIOL_INVALID_ARGUMENT;
	}

	/*
	 * If the file is in data mode, then switch it to define mode
	 */
	if (file->state == PNETCDF_DATA_MODE) {
		if ((ierr = ncmpi_redef(file->ncidp)) != NC_NOERR) {
			file->context->lib_type = SMIOL_LIBRARY_PNETCDF;
			file->context->lib_ierr = ierr;
			return SMIOL_LIBRARY_ERROR;
		}
		file->state = PNETCDF_DEFINE_MODE;
	}

	/*
	 * Define the variable
	 */
	if ((ierr = ncmpi_def_var(file->ncidp, varname, xtype, ndims, dimids, &varidp)) != NC_NOERR) {
		free(dimids);
		file->context->lib_type = SMIOL_LIBRARY_PNETCDF;
		file->context->lib_ierr = ierr;
		return SMIOL_LIBRARY_ERROR;
	}

	free(dimids);
#endif

	return SMIOL_SUCCESS;
}


/********************************************************************************
 *
 * SMIOL_inquire_var
 *
 * Inquires about an existing variable in a file.
 *
 * Inquires about a variable in a file, and optionally returns the type
 * of the variable, the dimensionality of the variable, and the names of
 * the dimensions of the variable. Which properties of the variable to return
 * (type, dimensionality, or dimension names) is indicated by the status of
 * the pointers for the corresponding properties: if the pointer is a non-NULL
 * pointer, the property will be set upon successful completion of this routine.
 *
 * If the names of a variable's dimensions are requested (by providing a non-NULL
 * actual argument for dimnames), the size of the dimnames array must be at least
 * the number of dimensions in the variable, and each character string pointed
 * to by an element of dimnames must be large enough to accommodate the corresponding
 * dimension name.
 *
 ********************************************************************************/
int SMIOL_inquire_var(struct SMIOL_file *file, const char *varname, int *vartype, int *ndims, char **dimnames)
{
#ifdef SMIOL_PNETCDF
	int *dimids;
	int varidp;
	int ierr;
	int i;
	int xtypep;
	int ndimsp;
#endif

	/*
	 * Check that file handle is valid
	 */
	if (file == NULL) {
		return SMIOL_INVALID_ARGUMENT;
	}

	/*
	 * Check that variable name is valid
	 */
	if (varname == NULL) {
		return SMIOL_INVALID_ARGUMENT;
	}

	/*
	 * If all output arguments are NULL, we can return early
	 */
	if (vartype == NULL && ndims == NULL && dimnames == NULL) {
		return SMIOL_SUCCESS;
	}

#ifdef SMIOL_PNETCDF
	/*
	 * Get variable ID
	 */
	if ((ierr = ncmpi_inq_varid(file->ncidp, varname, &varidp)) != NC_NOERR) {
		file->context->lib_type = SMIOL_LIBRARY_PNETCDF;
		file->context->lib_ierr = ierr;
		return SMIOL_LIBRARY_ERROR;
	}

	/*
	 * If requested, inquire about variable type
	 */
	if (vartype != NULL) {
		if ((ierr = ncmpi_inq_vartype(file->ncidp, varidp, &xtypep)) != NC_NOERR) {
			file->context->lib_type = SMIOL_LIBRARY_PNETCDF;
			file->context->lib_ierr = ierr;
			return SMIOL_LIBRARY_ERROR;
		}

		/* Convert parallel-netCDF variable type to SMIOL variable type */
		switch (xtypep) {
			case NC_FLOAT:
				*vartype = SMIOL_REAL32;
				break;
			case NC_DOUBLE:
				*vartype = SMIOL_REAL64;
				break;
			case NC_INT:
				*vartype = SMIOL_INT32;
				break;
			case NC_CHAR:
				*vartype = SMIOL_CHAR;
				break;
			default:
				*vartype = SMIOL_UNKNOWN_VAR_TYPE;
		}
	}

	/*
	 * All remaining properties will require the number of dimensions
	 */
	if (ndims != NULL || dimnames != NULL) {
		if ((ierr = ncmpi_inq_varndims(file->ncidp, varidp, &ndimsp)) != NC_NOERR) {
			file->context->lib_type = SMIOL_LIBRARY_PNETCDF;
			file->context->lib_ierr = ierr;
			return SMIOL_LIBRARY_ERROR;
		}
	}

	/*
	 * If requested, inquire about dimensionality
	 */
	if (ndims != NULL) {
		*ndims = ndimsp;
	}

	/*
	 * If requested, inquire about dimension names
	 */
	if (dimnames != NULL) {
		dimids = (int *)malloc(sizeof(int) * (size_t)ndimsp);
		if (dimids == NULL) {
			return SMIOL_MALLOC_FAILURE;
		}

		if ((ierr = ncmpi_inq_vardimid(file->ncidp, varidp, dimids)) != NC_NOERR) {
			file->context->lib_type = SMIOL_LIBRARY_PNETCDF;
			file->context->lib_ierr = ierr;
			free(dimids);
			return SMIOL_LIBRARY_ERROR;
		}

		for (i = 0; i < ndimsp; i++) {
			if (dimnames[i] == NULL) {
				return SMIOL_INVALID_ARGUMENT;
			}
			if ((ierr = ncmpi_inq_dimname(file->ncidp, dimids[i], dimnames[i])) != NC_NOERR) {
				file->context->lib_type = SMIOL_LIBRARY_PNETCDF;
				file->context->lib_ierr = ierr;
				free(dimids);
				return SMIOL_LIBRARY_ERROR;
			}
		}

		free(dimids);
	}
#endif

	return SMIOL_SUCCESS;
}


/********************************************************************************
 *
 * SMIOL_put_var
 *
 * Writes a variable to a file.
 *
 * Write values into a variable for an open SMIOL file. If varname is
 * dimensioned by a decomposed dimension, then a SMIOL_decomp struct should
 * be present that was created for that dimension with SMIOL_create_decomp.
 * 
 * If a variable is dimensioned by a decomposed dimension, then decomp should
 * be NULL. In such case, all values of buf should be the same across tasks, as
 * this routine does not guaranteed any task will write out its values over
 * another.
 *
 * Currently, this function cannot write out fields that are greater than 2GB
 * in size.
 *
 * This routine will need to be called by all tasks of an MPI communicator.
 *
 ********************************************************************************/
int SMIOL_put_var(struct SMIOL_file *file, struct SMIOL_decomp *decomp,
                  const char *varname, const void *buf)
{
#ifdef SMIOL_PNETCDF
	int i;
	int ierr;
	int varidp;
	int vartype;
	int ndims;
	int innerdim_size;
	int is_unlimited = 0;
	char **dimnames;
	size_t *start;
	size_t *count;
	size_t dsize;
	SMIOL_Offset dimsize;
	void *io_buf;
#endif
	if (file == NULL) {
		return SMIOL_INVALID_ARGUMENT;
	}

	if (varname == NULL) {
		return SMIOL_INVALID_ARGUMENT;
	}

	if (buf == NULL) {
		return SMIOL_INVALID_ARGUMENT;
	}
#ifdef SMIOL_PNETCDF
	/*
	 * Retrieve vartype and number of dims for this var
	 */
	ierr = SMIOL_inquire_var(file, varname, &vartype, &ndims, NULL);
	if (ierr != SMIOL_SUCCESS) {
		return ierr;
	}

	dimnames = (char **) malloc(sizeof(char *) * (size_t) ndims);
	for (i = 0; i < ndims; i++) {
		dimnames[i] = (char *)malloc(sizeof(char) * (size_t) 64);
	}
	
	ierr = SMIOL_inquire_var(file, varname, NULL, NULL, dimnames);
	if (ierr != SMIOL_SUCCESS) {
		return ierr;
	}

	start = malloc(sizeof(size_t) * (size_t) ndims);
	count = malloc(sizeof(size_t) * (size_t) ndims);

	innerdim_size = 1;
	for (i = 0; i < ndims; i++) {
		if (i == 0 && decomp != NULL) {
			start[i] = decomp->io_start;
			count[i] = decomp->io_count;
			continue;
		} else {
			ierr = SMIOL_inquire_dim(file, dimnames[i], &dimsize,
			                         &is_unlimited);
			if (ierr != SMIOL_SUCCESS) {
				return ierr;
			}

			/*
			 * If this var is dimensioned by a unlimited dimension,
			 * then write on the last set frame.
			 */
			if (is_unlimited) {
				start[i] = file->frame;
				count[i] = 1;
			} else {
				start[i] = 0;
				count[i] = dimsize;
			}
			innerdim_size *= dimsize;
		}
	}

	for (i = 0; i < ndims; i++) {
		free(dimnames[i]);
	}
	free(dimnames);

	switch (vartype) {
		case SMIOL_REAL32:
			dsize = sizeof(float);
			break;
		case SMIOL_REAL64:
			dsize = sizeof(double);
			break;
		case SMIOL_INT32:
			dsize = sizeof(int);
			break;
		case SMIOL_CHAR:
			dsize = sizeof(char);
			break;
	}

	/*
	 * Transfer fields
	 */
	if (decomp != NULL) {
		io_buf = malloc(dsize * count[0] * innerdim_size);
		ierr = transfer_field(decomp, SMIOL_COMP_TO_IO,
		                      dsize * innerdim_size,
				      buf,
				      io_buf);
		if (ierr != SMIOL_SUCCESS) {
			return ierr;
		}
	}

	/*
	 * Get variable ID
	 */
	ierr = ncmpi_inq_varid(file->ncidp, varname, &varidp);
	if (ierr != NC_NOERR) {
		file->context->lib_type = SMIOL_LIBRARY_PNETCDF;
		file->context->lib_ierr = ierr;
		return SMIOL_LIBRARY_ERROR;
	}

	/*
	 * If the file is in define mode, then switch it to data mode
	 */
	if (file->state == PNETCDF_DEFINE_MODE) {
		if ((ierr = ncmpi_enddef(file->ncidp)) != NC_NOERR) {
			file->context->lib_type = SMIOL_LIBRARY_PNETCDF;
			file->context->lib_ierr = ierr;
			return SMIOL_LIBRARY_ERROR;
		}
		file->state = PNETCDF_DATA_MODE;
	}

	/*
	 * Write the variable
	 */
	if (decomp) {
		/* If this is a decomposed variable, then use the result from
		 * transfer_field, io_buf.
		 */
		ierr = ncmpi_put_vara_all(file->ncidp, varidp,
					 (MPI_Offset *) start,
					 (MPI_Offset *) count,
					 io_buf,
					 (MPI_Offset) NULL, MPI_DATATYPE_NULL);
		free(io_buf);
	} else {
		/* If this is not a decomposed variable, then use buf to write
		 * out the this variable
		 */
		ierr = ncmpi_put_vara_all(file->ncidp, varidp,
					 (MPI_Offset *) start,
					 (MPI_Offset *) count,
					 buf,
					 (MPI_Offset) NULL, MPI_DATATYPE_NULL);
	}
	if (ierr != NC_NOERR) {
		file->context->lib_type = SMIOL_LIBRARY_PNETCDF;
		file->context->lib_ierr = ierr;
		return SMIOL_LIBRARY_ERROR;
	}

	free(start);
	free(count);
#endif

	return SMIOL_SUCCESS;
}


/********************************************************************************
 *
 * SMIOL_get_var
 *
 * Reads a variable from a file.
 *
 * Read a variable from an open SMIOL_file into buf.  If the variable is
 * dimensioned by a decomposed dimension then a decomp that was created with
 * SMIOL_create_decomp for the decomposed dimension should be present.
 *
 * If this variable is not defined by a decomposed dimension, then decomp
 * should be NULL.
 *
 * For decomposed variables, buf must be large enough to hold the number of
 * elements for a given task, based on the decomposition. For non-decomposed
 * variables, buf must be large enough to hold the variable for all tasks.
 *
 * Currently, this function cannot read in fields that are greater than 2GB
 * in size.
 *
 ********************************************************************************/
int SMIOL_get_var(struct SMIOL_file *file, struct SMIOL_decomp *decomp,
                  const char *varname, void *buf)
{
#ifdef SMIOL_PNETCDF
	int i;
	int ierr;
	int varidp;
	int vartype;
	int ndims;
	int is_unlimited;
	int innerdim_size;
	char **dimnames;
	size_t *start;
	size_t *count;
	size_t dsize;
	SMIOL_Offset dimsize;
	void *io_buf;
#endif
	if (file == NULL) {
		return SMIOL_INVALID_ARGUMENT;
	}

	if (varname == NULL) {
		return SMIOL_INVALID_ARGUMENT;
	}

	if (buf == NULL) {
		return SMIOL_INVALID_ARGUMENT;
	}
#ifdef SMIOL_PNETCDF
	/*
	 * Retrieve vartype and number of dims for this var
	 */
	ierr = SMIOL_inquire_var(file, varname, &vartype, &ndims, NULL);
	if (ierr != SMIOL_SUCCESS) {
		return ierr;
	}

	dimnames = (char **) malloc(sizeof(char *) * (size_t) ndims);
	for (i = 0; i < ndims; i++) {
		dimnames[i] = (char *)malloc(sizeof(char) * (size_t) 64);
	}

	ierr = SMIOL_inquire_var(file, varname, NULL, NULL, dimnames);
	if (ierr != SMIOL_SUCCESS) {
		return ierr;
	}

	start = malloc(sizeof(size_t) * (size_t) ndims);
	count = malloc(sizeof(size_t) * (size_t) ndims);

	/*
	 * Get the start and count values for the rest of the dimensions. If
	 * this is a decomposed field, then the first dimension's start and
	 * count values will come from the decomp.
	 */
	innerdim_size = 1;
	for (i = 0; i < ndims; i++) {
		if (i == 0 && decomp != NULL) {
			start[i] = decomp->io_start;
			count[i] = decomp->io_count;
			continue;
		} else {
			ierr = SMIOL_inquire_dim(file, dimnames[i], &dimsize,
			                         &is_unlimited);
			if (ierr != SMIOL_SUCCESS) {
				return ierr;
			}

			/*
			 * If this var is dimensioned by a unlimited dimension,
			 * then get the field from the last set frame.
			 */
			if (is_unlimited) {
				start[i] = file->frame;
				count[i] = 1;
			} else {
				start[i] = 0;
				count[i] = dimsize;
			}

			innerdim_size *= dimsize;
		}
	}

	for (i = 0; i < ndims; i++) {
		free(dimnames[i]);
	}
	free(dimnames);

	switch (vartype) {
		case SMIOL_REAL32:
			dsize = sizeof(float);
			break;
		case SMIOL_REAL64:
			dsize = sizeof(double);
			break;
		case SMIOL_INT32:
			dsize = sizeof(int);
			break;
		case SMIOL_CHAR:
			dsize = sizeof(char);
			break;
	}

	/*
	 * Get variable ID
	 */
	ierr = ncmpi_inq_varid(file->ncidp, varname, &varidp);
	if (ierr != NC_NOERR) {
		file->context->lib_type = SMIOL_LIBRARY_PNETCDF;
		file->context->lib_ierr = ierr;
		return SMIOL_LIBRARY_ERROR;
	}

	/*
	 * If the file is in define mode, then switch it to data mode
	 */
	if (file->state == PNETCDF_DEFINE_MODE) {
		if ((ierr = ncmpi_enddef(file->ncidp)) != NC_NOERR) {
			file->context->lib_type = SMIOL_LIBRARY_PNETCDF;
			file->context->lib_ierr = ierr;
			return SMIOL_LIBRARY_PNETCDF;
		}
		file->state = PNETCDF_DATA_MODE;
	}

	/*
	 * Read the variable
	 */
	if (decomp) {
		io_buf = malloc(dsize * count[0] * innerdim_size);
		ierr = ncmpi_get_vara_all(file->ncidp, varidp,
					 (MPI_Offset *) start,
					 (MPI_Offset *) count,
					 io_buf,
					 (MPI_Offset) NULL, MPI_DATATYPE_NULL);
	} else {
		ierr = ncmpi_get_vara_all(file->ncidp, varidp,
					 (MPI_Offset *) start,
					 (MPI_Offset *) count,
					 buf,
					 (MPI_Offset) NULL, MPI_DATATYPE_NULL);
	}

	if (ierr != NC_NOERR) {
		file->context->lib_type = SMIOL_LIBRARY_PNETCDF;
		file->context->lib_ierr = ierr;
		return SMIOL_LIBRARY_ERROR;
	}

	free(start);
	free(count);

	/*
	 * Transfer fields
	 */
	if (decomp) {
		ierr = transfer_field(decomp, SMIOL_IO_TO_COMP,
		                      dsize * innerdim_size,
		                      io_buf, buf);
		free(io_buf);
	}
#endif

	return SMIOL_SUCCESS;
}


/********************************************************************************
 *
 * SMIOL_define_att
 *
 * Defines a new attribute in a file.
 *
 * Detailed description.
 *
 ********************************************************************************/
int SMIOL_define_att(void)
{
	return SMIOL_SUCCESS;
}


/********************************************************************************
 *
 * SMIOL_inquire_att
 *
 * Inquires about an attribute in a file.
 *
 * Detailed description.
 *
 ********************************************************************************/
int SMIOL_inquire_att(void)
{
	return SMIOL_SUCCESS;
}


/********************************************************************************
 *
 * SMIOL_sync_file
 *
 * Forces all in-memory data to be flushed to disk.
 *
 * Upon success, all in-memory data for the file associatd with the file
 * handle will be flushed to the file system and SMIOL_SUCCESS will be
 * returned; otherwise, an error code is returned.
 *
 ********************************************************************************/
int SMIOL_sync_file(struct SMIOL_file *file)
{
#ifdef SMIOL_PNETCDF
	int ierr;
#endif

	/*
	 * Check that file is valid
	 */
	if (file == NULL) {
		return SMIOL_INVALID_ARGUMENT;
	}

#ifdef SMIOL_PNETCDF
	/*
	 * If the file is in define mode then switch it into data mode
	 */
	if (file->state == PNETCDF_DEFINE_MODE) {
		if ((ierr = ncmpi_enddef(file->ncidp)) != NC_NOERR) {
			file->context->lib_type = SMIOL_LIBRARY_PNETCDF;
			file->context->lib_ierr = ierr;
			return SMIOL_LIBRARY_ERROR;
		}
		file->state = PNETCDF_DATA_MODE;
	}

	if ((ierr = ncmpi_sync(file->ncidp)) != NC_NOERR) {
		file->context->lib_type = SMIOL_LIBRARY_PNETCDF;
		file->context->lib_ierr = ierr;
		return SMIOL_LIBRARY_ERROR;
	}
#endif

	return SMIOL_SUCCESS;
}


/********************************************************************************
 *
 * SMIOL_error_string
 *
 * Returns an error string for a specified error code.
 *
 * Returns an error string corresponding to a SMIOL error code. If the error code is
 * SMIOL_LIBRARY_ERROR and a valid SMIOL context is available, the SMIOL_lib_error_string
 * function should be called instead. The error string is null-terminated, but it
 * does not contain a newline character.
 *
 ********************************************************************************/
const char *SMIOL_error_string(int errno)
{
	switch (errno) {
	case SMIOL_SUCCESS:
		return "Success!";
	case SMIOL_MALLOC_FAILURE:
		return "malloc returned a null pointer";
	case SMIOL_INVALID_ARGUMENT:
		return "invalid subroutine argument";
	case SMIOL_MPI_ERROR:
		return "internal MPI call failed";
	case SMIOL_FORTRAN_ERROR:
		return "Fortran wrapper detected an inconsistency in C return values";
	case SMIOL_LIBRARY_ERROR:
		return "bad return code from a library call";
	default:
		return "Unknown error";
	}
}


/********************************************************************************
 *
 * SMIOL_lib_error_string
 *
 * Returns an error string for a third-party library called by SMIOL.
 *
 * Returns an error string corresponding to an error that was generated by
 * a third-party library that was called by SMIOL. The library that was the source
 * of the error, as well as the library-specific error code, are retrieved from
 * a SMIOL context. If successive library calls resulted in errors, only the error
 * string for the last of these errors will be returned. The error string is
 * null-terminated, but it does not contain a newline character.
 *
 ********************************************************************************/
const char *SMIOL_lib_error_string(struct SMIOL_context *context)
{
	if (context == NULL) {
		return "SMIOL_context argument is a NULL pointer";
	}

	switch (context->lib_type) {
#ifdef SMIOL_PNETCDF
	case SMIOL_LIBRARY_PNETCDF:
		return ncmpi_strerror(context->lib_ierr);
#endif
	default:
		return "Could not find matching library for the source of the error";
	}
}


/********************************************************************************
 *
 * SMIOL_set_option
 *
 * Sets an option for the SMIOL library.
 *
 * Detailed description.
 *
 ********************************************************************************/
int SMIOL_set_option(void)
{
	return SMIOL_SUCCESS;
}

/********************************************************************************
 *
 * SMIOL_set_frame
 *
 * Set the frame for the unlimited dimension for an open file
 *
 * For an open SMIOL file handle, set the frame for the unlimited dimension.
 * After setting the frame for a file, writing to a variable that is
 * dimensioned by the unlimited dimension will write to the last set frame,
 * overwriting any current data that maybe present in that frame.
 *
 * SMIOL_SUCCESS will be returned if the frame is successfully set otherwise an
 * error will return.
 *
 ********************************************************************************/
int SMIOL_set_frame(struct SMIOL_file *file, SMIOL_Offset frame)
{
	if (file == NULL) {
		return SMIOL_INVALID_ARGUMENT;
	}
	file->frame = frame;
	return SMIOL_SUCCESS;
}

/********************************************************************************
 *
 * SMIOL_get_frame
 *
 * Return the current frame of an open file
 *
 * Get the current frame of an open file. Upon success, SMIOL_SUCCESS will be
 * returned, otherwise an error will be returned.
 *
 ********************************************************************************/
int SMIOL_get_frame(struct SMIOL_file *file, SMIOL_Offset *frame)
{
	if (file == NULL || frame == NULL) {
		return SMIOL_INVALID_ARGUMENT;
	}
	*frame = file->frame;
	return SMIOL_SUCCESS;
}


/*******************************************************************************
 *
 * SMIOL_create_decomp
 *
 * Creates a mapping between compute elements and I/O elements.
 *
 * Given arrays of global element IDs that each task computes, the number of I/O
 * tasks, and the stride between I/O tasks, this routine works out a mapping of
 * elements between compute and I/O tasks.
 *
 * If all input arguments are determined to be valid and if the routine is
 * successful in working out a mapping, the decomp pointer is allocated and
 * given valid contents, and SMIOL_SUCCESS is returned; otherwise a non-success
 * error code is returned and the decomp pointer is NULL.
 *
 *******************************************************************************/
int SMIOL_create_decomp(struct SMIOL_context *context,
                        size_t n_compute_elements, SMIOL_Offset *compute_elements,
                        int num_io_tasks, int io_stride,
                        struct SMIOL_decomp **decomp)
{
	size_t i;
	size_t n_io_elements, n_io_elements_global;
	size_t io_start, io_count;
	SMIOL_Offset *io_elements;
	MPI_Comm comm;
	MPI_Datatype dtype;
	int ierr;


	/*
	 * Minimal check on the validity of arguments
	 */
	if (context == NULL) {
		return SMIOL_INVALID_ARGUMENT;
	}

	if (compute_elements == NULL && n_compute_elements != 0) {
		return SMIOL_INVALID_ARGUMENT;
	}

	comm = MPI_Comm_f2c(context->fcomm);

	/*
	 * Figure out MPI_Datatype for size_t... there must be a better way...
	 */
	switch (sizeof(size_t)) {
		case sizeof(uint64_t):
			dtype = MPI_UINT64_T;
			break;
		case sizeof(uint32_t):
			dtype = MPI_UINT32_T;
			break;
		case sizeof(uint16_t):
			dtype = MPI_UINT16_T;
			break;
		default:
			return SMIOL_MPI_ERROR;
	}

	/*
	 * Based on the number of compute elements for each task, determine
	 * the total number of elements across all tasks for I/O. The assumption
	 * is that the number of elements to read/write is equal to the size of
	 * the set of compute elements.
	 */
	n_io_elements = n_compute_elements;
	if (MPI_SUCCESS != MPI_Allreduce((const void *)&n_io_elements,
	                                 (void *)&n_io_elements_global,
	                                 1, dtype, MPI_SUM, comm)) {
		return SMIOL_MPI_ERROR;
	}

	/*
	 * Determine the contiguous range of elements to be read/written by
	 * this MPI task
	 */
	ierr = get_io_elements(context->comm_rank, num_io_tasks, io_stride,
	                       n_io_elements_global, &io_start, &io_count);

	/*
	 * Fill in io_elements from io_start through io_start + io_count - 1
	 */
	io_elements = NULL;
	if (io_count > 0) {
		io_elements = (SMIOL_Offset *)malloc(sizeof(SMIOL_Offset)
		                                     * n_io_elements_global);
		if (io_elements == NULL) {
			return SMIOL_MALLOC_FAILURE;
		}
		for (i = 0; i < io_count; i++) {
			io_elements[i] = (SMIOL_Offset)(io_start + i);
		}
	}

	/*
	 * Build the mapping between compute tasks and I/O tasks
	 */
	ierr = build_exchange(context,
	                      n_compute_elements, compute_elements,
	                      io_count, io_elements,
	                      decomp);

	free(io_elements);

	/*
	 * If decomp was successfully created, add io_start and io_count values
	 * to the decomp before returning
	 */
	if (ierr == SMIOL_SUCCESS) {
		(*decomp)->io_start = io_start;
		(*decomp)->io_count = io_count;
	}

	return ierr;
}


/********************************************************************************
 *
 * SMIOL_free_decomp
 *
 * Frees a mapping between compute elements and I/O elements.
 *
 * Free all memory of a SMIOL_decomp and returns SMIOL_SUCCESS. If decomp
 * points to NULL, then do nothing and return SMIOL_SUCCESS. After this routine
 * is called, no other SMIOL routines should use the freed SMIOL_decomp.
 *
 ********************************************************************************/
int SMIOL_free_decomp(struct SMIOL_decomp **decomp)
{
	if ((*decomp) == NULL) {
		return SMIOL_SUCCESS;
	}

	free((*decomp)->comp_list);
	free((*decomp)->io_list);
	free((*decomp));
	*decomp = NULL;

	return SMIOL_SUCCESS;
}
