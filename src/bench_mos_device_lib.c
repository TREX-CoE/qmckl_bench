#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <qmckl.h>
#include <qmckl_gpu.h>
#include <trexio.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <unistd.h>
#include <sys/time.h>

#include <omp.h>

#include <time.h>

#define ITERMAX 10
#define DEVICE_ID 0

int main(int argc, char** argv)
{
	long start, end;
	struct timeval timecheck;

	qmckl_context context;
	context = qmckl_context_create_device(DEVICE_ID);
	qmckl_exit_code rc;
	if (argc < 2) {
		fprintf(stderr,"Syntax: %s FILE.    (FILEs are data/*.h5)\n", argv[0]);
		exit(-1);
	}
	char* file_name = argv[1];
	trexio_t* trexio_file = trexio_open(file_name, 'r', TREXIO_HDF5, &rc);
	assert (rc == TREXIO_SUCCESS);

	int walk_num, elec_num;
	rc = trexio_read_electron_num(trexio_file, &elec_num);
	assert (rc == TREXIO_SUCCESS);
	rc = trexio_read_qmc_num(trexio_file, &walk_num);
	assert (rc == TREXIO_SUCCESS);


	double* elec_coord = malloc(sizeof(double)*walk_num*elec_num*3);
	double* elec_coord_device = omp_target_alloc(sizeof(double)*walk_num*elec_num*3, DEVICE_ID);

	assert (elec_coord != NULL);
	rc = trexio_read_qmc_point(trexio_file, elec_coord);
	assert (rc == TREXIO_SUCCESS);
	trexio_close(trexio_file);

	omp_target_memcpy(
			elec_coord_device, elec_coord,
			sizeof(double)*walk_num*elec_num*3,
			0, 0,
			DEVICE_ID, omp_get_initial_device()
			);


	printf("Reading %s.\n", file_name);
	rc = qmckl_trexio_read_device(context, file_name, 255);
	if (rc != QMCKL_SUCCESS) {
		printf("%s\n", qmckl_string_of_error(rc));
	}
	assert (rc == QMCKL_SUCCESS);

	int64_t mo_num;
	rc = qmckl_get_mo_basis_mo_num(context, &mo_num);
	assert (rc == QMCKL_SUCCESS);

	const int64_t size_max = 5*walk_num*elec_num*mo_num;
	double* mo_vgl = omp_target_alloc(size_max * sizeof(double), DEVICE_ID);
	assert (mo_vgl != NULL);

	rc = qmckl_set_electron_coord_device(context, 'N', walk_num, elec_coord_device, walk_num*elec_num*3);
	assert (rc == QMCKL_SUCCESS);

	rc = qmckl_get_mo_basis_mo_vgl_device(context, mo_vgl, size_max);
	assert (rc == QMCKL_SUCCESS);

	gettimeofday(&timecheck, NULL);
	start = (long)timecheck.tv_sec * 1000 + (long)timecheck.tv_usec / 1000;

	for (int i=0 ; i<ITERMAX ; ++i) {
		rc = qmckl_get_mo_basis_mo_vgl_inplace_device(context, mo_vgl, size_max);
	}
	gettimeofday(&timecheck, NULL);
	end = (long)timecheck.tv_sec * 1000 + (long)timecheck.tv_usec / 1000;

	printf("Time for the calculation of 1 step (ms): %10.1f\n", (double) (end-start) / (double) ITERMAX);
	rc = qmckl_context_destroy(context);
	omp_target_free(elec_coord, DEVICE_ID);
	omp_target_free(mo_vgl, DEVICE_ID);

}