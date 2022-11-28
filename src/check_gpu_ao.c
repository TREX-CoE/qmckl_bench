#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <qmckl.h>
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


	qmckl_context context;
	context = qmckl_context_create();
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

	//GPU PART

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
	rc = qmckl_trexio_read_device(context, file_name, 255, DEVICE_ID);
	if (rc != QMCKL_SUCCESS) {
		printf("%s\n", qmckl_string_of_error(rc));
	}
	assert (rc == QMCKL_SUCCESS);

	int64_t ao_num;
	rc = qmckl_get_ao_basis_ao_num(context, &ao_num);
	assert (rc == QMCKL_SUCCESS);

	const int64_t size_max = 5*walk_num*elec_num*ao_num;
	double * ao_vgl_device = omp_target_alloc (size_max * sizeof(double), DEVICE_ID);
	double* ao_vgl_device_return = malloc(size_max * sizeof(double));
	assert (ao_vgl_device != NULL);


	rc = qmckl_set_electron_coord_device(context, 'N', walk_num, elec_coord_device, walk_num*elec_num*3, DEVICE_ID);
	assert (rc == QMCKL_SUCCESS);

	rc = qmckl_get_ao_basis_ao_vgl_device(context, ao_vgl_device, size_max, DEVICE_ID);
	assert (rc == QMCKL_SUCCESS);

	omp_target_memcpy(
			ao_vgl_device_return, ao_vgl_device,
			sizeof(double)*walk_num*elec_num*3,
			0, 0,
			omp_get_initial_device, DEVICE_ID);

	qmckl_context_touch(context);

	free(elec_coord);
	omp_target_free(elec_coord_device);




 	double* elec_coord = malloc(sizeof(double)*walk_num*elec_num*3);
	

	int64_t ao_num_CPU;
	rc = qmckl_get_ao_basis_ao_num(context, &ao_num);
	assert (rc == QMCKL_SUCCESS);

	const int64_t size_max = 5*walk_num*elec_num*ao_num;
	double * ao_vgl_CPU = malloc (size_max * sizeof(double));
	assert (ao_vgl_CPU != NULL);

	rc = qmckl_set_electron_coord(context, 'N', walk_num, elec_coord, walk_num*elec_num*3);
	assert (rc == QMCKL_SUCCESS);

	rc = qmckl_get_ao_basis_ao_vgl(context, ao_vgl_CPU, size_max);
	assert (rc == QMCKL_SUCCESS);


	//Compare result	
	//
	
	for(int i = 0; i < size_max; i++){
		if(ao_vgl[i] != ao_vgl_device_return[i]){
			printf("Fail test");
			return 0;
		}
	}

}
