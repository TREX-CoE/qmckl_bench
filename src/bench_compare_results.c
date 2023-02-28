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

#define ITERMAX 5
#define DEVICE_ID 0

int main(int argc, char** argv)
{

	qmckl_context_device context;
	context = qmckl_context_create_device(DEVICE_ID);
	qmckl_exit_code rc;

	if (argc < 2) {
		fprintf(stderr,"Syntax: %s FILE.    (FILEs are data/*.h5)\n", argv[0]);
		exit(-1);
	}
	char* file_name = argv[1];
	trexio_t* trexio_file = trexio_open(file_name, 'r', TREXIO_HDF5, &rc);

	int walk_num, elec_num;
	rc = trexio_read_electron_num(trexio_file, &elec_num);
	rc = trexio_read_qmc_num(trexio_file, &walk_num);


	double* elec_coord = malloc(sizeof(double)*walk_num*elec_num*3);
	double* elec_coord_device = omp_target_alloc(sizeof(double)*walk_num*elec_num*3, DEVICE_ID);

	rc = trexio_read_qmc_point(trexio_file, elec_coord);
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

	int64_t ao_num;
	rc = qmckl_get_ao_basis_ao_num(context, &ao_num);

	const int64_t size_max = 5*walk_num*elec_num*ao_num;
	double * ao_vgl_device = omp_target_alloc (size_max * sizeof(double), DEVICE_ID);

	rc = qmckl_set_electron_coord_device(context, 'N', walk_num, elec_coord_device, walk_num*elec_num*3);

	rc = qmckl_get_ao_basis_ao_vgl_device(context, ao_vgl_device, size_max);

	qmckl_context_touch_device(context);

	printf("End of GPU version\n");










	qmckl_context context_CPU;
	context_CPU = qmckl_context_create();

	if (argc < 2) {
		fprintf(stderr,"Syntax: %s FILE.    (FILEs are data/*.h5)\n", argv[0]);
		exit(-1);
	}
	char* file_name_CPU = argv[1];
	trexio_t* trexio_file_CPU = trexio_open(file_name_CPU, 'r', TREXIO_HDF5, &rc);

	int walk_num_CPU, elec_num_CPU;
	rc = trexio_read_electron_num(trexio_file_CPU, &elec_num_CPU);
	rc = trexio_read_qmc_num(trexio_file_CPU, &walk_num_CPU);
	double* elec_coord_CPU = malloc(sizeof(double)*walk_num_CPU*elec_num_CPU*3);
	rc = trexio_read_qmc_point(trexio_file_CPU, elec_coord_CPU);
	trexio_close(trexio_file_CPU);

	printf("Reading %s.\n", file_name_CPU);
	rc = qmckl_trexio_read(context_CPU, file_name_CPU, 255);
	if (rc != QMCKL_SUCCESS) {
		printf("%s\n", qmckl_string_of_error(rc));
	}

	int64_t ao_num_CPU;
	rc = qmckl_get_ao_basis_ao_num(context, &ao_num_CPU);

	const int64_t size_max_CPU = 5*walk_num_CPU*elec_num_CPU*ao_num_CPU;
	double * ao_vgl_CPU = malloc (size_max * sizeof(double));

	rc = qmckl_set_electron_coord(context_CPU, 'N', walk_num_CPU, elec_coord_CPU, walk_num_CPU*elec_num_CPU*3);

	rc = qmckl_get_ao_basis_ao_vgl(context_CPU, ao_vgl_CPU, size_max_CPU);


	printf("End CPU\n");
	//Compare two output
	
	double* ao_vgl_back;
	ao_vgl_back = malloc(size_max * sizeof(double));

	omp_target_memcpy(
			ao_vgl_back,
			ao_vgl_device,
			sizeof(double)*size_max,
			0, 0,
			omp_get_initial_device(), DEVICE_ID);

	for(int i =0; i < size_max; i ++)
	{
		if(fabs((ao_vgl_back[i]- ao_vgl_CPU[i])) >  0.001)  
		{
			printf("BAD VALUES: %lf\n",fabs((ao_vgl_back[i] - ao_vgl_CPU[i])));
			printf("device: %lf, CPU: %lf, i= %d\n", ao_vgl_back[i], ao_vgl_CPU[i], i);
			return -1;
		}
	}

	printf("ok\n");

	// FREE TIME
	rc = qmckl_context_destroy(context);
	free(ao_vgl_CPU);
	free(ao_vgl_back);
	free(elec_coord);
	omp_target_free(elec_coord_device, DEVICE_ID);
	omp_target_free(ao_vgl_device, DEVICE_ID);

}
