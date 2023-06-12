#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <qmckl_gpu.h>
#include <trexio.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <unistd.h>
#include <sys/time.h>

#include <time.h>

#include "../data/jast_coef.h"

#define ITERMAX 3
#define DEVICE_ID 0

int main(int argc, char **argv) {
	long start, end;
	struct timeval timecheck;

	qmckl_context_device context;
	context = qmckl_context_create_device(DEVICE_ID);
	qmckl_exit_code_device rc;

	if (argc < 2) {
		fprintf(stderr, "Syntax: %s FILE.    (FILEs are data/*.h5)\n", argv[0]);
		exit(-1);
	}
	char *file_name = argv[1];
	trexio_t *trexio_file = trexio_open(file_name, 'r', TREXIO_HDF5, &rc);
	assert(rc == TREXIO_SUCCESS);

	int walk_num;
	rc = trexio_read_qmc_num(trexio_file, &walk_num);
	assert(rc == TREXIO_SUCCESS);

	int elec_num;
	rc = trexio_read_electron_num(trexio_file, &elec_num);
	assert(rc == TREXIO_SUCCESS);

	double *elec_coord = malloc(sizeof(double) * walk_num * elec_num * 3);
	double *elec_coord_device =
		qmckl_malloc_device(context, sizeof(double) * walk_num * elec_num * 3);
	assert(elec_coord != NULL);
	rc = trexio_read_qmc_point(trexio_file, elec_coord);
	assert(rc == TREXIO_SUCCESS);
	trexio_close(trexio_file);
	qmckl_memcpy_H2D(
		context, elec_coord_device, elec_coord,
		sizeof(double) * walk_num * elec_num * 3
	);

	printf("Reading %s.\n", file_name);
	rc = qmckl_trexio_read_device(context, file_name, 255);
	assert(rc == QMCKL_SUCCESS_DEVICE);

	int64_t nucl_num;
	rc = qmckl_get_nucleus_num_device(context, &nucl_num);
	assert(rc == TREXIO_SUCCESS);

	int64_t ao_num;
	rc = qmckl_get_ao_basis_ao_num_device(context, &ao_num);
	assert(rc == QMCKL_SUCCESS_DEVICE);

	rc = qmckl_set_electron_coord_device(context, 'N', walk_num, elec_coord_device,
								  walk_num * elec_num * 3);
	assert(rc == QMCKL_SUCCESS_DEVICE);

	const int64_t aord_num = 5;
	const int64_t bord_num = 5;
	const int64_t cord_num = 5;
	const int64_t type_nucl_num = 1;

	assert(rc == QMCKL_SUCCESS_DEVICE);

	rc = qmckl_set_jastrow_aord_num_device(context, aord_num);
	assert(rc == QMCKL_SUCCESS_DEVICE);

	rc = qmckl_set_jastrow_bord_num_device(context, bord_num);
	assert(rc == QMCKL_SUCCESS_DEVICE);

	rc = qmckl_set_jastrow_cord_num_device(context, cord_num);
	assert(rc == QMCKL_SUCCESS_DEVICE);

	rc = qmckl_set_jastrow_type_nucl_num_device(context, type_nucl_num);
	assert(rc == QMCKL_SUCCESS_DEVICE);

	double nucl_rescale[type_nucl_num];
	for (int64_t i = 0; i < type_nucl_num; ++i) {
		nucl_rescale[i] = 2.0;
	}
	double * nucl_rescale_device = qmckl_malloc_device(context, type_nucl_num * sizeof(double));
	qmckl_memcpy_H2D(context, nucl_rescale_device, nucl_rescale, type_nucl_num*sizeof(double));
	rc = qmckl_set_jastrow_rescale_factor_en_device(context, nucl_rescale_device,
												   type_nucl_num);
	assert(rc == QMCKL_SUCCESS_DEVICE);

	rc = qmckl_set_jastrow_rescale_factor_ee_device(context, 2.0);
	assert(rc == QMCKL_SUCCESS_DEVICE);

	int64_t *type_nucl_vector = malloc(nucl_num * sizeof(int64_t));
	for (int i = 0; i < nucl_num; ++i)
		type_nucl_vector[i] = 1;
	int64_t *type_nucl_vector_device = qmckl_malloc_device(context, nucl_num * sizeof(int64_t));
	qmckl_memcpy_H2D(
		context, type_nucl_vector_device, type_nucl_vector,
		nucl_num * sizeof(int64_t)
	);
	rc = qmckl_set_jastrow_type_nucl_vector_device(context, type_nucl_vector,
												  nucl_num);
	assert(rc == QMCKL_SUCCESS_DEVICE);
	free(type_nucl_vector);

	const double *aord_vector = &(jast_coef[0]);
	double * aord_vector_device = qmckl_malloc_device(context, 11080 * sizeof(double));
	qmckl_memcpy_H2D(
		context, aord_vector_device, aord_vector,
		11080 * sizeof(double)
	);
	rc = qmckl_set_jastrow_a_vector_device(context, aord_vector, 11080);
	assert(rc == QMCKL_SUCCESS_DEVICE);

	const double *bord_vector = &(jast_coef[1]);
	double * bord_vector_device = aord_vector_device+1;
	rc = qmckl_set_jastrow_b_vector_device(context, bord_vector_device, 11079);
	assert(rc == QMCKL_SUCCESS_DEVICE);

	const double *cord_vector = &(jast_coef[2]);
	double * cord_vector_device = aord_vector_device+2;
	rc = qmckl_set_jastrow_c_vector_device(context, cord_vector_device, 11078);
	assert(rc == QMCKL_SUCCESS_DEVICE);

	double jast_gl[walk_num][4][elec_num];
	double * jast_gl_device = qmckl_malloc_device(context, walk_num*4*elec_num);

	gettimeofday(&timecheck, NULL);
	start = (long)timecheck.tv_sec * 1000 + (long)timecheck.tv_usec / 1000;
	for (int i = 0; i < ITERMAX; ++i) {
		printf("%3d / %3d\n", i, ITERMAX);
		rc = qmckl_context_touch_device(context);
		rc = qmckl_get_jastrow_value_device(context, jast_gl_device,
										   walk_num);
	}
	gettimeofday(&timecheck, NULL);
	end = (long)timecheck.tv_sec * 1000 + (long)timecheck.tv_usec / 1000;

	printf("Value: Time for the calculation of 1 step (ms): %10.1f\n",
		   (double)(end - start) / (double)ITERMAX);

	gettimeofday(&timecheck, NULL);
	start = (long)timecheck.tv_sec * 1000 + (long)timecheck.tv_usec / 1000;
	for (int i = 0; i < ITERMAX; ++i) {
		printf("%3d / %3d\n", i, ITERMAX);
		rc = qmckl_context_touch_device(context);
		rc = qmckl_get_jastrow_gl_device(context, &(jast_gl[0][0][0]),
										walk_num * elec_num * 4);
	}
	gettimeofday(&timecheck, NULL);
	end = (long)timecheck.tv_sec * 1000 + (long)timecheck.tv_usec / 1000;

	printf("GL: Time for the calculation of 1 step (ms): %10.1f\n",
		   (double)(end - start) / (double)ITERMAX);

	rc = qmckl_context_destroy_device(context);

	free(elec_coord);
}
