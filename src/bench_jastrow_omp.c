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

#include "../data/jast_coef.h"

#define ITERMAX 3
#define DEVICE_ID 0



int main(int argc, char** argv)
{
  long start, end;
  struct timeval timecheck;

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
  assert (ao_vgl_device != NULL);

  rc = qmckl_set_electron_coord_device(context, 'N', walk_num, elec_coord_device, walk_num*elec_num*3, DEVICE_ID);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_ao_basis_ao_vgl_device(context, ao_vgl_device, size_max, DEVICE_ID);
  assert (rc == QMCKL_SUCCESS);


  rc = qmckl_set_electron_rescale_factor_en(context, 2.0);
  assert(rc == QMCKL_SUCCESS);

  rc = qmckl_set_electron_rescale_factor_ee(context, 2.0);
  assert(rc == QMCKL_SUCCESS);

  rc = qmckl_set_nucleus_rescale_factor (context, 2.0);
  assert(rc == QMCKL_SUCCESS);

  int64_t nucl_num;
  const int64_t aord_num = 5;
  const int64_t bord_num = 5;
  const int64_t cord_num = 5;
  const int64_t type_nucl_num = 1;

  rc = qmckl_get_nucleus_num(context, &nucl_num);

  assert(rc == QMCKL_SUCCESS);

  qmckl_init_jastrow(context);

  rc = qmckl_set_jastrow_ord_num(context, aord_num, bord_num, cord_num);
  assert(rc == QMCKL_SUCCESS);

  rc = qmckl_set_jastrow_type_nucl_num(context, type_nucl_num);
  assert(rc == QMCKL_SUCCESS);

  // int64_t* type_nucl_vector = malloc(nucl_num * sizeof(int64_t));
  int64_t* type_nucl_vector = omp_target_alloc(nucl_num * sizeof(int64_t), DEVICE_ID);
  #pragma omp target is_device_ptr(type_nucl_vector)
  {
  #pragma omp teams distribute parallel for
  for (int i=0 ; i<nucl_num; ++i) type_nucl_vector[i] = 1;
  }

  rc = qmckl_set_jastrow_type_nucl_vector_device(context, type_nucl_vector, nucl_num);
  assert(rc == QMCKL_SUCCESS);
  omp_target_free(type_nucl_vector, DEVICE_ID);


  // Put jast_coef on device
  double * jast_coef_device = omp_target_alloc(11080 * sizeof(double), DEVICE_ID);
  omp_target_memcpy(
    jast_coef_device, jast_coef,
    sizeof(double)*11080,
    0, 0,
    DEVICE_ID, omp_get_initial_device()
  );


  const double* aord_vector = jast_coef_device;

  rc = qmckl_set_jastrow_aord_vector_device(context, aord_vector, 11080, DEVICE_ID);
  assert(rc == QMCKL_SUCCESS);

  const double* bord_vector = jast_coef_device+1;
  rc = qmckl_set_jastrow_bord_vector_device(context, bord_vector, 11079, DEVICE_ID);
  assert(rc == QMCKL_SUCCESS);

  const double* cord_vector = jast_coef_device+2;
  rc = qmckl_set_jastrow_cord_vector_device(context, cord_vector, 11078, DEVICE_ID);
  assert(rc == QMCKL_SUCCESS);


  double * factor_een_deriv_e_device = omp_target_alloc(walk_num*4*elec_num*sizeof(double), DEVICE_ID);


  printf("[bench] 1\n");

  gettimeofday(&timecheck, NULL);
  start = (long)timecheck.tv_sec * 1000 + (long)timecheck.tv_usec / 1000;
  for (int i=0 ; i<ITERMAX ; ++i) {
    rc = qmckl_context_touch(context);
    rc = qmckl_get_jastrow_factor_een_deriv_e_device(context, factor_een_deriv_e_device, walk_num*elec_num*4, DEVICE_ID);
  }
  gettimeofday(&timecheck, NULL);
  end = (long)timecheck.tv_sec * 1000 + (long)timecheck.tv_usec / 1000;

  printf("[bench] 2\n");
  printf("Time for the calculation of 1 step (ms): %10.1f\n", (double) (end-start) / (double) ITERMAX);

}
