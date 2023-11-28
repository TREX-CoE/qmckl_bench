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

#include <time.h>

#define ITERMAX_AO 100
#define ITERMAX_MO 100

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
  assert (elec_coord != NULL);
  rc = trexio_read_qmc_point(trexio_file, elec_coord);
  assert (rc == TREXIO_SUCCESS);
  trexio_close(trexio_file);

  printf("Reading %s.\n", file_name);
  rc = qmckl_trexio_read(context, file_name, 255);
  if (rc != QMCKL_SUCCESS) {
    printf("%s\n", qmckl_string_of_error(rc));
  }
  assert (rc == QMCKL_SUCCESS);

  int64_t ao_num;
  rc = qmckl_get_ao_basis_ao_num(context, &ao_num);
  assert (rc == QMCKL_SUCCESS);

  const int64_t size_max_ao = 5*walk_num*elec_num*ao_num;
  double * ao_vgl = malloc (size_max_ao * sizeof(double));
  assert (ao_vgl != NULL);

  int64_t mo_num;
  rc = qmckl_get_mo_basis_mo_num(context, &mo_num);
  assert (rc == QMCKL_SUCCESS);

  const int64_t size_max_mo = 5*walk_num*elec_num*mo_num;
  double * mo_vgl = malloc (size_max_mo * sizeof(double));
  assert (mo_vgl != NULL);

//  { int precision = 20;
  for (int precision=2 ; precision < 54 ; precision++) {

    qmckl_set_numprec_precision(context, precision);

    rc = qmckl_set_electron_coord(context, 'N', walk_num, elec_coord, walk_num*elec_num*3);
    assert (rc == QMCKL_SUCCESS);

    rc = qmckl_get_ao_basis_ao_vgl(context, ao_vgl, size_max_ao);
    assert (rc == QMCKL_SUCCESS);

    rc = qmckl_get_mo_basis_mo_vgl(context, mo_vgl, size_max_mo);
    assert (rc == QMCKL_SUCCESS);

    gettimeofday(&timecheck, NULL);
    start = (long)timecheck.tv_sec * 1000 + (long)timecheck.tv_usec / 1000;

    for (int i=0 ; i<ITERMAX_AO ; ++i) {
      rc = qmckl_get_ao_basis_ao_vgl_inplace(context, ao_vgl, size_max_ao);
    }
    gettimeofday(&timecheck, NULL);
    end = (long)timecheck.tv_sec * 1000 + (long)timecheck.tv_usec / 1000;

    printf("%2d : %10.1f  ", precision, (double) (end-start) / (double) ITERMAX_AO);

    gettimeofday(&timecheck, NULL);
    start = (long)timecheck.tv_sec * 1000 + (long)timecheck.tv_usec / 1000;

    rc = qmckl_get_ao_basis_ao_value(context, ao_vgl, size_max_ao);
    assert (rc == QMCKL_SUCCESS);

    for (int i=0 ; i<ITERMAX_AO ; ++i) {
      rc = qmckl_get_ao_basis_ao_value_inplace(context, ao_vgl, size_max_ao);
    }
    gettimeofday(&timecheck, NULL);
    end = (long)timecheck.tv_sec * 1000 + (long)timecheck.tv_usec / 1000;

    printf("%10.1f  ", (double) (end-start) / (double) ITERMAX_AO);

    for (int i=0 ; i<ITERMAX_MO ; ++i) {
      rc = qmckl_get_mo_basis_mo_vgl_inplace(context, mo_vgl, size_max_mo);
    }
    gettimeofday(&timecheck, NULL);
    end = (long)timecheck.tv_sec * 1000 + (long)timecheck.tv_usec / 1000;

    printf("%10.1f  ", (double) (end-start) / (double) ITERMAX_MO);

    gettimeofday(&timecheck, NULL);
    start = (long)timecheck.tv_sec * 1000 + (long)timecheck.tv_usec / 1000;

    rc = qmckl_get_mo_basis_mo_value(context, mo_vgl, size_max_mo);
    assert (rc == QMCKL_SUCCESS);

    for (int i=0 ; i<ITERMAX_MO ; ++i) {
      rc = qmckl_get_mo_basis_mo_value_inplace(context, mo_vgl, size_max_mo);
    }
    gettimeofday(&timecheck, NULL);
    end = (long)timecheck.tv_sec * 1000 + (long)timecheck.tv_usec / 1000;

    printf("%10.1f  ", (double) (end-start) / (double) ITERMAX_MO);
    printf("\n");
  }

  rc = qmckl_context_destroy(context);
  free(elec_coord);
  free(ao_vgl);

}
